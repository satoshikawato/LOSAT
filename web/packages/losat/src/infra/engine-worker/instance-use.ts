// Whether the threaded instance of a runtime generation is in use: one Int32 on a
// SharedArrayBuffer that the Gateway (main thread), the Engine worker and its thread workers
// read and change synchronously (Atomics). It counts the calls of the Engine worker into the
// threaded instance, the threads of the search pool that run, and the grow with which the
// ThreadHost ends thread workers whose threads run (ThreadHost.terminate); a bit marks that
// the Gateway ends the generation (retire). Not counted: the creation of a memory and the
// instantiation of its main instance (runtime.ts, load), which come before any thread worker
// holds an instance on that memory, and before its keeper listens.
//
// While the instance is in use, a thread may grow its shared memory, and WebKit crashes the
// page when workers with instances on that memory end while a thread grows it (gateway.ts,
// retire; growToMaximum). The Gateway therefore ends the Engine worker at once only when the
// instance is not in use, and otherwise after the keeper has grown the memory to its maximum.
// A message from the Engine worker ("a threaded search starts") could not decide it: it is
// in flight for a moment after the call has started, and a warm-up of the threaded instance
// sent none. Here the Gateway's mark and the Engine worker's entry are one atomic operation
// each, so one of them sees the other: either the Gateway sees the instance in use, or the
// Engine worker sees the mark and does not call into the instance.
//
// The buffer holds no part of the engine's shared memory: the main thread must not hold that
// memory (memory-keeper.ts).
import { SLOT_RUNNING } from './protocol';

/** The bit that the Gateway sets when it ends the generation; the count of uses is below it. */
export const USE_RETIRING = 0x4000_0000;
const USE_COUNT = USE_RETIRING - 1;
/** Nothing notifies it: `Atomics.wait` on it pauses for its timeout (leaveLast). */
const PAUSE = typeof SharedArrayBuffer === 'function' ? new Int32Array(new SharedArrayBuffer(4)) : new Int32Array(1);

export class InstanceUse {
  private readonly word: Int32Array;

  /** `buffer`: a SharedArrayBuffer of at least 4 bytes, all sides of one generation the same. */
  constructor(readonly buffer: SharedArrayBuffer) {
    this.word = new Int32Array(buffer, 0, 1);
  }

  /** A new buffer for a generation: zero uses, not retiring. */
  static create(): InstanceUse {
    return new InstanceUse(new SharedArrayBuffer(4));
  }

  /**
   * The Engine worker, before it calls into the threaded instance: counts the call, unless the
   * Gateway ends the generation (false: the caller must not call into the instance).
   */
  enter(): boolean {
    for (;;) {
      const value = Atomics.load(this.word, 0);
      if ((value & USE_RETIRING) !== 0) return false;
      if (Atomics.compareExchange(this.word, 0, value, value + 1) === value) return true;
    }
  }

  /** The end of a call counted by `enter`. */
  leave(): void {
    Atomics.sub(this.word, 0, 1);
  }

  /**
   * The end of a use counted by `enter` that must be the last one (ThreadHost.terminate, once
   * it has counted off every thread that ran): a thread worker counts its thread off in two
   * steps (`endThread`), so this waits until the count holds this use alone, then leaves.
   * False, the use staying counted, as soon as the Gateway has marked the generation retiring:
   * it then saw a use and waits for the keeper's answer, so the caller ends nothing of what the
   * keeper needs.
   */
  leaveLast(): boolean {
    for (;;) {
      const value = Atomics.load(this.word, 0);
      if ((value & USE_RETIRING) !== 0) return false;
      const count = value & USE_COUNT;
      if (count === 0) throw new Error('the use count of the threaded instance lost a use');
      if (count === 1) {
        if (Atomics.compareExchange(this.word, 0, value, 0) === value) return true;
      } else {
        Atomics.wait(PAUSE, 0, 0, 1);
      }
    }
  }

  /**
   * A thread of the search pool that the ThreadHost starts in a slot it has set to
   * SLOT_RUNNING. It starts during a call of the Engine worker, which is counted, so the
   * Gateway already treats the instance as in use: it is counted without the check of `enter`.
   */
  startThread(): void {
    Atomics.add(this.word, 0, 1);
  }

  /**
   * A thread that ends: its slot (thread-host.ts) goes from SLOT_RUNNING to `next`, and then
   * the thread is counted off. The thread worker calls it when its thread returns or stops, and
   * the ThreadHost when the thread did not start or when it ends the worker of a thread that
   * runs: only the first of them changes the slot and counts the thread off (true).
   */
  endThread(slot: Int32Array, next: number): boolean {
    if (Atomics.compareExchange(slot, 0, SLOT_RUNNING, next) !== SLOT_RUNNING) return false;
    Atomics.sub(this.word, 0, 1);
    return true;
  }

  /** Whether the Gateway ends the generation (the Engine worker then leaves its workers as they are). */
  retiring(): boolean {
    return (Atomics.load(this.word, 0) & USE_RETIRING) !== 0;
  }

  /**
   * The Gateway, when it ends the generation: marks it retiring and returns whether the
   * instance was in use at that moment. If it was not, no call or thread can start after it.
   */
  retire(): boolean {
    return (Atomics.or(this.word, 0, USE_RETIRING) & USE_COUNT) !== 0;
  }

  /** The number of counted calls and threads (tests). */
  uses(): number {
    return Atomics.load(this.word, 0) & USE_COUNT;
  }
}
