// ThreadHost (plan §3.1, §4.7): the `wasi.thread-spawn` of the threaded module in the
// Engine worker. The engine starts N-1 threads for a search of N threads, on the thread
// that calls `losat_web2_run` (LOSAT/src/utils/threading.rs), and that call is
// synchronous, so a thread must start without the Engine worker returning to its event
// loop. Each thread therefore runs in a thread worker that was prepared beforehand: it has
// instantiated the module on the shared memory and waits for a start message. A thread
// worker is reused: when its thread ends, it instantiates the module again and becomes
// ready for the next thread. The Node host of the same ABI (LOSAT/tests/wasi_thread_host.js)
// starts a new Node worker for each thread instead, which Node can do while the spawning
// thread waits; a browser cannot be relied on to start a worker then.
//
// A search can build one thread pool after another (TBLASTN builds one for each batch of
// queries). The threads of a pool return only after the spawns of the next pool have
// returned: the engine's thread library holds a lock across `thread-spawn` that an ending
// thread needs (S09). So a search of N threads has two sets of N-1 thread workers: the next
// pool starts in the set that the previous pool did not use.
import {
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  SLOT_RUNNING,
  START_ABANDONED,
  START_PENDING,
  START_TAKEN,
  type ThreadCommand,
  type ThreadFault,
} from './protocol';

export interface ThreadHostOptions {
  readonly module: WebAssembly.Module;
  readonly memory: WebAssembly.Memory;
  /** BroadcastChannel name on which a thread that stops reports it (to the main thread). */
  readonly faultChannel: string;
  readonly createWorker?: () => Worker;
  /** How long a thread worker may take to become ready; 30 s by default. */
  readonly readyTimeoutMs?: number;
  /** How long a prepared thread worker may take to start a thread; 30 s by default. */
  readonly startTimeoutMs?: number;
}

export class ThreadHostError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ThreadHostError';
  }
}

interface Slot {
  readonly worker: Worker;
  readonly state: Int32Array;
  error: string | undefined;
}

const DEFAULT_TIMEOUT_MS = 30_000;

export class ThreadHost {
  private readonly slots: Slot[] = [];
  private nextTid = 1;
  private faults: BroadcastChannel | undefined;

  constructor(private readonly options: ThreadHostOptions) {}

  /**
   * Makes two sets of `count` thread workers (the threads of a search of `count` + 1
   * threads, and a spare set for its next pool) and waits until all of them are ready.
   */
  async prepare(count: number): Promise<void> {
    const total = 2 * count;
    while (this.slots.length < total) this.slots.push(this.createSlot());
    const deadline = performance.now() + (this.options.readyTimeoutMs ?? DEFAULT_TIMEOUT_MS);
    for (;;) {
      const failed = this.slots.find((slot) => Atomics.load(slot.state, 0) === SLOT_FAILED);
      if (failed !== undefined) throw new ThreadHostError(failed.error ?? 'a thread worker stopped');
      if (this.slots.filter((slot) => Atomics.load(slot.state, 0) === SLOT_READY).length >= total) return;
      if (performance.now() > deadline) throw new ThreadHostError('the thread workers did not become ready in time');
      await new Promise((resolve) => setTimeout(resolve, 2));
    }
  }

  /**
   * The `wasi.thread-spawn` import: starts a thread in a ready thread worker and returns
   * its thread ID, or -1 if none can start it (the engine then fails the search).
   */
  readonly spawn = (startArg: number): number => {
    const slot = this.claim();
    if (slot === undefined) return -1;
    const tid = this.nextTid++;
    const started = new Int32Array(new SharedArrayBuffer(4));
    const command: ThreadCommand = { type: 'start', tid, startArg, started: started.buffer as SharedArrayBuffer };
    slot.worker.postMessage(command);
    Atomics.wait(started, 0, START_PENDING, this.options.startTimeoutMs ?? DEFAULT_TIMEOUT_MS);
    // The thread worker takes the start (pending -> taken) unless the host gives up first
    // (pending -> abandoned): a thread must never start after its spawn returned -1, because
    // the engine frees its start argument then.
    const outcome = Atomics.compareExchange(started, 0, START_PENDING, START_ABANDONED);
    if (outcome === START_TAKEN) return tid;
    slot.error = `thread ${tid} did not start`;
    Atomics.store(slot.state, 0, SLOT_FAILED);
    // The search cannot go on as the engine planned it; the main thread ends this runtime.
    this.fault({ type: 'fault', tid, message: `The engine stopped: thread ${tid} did not start in time` });
    return -1;
  };

  /** Ends every thread worker. */
  terminate(): void {
    for (const slot of this.slots.splice(0)) slot.worker.terminate();
    this.faults?.close();
    this.faults = undefined;
  }

  private fault(message: ThreadFault): void {
    this.faults ??= new BroadcastChannel(this.options.faultChannel);
    this.faults.postMessage(message);
  }

  /**
   * Takes a ready slot, waiting for one that is being prepared again (its thread has
   * returned). A slot whose thread is still running is not waited for: it may be a thread of
   * the previous pool, which cannot return while this spawn waits.
   */
  private claim(): Slot | undefined {
    const deadline = performance.now() + (this.options.readyTimeoutMs ?? DEFAULT_TIMEOUT_MS);
    for (;;) {
      for (const slot of this.slots) {
        if (Atomics.compareExchange(slot.state, 0, SLOT_READY, SLOT_RUNNING) === SLOT_READY) return slot;
      }
      const preparing = this.slots.find((slot) => Atomics.load(slot.state, 0) === SLOT_PREPARING);
      if (preparing === undefined || performance.now() > deadline) return undefined;
      Atomics.wait(preparing.state, 0, SLOT_PREPARING, 50);
    }
  }

  private createSlot(): Slot {
    const worker =
      this.options.createWorker?.() ??
      new Worker(new URL('./thread-worker.ts', import.meta.url), { type: 'module', name: 'losat-thread' });
    const state = new Int32Array(new SharedArrayBuffer(4));
    Atomics.store(state, 0, SLOT_PREPARING);
    const slot: Slot = { worker, state, error: undefined };
    worker.onmessage = (event: MessageEvent<{ type?: string; message?: string }>) => {
      if (event.data?.type === 'failed') slot.error = event.data.message;
    };
    worker.onerror = (event) => {
      slot.error = `a thread worker failed: ${event.message}`;
      Atomics.store(state, 0, SLOT_FAILED);
    };
    const command: ThreadCommand = {
      type: 'prepare',
      module: this.options.module,
      memory: this.options.memory,
      slot: state.buffer as SharedArrayBuffer,
      faultChannel: this.options.faultChannel,
    };
    worker.postMessage(command);
    return slot;
  }
}
