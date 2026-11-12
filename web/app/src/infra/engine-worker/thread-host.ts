// ThreadHost (plan §3.1, §4.7): the `wasi.thread-spawn` of the threaded module in the
// Engine worker. The engine starts N-1 threads for a search of N threads, on the thread
// that calls `losat_web2_run` (LOSAT/src/utils/threading.rs), and that call is
// synchronous, so a thread must start without the Engine worker returning to its event
// loop. Each thread therefore runs in a thread worker that was prepared beforehand: it has
// instantiated the module on the shared memory and waits for a start message. A thread
// worker is reused: when its thread ends, it instantiates the module again and becomes
// ready for the next thread. The Node host of the same ABI (LOSAT/tests/wasi_thread_host.js)
// can also start a new Node worker for a thread while the spawning thread waits
// (LOSAT_WASI_PRESPAWN=0), though by default it prepares standby workers too; a browser
// cannot be relied on to start a worker then.
//
// A run of the engine builds one thread pool for all of its searches (`with_run_pool`,
// LOSAT/src/utils/threading.rs: a run searches once per hit list size of its formats, and
// BLASTP's block search can search again), so its N-1 threads start once. A pool after another
// would need thread workers of its own: the threads of a pool may not have returned when the next
// pool spawns, and they cannot return while it waits, because the engine's thread library holds a
// lock across `thread-spawn` that an ending thread needs (S09). A search of N threads therefore
// has one set of N-1 thread workers, all ready before the search starts (`prepare`).
import { growToMaximum } from '../reactor/shared-memory';
import type { InstanceUse } from './instance-use';
import {
  SLOT_ENDED,
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  SLOT_RUNNING,
  START_ABANDONED,
  START_PENDING,
  START_TAKEN,
  type KeeperCommand,
  type ThreadCommand,
  type ThreadFault,
} from './protocol';

export interface ThreadHostOptions {
  readonly module: WebAssembly.Module;
  readonly memory: WebAssembly.Memory;
  /** The maximum of `memory` in pages (shared-memory.ts). */
  readonly maximumPages: number;
  /** BroadcastChannel name on which a thread that stops reports it (to the main thread). */
  readonly faultChannel: string;
  /** BroadcastChannel name on which the main thread asks the keeper to grow the memory (memory-keeper.ts). */
  readonly keeperChannel: string;
  /**
   * The use count of the instance (instance-use.ts): each thread counts from its spawn until
   * it returns, so that the Gateway does not end the runtime at once while a thread may grow
   * the memory.
   */
  readonly use: InstanceUse;
  readonly createWorker?: () => Worker;
  /** Starts the keeper of the memory (tests may start a fake one). */
  readonly createKeeper?: () => Worker;
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
/** The longest wait of this worker's thread at a time (spawn), as the module's own waits. */
const WAIT_SLICE_MS = 100;

export class ThreadHost {
  private readonly slots: Slot[] = [];
  private nextTid = 1;
  private faults: BroadcastChannel | undefined;
  /** The keeper of the shared memory (memory-keeper.ts): it ends with the thread workers. */
  private readonly keeper: Worker;
  private keeperReady = false;
  private keeperError: string | undefined;

  constructor(private readonly options: ThreadHostOptions) {
    this.keeper =
      options.createKeeper?.() ??
      new Worker(new URL('./memory-keeper-worker.ts', import.meta.url), { type: 'module', name: 'losat-memory' });
    this.keeper.onmessage = (event: MessageEvent<{ type?: string }>) => {
      if (event.data?.type === 'ready') this.keeperReady = true;
    };
    this.keeper.onerror = (event) => {
      // A script that does not load fires a plain Event, without a message.
      this.keeperError = `the memory keeper failed: ${event.message || 'its script did not load'}`;
    };
    const keep: KeeperCommand = {
      type: 'keep',
      memory: options.memory,
      maximumPages: options.maximumPages,
      channel: options.keeperChannel,
    };
    this.keeper.postMessage(keep);
  }

  /**
   * Makes `count` thread workers (the threads of a search of `count` + 1 threads), if there
   * are fewer, and waits until `count` of them are ready: also those whose thread of the
   * previous search is still returning. The keeper must hold the memory too: without it, the
   * Gateway could not have the memory grown before it ends the runtime during the search.
   */
  async prepare(count: number): Promise<void> {
    // A runtime that the Gateway ends starts no thread worker (runtime.ts, load).
    if (this.options.use.retiring()) throw new ThreadHostError('the runtime is ending');
    const total = count;
    while (this.slots.length < total) this.slots.push(this.createSlot());
    const deadline = performance.now() + (this.options.readyTimeoutMs ?? DEFAULT_TIMEOUT_MS);
    for (;;) {
      const failed = this.slots.find((slot) => Atomics.load(slot.state, 0) === SLOT_FAILED);
      if (failed !== undefined) throw new ThreadHostError(failed.error ?? 'a thread worker stopped');
      if (this.keeperError !== undefined) throw new ThreadHostError(this.keeperError);
      const ready = this.slots.filter((slot) => Atomics.load(slot.state, 0) === SLOT_READY).length >= total;
      if (ready && this.keeperReady) return;
      if (performance.now() > deadline) {
        throw new ThreadHostError(ready ? 'the memory keeper did not start in time' : 'the thread workers did not become ready in time');
      }
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
    // Counted before it can start; the thread worker counts it off when the thread returns.
    this.options.use.startThread();
    const tid = this.nextTid++;
    const started = new Int32Array(new SharedArrayBuffer(4));
    const command: ThreadCommand = { type: 'start', tid, startArg, started: started.buffer as SharedArrayBuffer };
    slot.worker.postMessage(command);
    // In slices: WebKit never ends a worker whose thread waits without returning, and keeps the
    // shared memory with it (LOSAT/tests/wasi_shared_memory.js, the module's waits).
    const deadline = performance.now() + (this.options.startTimeoutMs ?? DEFAULT_TIMEOUT_MS);
    for (let left = deadline - performance.now(); left > 0; left = deadline - performance.now()) {
      if (Atomics.wait(started, 0, START_PENDING, Math.min(left, WAIT_SLICE_MS)) !== 'timed-out') break;
    }
    // The thread worker takes the start (pending -> taken) unless the host gives up first
    // (pending -> abandoned): a thread must never start after its spawn returned -1, because
    // the engine frees its start argument then.
    const outcome = Atomics.compareExchange(started, 0, START_PENDING, START_ABANDONED);
    if (outcome === START_TAKEN) return tid;
    slot.error = `thread ${tid} did not start`;
    this.options.use.endThread(slot.state, SLOT_FAILED);
    // The search cannot go on as the engine planned it; the main thread ends this runtime.
    this.fault({ type: 'fault', tid, message: `The engine stopped: thread ${tid} did not start in time` });
    return -1;
  };

  /**
   * Ends every thread worker and the keeper, and returns true; false, ending nothing, once the
   * Gateway ends the runtime, which may be waiting for the keeper's answer (instance-use.ts).
   * A thread that still runs may be growing the memory, and WebKit crashes the page when
   * workers with instances on it end while a thread grows it (growToMaximum). The memory is
   * then first grown to its maximum, where no thread can grow it, as the Gateway has the
   * keeper do before it ends a runtime whose threaded instance is in use. A slot runs from its
   * spawn until its thread returns, and a thread of a pool may not have returned when the
   * search did (see above), for example when the next search of the runtime runs on the other
   * module (runtime.ts, load). Such a thread cannot count itself off once its worker ends, so
   * it is counted off here, after the grow. Once a thread has started here, the ending is
   * itself a use, until every thread is counted off: a Gateway that ends the runtime meanwhile
   * asks the keeper rather than ending it at once, and then the keeper stays. A thread starts
   * only after `prepare`, so the keeper listens then; before, no thread is counted.
   */
  terminate(): boolean {
    const { use } = this.options;
    if (this.nextTid === 1) {
      if (use.retiring()) return false;
    } else {
      if (!use.enter()) return false;
      if (this.slots.some((slot) => Atomics.load(slot.state, 0) === SLOT_RUNNING)) {
        try {
          growToMaximum({ memory: this.options.memory, maximumPages: this.options.maximumPages });
        } catch {
          // The workers end all the same; the memory stays as it is (memory-keeper.ts).
        }
      }
      for (const slot of this.slots) use.endThread(slot.state, SLOT_ENDED);
      if (!use.leaveLast()) return false;
    }
    for (const slot of this.slots.splice(0)) slot.worker.terminate();
    this.keeper.terminate();
    this.faults?.close();
    this.faults = undefined;
    return true;
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
      // An error that ended a thread that ran counts it off.
      if (!this.options.use.endThread(state, SLOT_FAILED)) Atomics.store(state, 0, SLOT_FAILED);
    };
    const command: ThreadCommand = {
      type: 'prepare',
      module: this.options.module,
      memory: this.options.memory,
      slot: state.buffer as SharedArrayBuffer,
      faultChannel: this.options.faultChannel,
      use: this.options.use.buffer,
    };
    worker.postMessage(command);
    return slot;
  }
}
