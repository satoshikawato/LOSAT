// ThreadHost (src/infra/engine-worker/thread-host.ts) with fake thread workers: preparing
// them, starting a thread, a thread worker that does not start in time, the keeper of the
// shared memory, and the use count of the threads. The real thread workers run in the
// browsers (tests/e2e/engine.spec.ts).
import { Worker as NodeWorker } from 'node:worker_threads';
import { describe, expect, it, vi } from 'vitest';
import { InstanceUse } from '../../src/infra/engine-worker/instance-use';
import {
  SLOT_ENDED,
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  START_ABANDONED,
  START_PENDING,
  START_TAKEN,
  type KeeperCommand,
  type ThreadCommand,
  type ThreadFault,
} from '../../src/infra/engine-worker/protocol';
import { ThreadHost, ThreadHostError } from '../../src/infra/engine-worker/thread-host';

type Behaviour = 'ready-and-start' | 'ready-never-start' | 'never-ready';

class FakeThreadWorker {
  state: Int32Array | undefined;
  use: InstanceUse | undefined;
  pendingStart: Int32Array | undefined;
  starts: number[] = [];
  terminated = false;
  onmessage: ((event: MessageEvent) => void) | null = null;
  onerror: ((event: ErrorEvent) => void) | null = null;

  constructor(private readonly behaviour: Behaviour) {}

  postMessage(command: ThreadCommand): void {
    if (command.type === 'prepare') {
      this.state = new Int32Array(command.slot);
      this.use = new InstanceUse(command.use);
      if (this.behaviour !== 'never-ready') Atomics.store(this.state, 0, SLOT_READY);
      return;
    }
    const started = new Int32Array(command.started);
    if (this.behaviour === 'ready-and-start') {
      // As thread-worker.ts takes a start.
      if (Atomics.compareExchange(started, 0, START_PENDING, START_TAKEN) === START_PENDING) this.starts.push(command.tid);
      Atomics.notify(started, 0);
    } else {
      this.pendingStart = started;
    }
  }

  /** Its thread returns, as in thread-worker.ts: whether it counted the thread off. */
  returns(): boolean {
    return this.use!.endThread(this.state!, SLOT_PREPARING);
  }

  /** The late start of a worker whose host stopped waiting. */
  startLate(): boolean {
    return Atomics.compareExchange(this.pendingStart!, 0, START_PENDING, START_TAKEN) === START_PENDING;
  }

  terminate(): void {
    this.terminated = true;
  }
}

type KeeperBehaviour = 'ready' | 'never-ready' | 'fails';

class FakeKeeper {
  readonly received: KeeperCommand[] = [];
  terminated = false;
  onmessage: ((event: MessageEvent) => void) | null = null;
  onerror: ((event: ErrorEvent) => void) | null = null;

  constructor(private readonly behaviour: KeeperBehaviour) {}

  postMessage(command: KeeperCommand): void {
    this.received.push(command);
    // As memory-keeper-worker.ts answers, after the message loop turns.
    if (this.behaviour === 'ready') setTimeout(() => this.onmessage?.({ data: { type: 'ready' } } as MessageEvent), 1);
    // A script that does not load fires a plain Event, without a message.
    if (this.behaviour === 'fails') setTimeout(() => this.onerror?.({} as ErrorEvent), 1);
  }

  terminate(): void {
    this.terminated = true;
  }
}

/** A shared memory that records each grow and which workers had ended by then. */
class FakeMemory {
  pages = 2;
  /** A memory whose size cannot be read: growToMaximum throws. */
  unreadable = false;
  /** Called at each grow, before it takes effect (another thread acting meanwhile). */
  onGrow: (() => void) | undefined;
  readonly grows: Array<{ pages: number; threadWorkersEnded: boolean[]; keeperEnded: boolean[] }> = [];

  constructor(private readonly ended: () => { threadWorkers: boolean[]; keepers: boolean[] }) {}

  get buffer(): { byteLength: number } {
    if (this.unreadable) throw new Error('not readable');
    return { byteLength: this.pages * 65_536 };
  }

  grow(delta: number): number {
    this.onGrow?.();
    const { threadWorkers, keepers } = this.ended();
    this.grows.push({ pages: delta, threadWorkersEnded: threadWorkers, keeperEnded: keepers });
    const old = this.pages;
    this.pages += delta;
    return old;
  }
}

/** `readyTimeoutMs` bounds only the cases that fail; a test whose worker must become ready in time raises it. */
function host(behaviour: Behaviour, readyTimeoutMs = 100, keeperBehaviour: KeeperBehaviour = 'ready', startTimeoutMs = 50) {
  const workers: FakeThreadWorker[] = [];
  const keepers: FakeKeeper[] = [];
  const faultChannel = `test-faults-${Math.random()}`;
  const keeperChannel = `test-keeper-${Math.random()}`;
  const use = InstanceUse.create();
  const memory = new FakeMemory(() => ({
    threadWorkers: workers.map((worker) => worker.terminated),
    keepers: keepers.map((keeper) => keeper.terminated),
  }));
  /** The calls of endWaits: the grows before it, and whether the thread workers had ended. */
  const endedWaits: Array<{ grows: number; threadWorkersEnded: boolean[] }> = [];
  const threads = new ThreadHost({
    module: {} as WebAssembly.Module,
    memory: memory as unknown as WebAssembly.Memory,
    maximumPages: 65_535,
    faultChannel,
    keeperChannel,
    use,
    readyTimeoutMs,
    startTimeoutMs,
    createWorker: () => {
      const worker = new FakeThreadWorker(behaviour);
      workers.push(worker);
      return worker as unknown as Worker;
    },
    createKeeper: () => {
      const keeper = new FakeKeeper(keeperBehaviour);
      keepers.push(keeper);
      return keeper as unknown as Worker;
    },
    endWaits: () => {
      endedWaits.push({ grows: memory.grows.length, threadWorkersEnded: workers.map((worker) => worker.terminated) });
    },
  });
  return { threads, workers, keepers, memory, faultChannel, keeperChannel, use, endedWaits };
}

describe('ThreadHost', () => {
  it('prepares one thread worker for each thread and starts threads in them with increasing thread IDs', async () => {
    const { threads, workers } = host('ready-and-start');
    await threads.prepare(2);
    expect(workers).toHaveLength(2);
    expect(threads.spawn(1000)).toBe(1);
    expect(threads.spawn(2000)).toBe(2);
    expect(workers.map((worker) => worker.starts)).toEqual([[1], [2]]);
    threads.terminate();
    expect(workers.every((worker) => worker.terminated)).toBe(true);
  });

  it('fails a spawn at once when every thread worker runs a thread', async () => {
    const { threads } = host('ready-and-start');
    await threads.prepare(1);
    expect(threads.spawn(1000)).toBe(1);
    // The thread has not returned (its slot is running): the threads a spawn would wait for
    // may need this spawn to return before they can, so it does not wait.
    const failing = performance.now();
    expect(threads.spawn(2000)).toBe(-1);
    expect(performance.now() - failing).toBeLessThan(50);
    threads.terminate();
  });

  it('waits for a thread worker that is being prepared again after its thread returned', async () => {
    // A Node worker can take more than 100 ms to start while the other test files run.
    const { threads, workers } = host('ready-and-start', 5_000);
    await threads.prepare(1);
    expect(threads.spawn(1000)).toBe(1);
    // The thread returned; its worker becomes ready 30 ms later, while spawn waits.
    expect(workers[0]!.returns()).toBe(true);
    const state = workers[0]!.state!;
    const ready = new NodeWorker(
      `const { workerData } = require('node:worker_threads');
       const state = new Int32Array(workerData);
       setTimeout(() => { Atomics.store(state, 0, ${SLOT_READY}); Atomics.notify(state, 0); }, 30);`,
      { eval: true, workerData: state.buffer },
    );
    expect(threads.spawn(3000)).toBe(2);
    expect(workers[0]!.starts).toEqual([1, 2]);
    await ready.terminate();
    threads.terminate();
  });

  it('gives the shared memory and its channel to a keeper, waits for it, and ends it with the thread workers', async () => {
    const { threads, workers, keepers, memory, keeperChannel } = host('ready-and-start');
    expect(keepers).toHaveLength(1);
    const keeper = keepers[0]!;
    expect(keeper.received).toEqual([{ type: 'keep', memory, maximumPages: 65_535, channel: keeperChannel }]);
    await threads.prepare(1);
    threads.terminate();
    expect(keeper.terminated).toBe(true);
    expect(workers.every((worker) => worker.terminated)).toBe(true);
  });

  it('grows the memory to its maximum and ends its waits before it ends the workers while a thread runs, and not when none runs', async () => {
    // A thread that runs may be growing the memory: WebKit crashes the page when workers
    // with instances on it end then (a thread of a pool that has not returned when the
    // runtime moves to the other module, runtime.ts load).
    const running = host('ready-and-start');
    await running.threads.prepare(2);
    expect(running.threads.spawn(1000)).toBe(1);
    running.threads.terminate();
    expect(running.memory.grows).toEqual([{ pages: 65_533, threadWorkersEnded: [false, false], keeperEnded: [false] }]);
    // WebKit ends a worker whose thread waits only once the waits are ended (shared-memory.ts, endWaits).
    expect(running.endedWaits).toEqual([{ grows: 1, threadWorkersEnded: [false, false] }]);
    expect(running.workers.every((worker) => worker.terminated)).toBe(true);
    expect(running.keepers[0]!.terminated).toBe(true);

    const idle = host('ready-and-start');
    await idle.threads.prepare(2);
    idle.threads.terminate();
    expect(idle.memory.grows).toEqual([]);
    expect(idle.endedWaits).toEqual([]);
    expect(idle.workers.every((worker) => worker.terminated)).toBe(true);

    // A grow that fails does not keep the workers from ending.
    const broken = host('ready-and-start');
    await broken.threads.prepare(1);
    expect(broken.threads.spawn(1000)).toBe(1);
    broken.memory.unreadable = true;
    expect(() => broken.threads.terminate()).not.toThrow();
    expect(broken.workers.every((worker) => worker.terminated)).toBe(true);
    expect(broken.keepers[0]!.terminated).toBe(true);
  });

  it('counts each thread from its spawn until it returns, once, also when the workers end first', async () => {
    const { threads, workers, use } = host('ready-and-start');
    await threads.prepare(3);
    // The spawns come from a call of the Engine worker into the instance, which it counts.
    expect(use.enter()).toBe(true);
    expect([threads.spawn(1000), threads.spawn(2000), threads.spawn(3000)]).toEqual([1, 2, 3]);
    expect(use.uses()).toBe(4);
    use.leave();
    // A pool thread may return after the search did: until then the instance stays in use.
    expect(workers[0]!.returns()).toBe(true);
    expect(workers[0]!.returns()).toBe(false);
    expect(use.uses()).toBe(2);
    // The host ends the workers of the threads that still run: they cannot count themselves off.
    threads.terminate();
    expect(use.uses()).toBe(0);
    expect(workers.map((worker) => Atomics.load(worker.state!, 0))).toEqual([SLOT_PREPARING, SLOT_ENDED, SLOT_ENDED]);
    // A thread that returns while its worker ends is not counted off a second time.
    expect(workers[1]!.returns()).toBe(false);
    expect(use.uses()).toBe(0);
  });

  it('its grow before it ends the workers of running threads is a use: a Gateway that ends the runtime then keeps every worker', async () => {
    // The Gateway ends the runtime while the Engine worker grows the memory to end a thread
    // that still runs: it must see the instance in use, and then the keeper must stay to answer.
    const growing = host('ready-and-start');
    await growing.threads.prepare(1);
    expect(growing.use.enter()).toBe(true);
    expect(growing.threads.spawn(1000)).toBe(1);
    growing.use.leave();
    let inUse: boolean | undefined;
    growing.memory.onGrow = () => {
      // The thread returns meanwhile: the grow alone keeps the instance in use.
      expect(growing.workers[0]!.returns()).toBe(true);
      inUse = growing.use.retire();
    };
    expect(growing.threads.terminate()).toBe(false);
    expect(inUse).toBe(true);
    expect(growing.workers.some((worker) => worker.terminated)).toBe(false);
    expect(growing.keepers[0]!.terminated).toBe(false);

    // Once the runtime is retiring, it ends nothing and starts no thread worker, with threads or without.
    for (const withThread of [true, false]) {
      const retiring = host('ready-and-start');
      await retiring.threads.prepare(1);
      if (withThread) {
        expect(retiring.use.enter()).toBe(true);
        expect(retiring.threads.spawn(1000)).toBe(1);
        retiring.use.leave();
      }
      expect(retiring.use.retire()).toBe(withThread);
      expect(retiring.threads.terminate()).toBe(false);
      expect(retiring.memory.grows).toEqual([]);
      expect(retiring.workers.some((worker) => worker.terminated)).toBe(false);
      expect(retiring.keepers[0]!.terminated).toBe(false);
      await expect(retiring.threads.prepare(2)).rejects.toThrow('the runtime is ending');
      expect(retiring.workers).toHaveLength(1);
    }
  });

  it('a thread worker that fails while its thread runs counts the thread off', async () => {
    const { threads, workers, use } = host('ready-and-start');
    await threads.prepare(1);
    expect(use.enter()).toBe(true);
    expect(threads.spawn(1000)).toBe(1);
    use.leave();
    expect(use.uses()).toBe(1);
    workers[0]!.onerror?.({ message: 'boom' } as ErrorEvent);
    expect(use.uses()).toBe(0);
    expect(Atomics.load(workers[0]!.state!, 0)).toBe(SLOT_FAILED);
    await expect(threads.prepare(1)).rejects.toThrow('a thread worker failed: boom');
    expect(threads.terminate()).toBe(true);
  });

  it('fails to prepare when the keeper does not start in time, or fails', async () => {
    const late = host('ready-and-start', 100, 'never-ready');
    await expect(late.threads.prepare(1)).rejects.toThrow('the memory keeper did not start in time');
    late.threads.terminate();
    const failing = host('ready-and-start', 5_000, 'fails');
    await expect(failing.threads.prepare(1)).rejects.toThrow('the memory keeper failed: its script did not load');
    failing.threads.terminate();
  });

  it('fails to prepare when a thread worker does not become ready in time', async () => {
    const { threads } = host('never-ready');
    await expect(threads.prepare(1)).rejects.toBeInstanceOf(ThreadHostError);
    threads.terminate();
  });

  it('waits for a thread to start in slices of at most 100 ms, between which WebKit can end the worker', async () => {
    const { threads } = host('ready-never-start', 100, 'ready', 350);
    await threads.prepare(1);
    const wait = vi.spyOn(Atomics, 'wait');
    const start = performance.now();
    let timeouts: number[];
    try {
      expect(threads.spawn(1000)).toBe(-1);
      timeouts = wait.mock.calls.map((call) => call[3] as number);
    } finally {
      wait.mockRestore();
    }
    expect(performance.now() - start).toBeGreaterThanOrEqual(349);
    expect(timeouts.length).toBeGreaterThanOrEqual(4);
    expect(Math.max(...timeouts)).toBeLessThanOrEqual(100);
    threads.terminate();
  });

  it('abandons a thread that does not start in time, reports a fault, and the thread can no longer start', async () => {
    const { threads, workers, faultChannel, use } = host('ready-never-start');
    const faults = new BroadcastChannel(faultChannel);
    const received = new Promise<ThreadFault>((resolve) => {
      faults.onmessage = (event: MessageEvent<ThreadFault>) => resolve(event.data);
    });
    await threads.prepare(1);
    expect(threads.spawn(1000)).toBe(-1);
    const worker = workers[0]!;
    expect(Atomics.load(worker.pendingStart!, 0)).toBe(START_ABANDONED);
    expect(Atomics.load(worker.state!, 0)).toBe(SLOT_FAILED);
    // The thread that did not start is not counted.
    expect(use.uses()).toBe(0);
    expect(worker.startLate()).toBe(false);
    expect(await received).toEqual({ type: 'fault', tid: 1, message: 'The engine stopped: thread 1 did not start in time' });
    await expect(threads.prepare(1)).rejects.toThrow('thread 1 did not start');
    faults.close();
    threads.terminate();
  });
});
