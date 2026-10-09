// ThreadHost (src/infra/engine-worker/thread-host.ts) with fake thread workers: preparing
// them, starting a thread, and a thread worker that does not start in time. The real
// thread workers run in the browsers (tests/e2e/engine.spec.ts).
import { Worker as NodeWorker } from 'node:worker_threads';
import { describe, expect, it } from 'vitest';
import {
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  START_ABANDONED,
  START_PENDING,
  START_TAKEN,
  type ThreadCommand,
  type ThreadFault,
} from '../../src/infra/engine-worker/protocol';
import { ThreadHost, ThreadHostError } from '../../src/infra/engine-worker/thread-host';

type Behaviour = 'ready-and-start' | 'ready-never-start' | 'never-ready';

class FakeThreadWorker {
  state: Int32Array | undefined;
  pendingStart: Int32Array | undefined;
  starts: number[] = [];
  terminated = false;
  onmessage: ((event: MessageEvent) => void) | null = null;
  onerror: ((event: ErrorEvent) => void) | null = null;

  constructor(private readonly behaviour: Behaviour) {}

  postMessage(command: ThreadCommand): void {
    if (command.type === 'prepare') {
      this.state = new Int32Array(command.slot);
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

  /** The late start of a worker whose host stopped waiting. */
  startLate(): boolean {
    return Atomics.compareExchange(this.pendingStart!, 0, START_PENDING, START_TAKEN) === START_PENDING;
  }

  terminate(): void {
    this.terminated = true;
  }
}

/** `readyTimeoutMs` bounds only the cases that fail; a test whose worker must become ready in time raises it. */
function host(behaviour: Behaviour, readyTimeoutMs = 100) {
  const workers: FakeThreadWorker[] = [];
  const faultChannel = `test-faults-${Math.random()}`;
  const threads = new ThreadHost({
    module: {} as WebAssembly.Module,
    memory: {} as WebAssembly.Memory,
    faultChannel,
    readyTimeoutMs,
    startTimeoutMs: 50,
    createWorker: () => {
      const worker = new FakeThreadWorker(behaviour);
      workers.push(worker);
      return worker as unknown as Worker;
    },
  });
  return { threads, workers, faultChannel };
}

describe('ThreadHost', () => {
  it('prepares two sets of thread workers and starts threads in them with increasing thread IDs', async () => {
    const { threads, workers } = host('ready-and-start');
    await threads.prepare(2);
    expect(workers).toHaveLength(4);
    expect(threads.spawn(1000)).toBe(1);
    expect(threads.spawn(2000)).toBe(2);
    threads.terminate();
    expect(workers.every((worker) => worker.terminated)).toBe(true);
  });

  it('starts the next pool of a search in the spare set while the threads of the previous pool still run', async () => {
    const { threads, workers } = host('ready-and-start');
    await threads.prepare(1);
    expect(threads.spawn(1000)).toBe(1);
    // The thread of the first pool has not returned (its slot is running): the next pool's
    // spawn takes the spare thread worker at once instead of waiting for it.
    const started = performance.now();
    expect(threads.spawn(2000)).toBe(2);
    expect(performance.now() - started).toBeLessThan(50);
    expect(workers.map((worker) => worker.starts)).toEqual([[1], [2]]);
    // With every thread worker running, a spawn fails at once: the threads it would wait
    // for may need this spawn to return before they can.
    const failing = performance.now();
    expect(threads.spawn(3000)).toBe(-1);
    expect(performance.now() - failing).toBeLessThan(50);
    threads.terminate();
  });

  it('waits for a thread worker that is being prepared again after its thread returned', async () => {
    // A Node worker can take more than 100 ms to start while the other test files run.
    const { threads, workers } = host('ready-and-start', 5_000);
    await threads.prepare(1);
    expect(threads.spawn(1000)).toBe(1);
    expect(threads.spawn(2000)).toBe(2);
    // The first thread returned; its worker becomes ready 30 ms later, while spawn waits.
    const state = workers[0]!.state!;
    Atomics.store(state, 0, SLOT_PREPARING);
    const ready = new NodeWorker(
      `const { workerData } = require('node:worker_threads');
       const state = new Int32Array(workerData);
       setTimeout(() => { Atomics.store(state, 0, ${SLOT_READY}); Atomics.notify(state, 0); }, 30);`,
      { eval: true, workerData: state.buffer },
    );
    expect(threads.spawn(3000)).toBe(3);
    expect(workers[0]!.starts).toEqual([1, 3]);
    await ready.terminate();
    threads.terminate();
  });

  it('fails to prepare when a thread worker does not become ready in time', async () => {
    const { threads } = host('never-ready');
    await expect(threads.prepare(1)).rejects.toBeInstanceOf(ThreadHostError);
    threads.terminate();
  });

  it('abandons a thread that does not start in time, reports a fault, and the thread can no longer start', async () => {
    const { threads, workers, faultChannel } = host('ready-never-start');
    const faults = new BroadcastChannel(faultChannel);
    const received = new Promise<ThreadFault>((resolve) => {
      faults.onmessage = (event: MessageEvent<ThreadFault>) => resolve(event.data);
    });
    await threads.prepare(1);
    expect(threads.spawn(1000)).toBe(-1);
    const worker = workers[0]!;
    expect(Atomics.load(worker.pendingStart!, 0)).toBe(START_ABANDONED);
    expect(Atomics.load(worker.state!, 0)).toBe(SLOT_FAILED);
    expect(worker.startLate()).toBe(false);
    expect(await received).toEqual({ type: 'fault', tid: 1, message: 'The engine stopped: thread 1 did not start in time' });
    await expect(threads.prepare(1)).rejects.toThrow('thread 1 did not start');
    faults.close();
    threads.terminate();
  });
});
