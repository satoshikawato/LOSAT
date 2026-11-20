// The use count of a threaded instance (web/packages/losat/src/infra/engine-worker/instance-use.ts), which the
// Gateway, the Engine worker and its thread workers share; their use of it runs in the browsers
// (tests/e2e/engine.spec.ts).
import { Worker as NodeWorker } from 'node:worker_threads';
import { describe, expect, it } from 'vitest';
import { InstanceUse, USE_RETIRING } from '../../../packages/losat/src/infra/engine-worker/instance-use';
import {
  SLOT_ENDED,
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  SLOT_RUNNING,
} from '../../../packages/losat/src/infra/engine-worker/protocol';

function runningSlot(): Int32Array {
  const slot = new Int32Array(new SharedArrayBuffer(4));
  Atomics.store(slot, 0, SLOT_RUNNING);
  return slot;
}

describe('InstanceUse', () => {
  it('ends at once when nothing uses the instance, and then lets no call start', () => {
    const use = InstanceUse.create();
    expect(use.enter()).toBe(true);
    use.leave();
    expect(use.retire()).toBe(false);
    expect(use.retiring()).toBe(true);
    expect(use.enter()).toBe(false);
    expect(use.uses()).toBe(0);
  });

  it('is in use while a call or a thread runs; the calls and threads count off after the mark', () => {
    const use = InstanceUse.create();
    expect(use.enter()).toBe(true);
    const slot = runningSlot();
    use.startThread();
    expect(use.uses()).toBe(2);
    expect(use.retire()).toBe(true);
    // A second call does not start; the ones that run end as before.
    expect(use.enter()).toBe(false);
    use.leave();
    expect(use.endThread(slot, SLOT_PREPARING)).toBe(true);
    expect(use.uses()).toBe(0);
    expect(use.retiring()).toBe(true);
  });

  it('counts a thread off once, by whichever side moves its slot from running first', () => {
    const use = InstanceUse.create();
    expect(use.enter()).toBe(true);
    const slots = [runningSlot(), runningSlot()];
    use.startThread();
    use.startThread();
    use.leave();
    // The thread worker's thread returned before the ThreadHost ended its worker ...
    expect(use.endThread(slots[0]!, SLOT_PREPARING)).toBe(true);
    expect(use.endThread(slots[0]!, SLOT_ENDED)).toBe(false);
    // ... and the ThreadHost ended this one before its thread returned.
    expect(use.endThread(slots[1]!, SLOT_ENDED)).toBe(true);
    expect(use.endThread(slots[1]!, SLOT_FAILED)).toBe(false);
    expect(slots.map((slot) => Atomics.load(slot, 0))).toEqual([SLOT_PREPARING, SLOT_ENDED]);
    expect(use.uses()).toBe(0);
    // A slot that does not run has no thread to count off.
    const ready = new Int32Array(new SharedArrayBuffer(4));
    Atomics.store(ready, 0, SLOT_READY);
    expect(use.endThread(ready, SLOT_ENDED)).toBe(false);
    expect(use.uses()).toBe(0);
  });

  it('a last use that ends after the mark stays counted, so that its caller ends nothing the keeper needs', () => {
    const use = InstanceUse.create();
    expect(use.enter()).toBe(true);
    expect(use.leaveLast()).toBe(true);
    expect(use.uses()).toBe(0);
    expect(use.enter()).toBe(true);
    expect(use.retire()).toBe(true);
    expect(use.leaveLast()).toBe(false);
    expect(use.uses()).toBe(1);
  });

  it('a last use waits for a thread that another thread is counting off', async () => {
    const use = InstanceUse.create();
    expect(use.enter()).toBe(true);
    const slot = runningSlot();
    use.startThread();
    // The thread worker has moved the slot, and lowers the count 30 ms later (endThread's second step).
    expect(Atomics.compareExchange(slot, 0, SLOT_RUNNING, SLOT_PREPARING)).toBe(SLOT_RUNNING);
    const done = new Int32Array(new SharedArrayBuffer(4));
    const late = new NodeWorker(
      `const { workerData } = require('node:worker_threads');
       const [word, done] = [new Int32Array(workerData.word, 0, 1), new Int32Array(workerData.done)];
       setTimeout(() => { Atomics.store(done, 0, 1); Atomics.sub(word, 0, 1); }, 30);`,
      { eval: true, workerData: { word: use.buffer, done: done.buffer } },
    );
    expect(use.leaveLast()).toBe(true);
    expect(Atomics.load(done, 0)).toBe(1);
    expect(use.uses()).toBe(0);
    await late.terminate();
  });

  it('refuses to end a last use that the count no longer holds', () => {
    const use = InstanceUse.create();
    expect(() => use.leaveLast()).toThrow('lost a use');
  });

  it('shares the count through the buffer, and keeps the mark out of it', () => {
    const gateway = InstanceUse.create();
    const engine = new InstanceUse(gateway.buffer);
    expect(engine.enter()).toBe(true);
    expect(gateway.uses()).toBe(1);
    expect(gateway.retire()).toBe(true);
    expect(new Int32Array(gateway.buffer)[0]).toBe(USE_RETIRING + 1);
    expect(engine.retiring()).toBe(true);
    engine.leave();
    expect(gateway.uses()).toBe(0);
  });
});
