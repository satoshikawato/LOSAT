// Backpressure of a tool run (D-WW3-11, M4; src/infra/engine-worker/tool-credit.ts, tool-posts.ts): the Engine worker waits in
// the engine's call while the bytes that it posted and the page has not taken are above the limit. The wait is real, and it
// is the code of the Engine worker: a worker thread of Node runs `toolPosts` (the function that engine-worker.ts gives the
// engine) with a small limit and stays blocked in `Atomics.wait` until the main thread takes bytes (a slice of 10 minutes, so
// that only the main thread's notify can wake it within the test's time) or closes the run. A test fails if either `wait()` of
// tool-posts.ts is removed (the chunk's or the file frame's), or if the run's credit is not given to the posts.
import { fileURLToPath } from 'node:url';
import { Worker } from 'node:worker_threads';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { createToolCredit, postedBytes, TOOL_CREDIT_LIMIT, ToolCreditPage, ToolCreditWorker } from '../src/infra/engine-worker/tool-credit.ts';
import { toolPosts } from '../src/infra/engine-worker/tool-posts.ts';
import type { EngineEvent, WorkerTool } from '../src/infra/engine-worker/protocol.ts';

const workers: Worker[] = [];
afterEach(async () => {
  await Promise.all(workers.splice(0).map((worker) => worker.terminate()));
});

/** Starts the fixture worker thread and collects the steps that it reports. */
function start(buffer: SharedArrayBuffer, options: { limit?: number; sliceMs?: number; steps: { kind: 'chunk' | 'file'; size: number }[] }) {
  const steps: string[] = [];
  const worker = new Worker(new URL('./support/credit-worker.ts', import.meta.url), {
    workerData: { buffer, ...options },
    execArgv: ['--import', fileURLToPath(new URL('./support/register-ts.mjs', import.meta.url))],
  });
  workers.push(worker);
  worker.on('message', (step: string) => steps.push(step));
  worker.on('error', (error: Error) => steps.push(`error ${error.message}`));
  return { steps };
}

const pause = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));
/** Waits for a step of the worker thread, which starts slowly on a loaded machine. */
const until = (check: () => void) => vi.waitFor(check, { timeout: 30_000 });

describe('the credit of a tool run', () => {
  it('blocks the worker over the limit until the page takes bytes, and a close lets it go (a chunk, then a file frame, then a chunk)', async () => {
    const buffer = new SharedArrayBuffer(8);
    const page = new ToolCreditPage(buffer);
    const { steps } = start(buffer, {
      limit: 10,
      sliceMs: 600_000,
      steps: [
        { kind: 'chunk', size: 8 },
        { kind: 'file', size: 8 },
        { kind: 'chunk', size: 8 },
      ],
    });
    await until(() => expect(steps).toContain('after 1'));
    // 16 bytes over a limit of 10: the file frame is posted, and the worker waits for the page in the post of the frame.
    await until(() => expect(steps).toContain('posted 2'));
    await pause(300);
    expect(steps).toEqual(['posted 1', 'after 1', 'posted 2']);
    page.take(8);
    await until(() => expect(steps).toContain('after 2'));
    // The next chunk makes 16 again: the worker waits in the post of the chunk.
    await until(() => expect(steps).toContain('posted 3'));
    await pause(300);
    expect(steps.at(-1)).toBe('posted 3');
    // A run that ends (a cancel, a failure) closes the cell: the worker goes on without the bytes.
    page.close();
    await until(() => expect(steps).toContain('done'));
    expect(steps).toEqual(['posted 1', 'after 1', 'posted 2', 'after 2', 'posted 3', 'after 3', 'done']);
  });

  it('waits above 64 MiB and not at it, with the limit that the Engine worker has', async () => {
    expect(TOOL_CREDIT_LIMIT).toBe(64 * 1024 * 1024);
    const buffer = new SharedArrayBuffer(8);
    const page = new ToolCreditPage(buffer);
    // No limit is given: the posts take the constant, as in the Engine worker.
    const { steps } = start(buffer, {
      steps: [
        { kind: 'chunk', size: TOOL_CREDIT_LIMIT },
        { kind: 'file', size: 1 },
      ],
    });
    await until(() => expect(steps).toContain('posted 2'));
    await pause(300);
    // The limit itself does not block (after 1); one byte over does (no after 2).
    expect(steps).toEqual(['posted 1', 'after 1', 'posted 2']);
    page.take(1);
    await until(() => expect(steps).toContain('done'));
  });

  it('does not wait when the page takes what it is given', async () => {
    const buffer = new SharedArrayBuffer(8);
    const worker = new ToolCreditWorker(buffer, 10, 10);
    const page = new ToolCreditPage(buffer);
    for (let i = 0; i < 1000; i++) {
      worker.add(8);
      page.take(8);
      worker.wait();
    }
    expect(new Int32Array(buffer)[0]).toBe(0);
  });

  it('has a buffer where the page can share memory', () => {
    expect(createToolCredit()).toBeInstanceOf(SharedArrayBuffer);
  });
});

describe('the posts of a tool run in the Engine worker', () => {
  const command = (credit?: SharedArrayBuffer): WorkerTool => ({ type: 'tool', generation: 4, runId: 'r', argv: ['windowmasker'], ...(credit === undefined ? {} : { credit }) });

  it('counts the bytes of a chunk or a frame of bytes before the post (the post transfers the buffer), and counts no other frame', () => {
    const buffer = new SharedArrayBuffer(8);
    const cell = new Int32Array(buffer);
    const posted: EngineEvent[] = [];
    const posts = toolPosts(command(buffer), 4, (event, transfer) => {
      // A transfer empties the sender's buffer: the length is read before.
      for (const item of transfer) structuredClone(item, { transfer: [item] });
      posted.push(event);
    });
    posts.chunk('stdout', new Uint8Array(100));
    posts.chunk('stderr', new Uint8Array(20));
    posts.file({ number: 0, kind: 'begin', name: 'x' });
    posts.file({ number: 0, kind: 'bytes', offset: 0, bytes: new Uint8Array(1000) });
    posts.file({ number: 0, kind: 'truncate', length: 0 });
    posts.file({ number: 0, kind: 'end' });
    expect(cell[0]).toBe(1120);
    expect(posted.map((event) => event.type)).toEqual(['tool-chunk', 'tool-chunk', 'tool-file', 'tool-file', 'tool-file', 'tool-file']);
    expect(posted[0]).toMatchObject({ generation: 4, runId: 'r', stream: 'stdout' });
    // The page counts the same bytes of the same frames when it takes them (the buffer survived the transfer there).
    expect(postedBytes({ number: 0, kind: 'bytes', offset: 0, bytes: new Uint8Array(1000) })).toBe(1000);
    expect(postedBytes({ number: 0, kind: 'end' })).toBe(0);
  });

  it('counts nothing and waits for nothing without a credit', () => {
    const posts = toolPosts(command(), 1, () => undefined);
    posts.chunk('stdout', new Uint8Array(10));
    posts.file({ number: 0, kind: 'bytes', offset: 0, bytes: new Uint8Array(10) });
  });
});
