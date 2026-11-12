// A worker that the page ends while its thread waits on a shared memory gives that memory
// back. WebKit (Playwright's 2311) never ended such a worker when the wait had no timeout: it
// kept the worker's instance and its shared memory until the page closed, and a page can reserve
// a memory of 4 GiB about 15 times, so each cancel of a threaded search whose threads were
// waiting left one, until no threaded search could start. The shared-memory guard records each
// wait of the threaded module, and the keeper of the memory ends those waits before the Gateway
// ends the runtime: they then wait in slices with a call of the module's import sched_yield
// between them, where the browser ends the worker (LOSAT/tests/wasi_shared_memory.js,
// shared-memory.ts endWaits). This test ends more workers than WebKit could keep memories for,
// each with its keeper and waiting in a module that has the guard, as the Engine worker does;
// without the guard's waits or without the keeper's endWaits, WebKit cannot give the 16th its memory.
import { createRequire } from 'node:module';
import { expect } from '@playwright/test';
import { test } from './support/test';
import { buildHarness } from './support/browser';
import { openHarness, startHarnessServer, type HarnessServer } from './support/harness-server';

const guard = createRequire(import.meta.url)('../../../../LOSAT/tests/wasi_shared_memory.js') as {
  guardSharedMemory(bytes: Uint8Array): { bytes: Uint8Array; atomic: number };
};

/** The maximum that the threaded module declares: 4 GiB less one page. */
const MAXIMUM_PAGES = 65_535;
const ROUNDS = 20;

const uleb = (value: number): number[] => {
  const out: number[] = [];
  do {
    const byte = value & 127;
    value >>>= 7;
    out.push(byte | (value ? 128 : 0));
  } while (value);
  return out;
};
const vector = (items: number[][]): number[] => [...uleb(items.length), ...items.flat()];
const text = (value: string): number[] => [...uleb(value.length), ...Buffer.from(value)];
const section = (id: number, bytes: number[]): number[] => [id, ...uleb(bytes.length), ...bytes];

/**
 * Imports env.memory (shared, 1 to MAXIMUM_PAGES pages) and wasi_snapshot_preview1.sched_yield;
 * exports `park`, which waits at address 0 for the value 0 without a timeout.
 */
function waitModule(): Uint8Array {
  const park = [0, 0x41, 0, 0x41, 0, 0x42, 0x7f, 0xfe, 0x01, 2, 0, 0x1a, 0x0b];
  return Uint8Array.from([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    ...section(1, vector([[0x60, 0, 1, 0x7f], [0x60, 0, 0]])),
    ...section(2, vector([
      [...text('env'), ...text('memory'), 2, 3, 1, ...uleb(MAXIMUM_PAGES)],
      [...text('wasi_snapshot_preview1'), ...text('sched_yield'), 0, 0],
    ])),
    ...section(3, vector([[1]])),
    ...section(7, vector([[...text('park'), 0, 1]])),
    ...section(10, vector([[...uleb(park.length), ...park]])),
  ]);
}

let server: HarnessServer;
test.beforeAll(async () => {
  const guarded = guard.guardSharedMemory(waitModule());
  expect(guarded.atomic).toBe(1);
  server = await startHarnessServer(await buildHarness(), { extra: new Map([['wait.wasm', guarded.bytes]]) });
});
test.afterAll(async () => {
  await server?.close();
});

test('a worker that the page ends while it waits on a shared memory gives the memory back', async ({ context }) => {
  const page = await openHarness(context, server);
  const answers = await page.evaluate(
    ([url, rounds, pages]) => window.losatHarness!.engine.endWaitingWorkers(url, rounds, pages),
    ['/__extra/wait.wasm', ROUNDS, MAXIMUM_PAGES] as const,
  );
  expect(answers).toEqual(Array.from({ length: ROUNDS }, () => 'waiting'));
});
