// The threaded engine module that the site serves carries the shared-memory guard (release
// v0.2.0 DECISIONS 52): without it, Chromium can run a memory.fill or memory.copy of one
// thread against bounds that are stale after another thread grew the shared memory, and the
// search traps. The test reads the module from the production build (`vite preview`) and
// checks it with the guard's own verifier (LOSAT/tests/wasi_shared_memory.js), which accepts
// only the guard's helpers, call sites and counts. The serial module needs no guard.
import { createRequire } from 'node:module';
import { expect, test } from '@playwright/test';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';

test.skip(ENGINE === undefined, NO_ENGINE_REASON);

const guard = createRequire(import.meta.url)('../../../../LOSAT/tests/wasi_shared_memory.js') as {
  verifyGuarded(bytes: Uint8Array): { version: number; fill: number; copy: number; atomic: number };
};

test('the site serves the threaded engine module with the shared-memory guard', async ({ request, browserName }) => {
  // The bytes are the same for every browser.
  test.skip(browserName !== 'chromium', 'checked once, in the Chromium project');
  const response = await request.get(`/${ENGINE!.threads.name}`);
  expect(response.ok()).toBe(true);
  const bytes = await response.body();
  expect(bytes.equals(Buffer.from(ENGINE!.threads.bytes))).toBe(true);
  const record = guard.verifyGuarded(bytes);
  expect(record.fill + record.copy).toBeGreaterThan(0);
  expect(record.grow).toBeGreaterThan(0);
  expect(record.atomic).toBeGreaterThan(0);
});
