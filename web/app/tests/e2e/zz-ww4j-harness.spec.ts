// WW4j throwaway probe (never committed): the harness build in a beforeAll hook, as makeblastdb.spec.ts does it (WW4g's
// Firefox sighting: "beforeAll" hook timeout of 30000ms exceeded). Each repetition runs in a new worker, so each builds.
// The build's time is logged by the instrumentation (WW4J_LOG_DIR/harness-build.log).
import { expect } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, buildHarness } from './support/browser';
import { startHarnessServer, type HarnessServer } from './support/harness-server';

test.skip(!BUILD_HAS_ENGINE, 'the real engine only');
test.setTimeout(600_000);

let server: HarnessServer | undefined;
const began = Date.now();
let built = 0;

test.beforeAll(async () => {
  server = await startHarnessServer(await buildHarness());
  built = Date.now() - began;
});

test.afterAll(async () => {
  await server?.close();
});

test('ww4j harness probe: the beforeAll hook builds and serves the harness', () => {
  console.log(`[ww4j-harness-probe] ${test.info().project.name} worker ${test.info().workerIndex}: hook ended ${built} ms after the file loaded`);
  expect(server).toBeDefined();
});
