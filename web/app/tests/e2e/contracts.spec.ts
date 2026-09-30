// Runs the contract suites of tests/contract in real browser workers (plan §2.1 L; S10):
// the BlockStore contract against OPFS and memory, and the run output contract with the
// writer in a separate worker and the application's real Data worker as the receiver.
// The storage-full cases set the quota through CDP, which OPFS honours only in a profile
// on disk (support/profile.ts).
import { BLOCK_STORE_CASES } from '../contract/block-store.contract';
import { RUN_OUTPUT_CASES } from '../contract/run-output.contract';
import { buildHarness, serveHarness, setQuota, type HarnessFiles } from './support/browser';
import { expect, test, type Profile } from './support/profile';

test.skip(({ browserName }) => browserName !== 'chromium', 'uses a Chromium profile and the Chrome DevTools Protocol');
test.setTimeout(120_000);

let harness: HarnessFiles;

test.beforeAll(async () => {
  harness = await buildHarness();
});

async function openHarness(profile: Profile) {
  await serveHarness(profile.context, harness);
  const page = await profile.context.newPage();
  // The harness logs each case as it finishes.
  page.on('console', (message) => console.log(message.text()));
  await page.exposeFunction('setStorageQuota', (bytes: number | null) => setQuota(page, bytes));
  await page.goto('/__harness/index.html');
  await page.waitForFunction(() => window.losatHarness !== undefined);
  return page;
}

for (const backend of ['opfs', 'memory'] as const) {
  test(`BlockStore contract: ${backend} store in a dedicated worker`, async ({ profile }) => {
    const page = await openHarness(profile);
    const results = await page.evaluate((kind) => window.losatHarness!.blockStore(kind), backend);
    expect(results.filter((result) => !result.ok)).toEqual([]);
    expect(results.map((result) => result.name)).toEqual(BLOCK_STORE_CASES.map((c) => c.name));
  });
}

test('run output contract: a writer in its own worker and the real Data worker (OPFS)', async ({ profile }) => {
  const page = await openHarness(profile);
  const { backend, results } = await page.evaluate(() => window.losatHarness!.runOutput());
  expect(backend).toBe('opfs');
  expect(results.filter((result) => !result.ok)).toEqual([]);
  expect(results.map((result) => result.name)).toEqual(RUN_OUTPUT_CASES.map((c) => c.name));
});
