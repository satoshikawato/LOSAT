// The S10 data-layer contracts in Firefox and WebKit (S09): the BlockStore contract against
// OPFS and memory in a dedicated worker, the run output contract with the writer in its
// own worker and the application's real Data worker as the receiver, and the storage that
// the application chooses. Chromium runs the contracts in contracts.spec.ts, with a profile
// on disk and the storage-full cases, which need the Chrome DevTools Protocol to set the
// quota; here those cases are left out.
import { writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/profile';
import { BLOCK_STORE_CASES } from '../contract/block-store.contract';
import { RUN_OUTPUT_CASES } from '../contract/run-output.contract';
import { buildHarness, openApp, type HarnessFiles } from './support/browser';
import { openHarness, startHarnessServer, type HarnessServer } from './support/harness-server';

test.skip(({ browserName }) => browserName === 'chromium', 'Chromium runs these contracts in contracts.spec.ts');
test.setTimeout(120_000);

const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const withoutStorageFull = <T extends { name: string }>(cases: readonly T[]) => cases.filter((c) => !c.name.startsWith('storage full'));

let harness: HarnessFiles;
let server: HarnessServer;

test.beforeAll(async () => {
  harness = await buildHarness();
  server = await startHarnessServer(harness);
});

test.afterAll(async () => {
  await server?.close();
});

/** Runs the BlockStore contract in a dedicated worker of the harness page, without the storage-full cases. */
async function blockStoreContract(page: Page, backend: 'opfs' | 'memory'): Promise<void> {
  if (backend === 'opfs') {
    const opfs = await page.evaluate(() => typeof navigator.storage?.getDirectory === 'function');
    test.skip(!opfs, 'this browser has no Origin Private File System');
  }
  await page.exposeFunction('setStorageQuota', async () => undefined);
  const results = withoutStorageFull(await page.evaluate((kind) => window.losatHarness!.blockStore(kind), backend));
  expect(results.filter((result) => !result.ok)).toEqual([]);
  expect(results.map((result) => result.name)).toEqual(withoutStorageFull(BLOCK_STORE_CASES).map((c) => c.name));
}

// The OPFS store in a profile on disk, as users have it (support/profile.ts): WebKit refuses OPFS
// in Playwright's default, incognito contexts, where the application keeps results in memory and
// says so (the last test).
test('BlockStore contract: opfs store in a dedicated worker (without storage full)', async ({ profile }) => {
  await blockStoreContract(await openHarness(profile.context, server), 'opfs');
});

test('BlockStore contract: memory store in a dedicated worker (without storage full)', async ({ context }) => {
  await blockStoreContract(await openHarness(context, server), 'memory');
});

test('run output contract: a writer in its own worker and the real Data worker (without storage full)', async ({ context }) => {
  const page = await openHarness(context, server);
  await page.exposeFunction('setStorageQuota', async () => undefined);
  const { backend, results } = await page.evaluate(() => window.losatHarness!.runOutput());
  const kept = withoutStorageFull(results);
  expect(kept.filter((result) => !result.ok)).toEqual([]);
  expect(kept.map((result) => result.name)).toEqual(withoutStorageFull(RUN_OUTPUT_CASES).map((c) => c.name));
  console.log(`run output contract on the ${backend} store`);
});

test('the application keeps results in OPFS, or in memory with the reason shown', async ({ context, browserName }) => {
  const page = await openApp(context);
  const status = page.getByTestId('storage-status');
  const backend = await status.getAttribute('data-backend');
  expect(['opfs', 'memory']).toContain(backend);
  if (backend === 'memory') {
    await expect(status.getByTestId('storage-memory-note')).toContainText(
      'Results are kept in memory only and are lost when this tab is closed or reloaded (',
    );
  }
  const text = await status.innerText();
  if (EVIDENCE !== undefined) writeFileSync(join(EVIDENCE, `storage-${browserName}.json`), `${JSON.stringify({ backend, text }, null, 2)}\n`);
});
