// A measurement, not a test (line NP WW1, DESIGN_ww1 §7): a large WindowMasker counts file (NCBI's
// windowmasker_files/9606/wmasker.obinary, 724 MB) in the application, on the serial and the
// threaded engine of each browser. It runs only with LOSAT_WEB_WM_HUMAN, a directory with that file
// (wmasker.obinary), a query (q100k.fa) and a subject (s1m.fa) of the same organism, and the engine
// build (LOSAT_WEB_REACTORS); it writes its records to LOSAT_WEB_EVIDENCE.
// Each engine: a search that reads the counts and one that keeps them, with the run's details as
// the application shows them (the time to give the file to the engine, the engine's memory after
// the search), the time from Run to the end, and the SHA-256 of the stored outfmt 6.
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, openFiles, openParameters, openSearch, program, submit, waitStatus } from './support/search';

const HUMAN = process.env['LOSAT_WEB_WM_HUMAN'] || undefined;
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const FLAG = 'window_masker_db';

test.skip(HUMAN === undefined, 'a measurement: runs only with LOSAT_WEB_WM_HUMAN');
test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
test.setTimeout(1_800_000);

const sha256 = (bytes: Uint8Array): string => createHash('sha256').update(bytes).digest('hex');

interface Searched {
  readonly run: number;
  readonly threads: string;
  readonly wallMs: number;
  readonly counts: string;
  readonly engine: string;
  readonly memory: string;
  readonly outfmt6Sha256: string;
}

async function search(page: Page, run: number, threads: string): Promise<Searched> {
  await page.getByTestId('threads').selectOption(threads);
  const started = Date.now();
  await submit(page);
  await waitStatus(page, run, 'completed', 1_500_000);
  const wallMs = Date.now() - started;
  await chooseRun(page, run);
  await page.getByTestId('results-download-all').click();
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('results-download-menu').getByTestId('results-download-6').click()]);
  const outfmt6Sha256 = sha256(readFileSync((await download.path())!));
  await page.getByTestId('results-view-details').click();
  const details = page.getByTestId('run-record');
  await expect(details).toBeVisible();
  const row = async (label: string): Promise<string> =>
    ((await details.locator('dt').filter({ hasText: new RegExp(`^\\s*${label}\\s*$`) }).locator('xpath=following-sibling::dd[1]').first().textContent()) ?? '').replace(/\s+/g, ' ').trim();
  const result = { run, threads, wallMs, counts: await row('WindowMasker counts'), engine: await row('Engine'), memory: await row('Engine memory'), outfmt6Sha256 };
  await page.getByTestId('tab-search').click();
  return result;
}

test('the human counts in the application', async ({ page, browserName }) => {
  const dir = HUMAN!;
  await openSearch(page);
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'q100k.fa', text: readFileSync(join(dir, 'q100k.fa')) }]);
  await openFiles(page, 'subject', [{ name: 's1m.fa', text: readFileSync(join(dir, 's1m.fa')) }]);
  await openParameters(page);
  await page.getByTestId(`param-${FLAG}`).check();
  await page.getByTestId(`param-${FLAG}-file`).setInputFiles(join(dir, 'wmasker.obinary'));
  await expect(page.getByTestId(`param-${FLAG}-name`)).toContainText('wmasker.obinary');
  const runs: Searched[] = [];
  // Serial (1 thread): read, then kept; threaded (4 threads), a new instance: read, then kept.
  for (const [index, threads] of ['1', '1', '4', '4'].entries()) runs.push(await search(page, index + 1, threads));
  const record = { browser: browserName, at: new Date().toISOString(), runs };
  console.log(JSON.stringify(record));
  if (EVIDENCE !== undefined) writeFileSync(join(EVIDENCE, `window-masker-human-${browserName}.json`), `${JSON.stringify(record, null, 2)}\n`);
  for (const searched of runs) expect(searched.outfmt6Sha256, `run ${searched.run}`).toBe(runs[0]!.outfmt6Sha256);
  expect(runs.map((searched) => /kept from the previous search/.test(searched.counts))).toEqual([false, true, false, true]);
});
