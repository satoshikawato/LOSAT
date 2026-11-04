import { expect, test } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { navigate } from './support/navigate';
import { chooseRun, openSearch, paste, program, waitStatus } from './support/search';

test('the page is cross-origin isolated, so the threaded engine can use shared memory', async ({ page }) => {
  await navigate(page, '/');
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(true);
  expect(await page.evaluate(() => typeof SharedArrayBuffer)).toBe('function');
});

test(`paste, run, view and export with ${BUILD_HAS_ENGINE ? 'the engine' : 'the fake engine'}`, async ({ page }) => {
  // With the engine, more than the default 30 s, as in the other specs that search with it: on
  // CI, WebKit took 10-23 s of the 30 s to open the page alone, and the run was still preparing
  // when the test ended (losat-baselines/web-flaky-20261103/runs).
  if (BUILD_HAS_ENGINE) test.setTimeout(120_000);
  await openSearch(page);
  // The banner shows exactly when the build has no engine (web/AGENTS.md rule 7).
  await expect(page.getByTestId('fake-engine-banner')).toHaveCount(BUILD_HAS_ENGINE ? 0 : 1);

  // Each step waits until the form stops moving, or the click on "Run LOSAT" can be lost (WF1).
  await program(page, 'tblastx');
  await paste(page, 'query', '>q1\nACGTACGTACGT\n');
  await paste(page, 'subject', '>s1\nACGTACGTACGT\n');
  await page.getByTestId('run-losat').click();

  // "Run LOSAT" opens the results page of the run: its status, then its results and, in the
  // Outputs view, the outputs as written.
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  await page.getByTestId('results-view-outputs').click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:6');
  await expect(page.getByTestId('result-command')).toHaveText(
    'losat tblastx -query query.fa -subject subject.fa -outfmt 6',
  );
  await page.getByTestId('format-0').click();
  await expect(page.getByTestId('result-command')).toHaveText(
    'losat tblastx -query query.fa -subject subject.fa -outfmt 0',
  );
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:0');
  await expect(page.getByTestId('result-output')).toContainText(BUILD_HAS_ENGINE ? 'TBLASTX 2.17.0+' : 'FAKE ENGINE OUTPUT');

  const download = page.waitForEvent('download');
  await page.getByTestId('export-output').click();
  const file = await download;
  expect(file.suggestedFilename()).toMatch(/^losat-run1-tblastx-\d{8}T\d{6}Z\.outfmt0\.txt$/);
});

test('a BLASTN query of white space only has no record, and its run ends as the CLI\'s does', async ({ page }) => {
  // The index scan reads with the program's reader (NCBI's, DW-38), which reads no record from
  // it; the CLI then warns "Query is Empty!". (Before, bio's reader refused it before the run.)
  await openSearch(page);
  await program(page, 'blastn');
  await paste(page, 'query', ' \n');
  await paste(page, 'subject', '>s1\nACGTACGTACGT\n');
  await expect(page.getByTestId('query-source-0-summary')).toContainText('0 records');
  await page.getByTestId('run-losat').click();
  await waitStatus(page, 1, 'completed');
  if (BUILD_HAS_ENGINE) {
    await chooseRun(page, 1);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-diagnostics')).toContainText('Query is Empty!');
  }
});

test('the empty Results tab offers "Run a search", which opens the Search tab', async ({ page }) => {
  await openSearch(page);
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('results-empty')).toBeVisible();
  await page.getByTestId('results-run-a-search').click();
  await expect(page.getByTestId('tab-search')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'false');
});
