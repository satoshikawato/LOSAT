import { expect, test } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { openSearch, paste, program } from './support/search';

test('the page is cross-origin isolated, so the threaded engine can use shared memory', async ({ page }) => {
  await page.goto('/');
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(true);
  expect(await page.evaluate(() => typeof SharedArrayBuffer)).toBe('function');
});

test(`paste, run, view and export with ${BUILD_HAS_ENGINE ? 'the engine' : 'the fake engine'}`, async ({ page }) => {
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
  expect(file.suggestedFilename()).toBe('losat-run1-tblastx.outfmt0.txt');
});

test('a BLASTN query of white space only is refused with the index scan error, before it is started', async ({ page }) => {
  // The adapter's register reads such a file as no records and the CLI then warns
  // "Query is Empty!", but the index scan (bio's reader) refuses it; the application follows
  // the scan and does not rebuild BLASTN's rules (S09, docs/evidence/losat_web_w1/README.md).
  await openSearch(page);
  await program(page, 'blastn');
  await paste(page, 'query', ' \n');
  await paste(page, 'subject', '>s1\nACGTACGTACGT\n');
  await page.getByTestId('run-losat').click();
  await expect(page.getByTestId('search-message')).toHaveText('Query (pasted): Expected > at record start.');
  await expect(page.getByTestId('query-source-0-error')).toHaveText('This input cannot be read: Expected > at record start.');
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);
  await expect(page.getByTestId('tab-search')).toHaveAttribute('aria-pressed', 'true');
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
