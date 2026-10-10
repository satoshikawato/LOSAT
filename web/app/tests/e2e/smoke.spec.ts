import { expect, test } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';

test('the page is cross-origin isolated, so the threaded engine can use shared memory', async ({ page }) => {
  await page.goto('/');
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(true);
  expect(await page.evaluate(() => typeof SharedArrayBuffer)).toBe('function');
});

test(`paste, queue, run, view and export with ${BUILD_HAS_ENGINE ? 'the engine' : 'the fake engine'}`, async ({ page }) => {
  await page.goto('/');
  // The banner shows exactly when the build has no engine (web/AGENTS.md rule 7).
  await expect(page.getByTestId('fake-engine-banner')).toHaveCount(BUILD_HAS_ENGINE ? 0 : 1);

  await page.getByTestId('program-tblastx').check();
  await page.getByTestId('query-input').fill('>q1\nACGTACGTACGT\n');
  await page.getByTestId('subject-input').fill('>s1\nACGTACGTACGT\n');
  await page.getByTestId('add-to-queue').click();

  await expect(page.getByTestId('run-1-status')).toHaveText('completed');
  // The results tab shows the newest completed run; the outputs as written are in its Outputs view.
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  await page.getByTestId('results-view-outputs').click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:6');
  await expect(page.getByTestId('result-command')).toHaveText(
    'LOSAT tblastx -query query.fa -subject subject.fa -outfmt 6',
  );
  await page.getByTestId('format-0').click();
  await expect(page.getByTestId('result-command')).toHaveText(
    'LOSAT tblastx -query query.fa -subject subject.fa -outfmt 0',
  );
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:0');
  await expect(page.getByTestId('result-output')).toContainText(BUILD_HAS_ENGINE ? 'TBLASTX 2.17.0+' : 'FAKE ENGINE OUTPUT');

  const download = page.waitForEvent('download');
  await page.getByTestId('export-output').click();
  const file = await download;
  expect(file.suggestedFilename()).toBe('losat-run1-tblastx.outfmt0.txt');
});

test('a BLASTN query of white space only has no record: the engine reads it, and the run ends without a search', async ({ page }) => {
  // The index scan reads as the engine's reader does (kind 1 for a BLASTN query, S14): white
  // space only has no record, `register` accepts it without records, and a run gives NCBI's
  // "Query is Empty!" warning (docs/web/abi_v2.md §4).
  await page.goto('/');
  await page.getByTestId('program-blastn').check();
  await page.getByTestId('query-input').fill(' \n');
  await page.getByTestId('subject-input').fill('>s1\nACGTACGTACGT\n');
  await expect(page.getByTestId('query-source-0-summary')).toHaveText('0 records · 0 nt');
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'ok');
  await page.getByTestId('add-to-queue').click();
  await expect(page.getByTestId('run-1-status')).toHaveText('completed');
  if (BUILD_HAS_ENGINE) {
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-diagnostics')).toContainText('Query is Empty!');
  }
});
