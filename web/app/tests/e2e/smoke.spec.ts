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
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('result-command')).toHaveText(
    'LOSAT tblastx -query query.fa -subject subject.fa -outfmt 6',
  );
  await page.getByTestId('format-0').click();
  await expect(page.getByTestId('result-command')).toHaveText(
    'LOSAT tblastx -query query.fa -subject subject.fa -outfmt 0',
  );
  await expect(page.getByTestId('result-output')).toContainText(BUILD_HAS_ENGINE ? 'TBLASTX 2.17.0+' : 'FAKE ENGINE OUTPUT');

  const download = page.waitForEvent('download');
  await page.getByTestId('export-output').click();
  const file = await download;
  expect(file.suggestedFilename()).toBe('losat-run1-tblastx.outfmt0.txt');
});

test('a BLASTN query of white space only is refused with the index scan error, before it is queued', async ({ page }) => {
  // The adapter's register reads such a file as no records and the CLI then warns
  // "Query is Empty!", but the index scan (bio's reader) refuses it; the application follows
  // the scan and does not rebuild BLASTN's rules (S09, docs/evidence/losat_web_w1/README.md).
  await page.goto('/');
  await page.getByTestId('program-blastn').check();
  await page.getByTestId('query-input').fill(' \n');
  await page.getByTestId('subject-input').fill('>s1\nACGTACGTACGT\n');
  await page.getByTestId('add-to-queue').click();
  await expect(page.getByTestId('search-message')).toHaveText('Query FASTA: Expected > at record start.');
  await expect(page.getByTestId('run-1')).toHaveCount(0);
});
