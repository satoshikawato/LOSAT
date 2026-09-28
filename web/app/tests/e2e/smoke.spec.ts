import { expect, test } from '@playwright/test';

test('the page is cross-origin isolated, so the threaded engine can use shared memory', async ({ page }) => {
  await page.goto('/');
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(true);
  expect(await page.evaluate(() => typeof SharedArrayBuffer)).toBe('function');
});

test('paste, queue, run, view and export with the fake engine', async ({ page }) => {
  await page.goto('/');
  await expect(page.getByTestId('fake-engine-banner')).toBeVisible();

  await page.getByTestId('program-tblastx').check();
  await page.getByTestId('query-input').fill('>q1\nACGTACGTACGT\n');
  await page.getByTestId('subject-input').fill('>s1\nACGTACGTACGT\n');
  await page.getByTestId('add-to-queue').click();

  await expect(page.getByTestId('run-1-status')).toHaveText('completed');
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('result-command')).toHaveText(
    'LOSAT tblastx -query query.fa -subject subject.fa',
  );
  await page.getByTestId('format-0').click();
  await expect(page.getByTestId('result-output')).toContainText('FAKE ENGINE OUTPUT');

  const download = page.waitForEvent('download');
  await page.getByTestId('export-output').click();
  const file = await download;
  expect(file.suggestedFilename()).toBe('losat-run1-tblastx.outfmt0.txt');
});
