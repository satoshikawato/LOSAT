// Line WQA's D20 (2026-11-23), two of the three low defects of sweep 3 (S3-3 is D21), through the real application (the engine build, LOSAT_WEB_REACTORS):
//  - S3-1: a run that printed the command line's help or version says nothing of a comparison of options, outfmt command
//    lines or a search's memory;
//  - S3-2: the button's "(N runs)" is the number of runs that Run LOSAT starts, also when the Additional options hold a
//    help or version word and when a query file has no record left;
import { expect } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, openFiles, openParameters, openSearch, paste, program, showRecords, submit, waitStatus } from './support/search';

test.setTimeout(240_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine prints no help or version: the engine build runs these tests');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** Words that Run details and the results header say only of a search (S3-1). */
const SEARCH_WORDS = [/Support status/, /Supported/, /not compared/, /outfmt \d+:/, /Engine memory/, /after the search/, /read for this search/, /kept from the previous search/];

test('a run that printed help or version says nothing of options compared, outfmt command lines or a search (S3-1)', async ({ page }) => {
  for (const [index, word] of ['-version', '-help'].entries()) {
    const run = index + 1;
    if (run > 1) await page.getByTestId('tab-search').click();
    await program(page, 'blastn');
    await openParameters(page);
    await page.getByTestId('param-extra').fill(word);
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    await waitStatus(page, run, 'completed');
    await chooseRun(page, run);
    // The header: the program, with no support badge.
    await expect(page.getByTestId('results-program')).toBeVisible();
    await expect(page.getByTestId('verification-badge')).toHaveCount(0);
    await expect(page.getByTestId('results-subjects')).toHaveCount(0);
    await page.getByTestId('results-view-details').click();
    const details = page.getByTestId('run-details');
    await expect(details).toBeVisible();
    await expect(page.getByTestId('run-printed')).toBeVisible();
    await expect(page.getByTestId('verification-details')).toHaveCount(0);
    const text = (await details.textContent()) ?? '';
    for (const words of SEARCH_WORDS) expect(text, `${word}: ${words}`).not.toMatch(words);
    await expect(page.locator('[data-testid^="run-command-"]')).toHaveCount(0);
    await expect(page.getByTestId('run-input-query')).toHaveCount(0);
    // What is true of the run stays.
    await expect(page.getByTestId('run-record')).toContainText('Took');
  }
  // A search of the same page still says what it compared: the badge, the support status and the command lines.
  await page.getByTestId('tab-search').click();
  await openParameters(page);
  await page.getByTestId('param-extra').fill('');
  await paste(page, 'query', '>q\nACGTACGTACGTACGTACGTAAAC\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTACGTAAAC\n');
  await submit(page);
  await waitStatus(page, 3, 'completed');
  await chooseRun(page, 3);
  await expect(page.getByTestId('results-summary').getByTestId('verification-badge')).toBeVisible();
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('verification-details')).toBeVisible();
  await expect(page.getByTestId('run-command-6')).toBeVisible();
  await expect(page.getByTestId('run-input-query')).toBeVisible();
});

test('the button counts the runs that Run LOSAT starts: a help or version word, and a query file with no record left (S3-2)', async ({ page }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [
    { name: 'a.fa', text: '>a\nACGTACGTACGTACGTACGTAAAC\n' },
    { name: 'b.fa', text: '>b\nTTTACGTACGTACGTACGTACGT\n' },
  ]);
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTACGTAAAC\n');
  await page.getByTestId('query-mode-separate').check();
  await openParameters(page);
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT (2 runs)');
  await expect(page.getByTestId('run-losat-bottom')).toHaveText('Run LOSAT (2 runs)');

  // The word prints once, as the CLI prints it: one run, and the button says so.
  await page.getByTestId('param-extra').fill('-version-full');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT');
  await expect(page.getByTestId('run-losat-bottom')).toHaveText('Run LOSAT');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);

  // Without the word the two files are two runs again; leave out the only record of one, and the button says one run
  // (the file adds nothing to the search) and one run starts.
  await page.getByTestId('tab-search').click();
  await openParameters(page);
  await page.getByTestId('param-extra').fill('');
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT (2 runs)');
  await showRecords(page, 'query', 1);
  await page.getByTestId('query-source-1-record-0').uncheck();
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT');
  await expect(page.getByTestId('run-losat-bottom')).toHaveText('Run LOSAT');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await expect(page.getByTestId('run-3-option')).toHaveCount(0);
});
