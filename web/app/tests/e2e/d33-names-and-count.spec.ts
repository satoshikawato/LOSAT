// Line WQA's D33 (2026-11-28), two low defects of sweep 4, through the real application (the engine build, LOSAT_WEB_REACTORS):
//  - S4-3: a list file of the Database filters whose name starts or ends with a space is run, and shown, under that name,
//    quoted as a query file is (the shown command runs on the command line with the same bytes);
//  - S4-8: Run LOSAT's "(N runs)" never shows the count of an older argv after a help or version word is typed or taken out.
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { openFiles, openParameters, openSearch, paste, program, submit, waitStatus } from './support/search';

test.setTimeout(240_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine runs no database search and prints no help or version: the engine build runs these tests');

const TESTS = join(REPOSITORY, 'LOSAT/tests/db_search');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/**
 * Two list files of the Database filters, each run by itself (the command line refuses -seqidlist with -negative_seqidlist,
 * -gilist and the others together; the unit tests name every list file): the flag, the box's test id, and a name that starts
 * or ends with a space.
 */
const LISTS = [
  ['-seqidlist', 'seqidlist', 'ids.txt ', 'seqidlist/sl_tfp.txt'],
  ['-negative_seqidlist', 'negative_seqidlist', ' neg ids.txt', 'seqidlist/sl_tfp_neg.txt'],
] as const;

test('a list file whose name starts or ends with a space is run and shown under its name, quoted (S4-3)', async ({ page }) => {
  const database = readdirSync(join(TESTS, 'taxfilter'))
    .filter((name) => name.startsWith('tfp.'))
    .map((name) => ({ name, mimeType: 'application/octet-stream', buffer: readFileSync(join(TESTS, 'taxfilter', name)) }));
  await program(page, 'blastp');
  await paste(page, 'query', readFileSync(join(TESTS, 'taxfilter/q.faa'), 'utf8'));
  await page.getByTestId('subject-source-database').check();
  await page.getByTestId('subject-db-files-input').setInputFiles(database);
  await expect(page.getByTestId('subject-db-selected')).toContainText('tfp');
  await page.getByTestId('subject-db-filters-toggle').click();
  for (const [run, [flag, id, name, source]] of LISTS.entries()) {
    if (run > 0) {
      await page.getByTestId('tab-search').click();
      await page.getByTestId(`db-filter-${LISTS[run - 1]![1]}-remove`).click();
    }
    await page.getByTestId(`db-filter-${id}-input`).setInputFiles({ name, mimeType: 'text/plain', buffer: readFileSync(join(TESTS, source)) });
    await expect(page.getByTestId(`db-filter-${id}-name`)).toContainText(name.trim());
    await submit(page);
    await waitStatus(page, run + 1, 'completed');
    await page.getByTestId('tab-results').click();
    await page.getByTestId('results-view-details').click();
    // The command names the file as the user's file is called, quoted: it runs on the command line next to the file.
    const shown = (await page.getByTestId('run-argv').textContent()) ?? '';
    expect(shown, flag).toContain(`${flag} '${name}'`);
  }
});

test('the button names no count of an older argv after a help or version word is typed or taken out (S4-8)', async ({ page }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'a.fa', text: '>a\nACGTACGTACGTACGTACGTAAAC\n' }]);
  await openFiles(page, 'subject', [
    { name: 's1.fa', text: '>s1\nACGTACGTACGTACGTACGTAAAC\n' },
    { name: 's2.fa', text: '>s2\nTTTACGTACGTACGTACGTACGT\n' },
  ]);
  await page.getByTestId('subject-mode-separate').check();
  await openParameters(page);
  const button = page.getByTestId('run-losat');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  // A sampler of the button's text in every frame, from the edit to the engine's answer.
  // `from` is where the samples start that come after the edit's input event.
  await page.evaluate(() => {
    const w = window as unknown as { __labels: string[]; __from: number };
    w.__labels = [];
    w.__from = 0;
    document.addEventListener('input', () => void (w.__from = w.__labels.length), true);
    const sample = () => {
      w.__labels.push((document.querySelector('[data-testid="run-losat"]')?.textContent ?? '').trim());
      requestAnimationFrame(sample);
    };
    sample();
  });
  const sampled = async (edit: () => Promise<void>, settle: string): Promise<string[]> => {
    await edit();
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await expect(button).toHaveText(settle);
    return page.evaluate(() => {
      const w = window as unknown as { __labels: string[]; __from: number };
      return w.__labels.slice(w.__from);
    });
  };
  for (const word of ['-version', '-h']) {
    const typed = await sampled(() => page.getByTestId('param-extra').fill(word), 'Run LOSAT');
    expect(typed.filter((label) => label.includes('runs')), `after ${word}`).toEqual([]);
    // Taken out: the label says no count until the answer, then the 2 runs again.
    await sampled(() => page.getByTestId('param-extra').fill(''), 'Run LOSAT (2 runs)');
  }
  // A press at once after the word starts one run, as the label that waits for the answer says.
  await page.getByTestId('param-extra').fill('-version');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
});
