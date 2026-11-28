// Line WQA's D33 (2026-11-28), two low defects of sweep 4, through the real application (the engine build, LOSAT_WEB_REACTORS):
//  - S4-3: a list file of the Database filters whose name starts or ends with a space is run, and shown, under that name,
//    quoted as a query file is (the shown command runs on the command line with the same bytes);
//  - S4-8: Run LOSAT's "(N runs)" never shows the count of an older argv after a help or version word is typed or taken out;
//    the WindowMasker tab's button says one run for a help or version word, as the press starts one (R1, M3).
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
    // The line of the database's filter lists names the file the same way (R1, L1).
    await expect(page.getByTestId('run-database-filters')).toContainText(`${flag} '${name}'`);
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
    const w = window as unknown as { __labels: string[]; __widths: number[]; __from: number };
    w.__labels = [];
    w.__widths = [];
    w.__from = 0;
    document.addEventListener('input', () => void (w.__from = w.__labels.length), true);
    const sample = () => {
      const button = document.querySelector('[data-testid="run-losat"]');
      w.__labels.push((button?.textContent ?? '').trim());
      w.__widths.push(button?.getBoundingClientRect().width ?? 0);
      requestAnimationFrame(sample);
    };
    sample();
  });
  const sampled = async (edit: () => Promise<void>, settle: string): Promise<string[]> => {
    await edit();
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await expect(button).toHaveText(settle);
    return page.evaluate(() => {
      const w = window as unknown as { __labels: string[]; __widths: number[]; __from: number };
      // Nothing moves under a press: the button is as wide with every label.
      const widths = new Set(w.__widths.slice(Math.max(0, w.__from - 1)).map((width) => Math.round(width)));
      return [...w.__labels.slice(w.__from), ...(widths.size === 1 ? [] : [`WIDTHS ${[...widths].join(',')}`])];
    });
  };
  for (const word of ['-version', '-h']) {
    const typed = await sampled(() => page.getByTestId('param-extra').fill(word), 'Run LOSAT');
    expect(typed.filter((label) => label.includes('runs') || label.startsWith('WIDTHS')), `after ${word}`).toEqual([]);
    // Taken out: the label says no count until the answer, then the 2 runs again, in a button as wide as before.
    const taken = await sampled(() => page.getByTestId('param-extra').fill(''), 'Run LOSAT (2 runs)');
    expect(taken.filter((label) => label.startsWith('WIDTHS')), `after ${word} was taken out`).toEqual([]);
  }
  // An edit that does not print: the count is back after a round trip, and the button never changes its width.
  const ordinary = await sampled(() => page.getByTestId('param-extra').fill('-outfmt 17'), 'Run LOSAT (2 runs)');
  expect(ordinary.filter((label) => label.startsWith('WIDTHS')), 'after -outfmt 17').toEqual([]);
  await page.getByTestId('param-extra').fill('');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  // A press at once after the word starts one run, as the label that waits for the answer says.
  await page.getByTestId('param-extra').fill('-version');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
});

test('the WindowMasker button names no count for a help or version word, which starts one print request (S4-8, R1 M3)', async ({ page }) => {
  await page.getByTestId('tab-windowmasker').click();
  await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
  const files = ['a', 'b'].map((name) => ({
    name: `${name}.fa`,
    mimeType: 'application/octet-stream',
    buffer: Buffer.from(`>${name}\nACGTACGTACGTACGTACGTAAAC\n`),
  }));
  await page.getByTestId('wm-fasta-files').setInputFiles(files);
  await expect(page.locator('[data-testid="wm-fasta-list"] > li')).toHaveCount(2);
  await page.getByTestId('wm-fasta-mode-separate').check();
  const button = page.getByTestId('wm-run-losat');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  const toggle = page.getByTestId('wm-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await page.getByTestId('wm-param-extra').fill('-version');
  await expect(button).toHaveText('Run LOSAT');
  await page.getByTestId('wm-param-extra').fill('');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  await page.getByTestId('wm-param-extra').fill('-version');
  await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await button.click();
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
});
