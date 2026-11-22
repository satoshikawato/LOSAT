// Line WQA's D20 (2026-11-23), the three low defects of sweep 3, through the real application (the engine build,
// LOSAT_WEB_REACTORS; with LOSAT_WEB_NATIVE the bytes and the shown command are the native CLI's):
//  - S3-1: a run that printed the command line's help or version says nothing of a comparison of options, outfmt command
//    lines or a search's memory;
//  - S3-2: the button's "(N runs)" is the number of runs that Run LOSAT starts, also when the Additional options hold a
//    help or version word and when a query file has no record left;
//  - S3-3: the command shown for a compressed query or subject names the file as the user gave it, and run from a folder
//    with the user's files the native CLI writes the bytes that the application stored.
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync, zstdCompressSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, fasta, openFiles, openParameters, openSearch, paste, program, showRecords, storedOutputs, submit, waitStatus } from './support/search';

test.setTimeout(240_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine prints no help or version and inflates no file: the engine build runs these tests');

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

/** An input of a run: its bytes under the name that the user's file has. */
interface UserFile {
  readonly name: string;
  readonly bytes: Buffer;
}

/** The native CLI's outfmt 6 for the command that the application showed, run from a folder with the user's files. */
function nativeShown(command: string, files: readonly UserFile[]): Buffer {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'user-'));
  for (const file of files) writeFileSync(join(directory, file.name), file.bytes);
  const words = command.split(' ');
  expect(words[0], command).toBe('losat');
  const run = spawnSync(NATIVE!, [...words.slice(1), '-num_threads', '1'], { cwd: directory, maxBuffer: 1 << 28 });
  expect(run.status, `${command}: ${run.stderr.toString('utf8')}`).toBe(0);
  return run.stdout;
}

async function nameOfCommand(page: Page, run: number): Promise<string> {
  await chooseRun(page, run);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-6').click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${run}:6`);
  return ((await page.getByTestId('result-command').textContent()) ?? '').trim();
}

test('the command shown for a compressed query or subject names the file as the user gave it, and runs on the CLI (S3-3)', async ({ page }) => {
  const query = fasta('LC738884.fasta');
  const bz2 = spawnSync('bzip2', ['-9c'], { input: query, maxBuffer: 1 << 26 });
  expect(bz2.status, String(bz2.error ?? bz2.stderr)).toBe(0);
  const CASES: ReadonlyArray<{ readonly query: UserFile; readonly subject: UserFile }> = [
    { query: { name: 'q.fa.gz', bytes: gzipSync(query) }, subject: { name: 's.fa.zst', bytes: zstdCompressSync(query) } },
    { query: { name: 'q.fa.bz2', bytes: bz2.stdout }, subject: { name: 's.fa', bytes: query } },
  ];
  await program(page, 'blastn');
  for (const [index, c] of CASES.entries()) {
    const run = index + 1;
    if (run > 1) await page.getByTestId('tab-search').click();
    // The previous run's files go, the next run's come.
    if (run > 1) {
      await page.getByTestId('query-source-0-remove').click();
      await page.getByTestId('subject-source-0-remove').click();
    }
    await openFiles(page, 'query', [{ name: c.query.name, text: c.query.bytes }]);
    await openFiles(page, 'subject', [{ name: c.subject.name, text: c.subject.bytes }]);
    await submit(page);
    await waitStatus(page, run, 'completed');
    const expected = `losat blastn -query ${c.query.name} -subject ${c.subject.name} -outfmt 6`;
    // Outputs' command line, Run details' argument and command lines, and the printed argv of the run all name the user's files.
    expect(await nameOfCommand(page, run)).toBe(expected);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-argv')).toHaveText(`blastn -query ${c.query.name} -subject ${c.subject.name}`);
    await expect(page.getByTestId('run-command-6')).toHaveText(expected);
    await expect(page.getByTestId('run-command-0')).toHaveText(expected.replace('-outfmt 6', '-outfmt 0'));
    await expect(page.getByTestId('run-command-7')).toHaveText(expected.replace('-outfmt 6', '-outfmt 7'));
    if (NATIVE === undefined) test.info().annotations.push({ type: 'skip', description: `The CLI's bytes: ${NO_NATIVE_REASON}` });
    else {
      const stored = (await storedOutputs(page, run)).outputs[6];
      expect(stored.length).toBeGreaterThan(0);
      expect(nativeShown(expected, [c.query, c.subject]).equals(stored), `${expected}: the native CLI's standard output`).toBe(true);
    }
  }
});
