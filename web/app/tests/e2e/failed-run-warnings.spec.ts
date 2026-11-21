// A failed run keeps the warnings that the command line writes to its standard error before its
// error, and the results page shows them under the failure, as the engine wrote them (line WQA's
// sweep, D3, 2026-11-19; docs/web/abi_v2.md §5, stream 3). In the engine build, a BLASTN subject
// whose one record has no residues gives NCBI's `Subject sequence contains no data` warning and then
// its error `The average subject length is too short` (exit status 3); with LOSAT_WEB_NATIVE the
// warnings shown and the error are the native CLI's standard error for the run's own command line,
// line for line. The FakeEngine fails a run of an empty input with a marked warning. A run of two
// searches (ABI §7) has the outcome of the search of its own -outfmt: it fails as the command line
// fails, with that search's warnings only (D3 review, 2026-11-19), and the formatting warnings of its
// own format only (D3 review r2, 2026-11-19). A completed run whose -outfmt is structured shows the
// queries' warnings that the command line writes for the other formats that the run writes (D3
// review r3, 2026-11-19).
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page, type TestInfo } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openFiles, openParameters, openSearch, paste, pickRun, program, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** Pseudo-random nucleotides (a linear congruential generator): the same letters in every run. */
function dna(seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += 'ACGT'[(state >>> 16) % 4];
  }
  return out;
}

const QUERY = `>q1 first query\n${dna(7, 300)}\n`;
const EMPTY_RECORD = '>empty record\n';
const ERROR = 'BLAST engine error: The average subject length is too short';
const PROTEIN_QUERY = '>p\nMKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMAWRSTY\n';

/** The options of the run shown (`results-run-options`): its command line after the program and the inputs. */
async function runOptions(page: Page): Promise<string[]> {
  const text = (await page.getByTestId('results-run-options').innerText()).trim();
  return text === 'defaults' ? [] : text.split(/\s+/);
}

/**
 * The native CLI's run of the run shown: the same program, inputs and options (`runOptions`), one
 * thread; undefined, with a note, without LOSAT_WEB_NATIVE.
 */
async function nativeRun(page: Page, testInfo: TestInfo, mode: string, query: string, subject: string) {
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The native standard error: ${NO_NATIVE_REASON}` });
    return undefined;
  }
  const options = await runOptions(page);
  const directory = mkdtempSync(join(testInfo.outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  const cli = spawnSync(NATIVE, [mode, '-query', 'query.fa', '-subject', 'subject.fa', ...options, '-num_threads', '1'], { cwd: directory });
  return { status: cli.status, stderr: String(cli.stderr) };
}

test('a failed run shows the warnings written before its error, as the command line writes them', async ({ page }, testInfo) => {
  await program(page, 'blastn');
  await paste(page, 'query', QUERY);
  if (BUILD_HAS_ENGINE) await paste(page, 'subject', EMPTY_RECORD);
  else await openFiles(page, 'subject', [{ name: 'empty.fa', text: '' }]);
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  const status = page.getByTestId('results-status');
  await expect(status).toHaveAttribute('data-run-status', 'failed');
  // The command line's outfmt 0 (the run's own) holds its prolog before the error (D15, 2026-11-20).
  await expect(status).toContainText('It keeps the outputs written before the error. It keeps the warnings written before the error.');
  const warnings = page.getByTestId('results-failed-diagnostics');
  await expect(page.getByTestId('results-failed-warnings').getByRole('heading', { name: 'Warnings' })).toBeVisible();
  if (!BUILD_HAS_ENGINE) {
    await expect(status).toContainText('Run 1 failed: No sequence input was provided.');
    await expect(warnings).toHaveText('Warning: FAKE ENGINE OUTPUT - not a LOSAT search result');
    return;
  }
  await expect(status).toContainText(`Run 1 failed: ${ERROR}.`);
  const shown = await warnings.evaluate((element) => element.textContent ?? '');
  expect(shown).toBe('Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\n');

  // The standard error of the command line that ran: the warnings shown, then the error.
  const cli = await nativeRun(page, testInfo, 'blastn', QUERY, EMPTY_RECORD);
  if (cli !== undefined) {
    expect(cli.status).toBe(3);
    expect(cli.stderr).toBe(`${shown}${ERROR}\n`);
  }

  // The run keeps them: they come back with the run after another run was shown.
  await page.getByTestId('tab-search').click();
  await paste(page, 'subject', `>s1\n${dna(7, 300).slice(20, 280)}\n`);
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await pickRun(page, 2);
  await expect(page.getByTestId('results-failed-warnings')).toHaveCount(0);
  await pickRun(page, 1);
  await expect(warnings).toHaveText('Warning: [blastn] Subject_1 empty record: Subject sequence contains no data');
});

// BLASTN with -num_descriptions 5 makes two searches (docs/web/abi_v2.md §7), and both fail on the
// subject whose one record has no residues. The command line writes the formatting handler's
// `-num_descriptions ... ignored` warning for -outfmt 6 and 7 only: the run, whose own format is
// outfmt 0, shows the subject's warning alone, as the command line writes it (D3 review, 2026-11-19).
test('a failed run of two searches shows the warnings of its own command line only', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'two searches of the engine (engine build)');
  await program(page, 'blastn');
  await openParameters(page);
  await page.getByTestId('param-num_descriptions').fill('5');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', EMPTY_RECORD);
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  await expect(page.getByTestId('results-status')).toContainText(`Run 1 failed: ${ERROR}.`);
  const shown = await page.getByTestId('results-failed-diagnostics').evaluate((element) => element.textContent ?? '');
  expect(shown).toBe('Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\n');
  expect(await runOptions(page)).toEqual(['-num_descriptions', '5']);
  const cli = await nativeRun(page, testInfo, 'blastn', QUERY, EMPTY_RECORD);
  if (cli === undefined) return;
  expect(cli.status).toBe(3);
  expect(cli.stderr).toBe(`${shown}${ERROR}\n`);
});

// BLASTP with -num_alignments 0 makes two searches (docs/web/abi_v2.md §7): the tabular one fails
// alone before it reads the subject (NCBI's `No hits are being saved`, which the command line writes
// for -outfmt 6 and 7 only), and the format 0 search, that of the run's own outfmt 0, warns of the
// subject's empty first record and fails the run on the title of its second (`, ,`, where NCBI
// BLAST+ 2.17.0 crashes and LOSAT refuses). The run shows the warnings of the format 0 search alone,
// as the command line writes them for its command line (D3 review, 2026-11-19).
test('a run of two searches that its format 0 search fails shows the warnings of its own command line', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the engine fails the run on the title (engine build)');
  const subject = '>empty record\n>, ,\nMKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMA\n';
  await program(page, 'blastp');
  await openParameters(page);
  await page.getByTestId('param-num_alignments').fill('0');
  await paste(page, 'query', PROTEIN_QUERY);
  await paste(page, 'subject', subject);
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  const status = page.getByTestId('results-status');
  await expect(status).toContainText('Run 1 failed: subject record 2 has a title of punctuation');
  await expect(status).toContainText('It keeps the warnings written before the error.');
  const shown = await page.getByTestId('results-failed-diagnostics').evaluate((element) => element.textContent ?? '');
  expect(shown).toBe('Warning: [blastp] Subject_1 empty record: Subject sequence contains no data\n');
  expect(await runOptions(page)).toEqual(['-num_alignments', '0']);
  const cli = await nativeRun(page, testInfo, 'blastp', PROTEIN_QUERY, subject);
  if (cli === undefined) return;
  // The warning shown, then the error that ends the run.
  expect(cli.status).toBe(1);
  const lines = cli.stderr.split(/(?<=\n)/);
  expect(lines.at(-1)).toMatch(/^Error: subject record 2 has a title of punctuation/);
  expect(cli.stderr).toBe(`${shown}${lines.at(-1)}`);
});

// The run's own -outfmt 6 with -num_alignments 0 (a hit list of 0 for the tabular formats): the
// command line fails with NCBI's error after the formatting handler's warning (exit status 1), though
// outfmt 0, which the run also writes, has a hit list of its own. The form refuses the command line
// with that error before any run, as it refuses BLASTN's and BLASTX's (D3 review r2, 2026-11-19:
// BLASTP's validation let it run, and with -evalue 0 gave the e-value's error, which the command
// line gives only after the hit list's).
test('a command line whose own -outfmt has no hit list is refused with the command line\'s error', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the engine checks the hit list (engine build)');
  const subject = '>s1 ok\nMKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMA\n';
  const error = 'BLAST query/options error: No hits are being saved';
  await program(page, 'blastp');
  await paste(page, 'query', PROTEIN_QUERY);
  await paste(page, 'subject', subject);
  await openParameters(page);
  await page.getByTestId('param-num_alignments').fill('0');
  await page.getByTestId('param-evalue').fill('0');
  await page.getByTestId('param-extra').fill('-outfmt 6');
  const validation = page.getByTestId('argv-validation');
  await expect(validation).toHaveAttribute('data-state', 'invalid');
  await expect(validation).toContainText(error);
  await expect(validation).not.toContainText('expect value');
  await submit(page);
  await expect(page.getByTestId('search-message')).toContainText(error);
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The native standard error: ${NO_NATIVE_REASON}` });
    return;
  }
  const directory = mkdtempSync(join(testInfo.outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), PROTEIN_QUERY);
  writeFileSync(join(directory, 'subject.fa'), subject);
  const cli = spawnSync(NATIVE, ['blastp', '-query', 'query.fa', '-subject', 'subject.fa', '-num_alignments', '0', '-evalue', '0', '-outfmt', '6', '-num_threads', '1'], { cwd: directory });
  expect(cli.status).toBe(1);
  expect(String(cli.stderr)).toBe(`Warning: [blastp] Examining 5 or more matches is recommended\n${error}\nPlease refer to the LOSAT user manual.\n`);
});

// The formatting handler warns of -line_length for output formats above 4 only (NCBI's
// blast_args.cpp): a failed run whose own format is outfmt 0 shows the subject's warning alone, as
// the command line writes it, not the -line_length warning of the tabular formats that the run also
// writes (D3 review r2, 2026-11-19).
test('a failed run shows no formatting warning of the formats that its command line does not ask for', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the formatting warnings of the engine (engine build)');
  await program(page, 'blastn');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', EMPTY_RECORD);
  await openParameters(page);
  await page.getByTestId('param-line_length').fill('80');
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  await expect(page.getByTestId('results-status')).toContainText(`Run 1 failed: ${ERROR}.`);
  const shown = await page.getByTestId('results-failed-diagnostics').evaluate((element) => element.textContent ?? '');
  expect(shown).toBe('Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\n');
  expect(await runOptions(page)).toEqual(['-line_length', '80']);
  const cli = await nativeRun(page, testInfo, 'blastn', QUERY, EMPTY_RECORD);
  if (cli === undefined) return;
  expect(cli.status).toBe(3);
  expect(cli.stderr).toBe(`${shown}${ERROR}\n`);
});

// A structured -outfmt (5) writes none of the queries' warnings (NCBI's blast_format.cpp:1431-1452),
// and the run also writes outfmt 0, 6 and 7, whose command lines warn of a query without residues:
// the run's warnings are the CLI's standard error for -outfmt 5, then that warning, which the CLI's
// standard error for outfmt 6 adds (line WQA's sweep, D3 review r3 M2, 2026-11-19: it was lost).
test('a run with a structured -outfmt shows the query warnings that its other formats write', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the warnings of the engine (engine build)');
  const query = `>qempty no residues\n${QUERY}`;
  const subject = `>s1 first subject\n${dna(9, 120)}${dna(7, 300).slice(30, 270)}${dna(10, 120)}\n`;
  const warning = 'Warning: [blastn] Query_1 qempty no residues: Sequence contains no data \n';
  await program(page, 'blastn');
  await paste(page, 'query', query);
  await paste(page, 'subject', subject);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-outfmt 5');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await pickRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-diagnostics')).toBeVisible();
  const shown = await page.getByTestId('run-diagnostics').evaluate((element) => element.textContent ?? '');
  expect(shown).toBe(warning);
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The native standard error: ${NO_NATIVE_REASON}` });
    return;
  }
  const directory = mkdtempSync(join(testInfo.outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  const stderr = (outfmt: string) =>
    String(spawnSync(NATIVE!, ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-outfmt', outfmt, '-num_threads', '1'], { cwd: directory }).stderr);
  expect(stderr('5')).toBe('');
  expect(stderr('6')).toBe(warning);
});
