// A query whose reader fails on a later record runs as on the command line (line WQA's D10,
// 2026-11-19; docs/web/abi_v2.md §4 `losat_web2_register`): the command line reads its queries one
// batch at a time and ends with the reader's error in the batch that reads the record, after the
// reports of the batches before, so the input is not refused before the run. Its card lists the records
// before the error and warns with the command line's message; a run that reads that far ends with that
// message and keeps the outputs written before it, the reports of the batches before (line WQA's D28,
// 2026-11-27): the records of the error's batch have none, so with the error in the first batch no record has
// results (PSI-BLAST searches one query at a time: the records before the error keep theirs), outfmt 0 keeps its prolog
// (-html its page's head, outfmt 15 the start of its JSON), and outfmt 5, 6 and 7 are empty. Under -parse_deflines an
// error of an earlier batch (two records of one id) comes first, as on the command line. A query of comment lines only is a batch without a query: the
// card says that it has no sequence, and the run ends with NCBI's `Empty CBlastQueryVector`.
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openFiles, openParameters, openSearch, paste, pickRun, program, settled, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** Pseudo-random letters of `alphabet` (a linear congruential generator): the same letters in every run. */
function letters(alphabet: string, seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += alphabet[(state >>> 16) % alphabet.length];
  }
  return out;
}

const DNA = letters('ACGT', 7, 300);
const SUBJECT = `>s1 subject\n${DNA}\n`;
const NOT_PLAUSIBLE = (line: number) =>
  `LOSAT query error: CFastaReader: Near line ${line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`;
/**
 * The warning (review M-a, 2026-11-20; line WQA's D28, 2026-11-27): a run that reaches the error fails there; an error before
 * it (an ID repeated in an earlier batch under -parse_deflines) ends the run first.
 */
const WARNING = 'A run that reaches this error fails there, as the command line does: ';
/**
 * What the warning says that the run keeps (line WQA's D28 and its reviews r1 L1 and r2 M1, 2026-11-27: the reports of the
 * query batches before the error's; the records of the error's batch have none, one record in PSI-BLAST) and what to do (review L2), with the input's own line where the
 * run's differs (review L-b, 2026-11-20).
 */
const REMEDY = (own?: string) =>
  ` The records searched before it keep their results; the records read together with it have none (PSI-BLAST reads one record at a time). Correct that line${own === undefined ? '' : ` (${own})`} to search this input.`;

/** The native CLI's run of `words` (program first) on the inputs as query.fa and subject.fa, one thread. */
function nativeRun(words: readonly string[], query: string, subject: string): { status: number | null; stdout: string; stderr: string } {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  const args = [words[0]!, '-query', 'query.fa', '-subject', 'subject.fa', ...words.slice(1), '-num_threads', '1'];
  const result = spawnSync(NATIVE!, args, { cwd: directory, maxBuffer: 1 << 28 });
  return { status: result.status, stdout: result.stdout.toString('utf8'), stderr: result.stderr.toString('utf8') };
}

/** Runs the form and returns the results' status line of run `run`, which failed. */
async function failedRun(page: Page, run: number): Promise<ReturnType<Page['getByTestId']>> {
  await submit(page);
  await waitStatus(page, run, 'failed');
  await pickRun(page, run);
  const status = page.getByTestId('results-status');
  await expect(status).toContainText(`Run ${run} failed: `);
  return status;
}

test("a query whose reader fails on a later record runs, and the run ends with the reader's error", async ({ page }, testInfo) => {
  const query = `>q1 first\n${DNA.slice(0, 120)}\n>q2 second\n${DNA.slice(100, 260)}\n>bad\nACGTAC%%%%\n>q4\nACGT\n`;
  await program(page, 'blastn');
  await paste(page, 'subject', SUBJECT);
  await paste(page, 'query', query);
  await expect(page.getByTestId('query-source-0-summary')).toHaveText('2 records · 280 nt');
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-check', 'ok');
  await expect(check).toHaveAttribute('data-reader-error', '');
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(6)}${REMEDY()}`);
  await expect(page.getByTestId('query-source-0-error')).toHaveCount(0);
  const status = await failedRun(page, 1);
  if (!BUILD_HAS_ENGINE) return;
  await expect(status).toContainText(NOT_PLAUSIBLE(6));
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  // The error's batch is the first, so the command line writes no row either.
  const cli = nativeRun(['blastn', '-outfmt', '6'], query, SUBJECT);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`${NOT_PLAUSIBLE(6)}\n`);
  expect(cli.stdout).toBe('');
});

// BLASTP reads 10000 letters a batch: ten records of 1100 fill the first. Two records named q1 in the
// first batch end the command line there under -parse_deflines, before it reads the record of the
// second batch whose Seq-id is a database-ordinal id (`gnl|BL_ORD_ID|1`) that the reader refuses; the
// Web read every record before the run and failed with the second error (review of line WQA's D7, L1).
test('under -parse_deflines an error of an earlier query batch comes before the reader\'s error of a later one', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine reads no ids of the deflines: the engine build runs this test');
  const protein = letters('ACDEFGHIKLMNPQRSTVWY', 11, 1100);
  const names = ['q1', 'q1', ...Array.from({ length: 9 }, (_, i) => `q${i + 3}`), 'gnl|BL_ORD_ID|1'];
  const query = names.map((name) => `>${name}\n${protein}\n`).join('');
  const subject = `>s1\n${protein.slice(0, 300)}\n`;
  const message =
    "two query records have the Seq-id lcl|q1 under -parse_deflines in one query batch, or in any with -subject or outfmt 5 (NCBI BLAST+'s scope holds both and searches the letters of the first for both); this is not supported by LOSAT's BLASTP";
  await program(page, 'blastp');
  await paste(page, 'subject', subject);
  await paste(page, 'query', query);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-parse_deflines');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  // The validation of -parse_deflines comes after a pause in typing (the state above is still the earlier argv's), and
  // the engine then reads both inputs again; the query's check under it finds the reader's error of the last record,
  // which the card shows (it does not without the option). Run LOSAT is pressed once the page shows that, as a user
  // sees it: a press while the card's warning came was a press while the button moved (line WQA's D25, 2026-11-27).
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-reader-error', '', { timeout: 30_000 });
  // The card's warning holds only if the run reads that far: it does not say that the run ends with it (line WQA's D28).
  await expect(page.getByTestId('query-source-0-check')).toContainText(WARNING);
  await settled(page, 'subject');
  const status = await failedRun(page, 1);
  await expect(status).toContainText(message);
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastp', '-parse_deflines', '-outfmt', '6'], query, subject);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`Error: ${message}\n`);
});

// A reader's error in a later query batch (line WQA's D28, 2026-11-27): BLASTP reads 10000 letters a batch, so ten records
// of 1100 fill the first and the record that the reader refuses is in the second. The command line writes the first
// batch's hits and then fails with the reader's error; the run does the same and keeps those outputs, as its card says.
test("a query whose reader fails in a later batch: the run fails with the reader's error and keeps the outputs written before it", async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine writes the first batch's hits: the engine build runs this test");
  const protein = letters('ACDEFGHIKLMNPQRSTVWY', 11, 1100);
  const query = `${Array.from({ length: 10 }, (_, i) => `>q${i + 1}\n${protein}\n`).join('')}>bad\nACDEF%%%%\n>q12\nACDEFG\n`;
  const subject = `>s1\n${protein.slice(0, 300)}\n`;
  await program(page, 'blastp');
  await paste(page, 'subject', subject);
  await paste(page, 'query', query);
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-reader-error', '');
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(22)}${REMEDY()}`);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-outfmt 6');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  const status = await failedRun(page, 1);
  await expect(status).toContainText(NOT_PLAUSIBLE(22));
  await expect(status).toContainText('It keeps the outputs written before the error.');
  const shown = page.getByTestId('results-failed-outputs').getByTestId('result-output');
  await expect(shown).toContainText('q1\ts1\t');
  await expect(shown).toContainText('q10\ts1\t');
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's output: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastp', '-outfmt', '6'], query, subject);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`${NOT_PLAUSIBLE(22)}\n`);
  expect((await shown.textContent()) ?? '').toBe(cli.stdout);
});

// PSI-BLAST searches one query at a time (line WQA's D28, review r2 M1, 2026-11-27): its first batch is one record, so the
// records before the error keep their results, as the command line writes them, where BLASTP's short query keeps none.
test("a PSI-BLAST query whose reader fails on a later record: the run keeps the hits of the records before it", async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine writes the records' hits: the engine build runs this test");
  const protein = letters('ACDEFGHIKLMNPQRSTVWY', 11, 200);
  const query = `${Array.from({ length: 3 }, (_, i) => `>q${i + 1}\n${protein}\n`).join('')}>bad\nACDEF%%%%\n`;
  const subject = `>s1\n${protein}\n`;
  await program(page, 'psiblast');
  await paste(page, 'subject', subject);
  await paste(page, 'query', query);
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-reader-error', '');
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(8)}${REMEDY()}`);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-outfmt 6');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  const status = await failedRun(page, 1);
  await expect(status).toContainText(NOT_PLAUSIBLE(8));
  await expect(status).toContainText('It keeps the outputs written before the error.');
  const shown = page.getByTestId('results-failed-outputs').getByTestId('result-output');
  await expect(shown).toContainText('q1\ts1\t');
  await expect(shown).toContainText('q3\ts1\t');
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's output: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['psiblast', '-outfmt', '6'], query, subject);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`${NOT_PLAUSIBLE(8)}\n`);
  expect(cli.stdout).toContain('q3\ts1\t');
  expect((await shown.textContent()) ?? '').toBe(cli.stdout);
});

// One Combined query input reads its sources in turn (review of D10, 2026-11-20): the warning names the
// run's line, after the input before it (M2/L2), and the input's own line (L-b); the input after it is
// not read, and its card says so in place of "Accepted"; a source whose error is on a first line without
// a defline (a byte order mark) cannot follow another input, whose last record would read that line (M2):
// its card says only that, and the card after it does not say that the search stops there (L-c).
// Separate runs each alone.
test("a Combined query input reads up to a source's reader's error, and the warning names the run's line", async ({ page }, testInfo) => {
  const pasted = `>q1 first\n${DNA.slice(0, 120)}\n`;
  const bad = `>a\n${DNA.slice(0, 90)}\n>b\nACGTAC%%%%\n`;
  const after = `>c\n${DNA.slice(150, 240)}\n`;
  await program(page, 'blastn');
  await paste(page, 'subject', SUBJECT);
  await paste(page, 'query', pasted);
  await openFiles(page, 'query', [
    { name: 'bad.fa', text: bad },
    { name: 'after.fa', text: after },
  ]);
  const check = page.getByTestId('query-source-1-check');
  await expect(check).toHaveAttribute('data-reader-error', '');
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(6)}${REMEDY('line 4 of bad.fa')} The run does not read the inputs after this one.`);
  const later = page.getByTestId('query-source-2-check');
  await expect(later).toHaveAttribute('data-not-searched', '');
  await expect(later).toHaveText(
    'Not searched: a run does not read past the error in bad.fa, before this input. Correct that error, or choose Separate to search each input alone.',
  );
  await page.getByTestId('query-mode-separate').check();
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(4)}${REMEDY()}`);
  await expect(later).not.toHaveAttribute('data-not-searched', '');
  await page.getByTestId('query-mode-combined').check();
  await expect(check).toContainText(NOT_PLAUSIBLE(6));
  const status = await failedRun(page, 1);
  if (BUILD_HAS_ENGINE) await expect(status).toContainText(NOT_PLAUSIBLE(6));
  // A source whose reader's error is on its first line, read without a defline.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('query-source-2-remove').click();
  await page.getByTestId('query-source-1-remove').click();
  await openFiles(page, 'query', [
    { name: 'mark.fa', text: Buffer.from(`\uFEFF>q2\n${DNA.slice(0, 60)}\n`) },
    { name: 'after.fa', text: after },
  ]);
  await expect(page.getByTestId('query-source-1-mark')).toBeVisible();
  const mark = page.getByTestId('query-source-1-check');
  await expect(mark).toHaveAttribute('data-cannot-follow', '');
  await expect(mark).toHaveText(
    'This input starts with a line that is not a defline (">"), which cannot follow another input. Correct that line or remove this input.',
  );
  await expect(later).toHaveAttribute('data-check', 'ok');
  await expect(later).not.toHaveAttribute('data-not-searched', '');
  await expect(later).not.toHaveAttribute('data-cannot-follow', '');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText(
    'Query (file mark.fa) starts with a line that is not a defline (">"), which cannot follow another input. Correct that line or remove this input.',
  );
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastn', '-outfmt', '6'], `${pasted}${bad}${after}`, SUBJECT);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`${NOT_PLAUSIBLE(6)}\n`);
});

// The run passes every input of a Combined query as it is, one of comment and blank lines only too: the
// warning names the run's line, which counts that input's lines, and this input's own line (review L-a and
// L-b, 2026-11-20).
test("a Combined query's warning counts the lines of an input without records before it", async ({ page }, testInfo) => {
  const pasted = `>q\n${DNA.slice(0, 60)}\n`;
  const notes = '; notes\n\n# more\n';
  const bad = `>a\n${DNA.slice(0, 90)}\n>b\nACGTAC%%%%\n`;
  await program(page, 'blastn');
  await paste(page, 'subject', SUBJECT);
  await paste(page, 'query', pasted);
  await openFiles(page, 'query', [
    { name: 'notes.fa', text: notes },
    { name: 'bad.fa', text: bad },
  ]);
  const check = page.getByTestId('query-source-2-check');
  await expect(check).toHaveAttribute('data-reader-error', '');
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(9)}${REMEDY('line 4 of bad.fa')}`);
  const status = await failedRun(page, 1);
  if (BUILD_HAS_ENGINE) await expect(status).toContainText(NOT_PLAUSIBLE(9));
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastn', '-outfmt', '6'], `${pasted}${notes}${bad}`, SUBJECT);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`${NOT_PLAUSIBLE(9)}\n`);
});

test('a query of comment lines only has no sequence, and the run ends as the command line ends it', async ({ page }, testInfo) => {
  const query = '; a comment\n\n# another\n';
  await program(page, 'blastn');
  await paste(page, 'subject', SUBJECT);
  await paste(page, 'query', query);
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-check', 'ok');
  await expect(check).toHaveAttribute('data-sequence', 'none');
  if (!BUILD_HAS_ENGINE) return;
  const status = await failedRun(page, 1);
  await expect(status).toContainText('LOSAT engine error: Empty CBlastQueryVector');
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastn', '-outfmt', '6'], query, SUBJECT);
  expect(cli.status).toBe(3);
  expect(cli.stderr).toBe('LOSAT engine error: Empty CBlastQueryVector\n');
});
