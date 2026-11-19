// A query whose reader fails on a later record runs as on the command line (line WQA's D10,
// 2026-11-19; docs/web/abi_v2.md §4 `losat_web2_register`): the command line reads its queries one
// batch at a time and ends with the reader's error in the batch that reads the record, after the
// reports of the batches before, so the input is not refused before the run. Its card lists the records
// before the error and warns with the command line's message; the run ends with that message (a failed
// run keeps no results). Under -parse_deflines an error of an earlier batch (two records of one id)
// comes first, as on the command line. A query of comment lines only is a batch without a query: the
// card says that it has no sequence, and the run ends with NCBI's `Empty CBlastQueryVector`.
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openParameters, openSearch, paste, pickRun, program, submit, waitStatus } from './support/search';

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
  `BLAST query error: CFastaReader: Near line ${line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`;
const WARNING = 'A run reads the records before this error and then stops with it, as the command line does: ';

/** The native CLI's run of `words` (program first) on the inputs as query.fa and subject.fa, one thread. */
function nativeRun(words: readonly string[], query: string, subject: string): { status: number | null; stderr: string } {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  const args = [words[0]!, '-query', 'query.fa', '-subject', 'subject.fa', ...words.slice(1), '-num_threads', '1'];
  const result = spawnSync(NATIVE!, args, { cwd: directory, maxBuffer: 1 << 28 });
  return { status: result.status, stderr: result.stderr.toString('utf8') };
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
  await expect(check).toHaveText(`${WARNING}${NOT_PLAUSIBLE(6)}`);
  await expect(page.getByTestId('query-source-0-error')).toHaveCount(0);
  const status = await failedRun(page, 1);
  if (!BUILD_HAS_ENGINE) return;
  await expect(status).toContainText(NOT_PLAUSIBLE(6));
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastn', '-outfmt', '6'], query, SUBJECT);
  expect(cli.status).toBe(1);
  expect(cli.stderr).toBe(`${NOT_PLAUSIBLE(6)}\n`);
});

// BLASTP reads 10000 letters a batch: ten records of 1100 fill the first. Two records named q1 in the
// first batch end the command line there under -parse_deflines, before it reads the record of the
// second batch whose Seq-id names a nucleotide sequence (`ref|NM_…`); the Web read every record before
// the run and failed with the second error (review of line WQA's D7, L1).
test('under -parse_deflines an error of an earlier query batch comes before the reader\'s error of a later one', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine reads no ids of the deflines: the engine build runs this test');
  const protein = letters('ACDEFGHIKLMNPQRSTVWY', 11, 1100);
  const names = ['q1', 'q1', ...Array.from({ length: 9 }, (_, i) => `q${i + 3}`), 'ref|NM_000001.1|'];
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
  await expect(status).toContainText('BLAST engine error: Empty CBlastQueryVector');
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastn', '-outfmt', '6'], query, SUBJECT);
  expect(cli.status).toBe(3);
  expect(cli.stderr).toBe('BLAST engine error: Empty CBlastQueryVector\n');
});
