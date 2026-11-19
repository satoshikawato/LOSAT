// A failed run keeps the warnings that the command line writes to its standard error before its
// error, and the results page shows them under the failure, as the engine wrote them (line WQA's
// sweep, D3, 2026-11-19; docs/web/abi_v2.md §5, stream 3). In the engine build, a BLASTN subject
// whose one record has no residues gives NCBI's `Subject sequence contains no data` warning and then
// its error `The average subject length is too short` (exit status 3); with LOSAT_WEB_NATIVE the
// warnings shown and the error are the native CLI's standard error, line for line. The FakeEngine
// fails a run of an empty input with a marked warning. A run of two searches (ABI §7) that its
// format 0 search fails shows that search's warnings too (D3 review, 2026-11-19).
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect } from '@playwright/test';
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
  await expect(status).toContainText('A failed run keeps no results. It keeps the warnings written before the error.');
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

  // The command line's standard error: the warnings shown, then the error.
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The native standard error: ${NO_NATIVE_REASON}` });
  } else {
    const directory = mkdtempSync(join(testInfo.outputPath(), 'native-'));
    writeFileSync(join(directory, 'query.fa'), QUERY);
    writeFileSync(join(directory, 'subject.fa'), EMPTY_RECORD);
    const cli = spawnSync(NATIVE, ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-outfmt', '6', '-num_threads', '1'], { cwd: directory });
    expect(cli.status).toBe(3);
    expect(String(cli.stderr)).toBe(`${shown}${ERROR}\n`);
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

// BLASTP with -num_alignments 0 makes two searches (docs/web/abi_v2.md §7): the tabular one fails
// alone before it reads the subject (the CLI's standard error for -outfmt 6: the handler's warning and
// NCBI's error), and the format 0 search, the argv's own -outfmt 0, warns of the subject's empty first
// record and fails the run on the title of its second (`, ,`, where NCBI BLAST+ 2.17.0 crashes and
// LOSAT refuses). The warnings shown are the standard error for -outfmt 6 followed by what that for
// -outfmt 0 writes before its error.
test('a run that its format 0 search fails shows that search\'s warnings after those of the tabular formats', async ({ page }, testInfo) => {
  test.skip(!BUILD_HAS_ENGINE, 'the engine fails the run on the title (engine build)');
  const query = '>p\nMKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMAWRSTY\n';
  const subject = '>empty record\n>, ,\nMKVLAAGIVGLLLAAQPAMAMKVLAAGIVGLLLAAQPAMA\n';
  await program(page, 'blastp');
  await openParameters(page);
  await page.getByTestId('param-num_alignments').fill('0');
  await paste(page, 'query', query);
  await paste(page, 'subject', subject);
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  const status = page.getByTestId('results-status');
  await expect(status).toContainText('subject record 2 has a title of punctuation');
  await expect(status).toContainText('It keeps the warnings written before the error.');
  const shown = await page.getByTestId('results-failed-diagnostics').evaluate((element) => element.textContent ?? '');
  expect(shown).toBe(
    'Warning: [blastp] Examining 5 or more matches is recommended\nBLAST query/options error: No hits are being saved\nPlease refer to the BLAST+ user manual.\nWarning: [blastp] Subject_1 empty record: Subject sequence contains no data\n',
  );
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The native standard error: ${NO_NATIVE_REASON}` });
    return;
  }
  const directory = mkdtempSync(join(testInfo.outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  const cli = (outfmt: string) =>
    spawnSync(NATIVE!, ['blastp', '-query', 'query.fa', '-subject', 'subject.fa', '-num_alignments', '0', '-outfmt', outfmt, '-num_threads', '1'], {
      cwd: directory,
    });
  const tabular = cli('6');
  const own = cli('0');
  expect(tabular.status).toBe(1);
  expect(own.status).toBe(1);
  // The argv's own -outfmt 0: its warning, then the error that ends the run.
  const ownLines = String(own.stderr).split(/(?<=\n)/);
  expect(ownLines.at(-1)).toMatch(/^Error: subject record 2 has a title of punctuation/);
  expect(shown).toBe(`${String(tabular.stderr)}${ownLines.slice(0, -1).join('')}`);
});
