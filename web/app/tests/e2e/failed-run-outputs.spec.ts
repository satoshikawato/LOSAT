// A failed run keeps the outputs that the command line writes before its error, as it keeps the
// warnings (line WQA's sweep, D15, 2026-11-20; Owner 2026-11-20 17:0x, Q5=A: "a failed run keeps the
// outputs written before the error, the same bytes as the CLI, and shows them with the error's
// words"; docs/web/abi_v2.md §5). A BLASTP query of two batches whose second has no residues makes
// the engine write the first batch's hits and then fail with NCBI's `Sequence contains no data`
// (exit status 3), as the native CLI does in each format. The results page shows the output of the
// argv's own -outfmt under the failure, as the Outputs of a completed run show an output, and
// exports it; with LOSAT_WEB_NATIVE what is shown and what is exported are the bytes of the CLI's
// file for the same command line. The other outputs of the run (outfmt 0, 6 and 7 besides the
// argv's own) are not what the CLI writes for the argv, and are not offered. The FakeEngine writes a
// marked line to the output of the argv.
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Download, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openFiles, openParameters, openSearch, paste, pickRun, program, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

const FIXTURES = join(REPOSITORY, 'web/adapter/tests/fixtures/failed_outputs');
/** Two query batches: a protein of 12,000 letters that hits the subject, then a record without residues. */
const QUERY = readFileSync(join(FIXTURES, 'later_batch.q.faa'), 'utf8');
const SUBJECT = readFileSync(join(FIXTURES, 'later_batch.s.faa'), 'utf8');
const ERROR = 'BLAST engine error: Warning: Sequence contains no data';

/** The bytes of the file that `click` downloads. */
async function downloaded(page: Page, click: () => Promise<void>): Promise<Buffer> {
  const [download]: [Download, void] = await Promise.all([page.waitForEvent('download'), click()]);
  return readFileSync((await download.path())!);
}

/** The file that the native CLI writes with `-out` for the run's command line (`outfmt` in the options), and its status. */
function nativeFile(options: readonly string[]): { readonly file: Buffer; readonly status: number | null } | undefined {
  if (NATIVE === undefined) return undefined;
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), QUERY);
  writeFileSync(join(directory, 'subject.fa'), SUBJECT);
  const cli = spawnSync(NATIVE, ['blastp', '-query', 'query.fa', '-subject', 'subject.fa', ...options, '-num_threads', '1', '-out', 'out.txt'], { cwd: directory });
  return { file: readFileSync(join(directory, 'out.txt')), status: cli.status };
}

// The output of the argv's own -outfmt is the CLI's file for the argv: the hits of the first batch, and
// outfmt 0's prolog before them. Where the argv gives no -outfmt it is outfmt 0 (the page's tab
// "outfmt 0"), else 6, 7 or the output of its own (the tab "outfmt 6 qseqid sseqid").
for (const [outfmt, tab] of [
  [undefined, 'format-0'],
  ['6', 'format-6'],
  ['7', 'format-7'],
  ['6 qseqid sseqid', 'format-custom'],
] as const) {
  test(`a run that fails after its first query batch keeps what the CLI writes before the error: ${outfmt === undefined ? 'no -outfmt' : `-outfmt ${outfmt}`}`, async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, 'the engine fails the run after its first batch (engine build)');
    await program(page, 'blastp');
    await paste(page, 'query', QUERY);
    await paste(page, 'subject', SUBJECT);
    if (outfmt !== undefined) {
      await openParameters(page);
      await page.getByTestId('param-extra').fill(outfmt.includes(' ') ? `-outfmt "${outfmt}"` : `-outfmt ${outfmt}`);
    }
    await submit(page);
    await waitStatus(page, 1, 'failed');
    await pickRun(page, 1);
    const status = page.getByTestId('results-status');
    await expect(status).toHaveAttribute('data-run-status', 'failed');
    await expect(status).toContainText(`Run 1 failed: ${ERROR}.`);
    await expect(status).toContainText('It keeps the outputs written before the error.');
    await expect(status).not.toContainText('keeps no results');
    // Only the output that the CLI writes for the argv is offered, with its export.
    const outputs = page.getByTestId('results-failed-outputs');
    await expect(outputs.getByRole('heading', { name: 'Outputs' })).toBeVisible();
    await expect(outputs.getByTestId(tab)).toBeVisible();
    await expect(outputs.getByRole('tab')).toHaveCount(1);
    const shown = outputs.getByTestId('result-output');
    await expect(shown).toBeVisible();
    const text = (await shown.textContent()) ?? '';
    expect(text.length).toBeGreaterThan(0);
    const saved = await downloaded(page, () => outputs.getByTestId('export-output').click());
    expect(saved.toString('utf8')).toBe(text);
    const cli = nativeFile(outfmt === undefined ? [] : ['-outfmt', outfmt]);
    if (cli === undefined) test.info().annotations.push({ type: 'skip', description: `The bytes of the CLI's file: ${NO_NATIVE_REASON}` });
    else {
      expect(cli.status, 'the native CLI fails the same command line').toBe(3);
      expect(saved.equals(cli.file), 'the exported output is the CLI\'s file, byte for byte').toBe(true);
      expect(text).toBe(cli.file.toString('utf8'));
    }
    // The run's other runs and views are as they were: no hits, no other output, the warnings as before.
    await expect(page.getByTestId('results-hits')).toHaveCount(0);
  });
}

// A run that fails before its search wrote anything (here a subject without residues, where outfmt 6 is
// an empty file) keeps no output, and says so.
test('a failed run that wrote no output keeps none, and says so', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, 'the engine fails the run on the subject (engine build)');
  await program(page, 'blastn');
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGACTGACTAGCTAGCTAGCATGCATCGATCGATC\n');
  await paste(page, 'subject', '>empty record\n');
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-outfmt 6');
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  const status = page.getByTestId('results-status');
  await expect(status).toContainText('Run 1 failed: BLAST engine error: The average subject length is too short.');
  await expect(status).toContainText('It keeps the warnings written before the error.');
  await expect(status).not.toContainText('It keeps the outputs written before the error.');
  await expect(page.getByTestId('results-failed-outputs')).toHaveCount(0);
});

// The FakeEngine writes a marked line to the output of the failed run's argv and drops the others
// (its run of an empty input): the same page, without the engine.
test('a failed run of the FakeEngine keeps its marked output', async ({ page }) => {
  test.skip(BUILD_HAS_ENGINE, 'the FakeEngine build');
  await program(page, 'blastn');
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGACTGACTAGCTAGCTAGCATGCATCGATCGATC\n');
  await openFiles(page, 'subject', [{ name: 'empty.fa', text: '' }]);
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  await expect(page.getByTestId('results-status')).toContainText('It keeps the outputs written before the error.');
  const outputs = page.getByTestId('results-failed-outputs');
  await expect(outputs.getByTestId('format-0')).toBeVisible();
  await expect(outputs.getByRole('tab')).toHaveCount(1);
  await expect(outputs.getByTestId('result-output')).toHaveText('# FAKE ENGINE OUTPUT - not a LOSAT search result\n');
});
