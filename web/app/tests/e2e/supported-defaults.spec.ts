// The mark of a run with the defaults (WA-RS-03): in the engine build, a run of every mode and
// task with its default options is Supported, with one thread and with Auto. The defaults'
// option set has NCBI's records in outfmt 0, 6 and 7 (build/verification.ts; the def7.* rows of
// LOSAT/tests/outfmt0_manifest.tsv), and Auto's thread count, half of the processors up to 32,
// is at most the largest count that V-BR checked in the browsers (engine.spec.ts,
// domain/verification.ts). The inputs are larger than AUTO_SERIAL_BELOW_BYTES, so Auto runs the
// threaded module wherever the browser reports four processors or more.
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, fasta, openFiles, openSearch, program, submit, task, waitStatus } from './support/search';

test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine build marks every run "Development build"');
// Engine searches are slower in Firefox and WebKit than in Chromium (results.spec.ts).
test.setTimeout(600_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

interface Mode {
  readonly id: 'blastn' | 'blastp' | 'blastx' | 'tblastn' | 'tblastx';
  /** The tasks of its Task Selection, the default first (none for TBLASTX). */
  readonly tasks: readonly string[];
  /** Files of LOSAT/tests/fasta; together more than 20,000 bytes. */
  readonly query: string;
  readonly subject: readonly string[];
}

const MODES: readonly Mode[] = [
  { id: 'blastn', tasks: ['megablast', 'dc-megablast', 'blastn', 'blastn-short'], query: 'outfmt0/multi_query.fasta', subject: ['LC738870.fasta'] },
  { id: 'blastp', tasks: ['blastp', 'blastp-fast'], query: 'outfmt0/e2e_protein_query.faa', subject: ['HtNMV.faa'] },
  { id: 'blastx', tasks: ['blastx', 'blastx-fast'], query: 'outfmt0/pd_sub_q.fna', subject: ['HtNMV.faa'] },
  { id: 'tblastn', tasks: ['tblastn', 'tblastn-fast'], query: 'outfmt0/e2e_protein_query.faa', subject: ['LC738870.fasta'] },
  { id: 'tblastx', tasks: [], query: 'outfmt0/multi_query.fasta', subject: ['LC738870.fasta'] },
];

/** Runs the search as the form is, and returns the mark of the run with its reasons and the threads it used. */
async function runAndMark(page: Page, run: number): Promise<{ mark: string; why: string; engine: string }> {
  await submit(page);
  await waitStatus(page, run, 'completed');
  await chooseRun(page, run);
  const badge = page.getByTestId('results-summary').getByTestId('verification-badge');
  const mark = ((await badge.textContent()) ?? '').trim();
  const why = (await badge.getAttribute('title')) ?? '';
  await page.getByTestId('results-view-details').click();
  const engine = ((await page.getByTestId('run-record').locator('[data-detail="path"]').textContent()) ?? '').trim();
  await page.getByTestId('tab-search').click();
  return { mark, why, engine };
}

for (const mode of MODES) {
  test(`${mode.id}: a run with the defaults of each task is Supported with 1 thread and with Auto`, async ({ page }) => {
    await program(page, mode.id);
    await openFiles(page, 'query', [{ name: mode.query.split('/').at(-1)!, text: fasta(mode.query) }]);
    await openFiles(page, 'subject', mode.subject.map((path) => ({ name: path, text: fasta(path) })));
    let run = 0;
    for (const value of mode.tasks.length > 0 ? mode.tasks : [undefined]) {
      if (value !== undefined) await task(page, value);
      for (const threads of ['1', 'auto']) {
        await page.getByTestId('threads').selectOption(threads);
        const { mark, why, engine } = await runAndMark(page, ++run);
        const name = `${mode.id}${value === undefined ? '' : ` -task ${value}`}, Threads ${threads} (${engine})`;
        expect(mark, `${name}: ${why}`).toBe('Supported');
        if (threads === '1') expect(engine, name).toBe('serial, 1 thread');
        else if ((await page.evaluate(() => navigator.hardwareConcurrency)) >= 4) expect(engine, name).toMatch(/^threaded, \d+ threads$/);
        test.info().annotations.push({ type: 'mark', description: `${name}: ${mark}` });
      }
    }
  });
}

test('blastn: Auto with an input of 1 MB or more (up to 32 threads) is Supported', async ({ page }) => {
  // Four genomes of about 300 kB: Auto's larger limit (AUTO_MANY_THREADS_FROM_BYTES).
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'multi_query.fasta', text: fasta('outfmt0/multi_query.fasta') }]);
  const genomes = ['LC738868.fasta', 'LC738869.fasta', 'LC738870.fasta', 'LC738871.fasta'];
  await openFiles(page, 'subject', genomes.map((path) => ({ name: path, text: fasta(path) })));
  await page.getByTestId('threads').selectOption('auto');
  const { mark, why, engine } = await runAndMark(page, 1);
  expect(mark, `${engine}: ${why}`).toBe('Supported');
  test.info().annotations.push({ type: 'mark', description: `blastn, 4 genomes, Threads auto (${engine}): ${mark}` });
});
