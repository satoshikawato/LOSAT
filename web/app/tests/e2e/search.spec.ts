// The search screen (S12, W3; plan §7): research work and boundaries through the real
// application. The tests run with the FakeEngine build and with the engine build
// (LOSAT_WEB_REACTORS); what only the engine can show (its messages, the kept subject, a
// search long enough to edit or cancel while it runs) is checked in the engine build.
import { expect, test, type Page } from '@playwright/test';
import { gzipSync, zstdCompressSync } from 'node:zlib';
import { BUILD_HAS_ENGINE } from './support/browser';
import {
  chooseRun,
  fasta,
  openFiles,
  openParameters,
  openSearch,
  paste,
  pickRun,
  program,
  result,
  runMessage,
  settled,
  showRecords,
  submit,
  task,
  waitStatus,
} from './support/search';

/**
 * A BLASTP search of a proteome against itself (about 500,000 residues): many seconds in
 * every browser (S09's "Auto" measurements), long enough to edit or cancel during it.
 */
const SLOW = { program: 'blastp', query: 'NZ_CP006932.faa', subject: 'NZ_CP006932.faa' } as const;

async function startSlowSearch(page: Page): Promise<void> {
  await program(page, SLOW.program);
  await openFiles(page, 'query', [{ name: SLOW.query, text: fasta(SLOW.query) }]);
  await openFiles(page, 'subject', [{ name: SLOW.subject, text: fasta(SLOW.subject) }]);
  await submit(page);
  await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
}

/** Back to the search form (it stays mounted with its edits) after "Run LOSAT" went to the results. */
async function backToSearch(page: Page): Promise<void> {
  await page.getByTestId('tab-search').click();
}

// Engine searches are slower in Firefox and WebKit than in Chromium (W1 README).
test.setTimeout(300_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

test('a file is summarized with its records, not put in the paste box; files can also be dropped', async ({ page }) => {
  const records = Array.from({ length: 300 }, (_, i) => `>rec${i + 1} record number ${i + 1}\n${'ACGT'.repeat(25)}\n`).join('');
  await openFiles(page, 'query', [{ name: 'many.fa', text: records }]);
  await expect(page.getByTestId('query-input')).toHaveValue('');
  await expect(page.getByTestId('query-source-0-summary')).toHaveText('300 records · 30,000 nt');
  await expect(page.getByTestId('query-source-0')).toContainText('many.fa');
  await page.getByTestId('query-source-0').locator('summary', { hasText: 'First lines' }).click();
  await expect(page.getByTestId('query-source-0-head')).toContainText('>rec1 record number 1');
  // The record list draws only the rows in view.
  await showRecords(page, 'query');
  const drawn = await page.locator('[data-testid^="query-source-0-record-"]').count();
  expect(drawn).toBeGreaterThan(5);
  expect(drawn).toBeLessThan(40);
  await page.getByTestId('query-source-0-filter').fill('rec29');
  await expect(page.locator('[data-testid^="query-source-0-record-"]')).toHaveCount(11);

  // Drag and drop of a file onto the subject.
  const dataTransfer = await page.evaluateHandle(() => {
    const transfer = new DataTransfer();
    transfer.items.add(new File(['>dropped\nACGTACGTAC\n'], 'dropped.fa', { type: 'text/plain' }));
    return transfer;
  });
  await page.getByTestId('subject-dropzone').dispatchEvent('drop', { dataTransfer });
  await settled(page, 'subject');
  await expect(page.getByTestId('subject-source-0')).toContainText('dropped.fa');
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('1 record · 10 nt');
});

// A record that the checker refuses. The FakeEngine refuses a record with '!' in its sequence and
// names it in its message ("query record 3 (bad) ..."), so the screen offers to exclude it. The
// engine's register reads with the index scan's reader kind (DW-38), so it refuses no record that
// the index read; NCBI's reader skips the '!', and the record is left out with its checkbox.
const REFUSED_RECORD = '>bad\nACGT!ACGT\n';

test('a file named .gz is inflated when it is opened and shows its plain name', async ({ page }) => {
  await openFiles(page, 'subject', [{ name: 'LC738884.fasta.gz', text: gzipSync(fasta('LC738884.fasta')) }]);
  await expect(page.getByTestId('subject-source-0')).toHaveAttribute('data-status', 'ready');
  await expect(page.getByTestId('subject-source-0').locator('.source-name')).toHaveText('LC738884.fasta');
  await openFiles(page, 'subject', [{ name: 'LC738884.fasta', text: fasta('LC738884.fasta') }]);
  const summary = await page.getByTestId('subject-source-0-summary').textContent();
  expect(summary).toMatch(/record/);
  await expect(page.getByTestId('subject-source-1-summary')).toHaveText(summary!.trim());
});

test('a file named .zst is inflated when it is opened and shows its plain name', async ({ page }) => {
  await openFiles(page, 'subject', [{ name: 'LC738884.fasta.zst', text: zstdCompressSync(fasta('LC738884.fasta')) }]);
  await expect(page.getByTestId('subject-source-0')).toHaveAttribute('data-status', 'ready');
  await expect(page.getByTestId('subject-source-0').locator('.source-name')).toHaveText('LC738884.fasta');
  await openFiles(page, 'subject', [{ name: 'LC738884.fasta', text: fasta('LC738884.fasta') }]);
  const summary = await page.getByTestId('subject-source-0-summary').textContent();
  expect(summary).toMatch(/record/);
  await expect(page.getByTestId('subject-source-1-summary')).toHaveText(summary!.trim());
});

test('a file named .bz2 is refused when it is opened and shows the error on its card', async ({ page }) => {
  await openFiles(page, 'subject', [{ name: 'two.fa.bz2', text: '>a\nACGT\n' }]);
  await expect(page.getByTestId('subject-source-0')).toHaveAttribute('data-status', 'failed');
  await expect(page.getByTestId('subject-source-0-error')).toHaveText(
    "This input cannot be read: 'two.fa.bz2' is named .bz2; the web app reads .gz and .zst, not bzip2 (the command line does)",
  );
});

test('records: duplicate IDs by number, a record that the input check refuses (FakeEngine), exclusion, and a run after exclusion', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', `>dup first\nACGTACGTACGTAAACCCGGGTTT\n>dup second\nACGTACGTACGTAAACCCGGGTTA\n${REFUSED_RECORD}`);
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await expect(page.getByTestId('query-source-0-duplicates')).toContainText('2 records share an ID');
  const check = page.getByTestId('query-source-0-check');
  if (BUILD_HAS_ENGINE) {
    await expect(check).toHaveAttribute('data-check', 'ok');
    await showRecords(page, 'query');
    await page.getByTestId('query-source-0-record-2').uncheck();
  } else {
    await expect(check).toHaveAttribute('data-check', 'refused');
    await expect(check).toContainText('query record 3 (bad)');
    await expect(check).toContainText('FAKE ENGINE check');
    // The refused input cannot be started.
    await submit(page);
    await expect(page.getByTestId('search-message')).toContainText('Query (pasted): query record 3 (bad)');
    await expect(page.getByTestId('run-1-option')).toHaveCount(0);
    await page.getByTestId('query-source-0-exclude-refused').click();
  }
  await expect(check).toHaveAttribute('data-check', 'ok');
  await expect(page.getByTestId('query-source-0-summary')).toContainText('(2 included)');
  await submit(page);
  await waitStatus(page, 1, 'completed');

  // Leave out the first "dup" (by its number) and search again.
  await backToSearch(page);
  await showRecords(page, 'query');
  await page.getByTestId('query-source-0-record-0').uncheck();
  await expect(check).toHaveAttribute('data-check', 'ok');
  await expect(page.getByTestId('query-source-0-summary')).toContainText('(1 included)');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  if (BUILD_HAS_ENGINE) {
    const first = await result(page, 1, 7);
    const second = await result(page, 2, 7);
    expect(first.output.match(/^# Query: dup/gm)).toHaveLength(2);
    expect(second.output.match(/^# Query: dup second$/gm)).toHaveLength(1);
    expect(second.output).not.toContain('# Query: dup first');
  }
});

test('the same subject, another query: the engine keeps the subject (R1)', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await paste(page, 'query', '>qa\nACGTACGTACGTAAACCCGGGTTT\n');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await backToSearch(page);
  await paste(page, 'query', '>qb\nACGTACGTACGTAAACCCGGGTTA\n');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  if (BUILD_HAS_ENGINE) {
    await chooseRun(page, 2);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-record')).toContainText('kept from the previous search');
    expect((await result(page, 2, 7)).output).toContain('# Query: qb');
  }
});

test('several runs; the next search is edited while one runs; separate searches as a group', async ({ page }) => {
  if (BUILD_HAS_ENGINE) {
    await startSlowSearch(page);
  } else {
    await program(page, 'blastp');
    await openFiles(page, 'query', [{ name: 'q1.fa', text: '>q1\nMKLVVLAAGGHHKL\n' }]);
    await openFiles(page, 'subject', [{ name: 's1.fa', text: '>s1\nMKLVVLAAGGHHKL\n' }]);
    await submit(page);
    await expect(page.getByTestId('run-1-option')).toHaveCount(1);
  }

  // Edit the next search while run 1 runs: another program, query and option, two subjects searched separately.
  await backToSearch(page);
  await program(page, 'tblastx');
  await page.getByTestId('query-source-0-remove').click();
  await paste(page, 'query', '>next\nACGTACGTACGTACGTAC\n');
  await openParameters(page);
  await page.getByTestId('param-evalue').fill('1e-3');
  await page.getByTestId('subject-source-0-remove').click();
  await openFiles(page, 'subject', [
    { name: 'a.fa', text: '>a\nACGTACGTACGTACGTACGTAAA\n' },
    { name: 'b.fa', text: '>b\nTTTACGTACGTACGTACGTACGT\n' },
  ]);
  await page.getByTestId('subject-mode-separate').check();
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT (2 runs)');
  await expect(page.getByTestId('run-losat-bottom')).toHaveText('Run LOSAT (2 runs)');
  await submit(page);
  // The page goes to the results of the first run of the group.
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('run-3-option')).toHaveCount(1);
  await expect(page.getByTestId('run-2-option')).toContainText('query.fa vs a.fa');
  await expect(page.getByTestId('run-3-option')).toContainText('query.fa vs b.fa');
  await expect(page.getByTestId('results-run-options')).toContainText('group run 1 of 2');
  await pickRun(page, 3);
  await expect(page.getByTestId('results-run-options')).toContainText('-evalue 1e-3');
  await expect(page.getByTestId('results-run-options')).toContainText('group run 2 of 2');
  // Run 1 kept its own program, inputs and options.
  await expect(page.getByTestId('run-1-option')).toContainText('BLASTP');
  await pickRun(page, 1);
  await expect(page.getByTestId('results-run-options')).toHaveText('defaults');

  if (BUILD_HAS_ENGINE) {
    await expect(page.getByTestId('run-1-option')).toContainText(`${SLOW.query} vs ${SLOW.subject}`);
    // The group waits behind run 1; cancelling the group leaves run 1 running.
    await pickRun(page, 2);
    await expect(page.getByTestId('run-2-option')).toHaveAttribute('data-status', 'queued');
    await expect(page.getByTestId('run-status')).toHaveAttribute('data-status', 'queued');
    await expect(page.getByTestId('run-status')).toHaveAttribute('data-run', '2');
    await expect(page.getByTestId('run-status-phase')).toHaveText('Waiting for run 1');
    await expect(page.getByTestId('run-status-tab-warning')).toHaveText('Closing or reloading this tab ends the search');
    await page.getByTestId('run-status-cancel-group').click();
    await expect(page.getByTestId('run-2-option')).toHaveAttribute('data-status', 'cancelled');
    await expect(page.getByTestId('run-3-option')).toHaveAttribute('data-status', 'cancelled');
    await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', 'running');
    await pickRun(page, 1);
    await expect(page.getByTestId('run-status-phase')).toHaveText('Searching');
    await expect(page.getByTestId('run-status-phase').getByTestId('losat-loader')).toBeVisible(); // the loader before the words (DW-43)
    await expect(page.getByTestId('run-status-elapsed')).toHaveText(/^\d+:\d\d$/);
    await backToSearch(page);
    await submit(page); // the same group again, behind run 1
    await expect(page.getByTestId('run-5-option')).toHaveAttribute('data-status', 'queued');
    await pickRun(page, 1);
    await page.getByTestId('run-status-cancel').click();
    await waitStatus(page, 1, 'cancelled');
    // A run that ended leaves the status block, and its loader, for the results message.
    await expect(page.getByTestId('results-status')).toHaveAttribute('data-run-status', 'cancelled');
    await expect(page.getByTestId('run-status')).toHaveCount(0);
    await waitStatus(page, 4, 'completed');
    await waitStatus(page, 5, 'completed');
  } else {
    await expect(page.getByTestId('run-1-option')).toContainText('q1.fa vs s1.fa');
    await waitStatus(page, 1, 'completed');
    await waitStatus(page, 2, 'completed');
    await waitStatus(page, 3, 'completed');
  }
  const last = await result(page, BUILD_HAS_ENGINE ? 5 : 3, 6);
  expect(last.command).toBe('losat tblastx -query query.fa -subject b.fa -evalue 1e-3 -outfmt 6');
});

test('cancel a running search; the next one runs', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, 'needs a search that runs long enough to cancel (engine build)');
  await startSlowSearch(page);
  await page.getByTestId('run-status-cancel').click();
  await waitStatus(page, 1, 'cancelled');
  // The same subject, another query: the new runtime registers the subject again.
  await backToSearch(page);
  await page.getByTestId('query-source-0-remove').click();
  await paste(page, 'query', '>small\nMKLVVLAAGGHHKLMKLVVLAAGGHHKL\n');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await chooseRun(page, 2);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-record')).toContainText('read for this search');
});

test('the region of a role with one record: the drag, the fields, the limits and the argv', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', '>q1\nACGTACGTACGTAAACCCGGGTTT\n>q2\nACGTACGTACGTAAACCCGGGTTA\n');
  await paste(page, 'subject', `>s1\n${'ACGTACGTAC'.repeat(10)}\n`);
  // With two query records included, the query shows no subrange at all.
  await expect(page.getByTestId('query-region')).toHaveCount(0);
  await expect(page.getByTestId('query-region-unavailable')).toHaveCount(0);
  await expect(page.getByTestId('subject-region').locator('legend')).toHaveText('Subject subrange');
  await expect(page.getByTestId('subject-region-record')).toHaveText('Record s1');

  // Drag across the middle of the bar.
  const bar = page.getByTestId('subject-region-bar');
  await bar.scrollIntoViewIfNeeded();
  const box = (await bar.boundingBox())!;
  await page.mouse.move(box.x + box.width * 0.25, box.y + box.height / 2);
  await page.mouse.down();
  await page.mouse.move(box.x + box.width * 0.75, box.y + box.height / 2, { steps: 5 });
  await page.mouse.up();
  const start = Number(await page.getByTestId('subject-region-start').inputValue());
  const stop = Number(await page.getByTestId('subject-region-stop').inputValue());
  expect(start).toBeGreaterThan(20);
  expect(start).toBeLessThan(30);
  expect(stop).toBeGreaterThan(70);
  expect(stop).toBeLessThan(80);

  // Typed positions are limited to the record.
  await page.getByTestId('subject-region-start').fill('5');
  await page.getByTestId('subject-region-stop').fill('101');
  await expect(page.getByTestId('subject-region-problem')).toHaveText(
    'The stop must be between 1 and 100, the length of the record.',
  );
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText(
    'Subject region: The stop must be between 1 and 100, the length of the record.',
  );
  if (BUILD_HAS_ENGINE) {
    // A range of one letter is NCBI's error, which the engine's validate gives.
    await page.getByTestId('subject-region-stop').fill('5');
    await expect(page.getByTestId('argv-validation')).toHaveText(
      'BLAST engine error: Invalid specification of subject location (range cannot be empty)',
    );
  }
  await page.getByTestId('subject-region-stop').fill('060');
  await expect(page.getByTestId('subject-region-argument')).toHaveText('Positions in the results are those of the whole record.');
  // The query has one record once the other is excluded.
  await showRecords(page, 'query');
  await page.getByTestId('query-source-0-record-1').uncheck();
  await page.getByTestId('query-region-start').fill('3');
  await page.getByTestId('query-region-stop').fill('20');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('results-run-options')).toContainText('-subject_loc 5-60');
  expect((await result(page, 1, 6)).command).toBe(
    'losat blastn -query query.fa -subject subject.fa -query_loc 3-20 -subject_loc 5-60 -outfmt 6',
  );
});

test("the options form: the engine's defaults are not written; tasks, templates and genetic codes", async ({ page }) => {
  await program(page, 'tblastx');
  // TBLASTX has no task (no Task Selection); its genetic codes are in the query's and the subject's blocks.
  await expect(page.getByTestId('program-selection')).toHaveCount(0);
  await expect(page.getByTestId('query-panel').getByTestId('param-query_gencode')).toBeVisible();
  await expect(page.getByTestId('subject-panel').getByTestId('param-db_gencode')).toBeVisible();
  await openParameters(page);
  await expect(page.getByTestId('param-evalue')).toHaveAttribute('placeholder', '10');
  const queryCodes = page.getByTestId('param-query_gencode').locator('option');
  await expect(queryCodes).toHaveCount(27); // Default and the engine's 26 codes
  await expect(queryCodes.nth(0)).toHaveText('Default (1. Standard)');
  await expect(page.getByTestId('param-query_gencode')).toContainText('11. Bacterial, Archaeal and Plant Plastid');
  await page.getByTestId('param-evalue').fill('10');
  // A value equal to the engine's default is not written, and not marked.
  await expect(page.getByTestId('param-evalue-field')).not.toHaveClass(/\bchanged\b/);
  await page.getByTestId('param-db_gencode').selectOption('11');
  await expect(page.getByTestId('subject-gencode-note')).toContainText('Results with this code can differ from NCBI BLAST+, which translates local subjects otherwise (approved LOSAT exception).');
  await expect(page.getByTestId('param-db_gencode-field')).toHaveClass(/\bchanged\b/);
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAC\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTAC\n');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  expect((await result(page, 1, 7)).command).toBe('losat tblastx -query query.fa -subject subject.fa -db_gencode 11 -outfmt 7');
  // The run's mark is one of its words; the approved exception is in the run's details and in
  // the mark's title, not in the mark (the maintainer, 2026-11-04).
  await chooseRun(page, 1);
  const badge = page.getByTestId('results-summary').getByTestId('verification-badge');
  await expect(badge).toHaveText(BUILD_HAS_ENGINE ? /^(Supported|Not yet supported)$/ : 'Development build');
  if (BUILD_HAS_ENGINE) {
    await expect(badge).toHaveAttribute('title', /Approved exception: LOSAT translates the local subjects with -db_gencode 11\./);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('verification-exception')).toContainText('-db_gencode 11');
  }
  await page.getByTestId('tab-search').click();

  // TBLASTN's subject codes include 32 (the approved exception); its tasks (tblastn and
  // tblastn-fast) are radio buttons.
  await program(page, 'tblastn');
  await expect(page.getByTestId('param-db_gencode').locator('option[value="32"]')).toHaveCount(1);
  await expect(page.getByTestId('param-task-field')).toContainText('Algorithm');
  await expect(page.getByTestId('param-task').getByRole('radio')).toHaveCount(2);
  await expect(page.getByTestId('param-task-tblastn')).toBeChecked();

  // BLASTN: Discontiguous Word Options appear for discontiguous megablast, Template length
  // first; they stay while they hold a value; choosing the blastn task empties and hides them.
  await program(page, 'blastn');
  const legends = page.getByTestId('parameter-form').locator('.section-legend');
  const sections = ['General Parameters', 'Scoring Parameters', 'Filters and Masking', 'Other Parameters'];
  const withTemplates = [...sections.slice(0, 3), 'Discontiguous Word Options', 'Other Parameters'];
  await expect(legends).toHaveText(sections);
  await expect(page.getByTestId('param-template_type')).toHaveCount(0);
  await task(page, 'dc-megablast');
  await expect(legends).toHaveText(withTemplates);
  const length = (await page.getByTestId('param-template_length').boundingBox())!;
  const type = (await page.getByTestId('param-template_type').boundingBox())!;
  expect(length.y).toBeLessThan(type.y);
  await page.getByTestId('param-template_type').selectOption('optimal');
  await page.getByTestId('param-template_length').selectOption('21');
  await task(page, 'megablast');
  await expect(legends).toHaveText(withTemplates);
  await task(page, 'blastn');
  await expect(page.getByTestId('param-task-blastn')).toBeChecked();
  await expect(legends).toHaveText(sections);
  await task(page, 'dc-megablast');
  await expect(page.getByTestId('param-template_type')).toHaveValue('');
  await expect(page.getByTestId('param-template_length')).toHaveValue('');
  if (BUILD_HAS_ENGINE) {
    await task(page, 'blastn');
    await page.getByTestId('param-word_size').selectOption('other');
    await page.getByTestId('param-word_size-other').fill('3');
    await expect(page.getByTestId('argv-validation')).toContainText(
      "error: invalid value '3' for '-word_size <WORD_SIZE>': expected an integer >= 4",
    );
    await page.getByTestId('param-word_size-other').fill('');
    // Templates with the megablast task (the default, not written): the section stays, and the engine refuses them.
    await task(page, 'dc-megablast');
    await page.getByTestId('param-template_type').selectOption('coding');
    await page.getByTestId('param-template_length').selectOption('18');
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await task(page, 'megablast');
    await expect(page.getByTestId('argv-validation')).toContainText(
      'BLAST query/options error: Invalid discontiguous template parameters: word size must be either 11 or 12',
    );
  }
});

test("the engine's refusals of option values are shown as the engine writes them, and nothing is started", async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine's messages need the engine build");
  await program(page, 'tblastx');
  await openParameters(page);
  // A toolkit word in a free-text value (S08+a, D14): the arguments of the Low complexity regions filter.
  await page.getByTestId('param-seg-value').fill('-version');
  await expect(page.getByTestId('argv-validation')).toHaveText(
    "error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX",
  );
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAC\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTAC\n');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText(
    "error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX",
  );
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);
  // A value that LOSAT does not support, and one that NCBI refuses.
  await program(page, 'blastp');
  await page.getByTestId('param-comp_based_stats').selectOption('1');
  await expect(page.getByTestId('argv-validation')).toHaveText(
    "-comp_based_stats 1 is not supported by LOSAT's BLASTP",
  );
  await page.getByTestId('param-comp_based_stats').selectOption('2');
  await page.getByTestId('param-evalue').fill('0');
  await expect(page.getByTestId('argv-validation')).toContainText('BLAST query/options error:');
});

// WA-FM-01: scores without a Karlin-Altschul table row are not refused by the form's check. The CLI
// meets the tables only when a query has an ungapped block, so the run gives the CLI's result: a
// warning and no hit for 3,-1, and the table's error for `-task blastn -reward 7` (penalty -3).
test('BLASTN scores without a table row: the form accepts them and the run gives the engine result', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine's checks need the engine build");
  await program(page, 'blastn');
  await openParameters(page);
  const validation = page.getByTestId('argv-validation');
  await page.getByTestId('param-matchscores').selectOption('other');
  await page.getByTestId('param-reward').fill('3');
  await page.getByTestId('param-penalty').fill('-1');
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await paste(page, 'query', '>q\nACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA\n');
  await paste(page, 'subject', '>s\nACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA\n');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await backToSearch(page);
  await task(page, 'blastn');
  await page.getByTestId('param-matchscores').selectOption('other');
  await page.getByTestId('param-reward').fill('7');
  await page.getByTestId('param-penalty').fill('-3');
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 2, 'failed');
  expect(await runMessage(page, 2)).toContain('Substitution scores 7 and -3 are not supported');
});

// The index scan reads with the program's reader (SURVEY #3): NCBI's nucleotide reader stores 19 of
// these 24 letters (not E F I L P Q), so the record table and the engine's register agree on 19.
test('a protein-like record in a nucleotide input: the record table counts the letters the reader stores; the run completes', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', '>gene_x\nMKLVVLAAGGHHKLMKLVVLAAGG\n>q\nACGTACGTACGTAAACCCGGGTTT\n');
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'ok');
  await showRecords(page, 'query');
  const row = page.getByTestId('query-source-0-records').locator('.record-row').first();
  await expect(row.locator('.record-length')).toHaveText('19');
  await expect(row.locator('.badge.warn')).toHaveText('protein?');
  await submit(page);
  await waitStatus(page, 1, 'completed');
});

test('a query without a defline runs; the engine names it as NCBI does', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', 'ACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  // The results name it as the tabular outputs do.
  await chooseRun(page, 1);
  await expect(page.getByTestId('results-query-id')).toHaveText('Query_1');
  if (BUILD_HAS_ENGINE) expect((await result(page, 1, 6)).output).toMatch(/^Query_1\ts1\t/m);
});

test('empty inputs, every record excluded, a sequence without a defline, and the subrange of BLASTX', async ({ page }) => {
  await program(page, 'blastn');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Add the query sequences: paste them or open FASTA files.');
  // A gap line fails the index: NCBI's reader reads it as data of the record before it.
  await paste(page, 'query', '>q\nACGT\n>?10\nACGT\n');
  await expect(page.getByTestId('query-source-0-error')).toHaveText(
    "This input cannot be read: line 3 is a gap line ('>?'), whose residues have no bytes in the input for LOSAT's index to locate; gap lines are not supported in the browser",
  );
  if (BUILD_HAS_ENGINE) {
    // A first line that NCBI's data loaders would fetch as a Seq-id fails; a defline can go first.
    await paste(page, 'query', 'AB123456\nACGTACGTACGTACGT');
    await expect(page.getByTestId('query-source-0-error')).toContainText('may be a sequence identifier');
    await page.getByTestId('query-source-0-add-defline').click();
    await settled(page, 'query');
    await expect(page.getByTestId('query-input')).toHaveValue('>pasted_query\nAB123456\nACGTACGTACGTACGT');
  }
  // NCBI's reader reads a sequence without a defline as a record without a title (DW-38).
  await paste(page, 'query', 'ACGTACGTACGTACGT');
  await expect(page.getByTestId('query-source-0')).toHaveAttribute('data-status', 'ready');
  await expect(page.getByTestId('query-source-0-summary')).toContainText('1 record');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Add the subject sequences: paste them or open FASTA files.');
  await paste(page, 'subject', '>s\nACGTACGTACGT\n');
  await showRecords(page, 'subject');
  await page.getByTestId('subject-source-0-exclude-shown').click();
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Every subject record is excluded.');
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);

  // BLASTX takes -query_loc (DW-39): the query's only record shows its subrange.
  await expect(page.getByTestId('query-region')).toBeVisible();
  await program(page, 'blastx');
  await expect(page.getByTestId('param-query_gencode')).toBeVisible();
  await expect(page.getByTestId('query-region')).toBeVisible();
  await expect(page.getByTestId('run-losat')).toBeEnabled();
});

// A mode that reads a role with the other FASTA reader kind reads the source again; coming back to the mode that the
// source was set up in shows its table as it was left (WA-IN-01).
test("records left out and a subrange stay through a mode that reads the role with the other kind", async ({ page }) => {
  await program(page, 'blastn');
  const records = '>q1\nACGTACGTACGTAAACCCGGGTTT\n>q2\nACGTACGTACGTAAACCCGGGTTA\n>q3\nACGTACGTACGTAAACCCGGGTTC\n';
  await paste(page, 'query', records);
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await showRecords(page, 'query');
  await page.getByTestId('query-source-0-record-1').uncheck();
  await page.getByTestId('query-source-0-record-2').uncheck();
  await expect(page.getByTestId('query-source-0-summary')).toContainText('(1 included)');
  await page.getByTestId('query-region-start').fill('2');
  await page.getByTestId('query-region-stop').fill('12');

  // BLASTP reads the query with the protein flags: another table, with nothing left out.
  await program(page, 'blastp');
  await settled(page, 'query');
  await expect(page.getByTestId('query-source-0-summary')).not.toContainText('included');
  await program(page, 'blastn');
  await settled(page, 'query');
  await expect(page.getByTestId('query-source-0-summary')).toContainText('(1 included)');
  await showRecords(page, 'query');
  await expect(page.getByTestId('query-source-0-record-0')).toBeChecked();
  await expect(page.getByTestId('query-source-0-record-1')).not.toBeChecked();
  await expect(page.getByTestId('query-source-0-record-2')).not.toBeChecked();
  await expect(page.getByTestId('query-region-start')).toHaveValue('2');
  await expect(page.getByTestId('query-region-stop')).toHaveValue('12');
  // The search takes both: the first record only, with its subrange.
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await expect(page.getByTestId('results-run-options')).toContainText('-query_loc 2-12');
  expect((await result(page, 1, 6)).command).toBe('losat blastn -query query.fa -subject subject.fa -query_loc 2-12 -outfmt 6');
});

// An input without a sequence is run as the command line runs it (rc 0 and a warning), and the results say
// that the query has none, not that the query filters hid it (WA-IN-02).
const EMPTY_QUERY_MODES = [
  ['blastn', '  \n\n', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n'],
  ['blastp', '', '>s1\nMKVLAAGIVGLLLAQPSTEWMKVLAAGIVGLLL\n'],
  ['blastx', '\n', '>s1\nMKVLAAGIVGLLLAQPSTEWMKVLAAGIVGLLL\n'],
  ['tblastn', ' ', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n'],
  ['tblastx', '\n\n\n', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n'],
] as const;
for (const [id, query, subject] of EMPTY_QUERY_MODES) {
  test(`an empty query runs in ${id}; its results say that the query has no sequences, with the engine warning`, async ({ page }) => {
    test.skip(id !== 'blastn' && !BUILD_HAS_ENGINE, 'the other modes are checked with the engine');
    await program(page, id);
    if (query === '') {
      // An empty box has no source: an empty file does.
      await openFiles(page, 'query', [{ name: 'empty.fa', text: '' }]);
    } else {
      await paste(page, 'query', query);
    }
    await paste(page, 'subject', subject);
    await expect(page.locator('[data-testid^="query-source-"][data-status="ready"]')).toHaveCount(1);
    await expect(page.getByTestId('query-source-0')).toContainText('0 records');
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await pickRun(page, 1);
    const notice = page.locator('[data-testid="results-notice"][data-kind="empty-query"]');
    await expect(notice).toContainText('The query has no sequences.');
    await expect(page.locator('[data-testid="results-notice"][data-kind="no-queries"]')).toHaveCount(0);
    if (BUILD_HAS_ENGINE) await expect(page.getByTestId('results-notice-engine')).toHaveText(/Warning: \[\w+\] Query is Empty!/);
  });
}

// The browser can remove the origin's files under an open tab (clearing the site's data). The next run makes the tab's
// directory again; if it cannot, it says what happened, not the browser's own words (WA-IN-03). WebKit's build of
// Playwright has no OPFS, and keeps results in memory.
test("a run after the browser removed the tab's storage still completes", async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', '>q1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  test.skip((await page.getByTestId('storage-status').getAttribute('data-backend')) !== 'opfs', 'this browser keeps results in memory');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  // A synchronous access handle holds its file for a moment after it is closed: the removal is tried again.
  await expect
    .poll(
      () =>
        page.evaluate(async () => {
          try {
            await (await navigator.storage.getDirectory()).removeEntry('tmp', { recursive: true });
            return 'removed';
          } catch (error) {
            return (error as Error).name;
          }
        }),
      { timeout: 10_000 },
    )
    .toBe('removed');
  await backToSearch(page);
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await chooseRun(page, 2);
  // What the removal took with it says so when it is opened.
  await pickRun(page, 1);
  await expect(page.getByTestId('results-status')).toContainText('The browser removed the stored results of this run');
  await expect(page.getByTestId('results-status')).toContainText('Run it again.');
});

test('the sequence kind warning is an estimate and does not change the program', async ({ page }) => {
  await program(page, 'blastp');
  await paste(page, 'query', '>n\nACGTACGTACGTACGTACGTACGTACGTAC\n');
  await expect(page.getByTestId('query-source-0-kind-warning')).toContainText(
    '1 record looks like nucleotide sequences, but BLASTP reads the query as protein',
  );
  await expect(page.getByTestId('program-blastp')).toBeChecked();
});

test.describe('leaving and returning', () => {
  test.skip(!BUILD_HAS_ENGINE, 'needs a search that runs long enough (engine build)');

  test('the wake lock is held during a search; the page warns before it is left; the state is checked on return', async ({
    page,
    browserName,
  }) => {
    await page.addInitScript(() => {
      const log: string[] = [];
      (window as unknown as { __wakeLog: string[] }).__wakeLog = log;
      let visibility: DocumentVisibilityState = 'visible';
      Object.defineProperty(document, 'visibilityState', { get: () => visibility, configurable: true });
      (window as unknown as { __setVisible: (visible: boolean) => void }).__setVisible = (visible) => {
        visibility = visible ? 'visible' : 'hidden';
        document.dispatchEvent(new Event('visibilitychange'));
      };
      // A handheld device (coarse pointer): "Keep searching while you're away" is shown there.
      const matchMedia = window.matchMedia.bind(window);
      window.matchMedia = (query: string) =>
        query.includes('pointer: coarse') ? ({ matches: true, media: query } as MediaQueryList) : matchMedia(query);
      Object.defineProperty(navigator, 'wakeLock', {
        configurable: true,
        value: {
          request: async () => {
            log.push('request');
            const sentinel = new EventTarget() as EventTarget & { release(): Promise<void> };
            sentinel.release = async () => {
              log.push('release');
              sentinel.dispatchEvent(new Event('release'));
            };
            return sentinel;
          },
        },
      });
    });
    await openSearch(page);
    // On unless turned off (Owner 2026-11-03).
    await expect(page.getByTestId('keep-awake')).toBeChecked();
    await startSlowSearch(page);
    await expect(page.getByTestId('keep-awake').locator('xpath=..')).toHaveAttribute('data-wake-lock', 'held');
    await page.evaluate(() => (window as unknown as { __setVisible: (v: boolean) => void }).__setVisible(false));
    await page.waitForTimeout(1200);
    await page.evaluate(() => (window as unknown as { __setVisible: (v: boolean) => void }).__setVisible(true));
    await expect(page.getByTestId('resume-notice')).toContainText(/This tab was in the background for 0:0[1-9]/);
    await expect(page.getByTestId('resume-run-1')).toContainText('Run 1 was running when the tab was hidden');
    await expect(page.getByTestId('resume-data')).toHaveText('The stored results are available.');

    if (browserName === 'chromium') {
      // Leaving the page while a search runs asks first.
      const dialog = page.waitForEvent('dialog');
      void page.close({ runBeforeUnload: true });
      const shown = await dialog;
      expect(shown.type()).toBe('beforeunload');
      await shown.dismiss();
    }
    await page.getByTestId('run-status-cancel').click();
    await waitStatus(page, 1, 'cancelled');
    await expect(page.getByTestId('keep-awake').locator('xpath=..')).not.toHaveAttribute('data-wake-lock', 'held');
    await page.getByTestId('resume-dismiss').click();
    await expect(page.getByTestId('resume-notice')).toHaveCount(0);
    const log = await page.evaluate(() => (window as unknown as { __wakeLog: string[] }).__wakeLog);
    expect(log[0]).toBe('request');
    expect(log.at(-1)).toBe('release');
  });
});

test('the header, and the "?" help that shows on hover or keyboard focus and stays on click (DW-41)', async ({ page }) => {
  await openSearch(page);
  await expect(page).toHaveTitle('LOSAT');
  await expect(page.locator('.app-header h1')).toHaveText('LOSAT');
  await expect(page.locator('.app-header .tagline')).toHaveText('A miniaturized reimplementation of NCBI BLAST algorithm');
  // A desktop (fine pointer) is not offered the wake lock's choice; its screen does not turn off by itself.
  await expect(page.getByTestId('keep-awake')).toHaveCount(0);

  const toggle = page.getByTestId('param-task-help-toggle');
  const help = page.getByTestId('param-task-help');
  await expect(help).toHaveCount(0);
  // Pointed at: shown; pointed away: gone (as NCBI's helpbox).
  await toggle.hover();
  await expect(help).toHaveText('Task to execute: megablast, blastn, dc-megablast or blastn-short');
  await expect(help).toHaveAttribute('data-pinned', 'false');
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');
  await page.mouse.move(0, 0);
  await expect(help).toHaveCount(0);
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  // Clicked: stays while the pointer is elsewhere; a click elsewhere closes it.
  await toggle.click();
  await expect(help).toHaveAttribute('data-pinned', 'true');
  await page.mouse.move(0, 0);
  await expect(help).toBeVisible();
  await page.getByTestId('program-summary').click();
  await expect(help).toHaveCount(0);
  // Clicked again: the next click closes it, although the pointer is still over the "?".
  await toggle.click();
  await expect(help).toHaveAttribute('data-pinned', 'true');
  await toggle.click();
  await expect(help).toHaveCount(0);
  // Escape closes it too.
  await toggle.click();
  await expect(help).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(help).toHaveCount(0);
  // Reached with the keyboard (Tab from the task's radio): shown until the focus moves on.
  await page.getByTestId('param-task-megablast').click();
  await page.keyboard.press('Tab');
  await expect(toggle).toBeFocused();
  await expect(help).toBeVisible();
  await page.keyboard.press('Tab');
  await expect(help).toHaveCount(0);
});

test("the search screen in NCBI's order and words (W4b)", async ({ page }) => {
  await expect(page.getByTestId('program-summary')).toHaveText('BLASTN mode searches nucleotide subjects using a nucleotide query.');
  await expect(page.getByTestId('query-panel').locator('legend').first()).toHaveText('Enter Query Sequence');
  await expect(page.getByTestId('subject-panel').locator('legend').first()).toHaveText('Enter Subject Sequence');
  await expect(page.getByTestId('program-selection').locator('legend')).toHaveText('Task Selection');
  // The blocks one under the other, in NCBI's order.
  const order = ['program-tabs', 'program-summary', 'query-panel', 'subject-panel', 'program-selection', 'run-losat', 'algorithm-parameters'];
  const tops = [];
  for (const id of order) tops.push((await page.getByTestId(id).boundingBox())!.y);
  expect(tops).toEqual([...tops].sort((a, b) => a - b));
  // The rows of each block.
  for (const [role, title] of [
    ['query', 'Query'],
    ['subject', 'Subject'],
  ] as const) {
    const panel = page.getByTestId(`${role}-panel`);
    await expect(panel).toContainText('Enter accession number(s) or FASTA sequence(s)');
    // The subrange appears only once one record is included.
    await expect(page.getByTestId(`${role}-region`)).toHaveCount(0);
    await expect(panel).not.toContainText(`${title} subrange`);
    await expect(panel).toContainText('or drop files here');
    await expect(panel).toContainText('Or, upload file');
    await expect(page.getByTestId(`${role}-input`)).toHaveAccessibleName('Enter accession number(s) or FASTA sequence(s)');
  }
  await expect(page.getByTestId('query-panel')).toContainText('Job Title');
  await expect(page.getByTestId('query-panel')).toContainText('Enter a descriptive title for your search');
  await expect(page.getByTestId('subject-panel')).not.toContainText('Job Title');
  // "Clear" empties the paste box only.
  await openFiles(page, 'query', [{ name: 'q.fa', text: '>q\nACGTACGTACGT\n' }]);
  await paste(page, 'query', '>pasted\nACGTACGTACGT\n');
  await expect(page.locator('[data-testid^="query-source-"][data-status]')).toHaveCount(2);
  await page.getByTestId('query-clear').click();
  await expect(page.getByTestId('query-input')).toHaveValue('');
  await expect(page.locator('[data-testid^="query-source-"][data-status]')).toHaveCount(1);
  await expect(page.getByTestId('query-source-0')).toContainText('q.fa');

  // Task Selection: BLASTN's "Optimize for", the engine's default task checked.
  await expect(page.getByTestId('param-task-field')).toContainText('Optimize for');
  await expect(page.getByTestId('param-task').getByRole('radio')).toHaveCount(4);
  await expect(page.getByTestId('param-task')).toHaveAccessibleName(/^Optimize for -task/);
  for (const [value, label] of [
    ['megablast', 'Highly similar sequences (megablast)'],
    ['dc-megablast', 'More dissimilar sequences (discontiguous megablast)'],
    ['blastn', 'Somewhat similar sequences (blastn)'],
    ['blastn-short', 'Short sequences (blastn-short)'],
  ] as const) {
    await expect(page.getByTestId('param-task').locator('label', { has: page.getByTestId(`param-task-${value}`) })).toHaveText(label);
  }
  await expect(page.getByTestId('param-task-megablast')).toBeChecked();
  // The search button (NCBI's BLAST button) and the line that says what it searches.
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT');
  const line = page.getByTestId('search-summary-line');
  await expect(line).toHaveText(
    'Search nucleotide subjects with BLASTN, task megablast (optimized for highly similar sequences).',
  );
  await task(page, 'dc-megablast');
  await expect(line).toContainText('task dc-megablast (optimized for more dissimilar sequences)');

  // Algorithm parameters: NCBI's sections and names.
  await openParameters(page);
  const form = page.getByTestId('parameter-form');
  for (const label of [
    'Max target sequences',
    'Expect threshold',
    'Word size',
    'Match/Mismatch Scores',
    'Gap Costs',
    'Filter',
    'Low complexity regions',
    'Mask lower case letters',
    'Template length',
    'Template type',
  ]) {
    await expect(form.locator('.param-label', { hasText: label }).first()).toBeVisible();
  }
  // Each control shows the value that the search uses (discontiguous megablast's here), as NCBI's
  // form shows its values; NCBI's lists are selects, and NCBI's sentence is beside its control.
  await expect(page.getByTestId('param-word_size')).toHaveValue('11');
  await expect(page.getByTestId('param-word_size').locator('option')).toHaveText(['11', '12', 'Other']);
  // NCBI's pairs are one select each; the gap costs are those of the scores' list (no Linear: megablast only).
  await expect(page.getByTestId('param-matchscores')).toHaveValue('2,-3');
  await expect(page.getByTestId('param-matchscores').locator('option')).toHaveText(['1,-2', '1,-3', '1,-4', '2,-3', '4,-5', '1,-1', 'Other']);
  await expect(page.getByTestId('param-gapcosts')).toHaveValue('5 2');
  await expect(page.getByTestId('param-gapcosts').locator('option').first()).toHaveText('Existence: 4 Extension: 4');
  await expect(page.getByTestId('param-reward')).toHaveCount(0);
  await expect(page.getByTestId('param-dust')).toBeChecked();
  await expect(page.getByTestId('param-dust-value')).toHaveValue('20 64 1');
  await expect(page.getByTestId('param-soft_masking')).toBeChecked();
  await expect(page.getByTestId('param-max_target_seqs')).toHaveValue('500');
  await expect(page.getByTestId('param-max_target_seqs-field')).toContainText('Select the maximum number of aligned sequences to display');
  await expect(line.locator('strong')).toHaveText(['nucleotide subjects', 'BLASTN, task dc-megablast (optimized for more dissimilar sequences)']);
  await expect(page.getByTestId('param-matchscores')).toHaveAccessibleName('Match/Mismatch Scores -reward -penalty');
  await expect(page.getByTestId('param-gapcosts')).toHaveAccessibleName('Gap Costs -gapopen -gapextend');
  await expect(page.getByTestId('param-dust')).toHaveAccessibleName('Low complexity regions -dust');
  await expect(page.getByTestId('param-dust-value')).toHaveAccessibleName('Low complexity regions -dust arguments');

  await program(page, 'blastp');
  await expect(page.getByTestId('program-summary')).toHaveText('BLASTP mode searches protein subjects using a protein query.');
  await expect(page.getByTestId('param-task-field')).toContainText('Algorithm');
  await expect(page.getByTestId('param-task').locator('label', { has: page.getByTestId('param-task-blastp') })).toHaveText(
    'blastp (protein-protein BLAST)',
  );
  await expect(page.getByTestId('param-task').locator('label', { has: page.getByTestId('param-task-blastp-fast') })).toHaveText(
    'Quick BLASTP (blastp-fast)',
  );
  await expect(line).toHaveText('Search protein subjects with BLASTP, task blastp (protein-protein BLAST).');
  await expect(form.locator('.section-legend')).toHaveText(['General Parameters', 'Scoring Parameters', 'Filters and Masking', 'Other Parameters']);
  // The select shows the engine's default as its value, as NCBI's does.
  const cbs = page.getByTestId('param-comp_based_stats');
  await expect(cbs).toHaveValue('2');
  await expect(cbs.locator('option[value="2"]')).toHaveText('2: Conditional compositional score matrix adjustment');
  await expect(cbs.locator('option[value="0"]')).toHaveText('0: No adjustment');
  await expect(cbs.locator('option[value="3"]')).toHaveText('3: Universal compositional score matrix adjustment');
  await expect(page.getByTestId('param-seg-field')).toContainText('Low complexity regions');
  await expect(page.getByTestId('param-seg')).not.toBeChecked();
  await expect(page.getByTestId('param-seg-value')).toHaveCount(0);
  await expect(page.getByTestId('param-gapcosts')).toHaveValue('11 1');
  await program(page, 'tblastx');
  await expect(page.getByTestId('param-culling_limit-field')).toContainText('Max matches in a query range');
  await expect(line).toHaveText(
    'Search translated nucleotide subjects with TBLASTX (translated nucleotide query).',
  );
});

test('Algorithm parameters: closed at start, the changed values marked and counted, the defaults restored; the Job Title and the second button', async ({
  page,
}) => {
  const toggle = page.getByTestId('algorithm-parameters-toggle');
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  await expect(toggle).toHaveText(/^\+\s*Algorithm parameters\s*$/);
  await expect(page.getByTestId('param-evalue')).toHaveCount(0);
  await expect(page.getByTestId('run-losat-bottom')).toHaveCount(0);
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');
  await expect(toggle).toHaveText(/^−\s*Algorithm parameters\s*$/);
  await expect(page.getByTestId('algorithm-parameters')).toContainText(
    'Parameter values that differ from the default are highlighted in yellow and marked with ♦',
  );
  const evalue = page.getByTestId('param-evalue-field');
  await expect(evalue).not.toHaveClass(/\bchanged\b/);
  await page.getByTestId('param-evalue').fill('1e-3');
  await page.getByTestId('param-matchscores').selectOption('2,-3');
  await expect(evalue).toHaveClass(/\bchanged\b/);
  await expect(evalue).toContainText('♦');
  await expect(page.getByTestId('param-evalue')).toHaveAccessibleName(/^Expect threshold -evalue\s*, changed from the default$/);
  await expect(page.getByTestId('param-matchscores-field')).toHaveClass(/\bchanged\b/);
  await expect(page.getByTestId('param-matchscores')).toHaveAccessibleName(/^Match\/Mismatch Scores -reward -penalty\s*, changed from the default$/);
  await expect(page.getByTestId('param-gapcosts-field')).not.toHaveClass(/\bchanged\b/);
  // The count is that of the options written: -evalue, -reward and -penalty.
  await expect(toggle).toContainText('(3 changed)');
  // The task is chosen in Task Selection: it is not counted. 2,-3 is the blastn task's
  // default: it is neither counted nor written any more.
  await task(page, 'blastn');
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(1 changed)');
  await expect(page.getByTestId('param-matchscores-field')).not.toHaveClass(/\bchanged\b/);

  await page.getByTestId('job-title').fill('  Globin check  ');
  await paste(page, 'query', '>q\nACGTACGTACGTAAACCCGGGTTT\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTAAACCCGGGTTTACGT\n');
  await submit(page);
  await expect(page.getByTestId('results-run-options')).toHaveText('-task blastn -evalue 1e-3');
  await expect(page.getByTestId('results-job-title')).toHaveText('Globin check');
  await waitStatus(page, 1, 'completed');
  await backToSearch(page);

  // Closed, the bar keeps the count; a program with changed values opens it when it is chosen.
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  await expect(toggle).toContainText('(1 changed)');
  await program(page, 'blastp');
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  await program(page, 'blastn');
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');

  // "Restore default search parameters" clears the Algorithm parameters, not the task.
  await page.getByTestId('restore-defaults').click();
  // The controls show the blastn task's defaults again.
  await expect(page.getByTestId('param-evalue')).toHaveValue('10');
  await expect(page.getByTestId('param-matchscores')).toHaveValue('2,-3');
  await expect(evalue).not.toHaveClass(/\bchanged\b/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  await expect(page.getByTestId('param-task-blastn')).toBeChecked();
  await page.getByTestId('job-title').fill('');
  // The second button, below the parameters, starts the same search.
  const bottom = page.getByTestId('run-losat-bottom');
  await expect(bottom).toHaveText('Run LOSAT');
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await bottom.click();
  await expect(page.getByTestId('results-run-options')).toHaveText('-task blastn');
  await expect(page.getByTestId('results-job-title')).toHaveCount(0);
  // A started search shows no message: the page is on the run.
  await expect(page.getByTestId('search-message')).toHaveCount(0);
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'true');
  await waitStatus(page, 2, 'completed');
});

test("Algorithm parameters show the task's values; NCBI's lists with Other; a default chosen is not written, a given hit list size is", async ({
  page,
}) => {
  await program(page, 'blastn');
  await openParameters(page);
  const words = page.getByTestId('param-word_size');
  const wordsField = page.getByTestId('param-word_size-field');
  await expect(words).toHaveValue('28');
  await task(page, 'blastn');
  await expect(words).toHaveValue('11');
  await words.selectOption('15');
  await expect(wordsField).toHaveClass(/\bchanged\b/);
  await words.selectOption('11');
  await expect(wordsField).not.toHaveClass(/\bchanged\b/);
  await words.selectOption('other');
  await page.getByTestId('param-word_size-other').fill('9');
  await expect(wordsField).toHaveClass(/\bchanged\b/);
  await expect(words).toHaveValue('other');
  // A true/false option with a default is a checkbox in NCBI's "Mask" row.
  const mask = page.getByTestId('param-soft_masking');
  await expect(mask).toBeChecked();
  await mask.uncheck();
  await expect(page.getByTestId('param-soft_masking-field')).toHaveClass(/\bchanged\b/);
  // Max target sequences shows its default only while the other hit list options are not set
  // (BLAST+ tells them given from not given).
  await page.getByTestId('param-num_alignments').fill('3');
  await expect(page.getByTestId('param-max_target_seqs')).toHaveValue('');
  await paste(page, 'query', '>q\nACGTACGTACGTAAACCCGGGTTT\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTAAACCCGGGTTTACGT\n');
  await submit(page);
  await expect(page.getByTestId('results-run-options')).toHaveText('-task blastn -word_size 9 -soft_masking false -num_alignments 3');
  await waitStatus(page, 1, 'completed');
});

test("NCBI's selects of pairs follow the scores and the matrix; Other shows the two fields; Low complexity regions is a checkbox", async ({
  page,
}) => {
  await program(page, 'blastn');
  await openParameters(page);
  const scores = page.getByTestId('param-matchscores');
  const gaps = page.getByTestId('param-gapcosts');
  // Megablast's Linear (0 0) is its default.
  await expect(gaps).toHaveValue('0 0');
  await expect(gaps.locator('option:checked')).toHaveText('Linear');
  // "Other" with nothing typed stays while other state changes (review of SURVEY #13, 3).
  await gaps.selectOption('other');
  await page.getByTestId('job-title').fill('pairs');
  await expect(page.getByTestId('param-gapopen')).toHaveValue('0');
  await gaps.selectOption('0 0');
  await expect(page.getByTestId('param-gapopen')).toHaveCount(0);
  await page.getByTestId('job-title').fill('');
  await task(page, 'blastn');
  await expect(gaps).toHaveValue('5 2');
  // 4,-5 has its own list, whose default NCBI chooses (BLAST stops with 5 2 there).
  await scores.selectOption('4,-5');
  await expect(gaps).toHaveValue('12 8');
  await expect(gaps.locator('option:checked')).toHaveText('Existence: 12 Extension: 8');
  await expect(page.getByTestId('param-gapcosts-field')).toHaveClass(/\bchanged\b/);
  // Other: the two fields, with the values that the search uses.
  await gaps.selectOption('other');
  await expect(page.getByTestId('param-gapopen')).toHaveValue('12');
  await expect(page.getByTestId('param-gapopen')).toHaveAccessibleName(/^Existence -gapopen/);
  await page.getByTestId('param-gapopen').fill('6');
  await page.getByTestId('param-gapextend').fill('5');
  await expect(gaps).toHaveValue('other');
  // The filter: off writes no; on again is the default.
  const dust = page.getByTestId('param-dust');
  await dust.uncheck();
  await expect(page.getByTestId('param-dust-value')).toHaveCount(0);
  await expect(page.getByTestId('param-dust-field')).toHaveClass(/\bchanged\b/);
  await paste(page, 'query', '>q\nACGTACGTACGTAAACCCGGGTTT\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTAAACCCGGGTTTACGT\n');
  await submit(page);
  await expect(page.getByTestId('results-run-options')).toHaveText('-task blastn -reward 4 -penalty -5 -gapopen 6 -gapextend 5 -dust no');
  await waitStatus(page, 1, 'completed');
  await backToSearch(page);
  await dust.check();
  await expect(page.getByTestId('param-dust-value')).toHaveValue('');
  await expect(page.getByTestId('param-dust-value')).toHaveAttribute('placeholder', '20 64 1');
  await expect(page.getByTestId('param-dust-field')).not.toHaveClass(/\bchanged\b/);
  // blastn-short filters with no by default: on is yes, and the text box takes the arguments.
  await task(page, 'blastn-short');
  await expect(dust).not.toBeChecked();
  await dust.check();
  await expect(page.getByTestId('param-dust-value')).toHaveValue('yes');
  await page.getByTestId('param-dust-value').fill('20 32 1');
  await expect(dust).toBeChecked();

  // The protein programs' gap costs follow the matrix: its default, as NCBI's page and the engine have it.
  await program(page, 'blastp');
  const matrix = page.getByTestId('param-matrix');
  await expect(gaps).toHaveValue('11 1');
  await gaps.selectOption('10 1');
  await expect(page.getByTestId('param-gapcosts-field')).toHaveClass(/\bchanged\b/);
  await matrix.selectOption('PAM30');
  await expect(gaps).toHaveValue('9 1');
  await expect(page.getByTestId('param-gapcosts-field')).not.toHaveClass(/\bchanged\b/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(1 changed)');
  const seg = page.getByTestId('param-seg');
  await seg.check();
  await expect(page.getByTestId('param-seg-value')).toHaveValue('yes');
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(2 changed)');
  // BLASTX's list holds the one pair that LOSAT's BLASTX reads.
  await program(page, 'blastx');
  await expect(gaps.locator('option')).toHaveText(['Existence: 11 Extension: 1', 'Other']);
});

test("Additional options: the command line's other options after the form's, marked and counted, read by the engine; Max matches in a query range in every mode (DW-42)", async ({
  page,
}) => {
  await openParameters(page);
  for (const id of ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']) {
    await program(page, id);
    await expect(page.getByTestId('param-culling_limit-field')).toContainText('Max matches in a query range');
    // The last row of the form, in Other Parameters.
    await expect(page.getByTestId('parameter-form').locator('.parameter-section').last().locator('.param-row').last()).toHaveAttribute(
      'data-testid',
      'param-extra-field',
    );
  }
  await program(page, 'blastn');
  const field = page.getByTestId('param-extra-field');
  const extra = page.getByTestId('param-extra');
  const validation = page.getByTestId('argv-validation');
  await expect(field).not.toHaveClass(/\bchanged\b/);
  await expect(extra).toHaveAccessibleName('Additional options');
  // The engine reads the words as the command line does, and says what it does not read.
  await extra.fill('-qcov_hsp_perc 50 -foo');
  await expect(validation).toHaveText("error: unknown option or argument '-foo'; use -help for CLI v2 syntax");
  await extra.fill('-num_threads 2');
  await expect(validation).toHaveText('-num_threads is managed by LOSAT');
  await extra.fill("-strand 'plus");
  await expect(validation).toHaveText("Additional options: a quote (') is not closed");
  // An option of the form is set there.
  await extra.fill('-ungapped -evalue=1');
  await expect(validation).toHaveText('-evalue is set in the form: Expect threshold');
  await extra.fill("-ungapped  -strand 'plus'");
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await expect(field).toHaveClass(/\bchanged\b/);
  await expect(field).toContainText('♦');
  await expect(extra).toHaveAccessibleName(/^Additional options\s*, changed from the default$/);
  await page.getByTestId('param-culling_limit').fill('1');
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(2 changed)');
  await paste(page, 'query', '>q\nACGTACGTACGTAAACCCGGGTTT\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTAAACCCGGGTTTACGT\n');
  await submit(page);
  await expect(page.getByTestId('results-run-options')).toHaveText('-culling_limit 1 -ungapped -strand plus');
  await waitStatus(page, 1, 'completed');
  await backToSearch(page);
  // Each mode has its own; "Restore default search parameters" empties them.
  await program(page, 'blastp');
  await expect(extra).toHaveValue('');
  await program(page, 'blastn');
  await expect(extra).toHaveValue("-ungapped  -strand 'plus'");
  await page.getByTestId('restore-defaults').click();
  await expect(extra).toHaveValue('');
  await expect(field).not.toHaveClass(/\bchanged\b/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
});

test('narrow screens put the inputs one under the other, without horizontal scrolling', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await paste(page, 'query', '>NZ_CP006932.1 first\nACGTACGTACGTACGTAC\n>NZ_CP006932.1 second\nACGT!ACGT\n');
  // A record row puts its tags on a second line: the ID and the "refused" tag stay in the list.
  // Only the FakeEngine refuses the record (see REFUSED_RECORD), so only it shows the tag.
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', BUILD_HAS_ENGINE ? 'ok' : 'refused');
  await showRecords(page, 'query');
  const list = (await page.getByTestId('query-source-0-records').locator('.record-viewport').boundingBox())!;
  if (!BUILD_HAS_ENGINE) {
    const tag = (await page.getByTestId('query-source-0-refused-tag').boundingBox())!;
    expect(tag.x + tag.width).toBeLessThanOrEqual(list.x + list.width);
  }
  for (const id of await page.getByTestId('query-source-0-records').locator('.record-id').all()) {
    expect(await id.evaluate((element) => element.scrollWidth <= element.clientWidth)).toBe(true);
  }
  const query = (await page.getByTestId('query-panel').boundingBox())!;
  const subject = (await page.getByTestId('subject-panel').boundingBox())!;
  expect(subject.y).toBeGreaterThan(query.y + query.height - 1);
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  // The query block's paste box, drop zone and Job Title start and end at the same x (W4b screen review L13).
  const edges = await Promise.all(
    ['query-input', 'query-dropzone', 'job-title'].map(async (id) => {
      const box = (await page.getByTestId(id).boundingBox())!;
      return [Math.round(box.x), Math.round(box.x + box.width)];
    }),
  );
  expect(edges.slice(1)).toEqual([edges[0], edges[0]]);

  // With Algorithm parameters open, and the longest labels (W4b): still no horizontal scrolling.
  await openParameters(page);
  await page.getByTestId('param-evalue').fill('1e-5');
  for (const id of ['blastn', 'blastp', 'tblastn', 'tblastx']) {
    await program(page, id);
    if (id === 'tblastn') await page.getByTestId('param-db_gencode').selectOption('11');
    await expect(page.getByTestId('parameter-form')).toBeVisible();
    const width = await page.evaluate(() => [document.documentElement.scrollWidth, window.innerWidth]);
    expect(width[0], `${id}: the page is wider than the window`).toBeLessThanOrEqual(width[1]!);
  }
  // The labels are above their controls at this width.
  const label = (await page.getByTestId('param-evalue-field').locator('.param-label').boundingBox())!;
  const input = (await page.getByTestId('param-evalue').boundingBox())!;
  expect(input.y).toBeGreaterThanOrEqual(label.y + label.height - 1);
});
