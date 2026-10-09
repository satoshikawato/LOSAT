// The search screen (S12, W3; plan §7): research work and boundaries through the real
// application. The tests run with the FakeEngine build and with the engine build
// (LOSAT_WEB_REACTORS); what only the engine can show (its messages, the kept subject, a
// search long enough to edit or cancel while it runs) is checked in the engine build.
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import {
  fasta,
  openFiles,
  openParameters,
  paste,
  program,
  result,
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

/** Where a part of a run's card is: its top below the card's top, and its right edge from the card's right edge. */
async function placeInCard(page: Page, number: number, part: string): Promise<{ top: number; right: number }> {
  const card = (await page.getByTestId(`run-${number}`).boundingBox())!;
  const box = (await page.getByTestId(`run-${number}-${part}`).boundingBox())!;
  return { top: Math.round(box.y - card.y), right: Math.round(card.x + card.width - box.x - box.width) };
}

async function startSlowSearch(page: Page): Promise<void> {
  await program(page, SLOW.program);
  await openFiles(page, 'query', [{ name: SLOW.query, text: fasta(SLOW.query) }]);
  await openFiles(page, 'subject', [{ name: SLOW.subject, text: fasta(SLOW.subject) }]);
  await submit(page);
  await expect(page.getByTestId('run-1-status')).toHaveText('running', { timeout: 60_000 });
}

// Engine searches are slower in Firefox and WebKit than in Chromium (W1 README).
test.setTimeout(300_000);

test.beforeEach(async ({ page }) => {
  await page.goto('/');
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

// A record that the checker refuses. The FakeEngine names the record in its message ("query record 3
// (bad) ..."), so the screen offers to exclude it. The engine's register (S10) reads with NCBI's
// reader: it refuses a gap line ('>?') with a message about the line, which names no record, so the
// screen shows the refusal without the exclude button and the record is excluded with its checkbox.
const REFUSED_RECORD = BUILD_HAS_ENGINE ? '>?10\nACGT\n' : '>bad\nACGT!ACGT\n';

test('records: duplicate IDs by number, the engine refuses a record, exclusion, and a run after exclusion', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', `>dup first\nACGTACGTACGTAAACCCGGGTTT\n>dup second\nACGTACGTACGTAAACCCGGGTTA\n${REFUSED_RECORD}`);
  await paste(page, 'subject', '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n');
  await expect(page.getByTestId('query-source-0-duplicates')).toContainText('2 records share an ID');
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-check', 'refused');
  await expect(check).toContainText(BUILD_HAS_ENGINE ? 'not supported by LOSAT Web' : 'query record 3 (bad)');
  if (!BUILD_HAS_ENGINE) await expect(check).toContainText('FAKE ENGINE check');
  // The refused input cannot be queued.
  await submit(page);
  await expect(page.getByTestId('search-message')).toContainText(
    BUILD_HAS_ENGINE ? 'Query (pasted): ' : 'Query (pasted): query record 3 (bad)',
  );
  await expect(page.getByTestId('run-1')).toHaveCount(0);

  if (BUILD_HAS_ENGINE) {
    await expect(page.getByTestId('query-source-0-exclude-refused')).toHaveCount(0);
    await showRecords(page, 'query');
    await page.getByTestId('query-source-0-record-2').uncheck();
  } else {
    await page.getByTestId('query-source-0-exclude-refused').click();
  }
  await expect(check).toHaveAttribute('data-check', 'ok');
  await expect(page.getByTestId('query-source-0-summary')).toContainText('(2 included)');
  await submit(page);
  await waitStatus(page, 1, 'completed');

  // Leave out the first "dup" (by its number) and search again.
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
  await paste(page, 'query', '>qb\nACGTACGTACGTAAACCCGGGTTA\n');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  if (BUILD_HAS_ENGINE) {
    await page.getByTestId('run-2').locator('summary', { hasText: 'Details' }).click();
    await expect(page.getByTestId('run-2-details').locator('[data-detail="subject-retained"]')).toHaveText(
      'kept from the previous search',
    );
    expect((await result(page, 2, 7)).output).toContain('# Query: qb');
  }
});

test('several runs in the queue; the next job is edited while one runs; separate searches as a group', async ({ page }) => {
  if (BUILD_HAS_ENGINE) {
    await startSlowSearch(page);
  } else {
    await program(page, 'blastp');
    await openFiles(page, 'query', [{ name: 'q1.fa', text: '>q1\nMKLVVLAAGGHHKL\n' }]);
    await openFiles(page, 'subject', [{ name: 's1.fa', text: '>s1\nMKLVVLAAGGHHKL\n' }]);
    await submit(page);
    await expect(page.getByTestId('run-1')).toBeVisible();
  }

  // Edit the next job while run 1 runs: another program, query and option, two subjects searched separately.
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
  await expect(page.getByTestId('add-to-queue')).toHaveText(/^Run LOSAT \(2 runs(, add to queue)?\)$/);
  // While run 1 runs, the button says that the next runs wait in the queue.
  if (BUILD_HAS_ENGINE) await expect(page.getByTestId('add-to-queue')).toHaveText('Run LOSAT (2 runs, add to queue)');
  await expect(page.getByTestId('add-to-queue-bottom')).toHaveText(/^Run LOSAT \(2 runs(, add to queue)?\)$/);
  await submit(page);
  await expect(page.getByTestId('run-3')).toBeVisible();
  await expect(page.getByTestId('run-2')).toContainText('query.fa vs a.fa · group run 1 of 2');
  await expect(page.getByTestId('run-3')).toContainText('query.fa vs b.fa · group run 2 of 2');
  await expect(page.getByTestId('run-3-options')).toHaveText('Options: -evalue 1e-3');
  // Run 1 kept its own program, inputs and options.
  await expect(page.getByTestId('run-1')).toContainText('BLASTP');
  await expect(page.getByTestId('run-1-options')).toHaveText('Options: defaults');

  if (BUILD_HAS_ENGINE) {
    await expect(page.getByTestId('run-1')).toContainText(`${SLOW.query} vs ${SLOW.subject}`);
    // The group waits behind run 1; cancelling the group leaves run 1 running.
    await expect(page.getByTestId('run-2-status')).toHaveText('queued');
    // One card layout (W4b screen review L7): the phase line and Cancel are at the same place in
    // the running card and in the waiting one.
    await expect(page.getByTestId('run-2-phase')).toHaveText('Waiting');
    expect(await placeInCard(page, 2, 'progress')).toEqual(await placeInCard(page, 1, 'progress'));
    expect(await placeInCard(page, 2, 'cancel')).toEqual(await placeInCard(page, 1, 'cancel'));
    await page.getByTestId('run-2-cancel-group').click();
    await expect(page.getByTestId('run-2-status')).toHaveText('cancelled');
    await expect(page.getByTestId('run-3-status')).toHaveText('cancelled');
    await expect(page.getByTestId('run-1-status')).toHaveText('running');
    await expect(page.getByTestId('run-1-phase')).toHaveText('Searching');
    await expect(page.getByTestId('run-1-elapsed')).toHaveText(/^\d+:\d\d$/);
    await submit(page); // the same group again, behind run 1
    await expect(page.getByTestId('run-5-status')).toHaveText('queued');
    await page.getByTestId('run-1-cancel').click();
    await waitStatus(page, 1, 'cancelled');
    await waitStatus(page, 4, 'completed');
    await waitStatus(page, 5, 'completed');
  } else {
    await expect(page.getByTestId('run-1')).toContainText('q1.fa vs s1.fa');
    await waitStatus(page, 1, 'completed');
    await waitStatus(page, 2, 'completed');
    await waitStatus(page, 3, 'completed');
  }
  const last = await result(page, BUILD_HAS_ENGINE ? 5 : 3, 6);
  expect(last.command).toBe('LOSAT tblastx -query query.fa -subject b.fa -evalue 1e-3 -outfmt 6');
  // A finished card has the same lines: the run and its status, the time it took, "Open results"
  // at the right, then what it searched.
  const number = BUILD_HAS_ENGINE ? 5 : 3;
  await expect(page.getByTestId(`run-${number}-phase`)).toHaveText('Took');
  const box = async (part: string) => (await page.getByTestId(`run-${number}-${part}`).boundingBox())!;
  const [head, progress, open, options] = [await box('status'), await box('progress'), await box('open'), await box('options')];
  expect(progress.y).toBeGreaterThanOrEqual(head.y + head.height - 1);
  expect(open.y).toBeGreaterThanOrEqual(progress.y + progress.height - 1);
  expect(options.y).toBeGreaterThanOrEqual(open.y + open.height - 1);
  expect((await placeInCard(page, number, 'open')).right).toBeLessThanOrEqual(12);
});

test('cancel a running search; the next one runs', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, 'needs a search that runs long enough to cancel (engine build)');
  await startSlowSearch(page);
  await page.getByTestId('run-1-cancel').click();
  await waitStatus(page, 1, 'cancelled');
  // The same subject, another query: the new runtime registers the subject again.
  await page.getByTestId('query-source-0-remove').click();
  await paste(page, 'query', '>small\nMKLVVLAAGGHHKLMKLVVLAAGGHHKL\n');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await page.getByTestId('run-2').locator('summary', { hasText: 'Details' }).click();
  await expect(page.getByTestId('run-2-details').locator('[data-detail="subject-retained"]')).toHaveText(
    'read for this search',
  );
});

test('the region of a role with one record: the drag, the fields, the limits and the argv', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', '>q1\nACGTACGTACGTAAACCCGGGTTT\n>q2\nACGTACGTACGTAAACCCGGGTTA\n');
  await paste(page, 'subject', `>s1\n${'ACGTACGTAC'.repeat(10)}\n`);
  await expect(page.getByTestId('query-region-unavailable')).toContainText('(2 are)');
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
      'The engine refuses these options: BLAST engine error: Invalid specification of subject location (range cannot be empty)',
    );
  }
  await page.getByTestId('subject-region-stop').fill('060');
  await expect(page.getByTestId('subject-region-argument')).toContainText('-subject_loc 5-60');
  // The query has one record once the other is excluded.
  await showRecords(page, 'query');
  await page.getByTestId('query-source-0-record-1').uncheck();
  await page.getByTestId('query-region-start').fill('3');
  await page.getByTestId('query-region-stop').fill('20');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  expect((await result(page, 1, 6)).command).toBe(
    'LOSAT blastn -query query.fa -subject subject.fa -query_loc 3-20 -subject_loc 5-60 -outfmt 6',
  );
});

test("the options form: the engine's defaults are not written; tasks, templates and genetic codes", async ({ page }) => {
  await program(page, 'tblastx');
  // TBLASTX has no task (no Program Selection); its genetic codes are in the query's and the subject's blocks.
  await expect(page.getByTestId('program-selection')).toHaveCount(0);
  await expect(page.getByTestId('query-panel').getByTestId('param-query_gencode')).toBeVisible();
  await expect(page.getByTestId('subject-panel').getByTestId('param-db_gencode')).toBeVisible();
  await openParameters(page);
  await expect(page.getByTestId('param-evalue')).toHaveAttribute('placeholder', 'default: 10');
  const queryCodes = page.getByTestId('param-query_gencode').locator('option');
  await expect(queryCodes).toHaveCount(27); // Default and the engine's 26 codes
  await expect(queryCodes.nth(0)).toHaveText('Default (1. Standard)');
  await expect(page.getByTestId('param-query_gencode')).toContainText('11. Bacterial, Archaeal and Plant Plastid');
  await page.getByTestId('param-evalue').fill('10');
  // A value equal to the engine's default is not written, and not marked.
  await expect(page.getByTestId('param-evalue-field')).not.toHaveClass(/\bchanged\b/);
  await page.getByTestId('param-db_gencode').selectOption('11');
  await expect(page.getByTestId('subject-gencode-note')).toContainText('Approved LOSAT exception');
  await expect(page.getByTestId('param-db_gencode-field')).toHaveClass(/\bchanged\b/);
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAC\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTAC\n');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  expect((await result(page, 1, 7)).command).toBe('LOSAT tblastx -query query.fa -subject subject.fa -db_gencode 11 -outfmt 7');

  // TBLASTN's subject codes include 32 (the approved exception); its tasks are radio buttons.
  await program(page, 'tblastn');
  await expect(page.getByTestId('param-db_gencode').locator('option[value="32"]')).toHaveCount(1);
  await expect(page.getByTestId('param-task-field')).toContainText('Algorithm');
  await expect(page.getByTestId('param-task').getByRole('radio')).toHaveCount(2);
  await expect(page.getByTestId('param-task-tblastn')).toBeChecked();

  // BLASTN: Discontiguous Word Options appear for discontiguous megablast, Template length
  // first; they stay while they hold a value; choosing the blastn task empties and hides them.
  await program(page, 'blastn');
  const legends = page.getByTestId('parameter-form').locator('legend');
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
    await page.getByTestId('param-word_size').fill('3');
    await expect(page.getByTestId('argv-validation')).toContainText(
      "The engine refuses these options: error: invalid value '3' for '-word_size <WORD_SIZE>': expected an integer >= 4",
    );
    await page.getByTestId('param-word_size').fill('');
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

test("the engine's refusals of option values are shown as the engine writes them, and nothing is queued", async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine's messages need the engine build");
  await program(page, 'tblastx');
  await openParameters(page);
  // A toolkit word in a free-text value (S08+a, D14).
  await page.getByTestId('param-seg').fill('-version');
  await expect(page.getByTestId('argv-validation')).toHaveText(
    "The engine refuses these options: error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX",
  );
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAC\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTAC\n');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText(
    "error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX",
  );
  await expect(page.getByTestId('run-1')).toHaveCount(0);
  // A value that LOSAT does not support, and one that NCBI refuses.
  await program(page, 'blastp');
  await page.getByTestId('param-comp_based_stats').selectOption('1');
  await expect(page.getByTestId('argv-validation')).toHaveText(
    "The engine refuses these options: -comp_based_stats 1 is not supported by LOSAT's BLASTP",
  );
  await page.getByTestId('param-comp_based_stats').selectOption('');
  await page.getByTestId('param-evalue').fill('0');
  await expect(page.getByTestId('argv-validation')).toContainText('BLAST query/options error:');
});

test('empty inputs, every record excluded, a sequence without a defline, and BLASTX', async ({ page }) => {
  await program(page, 'blastn');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Add the query sequences: paste them or open FASTA files.');
  await paste(page, 'query', 'ACGTACGTACGTACGT');
  await expect(page.getByTestId('query-source-0-error')).toHaveText('This input cannot be read: Expected > at record start.');
  await page.getByTestId('query-source-0-add-defline').click();
  await settled(page, 'query');
  await expect(page.getByTestId('query-input')).toHaveValue('>pasted_query\nACGTACGTACGTACGT');
  await expect(page.getByTestId('query-source-0')).toHaveAttribute('data-status', 'ready');
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Add the subject sequences: paste them or open FASTA files.');
  await paste(page, 'subject', '>s\nACGTACGTACGT\n');
  await showRecords(page, 'subject');
  await page.getByTestId('subject-source-0-exclude-shown').click();
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Every subject record is excluded.');
  await expect(page.getByTestId('run-1')).toHaveCount(0);

  await program(page, 'blastx');
  await expect(page.getByTestId('program-unavailable')).toContainText('BLASTX joins LOSAT Web after its certification');
  await expect(page.getByTestId('add-to-queue')).toBeDisabled();
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
    await page.goto('/');
    await page.getByTestId('keep-awake').check();
    await startSlowSearch(page);
    await expect(page.getByTestId('wake-lock-held')).toBeVisible();
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
    await page.getByTestId('run-1-cancel').click();
    await waitStatus(page, 1, 'cancelled');
    await expect(page.getByTestId('wake-lock-held')).toHaveCount(0);
    await page.getByTestId('resume-dismiss').click();
    await expect(page.getByTestId('resume-notice')).toHaveCount(0);
    const log = await page.evaluate(() => (window as unknown as { __wakeLog: string[] }).__wakeLog);
    expect(log[0]).toBe('request');
    expect(log.at(-1)).toBe('release');
  });
});

test("the search screen in NCBI's order and words (W4b)", async ({ page }) => {
  await expect(page.getByTestId('program-summary')).toHaveText('BLASTN searches nucleotide subjects using a nucleotide query.');
  await expect(page.getByTestId('query-panel').locator('legend').first()).toHaveText('Enter Query Sequence');
  await expect(page.getByTestId('subject-panel').locator('legend').first()).toHaveText('Enter Subject Sequence');
  await expect(page.getByTestId('program-selection').locator('legend')).toHaveText('Program Selection');
  // The blocks one under the other, in NCBI's order.
  const order = ['program-tabs', 'program-summary', 'query-panel', 'subject-panel', 'program-selection', 'add-to-queue', 'algorithm-parameters'];
  const tops = [];
  for (const id of order) tops.push((await page.getByTestId(id).boundingBox())!.y);
  expect(tops).toEqual([...tops].sort((a, b) => a - b));
  // The rows of each block.
  for (const [role, title] of [
    ['query', 'Query'],
    ['subject', 'Subject'],
  ] as const) {
    const panel = page.getByTestId(`${role}-panel`);
    await expect(panel).toContainText('Enter FASTA sequence(s)');
    await expect(panel).toContainText(`${title} subrange`);
    await expect(panel).toContainText('Or, upload file');
    await expect(page.getByTestId(`${role}-input`)).toHaveAccessibleName('Enter FASTA sequence(s)');
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

  // Program Selection: BLASTN's "Optimize for", the engine's default task checked.
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
  await expect(page.getByTestId('add-to-queue')).toHaveText('Run LOSAT');
  const line = page.getByTestId('search-summary-line');
  await expect(line).toHaveText(
    'Search nucleotide subjects with BLASTN, task megablast (optimized for highly similar sequences). Runs in this browser.',
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
    'Match',
    'Mismatch',
    'Gap Costs',
    'Existence',
    'Extension',
    'Low complexity regions filter (DUST)',
    'Mask lower case letters',
    'Template length',
    'Template type',
  ]) {
    await expect(form.locator('.param-label', { hasText: label }).first()).toBeVisible();
  }
  await expect(page.getByTestId('param-reward')).toHaveAccessibleName('Match -reward');
  await expect(page.getByTestId('param-gapopen')).toHaveAccessibleName('Existence -gapopen');

  await program(page, 'blastp');
  await expect(page.getByTestId('program-summary')).toHaveText('BLASTP searches protein subjects using a protein query.');
  await expect(page.getByTestId('param-task-field')).toContainText('Algorithm');
  await expect(page.getByTestId('param-task').locator('label', { has: page.getByTestId('param-task-blastp') })).toHaveText(
    'blastp (protein-protein BLAST)',
  );
  await expect(page.getByTestId('param-task').locator('label', { has: page.getByTestId('param-task-blastp-fast') })).toHaveText(
    'Quick BLASTP (blastp-fast)',
  );
  await expect(line).toHaveText('Search protein subjects with BLASTP, task blastp (protein-protein BLAST). Runs in this browser.');
  await expect(form.locator('legend')).toHaveText(['General Parameters', 'Scoring Parameters', 'Filters and Masking', 'Other Parameters']);
  const cbs = page.getByTestId('param-comp_based_stats');
  await expect(cbs.locator('option').first()).toHaveText('Default (2: Conditional compositional score matrix adjustment)');
  await expect(cbs.locator('option[value="0"]')).toHaveText('0: No adjustment');
  await expect(cbs.locator('option[value="3"]')).toHaveText('3: Universal compositional score matrix adjustment');
  await expect(page.getByTestId('param-seg-field')).toContainText('Low complexity regions filter (SEG)');
  await program(page, 'tblastx');
  await expect(page.getByTestId('param-culling_limit-field')).toContainText('Max matches in a query range');
  await expect(line).toHaveText(
    'Search translated nucleotide subjects with TBLASTX (translated nucleotide query). Runs in this browser.',
  );
});

test('Algorithm parameters: closed at start, the changed values marked and counted, the defaults restored; the Job Title and the second button', async ({
  page,
}) => {
  const toggle = page.getByTestId('algorithm-parameters-toggle');
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  await expect(toggle).toHaveText(/^\+\s*Algorithm parameters\s*$/);
  await expect(page.getByTestId('param-evalue')).toHaveCount(0);
  await expect(page.getByTestId('add-to-queue-bottom')).toHaveCount(0);
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');
  await expect(toggle).toHaveText(/^−\s*Algorithm parameters\s*$/);
  await expect(page.getByTestId('algorithm-parameters')).toContainText(
    'Parameter values that differ from the default are highlighted in yellow and marked with ♦',
  );
  const evalue = page.getByTestId('param-evalue-field');
  await expect(evalue).not.toHaveClass(/\bchanged\b/);
  await page.getByTestId('param-evalue').fill('1e-3');
  await page.getByTestId('param-reward').fill('2');
  await expect(evalue).toHaveClass(/\bchanged\b/);
  await expect(evalue).toContainText('♦');
  await expect(page.getByTestId('param-evalue')).toHaveAccessibleName(/^Expect threshold -evalue\s*, changed from the default$/);
  await expect(page.getByTestId('param-reward-field')).toHaveClass(/\bchanged\b/);
  await expect(page.getByTestId('param-penalty-field')).not.toHaveClass(/\bchanged\b/);
  await expect(toggle).toContainText('(2 changed)');
  // The task is chosen in Program Selection: it is not counted.
  await task(page, 'blastn');
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(2 changed)');

  await page.getByTestId('job-title').fill('  Globin check  ');
  await paste(page, 'query', '>q\nACGTACGTACGTAAACCCGGGTTT\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTAAACCCGGGTTTACGT\n');
  await submit(page);
  await expect(page.getByTestId('run-1-options')).toHaveText('Options: -task blastn -evalue 1e-3 -reward 2');
  await expect(page.getByTestId('run-1-title')).toHaveText('Globin check');
  await waitStatus(page, 1, 'completed');
  // The title comes before the input names on the card (in one reading: the card changes as the run ends).
  const parts = await page.getByTestId('run-1').evaluate((card) => [...card.children].map((child) => child.className));
  expect(parts.indexOf('run-title')).toBeGreaterThan(-1);
  expect(parts.indexOf('run-title')).toBeLessThan(parts.findIndex((name) => name.startsWith('run-inputs')));

  // Closed, the bar keeps the count; a program with changed values opens it when it is chosen.
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  await expect(toggle).toContainText('(2 changed)');
  await program(page, 'blastp');
  await expect(toggle).toHaveAttribute('aria-expanded', 'false');
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  await program(page, 'blastn');
  await expect(toggle).toHaveAttribute('aria-expanded', 'true');

  // "Restore default search parameters" clears the Algorithm parameters, not the task.
  await page.getByTestId('restore-defaults').click();
  await expect(page.getByTestId('param-evalue')).toHaveValue('');
  await expect(page.getByTestId('param-reward')).toHaveValue('');
  await expect(evalue).not.toHaveClass(/\bchanged\b/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  await expect(page.getByTestId('param-task-blastn')).toBeChecked();
  await page.getByTestId('job-title').fill('');
  // The second button, below the parameters, queues the same search.
  const bottom = page.getByTestId('add-to-queue-bottom');
  await expect(bottom).toHaveText('Run LOSAT');
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await bottom.click();
  await expect(page.getByTestId('run-2-options')).toHaveText('Options: -task blastn');
  await expect(page.getByTestId('run-2-title')).toHaveCount(0);
  await expect(page.getByTestId('algorithm-parameters').getByTestId('search-message')).toHaveText('Added to the queue.');
  await waitStatus(page, 2, 'completed');
});

test('narrow screens put the inputs and the queue one under the other, without horizontal scrolling', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await paste(page, 'query', `>NZ_CP006932.1 first\nACGTACGTACGTACGTAC\n>NZ_CP006932.1 second\n${BUILD_HAS_ENGINE ? 'ACGT\n>?10\nACGT' : 'ACGT!ACGT'}\n`);
  // A record row puts its tags on a second line: the ID and the "refused" tag stay in the list.
  // The engine's refusals name no record (see REFUSED_RECORD), so only the FakeEngine shows the tag.
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'refused');
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
  const queue = (await page.getByTestId('queue').boundingBox())!;
  expect(subject.y).toBeGreaterThan(query.y + query.height - 1);
  expect(queue.y).toBeGreaterThan(subject.y);
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
