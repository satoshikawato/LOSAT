// The search screen (S12, W3; plan §7): research work and boundaries through the real
// application. The tests run with the FakeEngine build and with the engine build
// (LOSAT_WEB_REACTORS); what only the engine can show (its messages, the kept subject, a
// search long enough to edit or cancel while it runs) is checked in the engine build.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';

type Role = 'query' | 'subject';

const fasta = (path: string) => readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta', path));
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
  await expect(page.getByTestId('run-1-status')).toHaveText('running', { timeout: 60_000 });
}

async function settled(page: Page, role: Role, index = 0): Promise<void> {
  const source = page.getByTestId(`${role}-source-${index}`);
  await expect(source).toHaveAttribute('data-status', /ready|failed/, { timeout: 30_000 });
  if ((await source.getAttribute('data-status')) === 'ready') {
    await expect(page.getByTestId(`${role}-source-${index}-check`)).not.toHaveAttribute('data-check', 'pending', {
      timeout: 30_000,
    });
  }
}

async function paste(page: Page, role: Role, text: string): Promise<void> {
  await page.getByTestId(`${role}-input`).fill(text);
  await settled(page, role);
}

async function openFiles(page: Page, role: Role, files: ReadonlyArray<{ name: string; text: string | Buffer }>) {
  const before = await page.locator(`[data-testid^="${role}-source-"][data-status]`).count();
  await page.getByTestId(`${role}-files`).setInputFiles(
    files.map((file) => ({ name: file.name, mimeType: 'text/plain', buffer: Buffer.from(file.text) })),
  );
  for (let i = 0; i < files.length; i++) await settled(page, role, before + i);
}

/** Opens the record list of a source (it starts open for a few records). */
async function showRecords(page: Page, role: Role, index = 0): Promise<void> {
  await page
    .getByTestId(`${role}-source-${index}`)
    .locator('details.records')
    .evaluate((details) => ((details as HTMLDetailsElement).open = true));
}

async function program(page: Page, id: string): Promise<void> {
  await page.getByTestId(`program-${id}`).check();
}

async function submit(page: Page): Promise<void> {
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('add-to-queue').click();
}

/** Waits until run `run` ends, and fails with its error if it ends otherwise than `status`. */
async function waitStatus(page: Page, run: number, status: 'completed' | 'cancelled', timeout = 120_000): Promise<void> {
  const element = page.getByTestId(`run-${run}-status`);
  await expect(element).toHaveText(/^(completed|cancelled|failed)$/, { timeout });
  const actual = await element.textContent();
  if (actual !== status) {
    const error = (await page.getByTestId(`run-${run}-error`).textContent({ timeout: 1000 }).catch(() => null)) ?? '';
    throw new Error(`run ${run} ended ${actual}, not ${status}: ${error}`);
  }
}

/** The CLI command of a completed run, and its stored output of one format. */
async function result(page: Page, run: number, format: 0 | 6 | 7): Promise<{ command: string; output: string }> {
  await page.getByTestId('tab-results').click();
  const select = page.getByTestId('result-run');
  const label = await select.locator('option', { hasText: new RegExp(`Run ${run} ·`) }).textContent();
  await select.selectOption({ label: label!.trim() });
  await page.getByTestId(`format-${format}`).click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${run}:${format}`);
  const command = (await page.getByTestId('result-command').textContent()) ?? '';
  const output = (await page.getByTestId('result-output').textContent()) ?? '';
  await page.getByTestId('tab-search').click();
  return { command, output };
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
  await page.getByTestId('param-evalue').fill('1e-3');
  await page.getByTestId('subject-source-0-remove').click();
  await openFiles(page, 'subject', [
    { name: 'a.fa', text: '>a\nACGTACGTACGTACGTACGTAAA\n' },
    { name: 'b.fa', text: '>b\nTTTACGTACGTACGTACGTACGT\n' },
  ]);
  await page.getByTestId('subject-mode-separate').check();
  await expect(page.getByTestId('add-to-queue')).toContainText('(2 runs)');
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
  await expect(page.getByTestId('subject-region')).toContainText('Region of s1');

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
  await expect(page.getByTestId('param-evalue')).toHaveAttribute('placeholder', 'default: 10');
  const queryCodes = page.getByTestId('param-query_gencode').locator('option');
  await expect(queryCodes).toHaveCount(27); // Default and the engine's 26 codes
  await expect(queryCodes.nth(0)).toHaveText('Default (1. Standard)');
  await expect(page.getByTestId('param-query_gencode')).toContainText('11. Bacterial, Archaeal and Plant Plastid');
  await page.getByTestId('param-evalue').fill('10');
  await page.getByTestId('param-db_gencode').selectOption('11');
  await expect(page.getByTestId('subject-gencode-note')).toContainText('Approved LOSAT exception');
  await paste(page, 'query', '>q\nACGTACGTACGTACGTAC\n');
  await paste(page, 'subject', '>s\nACGTACGTACGTACGTAC\n');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  expect((await result(page, 1, 7)).command).toBe('LOSAT tblastx -query query.fa -subject subject.fa -db_gencode 11 -outfmt 7');

  // TBLASTN's subject codes include 32 (the approved exception).
  await program(page, 'tblastn');
  await expect(page.getByTestId('param-db_gencode').locator('option[value="32"]')).toHaveCount(1);

  // BLASTN: choosing the blastn task empties the discontiguous template.
  await program(page, 'blastn');
  await page.getByTestId('param-task').selectOption('dc-megablast');
  await page.getByTestId('param-template_type').selectOption('optimal');
  await page.getByTestId('param-template_length').selectOption('21');
  await page.getByTestId('param-task').selectOption('blastn');
  await expect(page.getByTestId('param-template_type')).toHaveValue('');
  await expect(page.getByTestId('param-template_length')).toHaveValue('');
  if (BUILD_HAS_ENGINE) {
    await page.getByTestId('param-word_size').fill('3');
    await expect(page.getByTestId('argv-validation')).toContainText(
      "The engine refuses these options: error: invalid value '3' for '-word_size <WORD_SIZE>': expected an integer >= 4",
    );
    await page.getByTestId('param-word_size').fill('');
    await page.getByTestId('param-template_type').selectOption('coding');
    await page.getByTestId('param-template_length').selectOption('18');
    await expect(page.getByTestId('argv-validation')).toContainText(
      'BLAST query/options error: Invalid lookup table type for discontiguous Mega BLAST',
    );
  }
});

test("the engine's refusals of option values are shown as the engine writes them, and nothing is queued", async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine's messages need the engine build");
  await program(page, 'tblastx');
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
});
