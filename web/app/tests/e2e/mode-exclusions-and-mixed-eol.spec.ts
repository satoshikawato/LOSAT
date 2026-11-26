// Line WQA's D16 (Owner Q2=B: exclusions survive a search mode change; only records newly parsed by
// the new reader enter as included) and D17 (combined inputs and record selections with mixed line
// endings no longer fail with InputMismatchError and match the native CLI's bytes on the exported
// input file).
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import {
  chooseRun,
  inputNames,
  openFiles,
  openSearch,
  paste,
  program,
  settled,
  showRecords,
  storedOutputs,
  submit,
  waitStatus,
  type Role,
} from './support/search';

test.setTimeout(180_000);

const FIXTURES = 'LOSAT/tests/fixtures/fasta_input';
const fixture = (name: string): Buffer => readFileSync(join(REPOSITORY, FIXTURES, name));

interface Held {
  readonly name: string;
  readonly bytes: Buffer;
}

const recordBox = (page: Page, role: Role, index: number) => page.getByTestId(`${role}-source-0-record-${index}`);

/** Downloads the made input file offered for `role` in Download All. */
async function downloadOfferedInput(page: Page, role: Role): Promise<Held> {
  await page.getByTestId('results-download-all').click();
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.getByTestId('results-download-menu').getByTestId(`results-download-${role}`).click(),
  ]);
  return { name: download.suggestedFilename(), bytes: readFileSync((await download.path())!) };
}

/** Runs the native CLI in a fresh folder with `files` and returns stdout and stderr. */
function nativeRun(
  words: readonly string[],
  files: readonly Held[],
  format: 0 | 6 | 7,
): { status: number | null; stdout: Buffer; stderr: Buffer } {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  for (const file of files) writeFileSync(join(directory, file.name), file.bytes);
  const result = spawnSync(NATIVE!, [...words, '-outfmt', String(format), '-num_threads', '1'], {
    cwd: directory,
    maxBuffer: 1 << 28,
    argv0: 'losat',
  });
  return { status: result.status, stdout: result.stdout, stderr: result.stderr };
}

/** Compares stored outfmt 0, 6, 7 and warnings of `run` against the native CLI using the offered input files. */
async function expectNativeParity(page: Page, run: number, programId: string, query: Held, subject: Held): Promise<void> {
  const { outputs, warnings } = await storedOutputs(page, run);
  const names = inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
  expect(names.query).toBe(query.name);
  expect(names.subject).toBe(subject.name);
  if (NATIVE === undefined) {
    test.info().annotations.push({ type: 'skip', description: `The CLI's bytes: ${NO_NATIVE_REASON}` });
    return;
  }
  const files: Held[] = [query, subject];
  for (const format of [0, 6, 7] as const) {
    const cli = nativeRun([programId, '-query', query.name, '-subject', subject.name], files, format);
    expect(cli.status, `-outfmt ${format}: ${cli.stderr.toString('utf8')}`).toBe(0);
    expect(outputs[format].equals(cli.stdout), `-outfmt ${format} bytes differ from native CLI`).toBe(true);
    expect(warnings, `-outfmt ${format} warnings differ from native CLI`).toBe(cli.stderr.toString('utf8'));
  }
}

/** Checks the summary text and each record checkbox state for source 0 of `role`. */
async function expectSelection(
  page: Page,
  role: Role,
  checked: readonly boolean[],
  summaryPatterns: readonly string[] = [],
): Promise<void> {
  const summary = page.getByTestId(`${role}-source-0-summary`);
  if (summaryPatterns.length === 0) await expect(summary).not.toContainText('included');
  for (const pattern of summaryPatterns) await expect(summary).toContainText(pattern);
  await showRecords(page, role);
  for (const [index, want] of checked.entries()) {
    const checkbox = recordBox(page, role, index);
    if (want) await expect(checkbox).toBeChecked();
    else await expect(checkbox).not.toBeChecked();
  }
}

/** Checks or unchecks a record in source 0 of `role`. */
async function toggleRecord(page: Page, role: Role, index: number, included: boolean): Promise<void> {
  await showRecords(page, role);
  const target = recordBox(page, role, index);
  if (included) await target.check();
  else await target.uncheck();
}

/** Removes the first `count` sources of `role` in reverse order. */
async function clearSources(page: Page, role: Role, count: number): Promise<void> {
  for (let i = count - 1; i >= 0; i--) {
    await page.getByTestId(`${role}-source-${i}-remove`).click();
  }
}

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

test('D16: exclusions survive mode switches, newly parsed records enter as included, and clear/replace leaks nothing', async ({ page }) => {
  await program(page, 'blastn');
  // Duplicate IDs + distinct records on query and subject.
  const queryText =
    '>dup first\nACGTACGTACGTAAACCCGGGTTT\n>dup second\nACGTACGTACGTAAACCCGGGTTA\n>q3\nACGTACGTACGTAAACCCGGGTTC\n';
  const subjectText =
    '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n>s2\nTTTTTTTTCCCCCCCCGGGGGGGGAAAAAAAA\n';
  await paste(page, 'query', queryText);
  await paste(page, 'subject', subjectText);

  await toggleRecord(page, 'query', 0, false);
  await toggleRecord(page, 'subject', 1, false);

  const queryFirstExcluded = [false, true, true] as const;
  // Verify initial state, same-kind switch (blastn -> tblastx), and cross-kind switches.
  for (const mode of [undefined, 'tblastx', 'blastp', 'blastx', 'tblastn'] as const) {
    if (mode !== undefined) {
      await program(page, mode);
      await settled(page, 'query');
      await settled(page, 'subject');
    }
    await expectSelection(page, 'query', queryFirstExcluded, ['(2 included)']);
    await expectSelection(page, 'subject', [true, false], ['(1 included)']);
  }

  // Edit exclusions in TBLASTN (protein query): re-include record 0 and exclude record 2, then
  // switch back to BLASTN and verify the edited selection persists.
  await toggleRecord(page, 'query', 0, true);
  await toggleRecord(page, 'query', 2, false);
  for (const mode of [undefined, 'blastn'] as const) {
    if (mode !== undefined) {
      await program(page, mode);
      await settled(page, 'query');
    }
    await expectSelection(page, 'query', [true, true, false], ['(2 included)']);
  }

  // Newly parsed records under another reader: under BLASTN (kind 1), r1 has only protein-only
  // letters before an invalid line, so kind 1 stops at r1 (only r0 parsed). Exclude r0, switch to
  // BLASTP (kind 2): r0 stays excluded while newly parsed r1 and r2 enter as included.
  await paste(page, 'query', '>r0\nACGTACGTACGT\n>r1\nEFFILP\n@@@@@@\n>r2\nACGTACGTACGT\n');
  await expectSelection(page, 'query', [true], ['1 record']);
  await recordBox(page, 'query', 0).uncheck();
  await expectSelection(page, 'query', [false], ['(0 included)']);

  await program(page, 'blastp');
  await settled(page, 'query');
  await expectSelection(page, 'query', queryFirstExcluded, ['3 records', '(2 included)']);

  // Clear and re-paste, and remove/re-add a file: zero stale exclusions leak.
  await page.getByTestId('query-input').fill('');
  await expect(page.getByTestId('query-source-0')).toHaveCount(0);
  await paste(page, 'query', `${queryText}`);
  await expectSelection(page, 'query', [true, true, true]);

  const subjectFile = { name: 's.fa', text: subjectText };
  await clearSources(page, 'subject', 1);
  await expect(page.getByTestId('subject-source-0')).toHaveCount(0);
  await openFiles(page, 'subject', [subjectFile]);
  await toggleRecord(page, 'subject', 0, false);
  await expectSelection(page, 'subject', [false, true], ['(1 included)']);
  await clearSources(page, 'subject', 1);
  await openFiles(page, 'subject', [{ ...subjectFile }]);
  await expectSelection(page, 'subject', [true, true]);
});

test('D17: combined inputs with mixed CR/CRLF/LF line endings complete and match native CLI bytes', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the FakeEngine has its own reader: the engine build exercises NCBI's stateful line reader");

  const subjectBytes = fixture('base.s.nuc.fa');
  const subjectHeld: Held = { name: 'base.s.nuc.fa', bytes: subjectBytes };
  const crNoFinal = fixture('eol_cr_nofinal.q.nuc.fa');
  const crFinal = fixture('eol_cr.q.nuc.fa');
  const lf = fixture('base.q.nuc.fa');

  await program(page, 'blastn');
  await openFiles(page, 'subject', [{ name: 'base.s.nuc.fa', text: subjectBytes }]);

  interface MixedCase {
    readonly run: number;
    readonly removeBefore: number;
    readonly setup: () => Promise<void>;
    readonly checkOffered?: (offered: Held) => Promise<void> | void;
  }

  const crFinalFile = { name: 'eol_cr.q.nuc.fa', text: crFinal };

  const expectQueryDetails = async (recordLabel: string, readDiff?: string, fromLine?: string) => {
    const resultsDiff = page.getByTestId('results-query-read-diff');
    if (readDiff === undefined) await expect(resultsDiff).toHaveCount(0);
    else await expect(resultsDiff).toHaveText(readDiff);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-input-query')).toContainText(recordLabel);
    if (fromLine !== undefined) await expect(page.getByTestId('run-input-query-from')).toHaveText(fromLine);
    const detailsDiff = page.getByTestId('run-input-query-read-diff');
    if (readDiff === undefined) await expect(detailsDiff).toHaveCount(0);
    else await expect(detailsDiff).toHaveText(readDiff);
  };

  const cases: readonly MixedCase[] = [
    // Case 1: pure-CR without trailing CR + pure-CR combined input.
    {
      run: 1,
      removeBefore: 0,
      setup: async () => {
        await openFiles(page, 'query', [{ name: 'eol_cr_nofinal.q.nuc.fa', text: crNoFinal }, crFinalFile]);
      },
      checkOffered: async (offered) => {
        // All 6 records stay intact (180 nt each) because pure-CR stays in CR mode: no diff notice.
        await expectQueryDetails('6 records');
        expect(offered.bytes.equals(Buffer.concat([crNoFinal, Buffer.from('\r'), crFinal]))).toBe(true);
      },
    },
    // Case 2: combined files whose record count differs when read as one file (table=5, run=6).
    {
      run: 2,
      removeBefore: 2,
      setup: async () => {
        await openFiles(page, 'query', [
          crFinalFile,
          { name: 'eol_mixed_lf_crlf_cr.q.nuc.fa', text: fixture('eol_mixed_lf_crlf_cr.q.nuc.fa') },
        ]);
      },
      checkOffered: async () => {
        await expectQueryDetails(
          '6 records',
          'Read as one file, as the command line reads it, these give 6 records (5 when each file is read alone).',
          'From eol_cr.q.nuc.fa (3 records) and eol_mixed_lf_crlf_cr.q.nuc.fa (2 records)',
        );
      },
    },
    // Case 3: pasted LF snippet + uploaded CRLF file combined input.
    {
      run: 3,
      removeBefore: 2,
      setup: async () => {
        await paste(page, 'query', lf.toString('utf8').trimEnd());
        await openFiles(page, 'query', [{ name: 'eol_crlf.q.nuc.fa', text: fixture('eol_crlf.q.nuc.fa') }]);
      },
    },
    // Case 4: record selection whose line-reader state changes a record's length when record #2 is excluded
    // (eol_lfcr_defline.q.nuc.fa without #2: table q3=180, engine reads q3=120).
    {
      run: 4,
      removeBefore: 2,
      setup: async () => {
        await openFiles(page, 'query', [{ name: 'eol_lfcr_defline.q.nuc.fa', text: fixture('eol_lfcr_defline.q.nuc.fa') }]);
        await showRecords(page, 'query');
        await page.getByTestId('query-source-0-record-1').uncheck();
        await expectSelection(page, 'query', [true, false, true], ['(2 included)']);
      },
      checkOffered: async () => {
        await expectQueryDetails(
          '2 records',
          'Read as one file, as the command line reads it, 1 record differs from the record table (record 2 (q3): length 120 instead of 180).',
        );
      },
    },
  ];

  for (const item of cases) {
    if (item.removeBefore > 0) {
      await page.getByTestId('tab-search').click();
      await clearSources(page, 'query', item.removeBefore);
    }
    await item.setup();
    await submit(page);
    await waitStatus(page, item.run, 'completed');
    await chooseRun(page, item.run);
    const offered = await downloadOfferedInput(page, 'query');
    await item.checkOffered?.(offered);
    await expectNativeParity(page, item.run, 'blastn', offered, subjectHeld);
  }
});
