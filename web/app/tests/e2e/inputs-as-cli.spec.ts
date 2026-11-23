// FASTA inputs that LOSAT Web reads as the command line reads them (line WQA's D13 and D16, 2026-11-20;
// docs/web/abi_v2.md §4 `register`, §9 *validation* and *scan*): a `>?` gap line (D13 (a)) and lines that NCBI's
// line reader joins (D13 (b)), which LOSAT Web refused ("not supported in the browser") while the CLI reads them;
// a first line that may be a Seq-id, which the CLI refuses too, refused with the CLI's words (D13 (c)); and that
// first line read as data when the file of `-conffile` sets `[BLAST] DATA_LOADERS=none`, as the CLI reads it (D16).
// The input of (a) is an NCBI-frozen fixture of LOSAT/tests/fixtures/fasta_input, that of (b) the record of
// web/adapter/tests/fixtures/d13/joined.nuc.fa (base.q.nuc.fa's lines, one ending in CR before a `!x` line). The engine reads them, so
// the tests run in the engine build (LOSAT_WEB_REACTORS); with LOSAT_WEB_NATIVE the saved outputs, the warnings and
// the refusal are the native CLI's of the same commit. With LOSAT_WEB_SCREENS=<directory> the input panels are
// saved at 1280x800 and 390x844 as <directory>/<browser>/<shot>-<width>.png once the expectations of the panel hold or
// fail, so that the unfixed build gives its pictures too (`shownThen`).
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { inputNames, openFiles, openParameters, openSearch, paste, program, showRecords, storedOutputs, submit, waitStatus, type Role } from './support/search';

const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;

test.setTimeout(180_000);
test.skip(!BUILD_HAS_ENGINE, "the FakeEngine has its own reader: the engine build reads the inputs as the CLI's reader does");

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

const FIXTURES = 'LOSAT/tests/fixtures/fasta_input';
const fixture = (name: string): Buffer => readFileSync(join(REPOSITORY, FIXTURES, name));

/** Saves the panel of `role` at 1280x800 and 390x844 (LOSAT_WEB_SCREENS), then restores the viewport. */
async function shoot(page: Page, role: Role, shot: string, browserName: string): Promise<void> {
  if (SCREENS === undefined) return;
  const directory = join(SCREENS, browserName);
  mkdirSync(directory, { recursive: true });
  const viewport = page.viewportSize();
  for (const [width, height] of [[1280, 800], [390, 844]] as const) {
    await page.setViewportSize({ width, height });
    const panel = page.getByTestId(`${role}-panel`);
    await panel.scrollIntoViewIfNeeded();
    await panel.screenshot({ path: join(directory, `${shot}-${width}.png`) });
  }
  if (viewport !== null) await page.setViewportSize(viewport);
}

/** Runs `expectations` on the panel of `role`, and saves the panel (`shoot`) whether they hold or not. */
async function shownThen(page: Page, role: Role, shot: string, browserName: string, expectations: () => Promise<void>): Promise<void> {
  try {
    await expectations();
  } finally {
    await shoot(page, role, shot, browserName);
  }
}

/**
 * An input of a search: its bytes, under the name that the application gives it in the argv (outfmt 0 shows it): the user's file name,
 * or for pasted text the name that the run gives (`searchAsTheCli` reads it from the run's details; DECISIONS 44).
 */
interface Input {
  readonly name: string;
  readonly bytes: Uint8Array | string;
}

/** The native CLI's run of `words` (program first) on the inputs under their names, with `files` beside them. */
function nativeRun(
  words: readonly string[],
  query: Input,
  subject: Input,
  format: 0 | 6 | 7,
  files: Readonly<Record<string, string>> = {},
): { status: number | null; stdout: Buffer; stderr: Buffer } {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, query.name), query.bytes);
  writeFileSync(join(directory, subject.name), subject.bytes);
  for (const [name, text] of Object.entries(files)) writeFileSync(join(directory, name), text);
  const args = [words[0]!, '-query', query.name, '-subject', subject.name, ...words.slice(1), '-outfmt', String(format), '-num_threads', '1'];
  const result = spawnSync(NATIVE!, args, { cwd: directory, maxBuffer: 1 << 28 });
  return { status: result.status, stdout: result.stdout, stderr: result.stderr };
}

/** Runs the search and compares its saved outputs and warnings with the native CLI's run of the same inputs. */
async function searchAsTheCli(
  page: Page,
  words: readonly string[],
  query: Input,
  subject: Input,
  files: Readonly<Record<string, string>> = {},
): Promise<void> {
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const { outputs, warnings } = await storedOutputs(page, 1);
  // The names that the run's command gives the inputs (storedOutputs leaves the page at the run's details).
  const names = inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
  query = { ...query, name: names.query };
  subject = { ...subject, name: names.subject };
  if (NATIVE === undefined) {
    test.info().annotations.push({ type: 'skip', description: `The bytes of the outputs: ${NO_NATIVE_REASON}` });
    return;
  }
  for (const format of [0, 6, 7] as const) {
    const cli = nativeRun(words, query, subject, format, files);
    expect(cli.status, `-outfmt ${format}: ${cli.stderr}`).toBe(0);
    expect(outputs[format].equals(cli.stdout), `-outfmt ${format}`).toBe(true);
    expect(warnings, `-outfmt ${format}: the warnings`).toBe(cli.stderr.toString('utf8'));
  }
}

// D13 (a): `>?100` in the second record is a gap of 100 N, which the CLI's reader puts in the record.
test('a >? gap line is read as the CLI reads it, with the gap\'s letters in its record', async ({ page, browserName }) => {
  const query = { name: 'gap_100.q.nuc.fa', bytes: fixture('gap_100.q.nuc.fa') };
  const subject = { name: 'base.s.nuc.fa', bytes: fixture('base.s.nuc.fa') };
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: query.name, text: query.bytes }]);
  await openFiles(page, 'subject', [{ name: subject.name, text: subject.bytes }]);
  await shownThen(page, 'query', 'gap-line', browserName, async () => {
    await expect(page.getByTestId('query-source-0-summary')).toHaveText('3 records · 640 nt');
    await showRecords(page, 'query');
    // The row of record #2 (its checkbox's parent): q2 with its 100 N.
    const row = page.getByTestId('query-source-0-record-1').locator('..');
    await expect(row.locator('.record-id')).toHaveText('q2');
    await expect(row.locator('.record-length')).toHaveText('280');
  });
  await searchAsTheCli(page, ['blastn'], query, subject);
});

// D13 (b): the file has LF line ends; NCBI's line reader joins the line that ends in CR alone with the `!x` after it,
// and the CLI reads 210 of the record's 270 letters. LOSAT Web refused it ("NCBI BLAST+ reads two of its lines as one").
test('lines that NCBI\'s line reader joins are read as the CLI reads them', async ({ page, browserName }) => {
  const query = { name: 'joined.nuc.fa', bytes: readFileSync(join(REPOSITORY, 'web/adapter/tests/fixtures/d13/joined.nuc.fa')) };
  const subject = { name: 'base.s.nuc.fa', bytes: fixture('base.s.nuc.fa') };
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: query.name, text: query.bytes }]);
  await openFiles(page, 'subject', [{ name: subject.name, text: subject.bytes }]);
  await shownThen(page, 'query', 'joined-lines', browserName, async () => {
    await expect(page.getByTestId('query-source-0-summary')).toHaveText('1 record · 210 nt');
    await showRecords(page, 'query');
  });
  await searchAsTheCli(page, ['blastn'], query, subject);
});

const PROTEIN = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ';
/** Pasted inputs: the run names them, and the test runs the CLI under those names. */
const SEQ_ID_QUERY = { name: 'pasted-query', bytes: `2HHB\n${PROTEIN}\n` };
const SEQ_ID_SUBJECT = { name: 'pasted-subject', bytes: `>s\n${PROTEIN}\n` };
/** The CLI's refusal of the query's first line `2HHB` with the data loaders on (the engine's reader, `seq_id.rs`). */
const SEQ_ID_REFUSAL =
  "the first line of the query (\"2HHB\") is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT's BLASTP (start the query with a '>' defline, or set DATA_LOADERS=none in the [BLAST] section of .ncbirc)";
/** The card's sentence after the CLI's, on where the browser reads DATA_LOADERS (line WQA's D13 round 2, L3, 2026-11-20). */
const LOADERS_HINT =
  'In the browser, LOSAT reads this setting from the file of -conffile, not from .ncbirc: add -conffile NAME to Additional options, and choose a file whose [BLAST] section has DATA_LOADERS=none.';

// D13 (c): the CLI refuses the input too (LOSAT does not fetch); the card gives the CLI's sentence.
test('a first line that may be a Seq-id is refused with the CLI\'s words', async ({ page, browserName }, testInfo) => {
  await program(page, 'blastp');
  await paste(page, 'query', SEQ_ID_QUERY.bytes);
  await paste(page, 'subject', SEQ_ID_SUBJECT.bytes);
  await shownThen(page, 'query', 'seq-id-words', browserName, async () => {
    await expect(page.getByTestId('query-source-0-error')).toHaveText(`This input cannot be read: ${SEQ_ID_REFUSAL}`);
    // The CLI's sentence names .ncbirc, which the browser has not: the card says where the setting goes here.
    await expect(page.getByTestId('query-source-0-loaders-hint')).toHaveText(LOADERS_HINT);
  });
  // The press that puts '>' before the line stays (the sentence's own remedy).
  await expect(page.getByTestId('query-source-0-put-defline-mark')).toBeVisible();
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastp'], SEQ_ID_QUERY, SEQ_ID_SUBJECT, 0);
  expect(cli.status).toBe(1);
  expect(cli.stderr.toString('utf8')).toBe(`Error: ${SEQ_ID_REFUSAL}\n`);
});

// D16: with `-conffile nol.ini` and a file of `[BLAST] DATA_LOADERS=none`, the CLI reads `2HHB` as data (H, H and B
// are residues; NCBI BLAST+ 2.17.0 gives the same outfmt 6, survey-d12 probe2 `sid.pdb`); the card reads the input
// again with the run's setting and the search runs.
test('with -conffile setting DATA_LOADERS=none, a first line that may be a Seq-id is read as data, as the CLI reads it', async ({
  page,
  browserName,
}) => {
  const settings = '[BLAST]\nDATA_LOADERS=none\n';
  await program(page, 'blastp');
  await paste(page, 'query', SEQ_ID_QUERY.bytes);
  await paste(page, 'subject', SEQ_ID_SUBJECT.bytes);
  await expect(page.getByTestId('query-source-0-error')).toContainText('may be a sequence identifier');
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-conffile nol.ini');
  await page.getByTestId('param-conffile-file').setInputFiles({ name: 'nol.ini', mimeType: 'text/plain', buffer: Buffer.from(settings) });
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await shownThen(page, 'query', 'data-loaders-none', browserName, async () => {
    await expect(page.getByTestId('query-source-0-summary')).toHaveText('1 record · 69 aa', { timeout: 30_000 });
    await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'ok', { timeout: 30_000 });
  });
  await searchAsTheCli(page, ['blastp', '-conffile', 'nol.ini'], SEQ_ID_QUERY, SEQ_ID_SUBJECT, { 'nol.ini': settings });
  // Without the file's setting the run reads with the data loaders on: the engine's check refuses the input again, as
  // the CLI refuses it.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('param-extra').fill('');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-check', 'refused', { timeout: 30_000 });
  await expect(check).toContainText(SEQ_ID_REFUSAL);
  // The check's refusal gets the press that puts '>' before the line, as the reader's does (line WQA's D13 round 2, L2).
  await expect(page.getByTestId('query-source-0-put-defline-mark')).toBeVisible();
});

// line WQA's D13 round 2, L2, 2026-11-20: an input read with the data loaders off of `-conffile`, which the engine's check
// refuses once the option is gone, gets the card's remedies of a Seq-id first line, as an input that the reader refuses
// gets them: the sentence on Fetch for an accession, and the press that puts '>' before the line, after which the
// engine reads the line as a defline and accepts the input. The CLI reads `A0A023GPI8` (an accession of UniProt) with
// the data loaders off as data (the letters A, A, G, P and I; LOSAT's CLI, probe of 2026-11-20) and refuses it with them on.
test("an input that the engine's check refuses for a first line that may be a Seq-id gets the card's remedies", async ({
  page,
  browserName,
}) => {
  const settings = '[BLAST]\nDATA_LOADERS=none\n';
  await program(page, 'blastp');
  await paste(page, 'query', `A0A023GPI8\n${PROTEIN}\n`);
  await paste(page, 'subject', SEQ_ID_SUBJECT.bytes);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-conffile nol.ini');
  await page.getByTestId('param-conffile-file').setInputFiles({ name: 'nol.ini', mimeType: 'text/plain', buffer: Buffer.from(settings) });
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  const check = page.getByTestId('query-source-0-check');
  await expect(check).toHaveAttribute('data-check', 'ok', { timeout: 30_000 });
  await page.getByTestId('param-extra').fill('');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await shownThen(page, 'query', 'seq-id-check-refused', browserName, async () => {
    await expect(check).toHaveAttribute('data-check', 'refused', { timeout: 30_000 });
    await expect(check).toContainText('the first line of the query ("A0A023GPI8") is not a defline and may be a sequence identifier');
    await expect(page.getByTestId('query-source-0-fetch-hint')).toHaveText('To fetch accessions, put only the accessions in the box, one per line.');
    await expect(page.getByTestId('query-source-0-put-defline-mark')).toHaveText("Put '>' before the first line");
    await expect(page.getByTestId('query-source-0-loaders-hint')).toHaveText(LOADERS_HINT);
  });
  await page.getByTestId('query-source-0-put-defline-mark').click();
  await expect(page.getByTestId('query-input')).toHaveValue(`>A0A023GPI8\n${PROTEIN}\n`);
  await expect(check).toHaveAttribute('data-check', 'ok', { timeout: 30_000 });
  await expect(page.getByTestId('query-source-0-summary')).toHaveText(`1 record · ${PROTEIN.length} aa`);
  await expect(page.getByTestId('query-source-0-loaders-hint')).toHaveCount(0);
});
