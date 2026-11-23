// Line WQA's D22 (2026-11-23): an input that is not the user's file as given (records left out, several files together, pasted
// text) is a file of its own that the run's command names: the command shows the name of a file whose bytes are exactly what the run
// searched, the page offers that file (Download All, Outputs), and says which records a file left out. Run on the shown command in a
// folder with the files that the page offers (and the user's own where the command names them), the native CLI writes the bytes
// that the application stored, for outfmt 0, 6 and 7 (the engine build, LOSAT_WEB_REACTORS and LOSAT_WEB_NATIVE). With
// LOSAT_WEB_SCREENS=<directory> the query card, Run details and Outputs are saved at 1280x800 and 390x844 as
// <directory>/<browser>/<shot>-<width>.png, whether the expectations hold or not (the unfixed build gives its pictures too).
import { spawnSync } from 'node:child_process';
import { gzipSync } from 'node:zlib';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, openFiles, openSearch, paste, program, showRecords, storedOutputs, submit, waitStatus, type Role } from './support/search';

const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;

test.setTimeout(240_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine searches nothing: the engine build runs these tests');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** A sequence of `length` letters from a fixed seed (the records share none, so a record hits only its own copy). */
function sequence(seed: number, length: number): string {
  let state = seed;
  let text = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    text += 'ACGT'[(state >>> 16) & 3];
  }
  return text;
}
const R1 = sequence(11, 300);
const R2 = sequence(22, 300);
const R3 = sequence(33, 300);
const record = (id: string, text: string): string => `>${id}\n${text}\n`;
/** Three query records; the subject holds a copy of r2 and of r3 (with more) and a record that hits nothing. */
const QUERY = record('r1', R1) + record('r2', R2) + record('r3', R3);
const SUBJECT = record('s1', R2 + sequence(44, 200)) + record('s2', R3 + sequence(55, 150)) + record('s3', sequence(66, 400));
const WITHOUT_R2 = record('r1', R1) + record('r3', R3);
const WITHOUT_S1 = record('s2', R3 + sequence(55, 150)) + record('s3', sequence(66, 400));

/** A file of the native CLI's folder. */
interface Held {
  readonly name: string;
  readonly bytes: Buffer | string;
}

/** The native CLI's standard output for the command that the application showed, run in a folder with `files`. */
function nativeShown(command: string, files: readonly Held[]): Buffer {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'folder-'));
  for (const file of files) writeFileSync(join(directory, file.name), file.bytes);
  const words = command.split(' ');
  expect(words[0], command).toBe('losat');
  const run = spawnSync(NATIVE!, [...words.slice(1), '-num_threads', '1'], { cwd: directory, maxBuffer: 1 << 28, argv0: 'losat' });
  expect(run.status, `${command}: ${run.stderr.toString('utf8')}`).toBe(0);
  return run.stdout;
}

/** Saves what a panel shows at 1280x800 and 390x844 (LOSAT_WEB_SCREENS), then restores the viewport. */
async function shoot(page: Page, testid: string, shot: string, browserName: string): Promise<void> {
  if (SCREENS === undefined) return;
  const directory = join(SCREENS, browserName);
  mkdirSync(directory, { recursive: true });
  const viewport = page.viewportSize();
  for (const [width, height] of [[1280, 800], [390, 844]] as const) {
    await page.setViewportSize({ width, height });
    const panel = page.getByTestId(testid);
    await panel.scrollIntoViewIfNeeded();
    await panel.screenshot({ path: join(directory, `${shot}-${width}.png`) });
  }
  if (viewport !== null) await page.setViewportSize(viewport);
}

/** Runs `expectations`, and saves the panel (`shoot`) whether they hold or not. */
async function shownThen(page: Page, testid: string, shot: string, browserName: string, expectations: () => Promise<void>): Promise<void> {
  try {
    await expectations();
  } finally {
    await shoot(page, testid, shot, browserName);
  }
}

/** The command that Outputs shows for run `run`'s outfmt 6. */
async function shownCommand(page: Page, run: number): Promise<string> {
  await chooseRun(page, run);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-6').click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${run}:6`);
  return ((await page.getByTestId('result-command').textContent()) ?? '').trim();
}

/** The file of `role` that Download All offers for the run shown: its name and its bytes. */
async function offered(page: Page, role: Role): Promise<Held> {
  await page.getByTestId('results-download-all').click();
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('results-download-menu').getByTestId(`results-download-${role}`).click()]);
  return { name: download.suggestedFilename(), bytes: readFileSync((await download.path())!) };
}

/** Leaves record `index` of the role's first file out of the search. */
async function leaveOut(page: Page, role: Role, index: number, source = 0): Promise<void> {
  await showRecords(page, role, source);
  await page.getByTestId(`${role}-source-${source}-record-${index}`).uncheck();
}

/** Every stored output of run `run` equals the native CLI's on the shown command, in a folder with `files`. */
async function expectNative(page: Page, run: number, shown: string, files: readonly Held[]): Promise<void> {
  if (NATIVE === undefined) {
    test.info().annotations.push({ type: 'skip', description: `The CLI's bytes: ${NO_NATIVE_REASON}` });
    return;
  }
  const stored = (await storedOutputs(page, run)).outputs;
  for (const format of [0, 6, 7] as const) {
    expect(stored[format].length, `outfmt ${format} is not empty`).toBeGreaterThan(0);
    expect(nativeShown(`${shown} -outfmt ${format}`, files).equals(stored[format]), `${shown} -outfmt ${format}: the native CLI's standard output`).toBe(true);
  }
}

test('a query with a record left out: the command names the file of the records searched, the page offers it and names the record (D22)', async ({ page, browserName }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'q.fa', text: QUERY }]);
  await openFiles(page, 'subject', [{ name: 's.fa', text: SUBJECT }]);
  await leaveOut(page, 'query', 1);
  await shownThen(page, 'query-panel', 'query-card', browserName, () => expect(page.getByTestId('query-source-0-summary')).toContainText('(2 included)'));
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const shown = 'losat blastn -query included_q.fa -subject s.fa';

  await chooseRun(page, 1);
  await expect(page.getByTestId('run-1-option')).toHaveText('Run 1 · BLASTN · included_q.fa vs s.fa');
  await expect(page.getByTestId('results-query-from')).toHaveText('q.fa without r2 (2 of 3 records)');
  await expect(page.getByTestId('results-subject-from')).toHaveCount(0);
  await page.getByTestId('results-view-details').click();
  await shownThen(page, 'run-details', 'run-details', browserName, async () => {
    await expect(page.getByTestId('run-argv')).toHaveText('blastn -query included_q.fa -subject s.fa');
    await expect(page.getByTestId('run-input-query')).toContainText('included_q.fa: 2 records, ');
    await expect(page.getByTestId('run-input-query-from')).toHaveText('From q.fa without r2 (2 of 3 records)');
    for (const format of [0, 6, 7]) await expect(page.getByTestId(`run-command-${format}`)).toHaveText(`${shown} -outfmt ${format}`);
    // The subject is the user's file as given: its row says nothing more.
    await expect(page.getByTestId('run-input-subject')).toContainText('s.fa: 3 records, ');
    await expect(page.getByTestId('run-input-subject-from')).toHaveCount(0);
  });
  await page.getByTestId('results-view-outputs').click();
  await shownThen(page, 'results-outputs', 'outputs', browserName, async () => {
    expect(await shownCommand(page, 1)).toBe(`${shown} -outfmt 6`);
    await expect(page.getByTestId('outputs-query-file')).toContainText('-query included_q.fa is q.fa without r2 (2 of 3 records).');
    await expect(page.getByTestId('outputs-subject-file')).toHaveCount(0);
  });

  // The file that the page offers holds the records searched, under the name of the command.
  const file = await offered(page, 'query');
  expect(file.name).toBe('included_q.fa');
  expect(file.bytes.toString('utf8')).toBe(WITHOUT_R2);
  await page.getByTestId('results-download-all').click();
  await expect(page.getByTestId('results-download-subject')).toHaveCount(0);
  await page.keyboard.press('Escape');
  // Outputs offers it too.
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('outputs-query-download').click()]);
  expect(download.suggestedFilename()).toBe('included_q.fa');
  expect(readFileSync((await download.path())!).toString('utf8')).toBe(WITHOUT_R2);

  await expectNative(page, 1, shown, [file, { name: 's.fa', bytes: SUBJECT }]);
  // The user's file as given searches r2 too: the shown command could not be that.
  if (NATIVE !== undefined) {
    const whole = nativeShown('losat blastn -query q.fa -subject s.fa -outfmt 6', [{ name: 'q.fa', bytes: QUERY }, { name: 's.fa', bytes: SUBJECT }]);
    expect(whole.toString('utf8')).toContain('r2\ts1\t');
    expect(nativeShown(`${shown} -outfmt 6`, [file, { name: 's.fa', bytes: SUBJECT }]).toString('utf8')).not.toContain('r2\t');
  }

  // r2 included again: the user's file as given, and the page says and offers nothing more.
  await page.getByTestId('tab-search').click();
  await showRecords(page, 'query');
  await page.getByTestId('query-source-0-record-1').check();
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await chooseRun(page, 2);
  await expect(page.getByTestId('run-2-option')).toHaveText('Run 2 · BLASTN · q.fa vs s.fa');
  await expect(page.getByTestId('results-query-from')).toHaveCount(0);
  await page.getByTestId('results-download-all').click();
  await expect(page.getByTestId('results-download-menu')).toBeVisible();
  await expect(page.getByTestId('results-download-query')).toHaveCount(0);
  await expect(page.getByTestId('results-download-subject')).toHaveCount(0);
  await page.keyboard.press('Escape');
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-input-query')).toContainText('q.fa: 3 records, ');
  await expect(page.getByTestId('run-input-query-from')).toHaveCount(0);
  await page.getByTestId('results-view-outputs').click();
  expect(await shownCommand(page, 2)).toBe('losat blastn -query q.fa -subject s.fa -outfmt 6');
  await expect(page.getByTestId('outputs-query-file')).toHaveCount(0);
  await expectNative(page, 2, 'losat blastn -query q.fa -subject s.fa', [{ name: 'q.fa', bytes: QUERY }, { name: 's.fa', bytes: SUBJECT }]);
});

test('a subject with a record left out is named and offered as the query is; several files and pasted text are offered too (D22)', async ({ page }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'q.fa', text: QUERY }]);
  await openFiles(page, 'subject', [{ name: 's.fa', text: SUBJECT }]);
  await leaveOut(page, 'subject', 0);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const shown = 'losat blastn -query q.fa -subject included_s.fa';
  await chooseRun(page, 1);
  await expect(page.getByTestId('results-subject-from')).toHaveText('s.fa without s1 (2 of 3 records)');
  await expect(page.getByTestId('results-query-from')).toHaveCount(0);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-input-subject-from')).toHaveText('From s.fa without s1 (2 of 3 records)');
  await expect(page.getByTestId('run-input-query-from')).toHaveCount(0);
  await page.getByTestId('results-view-outputs').click();
  expect(await shownCommand(page, 1)).toBe(`${shown} -outfmt 6`);
  await expect(page.getByTestId('outputs-subject-file')).toContainText('-subject included_s.fa is s.fa without s1 (2 of 3 records).');
  const file = await offered(page, 'subject');
  expect(file.name).toBe('included_s.fa');
  expect(file.bytes.toString('utf8')).toBe(WITHOUT_S1);
  // outfmt 0 and 7 write the subject's name too (`Input:`): they equal the CLI's on this file.
  await expectNative(page, 1, shown, [{ name: 'q.fa', bytes: QUERY }, file]);

  // Pasted text, and two files together with a record of one left out: files of LOSAT's, offered under the names of the command.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('subject-source-0-remove').click();
  await page.getByTestId('query-source-0-remove').click();
  await paste(page, 'query', QUERY);
  const a = record('a1', R2 + sequence(44, 200)) + record('a2', sequence(77, 300));
  const b = record('b1', R3 + sequence(55, 150)) + record('b2', sequence(88, 300));
  await openFiles(page, 'subject', [{ name: 'a.fa', text: a }, { name: 'b.fa', text: b }]);
  await leaveOut(page, 'subject', 1, 1);
  await submit(page);
  await waitStatus(page, 2, 'completed');
  const together = 'losat blastn -query query.fa -subject combined_subject.fa';
  await chooseRun(page, 2);
  await expect(page.getByTestId('results-subject-from')).toHaveText('a.fa (2 records) and b.fa without b2 (1 of 2 records)');
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-input-subject-from')).toHaveText('From a.fa (2 records) and b.fa without b2 (1 of 2 records)');
  // Pasted text with nothing left out is what was pasted: no line of its own, only the file in Download All.
  await expect(page.getByTestId('run-input-query-from')).toHaveCount(0);
  await page.getByTestId('results-view-outputs').click();
  expect(await shownCommand(page, 2)).toBe(`${together} -outfmt 6`);
  await expect(page.getByTestId('outputs-query-file')).toHaveCount(0);
  await expect(page.getByTestId('outputs-subject-file')).toContainText('-subject combined_subject.fa is a.fa (2 records) and b.fa without b2 (1 of 2 records).');
  const pasted = await offered(page, 'query');
  const combined = await offered(page, 'subject');
  expect(pasted.name).toBe('query.fa');
  expect(pasted.bytes.toString('utf8')).toBe(QUERY);
  expect(combined.name).toBe('combined_subject.fa');
  expect(combined.bytes.toString('utf8')).toBe(a + record('b1', R3 + sequence(55, 150)));
  await expectNative(page, 2, together, [pasted, combined]);
});

// Review M1: one name for two inputs could not be run on the offered files; and a compressed file is named as the user gave it.
test('the same file in both roles with other records left out gets two names; a compressed file with a record left out (D22)', async ({ page }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'q.fa', text: QUERY }]);
  await openFiles(page, 'subject', [{ name: 'q.fa', text: QUERY }]);
  await leaveOut(page, 'query', 1);
  await leaveOut(page, 'subject', 0);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const shown = 'losat blastn -query included_query_q.fa -subject included_subject_q.fa';
  await chooseRun(page, 1);
  await expect(page.getByTestId('run-1-option')).toHaveText('Run 1 · BLASTN · included_query_q.fa vs included_subject_q.fa');
  await page.getByTestId('results-view-outputs').click();
  expect(await shownCommand(page, 1)).toBe(`${shown} -outfmt 6`);
  const query = await offered(page, 'query');
  const subject = await offered(page, 'subject');
  expect([query.name, subject.name]).toEqual(['included_query_q.fa', 'included_subject_q.fa']);
  expect(query.bytes.toString('utf8')).toBe(WITHOUT_R2);
  expect(subject.bytes.toString('utf8')).toBe(record('r2', R2) + record('r3', R3));
  await expectNative(page, 1, shown, [query, subject]);

  // A compressed query with a record left out: the argv names the plain file of the records, the page the file given.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('query-source-0-remove').click();
  await page.getByTestId('subject-source-0-remove').click();
  await openFiles(page, 'query', [{ name: 'q.fa.gz', text: gzipSync(QUERY) }]);
  await openFiles(page, 'subject', [{ name: 's.fa', text: SUBJECT }]);
  await leaveOut(page, 'query', 1);
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await chooseRun(page, 2);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-argv')).toHaveText('blastn -query included_q.fa -subject s.fa');
  await expect(page.getByTestId('run-input-query-from')).toHaveText('From q.fa.gz without r2 (2 of 3 records)');
  const plain = await offered(page, 'query');
  expect(plain.name).toBe('included_q.fa');
  expect(plain.bytes.toString('utf8')).toBe(WITHOUT_R2);
  await expectNative(page, 2, 'losat blastn -query included_q.fa -subject s.fa', [plain, { name: 's.fa', bytes: SUBJECT }]);
});

test('a file name with a space: the made name is the one the browser saves, and the command runs on it (D22, round 5)', async ({ page }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'two genes.fna', text: QUERY }]);
  await openFiles(page, 'subject', [{ name: 's.fa', text: SUBJECT }]);
  await leaveOut(page, 'query', 1);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const shown = 'losat blastn -query included_two_genes.fna -subject s.fa';
  expect(await shownCommand(page, 1)).toBe(`${shown} -outfmt 6`);
  // Download All and the Outputs link save the file under the name of the command, whatever the browser does to a name.
  const file = await offered(page, 'query');
  expect(file.name).toBe('included_two_genes.fna');
  expect(file.bytes.toString('utf8')).toBe(WITHOUT_R2);
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('outputs-query-download').click()]);
  expect(download.suggestedFilename()).toBe('included_two_genes.fna');
  await expectNative(page, 1, shown, [file, { name: 's.fa', bytes: SUBJECT }]);
});
