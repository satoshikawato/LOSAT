// NCBI's application start through the real application (line RG item E1w-3b; docs/web/abi_v2.md
// §4 `losat_web2_start`, §5 stream 8): `-logfile NAME`, `-conffile NAME` and `-dryrun` typed in
// Additional options work as on the command line. The file of -logfile comes with the run, in the
// Outputs of a completed run, in "Download All" and under the status of a failed run (TBLASTX writes
// NCBI's error into it); the form asks for the file that -conffile names and the run reads it under
// that name; a dry run says that it searched nothing. With LOSAT_WEB_NATIVE, the downloaded file of
// -logfile and the outputs of the -conffile run are the bytes of the native CLI of the same commit.
// The FakeEngine does not imitate the start (its validate refuses the three words), so the tests run
// in the engine build (LOSAT_WEB_REACTORS). With LOSAT_WEB_SCREENS=<directory>, the states are saved
// as PNGs.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, mkdtempSync, readFileSync, statSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Download, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, openFiles, openParameters, openSearch, paste, pickRun, program, storedOutputs, submit, waitStatus } from './support/search';

const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;

test.setTimeout(180_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine does not read -logfile, -conffile or -dryrun: the engine build runs these tests');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

async function shoot(target: Locator | Page, name: string, browserName: string): Promise<void> {
  if (SCREENS === undefined) return;
  mkdirSync(join(SCREENS, browserName), { recursive: true });
  await target.screenshot({ path: join(SCREENS, browserName, `${name}.png`) });
}

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
const SUBJECT = `>s1 first subject\n${dna(9, 120)}${dna(7, 300).slice(30, 270)}${dna(10, 120)}\n`;
/** A TBLASTX subject record without residues: `register` and `validate` take it, and the run ends with NCBI's error (results.spec.ts). */
const EMPTY_SUBJECT = '>s1 alpha\n';
const AVERAGE_TOO_SHORT = 'The average subject length is too short';

/** The bytes of the file that `click` downloads. */
async function downloaded(page: Page, click: () => Promise<void>): Promise<{ readonly name: string; readonly bytes: Buffer }> {
  const [download]: [Download, void] = await Promise.all([page.waitForEvent('download'), click()]);
  return { name: download.suggestedFilename(), bytes: readFileSync((await download.path())!) };
}

/** A new directory of the test's output for a native run, with the inputs as query.fa and subject.fa. */
function nativeDirectory(query: string, subject: string | Uint8Array): string {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  return directory;
}

/**
 * The file of `-logfile run.log` and the standard error that the native CLI writes for the search `words` (program
 * first) on the inputs, with the run's one thread: undefined without LOSAT_WEB_NATIVE.
 */
function nativeLogRun(words: readonly string[], query: string, subject: string | Uint8Array): { file: Buffer; stderr: Buffer } | undefined {
  if (NATIVE === undefined) return undefined;
  const directory = nativeDirectory(query, subject);
  const run = spawnSync(NATIVE, [words[0]!, '-query', 'query.fa', '-subject', 'subject.fa', ...words.slice(1), '-logfile', 'run.log', '-num_threads', '1'], { cwd: directory });
  return { file: readFileSync(join(directory, 'run.log')), stderr: run.stderr };
}

/** The file of `-logfile run.log` of `nativeLogRun`. */
function nativeLogFile(words: readonly string[], query: string, subject: string | Uint8Array): Buffer | undefined {
  return nativeLogRun(words, query, subject)?.file;
}

// The file of -logfile of a completed run: TBLASTX writes the warning of -max_target_seqs 1 into it,
// BLASTN to its warnings (the file is made, empty), as NCBI BLAST+ 2.17.0 does. The Outputs show it
// after the formats, and it is saved unchanged under the run's name, from there and from "Download All".
test('the file of -logfile is a file of the run: shown after the formats and saved as the CLI writes it', async ({ page, browserName }, testInfo) => {
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the files: ${NO_NATIVE_REASON}` });
  const extra = page.getByTestId('param-extra');
  const cases = [
    ['blastn', ''],
    ['tblastx', 'Warning: [tblastx] Examining 5 or more matches is recommended\n'],
  ] as const;
  for (const [index, [id, text]] of cases.entries()) {
    const run = index + 1;
    if (run > 1) await page.getByTestId('tab-search').click();
    await program(page, id);
    await paste(page, 'query', QUERY);
    await paste(page, 'subject', SUBJECT);
    await openParameters(page);
    // "Max target sequences" is a field of the form.
    await page.getByTestId('param-max_target_seqs').selectOption('other');
    await page.getByTestId('param-max_target_seqs-other').fill('1');
    await extra.fill('-logfile run.log');
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    await waitStatus(page, run, 'completed');
    await chooseRun(page, run);
    await page.getByTestId('results-view-outputs').click();
    const tab = page.getByTestId('format-run-files');
    await expect(tab).toHaveText('-logfile run.log');
    await tab.click();
    const file = page.getByTestId('run-file-0');
    await expect(file.getByTestId('run-file-name')).toHaveText('run.log');
    if (text === '') await expect(file.getByTestId('run-file-empty')).toHaveText('The file is empty.');
    else await expect(file.getByTestId('run-file-text')).toHaveText(text);
    if (id === 'tblastx') await shoot(page.getByTestId('results-outputs'), '03-outputs-logfile', browserName);
    const saved = await downloaded(page, () => page.getByTestId('run-file-download-0').click());
    expect(saved.name).toMatch(new RegExp(`^losat-run${run}-${id}-.*Z\\.logfile\\.run\\.log$`));
    expect(saved.bytes.toString('utf8')).toBe(text);
    // "Download All" lists it after the outputs, and saves the same bytes under the same name.
    await page.getByTestId('results-download-all').click();
    const item = page.getByTestId('results-download-menu').getByTestId('results-download-run-file-0');
    await expect(item).toHaveText('run.log (-logfile)');
    const again = await downloaded(page, () => item.click());
    expect(again).toEqual(saved);
    const native = nativeLogFile([id, '-max_target_seqs', '1'], QUERY, SUBJECT);
    if (native !== undefined) expect(saved.bytes.equals(native), `${id} -logfile run.log`).toBe(true);
  }
});

/**
 * What the native CLI writes for the search `words` (program first) with `-logfile NAME` on the inputs, with the run's
 * one thread, in a new working directory where the name's directories exist: its standard error, and the file that the
 * name makes (undefined where it makes none), which is the name without the '/'s at its end, as NCBI reads it (`logs/`
 * makes `logs`). Undefined without LOSAT_WEB_NATIVE.
 */
function nativeNamedLog(
  words: readonly string[],
  name: string,
  query: string,
  subject: string,
): { readonly stderr: Buffer; readonly file: Buffer | undefined } | undefined {
  if (NATIVE === undefined) return undefined;
  const directory = nativeDirectory(query, subject);
  const made = name.replace(/\/+$/, '');
  const slash = made.lastIndexOf('/');
  if (slash > 0) mkdirSync(join(directory, made.slice(0, slash)), { recursive: true });
  const run = spawnSync(NATIVE, [words[0]!, '-query', 'query.fa', '-subject', 'subject.fa', ...words.slice(1), '-logfile', name, '-num_threads', '1'], { cwd: directory });
  expect(run.status, `native ${words.join(' ')} -logfile ${name}: ${run.stderr.toString('utf8')}`).toBe(0);
  const path = join(directory, made);
  return { stderr: run.stderr, file: existsSync(path) && statSync(path).isFile() ? readFileSync(path) : undefined };
}

// A -logfile name is read as the CLI on Linux reads it where the name's directories exist (ABI v2 §4 `losat_web2_start`;
// line WQA's sweep, D11, 2026-11-20): a name with a directory makes a file of the run under its whole name, as the CLI
// makes it there (TBLASTX writes the warning of -max_target_seqs 1 into it, BLASTN to its warnings, the file empty), and
// it is saved with each '/' as '_'; /dev/stderr is the standard error, so TBLASTX's warning is among its Warnings and the
// run has no file. Before, the press was refused: "the -logfile name logs/run.log has a directory, which LOSAT Web does
// not have: give a file name without one". `logs/` is the file `logs`, as NCBI reads the name and the CLI makes it, saved
// as `…logfile.logs` (D11 round 2, 2026-11-20; before, it was a directory and the warning was among the Warnings).
test('a -logfile name with a directory makes a file of that name, and /dev/stderr is the Warnings', async ({ page, browserName }, testInfo) => {
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the files: ${NO_NATIVE_REASON}` });
  const cases = [
    { id: 'tblastx', name: 'logs/run.log', saved: 'logs_run\\.log' },
    { id: 'blastn', name: 'a/b/run.log', saved: 'a_b_run\\.log' },
    { id: 'tblastx', name: 'logs/', saved: 'logs' },
    { id: 'tblastx', name: '/dev/stderr', saved: undefined },
  ] as const;
  for (const [index, { id, name, saved }] of cases.entries()) {
    // Each case on a page of its own: its run is run 1.
    if (index > 0) await openSearch(page);
    await program(page, id);
    await paste(page, 'query', QUERY);
    await paste(page, 'subject', SUBJECT);
    await openParameters(page);
    await page.getByTestId('param-max_target_seqs').selectOption('other');
    await page.getByTestId('param-max_target_seqs-other').fill('1');
    await page.getByTestId('param-extra').fill(`-logfile ${name}`);
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    // The press starts run 1 (before D11 it was refused with the sentence about the -logfile name).
    const message = page.getByTestId('search-message');
    await expect(page.getByTestId('run-1-option').or(message).first()).toBeAttached();
    const refusal = (await message.count()) > 0 ? await message.textContent() : '';
    expect(await page.getByTestId('run-1-option').count(), `-logfile ${name}: the press is refused: ${refusal}`).toBe(1);
    await waitStatus(page, 1, 'completed');
    await chooseRun(page, 1);
    const warning = `Warning: [${id}] Examining 5 or more matches is recommended\n`;
    const native = nativeNamedLog([id, '-max_target_seqs', '1'], name, QUERY, SUBJECT);
    await page.getByTestId('results-view-outputs').click();
    if (saved === undefined) {
      // The standard error: no file of the run, and the warning is among the run's Warnings, as the CLI writes it.
      await expect(page.getByTestId('format-run-files')).toHaveCount(0);
      await page.getByTestId('results-view-details').click();
      await expect(page.getByTestId('run-diagnostics')).toHaveText(warning.trimEnd());
      await expect(page.getByTestId('run-diagnostics-place')).toHaveCount(0);
      if (native !== undefined) {
        expect(native.stderr.toString('utf8'), `${id} -logfile ${name}: the CLI's standard error`).toBe(warning);
        expect(native.file, `${id} -logfile ${name}: no file`).toBeUndefined();
      }
      continue;
    }
    const text = id === 'tblastx' ? warning : '';
    const tab = page.getByTestId('format-run-files');
    await expect(tab).toHaveText(`-logfile ${name}`);
    await tab.click();
    const file = page.getByTestId('run-file-0');
    await expect(file.getByTestId('run-file-name')).toHaveText(name);
    if (text === '') await expect(file.getByTestId('run-file-empty')).toHaveText('The file is empty.');
    else await expect(file.getByTestId('run-file-text')).toHaveText(text);
    await shoot(page.getByTestId('results-outputs'), `08-outputs-logfile-${id}-${name.replace(/[^a-z]+/g, '-')}`, browserName);
    const download = await downloaded(page, () => page.getByTestId('run-file-download-0').click());
    expect(download.name).toMatch(new RegExp(`^losat-run1-${id}-.*Z\\.logfile\\.${saved}$`));
    expect(download.bytes.toString('utf8')).toBe(text);
    await page.getByTestId('results-download-all').click();
    await expect(page.getByTestId('results-download-menu').getByTestId('results-download-run-file-0')).toHaveText(`${name} (-logfile)`);
    await page.keyboard.press('Escape');
    if (native !== undefined) {
      expect(native.file?.equals(download.bytes), `${id} -logfile ${name}: the CLI's file`).toBe(true);
      expect(native.stderr.toString('utf8'), `${id} -logfile ${name}: the CLI's standard error`).toBe(id === 'tblastx' ? '' : warning);
    }
  }
});

// TBLASTX writes the error that ends its run into the file of -logfile (ABI v2 §5: stream 8 comes also
// when the run fails). The failed run's status keeps its sentence, says that the file is kept, and
// shows the file.
test('a failed TBLASTX run keeps the file of its -logfile, with NCBI\'s error, under its status', async ({ page, browserName }) => {
  await program(page, 'tblastx');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', EMPTY_SUBJECT);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-logfile run.log');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  const status = page.getByTestId('results-status');
  await expect(status).toHaveAttribute('data-run-status', 'failed');
  await expect(status).toContainText('Run 1 failed: ');
  await expect(status).toContainText(AVERAGE_TOO_SHORT);
  await expect(status).toContainText('It keeps the outputs written before the error. It keeps the file of its -logfile, run.log.');
  const file = page.getByTestId('run-files').getByTestId('run-file-0');
  await expect(file.getByTestId('run-file-name')).toHaveText('run.log');
  await expect(file.getByTestId('run-file-text')).toContainText(AVERAGE_TOO_SHORT);
  await shoot(page.getByTestId('results-panel'), '04-failed-tblastx-logfile', browserName);
  const shown = (await file.getByTestId('run-file-text').textContent()) ?? '';
  const saved = await downloaded(page, () => file.getByTestId('run-file-download-0').click());
  expect(saved.name).toMatch(/^losat-run1-tblastx-.*Z\.logfile\.run\.log$/);
  expect(saved.bytes.toString('utf8')).toBe(shown);
  const native = nativeLogFile(['tblastx'], QUERY, EMPTY_SUBJECT);
  if (native !== undefined) expect(saved.bytes.equals(native), 'tblastx -logfile run.log, failed').toBe(true);
});

// The file of -logfile holds what the CLI's file holds for the run's own command line (line WQA's
// sweep, D3 review r2, 2026-11-19). TBLASTX with -num_alignments 0 makes two searches (ABI §7): the
// tabular formats' fails alone with NCBI's `No hits are being saved`, and the run completes as the
// CLI's run of outfmt 0, whose file is empty. The other search's warning and error, which the file
// had, are the run's warnings (ABI §5, stream 3).
test('the file of -logfile of TBLASTX holds what the CLI writes for the run\'s own command line', async ({ page, browserName }) => {
  await program(page, 'tblastx');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT);
  await openParameters(page);
  await page.getByTestId('param-num_alignments').fill('0');
  await page.getByTestId('param-extra').fill('-logfile run.log');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-run-files').click();
  const file = page.getByTestId('run-file-0');
  await expect(file.getByTestId('run-file-name')).toHaveText('run.log');
  await expect(file.getByTestId('run-file-empty')).toHaveText('The file is empty.');
  await shoot(page.getByTestId('results-outputs'), '05-outputs-logfile-two-searches', browserName);
  const saved = await downloaded(page, () => file.getByTestId('run-file-download-0').click());
  expect(saved.bytes.length).toBe(0);
  const native = nativeLogFile(['tblastx', '-num_alignments', '0'], QUERY, SUBJECT);
  if (native !== undefined) expect(saved.bytes.equals(native), 'tblastx -num_alignments 0 -logfile run.log').toBe(true);
  // The tabular formats' warning and error are the run's warnings.
  await page.getByTestId('results-view-details').click();
  const warnings = page.getByTestId('run-diagnostics');
  await expect(warnings).toBeVisible();
  expect(await warnings.evaluate((element) => element.textContent ?? '')).toBe(
    'Warning: [tblastx] Examining 5 or more matches is recommended\nLOSAT query/options error: No hits are being saved\nPlease refer to the LOSAT user manual.\n',
  );
});

// Where TBLASTX's run messages go nowhere (-logfile /dev/null; ABI §4), the run shows no warning: the
// CLI writes none for any format, and the lines that the run's other formats add go nowhere too
// (line WQA's sweep, D3 review r3, 2026-11-19: the tabular formats' warning and error showed). Run
// details says that the run's warnings are not kept (ABI §5, stream 10; D7, 2026-11-20).
test('TBLASTX with -logfile /dev/null shows no warning, as the CLI writes none', async ({ page }) => {
  await program(page, 'tblastx');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT);
  await openParameters(page);
  await page.getByTestId('param-num_alignments').fill('0');
  await page.getByTestId('param-extra').fill('-logfile /dev/null');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-diagnostics-place')).toHaveText("The run's warnings are not kept, as the run's -logfile writes its messages nowhere.");
  await expect(page.getByTestId('run-diagnostics')).toHaveCount(0);
  await expect(page.getByTestId('run-diagnostics-none')).toHaveCount(0);
  if (NATIVE === undefined) return;
  for (const outfmt of ['0', '6']) {
    const cli = spawnSync(NATIVE, ['tblastx', '-query', 'query.fa', '-subject', 'subject.fa', '-num_alignments', '0', '-outfmt', outfmt, '-logfile', '/dev/null', '-num_threads', '1'], { cwd: nativeDirectory(QUERY, SUBJECT) });
    expect(String(cli.stderr), `-outfmt ${outfmt}`).toBe('');
  }
});

// Each output that a run could not write points to where its error is, as the engine says (ABI v2 §5, stream 9: the
// place where it wrote the error; line WQA's sweep, D7, 2026-11-20), and the notice, the Outputs and Run details say
// so. TBLASTX with -logfile run.log: -num_alignments 0, whose hit list of 0 fails the tabular formats, writes BLAST+'s
// error to the Warnings with the rest of what the tabular formats add, and the file of -logfile holds what the CLI's
// file of the argv's outfmt 0 holds, nothing (line WQA's sweep, D3 review r2, 2026-11-19); a subject title that NCBI
// cannot decode stops them with LOSAT's own refusal, which goes to the Warnings, as the CLI writes it to its standard
// error whatever -logfile is (the file stays empty, as the CLI's). With -logfile /dev/null BLAST+'s error is written
// nowhere, and LOSAT's still to the Warnings.
test("TBLASTX with -logfile: the outputs that could not be written point to where the run wrote their error", async ({ page }, testInfo) => {
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the CLI's file: ${NO_NATIVE_REASON}` });
  const notice = (kind: string) => page.locator(`[data-testid="results-notice"][data-kind="${kind}"]`);
  const bad = Buffer.concat([Buffer.from('>Subject_9 bad x'), Buffer.from([0x81]), Buffer.from(`y title\n${dna(7, 300)}\n`)]);
  const refusal = 'subject record 1 (Subject_9) has a title whose non-UTF-8 bytes';
  const nowhere = "not kept, as the run's -logfile writes its messages nowhere";
  const warnings = 'under Warnings in Run details';
  const cases = [
    {
      subject: SUBJECT,
      alignments: '0',
      extra: '-logfile run.log',
      kind: 'no-tabular-hits',
      text: `-num_alignments 0 leaves the hit list of outfmt 6 and 7 empty (the "No hits are being saved" error, ${warnings}), so this run has no HSP records; outfmt 0 is still written.`,
      place: warnings,
      file: '',
      diagnostics: 'No hits are being saved',
      details: undefined,
    },
    {
      subject: bad,
      alignments: undefined,
      extra: '-logfile run.log',
      kind: 'no-tabular-output',
      text: `outfmt 6 and 7 could not be written, so this run has no HSP records. Their error is ${warnings}.`,
      place: warnings,
      file: '',
      diagnostics: refusal,
      details: undefined,
    },
    {
      subject: SUBJECT,
      alignments: '0',
      extra: '-logfile /dev/null',
      kind: 'no-tabular-hits',
      text: `-num_alignments 0 leaves the hit list of outfmt 6 and 7 empty (the "No hits are being saved" error, ${nowhere}), so this run has no HSP records; outfmt 0 is still written.`,
      place: nowhere,
      file: undefined,
      diagnostics: undefined,
      details: `The run's warnings are ${nowhere}.`,
    },
    {
      subject: bad,
      alignments: undefined,
      extra: '-logfile /dev/null',
      kind: 'no-tabular-output',
      text: `outfmt 6 and 7 could not be written, so this run has no HSP records. Their error is ${warnings}.`,
      place: warnings,
      file: undefined,
      diagnostics: refusal,
      details: `The run's other warnings are ${nowhere}.`,
    },
  ] as const;
  for (const [index, { subject, alignments, extra, kind, text, place, file, diagnostics, details }] of cases.entries()) {
    // Each case on a page of its own: its run is run 1.
    if (index > 0) await openSearch(page);
    const run = 1;
    await program(page, 'tblastx');
    await paste(page, 'query', QUERY);
    if (typeof subject === 'string') await paste(page, 'subject', subject);
    else await openFiles(page, 'subject', [{ name: 'subject.fa', text: subject }]);
    await openParameters(page);
    // -num_alignments is a field of the form.
    if (alignments !== undefined) await page.getByTestId('param-num_alignments').fill(alignments);
    await page.getByTestId('param-extra').fill(extra);
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    await waitStatus(page, run, 'completed');
    await chooseRun(page, run);
    await expect(notice(kind)).toHaveText(text);
    await page.getByTestId('results-view-outputs').click();
    await page.getByTestId('format-6').click();
    await expect(page.getByTestId('outputs-not-written')).toHaveText(`outfmt 6 could not be written: its error is ${place}.`);
    await page.getByTestId('format-0').click();
    await expect(page.getByTestId('outputs-not-written')).toHaveCount(0);
    if (file === undefined) {
      await expect(page.getByTestId('format-run-files')).toHaveCount(0);
    } else {
      await page.getByTestId('format-run-files').click();
      const shown = page.getByTestId('run-file-0');
      if (file === '') {
        await expect(shown.getByTestId('run-file-empty')).toHaveText('The file is empty.');
        // The CLI's file of the same argv (outfmt 0, which writes the alignments) is empty too.
        const words = alignments === undefined ? ['tblastx'] : ['tblastx', '-num_alignments', alignments];
        const native = nativeLogFile(words, QUERY, subject);
        if (native !== undefined) expect(native.length, `${words.join(' ')} -logfile run.log: the file`).toBe(0);
        // The CLI's run of outfmt 6 writes LOSAT's refusal to its standard error, none to the file.
        const six = typeof subject === 'string' ? undefined : nativeLogRun(['tblastx', '-outfmt', '6'], QUERY, subject);
        if (six !== undefined) {
          expect(six.stderr.toString('utf8'), 'tblastx -outfmt 6 -logfile run.log: standard error').toContain(refusal);
          expect(six.file.length, 'tblastx -outfmt 6 -logfile run.log: the file').toBe(0);
        }
      } else {
        await expect(shown.getByTestId('run-file-text')).toContainText(file);
      }
    }
    // Run details: the Warnings (stream 3) and where the run's other messages are.
    await page.getByTestId('results-view-details').click();
    if (diagnostics === undefined) await expect(page.getByTestId('run-diagnostics')).toHaveCount(0);
    else await expect(page.getByTestId('run-diagnostics')).toContainText(diagnostics);
    if (details === undefined) await expect(page.getByTestId('run-diagnostics-place')).toHaveCount(0);
    else await expect(page.getByTestId('run-diagnostics-place')).toHaveText(details);
    await expect(page.getByTestId('run-diagnostics-none')).toHaveCount(0);
  }
});

// Run details says where the run's warnings are, as the engine says (ABI v2 §5, stream 10; line WQA's sweep, D7,
// 2026-11-20): TBLASTX writes the messages of its run into the file of -logfile (the warning of -max_target_seqs 1
// here, as the CLI's file holds it), or nowhere with -logfile /dev/null, and with -logfile /dev/fd/2, which the CLI
// opens on a file system without free space, where NCBI writes no message (D11 round 2; before, they were its
// Warnings); with -logfile -, and in BLASTN whatever the
// -logfile, they are its Warnings. "The engine wrote no warnings" is said where the warnings go and none came: a
// TBLASTX search without warnings leaves the file of -logfile empty, as the CLI's.
test("Run details says where the run's warnings are: the file of -logfile, nowhere, or its Warnings", async ({ page, browserName }) => {
  const cases = [
    { id: 'tblastx', extra: '-logfile run.log', warned: true, place: "The run's warnings are in the file of the run's -logfile, run.log, under Outputs.", file: 'Examining 5' },
    { id: 'tblastx', extra: '-logfile /dev/null', warned: true, place: "The run's warnings are not kept, as the run's -logfile writes its messages nowhere.", file: undefined },
    { id: 'tblastx', extra: '-logfile /dev/fd/2', warned: true, place: "The run's warnings are not kept, as the run's -logfile writes its messages nowhere.", file: undefined },
    { id: 'tblastx', extra: '-logfile -', warned: true, place: undefined, file: undefined },
    { id: 'blastn', extra: '-logfile run.log', warned: true, place: undefined, file: '' },
    { id: 'tblastx', extra: '-logfile run.log', warned: false, place: undefined, file: '' },
  ] as const;
  for (const [index, { id, extra, warned, place, file }] of cases.entries()) {
    // Each case on a page of its own: its run is run 1.
    if (index > 0) await openSearch(page);
    const run = 1;
    await program(page, id);
    await paste(page, 'query', QUERY);
    await paste(page, 'subject', SUBJECT);
    await openParameters(page);
    // "Max target sequences" is a field of the form.
    if (warned) {
      await page.getByTestId('param-max_target_seqs').selectOption('other');
      await page.getByTestId('param-max_target_seqs-other').fill('1');
    }
    await page.getByTestId('param-extra').fill(extra);
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    await waitStatus(page, run, 'completed');
    await chooseRun(page, run);
    await page.getByTestId('results-view-details').click();
    const warning = `Warning: [${id}] Examining 5 or more matches is recommended`;
    if (!warned) {
      await expect(page.getByTestId('run-diagnostics-none')).toHaveText('The engine wrote no warnings.');
      await expect(page.getByTestId('run-diagnostics')).toHaveCount(0);
      await expect(page.getByTestId('run-diagnostics-place')).toHaveCount(0);
    } else if (place === undefined) {
      await expect(page.getByTestId('run-diagnostics')).toContainText(warning);
      await expect(page.getByTestId('run-diagnostics-place')).toHaveCount(0);
    } else {
      await expect(page.getByTestId('run-diagnostics-place')).toHaveText(place);
      await expect(page.getByTestId('run-diagnostics')).toHaveCount(0);
    }
    if (warned) await expect(page.getByTestId('run-diagnostics-none')).toHaveCount(0);
    await shoot(page.getByTestId('run-details'), `07-run-details-${id}-${extra.replace(/[^a-z]+/g, '-')}${warned ? '' : '-no-warning'}`, browserName);
    await page.getByTestId('results-view-outputs').click();
    if (file === undefined) await expect(page.getByTestId('format-run-files')).toHaveCount(0);
    else {
      await page.getByTestId('format-run-files').click();
      const shown = page.getByTestId('run-file-0');
      if (file === '') await expect(shown.getByTestId('run-file-empty')).toHaveText('The file is empty.');
      else await expect(shown.getByTestId('run-file-text')).toContainText(warning);
    }
  }
});

// -conffile: the engine's start names the file, the form asks for it, and without it the engine's error
// (the CLI's) keeps the search from starting. A file of any name is read under the argv's name; a file
// that NCBI cannot read is refused with NCBI's report.
test('the file that -conffile names is chosen in the form and read under that name', async ({ page, browserName }, testInfo) => {
  const validation = page.getByTestId('argv-validation');
  const row = page.getByTestId('param-conffile');
  await program(page, 'blastn');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT);
  await openParameters(page);
  await expect(row).toHaveCount(0);
  await page.getByTestId('param-extra').fill('-conffile settings.ini');
  await expect(row).toBeVisible();
  await expect(page.getByTestId('param-conffile-argv')).toHaveText('settings.ini');
  await expect(page.getByTestId('param-conffile-missing')).toHaveText('Choose the file that the search reads as settings.ini.');
  const missing = 'Critical: (106.23) Registry file "settings.ini" not found, aborting.';
  await expect(validation).toHaveAttribute('data-state', 'invalid');
  await expect(validation).toHaveText(missing);
  await shoot(page.getByTestId('parameter-form'), '01-conffile-missing', browserName);
  await submit(page);
  await expect(page.getByTestId('search-message')).toContainText(missing);
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);

  // A file that NCBI's registry reader refuses: NCBI's report.
  await page.getByTestId('param-conffile-file').setInputFiles({ name: 'broken.ini', mimeType: 'text/plain', buffer: Buffer.from('[x\n') });
  await expect(page.getByTestId('param-conffile-name')).toHaveText('broken.ini (3 bytes)');
  await expect(validation).toHaveAttribute('data-state', 'invalid');
  await expect(validation).toContainText("Error: Invalid registry section (']' is missing)");

  // A good file of another name: the search reads it as settings.ini.
  const settings = '[BLAST]\n';
  await page.getByTestId('param-conffile-file').setInputFiles({ name: 'my settings.ini', mimeType: 'text/plain', buffer: Buffer.from(settings) });
  await expect(page.getByTestId('param-conffile-name')).toHaveText('my settings.ini (8 bytes)');
  await expect(row).toHaveClass(/changed/);
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await shoot(page.getByTestId('parameter-form'), '02-conffile-chosen', browserName);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-conffile')).toHaveText('settings.ini: my settings.ini, 8 bytes');
  await expect(page.getByTestId('run-conffile-sha256')).toHaveText(`SHA-256 ${createHash('sha256').update(settings).digest('hex')}`);
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the outputs: ${NO_NATIVE_REASON}` });
  else {
    const { outputs } = await storedOutputs(page, 1);
    const directory = nativeDirectory(QUERY, SUBJECT);
    writeFileSync(join(directory, 'settings.ini'), settings);
    const cli = spawnSync(NATIVE, ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-conffile', 'settings.ini', '-outfmt', '6', '-num_threads', '1'], { cwd: directory });
    expect(cli.status, String(cli.stderr)).toBe(0);
    expect(outputs[6].equals(cli.stdout)).toBe(true);
  }

  // Removed, the file is asked for again.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('param-conffile-remove').click();
  await expect(page.getByTestId('param-conffile-missing')).toBeVisible();
  await expect(validation).toHaveText(missing);
  // Without -conffile the row goes.
  await page.getByTestId('param-extra').fill('');
  await expect(row).toHaveCount(0);
  await expect(validation).toHaveAttribute('data-state', 'ok');
  // "Restore default search parameters" lets the chosen file go with Additional options.
  await page.getByTestId('param-extra').fill('-conffile settings.ini');
  await page.getByTestId('param-conffile-file').setInputFiles({ name: 'my settings.ini', mimeType: 'text/plain', buffer: Buffer.from(settings) });
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await page.getByTestId('restore-defaults').click();
  await expect(page.getByTestId('param-extra')).toHaveValue('');
  await expect(row).toHaveCount(0);
  await page.getByTestId('param-extra').fill('-conffile settings.ini');
  await expect(page.getByTestId('param-conffile-missing')).toBeVisible();
  await expect(validation).toHaveText(missing);
  // `-conffile=` names no file: no row asks for one, and the engine's message stands.
  await page.getByTestId('param-extra').fill('-conffile=');
  await expect(validation).toHaveText('Critical: (106.23) Registry file "" not found, aborting.');
  await expect(row).toHaveCount(0);
});

// -dryrun: the run reads the arguments and the settings and searches nothing (NCBI BLAST+ 2.17.0); the
// results say so in place of "no hits", and a -logfile file of the run is made, empty.
test('a dry run says that it searched nothing', async ({ page, browserName }) => {
  const dryRun = 'This run was a dry run (-dryrun): LOSAT read the arguments and the settings and ran no search.';
  await program(page, 'blastn');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-dryrun -logfile run.log');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  const notice = page.getByTestId('results-notice');
  await expect(notice).toHaveCount(1);
  await expect(notice).toHaveAttribute('data-kind', 'dry-run');
  await expect(notice).toHaveText(dryRun);
  await expect(page.getByTestId('section-graphic')).toHaveCount(0);
  await shoot(page.getByTestId('results-panel'), '05-dry-run', browserName);
  await page.getByTestId('pane-dotplot').click();
  await expect(page.getByTestId('results-dotplot-empty')).toHaveText(dryRun);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-dry-run')).toHaveText('a dry run: no search');
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-run-files').click();
  await expect(page.getByTestId('run-file-0').getByTestId('run-file-empty')).toHaveText('The file is empty.');
});

// The search strategy file of -export_search_strategy (line RG item F1w; ABI v2 §5 stream 8) is a file of the run as the
// file of -logfile is: shown under its option after the formats and saved as the CLI writes it. A name that the CLI
// cannot open is NCBI's error before the run, as on the command line.
test('the file of -export_search_strategy is a file of the run, saved as the CLI writes it', async ({ page, browserName }, testInfo) => {
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the file: ${NO_NATIVE_REASON}` });
  await program(page, 'blastn');
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT);
  await openParameters(page);
  const extra = page.getByTestId('param-extra');
  await extra.fill('-export_search_strategy /etc/s.asn');
  const validation = page.getByTestId('argv-validation');
  await expect(validation).toHaveAttribute('data-state', 'invalid');
  await expect(validation).toContainText('Argument "export_search_strategy". File is not accessible:  `/etc/s.asn\'');
  await extra.fill('-export_search_strategy out/s.asn -logfile run.log');
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-outputs').click();
  const tab = page.getByTestId('format-run-files');
  await expect(tab).toHaveText('-export_search_strategy out/s.asn, -logfile run.log');
  await tab.click();
  const file = page.getByTestId('run-file-0');
  await expect(file.getByTestId('run-file-name')).toHaveText('out/s.asn');
  await expect(file.getByTestId('run-file-text')).toContainText('Blast4-request ::= {');
  await shoot(page.getByTestId('results-outputs'), '09-outputs-strategy', browserName);
  const saved = await downloaded(page, () => page.getByTestId('run-file-download-0').click());
  expect(saved.name).toMatch(/^losat-run1-blastn-.*Z\.export_search_strategy\.out_s\.asn$/);
  await page.getByTestId('results-download-all').click();
  const item = page.getByTestId('results-download-menu').getByTestId('results-download-run-file-0');
  await expect(item).toHaveText('out/s.asn (-export_search_strategy)');
  expect(await downloaded(page, () => item.click())).toEqual(saved);
  if (NATIVE !== undefined) {
    const directory = nativeDirectory(QUERY, SUBJECT);
    mkdirSync(join(directory, 'out'));
    const cli = spawnSync(NATIVE, ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-export_search_strategy', 'out/s.asn', '-num_threads', '1'], { cwd: directory });
    expect(cli.status, cli.stderr.toString('utf8')).toBe(0);
    expect(saved.bytes.equals(readFileSync(join(directory, 'out/s.asn'))), "blastn -export_search_strategy out/s.asn: the CLI's file").toBe(true);
  }
});
