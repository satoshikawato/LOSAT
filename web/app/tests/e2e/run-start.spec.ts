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
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Download, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, openParameters, openSearch, paste, pickRun, program, storedOutputs, submit, waitStatus } from './support/search';

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
function nativeDirectory(query: string, subject: string): string {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, 'query.fa'), query);
  writeFileSync(join(directory, 'subject.fa'), subject);
  return directory;
}

/**
 * The file of `-logfile run.log` that the native CLI writes for the search `words` (program first) on
 * the inputs, with the run's one thread: undefined without LOSAT_WEB_NATIVE.
 */
function nativeLogFile(words: readonly string[], query: string, subject: string): Buffer | undefined {
  if (NATIVE === undefined) return undefined;
  const directory = nativeDirectory(query, subject);
  spawnSync(NATIVE, [words[0]!, '-query', 'query.fa', '-subject', 'subject.fa', ...words.slice(1), '-logfile', 'run.log', '-num_threads', '1'], { cwd: directory });
  return readFileSync(join(directory, 'run.log'));
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
  await expect(status).toContainText('A failed run keeps no results. It keeps the file of its -logfile, run.log.');
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
  await expect(validation).toContainText("Error: (CRegistryException::eSection) Invalid registry section (']' is missing)");

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
