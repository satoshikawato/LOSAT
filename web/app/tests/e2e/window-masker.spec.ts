// BLASTN's -window_masker_db through the real application (line NP WW1, DESIGN_ww1 §7): the
// "Species-specific repeats" field of "Filters and Masking" with its counts file (every build),
// and, in the engine build with the native CLI (LOSAT_WEB_REACTORS and LOSAT_WEB_NATIVE),
// searches that the application runs with the file: the stored outfmt 0, 6 and 7 and the
// warnings are the bytes of the native CLI of the same commit, and NCBI's frozen outfmt 6;
// the engine keeps the counts for the next BLASTN search of its instance and gives them up when
// the program changes. With LOSAT_WEB_SCREENS=<directory>, the form's states are saved as PNGs.
import { copyFileSync, mkdirSync, readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { expect, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, nativeExpectation, NO_NATIVE_REASON } from './support/native';
import { chooseRun, openFiles, openParameters, openSearch, program, storedOutputs, submit, task, waitStatus } from './support/search';

const COUNTS = 'LOSAT/tests/winmask/counts/mg1655_u9.obinary';
const QUERY = 'LOSAT/tests/fixtures/blastn_regression/inputs/wm_query.fa';
const SUBJECT = 'LOSAT/tests/fixtures/blastn_regression/inputs/s60k.fa';
const FLAG = 'window_masker_db';
/** NCBI BLAST+ 2.17.0's outfmt 6 of this search (row wm.blastn of LOSAT/tests/fixtures/blastn_regression/manifest.tsv). */
const NCBI_OUTFMT6 = '7774c0f582c486fd1f93d3412fd2b508937139e97005a158dd00a632b20f17cb';
const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;

test.setTimeout(300_000);

const read = (path: string): Buffer => readFileSync(join(REPOSITORY, path));
const baseName = (path: string): string => path.slice(path.lastIndexOf('/') + 1);
const sha256 = (bytes: Uint8Array | string): string => createHash('sha256').update(bytes).digest('hex');

async function shoot(target: Locator | Page, name: string, browserName: string): Promise<void> {
  if (SCREENS === undefined) return;
  mkdirSync(join(SCREENS, browserName), { recursive: true });
  await target.screenshot({ path: join(SCREENS, browserName, `${name}.png`) });
}

/** The inputs of the search: the query and the subject opened as files, the task blastn. */
async function setUpSearch(page: Page): Promise<void> {
  await openFiles(page, 'query', [{ name: baseName(QUERY), text: read(QUERY) }]);
  await openFiles(page, 'subject', [{ name: baseName(SUBJECT), text: read(SUBJECT) }]);
  await openParameters(page);
  await task(page, 'blastn');
}

/** Checks "Species-specific repeats filter for:", names the counts by a file (WW2: a taxonomy ID is the default) and chooses it. */
async function chooseCounts(page: Page): Promise<void> {
  await page.getByTestId(`param-${FLAG}`).check();
  await page.getByTestId(`param-${FLAG}-by`).selectOption('file');
  await page.getByTestId(`param-${FLAG}-file`).setInputFiles({ name: baseName(COUNTS), mimeType: 'application/octet-stream', buffer: read(COUNTS) });
  await expect(page.getByTestId(`param-${FLAG}-name`)).toContainText(baseName(COUNTS));
}

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

test('the field: a box and a file in the Filter row, an empty box keeps the search from starting, the file is an argument', async ({ page, browserName }) => {
  await program(page, 'blastn');
  await setUpSearch(page);
  const section = page.getByRole('group', { name: 'Filters and Masking' });
  const field = page.getByTestId(`param-${FLAG}-field`);
  const checkbox = page.getByTestId(`param-${FLAG}`);
  // NCBI's Filter row: the low complexity box and the species-specific repeats box, side by side.
  await expect(section.getByTestId('param-dust-field')).toBeVisible();
  await expect(field).toContainText('Species-specific repeats filter for:');
  await expect(checkbox).not.toBeChecked();
  await expect(page.getByTestId(`param-${FLAG}-choose`)).toHaveCount(0);
  await shoot(section, '01-unchecked', browserName);

  // Checked without a file: the sentence that asks for it, no argument, and Run does not start.
  await checkbox.check();
  await page.getByTestId(`param-${FLAG}-by`).selectOption('file');
  await expect(page.getByTestId(`param-${FLAG}-missing`)).toHaveText('Choose the WindowMasker counts file.');
  await expect(field).not.toHaveClass(/changed/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  await expect(page.getByTestId('param-extra')).toBeVisible();
  await shoot(section, '02-checked-no-file', browserName);
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Choose the WindowMasker counts file.');
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);
  await shoot(page.getByTestId('search-message'), '02b-run-refused', browserName);

  // The file chosen: its name and size, the field marked as changed, the name in the argument.
  await page.getByTestId(`param-${FLAG}-file`).setInputFiles({ name: baseName(COUNTS), mimeType: 'application/octet-stream', buffer: read(COUNTS) });
  await expect(page.getByTestId(`param-${FLAG}-name`)).toHaveText(`WindowMasker counts: ${baseName(COUNTS)} (568.6 kB)`);
  await expect(page.getByTestId(`param-${FLAG}-missing`)).toHaveCount(0);
  await expect(field).toHaveClass(/changed/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(1 changed)');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await shoot(section, '03-file-chosen', browserName);

  // "Additional options" does not take the option twice.
  await page.getByTestId('param-extra').fill('-window_masker_db other.obinary');
  await expect(page.getByTestId('argv-validation')).toContainText('-window_masker_db is set in the form: Species-specific repeats filter for:');
  await page.getByTestId('param-extra').fill('');

  // Restore default search parameters lets the file go; unchecking does too.
  await page.getByTestId('restore-defaults').click();
  await expect(checkbox).not.toBeChecked();
  await expect(page.getByTestId(`param-${FLAG}-name`)).toHaveCount(0);
  await task(page, 'blastn');
  await chooseCounts(page);
  await checkbox.uncheck();
  await expect(page.getByTestId(`param-${FLAG}-name`)).toHaveCount(0);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  // The same file can be chosen again.
  await chooseCounts(page);
});

test('a search with the counts file records the file and its SHA-256, and puts the argument in the command', async ({ page, browserName }) => {
  await program(page, 'blastn');
  await setUpSearch(page);
  await chooseCounts(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  const argv = page.getByTestId('run-argv');
  await expect(argv).toContainText(`-window_masker_db ${baseName(COUNTS)}`);
  const line = page.getByTestId('run-window-masker');
  await expect(line).toContainText(`${baseName(COUNTS)}, 568.6 kB`);
  await expect(line).toContainText(`SHA-256 ${sha256(read(COUNTS))}`);
  await expect(line).toContainText(BUILD_HAS_ENGINE ? 'read for this search in' : 'read for this search');
  await shoot(page.getByTestId('run-record'), '04-run-details', browserName);
});

test.describe('the engine and the native CLI', () => {
  test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);

  test('searches with the counts file give the bytes of the native CLI and of NCBI; the engine keeps the counts and gives them up with the program', async ({ page }, testInfo) => {
    // The native CLI runs the argv that the application builds, with the inputs named alike.
    const directory = testInfo.outputPath('native');
    mkdirSync(directory, { recursive: true });
    for (const path of [COUNTS, QUERY, SUBJECT]) copyFileSync(join(REPOSITORY, path), join(directory, baseName(path)));
    const native = nativeExpectation(
      ['blastn', '-query', baseName(QUERY), '-subject', baseName(SUBJECT), '-task', 'blastn', '-window_masker_db', baseName(COUNTS)],
      directory,
    );
    expect(native.sha256[6]).toBe(NCBI_OUTFMT6);

    await program(page, 'blastn');
    await setUpSearch(page);
    await chooseCounts(page);
    const threads = page.getByTestId('threads');
    const details = page.getByTestId('run-record');

    /** Searches once more and checks the run against the native CLI; returns its counts line, if it has one. */
    async function search(run: number, expected: 'read' | 'kept' | 'none'): Promise<void> {
      await submit(page);
      await waitStatus(page, run, 'completed');
      const { outputs, warnings } = await storedOutputs(page, run);
      if (expected !== 'none') {
        expect({ run, 0: sha256(outputs[0]), 6: sha256(outputs[6]), 7: sha256(outputs[7]) }).toEqual({ run, 0: native.sha256[0], 6: native.sha256[6], 7: native.sha256[7] });
        expect(sha256(warnings), `run ${run}: the warnings`).toBe(native.stderrSha256);
      }
      await page.getByTestId('results-view-details').click();
      const line = page.getByTestId('run-window-masker');
      if (expected === 'none') await expect(line).toHaveCount(0);
      else {
        await expect(line).toContainText(`SHA-256 ${sha256(read(COUNTS))}`);
        await expect(line).toContainText(expected === 'kept' ? 'kept from the previous search' : 'read for this search in');
      }
      await expect(details).toBeVisible();
      await page.getByTestId('tab-search').click();
    }

    await threads.selectOption('1');
    await search(1, 'read');
    // The engine holds the counts: the next search does not read the file again.
    await search(2, 'kept');
    // The threaded instance is another instance: it reads the file.
    await threads.selectOption('2');
    await search(3, 'read');
    await search(4, 'kept');
    // Another program: the counts are given up before it runs.
    await program(page, 'tblastx');
    await search(5, 'none');
    await program(page, 'blastn');
    await search(6, 'read');
    // And back to one thread, a new serial instance.
    await threads.selectOption('1');
    await search(7, 'read');
  });
});
