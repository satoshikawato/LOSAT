// BLASTN's -window_masker_taxid through the real application (line NP WW2, DESIGN_ww2 §4-§5): the
// "Species-specific repeats filter for:" row with a taxonomy ID, and the line under the box that
// downloads the counts of the ID from ftp.ncbi.nlm.nih.gov only when the user presses "Download",
// keeps them in this browser (OPFS windowmasker/<taxid>/), continues a download that stopped with a
// Range, and starts again when the file changed on the site. The site is stubbed with routes at
// its real address (taxonomy ID 562 serves LOSAT/tests/winmask/counts/mg1655_u9.obinary), so the
// page runs under the deployed CSP (public/_headers connect-src); the tests count every request
// that leaves the site and never contact NCBI. In the engine build with the native CLI
// (LOSAT_WEB_REACTORS and LOSAT_WEB_NATIVE), a search with the downloaded counts gives the bytes of
// the native CLI that reads <dir>/562/wmasker.obinary (WINDOW_MASKER_PATH=<dir>). With
// LOSAT_WEB_SCREENS=<directory>, the row's states are saved as PNGs.
import { copyFileSync, mkdirSync, readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { join } from 'node:path';
import { expect, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, nativeExpectation, NO_NATIVE_REASON } from './support/native';
import { E2E_ORIGIN } from './support/origin';
import { reload } from './support/navigate';
import { chooseRun, described, openFiles, openParameters, openSearch, program, storedOutputs, submit, task, waitStatus } from './support/search';

const COUNTS_PATH = 'LOSAT/tests/winmask/counts/mg1655_u9.obinary';
const QUERY = 'LOSAT/tests/fixtures/blastn_regression/inputs/wm_query.fa';
const SUBJECT = 'LOSAT/tests/fixtures/blastn_regression/inputs/s60k.fa';
const FLAG = 'window_masker_db';
const TAXID = 562;
const URL = `https://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/${TAXID}/wmasker.obinary`;
const LAST_MODIFIED = 'Thu, 01 Oct 2026 08:00:00 GMT';
const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;

test.setTimeout(300_000);

const read = (path: string): Buffer => readFileSync(join(REPOSITORY, path));
const baseName = (path: string): string => path.slice(path.lastIndexOf('/') + 1);
const sha256 = (bytes: Uint8Array | string): string => createHash('sha256').update(bytes).digest('hex');
const COUNTS = read(COUNTS_PATH);

async function shoot(target: Locator | Page, name: string, browserName: string): Promise<void> {
  if (SCREENS === undefined) return;
  mkdirSync(join(SCREENS, browserName), { recursive: true });
  await target.screenshot({ path: join(SCREENS, browserName, `${name}.png`) });
}

/** The file the stub serves, and how it answers: a cut ends the body after that many bytes; a gate holds the answer. */
interface Site {
  file: Buffer;
  lastModified: string;
  /** Ends the next answer's body early (the Content-Length says more than comes). */
  cutAt: number | undefined;
  /** The status the next answer has instead of the file. */
  status: number | undefined;
  /** Holds the next answer until it is opened. */
  gate: Promise<void> | undefined;
  /** The Range header of each request that reached the site. */
  readonly ranges: Array<string | undefined>;
  /** Every request that left the application's origin. */
  readonly away: string[];
}

/** Stubs the site of the counts at its real address (a Range is answered as NCBI's FTP does). */
async function stubSite(page: Page, options: Partial<Pick<Site, 'file' | 'lastModified'>> = {}): Promise<Site> {
  const site: Site = { file: options.file ?? COUNTS, lastModified: options.lastModified ?? LAST_MODIFIED, cutAt: undefined, status: undefined, gate: undefined, ranges: [], away: [] };
  page.on('request', (request) => {
    const url = request.url();
    if (/^https?:/.test(url) && !url.startsWith(E2E_ORIGIN)) site.away.push(url);
  });
  await page.route(/^https:\/\/ftp\.ncbi\.nlm\.nih\.gov\//, async (route) => {
    const request = route.request();
    const cors = { 'Access-Control-Allow-Origin': '*', 'Access-Control-Expose-Headers': 'Content-Range, Content-Length, Last-Modified' };
    if (request.url() !== URL) {
      await route.fulfill({ status: 404, headers: cors, body: 'Not Found' });
      return;
    }
    const range = request.headers()['range'];
    site.ranges.push(range);
    const gate = site.gate;
    site.gate = undefined;
    await gate;
    const status = site.status;
    site.status = undefined;
    if (status !== undefined) {
      await route.fulfill({ status, headers: cors, body: 'Error' });
      return;
    }
    const match = range === undefined ? null : /^bytes=(\d+)-$/.exec(range);
    const from = match === null ? 0 : Number(match[1]);
    const rest = site.file.subarray(from);
    const cut = site.cutAt;
    site.cutAt = undefined;
    const headers: Record<string, string> = {
      ...cors,
      'Content-Type': 'application/octet-stream',
      'Content-Length': String(rest.length),
      'Last-Modified': site.lastModified,
    };
    if (match !== null) headers['Content-Range'] = `bytes ${from}-${site.file.length - 1}/${site.file.length}`;
    await route.fulfill({ status: match === null ? 200 : 206, headers, body: cut === undefined ? rest : rest.subarray(0, cut) });
  });
  return site;
}

/** Checks "Species-specific repeats filter for:" (the taxonomy ID is the default) and types the ID. */
async function typeTaxid(page: Page, id: number | string = TAXID): Promise<void> {
  const checkbox = page.getByTestId(`param-${FLAG}`);
  if (!(await checkbox.isChecked())) await checkbox.check();
  await page.getByTestId(`param-${FLAG}-taxid`).fill(String(id));
}

const note = (page: Page) => page.getByTestId('counts-note');
/** The sentence of the saved counts: kept in the browser, or (where the browser keeps nothing across sessions) in this tab. */
const SAVED = /^(?:Saved in this browser|Kept in this tab only) \(2026-10-01, 568\.6 kB\)$/;

/** The inputs of the search: the query and the subject opened as files, the task blastn. */
async function setUpSearch(page: Page): Promise<void> {
  await openFiles(page, 'query', [{ name: baseName(QUERY), text: read(QUERY) }]);
  await openFiles(page, 'subject', [{ name: baseName(SUBJECT), text: read(SUBJECT) }]);
  await openParameters(page);
  await task(page, 'blastn');
}

async function openBlastn(page: Page): Promise<void> {
  await openSearch(page);
  await program(page, 'blastn');
  await openParameters(page);
}

test('the row: a taxonomy ID by default, a select for a counts file, no request before Download, and the option refused in Additional options', async ({ page, browserName }) => {
  const site = await stubSite(page);
  await openBlastn(page);
  const section = page.getByRole('group', { name: 'Filters and Masking' });
  const checkbox = page.getByTestId(`param-${FLAG}`);
  await expect(page.getByTestId(`param-${FLAG}-by`)).toHaveCount(0);
  await shoot(section, '00-row-unchecked', browserName);

  // Checked: the taxonomy ID is asked for; the search does not start.
  await checkbox.check();
  const by = page.getByTestId(`param-${FLAG}-by`);
  await expect(by).toHaveValue('taxid');
  await expect(by.locator('option')).toHaveText(['Taxonomy ID', 'Counts file']);
  await expect(page.getByTestId(`param-${FLAG}-missing`)).toHaveText('Enter a taxonomy ID.');
  await expect(page.getByTestId('counts-line')).toHaveCount(0);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveCount(0);
  await shoot(section, '01-taxid-empty', browserName);

  // Typed: the line that downloads, with its in-place sentence; nothing was asked of the site.
  await page.getByTestId(`param-${FLAG}-taxid`).fill('562');
  await expect(note(page)).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');
  await expect(page.getByTestId('counts-download')).toBeVisible();
  await expect(page.getByTestId(`param-${FLAG}-field`)).toHaveClass(/changed/);
  await expect(page.getByTestId('algorithm-parameters-changed')).toHaveText('(1 changed)');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await shoot(section, '02-not-downloaded', browserName);
  expect(site.away).toEqual([]);
  expect(site.ranges).toEqual([]);

  // The box takes digits only; an ID of 0 is not one.
  await page.getByTestId(`param-${FLAG}-taxid`).fill('5a6x2');
  await expect(page.getByTestId(`param-${FLAG}-taxid`)).toHaveValue('562');
  await page.getByTestId(`param-${FLAG}-taxid`).fill('0');
  await expect(page.getByTestId(`param-${FLAG}-missing`)).toHaveText('Enter a taxonomy ID.');
  await page.getByTestId(`param-${FLAG}-taxid`).fill('562');

  // A search needs the counts: the line above Run says so, and Run says it too.
  await openFiles(page, 'query', [{ name: baseName(QUERY), text: read(QUERY) }]);
  await openFiles(page, 'subject', [{ name: baseName(SUBJECT), text: read(SUBJECT) }]);
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Download the counts for 562.');
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);

  // "Additional options" does not take the option twice.
  await page.getByTestId('param-extra').fill('-window_masker_taxid 9606');
  await expect(page.getByTestId('argv-validation')).toContainText('-window_masker_taxid is set in the form: Species-specific repeats filter for:');
  await page.getByTestId('param-extra').fill('');

  // The other way: a counts file of the user's own (WW1), which lets the ID go.
  await by.selectOption('file');
  await expect(page.getByTestId(`param-${FLAG}-taxid`)).toHaveCount(0);
  await expect(page.getByTestId('counts-line')).toHaveCount(0);
  await expect(page.getByTestId(`param-${FLAG}-missing`)).toHaveText('Choose the WindowMasker counts file.');
  await shoot(section, '03-counts-file', browserName);
  await by.selectOption('taxid');
  await expect(page.getByTestId(`param-${FLAG}-missing`)).toHaveText('Enter a taxonomy ID.');

  // Restore default search parameters, and unchecking, let the ID go.
  await page.getByTestId(`param-${FLAG}-taxid`).fill('562');
  await page.getByTestId('restore-defaults').click();
  await expect(checkbox).not.toBeChecked();
  await expect(page.getByTestId('counts-line')).toHaveCount(0);
  expect(site.away).toEqual([]);
});

test('Download keeps the counts in this browser; a search reads them; they load again without a request; Delete removes them', async ({ page, browserName }) => {
  const site = await stubSite(page);
  await openBlastn(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  const section = page.getByRole('group', { name: 'Filters and Masking' });
  await typeTaxid(page);
  await expect(note(page)).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');

  // Nothing is asked until Download is pressed; the request is for the file's URL, with no Range.
  await page.getByTestId('counts-download').click();
  await expect(note(page)).toHaveText(SAVED);
  expect(site.away).toEqual([URL]);
  expect(site.ranges).toEqual([undefined]);
  await expect(page.getByTestId('counts-delete')).toBeVisible();
  await shoot(section, '05-saved', browserName);

  // The search reads the saved counts.
  await setUpSearch(page);
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-argv')).toContainText(`-window_masker_taxid ${TAXID}`);
  const line = page.getByTestId('run-window-masker');
  await expect(line).toContainText(`${TAXID}: wmasker.obinary from ftp.ncbi.nlm.nih.gov (2026-10-01), 568.6 kB,`);
  await expect(line).toContainText(`SHA-256 ${sha256(COUNTS)}`);
  await expect(page.getByTestId('run-window-masker-url')).toHaveText(URL);
  await shoot(page.getByTestId('run-record'), '06-run-details', browserName);
  expect(site.away).toEqual([URL]);

  if (backend !== 'opfs') return; // Nothing lasts without OPFS: the next load downloads again.

  // After a reload the counts are found in the browser: no request, and the search runs.
  await reload(page);
  await described(page);
  await program(page, 'blastn');
  await openParameters(page);
  await typeTaxid(page);
  await expect(note(page)).toHaveText(SAVED);
  await setUpSearch(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-window-masker')).toContainText(`SHA-256 ${sha256(COUNTS)}`);
  expect(site.away).toEqual([URL]);

  // Delete removes them: the sentence returns, and a search needs the download again.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('counts-delete').click();
  await expect(note(page)).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');
  await expect(page.getByTestId('counts-delete')).toHaveCount(0);
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Download the counts for 562.');
  const storage = await page.evaluate(async () => {
    const root = await navigator.storage.getDirectory();
    try {
      const directory = await root.getDirectoryHandle('windowmasker');
      const names: string[] = [];
      for await (const name of (directory as unknown as { keys(): AsyncIterable<string> }).keys()) names.push(name);
      return names;
    } catch {
      return [];
    }
  });
  expect(storage).toEqual([]);
});

/** Cuts the first download after `bytes` bytes of the body: the connection ended early. */
async function cutFirstDownload(page: Page, site: Site, bytes: number): Promise<void> {
  await typeTaxid(page);
  site.cutAt = bytes;
  await page.getByTestId('counts-download').click();
}

/** Where the browser keeps nothing across sessions, a download that stopped leaves nothing to continue: the error, and Download again. */
async function expectNothingKept(page: Page): Promise<void> {
  await expect(page.getByTestId('counts-error')).toContainText('ftp.ncbi.nlm.nih.gov could not be reached');
  await expect(page.getByTestId('counts-download')).toBeVisible();
  await expect(page.getByTestId('counts-continue')).toHaveCount(0);
  await expect(note(page)).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');
}

test('a download that stopped keeps what came; Continue asks for the rest with a Range and the file is the same', async ({ page, browserName }) => {
  const site = await stubSite(page);
  await openBlastn(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  const section = page.getByRole('group', { name: 'Filters and Masking' });
  await cutFirstDownload(page, site, 200_000);
  if (backend !== 'opfs') return expectNothingKept(page);
  await expect(page.getByTestId('counts-error')).toContainText('ftp.ncbi.nlm.nih.gov could not be reached');
  const partial = await note(page).textContent();
  expect(partial).toMatch(/^Downloaded \d+\.\d kB of 568\.6 kB$/);
  await expect(page.getByTestId('counts-continue')).toBeVisible();
  await expect(page.getByTestId('counts-delete')).toBeVisible();
  await shoot(section, '04-partial', browserName);

  // The partial copy outlives a reload.
  await reload(page);
  await described(page);
  await program(page, 'blastn');
  await openParameters(page);
  await typeTaxid(page);
  await expect(note(page)).toHaveText(partial!);
  expect(site.ranges).toEqual([undefined]);

  // Continue: held, to show the progress; then the rest comes by Range from the saved byte.
  const saved = Number(/^Downloaded ([\d.]+) kB/.exec(partial!)![1]) * 1000;
  let open!: () => void;
  site.gate = new Promise<void>((resolve) => (open = resolve));
  await page.getByTestId('counts-continue').click();
  await expect(note(page)).toHaveText(/^Downloading \d+\.\d kB of 568\.6 kB from ftp\.ncbi\.nlm\.nih\.gov…$/);
  await expect(page.getByTestId('counts-cancel')).toBeVisible();
  await shoot(section, '03b-downloading', browserName);
  open();
  await expect(note(page)).toHaveText('Saved in this browser (2026-10-01, 568.6 kB)');
  expect(site.ranges).toHaveLength(2);
  expect(site.ranges[1]).toMatch(/^bytes=\d+-$/);
  expect(Math.abs(Number(site.ranges[1]!.slice(6, -1)) - saved)).toBeLessThan(100);

  await setUpSearch(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-window-masker')).toContainText(`SHA-256 ${sha256(COUNTS)}`);
});

test('Cancel stops the download and keeps what came', async ({ page }) => {
  const site = await stubSite(page);
  await openBlastn(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  await cutFirstDownload(page, site, 200_000);
  if (backend !== 'opfs') return expectNothingKept(page);
  await expect(page.getByTestId('counts-continue')).toBeVisible();
  let open!: () => void;
  site.gate = new Promise<void>((resolve) => (open = resolve));
  await page.getByTestId('counts-continue').click();
  await page.getByTestId('counts-cancel').click();
  open();
  await expect(note(page)).toHaveText(/^Downloaded \d+\.\d kB of 568\.6 kB$/);
  await expect(page.getByTestId('counts-continue')).toBeVisible();
  await page.getByTestId('counts-delete').click();
  await expect(note(page)).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');
});

test('a file that changed on the site (another Last-Modified) is downloaded again from its first byte', async ({ page }) => {
  const site = await stubSite(page);
  await openBlastn(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  await cutFirstDownload(page, site, 200_000);
  if (backend !== 'opfs') return expectNothingKept(page);
  await expect(page.getByTestId('counts-continue')).toBeVisible();
  // The site has another file now.
  const changed = Buffer.concat([COUNTS, Buffer.from([1, 2, 3])]);
  site.file = changed;
  site.lastModified = 'Fri, 02 Oct 2026 08:00:00 GMT';
  await page.getByTestId('counts-continue').click();
  await expect(note(page)).toHaveText('Saved in this browser (2026-10-02, 568.6 kB)');
  // The Range was refused by the other file's date; the file was asked again whole.
  expect(site.ranges.length).toBe(3);
  expect(site.ranges[1]).toMatch(/^bytes=\d+-$/);
  expect(site.ranges[2]).toBeUndefined();
  await setUpSearch(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-window-masker')).toContainText(`SHA-256 ${sha256(changed)}`);
  await expect(page.getByTestId('run-window-masker')).toContainText('(2026-10-02)');
});

test('a site with no counts for the ID, and one that fails, are named and leave the form as it was', async ({ page }) => {
  const site = await stubSite(page);
  await openBlastn(page);
  await typeTaxid(page, 12345);
  await page.getByTestId('counts-download').click();
  await expect(page.getByTestId('counts-error')).toHaveText('ftp.ncbi.nlm.nih.gov has no WindowMasker counts for 12345.');
  await expect(page.getByTestId('counts-download')).toBeVisible();
  await typeTaxid(page, TAXID);
  site.status = 500;
  await page.getByTestId('counts-download').click();
  await expect(page.getByTestId('counts-error')).toHaveText('ftp.ncbi.nlm.nih.gov answered HTTP 500: try again in a moment.');
  await expect(note(page)).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');
  // Without a retry of its own: one request each.
  expect(site.away).toEqual(['https://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/12345/wmasker.obinary', URL]);
  await page.getByTestId('counts-download').click();
  await expect(note(page)).toHaveText(SAVED);
});

test.describe('the engine and the native CLI', () => {
  test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);

  test('a search with the downloaded counts gives the bytes of the native CLI that reads <dir>/562/wmasker.obinary; a Continue gives the same bytes', async ({ page }, testInfo) => {
    // The native CLI runs the argv that the application builds, with the counts where WINDOW_MASKER_PATH names.
    const directory = testInfo.outputPath('native');
    mkdirSync(join(directory, 'wm', String(TAXID)), { recursive: true });
    copyFileSync(join(REPOSITORY, COUNTS_PATH), join(directory, 'wm', String(TAXID), 'wmasker.obinary'));
    for (const path of [QUERY, SUBJECT]) copyFileSync(join(REPOSITORY, path), join(directory, baseName(path)));
    const native = nativeExpectation(
      ['blastn', '-query', baseName(QUERY), '-subject', baseName(SUBJECT), '-task', 'blastn', '-window_masker_taxid', String(TAXID)],
      directory,
      { WINDOW_MASKER_PATH: join(directory, 'wm') },
    );

    const site = await stubSite(page);
    await openBlastn(page);
    const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
    // A download that stopped and was continued (where the browser keeps it), or one download.
    if (backend === 'opfs') {
      await cutFirstDownload(page, site, 200_000);
      await expect(page.getByTestId('counts-continue')).toBeVisible();
      await page.getByTestId('counts-continue').click();
    } else {
      await typeTaxid(page);
      await page.getByTestId('counts-download').click();
    }
    await expect(note(page)).toHaveText(SAVED);
    await setUpSearch(page);
    const threads = page.getByTestId('threads');
    await threads.selectOption('1');
    for (const [run, expected] of [[1, 'read'], [2, 'kept']] as const) {
      await submit(page);
      await waitStatus(page, run, 'completed');
      const { outputs, warnings } = await storedOutputs(page, run);
      expect({ run, 0: sha256(outputs[0]), 6: sha256(outputs[6]), 7: sha256(outputs[7]) }).toEqual({ run, 0: native.sha256[0], 6: native.sha256[6], 7: native.sha256[7] });
      expect(sha256(warnings), `run ${run}: the warnings`).toBe(native.stderrSha256);
      const line = page.getByTestId('run-window-masker');
      await expect(line).toContainText(`SHA-256 ${sha256(COUNTS)}`);
      await expect(line).toContainText(expected === 'kept' ? 'kept from the previous search' : 'read for this search in');
      await page.getByTestId('tab-search').click();
    }
  });
});
