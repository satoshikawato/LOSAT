// Line WF-2 (docs/web/fetch_design.md §9): "Example" downloads a mode's subject set only when the
// user presses it, puts the set in the subject list and the mode's query in the query box, and
// keeps the set in this browser (OPFS examples/), so that the next time it loads without a
// request. The download hosts are stubbed with routes at their real addresses, so the page runs
// under the deployed CSP (public/_headers connect-src); the tests count every request that
// leaves the site. Line WF-3: the saved copy goes to a folder that the user picks (Chromium; the
// picker is stubbed with a folder of OPFS) or to the browser's download (Firefox). The engine's
// decoders inflate a gzip body in the Data worker (swebok-improve SW8 D2): the engine build gets the
// set as two gzip members (as `cat a.gz b.gz` writes it, which the browsers' DecompressionStream
// refused), and a build without the engine, which has no decoder, gets it already decoded on the way.
import { readFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { exampleById, exampleQuery } from '../../src/domain/examples';
import { BUILD_HAS_ENGINE } from './support/browser';
import { reload } from './support/navigate';
import { E2E_ORIGIN } from './support/origin';
import { openSearch, program, settled } from './support/search';

const SIXTEEN_S = exampleById('refseq-16s-bacteria');
const FIRST = '>NR_102783.2 Bacillus subtilis 16S\nACGTACGTACGTACGTACGT\n';
const SET = `${FIRST}>NR_114042.1 Escherichia coli 16S\nACGTTCGTACGTACGTTCGT\n`;
const HOSTS = /^https:\/\/(?:ftp\.ncbi\.nlm\.nih\.gov|rest\.uniprot\.org)\//;

/**
 * Stubs the download hosts (the 16S set as gzip, or the bytes given; UniProt with a 500), and
 * returns the requests that leave the site.
 */
async function stubHosts(page: Page, body?: Buffer): Promise<string[]> {
  const away: string[] = [];
  page.on('request', (request) => {
    const url = request.url();
    if (/^https?:/.test(url) && !url.startsWith(E2E_ORIGIN)) away.push(url);
  });
  const gz = body ?? (BUILD_HAS_ENGINE ? Buffer.concat([gzipSync(FIRST), gzipSync(SET.slice(FIRST.length))]) : Buffer.from(SET));
  await page.route(HOSTS, async (route) => {
    const url = route.request().url();
    if (url === SIXTEEN_S.url) {
      await route.fulfill({
        status: 200,
        contentType: 'application/x-gzip',
        headers: {
          'Access-Control-Allow-Origin': '*',
          'Content-Length': String(gz.length),
          'Last-Modified': 'Tue, 27 Oct 2026 12:04:14 GMT',
        },
        body: gz,
      });
    } else {
      await route.fulfill({ status: 500, headers: { 'Access-Control-Allow-Origin': '*' }, body: 'Error messages\nInternal server error' });
    }
  });
  return away;
}

test('an Example downloads its set only when asked; the saved copy loads again without a request; Delete removes it', async ({
  page,
}) => {
  const away = await stubHosts(page);
  await openSearch(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  console.log(`examples on the ${backend} store`);
  const note = page.getByTestId('example-refseq-16s-bacteria-note');
  await expect(note).toHaveText('Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov');
  await expect(page.getByTestId('example-refseq-its-fungi-note')).toHaveText('Downloads 3.0 MB from ftp.ncbi.nlm.nih.gov');
  await page.getByTestId('query-input').fill('>old\nACGT\n');
  expect(away).toEqual([]);

  await page.getByTestId('example-refseq-16s-bacteria').click();
  await expect(page.getByTestId('query-input')).toHaveValue(exampleQuery(SIXTEEN_S, 'blastn'));
  await expect(page.getByTestId('subject-source-0')).toContainText('bacteria.16SrRNA.fna');
  await settled(page, 'subject');
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('2 records · 40 nt');
  expect(away).toEqual([SIXTEEN_S.url]);
  if (backend !== 'opfs') {
    // Nothing lasts without OPFS: the sentence stays, and the next load downloads again.
    await expect(note).toHaveText('Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov');
    return;
  }
  await expect(note).toHaveText('Saved in this browser (2026-10-27, 111 bytes)');

  await reload(page);
  await expect(note).toHaveText('Saved in this browser (2026-10-27, 111 bytes)');
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await settled(page, 'subject');
  await expect(page.getByTestId('subject-source-0')).toContainText('bacteria.16SrRNA.fna');
  expect(away).toEqual([SIXTEEN_S.url]);

  await page.getByTestId('example-refseq-16s-bacteria-delete').click();
  await expect(note).toHaveText('Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov');
  await expect(page.getByTestId('subject-source-0')).toHaveCount(0);
  await expect(page.getByTestId('example-refseq-16s-bacteria-delete')).toHaveCount(0);
});

test('a failed download names the site and the status and leaves the boxes as they were', async ({ page }) => {
  const away = await stubHosts(page);
  await openSearch(page);
  await program(page, 'blastp');
  await expect(page.getByTestId('example-uniprot-ecoli-k12-note')).toHaveText('Downloads 981.2 kB from rest.uniprot.org');
  await expect(page.getByTestId('example-pdbaa-note')).toHaveText('Downloads 36.8 MB from ftp.ncbi.nlm.nih.gov');
  await page.getByTestId('query-input').fill('>old\nMKV\n');
  await page.getByTestId('example-uniprot-ecoli-k12').click();
  await expect(page.getByTestId('example-error')).toHaveText('rest.uniprot.org answered HTTP 500: try again in a moment.');
  await expect(page.getByTestId('query-input')).toHaveValue('>old\nMKV\n');
  await expect(page.getByTestId('subject-source-0')).toHaveCount(0);
  expect(away).toEqual([exampleById('uniprot-ecoli-k12').url]);
});

test("a damaged gzip body names the site and the engine's sentence, and leaves the boxes as they were", async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine's decoders inflate the set (engine build)");
  const gz = gzipSync(SET);
  const away = await stubHosts(page, gz.subarray(0, gz.length >> 1));
  await openSearch(page);
  await page.getByTestId('query-input').fill('>old\nACGT\n');
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await expect(page.getByTestId('example-error')).toHaveText(
    "The file from ftp.ncbi.nlm.nih.gov could not be decompressed ('bacteria.16SrRNA.fna.gz' is named .gz but is not gzip data (incomplete deflate stream); LOSAT reads a file named .gz through gzip).",
  );
  await expect(page.getByTestId('query-input')).toHaveValue('>old\nACGT\n');
  await expect(page.getByTestId('subject-source-0')).toHaveCount(0);
  expect(away).toEqual([SIXTEEN_S.url]);
});

test('Delete in one tab removes the subject source made from the copy in another tab', async ({ context }) => {
  const a = await context.newPage();
  const b = await context.newPage();
  await stubHosts(a);
  await stubHosts(b);
  await openSearch(a);
  test.skip((await a.getByTestId('storage-status').getAttribute('data-backend')) !== 'opfs', 'nothing is saved without OPFS');
  await a.getByTestId('example-refseq-16s-bacteria').click();
  await settled(a, 'subject');
  await openSearch(b);
  await expect(b.getByTestId('example-refseq-16s-bacteria-note')).toHaveText('Saved in this browser (2026-10-27, 111 bytes)');
  await b.getByTestId('example-refseq-16s-bacteria-delete').click();
  await expect(a.getByTestId('subject-source-0')).toHaveCount(0);
  await expect(a.getByTestId('example-refseq-16s-bacteria-note')).toHaveText('Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov');
});

test("the saved copy goes to a folder that the user picks, or to the browser's download, without a request", async ({ page }) => {
  const away = await stubHosts(page);
  // The folder picker cannot be driven: where the browser has one, it picks the folder "picked" of OPFS,
  // and refuses as the real one does when the click's activation has passed.
  await page.addInitScript(() => {
    const target = window as { showDirectoryPicker?: unknown };
    if (typeof target.showDirectoryPicker === 'function') {
      target.showDirectoryPicker = async () => {
        if (!navigator.userActivation.isActive) {
          throw new DOMException('Must be handling a user gesture to show a file picker.', 'SecurityError');
        }
        return (await navigator.storage.getDirectory()).getDirectoryHandle('picked', { create: true });
      };
    }
  });
  await openSearch(page);
  test.skip((await page.getByTestId('storage-status').getAttribute('data-backend')) !== 'opfs', 'nothing is saved without OPFS');
  const folders = await page.evaluate(() => typeof (window as { showDirectoryPicker?: unknown }).showDirectoryPicker === 'function');
  console.log(`the copy goes to ${folders ? 'a folder' : 'the download'}`);
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await settled(page, 'subject');
  const save = page.getByTestId('example-refseq-16s-bacteria-save');
  if (folders) {
    await expect(save).toHaveText('Save to folder…');
    await save.click();
    await expect(page.getByTestId('example-note')).toHaveText('Saved bacteria.16SrRNA.fna in picked.');
    const written = await page.evaluate(async () => {
      const folder = await (await navigator.storage.getDirectory()).getDirectoryHandle('picked');
      return (await (await folder.getFileHandle('bacteria.16SrRNA.fna')).getFile()).text();
    });
    expect(written).toBe(SET);
  } else {
    await expect(save).toHaveText('Download');
    const [download] = await Promise.all([page.waitForEvent('download'), save.click()]);
    expect(download.suggestedFilename()).toBe('bacteria.16SrRNA.fna');
    expect(readFileSync(await download.path(), 'utf8')).toBe(SET);
  }
  expect(away).toEqual([SIXTEEN_S.url]);
});
