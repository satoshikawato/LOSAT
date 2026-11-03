// Line WF-2 (docs/web/fetch_design.md §9): "Example" downloads a mode's subject set only when the
// user presses it, puts the set in the subject list and the mode's query in the query box, and
// keeps the set in this browser (OPFS examples/), so that the next time it loads without a
// request. The download hosts are stubbed with routes at their real addresses, so the page runs
// under the deployed CSP (public/_headers connect-src); the tests count every request that
// leaves the site.
import { gzipSync } from 'node:zlib';
import { expect, test, type Page } from '@playwright/test';
import { exampleById, exampleQuery } from '../../src/domain/examples';
import { E2E_ORIGIN } from './support/origin';
import { openSearch, program, settled } from './support/search';

const SIXTEEN_S = exampleById('refseq-16s-bacteria');
const SET = '>NR_102783.2 Bacillus subtilis 16S\nACGTACGTACGTACGTACGT\n>NR_114042.1 Escherichia coli 16S\nACGTTCGTACGTACGTTCGT\n';
const HOSTS = /^https:\/\/(?:ftp\.ncbi\.nlm\.nih\.gov|rest\.uniprot\.org)\//;

/** Stubs the download hosts (the 16S set as gzip, UniProt with a 500), and returns the requests that leave the site. */
async function stubHosts(page: Page): Promise<string[]> {
  const away: string[] = [];
  page.on('request', (request) => {
    const url = request.url();
    if (/^https?:/.test(url) && !url.startsWith(E2E_ORIGIN)) away.push(url);
  });
  const gz = gzipSync(SET);
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

  await page.reload();
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
