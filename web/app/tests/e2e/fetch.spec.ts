// Line WF-1 (docs/web/fetch_design.md §4): accessions in a paste box are fetched only when the
// user presses "Fetch", from the source chosen (NCBI, EBI, CNCB), and their FASTA is written
// into the box. The sources are stubbed with routes at their real addresses, so the page runs
// under the deployed CSP (public/_headers connect-src); the tests also count every request
// that leaves the site.
import { expect, test, type Page } from '@playwright/test';
import { E2E_ORIGIN } from './support/origin';
import { openSearch, paste, program, settled, submit } from './support/search';

const RECORDS: Readonly<Record<string, string>> = {
  'NR_074540.1': '>NR_074540.1 Bacillus cereus ATCC 14579 16S ribosomal RNA (rrnA), partial sequence\nTTATTGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAA\n',
  'AB000001.1': '>ENA|AB000001|AB000001.1 Rhizoctonia solani genes for 18S rRNA\nAATTTTAATGAAAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGC\n',
  P12345: '>sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial\nMALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYKRDTNSKK\n',
};
const SOURCES = /^https:\/\/(?:eutils\.ncbi\.nlm\.nih\.gov|www\.ebi\.ac\.uk|rest\.uniprot\.org|ngdc\.cncb\.ac\.cn)\//;

/** Stubs the sources (an unknown accession: 404), and returns the requests that leave the site, as they are made. */
async function stubSources(page: Page): Promise<string[]> {
  const away: string[] = [];
  page.on('request', (request) => {
    const url = request.url();
    if (/^https?:/.test(url) && !url.startsWith(E2E_ORIGIN)) away.push(url);
  });
  await page.route(SOURCES, async (route) => {
    const url = new URL(route.request().url());
    const accession =
      url.searchParams.get('id') ?? url.searchParams.get('acc') ?? url.pathname.split('/').pop()!.replace(/\.fasta$/, '');
    const body = RECORDS[accession];
    await route.fulfill({
      status: body === undefined ? 404 : 200,
      contentType: 'text/plain',
      headers: { 'Access-Control-Allow-Origin': '*' },
      body: body ?? 'not found',
    });
  });
  return away;
}

async function fetchInto(page: Page, role: 'query' | 'subject', accessions: string): Promise<void> {
  await page.getByTestId(`${role}-input`).fill(accessions);
  await page.getByTestId(`${role}-fetch`).click();
}

test('an accession is fetched from NCBI only when asked, and its FASTA is read as a pasted text', async ({ page }) => {
  const away = await stubSources(page);
  await openSearch(page);
  await page.getByTestId('query-input').fill('NR_074540.1');
  await expect(page.getByTestId('query-fetch-note')).toHaveText('Sends the accession to NCBI');
  await expect(page.getByTestId('query-source-0')).toHaveCount(0);
  await paste(page, 'subject', '>s\nGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCG\n');
  // Run does not fetch: it asks for the fetch.
  await submit(page);
  await expect(page.getByTestId('search-message')).toHaveText('Fetch the query accession first.');
  expect(away).toEqual([]);

  await page.getByTestId('query-fetch').click();
  await expect(page.getByTestId('query-input')).toHaveValue(RECORDS['NR_074540.1']!);
  await settled(page, 'query');
  await expect(page.getByTestId('query-source-0-summary')).toHaveText('1 record · 70 nt');
  await expect(page.getByTestId('query-fetch-row')).toHaveCount(0);
  expect(away).toEqual([
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_074540.1&rettype=fasta&retmode=text',
  ]);
});

test('the source chosen stays: EBI answers from ENA and UniProt, CNCB from GenBase; a missing record is reported', async ({
  page,
}) => {
  const away = await stubSources(page);
  await openSearch(page);
  await page.getByTestId('query-input').fill('AB000001.1');
  await page.getByTestId('query-fetch-source').selectOption('ebi');
  await expect(page.getByTestId('query-fetch-note')).toHaveText('Sends the accession to EBI');
  await page.getByTestId('query-fetch').click();
  await expect(page.getByTestId('query-input')).toHaveValue(RECORDS['AB000001.1']!);

  await program(page, 'blastp');
  await fetchInto(page, 'subject', 'P12345\nP99999');
  await expect(page.getByTestId('subject-fetch-error')).toHaveText('EBI has no protein record P99999.');
  await expect(page.getByTestId('subject-input')).toHaveValue('P12345\nP99999');
  await fetchInto(page, 'subject', 'P12345');
  await expect(page.getByTestId('subject-input')).toHaveValue(RECORDS['P12345']!);

  await program(page, 'blastn');
  await page.getByTestId('query-input').fill('NR_074540.1');
  await page.getByTestId('query-fetch-source').selectOption('cncb');
  await page.getByTestId('query-fetch').click();
  await expect(page.getByTestId('query-input')).toHaveValue(RECORDS['NR_074540.1']!);
  expect(away).toEqual([
    'https://www.ebi.ac.uk/ena/browser/api/fasta/AB000001.1',
    'https://rest.uniprot.org/uniprotkb/P12345.fasta',
    'https://rest.uniprot.org/uniprotkb/P99999.fasta',
    'https://rest.uniprot.org/uniprotkb/P12345.fasta',
    'https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=NR_074540.1',
  ]);
});
