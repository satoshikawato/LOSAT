// Line WF-1 (docs/web/fetch_design.md §4): accessions in a paste box are fetched only when the
// user presses "Fetch", from the source chosen (NCBI, EBI, CNCB), and their FASTA is written
// into the box. Line WF-3 (§10): "Fastest", the default, races a fixed public record at the
// three sources on the first Fetch and sends the accessions to the first that answers. Line WF-5
// (§11): a box that no one source serves has a Fetch per group of accessions, one source each. The
// sources are stubbed with routes at their real addresses, so the page runs under the deployed
// CSP (public/_headers connect-src); the tests also count every request that leaves the site.
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { E2E_ORIGIN } from './support/origin';
import { openSearch, paste, program, settled, submit } from './support/search';

const RECORDS: Readonly<Record<string, string>> = {
  'NR_074540.1': '>NR_074540.1 Bacillus cereus ATCC 14579 16S ribosomal RNA (rrnA), partial sequence\nTTATTGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAA\n',
  'AB000001.1': '>ENA|AB000001|AB000001.1 Rhizoctonia solani genes for 18S rRNA\nAATTTTAATGAAAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGC\n',
  P12345: '>sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial\nMALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYKRDTNSKK\n',
  'MH011443.1': '>MH011443.1 Homo sapiens TP53 (TP53) gene, exon 5 and partial cds\nTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCGTGGCCATCTAC\n',
  'C_AA004835.1': '>C_AA004835.1 GenBase protein\nMESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQ\n',
  'NP_414542.1': '>NP_414542.1 thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]\nMKRISTTITTTITITTGNGAG\n',
  P0AD86: '>sp|P0AD86|LPT_ECOLI thr operon leader peptide OS=Escherichia coli (strain K12)\nMKRISTTITTTITITTGNGAG\n',
};
const PROBES = [
  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MH011443.1&rettype=fasta&retmode=text',
  'https://www.ebi.ac.uk/ena/browser/api/fasta/MH011443.1',
  'https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=MH011443.1',
];
const SOURCES = /^https:\/\/(?:eutils\.ncbi\.nlm\.nih\.gov|www\.ebi\.ac\.uk|rest\.uniprot\.org|ngdc\.cncb\.ac\.cn)\//;

/**
 * Stubs the sources (an unknown accession: 404), each answering after its delay in milliseconds,
 * and returns the requests that leave the site, as they are made.
 */
async function stubSources(page: Page, delays: Readonly<Record<string, number>> = {}): Promise<string[]> {
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
    const delay = delays[url.host] ?? 0;
    if (delay > 0) await new Promise((resolve) => setTimeout(resolve, delay));
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
  await page.getByTestId('query-fetch-source').selectOption('ncbi');
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

test('Fastest races a fixed public record at the three sources on the first Fetch only, and sends the accession to the first that answers', async ({
  page,
}) => {
  const away = await stubSources(page, { 'eutils.ncbi.nlm.nih.gov': 1500, 'ngdc.cncb.ac.cn': 1000 });
  await openSearch(page);
  await page.getByTestId('query-input').fill('AB000001.1');
  await expect(page.getByTestId('query-fetch-source')).toHaveValue('fastest');
  await expect(page.getByTestId('query-fetch-note')).toHaveText('Sends the accession to the fastest of NCBI, EBI and CNCB');
  // EBI has no RefSeq records: it is not named for a RefSeq accession.
  await page.getByTestId('subject-input').fill('NR_074540.1');
  await expect(page.getByTestId('subject-fetch-note')).toHaveText('Sends the accession to the fastest of NCBI and CNCB');
  expect(away).toEqual([]);

  await page.getByTestId('query-fetch').click();
  await expect(page.getByTestId('query-input')).toHaveValue(RECORDS['AB000001.1']!);
  // The race asked for its own record only; the accession went to EBI, the first to answer.
  expect([...away.slice(0, 3)].sort()).toEqual([...PROBES].sort());
  expect(away.slice(3)).toEqual(['https://www.ebi.ac.uk/ena/browser/api/fasta/AB000001.1']);
  await expect(page.getByTestId('subject-fetch-source').locator('option[value="fastest"]')).toHaveText('Fastest (CNCB)');
  await expect(page.getByTestId('subject-fetch-note')).toHaveText('Sends the accession to CNCB');

  // The order is kept: no race again.
  await page.getByTestId('subject-fetch').click();
  await expect(page.getByTestId('subject-input')).toHaveValue(RECORDS['NR_074540.1']!);
  expect(away.slice(4)).toEqual(['https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=NR_074540.1']);
  await page.getByTestId('query-input').fill('AB000001.1');
  await expect(page.getByTestId('query-fetch-source').locator('option[value="fastest"]')).toHaveText('Fastest (EBI)');
  await expect(page.getByTestId('query-fetch-note')).toHaveText('Sends the accession to EBI');
  expect(away).toHaveLength(5);
});

test('a box that no one site serves has a Fetch per group, each sends its accessions to one site only, and the box is written once all have come', async ({
  page,
}) => {
  const away = await stubSources(page);
  await openSearch(page);
  await program(page, 'blastp');
  await page.getByTestId('subject-input').fill('C_AA004835.1\nP12345\nP99999');
  await expect(page.getByTestId('subject-fetch-note')).toHaveText('No one site has all of these accessions');
  await expect(page.getByTestId('subject-fetch')).toHaveCount(0);
  await expect(page.getByTestId('subject-fetch-group-0')).toHaveText(/^Fetch C_AA004835\.1\s*from\s*CNCB$/);
  await expect(page.getByRole('button', { name: 'Fetch C_AA004835.1 from CNCB' })).toBeVisible();
  await expect(page.getByTestId('subject-fetch-group-1-button')).toHaveText('Fetch P12345 and 1 more');
  await expect(page.getByRole('button', { name: 'Fetch P12345 and 1 more from NCBI' })).toHaveAttribute('title', 'P12345, P99999');
  await expect(page.getByTestId('subject-fetch-group-1-source')).toHaveValue('ncbi');
  expect(away).toEqual([]);

  await page.getByTestId('subject-fetch-group-0-button').click();
  await expect(page.getByTestId('subject-fetch-group-0-button')).toHaveText('Fetched C_AA004835.1');
  await expect(page.getByTestId('subject-fetch-group-0-button')).toBeDisabled();
  await expect(page.getByTestId('subject-input')).toHaveValue('C_AA004835.1\nP12345\nP99999');
  expect(away).toEqual(['https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=C_AA004835.1']);

  // NCBI has no P99999 here: the group fails alone, names the site and the accession, and keeps the other group.
  await page.getByTestId('subject-fetch-group-1-button').click();
  await expect(page.getByTestId('subject-fetch-error')).toHaveText('NCBI has no protein record P99999.');
  await expect(page.getByTestId('subject-fetch-group-0-button')).toHaveText('Fetched C_AA004835.1');
  // Another site for the group drops the message of the failed one.
  await page.getByTestId('subject-fetch-group-1-source').selectOption('ebi');
  await expect(page.getByRole('button', { name: 'Fetch P12345 and 1 more from EBI' })).toBeVisible();
  await expect(page.getByTestId('subject-fetch-error')).toHaveCount(0);
  await page.getByTestId('subject-input').fill('C_AA004835.1\nP12345');
  // An edit of the box drops what came.
  await expect(page.getByTestId('subject-fetch-group-0-button')).toHaveText('Fetch C_AA004835.1');
  await expect(page.getByTestId('subject-fetch-group-1-source')).toHaveValue('ncbi');
  await page.getByTestId('subject-fetch-group-1-source').selectOption('ebi');
  await page.getByTestId('subject-fetch-group-1-button').click();
  await expect(page.getByTestId('subject-fetch-group-1-button')).toHaveText('Fetched P12345');
  await page.getByTestId('subject-fetch-group-0-button').click();
  await expect(page.getByTestId('subject-input')).toHaveValue(RECORDS['C_AA004835.1']! + RECORDS['P12345']!);
  await settled(page, 'subject');
  await expect(page.getByTestId('subject-fetch-row')).toHaveCount(0);
  // One site a press, only the group's accessions, and no race.
  expect(away).toEqual([
    'https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=C_AA004835.1',
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=P12345&rettype=fasta&retmode=text',
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=P99999&rettype=fasta&retmode=text',
    'https://rest.uniprot.org/uniprotkb/P12345.fasta',
    'https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=C_AA004835.1',
  ]);
});

test('a site chosen by hand starts each group it serves; a box that one site serves keeps the one Fetch', async ({ page }) => {
  const away = await stubSources(page);
  await openSearch(page);
  await program(page, 'blastp');
  await page.getByTestId('subject-input').fill('P12345');
  await page.getByTestId('subject-fetch-source').selectOption('ebi');
  await page.getByTestId('subject-input').fill('C_AA004835.1\nP12345');
  await expect(page.getByTestId('subject-fetch-note')).toHaveText('No one site has all of these accessions');
  await expect(page.getByTestId('subject-fetch-group-1-source')).toHaveValue('ebi');
  await expect(page.getByTestId('subject-fetch-group-0')).toHaveText(/^Fetch C_AA004835\.1\s*from\s*CNCB$/);
  await page.getByTestId('subject-input').fill('C_AA004835.1');
  await expect(page.getByTestId('subject-fetch-groups')).toHaveCount(0);
  await expect(page.getByTestId('subject-fetch-source')).toHaveValue('ebi');
  await page.getByTestId('subject-fetch').click();
  await expect(page.getByTestId('subject-fetch-error')).toHaveText(
    'EBI (UniProt) has UniProt accessions only: fetch C_AA004835.1 from CNCB.',
  );
  expect(away).toEqual([]);
});

test('a group starts at the race order of the tab, and keeps the site its records came from', async ({ page }) => {
  const away = await stubSources(page, { 'eutils.ncbi.nlm.nih.gov': 1500, 'ngdc.cncb.ac.cn': 1000 });
  await openSearch(page);
  await program(page, 'blastp');
  await page.getByTestId('subject-input').fill('C_AA004835.1\nP12345');
  await expect(page.getByTestId('subject-fetch-group-1-source')).toHaveValue('ncbi');
  await page.getByTestId('subject-fetch-group-1-button').click();
  await expect(page.getByTestId('subject-fetch-status')).toHaveText('Fetching 1 of 1 from NCBI…');
  await expect(page.getByRole('button', { name: 'Fetched P12345 from NCBI' })).toBeDisabled();

  // The query's Fetch races the protein sites; EBI answers first.
  await page.getByTestId('query-input').fill('P12345');
  await page.getByTestId('query-fetch').click();
  await expect(page.getByTestId('query-input')).toHaveValue(RECORDS['P12345']!);
  // What came stays named by where it came from; a new box starts at the race's first.
  await expect(page.getByRole('button', { name: 'Fetched P12345 from NCBI' })).toBeDisabled();
  await expect(page.getByTestId('subject-fetch-group-1-source')).toHaveValue('ncbi');
  await page.getByTestId('subject-input').fill('C_AA004835.1\nP12345\n');
  await expect(page.getByTestId('subject-fetch-group-1-source')).toHaveValue('ebi');
  expect(away.slice(0, 1)).toEqual([
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=P12345&rettype=fasta&retmode=text',
  ]);
  // The race's fixed records, then the query's accession at EBI; nothing more of the subject.
  expect([...away.slice(1)].sort()).toEqual(
    [
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_414542.1&rettype=fasta&retmode=text',
      'https://rest.uniprot.org/uniprotkb/P0AD86.fasta',
      'https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=NP_414542.1',
      'https://rest.uniprot.org/uniprotkb/P12345.fasta',
    ].sort(),
  );
});

test('accessions on one line are put one per line by one press, for Fetch (WA2-OP-04)', async ({ page }) => {
  const away = await stubSources(page);
  await openSearch(page);
  // Accessions on one line, separated by white space or commas, cannot be read (NCBI's reader
  // error): one press puts each on its own line, where Fetch takes them.
  for (const text of ['AB000001.1 MH011443.1', 'AB000001.1,MH011443.1', 'ab000001.1, MH011443.1\n']) {
    await paste(page, 'query', text);
    await expect(page.getByTestId('query-source-0-error')).toContainText('This input cannot be read: ');
    await expect(page.getByTestId('query-source-0-split-accessions')).toHaveText('Put each accession on its own line');
    await expect(page.getByTestId('query-source-0-add-defline')).toHaveCount(0);
    await expect(page.getByTestId('query-fetch-row')).toHaveCount(0);
  }
  await page.getByTestId('query-source-0-split-accessions').click();
  await expect(page.getByTestId('query-input')).toHaveValue('ab000001.1\nMH011443.1');
  await expect(page.getByTestId('query-fetch-row')).toBeVisible();
  await expect(page.getByTestId('query-source-0')).toHaveCount(0);
  // The card goes; the focus goes to the Fetch, not to the page's body.
  await expect(page.getByTestId('query-fetch')).toBeFocused();
  // A box that no one site serves: the focus goes to the first group's Fetch.
  await program(page, 'blastp');
  await paste(page, 'subject', 'C_AA004835.1, P12345');
  await page.getByTestId('subject-source-0-split-accessions').press('Enter');
  await expect(page.getByTestId('subject-input')).toHaveValue('C_AA004835.1\nP12345');
  await expect(page.getByTestId('subject-fetch-group-0-button')).toBeFocused();
  await program(page, 'blastn');
  // A box with anything but accessions keeps the defline's button only.
  await paste(page, 'query', 'AB000001.1 MH011443.1 ACGT-1');
  await expect(page.getByTestId('query-source-0-split-accessions')).toHaveCount(0);
  await expect(page.getByTestId('query-source-0-add-defline')).toBeVisible();
  expect(away).toEqual([]);
});

test('a first line that is an accession and NCBI BLAST+ would fetch as a Seq-id points to Fetch; a GI, a PDB id and a tagged id do not (WA2-OP-04)', async ({
  page,
}) => {
  // The FakeEngine's reader never takes a first line for a Seq-id (infra/fake/fake-fasta.ts).
  test.skip(!BUILD_HAS_ENGINE, 'the reason is the adapter scan of the engine build');
  const away = await stubSources(page);
  await openSearch(page);
  // The adapter's reason is NCBI's text and the defline; the card adds the Fetch sentence for an accession.
  const reason = (line: string) =>
    `This input cannot be read: the first line ("${line}") is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported in the browser (start the input with a '>' defline)`;
  const fetchSentence = 'To fetch accessions, put only the accessions in the box, one per line.';
  for (const [role, text, line] of [
    ['query', 'AB000001.1\n>x\nACGTACGT', 'AB000001.1'],
    ['query', 'AB000001.1\nACGTACGT', 'AB000001.1'],
    ['subject', ' NR_074540.1 \nACGTACGT', 'NR_074540.1'],
  ] as const) {
    await paste(page, role, text);
    await expect(page.getByTestId(`${role}-source-0-error`)).toHaveText(reason(line));
    await expect(page.getByTestId(`${role}-source-0-fetch-hint`)).toHaveText(fetchSentence);
    await expect(page.getByTestId(`${role}-source-0-add-defline`)).toBeVisible();
    await expect(page.getByTestId(`${role}-source-0-split-accessions`)).toHaveCount(0);
  }
  // A GI number, a PDB id and a tagged id are Seq-ids to the engine, but not what Fetch takes: no sentence.
  for (const [role, program_, text, line] of [
    ['query', 'blastn', '12345678\nACGTACGT', '12345678'],
    ['query', 'blastp', '4HHB\nMKLV', '4HHB'],
    ['subject', 'blastp', '1ABC_A\nMKLV', '1ABC_A'],
    ['query', 'blastp', 'sp|P04637|P53_HUMAN\nMKLV', 'sp|P04637|P53_HUMAN'],
    ['query', 'blastn', 'ref|NM_000546.6|\nACGTACGT', 'ref|NM_000546.6|'],
  ] as const) {
    await program(page, program_);
    await paste(page, role, text);
    await expect(page.getByTestId(`${role}-source-0-error`)).toHaveText(reason(line));
    await expect(page.getByTestId(`${role}-source-0-fetch-hint`)).toHaveCount(0);
  }
  await program(page, 'blastn');
  // Accessions only, more than one on a line: the reason, the sentence and the press that splits them.
  await paste(page, 'subject', 'NR_074540.1\nAB000001.1 MH011443.1');
  await expect(page.getByTestId('subject-source-0-error')).toHaveText(reason('NR_074540.1'));
  await expect(page.getByTestId('subject-source-0-fetch-hint')).toHaveText(fetchSentence);
  await expect(page.getByTestId('subject-source-0-split-accessions')).toBeVisible();
  expect(away).toEqual([]);
});
