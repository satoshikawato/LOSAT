// LOSAT Web's own files of a run (S15 item 2; design §12.1, §12.3) through the real application: the
// Outputs tab's "LOSAT Web formats" section, apart from the compatibility outputs; CSV, JSON and the
// report of each scope (the whole run, after the view filters, the subjects marked in the
// Descriptions), read back as the browser saved them; the statement that each is a LOSAT Web
// format; the counts; IDs and a title made to break HTML, CSV and JSON, written as text; and the
// report opened with the network blocked: it renders and makes no request.
//
// The FakeEngine build knows the HSPs (src/infra/fake/fake-engine.ts): three queries on three
// subjects give 13 HSPs with E values of index * 1e-5. The engine build runs a small BLASTN search
// and checks the files against its stored outputs. The outfmt 6 fields of the files are compared
// with the stored outfmt 6 text as the Outputs view shows it, never with values computed here.
import { readFile } from 'node:fs/promises';
import { expect, test, type BrowserContext, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { paste, program, showOutput, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await page.goto('/');
});

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

/** IDs (a FASTA ID ends at the first space) and a title that must stay text in every file. */
const HOSTILE = {
  query: '</style><b>q1</b>',
  subjects: ['"><img/src=x/onerror=alert(1)>', '<script>alert(2)</script>,"q"', 'javascript:alert(3)'],
  title: '<i>Exports</i> & "quotes" \'here\'',
};
const QUERIES = `>${HOSTILE.query} first query\n${dna(1, 40)}\n>q2\n${dna(2, 30)}\n>q3\n${dna(3, 24)}\n`;
const SUBJECTS = HOSTILE.subjects.map((id, i) => `>${id} subject ${i + 1}\n${dna(10 + i, 50)}\n`).join('');
/** The FakeEngine's HSPs of these inputs: q1 4, q2 5, q3 4 (index 0-12, E value index * 1e-5). */
const HSPS = 13;
const CSP = "default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'";

/** The records of a CSV text (RFC 4180, CRLF), the header's first: a reader independent of the writer. */
function parseCsv(text: string): string[][] {
  const records: string[][] = [];
  let fields: string[] = [];
  let i = 0;
  while (i < text.length) {
    let field = '';
    if (text[i] === '"') {
      i++;
      for (;;) {
        const quote = text.indexOf('"', i);
        field += text.slice(i, quote);
        i = quote + 1;
        if (text[i] !== '"') break;
        field += '"';
        i++;
      }
    } else {
      while (i < text.length && text[i] !== ',' && text[i] !== '\r') field += text[i++];
    }
    fields.push(field);
    if (text[i] === ',') i++;
    else {
      expect(text.slice(i, i + 2)).toBe('\r\n');
      i += 2;
      records.push(fields);
      fields = [];
    }
  }
  return records;
}

/** The stored output of one format of run 1 as the Outputs view shows it. */
async function storedText(page: Page, format: 0 | 6): Promise<string> {
  await showOutput(page, 1, format, false);
  return (await page.getByTestId('result-output').textContent()) ?? '';
}

/** The text of the file that clicking a button saves (the browser's download). */
async function saved(page: Page, testid: string): Promise<{ name: string; text: string }> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await download;
  return { name: file.suggestedFilename(), text: await readFile((await file.path())!, 'utf8') };
}

const scope = (page: Page, name: 'all' | 'filtered' | 'marked') => page.getByTestId(`export-scope-${name}`);

/** Saves one format of the chosen scope and waits for the summary line. */
async function exportFile(page: Page, format: 'csv' | 'json' | 'report', hsps: number): Promise<{ name: string; text: string }> {
  const file = await saved(page, `export-${format}`);
  const summary = page.getByTestId('export-summary');
  await expect(summary).toHaveAttribute('data-format', format);
  await expect(summary).toHaveAttribute('data-hsps', String(hsps));
  await expect(summary).toContainText(`Saved ${file.name}: ${hsps} ${hsps === 1 ? 'HSP' : 'HSPs'}`);
  return file;
}

/** The stored outfmt 6 rows of run 1 as the Outputs view shows them, split into their fields. */
async function storedOutfmt6(page: Page): Promise<string[][]> {
  return (await storedText(page, 6))
    .split('\n')
    .filter((line) => line !== '' && !line.startsWith('#'))
    .map((line) => line.split('\t'));
}

/**
 * Opens a saved report in a page of its own with every request refused and counted: it must render
 * from its own text, run nothing and ask for nothing.
 */
async function openReport(context: BrowserContext, html: string) {
  const report = await context.newPage();
  const requests: string[] = [];
  await report.route('**/*', (route) => {
    requests.push(route.request().url());
    return route.abort();
  });
  const dialogs: string[] = [];
  report.on('dialog', (dialog) => {
    dialogs.push(dialog.message());
    void dialog.dismiss();
  });
  await report.setContent(html, { waitUntil: 'load' });
  await expect(report.locator('h1')).toHaveText('LOSAT Web report');
  const page = await report.evaluate(() => ({
    active: document.querySelectorAll('script, img, iframe, object, embed, link, a, form, svg, base, input, button').length,
    csp: document.querySelector('meta[http-equiv="Content-Security-Policy"]')?.getAttribute('content'),
    text: document.body.textContent ?? '',
    styles: document.querySelectorAll('style').length,
  }));
  return { report, requests, dialogs, ...page };
}

test.describe('FakeEngine build', () => {
  test.skip(BUILD_HAS_ENGINE, "the expected files follow the FakeEngine's HSPs");

  test('CSV, JSON and the report of each scope say what they are and keep hostile IDs as text', async ({ page, context }) => {
    await program(page, 'blastn');
    await paste(page, 'query', QUERIES);
    await paste(page, 'subject', SUBJECTS);
    await page.getByTestId('job-title').fill(HOSTILE.title);
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    const outfmt6 = await storedOutfmt6(page);
    expect(outfmt6).toHaveLength(HSPS);

    // The Outputs tab, which "Download All" opens: the compatibility outputs, then LOSAT Web's formats apart.
    await page.getByTestId('results-download-all').click();
    await expect(page.getByTestId('results-outputs')).toContainText('The outputs as the engine wrote them. The view filters do not change them.');
    await expect(page.getByTestId('export-formats-note')).toContainText('application formats of LOSAT Web, not NCBI BLAST formats');
    await expect(scope(page, 'all')).toHaveAttribute('data-count', String(HSPS));
    await expect(scope(page, 'all')).toBeChecked();
    await expect(scope(page, 'filtered')).toHaveAttribute('data-count', String(HSPS));
    await expect(scope(page, 'marked')).toHaveAttribute('data-count', '0');
    await expect(scope(page, 'marked')).toBeDisabled();
    await expect(page.getByTestId('export-json-aligned')).toBeChecked();

    // --- the whole run ---
    const csv = await exportFile(page, 'csv', HSPS);
    expect(csv.name).toBe('losat-run1-blastn-hsps.csv');
    const records = parseCsv(csv.text);
    expect(records[0]!.slice(0, 8)).toEqual(['run', 'query_record', 'query_id', 'subject_record', 'subject_id', 'hsp', 'index', 'rank']);
    expect(records.slice(1).map((fields) => fields.slice(8, 20))).toEqual(outfmt6);
    expect(records[1]!.slice(0, 8)).toEqual(['1', '1', HOSTILE.query, '1', HOSTILE.subjects[0], '1.1', '0', '0']);
    expect(records[2]!.slice(3, 6)).toEqual(['2', HOSTILE.subjects[1], '1.2']);
    // The ID with a comma and quotes is quoted, its quotes doubled; nothing is altered.
    expect(csv.text).toContain('"<script>alert(2)</script>,""q"""');

    const json = await exportFile(page, 'json', HSPS);
    expect(json.name).toBe('losat-run1-blastn-hsps.json');
    const doc = JSON.parse(json.text);
    expect(doc.format).toBe('LOSAT Web HSP export');
    expect(doc.note).toMatch(/not an NCBI BLAST output/);
    expect(doc.run).toMatchObject({ number: 1, title: HOSTILE.title, program: 'blastn', engine_build: 'fake-engine' });
    expect(doc.run.query).toMatchObject({ name: 'query.fa', records: 3 });
    expect(doc.scope).toMatchObject({ name: 'all', hsps: HSPS, filters: null, aligned_sequences: true });
    expect(doc.hsps.map((hsp: { outfmt6: Record<string, string> }) => Object.values(hsp.outfmt6))).toEqual(outfmt6);
    expect(doc.hsps.map((hsp: { record: { index: number } }) => hsp.record.index)).toEqual([...Array(HSPS).keys()]);
    expect(doc.hsps[0]).toMatchObject({ query_id: HOSTILE.query, subject_id: HOSTILE.subjects[0], hsp: '1.1' });
    expect(doc.hsps[0].record.query_aligned).toMatch(/^FAKE/);

    const all = await exportFile(page, 'report', HSPS);
    expect(all.name).toBe('losat-run1-blastn-report.html');
    const opened = await openReport(context, all.text);
    expect(opened.requests).toEqual([]);
    expect(opened.dialogs).toEqual([]);
    expect(opened.active).toBe(0);
    expect(opened.styles).toBe(1);
    expect(opened.csp).toBe(CSP);
    expect(opened.text).toContain('An application format of LOSAT Web.');
    expect(opened.text).toContain('It is not an NCBI BLAST report');
    for (const text of [HOSTILE.query, ...HOSTILE.subjects, HOSTILE.title]) expect(opened.text).toContain(text);
    expect(opened.text).toContain('LOSAT blastn -query query.fa -subject subject.fa -outfmt 0');
    expect(opened.text).toContain('Development build');
    expect(opened.text).toContain('Warning: FAKE ENGINE OUTPUT');
    await expect(opened.report.locator('section.query')).toHaveCount(3);
    await expect(opened.report.locator('pre.section')).toHaveCount(9);
    expect(all.text).toContain('<dt>Alignments</dt><dd>Included: the outfmt 0 headings and sections of the HSPs that outfmt 0 shows, as written.</dd>');
    await opened.report.close();

    // --- after the view filters: E values of at most 5e-5 keep HSPs 0-5 ---
    await page.getByTestId('filter-evalue').fill('0.00005');
    await page.getByTestId('filter-apply').click();
    await expect(scope(page, 'filtered')).toHaveAttribute('data-count', '6');
    await expect(scope(page, 'all')).toHaveAttribute('data-count', String(HSPS));
    await scope(page, 'filtered').check();
    const filteredCsv = parseCsv((await exportFile(page, 'csv', 6)).text);
    expect(filteredCsv.slice(1).map((fields) => fields[5])).toEqual(['1.1', '1.2', '1.3', '1.4', '2.1', '2.2']);
    const filteredJson = JSON.parse((await exportFile(page, 'json', 6)).text);
    expect(filteredJson.scope).toMatchObject({ name: 'filtered', hsps: 6, filters: { max_e_value: 0.00005 } });
    const filteredReport = await exportFile(page, 'report', 6);
    expect(filteredReport.text).toContain('<dt>View filters</dt><dd>E value at most 0.00005</dd>');
    expect(filteredReport.text.match(/<section class="query">/g)).toHaveLength(2);

    // --- the subjects marked in the Descriptions: q1's second subject (HSPs 1.2 and 1.3) ---
    await page.getByTestId('results-view-hits').click();
    await page.getByTestId('subject-mark-1').check();
    await page.getByTestId('results-view-outputs').click();
    await expect(scope(page, 'marked')).toHaveAttribute('data-count', '2');
    await expect(scope(page, 'marked')).toBeEnabled();
    await scope(page, 'marked').check();
    expect(parseCsv((await exportFile(page, 'csv', 2)).text).slice(1).map((fields) => fields[5])).toEqual(['1.2', '1.3']);
    await page.getByTestId('export-json-aligned').uncheck();
    const markedJson = JSON.parse((await exportFile(page, 'json', 2)).text);
    expect(markedJson.scope).toMatchObject({
      name: 'marked',
      hsps: 2,
      aligned_sequences: false,
      query: { record: 1, id: HOSTILE.query },
      marked_subjects: [{ record: 2, id: HOSTILE.subjects[1] }],
    });
    expect(markedJson.hsps[0].record).not.toHaveProperty('query_aligned');
    const marked = await exportFile(page, 'report', 2);
    expect(marked.name).toBe('losat-run1-blastn-report.html');
    const markedReport = await openReport(context, marked.text);
    expect(markedReport.requests).toEqual([]);
    expect(markedReport.dialogs).toEqual([]);
    expect(markedReport.active).toBe(0);
    expect(markedReport.text).toContain(`Marked subjects2: ${HOSTILE.subjects[1]}`);
    await expect(markedReport.report.locator('section.query')).toHaveCount(1);
    await markedReport.report.close();
    // The report without the alignments (fix round 2): the HSP tables only, and the page says so.
    await expect(page.getByTestId('export-report-alignments')).toBeChecked();
    await page.getByTestId('export-report-alignments').uncheck();
    const tablesOnly = await exportFile(page, 'report', 2);
    await expect(page.getByTestId('export-summary')).toContainText('(Marked subjects, without the alignments)');
    expect(tablesOnly.text).toContain('<dt>Alignments</dt><dd>Not included: this report was saved without the outfmt 0 text of the alignments.');
    expect(tablesOnly.text).not.toContain('<pre class="section">');
    expect(tablesOnly.text).toContain('<td>1.2</td>');
    expect(tablesOnly.text.length).toBeLessThan(marked.text.length);

    // A filter that hides every HSP empties the filtered scope (and the marked one), and the choice goes back to the whole run.
    await scope(page, 'filtered').check();
    await page.getByTestId('filter-bits').fill('1000');
    await page.getByTestId('filter-apply').click();
    await expect(scope(page, 'filtered')).toHaveAttribute('data-count', '0');
    await expect(scope(page, 'filtered')).toBeDisabled();
    await expect(scope(page, 'marked')).toBeDisabled();
    await expect(scope(page, 'all')).toBeChecked();
    await expect(page.getByTestId('export-csv')).toBeEnabled();
  });
});

test.describe('engine build', () => {
  test.skip(!BUILD_HAS_ENGINE, 'the engine build searches for real');

  test('the files of a BLASTN search hold its stored outputs as written', async ({ page, context }) => {
    const a = dna(11, 60);
    await program(page, 'blastn');
    await paste(page, 'query', `>rec1 query\n${a}\n`);
    await paste(page, 'subject', `>s1 first\n${a}${dna(12, 60)}${a}${dna(13, 60)}\n>s2 second\n${dna(14, 60)}${a}${dna(15, 60)}\n`);
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    const outfmt0 = await storedText(page, 0);
    const outfmt6 = await storedOutfmt6(page);
    expect(outfmt6.length).toBeGreaterThanOrEqual(3);
    await expect(scope(page, 'all')).toHaveAttribute('data-count', String(outfmt6.length));

    const csv = parseCsv((await exportFile(page, 'csv', outfmt6.length)).text);
    expect(csv.slice(1).map((fields) => fields.slice(8, 20))).toEqual(outfmt6);
    const doc = JSON.parse((await exportFile(page, 'json', outfmt6.length)).text);
    expect(doc.hsps.map((hsp: { outfmt6: Record<string, string> }) => Object.values(hsp.outfmt6))).toEqual(outfmt6);
    expect(doc.run.engine_build).not.toBe('');
    for (const hsp of doc.hsps) {
      expect(typeof hsp.record.bit_score).toBe('number');
      expect(hsp.record.query_aligned.length).toBe(hsp.record.subject_aligned.length);
    }

    const report = await openReport(context, (await exportFile(page, 'report', outfmt6.length)).text);
    expect(report.requests).toEqual([]);
    expect(report.active).toBe(0);
    const pres = await report.report.locator('pre.section, pre.heading').allTextContents();
    expect(pres.length).toBeGreaterThan(outfmt6.length);
    for (const text of pres) expect(outfmt0).toContain(text);
    await report.report.close();

    await page.getByTestId('export-report-alignments').uncheck();
    const tables = await openReport(context, (await exportFile(page, 'report', outfmt6.length)).text);
    expect(tables.requests).toEqual([]);
    await expect(tables.report.locator('pre.section, pre.heading')).toHaveCount(0);
    await expect(tables.report.locator('section.query tbody tr')).toHaveCount(outfmt6.length);
    expect(tables.text).toContain('Not included: this report was saved without the outfmt 0 text of the alignments.');
    await tables.report.close();
  });
});
