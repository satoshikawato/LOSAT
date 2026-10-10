// The candidate tray and extraction (S14, W5; plan §5.7, design §11.3-§11.4,
// docs/web/ncbi_ui_mapping.md "候補（S14）") through the real application: HSPs added from the
// Descriptions' marks, the Alignments and the dot plot's popup, of two runs; notes, orders, moves,
// removal and "Show in results"; and the two downloads, read back as the browser saved them.
//
// The FakeEngine build checks the files byte for byte: the FakeEngine's HSPs lie at known
// coordinates of the inputs below, so the expected FASTA is written here from the input letters.
// The engine build checks real searches (a BLASTN hit on the minus strand of one of two records
// with the same ID, a TBLASTN run): the saved residues are the input file's letters at the
// intervals that outfmt 6 reports, read from the File by the Data worker.
import { readFile } from 'node:fs/promises';
import { expect, test, type Locator, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, fasta, openFiles, paste, program, showOutput, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await page.goto('/');
});

// --- inputs and files ----------------------------------------------------------------------------

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

/** Letters `from` to `to` (1-based, both included) of a sequence. */
const slice = (letters: string, from: number, to: number) => letters.slice(from - 1, to);

/** A FASTA record with its letters in lines of `width`. */
const record = (title: string, letters: string, width = 30) => `>${title}\n${letters.match(new RegExp(`.{1,${width}}`, 'g'))!.join('\n')}\n`;

/** Subject letters with a lower-case stretch: the extraction keeps the case of the file. */
const S = {
  s1: dna(31, 10) + dna(32, 10).toLowerCase() + dna(33, 30),
  s2: dna(34, 37),
  s3: dna(35, 50) + dna(36, 30).toLowerCase(),
};
const Q = { q1: dna(41, 40), q2: dna(42, 30), q3: dna(43, 24) };
const SUBJECTS = record('s1 first subject', S.s1) + record('s2', S.s2) + record('s3 third subject', S.s3);
const QUERIES = record('q1 first query', Q.q1) + record('q2', Q.q2) + record('q3', Q.q3);
const PROTEIN = 'MKVLAAGIVGLLLAHHKKEEDDPPWWRRSS';

interface FastaRecord {
  readonly header: string;
  readonly letters: string;
  readonly lines: readonly string[];
}

/** The records of a FASTA text: the header line without `>`, the letters, and the letter lines. */
function parseFasta(text: string): FastaRecord[] {
  return text
    .split('>')
    .slice(1)
    .map((block) => {
      const [header = '', ...lines] = block.split('\n');
      const letterLines = lines.filter((line) => line !== '');
      return { header, letters: letterLines.join(''), lines: letterLines };
    });
}

/** The records of a FASTA file of the repository, by position: ID and letters. */
function fileRecords(text: string): { id: string; letters: string }[] {
  return parseFasta(text).map(({ header, lines }) => ({ id: header.trim().split(/\s+/)[0]!, letters: lines.join('').replace(/\s/g, '') }));
}

/** The text of the file that clicking a button saves (the browser's download). */
async function saved(page: Page, testid: string): Promise<{ name: string; text: string }> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await download;
  return { name: file.suggestedFilename(), text: await readFile((await file.path())!, 'utf8') };
}

/** The FakeEngine's aligned row of an HSP (src/infra/fake/fake-engine.ts `fakeAlignedRows`). */
const fakeRow = (length: number, width: number) => 'FAKE'.repeat(Math.ceil(length / 4)).slice(0, length) + '-'.repeat(width - length);

/** The rows of an outfmt 6 text, split into their fields (the FakeEngine writes a comment line first). */
function outfmt6Rows(text: string): string[][] {
  return text
    .split('\n')
    .filter((line) => line !== '' && !line.startsWith('#'))
    .map((line) => line.split('\t'));
}

const reverseComplement = (letters: string) =>
  [...letters]
    .reverse()
    .map((letter) => ({ A: 'T', C: 'G', G: 'C', T: 'A', a: 't', c: 'g', g: 'c', t: 'a' })[letter] ?? letter)
    .join('');

// --- the application -----------------------------------------------------------------------------

async function search(page: Page, number: number, title = ''): Promise<void> {
  await page.getByTestId('job-title').fill(title);
  await submit(page);
  await waitStatus(page, number, 'completed');
}

/** Run 1: BLASTN of the three queries against the three subjects. */
async function blastnRun(page: Page, title = 'First run'): Promise<void> {
  await program(page, 'blastn');
  await paste(page, 'query', QUERIES);
  await paste(page, 'subject', SUBJECTS);
  await search(page, 1, title);
}

/** A TBLASTN run of a protein query against the same subjects. */
async function tblastnRun(page: Page, number: number, title: string): Promise<void> {
  await page.getByTestId('tab-search').click();
  await program(page, 'tblastn');
  await paste(page, 'query', record('p1', PROTEIN));
  await search(page, number, title);
}

const confirmation = (page: Page) => page.getByTestId('candidates-added');
const tabCount = (page: Page) => page.getByTestId('tab-candidates-count');
const rows = (page: Page) => page.getByTestId('candidate-list').locator('[data-testid^="candidate-"][data-key]');
const row = (page: Page, n: number) => page.getByTestId(`candidate-${n}`);
const field = (locator: Locator, name: string) => locator.locator(`[data-field="${name}"]`);

/** Shows a tab of the results screen. */
async function show(page: Page, testid: string): Promise<void> {
  await page.getByTestId(testid).click();
  await expect(page.getByTestId(testid)).toHaveAttribute('aria-pressed', 'true');
}

/** Marks every subject of the selected query in the Descriptions and adds their HSPs. */
async function addAllOfQuery(page: Page, expected: string): Promise<void> {
  await show(page, 'results-view-hits');
  await page.getByTestId('descriptions-select-all').check();
  await page.getByTestId('descriptions-add-candidates').click();
  await expect(confirmation(page)).toHaveText(expected);
}

/** The data-hsp ("qIdx:rank") and run of each row of the tray, in order. */
async function trayOrder(page: Page): Promise<string[]> {
  return rows(page).evaluateAll((elements) => elements.map((element) => `${element.getAttribute('data-run')}/${element.getAttribute('data-hsp')}`));
}

// --- FakeEngine build ----------------------------------------------------------------------------

test.describe('FakeEngine build', () => {
  test.skip(BUILD_HAS_ENGINE, 'the expected files follow the FakeEngine\'s HSPs');

  test('an empty tray says how to add candidates', async ({ page }) => {
    await expect(tabCount(page)).toHaveText('0');
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidates-empty')).toContainText('Add to candidates');
    await expect(page.getByTestId('candidate-list')).toHaveCount(0);
  });

  test('candidates are added from the Descriptions\' marks, the Alignments and the dot plot popup, of two runs; Origins name the runs', async ({
    page,
  }) => {
    await blastnRun(page);
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');

    // The Descriptions of q1: s1 (1 HSP), s2 (2 HSPs), s3 (1 HSP). The marks are check boxes beside the rows.
    await expect(page.getByTestId('descriptions-selected')).toHaveText('0 sequences selected');
    await expect(page.getByTestId('descriptions-add-candidates')).toBeDisabled();
    await page.getByTestId('subject-mark-0').check();
    await page.getByTestId('subject-mark-1').check();
    await expect(page.getByTestId('descriptions-selected')).toHaveText('2 sequences selected');
    await expect(page.getByTestId('descriptions-select-all')).toHaveJSProperty('indeterminate', true);
    // Marking does not select the row: the selected subject is still the first.
    await expect(page.getByTestId('subject-row-0')).toHaveAttribute('aria-pressed', 'true');
    await page.getByTestId('subject-mark-0').uncheck();
    await expect(page.getByTestId('subject-row-0')).toHaveAttribute('aria-pressed', 'true');
    await page.getByTestId('subject-mark-0').check();
    await page.getByTestId('descriptions-add-candidates').click();
    await expect(confirmation(page)).toHaveText('3 HSPs added to Candidates.');
    await expect(tabCount(page)).toHaveText('3');
    // "select all" marks every listed subject; s1 alone again is already in the tray.
    await page.getByTestId('descriptions-select-all').check();
    await expect(page.getByTestId('descriptions-selected')).toHaveText('3 sequences selected');
    await page.getByTestId('descriptions-select-all').uncheck();
    await expect(page.getByTestId('descriptions-selected')).toHaveText('0 sequences selected');
    await page.getByTestId('subject-mark-0').check();
    await page.getByTestId('descriptions-add-candidates').click();
    await expect(confirmation(page)).toHaveText('This HSP is already in Candidates.');
    await expect(tabCount(page)).toHaveText('3');

    // The Alignments of s3, under the Descriptions on the one page: "Add to candidates" beside its Range, then "In candidates".
    await page.getByTestId('subject-row-2').click();
    await expect(page.getByTestId('results-alignments-heading')).toBeInViewport();
    const range3 = page.getByTestId('range-add-0-3');
    await expect(range3).toHaveText('Add to candidates');
    await expect(page.getByTestId('alignments-add-subject')).toHaveText('Add all matches to candidates');
    await range3.click();
    await expect(confirmation(page)).toHaveText('1 HSP added to Candidates.');
    await expect(tabCount(page)).toHaveText('4');
    await expect(range3).toHaveText('In candidates');
    await expect(range3).toHaveAttribute('aria-disabled', 'true');
    await expect(page.getByTestId('alignments-add-subject')).toHaveText('All matches in candidates');
    // s2, added from the Descriptions: both Ranges say so.
    await page.getByTestId('alignments-prev-subject').click();
    await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', '1');
    await expect(page.getByTestId('range-add-0-1')).toHaveText('In candidates');
    await expect(page.getByTestId('range-add-0-2')).toHaveText('In candidates');
    // q3's subject s1 has one HSP of one letter: "Add all matches to candidates".
    await page.getByTestId('query-row-2').click();
    await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', '0');
    await page.getByTestId('alignments-add-subject').click();
    await expect(confirmation(page)).toHaveText('1 HSP added to Candidates.');
    await expect(tabCount(page)).toHaveText('5');

    // The dot plot of q2 and s1 (2 HSPs): the popup adds the selected HSP, then the next one.
    await page.getByTestId('query-row-1').click();
    await show(page, 'pane-dotplot');
    const canvas = page.getByTestId('dotplot-canvas');
    await expect(canvas).toHaveAttribute('data-segments', '2');
    await canvas.focus();
    await page.keyboard.press('Enter');
    const popup = page.getByTestId('dotplot-popup');
    await expect(popup).toHaveAttribute('data-hsp', '1:0');
    await expect(popup.getByTestId('dotplot-popup-add')).toHaveText('Add to candidates');
    await popup.getByTestId('dotplot-popup-add').click();
    await expect(confirmation(page)).toHaveText('1 HSP added to Candidates.');
    await expect(popup.getByTestId('dotplot-popup-add')).toHaveText('In candidates');
    await expect(popup.getByTestId('dotplot-popup-add')).toHaveAttribute('aria-disabled', 'true');
    await page.keyboard.press('n');
    await expect(popup).toHaveAttribute('data-hsp', '1:1');
    await expect(popup.getByTestId('dotplot-popup-add')).toHaveText('Add to candidates');
    await popup.getByTestId('dotplot-popup-add').click();
    await expect(tabCount(page)).toHaveText('7');
    await page.keyboard.press('Escape');

    // Run 2, TBLASTN: every subject of its query, from "select all".
    await tblastnRun(page, 2, 'Second run');
    await chooseRun(page, 2);
    await addAllOfQuery(page, '4 HSPs added to Candidates.');
    await expect(tabCount(page)).toHaveText('11');

    // The tray, in the order of addition.
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidate-list')).toHaveAttribute('data-count', '11');
    expect(await trayOrder(page)).toEqual(['1/0:0', '1/0:1', '1/0:2', '1/0:3', '1/2:0', '1/1:0', '1/1:1', '2/0:0', '2/0:1', '2/0:2', '2/0:3']);
    const expected = [
      { n: 1, run: 'Run 1 First run', program: 'BLASTN', query: 'q1', subject: 's1', range: '1 to 25', strand: 'plus' },
      { n: 3, run: 'Run 1 First run', program: 'BLASTN', query: 'q1', subject: 's2', range: '2 to 19', strand: 'minus' },
      { n: 5, run: 'Run 1 First run', program: 'BLASTN', query: 'q3', subject: 's1', range: '1 to 1', strand: 'unknown' },
      { n: 8, run: 'Run 2 Second run', program: 'TBLASTN', query: 'p1', subject: 's1', range: '1 to 25', strand: 'plus' },
      { n: 10, run: 'Run 2 Second run', program: 'TBLASTN', query: 'p1', subject: 's2', range: '2 to 19', strand: 'minus' },
    ];
    for (const item of expected) {
      const line = row(page, item.n);
      await expect(field(line, 'n')).toHaveText(String(item.n));
      await expect(field(line, 'run')).toHaveText(item.run);
      await expect(field(line, 'program')).toHaveText(item.program);
      await expect(field(line, 'query')).toHaveText(item.query);
      await expect(field(line, 'subject')).toHaveText(item.subject);
      await expect(field(line, 'range')).toHaveText(item.range);
      await expect(field(line, 'strand')).toHaveText(item.strand);
      // The outfmt 6 row's E value and bit score as written (the FakeEngine writes FAKE).
      await expect(field(line, 'evalue')).toHaveText('FAKE');
      await expect(field(line, 'bitscore')).toHaveText('FAKE');
    }
    // New candidates are selected.
    await expect(page.getByTestId('candidates-selected')).toHaveText('11 of 11 candidates selected');

    // Origins: one entry per run, with what it searched.
    const origins = page.getByTestId('candidate-origins');
    await expect(origins.locator('li[data-testid^="candidate-origin-"]')).toHaveCount(2);
    const first = page.getByTestId('candidate-origin-1');
    await expect(first.locator('h4')).toContainText('Run 1 · First run');
    await expect(first.locator('h4')).toContainText('7 candidates');
    await expect(first.locator('[data-detail="program"]')).toHaveText('BLASTN');
    await expect(first.locator('[data-detail="query"]')).toContainText('query.fa');
    await expect(first.locator('[data-detail="query"]')).toContainText(/SHA-256 [0-9a-f]{64}/);
    await expect(first.locator('[data-detail="subject"]')).toContainText(/subject\.fa\s*SHA-256 [0-9a-f]{64}/);
    await expect(first.locator('[data-detail="build"]')).toHaveText('fake-engine');
    await expect(first.locator('[data-detail="ended"]')).toHaveText(/^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$/);
    const second = page.getByTestId('candidate-origin-2');
    await expect(second.locator('h4')).toContainText('Run 2 · Second run');
    await expect(second.locator('h4')).toContainText('4 candidates');
    await expect(second.locator('[data-detail="program"]')).toHaveText('TBLASTN');
  });

  test('notes, orders, moves and removal; "Show in results" lands on the HSP and says which view filters it cleared', async ({ page }) => {
    await blastnRun(page);
    await page.getByTestId('tab-results').click();
    await addAllOfQuery(page, '4 HSPs added to Candidates.');
    await tblastnRun(page, 2, 'Second run');
    await chooseRun(page, 2);
    await addAllOfQuery(page, '4 HSPs added to Candidates.');
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidate-list')).toHaveAttribute('data-count', '8');

    // A note stays with its candidate, also when the tray is shown again and sorted.
    await page.getByTestId('candidate-note-2').fill('check the second exon');
    await page.getByTestId('tab-search').click();
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidate-note-2')).toHaveValue('check the second exon');

    // Subject: the record's name, then its position on the record, then the run.
    await page.getByTestId('candidates-sort').selectOption('subject');
    expect(await trayOrder(page)).toEqual(['1/0:0', '2/0:0', '1/0:1', '2/0:1', '1/0:2', '2/0:2', '1/0:3', '2/0:3']);
    await expect(page.getByTestId('candidate-note-3')).toHaveValue('check the second exon');
    // Run: the run's number, then the engine's order.
    await page.getByTestId('candidates-sort').selectOption('run');
    expect(await trayOrder(page)).toEqual(['1/0:0', '1/0:1', '1/0:2', '1/0:3', '2/0:0', '2/0:1', '2/0:2', '2/0:3']);
    await page.getByTestId('candidates-sort').selectOption('added');

    // Up and Down move one place; the focus stays on the moved row's button; the order is the user's.
    await expect(page.getByTestId('candidate-up-1')).toHaveAttribute('aria-disabled', 'true');
    await page.getByTestId('candidate-down-1').click();
    expect((await trayOrder(page)).slice(0, 2)).toEqual(['1/0:1', '1/0:0']);
    await expect(page.getByTestId('candidate-down-2')).toBeFocused();
    await expect(page.getByTestId('candidates-sort')).toHaveValue('custom');
    await page.keyboard.press('Enter');
    expect((await trayOrder(page)).slice(0, 3)).toEqual(['1/0:1', '1/0:2', '1/0:0']);
    await page.getByTestId('candidate-up-3').click();
    await expect(page.getByTestId('candidate-up-2')).toBeFocused();
    expect((await trayOrder(page)).slice(0, 3)).toEqual(['1/0:1', '1/0:0', '1/0:2']);

    // Remove one; then the selected ones.
    await page.getByTestId('candidate-remove-3').click();
    await expect(page.getByTestId('candidate-list')).toHaveAttribute('data-count', '7');
    await expect(tabCount(page)).toHaveText('7');
    await expect(page.getByTestId('candidate-remove-3')).toBeFocused();
    expect(await trayOrder(page)).not.toContain('1/0:2');
    await page.getByTestId('candidates-select-all').uncheck();
    await expect(page.getByTestId('candidates-selected')).toHaveText('0 of 7 candidates selected');
    await expect(page.getByTestId('extract-download')).toBeDisabled();
    await page.getByTestId('candidate-mark-1').check();
    await page.getByTestId('candidate-mark-2').check();
    await page.getByTestId('candidates-remove-selected').click();
    await expect(tabCount(page)).toHaveText('5');
    expect(await trayOrder(page)).toEqual(['1/0:3', '2/0:0', '2/0:1', '2/0:2', '2/0:3']);

    // "Show in results" of a candidate of run 2 that the view filters hide: the filter is cleared and said so.
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '2');
    await page.getByTestId('filter-subject').fill('s3');
    await page.getByTestId('filter-apply').click();
    await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', '1');
    await page.getByTestId('tab-candidates').click();
    await page.getByTestId('candidate-reveal-3').click();
    await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('results-view-hits')).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', '0:1');
    await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', '1');
    await expect(page.getByTestId('revealed-message')).toHaveText(
      'The view filters hid HSP 1.2 of run 2, so these were cleared: Subject ID contains.',
    );
    await expect(page.getByTestId('range-0-1').getByTestId('range-label')).toBeFocused();
    await expect(page.getByTestId('filter-subject')).toHaveValue('');
    // Another HSP selected: the message goes.
    await page.getByTestId('range-0-1').getByTestId('range-next').click();
    await expect(page.getByTestId('revealed-message')).toHaveCount(0);

    // A candidate of run 1 opens run 1 (no filter hid it).
    await page.getByTestId('tab-candidates').click();
    await page.getByTestId('candidate-reveal-1').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', '0:3');
    await expect(page.getByTestId('revealed-message')).toHaveCount(0);
    await expect(page.getByTestId('range-0-3').getByTestId('range-label')).toBeFocused();
  });

  test('downloads: the residues of the input at the HSPs (regions, flanks cut at the ends, spanning, complete, queries) and the aligned rows of the HSP records', async ({
    page,
  }) => {
    await blastnRun(page);
    await page.getByTestId('tab-results').click();
    await addAllOfQuery(page, '4 HSPs added to Candidates.');
    // q3's HSP of one letter on s1.
    await page.getByTestId('query-row-2').click();
    await page.getByTestId('subject-mark-0').check();
    await page.getByTestId('descriptions-add-candidates').click();
    await expect(tabCount(page)).toHaveText('5');
    await page.getByTestId('tab-candidates').click();

    const header = (id: string, from: number, to: number, role: string, k: number, length: number, hsps: string, strand: string, requested = '') =>
      `${id}:${from}-${to} run=1 ${role}_record=${k} length=${length} unit=nt hsps=${hsps} hit_strand=${strand}${requested === '' ? '' : ` requested=${requested}`}`;
    const check = (records: FastaRecord[], expected: { header: string; letters: string }[]) => {
      expect(records.map((r) => r.header)).toEqual(expected.map((e) => e.header));
      records.forEach((r, i) => {
        expect(r.letters, r.header).toBe(expected[i]!.letters);
        // 60 letters per line.
        expect(r.lines.every((line, j) => line.length === 60 || (j === r.lines.length - 1 && line.length <= 60)), r.header).toBe(true);
      });
    };

    // The hit regions on the subjects (the default).
    let file = await saved(page, 'extract-download');
    expect(file.name).toBe('losat-candidates.fa');
    check(parseFasta(file.text), [
      { header: header('s1', 1, 25, 'subject', 1, 50, '1.1', 'plus'), letters: slice(S.s1, 1, 25) },
      { header: header('s2', 1, 18, 'subject', 2, 37, '1.2', 'plus'), letters: slice(S.s2, 1, 18) },
      { header: header('s2', 2, 19, 'subject', 2, 37, '1.3', 'minus'), letters: slice(S.s2, 2, 19) },
      { header: header('s3', 1, 40, 'subject', 3, 80, '1.4', 'plus'), letters: slice(S.s3, 1, 40) },
      { header: header('s1', 1, 1, 'subject', 1, 50, '3.1', 'unknown'), letters: slice(S.s1, 1, 1) },
    ]);
    const summary = page.getByTestId('extract-summary');
    await expect(summary).toHaveAttribute('data-output', 'sequences');
    await expect(summary).toContainText('Saved losat-candidates.fa: 5 sequences of subject records from 5 candidates');
    await expect(page.getByTestId('extract-clipped')).toHaveCount(0);
    await expect(page.getByTestId('extract-unknown-strand')).toHaveCount(1);
    await expect(page.getByTestId('extract-unknown-strand')).toContainText('HSP 3.1 of run 1 covers one letter of subject record 1');

    // With flanks: the records' ends cut them, and the header and the summary keep the request.
    await page.getByTestId('extract-region-flanked').check();
    await page.getByTestId('extract-flank-left').fill('x');
    await expect(page.getByTestId('extract-flank-error')).toBeVisible();
    await expect(page.getByTestId('extract-download')).toBeDisabled();
    await page.getByTestId('extract-flank-left').fill('5');
    await page.getByTestId('extract-flank-right').fill('30');
    file = await saved(page, 'extract-download');
    check(parseFasta(file.text), [
      { header: header('s1', 1, 50, 'subject', 1, 50, '1.1', 'plus', '-4-55'), letters: S.s1 },
      { header: header('s2', 1, 37, 'subject', 2, 37, '1.2', 'plus', '-4-48'), letters: S.s2 },
      { header: header('s2', 1, 37, 'subject', 2, 37, '1.3', 'minus', '-3-49'), letters: S.s2 },
      { header: header('s3', 1, 70, 'subject', 3, 80, '1.4', 'plus', '-4-70'), letters: slice(S.s3, 1, 70) },
      { header: header('s1', 1, 31, 'subject', 1, 50, '3.1', 'unknown', '-4-31'), letters: slice(S.s1, 1, 31) },
    ]);
    await expect(page.getByTestId('extract-clipped')).toHaveCount(5);
    await expect(page.getByTestId('extract-clipped').first()).toHaveText('s1 (run 1, subject record 1, HSP 1.1): requested -4 to 55 nt, written 1 to 50 nt of 50 nt');

    // One region spanning the HSPs of each record (the order of the records' first HSPs).
    await page.getByTestId('extract-region-hit').check();
    await page.getByTestId('extract-join-spanning').check();
    file = await saved(page, 'extract-download');
    check(parseFasta(file.text), [
      { header: header('s1', 1, 25, 'subject', 1, 50, '1.1,3.1', 'mixed'), letters: slice(S.s1, 1, 25) },
      { header: header('s2', 1, 19, 'subject', 2, 37, '1.2,1.3', 'mixed'), letters: slice(S.s2, 1, 19) },
      { header: header('s3', 1, 40, 'subject', 3, 80, '1.4', 'plus'), letters: slice(S.s3, 1, 40) },
    ]);

    // The complete sequences, once per record.
    await page.getByTestId('extract-region-whole').check();
    await expect(page.getByTestId('extract-join-spanning')).toBeDisabled();
    file = await saved(page, 'extract-download');
    check(parseFasta(file.text), [
      { header: header('s1', 1, 50, 'subject', 1, 50, '1.1,3.1', 'mixed'), letters: S.s1 },
      { header: header('s2', 1, 37, 'subject', 2, 37, '1.2,1.3', 'mixed'), letters: S.s2 },
      { header: header('s3', 1, 80, 'subject', 3, 80, '1.4', 'plus'), letters: S.s3 },
    ]);

    // The queries' hit regions.
    await page.getByTestId('extract-region-hit').check();
    await page.getByTestId('extract-join-separate').check();
    await page.getByTestId('extract-role-query').check();
    file = await saved(page, 'extract-download');
    check(parseFasta(file.text), [
      { header: header('q1', 1, 20, 'query', 1, 40, '1.1', 'plus'), letters: slice(Q.q1, 1, 20) },
      { header: header('q1', 1, 20, 'query', 1, 40, '1.2', 'plus'), letters: slice(Q.q1, 1, 20) },
      { header: header('q1', 2, 21, 'query', 1, 40, '1.3', 'plus'), letters: slice(Q.q1, 2, 21) },
      { header: header('q1', 1, 20, 'query', 1, 40, '1.4', 'plus'), letters: slice(Q.q1, 1, 20) },
      { header: header('q3', 1, 1, 'query', 3, 24, '3.1', 'unknown'), letters: slice(Q.q3, 1, 1) },
    ]);
    await expect(summary).toContainText('5 sequences of query records');

    // The aligned rows of the HSP records, as written, in a file of their own; HSP 1.3's record has none.
    file = await saved(page, 'extract-aligned');
    expect(file.name).toBe('losat-candidates-aligned.fa');
    const aligned = (q: [number, number], s: [number, number], hsp: string, qid: string, k: number, sid: string, sk: number) => {
      const width = Math.max(q[1] - q[0] + 1, s[1] - s[0] + 1);
      return [
        { header: `${qid}:${q[0]}-${q[1]} run=1 query_record=${k} hsp=${hsp} aligned`, letters: fakeRow(q[1] - q[0] + 1, width) },
        { header: `${sid}:${s[0]}-${s[1]} run=1 subject_record=${sk} hsp=${hsp} aligned`, letters: fakeRow(s[1] - s[0] + 1, width) },
      ];
    };
    check(parseFasta(file.text), [
      ...aligned([1, 20], [1, 25], '1.1', 'q1', 1, 's1', 1),
      ...aligned([1, 20], [1, 18], '1.2', 'q1', 1, 's2', 2),
      ...aligned([1, 20], [1, 40], '1.4', 'q1', 1, 's3', 3),
      ...aligned([1, 1], [1, 1], '3.1', 'q3', 3, 's1', 1),
    ]);
    await expect(summary).toHaveAttribute('data-output', 'alignments');
    await expect(summary).toContainText('Saved losat-candidates-aligned.fa: 4 alignments from 5 candidates');
    await expect(page.getByTestId('extract-missing')).toHaveText(
      'HSP 1.3 of run 1 has no aligned sequences in its HSP record, so it has no gapped alignment to write.',
    );
  });

  test('a phone: the Descriptions\' marks and the tray fit the screen; a row shows #, Subject, Range, the note and the actions first', async ({ page }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await blastnRun(page);
    await page.getByTestId('tab-results').click();
    await page.getByTestId('subject-mark-1').check();
    await page.getByTestId('descriptions-add-candidates').click();
    await expect(confirmation(page)).toHaveText('2 HSPs added to Candidates.');
    const noSideScroll = () => page.evaluate(() => document.documentElement.scrollWidth <= document.documentElement.clientWidth);
    expect(await noSideScroll()).toBe(true);
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidate-list')).toHaveAttribute('data-count', '2');
    expect(await noSideScroll()).toBe(true);
    const line = row(page, 2);
    const boxes = await Promise.all(
      [field(line, 'n'), field(line, 'subject'), field(line, 'range'), page.getByTestId('candidate-note-2'), page.getByTestId('candidate-remove-2')].map(
        (locator) => locator.boundingBox(),
      ),
    );
    for (const box of boxes) {
      expect(box).not.toBeNull();
      expect(box!.x).toBeGreaterThanOrEqual(0);
      expect(box!.x + box!.width).toBeLessThanOrEqual(390);
    }
    // #, Subject and Range on the first line; the note under them, the actions under the note, the other values last.
    type Box = NonNullable<(typeof boxes)[number]>;
    const [n, subject, range, note, remove] = boxes.map((box) => box!) as [Box, Box, Box, Box, Box];
    expect(Math.abs(subject.y - n.y)).toBeLessThan(4);
    expect(Math.abs(range.y - n.y)).toBeLessThan(4);
    expect(note.y).toBeGreaterThan(n.y + n.height - 1);
    expect(remove.y).toBeGreaterThan(note.y + note.height - 1);
    const rest = (await field(line, 'evalue').boundingBox())!;
    expect(rest.y).toBeGreaterThan(remove.y + remove.height - 1);
    await expect(field(line, 'subject')).toHaveText('s2');
    await expect(field(line, 'range')).toHaveText('2 to 19');
    // The touch targets keep at least 24 px.
    const mark = (await page.getByTestId('candidate-mark-2').locator('..').boundingBox())!;
    expect(Math.min(mark.width, mark.height)).toBeGreaterThanOrEqual(24);
  });
});

// --- engine build --------------------------------------------------------------------------------

test.describe('engine build', () => {
  test.skip(!BUILD_HAS_ENGINE, 'real searches need the engine (LOSAT_WEB_REACTORS)');

  /** The HSPs of a run's outfmt 6 (Outputs view): subject ID and interval, strand, and query interval. */
  async function outfmt6(page: Page, run: number): Promise<string[][]> {
    await showOutput(page, run, 6, false);
    const rows = outfmt6Rows((await page.getByTestId('result-output').textContent()) ?? '');
    await show(page, 'results-view-hits');
    return rows;
  }
  const span = (start: string, end: string) => [Math.min(Number(start), Number(end)), Math.max(Number(start), Number(end))] as const;
  /** `id:from-to run=N role_record=K ... hit_strand=S` of a sequence header. */
  function parseHeader(header: string) {
    const match = /^(\S+):(\d+)-(\d+) run=(\d+) (query|subject)_record=(\d+) length=(\d+) unit=(nt|aa) hsps=(\S+) hit_strand=(\w+)/.exec(header);
    expect(match, header).not.toBeNull();
    const [, id, from, to, run, role, k, length, unit, hsps, strand] = match!;
    return { id: id!, from: Number(from), to: Number(to), run: Number(run), role: role!, k: Number(k), length: Number(length), unit: unit!, hsps: hsps!, strand: strand! };
  }

  test('BLASTN: a hit on the minus strand of one of two subject records with the same ID; the saved residues are the file\'s letters at the outfmt 6 intervals', async ({
    page,
  }) => {
    // The query lies forward in the first "dup" record and reverse-complemented in the second;
    // the third record is unrelated. The flanks are lower case in the file.
    const query = dna(51, 150);
    const records = [
      { title: 'dup first copy', letters: dna(52, 90) + dna(53, 10).toLowerCase() + query + dna(54, 10).toLowerCase() + dna(55, 90) },
      { title: 'dup second copy', letters: dna(56, 70) + reverseComplement(query) + dna(57, 130) },
      { title: 'other', letters: dna(58, 300) },
    ];
    await program(page, 'blastn');
    await openFiles(page, 'query', [{ name: 'query.fa', text: record('q1', query, 70) }]);
    await openFiles(page, 'subject', [{ name: 'dups.fa', text: records.map((r) => record(r.title, r.letters, 70)).join('') }]);
    await search(page, 1, 'Duplicate IDs');
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    const hits = await outfmt6(page, 1);
    expect(hits.length).toBeGreaterThanOrEqual(2);
    await addAllOfQuery(page, `${hits.length} HSPs added to Candidates.`);
    await page.getByTestId('tab-candidates').click();

    // Hit regions: one sequence per HSP, at the outfmt 6 subject interval, of the right record.
    let file = await saved(page, 'extract-download');
    let pieces = parseFasta(file.text).map((piece) => ({ ...parseHeader(piece.header), letters: piece.letters }));
    expect(pieces.map((p) => [p.id, p.from, p.to].join(' ')).sort()).toEqual(hits.map((h) => [h[1], ...span(h[8]!, h[9]!)].join(' ')).sort());
    for (const piece of pieces) {
      expect(piece.letters, `${piece.id} record ${piece.k}`).toBe(slice(records[piece.k - 1]!.letters, piece.from, piece.to));
    }
    const forward = pieces.find((p) => p.k === 1)!;
    const reverse = pieces.find((p) => p.k === 2)!;
    expect(forward.strand).toBe('plus');
    expect(reverse.strand).toBe('minus');
    // The right record of the two with the ID "dup": the query's letters, forward and reverse-complemented.
    expect(query).toContain(forward.letters.toUpperCase());
    expect(query).toContain(reverseComplement(reverse.letters).toUpperCase());

    // With flanks of 10: the lower-case letters of the file come with the first record's hit.
    await page.getByTestId('extract-region-flanked').check();
    await page.getByTestId('extract-flank-left').fill('10');
    await page.getByTestId('extract-flank-right').fill('10');
    file = await saved(page, 'extract-download');
    pieces = parseFasta(file.text).map((piece) => ({ ...parseHeader(piece.header), letters: piece.letters }));
    for (const piece of pieces) expect(piece.letters).toBe(slice(records[piece.k - 1]!.letters, piece.from, piece.to));
    expect(pieces.find((p) => p.k === 1)!.letters).toMatch(/[acgt]{10}$/);

    // The aligned rows: the HSP records' rows, as the outfmt 0 section shows them (Query and Sbjct lines).
    file = await saved(page, 'extract-aligned');
    const alignedRecords = parseFasta(file.text);
    expect(alignedRecords).toHaveLength(2 * hits.length);
    for (let i = 0; i < alignedRecords.length; i += 2) {
      const queryRow = alignedRecords[i]!.letters;
      const subjectRow = alignedRecords[i + 1]!.letters;
      expect(queryRow.length).toBe(subjectRow.length);
      expect(alignedRecords[i + 1]!.header).toMatch(/^dup:\d+-\d+ run=1 subject_record=[12] hsp=1\.\d+ aligned$/);
    }
  });

  test('TBLASTN: the saved subject residues (nt) and query residues (aa) are the files\' letters at the outfmt 6 intervals', async ({ page }) => {
    const queryFile = fasta('outfmt0/e2e_protein_query.faa');
    const subjectFile = fasta('outfmt0/e2e_amb_subject.fna');
    await program(page, 'tblastn');
    await openFiles(page, 'query', [{ name: 'e2e_protein_query.faa', text: queryFile }]);
    await openFiles(page, 'subject', [{ name: 'e2e_amb_subject.fna', text: subjectFile }]);
    await search(page, 1);
    await page.getByTestId('tab-results').click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    const queries = fileRecords(queryFile.toString('utf8'));
    const subjects = fileRecords(subjectFile.toString('utf8'));
    const hits = await outfmt6(page, 1);
    // The query that the results screen selects first: its rows.
    const qid = (await page.getByTestId('query-list').locator('[aria-pressed="true"]').getAttribute('data-testid'))!;
    const qIdx = Number(qid.replace('query-row-', ''));
    const ofQuery = hits.filter((h) => h[0] === queries[qIdx]!.id);
    expect(ofQuery.length).toBeGreaterThan(0);
    await addAllOfQuery(page, `${ofQuery.length} ${ofQuery.length === 1 ? 'HSP' : 'HSPs'} added to Candidates.`);
    await page.getByTestId('tab-candidates').click();

    let file = await saved(page, 'extract-download');
    let pieces = parseFasta(file.text).map((piece) => ({ ...parseHeader(piece.header), letters: piece.letters }));
    expect(pieces.map((p) => [p.id, p.from, p.to].join(' ')).sort()).toEqual(ofQuery.map((h) => [h[1], ...span(h[8]!, h[9]!)].join(' ')).sort());
    for (const piece of pieces) {
      expect(piece.unit).toBe('nt');
      expect(piece.letters).toBe(slice(subjects[piece.k - 1]!.letters, piece.from, piece.to));
    }

    await page.getByTestId('extract-role-query').check();
    file = await saved(page, 'extract-download');
    pieces = parseFasta(file.text).map((piece) => ({ ...parseHeader(piece.header), letters: piece.letters }));
    expect(pieces.map((p) => [p.from, p.to].join(' ')).sort()).toEqual(ofQuery.map((h) => span(h[6]!, h[7]!).join(' ')).sort());
    for (const piece of pieces) {
      expect(piece.unit).toBe('aa');
      expect(piece.k).toBe(qIdx + 1);
      expect(piece.letters).toBe(slice(queries[qIdx]!.letters, piece.from, piece.to));
    }
  });
});
