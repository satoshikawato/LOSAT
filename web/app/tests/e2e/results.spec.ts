// The results screen (S13, W4; plan §5.7, design §11): Run -> Query -> Subject -> HSP through
// the real application, in NCBI BLAST's order and words since W4b (docs/web/ncbi_ui_mapping.md
// §2): the header block, "Results for", and the tabs Descriptions, Graphic Summary, Alignments,
// Dot Plot, Run details and Outputs. The tests run with the FakeEngine build and with the engine
// build (LOSAT_WEB_REACTORS); what only the engine can show (outfmt 0's frames, a run that the
// engine refuses, the verification badge of a real run, the bins of real bit scores) is checked in
// the engine build. The lists, the Alignments, the Graphic Summary and the dot plot are compared
// with the stored outputs as the Outputs view shows them (the outfmt 6 rows and the outfmt 0
// text), never with values computed here.
//
// Every program that LOSAT Web runs is covered: BLASTN, BLASTP, TBLASTN and TBLASTX. BLASTX
// is not in ABI v2 until session SX (docs/losat_web_gui_sessions/session_sx_blastx_integration.md;
// the search screen refuses it, search.spec.ts), so it has no results screen to test yet.
import { readFileSync } from 'node:fs';
import { basename, join } from 'node:path';
import { expect, test, type Locator, type Page } from '@playwright/test';
import { generateVerificationTable } from '../../build/verification';
import type { OutputFormat } from '../../src/domain/output-format';
import { optionKey } from '../../src/domain/verification';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { fasta, openFiles, openParameters, paste, program, showOutput, submit, task, waitStatus } from './support/search';
import { repeats } from './support/synthetic';

type ProgramId = 'blastn' | 'blastp' | 'tblastn' | 'tblastx';
type Unit = 'nt' | 'aa';
type Orientation = 'forward' | 'reverse' | 'unknown';
type Tab = 'hits' | 'graphic' | 'alignment' | 'dotplot' | 'details' | 'outputs';

const FORMATS = [0, 6, 7] as const satisfies readonly OutputFormat[];
const LABELS: Readonly<Record<ProgramId, string>> = { blastn: 'BLASTN', blastp: 'BLASTP', tblastn: 'TBLASTN', tblastx: 'TBLASTX' };
/** The tabs' test IDs (W4's: the hits view, the HSP's panes, the run's views). */
const TABS: Readonly<Record<Tab, string>> = {
  hits: 'results-view-hits',
  graphic: 'results-view-graphic',
  alignment: 'pane-alignment',
  dotplot: 'pane-dotplot',
  details: 'results-view-details',
  outputs: 'results-view-outputs',
};

// Engine searches are slower in Firefox and WebKit than in Chromium (W1 README).
test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await page.goto('/');
});

// --- inputs ----------------------------------------------------------------------------------

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

const A = dna(11, 60);
const S1 = `${A}${dna(12, 60)}${A}${dna(13, 60)}`;
/**
 * BLASTN subjects for the selection, filter and many-query tests: s1 holds the 60 letters A
 * twice, s2 once, s3 not at all.
 */
const PANEL_SUBJECTS = `>s1 first subject\n${S1}\n>s2 second subject\n${dna(14, 60)}${A}${dna(15, 60)}\n>s3 third subject\n${dna(16, 120)}\n`;
/**
 * 150 queries. Every fourth (#4, #8, ...; as the FakeEngine leaves every fourth query
 * without hits) is unrelated to the subjects and has no hits; the others are 60 letters of s1.
 * The first is A, so it has two HSPs on s1 and one on s2.
 */
const PANEL_QUERIES = Array.from({ length: 150 }, (_, i) => {
  const offset = (i * 7) % (S1.length - 59);
  return `>rec${i + 1}\n${i % 4 === 3 ? dna(1000 + i, 60) : S1.slice(offset, offset + 60)}\n`;
}).join('');

/**
 * Searches of the outfmt 0 fixtures of LOSAT/tests/outfmt0_manifest.tsv with the program's
 * default options: multi.megablast (both strands, a subject with an HSP on each, a query
 * without hits), method.blastp, e2e.tblastn.amb (minus frames; inputs large enough for Auto
 * to choose the threaded runtime) and tblastx.ambig.0 (frames of opposite signs).
 */
const PROGRAM_CASES: readonly {
  readonly id: ProgramId;
  readonly query: string;
  readonly subject: string;
  readonly units: { readonly query: Unit; readonly subject: Unit };
  /** The outfmt 0 text of the frames that the HSP list shows as `q/s` (translated programs). */
  readonly frameLine?: (query: string, subject: string) => string;
}[] = [
  { id: 'blastn', query: 'outfmt0/multi_query.fasta', subject: 'outfmt0/multi_subject.fasta', units: { query: 'nt', subject: 'nt' } },
  { id: 'blastp', query: 'outfmt0/e2e_protein_query.faa', subject: 'outfmt0/e2e_protein_subject.faa', units: { query: 'aa', subject: 'aa' } },
  {
    id: 'tblastn',
    query: 'outfmt0/e2e_protein_query.faa',
    subject: 'outfmt0/e2e_amb_subject.fna',
    units: { query: 'aa', subject: 'nt' },
    frameLine: (_query, subject) => ` Frame = ${subject}\n`,
  },
  {
    id: 'tblastx',
    query: 'outfmt0/tblastx_ambig_query.fasta',
    subject: 'outfmt0/tblastx_ambig_subject.fasta',
    units: { query: 'nt', subject: 'nt' },
    frameLine: (query, subject) => ` Frame = ${query}/${subject}\n`,
  },
];

// --- helpers ---------------------------------------------------------------------------------

const count = (n: number) => n.toLocaleString('en-US');
const plural = (n: number, one: string) => `${count(n)} ${n === 1 ? one : `${one}s`}`;

/** The rows of an outfmt 6 text, split into their fields (the FakeEngine writes a comment line first). */
function outfmt6Rows(text: string): string[][] {
  return text
    .split('\n')
    .filter((line) => line !== '' && !line.startsWith('#'))
    .map((line) => line.split('\t'));
}

/** The records of a FASTA file: the ID (the title's first word) and the number of residues. */
function fastaRecords(bytes: Buffer): { id: string; length: number }[] {
  return bytes
    .toString('utf8')
    .split('>')
    .slice(1)
    .map((block) => {
      const [title = '', ...lines] = block.split('\n');
      return { id: title.trim().split(/\s+/)[0]!, length: lines.join('').replace(/\s/g, '').length };
    });
}

/** NCBI's "Alignment Scores" bin of a bit score as outfmt 6 writes it: < 40, 40 - 50, 50 - 80, 80 - 200, >= 200. */
function ncbiBin(bitscore: string): number {
  const bits = Number(bitscore);
  return bits < 40 ? 0 : bits < 50 ? 1 : bits < 80 ? 2 : bits < 200 ? 3 : 4;
}

/**
 * Whether a bit score as outfmt 6 writes it lies within half a unit of its last digit of a bin
 * edge (0.05 for one decimal, 0.5 for a whole number): the record's value that the Graphic
 * Summary bins may then lie on the other side of the edge (39.96 is written 40.0). The tests
 * have no hook to the record's value, so such HSPs are left out of the bin check.
 */
function nearBinEdge(bitscore: string): boolean {
  const decimals = /\.(\d+)$/.exec(bitscore)?.[1]?.length ?? 0;
  const half = 0.5 * 10 ** -decimals;
  return [40, 50, 80, 200].some((edge) => Math.abs(Number(bitscore) - edge) <= half);
}

/** NCBI's "Range n: a to b" of an HSP: its position among the subject's HSPs, and its outfmt 6 subject coordinates in ascending order. */
function rangeLabel(n: number, fields: readonly string[]): string {
  const [start, end] = [Number(fields[8]), Number(fields[9])];
  return `Range ${n}: ${Math.min(start, end)} to ${Math.max(start, end)}`;
}

async function text(locator: Locator): Promise<string> {
  return (await locator.textContent()) ?? '';
}

/** A text matched as it is in a regular expression. */
const literal = (value: string) => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

/** Shows a tab of the results screen. */
async function show(page: Page, tab: Tab): Promise<void> {
  const button = page.getByTestId(TABS[tab]);
  await button.click();
  await expect(button).toHaveAttribute('aria-pressed', 'true');
}

/** Queues the search set in the form and waits until it, run `number`, completes. */
async function run(page: Page, number: number): Promise<void> {
  await submit(page);
  await waitStatus(page, number, 'completed');
}

/** Opens a completed run's results from the queue and waits until they are read. */
async function openFromQueue(page: Page, number: number): Promise<void> {
  await page.getByTestId(`run-${number}-open`).click();
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(number));
}

/** The commands and stored outputs of the run on the results screen, read in its Outputs view; the Descriptions are shown after. */
async function readOutputs(page: Page, number: number): Promise<{ command: Record<OutputFormat, string>; text: Record<OutputFormat, string> }> {
  const command = {} as Record<OutputFormat, string>;
  const output = {} as Record<OutputFormat, string>;
  for (const format of FORMATS) {
    await showOutput(page, number, format, false);
    command[format] = await text(page.getByTestId('result-command'));
    output[format] = await text(page.getByTestId('result-output'));
  }
  await show(page, 'hits');
  return { command, text: output };
}

const subjectRows = (page: Page) => page.getByTestId('subject-list').locator('[data-testid^="subject-row-"]');
const hspRows = (page: Page) => page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"]');
const rangeBlocks = (page: Page) => page.getByTestId('alignments').locator('section.range-block');

/**
 * Watches the reads of outfmt 0 byte ranges that the page asks of the Data worker (its requests
 * `readOutputRange` of format 0, counted as they are sent, answered or not), and can make the
 * next ones fail (the worker answers that it has no such method), as a read of storage may.
 */
async function watchReads(page: Page): Promise<{ asked: () => Promise<number>; failNext: (n: number) => Promise<void> }> {
  await page.evaluate(() => {
    const w = window as unknown as { __reads?: number; __fail?: number };
    if (w.__reads !== undefined) return;
    w.__reads = 0;
    w.__fail = 0;
    const post = Worker.prototype.postMessage as (this: Worker, message: unknown, transfer: Transferable[]) => void;
    Worker.prototype.postMessage = function (this: Worker, message: unknown, transfer: Transferable[]) {
      const request = message as { method?: unknown; args?: unknown[] } | null;
      if (typeof request === 'object' && request !== null && request.method === 'readOutputRange' && request.args?.[1] === 0) {
        w.__reads!++;
        if (w.__fail! > 0) {
          w.__fail!--;
          message = { ...request, method: 'readOutputRange (made to fail by the test)' };
        }
      }
      post.call(this, message, transfer);
    } as typeof Worker.prototype.postMessage;
  });
  return {
    asked: () => page.evaluate(() => (window as unknown as { __reads: number }).__reads),
    failNext: async (n) => {
      await page.evaluate((count) => ((window as unknown as { __fail: number }).__fail = count), n);
    },
  };
}

/** Waits until the Alignments ask for no more sections and every section asked for has arrived. */
async function settledSections(page: Page, asked: () => Promise<number>): Promise<void> {
  let last = '';
  const state = async () =>
    `${await asked()}/${await page.getByTestId('alignments').locator('[data-testid="range-section"]:not([data-state="pending"])').count()}`;
  await expect
    .poll(
      async () => {
        const now = await state();
        const same = now === last;
        last = now;
        return same;
      },
      { intervals: [700], timeout: 30_000 },
    )
    .toBe(true);
}
const subjectIndex = async (row: Locator) => Number((await row.getAttribute('data-testid'))!.replace('subject-row-', ''));

/** The kinds of the notices shown, without the outfmt 0 notice (which the FakeEngine's third subject adds). */
async function noticeKinds(page: Page): Promise<string[]> {
  const kinds = await page.getByTestId('results-notice').evaluateAll((elements) => elements.map((e) => (e as HTMLElement).dataset['kind'] ?? ''));
  return kinds.filter((kind) => kind !== 'outfmt0-partial');
}

/** The HSP id ("<query index>:<rank>") of an HSP row. */
async function hspId(row: Locator): Promise<string> {
  const [, qIdx, rank] = /^hsp-row-(\d+)-(\d+)$/.exec((await row.getAttribute('data-testid')) ?? '')!;
  return `${qIdx}:${rank}`;
}

/**
 * Checks that the HSP `id` of subject `sIdx` is the one selection of every tab: the Alignments'
 * subject, its HSP table and detail, and the Descriptions' row. The Descriptions are shown after.
 */
async function expectSelected(page: Page, id: string, sIdx: number): Promise<void> {
  await show(page, 'alignment');
  await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', String(sIdx));
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', id);
  await expect(page.getByTestId(`range-${id.replace(':', '-')}`).getByTestId('hsp-detail')).toHaveCount(1);
  await expect(page.getByTestId(`hsp-row-${id.replace(':', '-')}`)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('hsp-list').locator('[aria-pressed="true"]')).toHaveCount(1);
  await show(page, 'hits');
  await expect(page.getByTestId(`subject-row-${sIdx}`)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('subject-list').locator('[aria-pressed="true"]')).toHaveCount(1);
}

/**
 * At the desktop size a table fits its pane (S13 screen review H1, M1): it does not scroll
 * sideways or say that it does, a value that its column cuts has its full text in its title, the
 * coordinates, frames and orientation are not cut, and the header of a column of right-aligned
 * values ends where its values end. The table is in the tab shown.
 */
async function expectTableFits(page: Page, table: 'subject-table' | 'hsp-table'): Promise<void> {
  const scroll = page.getByTestId(table).locator('.table-scroll');
  const [scrollWidth, clientWidth] = await scroll.evaluate((element) => [element.scrollWidth, element.clientWidth]);
  expect(scrollWidth, `${table} scrolls sideways`).toBeLessThanOrEqual(clientWidth);
  await expect(page.getByTestId(`${table}-scroll-hint`)).toHaveCount(0);
  const report = await scroll.evaluate((element) => {
    const cells = [...element.querySelectorAll<HTMLElement>('.table-row [data-field]')];
    const untitled = cells
      .filter((cell) => cell.scrollWidth > cell.clientWidth && cell.title.replace(/\s/g, '') !== (cell.textContent ?? '').replace(/\s/g, ''))
      .map((cell) => `${cell.dataset['field']}: ${cell.textContent}`);
    const cut = cells
      .filter((cell) => ['query', 'subject', 'frames', 'orientation'].includes(cell.dataset['field']!) && cell.scrollWidth > cell.clientWidth)
      .map((cell) => `${cell.dataset['field']}: ${cell.textContent}`);
    return { untitled, cut };
  });
  expect(report, table).toEqual({ untitled: [], cut: [] });
  expect(await misalignedHeaders(page, table)).toEqual([]);
}

/**
 * The headers of right-aligned values that do not end where the first row's value ends (S13
 * screen review M1), also when the rows have a scroll bar. Header and row cells are in the
 * same order.
 */
function misalignedHeaders(page: Page, table: string): Promise<string[]> {
  return page
    .getByTestId(table)
    .locator('.table-scroll')
    .evaluate((element) => {
      const textRight = (node: Element) => {
        const range = document.createRange();
        range.selectNodeContents(node);
        return range.getBoundingClientRect().right;
      };
      const head = [...element.querySelector('.table-head')!.children];
      const row = [...element.querySelector('.table-row')!.children];
      return row.flatMap((cell, i) => {
        if (!cell.classList.contains('num') || head[i] === undefined) return [];
        const [headRight, cellRight] = [textRight(head[i]!), textRight(cell)];
        return Math.abs(headRight - cellRight) > 2 ? [`${head[i]!.textContent?.trim()}: header ends at ${headRight}, values at ${cellRight}`] : [];
      });
    });
}

/**
 * Checks the Range blocks of the selected subject in the Alignments (W4b): the window of the first
 * Ranges (the selected HSP is the subject's first), each label NCBI's "Range n: a to b" of the HSP's
 * outfmt 6 row, in the engine's order; and the sections of the first `sections` other blocks,
 * read as they come into view, are outfmt 0's text.
 */
async function expectRanges(page: Page, queryRows: readonly string[][], pairRanks: readonly number[], out0: string, sections: number): Promise<void> {
  const blocks = rangeBlocks(page);
  await expect(blocks).toHaveCount(Math.min(pairRanks.length, 26));
  let checked = 0;
  for (const [position, block] of (await blocks.all()).entries()) {
    const rank = Number((await block.getAttribute('data-testid'))!.split('-')[2]);
    expect(rank).toBe(pairRanks[position]);
    await expect(block.getByTestId('range-label')).toHaveText(rangeLabel(position + 1, queryRows[rank]!));
    const section = block.getByTestId('range-section');
    if (checked >= sections || (await section.count()) === 0) continue;
    await section.scrollIntoViewIfNeeded();
    await expect(section).toHaveAttribute('data-state', 'ready');
    // An HSP's section of outfmt 0 starts with its Score line (an empty text would pass toContain).
    const body = await text(section);
    expect(body).toMatch(/^ Score = /);
    expect(out0).toContain(body);
    checked++;
  }
}

/**
 * The Graphic Summary of the selected query `qIdx` (W4b): one row per subject drawn, each drawn
 * HSP's bin NCBI's bin of its outfmt 6 bit score (the engine's values only), the popover of a
 * hovered HSP the outfmt 6 strings; a click selects the HSP and shows its Range in the Alignments.
 */
async function expectGraphicSummary(page: Page, queryRows: readonly string[][], qIdx: number, subjects: number): Promise<void> {
  await show(page, 'graphic');
  const canvas = page.getByTestId('graphic-canvas');
  const rows = Math.min(subjects, 100);
  await expect(canvas).toHaveAttribute('data-rows', String(rows));
  const drawn = Number(await canvas.getAttribute('data-hsps'));
  await expect(page.getByTestId('graphic-title')).toHaveText(`Distribution of ${count(drawn)} HSPs on ${count(rows)} subject sequences`);
  await expect(page.getByTestId('graphic-legend').locator('[data-bin]')).toHaveText(['< 40', '40 - 50', '50 - 80', '80 - 200', '>= 200']);
  const targets = JSON.parse((await canvas.getAttribute('data-targets'))!) as { hsp: string; x: number; y: number; bin: number }[];
  expect(targets.length).toBeGreaterThan(0);
  const rowOf = new Map<string, number>();
  for (const target of targets) {
    const [q, rank] = target.hsp.split(':').map(Number);
    expect(q).toBe(qIdx);
    const fields = queryRows[rank!]!;
    if (BUILD_HAS_ENGINE && !nearBinEdge(fields[11]!)) expect(target.bin, `${target.hsp}: bit score ${fields[11]}`).toBe(ncbiBin(fields[11]!));
    else expect([0, 1, 2, 3, 4]).toContain(target.bin);
    // One row per subject: the HSPs of a subject share a row, and no other subject's are on it.
    expect(rowOf.get(fields[1]!) ?? target.y, `the row of ${fields[1]}`).toBe(target.y);
    rowOf.set(fields[1]!, target.y);
  }
  expect(new Set(rowOf.values()).size).toBe(rowOf.size);

  // The last HSP drawn: the popover of the hovered bar, then a click.
  const target = targets.at(-1)!;
  const rank = Number(target.hsp.split(':')[1]);
  const fields = queryRows[rank]!;
  await canvas.hover({ position: { x: target.x, y: target.y } });
  const popover = page.getByTestId('graphic-popover');
  await expect(popover).toHaveAttribute('data-hsp', target.hsp);
  await expect(popover.locator('.graphic-popover-title')).toHaveText(fields[1]!);
  await expect(popover).toContainText(`HSP ${rank + 1} · Bit score ${fields[11]} · E value ${fields[10]}`);
  await canvas.click({ position: { x: target.x, y: target.y } });
  await expect(page.getByTestId('pane-alignment')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', target.hsp);
  const block = page.getByTestId(`range-${target.hsp.replace(':', '-')}`);
  await expect(block).toBeInViewport();
  await expect(block.getByTestId('range-label')).toBeFocused();
}

/** Checks that the page is no wider than the window; the message names the elements that reach beyond it. */
async function expectNoSideScroll(page: Page, state: string): Promise<void> {
  const { scrollWidth, innerWidth, beyond } = await page.evaluate(() => {
    const edge = window.innerWidth + 0.5;
    const beyond = [...document.body.querySelectorAll<HTMLElement>('*')]
      .filter((element) => element.getBoundingClientRect().right > edge && (element.parentElement?.getBoundingClientRect().right ?? 0) <= edge)
      .slice(0, 10)
      .map((e) => `${e.tagName.toLowerCase()}.${[...e.classList].join('.')}[${e.dataset['testid'] ?? ''}] to ${Math.round(e.getBoundingClientRect().right)}`);
    return { scrollWidth: document.documentElement.scrollWidth, innerWidth: window.innerWidth, beyond };
  });
  expect(scrollWidth, `${state}: the page is wider than the window (${beyond.join(', ')})`).toBeLessThanOrEqual(innerWidth);
}

/** The columns of a table whose headers are inside the table's box, from left to right. */
function columnsInView(page: Page, table: string): Promise<string[]> {
  return page
    .getByTestId(table)
    .locator('.table-scroll')
    .evaluate((scroll) => {
      const box = scroll.getBoundingClientRect();
      return [...scroll.querySelectorAll<HTMLElement>('.table-head > *')]
        .map((cell) => ({ col: cell.dataset['col'] ?? cell.textContent?.trim() ?? '', rect: cell.getBoundingClientRect() }))
        .filter(({ rect }) => rect.left >= box.left - 0.5 && rect.right <= box.right + 0.5)
        .sort((a, b) => a.rect.left - b.rect.left)
        .map(({ col }) => col);
    });
}

/** The orientation that an HSP's outfmt 6 coordinates give (ABI v2 §8), or undefined when only its frames could. */
function orientationOf(fields: readonly string[], units: { readonly query: Unit; readonly subject: Unit }, translated: boolean): Orientation | undefined {
  const direction = (start: number, end: number, unit: Unit): number | undefined => {
    if (start !== end) return start < end ? 1 : -1;
    if (unit === 'aa') return 1;
    return translated ? undefined : 0;
  };
  const q = direction(Number(fields[6]), Number(fields[7]), units.query);
  const s = direction(Number(fields[8]), Number(fields[9]), units.subject);
  if (q === undefined || s === undefined) return undefined;
  if (q === 0 || s === 0) return 'unknown';
  return q === s ? 'forward' : 'reverse';
}

// The verification table that the engine build's badges come from (build/verification.ts),
// and the engine's option grammar (the describe snapshot that the table is generated with).
const TABLE = BUILD_HAS_ENGINE ? generateVerificationTable(REPOSITORY) : undefined;
const DESCRIBE = JSON.parse(readFileSync(join(REPOSITORY, 'web/app/src/infra/fake/describe.json'), 'utf8')) as Record<
  string,
  { parameters: { flag: string; takes_value: boolean; default?: string }[] }
>;

/**
 * The badge level that the table gives a run: certified only when its option set was compared
 * with NCBI for outfmt 0, 6 and 7 and its runtime path and thread count were checked in the
 * browser (plan §6.1). It reads the table, not the application's badge function.
 */
function tableLevel(programId: ProgramId, argv: readonly string[], path: string, threads: number): 'certified' | 'outside' {
  const parameters = DESCRIBE[programId]!.parameters;
  const takesValue = new Map(parameters.map((p) => [p.flag, p.takes_value]));
  const task = parameters.find((p) => p.flag === '-task')?.default;
  const key = optionKey(argv, { takesValue: (flag) => takesValue.get(flag) ?? false, ...(task === undefined ? {} : { defaultTask: task }) });
  const verification = TABLE!.programs[programId];
  const compared = key === undefined ? [] : (verification?.optionSets[key] ?? []);
  const browser = verification?.browser;
  const checked =
    browser !== undefined &&
    FORMATS.every((format) => compared.includes(format) && browser.formats.includes(format)) &&
    (browser.paths as readonly string[]).includes(path) &&
    browser.threads.includes(threads);
  return checked ? 'certified' : 'outside';
}

// --- every program -----------------------------------------------------------------------------

for (const c of PROGRAM_CASES) {
  test(`${LABELS[c.id]}: the header block is the run's; the Descriptions, Graphic Summary and HSP table are the outfmt 6 rows, the Alignments outfmt 0's text; units, frames and run details`, async ({
    page,
  }) => {
    const translated = c.frameLine !== undefined;
    // The desktop size of the screen review (S13): the tables fit their pane.
    await page.setViewportSize({ width: 1280, height: 900 });
    await program(page, c.id);
    const [queryFile, subjectFile] = [fasta(c.query), fasta(c.subject)];
    await openFiles(page, 'query', [{ name: basename(c.query), text: queryFile }]);
    await openFiles(page, 'subject', [{ name: basename(c.subject), text: subjectFile }]);
    // A Job Title for one program (W4b): the header block shows it; the others have none.
    const title = c.id === 'blastn' ? 'Fixtures of outfmt 0, BLASTN' : '';
    if (title !== '') await page.getByTestId('job-title').fill(title);
    await run(page, 1);
    await openFromQueue(page, 1);
    await expect(page.getByTestId('results-run').locator('option:checked')).toHaveText(
      `Run 1 · ${LABELS[c.id]} · ${basename(c.query)} vs ${basename(c.subject)}`,
    );
    const outputs = await readOutputs(page, 1);
    const rows = outfmt6Rows(outputs.text[6]);
    const out0 = outputs.text[0];

    // The header block (NCBI's Job Title, RID, Program, Query ID ...) holds the run's snapshot.
    const [queries, subjects] = [fastaRecords(queryFile), fastaRecords(subjectFile)];
    if (title === '') await expect(page.getByTestId('results-job-title')).toHaveCount(0);
    else await expect(page.getByTestId('results-job-title')).toHaveText(title);
    const defaultTask = DESCRIBE[c.id]!.parameters.find((p) => p.flag === '-task')?.default;
    await expect(page.getByTestId('results-program')).toHaveText(defaultTask === undefined ? LABELS[c.id] : `${LABELS[c.id]} (task ${defaultTask})`);
    await expect(page.getByTestId('results-run-options')).toHaveText('defaults');
    await expect(page.getByTestId('results-summary').getByTestId('verification-badge')).toBeVisible();
    if (queries.length === 1) {
      await expect(page.getByTestId('results-query-id')).toHaveText(queries[0]!.id);
      await expect(page.getByTestId('results-query-length')).toHaveText(`${count(queries[0]!.length)} ${c.units.query}`);
      // "Results for" is for runs of several queries.
      await expect(page.getByTestId('query-list')).toHaveCount(0);
    } else {
      await expect(page.getByTestId('results-query-id')).toHaveCount(0);
      await expect(page.locator('.query-picker h3')).toHaveText('Results for');
      await expect(page.getByTestId('query-list')).toHaveAttribute('data-count', String(queries.length));
    }
    if (subjects.length === 1) {
      await expect(page.getByTestId('results-subject-id')).toHaveText(subjects[0]!.id);
      await expect(page.getByTestId('results-subject-length')).toHaveText(`${count(subjects[0]!.length)} ${c.units.subject}`);
    } else {
      await expect(page.getByTestId('results-subject-id')).toHaveCount(0);
      await expect(page.getByTestId('results-subjects')).toHaveText(`${basename(c.subject)}, ${plural(subjects.length, 'record')}`);
    }
    // "Download All" opens the Outputs.
    await page.getByTestId('results-download-all').click();
    await expect(page.getByTestId(TABS.outputs)).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('results-outputs')).toBeVisible();
    await show(page, 'hits');

    // The Descriptions: NCBI's heading and order of columns, the Subject ID last.
    await expect(page.getByTestId('subject-table').locator('.tool-band h3')).toContainText('Sequences producing significant alignments');
    expect(await page.getByTestId('subject-table').locator('.table-head > *').evaluateAll((cells) => cells.map((cell) => (cell as HTMLElement).dataset['col']))).toEqual([
      'order',
      'description',
      'bitscore',
      'evalue',
      'hsps',
      'length',
      'sseqid',
    ]);
    // Units of the lists.
    await expect(page.getByTestId('subject-sort-length')).toHaveText(`Length (${c.units.subject})`);
    await show(page, 'alignment');
    await expect(page.getByTestId('hsp-sort-qStart')).toHaveText(`Query (${c.units.query})`);
    await expect(page.getByTestId('hsp-sort-sStart')).toHaveText(`Subject (${c.units.subject})`);

    // The first query with hits is selected when the run opens; the first three are checked.
    const detail = page.getByTestId('hsp-detail');
    await expect(detail).toHaveAttribute('data-state', 'ready');
    const selectedQuery = Number((await detail.getAttribute('data-hsp'))!.split(':')[0]);
    const withHits: number[] = [];
    if (queries.length === 1) withHits.push(0);
    else {
      await expect(page.getByTestId(`query-row-${selectedQuery}`)).toHaveAttribute('aria-pressed', 'true');
      for (const row of await page.getByTestId('query-list').locator('[data-testid^="query-row-"]').all()) {
        if (!(await text(row)).includes('no hits')) withHits.push(Number((await row.getAttribute('data-testid'))!.replace('query-row-', '')));
      }
    }
    expect(withHits[0]).toBe(selectedQuery);
    const seen = { queries: 0, subjects: 0, hsps: 0, orientations: new Set<string>(), frames: new Set<string>() };
    let last = { qIdx: selectedQuery, rows: [] as string[][], subjects: 0 };
    for (const qIdx of withHits.slice(0, 3)) {
      const queryRow = page.getByTestId(`query-row-${qIdx}`);
      if (queries.length > 1) {
        await queryRow.click();
        await expect(queryRow).toHaveAttribute('aria-pressed', 'true');
      }
      await show(page, 'alignment');
      await expect(detail).toHaveAttribute('data-hsp', new RegExp(`^${qIdx}:`));
      await expect(detail).toHaveAttribute('data-state', 'ready');
      // The query's HSPs are its outfmt 6 rows, by rank.
      const qseqid = (await text(page.getByTestId('detail-row'))).split('\t')[0]!;
      const queryRows = rows.filter((fields) => fields[0] === qseqid);
      const sseqids = [...new Set(queryRows.map((fields) => fields[1]!))];
      expect(queryRows.length).toBeGreaterThan(0);
      if (queries.length > 1) await expect(queryRow).toContainText(`${c.units.query} · ${count(sseqids.length)} subj., ${count(queryRows.length)} HSPs`);
      last = { qIdx, rows: queryRows, subjects: sseqids.length };

      // The Descriptions: each subject's first outfmt 6 row, in the engine's order.
      await show(page, 'hits');
      await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', String(sseqids.length));
      const drawnSubjects = await subjectRows(page).count();
      expect(drawnSubjects).toBe(Math.min(sseqids.length, 22));
      for (let i = 0; i < drawnSubjects; i++) {
        const row = subjectRows(page).nth(i);
        const order = Number(await row.getAttribute('data-order'));
        expect(order).toBe(i + 1);
        const pair = queryRows.filter((fields) => fields[1] === sseqids[order - 1]);
        expect(await text(row.locator('[data-field="sseqid"]'))).toBe(sseqids[order - 1]);
        expect(await text(row.locator('[data-field="bitscore"]'))).toBe(pair[0]![11]);
        expect(await text(row.locator('[data-field="evalue"]'))).toBe(pair[0]![10]);
        expect((await text(row.locator('[data-field="hsps"]'))).trim()).toMatch(new RegExp(`^${count(pair.length)}(max)?$`));
      }

      // The first subjects in the Alignments: the HSP table holds every field of the HSPs' outfmt 6
      // rows, the Ranges their subject coordinates, the sections and the heading outfmt 0's text.
      for (let i = 0; i < Math.min(drawnSubjects, 3); i++) {
        await show(page, 'hits');
        const subject = subjectRows(page).nth(i);
        await subject.click();
        await expect(subject).toHaveAttribute('aria-pressed', 'true');
        const sIdx = await subjectIndex(subject);
        const sseqid = sseqids[i]!;
        const pairRanks = queryRows.flatMap((fields, rank) => (fields[1] === sseqid ? [rank] : []));
        const description = subject.locator('[data-field="description"]');
        await expect(description).not.toHaveText('…');
        const inOutfmt0 = (await text(description)) !== 'not in outfmt 0';
        const descriptionText = await text(description);
        const length = (await text(subject.locator('[data-field="length"]'))).replace(/,/g, '');

        await show(page, 'alignment');
        await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', String(sIdx));
        await expect(page.getByTestId('alignments-summary')).toHaveText(
          new RegExp(`^\\s*Sequence ID: ${sseqid.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}\\s+Length: ${count(Number(length))}\\s+Number of Matches: ${count(pairRanks.length)}\\s*$`),
        );
        await expect(page.getByTestId('hsp-list')).toHaveAttribute('data-count', String(pairRanks.length));
        await expect(hspRows(page)).toHaveCount(Math.min(pairRanks.length, 20));
        for (const row of await hspRows(page).all()) {
          const [, rank] = (await hspId(row)).split(':').map(Number);
          expect(pairRanks).toContain(rank);
          const fields = queryRows[rank!]!;
          const field = (name: string) => text(row.locator(`[data-field="${name}"]`));
          expect(await field('bitscore')).toBe(fields[11]);
          expect(await field('evalue')).toBe(fields[10]);
          expect(await field('pident')).toBe(fields[2]);
          expect(await field('length')).toBe(fields[3]);
          expect(await field('mismatch')).toBe(fields[4]);
          expect(await field('gapopen')).toBe(fields[5]);
          expect(await field('query')).toBe(`${fields[6]}–${fields[7]}`);
          expect(await field('subject')).toBe(`${fields[8]}–${fields[9]}`);
          expect(await field('outfmt0')).toBe(inOutfmt0 ? 'shown' : 'not shown');
          const orientation = orientationOf(fields, c.units, translated);
          if (orientation !== undefined) expect(await row.getAttribute('data-orientation')).toBe(orientation);
          if (translated) {
            expect(await field('frames')).toMatch(c.id === 'tblastn' ? /^–\/[+-][123]$/ : /^[+-][123]\/[+-][123]$/);
            // The frames' signs agree with the coordinates' directions.
            const [q, s] = (await field('frames')).split('/');
            const sign = (frame: string | undefined) => (frame === '–' ? 1 : frame!.startsWith('-') ? -1 : 1);
            expect(sign(q) === sign(s) ? 'forward' : 'reverse').toBe(await row.getAttribute('data-orientation'));
            seen.frames.add(await field('frames'));
          } else {
            await expect(row.locator('[data-field="frames"]')).toHaveCount(0);
          }
          seen.hsps++;
          seen.orientations.add((await row.getAttribute('data-orientation')) ?? '');
        }

        // The detail of the subject's first HSP (selected with the subject), in its Range.
        const first = hspRows(page).first();
        await expect(first).toHaveAttribute('aria-pressed', 'true');
        const id = await hspId(first);
        await expect(detail).toHaveAttribute('data-hsp', id);
        await expect(detail).toHaveAttribute('data-state', 'ready');
        await expect(page.getByTestId(`range-${id.replace(':', '-')}`).getByTestId('hsp-detail')).toHaveCount(1);
        const fields = queryRows[Number(id.split(':')[1])]!;
        expect((await text(page.getByTestId('detail-row'))).replace(/\n$/, '')).toBe(fields.join('\t'));
        // At 1280 px the row shows all its fields in the Range block: it does not run past its
        // box (W4b screen review M1: the bit score was cut).
        expect(await page.getByTestId('detail-row').evaluate((row) => row.scrollWidth - row.clientWidth), 'the outfmt 6 row overflows').toBeLessThanOrEqual(0);
        await expectRanges(page, queryRows, pairRanks, out0, 3);
        seen.subjects++;
        if (!inOutfmt0) {
          await expect(page.getByTestId('detail-not-in-outfmt0')).toBeVisible();
          await expect(page.getByTestId('detail-heading')).toHaveCount(0);
          continue;
        }
        const heading = await text(page.getByTestId('detail-heading'));
        const section = await text(page.getByTestId('detail-section'));
        expect(heading.length).toBeGreaterThan(0);
        expect(section.length).toBeGreaterThan(0);
        expect(out0).toContain(heading);
        expect(out0).toContain(section);
        // The description is the heading's title; the length is the heading's Length=.
        expect(heading).toContain(`${descriptionText}\nLength=${length}\n`);
        const orientation = await first.getAttribute('data-orientation');
        if (c.id === 'blastn') expect(section).toContain(` Strand=Plus/${orientation === 'reverse' ? 'Minus' : 'Plus'}\n`);
        if (translated && BUILD_HAS_ENGINE) {
          const [q, s] = (await text(first.locator('[data-field="frames"]'))).split('/');
          expect(section).toContain(c.frameLine!(q!, s!));
        }
      }
      seen.queries++;
    }
    await show(page, 'alignment');
    await expectTableFits(page, 'hsp-table');
    await show(page, 'hits');
    await expectTableFits(page, 'subject-table');
    // The "With hits only" box sits next to its label (S13 screen review L1).
    if (queries.length > 1) expect((await page.getByTestId('filter-hits-only').boundingBox())!.width).toBeLessThan(30);
    test.info().annotations.push({
      type: 'coverage',
      description:
        `${c.id}: ${seen.queries} queries, ${seen.subjects} subjects opened, ${seen.hsps} HSP rows; ` +
        `orientations ${[...seen.orientations].sort().join(' ')}${translated ? `; frames ${[...seen.frames].sort().join(' ')}` : ''}`,
    });

    // The Graphic Summary of the last query checked; its click selects an HSP and shows its Range.
    await expectGraphicSummary(page, last.rows, last.qIdx, last.subjects);

    // The dot plot of the selected pair, in the units of the records.
    const fields = (await text(page.getByTestId('detail-row'))).replace(/\n$/, '').split('\t');
    const selected = await detail.getAttribute('data-hsp');
    const frames = translated ? await text(page.getByTestId(`hsp-row-${selected!.replace(':', '-')}`).locator('[data-field="frames"]')) : '';
    await show(page, 'dotplot');
    const canvas = page.getByTestId('dotplot-canvas');
    await expect(canvas).toHaveAttribute('data-segments', (await page.getByTestId('hsp-list').getAttribute('data-count'))!);
    await expect(page.getByTestId('dotplot')).toContainText(new RegExp(`\\([\\d,]+ ${c.units.query}\\) against subject .* \\([\\d,]+ ${c.units.subject}\\)\\.`));
    await expect(page.getByTestId('dotplot-selected')).toContainText(
      `query ${fields[6]}–${fields[7]} ${c.units.query}, subject ${fields[8]}–${fields[9]} ${c.units.subject}` +
        (translated ? `, frames ${frames.replace('/', ' / ')}.` : '.'),
    );
    // The same scale on both axes, a residue counting as 3 nt against nucleotides (S13b decision
    // 26: TBLASTN's query), unless the shorter side was raised to 120 px (then the note says so).
    const [queryLength, subjectLength] = /\(([\d,]+) \w+\) against subject .* \(([\d,]+) \w+\)\./
      .exec(await text(page.getByTestId('dotplot')))!
      .slice(1)
      .map((n) => Number(n.replace(/,/g, '')));
    const ratio = (queryLength! * (c.units.query === 'aa' && c.units.subject === 'nt' ? 3 : 1)) / subjectLength!;
    const [width, height] = (await canvas.getAttribute('data-plot'))!.split('x').map(Number);
    if ((await page.getByTestId('dotplot-scale-note').count()) === 0) {
      expect(Math.min(Math.abs(width! - height! * ratio), Math.abs(height! - width! / ratio))).toBeLessThanOrEqual(0.5);
    } else expect(Math.min(width!, height!)).toBe(120);
    // The axis titles as drawn: the IDs and the units of the visible spans.
    const titles = async () => JSON.parse((await canvas.getAttribute('data-titles'))!) as [string, string];
    expect((await titles())[1]).toMatch(new RegExp(`^Subject ${literal(fields[1]!)} \\((${c.units.subject === 'aa' ? 'aa|kaa' : 'bp|kbp'})\\)$`));
    if (c.id === 'tblastn') {
      // TBLASTN's pair is to scale at 1280 px, its protein axis drawn at 3 nt per aa, which its title
      // says; on a phone the 120 px minimum changes the scale, and the note says so instead (W4b
      // screen review finding 4).
      await expect(page.getByTestId('dotplot-scale-note')).toHaveCount(0);
      expect((await titles())[0]).toBe(`Query ${fields[0]} (aa; drawn at 3 nt per aa)`);
      await page.setViewportSize({ width: 390, height: 844 });
      await expect(page.getByTestId('dotplot-scale-note')).toBeVisible();
      await expect.poll(async () => (await titles())[0]).toBe(`Query ${fields[0]} (aa)`);
      await page.setViewportSize({ width: 1280, height: 900 });
      await expect(page.getByTestId('dotplot-scale-note')).toHaveCount(0);
    } else {
      expect((await titles())[0]).toMatch(new RegExp(`^Query ${literal(fields[0]!)} \\((${c.units.query === 'aa' ? 'aa|kaa' : 'bp|kbp'})\\)$`));
    }

    // Run details: the command of each format is the Outputs view's, and the verification badge.
    await show(page, 'details');
    for (const format of FORMATS) await expect(page.getByTestId(`run-command-${format}`)).toHaveText(outputs.command[format]);
    // Times in ISO 8601 form (local time), and the values of both lists start at the same x (S13 screen review L5).
    const iso = /^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$/;
    await expect(page.getByTestId('run-snapshot').locator('[data-detail="queued"]')).toHaveText(iso);
    await expect(page.getByTestId('run-record').locator('[data-detail="started"]')).toHaveText(iso);
    const valueStarts = await page
      .locator('[data-testid="run-snapshot"] > dd, [data-testid="run-record"] > dd')
      .evaluateAll((values) => [...new Set(values.map((dd) => Math.round(dd.getBoundingClientRect().left)))]);
    expect(valueStarts).toHaveLength(1);
    const badge = page.getByTestId('verification-badge');
    if (!BUILD_HAS_ENGINE) {
      await expect(badge).toHaveAttribute('data-level', 'development');
      await expect(badge).toHaveText('Development build');
    } else {
      const argv = (await text(page.getByTestId('run-argv'))).split(' ');
      const [, path, threads] = /^\s*(serial|threaded), (\d+) threads?\s*$/.exec(await text(page.getByTestId('run-record').locator('[data-detail="path"]')))!;
      const level = tableLevel(c.id, argv, path!, Number(threads));
      await expect(badge).toHaveAttribute('data-level', level);
      await expect(badge).toHaveText(level === 'certified' ? 'Certified profile' : 'Engine-supported, outside certified profile');
      test.info().annotations.push({ type: 'verification', description: `${c.id} ${argv.slice(5).join(' ') || 'defaults'}, ${path} ${threads}: ${level}` });
    }
    await expect(page.getByTestId('verification-details')).toHaveAttribute('data-level', (await badge.getAttribute('data-level'))!);
  });
}

// --- outfmt 0 and the limits -------------------------------------------------------------------

test('an HSP that outfmt 0 does not show, and a hit list that may have reached its limit', async ({ page }) => {
  // many.blastn of the outfmt 0 manifest: 260 subjects with hits, outfmt 0 shows the alignments
  // of the first 250 (BLAST+'s -num_alignments). The FakeEngine leaves its third subject out.
  const [total, shown] = BUILD_HAS_ENGINE ? [260, 250] : [3, 2];
  await program(page, 'blastn');
  await task(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'many_query.fasta', text: fasta('outfmt0/many_query.fasta') }]);
  await openFiles(page, 'subject', [{ name: 'many_subject.fasta', text: fasta('outfmt0/many_subject.fasta') }]);
  await run(page, 1);
  await openFromQueue(page, 1);
  await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', String(total));
  const partial = page.locator('[data-testid="results-notice"][data-kind="outfmt0-partial"]');
  await expect(partial).toContainText(`outfmt 0 shows the alignments of the first ${shown} subjects of this query`);
  // Fewer subjects than the default hit list (500): no limit notice.
  expect(await noticeKinds(page)).toEqual([]);
  // The headers stay over their values when the list has a scroll bar.
  expect(await misalignedHeaders(page, 'subject-table')).toEqual([]);

  // The last subject in the engine's order: the "#" column the other way round. The sort keeps
  // the selected (first) subject, at the list's end now, and shows the list's first rows.
  await show(page, 'alignment');
  const selected = await page.getByTestId('hsp-detail').getAttribute('data-hsp');
  await show(page, 'hits');
  await page.getByTestId('subject-sort-order').click();
  await expect(page.getByTestId('subject-sort-order').locator('..')).toHaveAttribute('aria-sort', 'descending');
  const subjectList = page.getByTestId('subject-list');
  await expect.poll(() => subjectList.evaluate((element) => element.scrollTop)).toBe(0);
  await subjectList.evaluate((element) => (element.scrollTop = element.scrollHeight));
  await expect(subjectRows(page).and(page.locator('[aria-pressed="true"]'))).toHaveAttribute('data-order', '1');
  await show(page, 'alignment');
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', selected!);
  await show(page, 'hits');
  await subjectList.evaluate((element) => (element.scrollTop = 0));
  const last = subjectRows(page).first();
  await expect(last).toHaveAttribute('data-order', String(total));
  await expect(last.locator('[data-field="description"]')).toHaveText('not in outfmt 0');
  const sseqid = await text(last.locator('[data-field="sseqid"]'));
  await last.click();
  await expect(last).toHaveAttribute('aria-pressed', 'true');
  await show(page, 'alignment');
  for (const row of await hspRows(page).all()) await expect(row.locator('[data-field="outfmt0"]')).toHaveText('not shown');
  await expect(page.getByTestId('detail-not-in-outfmt0')).toContainText(
    `outfmt 0 does not show this HSP. It shows the alignments of the first ${shown} subjects of this query`,
  );
  await expect(page.getByTestId('detail-heading')).toHaveCount(0);
  await expect(page.getByTestId('detail-section')).toHaveCount(0);
  await expect(page.getByTestId('range-section')).toHaveCount(0);
  const row = (await text(page.getByTestId('detail-row'))).replace(/\n$/, '');
  expect(row.split('\t')[1]).toBe(sseqid);
  // The Alignments are the tab of the next run opened, and of the results shown again (W4b).
  await page.getByTestId('tab-search').click();
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId(TABS.alignment)).toHaveAttribute('aria-pressed', 'true');
  // outfmt 6 has the HSP; outfmt 0 has no alignment heading for its subject.
  await showOutput(page, 1, 6, false);
  expect((await text(page.getByTestId('result-output'))).split('\n')).toContain(row);
  await showOutput(page, 1, 0, false);
  const escaped = sseqid.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  expect(await text(page.getByTestId('result-output'))).not.toMatch(new RegExp(`^> ?${escaped}(\\s|$)`, 'm'));
  await show(page, 'alignment');

  // An explicit -max_target_seqs that the query's subjects reach (many.mts255.blastn): the
  // notice says that more subjects may match, not that hits were lost.
  const limit = BUILD_HAS_ENGINE ? 255 : 3;
  await page.getByTestId('tab-search').click();
  await openParameters(page);
  await page.getByTestId('param-max_target_seqs').fill(String(limit));
  await run(page, 2);
  await openFromQueue(page, 2);
  await expect(page.getByTestId(TABS.alignment)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('results-run-options')).toHaveText(`-task blastn -max_target_seqs ${limit}`);
  await expect(page.getByTestId('results-program')).toHaveText('BLASTN (task blastn)');
  const notice = page.locator('[data-testid="results-notice"][data-kind="subject-limit"]');
  await expect(notice).toContainText(
    `This query has ${limit} subjects, the most that the search keeps (-max_target_seqs: ${limit}). More subjects may match; ` +
      'a search with a larger -max_target_seqs would show them.',
  );
  expect(await text(notice)).not.toMatch(/\b(lost|missing|dropped|discarded|omitted|truncated)\b/i);
  expect(await noticeKinds(page)).toEqual(['subject-limit']);
  if (BUILD_HAS_ENGINE) {
    // With -max_target_seqs, outfmt 0 shows the alignments of every subject kept.
    await expect(page.getByTestId('alignments-subject')).toBeVisible();
    await show(page, 'hits');
    await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', String(limit));
    await expect(partial).toHaveCount(0);
  }
});

test("BLASTN: an HSP of one letter has no orientation in its record; the note points to outfmt 0's Strand= line", async ({ page }) => {
  // The example of docs/evidence/losat_web_e2c/AUTHORITY.md §P: a one-letter HSP on the minus
  // strand. The FakeEngine writes a one-letter BLASTN HSP for its third query.
  const qIdx = BUILD_HAS_ENGINE ? 0 : 2;
  await program(page, 'blastn');
  await task(page, 'blastn');
  await openParameters(page);
  await page.getByTestId('param-word_size').fill('4');
  await paste(page, 'query', BUILD_HAS_ENGINE ? '>q\nTAGGACGG\n' : '>q1\nACGTACGTACGT\n>q2\nACGTACGTACGT\n>q\nTAGGACGG\n');
  await paste(page, 'subject', '>s\nYCAYAANTNCRGYACT\n');
  await run(page, 1);
  await openFromQueue(page, 1);
  // "Results for" only for the FakeEngine's three queries.
  if (BUILD_HAS_ENGINE) await expect(page.getByTestId('query-list')).toHaveCount(0);
  else await page.getByTestId(`query-row-${qIdx}`).click();
  await expect(page.getByTestId('subject-row-0')).toHaveAttribute('aria-pressed', 'true');
  await show(page, 'alignment');
  const single = page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"][data-orientation="unknown"]');
  await expect(single).toHaveCount(1);
  await single.click();
  await expect(single).toHaveAttribute('aria-pressed', 'true');
  expect(await text(single.locator('[data-field="query"]'))).toMatch(/^(\d+)–\1$/);
  expect(await text(single.locator('[data-field="subject"]'))).toMatch(/^(\d+)–\1$/);
  await expect(single.locator('[data-field="orientation"]')).toHaveText('Not in the record');
  const note = page.getByTestId('detail-strand-note');
  await expect(note).toContainText('so its coordinates do not show its strand');
  await expect(note.locator('code')).toHaveText('Strand=');
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready');
  const section = await text(page.getByTestId('detail-section'));
  expect(section).toMatch(/^ Strand=(Plus|Minus)\/(Plus|Minus)$/m);
  if (BUILD_HAS_ENGINE) expect(section).toContain(' Strand=Plus/Minus\n');
  await show(page, 'dotplot');
  await expect(page.getByTestId('dotplot')).toContainText('One letter: the HSP record does not say its strand');
  await showOutput(page, 1, 0, false);
  expect(await text(page.getByTestId('result-output'))).toContain(section);
});

// --- selection, the Alignments, filters, notices, many queries ----------------------------------------

async function panelRun(page: Page): Promise<void> {
  await program(page, 'blastn');
  await paste(page, 'query', PANEL_QUERIES);
  await paste(page, 'subject', PANEL_SUBJECTS);
  await run(page, 1);
  await openFromQueue(page, 1);
}

/** Selects the first subject of the selected query with two HSPs, and returns its index. */
async function subjectWithTwoHsps(page: Page): Promise<number> {
  await show(page, 'hits');
  for (const row of await subjectRows(page).all()) {
    if ((await text(row.locator('[data-field="hsps"]'))).trim() === '2') {
      await row.click();
      await expect(row).toHaveAttribute('aria-pressed', 'true');
      return subjectIndex(row);
    }
  }
  throw new Error('no subject with two HSPs');
}

test('selection by HSP identity: sorting keeps the selected HSP; another subject or query moves it', async ({ page }) => {
  await panelRun(page);
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'true');
  // A subject of the first query with two HSPs; choose its second HSP.
  const sIdx = await subjectWithTwoHsps(page);
  await show(page, 'alignment');
  await expect(hspRows(page)).toHaveCount(2);
  const second = hspRows(page).nth(1);
  const id = await hspId(second);
  await second.click();
  await expectSelected(page, id, sIdx);

  // The "#" column the other way round: the selected HSP is now first, and still selected.
  await show(page, 'alignment');
  await page.getByTestId('hsp-sort-rank').click();
  await expect(page.getByTestId('hsp-sort-rank').locator('..')).toHaveAttribute('aria-sort', 'descending');
  expect(await hspId(hspRows(page).first())).toBe(id);
  await expectSelected(page, id, sIdx);
  for (const sort of [
    'hsp-sort-bitScore',
    'hsp-sort-eValue',
    'hsp-sort-qStart',
    'hsp-sort-sStart',
    'subject-sort-bitScore',
    'subject-sort-eValue',
    'subject-sort-length',
    'subject-sort-hsps',
    'subject-sort-order',
  ]) {
    // Each column one way, then the other.
    for (let click = 0; click < 2; click++) {
      await show(page, sort.startsWith('hsp-') ? 'alignment' : 'hits');
      await page.getByTestId(sort).click();
      await expect(page.getByTestId(sort).locator('..')).toHaveAttribute('aria-sort', /^(ascending|descending)$/);
      await expectSelected(page, id, sIdx);
    }
  }

  // Another subject: the selection moves to its first HSP in the current order.
  const unselected = page.getByTestId('subject-list').locator('[data-testid^="subject-row-"][aria-pressed="false"]').first();
  const otherIdx = await subjectIndex(unselected);
  await page.getByTestId(`subject-row-${otherIdx}`).click();
  await expect(page.getByTestId(`subject-row-${otherIdx}`)).toHaveAttribute('aria-pressed', 'true');
  await show(page, 'alignment');
  const firstOfOther = await hspId(hspRows(page).first());
  expect(firstOfOther).not.toBe(id);
  await expectSelected(page, firstOfOther, otherIdx);

  // Another query: the selection moves to its first subject's first HSP.
  await page.getByTestId('query-row-1').click();
  await expect(page.getByTestId('query-row-1')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'false');
  const firstSubject = await subjectIndex(subjectRows(page).first());
  await show(page, 'alignment');
  const firstHsp = await hspId(hspRows(page).first());
  expect(firstHsp).toMatch(/^1:/);
  await expectSelected(page, firstHsp, firstSubject);
});

test('the Alignments: a block per Range; Next, Previous and First Match; the previous and next subject; Ranges far apart in a window', async ({
  page,
}) => {
  await panelRun(page);
  const sIdx = await subjectWithTwoHsps(page);
  const subjectOrder = await subjectRows(page).evaluateAll((rows) => rows.map((row) => Number((row as HTMLElement).dataset['testid']!.replace('subject-row-', ''))));
  await show(page, 'alignment');
  const [firstId, secondId] = [await hspId(hspRows(page).nth(0)), await hspId(hspRows(page).nth(1))];
  const [first, second] = [page.getByTestId(`range-${firstId.replace(':', '-')}`), page.getByTestId(`range-${secondId.replace(':', '-')}`)];
  await expect(rangeBlocks(page)).toHaveCount(2);
  await expect(first.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', firstId);
  // At the ends, the buttons that would leave the subject's Ranges are disabled.
  await expect(first.getByTestId('range-previous')).toBeDisabled();
  await expect(first.getByTestId('range-first')).toBeDisabled();
  await expect(second.getByTestId('range-next')).toBeDisabled();
  // Next Match: the selection, the HSP table and the focus move to the next Range.
  await first.getByTestId('range-next').click();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', secondId);
  await expect(second.getByTestId('hsp-detail')).toHaveCount(1);
  await expect(page.getByTestId(`hsp-row-${secondId.replace(':', '-')}`)).toHaveAttribute('aria-pressed', 'true');
  await expect(second.getByTestId('range-label')).toBeFocused();
  // The first Range shows its section, read when it came into view.
  await expect(first.getByTestId('range-section')).toHaveAttribute('data-state', 'ready');
  await second.getByTestId('range-previous').click();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', firstId);
  await expect(first.getByTestId('range-label')).toBeFocused();
  await first.getByTestId('range-next').click();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', secondId);
  await second.getByTestId('range-first').click();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', firstId);
  // A click in the HSP table selects the HSP and brings its Range into view.
  await page.getByTestId(`hsp-row-${secondId.replace(':', '-')}`).click();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', secondId);
  await expect(second.getByTestId('range-label')).toBeInViewport();

  // Previous and Next subject follow the Descriptions' order; "Descriptions" goes back to them.
  const at = subjectOrder.indexOf(sIdx);
  const previous = page.getByTestId('alignments-prev-subject');
  const next = page.getByTestId('alignments-next-subject');
  if (at === 0) await expect(previous).toBeDisabled();
  if (at < subjectOrder.length - 1) {
    await next.click();
    await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', String(subjectOrder[at + 1]));
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', /^0:/);
    await previous.click();
    await expect(page.getByTestId('alignments-subject')).toHaveAttribute('data-subject', String(sIdx));
  }
  await page.getByTestId('alignments-descriptions').click();
  await expect(page.getByTestId(TABS.hits)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId(`subject-row-${sIdx}`)).toHaveAttribute('aria-pressed', 'true');

  if (!BUILD_HAS_ENGINE) return;
  // A pair of many HSPs (a repeated unit searched against itself, about two HSPs a copy): the
  // selected Range and 25 on each side; more on request; sections read only where they are seen.
  await page.getByTestId('tab-search').click();
  for (const role of ['query', 'subject'] as const) {
    const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
    while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
    await openFiles(page, role, [{ name: 'repeat.fna', text: repeats(20, 60) }]);
  }
  await run(page, 2);
  await openFromQueue(page, 2);
  await expect(page.getByTestId(TABS.hits)).toHaveAttribute('aria-pressed', 'true');
  const { asked, failNext } = await watchReads(page);
  const before = await asked();
  await show(page, 'alignment');
  const total = Number(await page.getByTestId('hsp-list').getAttribute('data-count'));
  expect(total).toBeGreaterThan(80);
  await expect(rangeBlocks(page)).toHaveCount(26);
  await expect(page.getByTestId('alignments-show-earlier')).toHaveCount(0);
  // Each block that was seen asked for its section once (the reads counted as they were asked
  // for): every read asked for arrived in a block; the blocks never seen still wait.
  const sectionsIn = (state: string) => page.getByTestId('alignments').locator(`[data-testid="range-section"][data-state="${state}"]`);
  await settledSections(page, asked);
  expect((await asked()) - before).toBe(await sectionsIn('ready').count());
  // Further down: the blocks around the 13th ask for theirs. The first read asked for fails.
  await failNext(1);
  await rangeBlocks(page).nth(12).evaluate((block) => block.scrollIntoView({ block: 'start' }));
  await settledSections(page, asked);
  const readNow = await sectionsIn('ready').count();
  expect(readNow).toBeGreaterThan(0);
  expect((await asked()) - before).toBe(readNow + 1);
  expect(readNow).toBeLessThan(24);
  await expect(rangeBlocks(page).nth(1).getByTestId('range-section')).toHaveAttribute('data-state', 'pending');
  // The failed block says so; it is read again when it comes back into view.
  const failed = sectionsIn('failed');
  await expect(failed).toHaveCount(1);
  await expect(failed).toContainText('The alignment could not be read.');
  const failedRange = (await failed.getAttribute('data-range'))!;
  await page.evaluate(() => window.scrollTo(0, 0));
  await page.getByTestId(`range-${failedRange}`).evaluate((block) => block.scrollIntoView({ block: 'start' }));
  await expect(page.getByTestId(`range-${failedRange}`).getByTestId('range-section')).toHaveAttribute('data-state', 'ready');
  await expect(page.getByTestId(`range-${failedRange}`).getByTestId('range-section')).toHaveText(/^ Score = /);
  // The buttons are pressed without the mouse: sections read as they come into view move the blocks.
  await page.getByTestId('alignments-show-later').dispatchEvent('click');
  await expect(rangeBlocks(page)).toHaveCount(Math.min(total, 51));
  // Failed reads, read again with "Try again". Every read fails for a while, so that no section
  // arrives to move the blocks (a failed block that leaves the screen is read again on its return).
  await failNext(1000);
  await rangeBlocks(page).nth(45).evaluate((block) => block.scrollIntoView({ block: 'start' }));
  await settledSections(page, asked);
  await failNext(0);
  const retried = rangeBlocks(page).nth(45).getByTestId('range-section');
  await expect(retried).toHaveAttribute('data-state', 'failed');
  await expect(retried).toBeInViewport();
  await retried.getByTestId('range-retry').click();
  await expect(retried).toHaveAttribute('data-state', 'ready');
  await expect(retried).toHaveText(/^ Score = /);
  // The last Range, chosen in the HSP table: the window moves to it.
  await page.getByTestId('hsp-sort-rank').click();
  await expect(page.getByTestId('hsp-sort-rank').locator('..')).toHaveAttribute('aria-sort', 'descending');
  await hspRows(page).first().click();
  const lastId = await hspId(hspRows(page).first());
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', lastId);
  await expect(rangeBlocks(page).last()).toHaveAttribute('data-testid', `range-${lastId.replace(':', '-')}`);
  await expect(rangeBlocks(page).last()).toHaveAttribute('data-n', String(total));
  await expect(rangeBlocks(page)).toHaveCount(26);
  await expect(page.getByTestId('alignments-show-later')).toHaveCount(0);
  // Earlier Ranges are added above the blocks in view, and the screen stays where it was: it shows
  // the Ranges just before the one that was first (whose sections, read as they come into view,
  // grow below their labels), not the first Ranges added.
  const earlier = page.getByTestId('alignments-show-earlier');
  await earlier.scrollIntoViewIfNeeded();
  const firstShown = Number(await rangeBlocks(page).first().getAttribute('data-n'));
  await earlier.dispatchEvent('click');
  await expect(rangeBlocks(page)).toHaveCount(51);
  await expect(rangeBlocks(page).nth(25)).toHaveAttribute('data-n', String(firstShown));
  await page.waitForTimeout(500);
  await expect(rangeBlocks(page).first()).not.toBeInViewport();
  const nearInView = () =>
    rangeBlocks(page).evaluateAll((blocks) =>
      blocks.slice(20, 26).filter((block) => {
        const box = block.getBoundingClientRect();
        return box.bottom > 0 && box.top < window.innerHeight;
      }).length,
    );
  expect(await nearInView()).toBeGreaterThan(0);
});

test('many queries; view filters change the view, not the search; the notices tell filtered out, no hits and a failed run apart', async ({
  page,
}) => {
  await panelRun(page);
  const queue = page.getByTestId('queue').locator(':scope > li');
  await expect(queue).toHaveCount(1);
  // The queue's card has one layout, at the same place beside either tab (W4 screen review low 5).
  const card = async () => {
    const [run, open] = [(await page.getByTestId('run-1').boundingBox())!, (await page.getByTestId('run-1-open').boundingBox())!];
    return [run.x, run.width, run.height, open.x - run.x, open.y - run.y].map(Math.round);
  };
  const inResults = await card();
  await page.getByTestId('tab-search').click();
  expect(await card()).toEqual(inResults);
  await page.getByTestId('tab-results').click();
  // The selection is read in the Alignments; "Results for", "Filter Results" and the notices are above the tabs.
  await show(page, 'alignment');
  const detail = page.getByTestId('hsp-detail');

  // A query without hits: NCBI's words; the tabs stay.
  await page.getByTestId('query-row-3').click();
  await expect(page.getByTestId('query-row-3')).toContainText('no hits');
  await expect(page.locator('[data-testid="results-notice"][data-kind="no-hits"]')).toHaveText('No significant similarity found for this query.');
  expect(await noticeKinds(page)).toEqual(['no-hits']);
  await expect(page.getByTestId('alignments')).toHaveCount(0);
  await expect(page.getByTestId(TABS.alignment)).toBeVisible();
  await show(page, 'hits');
  await expect(page.getByTestId('subject-table')).toHaveCount(0);
  await show(page, 'graphic');
  await expect(page.getByTestId('graphic-summary')).toHaveCount(0);
  await show(page, 'alignment');

  // The query picker draws only the rows in view, and finds a query by its ID.
  const list = page.getByTestId('query-list');
  await expect(list).toHaveAttribute('data-count', '150');
  await expect(page.getByTestId('query-count')).toHaveText('150 of 150 queries');
  const drawn = await page.locator('[data-testid^="query-row-"]').count();
  expect(drawn).toBeGreaterThan(5);
  expect(drawn).toBeLessThan(40);
  await list.evaluate((element) => (element.scrollTop = element.scrollHeight));
  await expect(page.getByTestId('query-row-149')).toBeVisible();
  await expect(page.getByTestId('query-row-0')).toHaveCount(0);
  // The selected query (#4, no hits) leaves the list: the first query listed is selected.
  await page.getByTestId('filter-hits-only').check();
  await expect(page.getByTestId('query-count')).toHaveText('113 of 150 queries');
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'true');
  await expect(detail).toHaveAttribute('data-hsp', /^0:/);
  expect(await noticeKinds(page)).not.toContain('no-hits');
  await page.getByTestId('filter-hits-only').uncheck();
  await page.getByTestId('query-filter').fill('rec137');
  await expect(page.getByTestId('query-count')).toHaveText('1 of 150 queries');
  await expect(page.locator('[data-testid^="query-row-"]')).toHaveCount(1);
  await page.getByTestId('query-row-136').click();
  await expect(page.getByTestId('query-row-136')).toHaveAttribute('aria-pressed', 'true');
  await expect(detail).toHaveAttribute('data-hsp', /^136:/);
  await page.getByTestId('query-filter').fill('');
  await expect(page.getByTestId('query-count')).toHaveText('150 of 150 queries');
  await list.evaluate((element) => (element.scrollTop = 0));
  await page.getByTestId('query-row-0').click();
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'true');

  // A subject filter that matches nothing: every HSP of the query is hidden, and one click clears it.
  await show(page, 'hits');
  const [, subjects, hsps] = /([\d,]+) subj\., ([\d,]+) HSPs/.exec(await text(page.getByTestId('query-row-0')))!;
  await page.getByTestId('filter-subject').fill('no-such-subject');
  await page.getByTestId('filter-apply').click();
  const filteredOut = page.locator('[data-testid="results-notice"][data-kind="filtered-out"]');
  await expect(filteredOut).toContainText(
    `No HSPs of this query match the view filters (${plural(Number(hsps), 'HSP')} on ${plural(Number(subjects), 'subject')} hidden).`,
  );
  expect(await noticeKinds(page)).toEqual(['filtered-out']);
  await expect(page.getByTestId('subject-table')).toHaveCount(0);
  await show(page, 'alignment');
  await expect(page.getByTestId('hsp-table')).toHaveCount(0);
  await show(page, 'hits');
  await page.getByTestId('results-notice-clear').click();
  await expect(filteredOut).toHaveCount(0);
  await expect(page.getByTestId('filter-subject')).toHaveValue('');
  await expect(subjectRows(page)).toHaveCount(Number(subjects));

  // A filter that hides some HSPs: s1 only.
  await page.getByTestId('filter-subject').fill('s1');
  await page.getByTestId('filter-subject').press('Enter');
  await expect(subjectRows(page)).toHaveCount(1);
  await expect(subjectRows(page).locator('[data-field="sseqid"]')).toHaveText('s1');
  const kept = Number((await text(subjectRows(page).locator('[data-field="hsps"]'))).trim());
  await expect(page.locator('[data-testid="results-notice"][data-kind="filtered-some"]')).toContainText(
    `The view filters hide ${plural(Number(hsps) - kept, 'HSP')} and ${plural(Number(subjects) - 1, 'subject')} of this query.`,
  );
  // The boxes show the filters in force when the results tab is shown again (the form is drawn
  // again), and applying one box keeps the others.
  await page.getByTestId('filter-bits').fill('5000');
  await page.getByTestId('filter-apply').click();
  await expect(filteredOut).toHaveCount(1);
  await page.getByTestId('tab-search').click();
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('filter-bits')).toHaveValue('5000');
  await expect(page.getByTestId('filter-subject')).toHaveValue('s1');
  await page.getByTestId('filter-subject').press('Enter');
  await expect(filteredOut).toHaveCount(1);
  await page.getByTestId('filter-bits').fill('');
  await page.getByTestId('filter-apply').click();
  await expect(subjectRows(page)).toHaveCount(1);
  await expect(filteredOut).toHaveCount(0);
  // A value that is not a number is not applied.
  await page.getByTestId('filter-evalue').fill('abc');
  await page.getByTestId('filter-apply').click();
  await expect(page.getByTestId('filter-problem')).toHaveText('E value must be a number, such as 1e-5 or 50.');
  await page.getByTestId('filter-clear').click();
  await expect(page.getByTestId('filter-problem')).toHaveCount(0);
  await expect(subjectRows(page)).toHaveCount(Number(subjects));
  expect(await noticeKinds(page)).toEqual([]);
  // None of this was a search.
  await expect(queue).toHaveCount(1);
  await expect(page.getByTestId('results-run').locator('option:not([disabled])')).toHaveCount(1);
  await expect(page.getByTestId('run-1-status')).toHaveText('completed');

  if (BUILD_HAS_ENGINE) {
    // A run that fails: TBLASTX refuses at run time a subject title with an HTML character
    // reference, which NCBI decodes in outfmt 0 (docs/web/abi_v2.md §4; screens.spec.ts).
    const sequence = dna(7, 300);
    // A view filter of run 1 is not carried to the next run that opens.
    await page.getByTestId('filter-subject').fill('s1');
    await page.getByTestId('filter-subject').press('Enter');
    await expect(subjectRows(page)).toHaveCount(1);
    await page.getByTestId('tab-search').click();
    await program(page, 'tblastx');
    await paste(page, 'query', `>q1\n${sequence}\n`);
    await paste(page, 'subject', `>s1 alpha &amp; beta\n${sequence}\n`);
    await submit(page);
    await waitStatus(page, 2, 'failed');
    await expect(page.getByTestId('run-2-open')).toHaveCount(0);
    await page.getByTestId('tab-results').click();
    const select = page.getByTestId('results-run');
    await select.selectOption({ label: 'Run 2 · TBLASTX · query.fa vs subject.fa (failed)' });
    const status = page.getByTestId('results-status');
    await expect(status).toHaveAttribute('data-run-status', 'failed');
    await expect(status).toHaveAttribute('data-phase', 'unavailable');
    await expect(status).toContainText('Run 2 failed: ');
    await expect(status).toContainText('A failed run keeps no results.');
    // The header block still says which run this is.
    await expect(page.getByTestId('results-program')).toHaveText('TBLASTX');
    await expect(page.getByTestId('results-download-all')).toHaveCount(0);
    await expect(page.getByTestId('results-hits')).toHaveCount(0);
    await expect(page.getByTestId('results-notice')).toHaveCount(0);
    const run1 = await select.locator('option', { hasText: /^Run 1 · / }).getAttribute('value');
    await select.selectOption(run1!);
    await expect(page.getByTestId('filter-subject')).toHaveValue('');
    await expect(subjectRows(page)).toHaveCount(Number(subjects));
    expect(await noticeKinds(page)).toEqual([]);
  }
});

// --- the dot plot ------------------------------------------------------------------------------

/**
 * Counts the canvas's pixels of an HSP colour of blast2dotplot.py (#1f77b4 forward, #ff7f0e
 * reverse) at any of its opacities over the white plot (1, 0.8, 0.6, 0.4).
 */
function colouredPixels(canvas: Locator): Promise<number> {
  return canvas.evaluate((element) => {
    const c = element as HTMLCanvasElement;
    const { data } = c.getContext('2d')!.getImageData(0, 0, c.width, c.height);
    const colours = [
      [31, 119, 180],
      [255, 127, 14],
    ];
    const tints = colours.flatMap((rgb) => [1, 0.8, 0.6, 0.4].map((alpha) => rgb.map((v) => Math.round(255 - alpha * (255 - v)))));
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      const [r, g, b] = [data[i]!, data[i + 1]!, data[i + 2]!];
      if (tints.some(([tr, tg, tb]) => Math.abs(r - tr!) + Math.abs(g - tg!) + Math.abs(b - tb!) < 30)) n++;
    }
    return n;
  });
}

test('the dot plot: the HSPs of the pair on a canvas; zoom; choosing an HSP on it or in the list; its popup', async ({ page }) => {
  // q2 is A, which the subject holds twice: two HSPs far apart (the FakeEngine also writes
  // two HSPs for its second query).
  await program(page, 'blastn');
  await paste(page, 'query', `>q1\n${dna(21, 20)}\n>q2\n${A}\n`);
  await paste(page, 'subject', `>s1\n${A}${dna(12, 60)}${A}\n`);
  await run(page, 1);
  await openFromQueue(page, 1);
  const rows = outfmt6Rows((await readOutputs(page, 1)).text[6]);
  await page.getByTestId('query-row-1').click();
  // The Dot Plot tab has the HSP table under the plot.
  await show(page, 'dotplot');
  await expect(hspRows(page)).toHaveCount(2);
  const [firstId, secondId] = [await hspId(hspRows(page).nth(0)), await hspId(hspRows(page).nth(1))];
  const canvas = page.getByTestId('dotplot-canvas');
  await expect(canvas).toHaveAttribute('data-segments', '2');
  await expect(canvas).toHaveAttribute('data-selected', firstId);
  const full = '0,60,0,180';
  await expect(canvas).toHaveAttribute('data-view', full);

  // The canvas is drawn: pixels of the HSPs' colours.
  expect(await colouredPixels(canvas)).toBeGreaterThan(20);

  // Zoom in, in again, out, and back to the whole sequences.
  const span = async () => {
    const [x0, x1, y0, y1] = (await canvas.getAttribute('data-view'))!.split(',').map(Number);
    return { x: x1! - x0!, y: y1! - y0! };
  };
  await page.getByTestId('dotplot-zoom-in').click();
  await expect(canvas).not.toHaveAttribute('data-view', full);
  const once = await span();
  expect(once.x).toBeLessThan(60);
  // Both axes keep the same scale.
  expect(once.y / once.x).toBeCloseTo(3, 5);
  await page.getByTestId('dotplot-zoom-in').click();
  await expect.poll(async () => (await span()).x).toBeLessThan(once.x);
  const twice = await span();
  await page.getByTestId('dotplot-zoom-out').click();
  await expect.poll(async () => (await span()).x).toBeGreaterThan(twice.x);
  await page.getByTestId('dotplot-reset').click();
  await expect(canvas).toHaveAttribute('data-view', full);

  // The mouse wheel zooms only with Ctrl or ⌘ held; otherwise the page scrolls, as it does on a
  // phone over the plot (S13 screen review L4).
  const wheel = (keys: { ctrlKey?: boolean; metaKey?: boolean }) =>
    canvas.evaluate((element, held: { ctrlKey?: boolean; metaKey?: boolean }) => {
      const event = new WheelEvent('wheel', { deltaY: -100, bubbles: true, cancelable: true, ...held });
      element.dispatchEvent(event);
      return event.defaultPrevented;
    }, keys);
  expect(await wheel({})).toBe(false);
  await expect(canvas).toHaveAttribute('data-view', full);
  for (const modifier of [{ ctrlKey: true }, { metaKey: true }]) {
    expect(await wheel(modifier)).toBe(true);
    await expect(canvas).not.toHaveAttribute('data-view', full);
    await page.getByTestId('dotplot-reset').click();
    await expect(canvas).toHaveAttribute('data-view', full);
  }
  // Screen readers: the plot is an application (they pass its keys to it, as for the Graphic
  // Summary), whose name gives the keys.
  await expect(
    page.getByRole('application', {
      name: /^Dot plot of 2 HSPs\. .*Ctrl or ⌘ and the mouse wheel.* n and p to select the next or previous HSP; Enter to show the selected HSP, and Escape to close it\.$/,
    }),
  ).toHaveAttribute('data-testid', 'dotplot-canvas');
  await expect(canvas).toHaveCSS('touch-action', 'pan-y');
  await expect(page.getByTestId('dotplot-legend-selected')).toBeVisible();

  // An HSP chosen in the list is selected on the plot; "Zoom to HSP" frames it, from a view
  // zoomed in away from it (with the same scale on both axes, the frame of an HSP that spans the
  // whole query is the whole plot).
  await hspRows(page).nth(1).click();
  await expect(canvas).toHaveAttribute('data-selected', secondId);
  const [qStart, qEnd] = (await text(hspRows(page).nth(1).locator('[data-field="query"]'))).split('–').map(Number);
  const [sStart, sEnd] = (await text(hspRows(page).nth(1).locator('[data-field="subject"]'))).split('–').map(Number);
  await page.getByTestId('dotplot-zoom-in').click();
  await page.getByTestId('dotplot-zoom-in').click();
  await expect.poll(async () => (await span()).x).toBeLessThan(30);
  await page.getByTestId('dotplot-zoom-hsp').click();
  const frames = async () => {
    const [x0, x1, y0, y1] = (await canvas.getAttribute('data-view'))!.split(',').map(Number);
    return (
      x0! <= Math.min(qStart!, qEnd!) && x1! >= Math.max(qStart!, qEnd!) && y0! <= Math.min(sStart!, sEnd!) && y1! >= Math.max(sStart!, sEnd!)
    );
  };
  await expect.poll(frames).toBe(true);
  await page.getByTestId('dotplot-reset').click();

  // A click on the other HSP's line selects it, in the plot and in the list, and opens its popup:
  // the values of its outfmt 6 row as written.
  const targets = JSON.parse((await canvas.getAttribute('data-targets'))!) as { hsp: string; x: number; y: number }[];
  const target = targets.find((t) => t.hsp === firstId)!;
  await canvas.click({ position: { x: target.x, y: target.y } });
  await expect(canvas).toHaveAttribute('data-selected', firstId);
  await expect(hspRows(page).nth(0)).toHaveAttribute('aria-pressed', 'true');
  await expect(hspRows(page).nth(1)).toHaveAttribute('aria-pressed', 'false');
  const firstQuery = await text(hspRows(page).nth(0).locator('[data-field="query"]'));
  await expect(page.getByTestId('dotplot-selected')).toContainText(`query ${firstQuery} nt`);
  const popup = page.getByTestId('dotplot-popup');
  await expect(popup).toHaveAttribute('data-hsp', firstId);
  // Beside the line, off one of its ends: it does not cover the line's midpoint (W4b screen review L10).
  await expect(popup).toHaveAttribute('data-place', 'beside');
  const covers = async (point: { x: number; y: number }) => {
    const [at, plot] = [(await popup.boundingBox())!, (await canvas.boundingBox())!];
    const [x, y] = [plot.x + point.x, plot.y + point.y];
    return x >= at.x && x <= at.x + at.width && y >= at.y && y <= at.y + at.height;
  };
  expect(await covers(target)).toBe(false);
  const fields = rows.filter((row) => row[0] === 'q2')[Number(firstId.split(':')[1])]!;
  const value = (name: string) => popup.locator(`dd[data-field="${name}"]`);
  await expect(value('bitscore')).toHaveText(fields[11]!);
  await expect(value('evalue')).toHaveText(fields[10]!);
  await expect(value('pident')).toHaveText(fields[2]!);
  await expect(value('query')).toHaveText(`${fields[6]}–${fields[7]} nt`);
  await expect(value('subject')).toHaveText(`${fields[8]}–${fields[9]} nt`);
  await expect(value('outfmt0')).toHaveText(BUILD_HAS_ENGINE ? 'shown' : 'shown');
  // Escape closes it, and the focus is on the plot.
  await page.keyboard.press('Escape');
  await expect(popup).toHaveCount(0);
  await expect(canvas).toBeFocused();
  // The keyboard: n selects the next HSP; Enter opens its popup, with the focus in it; Escape
  // closes it and gives the focus back to the plot.
  await page.keyboard.press('n');
  await expect(canvas).toHaveAttribute('data-selected', secondId);
  await expect(hspRows(page).nth(1)).toHaveAttribute('aria-pressed', 'true');
  await page.keyboard.press('Enter');
  await expect(popup).toHaveAttribute('data-hsp', secondId);
  await expect(popup).toBeFocused();
  expect(await covers(targets.find((t) => t.hsp === secondId)!)).toBe(false);
  await expect(page.getByRole('dialog', { name: `HSP ${Number(secondId.split(':')[1]) + 1}` })).toHaveAttribute('data-testid', 'dotplot-popup');
  await page.keyboard.press('Escape');
  await expect(popup).toHaveCount(0);
  await expect(canvas).toBeFocused();
  // Reached with Tab, the plot opens the popup with Enter, and the popup shows its focus in every
  // browser (W4b screen review L11: Firefox showed none).
  await page.getByTestId('dotplot-reset').focus();
  await page.keyboard.press('Tab');
  await expect(canvas).toBeFocused();
  await page.keyboard.press('Enter');
  await expect(popup).toBeFocused();
  await expect(popup).toHaveCSS('outline-style', 'solid');
  await expect(popup).toHaveCSS('outline-width', '2px');
  // "Show alignment": the Alignments, with the HSP's Range in view and the focus on its label.
  await popup.getByTestId('dotplot-popup-alignment').click();
  await expect(page.getByTestId(TABS.alignment)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', secondId);
  const range = page.getByTestId(`range-${secondId.replace(':', '-')}`);
  await expect(range.getByTestId('range-label')).toBeFocused();
  await expect(range).toBeInViewport();
});

// --- the phone size (S13 screen review) ------------------------------------------------------------

test('narrow screens: the results are no wider than the screen; "Open results" shows them; the tables show the key values first', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  // A translated search (frames) of files with long names: the run's label in the run picker is long.
  const c = PROGRAM_CASES.find((programCase) => programCase.id === 'tblastx')!;
  await program(page, c.id);
  await openFiles(page, 'query', [{ name: basename(c.query), text: fasta(c.query) }]);
  await openFiles(page, 'subject', [{ name: basename(c.subject), text: fasta(c.subject) }]);
  await run(page, 1);

  // The queue is under the search form: "Open results" brings the results' heading into view and focuses it.
  await page.getByTestId('run-1-open').click();
  const heading = page.getByTestId('results-heading');
  await expect(heading).toBeFocused();
  await expect(heading).toBeInViewport();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  await expect(page.getByTestId('subject-table')).toBeVisible();
  await expectNoSideScroll(page, 'descriptions');

  // The key values first; the tables say that they scroll sideways.
  expect((await columnsInView(page, 'subject-table')).slice(0, 4)).toEqual(['order', 'sseqid', 'bitscore', 'evalue']);
  // The Subject column is at least 12ch wide, so that IDs with a common prefix stay apart.
  const subjectWidth = await page
    .getByTestId('subject-table')
    .locator('.table-row')
    .first()
    .evaluate((row) => {
      const probe = document.createElement('span');
      probe.style.cssText = 'position:absolute;visibility:hidden;width:12ch';
      row.append(probe);
      const twelve = probe.getBoundingClientRect().width;
      probe.remove();
      return row.querySelector<HTMLElement>('[data-field="sseqid"]')!.getBoundingClientRect().width - twelve;
    });
  expect(subjectWidth).toBeGreaterThanOrEqual(-0.5);
  /**
   * The table's sideways scroll after a click on a row (W4 screen review middle 2: Firefox scrolled
   * it). The mouse presses the row where it is in view: Playwright's click() would first scroll the
   * whole row's button into view itself.
   */
  const clickKeepsScroll = async (table: string, row: Locator) => {
    const scroll = page.getByTestId(table).locator('.table-scroll');
    await scroll.evaluate((element) => (element.scrollLeft = 0));
    await row.scrollIntoViewIfNeeded();
    await scroll.evaluate((element) => (element.scrollLeft = 0));
    const box = (await row.boundingBox())!;
    await page.mouse.click(box.x + 24, box.y + box.height / 2);
    await expect(row).toHaveAttribute('aria-pressed', 'true');
    await page.waitForTimeout(100);
    expect(await scroll.evaluate((element) => element.scrollLeft), `${table}: scrolled sideways by a click`).toBe(0);
  };
  await clickKeepsScroll('subject-table', subjectRows(page).first());

  await show(page, 'alignment');
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready');
  await expectNoSideScroll(page, 'alignments');
  // The HSP table: the frames before the ranges, so that they are in view (W4 screen review middle 2).
  expect((await columnsInView(page, 'hsp-table')).slice(0, 4)).toEqual(['rank', 'bitscore', 'evalue', 'frames']);
  await clickKeepsScroll('hsp-table', hspRows(page).nth(Math.min(1, (await hspRows(page).count()) - 1)));
  for (const table of ['subject-table', 'hsp-table']) {
    await show(page, table === 'subject-table' ? 'hits' : 'alignment');
    await expect(page.getByTestId(`${table}-scroll-hint`)).toBeVisible();
    // Each value is under its header.
    const offsets = await page
      .getByTestId(table)
      .locator('.table-scroll')
      .evaluate((scroll) => {
        const row = scroll.querySelector('.table-row')!;
        return [...scroll.querySelectorAll<HTMLElement>('.table-head [data-col]')].map((head) => {
          const cell = row.querySelector<HTMLElement>(`[data-field="${head.dataset['col']}"]`);
          return `${head.dataset['col']} ${cell === null ? 'missing' : Math.round(cell.getBoundingClientRect().left - head.getBoundingClientRect().left)}`;
        });
      });
    expect(offsets.filter((offset) => !offset.endsWith(' 0'))).toEqual([]);
    // Sort headers are touch targets of at least 24 px (WCAG 2.5.8).
    for (const button of await page.locator('.sort-button').all()) expect((await button.boundingBox())!.height).toBeGreaterThanOrEqual(24);
  }

  // The Graphic Summary, the dot plot, the run details and the outputs.
  await show(page, 'graphic');
  await expect(page.getByTestId('graphic-canvas')).toHaveAttribute('data-rows', /^[1-9]/);
  await expectNoSideScroll(page, 'graphic summary');
  await show(page, 'dotplot');
  await expect(page.getByTestId('dotplot-canvas')).toHaveAttribute('data-segments', /^[1-9]/);
  await expectNoSideScroll(page, 'dot plot');
  // The HSP popup is under the plot on a phone, not over it (W4b screen review L10).
  await page.getByTestId('dotplot-canvas').focus();
  await page.keyboard.press('Enter');
  const popup = page.getByTestId('dotplot-popup');
  await expect(popup).toHaveAttribute('data-place', 'below');
  await expect(popup).toBeFocused();
  expect((await popup.boundingBox())!.y).toBeGreaterThanOrEqual((await page.getByTestId('dotplot-canvas').boundingBox())!.y + (await page.getByTestId('dotplot-canvas').boundingBox())!.height);
  await expectNoSideScroll(page, 'dot plot popup');
  await page.keyboard.press('Escape');
  await expect(popup).toHaveCount(0);
  await show(page, 'details');
  await expect(page.getByTestId('run-details')).toBeVisible();
  await expectNoSideScroll(page, 'run details');
  await showOutput(page, 1, 0, false);
  await expectNoSideScroll(page, 'outputs');

  // "Open results" again from the search tab, the page scrolled to the queue at its end, while
  // the results are read: the heading stays in view.
  await page.getByTestId('tab-search').click();
  await page.evaluate(() => window.scrollTo(0, document.documentElement.scrollHeight));
  await page.getByTestId('run-1-open').click();
  await expect(heading).toBeFocused();
  await expect(page.getByTestId('results-outputs')).toBeVisible();
  await expect(heading).toBeInViewport({ ratio: 1 });
});
