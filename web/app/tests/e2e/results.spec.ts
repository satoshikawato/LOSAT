// The results screen (S13, W4; plan §5.7, design §11): Run -> Query -> Subject -> HSP through
// the real application. The tests run with the FakeEngine build and with the engine build
// (LOSAT_WEB_REACTORS); what only the engine can show (outfmt 0's frames, a run that the
// engine refuses, the verification badge of a real run) is checked in the engine build. The lists, the HSP detail and the dot plot are compared with the stored outputs as
// the Outputs view shows them (the outfmt 6 rows and the outfmt 0 text), never with values
// computed here.
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
import { fasta, openFiles, paste, program, showOutput, submit, waitStatus } from './support/search';

type ProgramId = 'blastn' | 'blastp' | 'tblastn' | 'tblastx';
type Unit = 'nt' | 'aa';
type Orientation = 'forward' | 'reverse' | 'unknown';

const FORMATS = [0, 6, 7] as const satisfies readonly OutputFormat[];
const LABELS: Readonly<Record<ProgramId, string>> = { blastn: 'BLASTN', blastp: 'BLASTP', tblastn: 'TBLASTN', tblastx: 'TBLASTX' };

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

async function text(locator: Locator): Promise<string> {
  return (await locator.textContent()) ?? '';
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

/** The commands and stored outputs of the run on the results screen, read in its Outputs view. */
async function readOutputs(page: Page, number: number): Promise<{ command: Record<OutputFormat, string>; text: Record<OutputFormat, string> }> {
  const command = {} as Record<OutputFormat, string>;
  const output = {} as Record<OutputFormat, string>;
  for (const format of FORMATS) {
    await showOutput(page, number, format, false);
    command[format] = await text(page.getByTestId('result-command'));
    output[format] = await text(page.getByTestId('result-output'));
  }
  await page.getByTestId('results-view-hits').click();
  return { command, text: output };
}

const subjectRows = (page: Page) => page.getByTestId('subject-list').locator('[data-testid^="subject-row-"]');
const hspRows = (page: Page) => page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"]');

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

/** Checks that the HSP `id` of subject `sIdx` is the one selection of the lists and the detail. */
async function expectSelected(page: Page, id: string, sIdx: number): Promise<void> {
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', id);
  await expect(page.getByTestId(`hsp-row-${id.replace(':', '-')}`)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('hsp-list').locator('[aria-pressed="true"]')).toHaveCount(1);
  await expect(page.getByTestId(`subject-row-${sIdx}`)).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('subject-list').locator('[aria-pressed="true"]')).toHaveCount(1);
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
  test(`${LABELS[c.id]}: the lists are the outfmt 6 rows, the detail is outfmt 0's text; units, frames and run details`, async ({ page }) => {
    const translated = c.frameLine !== undefined;
    await program(page, c.id);
    await openFiles(page, 'query', [{ name: basename(c.query), text: fasta(c.query) }]);
    await openFiles(page, 'subject', [{ name: basename(c.subject), text: fasta(c.subject) }]);
    await run(page, 1);
    await openFromQueue(page, 1);
    await expect(page.getByTestId('results-run').locator('option:checked')).toHaveText(
      `Run 1 · ${LABELS[c.id]} · ${basename(c.query)} vs ${basename(c.subject)}`,
    );
    const outputs = await readOutputs(page, 1);
    const rows = outfmt6Rows(outputs.text[6]);
    const out0 = outputs.text[0];

    // Units of the lists.
    await expect(page.getByTestId('subject-sort-length')).toHaveText(`Length (${c.units.subject})`);
    await expect(page.getByTestId('hsp-sort-qStart')).toHaveText(`Query (${c.units.query})`);
    await expect(page.getByTestId('hsp-sort-sStart')).toHaveText(`Subject (${c.units.subject})`);

    // The first query with hits is selected when the run opens; the first three are checked.
    const detail = page.getByTestId('hsp-detail');
    await expect(detail).toHaveAttribute('data-state', 'ready');
    const selectedQuery = Number((await detail.getAttribute('data-hsp'))!.split(':')[0]);
    await expect(page.getByTestId(`query-row-${selectedQuery}`)).toHaveAttribute('aria-pressed', 'true');
    const withHits: number[] = [];
    for (const row of await page.getByTestId('query-list').locator('[data-testid^="query-row-"]').all()) {
      if (!(await text(row)).includes('no hits')) withHits.push(Number((await row.getAttribute('data-testid'))!.replace('query-row-', '')));
    }
    expect(withHits[0]).toBe(selectedQuery);
    const seen = { queries: 0, subjects: 0, hsps: 0, orientations: new Set<string>(), frames: new Set<string>() };
    for (const qIdx of withHits.slice(0, 3)) {
      const queryRow = page.getByTestId(`query-row-${qIdx}`);
      await queryRow.click();
      await expect(queryRow).toHaveAttribute('aria-pressed', 'true');
      await expect(detail).toHaveAttribute('data-hsp', new RegExp(`^${qIdx}:`));
      // The query's HSPs are its outfmt 6 rows, by rank.
      const qseqid = (await text(page.getByTestId('detail-row'))).split('\t')[0]!;
      const queryRows = rows.filter((fields) => fields[0] === qseqid);
      const sseqids = [...new Set(queryRows.map((fields) => fields[1]!))];
      expect(queryRows.length).toBeGreaterThan(0);
      await expect(queryRow).toContainText(`${c.units.query} · ${count(sseqids.length)} subj., ${count(queryRows.length)} HSPs`);

      // The subject list: each subject's first outfmt 6 row, in the engine's order.
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

      // The HSPs of the first subjects: every field of the HSP's outfmt 6 row; the detail is outfmt 0's text.
      for (let i = 0; i < Math.min(drawnSubjects, 3); i++) {
        const subject = subjectRows(page).nth(i);
        await subject.click();
        await expect(subject).toHaveAttribute('aria-pressed', 'true');
        const sseqid = sseqids[i]!;
        const pairRanks = queryRows.flatMap((fields, rank) => (fields[1] === sseqid ? [rank] : []));
        await expect(page.getByTestId('hsp-list')).toHaveAttribute('data-count', String(pairRanks.length));
        await expect(hspRows(page)).toHaveCount(Math.min(pairRanks.length, 20));
        const description = subject.locator('[data-field="description"]');
        await expect(description).not.toHaveText('…');
        const inOutfmt0 = (await text(description)) !== 'not in outfmt 0';
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

        // The detail of the subject's first HSP (selected with the subject).
        const first = hspRows(page).first();
        await expect(first).toHaveAttribute('aria-pressed', 'true');
        const id = await hspId(first);
        await expect(detail).toHaveAttribute('data-hsp', id);
        await expect(detail).toHaveAttribute('data-state', 'ready');
        const fields = queryRows[Number(id.split(':')[1])]!;
        expect((await text(page.getByTestId('detail-row'))).replace(/\n$/, '')).toBe(fields.join('\t'));
        seen.subjects++;
        if (!inOutfmt0) {
          await expect(page.getByTestId('detail-not-in-outfmt0')).toBeVisible();
          continue;
        }
        const heading = await text(page.getByTestId('detail-heading'));
        const section = await text(page.getByTestId('detail-section'));
        expect(heading.length).toBeGreaterThan(0);
        expect(section.length).toBeGreaterThan(0);
        expect(out0).toContain(heading);
        expect(out0).toContain(section);
        // The description is the heading's title; the length is the heading's Length=.
        const length = (await text(subject.locator('[data-field="length"]'))).replace(/,/g, '');
        expect(heading).toContain(`${await text(description)}\nLength=${length}\n`);
        const orientation = await first.getAttribute('data-orientation');
        if (c.id === 'blastn') expect(section).toContain(` Strand=Plus/${orientation === 'reverse' ? 'Minus' : 'Plus'}\n`);
        if (translated && BUILD_HAS_ENGINE) {
          const [q, s] = (await text(first.locator('[data-field="frames"]'))).split('/');
          expect(section).toContain(c.frameLine!(q!, s!));
        }
      }
      seen.queries++;
    }
    test.info().annotations.push({
      type: 'coverage',
      description:
        `${c.id}: ${seen.queries} queries, ${seen.subjects} subjects opened, ${seen.hsps} HSP rows; ` +
        `orientations ${[...seen.orientations].sort().join(' ')}${translated ? `; frames ${[...seen.frames].sort().join(' ')}` : ''}`,
    });

    // The dot plot of the selected pair, in the units of the records.
    const fields = (await text(page.getByTestId('detail-row'))).replace(/\n$/, '').split('\t');
    const frames = translated ? await text(hspRows(page).first().locator('[data-field="frames"]')) : '';
    await page.getByTestId('pane-dotplot').click();
    const canvas = page.getByTestId('dotplot-canvas');
    await expect(canvas).toHaveAttribute('data-segments', (await page.getByTestId('hsp-list').getAttribute('data-count'))!);
    await expect(page.getByTestId('dotplot')).toContainText(new RegExp(`\\([\\d,]+ ${c.units.query}\\) against subject .* \\([\\d,]+ ${c.units.subject}\\)\\.`));
    await expect(page.getByTestId('dotplot-selected')).toContainText(
      `query ${fields[6]}–${fields[7]} ${c.units.query}, subject ${fields[8]}–${fields[9]} ${c.units.subject}` +
        (translated ? `, frames ${frames.replace('/', ' / ')}.` : '.'),
    );
    await page.getByTestId('pane-alignment').click();

    // Run details: the command of each format is the Outputs view's, and the verification badge.
    await page.getByTestId('results-view-details').click();
    for (const format of FORMATS) await expect(page.getByTestId(`run-command-${format}`)).toHaveText(outputs.command[format]);
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
  await page.getByTestId('param-task').selectOption('blastn');
  await openFiles(page, 'query', [{ name: 'many_query.fasta', text: fasta('outfmt0/many_query.fasta') }]);
  await openFiles(page, 'subject', [{ name: 'many_subject.fasta', text: fasta('outfmt0/many_subject.fasta') }]);
  await run(page, 1);
  await openFromQueue(page, 1);
  await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', String(total));
  const partial = page.locator('[data-testid="results-notice"][data-kind="outfmt0-partial"]');
  await expect(partial).toContainText(`outfmt 0 shows the alignments of the first ${shown} subjects of this query`);
  // Fewer subjects than the default hit list (500): no limit notice.
  expect(await noticeKinds(page)).toEqual([]);

  // The last subject in the engine's order: the "#" column the other way round. The sort keeps
  // the selected (first) subject, at the list's end now, and shows the list's first rows.
  const selected = await page.getByTestId('hsp-detail').getAttribute('data-hsp');
  await page.getByTestId('subject-sort-order').click();
  await expect(page.getByTestId('subject-sort-order').locator('..')).toHaveAttribute('aria-sort', 'descending');
  const subjectList = page.getByTestId('subject-list');
  await expect.poll(() => subjectList.evaluate((element) => element.scrollTop)).toBe(0);
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', selected!);
  await subjectList.evaluate((element) => (element.scrollTop = element.scrollHeight));
  await expect(subjectRows(page).and(page.locator('[aria-pressed="true"]'))).toHaveAttribute('data-order', '1');
  await subjectList.evaluate((element) => (element.scrollTop = 0));
  const last = subjectRows(page).first();
  await expect(last).toHaveAttribute('data-order', String(total));
  await expect(last.locator('[data-field="description"]')).toHaveText('not in outfmt 0');
  const sseqid = await text(last.locator('[data-field="sseqid"]'));
  await last.click();
  await expect(last).toHaveAttribute('aria-pressed', 'true');
  for (const row of await hspRows(page).all()) await expect(row.locator('[data-field="outfmt0"]')).toHaveText('not shown');
  await expect(page.getByTestId('detail-not-in-outfmt0')).toContainText(
    `outfmt 0 does not show this HSP. It shows the alignments of the first ${shown} subjects of this query`,
  );
  await expect(page.getByTestId('detail-heading')).toHaveCount(0);
  await expect(page.getByTestId('detail-section')).toHaveCount(0);
  const row = (await text(page.getByTestId('detail-row'))).replace(/\n$/, '');
  expect(row.split('\t')[1]).toBe(sseqid);
  // outfmt 6 has the HSP; outfmt 0 has no alignment heading for its subject.
  await showOutput(page, 1, 6, false);
  expect((await text(page.getByTestId('result-output'))).split('\n')).toContain(row);
  await showOutput(page, 1, 0, false);
  const escaped = sseqid.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  expect(await text(page.getByTestId('result-output'))).not.toMatch(new RegExp(`^> ?${escaped}(\\s|$)`, 'm'));

  // An explicit -max_target_seqs that the query's subjects reach (many.mts255.blastn): the
  // notice says that more subjects may match, not that hits were lost.
  const limit = BUILD_HAS_ENGINE ? 255 : 3;
  await page.getByTestId('tab-search').click();
  await page.getByTestId('param-max_target_seqs').fill(String(limit));
  await run(page, 2);
  await openFromQueue(page, 2);
  const notice = page.locator('[data-testid="results-notice"][data-kind="subject-limit"]');
  await expect(notice).toContainText(
    `This query has ${limit} subjects, the most that the search keeps (-max_target_seqs: ${limit}). More subjects may match; ` +
      'a search with a larger -max_target_seqs would show them.',
  );
  expect(await text(notice)).not.toMatch(/\b(lost|missing|dropped|discarded|omitted|truncated)\b/i);
  expect(await noticeKinds(page)).toEqual(['subject-limit']);
  if (BUILD_HAS_ENGINE) {
    // With -max_target_seqs, outfmt 0 shows the alignments of every subject kept.
    await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', String(limit));
    await expect(partial).toHaveCount(0);
  }
});

test("BLASTN: an HSP of one letter has no orientation in its record; the note points to outfmt 0's Strand= line", async ({ page }) => {
  // The example of docs/evidence/losat_web_e2c/AUTHORITY.md §P: a one-letter HSP on the minus
  // strand. The FakeEngine writes a one-letter BLASTN HSP for its third query.
  const qIdx = BUILD_HAS_ENGINE ? 0 : 2;
  await program(page, 'blastn');
  await page.getByTestId('param-task').selectOption('blastn');
  await page.getByTestId('param-word_size').fill('4');
  await paste(page, 'query', BUILD_HAS_ENGINE ? '>q\nTAGGACGG\n' : '>q1\nACGTACGTACGT\n>q2\nACGTACGTACGT\n>q\nTAGGACGG\n');
  await paste(page, 'subject', '>s\nYCAYAANTNCRGYACT\n');
  await run(page, 1);
  await openFromQueue(page, 1);
  await page.getByTestId(`query-row-${qIdx}`).click();
  await expect(page.getByTestId('subject-row-0')).toHaveAttribute('aria-pressed', 'true');
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
  await page.getByTestId('pane-dotplot').click();
  await expect(page.getByTestId('dotplot')).toContainText('One letter: the HSP record does not say its strand');
  await showOutput(page, 1, 0, false);
  expect(await text(page.getByTestId('result-output'))).toContain(section);
});

// --- selection, filters, notices, many queries ------------------------------------------------------

async function panelRun(page: Page): Promise<void> {
  await program(page, 'blastn');
  await paste(page, 'query', PANEL_QUERIES);
  await paste(page, 'subject', PANEL_SUBJECTS);
  await run(page, 1);
  await openFromQueue(page, 1);
}

test('selection by HSP identity: sorting keeps the selected HSP; another subject or query moves it', async ({ page }) => {
  await panelRun(page);
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'true');
  // A subject of the first query with two HSPs; choose its second HSP.
  const pairs = await subjectRows(page).all();
  let sIdx = -1;
  for (const row of pairs) {
    if ((await text(row.locator('[data-field="hsps"]'))).trim() === '2') {
      sIdx = Number((await row.getAttribute('data-testid'))!.replace('subject-row-', ''));
      await row.click();
      break;
    }
  }
  expect(sIdx).toBeGreaterThanOrEqual(0);
  await expect(hspRows(page)).toHaveCount(2);
  const second = hspRows(page).nth(1);
  const id = await hspId(second);
  await second.click();
  await expectSelected(page, id, sIdx);

  // The "#" column the other way round: the selected HSP is now first, and still selected.
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
      await page.getByTestId(sort).click();
      await expect(page.getByTestId(sort).locator('..')).toHaveAttribute('aria-sort', /^(ascending|descending)$/);
      await expectSelected(page, id, sIdx);
    }
  }

  // Another subject: the selection moves to its first HSP in the current order.
  const unselected = page.getByTestId('subject-list').locator('[data-testid^="subject-row-"][aria-pressed="false"]').first();
  const otherIdx = Number((await unselected.getAttribute('data-testid'))!.replace('subject-row-', ''));
  await page.getByTestId(`subject-row-${otherIdx}`).click();
  await expect(page.getByTestId(`subject-row-${otherIdx}`)).toHaveAttribute('aria-pressed', 'true');
  const firstOfOther = await hspId(hspRows(page).first());
  expect(firstOfOther).not.toBe(id);
  await expectSelected(page, firstOfOther, otherIdx);

  // Another query: the selection moves to its first subject's first HSP.
  await page.getByTestId('query-row-1').click();
  await expect(page.getByTestId('query-row-1')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'false');
  const firstSubject = subjectRows(page).first();
  const firstHsp = await hspId(hspRows(page).first());
  expect(firstHsp).toMatch(/^1:/);
  await expectSelected(page, firstHsp, Number((await firstSubject.getAttribute('data-testid'))!.replace('subject-row-', '')));
});

test('many queries; view filters change the view, not the search; the notices tell filtered out, no hits and a failed run apart', async ({
  page,
}) => {
  await panelRun(page);
  const queue = page.getByTestId('queue').locator(':scope > li');
  await expect(queue).toHaveCount(1);

  // A query without hits.
  await page.getByTestId('query-row-3').click();
  await expect(page.getByTestId('query-row-3')).toContainText('no hits');
  await expect(page.locator('[data-testid="results-notice"][data-kind="no-hits"]')).toHaveText('No hits for this query.');
  expect(await noticeKinds(page)).toEqual(['no-hits']);
  await expect(page.getByTestId('subject-table')).toHaveCount(0);
  await expect(page.getByTestId('hsp-table')).toHaveCount(0);

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
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', /^0:/);
  expect(await noticeKinds(page)).not.toContain('no-hits');
  await page.getByTestId('filter-hits-only').uncheck();
  await page.getByTestId('query-filter').fill('rec137');
  await expect(page.getByTestId('query-count')).toHaveText('1 of 150 queries');
  await expect(page.locator('[data-testid^="query-row-"]')).toHaveCount(1);
  await page.getByTestId('query-row-136').click();
  await expect(page.getByTestId('query-row-136')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', /^136:/);
  await page.getByTestId('query-filter').fill('');
  await expect(page.getByTestId('query-count')).toHaveText('150 of 150 queries');
  await list.evaluate((element) => (element.scrollTop = 0));
  await page.getByTestId('query-row-0').click();
  await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'true');

  // A subject filter that matches nothing: every HSP of the query is hidden, and one click clears it.
  const [, subjects, hsps] = /([\d,]+) subj\., ([\d,]+) HSPs/.exec(await text(page.getByTestId('query-row-0')))!;
  await page.getByTestId('filter-subject').fill('no-such-subject');
  await page.getByTestId('filter-apply').click();
  const filteredOut = page.locator('[data-testid="results-notice"][data-kind="filtered-out"]');
  await expect(filteredOut).toContainText(
    `No HSPs of this query match the view filters (${plural(Number(hsps), 'HSP')} on ${plural(Number(subjects), 'subject')} hidden).`,
  );
  expect(await noticeKinds(page)).toEqual(['filtered-out']);
  await expect(page.getByTestId('subject-table')).toHaveCount(0);
  await expect(page.getByTestId('hsp-table')).toHaveCount(0);
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

test('the dot plot: the HSPs of the pair on a canvas; zoom; choosing an HSP on it or in the list', async ({ page }) => {
  // q2 is A, which the subject holds twice: two HSPs far apart (the FakeEngine also writes
  // two HSPs for its second query).
  await program(page, 'blastn');
  await paste(page, 'query', `>q1\n${dna(21, 20)}\n>q2\n${A}\n`);
  await paste(page, 'subject', `>s1\n${A}${dna(12, 60)}${A}\n`);
  await run(page, 1);
  await openFromQueue(page, 1);
  await page.getByTestId('query-row-1').click();
  await expect(hspRows(page)).toHaveCount(2);
  const [firstId, secondId] = [await hspId(hspRows(page).nth(0)), await hspId(hspRows(page).nth(1))];
  await page.getByTestId('pane-dotplot').click();
  const canvas = page.getByTestId('dotplot-canvas');
  await expect(canvas).toHaveAttribute('data-segments', '2');
  await expect(canvas).toHaveAttribute('data-selected', firstId);
  const full = '0,60,0,180';
  await expect(canvas).toHaveAttribute('data-view', full);

  // The canvas is drawn: pixels of the HSPs' colours (forward blue, reverse orange).
  const coloured = await canvas.evaluate((element) => {
    const canvasElement = element as HTMLCanvasElement;
    const { data } = canvasElement.getContext('2d')!.getImageData(0, 0, canvasElement.width, canvasElement.height);
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      const [r, g, b] = [data[i]!, data[i + 1]!, data[i + 2]!];
      if (Math.abs(r - 31) + Math.abs(g - 95) + Math.abs(b - 191) < 60 || Math.abs(r - 194) + Math.abs(g - 65) + Math.abs(b - 12) < 60) n++;
    }
    return n;
  });
  expect(coloured).toBeGreaterThan(20);

  // Zoom in, in again, out, and back to the whole sequences.
  const span = async () => {
    const [x0, x1, y0, y1] = (await canvas.getAttribute('data-view'))!.split(',').map(Number);
    return { x: x1! - x0!, y: y1! - y0! };
  };
  await page.getByTestId('dotplot-zoom-in').click();
  await expect(canvas).not.toHaveAttribute('data-view', full);
  const once = await span();
  expect(once.x).toBeLessThan(60);
  await page.getByTestId('dotplot-zoom-in').click();
  await expect.poll(async () => (await span()).x).toBeLessThan(once.x);
  const twice = await span();
  await page.getByTestId('dotplot-zoom-out').click();
  await expect.poll(async () => (await span()).x).toBeGreaterThan(twice.x);
  await page.getByTestId('dotplot-reset').click();
  await expect(canvas).toHaveAttribute('data-view', full);

  // An HSP chosen in the list is selected on the plot; "Zoom to HSP" frames it.
  await hspRows(page).nth(1).click();
  await expect(canvas).toHaveAttribute('data-selected', secondId);
  const [qStart, qEnd] = (await text(hspRows(page).nth(1).locator('[data-field="query"]'))).split('–').map(Number);
  const [sStart, sEnd] = (await text(hspRows(page).nth(1).locator('[data-field="subject"]'))).split('–').map(Number);
  await page.getByTestId('dotplot-zoom-hsp').click();
  await expect(canvas).not.toHaveAttribute('data-view', full);
  const [x0, x1, y0, y1] = (await canvas.getAttribute('data-view'))!.split(',').map(Number);
  expect(x0).toBeLessThanOrEqual(Math.min(qStart!, qEnd!));
  expect(x1).toBeGreaterThanOrEqual(Math.max(qStart!, qEnd!));
  expect(y0).toBeLessThanOrEqual(Math.min(sStart!, sEnd!));
  expect(y1).toBeGreaterThanOrEqual(Math.max(sStart!, sEnd!));
  await page.getByTestId('dotplot-reset').click();

  // A click on the other HSP's line selects it, in the plot and in the list.
  const targets = JSON.parse((await canvas.getAttribute('data-targets'))!) as { hsp: string; x: number; y: number }[];
  const target = targets.find((t) => t.hsp === firstId)!;
  await canvas.click({ position: { x: target.x, y: target.y } });
  await expect(canvas).toHaveAttribute('data-selected', firstId);
  await expect(hspRows(page).nth(0)).toHaveAttribute('aria-pressed', 'true');
  await expect(hspRows(page).nth(1)).toHaveAttribute('aria-pressed', 'false');
  const firstQuery = await text(hspRows(page).nth(0).locator('[data-field="query"]'));
  await expect(page.getByTestId('dotplot-selected')).toContainText(`query ${firstQuery} nt`);
  // The keyboard: n selects the next HSP.
  await canvas.focus();
  await page.keyboard.press('n');
  await expect(canvas).toHaveAttribute('data-selected', secondId);
  await expect(hspRows(page).nth(1)).toHaveAttribute('aria-pressed', 'true');
});
