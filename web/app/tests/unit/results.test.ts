// The results screen's domain and application parts (S13): the HSP table and its index,
// sorting and filtering with engine values, the limits and the orientation, the
// verification badge, and the ResultsBrowser's selection by HSP identity.
import { describe, expect, it } from 'vitest';
import type { AppState, RunView } from '../../src/application/coordinator';
import { ResultsBrowser, type ResultsDeps } from '../../src/application/results';
import { Store } from '../../src/application/store';
import { hspTable, out0Range, out6Range, frame, subjectSpan } from '../../src/domain/hsp-table';
import { headingTitle } from '../../src/domain/outfmt0';
import { splitOutfmt6Row } from '../../src/domain/outfmt6';
import {
  buildResultIndex,
  filterQuery,
  hitLimits,
  orientation,
  sortHsps,
  sortSubjects,
  windowAround,
  type SubjectGroup,
} from '../../src/domain/result-index';
import {
  approvedExceptions,
  optionKey,
  verificationBadge,
  RANGE_VALUE,
  type BadgeInput,
  type VerificationTable,
} from '../../src/domain/verification';
import type { RunSnapshot } from '../../src/domain/run';
import type { HspRecord, ProgramDescription } from '../../src/ports/engine';

// --- a small run: two queries (the second without hits), three subjects ----------------------

interface Spec {
  readonly q: number;
  readonly s: number;
  readonly bits: number;
  readonly e: number;
  readonly coords: readonly [number, number, number, number];
  readonly shown?: boolean;
  readonly frames?: readonly [number | null, number | null];
}

const SPECS: readonly Spec[] = [
  { q: 0, s: 1, bits: 90, e: 1e-20, coords: [1, 50, 101, 150] },
  { q: 0, s: 1, bits: 40, e: 1e-5, coords: [60, 80, 300, 280] },
  { q: 0, s: 0, bits: 70, e: 1e-12, coords: [5, 45, 1, 41] },
  { q: 0, s: 2, bits: 70, e: 1e-12, coords: [7, 7, 9, 9], shown: false },
];

function makeRun(specs: readonly Spec[] = SPECS) {
  let out6 = '';
  let out0 = 'header\n';
  const records: HspRecord[] = [];
  const ranks = new Map<number, number>();
  const headings = new Map<string, [number, number]>();
  specs.forEach((spec, index) => {
    const [qs, qe, ss, se] = spec.coords;
    const row = [`q${spec.q}`, `s${spec.s}`, '99.000', '50', '0', '0', qs, qe, ss, se, String(spec.e), String(spec.bits)].join('\t') + '\n';
    const out6Start = out6.length;
    out6 += row;
    const shown = spec.shown ?? true;
    const key = `${spec.q}:${spec.s}`;
    if (shown && !headings.has(key)) {
      const heading = `> s${spec.s} subject number ${spec.s}\nLength=500\n\n`;
      headings.set(key, [out0.length, out0.length + heading.length]);
      out0 += heading;
    }
    const sectionStart = out0.length;
    if (shown) out0 += ` Score = ${spec.bits} bits\n\nQuery  ${qs}  ACGT  ${qe}\n`;
    const rank = ranks.get(spec.q) ?? 0;
    ranks.set(spec.q, rank + 1);
    records.push({
      index,
      q_idx: spec.q,
      s_idx: spec.s,
      rank,
      raw_score: spec.bits * 2,
      bit_score: spec.bits,
      e_value: spec.e,
      q_start: qs,
      q_end: qe,
      s_start: ss,
      s_end: se,
      query_frame: spec.frames?.[0] ?? null,
      subject_frame: spec.frames?.[1] ?? null,
      subject_length: 500,
      query_aligned: null,
      subject_aligned: null,
      out6: [out6Start, out6.length],
      out0: shown ? [sectionStart, out0.length] : null,
      out0_subject: shown ? headings.get(key)! : null,
    });
  });
  return { records, out6: new TextEncoder().encode(out6), out0: new TextEncoder().encode(out0) };
}

describe('outfmt 6 rows and outfmt 0 headings', () => {
  it('splits a row into its twelve fields, and refuses another number of fields', () => {
    const row = splitOutfmt6Row('q1\ts1\t98.684\t76\t1\t0\t1\t76\t1\t76\t2.37e-38\t138\n');
    expect(row).toMatchObject({ qseqid: 'q1', sseqid: 's1', pident: '98.684', evalue: '2.37e-38', bitscore: '138', send: '76' });
    expect(() => splitOutfmt6Row('q1\ts1\t1\n')).toThrow(/12 fields, this one 3/);
  });

  it('takes the title of a subject heading as written', () => {
    expect(headingTitle('> s1 a long title that NCBI\nwrapped here\nLength=500\n\n')).toBe('s1 a long title that NCBI\nwrapped here');
  });
});

describe('the HSP table and its index', () => {
  const { records } = makeRun();
  // The records in another order: the index follows the HSP index, not the record order.
  const table = hspTable([records[2]!, records[0]!, records[3]!, records[1]!]);
  const index = buildResultIndex(table);

  it('keeps the small fields of each record, with null as -1 or 0', () => {
    expect(table.count).toBe(4);
    const row = [...table.index].indexOf(3);
    expect(out0Range(table, row)).toBeUndefined();
    expect(out6Range(table, row)).toEqual(records[3]!.out6);
    expect(frame(table.queryFrame, row)).toBeUndefined();
  });

  it('groups the HSPs by query and subject in the engine order', () => {
    expect([...index.queries.keys()]).toEqual([0]);
    const query = index.queries.get(0)!;
    expect(query.hspCount).toBe(4);
    expect(query.subjects.map((s) => s.sIdx)).toEqual([1, 0, 2]);
    expect(query.subjects[0]!.rows.map((row) => table.index[row])).toEqual([0, 1]);
  });

  it('sorts subjects and HSPs by engine values, ties in the engine order', () => {
    const subjects = index.queries.get(0)!.subjects;
    const order = (list: readonly SubjectGroup[]) => list.map((s) => s.sIdx);
    expect(order(sortSubjects(subjects, { key: 'bitScore', descending: true }, table, () => 0))).toEqual([1, 0, 2]);
    expect(order(sortSubjects(subjects, { key: 'bitScore', descending: false }, table, () => 0))).toEqual([0, 2, 1]);
    expect(order(sortSubjects(subjects, { key: 'eValue', descending: false }, table, () => 0))).toEqual([1, 0, 2]);
    expect(order(sortSubjects(subjects, { key: 'hsps', descending: true }, table, () => 0))).toEqual([1, 0, 2]);
    expect(order(sortSubjects(subjects, { key: 'length', descending: true }, table, (s) => [10, 5, 20][s]!))).toEqual([2, 0, 1]);
    const rows = subjects[0]!.rows;
    expect(sortHsps(rows, { key: 'bitScore', descending: false }, table).map((r) => table.index[r])).toEqual([1, 0]);
    expect(sortHsps(rows, { key: 'sStart', descending: true }, table).map((r) => table.index[r])).toEqual([1, 0]);
  });

  it('filters HSPs by E value, bit score and subject text, and counts what it hides', () => {
    const query = index.queries.get(0)!;
    const ids = (s: SubjectGroup) => [`s${s.sIdx}`];
    const byE = filterQuery(query, { maxEValue: 1e-10 }, table, ids);
    expect(byE.subjects.map((s) => s.sIdx)).toEqual([1, 0, 2]);
    expect(byE.subjects[0]!.rows).toHaveLength(1);
    expect(byE).toMatchObject({ hiddenHsps: 1, hiddenSubjects: 0 });
    const byBits = filterQuery(query, { minBitScore: 80 }, table, ids);
    expect(byBits.subjects.map((s) => s.sIdx)).toEqual([1]);
    expect(byBits).toMatchObject({ hiddenHsps: 3, hiddenSubjects: 2 });
    expect(filterQuery(query, { subjectText: 'S2' }, table, ids).subjects.map((s) => s.sIdx)).toEqual([2]);
    expect(filterQuery(query, { minBitScore: 1000 }, table, ids)).toMatchObject({ subjects: [], hiddenHsps: 4, hiddenSubjects: 3 });
    // Query filters do not hide HSPs.
    expect(filterQuery(query, { queriesWithHitsOnly: true, queryText: 'x' }, table, ids).subjects).toBe(query.subjects);
  });

  it("gives an HSP's subject coordinates in ascending order (NCBI's Range)", () => {
    const rowOf = (i: number) => [...table.index].indexOf(i);
    expect(subjectSpan(table, rowOf(0))).toEqual([101, 150]);
    expect(subjectSpan(table, rowOf(1))).toEqual([280, 300]);
    expect(subjectSpan(table, rowOf(3))).toEqual([9, 9]);
  });

  it('puts a window of Ranges around a position, within the list', () => {
    expect(windowAround(30, 100, 25)).toEqual({ start: 5, end: 56 });
    expect(windowAround(0, 100, 25)).toEqual({ start: 0, end: 26 });
    expect(windowAround(99, 100, 25)).toEqual({ start: 74, end: 100 });
    expect(windowAround(3, 5, 25)).toEqual({ start: 0, end: 5 });
    expect(windowAround(-1, 5, 25)).toEqual({ start: 0, end: 5 });
    expect(windowAround(7, 5, 1)).toEqual({ start: 3, end: 5 });
    expect(windowAround(0, 0, 25)).toEqual({ start: 0, end: 0 });
  });

  it('gives the orientation from coordinates, frames, or not at all', () => {
    const nucleotide = { query: 'nucleotide', subject: 'nucleotide' } as const;
    const rowOf = (i: number) => [...table.index].indexOf(i);
    expect(orientation(table, rowOf(0), nucleotide)).toBe('forward');
    expect(orientation(table, rowOf(1), nucleotide)).toBe('reverse');
    expect(orientation(table, rowOf(3), nucleotide)).toBe('unknown');
    expect(orientation(table, rowOf(3), { query: 'protein', subject: 'protein' })).toBe('forward');
    const framed = hspTable(makeRun([{ q: 0, s: 0, bits: 1, e: 1, coords: [3, 3, 9, 9], frames: [1, -2] }]).records);
    expect(orientation(framed, 0, nucleotide)).toBe('reverse');
    const both = hspTable(makeRun([{ q: 0, s: 0, bits: 1, e: 1, coords: [30, 1, 90, 61], frames: [-1, -3] }]).records);
    expect(orientation(both, 0, nucleotide)).toBe('forward');
  });
});

describe('hitLimits', () => {
  const described = [
    { flag: '-max_target_seqs', help: 'Maximum number of aligned sequences to keep [default: 500]' },
    { flag: '-max_hsps', help: 'Maximum number of HSPs per subject' },
  ];

  it('reads the limits of the argv, or the default that the engine describes', () => {
    expect(hitLimits(['blastp', '-query', 'q', '-subject', 's'], described)).toEqual({ maxTargetSeqs: 500 });
    expect(hitLimits(['blastp', '-max_target_seqs', '5', '-max_hsps', '2'], described)).toEqual({ maxTargetSeqs: 5, maxHsps: 2 });
    expect(hitLimits(['blastp'], [{ flag: '-max_target_seqs', help: 'x', defaultValue: '100' }])).toEqual({ maxTargetSeqs: 100 });
  });

  it('claims no limit that it cannot read', () => {
    expect(hitLimits(['blastp'], [{ flag: '-max_target_seqs', help: 'Maximum number' }])).toEqual({});
    expect(hitLimits(['blastp', '-max_target_seqs', '1e3'], [])).toEqual({});
  });
});

// --- verification ---------------------------------------------------------------------------

const grammar = {
  takesValue: (flag: string) => ['-task', '-evalue', '-query_loc', '-subject_loc', '-db_gencode', '-penalty'].includes(flag),
  defaultTask: 'megablast',
};

describe('optionKey', () => {
  it('is the same for the same options in any order, with the default task written or not', () => {
    const a = optionKey(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '1e-5', '-lcase_masking'], grammar);
    const b = optionKey(['-lcase_masking', '-task', 'megablast', '-outfmt', '7', '-evalue', '1e-5', '-num_threads', '4'], grammar);
    expect(a).toBe(b);
    expect(a).toBe(JSON.stringify([['-evalue', '1e-5'], ['-lcase_masking', ''], ['-task', 'megablast']]));
  });

  it('keeps a value that starts with "-", and names a range without its coordinates', () => {
    expect(optionKey(['-penalty', '-3', '-query_loc', '10-200'], grammar)).toBe(
      JSON.stringify([['-penalty', '-3'], ['-query_loc', RANGE_VALUE], ['-task', 'megablast']]),
    );
  });

  it('reads nothing from words that are not options', () => {
    expect(optionKey(['-evalue'], grammar)).toBeUndefined();
    expect(optionKey(['-lcase_masking', 'stray'], grammar)).toBeUndefined();
  });
});

describe('verificationBadge', () => {
  const key = optionKey(['blastn'], grammar)!;
  const table: VerificationTable = {
    ncbi: '2.17.0',
    sources: [],
    programs: {
      blastn: {
        optionSets: { [key]: [0, 6, 7], [optionKey(['-evalue', '1'], grammar)!]: [0, 6] },
        browser: { formats: [0, 6, 7], paths: ['serial', 'threaded'], threads: [1, 2, 4], browsers: ['Chromium', 'Firefox', 'WebKit'] },
      },
    },
  };
  const input: BadgeInput = {
    program: 'blastn',
    argv: ['blastn', '-query', 'q', '-subject', 's'],
    formats: [0, 6, 7],
    runtimePath: 'threaded',
    threads: 4,
    grammar,
  };

  it('certifies a run inside a profile, and says what that means', () => {
    const badge = verificationBadge(input, table);
    expect(badge.level).toBe('certified');
    expect(badge.label).toBe('Certified profile');
    expect(badge.details.join(' ')).toMatch(/not a comparison of your inputs with NCBI/);
  });

  it('says why a run is outside the certified profiles', () => {
    const reasons = (change: Partial<BadgeInput>) => verificationBadge({ ...input, ...change }, table);
    expect(reasons({ argv: ['blastn', '-evalue', '1'] }).details.join(' ')).toMatch(/for outfmt 0 and 6, not for outfmt 7/);
    expect(reasons({ argv: ['blastn', '-evalue', '2'] }).details.join(' ')).toMatch(/were not compared with NCBI BLAST\+ 2\.17\.0/);
    expect(reasons({ threads: 3 }).details.join(' ')).toMatch(/used 3 threads; the browser runtime was checked with 1, 2 and 4 threads/);
    expect(reasons({ program: 'tblastx', argv: ['tblastx'] }).details.join(' ')).toMatch(/browser runtime of tblastx was not checked/);
    expect(reasons({ threads: 3 }).label).toBe('Engine-supported, outside certified profile');
    expect(reasons({ runtimePath: 'fake' }).level).toBe('development');
  });

  it('names the approved exceptions of a non-default subject genetic code', () => {
    expect(approvedExceptions('tblastn', ['tblastn', '-db_gencode', '4'])[0]).toMatch(/PD-TLOSAN-LOCAL-GENCODE-32.*-db_gencode 4/);
    expect(approvedExceptions('tblastx', ['tblastx', '-db_gencode', '11'])).toHaveLength(1);
    expect(approvedExceptions('tblastx', ['tblastx', '-db_gencode', '1'])).toEqual([]);
    expect(approvedExceptions('blastn', ['blastn', '-db_gencode', '4'])).toEqual([]);
  });
});

// --- ResultsBrowser ---------------------------------------------------------------------------

const DESCRIPTION: ProgramDescription = {
  program: 'blastn',
  formats: [0, 6, 7],
  parameters: [
    { flag: '-task', help: 'Task', takesValue: true, defaultValue: 'megablast' },
    { flag: '-max_target_seqs', help: 'Maximum number of aligned sequences to keep (default: 500)', takesValue: true },
    { flag: '-max_hsps', help: 'Maximum number of HSPs per subject', takesValue: true },
  ],
};

function snapshot(runId: string, argv: readonly string[] = ['blastn', '-query', 'q.fa', '-subject', 's.fa']): RunSnapshot {
  const input = (records: readonly string[]) => ({
    name: 'x.fa',
    bytes: new Uint8Array(),
    sha256: '0',
    revisionIds: [],
    records: records.map((id) => ({ id, length: 500 })),
  });
  return { runId, number: 1, program: 'blastn', argv, query: input(['q0', 'q1']), subject: input(['s0', 's1', 's2']), requestedThreads: 1, queuedAt: 0 };
}

function setup(options: { status?: RunView['status']; specs?: readonly Spec[]; argv?: readonly string[]; unreadable?: boolean } = {}) {
  const run = makeRun(options.specs);
  const view: RunView = {
    snapshot: snapshot('r1', options.argv),
    status: options.status ?? 'completed',
    record: { runtimePath: 'serial', threads: 1, ...(options.status === 'failed' ? { error: 'boom' } : {}) },
  };
  const runs = new Store<AppState>({ runs: [view] });
  const pending: Array<() => void> = [];
  let hold = false;
  let tableReads = 0;
  let rangeReads = 0;
  let failRanges = false;
  const deps: ResultsDeps = {
    data: {
      readHitTable: async () => {
        tableReads++;
        if (options.unreadable === true) throw new Error('the stored records are gone');
        return hspTable(run.records);
      },
      readOutput: async (_id, format) => (format === 6 ? run.out6.slice() : run.out0.slice()),
      readOutputRange: (_id, _format, start, end) =>
        new Promise((resolve, reject) => {
          rangeReads++;
          const answer = () => (failRanges ? reject(new Error('the stored output is gone')) : resolve(run.out0.slice(start, end)));
          if (hold) pending.push(answer);
          else answer();
        }),
      readDiagnostics: async () => '',
    },
    describe: async () => DESCRIPTION,
    runs,
    verification: { ncbi: '2.17.0', sources: [], programs: {} },
  };
  const results = new ResultsBrowser(deps);
  return {
    results,
    runs,
    view,
    hold: (on: boolean) => (hold = on),
    /** Answers the held reads in the order they were made, or the last first. */
    release: (lastFirst = false) => {
      const answers = pending.splice(0);
      (lastFirst ? answers.reverse() : answers).forEach((answer) => answer());
    },
    tableReads: () => tableReads,
    rangeReads: () => rangeReads,
    failRanges: (on: boolean) => (failRanges = on),
  };
}

const settle = () => new Promise((resolve) => setTimeout(resolve, 0));

describe('ResultsBrowser', () => {
  it('opens a run on its first query with hits, its first subject and HSP, and reads the HSP as written', async () => {
    const { results } = setup();
    await results.open('r1');
    await settle();
    const state = results.state.get();
    expect(state.phase).toBe('ready');
    expect(state.queries.map((q) => [q.id, q.subjects, q.hsps])).toEqual([
      ['q0', 3, 4],
      ['q1', 0, 0],
    ]);
    expect(state.subjects.map((s) => [s.sIdx, s.first.sseqid, s.first.bitscore, s.hspCount, s.inOutfmt0])).toEqual([
      [1, 's1', '90', 2, true],
      [0, 's0', '70', 1, true],
      [2, 's2', '70', 1, false],
    ]);
    expect(state.hsp).toEqual({ runId: 'r1', qIdx: 0, rank: 0 });
    expect(state.hsps.map((h) => [h.id.rank, h.orientation])).toEqual([
      [0, 'forward'],
      [1, 'reverse'],
    ]);
    expect(state.detail).toMatchObject({ state: 'ready', heading: '> s1 subject number 1\nLength=500\n\n' });
    expect(state.detail!.section).toBe(' Score = 90 bits\n\nQuery  1  ACGT  50\n');
    expect(state.detail!.row).toBe('q0\ts1\t99.000\t50\t0\t0\t1\t50\t101\t150\t1e-20\t90\n');
    expect(state.outfmt0Subjects).toBe(2);
    expect(state.atSubjectLimit).toBe(false);
    expect(state.loaded!.badge.level).toBe('outside');
  });

  it('follows an HSP of another subject, and shows a query without hits', async () => {
    const { results } = setup();
    await results.open('r1');
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 3 });
    await settle();
    let state = results.state.get();
    expect(state.sIdx).toBe(2);
    expect(state.hsps[0]!.orientation).toBe('unknown');
    expect(state.detail).toMatchObject({ state: 'ready' });
    expect(state.detail!.section).toBeUndefined();
    results.selectQuery(1);
    state = results.state.get();
    expect(state.queryTotals).toEqual({ subjects: 0, hsps: 0 });
    expect(state.subjects).toEqual([]);
    expect(state.hsp).toBeUndefined();
  });

  it('drops a detail that arrives after the selection moved on', async () => {
    const { results, hold, release } = setup();
    await results.open('r1');
    await settle();
    hold(true);
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 1 });
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 2 });
    // The reads for rank 2 answer first: rank 1's detail arrives last, and is dropped.
    release(true);
    await settle();
    hold(false);
    const { detail, hsp } = results.state.get();
    expect(hsp).toEqual({ runId: 'r1', qIdx: 0, rank: 2 });
    expect(detail!.id).toEqual(hsp);
    expect(detail!.section).toBe(' Score = 70 bits\n\nQuery  5  ACGT  45\n');
  });

  it('keeps the selected HSP through a filter that shows it, and moves off one that hides it', async () => {
    const { results } = setup();
    await results.open('r1');
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 1 });
    results.setFilters({ maxEValue: 1e-4 });
    expect(results.state.get().hsp).toEqual({ runId: 'r1', qIdx: 0, rank: 1 });
    results.setFilters({ maxEValue: 1e-10 });
    expect(results.state.get().hsp).toEqual({ runId: 'r1', qIdx: 0, rank: 0 });
    expect(results.state.get().hidden).toEqual({ subjects: 0, hsps: 1 });
    results.setFilters({ minBitScore: 1000 });
    const state = results.state.get();
    expect(state.subjects).toEqual([]);
    expect(state.hsp).toBeUndefined();
    expect(state.hidden).toEqual({ subjects: 3, hsps: 4 });
    results.setFilters({ queriesWithHitsOnly: true });
    expect(results.state.get().queries.map((q) => q.id)).toEqual(['q0']);
  });

  it('moves off a query that a query filter hides, and selects none when no query is listed', async () => {
    const { results } = setup();
    await results.open('r1');
    results.selectQuery(1);
    results.setFilters({ queriesWithHitsOnly: true });
    let state = results.state.get();
    expect(state.queries.map((q) => q.id)).toEqual(['q0']);
    expect(state.qIdx).toBe(0);
    expect(state.hsp).toEqual({ runId: 'r1', qIdx: 0, rank: 0 });
    expect(state.queryTotals).toEqual({ subjects: 3, hsps: 4 });

    results.setFilters({ queryText: 'none' });
    state = results.state.get();
    expect(state.queries).toEqual([]);
    expect([state.qIdx, state.sIdx, state.hsp, state.queryTotals]).toEqual([undefined, undefined, undefined, undefined]);
    expect(state.subjects).toEqual([]);

    results.setFilters({ queryText: 'q1' });
    state = results.state.get();
    expect(state.qIdx).toBe(1);
    expect(state.queryTotals).toEqual({ subjects: 0, hsps: 0 });
    results.setFilters({ queryText: 'q' });
    expect(results.state.get().qIdx).toBe(1);
  });

  it('opens another run without view filters, and keeps them for the same run', async () => {
    const { results, runs, view } = setup();
    runs.set({ runs: [view, { ...view, snapshot: snapshot('r2') }] });
    await results.open('r1');
    results.setFilters({ subjectText: 's0' });
    await results.open('r1');
    expect(results.state.get().filters).toEqual({ subjectText: 's0' });
    expect(results.state.get().subjects.map((s) => s.sIdx)).toEqual([0]);
    await results.open('r2');
    const state = results.state.get();
    expect(state.runId).toBe('r2');
    expect(state.filters).toEqual({});
    expect(state.subjects.map((s) => s.sIdx)).toEqual([1, 0, 2]);
  });

  it('says that a query reached the hit list limit, and a subject the HSP limit', async () => {
    const { results } = setup({ argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-max_target_seqs', '3', '-max_hsps', '2'] });
    await results.open('r1');
    const state = results.state.get();
    expect(state.atSubjectLimit).toBe(true);
    expect(state.subjects.map((s) => s.atHspLimit)).toEqual([true, false, false]);
  });

  it('shows why a run has no results, and opens it when it completes', async () => {
    const failed = setup({ status: 'failed' });
    await failed.results.open('r1');
    expect(failed.results.state.get()).toMatchObject({ phase: 'unavailable', message: 'Run 1 failed: boom. A failed run keeps no results.' });

    const running = setup({ status: 'running' });
    await running.results.open('r1');
    expect(running.results.state.get().phase).toBe('unavailable');
    running.runs.set({ runs: [{ ...running.view, status: 'completed' }] });
    await settle();
    await settle();
    expect(running.results.state.get().phase).toBe('ready');
  });

  it('reads a run whose results could not be read once, not at each change of the runs', async () => {
    const { results, runs, view, tableReads } = setup({ unreadable: true });
    await results.open('r1');
    expect(results.state.get().phase).toBe('failed');
    expect(results.state.get().message).toContain('The results of run 1 could not be read: ');
    runs.set({ runs: [view, { ...view, snapshot: snapshot('r2'), status: 'running' }] });
    await settle();
    expect(tableReads()).toBe(1);
    expect(results.state.get().phase).toBe('failed');
  });

  it('reads the headings of a run opened while those of the run before were being read', async () => {
    const { results, runs, view, hold, release } = setup();
    runs.set({ runs: [view, { ...view, snapshot: snapshot('r2') }] });
    await results.open('r1');
    await settle();
    hold(true);
    results.requestHeadings([0, 1]);
    hold(false);
    await results.open('r2');
    await settle();
    results.requestHeadings([1]);
    await settle();
    release();
    await settle();
    const state = results.state.get();
    expect(state.runId).toBe('r2');
    expect([...state.headings]).toEqual([[1, '> s1 subject number 1\nLength=500\n\n']]);
  });

  it('sorts subjects by the values that the list shows: the first HSP and the HSP count before the filters', async () => {
    const { results } = setup({
      specs: [
        { q: 0, s: 1, bits: 60, e: 1e-10, coords: [1, 50, 101, 150] },
        { q: 0, s: 1, bits: 40, e: 1e-20, coords: [60, 80, 300, 280] },
        { q: 0, s: 0, bits: 50, e: 1e-16, coords: [5, 45, 1, 41] },
      ],
    });
    await results.open('r1');
    // The filter hides the first HSP of s1 and keeps its second.
    results.setFilters({ maxEValue: 1e-15 });
    const shown = () => results.state.get().subjects.map((s) => [s.sIdx, s.first.bitscore, s.first.evalue, s.hspCount, s.rows.length]);
    results.setSubjectSort({ key: 'eValue', descending: false });
    expect(shown()).toEqual([
      [0, '50', '1e-16', 1, 1],
      [1, '60', '1e-10', 2, 1],
    ]);
    results.setSubjectSort({ key: 'bitScore', descending: true });
    expect(shown().map(([sIdx]) => sIdx)).toEqual([1, 0]);
    results.setSubjectSort({ key: 'hsps', descending: false });
    expect(shown().map(([sIdx]) => sIdx)).toEqual([0, 1]);
  });

  it("lists the selected subject's Ranges in the engine's order, numbered among all its HSPs, as the filters show them", async () => {
    const { results } = setup({
      specs: [
        { q: 0, s: 1, bits: 60, e: 1e-10, coords: [1, 50, 101, 150] },
        { q: 0, s: 0, bits: 50, e: 1e-16, coords: [5, 45, 1, 41] },
        { q: 0, s: 1, bits: 40, e: 1e-20, coords: [60, 80, 300, 280] },
      ],
    });
    await results.open('r1');
    const ranges = () => results.state.get().ranges.map((r) => [r.id.rank, r.n, r.from, r.to, r.inOutfmt0]);
    expect(results.state.get().sIdx).toBe(1);
    expect(ranges()).toEqual([
      [0, 1, 101, 150, true],
      [2, 2, 280, 300, true],
    ]);
    // The HSP table's order does not change the Ranges'.
    results.setHspSort({ key: 'bitScore', descending: false });
    expect(ranges()).toEqual([
      [0, 1, 101, 150, true],
      [2, 2, 280, 300, true],
    ]);
    // A filter that hides the first HSP keeps the second's number.
    results.setFilters({ maxEValue: 1e-15 });
    expect(ranges()).toEqual([[2, 2, 280, 300, true]]);
    results.selectSubject(0);
    expect(ranges()).toEqual([[1, 1, 1, 41, true]]);
    results.setFilters({ minBitScore: 1000 });
    expect(results.state.get().ranges).toEqual([]);
  });

  it('reads a section of outfmt 0 once per run, none for an HSP that outfmt 0 does not show, and again after a failed read', async () => {
    const { results, runs, view, rangeReads, failRanges } = setup();
    runs.set({ runs: [view, { ...view, snapshot: snapshot('r2') }] });
    await results.open('r1');
    await settle();
    // The selected HSP's detail has read its section: the Range reads none.
    const reads = rangeReads();
    expect(await results.readSection({ runId: 'r1', qIdx: 0, rank: 0 })).toBe(' Score = 90 bits\n\nQuery  1  ACGT  50\n');
    expect(rangeReads()).toBe(reads);
    expect(await results.readSection({ runId: 'r1', qIdx: 0, rank: 1 })).toBe(' Score = 40 bits\n\nQuery  60  ACGT  80\n');
    expect(await results.readSection({ runId: 'r1', qIdx: 0, rank: 1 })).toBe(' Score = 40 bits\n\nQuery  60  ACGT  80\n');
    expect(rangeReads()).toBe(reads + 1);
    expect(await results.readSection({ runId: 'r1', qIdx: 0, rank: 3 })).toBeUndefined();
    expect(await results.readSection({ runId: 'r2', qIdx: 0, rank: 1 })).toBeUndefined();
    expect(rangeReads()).toBe(reads + 1);
    // Another run reads its own sections.
    await results.open('r2');
    await settle();
    const before = rangeReads();
    failRanges(true);
    await expect(results.readSection({ runId: 'r2', qIdx: 0, rank: 2 })).rejects.toThrow('the stored output is gone');
    failRanges(false);
    expect(await results.readSection({ runId: 'r2', qIdx: 0, rank: 2 })).toBe(' Score = 70 bits\n\nQuery  5  ACGT  45\n');
    expect(rangeReads()).toBe(before + 2);
  });

  it("reads a subject's heading once per run for the selected HSP's detail and the lists, and finds HSPs by rank", async () => {
    const { results, rangeReads } = setup();
    await results.open('r1');
    await settle();
    // The detail of the first HSP read its subject's heading and its section.
    expect(rangeReads()).toBe(2);
    results.requestHeadings([1]);
    await settle();
    expect([...results.state.get().headings]).toEqual([[1, '> s1 subject number 1\nLength=500\n\n']]);
    expect(rangeReads()).toBe(2);
    // Another HSP of the subject reads only its section.
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 1 });
    await settle();
    expect(results.state.get().detail).toMatchObject({ state: 'ready', heading: '> s1 subject number 1\nLength=500\n\n' });
    expect(rangeReads()).toBe(3);
    // An HSP of another subject: its query's rank finds its row.
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 2 });
    await settle();
    expect(results.state.get()).toMatchObject({ sIdx: 0, detail: { state: 'ready', heading: '> s0 subject number 0\nLength=500\n\n' } });
    // A rank that the query does not have selects nothing.
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 9 });
    results.selectHsp({ runId: 'r1', qIdx: 1, rank: 0 });
    expect(results.state.get().hsp).toEqual({ runId: 'r1', qIdx: 0, rank: 2 });
  });

  it("names the task of a run: the argv's, else the engine's default", async () => {
    const byDefault = setup();
    await byDefault.results.open('r1');
    expect(byDefault.results.state.get().loaded!.task).toBe('megablast');
    const given = setup({ argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-task', 'blastn'] });
    await given.results.open('r1');
    expect(given.results.state.get().loaded!.task).toBe('blastn');
  });
});
