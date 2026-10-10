// The candidate tray (application/candidates.ts) and what the results browser gives it
// (application/results.ts: candidate sources, every HSP of a subject, the row marks of the
// subject list, and going back to a candidate's result). Runs, records and residues are made
// here; the engine tests (extraction-engine.test.ts) take candidates from real searches.
import { describe, expect, it } from 'vitest';
import {
  ALIGNMENTS_FILE,
  CandidateTray,
  candidateKey,
  FASTA_MIME,
  SEQUENCES_FILE,
  type CandidateSource,
  type CandidateTrayDeps,
} from '../../src/application/candidates';
import type { AppState, RunView } from '../../src/application/coordinator';
import { ResultsBrowser, type HspId } from '../../src/application/results';
import { Store } from '../../src/application/store';
import type { Interval } from '../../src/domain/coordinates';
import { hspTable } from '../../src/domain/hsp-table';
import { filtersHiding } from '../../src/domain/result-index';
import type { ProgramId } from '../../src/domain/programs';
import type { RunStatus } from '../../src/domain/run';
import type { RecordResidues } from '../../src/ports/data';
import type { HspRecord } from '../../src/ports/engine';

const encoder = new TextEncoder();
const decoder = new TextDecoder();

// --- runs and their HSPs ----------------------------------------------------------------------------

interface HspSpec {
  readonly q: number;
  readonly s: number;
  /** q_start, q_end, s_start, s_end */
  readonly coords: readonly [number, number, number, number];
  readonly bits?: number;
  readonly e?: number;
  readonly frames?: readonly [number | null, number | null];
  readonly aligned?: readonly [string, string];
}

interface RunSpec {
  readonly runId: string;
  readonly number: number;
  readonly status?: RunStatus;
  readonly program?: ProgramId;
  readonly title?: string;
  readonly queries?: ReadonlyArray<readonly [string, number]>;
  readonly subjects?: ReadonlyArray<readonly [string, number]>;
  readonly hsps?: readonly HspSpec[];
}

const QUERIES: ReadonlyArray<readonly [string, number]> = [
  ['q0', 300],
  ['q1', 200],
];
const SUBJECTS: ReadonlyArray<readonly [string, number]> = [
  ['s0', 100],
  ['', 50],
  ['s2', 400],
];

function runView(spec: RunSpec): RunView {
  const input = (name: string, records: ReadonlyArray<readonly [string, number]>) => ({
    name,
    bytes: new Uint8Array(),
    sha256: `sha-${spec.runId}-${name}`,
    revisionIds: [`${spec.runId}:${name}`],
    records: records.map(([id, length]) => ({ id, length })),
  });
  const program = spec.program ?? 'blastn';
  return {
    snapshot: {
      runId: spec.runId,
      number: spec.number,
      program,
      ...(spec.title === undefined ? {} : { title: spec.title }),
      argv: [program, '-query', 'q.fa', '-subject', 's.fa', '-evalue', '1e-5'],
      query: input('q.fa', spec.queries ?? QUERIES),
      subject: input('s.fa', spec.subjects ?? SUBJECTS),
      requestedThreads: 1,
      queuedAt: 0,
    },
    status: spec.status ?? 'completed',
    record: { runtimePath: 'serial', threads: 1, engineBuild: 'build-7', endedAt: 1000 + spec.number },
  };
}

/** The HSP records and outfmt 6 text of a run; ranks follow the order of the specs per query. */
function outputs(specs: readonly HspSpec[], subjects: ReadonlyArray<readonly [string, number]> = SUBJECTS) {
  let out6 = '';
  const ranks = new Map<number, number>();
  const records: HspRecord[] = specs.map((spec, index) => {
    const [qs, qe, ss, se] = spec.coords;
    const start = out6.length;
    const sseqid = subjects[spec.s]![0] || `Subject_${spec.s + 1}`;
    out6 += [`q${spec.q}`, sseqid, '99.000', '50', '0', '0', qs, qe, ss, se, String(spec.e ?? 1e-10), String(spec.bits ?? 50)].join('\t') + '\n';
    const rank = ranks.get(spec.q) ?? 0;
    ranks.set(spec.q, rank + 1);
    return {
      index,
      q_idx: spec.q,
      s_idx: spec.s,
      rank,
      raw_score: 1,
      bit_score: spec.bits ?? 50,
      e_value: spec.e ?? 1e-10,
      q_start: qs,
      q_end: qe,
      s_start: ss,
      s_end: se,
      query_frame: spec.frames?.[0] ?? null,
      subject_frame: spec.frames?.[1] ?? null,
      subject_length: subjects[spec.s]![1],
      query_aligned: spec.aligned?.[0] ?? null,
      subject_aligned: spec.aligned?.[1] ?? null,
      out6: [start, out6.length],
      out0: null,
      out0_subject: null,
    };
  });
  return { records, out6: encoder.encode(out6) };
}

const SPECS: readonly HspSpec[] = [
  { q: 0, s: 2, coords: [1, 60, 301, 360], bits: 90, e: 1e-30, aligned: ['ACGT-A', 'ACGTTA'] },
  { q: 0, s: 2, coords: [70, 100, 120, 90], bits: 40, e: 1e-5, aligned: ['acg', 'ACG'] },
  { q: 0, s: 0, coords: [5, 45, 1, 41], bits: 70, e: 1e-12 },
  { q: 0, s: 1, coords: [7, 7, 9, 9], bits: 20, e: 1e-2 },
  { q: 1, s: 0, coords: [10, 30, 80, 100], bits: 30, e: 1e-4 },
];

/** A results browser over runs whose HSPs the specs give. */
function browser(runs: Store<AppState>, hsps: ReadonlyMap<string, readonly HspSpec[]>) {
  const subjectsOf = (runId: string) =>
    runs
      .get()
      .runs.find((run) => run.snapshot.runId === runId)!
      .snapshot.subject.records.map(({ id, length }) => [id, length] as const);
  const made = new Map([...hsps].map(([runId, specs]) => [runId, outputs(specs, subjectsOf(runId))]));
  return new ResultsBrowser({
    data: {
      readHitTable: async (runId) => hspTable(made.get(runId)!.records),
      readOutput: async (runId, format) => (format === 6 ? made.get(runId)!.out6.slice() : new Uint8Array()),
      readOutputRange: async () => new Uint8Array(),
      readDiagnostics: async () => '',
    },
    describe: async (program) => ({ program, formats: [0, 6, 7], parameters: [] }),
    runs,
    verification: { ncbi: '2.17.0', sources: [], programs: {} },
  });
}

// --- the tray's dependencies --------------------------------------------------------------------------

/** Letters of each record of a run's inputs, by revision and position (as `readResidues` reads them). */
function letters(length: number, seed: number): string {
  let out = '';
  for (let i = 0; i < length; i++) out += 'ACGTacgu'[(i * 7 + seed * 13 + ((i * i) >> 3)) % 8];
  return out;
}

function trayDeps(runs: Store<AppState>, options: { fail?: string } = {}) {
  const saved: Array<{ name: string; text: string; mime: string }> = [];
  const residueCalls: Array<{ revisionIds: readonly string[]; position: number; intervals: readonly Interval[] }> = [];
  const hspCalls: Array<{ runId: string; indices: readonly number[] }> = [];
  const records = new Map<string, readonly HspRecord[]>();
  let clock = 5000;
  const deps: CandidateTrayDeps = {
    runs,
    data: {
      readResidues: async (revisionIds, position, intervals): Promise<RecordResidues> => {
        residueCalls.push({ revisionIds, position, intervals });
        if (options.fail !== undefined) throw new Error(options.fail);
        const view = runs.get().runs.find((run) => run.snapshot.query.revisionIds[0] === revisionIds[0] || run.snapshot.subject.revisionIds[0] === revisionIds[0])!;
        const role = view.snapshot.query.revisionIds[0] === revisionIds[0] ? 'query' : 'subject';
        const record = view.snapshot[role].records[position]!;
        const text = letters(record.length, position + (role === 'query' ? 100 : 0));
        return {
          origin: { sourceId: 'src', sourceName: 'x.fa', revisionId: revisionIds[0]!, recordIndex: position, id: record.id, length: record.length, sha256: '0' },
          residues: intervals.map(({ from, to }) => encoder.encode(text.slice(from - 1, to))),
        };
      },
      readHspRecords: async (runId, indices) => {
        hspCalls.push({ runId, indices });
        return indices.map((index) => records.get(runId)![index]!);
      },
    },
    downloader: { save: (name, bytes, mime) => saved.push({ name, text: decoder.decode(bytes), mime }) },
    now: () => clock++,
  };
  return { deps, saved, residueCalls, hspCalls, records };
}

/** Candidate sources of runs made here, through the results browser of each run. */
async function sourcesOf(runs: Store<AppState>, hsps: ReadonlyMap<string, readonly HspSpec[]>, runId: string, ids?: readonly HspId[]) {
  const results = browser(runs, hsps);
  await results.open(runId);
  const all = hsps.get(runId)!.map((_, index) => index);
  const table = results.state.get().loaded!.index.table;
  return results.candidateSources(ids ?? all.map((row) => ({ runId, qIdx: table.qIdx[row]!, rank: table.rank[row]! })));
}

function setup(specs: readonly RunSpec[] = [{ runId: 'r1', number: 1, hsps: SPECS }]) {
  const runs = new Store<AppState>({ runs: specs.map(runView) });
  const hsps = new Map(specs.map((spec) => [spec.runId, spec.hsps ?? SPECS]));
  const tray = trayDeps(runs);
  for (const spec of specs) tray.records.set(spec.runId, outputs(spec.hsps ?? SPECS, spec.subjects ?? SUBJECTS).records);
  return { runs, hsps, tray, candidates: new CandidateTray(tray.deps) };
}

const keys = (tray: CandidateTray) => tray.state.get().candidates.map((candidate) => candidate.key);

// --- the results browser's part ---------------------------------------------------------------------

describe('what the results browser gives the tray', () => {
  it('makes candidate sources from the HSP table, the run snapshot and the outfmt 6 rows, without reading', async () => {
    const { runs, hsps } = setup([{ runId: 'r1', number: 4, title: 'mito', hsps: SPECS }]);
    const results = browser(runs, hsps);
    await results.open('r1');
    const [first, onePosition] = results.candidateSources([
      { runId: 'r1', qIdx: 0, rank: 1 },
      { runId: 'r1', qIdx: 0, rank: 3 },
    ]);
    expect(first).toEqual({
      id: { runId: 'r1', qIdx: 0, rank: 1 },
      index: 1,
      run: { runId: 'r1', number: 4, title: 'mito', program: 'blastn' },
      query: { position: 0, id: 'q0', length: 300, kind: 'nucleotide', unit: 'nt' },
      subject: { position: 2, id: 's2', length: 400, kind: 'nucleotide', unit: 'nt' },
      coordinates: { q_start: 70, q_end: 100, s_start: 120, s_end: 90, query_frame: null, subject_frame: null },
      row: { qseqid: 'q0', sseqid: 's2', pident: '99.000', length: '50', mismatch: '0', gapopen: '0', qstart: '70', qend: '100', sstart: '120', send: '90', evalue: '0.00001', bitscore: '40' },
    });
    expect(onePosition!.subject).toEqual({ position: 1, id: '', length: 50, kind: 'nucleotide', unit: 'nt' });
    expect(onePosition!.row.sseqid).toBe('Subject_2');
    expect(() => results.candidateSources([{ runId: 'r1', qIdx: 0, rank: 9 }])).toThrow('run 4 has no HSP 1.10');
    expect(() => results.candidateSources([{ runId: 'r2', qIdx: 0, rank: 0 }])).toThrow('(it is of another run)');
  });

  it('keeps the frames of translated records', async () => {
    const specs: HspSpec[] = [{ q: 0, s: 0, coords: [30, 1, 10, 19], frames: [-2, 1] }];
    const { runs, hsps } = setup([{ runId: 'x', number: 1, program: 'tblastx', hsps: specs }]);
    const [source] = await sourcesOf(runs, hsps, 'x');
    expect(source!.coordinates).toMatchObject({ query_frame: -2, subject_frame: 1 });
    expect(source!.query.unit).toBe('nt');
  });

  it('lists every HSP of a subject in the engine order, whatever the filters show', async () => {
    const { runs, hsps } = setup();
    const results = browser(runs, hsps);
    await results.open('r1');
    results.setFilters({ maxEValue: 1e-20 });
    expect(results.state.get().hsps.map((hsp) => hsp.id.rank)).toEqual([0]);
    expect(results.hspIdsOfSubject(0, 2)).toEqual([
      { runId: 'r1', qIdx: 0, rank: 0 },
      { runId: 'r1', qIdx: 0, rank: 1 },
    ]);
    expect(results.hspIdsOfSubject(1, 2)).toEqual([]);
    expect(results.hspIdsOfSubject(5, 0)).toEqual([]);
  });

  it('marks subjects that the list shows, keeps them through sorting, and drops them with the query, the run or a filter', async () => {
    const { runs, hsps } = setup([
      { runId: 'r1', number: 1 },
      { runId: 'r2', number: 2 },
    ]);
    const results = browser(runs, hsps);
    await results.open('r1');
    expect(results.state.get().subjects.map((s) => s.sIdx)).toEqual([2, 0, 1]);
    results.markSubjects([0, 2, 7], true);
    expect([...results.state.get().marked].sort()).toEqual([0, 2]);
    expect(results.markedHspIds().map((id) => id.rank)).toEqual([0, 1, 2]);
    results.setSubjectSort({ key: 'bitScore', descending: false });
    expect([...results.state.get().marked].sort()).toEqual([0, 2]);
    // The list's order: subject 0 (bits 70) before subject 2 (bits 90) ascending.
    expect(results.markedHspIds().map((id) => id.rank)).toEqual([2, 0, 1]);
    results.selectSubject(0);
    expect(results.state.get().marked.size).toBe(2);
    // A filter that hides subject 0 drops its mark.
    results.setFilters({ subjectText: 's2' });
    expect([...results.state.get().marked]).toEqual([2]);
    results.setFilters({});
    results.markAll(true);
    expect([...results.state.get().marked].sort()).toEqual([0, 1, 2]);
    results.markAll(false);
    expect(results.state.get().marked.size).toBe(0);
    results.markAll(true);
    results.selectQuery(1);
    expect(results.state.get().marked.size).toBe(0);
    results.markAll(true);
    expect(results.markedHspIds()).toEqual([{ runId: 'r1', qIdx: 1, rank: 0 }]);
    // An HSP of another query moves the query: its marks go.
    results.selectHsp({ runId: 'r1', qIdx: 0, rank: 2 });
    expect(results.state.get().marked.size).toBe(0);
    results.markAll(true);
    await results.open('r2');
    expect(results.state.get().marked.size).toBe(0);
    expect(results.markedHspIds()).toEqual([]);
  });

  it('goes back to a candidate: opens its run, selects it, and clears only the filters that hide it', async () => {
    const { runs, hsps } = setup([
      { runId: 'r1', number: 1 },
      { runId: 'r2', number: 2 },
    ]);
    const results = browser(runs, hsps);
    await results.open('r1');
    expect(await results.reveal({ runId: 'r2', qIdx: 0, rank: 2 })).toBe(true);
    let state = results.state.get();
    expect([state.runId, state.qIdx, state.sIdx, state.hsp]).toEqual(['r2', 0, 0, { runId: 'r2', qIdx: 0, rank: 2 }]);
    expect(state.revealed).toEqual({ id: { runId: 'r2', qIdx: 0, rank: 2 }, cleared: [] });

    // Rank 3 (E value 1e-2, bits 20, subject without an ID) is hidden by three of the four filters.
    results.setFilters({ maxEValue: 1e-3, minBitScore: 10, subjectText: 'subject_', queryText: 'q1' });
    expect(results.state.get().qIdx).toBe(1);
    expect(await results.reveal({ runId: 'r2', qIdx: 0, rank: 3 })).toBe(true);
    state = results.state.get();
    expect(state.filters).toEqual({ minBitScore: 10, subjectText: 'subject_' });
    expect([state.qIdx, state.sIdx, state.hsp]).toEqual([0, 1, { runId: 'r2', qIdx: 0, rank: 3 }]);
    expect(state.hsps.map((hsp) => hsp.id.rank)).toEqual([3]);
    expect(state.revealed).toEqual({
      id: { runId: 'r2', qIdx: 0, rank: 3 },
      cleared: ['queryText', 'maxEValue'],
      message: 'The view filters hid HSP 1.4 of run 2, so these were cleared: Find queries by ID, E value.',
    });

    expect(await results.reveal({ runId: 'r2', qIdx: 1, rank: 4 })).toBe(false);
    expect(results.state.get().revealed).toEqual({ id: { runId: 'r2', qIdx: 1, rank: 4 }, cleared: [], message: 'Run 2 has no HSP 2.5.' });
    expect(results.state.get().hsp).toEqual({ runId: 'r2', qIdx: 0, rank: 3 });
  });

  it('does not go back to a run without results', async () => {
    const { runs, hsps } = setup([{ runId: 'r1', number: 1, status: 'cancelled' }]);
    const results = browser(runs, hsps);
    expect(await results.reveal({ runId: 'r1', qIdx: 0, rank: 0 })).toBe(false);
    expect(results.state.get()).toMatchObject({ phase: 'unavailable', message: 'Run 1 was cancelled. A cancelled run keeps no results.' });
  });

  it('names the filters that hide an HSP as filterQuery and the query list apply them', () => {
    const hsp = { queryId: 'Contig_7', subjectIds: ['', 'Subject_3'], eValue: 1e-5, bitScore: 42 };
    expect(filtersHiding({}, hsp)).toEqual([]);
    expect(filtersHiding({ maxEValue: 1e-5, minBitScore: 42, subjectText: 'SUBJECT', queryText: 'contig', queriesWithHitsOnly: true }, hsp)).toEqual([]);
    expect(filtersHiding({ maxEValue: 1e-6, minBitScore: 42.5, subjectText: 'x', queryText: 'y' }, hsp)).toEqual(['queryText', 'subjectText', 'maxEValue', 'minBitScore']);
    expect(filtersHiding({ maxEValue: 1 }, { ...hsp, eValue: Number.NaN })).toEqual(['maxEValue']);
  });
});

// --- the tray -------------------------------------------------------------------------------------------

describe('CandidateTray', () => {
  it('adds HSPs once each, in order, selected, with an empty note and when they were added', async () => {
    const { runs, hsps, candidates } = setup();
    const sources = await sourcesOf(runs, hsps, 'r1');
    expect(candidates.add(sources.slice(0, 3))).toEqual({ ok: true, added: 3, already: 0 });
    expect(candidates.add([sources[1]!, sources[4]!, sources[4]!])).toEqual({ ok: true, added: 1, already: 2 });
    const state = candidates.state.get();
    expect(keys(candidates)).toEqual(['r1/0/0', 'r1/0/1', 'r1/0/2', 'r1/1/0']);
    expect(state.candidates.map((c) => [c.note, c.addedAt, c.serial])).toEqual([
      ['', 5000, 0],
      ['', 5000, 1],
      ['', 5000, 2],
      ['', 5001, 3],
    ]);
    expect([...state.selected]).toEqual(keys(candidates));
    expect(candidateKey({ runId: 'r1', qIdx: 1, rank: 0 })).toBe('r1/1/0');
    expect(state.candidates[3]!.row.evalue).toBe('0.0001');
  });

  it.each<[RunStatus | 'unknown', string]>([
    ['queued', 'Run 2 has not completed yet (it is queued).'],
    ['preparing', 'Run 2 has not completed yet (it is preparing).'],
    ['running', 'Run 2 has not completed yet (it is running).'],
    ['finalizing', 'Run 2 has not completed yet (it is finalizing).'],
    ['cancelled', 'Run 2 was cancelled, and a cancelled run keeps no results to take candidates from.'],
    ['failed', 'Run 2 failed, and a failed run keeps no results to take candidates from.'],
    ['unknown', 'Run 2 is not in this working session.'],
  ])('refuses HSPs, extraction and export of a run that is %s (REQ-10)', async (status, why) => {
    const { runs, hsps, candidates, tray } = setup([
      { runId: 'r1', number: 1 },
      { runId: 'r2', number: 2 },
    ]);
    const good = await sourcesOf(runs, hsps, 'r1');
    const other = await sourcesOf(runs, hsps, 'r2');
    const withStatus = (list: readonly RunView[]) =>
      status === 'unknown' ? list.filter((run) => run.snapshot.runId !== 'r2') : list.map((run) => (run.snapshot.runId === 'r2' ? { ...run, status } : run));
    candidates.add(good.slice(0, 2));
    runs.set({ runs: withStatus(runs.get().runs) });
    const before = candidates.state.get().candidates;
    expect(candidates.add([good[2]!, other[0]!])).toEqual({ ok: false, message: `Only HSPs of completed runs can be added to the candidates. ${why}` });
    expect(candidates.state.get().candidates).toBe(before);
    expect(candidates.state.get().message).toContain(why);

    // Candidates taken while the run was completed (a state that the coordinator never goes back
    // from, made here to test the refusal) are refused too.
    runs.set({ runs: runs.get().runs.map((run): RunView => ({ ...run, status: 'completed' })).concat(status === 'unknown' ? [runView({ runId: 'r2', number: 2 })] : []) });
    candidates.add([other[0]!]);
    runs.set({ runs: withStatus(runs.get().runs) });
    expect(await candidates.extract()).toEqual({ ok: false, message: `Only candidates of completed runs can be extracted. ${why}` });
    expect(await candidates.exportAlignments()).toEqual({ ok: false, message: `Only candidates of completed runs can be exported. ${why}` });
    expect(tray.saved).toEqual([]);
    expect(tray.residueCalls).toEqual([]);
    // Without the candidate of that run, the others are written.
    candidates.select([candidateKey(other[0]!.id)], false);
    expect((await candidates.extract()).ok).toBe(true);
    expect(candidates.state.get().message).toBeUndefined();
  });

  it('removes, moves, notes, selects and clears', async () => {
    const { runs, hsps, candidates } = setup();
    candidates.add(await sourcesOf(runs, hsps, 'r1'));
    candidates.move('r1/0/0', 3);
    expect(keys(candidates)).toEqual(['r1/0/1', 'r1/0/2', 'r1/0/3', 'r1/0/0', 'r1/1/0']);
    candidates.move('r1/1/0', -5);
    candidates.move('r1/0/2', 99);
    candidates.move('nothing', 0);
    expect(keys(candidates)).toEqual(['r1/1/0', 'r1/0/1', 'r1/0/3', 'r1/0/0', 'r1/0/2']);
    candidates.setNote('r1/0/3', 'check the repeat');
    candidates.setNote('nothing', 'x');
    expect(candidates.state.get().candidates.map((c) => c.note)).toEqual(['', '', 'check the repeat', '', '']);
    candidates.select(['r1/0/1', 'r1/0/3', 'nothing'], false);
    expect([...candidates.state.get().selected].sort()).toEqual(['r1/0/0', 'r1/0/2', 'r1/1/0']);
    candidates.remove(['r1/0/0', 'r1/0/1']);
    expect(keys(candidates)).toEqual(['r1/1/0', 'r1/0/3', 'r1/0/2']);
    expect([...candidates.state.get().selected].sort()).toEqual(['r1/0/2', 'r1/1/0']);
    candidates.selectAll(false);
    expect(candidates.state.get().selected.size).toBe(0);
    candidates.select(['nothing', 'r1/0/3'], true);
    expect([...candidates.state.get().selected]).toEqual(['r1/0/3']);
    candidates.selectAll(true);
    expect(candidates.state.get().selected.size).toBe(3);
    // An HSP removed and added again comes back at the end, without its note.
    candidates.remove(['r1/0/3']);
    candidates.add(await sourcesOf(runs, hsps, 'r1', [{ runId: 'r1', qIdx: 0, rank: 3 }]));
    expect(keys(candidates)).toEqual(['r1/1/0', 'r1/0/2', 'r1/0/3']);
    expect(candidates.state.get().candidates[2]!.note).toBe('');
    candidates.clear();
    expect(candidates.state.get()).toEqual({ candidates: [], selected: new Set() });
  });

  it('sorts by addition, by run and engine order, and by subject and position; sorting rewrites the order', async () => {
    const { runs, hsps, candidates } = setup([
      { runId: 'r1', number: 1 },
      { runId: 'r2', number: 2 },
    ]);
    const one = await sourcesOf(runs, hsps, 'r1');
    const two = await sourcesOf(runs, hsps, 'r2');
    candidates.add([two[4]!, one[3]!, two[0]!]);
    candidates.add([one[0]!, two[2]!, one[1]!]);
    const added = ['r2/1/0', 'r1/0/3', 'r2/0/0', 'r1/0/0', 'r2/0/2', 'r1/0/1'];
    expect(keys(candidates)).toEqual(added);
    candidates.sortBy('run');
    expect(keys(candidates)).toEqual(['r1/0/0', 'r1/0/1', 'r1/0/3', 'r2/0/0', 'r2/0/2', 'r2/1/0']);
    // Subject IDs: "Subject_2" (no ID) < "s0" < "s2"; on s0, 1-41 before 80-100; on s2, 90-120 before 301-360.
    candidates.sortBy('subject');
    expect(keys(candidates)).toEqual(['r1/0/3', 'r2/0/2', 'r2/1/0', 'r1/0/1', 'r1/0/0', 'r2/0/0']);
    candidates.move('r2/0/0', 0);
    expect(keys(candidates)[0]).toBe('r2/0/0');
    candidates.sortBy('added');
    expect(keys(candidates)).toEqual(added);
  });

  it('lists where the candidates come from, one entry per run, in run order', async () => {
    const { runs, hsps, candidates } = setup([
      { runId: 'r1', number: 1, title: 'first' },
      { runId: 'r2', number: 2 },
      { runId: 'r3', number: 3 },
    ]);
    candidates.add((await sourcesOf(runs, hsps, 'r2')).slice(0, 1));
    candidates.add((await sourcesOf(runs, hsps, 'r1')).slice(0, 3));
    expect(candidates.origins()).toEqual([
      {
        runId: 'r1',
        number: 1,
        title: 'first',
        program: 'blastn',
        options: ['-evalue', '1e-5'],
        query: { name: 'q.fa', sha256: 'sha-r1-q.fa' },
        subject: { name: 's.fa', sha256: 'sha-r1-s.fa' },
        engineBuild: 'build-7',
        endedAt: 1001,
        candidates: 3,
      },
      expect.objectContaining({ runId: 'r2', number: 2, candidates: 1 }),
    ]);
    expect(candidates.origins()[1]).not.toHaveProperty('title');
  });

  it('extracts the selected candidates in tray order: one read per record, clipped ends and unknown strands reported', async () => {
    const { runs, hsps, candidates, tray } = setup();
    candidates.add(await sourcesOf(runs, hsps, 'r1'));
    candidates.select(['r1/1/0'], false);
    candidates.move('r1/0/2', 0);
    let busy: unknown;
    const unsubscribe = candidates.state.subscribe((state) => {
      if (state.busy !== undefined) busy = state.busy;
    });
    const result = await candidates.extract({ region: { kind: 'flanked', flanks: { left: 5, right: 20 } } });
    unsubscribe();
    expect(busy).toBe('sequences');
    expect(candidates.state.get().busy).toBeUndefined();
    expect(tray.residueCalls.map((call) => [call.revisionIds, call.position, call.intervals])).toEqual([
      [['r1:s.fa'], 0, [{ from: 1, to: 61 }]],
      [['r1:s.fa'], 2, [{ from: 296, to: 380 }, { from: 85, to: 140 }]],
      [['r1:s.fa'], 1, [{ from: 4, to: 29 }]],
    ]);
    expect(tray.saved.map(({ name, mime }) => [name, mime])).toEqual([[SEQUENCES_FILE, FASTA_MIME]]);
    const text = tray.saved[0]!.text;
    const s0 = letters(100, 0);
    expect(text.split('\n').filter((line) => line.startsWith('>'))).toEqual([
      '>s0:1-61 run=1 subject_record=1 length=100 unit=nt hsps=1.3 hit_strand=plus requested=-4-61',
      '>s2:296-380 run=1 subject_record=3 length=400 unit=nt hsps=1.1 hit_strand=plus',
      '>s2:85-140 run=1 subject_record=3 length=400 unit=nt hsps=1.2 hit_strand=minus',
      '>Subject_2:4-29 run=1 subject_record=2 length=50 unit=nt hsps=1.4 hit_strand=unknown',
    ]);
    expect(text.startsWith(`>s0:1-61 run=1 subject_record=1 length=100 unit=nt hsps=1.3 hit_strand=plus requested=-4-61\n${s0.slice(0, 60)}\n${s0[60]}\n`)).toBe(true);
    expect(result).toEqual({ ok: true, summary: candidates.state.get().last });
    expect(candidates.state.get().last).toEqual({
      output: 'sequences',
      fileName: SEQUENCES_FILE,
      role: 'subject',
      candidates: 4,
      sequences: 4,
      bytes: encoder.encode(text).length,
      clipped: [{ name: 's0', runNumber: 1, role: 'subject', position: 0, hsps: ['1.3'], requested: { from: -4, to: 61 }, actual: { from: 1, to: 61 }, recordLength: 100, unit: 'nt' }],
      unknownStrand: [expect.stringContaining('HSP 1.4 of run 1 covers one letter of subject record 2')],
    });

    // The queries, whole and spanning: one sequence per record.
    tray.saved.length = 0;
    const whole = await candidates.extract({ role: 'query', region: { kind: 'whole' }, join: 'spanning' });
    expect(whole.ok && whole.summary.sequences).toBe(1);
    expect(tray.saved[0]!.text.split('\n')[0]).toBe('>q0:1-300 run=1 query_record=1 length=300 unit=nt hsps=1.3,1.1,1.2,1.4 hit_strand=mixed');
  });

  it('saves nothing and leaves the tray as it was when a read fails, and refuses without a selection or while busy', async () => {
    const runs = new Store<AppState>({ runs: [runView({ runId: 'r1', number: 1 })] });
    const hsps = new Map([['r1', SPECS]]);
    const failing = trayDeps(runs, { fail: 'the file is gone' });
    const candidates = new CandidateTray(failing.deps);
    candidates.add(await sourcesOf(runs, hsps, 'r1'));
    candidates.setNote('r1/0/0', 'kept');
    const before = candidates.state.get();
    const result = await candidates.extract();
    expect(result).toEqual({ ok: false, message: 'The sequences could not be extracted: the file is gone' });
    const after = candidates.state.get();
    expect([after.candidates, after.selected, after.busy, after.message]).toEqual([before.candidates, before.selected, undefined, result.ok ? '' : result.message]);
    expect(failing.saved).toEqual([]);

    // A flank that is not a whole number is refused before anything is read.
    failing.residueCalls.length = 0;
    expect((await candidates.extract({ region: { kind: 'flanked', flanks: { left: -1, right: 0 } } })).ok).toBe(false);
    expect(failing.residueCalls).toEqual([]);

    const { candidates: slow, tray } = setup();
    slow.add(await sourcesOf(runs, hsps, 'r1'));
    const first = slow.extract();
    expect(await slow.exportAlignments()).toEqual({ ok: false, message: 'The candidates are being written to a file; wait until that ends.' });
    expect((await first).ok).toBe(true);
    slow.selectAll(false);
    expect(await slow.extract()).toEqual({ ok: false, message: 'Select the candidates to extract.' });
    expect(await slow.exportAlignments()).toEqual({ ok: false, message: 'Select the candidates to export.' });
    expect(tray.saved).toHaveLength(1);
  });

  it('exports the gapped alignments of the selected candidates from their HSP records: one read per run, apart from the sequences', async () => {
    const { runs, hsps, candidates, tray } = setup([
      { runId: 'r1', number: 1 },
      { runId: 'r2', number: 2, hsps: [{ q: 1, s: 2, coords: [50, 41, 200, 209], aligned: ['GGTTAACCAA', 'GGTTAACCTA'] }] },
    ]);
    candidates.add([...(await sourcesOf(runs, hsps, 'r1')), ...(await sourcesOf(runs, hsps, 'r2'))]);
    candidates.move('r2/1/0', 1);
    const result = await candidates.exportAlignments();
    expect(tray.hspCalls).toEqual([
      { runId: 'r1', indices: [0, 1, 2, 3, 4] },
      { runId: 'r2', indices: [0] },
    ]);
    expect(tray.residueCalls).toEqual([]);
    expect(tray.saved.map(({ name, mime }) => [name, mime])).toEqual([[ALIGNMENTS_FILE, FASTA_MIME]]);
    expect(tray.saved[0]!.text).toBe(
      [
        '>q0:1-60 run=1 query_record=1 hsp=1.1 aligned',
        'ACGT-A',
        '>s2:301-360 run=1 subject_record=3 hsp=1.1 aligned',
        'ACGTTA',
        '>q1:50-41 run=2 query_record=2 hsp=2.1 aligned',
        'GGTTAACCAA',
        '>s2:200-209 run=2 subject_record=3 hsp=2.1 aligned',
        'GGTTAACCTA',
        '>q0:70-100 run=1 query_record=1 hsp=1.2 aligned',
        'acg',
        '>s2:120-90 run=1 subject_record=3 hsp=1.2 aligned',
        'ACG',
        '',
      ].join('\n'),
    );
    expect(result.ok && result.summary).toEqual({
      output: 'alignments',
      fileName: ALIGNMENTS_FILE,
      candidates: 6,
      alignments: 3,
      bytes: encoder.encode(tray.saved[0]!.text).length,
      missing: [
        'HSP 1.3 of run 1 has no aligned sequences in its HSP record, so it has no gapped alignment to write.',
        'HSP 1.4 of run 1 has no aligned sequences in its HSP record, so it has no gapped alignment to write.',
        'HSP 2.1 of run 1 has no aligned sequences in its HSP record, so it has no gapped alignment to write.',
      ],
    });

    // A record that is not the candidate's HSP, and candidates without aligned rows, save nothing.
    tray.saved.length = 0;
    tray.records.set('r2', [{ ...tray.records.get('r2')![0]!, rank: 5 }]);
    expect(await candidates.exportAlignments()).toEqual({ ok: false, message: 'The alignments could not be exported: HSP record 0 of run 2 is not HSP 2.1' });
    candidates.selectAll(false);
    candidates.select(['r1/0/2'], true);
    const none = await candidates.exportAlignments();
    expect(none.ok).toBe(false);
    expect(none.ok ? '' : none.message).toContain('no selected candidate has aligned sequences in its HSP record');
    expect(tray.saved).toEqual([]);
  });
});

// --- size ---------------------------------------------------------------------------------------------------

describe('the tray at the size of a large run', () => {
  /** Sources of `n` HSPs on 50 subjects of 3 runs, made without a results browser. */
  function many(n: number): CandidateSource[] {
    const runs = ['ra', 'rb', 'rc'];
    return Array.from({ length: n }, (_, i) => {
      const runId = runs[i % 3]!;
      const subject = (i * 7919) % 50;
      const start = (i * 104_729) % 1_000_000;
      return {
        id: { runId, qIdx: i % 11, rank: i },
        index: i,
        run: { runId, number: (i % 3) + 1, program: 'blastn' },
        query: { position: i % 11, id: `q${i % 11}`, length: 5000, kind: 'nucleotide', unit: 'nt' },
        subject: { position: subject, id: `chr${subject}`, length: 2_000_000, kind: 'nucleotide', unit: 'nt' },
        coordinates: { q_start: 1, q_end: 100, s_start: start + 1, s_end: start + 100, query_frame: null, subject_frame: null },
        row: { qseqid: 'q', sseqid: 's', pident: '100.000', length: '100', mismatch: '0', gapopen: '0', qstart: '1', qend: '100', sstart: '1', send: '100', evalue: '1e-50', bitscore: '185' },
      };
    });
  }

  it('takes every HSP of a subject of 3000 copies, and of 200 marked subjects, from the results browser', async () => {
    const subjects = Array.from({ length: 201 }, (_, i) => [`chr${i}`, 1_000_000] as const);
    const specs: HspSpec[] = [];
    for (let i = 0; i < 3000; i++) specs.push({ q: 0, s: 0, coords: [1, 500, 1 + i * 300, 500 + i * 300], bits: 900 });
    for (let s = 1; s <= 200; s++) for (let i = 0; i < 5; i++) specs.push({ q: 0, s, coords: [1, 100, 1 + i * 1000, 100 + i * 1000], bits: 100 });
    const { runs, hsps, candidates } = setup([{ runId: 'big', number: 1, subjects, hsps: specs }]);
    const results = browser(runs, hsps);
    await results.open('big');
    let started = performance.now();
    const copies = results.candidateSources(results.hspIdsOfSubject(0, 0));
    expect(candidates.add(copies)).toEqual({ ok: true, added: 3000, already: 0 });
    const oneSubject = performance.now() - started;
    results.markAll(true);
    results.markSubjects([0], false);
    started = performance.now();
    const marked = results.candidateSources(results.markedHspIds());
    expect(candidates.add(marked)).toEqual({ ok: true, added: 1000, already: 0 });
    const markedSubjects = performance.now() - started;
    expect(candidates.state.get().candidates).toHaveLength(4000);
    expect(new Set(marked.map((source) => source.subject.position)).size).toBe(200);
    // Both take a few milliseconds; the bounds leave room for a busy machine.
    expect(oneSubject).toBeLessThan(500);
    expect(markedSubjects).toBeLessThan(500);
  });

  function exercise(n: number): number {
    const runs = new Store<AppState>({ runs: ['ra', 'rb', 'rc'].map((runId, i) => runView({ runId, number: i + 1 })) });
    const tray = new CandidateTray(trayDeps(runs).deps);
    const sources = many(n);
    const started = performance.now();
    tray.add(sources);
    tray.add(sources.slice(0, n / 2));
    tray.sortBy('subject');
    tray.sortBy('run');
    for (let i = 0; i < 20; i++) tray.move(tray.state.get().candidates[(i * 997) % n]!.key, (i * 7) % n);
    tray.setNote(tray.state.get().candidates[n - 1]!.key, 'last');
    tray.select(sources.slice(0, n / 4).map((source) => candidateKey(source.id)), false);
    tray.selectAll(true);
    tray.origins();
    tray.sortBy('added');
    tray.remove(sources.filter((_, i) => i % 2 === 0).map((source) => candidateKey(source.id)));
    const elapsed = performance.now() - started;
    expect(tray.state.get().candidates).toHaveLength(n / 2);
    return elapsed;
  }

  it('adds 20 000 candidates, and sorts, moves, notes, selects and removes them in time linear in their number', () => {
    const best = (n: number) => Math.min(exercise(n), exercise(n), exercise(n));
    const small = best(5000);
    const large = best(20_000);
    // Four times the candidates: linear work takes about 4 times as long (n log n a little more),
    // quadratic work 16 times.
    expect(large / Math.max(small, 1)).toBeLessThan(9);
    expect(large).toBeLessThan(2000);
  });
});
