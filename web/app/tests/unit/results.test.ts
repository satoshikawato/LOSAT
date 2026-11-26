// The results screen's domain and application parts (S13): the HSP table and its index,
// sorting and filtering with engine values, the limits and the orientation, the
// verification badge, and the ResultsBrowser's selection by HSP identity.
import { describe, expect, it } from 'vitest';
import type { MessagePlace } from 'losat/domain';
import type { HspRecord, ProgramDescription } from 'losat/ports';
import type { AppState, RunView, ToolRunView } from '../../src/application/coordinator';
import { keptOutputs, noHspReason, ResultsBrowser, type ResultsDeps } from '../../src/application/results';
import { Store } from '../../src/application/store';
import { hspTable, out0Range, out6Range, frame, subjectSpan } from '../../src/domain/hsp-table';
import { headingTitle } from '../../src/domain/outfmt0';
import { shownScores, splitOutfmt6Row } from '../../src/domain/outfmt6';
import {
  buildResultIndex,
  filterQuery,
  flagGiven,
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
  toolOptionGrammar,
  toolOptionKey,
  verificationBadge,
  RANGE_VALUE,
  type BadgeInput,
  type VerificationTable,
} from '../../src/domain/verification';
import type { RunRecord, RunSnapshot } from '../../src/domain/run';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import type { FailedOutput } from '../../src/ports/data';
import { NO_HSP_TEXT } from '../../src/ui/no-hsp-text';
import { NO_WARNINGS, placeParts, placeText, warningsPlaceParts, warningsPlaceText } from '../../src/ui/message-place-text';

// --- a small run: two queries (the second without hits), three subjects ----------------------

interface Spec {
  readonly q: number;
  readonly s: number;
  readonly round?: number;
  readonly bits: number;
  readonly e: number;
  readonly coords: readonly [number, number, number, number];
  readonly shown?: boolean;
  readonly frames?: readonly [number | null, number | null];
  /** The outfmt 6 texts of the E value and the bit score, where they are not the numbers as JavaScript writes them. */
  readonly text?: { readonly e?: string; readonly bits?: string };
  /** The outfmt 0 text of the subject's heading and of the HSP's section, where they are not the plain ones. */
  readonly out0?: { readonly heading?: string; readonly section?: string };
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
  const ranks = new Map<string, number>();
  const headings = new Map<string, [number, number]>();
  specs.forEach((spec, index) => {
    const [qs, qe, ss, se] = spec.coords;
    const row = [`q${spec.q}`, `s${spec.s}`, '99.000', '50', '0', '0', qs, qe, ss, se, spec.text?.e ?? String(spec.e), spec.text?.bits ?? String(spec.bits)].join('\t') + '\n';
    const out6Start = out6.length;
    out6 += row;
    const shown = spec.shown ?? true;
    const round = spec.round ?? 0;
    const key = `${spec.q}:${round}:${spec.s}`;
    if (shown && !headings.has(key)) {
      const heading = spec.out0?.heading ?? `> s${spec.s} subject number ${spec.s}\nLength=500\n\n`;
      headings.set(key, [out0.length, out0.length + heading.length]);
      out0 += heading;
    }
    const sectionStart = out0.length;
    if (shown) out0 += spec.out0?.section ?? ` Score = ${spec.bits} bits\n\nQuery  ${qs}  ACGT  ${qe}\n`;
    const rankKey = `${spec.q}:${round}`;
    const rank = ranks.get(rankKey) ?? 0;
    ranks.set(rankKey, rank + 1);
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
      ...(spec.round === undefined ? {} : { round: spec.round }),
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

  it('reads the E value and bit score that each row shows, and keeps the record value where there is no row', () => {
    const rows = ['q\ts\t100.000\t9\t0\t0\t1\t9\t1\t9\t0.0\t400\n', 'q\ts\t99.0\t9\t0\t0\t1\t9\t1\t9\t2e-50\t 45.4\r\n'];
    const out6 = new TextEncoder().encode(rows.join(''));
    const record = (index: number, e: number, bits: number, range: readonly [number, number] | null): HspRecord => ({
      ...{ index, q_idx: 0, s_idx: 0, rank: index, raw_score: 1, bit_score: bits, e_value: e },
      ...{ q_start: 1, q_end: 9, s_start: 1, s_end: 9, query_frame: null, subject_frame: null, subject_length: 9 },
      ...{ query_aligned: null, subject_aligned: null, out6: range, out0: null, out0_subject: null },
    });
    const table = hspTable([
      record(0, 3e-200, 400.4, [0, rows[0]!.length]),
      // A range that ends past the stored output reads what there is.
      record(1, 2.4e-50, 45.42, [rows[0]!.length, out6.length + 10]),
      record(2, 0.5, 20.1, null),
    ]);
    const scores = shownScores(table, out6);
    expect([...scores.eValue]).toEqual([0, 2e-50, 0.5]);
    expect([...scores.bitScore]).toEqual([400, 45.4, 20.1]);
    // A field that is empty or not a number (the FakeEngine writes FAKE) keeps the record value.
    const odd = new TextEncoder().encode('q\ts\t1\t9\t0\t0\t1\t9\t1\t9\tFAKE\t\n');
    const kept = shownScores(hspTable([record(0, 1e-5, 33.3, [0, odd.length])]), odd);
    expect([...kept.eValue, ...kept.bitScore]).toEqual([1e-5, 33.3]);
    // The record's values stay as they are.
    expect([...table.eValue]).toEqual([3e-200, 2.4e-50, 0.5]);
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
  // The hit list options have no default in NCBI's argument table (BLAST+ tells them given from not
  // given): describe's `defaults` give the values that the engine resolves, by task.
  const resolved = (value: string) => ({ resolvedDefault: { byTask: { blastp: value, 'blastp-short': value } } });
  const described = [
    { flag: '-task', defaultValue: 'blastp' },
    { flag: '-max_target_seqs', ...resolved('500') },
    { flag: '-max_hsps' },
  ];

  it('reads the limits of the argv, or the default that the engine describes', () => {
    expect(hitLimits(['blastp', '-query', 'q', '-subject', 's'], described)).toEqual({ maxTargetSeqs: 500 });
    expect(hitLimits(['blastp', '-max_target_seqs', '5', '-max_hsps', '2'], described)).toEqual({ maxTargetSeqs: 5, maxHsps: 2 });
    expect(hitLimits(['blastp'], [{ flag: '-max_target_seqs', defaultValue: '100' }])).toEqual({ maxTargetSeqs: 100 });
    // The run's task chooses the resolved value: its -task, else the default task; '' for a program without -task.
    const byTask = [{ flag: '-task', defaultValue: 'blastp' }, { flag: '-max_target_seqs', resolvedDefault: { byTask: { blastp: '500', 'blastp-short': '7' } } }];
    expect(hitLimits(['blastp', '-task', 'blastp-short'], byTask)).toEqual({ maxTargetSeqs: 7 });
    expect(hitLimits(['blastp'], byTask)).toEqual({ maxTargetSeqs: 500 });
    expect(hitLimits(['tblastx'], [{ flag: '-max_target_seqs', resolvedDefault: { byTask: { '': '500' } } }])).toEqual({ maxTargetSeqs: 500 });
    // -task=blastp-short is the engine's blastp-short (runTask), and the engine's value for the task comes before
    // NCBI's default, as in the search form (taskDefault).
    expect(hitLimits(['blastp', '-task=blastp-short'], byTask)).toEqual({ maxTargetSeqs: 7 });
    const both = [{ flag: '-task', defaultValue: 'blastp' }, { flag: '-max_target_seqs', defaultValue: '500', resolvedDefault: { byTask: { 'blastp-short': '7' } } }];
    expect(hitLimits(['blastp', '-task', 'blastp-short'], both)).toEqual({ maxTargetSeqs: 7 });
    expect(hitLimits(['blastp'], both)).toEqual({ maxTargetSeqs: 500 });
  });

  // The tabular formats keep -max_target_seqs, else -num_alignments, else the default of -max_target_seqs;
  // outfmt 0 the larger of -num_descriptions (500) and -num_alignments (250) where either is given
  // (BLAST+ blast_args.cpp:2894-2973), and the run searches twice when the sizes differ (abi_v2.md §7).
  it('reads the hit list sizes of -num_descriptions and -num_alignments', () => {
    const sized = [
      ...described,
      { flag: '-num_descriptions', ...resolved('500') },
      { flag: '-num_alignments', ...resolved('250') },
    ];
    const limits = (...words: string[]) => hitLimits(['blastp', '-query', 'q', '-subject', 's', ...words], sized);
    // -num_alignments is the tabular formats' hit list; outfmt 0 keeps 500 (the larger of 500 and 3).
    expect(limits('-num_alignments', '3')).toEqual({
      maxTargetSeqs: 3,
      subjectFlag: '-num_alignments',
      outfmt0: { alignments: 3, hitList: 500, separate: true },
    });
    // -num_descriptions 3 alone: outfmt 0 keeps 250 and the tabular formats 500.
    expect(limits('-num_descriptions', '3')).toEqual({ maxTargetSeqs: 500, outfmt0: { alignments: 250, hitList: 250, separate: true } });
    // The numbers that make the sizes equal are one search.
    expect(limits('-num_descriptions', '2', '-num_alignments', '4')).toEqual({
      maxTargetSeqs: 4,
      subjectFlag: '-num_alignments',
      outfmt0: { alignments: 4, hitList: 4, separate: false },
    });
    expect(limits('-num_descriptions', '500').outfmt0).toEqual({ alignments: 250, hitList: 500, separate: false });
    expect(limits('-num_alignments', '600').outfmt0).toEqual({ alignments: 600, hitList: 600, separate: false });
    expect(limits('-num_descriptions', '4', '-num_alignments', '2').outfmt0).toEqual({ alignments: 2, hitList: 4, separate: true });
    // The adapter accepts `-flag=value` as the CLI does: the same limits as the two words, and flagGiven sees it.
    expect(limits('-num_alignments=3')).toEqual(limits('-num_alignments', '3'));
    expect(limits('-num_descriptions=4', '-num_alignments=2').outfmt0).toEqual({ alignments: 2, hitList: 4, separate: true });
    expect(limits('-max_target_seqs=7')).toEqual({ maxTargetSeqs: 7 });
    // -sorthits (any value, also 0) orders outfmt 0's hit list by its key, in either spelling.
    expect(limits('-sorthits', '3')).toEqual({ maxTargetSeqs: 500, sorted: true });
    expect(limits('-sorthits=0', '-num_alignments', '3').sorted).toBe(true);
    expect(limits('-sorthsps', '1').sorted).toBeUndefined();
    expect(flagGiven(['blastp', '-num_alignments=3'], '-num_alignments')).toBe(true);
    expect(flagGiven(['blastp', '-num_alignments=3'], '-num_descriptions')).toBe(false);
    expect(flagGiven(['blastp', '-num_alignments', '3'], '-num_alignments')).toBe(true);
    // -max_target_seqs excludes both options: the same hit list in every format, no outfmt 0 entry.
    expect(limits('-max_target_seqs', '7')).toEqual({ maxTargetSeqs: 7 });
    // Without the described defaults, outfmt 0's size is not claimed (and a separate search is assumed).
    expect(hitLimits(['blastp', '-num_descriptions', '3'], described).outfmt0).toEqual({ separate: true });
  });

  it('claims no limit that it cannot read', () => {
    // NCBI's help text states a default (Default = `500'), but the limits read no help text.
    expect(hitLimits(['blastp'], [{ flag: '-max_target_seqs' }])).toEqual({});
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

  // Line RG E1w-3b: -logfile names a file of the run's messages, as -out names its output; the
  // outputs do not change with it. -conffile and -dryrun change the run: they stay options, which no
  // record compares (the engine's describe does not list them).
  it('leaves -logfile out of the options, as -out, and keeps -conffile and -dryrun', () => {
    const plain = optionKey(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '1e-5'], grammar);
    expect(optionKey(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '1e-5', '-logfile', 'run.log'], grammar)).toBe(plain);
    expect(optionKey(['blastn', '-logfile=run.log', '-evalue', '1e-5'], grammar)).toBe(plain);
    expect(optionKey(['blastn', '-evalue', '1e-5', '-dryrun'], grammar)).not.toBe(plain);
    expect(optionKey(['blastn', '-evalue', '1e-5', '-conffile', 'settings.ini'], grammar)).toBeUndefined();
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

  it('supports a run inside a profile, and says what that means', () => {
    const badge = verificationBadge(input, table);
    expect(badge.level).toBe('certified');
    expect(badge.label).toBe('Supported');
    expect(badge.details.join(' ')).toMatch(/not a comparison of your inputs with NCBI/);
    // No word of a certification or a verification (DW-36; the maintainer, 2026-11-04): the
    // badge says that the run is supported.
    expect(`${badge.label} ${badge.details.join(' ')}`).not.toMatch(/verif|certif/i);
  });

  it('says why a run is not supported', () => {
    const reasons = (change: Partial<BadgeInput>) => verificationBadge({ ...input, ...change }, table);
    expect(reasons({ argv: ['blastn', '-evalue', '1'] }).details.join(' ')).toMatch(/for outfmt 0 and 6, not for outfmt 7/);
    expect(reasons({ argv: ['blastn', '-evalue', '2'] }).details.join(' ')).toMatch(/were not compared with NCBI BLAST\+ 2\.17\.0/);
    expect(reasons({ threads: 3 }).details.join(' ')).toMatch(/used 3 threads; the browser runtime was checked with 1, 2 and 4 threads/);
    expect(reasons({ program: 'tblastx', argv: ['tblastx'] }).details.join(' ')).toMatch(/browser runtime of tblastx was not checked/);
    expect(reasons({ threads: 3 }).label).toBe('Not yet supported');
    for (const badge of [reasons({ argv: ['blastn', '-evalue', '1'] }), reasons({ argv: ['blastn', '-evalue', '2'] }), reasons({ threads: 3 }), reasons({ runtimePath: 'fake' })]) {
      expect(`${badge.label} ${badge.details.join(' ')}`).not.toMatch(/verif|certif/i);
    }
    expect(reasons({ runtimePath: 'fake' }).level).toBe('development');
  });

  // The records compare the formats that every run writes; the output of the run's own -outfmt
  // with fields is not among them, so such a run is not yet supported, and the reasons say so (W10).
  it("does not support a run with the output of its -outfmt, and names that output", () => {
    const badge = verificationBadge({ ...input, argv: [...input.argv, '-outfmt', '6 qseqid sseqid'], custom: '6 qseqid sseqid' }, table);
    expect(badge.level).toBe('outside');
    expect(badge.label).toBe('Not yet supported');
    expect(badge.details).toEqual([
      'These options were compared with NCBI BLAST+ 2.17.0 for outfmt 0, 6 and 7, not for outfmt 6 qseqid sseqid.',
      'The browser runtime was not checked for outfmt 6 qseqid sseqid.',
    ]);
    // -outfmt itself is not an option of the profile: the same options without it are supported, in either form.
    expect(verificationBadge({ ...input, argv: [...input.argv, '-outfmt', '7'] }, table).level).toBe('certified');
    expect(verificationBadge({ ...input, argv: [...input.argv, '-outfmt=7'] }, table).level).toBe('certified');
  });

  it('supports a thread count only when it was checked in the browser (WA-RS-03)', () => {
    // Auto picks a count from the user's processors; a count that V-BR did not run is not supported.
    const level = (change: Partial<BadgeInput>) => verificationBadge({ ...input, ...change }, table).level;
    expect(level({ threads: 2 })).toBe('certified');
    expect(level({ runtimePath: 'serial', threads: 1 })).toBe('certified');
    expect(level({ threads: 3 })).toBe('outside');
    expect(level({ threads: 8 })).toBe('outside');
    expect(level({ threads: undefined })).toBe('outside');
  });

  it('says the thread counts checked as ranges, for the path of the run', () => {
    const every = Array.from({ length: 32 }, (_, i) => i + 1);
    const wide: VerificationTable = { ...table, programs: { blastn: { ...table.programs.blastn!, browser: { ...table.programs.blastn!.browser!, threads: every } } } };
    const details = (change: Partial<BadgeInput>) => verificationBadge({ ...input, ...change }, wide).details.join(' ');
    expect(details({ threads: 6 })).toMatch(/threaded runtime in the browser was checked against the native engine with 2 to 32 threads in/);
    expect(details({ runtimePath: 'serial', threads: 1 })).toMatch(/serial runtime in the browser was checked against the native engine with 1 thread in/);
    expect(details({ threads: 40 })).toMatch(/used 40 threads; the browser runtime was checked with 1 to 32 threads\./);
    expect(details({ threads: 4 })).not.toMatch(/\b3, 4\b/);
  });

  it('names the approved exceptions of a non-default subject genetic code', () => {
    expect(approvedExceptions('tblastn', ['tblastn', '-db_gencode', '4'])[0]).toMatch(/PD-TLOSAN-LOCAL-GENCODE-32.*-db_gencode 4/);
    expect(approvedExceptions('tblastx', ['tblastx', '-db_gencode', '11'])).toHaveLength(1);
    expect(approvedExceptions('tblastx', ['tblastx', '-db_gencode', '1'])).toEqual([]);
    expect(approvedExceptions('blastn', ['blastn', '-db_gencode', '4'])).toEqual([]);
    expect(approvedExceptions('tblastn', ['tblastn', '-db_gencode', '4']).join(' ')).not.toMatch(/verif|certif/i);
  });
});

// --- ResultsBrowser ---------------------------------------------------------------------------

const DESCRIPTION: ProgramDescription = {
  program: 'blastn',
  formats: [0, 6, 7],
  parameters: [
    { flag: '-task', help: 'Task', takesValue: true, defaultValue: 'megablast' },
    {
      flag: '-max_target_seqs',
      help: "Maximum number of aligned sequences to keep \n(value of 5 or more is recommended)\nDefault = `500'",
      takesValue: true,
      resolvedDefault: { byTask: { megablast: '500' } },
    },
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

function setup(
  options: {
    status?: RunView['status'];
    specs?: readonly Spec[];
    argv?: readonly string[];
    unreadable?: boolean;
    error?: string;
    snapshot?: Partial<RunSnapshot>;
    record?: Partial<RunRecord>;
    result?: RunView['result'];
  } = {},
) {
  const run = makeRun(options.specs);
  const view: RunView = {
    snapshot: { ...snapshot('r1', options.argv), ...options.snapshot },
    status: options.status ?? 'completed',
    record: { runtimePath: 'serial', threads: 1, messages: 'diagnostics', ...(options.status === 'failed' ? { error: options.error ?? 'boom' } : {}), ...options.record },
    ...(options.result === undefined ? {} : { result: options.result }),
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
  it('uses RunRecord.readRecords for query and subject metadata when an assembled input read differently (D17)', async () => {
    const { results } = setup({
      record: {
        readRecords: {
          query: [
            { id: 'q0', length: 120 },
            { id: 'q1', length: 180 },
            { id: 'q2', length: 90 },
          ],
          subject: [
            { id: 's0', length: 300 },
            { id: 's1', length: 310 },
            { id: 's2', length: 320 },
          ],
        },
      },
    });
    await results.open('r1');
    await settle();
    const state = results.state.get();
    expect(state.queries.map((q) => [q.id, q.length])).toEqual([
      ['q0', 120],
      ['q1', 180],
      ['q2', 90],
    ]);
    expect(state.subjects.map((s) => [s.recordId, s.length])).toEqual([
      ['s1', 310],
      ['s0', 300],
      ['s2', 320],
    ]);
  });

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
    expect(state.hsp).toEqual({ runId: 'r1', qIdx: 0, round: 0, rank: 0 });
    expect(state.hsps.map((h) => [h.id.rank, h.orientation])).toEqual([
      [0, 'forward'],
      [1, 'reverse'],
    ]);
    expect(state.detail).toMatchObject({ state: 'ready', heading: '> s1 subject number 1\nLength=500\n\n' });
    expect(state.detail!.section).toBe(' Score = 90 bits\n\nQuery  1  ACGT  50\n');
    expect(state.detail!.row).toBe('q0\ts1\t99.000\t50\t0\t0\t1\t50\t101\t150\t1e-20\t90\n');
    expect(state.outfmt0Subjects).toBe(2);
    expect(state.atSubjectLimit).toBe(false);
    expect(state.loaded!.badge!.level).toBe('outside');
  });

  it('follows an HSP of another subject, and shows a query without hits', async () => {
    const { results } = setup();
    await results.open('r1');
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 3 });
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
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 1 });
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 2 });
    // The reads for rank 2 answer first: rank 1's detail arrives last, and is dropped.
    release(true);
    await settle();
    hold(false);
    const { detail, hsp } = results.state.get();
    expect(hsp).toEqual({ runId: 'r1', qIdx: 0, round: 0, rank: 2 });
    expect(detail!.id).toEqual(hsp);
    expect(detail!.section).toBe(' Score = 70 bits\n\nQuery  5  ACGT  45\n');
  });

  it('filters by the E value and bit score that the lists show: 0.0 is 0, 123 is 123 (WA-RS-01)', async () => {
    // NCBI writes an E value below 1e-180 as 0.0, and a bit score of 100 or more without decimals.
    const specs: readonly Spec[] = [
      { q: 0, s: 0, bits: 400.4, e: 3e-200, coords: [1, 50, 1, 50], text: { e: '0.0', bits: '400' } },
      { q: 0, s: 1, bits: 122.6, e: 2.4e-50, coords: [1, 50, 1, 50], text: { e: '2e-50', bits: '123' } },
      { q: 0, s: 2, bits: 40.2, e: 1e-5, coords: [1, 50, 1, 50], text: { e: '1e-05', bits: '40.2' } },
    ];
    const { results } = setup({ specs });
    await results.open('r1');
    const shown = () => results.state.get().subjects.map((s) => [s.sIdx, s.first.evalue, s.first.bitscore]);
    results.setFilters({ maxEValue: 0 });
    expect(shown()).toEqual([[0, '0.0', '400']]);
    results.setFilters({ maxEValue: 2e-50 });
    expect(shown()).toEqual([
      [0, '0.0', '400'],
      [1, '2e-50', '123'],
    ]);
    results.setFilters({ minBitScore: 123 });
    expect(shown()).toEqual([
      [0, '0.0', '400'],
      [1, '2e-50', '123'],
    ]);
    results.setFilters({ minBitScore: 400.4 });
    expect(shown()).toEqual([]);
    expect(results.state.get().hidden).toEqual({ subjects: 3, hsps: 3 });
  });

  it('keeps the selected HSP through a filter that shows it, and moves off one that hides it', async () => {
    const { results } = setup();
    await results.open('r1');
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 1 });
    results.setFilters({ maxEValue: 1e-4 });
    expect(results.state.get().hsp).toEqual({ runId: 'r1', qIdx: 0, round: 0, rank: 1 });
    results.setFilters({ maxEValue: 1e-10 });
    expect(results.state.get().hsp).toEqual({ runId: 'r1', qIdx: 0, round: 0, rank: 0 });
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
    expect(state.hsp).toEqual({ runId: 'r1', qIdx: 0, round: 0, rank: 0 });
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

  // line WQA's sweep, D7, 2026-11-19: a run of -parse_deflines names its records by the IDs that its
  // outputs show (the run record's shownIds), in the query list and in the filters.
  it('lists and filters the records of a -parse_deflines run by the IDs that the outputs show', async () => {
    const { results } = setup({ record: { shownIds: { query: ['AB000001.1', 'AB000002.1'], subject: ['NP_1.1', 'NP_2.1', 'NP_3.1'] } } });
    await results.open('r1');
    expect(results.state.get().queries.map((q) => q.id)).toEqual(['AB000001.1', 'AB000002.1']);
    results.setFilters({ queryText: 'ab000002' });
    expect(results.state.get().queries.map((q) => q.id)).toEqual(['AB000002.1']);
    results.setFilters({ queryText: 'AB000001', subjectText: 'np_1.1' });
    expect(results.state.get().subjects.map((s) => s.sIdx)).toEqual([0]);
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
    expect(failed.results.state.get()).toMatchObject({ phase: 'unavailable', message: 'Run 1 failed: boom. It keeps no outputs.' });

    // One full stop at the end of the reason, whatever it ends with (WA-RN-02).
    for (const [error, said] of [
      ['Earlier results are kept.', 'Earlier results are kept.'],
      ['LOSAT engine error: Warning: Sequence contains no data \n', 'LOSAT engine error: Warning: Sequence contains no data.'],
      ['Is it full?', 'Is it full?'],
      ['The file (see the log.)', 'The file (see the log.)'],
      ['', 'no reason was recorded.'],
    ] as const) {
      const run = setup({ status: 'failed', error });
      await run.results.open('r1');
      expect(run.results.state.get().message).toBe(`Run 1 failed: ${said} It keeps no outputs.`);
    }

    // A failed run that made the file of its -logfile keeps it (line RG E1w-3b).
    const logged = setup({ status: 'failed', error: 'LOSAT engine error: The average subject length is too short', record: { runFiles: [{ name: 'run.log', bytes: new Uint8Array([65]) }] } });
    await logged.results.open('r1');
    expect(logged.results.state.get().message).toBe(
      'Run 1 failed: LOSAT engine error: The average subject length is too short. It keeps the file of its -logfile, run.log.',
    );
    // ...and its search strategy file, each under its option (line RG item F1w).
    const both = setup({
      status: 'failed',
      error: 'LOSAT engine error: x',
      record: {
        runFiles: [
          { name: 'out/s.asn', bytes: new Uint8Array(), option: 'export_search_strategy' },
          { name: 'run.log', bytes: new Uint8Array([65]) },
        ],
      },
    });
    await both.results.open('r1');
    expect(both.results.state.get().message).toBe(
      'Run 1 failed: LOSAT engine error: x. It keeps the file of its -export_search_strategy, out/s.asn. It keeps the file of its -logfile, run.log.',
    );
    // A failed run keeps the warnings written before its error (line WQA's sweep, D3, 2026-11-19); white space alone is none.
    const warned = setup({
      status: 'failed',
      error: 'LOSAT engine error: The average subject length is too short',
      record: { diagnostics: 'Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\n' },
    });
    await warned.results.open('r1');
    expect(warned.results.state.get().message).toBe(
      'Run 1 failed: LOSAT engine error: The average subject length is too short. It keeps the warnings written before the error.',
    );
    const blank = setup({ status: 'failed', error: 'boom', record: { diagnostics: '\n' } });
    await blank.results.open('r1');
    expect(blank.results.state.get().message).toBe('Run 1 failed: boom. It keeps no outputs.');

    // A failed run keeps the outputs written before its error (line WQA's sweep, D15, 2026-11-20): said when one holds bytes,
    // with the others' words after it, and not said for an output of no bytes. The outputs that the engine dropped (place
    // `run`) are not the run's.
    const result = (lengths: { 0: number; 6: number; 7: number }, failed: readonly FailedOutput[]) => ({
      runId: 'r1',
      byteLengths: lengths,
      hitCount: 0,
      failed,
    });
    const dropped: readonly FailedOutput[] = [
      { output: 6, error: 'run' },
      { output: 7, error: 'run' },
    ];
    const wrote = setup({ status: 'failed', error: 'boom', result: result({ 0: 282, 6: 0, 7: 0 }, dropped), record: { diagnostics: 'Warning: x\n' } });
    await wrote.results.open('r1');
    expect(wrote.results.state.get().message).toBe('Run 1 failed: boom. It keeps the outputs written before the error. It keeps the warnings written before the error.');
    expect(keptOutputs(wrote.view)).toEqual([0]);
    const empty = setup({ status: 'failed', error: 'boom', result: result({ 0: 0, 6: 0, 7: 0 }, dropped) });
    await empty.results.open('r1');
    expect(empty.results.state.get().message).toBe('Run 1 failed: boom. It keeps no outputs.');
    // The bytes of an output that the engine dropped do not count (they are none, but the rule is the engine's).
    const lost = setup({ status: 'failed', error: 'boom', result: result({ 0: 0, 6: 9, 7: 0 }, dropped) });
    await lost.results.open('r1');
    expect(lost.results.state.get().message).toBe('Run 1 failed: boom. It keeps no outputs.');
    // A completed run keeps every output of its argv.
    expect(keptOutputs({ ...wrote.view, status: 'completed', result: result({ 0: 1, 6: 1, 7: 1 }, []) })).toEqual([0, 6, 7]);

    const running = setup({ status: 'running' });
    await running.results.open('r1');
    expect(running.results.state.get().phase).toBe('unavailable');
    running.runs.set({ runs: [{ ...running.view, status: 'completed' }] });
    await settle();
    await settle();
    expect(running.results.state.get().phase).toBe('ready');
  });

  // The outputs that the run could not write, and where their errors are, are the engine's word (its signal of
  // stream 9, `ResultSetRef.failed`), not a guess from their lengths or the program: outfmt 0 that failed alone in the
  // search of the other formats, or whose search of its own (-num_alignments) failed; outfmt 6 and 7 that failed,
  // which leave the run without HSP records; their errors on stream 3 or in the file of TBLASTX's -logfile, or nowhere
  // (line WQA's sweep, D7, 2026-11-19 and 2026-11-20).
  it('takes the outputs that could not be written, and where their errors are, from the engine', async () => {
    const result = (failed: readonly FailedOutput[], out0 = 0) => ({ runId: 'r1', byteLengths: { 0: out0, 6: 10, 7: 10 }, hitCount: 4, failed });
    const one = ['blastn', '-query', 'q.fa', '-subject', 's.fa'];
    const separate = [...one, '-num_alignments', '3'];
    const failedIn = (error: MessagePlace, ...outputs: readonly (0 | 6 | 7)[]): FailedOutput[] => outputs.map((output) => ({ output, error }));
    const cases: ReadonlyArray<readonly [readonly string[], readonly FailedOutput[], number, MessagePlace | undefined, MessagePlace | undefined]> = [
      [one, failedIn('diagnostics', 0), 0, 'diagnostics', undefined],
      [separate, failedIn('diagnostics', 0), 0, 'diagnostics', undefined],
      [separate, failedIn('logfile', 0), 0, 'logfile', undefined],
      // An empty outfmt 0 that the engine did not name: written, as the engine said.
      [one, [], 0, undefined, undefined],
      [separate, [], 0, undefined, undefined],
      [one, failedIn('diagnostics', 6, 7), 100, undefined, 'diagnostics'],
      [one, failedIn('nowhere', 6, 7), 100, undefined, 'nowhere'],
    ];
    for (const [argv, failed, out0, outfmt0Error, tabularError] of cases) {
      const { results } = setup({ argv, result: result(failed, out0) });
      await results.open('r1');
      const loaded = results.state.get().loaded!;
      expect(loaded.outfmt0Error, JSON.stringify([argv, failed])).toBe(outfmt0Error);
      expect(loaded.tabularError, JSON.stringify([argv, failed])).toBe(tabularError);
    }
    const tabular = setup({ specs: [], result: result(failedIn('logfile', 6, 7), 100) });
    await tabular.results.open('r1');
    expect(noHspReason(tabular.results.state.get())).toBe('no-tabular-output');
    expect(NO_HSP_TEXT['no-tabular-output']).toBe('outfmt 6 and 7 could not be written, so this run has no HSP records.');
    // The words of each place: the run's Warnings, the file of its -logfile under Outputs (by the name of the run's
    // file), or nowhere.
    const run = tabular.results.state.get().loaded!.run;
    const logged = { ...run, record: { ...run.record, runFiles: [{ name: 'run.log', bytes: new Uint8Array(0) }] } };
    expect(placeText('diagnostics', run)).toBe('under Warnings in Run details');
    expect(placeText('logfile', logged)).toBe("in the file of the run's -logfile, run.log, under Outputs");
    expect(placeText('nowhere', run)).toBe("not kept, as the run's -logfile writes its messages nowhere");
    // -num_alignments 0 keeps its own notice: the tabular search fails because of that hit list.
    const zero = setup({
      specs: [],
      argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-num_alignments', '0'],
      result: result(failedIn('diagnostics', 6, 7), 100),
    });
    await zero.results.open('r1');
    expect(noHspReason(zero.results.state.get())).toBe('no-tabular-hits');
  });

  // Where the run's messages went is the engine's word (RunRecord.messages, its signal of stream 10), not a guess from
  // the program or the argv: empty diagnostics mean no warnings where the messages go to them; for TBLASTX with
  // -logfile they are in that file (under Outputs) or nowhere, also when the diagnostics hold lines, and an empty file
  // holds none (line WQA's sweep, D7, 2026-11-20).
  it("takes where the run's messages went from the engine, and says where its warnings are", async () => {
    for (const messages of ['diagnostics', 'logfile', 'nowhere'] as const) {
      const { results } = setup({ record: { messages } });
      await results.open('r1');
      expect(results.state.get().loaded!.messages).toBe(messages);
    }
    const warning = new TextEncoder().encode('Warning: [tblastx] Examining 5 or more matches is recommended\n');
    const { results } = setup({ record: { messages: 'logfile', runFiles: [{ name: 'run.log', bytes: warning }] } });
    await results.open('r1');
    const run = results.state.get().loaded!.run;
    expect(NO_WARNINGS).toBe('The engine wrote no warnings.');
    expect(warningsPlaceText('diagnostics', '', run)).toBe(NO_WARNINGS);
    expect(warningsPlaceText('diagnostics', 'Warning: x\n', run)).toBeUndefined();
    expect(warningsPlaceText('logfile', '', run)).toBe("The run's warnings are in the file of the run's -logfile, run.log, under Outputs.");
    expect(warningsPlaceText('logfile', 'Warning: x\n', run)).toBe("The run's other warnings are in the file of the run's -logfile, run.log, under Outputs.");
    expect(warningsPlaceText('nowhere', '', run)).toBe("The run's warnings are not kept, as the run's -logfile writes its messages nowhere.");
    expect(warningsPlaceText('nowhere', 'Warning: x\n', run)).toBe("The run's other warnings are not kept, as the run's -logfile writes its messages nowhere.");
    // TBLASTX with -logfile run.log on a search without warnings: the file is empty, as the CLI's, and Run details says
    // that the engine wrote none (review 3 L1), or nothing besides the diagnostics where they hold lines.
    const empty = setup({ record: { messages: 'logfile', runFiles: [{ name: 'run.log', bytes: new Uint8Array(0) }] } });
    await empty.results.open('r1');
    const emptyRun = empty.results.state.get().loaded!.run;
    expect(warningsPlaceText('logfile', '', emptyRun)).toBe(NO_WARNINGS);
    expect(warningsPlaceText('logfile', 'Error: x\n', emptyRun)).toBeUndefined();
    // The search strategy file is not the file of -logfile: an empty one beside a log that holds warnings says nothing,
    // and its name is not the log's (line RG item F1w).
    const strategy = { name: 's.asn', bytes: new Uint8Array(0), option: 'export_search_strategy' } as const;
    const withStrategy = setup({ record: { messages: 'logfile', runFiles: [strategy, { name: 'run.log', bytes: warning }] } });
    await withStrategy.results.open('r1');
    const strategyRun = withStrategy.results.state.get().loaded!.run;
    expect(warningsPlaceText('logfile', '', strategyRun)).toBe("The run's warnings are in the file of the run's -logfile, run.log, under Outputs.");
    const strategyOnly = setup({ record: { messages: 'logfile', runFiles: [strategy] } });
    await strategyOnly.results.open('r1');
    expect(warningsPlaceText('logfile', '', strategyOnly.results.state.get().loaded!.run)).toBe(
      "The run's warnings are in the file of the run's -logfile under Outputs.",
    );
    // Files that are not kept say nothing of what they held.
    const unkept = setup({ record: { messages: 'logfile' } });
    await unkept.results.open('r1');
    expect(warningsPlaceText('logfile', '', unkept.results.state.get().loaded!.run)).toBe("The run's warnings are in the file of the run's -logfile under Outputs.");
    // The page shows the option name unbroken (`option-name`): a narrow screen broke "-logfile" after its hyphen.
    for (const place of ['logfile', 'nowhere'] as const) {
      expect(placeParts(place, run).filter((part) => part.option === true).map((part) => part.text)).toEqual(['-logfile']);
      expect(warningsPlaceParts(place, '', run)!.filter((part) => part.option === true).map((part) => part.text)).toEqual(['-logfile']);
    }
    expect(placeParts('diagnostics', run).some((part) => part.option === true)).toBe(false);
    // A completed run whose record does not say where its messages went is not read as one without warnings.
    const unsaid = setup();
    const record: RunRecord = Object.fromEntries(Object.entries(unsaid.view.record).filter(([key]) => key !== 'messages'));
    unsaid.runs.set({ runs: [{ ...unsaid.view, record }] });
    await unsaid.results.open('r1');
    expect(unsaid.results.state.get().phase).toBe('failed');
  });

  // A dry run (-dryrun; line RG E1w-3b) completes without a search: its notice says so in place of "no hits".
  it('says that a completed dry run searched nothing', async () => {
    const { results } = setup({ specs: [], snapshot: { dryRun: true } });
    await results.open('r1');
    expect(results.state.get().phase).toBe('ready');
    expect(noHspReason(results.state.get())).toBe('dry-run');
    expect(NO_HSP_TEXT['dry-run']).toBe('This run was a dry run (-dryrun): LOSAT read the arguments and the settings and ran no search.');
    const searched = setup({ specs: [] });
    await searched.results.open('r1');
    expect(noHspReason(searched.results.state.get())).toBe('no-hits');
  });

  // Line WQA's D14, 2026-11-20: a run of a help or version word completes without a search: its notice says so and where
  // the text is, in place of "no hits", also when the argv has -dryrun (the command line prints the text first).
  it("says that a completed run printed the command line's help or version and searched nothing", async () => {
    for (const prints of ['help', 'version'] as const) {
      const { results } = setup({ specs: [], snapshot: { prints, dryRun: true } });
      await results.open('r1');
      expect(results.state.get().phase).toBe('ready');
      expect(noHspReason(results.state.get())).toBe(`printed-${prints}`);
    }
    expect(NO_HSP_TEXT['printed-help']).toBe("This run printed the command line's help and ran no search. The text is in Outputs.");
    expect(NO_HSP_TEXT['printed-version']).toBe("This run printed the command line's version and ran no search. The text is in Outputs.");
  });

  // Line WQA's D20, S3-1, 2026-11-23: a run that printed the command line's help or version compared no options: no badge.
  it('gives no support badge to a run that printed the command line\'s help or version', async () => {
    for (const prints of ['help', 'version'] as const) {
      const { results } = setup({ specs: [], snapshot: { prints } });
      await results.open('r1');
      expect(results.state.get().phase).toBe('ready');
      expect(results.state.get().loaded!.badge).toBeUndefined();
    }
    const { results } = setup({ specs: [] });
    await results.open('r1');
    expect(results.state.get().loaded!.badge).toBeDefined();
  });

  // WA-FM-07: the command of a run without results is read with the grammar too.
  it("gives the grammar of the run's mode, for its command, also for a run without results", async () => {
    for (const options of [{ status: 'failed' }, { status: 'completed' }, { unreadable: true }] as const) {
      const { results } = setup(options);
      await results.open('r1');
      await settle();
      const grammar = results.state.get().grammar;
      expect(grammar?.takesValue('-max_hsps'), JSON.stringify(options)).toBe(true);
      expect(grammar?.takesValue('-lcase_masking'), JSON.stringify(options)).toBe(false);
    }
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
    expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 0, rank: 0 })).toBe(' Score = 90 bits\n\nQuery  1  ACGT  50\n');
    expect(rangeReads()).toBe(reads);
    expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 0, rank: 1 })).toBe(' Score = 40 bits\n\nQuery  60  ACGT  80\n');
    expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 0, rank: 1 })).toBe(' Score = 40 bits\n\nQuery  60  ACGT  80\n');
    expect(rangeReads()).toBe(reads + 1);
    expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 0, rank: 3 })).toBeUndefined();
    expect(await results.readSection({ runId: 'r2', qIdx: 0, round: 0, rank: 1 })).toBeUndefined();
    expect(rangeReads()).toBe(reads + 1);
    // Another run reads its own sections.
    await results.open('r2');
    await settle();
    const before = rangeReads();
    failRanges(true);
    await expect(results.readSection({ runId: 'r2', qIdx: 0, round: 0, rank: 2 })).rejects.toThrow('the stored output is gone');
    failRanges(false);
    expect(await results.readSection({ runId: 'r2', qIdx: 0, round: 0, rank: 2 })).toBe(' Score = 70 bits\n\nQuery  5  ACGT  45\n');
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
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 1 });
    await settle();
    expect(results.state.get().detail).toMatchObject({ state: 'ready', heading: '> s1 subject number 1\nLength=500\n\n' });
    expect(rangeReads()).toBe(3);
    // An HSP of another subject: its query's rank finds its row.
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 2 });
    await settle();
    expect(results.state.get()).toMatchObject({ sIdx: 0, detail: { state: 'ready', heading: '> s0 subject number 0\nLength=500\n\n' } });
    // A rank that the query does not have selects nothing.
    results.selectHsp({ runId: 'r1', qIdx: 0, round: 0, rank: 9 });
    results.selectHsp({ runId: 'r1', qIdx: 1, round: 0, rank: 0 });
    expect(results.state.get().hsp).toEqual({ runId: 'r1', qIdx: 0, round: 0, rank: 2 });
  });

  it('separates multi-round PSI-BLAST HSPs by round, defaults to the last round, and switches rounds cleanly', async () => {
    const multiRoundSpecs: readonly Spec[] = [
      // Round 1: s0 is top hit (e=1.35e-19), s1 is second (e=2.76e-18)
      { q: 0, s: 0, round: 1, bits: 63.9, e: 1.35e-19, coords: [1, 30, 1, 30], text: { e: '1.35e-19', bits: '63.9' }, out0: { section: ' Score = 63.9 bits (R1 s0)\n' } },
      { q: 0, s: 1, round: 1, bits: 60.5, e: 2.76e-18, coords: [1, 30, 1, 30], text: { e: '2.76e-18', bits: '60.5' }, out0: { section: ' Score = 60.5 bits (R1 s1)\n' } },
      // Round 2: s1 becomes top hit (e=5.59e-18), s0 is second (e=8.24e-18)
      { q: 0, s: 1, round: 2, bits: 59.7, e: 5.59e-18, coords: [1, 30, 1, 30], text: { e: '5.59e-18', bits: '59.7' }, out0: { section: ' Score = 59.7 bits (R2 s1)\n' } },
      { q: 0, s: 0, round: 2, bits: 59.3, e: 8.24e-18, coords: [1, 30, 1, 30], text: { e: '8.24e-18', bits: '59.3' }, out0: { section: ' Score = 59.3 bits (R2 s0)\n' } },
    ];
    const index = buildResultIndex(hspTable(makeRun(multiRoundSpecs).records));
    const q0 = index.queries.get(0)!;
    expect(q0.rounds.map((r) => [r.round, r.subjects.map((s) => s.sIdx), r.hspCount])).toEqual([
      [1, [0, 1], 2],
      [2, [1, 0], 2],
    ]);
    expect(q0.subjects.map((s) => s.sIdx)).toEqual([1, 0]);
    expect(q0.hspCount).toBe(2);

    const { results } = setup({ specs: multiRoundSpecs, argv: ['psiblast', '-query', 'q.fa', '-subject', 's.fa', '-num_iterations', '2'] });
    await results.open('r1');
    await settle();
    let state = results.state.get();
    expect(state.rounds).toEqual([1, 2]);
    expect(state.round).toBe(2);
    expect(state.queryTotals).toEqual({ subjects: 2, hsps: 2 });
    expect(state.subjects.map((s) => [s.sIdx, s.hspCount, s.first.evalue])).toEqual([
      [1, 1, '5.59e-18'],
      [0, 1, '8.24e-18'],
    ]);
    expect(state.hsp).toEqual({ runId: 'r1', qIdx: 0, round: 2, rank: 0 });
    expect(state.detail?.section).toBe(' Score = 59.7 bits (R2 s1)\n');

    results.selectRound(1);
    await settle();
    state = results.state.get();
    expect(state.round).toBe(1);
    expect(state.queryTotals).toEqual({ subjects: 2, hsps: 2 });
    expect(state.subjects.map((s) => [s.sIdx, s.hspCount, s.first.evalue])).toEqual([
      [0, 1, '1.35e-19'],
      [1, 1, '2.76e-18'],
    ]);
    expect(state.hsp).toEqual({ runId: 'r1', qIdx: 0, round: 1, rank: 0 });
    expect(state.detail?.section).toBe(' Score = 63.9 bits (R1 s0)\n');
    expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 2, rank: 0 })).toBe(' Score = 59.7 bits (R2 s1)\n');
    expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 1, rank: 0 })).toBe(' Score = 63.9 bits (R1 s0)\n');
  });

  describe('outfmt 0 of a run made with -html', () => {
    // NCBI's HTML: a report link and an anchor in the heading, a font tag on a masked row, and the title
    // encoded with &quot; &amp; &lt; &gt; only.
    const HEADING = '><a title="Show report for s1" href="https://www.ncbi.nlm.nih.gov/protein/s1?a=1&b=2" >s1</a><a name=s1></a> subject &lt;one&gt; &amp;lt; &quot;q&quot;\nLength=500\n\n';
    const SECTION = ' Score = 90 bits\n\nQuery  1  <font color="#000000">acgt</font>  50\n';
    const specs: readonly Spec[] = [{ q: 0, s: 1, bits: 90, e: 1e-20, coords: [1, 50, 101, 150], out0: { heading: HEADING, section: SECTION } }];
    const PLAIN_HEADING = '>s1 subject <one> &lt; "q"\nLength=500\n\n';

    it('shows the text that the HTML displays in the detail, the lists and the Range read, and keeps the stored bytes', async () => {
      const { results } = setup({ specs, argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-html'] });
      await results.open('r1');
      // The first range read of such a run also loads the module that reads its HTML (domain/outfmt0-html.ts).
      await expect.poll(() => results.state.get().detail?.state, { timeout: 10_000 }).toBe('ready');
      const state = results.state.get();
      expect(state.detail).toMatchObject({ state: 'ready', heading: PLAIN_HEADING });
      expect(state.detail!.section).toBe(' Score = 90 bits\n\nQuery  1  acgt  50\n');
      // The subject list's title and the readers of the Alignments: one decoding for all of them.
      expect(headingTitle(state.detail!.heading!)).toBe('s1 subject <one> &lt; "q"');
      expect(await results.readSection({ runId: 'r1', qIdx: 0, round: 0, rank: 0 })).toBe(' Score = 90 bits\n\nQuery  1  acgt  50\n');
      results.requestHeadings([1]);
      await expect.poll(() => results.state.get().headings.get(1)).toBe(PLAIN_HEADING);
    });

    it('shows a run without -html as written, a plain title with a literal &lt; included', async () => {
      const plain: readonly Spec[] = [
        { q: 0, s: 1, bits: 90, e: 1e-20, coords: [1, 50, 101, 150], out0: { heading: '> s1 a &lt;b&gt; <i>\nLength=500\n\n', section: ' Score = 90 bits\n\nQuery  1  <b>  50\n' } },
      ];
      const { results } = setup({ specs: plain });
      await results.open('r1');
      await settle();
      const state = results.state.get();
      expect(state.detail).toMatchObject({ heading: '> s1 a &lt;b&gt; <i>\nLength=500\n\n' });
      expect(state.detail!.section).toBe(' Score = 90 bits\n\nQuery  1  <b>  50\n');
    });
  });

  it("names the task of a run: the argv's, else the engine's default", async () => {
    const byDefault = setup();
    await byDefault.results.open('r1');
    expect(byDefault.results.state.get().loaded!.task).toBe('megablast');
    const given = setup({ argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-task', 'blastn'] });
    await given.results.open('r1');
    expect(given.results.state.get().loaded!.task).toBe('blastn');
  });

  it('gives no tool badge when describeTool rejects or is absent, and shows Supported once describeTool succeeds', async () => {
    // cases.tsv k_obinary_file argv: Supported with -xmlhelp grammar, no key with fallback grammar.
    const argv = ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-unit', '5', '-sformat', 'obinary', '-out', 'c.counts'];
    const description = await new FakeEngine().describeTool('windowmasker');
    const key = toolOptionKey(argv, toolOptionGrammar(description))!;
    const toolRun: ToolRunView = {
      snapshot: {
        kind: 'tool',
        runId: 't1',
        number: 1,
        tool: 'windowmasker',
        stage: 'mk_counts',
        argv,
        directory: [{ name: 'genome.fa', size: 100, origin: { kind: 'pasted' } }],
        queuedAt: 0,
      },
      status: 'completed',
      record: { runtimePath: 'serial', exitStatus: 0 },
    };
    const runs = new Store<AppState>({ runs: [toolRun] });
    const verification: VerificationTable = {
      ncbi: '2.17.0',
      sources: [],
      programs: {},
      tools: { windowmasker: { optionSets: [key] } },
    };
    const data: ResultsDeps['data'] = {
      readHitTable: async () => hspTable([]),
      readOutput: async () => new Uint8Array(),
      readOutputRange: async () => new Uint8Array(),
      readDiagnostics: async () => '',
    };

    let describeCalls = 0;
    const results = new ResultsBrowser({
      data,
      describe: async () => DESCRIPTION,
      describeTool: async () => {
        describeCalls++;
        if (describeCalls === 1) throw new Error('engine busy');
        return description;
      },
      runs,
      verification,
    });

    // Rejected describeTool: fallback grammar formats CommandText, but no badge is shown and the rejection is not cached.
    await results.open('t1');
    await settle();
    expect(results.state.get().phase).toBe('ready');
    expect(results.state.get().grammar).toBeDefined();
    expect(results.state.get().toolBadge).toBeUndefined();

    // Opening again asks describeTool again; once it resolves with -xmlhelp, the badge is Supported.
    await results.open('t1');
    await settle();
    expect(describeCalls).toBe(2);
    expect(results.state.get().toolBadge).toMatchObject({ level: 'certified', label: 'Supported' });

    // Without describeTool at all, no badge is shown either.
    const noDescribe = new ResultsBrowser({
      data,
      describe: async () => DESCRIPTION,
      runs,
      verification,
    });
    await noDescribe.open('t1');
    await settle();
    expect(noDescribe.state.get().grammar).toBeDefined();
    expect(noDescribe.state.get().toolBadge).toBeUndefined();
  });
});
