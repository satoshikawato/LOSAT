// The results screen's view of a run (plan §5.7, docs/web/results_columns.md): the HSPs of
// each query grouped by subject in the engine's order, sorting and filtering with the
// engine's values, the limits that a query may have reached, and the orientation of an
// HSP. Nothing here computes or formats a BLAST value.
import type { HspTable } from './hsp-table';
import type { SequenceKind } from './programs';

export interface SubjectGroup {
  /** 0-based record index of the subject in the run's subject input. */
  readonly sIdx: number;
  /** Table rows of the subject's HSPs in this query, in the engine's order (rank). */
  readonly rows: readonly number[];
}

export interface QueryGroup {
  /** 0-based record index of the query in the run's query input. */
  readonly qIdx: number;
  /** Subjects in the engine's order: the order of their first HSPs. */
  readonly subjects: readonly SubjectGroup[];
  readonly hspCount: number;
}

export interface ResultIndex {
  readonly table: HspTable;
  /** The queries that have HSPs; a query without HSPs has no group. */
  readonly queries: ReadonlyMap<number, QueryGroup>;
}

/** Groups the rows of the table by query, then by subject, in the order of the HSP index. */
export function buildResultIndex(table: HspTable): ResultIndex {
  const order = Array.from({ length: table.count }, (_, row) => row).sort((a, b) => table.index[a]! - table.index[b]!);
  const queries = new Map<number, { qIdx: number; subjects: Map<number, { sIdx: number; rows: number[] }>; hspCount: number }>();
  for (const row of order) {
    const qIdx = table.qIdx[row]!;
    let query = queries.get(qIdx);
    if (query === undefined) {
      query = { qIdx, subjects: new Map(), hspCount: 0 };
      queries.set(qIdx, query);
    }
    const sIdx = table.sIdx[row]!;
    let subject = query.subjects.get(sIdx);
    if (subject === undefined) {
      subject = { sIdx, rows: [] };
      query.subjects.set(sIdx, subject);
    }
    subject.rows.push(row);
    query.hspCount++;
  }
  const frozen = new Map<number, QueryGroup>();
  for (const [qIdx, query] of queries) {
    frozen.set(
      qIdx,
      Object.freeze({
        qIdx,
        subjects: Object.freeze([...query.subjects.values()].map((s) => Object.freeze({ sIdx: s.sIdx, rows: Object.freeze(s.rows) }))),
        hspCount: query.hspCount,
      }),
    );
  }
  return Object.freeze({ table, queries: frozen });
}

// --- sorting with engine values ---------------------------------------------------------------

export type SubjectSortKey = 'order' | 'bitScore' | 'eValue' | 'length' | 'hsps';
export type HspSortKey = 'rank' | 'bitScore' | 'eValue' | 'qStart' | 'sStart';

export interface SortSpec<K extends string> {
  readonly key: K;
  readonly descending: boolean;
}

/**
 * Sorts subjects by a value of their first HSP (the HSP whose outfmt 6 values the subject
 * list shows), their record length, or their HSP count. Ties keep the engine's order.
 */
export function sortSubjects(
  subjects: readonly SubjectGroup[],
  sort: SortSpec<SubjectSortKey>,
  table: HspTable,
  subjectLength: (sIdx: number) => number,
): readonly SubjectGroup[] {
  if (sort.key === 'order' && !sort.descending) return subjects;
  const value = (subject: SubjectGroup): number => {
    const first = subject.rows[0]!;
    switch (sort.key) {
      case 'order':
        return table.index[first]!;
      case 'bitScore':
        return table.bitScore[first]!;
      case 'eValue':
        return table.eValue[first]!;
      case 'length':
        return subjectLength(subject.sIdx);
      case 'hsps':
        return subject.rows.length;
    }
  };
  return stableSort(subjects, value, sort.descending, (s) => table.index[s.rows[0]!]!);
}

/** Sorts the rows of HSPs by an engine value. Ties keep the engine's order. */
export function sortHsps(rows: readonly number[], sort: SortSpec<HspSortKey>, table: HspTable): readonly number[] {
  if (sort.key === 'rank' && !sort.descending) return rows;
  const column = {
    rank: table.rank,
    bitScore: table.bitScore,
    eValue: table.eValue,
    qStart: table.qStart,
    sStart: table.sStart,
  }[sort.key];
  return stableSort(rows, (row) => column[row]!, sort.descending, (row) => table.index[row]!);
}

function stableSort<T>(items: readonly T[], value: (item: T) => number, descending: boolean, order: (item: T) => number): T[] {
  const sign = descending ? -1 : 1;
  return [...items].sort((a, b) => {
    const va = value(a);
    const vb = value(b);
    if (va !== vb) return va < vb ? -sign : sign;
    return order(a) - order(b);
  });
}

// --- view filters (ViewState) -------------------------------------------------------------------

/** Filters of the view (design §11.2). They change what is shown, never the search. */
export interface ViewFilters {
  /** Shows HSPs whose E value is at most this. */
  readonly maxEValue?: number;
  /** Shows HSPs whose bit score is at least this. */
  readonly minBitScore?: number;
  /** Shows subjects whose ID contains this text (case-insensitive). */
  readonly subjectText?: string;
  /** Lists only the queries that have HSPs. */
  readonly queriesWithHitsOnly?: boolean;
  /** Lists only the queries whose ID contains this text (case-insensitive). */
  readonly queryText?: string;
}

export const NO_FILTERS: ViewFilters = Object.freeze({});

export function hasHspFilters(filters: ViewFilters): boolean {
  return filters.maxEValue !== undefined || filters.minBitScore !== undefined || (filters.subjectText ?? '') !== '';
}

export interface FilteredQuery {
  /** Subjects with at least one HSP that passes, with only the HSPs that pass. */
  readonly subjects: readonly SubjectGroup[];
  readonly hiddenSubjects: number;
  readonly hiddenHsps: number;
}

/**
 * Applies the HSP filters to one query. `subjectIds(sIdx)` gives the texts that the text
 * filter searches: the subject's ID in the record table and its outfmt 6 `sseqid`.
 */
export function filterQuery(
  query: QueryGroup,
  filters: ViewFilters,
  table: HspTable,
  subjectIds: (subject: SubjectGroup) => readonly string[],
): FilteredQuery {
  if (!hasHspFilters(filters)) return { subjects: query.subjects, hiddenSubjects: 0, hiddenHsps: 0 };
  const text = (filters.subjectText ?? '').toLowerCase();
  const subjects: SubjectGroup[] = [];
  let hiddenHsps = 0;
  for (const subject of query.subjects) {
    if (text !== '' && !subjectIds(subject).some((id) => id.toLowerCase().includes(text))) {
      hiddenHsps += subject.rows.length;
      continue;
    }
    const rows = subject.rows.filter(
      (row) =>
        (filters.maxEValue === undefined || table.eValue[row]! <= filters.maxEValue) &&
        (filters.minBitScore === undefined || table.bitScore[row]! >= filters.minBitScore),
    );
    hiddenHsps += subject.rows.length - rows.length;
    if (rows.length > 0) subjects.push(rows.length === subject.rows.length ? subject : { sIdx: subject.sIdx, rows });
  }
  return { subjects, hiddenSubjects: query.subjects.length - subjects.length, hiddenHsps };
}

// --- limits that a query may have reached -----------------------------------------------------------

export interface HitLimits {
  /** Subjects kept per query (-max_target_seqs, or the engine's described default). */
  readonly maxTargetSeqs?: number;
  /** HSPs kept per subject (-max_hsps), if set. */
  readonly maxHsps?: number;
}

/** The parts of the engine's description of an option that the limits read. */
export interface DescribedOption {
  readonly flag: string;
  readonly help: string;
  readonly defaultValue?: string;
}

/**
 * The hit limits of a run from its argv; where the argv does not set -max_target_seqs,
 * from the default that the engine describes (its `default`, or the "default: N" of its
 * help text). A limit that cannot be read is left out, so nothing is claimed about it.
 */
export function hitLimits(argv: readonly string[], described: readonly DescribedOption[]): HitLimits {
  const given = (flag: string): string | undefined => {
    const at = argv.lastIndexOf(flag);
    return at >= 0 ? argv[at + 1] : undefined;
  };
  const count = (text: string | undefined): number | undefined =>
    text !== undefined && /^\s*\+?\d+\s*$/.test(text) ? Number.parseInt(text, 10) : undefined;
  const option = described.find((o) => o.flag === '-max_target_seqs');
  const describedDefault = option?.defaultValue ?? option?.help.match(/default:\s*(\d+)/i)?.[1];
  const maxTargetSeqs = count(given('-max_target_seqs')) ?? count(describedDefault);
  const maxHsps = count(given('-max_hsps'));
  return {
    ...(maxTargetSeqs === undefined ? {} : { maxTargetSeqs }),
    ...(maxHsps === undefined ? {} : { maxHsps }),
  };
}

// --- orientation ----------------------------------------------------------------------------------

export type Orientation = 'forward' | 'reverse' | 'unknown';

/**
 * Whether the HSP runs along both sequences in the same direction. Coordinates tell it
 * (start > end is the minus strand, ABI v2 §8); for one letter (start = end) the frame's
 * sign tells it, and a protein runs forward. A BLASTN HSP of one letter has neither, so its
 * orientation is unknown from the record (only the outfmt 0 section's Strand= line shows it).
 */
export function orientation(table: HspTable, row: number, kinds: { readonly query: SequenceKind; readonly subject: SequenceKind }): Orientation {
  const direction = (start: number, end: number, frame: number, kind: SequenceKind): number => {
    if (start !== end) return start < end ? 1 : -1;
    if (frame !== 0) return Math.sign(frame);
    return kind === 'protein' ? 1 : 0;
  };
  const q = direction(table.qStart[row]!, table.qEnd[row]!, table.queryFrame[row]!, kinds.query);
  const s = direction(table.sStart[row]!, table.sEnd[row]!, table.subjectFrame[row]!, kinds.subject);
  if (q === 0 || s === 0) return 'unknown';
  return q === s ? 'forward' : 'reverse';
}
