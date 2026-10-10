// LOSAT Web's own files of a run's HSPs (S15 item 2; design §12.1, §12.3): CSV and JSON made from
// the HSP records, for the whole run, after the view filters, or for the subjects marked in the
// Descriptions. They are application formats, not NCBI BLAST outputs, and say so where the format
// has room for it (CSV has none: the export screen says it). Nothing here computes or formats a
// BLAST value (web/AGENTS.md rule 1): the outfmt 6 fields are the text of the HSP's outfmt 6 row
// as the engine wrote it, and the record's numbers are the engine's. The functions make the parts
// of a file; application/result-export.ts writes them in order, in blocks (the Writer contract).
import { hspLabel } from './extraction';
import type { HspTable } from './hsp-table';
import { OUTFMT6_FIELDS, type Outfmt6Row } from './outfmt6';
import { filterQuery, type ResultIndex, type SubjectGroup, type ViewFilters } from './result-index';

// --- scopes ---------------------------------------------------------------------------------------

/** Which HSPs a file holds (design §12.1: the run's whole result, after the filters, the selection). */
export type ExportScope = 'all' | 'filtered' | 'marked';
export const EXPORT_SCOPES: readonly ExportScope[] = Object.freeze(['all', 'filtered', 'marked']);

export const SCOPE_LABELS: Readonly<Record<ExportScope, string>> = Object.freeze({
  all: 'Whole run',
  filtered: 'After the view filters',
  marked: 'Marked subjects',
});

/** Table rows in the engine's order (the HSP `index`); the table's rows usually are already. */
export function engineOrder(table: HspTable, rows: ArrayLike<number>): Int32Array {
  const ordered = Int32Array.from(rows);
  for (let i = 1; i < ordered.length; i++) {
    if (table.index[ordered[i - 1]!]! > table.index[ordered[i]!]!) {
      return ordered.sort((a, b) => table.index[a]! - table.index[b]!);
    }
  }
  return ordered;
}

/** Every HSP of the run, in the engine's order. */
export function allRows(table: HspTable): Int32Array {
  return engineOrder(table, Int32Array.from({ length: table.count }, (_, row) => row));
}

/**
 * The HSPs that the view's HSP filters keep (E value, bit score, subject text: `filterQuery`, as
 * the Descriptions apply them) in the queries that the query filters list, in the engine's order.
 * `subjectIds` gives the texts that the subject text filter searches, as for the Descriptions.
 */
export function filteredRows(
  index: ResultIndex,
  queries: Iterable<number>,
  filters: ViewFilters,
  subjectIds: (subject: SubjectGroup) => readonly string[],
): Int32Array {
  const rows: number[] = [];
  for (const qIdx of queries) {
    const query = index.queries.get(qIdx);
    if (query === undefined) continue;
    for (const subject of filterQuery(query, filters, index.table, subjectIds).subjects) {
      for (const row of subject.rows) rows.push(row);
    }
  }
  return engineOrder(index.table, rows);
}

/**
 * The HSPs of the marked subjects that the view filters keep: `subjects` are the subjects that the
 * Descriptions list (their rows already filtered), in the engine's order whatever the list's sort.
 */
export function markedRows(table: HspTable, subjects: readonly SubjectGroup[], marked: ReadonlySet<number>): Int32Array {
  const rows: number[] = [];
  for (const subject of subjects) {
    if (!marked.has(subject.sIdx)) continue;
    for (const row of subject.rows) rows.push(row);
  }
  return engineOrder(table, rows);
}

// --- one HSP ----------------------------------------------------------------------------------------

/** What the files say about one HSP besides its record's numbers. */
export interface ExportedHsp {
  /** The run's number in the working session. */
  readonly run: number;
  /** 0-based positions of the query and subject records in the run's inputs (the record's q_idx, s_idx). */
  readonly qIdx: number;
  readonly sIdx: number;
  /** The records' IDs in the run's record tables (empty for a record without a title). */
  readonly queryId: string;
  readonly subjectId: string;
  readonly index: number;
  readonly rank: number;
  /** The fields of the HSP's outfmt 6 row as written; undefined where outfmt 6 does not show the HSP. */
  readonly outfmt6: Outfmt6Row | undefined;
  /** The record's frames (null where it has none). */
  readonly queryFrame: number | null;
  readonly subjectFrame: number | null;
  /** Whether outfmt 0 shows the HSP's alignment (the record has an `out0` range). */
  readonly inOutfmt0: boolean;
}

// --- CSV --------------------------------------------------------------------------------------------

/** The columns of the CSV: the HSP's identity, its outfmt 6 fields as written, its frames and outfmt 0. */
export const CSV_COLUMNS: readonly string[] = Object.freeze([
  'run',
  'query_record',
  'query_id',
  'subject_record',
  'subject_id',
  'hsp',
  'index',
  'rank',
  ...OUTFMT6_FIELDS,
  'query_frame',
  'subject_frame',
  'in_outfmt0',
]);

/** RFC 4180 ends every record, the header's too, with CR LF. */
export const CSV_LINE_END = '\r\n';

/** A field as RFC 4180 writes it: in quotes, its quotes doubled, when it holds `,` `"` CR or LF. */
export function csvField(text: string): string {
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export function csvRecord(fields: readonly string[]): string {
  return fields.map(csvField).join(',') + CSV_LINE_END;
}

export const csvHeader = (): string => csvRecord(CSV_COLUMNS);

/** The CSV record of an HSP. IDs are written as they are, never altered (a spreadsheet may read `=…` as a formula). */
export function csvLine(hsp: ExportedHsp): string {
  return csvRecord([
    String(hsp.run),
    String(hsp.qIdx + 1),
    hsp.queryId,
    String(hsp.sIdx + 1),
    hsp.subjectId,
    hspLabel(hsp.qIdx, hsp.rank),
    String(hsp.index),
    String(hsp.rank),
    ...OUTFMT6_FIELDS.map((name) => hsp.outfmt6?.[name] ?? ''),
    hsp.queryFrame === null ? '' : String(hsp.queryFrame),
    hsp.subjectFrame === null ? '' : String(hsp.subjectFrame),
    String(hsp.inOutfmt0),
  ]);
}

// --- what the JSON and the report say about the run and the scope -------------------------------------

export interface ExportInput {
  /** The name passed as -query / -subject. */
  readonly name: string;
  readonly records: number;
  /** Length and SHA-256 of the FASTA bytes given to the engine (the included records, joined). */
  readonly bytes: number;
  readonly sha256: string;
}

export interface ExportRun {
  readonly number: number;
  readonly title?: string;
  readonly program: string;
  /** The program's name as the screen shows it, e.g. BLASTN. */
  readonly programLabel: string;
  /** The task that the search used (the argv's, or the engine's default). */
  readonly task?: string;
  readonly argv: readonly string[];
  readonly requestedThreads: number | 'auto';
  readonly engineBuild?: string;
  readonly runtimePath?: string;
  readonly threads?: number;
  /** Times in milliseconds since the epoch. */
  readonly queuedAt: number;
  readonly startedAt?: number;
  readonly endedAt?: number;
  readonly group?: { readonly position: number; readonly size: number };
  readonly query: ExportInput;
  readonly subject: ExportInput;
  /** The run's verification badge (domain/verification.ts): about the engine and the options, not about the file. */
  readonly verification: {
    readonly level: string;
    readonly label: string;
    readonly details: readonly string[];
    readonly exceptions: readonly string[];
  };
}

export interface ExportRecordRef {
  /** 0-based position of the record in the run's input. */
  readonly position: number;
  readonly id: string;
}

export interface ExportScopeInfo {
  readonly scope: ExportScope;
  readonly hsps: number;
  /** The view filters that made the scope (the filtered and marked scopes). */
  readonly filters?: ViewFilters;
  /** The query whose subjects are marked (the marked scope). */
  readonly query?: ExportRecordRef;
  /** The marked subjects, in the order that the Descriptions list them (the marked scope). */
  readonly markedSubjects?: readonly ExportRecordRef[];
}

/** The view filters, as the JSON names them; only those set. */
export function filtersJson(filters: ViewFilters): Record<string, number | string | boolean> {
  return {
    ...(filters.maxEValue === undefined ? {} : { max_e_value: filters.maxEValue }),
    ...(filters.minBitScore === undefined ? {} : { min_bit_score: filters.minBitScore }),
    ...((filters.subjectText ?? '') === '' ? {} : { subject_id_contains: filters.subjectText! }),
    ...(filters.queriesWithHitsOnly === undefined ? {} : { queries_with_hits_only: filters.queriesWithHitsOnly }),
    ...((filters.queryText ?? '') === '' ? {} : { query_id_contains: filters.queryText! }),
  };
}

/** The view filters in words (the screen's labels), one per filter set; none for no filter. */
export function filterWords(filters: ViewFilters): readonly string[] {
  return [
    ...(filters.maxEValue === undefined ? [] : [`E value at most ${filters.maxEValue}`]),
    ...(filters.minBitScore === undefined ? [] : [`Bit score at least ${filters.minBitScore}`]),
    ...((filters.subjectText ?? '') === '' ? [] : [`Subject ID contains "${filters.subjectText}"`]),
    ...(filters.queriesWithHitsOnly === true ? ['Queries with hits only'] : []),
    ...((filters.queryText ?? '') === '' ? [] : [`Query ID contains "${filters.queryText}"`]),
  ];
}

const iso = (ms: number): string => new Date(ms).toISOString();

// --- JSON -------------------------------------------------------------------------------------------

export const JSON_FORMAT = 'LOSAT Web HSP export';
export const JSON_SCHEMA = 1;
export const JSON_NOTE =
  'Made by LOSAT Web from the HSP records of one run: an application format of LOSAT Web, not an NCBI BLAST output, ' +
  'and not compared with NCBI BLAST+. "outfmt6" holds the fields of the HSP\'s outfmt 6 row as the engine wrote them ' +
  '(the formatted values; null where outfmt 6 does not show the HSP). "record" holds the numbers of the HSP record as ' +
  'the engine wrote them (not rounded), and its aligned rows when they are included. query_record and subject_record ' +
  'are 1-based positions in the run\'s inputs; index, rank, q_idx and s_idx are 0-based.';

/** The numbers and aligned rows of an HSP record (docs/web/abi_v2.md §8) that the JSON writes. */
export interface ExportedRecord {
  readonly index: number;
  readonly q_idx: number;
  readonly s_idx: number;
  readonly rank: number;
  readonly raw_score: number;
  readonly bit_score: number;
  readonly e_value: number;
  readonly q_start: number;
  readonly q_end: number;
  readonly s_start: number;
  readonly s_end: number;
  readonly query_frame: number | null;
  readonly subject_frame: number | null;
  readonly subject_length: number | null;
  readonly query_aligned: string | null;
  readonly subject_aligned: string | null;
}

const RECORD_NUMBERS = [
  'index',
  'q_idx',
  's_idx',
  'rank',
  'raw_score',
  'bit_score',
  'e_value',
  'q_start',
  'q_end',
  's_start',
  's_end',
  'query_frame',
  'subject_frame',
  'subject_length',
] as const;

const inputJson = (input: ExportInput) => ({
  name: input.name,
  records: input.records,
  engine_input_bytes: input.bytes,
  engine_input_sha256: input.sha256,
});

export function runJson(run: ExportRun): Record<string, unknown> {
  return {
    number: run.number,
    ...(run.title === undefined ? {} : { title: run.title }),
    program: run.program,
    ...(run.task === undefined ? {} : { task: run.task }),
    argv: [...run.argv],
    requested_threads: run.requestedThreads,
    ...(run.engineBuild === undefined ? {} : { engine_build: run.engineBuild }),
    ...(run.runtimePath === undefined ? {} : { runtime_path: run.runtimePath }),
    ...(run.threads === undefined ? {} : { threads: run.threads }),
    queued_at: iso(run.queuedAt),
    ...(run.startedAt === undefined ? {} : { started_at: iso(run.startedAt) }),
    ...(run.endedAt === undefined ? {} : { ended_at: iso(run.endedAt) }),
    ...(run.group === undefined ? {} : { group: { position: run.group.position, size: run.group.size } }),
    query: inputJson(run.query),
    subject: inputJson(run.subject),
    verification: {
      level: run.verification.level,
      label: run.verification.label,
      details: [...run.verification.details],
      exceptions: [...run.verification.exceptions],
    },
  };
}

const recordRef = (ref: ExportRecordRef) => ({ record: ref.position + 1, id: ref.id });

export function scopeJson(scope: ExportScopeInfo, aligned: boolean): Record<string, unknown> {
  return {
    name: scope.scope,
    label: SCOPE_LABELS[scope.scope],
    hsps: scope.hsps,
    filters: scope.filters === undefined ? null : filtersJson(scope.filters),
    ...(scope.query === undefined ? {} : { query: recordRef(scope.query) }),
    ...(scope.markedSubjects === undefined ? {} : { marked_subjects: scope.markedSubjects.map(recordRef) }),
    aligned_sequences: aligned,
  };
}

/** The JSON document up to the opening of its `hsps` array. */
export function jsonHead(run: ExportRun, scope: ExportScopeInfo, options: { readonly aligned: boolean; readonly exportedAt: number }): string {
  const head = JSON.stringify({
    format: JSON_FORMAT,
    schema: JSON_SCHEMA,
    note: JSON_NOTE,
    exported_at: iso(options.exportedAt),
    run: runJson(run),
    scope: scopeJson(scope, options.aligned),
  });
  return `${head.slice(0, -1)},"hsps":[`;
}

/** One HSP of the `hsps` array, on its own line; `first` for the first, which no comma precedes. */
export function jsonHsp(hsp: ExportedHsp, record: ExportedRecord, aligned: boolean, first: boolean): string {
  const numbers: Record<string, number | null | string> = {};
  for (const name of RECORD_NUMBERS) numbers[name] = record[name];
  if (aligned) {
    numbers.query_aligned = record.query_aligned;
    numbers.subject_aligned = record.subject_aligned;
  }
  const item = {
    run: hsp.run,
    query_record: hsp.qIdx + 1,
    query_id: hsp.queryId,
    subject_record: hsp.sIdx + 1,
    subject_id: hsp.subjectId,
    hsp: hspLabel(hsp.qIdx, hsp.rank),
    index: hsp.index,
    rank: hsp.rank,
    in_outfmt0: hsp.inOutfmt0,
    outfmt6: hsp.outfmt6 === undefined ? null : hsp.outfmt6,
    record: numbers,
  };
  return `${first ? '\n' : ',\n'}${JSON.stringify(item)}`;
}

export const JSON_TAIL = '\n]}\n';
