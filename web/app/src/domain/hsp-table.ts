// The HSP records of a run as columns (plan §5.7, design §10.2): the small fields of every
// HSP record (docs/web/abi_v2.md §8) in typed arrays, without the aligned sequences. The
// results screen sorts and filters with these engine values and reads the text it shows
// from the byte ranges (`out6`, `out0`, `out0_subject`) of the stored outputs. The Data
// worker builds the table and transfers its buffers, so a run with many HSPs never becomes
// an array of objects on the UI thread.

/** The fields of an HSP record that the table keeps (the names of `HspRecord`). */
export interface HspSummary {
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
  readonly out6: readonly [number, number] | null;
  readonly out0: readonly [number, number] | null;
  readonly out0_subject: readonly [number, number] | null;
}

/**
 * Row `i` of every column is the HSP with `index[i]`, in the order of the records (the
 * final hit list of the run). Frames are 0 where the record has none; ranges are -1 where
 * the record has none. Byte offsets may exceed 2^32 and stay below 2^53 (Float64Array).
 */
export interface HspTable {
  readonly count: number;
  readonly index: Int32Array;
  readonly qIdx: Int32Array;
  readonly sIdx: Int32Array;
  readonly rank: Int32Array;
  readonly rawScore: Float64Array;
  readonly bitScore: Float64Array;
  readonly eValue: Float64Array;
  readonly qStart: Float64Array;
  readonly qEnd: Float64Array;
  readonly sStart: Float64Array;
  readonly sEnd: Float64Array;
  readonly queryFrame: Int8Array;
  readonly subjectFrame: Int8Array;
  readonly out6Start: Float64Array;
  readonly out6End: Float64Array;
  readonly out0Start: Float64Array;
  readonly out0End: Float64Array;
  readonly out0SubjectStart: Float64Array;
  readonly out0SubjectEnd: Float64Array;
}

export function hspTable(records: readonly HspSummary[]): HspTable {
  const count = records.length;
  const table = {
    count,
    index: new Int32Array(count),
    qIdx: new Int32Array(count),
    sIdx: new Int32Array(count),
    rank: new Int32Array(count),
    rawScore: new Float64Array(count),
    bitScore: new Float64Array(count),
    eValue: new Float64Array(count),
    qStart: new Float64Array(count),
    qEnd: new Float64Array(count),
    sStart: new Float64Array(count),
    sEnd: new Float64Array(count),
    queryFrame: new Int8Array(count),
    subjectFrame: new Int8Array(count),
    out6Start: new Float64Array(count),
    out6End: new Float64Array(count),
    out0Start: new Float64Array(count),
    out0End: new Float64Array(count),
    out0SubjectStart: new Float64Array(count),
    out0SubjectEnd: new Float64Array(count),
  };
  records.forEach((record, i) => {
    table.index[i] = record.index;
    table.qIdx[i] = record.q_idx;
    table.sIdx[i] = record.s_idx;
    table.rank[i] = record.rank;
    table.rawScore[i] = record.raw_score;
    table.bitScore[i] = record.bit_score;
    table.eValue[i] = record.e_value;
    table.qStart[i] = record.q_start;
    table.qEnd[i] = record.q_end;
    table.sStart[i] = record.s_start;
    table.sEnd[i] = record.s_end;
    table.queryFrame[i] = record.query_frame ?? 0;
    table.subjectFrame[i] = record.subject_frame ?? 0;
    table.out6Start[i] = record.out6?.[0] ?? -1;
    table.out6End[i] = record.out6?.[1] ?? -1;
    table.out0Start[i] = record.out0?.[0] ?? -1;
    table.out0End[i] = record.out0?.[1] ?? -1;
    table.out0SubjectStart[i] = record.out0_subject?.[0] ?? -1;
    table.out0SubjectEnd[i] = record.out0_subject?.[1] ?? -1;
  });
  return table;
}

/** A byte range [start, end) of a stored output, or undefined where the record has none. */
export type ByteRange = readonly [number, number];

const range = (start: Float64Array, end: Float64Array, row: number): ByteRange | undefined =>
  start[row]! < 0 ? undefined : [start[row]!, end[row]!];

export const out6Range = (table: HspTable, row: number) => range(table.out6Start, table.out6End, row);
export const out0Range = (table: HspTable, row: number) => range(table.out0Start, table.out0End, row);
export const out0SubjectRange = (table: HspTable, row: number) =>
  range(table.out0SubjectStart, table.out0SubjectEnd, row);

/** A frame of the record, or undefined where it has none. */
export const frame = (column: Int8Array, row: number): number | undefined => (column[row] === 0 ? undefined : column[row]);

/**
 * The smaller and the larger subject coordinate of an HSP record: NCBI's "Range n: a to b" of an
 * alignment names the subject's positions in ascending order (docs/web/ncbi_ui_mapping.md
 * "Alignments"). The coordinates are the record's; nothing is computed from them.
 */
export function subjectSpan(table: HspTable, row: number): readonly [number, number] {
  const [start, end] = [table.sStart[row]!, table.sEnd[row]!];
  return start <= end ? [start, end] : [end, start];
}
