// Datasets (plan §5.2, §5.4): the record table of one source FASTA, and immutable
// revisions of it that say which records a run uses. The record fields are those of the
// ABI v2 *scan* response (docs/web/abi_v2.md §9), so no mapping layer is needed.

/** The FASTA reader whose rules the index scan follows (ABI v2 `scan_begin`, plan TD-8). */
export type FastaParserKind = 0 | 1;

/** How to find the byte offset of residue `i` of a record (ABI v2 §9 `line_layout`). */
export type LineLayout =
  | { readonly kind: 'uniform'; readonly width: number; readonly eol: number }
  | { readonly kind: 'checkpoints'; readonly every: number; readonly offsets: readonly number[] };

/** One record of the index scan. Offsets are byte offsets in the source and stay below 2^53. */
export interface IndexedRecord {
  readonly index: number;
  readonly id: string;
  /** Offset of the `>` of the header line. */
  readonly header_offset: number;
  /** Offset of the first byte after the header line. */
  readonly sequence_offset: number;
  /** Offset of the next header line, or the end of the input. */
  readonly end_offset: number;
  /** Sequence length in bytes, as the parser reports it. */
  readonly length: number;
  readonly line_layout: LineLayout;
  /** Count of each byte value of the sequence: 0x21-0x7E as the character, others as "0xNN". */
  readonly residue_counts: Readonly<Record<string, number>>;
}

export interface DatasetRecord extends IndexedRecord {
  /** Lower-case hex SHA-256 of the record's original bytes [header_offset, end_offset). */
  readonly sha256: string;
}

/**
 * An immutable view of one source: its record table and the records that runs leave out.
 * Changing the selection creates a new revision (plan §5.2).
 */
export interface DatasetRevision {
  readonly revisionId: string;
  readonly sourceId: string;
  readonly parser: FastaParserKind;
  readonly records: readonly DatasetRecord[];
  /** Indices of the records that runs leave out, ascending and without duplicates. */
  readonly excluded: readonly number[];
}

/** What the engine's `register` reports for each record (ABI v2 §9 *register*). */
export interface RecordKey {
  readonly id: string;
  readonly length: number;
}

export function includedRecords(revision: DatasetRevision): readonly DatasetRecord[] {
  const excluded = new Set(revision.excluded);
  return revision.records.filter((record) => !excluded.has(record.index));
}

export function recordKey(record: RecordKey): RecordKey {
  return { id: record.id, length: record.length };
}

/** Sorts and deduplicates record indices; rejects an index that is not in the table. */
export function normalizeExclusion(indices: readonly number[], recordCount: number): readonly number[] {
  for (const index of indices) {
    if (!Number.isInteger(index) || index < 0 || index >= recordCount) {
      throw new RangeError(`record index ${index} is not in the record table (${recordCount} records)`);
    }
  }
  return Object.freeze([...new Set(indices)].sort((a, b) => a - b));
}

/**
 * Compares the record table of a run input with the records that the engine read from the
 * same bytes (plan §5.4). Returns a description of the first difference, or undefined if
 * they agree.
 */
export function recordMismatch(expected: readonly RecordKey[], actual: readonly RecordKey[]): string | undefined {
  const count = Math.min(expected.length, actual.length);
  for (let i = 0; i < count; i++) {
    const want = expected[i]!;
    const got = actual[i]!;
    if (want.id !== got.id || want.length !== got.length) {
      return (
        `record ${i + 1} is "${want.id}" (length ${want.length}) in the record table, ` +
        `but the engine read "${got.id}" (length ${got.length})`
      );
    }
  }
  if (expected.length !== actual.length) {
    return `the record table has ${expected.length} records, but the engine read ${actual.length}`;
  }
  return undefined;
}
