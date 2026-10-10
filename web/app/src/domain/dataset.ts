// Datasets (plan §5.2, §5.4): the record table of one source FASTA, and immutable
// revisions of it that say which records a run uses. The record fields are those of the
// ABI v2 *scan* response (docs/web/abi_v2.md §9), so no mapping layer is needed.

/**
 * The FASTA reader whose rules the index scan follows (ABI v2 `scan_begin`, abi_v2.md §9;
 * plan TD-8): the engine's port of NCBI BLAST+'s reader, with the flags of nucleotide input
 * (kind 1) or of protein input (kind 2), so that the record table has the records that the
 * engine searches (plan DW-23 (6)). The program and the role choose it (`indexParser`).
 */
export type FastaParserKind = 1 | 2;

/** How to find the byte offset of residue `i` of a record (ABI v2 §9 `line_layout`). */
export type LineLayout =
  | { readonly kind: 'uniform'; readonly width: number; readonly eol: number }
  | { readonly kind: 'checkpoints'; readonly every: number; readonly offsets: readonly number[] };

/** One record of the index scan. Offsets are byte offsets in the source and stay below 2^53. */
export interface IndexedRecord {
  readonly index: number;
  readonly id: string;
  /**
   * Offset of the defline's `>`. The input's first record without a defline (residues
   * before the first `>`) has 0 here and in `sequence_offset`, so its bytes start the input.
   */
  readonly header_offset: number;
  /** Offset of the first byte after the defline's end of line (0 for a first record without a defline). */
  readonly sequence_offset: number;
  /** Offset of the next record's `>`, or the end of the input. */
  readonly end_offset: number;
  /** The number of residues that the reader stores. */
  readonly length: number;
  readonly line_layout: LineLayout;
  /**
   * Count of each stored residue, upper-cased (kind 1 counts `U` as `T`), keyed by the
   * character; a byte outside 0x21-0x7E would be "0xNN", but kinds 1 and 2 store letters and `*` only.
   */
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

/**
 * LOSAT Web's refusal of an input whose first line is not a defline and may be a sequence
 * identifier that NCBI BLAST+ fetches through a data loader (the adapter's scan: `the first line
 * ("…") is not a defline and may be ... (start the input with a '>' defline)`). A defline in
 * front of the text is the fix for this refusal only, not for a gap line or a `Near line N` one.
 */
export function isFirstLineRefusal(message: string): boolean {
  return /\bthe first line \(.*\) is not a defline\b/s.test(message);
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

/**
 * The IDs that more than one record of the list has. Records are told apart by their
 * index in the record table, never by ID (design §2.1: duplicate IDs keep an internal ID).
 */
export function duplicateIds(records: ReadonlyArray<{ readonly id: string }>): ReadonlySet<string> {
  const seen = new Set<string>();
  const duplicates = new Set<string>();
  for (const { id } of records) {
    if (seen.has(id)) duplicates.add(id);
    else seen.add(id);
  }
  return duplicates;
}
