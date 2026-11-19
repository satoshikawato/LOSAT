import type { FastaParserKind, IndexedRecord, RecordKey, ShownIds } from 'losat/domain';
// Datasets (plan §5.2, §5.4): the record table of one source FASTA, and immutable
// revisions of it that say which records a run uses. The record fields are those of the
// ABI v2 *scan* response (docs/web/abi_v2.md §9), so no mapping layer is needed.

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
  /**
   * The reader's error that ends the records (the *scan* response's `reader_error`): the bytes from
   * `offset` on, the record that the reader fails on and what follows it, belong to no record, and
   * a run of the source reads its records before them and ends with the error, as the command
   * line's query batches do (line WQA's D10, 2026-11-19). `defline`: whether that record starts
   * with a defline; one without cannot follow another input (review M2, 2026-11-20).
   */
  readonly readerError?: { readonly offset: number; readonly defline: boolean; readonly message: string };
}

export function includedRecords(revision: DatasetRevision): readonly DatasetRecord[] {
  const excluded = new Set(revision.excluded);
  return revision.records.filter((record) => !excluded.has(record.index));
}

/**
 * The name of a run input's record in the results: the ID that the run's outputs show for it where
 * the engine read the IDs of the deflines (`shown`, -parse_deflines: `AB000001.1` of
 * `>gb|AB000001.1| ...`; line WQA's sweep, D7, 2026-11-19), else its ID, or for a record without a
 * title the name of NCBI's tabular outputs, `Query_N` or `Subject_N` (N its number in the run input).
 */
export function recordName(id: string, role: 'query' | 'subject', index: number, shown?: ShownIds): string {
  const parsed = shown?.[role][index];
  if (parsed !== undefined) return parsed;
  return id !== '' ? id : `${role === 'query' ? 'Query' : 'Subject'}_${index + 1}`;
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
