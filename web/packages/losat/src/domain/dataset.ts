// The records of an engine input, as the engine's index scan and `register` report them
// (docs/web/abi_v2.md §9). The fields are those of the ABI v2 responses, so no mapping layer is needed.

/**
 * The FASTA reader whose rules the index scan follows (ABI v2 `scan_begin`, plan TD-8): NCBI
 * BLAST+'s reader with the nucleotide flags (1) or the protein flags (2).
 */
export type FastaParserKind = 1 | 2;

/** How to find the byte offset of residue `i` of a record (ABI v2 §9 `line_layout`). */
export type LineLayout =
  | { readonly kind: 'uniform'; readonly width: number; readonly eol: number }
  | { readonly kind: 'checkpoints'; readonly every: number; readonly offsets: readonly number[] };

/** One record of the index scan. Offsets are byte offsets in the source and stay below 2^53. */
export interface IndexedRecord {
  readonly index: number;
  /** The title up to its first space; empty for a record without a title. */
  readonly id: string;
  /** Offset of the defline's `>`; 0 for a first record without a defline. */
  readonly header_offset: number;
  /** Offset of the first byte after the defline; 0 for a first record without a defline. */
  readonly sequence_offset: number;
  /** Offset of the next defline, or the end of the input. */
  readonly end_offset: number;
  /** The number of residues that the reader stores. */
  readonly length: number;
  readonly line_layout: LineLayout;
  /** Count of each stored residue, upper case (kind 1: `U` as `T`): `A`-`Z` and `*`. */
  readonly residue_counts: Readonly<Record<string, number>>;
}

/** What the engine's `register` reports for each record (ABI v2 §9 *register*). */
export interface RecordKey {
  readonly id: string;
  readonly length: number;
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
