// The records of an engine input, as the engine's index scan and `register` report them
// (docs/web/abi_v2.md §9). The fields are those of the ABI v2 responses, so no mapping layer is needed.
import type { InputRole, ProgramId } from './programs.js';

/**
 * The kind of FASTA reader of a program's input (plan TD-8): NCBI BLAST+'s reader with the
 * nucleotide flags (1) or the protein flags (2). Two readers of one kind read the same records
 * from the same bytes; only the words of their refusals differ (`FastaReader`).
 */
export type FastaParserKind = 1 | 2;

/**
 * The engine's reader of a program's input of a role (ABI v2 `scan_begin` and `register`, version 14):
 * the reader of the CLI for that input, with NCBI's data loaders on or off as the run's registry sets
 * them (`[BLAST] DATA_LOADERS` of `-conffile`; the *validation* response's `data_loaders`) and with
 * `-parse_deflines` on or off (`parse_deflines`). With the loaders on, the reader refuses a first line
 * that may be a Seq-id, in the CLI's words; with them off it reads that line as data (line WQA's D13
 * and D16, 2026-11-20). With `parseDeflines` on, the scan reads each record's Seq-id molecule and
 * letters as the engine's reader reads `-parse_deflines` (line RG-XT item X7-5, 2026-11-25).
 */
export interface FastaReader {
  readonly program: ProgramId;
  readonly role: InputRole;
  readonly dataLoaders: boolean;
  readonly parseDeflines: boolean;
}

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
  /** The number of residues that the reader stores, the letters of a `>?` gap line included. */
  readonly length: number;
  /** Count of each stored residue, upper case (kind 1: `U` as `T`): `A`-`Z` and `*`; a gap's letters are `N` (kind 1) or `X` (kind 2). */
  readonly residue_counts: Readonly<Record<string, number>>;
}

/** What the engine's `register` reports for each record (ABI v2 §9 *register*). */
export interface RecordKey {
  readonly id: string;
  readonly length: number;
}

/**
 * The IDs that the outputs of a run show for the records of its query and its subject, in the
 * order of the records, where the engine read them from the deflines (`-parse_deflines`: ABI v2 §9
 * *inputs*, `shown_id`): `AB000001.1` of `>gb|AB000001.1| ...`, the first two columns of outfmt 6
 * and 7. A run without them shows each record's ID (`RecordKey.id`). A run of `-db` has no subject records:
 * `subject` is empty, and a hit's `s_idx` is its subject's OID in the database, whose IDs the outputs show
 * (ABI v2 §8).
 */
export interface ShownIds {
  readonly query: readonly string[];
  readonly subject: readonly string[];
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
