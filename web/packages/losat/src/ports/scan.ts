// Index scan port: ABI v2 `scan_begin` / `scan_chunk` / `scan_end` (docs/web/abi_v2.md §4,
// §9; plan §5.4, TD-8). The Data worker builds the record table of a source with it, to
// extract original records later. The search itself always reads its input with the
// program's own parser, and the engine checks the table at `register`.
//
// The real implementation is the adapter's serial reactor inside the Data worker (S09);
// the development build and the tests use the FakeScanner. Every implementation must pass
// web/app/tests/contract/record-scanner.contract.ts.
import type { FastaReader, IndexedRecord } from '../domain/dataset.js';

/**
 * The reader's own error on the record after the records of a scan (NCBI's `CheckDataLine`): the
 * reading stops there, as the command line's query batches stop (line WQA's D10, 2026-11-19).
 */
export interface ReaderStop {
  /** The offset of the record that the reader fails on (its `>`; 0 for a first record without a defline). */
  readonly offset: number;
  /**
   * Whether that record starts with a defline: `false` for the input's first record read without one
   * (a byte-order mark, a first line of other text), which cannot follow another input's records
   * (review M2, 2026-11-20).
   */
  readonly defline: boolean;
  /** The reader's message (`CFastaReader: Near line …`). */
  readonly message: string;
}

/** The ABI v2 *scan* response. */
export interface ScanResponse {
  readonly records: readonly IndexedRecord[];
  /** The reader's error that ends the records, if one does. */
  readonly reader_error?: ReaderStop;
}

export interface RecordScanner {
  /**
   * Scans one input with the engine's reader of it (`reader`), given as chunks of any size in
   * order. The reader's own error ends the records (`reader_error`); its refusal in the CLI's words
   * (a first line that may be a Seq-id with the data loaders on, a record longer than LOSAT reads)
   * rejects the scan.
   */
  scan(reader: FastaReader, chunks: AsyncIterable<Uint8Array>): Promise<ScanResponse>;
}
