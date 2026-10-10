// Index scan port: ABI v2 `scan_begin` / `scan_chunk` / `scan_end` (docs/web/abi_v2.md §4,
// §9; plan §5.4, TD-8). The Data worker builds the record table of a source with it, to
// extract original records later. It reads as the engine's reader does (kinds 1 and 2), so
// the table has the records that the search reads; the engine checks the table at
// `register`.
//
// The real implementation is the adapter's serial reactor inside the Data worker (S09).
// The FakeEngine build and the tests use the FakeScanner. Every implementation must pass
// tests/contract/record-scanner.contract.ts.
import type { FastaParserKind, IndexedRecord } from '../domain/dataset';

/** The ABI v2 *scan* response. */
export interface ScanResponse {
  readonly records: readonly IndexedRecord[];
}

export interface RecordScanner {
  /**
   * Scans one input, given as chunks of any size in order. Rejects with the reader's error
   * message (for example "CFastaReader: Near line 2, there's a line that doesn't look like
   * plausible data, ...") or LOSAT Web's rejection ("... not supported by LOSAT Web").
   */
  scan(parser: FastaParserKind, chunks: AsyncIterable<Uint8Array>): Promise<ScanResponse>;
}
