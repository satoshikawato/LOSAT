// The index scan of the adapter's serial reactor (ABI v2 `scan_begin` / `scan_chunk` /
// `scan_end`, docs/web/abi_v2.md §4, §9; plan TD-8): the RecordScanner of the Data worker.
// The response is the *scan* JSON as the adapter writes it, and a failure is the parser's
// error text.
import type { RecordScanner, ScanResponse } from '../../ports/scan';
import type { ReactorAbi } from './abi';

export class ReactorScanner implements RecordScanner {
  constructor(private readonly reactor: () => Promise<ReactorAbi>) {}

  async scan(parser: number, chunks: AsyncIterable<Uint8Array>): Promise<ScanResponse> {
    const abi = await this.reactor();
    const scanner = abi.scanBegin(parser);
    let ended = false;
    try {
      for await (const chunk of chunks) abi.scanChunk(scanner, chunk);
      ended = true;
      return abi.scanEnd(scanner) as ScanResponse;
    } finally {
      // A scan that stops early still frees its scanner in the engine.
      if (!ended && abi.stoppedBy === undefined) {
        try {
          abi.scanEnd(scanner);
        } catch {
          // The scan's own error is the one to report.
        }
      }
    }
  }
}

/**
 * Opens a reactor on first use, and again after it failed to open or stopped (a trap
 * leaves it unusable). The returned function gives the current one.
 */
export function reopening(open: () => Promise<ReactorAbi>): () => Promise<ReactorAbi> {
  let current: Promise<ReactorAbi> | undefined;
  return async () => {
    if (current === undefined) {
      current = open();
      return current;
    }
    const pending = current;
    const abi = await pending.catch(() => undefined);
    if (abi !== undefined && abi.stoppedBy === undefined) return abi;
    // Only the first caller that sees the failure opens a new one.
    if (current === pending) current = open();
    return current;
  };
}
