// Run output channel (plan §3.1): the engine sends the outputs of one run straight to the
// Data worker over the MessagePort that `DataGateway.openRun` returns, so the outputs do
// not pass through the main thread. The engine side writes with RunOutputWriter
// (src/infra/run-output/writer.ts); tests/contract/run-output.contract.ts is the contract.
import type { OutputFormat } from '../domain/output-format';

export const HITS_STREAM = 1;
export const DIAGNOSTICS_STREAM = 3;

/**
 * The ABI v2 `emit` streams of a run (docs/web/abi_v2.md §5): 0, 6 and 7 are the output
 * formats, 1 the HSP records (JSON Lines) and 3 the diagnostics (UTF-8).
 */
export type OutputStream = OutputFormat | typeof HITS_STREAM | typeof DIAGNOSTICS_STREAM;
export const OUTPUT_STREAMS: readonly OutputStream[] = Object.freeze([0, 6, 7, HITS_STREAM, DIAGNOSTICS_STREAM]);

/**
 * The messages of one run, in order: any number of chunks, then one `end` that states how
 * many chunks and bytes were sent. A run can be committed only after its `end` arrived and
 * the totals agree; a run that fails or is cancelled sends no `end` and is discarded.
 */
export type RunOutputMessage =
  | { readonly type: 'chunk'; readonly stream: OutputStream; readonly bytes: Uint8Array }
  | { readonly type: 'end'; readonly chunks: number; readonly bytes: number };
