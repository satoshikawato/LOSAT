// Run output channel (plan §3.1): the engine sends the outputs of one run straight to the
// Data worker over the MessagePort that `DataGateway.openRun` returns, so the outputs do
// not pass through the main thread. The engine side writes with RunOutputWriter
// (src/infra/run-output/writer.ts); web/app/tests/contract/run-output.contract.ts is the contract.
import type { RunOutput } from '../domain/output-format.js';

export const HITS_STREAM = 1;
export const DIAGNOSTICS_STREAM = 3;

/**
 * The ABI v2 `emit` streams of a run (docs/web/abi_v2.md §5): 0, 6 and 7 are the output
 * formats, 4 the output of the argv's -outfmt (another format, or 6 or 7 with fields; only such a
 * run writes it),
 * 1 the HSP records (JSON Lines) and 3 the diagnostics (UTF-8).
 */
export type OutputStream = RunOutput | typeof HITS_STREAM | typeof DIAGNOSTICS_STREAM;
export const OUTPUT_STREAMS: readonly OutputStream[] = Object.freeze([0, 6, 7, 4, HITS_STREAM, DIAGNOSTICS_STREAM]);

/**
 * Where the error of an output that a run could not write is, as the engine's signal of stream 9 says (docs/web/abi_v2.md
 * §5): `diagnostics`, the run's diagnostics (stream 3, the CLI's standard error); `logfile`, the file of the run's
 * -logfile (stream 8: TBLASTX writes the messages of its run there, as the CLI does); `nowhere`, where -logfile names no
 * file (TBLASTX's `-logfile ""` or `-logfile /dev/null`, which write them nowhere, as the CLI does). Line WQA's sweep,
 * D7, 2026-11-20.
 */
export type OutputErrorPlace = 'diagnostics' | 'logfile' | 'nowhere';
export const OUTPUT_ERROR_PLACES: readonly OutputErrorPlace[] = Object.freeze(['diagnostics', 'logfile', 'nowhere']);

/** Whether `value` is a place of the error of an output that a run could not write. */
export function isOutputErrorPlace(value: unknown): value is OutputErrorPlace {
  return OUTPUT_ERROR_PLACES.includes(value as OutputErrorPlace);
}

/**
 * The messages of one run, in order: any number of chunks and discards, then one `end` that
 * states how many chunks and bytes were sent. A `discard` says that an output failed and where its
 * error is (the engine's signal of stream 9, docs/web/abi_v2.md §5: it failed alone, or its search
 * failed in a run of two searches) and drops what arrived of it, which then stays empty; the totals
 * of `end` count every chunk sent, discarded or not. A run can be committed only after its `end`
 * arrived and the totals agree; a run that fails or is cancelled sends no `end` and is discarded.
 */
export type RunOutputMessage =
  | { readonly type: 'chunk'; readonly stream: OutputStream; readonly bytes: Uint8Array }
  | { readonly type: 'discard'; readonly stream: RunOutput; readonly error: OutputErrorPlace }
  | { readonly type: 'end'; readonly chunks: number; readonly bytes: number };

/** Whether `stream` is an output of a run (not the HSP records or the diagnostics): what a `discard` may name. */
export function isRunOutput(stream: unknown): stream is RunOutput {
  return OUTPUT_STREAMS.includes(stream as OutputStream) && stream !== HITS_STREAM && stream !== DIAGNOSTICS_STREAM;
}
