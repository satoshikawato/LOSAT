// Runtime choices of the Engine worker gateway (plan §4.7, §5.5). They choose how a
// search runs, never what it computes: the thread count does not change the output
// (LOSAT reduces parallel work back into the NCBI order), and a renewed instance gives
// the same bytes as a warm one (plan §4.6 gate). The values come from the S09
// measurements in the browsers (docs/evidence/losat_web_w1/README.md, "Auto" and
// "Memory").

/**
 * Auto runs the serial module when the query and subject FASTA together are smaller than
 * this. gbdraw's default (500,000 characters) was the starting point. In S09, 4 threads
 * made the search phase about three times shorter wherever the work spreads over many
 * records (BLASTP of a proteome against itself, from 20,000 residues on), and changed it
 * little in absolute terms where it does not (one genome against another: at most 69 ms
 * longer); below this size every search took less than half a second.
 */
export const AUTO_SERIAL_BELOW_BYTES = 20_000;

/** Auto uses at most this many threads: the counts that V-ABI and V-BR check are 1, 2 and 4. */
export const AUTO_MAX_THREADS = 4;

/** The threads of a search: an explicit request as it is, Auto from the input size. */
export function chooseThreads(requested: number | 'auto', inputBytes: number, hardwareConcurrency: number): number {
  if (requested !== 'auto') return requested;
  if (inputBytes < AUTO_SERIAL_BELOW_BYTES) return 1;
  const hardware = Number.isInteger(hardwareConcurrency) && hardwareConcurrency > 0 ? hardwareConcurrency : 2;
  return Math.max(1, Math.min(AUTO_MAX_THREADS, Math.floor(hardware / 2)));
}

/**
 * When the Engine worker is renewed before the next search (plan §5.5, G12): the reactor's
 * linear memory never shrinks, and repeated calls could make it grow
 * (docs/wasm_reactor_memory_followup_20260914.md).
 */
export interface RenewalLimits {
  /** Searches of one instance. */
  readonly maxRuns: number;
  /** Linear memory of the instance after a search, in bytes. */
  readonly highWaterBytes: number;
}

/**
 * In S09 the linear memory reached its level in the first one to three searches and then
 * stayed there, over 20 repeats of each program and 30 searches of every program in turn;
 * the largest search measured (two E. coli genomes) left about 240 MiB. So the count of searches
 * renews nothing; the high-water mark, half the threaded module's maximum (1 GiB, plan
 * TD-7), gives the memory of an unusually large search back to the browser.
 */
export const DEFAULT_RENEWAL: RenewalLimits = Object.freeze({
  maxRuns: Number.POSITIVE_INFINITY,
  highWaterBytes: 512 * 1024 * 1024,
});

/** Why the runtime must be renewed after a search, or undefined. */
export function renewalReason(
  memory: { readonly linearBytesAfter: number; readonly instanceRuns: number },
  limits: RenewalLimits,
): string | undefined {
  if (memory.linearBytesAfter >= limits.highWaterBytes) {
    return `the engine memory reached ${memory.linearBytesAfter} bytes (limit ${limits.highWaterBytes})`;
  }
  if (memory.instanceRuns >= limits.maxRuns) return `the engine instance ran ${memory.instanceRuns} searches`;
  return undefined;
}
