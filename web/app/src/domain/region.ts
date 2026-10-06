// The region of an input (plan §5.3, DW-9): one range on the only record of a role, which
// becomes `-query_loc start-stop` or `-subject_loc start-stop`. Positions are 1-based and
// inclusive, in the record's letters (nucleotides or residues). The form limits both
// positions to the record (1 to its length), as the preview's selection does (S12
// instructions, item 3); every other rule (for example start < stop) is the engine's,
// whose `validate` gives NCBI's message.
import type { InputRole } from './programs';

/** The two fields as the user typed them. */
export interface RegionText {
  readonly start: string;
  readonly stop: string;
}

export const REGION_FLAG: Readonly<Record<InputRole, string>> = { query: '-query_loc', subject: '-subject_loc' };

const DIGITS = /^[0-9]+$/;

/** Why the fields do not make a range within a record of `length` letters, or undefined. */
export function regionProblem(region: RegionText, length: number): string | undefined {
  const start = region.start.trim();
  const stop = region.stop.trim();
  if (start === '' || stop === '') return 'Enter both the start and the stop.';
  if (!DIGITS.test(start) || !DIGITS.test(stop)) return 'The start and the stop are whole numbers.';
  for (const [name, text] of [
    ['start', start],
    ['stop', stop],
  ] as const) {
    const value = Number(text);
    if (value < 1 || value > length) return `The ${name} must be between 1 and ${length}, the length of the record.`;
  }
  return undefined;
}

/** The option's value, `start-stop`, without white space or leading zeros. */
export function regionValue(region: RegionText): string {
  return `${Number(region.start.trim())}-${Number(region.stop.trim())}`;
}

/** The region as a range of positions, for the preview, or undefined while it is incomplete. */
export function regionRange(region: RegionText, length: number): { readonly start: number; readonly stop: number } | undefined {
  if (regionProblem(region, length) !== undefined) return undefined;
  return { start: Number(region.start.trim()), stop: Number(region.stop.trim()) };
}

/** A region from two positions of the preview, ordered and limited to the record. */
export function regionFromPositions(a: number, b: number, length: number): RegionText {
  const clamp = (value: number) => Math.min(length, Math.max(1, Math.round(value)));
  const [start, stop] = [clamp(a), clamp(b)].sort((x, y) => x - y) as [number, number];
  return { start: String(start), stop: String(stop) };
}
