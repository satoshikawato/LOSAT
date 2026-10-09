// Scales and display classes of the results' plots (S13b; docs/web/ncbi_ui_mapping.md, "Dot Plot"
// and "Graphic Summary"). The dot plot's ticks follow the Owner's blast2dotplot.py (its tick_size
// table and its bp / kbp / Mbp labels), applied to the visible span so that zooming chooses them
// again. The classes only sort an engine value into the bins of a fixed table (the identity
// classes of the script, the "Alignment Scores" bins of NCBI's Graphic Summary): they compute no
// BLAST value, and the values shown next to them stay the engine's strings.

// --- ticks of the dot plot -------------------------------------------------------------------------

export interface TickSteps {
  readonly major: number;
  readonly minor: number;
}

/** The major and minor tick steps for a span of letters: the `tick_size` table of blast2dotplot.py. */
export function tickSteps(span: number): TickSteps {
  if (span <= 1000) return { major: 100, minor: 10 };
  if (span <= 5000) return { major: 500, minor: 100 };
  if (span <= 10000) return { major: 1000, minor: 500 };
  if (span <= 20000) return { major: 2000, minor: 1000 };
  if (span <= 50000) return { major: 5000, minor: 1000 };
  if (span <= 100000) return { major: 10000, minor: 1000 };
  if (span <= 1000000) return { major: 50000, minor: 10000 };
  if (span <= 5000000) return { major: 1000000, minor: 200000 };
  return { major: 1000000, minor: 500000 };
}

export interface AxisUnit {
  /** bp, kbp or Mbp on a nucleotide axis; aa, kaa or Maa on a protein axis. */
  readonly name: string;
  readonly divisor: number;
}

/**
 * The unit of an axis's labels by its visible span, as blast2dotplot.py's base_to_kbp chooses it
 * by the length (under 5000 letters, under 1,000,000, and above). `residue` is the run's unit of
 * the axis's records: "aa" for proteins, anything else ("nt") for nucleotides.
 */
export function axisUnit(span: number, residue: string): AxisUnit {
  const names = residue === 'aa' ? (['aa', 'kaa', 'Maa'] as const) : (['bp', 'kbp', 'Mbp'] as const);
  if (span < 5000) return { name: names[0], divisor: 1 };
  if (span < 1000000) return { name: names[1], divisor: 1000 };
  return { name: names[2], divisor: 1000000 };
}

/** A position in its axis's unit: whole values without a decimal point, others with the decimals they need. */
export function tickLabel(position: number, unit: AxisUnit): string {
  return (position / unit.divisor).toLocaleString('en-US', { maximumFractionDigits: 3 });
}

/** Minor ticks are left out when more than this many would be drawn. */
export const MAX_MINOR_TICKS = 50;

export interface AxisTicks {
  readonly unit: AxisUnit;
  readonly steps: TickSteps;
  /** Positions of the major ticks (the labelled ones), in letters. */
  readonly major: readonly number[];
  /** Positions of the minor ticks that are not major ticks; none when there would be too many. */
  readonly minor: readonly number[];
}

/** The ticks of an axis whose visible range is `from` to `to` (positions 0 to the record's length). */
export function axisTicks(from: number, to: number, residue: string): AxisTicks {
  const span = to - from;
  const steps = tickSteps(span);
  const unit = axisUnit(span, residue);
  if (!(span > 0)) return { unit, steps, major: [], minor: [] };
  const major = multiples(from, to, steps.major);
  const minorCount = count(from, to, steps.minor) - major.length;
  const minor = minorCount > MAX_MINOR_TICKS ? [] : multiples(from, to, steps.minor).filter((t) => t % steps.major !== 0);
  return { unit, steps, major, minor };
}

const count = (from: number, to: number, step: number): number => Math.max(0, Math.floor(to / step) - Math.ceil(from / step) + 1);

function multiples(from: number, to: number, step: number): number[] {
  const out: number[] = [];
  for (let k = Math.ceil(from / step), last = Math.floor(to / step); k <= last; k++) out.push(k * step);
  return out;
}

// --- the Graphic Summary's ruler -------------------------------------------------------------------

/**
 * Labelled positions of a ruler along a sequence of `length` letters: 1, the length, and round
 * positions (1, 2 or 5 times a power of ten) between them, at most about `labels` in all. Round
 * positions closer than half a step to either end are left out, so that labels do not collide.
 */
export function rulerTicks(length: number, labels: number): number[] {
  if (length <= 1) return [1];
  const raw = (length - 1) / Math.max(1, labels - 1);
  const power = 10 ** Math.floor(Math.log10(raw));
  const step = Math.max(1, [1, 2, 5, 10].map((m) => m * power).find((s) => s >= raw) ?? 10 * power);
  const out = [1];
  for (let p = step; p < length; p += step) if (p - 1 >= step / 2 && length - p >= step / 2) out.push(p);
  out.push(length);
  return out;
}

// --- display classes ---------------------------------------------------------------------------------

export interface IdentityClass {
  readonly label: string;
  readonly opacity: number;
}

/** blast2dotplot.py's opacity classes of the percent identity, its whole part (`int()`) compared. */
export const IDENTITY_CLASSES: readonly IdentityClass[] = Object.freeze([
  { label: '≤ 60', opacity: 0.4 },
  { label: '61–70', opacity: 0.6 },
  { label: '71–80', opacity: 0.8 },
  { label: '> 80', opacity: 1 },
]);

/**
 * The opacity class (an index of IDENTITY_CLASSES) of an outfmt 6 `pident` string as written: the
 * string read as a number and truncated to its whole part, as the script's `int()` does, then
 * ≤ 60, ≤ 70, ≤ 80 or above. Truncating the rounded `pident` (three decimals) gives the class that
 * the script's truncation gives. A string that is not a number (the FakeEngine's) has the last
 * class, drawn opaque.
 */
export function identityClass(pident: string): number {
  const value = pident.trim() === '' ? NaN : Math.trunc(Number(pident));
  if (!Number.isFinite(value)) return IDENTITY_CLASSES.length - 1;
  if (value <= 60) return 0;
  if (value <= 70) return 1;
  if (value <= 80) return 2;
  return 3;
}

export interface ScoreBin {
  readonly label: string;
  readonly color: string;
}

/** The "Alignment Scores" bins of NCBI's Graphic Summary, half-open: < 40, [40, 50), [50, 80), [80, 200), >= 200. */
export const SCORE_BINS: readonly ScoreBin[] = Object.freeze([
  { label: '< 40', color: '#000000' },
  { label: '40 - 50', color: '#0020e9' },
  { label: '50 - 80', color: '#75ea4c' },
  { label: '80 - 200', color: '#db3de9' },
  { label: '>= 200', color: '#db3324' },
]);

/** The bin (an index of SCORE_BINS) of an HSP record's bit score; a value that is not a number falls in the first. */
export function scoreBin(bitScore: number): number {
  if (!(bitScore >= 40)) return 0;
  if (bitScore < 50) return 1;
  if (bitScore < 80) return 2;
  if (bitScore < 200) return 3;
  return 4;
}
