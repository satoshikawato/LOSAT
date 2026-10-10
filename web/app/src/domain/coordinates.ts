// Coordinates of records and HSPs (design §11.3, plan §5.7): the one place that turns the
// coordinates of an HSP record into intervals of the original record, with their strand and
// unit, and that applies flanks and the record's ends. The HSP record's coordinates are those
// that outfmt 6 writes (docs/web/abi_v2.md §8): 1-based, both ends included, start > end on the
// minus strand, in the letters of the record - the nucleotides of a translated sequence (the
// BLASTX query, the TBLASTN subject, both of TBLASTX; docs/web/results_columns.md), so a frame
// changes no coordinate here. Nothing here computes a BLAST value. Drawing transforms (the dot
// plot's view, plot-geometry.ts) stay with the drawing; no other module turns an HSP's start and
// end into an interval or a strand.
import type { InputRole, ProgramId, SequenceKind } from './programs';

/** The unit of a record's letters and of every coordinate on it. */
export type Unit = 'nt' | 'aa';

/** How an HSP runs along a sequence; unknown where the record does not say (a BLASTN HSP of one letter). */
export type Strand = 'plus' | 'minus' | 'unknown';

/** An interval of a record: 1-based positions, both ends included, `from` <= `to`, on the record's forward strand. */
export interface Interval {
  readonly from: number;
  readonly to: number;
}

/** The interval of an HSP on one of its sequences, with the strand that its coordinates give. */
export interface HspSpan extends Interval {
  readonly strand: Strand;
}

/** The coordinates and frames of an HSP record (docs/web/abi_v2.md §8; `HspRecord`, `HspSummary`). */
export interface HspCoordinates {
  readonly q_start: number;
  readonly q_end: number;
  readonly s_start: number;
  readonly s_end: number;
  readonly query_frame: number | null;
  readonly subject_frame: number | null;
}

export const unitOf = (kind: SequenceKind): Unit => (kind === 'nucleotide' ? 'nt' : 'aa');

export const intervalLength = (interval: Interval): number => interval.to - interval.from + 1;

/** "101-250": the record's positions, the smaller first. */
export const intervalText = (interval: Interval): string => `${interval.from}-${interval.to}`;

/** The positions of an HSP record on one sequence, the smaller first; the strand is dropped. */
export function interval(start: number, end: number): Interval {
  return start <= end ? { from: start, to: end } : { from: end, to: start };
}

/**
 * Whether the program translates the sequences of this kind, so that its HSP records' frames
 * mean something: the nucleotide sequences of every program but BLASTN (docs/web/results_columns.md
 * "Frames"). The engine's BLASTP records carry frame 1 for both sequences, which is not a frame.
 */
export const translates = (program: ProgramId, kind: SequenceKind): boolean => kind === 'nucleotide' && program !== 'blastn';

/**
 * The strand of an HSP on a sequence: the coordinates tell it (start > end is the minus strand);
 * for one letter (start = end) the frame's sign tells it, and a protein runs forward. A BLASTN
 * HSP of one letter has neither, so its strand is unknown from the record (only the `Strand=`
 * line of its outfmt 0 section shows it). A frame of 0, null or undefined is no frame.
 */
export function strandOf(start: number, end: number, frame: number | null | undefined, kind: SequenceKind): Strand {
  if (start !== end) return start < end ? 'plus' : 'minus';
  if (frame !== undefined && frame !== null && frame !== 0) return frame > 0 ? 'plus' : 'minus';
  return kind === 'protein' ? 'plus' : 'unknown';
}

/** The span of an HSP on one sequence from its record's coordinates. */
export function hspSpan(start: number, end: number, frame: number | null | undefined, kind: SequenceKind): HspSpan {
  return { ...interval(start, end), strand: strandOf(start, end, frame, kind) };
}

/** The span of an HSP record on its query or its subject, a sequence of `kind`. */
export function spanOn(hsp: HspCoordinates, role: InputRole, kind: SequenceKind): HspSpan {
  return role === 'query'
    ? hspSpan(hsp.q_start, hsp.q_end, hsp.query_frame, kind)
    : hspSpan(hsp.s_start, hsp.s_end, hsp.subject_frame, kind);
}

/**
 * Letters added on each side of an interval, on the record's axis whatever the strand: `left`
 * toward position 1, `right` toward the record's end (design §11.4).
 */
export interface Flanks {
  readonly left: number;
  readonly right: number;
}

export const NO_FLANKS: Flanks = Object.freeze({ left: 0, right: 0 });

/** An interval cut to its record: what was asked for, and what the record has (design §11.4). */
export interface ClippedInterval {
  readonly requested: Interval;
  /** `requested` within 1 and the record's length. */
  readonly actual: Interval;
  readonly clippedLeft: boolean;
  readonly clippedRight: boolean;
}

export const isClipped = (clipped: ClippedInterval): boolean => clipped.clippedLeft || clipped.clippedRight;

/** Whether an interval is whole positions within a record of `length` letters. */
export function withinRecord(interval: Interval, length: number): boolean {
  const { from, to } = interval;
  return Number.isSafeInteger(from) && Number.isSafeInteger(to) && from >= 1 && from <= to && to <= length;
}

/**
 * Cuts an interval to a record of `length` letters. The interval must reach into the record
 * (an HSP's interval always does); an interval outside it is an error of the caller's data.
 */
export function clipToRecord(requested: Interval, length: number): ClippedInterval {
  if (!Number.isSafeInteger(length) || length < 1) throw new RangeError(`a record of ${length} letters has no interval`);
  if (!Number.isSafeInteger(requested.from) || !Number.isSafeInteger(requested.to) || requested.from > requested.to) {
    throw new RangeError(`${requested.from}-${requested.to} is not an interval`);
  }
  if (requested.to < 1 || requested.from > length) {
    throw new RangeError(`the interval ${intervalText(requested)} is outside the record of ${length} letters`);
  }
  return {
    requested: { from: requested.from, to: requested.to },
    actual: { from: Math.max(1, requested.from), to: Math.min(length, requested.to) },
    clippedLeft: requested.from < 1,
    clippedRight: requested.to > length,
  };
}

/** The interval with its flanks, cut to the record; `requested` keeps the flanks that the record does not have. */
export function flanked(interval: Interval, flanks: Flanks, length: number): ClippedInterval {
  for (const side of ['left', 'right'] as const) {
    const value = flanks[side];
    if (!Number.isSafeInteger(value) || value < 0) throw new RangeError(`the ${side} flank must be a whole number of letters, not ${value}`);
  }
  return clipToRecord({ from: interval.from - flanks.left, to: interval.to + flanks.right }, length);
}

/**
 * The one interval from the smallest to the largest position of several (design §11.4's "one
 * region that includes what lies between the HSPs"). Overlapping or touching intervals are never
 * joined otherwise: the caller chooses this over separate intervals.
 */
export function spanning(intervals: readonly Interval[]): Interval {
  if (intervals.length === 0) throw new RangeError('no interval to span');
  let from = Number.POSITIVE_INFINITY;
  let to = Number.NEGATIVE_INFINITY;
  for (const each of intervals) {
    from = Math.min(from, each.from);
    to = Math.max(to, each.to);
  }
  return { from, to };
}

/** The strand that several spans share, or 'mixed' where they differ. */
export function commonStrand(spans: readonly HspSpan[]): Strand | 'mixed' {
  if (spans.length === 0) throw new RangeError('no span');
  const strands = new Set(spans.map((span) => span.strand));
  return strands.size === 1 ? spans[0]!.strand : 'mixed';
}
