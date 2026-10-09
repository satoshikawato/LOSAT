// Geometry of the dot plot (S13b; docs/web/ncbi_ui_mapping.md, "Dot Plot"), after the Owner's
// blast2dotplot.py: X is the query and Y the subject, the origin at the top left with subject
// positions growing downward, and the same scale on both axes. Positions run from 0 to the
// record's length, and an HSP is drawn from (qstart, sstart) to (qend, send) of its record, as the
// script draws it. The view is the visible range; zoom keeps the two axes' spans in the ratio of
// the whole plot, so that the scale stays the same on both axes at every zoom. Where one axis is
// protein and the other nucleotide (TBLASTN, BLASTX), a residue counts as 3 letters for the
// proportions only (S13b decision 26), so that the HSPs run at 45°; positions stay the records'.
import type { SequenceKind } from './programs';

/** The lengths of the query (x) and the subject (y) records, in their letters. */
export interface Extent {
  readonly x: number;
  readonly y: number;
}

/** The visible range: query positions x0 to x1, subject positions y0 to y1. */
export interface View {
  readonly x0: number;
  readonly x1: number;
  readonly y0: number;
  readonly y1: number;
}

/** A rectangle in CSS pixels. */
export interface Box {
  readonly left: number;
  readonly top: number;
  readonly width: number;
  readonly height: number;
}

/** The shortest side of the plot in CSS pixels: a much shorter sequence is drawn this long, not to scale. */
export const MIN_SIDE = 120;
/** The longest side of the plot in CSS pixels (the script's figure size). */
export const MAX_SIDE = 1000;
/** The fewest letters an axis shows when zoomed in. */
export const MIN_SPAN = 4;

export const fullView = (extent: Extent): View => ({ x0: 0, x1: extent.x, y0: 0, y1: extent.y });

export interface PlotSize {
  readonly width: number;
  readonly height: number;
  /** Whether both axes have the same scale (false when the shorter side was raised to MIN_SIDE). */
  readonly toScale: boolean;
}

/** The letters that one position of each axis counts as for the plot's proportions. */
export interface Weights {
  readonly x: number;
  readonly y: number;
}

export const EQUAL_WEIGHTS: Weights = { x: 1, y: 1 };

/**
 * The weights of the axes of a query and a subject of these kinds: a protein axis against a
 * nucleotide one counts 3 per residue (a codon; TBLASTN's query, BLASTX's subject), else 1 each.
 */
export function codonWeights(query: SequenceKind, subject: SequenceKind): Weights {
  if (query === subject) return EQUAL_WEIGHTS;
  return query === 'protein' ? { x: 3, y: 1 } : { x: 1, y: 3 };
}

/**
 * The size of the plot area (inside the frame): the longer sequence (its length times its
 * weight) spans `available` pixels (at most MAX_SIDE), the shorter one in proportion, but at least
 * MIN_SIDE.
 */
export function plotSize(extent: Extent, available: number, weights: Weights = EQUAL_WEIGHTS): PlotSize {
  const longest = Math.max(1, Math.floor(Math.min(MAX_SIDE, available)));
  const [x, y] = [Math.max(extent.x, 0) * weights.x, Math.max(extent.y, 0) * weights.y];
  const longer = Math.max(x, y, 1);
  const side = (length: number) => (longest * length) / longer;
  const minimum = Math.min(MIN_SIDE, longest);
  const width = side(x);
  const height = side(y);
  return {
    width: Math.max(minimum, Math.round(width)),
    height: Math.max(minimum, Math.round(height)),
    toScale: width >= minimum && height >= minimum,
  };
}

// --- sequence positions and pixels -------------------------------------------------------------------

export const toPixelX = (x: number, view: View, box: Box): number => box.left + ((x - view.x0) / (view.x1 - view.x0)) * box.width;
export const toPixelY = (y: number, view: View, box: Box): number => box.top + ((y - view.y0) / (view.y1 - view.y0)) * box.height;
export const fromPixelX = (px: number, view: View, box: Box): number => view.x0 + ((px - box.left) / box.width) * (view.x1 - view.x0);
export const fromPixelY = (py: number, view: View, box: Box): number => view.y0 + ((py - box.top) / box.height) * (view.y1 - view.y0);

// --- zoom and pan ------------------------------------------------------------------------------------

/** The view's scale: the fraction of the query that it shows (1 for the whole plot). */
const scaleOf = (view: View, extent: Extent): number => (view.x1 - view.x0) / extent.x;

/** The smallest scale: each axis shows at least MIN_SPAN letters (or its whole record if shorter). */
const minScale = (extent: Extent): number => Math.max(Math.min(MIN_SPAN, extent.x) / extent.x, Math.min(MIN_SPAN, extent.y) / extent.y);

/** Moves a view (without changing its spans) so that it lies within the records. */
export function clampView(view: View, extent: Extent): View {
  const fit = (a: number, b: number, max: number): [number, number] => {
    const span = Math.min(b - a, max);
    const start = Math.min(Math.max(0, a), max - span);
    return [start, start + span];
  };
  const [x0, x1] = fit(view.x0, view.x1, extent.x);
  const [y0, y1] = fit(view.y0, view.y1, extent.y);
  return { x0, x1, y0, y1 };
}

/** A view of `scale` centred on a point (sequence positions), moved into the records. */
function viewAt(scale: number, cx: number, cy: number, extent: Extent): View {
  const s = Math.min(1, Math.max(minScale(extent), scale));
  const w = s * extent.x;
  const h = s * extent.y;
  return clampView({ x0: cx - w / 2, x1: cx + w / 2, y0: cy - h / 2, y1: cy + h / 2 }, extent);
}

/**
 * Zooms by `factor` (above 1 in, below 1 out) about a point given in sequence positions, which
 * stays where it is on the screen unless the view meets an end of a record. Both spans change by
 * the same factor, so the scale stays the same on both axes; the zoom stops at the whole plot and
 * at MIN_SPAN letters.
 */
export function zoomView(view: View, factor: number, at: { readonly x: number; readonly y: number }, extent: Extent): View {
  const scale = scaleOf(view, extent);
  const next = Math.min(1, Math.max(minScale(extent), scale / factor));
  const f = scale / next;
  const w = next * extent.x;
  const h = next * extent.y;
  const x0 = at.x - (at.x - view.x0) / f;
  const y0 = at.y - (at.y - view.y0) / f;
  return clampView({ x0, x1: x0 + w, y0, y1: y0 + h }, extent);
}

/** Moves the view by `dx` query and `dy` subject letters, within the records. */
export function panView(view: View, dx: number, dy: number, extent: Extent): View {
  return clampView({ x0: view.x0 + dx, x1: view.x1 + dx, y0: view.y0 + dy, y1: view.y1 + dy }, extent);
}

/**
 * The view that frames a range (an HSP from its start to its end on both sequences) with a margin
 * of a quarter of its length (at least 10 letters) on each side, at the same scale on both axes.
 */
export function frameRange(range: View, extent: Extent): View {
  const [qa, qb] = [Math.min(range.x0, range.x1), Math.max(range.x0, range.x1)];
  const [sa, sb] = [Math.min(range.y0, range.y1), Math.max(range.y0, range.y1)];
  const padX = Math.max(10, (qb - qa) * 0.25);
  const padY = Math.max(10, (sb - sa) * 0.25);
  const scale = Math.max((qb - qa + 2 * padX) / extent.x, (sb - sa + 2 * padY) / extent.y);
  return viewAt(scale, (qa + qb) / 2, (sa + sb) / 2, extent);
}

// --- picking and placing -----------------------------------------------------------------------------

/** The distance in pixels from a point to the segment from (ax, ay) to (bx, by). */
export function segmentDistance(x: number, y: number, ax: number, ay: number, bx: number, by: number): number {
  const lx = bx - ax;
  const ly = by - ay;
  const length = lx * lx + ly * ly;
  const t = length === 0 ? 0 : Math.max(0, Math.min(1, ((x - ax) * lx + (y - ay) * ly) / length));
  return Math.hypot(x - (ax + t * lx), y - (ay + t * ly));
}

/**
 * Whether a line from (ax, ay) to (bx, by) can show in the view: its bounding box, widened by
 * `padX` and `padY` (the reach of its stroke beyond the line, in the axes' positions), meets the
 * view. A line that is left out leaves no pixel in the frame.
 */
export function lineNearView(ax: number, ay: number, bx: number, by: number, view: View, padX: number, padY: number): boolean {
  return (
    (ax > bx ? ax : bx) >= view.x0 - padX &&
    (ax < bx ? ax : bx) <= view.x1 + padX &&
    (ay > by ? ay : by) >= view.y0 - padY &&
    (ay < by ? ay : by) <= view.y1 + padY
  );
}

/** The part of the segment from (ax, ay) to (bx, by) inside a box (Liang–Barsky), or undefined if none. */
export function clipToBox(ax: number, ay: number, bx: number, by: number, box: Box): [number, number, number, number] | undefined {
  const [dx, dy] = [bx - ax, by - ay];
  let [t0, t1] = [0, 1];
  // For each edge: the segment's motion towards the outside, and the room inside at its start.
  const edges: readonly (readonly [number, number])[] = [
    [-dx, ax - box.left],
    [dx, box.left + box.width - ax],
    [-dy, ay - box.top],
    [dy, box.top + box.height - ay],
  ];
  for (const [p, q] of edges) {
    if (p === 0) {
      if (q < 0) return undefined;
      continue;
    }
    const r = q / p;
    if (p < 0) {
      if (r > t1) return undefined;
      if (r > t0) t0 = r;
    } else {
      if (r < t0) return undefined;
      if (r < t1) t1 = r;
    }
  }
  return [ax + t0 * dx, ay + t0 * dy, ax + t1 * dx, ay + t1 * dy];
}

/**
 * The line nearest to a point within `reach` pixels, or -1: lines `0 .. count - 1` with their
 * ends in pixels (`ends`). Only what the frame `box` shows counts: a point outside the frame picks
 * nothing, and a line is measured by its part inside the frame, so that a line drawn outside the
 * view (or its part outside it) is never picked. Of lines at the same distance, the last wins.
 */
export function pickLine(
  x: number,
  y: number,
  count: number,
  ends: (i: number) => readonly [number, number, number, number],
  box: Box,
  reach: number,
): number {
  if (x < box.left || x > box.left + box.width || y < box.top || y > box.top + box.height) return -1;
  let best = -1;
  let bestDistance = reach;
  for (let i = 0; i < count; i++) {
    const [ax, ay, bx, by] = ends(i);
    if (Math.max(ax, bx) < x - reach || Math.min(ax, bx) > x + reach || Math.max(ay, by) < y - reach || Math.min(ay, by) > y + reach) continue;
    const seen = clipToBox(ax, ay, bx, by, box);
    if (seen === undefined) continue;
    const d = segmentDistance(x, y, ...seen);
    if (d <= bestDistance) {
      best = i;
      bestDistance = d;
    }
  }
  return best;
}

/**
 * Where to put a box of `size` next to an anchor point within `bounds` (both in CSS pixels from
 * the bounds' top left): to the right of and below the point, `gap` pixels away; to the left or
 * above where it would leave the bounds; and pushed inside the bounds where neither side has room.
 */
export function placeBox(
  anchor: { readonly x: number; readonly y: number },
  size: { readonly width: number; readonly height: number },
  bounds: { readonly width: number; readonly height: number },
  gap = 12,
): { left: number; top: number } {
  const side = (at: number, length: number, room: number): number => {
    if (at + gap + length <= room) return at + gap;
    if (at - gap - length >= 0) return at - gap - length;
    return Math.max(0, Math.min(room - length, at - length / 2));
  };
  return { left: Math.round(side(anchor.x, size.width, bounds.width)), top: Math.round(side(anchor.y, size.height, bounds.height)) };
}

/**
 * Where to put a box of `size` beside a line from (ax, ay) to (bx, by), within `bounds` (CSS
 * pixels from the bounds' top left): the dot plot's popup of the selected HSP, which hid half of
 * the line when it was put at the line's midpoint (W4b screen review L10). The box goes `gap`
 * pixels off one of the line's ends, towards one of the four corners; of these eight places, moved
 * inside the bounds, the one that covers the least of the line (with a margin of half the gap) and
 * then needs the least move. None that covers the line's midpoint is taken: undefined where every
 * one does, or where the box is larger than the bounds (the dot plot then shows it below itself).
 */
export function placeBeside(
  line: readonly [number, number, number, number],
  size: { readonly width: number; readonly height: number },
  bounds: { readonly width: number; readonly height: number },
  gap = 12,
): { left: number; top: number } | undefined {
  if (size.width > bounds.width || size.height > bounds.height) return undefined;
  const [ax, ay, bx, by] = line;
  const [mx, my] = [(ax + bx) / 2, (ay + by) / 2];
  const margin = gap / 2;
  let best: { left: number; top: number } | undefined;
  let bestCost = Infinity;
  for (const [ex, ey] of [
    [ax, ay],
    [bx, by],
  ] as const) {
    for (const [dx, dy] of [
      [1, 1],
      [1, -1],
      [-1, 1],
      [-1, -1],
    ] as const) {
      const wanted = { left: dx > 0 ? ex + gap : ex - gap - size.width, top: dy > 0 ? ey + gap : ey - gap - size.height };
      const left = Math.max(0, Math.min(bounds.width - size.width, wanted.left));
      const top = Math.max(0, Math.min(bounds.height - size.height, wanted.top));
      const area: Box = { left: left - margin, top: top - margin, width: size.width + 2 * margin, height: size.height + 2 * margin };
      if (mx >= area.left && mx <= area.left + area.width && my >= area.top && my <= area.top + area.height) continue;
      const covered = clipToBox(ax, ay, bx, by, area);
      const cost = (covered === undefined ? 0 : 1000 * Math.hypot(covered[2] - covered[0], covered[3] - covered[1])) + Math.abs(left - wanted.left) + Math.abs(top - wanted.top);
      if (cost < bestCost) {
        bestCost = cost;
        best = { left: Math.round(left), top: Math.round(top) };
      }
    }
  }
  return best;
}
