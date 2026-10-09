// Geometry of the dot plot (S13b; docs/web/ncbi_ui_mapping.md, "Dot Plot"), after the Owner's
// blast2dotplot.py: X is the query and Y the subject, the origin at the top left with subject
// positions growing downward, and the same scale on both axes. Positions run from 0 to the
// record's length, and an HSP is drawn from (qstart, sstart) to (qend, send) of its record, as the
// script draws it. The view is the visible range; zoom keeps the two axes' spans in the ratio of
// the whole plot, so that the scale stays the same on both axes at every zoom.

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

/**
 * The size of the plot area (inside the frame): the longer sequence spans `available` pixels (at
 * most MAX_SIDE), the shorter one in proportion, but at least MIN_SIDE.
 */
export function plotSize(extent: Extent, available: number): PlotSize {
  const longest = Math.max(1, Math.floor(Math.min(MAX_SIDE, available)));
  const longer = Math.max(extent.x, extent.y, 1);
  const side = (length: number) => (longest * Math.max(length, 0)) / longer;
  const minimum = Math.min(MIN_SIDE, longest);
  const width = side(extent.x);
  const height = side(extent.y);
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
