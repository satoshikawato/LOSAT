// The dot plot's geometry (src/domain/plot-geometry.ts): the same scale on both axes (3 letters a
// residue on a protein axis against a nucleotide one), the mapping between positions and pixels
// with the origin at the top left, zoom and pan within the records, the lines that can show in the
// view and the line under a point, and the placing of the HSP popup.
import { describe, expect, it } from 'vitest';
import {
  clampView,
  clipToBox,
  codonWeights,
  EQUAL_WEIGHTS,
  frameRange,
  fromPixelX,
  fromPixelY,
  fullView,
  lineNearView,
  MIN_SIDE,
  panView,
  pickLine,
  placeBeside,
  placeBox,
  plotSize,
  segmentDistance,
  toPixelX,
  toPixelY,
  zoomView,
  type View,
} from '../../src/domain/plot-geometry';

const ratio = (view: View) => (view.x1 - view.x0) / (view.y1 - view.y0);

describe('plotSize', () => {
  it('spans the available width with the longer sequence, the shorter one in proportion', () => {
    expect(plotSize({ x: 2000, y: 1000 }, 800)).toEqual({ width: 800, height: 400, toScale: true });
    // The subject is the longer: it spans the available pixels downward.
    expect(plotSize({ x: 500, y: 1000 }, 600)).toEqual({ width: 300, height: 600, toScale: true });
  });

  it('is at most 1000 px and the shorter side at least 120 px (then not to scale)', () => {
    expect(plotSize({ x: 5000, y: 5000 }, 1400)).toEqual({ width: 1000, height: 1000, toScale: true });
    const thin = plotSize({ x: 100000, y: 500 }, 900);
    expect(thin).toEqual({ width: 900, height: MIN_SIDE, toScale: false });
  });

  it('counts 3 letters a residue on a protein axis against a nucleotide one (TBLASTN, BLASTX; decision 26)', () => {
    expect(codonWeights('protein', 'nucleotide')).toEqual({ x: 3, y: 1 });
    expect(codonWeights('nucleotide', 'protein')).toEqual({ x: 1, y: 3 });
    expect(codonWeights('nucleotide', 'nucleotide')).toEqual(EQUAL_WEIGHTS);
    expect(codonWeights('protein', 'protein')).toEqual(EQUAL_WEIGHTS);
    // TBLASTN: a 300 aa query against 3,000 nt is 900 against 3,000, to scale (not 60 px raised to 120).
    expect(plotSize({ x: 300, y: 3000 }, 600)).toEqual({ width: MIN_SIDE, height: 600, toScale: false });
    expect(plotSize({ x: 300, y: 3000 }, 600, codonWeights('protein', 'nucleotide'))).toEqual({ width: 180, height: 600, toScale: true });
    // BLASTX: 900 nt against 300 aa is square.
    expect(plotSize({ x: 900, y: 300 }, 600, codonWeights('nucleotide', 'protein'))).toEqual({ width: 600, height: 600, toScale: true });
    // A short protein against a genome still meets the 120 px minimum, and says so.
    expect(plotSize({ x: 100, y: 1_000_000 }, 800, codonWeights('protein', 'nucleotide'))).toEqual({ width: MIN_SIDE, height: 800, toScale: false });
  });
});

describe('positions and pixels', () => {
  const view = { x0: 100, x1: 300, y0: 0, y1: 50 };
  const box = { left: 50, top: 40, width: 400, height: 100 };

  it('maps query positions to x and subject positions downward to y, both ways', () => {
    expect(toPixelX(100, view, box)).toBe(50);
    expect(toPixelX(300, view, box)).toBe(450);
    expect(toPixelY(0, view, box)).toBe(40);
    expect(toPixelY(50, view, box)).toBe(140);
    expect(fromPixelX(toPixelX(217, view, box), view, box)).toBeCloseTo(217);
    expect(fromPixelY(toPixelY(12.5, view, box), view, box)).toBeCloseTo(12.5);
  });

  it('measures the distance from a point to a segment', () => {
    expect(segmentDistance(5, 5, 0, 0, 10, 0)).toBe(5);
    expect(segmentDistance(-3, 4, 0, 0, 10, 0)).toBe(5);
    expect(segmentDistance(2, 2, 2, 2, 2, 2)).toBe(0);
  });
});

describe('zoom and pan', () => {
  const extent = { x: 60, y: 180 };
  const full = fullView(extent);

  it('zooms about a point, keeping the two spans in the plot\'s ratio, and back to the whole plot', () => {
    const zoomed = zoomView(full, 1.5, { x: 30, y: 90 }, extent);
    expect(zoomed.x1 - zoomed.x0).toBeCloseTo(40);
    expect(ratio(zoomed)).toBeCloseTo(ratio(full));
    expect(zoomed.x0).toBeCloseTo(10);
    const back = zoomView(zoomed, 1 / 1.5, { x: 30, y: 90 }, extent);
    expect([back.x0, back.x1, back.y0, back.y1].map(Math.round)).toEqual([0, 60, 0, 180]);
  });

  it('keeps the point under the pointer in place', () => {
    const zoomed = zoomView(full, 2, { x: 15, y: 45 }, extent);
    // 15 was a quarter of the way along x; it still is.
    expect((15 - zoomed.x0) / (zoomed.x1 - zoomed.x0)).toBeCloseTo(0.25);
  });

  it('stops at the whole plot and at 4 letters', () => {
    expect(zoomView(full, 0.5, { x: 30, y: 90 }, extent)).toEqual(full);
    let view = full;
    for (let i = 0; i < 40; i++) view = zoomView(view, 2, { x: 30, y: 90 }, extent);
    expect(view.x1 - view.x0).toBeCloseTo(4);
    expect(ratio(view)).toBeCloseTo(ratio(full));
  });

  it('pans within the records', () => {
    const zoomed = zoomView(full, 2, { x: 30, y: 90 }, extent);
    expect(panView(zoomed, -100, 0, extent).x0).toBe(0);
    const moved = panView(zoomed, 5, 1000, extent);
    expect(moved.x0).toBeCloseTo(zoomed.x0 + 5);
    expect(moved.y1).toBe(180);
    expect(clampView({ x0: -5, x1: 25, y0: 170, y1: 260 }, extent)).toEqual({ x0: 0, x1: 30, y0: 90, y1: 180 });
  });

  it('frames an HSP (reverse ones too) at the same scale on both axes', () => {
    const big = { x: 10000, y: 20000 };
    const view = frameRange({ x0: 5200, x1: 5000, y0: 900, y1: 1000 }, big);
    expect(view.x0).toBeLessThanOrEqual(5000);
    expect(view.x1).toBeGreaterThanOrEqual(5200);
    expect(view.y0).toBeLessThanOrEqual(900);
    expect(view.y1).toBeGreaterThanOrEqual(1000);
    expect(ratio(view)).toBeCloseTo(0.5);
    // At an end of the records, the view moves inside and still holds the HSP.
    const edge = frameRange({ x0: 1, x1: 30, y0: 170, y1: 180 }, extent);
    expect(edge.x0).toBe(0);
    expect(edge.y1).toBe(180);
    expect(edge.x1).toBeGreaterThanOrEqual(30);
  });
});

describe('the lines that show in the view, and the line under a point', () => {
  it('keeps a line whose stroke reaches into the view, and leaves out one that cannot', () => {
    const view = { x0: 1000, x1: 2000, y0: 0, y1: 1000 };
    // 0.002 px a letter (a 1 Mbp pair zoomed 2x): a stroke of 1.5 px reaches 750 letters.
    const pad = 1.5 / 0.002;
    expect(lineNearView(100, 10, 600, 20, view, 0, 0)).toBe(false);
    expect(lineNearView(100, 10, 600, 20, view, pad, pad)).toBe(true);
    expect(lineNearView(600, 20, 100, 10, view, pad, pad)).toBe(true);
    expect(lineNearView(0, 10, 200, 20, view, pad, pad)).toBe(false);
    expect(lineNearView(1500, 1200, 1600, 1700, view, pad, pad)).toBe(true);
    expect(lineNearView(1500, 1800, 1600, 1900, view, pad, pad)).toBe(false);
    expect(lineNearView(1200, 300, 1300, 400, view, 0, 0)).toBe(true);
  });

  const box = { left: 10, top: 10, width: 100, height: 100 };

  it('clips a segment to a box', () => {
    expect(clipToBox(20, 20, 80, 80, box)).toEqual([20, 20, 80, 80]);
    expect(clipToBox(0, 60, 50, 60, box)).toEqual([10, 60, 50, 60]);
    expect(clipToBox(0, 0, 120, 120, box)).toEqual([10, 10, 110, 110]);
    expect(clipToBox(0, 40, 9, 52, box)).toBeUndefined();
    expect(clipToBox(5, 5, 5, 5, box)).toBeUndefined();
    expect(clipToBox(50, 50, 50, 50, box)).toEqual([50, 50, 50, 50]);
  });

  it('picks the nearest line within reach by its part inside the frame, and nothing outside the frame', () => {
    const lines: [number, number, number, number][] = [
      [0, 40, 9, 52], // wholly left of the frame, 3.6 px from (12, 50)
      [5, 15, 60, 120], // enters the frame at (10, 24.5): its part 6.7 px from (11, 12) lies outside
      [20, 20, 80, 80],
      [0, 60, 50, 60], // crosses the left edge
      [60, 90, 100, 90],
      [60, 92, 100, 92],
    ];
    const pick = (x: number, y: number) => pickLine(x, y, lines.length, (i) => lines[i]!, box, 8);
    expect(pick(12, 50)).toBe(-1);
    expect(pick(11, 12)).toBe(-1);
    expect(pick(52, 48)).toBe(2);
    expect(pick(12, 62)).toBe(3);
    // Outside the frame, nothing; the line at x 12 is 4 px away.
    expect(pick(8, 60)).toBe(-1);
    // Of two lines at the same distance, the last.
    expect(pick(80, 91)).toBe(5);
    expect(pick(80, 30)).toBe(-1);
  });
});

describe('placeBox', () => {
  const bounds = { width: 400, height: 300 };
  const size = { width: 160, height: 100 };

  it('puts the box right of and below the anchor', () => {
    expect(placeBox({ x: 50, y: 40 }, size, bounds)).toEqual({ left: 62, top: 52 });
  });

  it('flips to the left and above near the far edges', () => {
    expect(placeBox({ x: 350, y: 280 }, size, bounds)).toEqual({ left: 178, top: 168 });
  });

  it('stays inside bounds that have room on neither side', () => {
    expect(placeBox({ x: 100, y: 50 }, { width: 180, height: 100 }, { width: 220, height: 120 })).toEqual({ left: 10, top: 0 });
  });
});

describe('placeBeside', () => {
  /** Whether a box of `size` at `at` covers a point, or (with the margin of half the gap) any of a line. */
  const covers = (at: { left: number; top: number }, size: { width: number; height: number }, x: number, y: number) =>
    x >= at.left && x <= at.left + size.width && y >= at.top && y <= at.top + size.height;
  const touches = (at: { left: number; top: number }, size: { width: number; height: number }, line: readonly [number, number, number, number]) =>
    clipToBox(...line, { left: at.left - 6, top: at.top - 6, width: size.width + 12, height: size.height + 12 }) !== undefined;

  it('puts the popup off an end of the line, leaving the line uncovered (W4b screen review L10)', () => {
    // The reverse HSP of the review's capture: its popup hid the upper half of the line.
    const line = [481, 717, 767, 432] as const;
    const size = { width: 228, height: 245 };
    const at = placeBeside(line, size, { width: 924, height: 800 })!;
    expect(covers(at, size, 624, 574.5)).toBe(false);
    expect(touches(at, size, line)).toBe(false);
    // Beside one of the ends: the box's nearest corner is the gap away from it.
    const near = (x: number, y: number) =>
      Math.hypot(Math.max(at.left - x, 0, x - at.left - size.width), Math.max(at.top - y, 0, y - at.top - size.height)) <= 12 * Math.SQRT2 + 0.5;
    expect(near(481, 717) || near(767, 432)).toBe(true);
  });

  it('takes the first place that covers nothing of the line and needs no move', () => {
    // A forward HSP from the top left: the boxes beside its first end cover it or leave the bounds.
    expect(placeBeside([60, 60, 300, 300], { width: 160, height: 100 }, { width: 600, height: 500 })).toEqual({ left: 312, top: 312 });
    // Near the far corner, the box goes back over the line's side that is free.
    const at = placeBeside([300, 300, 590, 490], { width: 160, height: 100 }, { width: 600, height: 500 })!;
    expect(touches(at, { width: 160, height: 100 }, [300, 300, 590, 490])).toBe(false);
  });

  it('places a box beside a point (a line out of view, or of one pixel)', () => {
    expect(placeBeside([100, 100, 100, 100], { width: 50, height: 40 }, { width: 400, height: 300 })).toEqual({ left: 112, top: 112 });
  });

  it('gives no place where the box would cover the midpoint wherever it goes, or does not fit', () => {
    expect(placeBeside([10, 100, 290, 100], { width: 280, height: 190 }, { width: 300, height: 200 })).toBeUndefined();
    expect(placeBeside([10, 10, 20, 20], { width: 320, height: 100 }, { width: 300, height: 200 })).toBeUndefined();
    expect(placeBeside([10, 10, 20, 20], { width: 100, height: 210 }, { width: 300, height: 200 })).toBeUndefined();
  });
});
