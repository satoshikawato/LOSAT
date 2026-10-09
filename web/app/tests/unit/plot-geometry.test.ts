// The dot plot's geometry (src/domain/plot-geometry.ts): the same scale on both axes, the mapping
// between positions and pixels with the origin at the top left, zoom and pan within the records,
// and the placing of the HSP popup.
import { describe, expect, it } from 'vitest';
import {
  clampView,
  frameRange,
  fromPixelX,
  fromPixelY,
  fullView,
  MIN_SIDE,
  panView,
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
