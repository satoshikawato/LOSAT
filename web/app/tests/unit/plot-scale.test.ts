// The plots' scales and display classes (src/domain/plot-scale.ts): blast2dotplot.py's tick table
// and labels on the visible span, the Graphic Summary's ruler, the identity opacity classes and
// NCBI's "Alignment Scores" bins.
import { describe, expect, it } from 'vitest';
import {
  axisTicks,
  axisUnit,
  identityClass,
  IDENTITY_CLASSES,
  MAX_MINOR_TICKS,
  rulerTicks,
  scoreBin,
  SCORE_BINS,
  tickLabel,
  tickSteps,
} from '../../src/domain/plot-scale';

describe('tickSteps', () => {
  it('is the tick_size table of blast2dotplot.py, its bounds inclusive', () => {
    const table: [number, number, number][] = [
      [1, 100, 10],
      [1000, 100, 10],
      [1001, 500, 100],
      [5000, 500, 100],
      [5001, 1000, 500],
      [10000, 1000, 500],
      [10001, 2000, 1000],
      [20000, 2000, 1000],
      [20001, 5000, 1000],
      [50000, 5000, 1000],
      [50001, 10000, 1000],
      [100000, 10000, 1000],
      [100001, 50000, 10000],
      [1000000, 50000, 10000],
      [1000001, 1000000, 200000],
      [5000000, 1000000, 200000],
      [5000001, 1000000, 500000],
      [3e9, 1000000, 500000],
    ];
    for (const [span, major, minor] of table) expect(tickSteps(span), `span ${span}`).toEqual({ major, minor });
  });
});

describe('axisUnit and tickLabel', () => {
  it('chooses bp, kbp or Mbp by the span (aa, kaa or Maa on a protein axis)', () => {
    expect(axisUnit(4999, 'nt')).toEqual({ name: 'bp', divisor: 1 });
    expect(axisUnit(5000, 'nt')).toEqual({ name: 'kbp', divisor: 1000 });
    expect(axisUnit(999999, 'nt')).toEqual({ name: 'kbp', divisor: 1000 });
    expect(axisUnit(1000000, 'nt')).toEqual({ name: 'Mbp', divisor: 1000000 });
    expect(axisUnit(300, 'aa').name).toBe('aa');
    expect(axisUnit(6000, 'aa').name).toBe('kaa');
    expect(axisUnit(2e6, 'aa').name).toBe('Maa');
  });

  it('writes a decimal only where the value is not whole in the unit', () => {
    const kbp = axisUnit(5000, 'nt');
    expect(tickLabel(3000, kbp)).toBe('3');
    expect(tickLabel(2500, kbp)).toBe('2.5');
    expect(tickLabel(0, kbp)).toBe('0');
    const mbp = axisUnit(1e6, 'nt');
    expect(tickLabel(1050000, mbp)).toBe('1.05');
    expect(tickLabel(12000000, mbp)).toBe('12');
    // A bp label far along a long record keeps its thousands separators.
    expect(tickLabel(2345500, axisUnit(3000, 'nt'))).toBe('2,345,500');
  });
});

describe('axisTicks', () => {
  it('places the major and minor ticks of the visible span, minor ticks off the major ones', () => {
    const ticks = axisTicks(0, 400, 'nt');
    expect(ticks.major).toEqual([0, 100, 200, 300, 400]);
    expect(ticks.minor).toHaveLength(36);
    expect(ticks.minor).not.toContain(100);
    expect(ticks.unit.name).toBe('bp');
  });

  it('chooses the steps again for a zoomed view and starts at the first multiple in view', () => {
    const ticks = axisTicks(123456.7, 131000, 'nt');
    expect(ticks.steps).toEqual({ major: 1000, minor: 500 });
    expect(ticks.major[0]).toBe(124000);
    expect(ticks.major.at(-1)).toBe(131000);
    expect(ticks.minor[0]).toBe(123500);
    expect(ticks.unit.name).toBe('kbp');
  });

  it('leaves out the minor ticks when more than 50 would be drawn', () => {
    expect(axisTicks(0, 500, 'nt').minor).toHaveLength(45);
    expect(axisTicks(0, 600, 'nt').minor).toEqual([]);
    expect(axisTicks(0, 600, 'nt').major).toHaveLength(7);
    expect(MAX_MINOR_TICKS).toBe(50);
  });

  it('has no ticks for an empty span', () => {
    expect(axisTicks(10, 10, 'nt').major).toEqual([]);
  });
});

describe('rulerTicks', () => {
  it('labels 1, the length and round positions between', () => {
    expect(rulerTicks(628, 8)).toEqual([1, 100, 200, 300, 400, 500, 628]);
    expect(rulerTicks(10000, 6)).toEqual([1, 2000, 4000, 6000, 8000, 10000]);
  });

  it('leaves out round positions too close to an end, and copes with tiny records', () => {
    expect(rulerTicks(1020, 11)).toEqual([1, 200, 400, 600, 800, 1020]);
    expect(rulerTicks(3, 10)).toEqual([1, 2, 3]);
    expect(rulerTicks(1, 10)).toEqual([1]);
  });
});

describe('identityClass', () => {
  it('truncates the pident string and classes it as blast2dotplot.py does', () => {
    const opacity = (pident: string) => IDENTITY_CLASSES[identityClass(pident)]!.opacity;
    expect(opacity('60.000')).toBe(0.4);
    expect(opacity('60.999')).toBe(0.4);
    expect(opacity('61.000')).toBe(0.6);
    expect(opacity('70.990')).toBe(0.6);
    expect(opacity('71.0')).toBe(0.8);
    expect(opacity('80.952')).toBe(0.8);
    expect(opacity('81.000')).toBe(1);
    expect(opacity('100.000')).toBe(1);
    expect(opacity('12.5')).toBe(0.4);
  });

  it('draws a pident that is not a number (the FakeEngine writes FAKE) opaque', () => {
    expect(identityClass('FAKE')).toBe(3);
    expect(identityClass('')).toBe(3);
  });
});

describe('scoreBin', () => {
  it('uses NCBI\'s half-open bins and colours', () => {
    expect([0, 39.99, 40, 49.9, 50, 79.9, 80, 199.9, 200, 1e6].map(scoreBin)).toEqual([0, 0, 1, 1, 2, 2, 3, 3, 4, 4]);
    expect(SCORE_BINS.map((bin) => bin.color)).toEqual(['#000000', '#0020e9', '#75ea4c', '#db3de9', '#db3324']);
    expect(SCORE_BINS.map((bin) => bin.label)).toEqual(['< 40', '40 - 50', '50 - 80', '80 - 200', '>= 200']);
    expect(scoreBin(Number.NaN)).toBe(0);
  });
});
