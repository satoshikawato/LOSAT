// The one coordinate module (domain/coordinates.ts, instruction item 1, design §11.3): HSP record
// coordinates (outfmt 6's: 1-based, both ends included, start > end on the minus strand, in the
// record's own letters) to record intervals, strands and units; flanks, ends and spanning.
import { describe, expect, it } from 'vitest';
import {
  clipToRecord,
  commonStrand,
  flanked,
  hspSpan,
  interval,
  intervalLength,
  intervalText,
  isClipped,
  NO_FLANKS,
  spanning,
  spanOn,
  strandOf,
  translates,
  unitOf,
  withinRecord,
  type HspCoordinates,
} from '../../src/domain/coordinates';
import { residueUnit } from '../../src/domain/programs';

const hsp = (q: [number, number], s: [number, number], frames: [number | null, number | null] = [null, null]): HspCoordinates => ({
  q_start: q[0],
  q_end: q[1],
  s_start: s[0],
  s_end: s[1],
  query_frame: frames[0],
  subject_frame: frames[1],
});

describe('record intervals and strands of HSP coordinates', () => {
  it('orders the two coordinates and keeps the strand they give', () => {
    expect(interval(10, 50)).toEqual({ from: 10, to: 50 });
    expect(interval(50, 10)).toEqual({ from: 10, to: 50 });
    expect(hspSpan(10, 50, null, 'nucleotide')).toEqual({ from: 10, to: 50, strand: 'plus' });
    expect(hspSpan(50, 10, null, 'nucleotide')).toEqual({ from: 10, to: 50, strand: 'minus' });
    expect(intervalLength({ from: 10, to: 50 })).toBe(41);
    expect(intervalText({ from: 10, to: 50 })).toBe('10-50');
  });

  it('BLASTN: the query on its plus strand, a minus-strand subject, and a one-letter HSP whose strand the record does not give', () => {
    const minus = hsp([1, 100], [500, 401]);
    expect(spanOn(minus, 'query', 'nucleotide')).toEqual({ from: 1, to: 100, strand: 'plus' });
    expect(spanOn(minus, 'subject', 'nucleotide')).toEqual({ from: 401, to: 500, strand: 'minus' });
    const one = hsp([7, 7], [30, 30]);
    expect(spanOn(one, 'query', 'nucleotide')).toEqual({ from: 7, to: 7, strand: 'unknown' });
    expect(spanOn(one, 'subject', 'nucleotide')).toEqual({ from: 30, to: 30, strand: 'unknown' });
  });

  it('frames change no coordinate: TBLASTN, TBLASTX and BLASTX positions are the nucleotides of the record', () => {
    // TBLASTN: protein query 5-40, subject frame -2 on nucleotides 908 down to 801 (36 codons).
    const tblastn = hsp([5, 40], [908, 801], [null, -2]);
    expect(spanOn(tblastn, 'query', 'protein')).toEqual({ from: 5, to: 40, strand: 'plus' });
    expect(spanOn(tblastn, 'subject', 'nucleotide')).toEqual({ from: 801, to: 908, strand: 'minus' });
    // BLASTX (synthetic: it runs in the app from SX): the query is translated, frame -3.
    const blastx = hsp([300, 181], [11, 50], [-3, null]);
    expect(spanOn(blastx, 'query', 'nucleotide')).toEqual({ from: 181, to: 300, strand: 'minus' });
    expect(spanOn(blastx, 'subject', 'protein')).toEqual({ from: 11, to: 50, strand: 'plus' });
    const blastxPlus = hsp([4, 63], [1, 20], [1, null]);
    expect(spanOn(blastxPlus, 'query', 'nucleotide')).toEqual({ from: 4, to: 63, strand: 'plus' });
    // TBLASTX: both translated; the coordinates say the strand, the frame only where start = end.
    const tblastx = hsp([12, 101], [400, 311], [3, -1]);
    expect(spanOn(tblastx, 'query', 'nucleotide').strand).toBe('plus');
    expect(spanOn(tblastx, 'subject', 'nucleotide').strand).toBe('minus');
    expect(strandOf(9, 9, -1, 'nucleotide')).toBe('minus');
    expect(strandOf(9, 9, 2, 'nucleotide')).toBe('plus');
    expect(strandOf(9, 9, 0, 'nucleotide')).toBe('unknown');
    expect(strandOf(9, 9, undefined, 'protein')).toBe('plus');
  });

  it('knows which sequences a program translates and their units', () => {
    expect(translates('blastn', 'nucleotide')).toBe(false);
    expect(translates('blastp', 'protein')).toBe(false);
    expect(translates('tblastn', 'nucleotide')).toBe(true);
    expect(translates('tblastn', 'protein')).toBe(false);
    expect(translates('tblastx', 'nucleotide')).toBe(true);
    expect(translates('blastx', 'nucleotide')).toBe(true);
    expect(unitOf('nucleotide')).toBe('nt');
    expect(unitOf('protein')).toBe('aa');
    expect(residueUnit('nucleotide')).toBe(unitOf('nucleotide'));
    expect(residueUnit('protein')).toBe(unitOf('protein'));
  });
});

describe('flanks, the ends of a record and spanning', () => {
  it('adds left and right flanks on the record axis whatever the strand', () => {
    const minus = hspSpan(500, 401, null, 'nucleotide');
    expect(flanked(minus, { left: 10, right: 30 }, 1000)).toEqual({
      requested: { from: 391, to: 530 },
      actual: { from: 391, to: 530 },
      clippedLeft: false,
      clippedRight: false,
    });
    expect(flanked(minus, NO_FLANKS, 1000).actual).toEqual({ from: 401, to: 500 });
  });

  it('cuts at both ends and keeps the requested interval', () => {
    const clipped = flanked({ from: 5, to: 95 }, { left: 20, right: 20 }, 100);
    expect(clipped).toEqual({ requested: { from: -15, to: 115 }, actual: { from: 1, to: 100 }, clippedLeft: true, clippedRight: true });
    expect(isClipped(clipped)).toBe(true);
    expect(isClipped(clipToRecord({ from: 1, to: 100 }, 100))).toBe(false);
    expect(clipToRecord({ from: 0, to: 3 }, 100)).toMatchObject({ actual: { from: 1, to: 3 }, clippedLeft: true, clippedRight: false });
  });

  it('refuses intervals outside the record, negative or fractional flanks', () => {
    expect(() => clipToRecord({ from: 101, to: 120 }, 100)).toThrow(/outside the record/);
    expect(() => clipToRecord({ from: -5, to: 0 }, 100)).toThrow(/outside the record/);
    expect(() => clipToRecord({ from: 5, to: 4 }, 100)).toThrow(/not an interval/);
    expect(() => clipToRecord({ from: 1, to: 1 }, 0)).toThrow(RangeError);
    expect(() => flanked({ from: 5, to: 9 }, { left: -1, right: 0 }, 100)).toThrow(/left flank/);
    expect(() => flanked({ from: 5, to: 9 }, { left: 0, right: 1.5 }, 100)).toThrow(/right flank/);
  });

  it('spans several intervals without joining them otherwise', () => {
    expect(spanning([{ from: 40, to: 60 }, { from: 10, to: 20 }, { from: 50, to: 55 }])).toEqual({ from: 10, to: 60 });
    expect(() => spanning([])).toThrow(RangeError);
  });

  it('gives the strand that spans share, or mixed', () => {
    const plus = hspSpan(1, 5, null, 'nucleotide');
    const minus = hspSpan(9, 6, null, 'nucleotide');
    const one = hspSpan(3, 3, null, 'nucleotide');
    expect(commonStrand([plus, plus])).toBe('plus');
    expect(commonStrand([minus])).toBe('minus');
    expect(commonStrand([one])).toBe('unknown');
    expect(commonStrand([plus, minus])).toBe('mixed');
    expect(commonStrand([plus, one])).toBe('mixed');
  });

  it('checks that an interval lies in its record', () => {
    expect(withinRecord({ from: 1, to: 100 }, 100)).toBe(true);
    expect(withinRecord({ from: 0, to: 10 }, 100)).toBe(false);
    expect(withinRecord({ from: 90, to: 101 }, 100)).toBe(false);
    expect(withinRecord({ from: 5, to: 4 }, 100)).toBe(false);
    expect(withinRecord({ from: 1.5, to: 4 }, 100)).toBe(false);
  });
});
