// Extraction written as it is read (domain/extraction.ts `extractionSteps`, S15 item 1, design
// §12.1, code review L1): the steps of a plan read at most a bounded number of residues at once,
// in the order of the file, and a sequence written in parts has the bytes of one written whole.
import { describe, expect, it } from 'vitest';
import {
  EXTRACTION_READ_RESIDUES,
  extractionSteps,
  extractionTarget,
  FASTA_LINE_WIDTH,
  fastaLines,
  planExtraction,
  readRequests,
  sequenceFasta,
  sequenceHeaderLine,
  sequencePartLines,
  type ExtractedHsp,
  type ExtractionPlan,
  type RunRef,
} from '../../src/domain/extraction';

const RUN: RunRef = { runId: 'run-a', number: 2, program: 'blastn' };
const decoder = new TextDecoder();
const text = (bytes: Uint8Array) => decoder.decode(bytes);

const hsp = (s_idx: number, rank: number, s_start: number, s_end: number): ExtractedHsp => ({
  q_idx: 0,
  s_idx,
  rank,
  q_start: 1,
  q_end: 10,
  s_start,
  s_end,
  query_frame: null,
  subject_frame: null,
});

/** Records of 10,000, 500 and 2,000 letters; HSPs given as [record, start, end]. */
const LENGTHS = [10_000, 500, 2_000];
function plan(hsps: ReadonlyArray<readonly [number, number, number]>, region: 'hit' | 'whole' = 'hit'): ExtractionPlan {
  const targets = hsps.map(([s, from, to], rank) => extractionTarget(RUN, hsp(s, rank, from, to), 'subject', { id: `r${s}`, length: LENGTHS[s]! }));
  return planExtraction(targets, { region: { kind: region }, join: 'separate' });
}

/** Letters of a record, mixed case, the same for the same record. */
function letters(length: number, seed: number): Uint8Array {
  return Uint8Array.from({ length }, (_, i) => 'ACGTacgtNu'.charCodeAt((i * 31 + seed * 7 + ((i * i) >> 4)) % 10));
}

describe('extraction steps', () => {
  it('reads pieces that fit in one read together, one request per record, as readRequests does', () => {
    const small = plan([
      [0, 100, 149],
      [1, 5, 1],
      [0, 700, 709],
      [2, 1, 60],
    ]);
    expect(extractionSteps(small)).toEqual([{ kind: 'pieces', first: 0, end: 4, requests: readRequests(small) }]);
    expect(EXTRACTION_READ_RESIDUES % FASTA_LINE_WIDTH).toBe(0);
    expect(EXTRACTION_READ_RESIDUES).toBeGreaterThan(4_000_000);
  });

  it('starts a new step where the residues or the pieces of a step would pass their bound, in tray order', () => {
    const p = plan([
      [0, 1, 600],
      [1, 1, 300],
      [0, 601, 900],
      [2, 1, 60],
      [2, 61, 120],
      [2, 121, 180],
    ]);
    const steps = extractionSteps(p, 960, 2);
    expect(steps.map((step) => (step.kind === 'pieces' ? [step.first, step.end, step.requests.map((r) => [r.position, r.intervals.length])] : step))).toEqual([
      [0, 2, [[0, 1], [1, 1]]],
      [2, 4, [[0, 1], [2, 1]]],
      [4, 6, [[2, 2]]],
    ]);
    // By residues alone: 600 + 300 fit in 960, the next 300 does not.
    expect(extractionSteps(p, 960, 1000).map((step) => step.kind === 'pieces' && [step.first, step.end])).toEqual([
      [0, 2],
      [2, 6],
    ]);
  });

  it('reads a piece longer than one read alone, in parts of whole lines, and marks a whole record', () => {
    const p = plan(
      [
        [1, 10, 20],
        [0, 1, 10_000],
        [2, 300, 1],
      ],
      'whole',
    );
    const steps = extractionSteps(p, 960);
    expect(steps.map((step) => step.kind)).toEqual(['pieces', 'parts', 'parts']);
    const whole = steps[1]!;
    expect(whole.kind === 'parts' && [whole.piece, whole.wholeRecord, whole.parts.length]).toEqual([1, true, 11]);
    if (whole.kind !== 'parts') throw new Error('parts expected');
    expect(whole.parts[0]).toEqual({ from: 1, to: 960 });
    expect(whole.parts[10]).toEqual({ from: 9601, to: 10_000 });
    for (const [i, part] of whole.parts.entries()) {
      expect(part.to - part.from + 1).toBeLessThanOrEqual(960);
      if (i > 0) expect(part.from).toBe(whole.parts[i - 1]!.to + 1);
    }
    // A long interval that is not the whole record is read in parts too, but not checked as a whole.
    const partial = extractionSteps(plan([[0, 20, 2000]]), 960);
    expect(partial).toEqual([
      {
        kind: 'parts',
        piece: 0,
        parts: [
          { from: 20, to: 979 },
          { from: 980, to: 1939 },
          { from: 1940, to: 2000 },
        ],
        wholeRecord: false,
      },
    ]);
  });

  it('refuses a read that is not a whole number of lines', () => {
    const p = plan([[0, 1, 10]]);
    for (const bad of [1000, 59, 0, -60, 60.5]) expect(() => extractionSteps(p, bad)).toThrow(RangeError);
    expect(() => extractionSteps(p, 960, 0)).toThrow(RangeError);
    expect(extractionSteps(p, 60)).toHaveLength(1);
  });

  it('writes a sequence in parts with the same bytes as written whole', () => {
    for (const [from, to, read] of [
      [1, 10_000, 960],
      [1, 10_000, 60],
      [7, 9_999, 120],
      [1, 961, 960],
      [1, 960, 60],
    ] as const) {
      const p = plan([[0, from, to]]);
      const piece = p.pieces[0]!;
      const record = letters(10_000, 3);
      const expected = text(sequenceFasta(piece, record.subarray(from - 1, to)));
      const steps = extractionSteps(p, read);
      let got = '';
      for (const step of steps) {
        if (step.kind === 'pieces') {
          got += text(sequenceFasta(piece, record.subarray(from - 1, to)));
        } else {
          got += text(sequenceHeaderLine(piece));
          for (const part of step.parts) got += text(sequencePartLines(piece, part, record.subarray(part.from - 1, part.to)));
        }
      }
      expect(got).toBe(expected);
    }
    const piece = plan([[0, 1, 200]]).pieces[0]!;
    expect(() => sequencePartLines(piece, { from: 1, to: 120 }, letters(119, 1))).toThrow(/119 residues were given for 1-120 \(120 nt\)/);
    expect(text(fastaLines(new TextEncoder().encode('ABCDE'), 2))).toBe('AB\nCD\nE\n');
    expect(fastaLines(new Uint8Array())).toEqual(new Uint8Array());
  });
});
