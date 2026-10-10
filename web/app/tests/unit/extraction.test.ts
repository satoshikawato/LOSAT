// Extraction planning and its two FASTA formats (domain/extraction.ts, instruction items 2-4,
// design §11.4): regions, joins, record ends, strands, the sequence header and the gapped
// alignment records written from the HSP record's aligned rows.
import { describe, expect, it } from 'vitest';
import {
  alignmentFasta,
  alignmentHeaders,
  extractionTarget,
  fastaRecord,
  hspLabel,
  missingAlignmentNote,
  planExtraction,
  readRequests,
  recordLabel,
  sequenceFasta,
  sequenceHeader,
  type AlignedHsp,
  type ExtractionTarget,
  type RunRef,
} from '../../src/domain/extraction';

const decoder = new TextDecoder();
const text = (bytes: Uint8Array | undefined) => decoder.decode(bytes);

const BLASTN: RunRef = { runId: 'run-a', number: 3, program: 'blastn' };
const TBLASTN: RunRef = { runId: 'run-b', number: 4, program: 'tblastn' };
const TBLASTX: RunRef = { runId: 'run-c', number: 5, program: 'tblastx' };
const BLASTX: RunRef = { runId: 'run-x', number: 6, program: 'blastx' };

function record(fields: Partial<AlignedHsp> & Pick<AlignedHsp, 'q_start' | 'q_end' | 's_start' | 's_end'>): AlignedHsp {
  return {
    q_idx: 0,
    s_idx: 0,
    rank: 0,
    query_frame: null,
    subject_frame: null,
    query_aligned: null,
    subject_aligned: null,
    ...fields,
  };
}

const SUBJECT = { id: 'chr1', length: 1000 };
const target = (hsp: AlignedHsp, run: RunRef = BLASTN, role: 'query' | 'subject' = 'subject', rec = SUBJECT): ExtractionTarget =>
  extractionTarget(run, hsp, role, rec);

describe('extraction targets', () => {
  it('takes the record, unit, span and label of an HSP on either role', () => {
    const hsp = record({ q_idx: 2, s_idx: 5, rank: 1, q_start: 1, q_end: 100, s_start: 600, s_end: 501 });
    expect(target(hsp)).toEqual({
      runId: 'run-a',
      runNumber: 3,
      role: 'subject',
      position: 5,
      recordId: 'chr1',
      recordLength: 1000,
      unit: 'nt',
      span: { from: 501, to: 600, strand: 'minus' },
      hsp: '3.2',
    });
    expect(target(hsp, BLASTN, 'query', { id: 'q', length: 120 })).toMatchObject({ position: 2, span: { from: 1, to: 100, strand: 'plus' } });
    expect(hspLabel(0, 0)).toBe('1.1');
  });

  it('BLASTX (synthetic records until SX): the query positions are nucleotides whatever the frame', () => {
    const minus = record({ q_start: 330, q_end: 211, s_start: 1, s_end: 40, query_frame: -2 });
    expect(target(minus, BLASTX, 'query', { id: 'contig', length: 400 })).toMatchObject({
      unit: 'nt',
      span: { from: 211, to: 330, strand: 'minus' },
    });
    expect(target(minus, BLASTX, 'subject', { id: 'prot', length: 50 })).toMatchObject({ unit: 'aa', span: { from: 1, to: 40, strand: 'plus' } });
    const plan = planExtraction([target(minus, BLASTX, 'query', { id: 'contig', length: 400 })], { region: { kind: 'flanked', flanks: { left: 5, right: 100 } }, join: 'separate' });
    expect(sequenceHeader(plan.pieces[0]!)).toBe(
      'contig:206-400 run=6 query_record=1 length=400 unit=nt hsps=1.1 hit_strand=minus requested=206-430',
    );
  });

  it('TBLASTN: a minus frame subject in nucleotides and the protein query in residues', () => {
    const hsp = record({ q_start: 2, q_end: 30, s_start: 900, s_end: 814, subject_frame: -1 });
    expect(target(hsp, TBLASTN, 'subject', { id: 'genome', length: 5000 })).toMatchObject({ unit: 'nt', span: { from: 814, to: 900, strand: 'minus' } });
    expect(target(hsp, TBLASTN, 'query', { id: 'p', length: 30 })).toMatchObject({ unit: 'aa', span: { from: 2, to: 30, strand: 'plus' } });
  });

  it('refuses an HSP that is not within its record', () => {
    const hsp = record({ q_start: 1, q_end: 10, s_start: 990, s_end: 1010 });
    expect(() => planExtraction([target(hsp)], { region: { kind: 'hit' }, join: 'separate' })).toThrow(/not within subject record 1/);
  });
});

describe('planning', () => {
  const a = record({ rank: 0, q_start: 1, q_end: 50, s_start: 100, s_end: 149 });
  const b = record({ rank: 1, q_start: 60, q_end: 80, s_start: 130, s_end: 110 });
  const c = record({ rank: 2, q_start: 90, q_end: 99, s_start: 700, s_end: 709 });
  const other = record({ q_idx: 1, s_idx: 1, rank: 0, q_start: 5, q_end: 9, s_start: 5, s_end: 1 });
  const targets = [target(a), target(other, BLASTN, 'subject', { id: 'chr1', length: 40 }), target(b), target(c)];

  it('separate: one sequence per HSP in tray order, overlaps kept', () => {
    const plan = planExtraction(targets, { region: { kind: 'hit' }, join: 'separate' });
    expect(plan.pieces.map((p) => [p.position, p.interval.actual, p.hsps, p.strand])).toEqual([
      [0, { from: 100, to: 149 }, ['1.1'], 'plus'],
      [1, { from: 1, to: 5 }, ['2.1'], 'minus'],
      [0, { from: 110, to: 130 }, ['1.2'], 'minus'],
      [0, { from: 700, to: 709 }, ['1.3'], 'plus'],
    ]);
    expect(plan.notes).toEqual([]);
    expect(readRequests(plan)).toEqual([
      { runId: 'run-a', role: 'subject', position: 0, intervals: [{ from: 100, to: 149 }, { from: 110, to: 130 }, { from: 700, to: 709 }], pieces: [0, 2, 3] },
      { runId: 'run-a', role: 'subject', position: 1, intervals: [{ from: 1, to: 5 }], pieces: [1] },
    ]);
  });

  it('spanning: one interval from the smallest to the largest position on each record, then the flanks', () => {
    const plan = planExtraction(targets, { region: { kind: 'flanked', flanks: { left: 10, right: 400 } }, join: 'spanning' });
    expect(plan.pieces.map((p) => [p.position, p.interval, p.hsps, p.strand])).toEqual([
      [0, { requested: { from: 90, to: 1109 }, actual: { from: 90, to: 1000 }, clippedLeft: false, clippedRight: true }, ['1.1', '1.2', '1.3'], 'mixed'],
      [1, { requested: { from: -9, to: 405 }, actual: { from: 1, to: 40 }, clippedLeft: true, clippedRight: true }, ['2.1'], 'minus'],
    ]);
    expect(sequenceHeader(plan.pieces[0]!)).toBe(
      'chr1:90-1000 run=3 subject_record=1 length=1000 unit=nt hsps=1.1,1.2,1.3 hit_strand=mixed requested=90-1109',
    );
    expect(sequenceHeader(plan.pieces[1]!)).toBe('chr1:1-40 run=3 subject_record=2 length=40 unit=nt hsps=2.1 hit_strand=minus requested=-9-405');
  });

  it('whole: one sequence per record listing every HSP of the tray on it, whatever the join', () => {
    for (const join of ['separate', 'spanning'] as const) {
      const plan = planExtraction(targets, { region: { kind: 'whole' }, join });
      expect(plan.pieces.map((p) => [p.position, p.interval.actual, p.hsps])).toEqual([
        [0, { from: 1, to: 1000 }, ['1.1', '1.2', '1.3']],
        [1, { from: 1, to: 40 }, ['2.1']],
      ]);
      expect(sequenceHeader(plan.pieces[0]!)).toBe('chr1:1-1000 run=3 subject_record=1 length=1000 unit=nt hsps=1.1,1.2,1.3 hit_strand=mixed');
    }
  });

  it('never joins records of different runs or roles, even with the same ID, length and position', () => {
    const run2: RunRef = { ...BLASTN, runId: 'run-other', number: 9 };
    const plan = planExtraction([target(a), target(a, run2), target(a, BLASTN, 'query', SUBJECT)], { region: { kind: 'hit' }, join: 'spanning' });
    expect(plan.pieces.map((p) => [p.runNumber, p.role, p.interval.actual])).toEqual([
      [3, 'subject', { from: 100, to: 149 }],
      [9, 'subject', { from: 100, to: 149 }],
      [3, 'query', { from: 1, to: 50 }],
    ]);
  });

  it('says that it does not decide the strand of a one-letter BLASTN HSP, and why', () => {
    const one = record({ rank: 4, q_start: 10, q_end: 10, s_start: 20, s_end: 20 });
    const plan = planExtraction([target(one)], { region: { kind: 'hit' }, join: 'separate' });
    expect(plan.pieces[0]!.strand).toBe('unknown');
    expect(sequenceHeader(plan.pieces[0]!)).toMatch(/hit_strand=unknown$/);
    expect(plan.notes).toEqual([
      'HSP 1.5 of run 3 covers one letter of subject record 1, so its HSP record does not give its strand and the extraction does not decide it (hit_strand=unknown); the Strand= line of its outfmt 0 section shows it.',
    ]);
  });

  it('refuses HSPs that disagree about one record', () => {
    expect(() => planExtraction([target(a), target(b, BLASTN, 'subject', { id: 'chr1', length: 999 })], { region: { kind: 'hit' }, join: 'separate' })).toThrow(
      /disagree about the record/,
    );
  });
});

describe('sequence FASTA', () => {
  it('writes the header and the residues as read, 60 per line, case kept', () => {
    const plan = planExtraction([target(record({ q_start: 1, q_end: 9, s_start: 1, s_end: 130 }), BLASTN, 'subject', { id: 'x', length: 130 })], {
      region: { kind: 'hit' },
      join: 'separate',
    });
    const residues = new TextEncoder().encode('acgU'.repeat(32) + 'NN');
    const fasta = text(sequenceFasta(plan.pieces[0]!, residues));
    const lines = fasta.split('\n');
    expect(lines[0]).toBe('>x:1-130 run=3 subject_record=1 length=130 unit=nt hsps=1.1 hit_strand=plus');
    expect(lines.slice(1).map((line) => line.length)).toEqual([60, 60, 10, 0]);
    expect(lines.slice(1).join('')).toBe('acgU'.repeat(32) + 'NN');
    expect(() => sequenceFasta(plan.pieces[0]!, residues.subarray(1))).toThrow(/129 residues were given/);
  });

  it('names a record without an ID as the engine does in outfmt 6', () => {
    expect(recordLabel('', 'query', 0)).toBe('Query_1');
    expect(recordLabel('', 'subject', 6)).toBe('Subject_7');
    expect(recordLabel('seq1', 'subject', 6)).toBe('seq1');
    const plan = planExtraction([target(record({ q_start: 1, q_end: 9, s_start: 3, s_end: 9 }), BLASTN, 'subject', { id: '', length: 10 })], {
      region: { kind: 'hit' },
      join: 'separate',
    });
    expect(sequenceHeader(plan.pieces[0]!)).toBe('Subject_1:3-9 run=3 subject_record=1 length=10 unit=nt hsps=1.1 hit_strand=plus');
  });

  it('wraps any letters at the given width', () => {
    expect(text(fastaRecord('h', new TextEncoder().encode('ABCDE'), 2))).toBe('>h\nAB\nCD\nE\n');
    expect(text(fastaRecord('h é', new TextEncoder().encode('AB'), 2))).toBe('>h é\nAB\n');
  });
});

describe('gapped alignment FASTA', () => {
  it('writes both rows as the engine wrote them, with the search direction of the coordinates', () => {
    const hsp = record({
      q_idx: 1,
      s_idx: 2,
      rank: 0,
      q_start: 1,
      q_end: 70,
      s_start: 500,
      s_end: 430,
      query_aligned: 'ACGT-acgt'.repeat(8),
      subject_aligned: 'ACGTTacgt'.repeat(8),
    });
    const fasta = text(alignmentFasta(BLASTN, hsp, { query: 'q2', subject: 'chr9' }));
    const lines = fasta.split('\n');
    expect(lines[0]).toBe('>q2:1-70 run=3 query_record=2 hsp=2.1 aligned');
    expect(lines[1]).toBe(('ACGT-acgt'.repeat(8)).slice(0, 60));
    expect(lines[2]).toBe(('ACGT-acgt'.repeat(8)).slice(60));
    expect(lines[3]).toBe('>chr9:500-430 run=3 subject_record=3 hsp=2.1 aligned');
    expect(lines.slice(4).join('')).toBe('ACGTTacgt'.repeat(8));
  });

  it('shows frames for translated sequences only', () => {
    const hsp = record({ q_start: 3, q_end: 92, s_start: 400, s_end: 311, query_frame: 3, subject_frame: -1, query_aligned: 'MK', subject_aligned: 'MR' });
    expect(alignmentHeaders(TBLASTX, hsp, { query: 'a', subject: '' })).toEqual({
      query: 'a:3-92 run=5 query_record=1 hsp=1.1 aligned frame=3',
      subject: 'Subject_1:400-311 run=5 subject_record=1 hsp=1.1 aligned frame=-1',
    });
    const tblastn = record({ q_start: 1, q_end: 30, s_start: 90, s_end: 1, query_frame: 1, subject_frame: -2 });
    expect(alignmentHeaders(TBLASTN, tblastn, { query: 'p', subject: 'g' })).toEqual({
      query: 'p:1-30 run=4 query_record=1 hsp=1.1 aligned',
      subject: 'g:90-1 run=4 subject_record=1 hsp=1.1 aligned frame=-2',
    });
    const blastp = record({ q_start: 1, q_end: 30, s_start: 1, s_end: 30, query_frame: 1, subject_frame: 1 });
    expect(alignmentHeaders({ runId: 'p', number: 1, program: 'blastp' }, blastp, { query: 'p', subject: 's' }).subject).toBe('s:1-30 run=1 subject_record=1 hsp=1.1 aligned');
    const blastx = record({ q_start: 300, q_end: 211, s_start: 1, s_end: 30, query_frame: -1 });
    expect(alignmentHeaders(BLASTX, blastx, { query: 'n', subject: 'p' }).query).toBe('n:300-211 run=6 query_record=1 hsp=1.1 aligned frame=-1');
  });

  it('reports an HSP without aligned rows instead of inventing them', () => {
    const hsp = record({ rank: 2, q_start: 1, q_end: 9, s_start: 1, s_end: 9 });
    expect(alignmentFasta(BLASTN, hsp, { query: 'q', subject: 's' })).toBeUndefined();
    expect(missingAlignmentNote(BLASTN, hsp)).toBe('HSP 1.3 of run 3 has no aligned sequences in its HSP record, so it has no gapped alignment to write.');
  });
});
