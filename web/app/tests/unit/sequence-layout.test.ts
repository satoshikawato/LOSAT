// Reading a record's residues from its source (domain/sequence-layout.ts, abi_v2.md §9): records
// that the test writes itself, with the layout it computes from the bytes it wrote
// (support/fasta-writer.ts), read back for random intervals and every chunk size.
import { describe, expect, it } from 'vitest';
import type { IndexedRecord } from '../../src/domain/dataset';
import { ForwardReader, readerKind, readStart, residueCounts, stores, uniformOffset } from '../../src/domain/sequence-layout';
import { between, counts, FastaWriter, indexedRecord, residues, seeded, type Lines, type ReaderKind, type WrittenRecord } from './support/fasta-writer';

const NUCLEOTIDES = 'ACGTUNRYKMSWBDHV';
const AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWYUOBZX*';

/** Reads residues `from` to `to` (1-based) of a record from `bytes`, fed `chunk` bytes at a time. */
function read(bytes: Uint8Array, record: IndexedRecord, kind: ReaderKind, from: number, to: number, chunk: number): string {
  const start = readStart(record.line_layout, record.sequence_offset, from - 1);
  const reader = new ForwardReader(kind, start.skip, to - from + 1);
  for (let at = start.offset; at < record.end_offset && reader.feed(bytes.subarray(at, Math.min(record.end_offset, at + chunk))); at += chunk) {
    // feed() returns false once the interval is read.
  }
  return latin1.decode(reader.residues());
}

const latin1 = new TextDecoder('latin1');

interface Case {
  readonly name: string;
  readonly kind: ReaderKind;
  readonly lines: Lines;
  readonly length: number;
  readonly every?: number;
}

const CASES: readonly Case[] = [
  { name: 'kind 1, uniform, LF', kind: 1, lines: { kind: 'uniform', width: 7, eol: '\n' }, length: 40 },
  { name: 'kind 1, uniform, CR LF', kind: 1, lines: { kind: 'uniform', width: 5, eol: '\r\n' }, length: 33 },
  { name: 'kind 1, uniform, eol 3 from trailing spaces', kind: 1, lines: { kind: 'uniform', width: 6, eol: '  \n' }, length: 31 },
  { name: 'kind 2, uniform, eol 4 from a tab and CR LF', kind: 2, lines: { kind: 'uniform', width: 4, eol: ' \t\r\n' }, length: 18 },
  { name: 'kind 1, one line', kind: 1, lines: { kind: 'uniform', width: 100, eol: '\n' }, length: 23 },
  { name: 'kind 1, ragged LF with noise', kind: 1, lines: { kind: 'ragged', minWidth: 1, maxWidth: 9, eol: '\n', noise: true }, length: 45, every: 7 },
  { name: 'kind 1, ragged CR LF with noise', kind: 1, lines: { kind: 'ragged', minWidth: 2, maxWidth: 8, eol: '\r\n', noise: true }, length: 40, every: 6 },
  { name: 'kind 2, ragged LF with noise', kind: 2, lines: { kind: 'ragged', minWidth: 1, maxWidth: 11, eol: '\n', noise: true }, length: 38, every: 5 },
  { name: 'kind 2, ragged CR LF, comments only', kind: 2, lines: { kind: 'ragged', minWidth: 3, maxWidth: 6, eol: '\r\n', comments: true }, length: 30, every: 4 },
];

describe('reading residues forward from the layout of a record', () => {
  for (const c of CASES) {
    it(`${c.name}: every interval start, random ends, every chunk size`, () => {
      const random = seeded(c.length * 31 + c.kind);
      const writer = new FastaWriter(random, c.kind).raw('; a comment before the first record\n\n');
      const alphabet = c.kind === 1 ? NUCLEOTIDES : AMINO_ACIDS;
      writer.record('first one', residues(random, 3, alphabet, 0.5), { kind: 'uniform', width: 60, eol: '\n' });
      const target = writer.record('target with a description', residues(random, c.length, alphabet, 0.4), c.lines);
      writer.record('last', residues(random, 4, alphabet), { kind: 'uniform', width: 60, eol: '\n' });
      const bytes = writer.bytes();
      const record = indexedRecord(target, 1, c.kind, c.every);
      expect(record.line_layout.kind).toBe(c.lines.kind === 'uniform' ? 'uniform' : 'checkpoints');
      for (let from = 1; from <= c.length; from++) {
        const to = between(random, from, c.length);
        const expected = target.letters.slice(from - 1, to);
        const span = record.end_offset - readStart(record.line_layout, record.sequence_offset, from - 1).offset;
        for (let chunk = 1; chunk <= span; chunk++) {
          expect(read(bytes, record, c.kind, from, to, chunk), `${from}-${to} in chunks of ${chunk}`).toBe(expected);
        }
      }
      expect(read(bytes, record, c.kind, 1, c.length, 3)).toBe(target.letters);
    });
  }

  it('reads records of more than 65536 residues through their checkpoints, and a uniform one with CR LF', () => {
    const random = seeded(7);
    const writer = new FastaWriter(random, 1);
    const ragged = writer.record('ragged', residues(random, 140_000, 'ACGTUN', 0.3), { kind: 'ragged', minWidth: 50, maxWidth: 90, eol: '\n', noise: true });
    const uniform = writer.record('uniform', residues(random, 70_001, 'acgtACGT'), { kind: 'uniform', width: 60, eol: '\r\n' });
    const bytes = writer.bytes();
    const raggedRecord = indexedRecord(ragged, 0, 1);
    const uniformRecord = indexedRecord(uniform, 1, 1);
    expect(raggedRecord.line_layout).toMatchObject({ kind: 'checkpoints', every: 65_536 });
    expect(raggedRecord.line_layout.kind === 'checkpoints' && raggedRecord.line_layout.offsets).toHaveLength(3);
    for (const [from, to] of [
      [1, 140_000],
      [65_530, 65_545],
      [65_536, 65_537],
      [65_537, 131_073],
      [131_072, 131_074],
      [139_990, 140_000],
    ] as const) {
      for (const chunk of [1, 7, 4096, 1 << 20]) {
        expect(read(bytes, raggedRecord, 1, from, to, chunk), `${from}-${to} in chunks of ${chunk}`).toBe(ragged.letters.slice(from - 1, to));
      }
    }
    for (const [from, to] of [
      [1, 70_001],
      [59, 62],
      [65_536, 70_001],
    ] as const) {
      expect(read(bytes, uniformRecord, 1, from, to, 1000)).toBe(uniform.letters.slice(from - 1, to));
    }
  });

  it('places residue i of a uniform record with the formula, for any positive eol', () => {
    for (const eol of ['\n', '\r\n', ' \n', '\t \r\n', '     \n']) {
      const writer = new FastaWriter(seeded(eol.length), 1);
      const written: WrittenRecord = writer.record('u', 'ACGT'.repeat(9), { kind: 'uniform', width: 8, eol });
      writer.bytes();
      written.offsets.forEach((offset, i) => expect(uniformOffset({ width: 8, eol: eol.length }, written.sequenceOffset, i)).toBe(offset));
    }
    expect(() => uniformOffset({ width: 0, eol: 1 }, 0, 0)).toThrow(RangeError);
  });

  it('keeps the bytes of the file: case and U as written', () => {
    const writer = new FastaWriter(seeded(1), 1);
    const target = writer.record('x', 'acgUuTNnrY', { kind: 'ragged', minWidth: 3, maxWidth: 3, eol: '\n' });
    const bytes = writer.bytes();
    expect(read(bytes, indexedRecord(target, 0, 1), 1, 1, 10, 2)).toBe('acgUuTNnrY');
  });

  it('follows the line rules of abi_v2.md §9 on hand-written lines', () => {
    const text = '>r\nAC GT\n  ;GG comment\n\t# CC\n ! TT\nA;CC rest\n\v\fG-1.T*X\r\nUa\rc\n';
    const bytes = new TextEncoder().encode(text);
    const kind1 = (from: number, to: number) => {
      const reader = new ForwardReader(1, from, to - from);
      reader.feed(bytes.subarray(3));
      return latin1.decode(reader.residues());
    };
    expect(kind1(0, 10)).toBe('ACGTAGTUac');
    const protein = new ForwardReader(2, 0, 20);
    protein.feed(bytes.subarray(3));
    expect(latin1.decode(protein.residues())).toBe('ACGTAGT*XUac');
  });

  it('counts residues as the record table does and knows the stored letters of each kind', () => {
    expect(residueCounts(1, new TextEncoder().encode('acgUuTN'))).toEqual(counts('acgUuTN', 1));
    expect(residueCounts(1, new TextEncoder().encode('Uu'))).toEqual({ T: 2 });
    expect(residueCounts(2, new TextEncoder().encode('Uu*m'))).toEqual({ U: 2, '*': 1, M: 1 });
    expect([...'ABCDGHKMNRSTUVWYabcdghkmnrstuvwy'].every((c) => stores(1, c.charCodeAt(0)))).toBe(true);
    expect([...'EFIJLOPQXZ*-0 ;'].some((c) => stores(1, c.charCodeAt(0)))).toBe(false);
    expect([...'AZaz*'].every((c) => stores(2, c.charCodeAt(0)))).toBe(true);
    expect([...'-0 ;@[`{'].some((c) => stores(2, c.charCodeAt(0)))).toBe(false);
    // Bytes 0x80-0xFF (UTF-8, a byte order mark) are never residues.
    for (let byte = 0x80; byte < 0x100; byte++) expect(stores(1, byte) || stores(2, byte), `0x${byte.toString(16)}`).toBe(false);
  });

  it('accepts reader kinds 1 and 2 only', () => {
    expect(readerKind(1)).toBe(1);
    expect(readerKind(2)).toBe(2);
    expect(() => readerKind(0)).toThrow(/kind 0 cannot be extracted/);
    expect(() => readerKind(3)).toThrow(RangeError);
  });

  it('finds the checkpoint at or before a residue', () => {
    const layout = { kind: 'checkpoints', every: 10, offsets: [5, 30, 52] } as const;
    expect(readStart(layout, 5, 0)).toEqual({ offset: 5, skip: 0 });
    expect(readStart(layout, 5, 19)).toEqual({ offset: 30, skip: 9 });
    expect(readStart(layout, 5, 20)).toEqual({ offset: 52, skip: 0 });
    expect(() => readStart(layout, 5, 30)).toThrow(/no checkpoint for residue 31/);
    expect(readStart({ kind: 'uniform', width: 4, eol: 2 }, 10, 5)).toEqual({ offset: 17, skip: 0 });
  });
});
