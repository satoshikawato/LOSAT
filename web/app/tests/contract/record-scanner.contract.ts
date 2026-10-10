// RecordScanner contract: the ABI v2 index scan (docs/web/abi_v2.md §4, §9), parser kinds 1
// (NCBI BLAST+'s reader with the nucleotide flags) and 2 (with the protein flags). The
// expected values are what the adapter's scan (web/adapter/src/scan/ncbi.rs, checked against
// the engine's reader by web/adapter/tests/scan_ncbi_properties.rs) reports; the adapter is
// the authority (plan TD-8). The FakeScanner runs these cases in Vitest, and the serial
// reactor runs them in Node (tests/unit/engine-runtime.test.ts) and in a browser worker. The
// cases are those on which the two agree: the reader's errors and LOSAT Web's rejections are
// the reactor's only (engine-runtime.test.ts). The layout kind is not prescribed: every
// record's residues are read through its layout by a walk of this file's own.
import type { FastaParserKind, IndexedRecord } from '../../src/domain/dataset';
import type { RecordScanner, ScanResponse } from '../../src/ports/scan';
import { check, rejects, same, sameValue, type ContractCase } from './contract';

export interface RecordScannerEnv {
  readonly scanner: RecordScanner;
}

interface ExpectedRecord {
  readonly id: string;
  readonly header_offset: number;
  readonly sequence_offset: number;
  readonly end_offset: number;
  /** The stored residues, upper-cased (`U` as `T` in kind 1); one string, or one for each kind. */
  readonly residues: string | Readonly<Record<FastaParserKind, string>>;
}

interface Corpus {
  readonly name: string;
  readonly kinds: readonly FastaParserKind[];
  readonly input: () => Uint8Array;
  readonly expected: () => readonly ExpectedRecord[];
}

const BOTH: readonly FastaParserKind[] = [1, 2];
const encoder = new TextEncoder();
const text = (value: string) => () => encoder.encode(value);
const records = (...list: ExpectedRecord[]) => () => list;
const CHUNKINGS = [0, 1, 3, 4096] as const;
const CR = 0x0d;
const LF = 0x0a;

/** A long record whose lines alternate between two widths (or have one width). */
function longRecord(residues: number, widths: readonly number[]): Pick<Corpus, 'input' | 'expected'> {
  const letters = 'ACGTNacgtn';
  let sequence = '';
  for (let i = 0; i < residues; i++) sequence += letters[(i * 7 + (i >>> 5)) % letters.length];
  const header = '>long record\n';
  let body = '';
  for (let at = 0, line = 0; at < residues; line++) {
    const width = widths[line % widths.length]!;
    body += `${sequence.slice(at, at + width)}\n`;
    at += width;
  }
  const input = header + body;
  return {
    input: text(input),
    expected: records({
      id: 'long',
      header_offset: 0,
      sequence_offset: header.length,
      end_offset: input.length,
      residues: sequence.toUpperCase(),
    }),
  };
}

const CORPUS: readonly Corpus[] = [
  {
    name: 'several records with regular lines (LF)',
    kinds: BOTH,
    input: text('>a desc\nACGT\nAC\n>b\nGGCC\n'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 8, end_offset: 16, residues: 'ACGTAC' },
      { id: 'b', header_offset: 16, sequence_offset: 19, end_offset: 24, residues: 'GGCC' },
    ),
  },
  {
    name: 'CR LF line ends',
    kinds: BOTH,
    input: text('>a\r\nACGT\r\nAC\r\n>b x\r\nGG\r\n'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 4, end_offset: 14, residues: 'ACGTAC' },
      { id: 'b', header_offset: 14, sequence_offset: 20, end_offset: 24, residues: 'GG' },
    ),
  },
  {
    name: 'CR line ends',
    kinds: BOTH,
    input: text('>a\rACGT\rAC\r>b x\rGG\r'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 11, residues: 'ACGTAC' },
      { id: 'b', header_offset: 11, sequence_offset: 16, end_offset: 19, residues: 'GG' },
    ),
  },
  {
    name: 'LF and CR LF line ends in one file',
    kinds: BOTH,
    input: text('>a\nACGT\r\nAC\n>b\r\nGG\n'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 12, residues: 'ACGTAC' },
      { id: 'b', header_offset: 12, sequence_offset: 16, end_offset: 19, residues: 'GG' },
    ),
  },
  {
    name: 'lower case residues are counted upper-cased',
    kinds: BOTH,
    input: text('>a\nacgtNn\nAcGt\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 15, residues: 'ACGTNNACGT' }),
  },
  {
    name: 'comment lines (!, # and ;) before and inside records',
    kinds: BOTH,
    input: text('#c\n>a\n!x\nACGT\n#y\n ;z\nGGCC\n'),
    expected: records({ id: 'a', header_offset: 3, sequence_offset: 6, end_offset: 26, residues: 'ACGTGGCC' }),
  },
  {
    name: '; ends the data of a line',
    kinds: BOTH,
    input: text('>a\nACGTAC;GT\n\tGGCC ;TT\nAA\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 26, residues: 'ACGTACGGCCAA' }),
  },
  {
    name: 'blank lines and white space',
    kinds: BOTH,
    input: text('>a\n\nAC GT\n \t\nGG\n\n>b\n\n'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 17, residues: 'ACGTGG' },
      { id: 'b', header_offset: 17, sequence_offset: 20, end_offset: 21, residues: '' },
    ),
  },
  {
    name: 'a first record without a defline',
    kinds: BOTH,
    input: text('ACGTACGTAC\nGG\n>b\nTT\n'),
    expected: records(
      { id: '', header_offset: 0, sequence_offset: 0, end_offset: 14, residues: 'ACGTACGTACGG' },
      { id: 'b', header_offset: 14, sequence_offset: 17, end_offset: 20, residues: 'TT' },
    ),
  },
  {
    name: 'a first record without a defline after a comment line',
    kinds: BOTH,
    input: text('#c\nACGT\n>b\nAC\n'),
    expected: records(
      { id: '', header_offset: 0, sequence_offset: 0, end_offset: 8, residues: 'ACGT' },
      { id: 'b', header_offset: 8, sequence_offset: 11, end_offset: 14, residues: 'AC' },
    ),
  },
  {
    name: 'records without residues, the last a defline without a newline',
    kinds: BOTH,
    input: text('>a\n>b\nAC\n>c'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 3, residues: '' },
      { id: 'b', header_offset: 3, sequence_offset: 6, end_offset: 9, residues: 'AC' },
      { id: 'c', header_offset: 9, sequence_offset: 11, end_offset: 11, residues: '' },
    ),
  },
  {
    name: 'U is stored as T in kind 1',
    kinds: BOTH,
    input: text('>a\nACGUuacgu\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 13, residues: { 1: 'ACGTTACGT', 2: 'ACGUUACGU' } }),
  },
  {
    name: '* and the protein letters are stored in kind 2 only',
    kinds: BOTH,
    input: text('>p\nMKV*LLe*\n'),
    expected: records({ id: 'p', header_offset: 0, sequence_offset: 3, end_offset: 12, residues: { 1: 'MKV', 2: 'MKV*LLE*' } }),
  },
  {
    name: 'bytes that the kind does not store are skipped',
    kinds: BOTH,
    input: text('>a\nAC-GT.12xx\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 14, residues: { 1: 'ACGT', 2: 'ACGTXX' } }),
  },
  {
    name: 'the ID is the title up to its first space, after the white space that follows >',
    kinds: BOTH,
    input: text('>  id1 some description\nAC\n'),
    expected: records({ id: 'id1', header_offset: 0, sequence_offset: 24, end_offset: 27, residues: 'AC' }),
  },
  {
    name: 'a non-ASCII title',
    kinds: BOTH,
    input: text('>é1 d\nAC\n'),
    expected: records({ id: 'é1', header_offset: 0, sequence_offset: 7, end_offset: 10, residues: 'AC' }),
  },
  { name: 'an empty input has no records', kinds: BOTH, input: text(''), expected: records() },
  {
    name: 'white space, blank and comment lines only have no records',
    kinds: BOTH,
    input: text(' \n\n#c\n;d\n!e\n'),
    expected: records(),
  },
  { name: 'a long record (4 checkpoints) with lines of 60 and 61 residues', kinds: BOTH, ...longRecord(3 * 65_536 + 100, [60, 61]) },
  { name: 'a long record with regular 80-residue lines', kinds: BOTH, ...longRecord(200_000, [80]) },
];

async function* chunked(bytes: Uint8Array, size: number): AsyncGenerator<Uint8Array> {
  if (size === 0) {
    yield bytes;
    return;
  }
  for (let at = 0; at < bytes.length; at += size) yield bytes.slice(at, at + size);
}

const NUCLEOTIDE_LETTERS = 'ABCDGHKMNRSTUVWY';

/** The residue that `kind` stores for a byte (upper-cased, `U` as `T` in kind 1), or undefined. */
function residueOf(byte: number, kind: FastaParserKind): string | undefined {
  const character = String.fromCharCode(byte);
  if (!/^[A-Za-z*]$/.test(character)) return undefined;
  const letter = character.toUpperCase();
  if (kind === 2) return letter;
  if (!NUCLEOTIDE_LETTERS.includes(letter)) return undefined;
  return letter === 'U' ? 'T' : letter;
}

/**
 * Reads `count` residues forward from `offset` under the rules of `kind` (abi_v2.md §9): CR
 * and LF end a line; at the start of a line space, tab, VT and FF are skipped, and the whole
 * line when its first other byte is `!`, `#` or `;`; elsewhere `;` skips the rest of the line;
 * every byte that the kind does not store is skipped.
 */
function readForward(bytes: Uint8Array, offset: number, end: number, count: number, kind: FastaParserKind): string[] {
  const out: string[] = [];
  let lineStart = false;
  let skipLine = false;
  for (let at = offset; at < end && out.length < count; at++) {
    const byte = bytes[at]!;
    if (byte === CR || byte === LF) {
      lineStart = true;
      skipLine = false;
      continue;
    }
    if (skipLine) continue;
    if (lineStart) {
      if (byte === 0x20 || byte === 0x09 || byte === 0x0b || byte === 0x0c) continue;
      lineStart = false;
      if (byte === 0x21 || byte === 0x23 || byte === 0x3b) {
        skipLine = true;
        continue;
      }
    }
    if (byte === 0x3b) {
      skipLine = true;
      continue;
    }
    const residue = residueOf(byte, kind);
    if (residue !== undefined) out.push(residue);
  }
  return out;
}

/** The residues of a record, found through its line layout. */
function locate(bytes: Uint8Array, record: IndexedRecord, kind: FastaParserKind): string[] {
  const layout = record.line_layout;
  const where = `record ${record.index + 1} (${record.id})`;
  const out: string[] = [];
  if (layout.kind === 'uniform') {
    check(layout.width > 0 || record.length === 0, `${where}: a uniform layout needs a width`);
    check(layout.eol > 0, `${where}: a uniform layout needs a positive eol`);
    for (let i = 0; i < record.length; i++) {
      const at = record.sequence_offset + Math.floor(i / layout.width) * (layout.width + layout.eol) + (i % layout.width);
      check(at < record.end_offset, `${where}: residue ${i} lies outside the record`);
      const residue = residueOf(bytes[at]!, kind);
      check(residue !== undefined, `${where}: residue ${i} is at a byte that kind ${kind} does not store`);
      out.push(residue);
    }
    return out;
  }
  check(layout.every > 0, `${where}: checkpoints need a positive spacing`);
  same(layout.offsets.length, Math.ceil(record.length / layout.every), `${where}: checkpoint count`);
  layout.offsets.forEach((offset, k) => {
    check(residueOf(bytes[offset]!, kind) !== undefined, `${where}: checkpoint ${k} is not at a residue`);
    const count = Math.min(layout.every, record.length - k * layout.every);
    for (const residue of readForward(bytes, offset, record.end_offset, count, kind)) out.push(residue);
  });
  return out;
}

function countResidues(residues: readonly string[]): Record<string, number> {
  const counts: Record<string, number> = {};
  for (const residue of residues) counts[residue] = (counts[residue] ?? 0) + 1;
  return counts;
}

function checkResponse(bytes: Uint8Array, response: ScanResponse, expected: readonly ExpectedRecord[], kind: FastaParserKind): void {
  same(response.records.length, expected.length, 'record count');
  response.records.forEach((record, i) => {
    const want = expected[i]!;
    const residues = [...(typeof want.residues === 'string' ? want.residues : want.residues[kind])];
    same(record.index, i, `record ${i + 1} index`);
    same(record.id, want.id, `record ${i + 1} id`);
    same(record.header_offset, want.header_offset, `record ${want.id} header_offset`);
    same(record.sequence_offset, want.sequence_offset, `record ${want.id} sequence_offset`);
    same(record.end_offset, want.end_offset, `record ${want.id} end_offset`);
    same(record.length, residues.length, `record ${want.id} length`);
    sameValue(record.residue_counts, countResidues(residues), `record ${want.id} residue_counts`);
    // The walk of this file finds `length` residues whose counts are `residue_counts`.
    const located = locate(bytes, record, kind);
    same(located.length, record.length, `record ${want.id}: the number of residues found through the line layout`);
    sameValue(countResidues(located), record.residue_counts, `record ${want.id}: the counts of the residues found through the line layout`);
    same(located.join(''), residues.join(''), `record ${want.id}: the residues found through the line layout`);
  });
}

export const RECORD_SCANNER_CASES: readonly ContractCase<RecordScannerEnv>[] = [
  ...CORPUS.flatMap((corpus) =>
    corpus.kinds.map(
      (kind): ContractCase<RecordScannerEnv> => ({
        name: `scan kind ${kind}: ${corpus.name}`,
        async run({ scanner }) {
          const bytes = corpus.input();
          const whole = await scanner.scan(kind, chunked(bytes, 0));
          checkResponse(bytes, whole, corpus.expected(), kind);
          for (const size of CHUNKINGS.slice(1)) {
            sameValue(await scanner.scan(kind, chunked(bytes, size)), whole, `the response for chunks of ${size} bytes`);
          }
        },
      }),
    ),
  ),
  {
    name: 'scan: a multi-byte character split across chunks',
    async run({ scanner }) {
      const bytes = encoder.encode('>éé x\nAC\n');
      for (const size of [1, 2, 3]) {
        const response = await scanner.scan(1, chunked(bytes, size));
        same(response.records[0]?.id, 'éé', `id with chunks of ${size} bytes`);
      }
    },
  },
  {
    // Kind 0 (the `bio::io::fasta` reader before session SF) stays in the adapter for itself;
    // the application uses kinds 1 and 2 only.
    name: 'scan: an unknown parser kind is refused',
    async run({ scanner }) {
      await rejects(scanner.scan(3 as unknown as FastaParserKind, chunked(encoder.encode('>a\nAC\n'), 0)), /parser kind 3/, 'kind 3');
    },
  },
];
