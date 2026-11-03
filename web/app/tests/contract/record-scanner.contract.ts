// RecordScanner contract: the ABI v2 index scan (docs/web/abi_v2.md §4, §9), parser kinds 1
// (NCBI BLAST+'s reader with the nucleotide flags) and 2 (with the protein flags). The
// expected values are what the adapter's scan (web/adapter/src/scan/ncbi.rs, checked against
// the engine's reader by web/adapter/tests/scan_ncbi_properties.rs) reports; the adapter is
// the authority (plan TD-8). The FakeScanner runs these cases in Vitest, and
// tests/unit/engine-runtime.test.ts runs them against the adapter's serial reactor. The
// cases are those that the FakeScanner follows too (its simplifications are listed in
// src/infra/fake/fake-fasta.ts). The layout kind is not prescribed: the cases check that the
// layout finds every residue.
import type { FastaParserKind, IndexedRecord } from '../../src/domain/dataset';
import type { RecordScanner, ScanResponse } from '../../src/ports/scan';
import { check, rejects, same, sameBytes, sameValue, type ContractCase } from './contract';

export interface RecordScannerEnv {
  readonly scanner: RecordScanner;
}

interface ExpectedRecord {
  readonly id: string;
  readonly header_offset: number;
  readonly sequence_offset: number;
  readonly end_offset: number;
  /** The residues that the kind stores, upper case (`U` as `T` in kind 1). */
  readonly residues: string;
}

interface Corpus {
  readonly name: string;
  readonly input: () => Uint8Array;
  readonly expected: (kind: FastaParserKind) => readonly ExpectedRecord[];
}

const KINDS: readonly FastaParserKind[] = [1, 2];
const encoder = new TextEncoder();
const text = (value: string) => () => encoder.encode(value);
/** The same records for both kinds. */
const records =
  (...list: ExpectedRecord[]) =>
  () =>
    list;
/** The records of kind 1 and of kind 2. */
const byKind =
  (kind1: readonly ExpectedRecord[], kind2: readonly ExpectedRecord[]) =>
  (kind: FastaParserKind) =>
    kind === 1 ? kind1 : kind2;
const CHUNKINGS = [0, 1, 3, 4096] as const;

/** A long record whose lines alternate between two widths (or have one width). */
function longRecord(residues: number, widths: readonly number[]): { input: () => Uint8Array; expected: () => ExpectedRecord[] } {
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
    expected: () => [
      {
        id: 'long',
        header_offset: 0,
        sequence_offset: header.length,
        end_offset: input.length,
        residues: sequence.toUpperCase(),
      },
    ],
  };
}

const CORPUS: readonly Corpus[] = [
  {
    name: 'one record with regular lines (LF)',
    input: text('>a desc\nACGT\nAC\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 8, end_offset: 16, residues: 'ACGTAC' }),
  },
  {
    name: 'CRLF line ends',
    input: text('>a\r\nACGT\r\nAC\r\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 4, end_offset: 14, residues: 'ACGTAC' }),
  },
  {
    name: 'CR line ends',
    input: text('>a\rACGT\rAC\r'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 11, residues: 'ACGTAC' }),
  },
  {
    name: 'lines of different lengths',
    input: text('>a\nAC\nACGT\nA\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 13, residues: 'ACACGTA' }),
  },
  {
    name: 'several records, the last without a final newline',
    input: text('>a\nAAAA\n>b x\nCC'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 8, residues: 'AAAA' },
      { id: 'b', header_offset: 8, sequence_offset: 13, end_offset: 15, residues: 'CC' },
    ),
  },
  {
    name: 'a blank line and a line of white space inside a sequence',
    input: text('>a\nAC\n\n  \nGT\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 13, residues: 'ACGT' }),
  },
  {
    name: 'white space, digits and hyphens around and inside a sequence line are not residues',
    input: text('>a\n\t1 AC-GT 5\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 14, residues: 'ACGT' }),
  },
  {
    name: 'a non-ASCII title',
    input: text('>é1 d\nAC\n'),
    expected: records({ id: 'é1', header_offset: 0, sequence_offset: 7, end_offset: 10, residues: 'AC' }),
  },
  {
    name: 'a title byte that is not UTF-8 is U+FFFD in the ID',
    input: () => new Uint8Array([...encoder.encode('>a'), 0xff, ...encoder.encode('b c\nAC\n')]),
    expected: records({ id: 'a�b', header_offset: 0, sequence_offset: 7, end_offset: 10, residues: 'AC' }),
  },
  {
    name: 'a U+FEFF right after > is part of the ID',
    input: () => new Uint8Array([0x3e, 0xef, 0xbb, 0xbf, ...encoder.encode('abc d\nAC\n')]),
    expected: records({ id: '\ufeffabc', header_offset: 0, sequence_offset: 10, end_offset: 13, residues: 'AC' }),
  },
  {
    name: 'white space after > is not part of the title, and a tab ends it',
    input: text('>  a b\nAC\n>c\td e\nGT\n'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 7, end_offset: 10, residues: 'AC' },
      { id: 'c', header_offset: 10, sequence_offset: 17, end_offset: 20, residues: 'GT' },
    ),
  },
  {
    name: 'a defline that starts with >?_ is read as one that starts with >',
    input: text('>?_a\nAC\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 5, end_offset: 8, residues: 'AC' }),
  },
  {
    name: 'a defline without a title',
    input: text('>\nAC\n'),
    expected: records({ id: '', header_offset: 0, sequence_offset: 2, end_offset: 5, residues: 'AC' }),
  },
  {
    name: 'residues are stored upper case; kind 1 stores U as T',
    input: text('>a\nacgu\n'),
    expected: byKind(
      [{ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 8, residues: 'ACGT' }],
      [{ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 8, residues: 'ACGU' }],
    ),
  },
  {
    name: 'kind 1 stores the nucleotide letters only, kind 2 every letter and *',
    input: text('>gene_x\nMKLVVLAAGGHHKLMKLVVLAAGG\n>p\nMK*LV\n'),
    expected: byKind(
      [
        { id: 'gene_x', header_offset: 0, sequence_offset: 8, end_offset: 33, residues: 'MKVVAAGGHHKMKVVAAGG' },
        { id: 'p', header_offset: 33, sequence_offset: 36, end_offset: 42, residues: 'MKV' },
      ],
      [
        { id: 'gene_x', header_offset: 0, sequence_offset: 8, end_offset: 33, residues: 'MKLVVLAAGGHHKLMKLVVLAAGG' },
        { id: 'p', header_offset: 33, sequence_offset: 36, end_offset: 42, residues: 'MK*LV' },
      ],
    ),
  },
  {
    name: 'a ! inside a sequence line is skipped',
    input: text('>a\nACGT!ACGT\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 13, residues: 'ACGTACGT' }),
  },
  {
    name: 'comment lines and a ; after a residue skip the rest of their line',
    input: text('>a\nAC\n# GG\n ;TT\n!CC\nGT;AA\nC\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 28, residues: 'ACGTC' }),
  },
  {
    name: 'a record without residues',
    input: text('>a\n>b\nAC\n'),
    expected: records(
      { id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 3, residues: '' },
      { id: 'b', header_offset: 3, sequence_offset: 6, end_offset: 9, residues: 'AC' },
    ),
  },
  {
    name: 'a defline without a newline at the end of the input',
    input: text('>a'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 2, end_offset: 2, residues: '' }),
  },
  { name: 'an empty input has no records', input: text(''), expected: records() },
  { name: 'white space and blank lines only have no records', input: text(' \n\n\t\n'), expected: records() },
  { name: 'comment lines only have no records', input: text('#c\n;c\n!c\n'), expected: records() },
  {
    name: 'a comment line and a blank line before the first defline',
    input: text('#comment\n\n>a\nAC\n'),
    expected: records({ id: 'a', header_offset: 10, sequence_offset: 13, end_offset: 16, residues: 'AC' }),
  },
  {
    name: 'a first record without a defline',
    input: text('ACGT\n>b\nAC\n'),
    expected: records(
      { id: '', header_offset: 0, sequence_offset: 0, end_offset: 5, residues: 'ACGT' },
      { id: 'b', header_offset: 5, sequence_offset: 8, end_offset: 11, residues: 'AC' },
    ),
  },
  {
    name: 'a sequence without a defline',
    input: text('ACGTACGTACGTACGT'),
    expected: records({ id: '', header_offset: 0, sequence_offset: 0, end_offset: 16, residues: 'ACGTACGTACGTACGT' }),
  },
  {
    name: 'a first line of letters that kind 1 does not store',
    input: text('x\n>a\nAC\n'),
    expected: byKind(
      [
        { id: '', header_offset: 0, sequence_offset: 0, end_offset: 2, residues: '' },
        { id: 'a', header_offset: 2, sequence_offset: 5, end_offset: 8, residues: 'AC' },
      ],
      [
        { id: '', header_offset: 0, sequence_offset: 0, end_offset: 2, residues: 'X' },
        { id: 'a', header_offset: 2, sequence_offset: 5, end_offset: 8, residues: 'AC' },
      ],
    ),
  },
  {
    name: 'a line that starts with white space and > is a data line',
    input: text(' >a\nAC\n'),
    expected: records({ id: '', header_offset: 0, sequence_offset: 0, end_offset: 7, residues: 'AAC' }),
  },
  {
    name: 'only the data lines before a record has a residue are checked',
    input: text('>a\nA$\n$$$$$$\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 13, residues: 'A' }),
  },
  { name: 'a long record (4 checkpoints) with lines of 60 and 61 residues', ...longRecord(3 * 65_536 + 100, [60, 61]) },
  { name: 'a long record with regular 80-residue lines', ...longRecord(200_000, [80]) },
];

const NOT_PLAUSIBLE = (line: number) =>
  new RegExp(`CFastaReader: Near line ${line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment\\.$`);
const GAP = (line: number) => new RegExp(`line ${line} is a gap line \\('>\\?'\\), .* gap lines are not supported by LOSAT Web$`);

const ERRORS: ReadonlyArray<{ readonly name: string; readonly input: Uint8Array; readonly error: RegExp }> = [
  { name: 'a gap line', input: encoder.encode('>a\nACGT\n>?10\nAC\n'), error: GAP(3) },
  { name: 'a gap line after white space', input: encoder.encode('>a\nACGT\n >?10\n'), error: GAP(3) },
  { name: 'a gap line at the end of the input', input: encoder.encode('>a\r\nAC\r\n>?'), error: GAP(3) },
  { name: 'a first data line that is not plausible data', input: encoder.encode('>a\n$$$$$$\nACGT\n'), error: NOT_PLAUSIBLE(2) },
  { name: 'a ; early in a first data line', input: encoder.encode('>a\nAC;x\nGT\n'), error: NOT_PLAUSIBLE(2) },
  { name: 'a later record whose first data line is not plausible data', input: encoder.encode('>a\nACGT\n>b\n$$$$$$\n'), error: NOT_PLAUSIBLE(4) },
  {
    name: 'a byte order mark',
    input: new Uint8Array([0xef, 0xbb, 0xbf, ...encoder.encode('>a\nAC\n')]),
    error: NOT_PLAUSIBLE(1),
  },
  {
    name: 'a byte that is not UTF-8 in a first data line',
    input: new Uint8Array([...encoder.encode('>a\nAC'), 0xff, ...encoder.encode('GT\n')]),
    error: NOT_PLAUSIBLE(2),
  },
];

async function* chunked(bytes: Uint8Array, size: number): AsyncGenerator<Uint8Array> {
  if (size === 0) {
    yield bytes;
    return;
  }
  for (let at = 0; at < bytes.length; at += size) yield bytes.slice(at, at + size);
}

/** The residue that the kind stores for a byte (upper case; `U` as `T` in kind 1), or undefined. */
function stored(byte: number, kind: FastaParserKind): number | undefined {
  const upper = byte >= 0x61 && byte <= 0x7a ? byte - 0x20 : byte;
  if (kind === 2) return (upper >= 0x41 && upper <= 0x5a) || upper === 0x2a ? upper : undefined;
  if (!'ABCDGHKMNRSTUVWY'.includes(String.fromCharCode(upper))) return undefined;
  return upper === 0x55 ? 0x54 : upper;
}

/**
 * Reads residues forward from a checkpoint (whose byte is a residue) under the rules of
 * kinds 1 and 2 (abi_v2.md §9): CR and LF end a line; a line whose first byte other than
 * space, tab, VT and FF is `!`, `#` or `;` is skipped; elsewhere a `;` skips the rest of the
 * line; a byte that the kind stores is the next residue; every other byte is skipped.
 */
function readForward(bytes: Uint8Array, offset: number, end: number, count: number, kind: FastaParserKind): number[] {
  const out: number[] = [];
  let lineStart = false;
  let skipLine = false;
  for (let i = offset; i < end && out.length < count; i++) {
    const byte = bytes[i]!;
    if (byte === 0x0a || byte === 0x0d) {
      lineStart = true;
      skipLine = false;
      continue;
    }
    if (skipLine) continue;
    if (lineStart) {
      if (byte === 0x20 || byte === 0x09 || byte === 0x0b || byte === 0x0c) continue;
      lineStart = false;
      if (byte === 0x21 || byte === 0x23) {
        skipLine = true;
        continue;
      }
    }
    if (byte === 0x3b) {
      skipLine = true;
      continue;
    }
    const residue = stored(byte, kind);
    if (residue !== undefined) out.push(residue);
  }
  return out;
}

/** The residues of a record, found through its line layout. */
function locate(bytes: Uint8Array, record: IndexedRecord, kind: FastaParserKind): Uint8Array {
  const out = new Uint8Array(record.length);
  const layout = record.line_layout;
  if (layout.kind === 'uniform') {
    check(layout.width > 0 || record.length === 0, `record ${record.id}: a uniform layout needs a width`);
    for (let i = 0; i < record.length; i++) {
      const at = record.sequence_offset + Math.floor(i / layout.width) * (layout.width + layout.eol) + (i % layout.width);
      check(at < record.end_offset, `record ${record.id}: residue ${i} lies outside the record`);
      out[i] = stored(bytes[at]!, kind) ?? 0;
    }
    return out;
  }
  check(layout.every > 0, `record ${record.id}: checkpoints need a positive spacing`);
  same(layout.offsets.length, Math.ceil(record.length / layout.every), `record ${record.id}: checkpoint count`);
  layout.offsets.forEach((offset, k) => {
    const first = k * layout.every;
    out.set(readForward(bytes, offset, record.end_offset, Math.min(layout.every, record.length - first), kind), first);
  });
  return out;
}

function countBytes(bytes: Uint8Array): Record<string, number> {
  const counts: Record<string, number> = {};
  for (const byte of bytes) {
    const key = String.fromCharCode(byte);
    counts[key] = (counts[key] ?? 0) + 1;
  }
  return counts;
}

function checkResponse(bytes: Uint8Array, response: ScanResponse, expected: readonly ExpectedRecord[], kind: FastaParserKind): void {
  same(response.records.length, expected.length, 'record count');
  response.records.forEach((record, i) => {
    const want = expected[i]!;
    const residues = encoder.encode(want.residues);
    same(record.index, i, `record ${i + 1} index`);
    same(record.id, want.id, `record ${i + 1} id`);
    same(record.header_offset, want.header_offset, `record ${want.id} header_offset`);
    same(record.sequence_offset, want.sequence_offset, `record ${want.id} sequence_offset`);
    same(record.end_offset, want.end_offset, `record ${want.id} end_offset`);
    same(record.length, residues.length, `record ${want.id} length`);
    sameValue(record.residue_counts, countBytes(residues), `record ${want.id} residue_counts`);
    sameBytes(locate(bytes, record, kind), residues, `record ${want.id} residues found through the line layout`);
  });
}

export const RECORD_SCANNER_CASES: readonly ContractCase<RecordScannerEnv>[] = [
  ...CORPUS.map(
    (corpus): ContractCase<RecordScannerEnv> => ({
      name: `scan: ${corpus.name}`,
      async run({ scanner }) {
        const bytes = corpus.input();
        for (const kind of KINDS) {
          const whole = await scanner.scan(kind, chunked(bytes, 0));
          checkResponse(bytes, whole, corpus.expected(kind), kind);
          for (const size of CHUNKINGS.slice(1)) {
            sameValue(await scanner.scan(kind, chunked(bytes, size)), whole, `kind ${kind}: the response for chunks of ${size} bytes`);
          }
        }
      },
    }),
  ),
  ...ERRORS.map(
    (error): ContractCase<RecordScannerEnv> => ({
      name: `scan: the reader's error for ${error.name}`,
      async run({ scanner }) {
        for (const kind of KINDS) {
          for (const size of CHUNKINGS) {
            await rejects(scanner.scan(kind, chunked(error.input, size)), error.error, `kind ${kind}, chunks of ${size} bytes`);
          }
        }
      },
    }),
  ),
  {
    name: 'scan: a multi-byte character split across chunks',
    async run({ scanner }) {
      const bytes = encoder.encode('>éé x\nAC\n');
      for (const kind of KINDS) {
        for (const size of [1, 2, 3]) {
          const response = await scanner.scan(kind, chunked(bytes, size));
          same(response.records[0]?.id, 'éé', `kind ${kind}: id with chunks of ${size} bytes`);
        }
      }
    },
  },
  {
    name: 'scan: an unknown parser kind is refused, kind 0 (bio\'s reader, removed) too',
    async run({ scanner }) {
      for (const kind of [0, 3]) {
        await rejects(scanner.scan(kind as FastaParserKind, chunked(encoder.encode('>a\nAC\n'), 0)), new RegExp(`parser kind ${kind}$`), `kind ${kind}`);
      }
    },
  },
];
