// RecordScanner contract: the ABI v2 index scan (docs/web/abi_v2.md §4, §9), parser kind
// 0. The expected values are what the adapter's scan (web/adapter/src/scan.rs, which
// reproduces `bio::io::fasta` 1.6.0 and is checked by web/adapter/tests/scan_properties.rs)
// reports; the adapter is the authority (plan TD-8). The FakeScanner runs these cases in
// Vitest; S09 runs them against the adapter's serial reactor. The layout kind is not
// prescribed: the cases check that the layout finds every residue.
import type { IndexedRecord } from '../../src/domain/dataset';
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
  /** The sequence as the parser reports it. */
  readonly residues: string;
}

interface Corpus {
  readonly name: string;
  readonly input: () => Uint8Array;
  readonly expected: () => readonly ExpectedRecord[];
}

const encoder = new TextEncoder();
const text = (value: string) => () => encoder.encode(value);
const records = (...list: ExpectedRecord[]) => () => list;
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
    expected: records({
      id: 'long',
      header_offset: 0,
      sequence_offset: header.length,
      end_offset: input.length,
      residues: sequence,
    }),
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
    name: 'a blank line inside a sequence',
    input: text('>a\nAC\n\nGT\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 10, residues: 'ACGT' }),
  },
  {
    name: 'interior white space belongs to the sequence, trailing white space does not',
    input: text('>a\nAC GT  \nA\t\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 14, residues: 'AC GTA' }),
  },
  {
    name: 'a non-ASCII header',
    input: text('>é1 d\nAC\n'),
    expected: records({ id: 'é1', header_offset: 0, sequence_offset: 7, end_offset: 10, residues: 'AC' }),
  },
  {
    name: 'residues keep their case',
    input: text('>a\nacgtN\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 9, residues: 'acgtN' }),
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
    name: 'an empty record ends the input (bio stops reading there)',
    input: text('>a\nAC\n>\n>b\nGG\n'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 6, residues: 'AC' }),
  },
  {
    name: 'a header line without a newline at the end of the input',
    input: text('>a'),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 2, end_offset: 2, residues: '' }),
  },
  { name: 'an empty input has no records', input: text(''), expected: records() },
  {
    name: 'a U+FEFF at the start of a header is part of the ID',
    input: text(`>${String.fromCharCode(0xfeff)}a\nAC\n`),
    expected: records({ id: `${String.fromCharCode(0xfeff)}a`, header_offset: 0, sequence_offset: 6, end_offset: 9, residues: 'AC' }),
  },
  {
    name: 'bytes after the empty record that ends the input are not read',
    input: () => new Uint8Array([...encoder.encode('>a\nAC\n>\n>b\nGG\n'), 0xff]),
    expected: records({ id: 'a', header_offset: 0, sequence_offset: 3, end_offset: 6, residues: 'AC' }),
  },
  { name: 'a long record (4 checkpoints) with lines of 60 and 61 residues', ...longRecord(3 * 65_536 + 100, [60, 61]) },
  { name: 'a long record with regular 80-residue lines', ...longRecord(200_000, [80]) },
];

const ERRORS: ReadonlyArray<{ readonly name: string; readonly input: Uint8Array; readonly error: RegExp }> = [
  { name: 'text before the first header', input: encoder.encode('x\n>a\nAC\n'), error: /Expected > at record start\./ },
  { name: 'a blank line before the first header', input: encoder.encode('\n>a\nAC\n'), error: /Expected > at record start\./ },
  { name: 'white space only', input: encoder.encode(' \n'), error: /Expected > at record start\./ },
  {
    name: 'a first line that is not a header, before bytes that are not UTF-8',
    input: new Uint8Array([0x78, 0x0a, 0xff]),
    error: /Expected > at record start\./,
  },
  {
    name: 'bytes that are not UTF-8',
    input: new Uint8Array([...encoder.encode('>a\nAC'), 0xff, ...encoder.encode('GT\n')]),
    error: /stream did not contain valid UTF-8/,
  },
];

async function* chunked(bytes: Uint8Array, size: number): AsyncGenerator<Uint8Array> {
  if (size === 0) {
    yield bytes;
    return;
  }
  for (let at = 0; at < bytes.length; at += size) yield bytes.slice(at, at + size);
}

/** Reads residues forward from `offset` under parser kind 0: each line without its trailing white space. */
function readForward(bytes: Uint8Array, offset: number, end: number, count: number): number[] {
  const out: number[] = [];
  let position = offset;
  while (out.length < count && position < end) {
    let lineEnd = bytes.indexOf(0x0a, position);
    if (lineEnd < 0 || lineEnd > end) lineEnd = end;
    let keep = lineEnd;
    while (keep > position && (bytes[keep - 1] === 0x20 || (bytes[keep - 1]! >= 0x09 && bytes[keep - 1]! <= 0x0d))) keep--;
    for (let i = position; i < keep && out.length < count; i++) out.push(bytes[i]!);
    position = lineEnd + 1;
  }
  return out;
}

/** The residues of a record, found through its line layout. */
function locate(bytes: Uint8Array, record: IndexedRecord): Uint8Array {
  const out = new Uint8Array(record.length);
  const layout = record.line_layout;
  if (layout.kind === 'uniform') {
    check(layout.width > 0 || record.length === 0, `record ${record.id}: a uniform layout needs a width`);
    for (let i = 0; i < record.length; i++) {
      const at = record.sequence_offset + Math.floor(i / layout.width) * (layout.width + layout.eol) + (i % layout.width);
      check(at < record.end_offset, `record ${record.id}: residue ${i} lies outside the record`);
      out[i] = bytes[at]!;
    }
    return out;
  }
  check(layout.every > 0, `record ${record.id}: checkpoints need a positive spacing`);
  same(layout.offsets.length, Math.ceil(record.length / layout.every), `record ${record.id}: checkpoint count`);
  layout.offsets.forEach((offset, k) => {
    const first = k * layout.every;
    out.set(readForward(bytes, offset, record.end_offset, Math.min(layout.every, record.length - first)), first);
  });
  return out;
}

function countBytes(bytes: Uint8Array): Record<string, number> {
  const counts: Record<string, number> = {};
  for (const byte of bytes) {
    const key = byte >= 0x21 && byte <= 0x7e ? String.fromCharCode(byte) : `0x${byte.toString(16).padStart(2, '0')}`;
    counts[key] = (counts[key] ?? 0) + 1;
  }
  return counts;
}

function checkResponse(bytes: Uint8Array, response: ScanResponse, expected: readonly ExpectedRecord[]): void {
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
    sameBytes(locate(bytes, record), residues, `record ${want.id} residues found through the line layout`);
  });
}

export const RECORD_SCANNER_CASES: readonly ContractCase<RecordScannerEnv>[] = [
  ...CORPUS.map(
    (corpus): ContractCase<RecordScannerEnv> => ({
      name: `scan: ${corpus.name}`,
      async run({ scanner }) {
        const bytes = corpus.input();
        const whole = await scanner.scan(0, chunked(bytes, 0));
        checkResponse(bytes, whole, corpus.expected());
        for (const size of CHUNKINGS.slice(1)) {
          sameValue(await scanner.scan(0, chunked(bytes, size)), whole, `the response for chunks of ${size} bytes`);
        }
      },
    }),
  ),
  ...ERRORS.map(
    (error): ContractCase<RecordScannerEnv> => ({
      name: `scan: the parser's error for ${error.name}`,
      async run({ scanner }) {
        for (const size of CHUNKINGS) {
          await rejects(scanner.scan(0, chunked(error.input, size)), error.error, `chunks of ${size} bytes`);
        }
      },
    }),
  ),
  {
    name: 'scan: a multi-byte character split across chunks',
    async run({ scanner }) {
      const bytes = encoder.encode('>éé x\nAC\n');
      for (const size of [1, 2, 3]) {
        const response = await scanner.scan(0, chunked(bytes, size));
        same(response.records[0]?.id, 'éé', `id with chunks of ${size} bytes`);
      }
    },
  },
  {
    // Parser kind 1 (BLASTX's reader) joins in SX; SX replaces this case with its corpus.
    name: 'scan: parser kind 1 is refused until SX',
    async run({ scanner }) {
      await rejects(scanner.scan(1, chunked(encoder.encode('>a\nAC\n'), 0)), /parser kind 1/, 'kind 1');
    },
  },
];
