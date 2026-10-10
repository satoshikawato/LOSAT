// A fake of the ABI v2 index scan (docs/web/abi_v2.md §9) for tests and for the FakeEngine
// build. It is not the engine's FASTA reader and must not be used for anything else. It
// follows the rules of parser kinds 1 and 2 (NCBI BLAST+'s reader with the nucleotide or
// the protein flags) closely enough to pass tests/contract/record-scanner.contract.ts:
// - a record starts at a line whose first byte is `>`; residues before the first `>` line
//   make a first record without a defline (offsets 0 and 0);
// - the ID is the title up to its first space; the title follows `>` and its white space,
//   and ends at the first byte below 0x20, without trailing white space;
// - at the start of a line, space, tab, VT and FF are skipped, and the whole line when its
//   first other byte is `!`, `#` or `;`; elsewhere `;` skips the rest of the line;
// - the residues are the letters that the kind stores (kind 1: `ABCDGHKMNRSTUVWY` in either
//   case; kind 2: every ASCII letter and `*`), counted upper-cased with `U` as `T` in kind 1;
//   every other byte is skipped;
// - an input of white space, blank and comment lines has no record.
// Its simplifications: LF, CR LF and a lone CR each end one line (NCBI's reader joins two
// lines of a file that mixes CR with LF line ends, which the engine's scan then rejects); it
// never fails (NCBI's reader errors and LOSAT Web's rejections, a `>?` gap line among them,
// are the engine's); and it always reports the checkpoints layout, every 65536 residues.
import { recordKey, type FastaParserKind, type IndexedRecord, type RecordKey } from '../../domain/dataset';
import { indexParser, type InputRole, type ProgramId } from '../../domain/programs';
import type { InputCheck, InputChecker } from '../../ports/input-check';
import type { RecordScanner, ScanResponse } from '../../ports/scan';
import { concatBytes } from '../bytes';

const GT = 0x3e;
const LF = 0x0a;
const CR = 0x0d;
const SPACE = 0x20;
const BANG = 0x21;
const HASH = 0x23;
const STAR = 0x2a;
const SEMICOLON = 0x3b;
const CHECKPOINT_EVERY = 65_536;
/** The letters that kind 1 stores, upper-cased. */
const NUCLEOTIDE = new Set([...'ABCDGHKMNRSTUVWY'].map((letter) => letter.charCodeAt(0)));

export class FakeScanner implements RecordScanner {
  async scan(parser: number, chunks: AsyncIterable<Uint8Array>): Promise<ScanResponse> {
    if (parser !== 1 && parser !== 2) throw new Error(`unknown or unavailable FASTA parser kind ${parser}`);
    const parts: Uint8Array[] = [];
    for await (const chunk of chunks) parts.push(chunk.slice());
    return { records: fakeScan(concatBytes(parts), parser) };
  }
}

/** The records that the fake reads with `kind`, as the FakeEngine's `register` reports them. */
export function fakeRecordKeys(bytes: Uint8Array, kind: FastaParserKind): RecordKey[] {
  return fakeScan(bytes, kind).map(recordKey);
}

export function fakeScan(bytes: Uint8Array, kind: FastaParserKind): IndexedRecord[] {
  return scanRecords(bytes, kind).map(({ record }) => record);
}

interface Builder {
  readonly index: number;
  readonly id: string;
  readonly headerOffset: number;
  readonly sequenceOffset: number;
  length: number;
  readonly counts: Map<number, number>;
  readonly checkpoints: number[];
  /** A data line has `!` (the FakeInputChecker's refusal). */
  bang: boolean;
}

interface Scanned {
  readonly record: IndexedRecord;
  readonly bang: boolean;
}

function scanRecords(bytes: Uint8Array, kind: FastaParserKind): Scanned[] {
  const records: Scanned[] = [];
  let current: Builder | undefined;
  let offset = 0;
  while (offset < bytes.length) {
    let end = offset;
    while (end < bytes.length && bytes[end] !== LF && bytes[end] !== CR) end++;
    const next = end === bytes.length ? end : end + (bytes[end] === CR && bytes[end + 1] === LF ? 2 : 1);
    if (bytes[offset] === GT) {
      if (current !== undefined) records.push(finish(current, offset));
      current = startRecord(bytes, offset, end, next, records.length);
    } else {
      let at = offset;
      while (at < end && isLeadingSpace(bytes[at]!)) at++;
      const first = bytes[at];
      if (at < end && first !== BANG && first !== HASH && first !== SEMICOLON) {
        // Data before the first defline: the input's first record, without a defline.
        current ??= newBuilder(0, '', 0, 0);
        for (; at < end && bytes[at] !== SEMICOLON; at++) addByte(current, bytes[at]!, at, kind);
      }
    }
    offset = next;
  }
  if (current !== undefined) records.push(finish(current, bytes.length));
  return records;
}

function newBuilder(index: number, id: string, headerOffset: number, sequenceOffset: number): Builder {
  return { index, id, headerOffset, sequenceOffset, length: 0, counts: new Map(), checkpoints: [], bang: false };
}

function startRecord(bytes: Uint8Array, start: number, lineEnd: number, next: number, index: number): Builder {
  let from = start + 1;
  while (from < lineEnd && isLeadingSpace(bytes[from]!)) from++;
  let to = Math.min(from + 1, lineEnd);
  while (to < lineEnd && bytes[to]! >= SPACE) to++;
  while (to > from && bytes[to - 1] === SPACE) to--;
  const title = bytes.subarray(from, to);
  const space = title.indexOf(SPACE);
  const id = new TextDecoder('utf-8', { fatal: false }).decode(space < 0 ? title : title.subarray(0, space));
  return newBuilder(index, id, start, next);
}

function addByte(builder: Builder, byte: number, offset: number, kind: FastaParserKind): void {
  if (byte === BANG) builder.bang = true;
  const residue = storedResidue(byte, kind);
  if (residue === undefined) return;
  if (builder.length % CHECKPOINT_EVERY === 0) builder.checkpoints.push(offset);
  builder.counts.set(residue, (builder.counts.get(residue) ?? 0) + 1);
  builder.length++;
}

/** The residue that the kind stores for a byte, upper-cased, or undefined for a byte it skips. */
function storedResidue(byte: number, kind: FastaParserKind): number | undefined {
  const upper = byte >= 0x61 && byte <= 0x7a ? byte - 0x20 : byte;
  if (kind === 1) {
    if (!NUCLEOTIDE.has(upper)) return undefined;
    return upper === 0x55 ? 0x54 : upper; // U is stored as T
  }
  return (upper >= 0x41 && upper <= 0x5a) || byte === STAR ? upper : undefined;
}

function finish(builder: Builder, end: number): Scanned {
  const residueCounts: Record<string, number> = {};
  for (const [residue, count] of [...builder.counts].sort(([a], [b]) => a - b)) {
    residueCounts[String.fromCharCode(residue)] = count;
  }
  return {
    record: {
      index: builder.index,
      id: builder.id,
      header_offset: builder.headerOffset,
      sequence_offset: builder.sequenceOffset,
      end_offset: end,
      length: builder.length,
      line_layout: { kind: 'checkpoints', every: CHECKPOINT_EVERY, offsets: builder.checkpoints },
      residue_counts: residueCounts,
    },
    bang: builder.bang,
  };
}

/** Space, tab, VT and FF: skipped at the start of a line and of a title. */
function isLeadingSpace(byte: number): boolean {
  return byte === SPACE || byte === 0x09 || byte === 0x0b || byte === 0x0c;
}

/**
 * The input check of the development build (ports/input-check.ts). It is not the engine's:
 * it accepts what the FakeScanner reads with the role's kind, except a record with `!` in a
 * data line, which it refuses with a message shaped as the engine's (the record's number and
 * ID), so that the screen's handling of a refused record can be tried without the engine.
 */
export class FakeInputChecker implements InputChecker {
  async check(program: ProgramId, role: InputRole, bytes: Uint8Array): Promise<InputCheck> {
    if (program === 'blastx') return { ok: false, message: 'blastx is not available in LOSAT Web ABI v2 yet' };
    const records = scanRecords(bytes, indexParser(program, role));
    const refused = records.find(({ bang }) => bang)?.record;
    if (refused !== undefined) {
      return {
        ok: false,
        message: `${role} record ${refused.index + 1} (${refused.id}) has '!' in its sequence (FAKE ENGINE check, not LOSAT's)`,
      };
    }
    return { ok: true, records: records.map(({ record }) => recordKey(record)) };
  }
}
