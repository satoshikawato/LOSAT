// A fake of the ABI v2 index scan (docs/web/abi_v2.md §9) for tests and for the
// development build. It is not the engine's FASTA reader and must not be used for anything
// else. It follows the rules of parser kinds 1 and 2 in the ABI document (NCBI BLAST+'s
// reader with the nucleotide or the protein flags) closely enough to pass
// tests/contract/record-scanner.contract.ts, with these simplifications: CR, LF and CR LF
// end a line, without the quirks of NCBI's line reader (which joins the rest of a line to the
// next after a CR in a file of LF line ends, or the reverse, and has buffer ends); it never
// takes a first line for a Seq-id that NCBI's data loaders would fetch, nor refuses a record
// over 2147483647 letters; and it always reports the checkpoints layout.
import { recordKey, type FastaParserKind, type IndexedRecord, type RecordKey } from '../../domain/dataset';
import { fastaParser, programById, type InputRole, type ProgramId } from '../../domain/programs';
import type { InputCheck, InputChecker } from '../../ports/input-check';
import type { RecordScanner, ScanResponse } from '../../ports/scan';
import { concatBytes } from '../bytes';

const GT = 0x3e;
const LF = 0x0a;
const CR = 0x0d;
const QUESTION = 0x3f;
const UNDERSCORE = 0x5f;
const SEMICOLON = 0x3b;
const CHECKPOINT_EVERY = 65_536;
/** `CheckDataLine`'s window. */
const CHECK_WINDOW = 70;
/** The letters that kind 1 stores, in either case. */
const NUCLEOTIDES = new Set([...'ABCDGHKMNRSTUVWY'].map((letter) => letter.charCodeAt(0)));

export class FakeScanner implements RecordScanner {
  async scan(parser: FastaParserKind, chunks: AsyncIterable<Uint8Array>): Promise<ScanResponse> {
    if (parser !== 1 && parser !== 2) throw new Error(`unknown or unavailable FASTA parser kind ${String(parser)}`);
    const parts: Uint8Array[] = [];
    for await (const chunk of chunks) parts.push(chunk.slice());
    return { records: fakeScan(concatBytes(parts), parser) };
  }
}

/** The records that the fake reads for a program's input, as the FakeEngine's `register` reports them. */
export function fakeRecordKeys(bytes: Uint8Array, program: ProgramId, role: InputRole): RecordKey[] {
  return fakeScan(bytes, fastaParser(programById(program), role)).map(recordKey);
}

interface Builder {
  readonly id: string;
  readonly headerOffset: number;
  readonly sequenceOffset: number;
  length: number;
  readonly counts: number[];
  readonly checkpoints: number[];
}

/** `isspace` of the C locale. */
function space(byte: number): boolean {
  return byte === 0x20 || (byte >= 0x09 && byte <= 0x0d);
}

/** The residue that the kind stores for a byte (upper case; `U` as `T` in kind 1), or undefined. */
function stored(byte: number, parser: FastaParserKind): number | undefined {
  const upper = byte >= 0x61 && byte <= 0x7a ? byte - 0x20 : byte;
  if (parser === 2) return (upper >= 0x41 && upper <= 0x5a) || upper === 0x2a ? upper : undefined;
  if (!NUCLEOTIDES.has(upper)) return undefined;
  return upper === 0x55 ? 0x54 : upper;
}

/**
 * How a line that starts with `>` at `at` begins: undefined for a gap line (`>?`), else
 * the number of bytes before its title (3 for `>?_`, which NCBI reads as `>`).
 */
function angle(bytes: Uint8Array, at: number, end: number): number | undefined {
  const byte = (i: number) => (at + i < end ? bytes[at + i] : undefined);
  if (byte(1) !== QUESTION) return 1;
  return byte(2) === UNDERSCORE && byte(3) !== QUESTION ? 3 : undefined;
}

function gapError(line: number): Error {
  return new Error(
    `line ${line} is a gap line ('>?'), whose residues have no bytes in the input for LOSAT's index to locate; gap lines are not supported in the browser`,
  );
}

/** The ID of a defline: its title (NCBI's `ParseDefLine`) up to the first space. */
function titleId(bytes: Uint8Array, start: number, end: number): string {
  let at = start;
  while (at < end && space(bytes[at]!)) at++;
  if (at === end) return '';
  let stop = at + 1;
  while (stop < end && bytes[stop]! >= 0x20) stop++;
  while (stop > at && space(bytes[stop - 1]!)) stop--;
  const title = bytes.subarray(at, stop);
  const firstSpace = title.indexOf(0x20);
  return new TextDecoder('utf-8', { ignoreBOM: true }).decode(firstSpace < 0 ? title : title.subarray(0, firstSpace));
}

/** NCBI's `CheckDataLine` of the first bytes of a data line (from its first byte other than white space). */
function plausible(line: readonly number[]): boolean {
  let end = 0;
  line.forEach((byte, i) => {
    if (!space(byte)) end = i < CHECK_WINDOW ? i + 1 : CHECK_WINDOW + 1;
  });
  const length = Math.min(end, CHECK_WINDOW);
  let good = 0;
  let bad = 0;
  for (const c of line.slice(0, length)) {
    if ((c >= 0x41 && c <= 0x5a) || (c >= 0x61 && c <= 0x7a) || c === 0x2a) good++;
    else if (c === 0x2d || space(c) || (c >= 0x30 && c <= 0x39)) continue;
    else if (c === SEMICOLON) break;
    else bad++;
  }
  return !(bad >= Math.floor(good / 3) && (length > 3 || good === 0 || bad > good));
}

export function fakeScan(bytes: Uint8Array, parser: FastaParserKind): IndexedRecord[] {
  const records: IndexedRecord[] = [];
  let current: Builder | undefined;
  let needDefline = true;
  let lineNumber = 0;
  const close = (end: number) => {
    if (current !== undefined) records.push(finish(current, records.length, end));
  };
  for (let start = 0; start < bytes.length; ) {
    let end = start;
    while (end < bytes.length && bytes[end] !== LF && bytes[end] !== CR) end++;
    let next = end;
    if (next < bytes.length) next += bytes[next] === CR && bytes[next + 1] === LF ? 2 : 1;
    lineNumber++;
    if (bytes[start] === GT) {
      const skip = angle(bytes, start, end);
      if (skip === undefined) throw gapError(lineNumber);
      close(start);
      needDefline = false;
      current = newBuilder(titleId(bytes, start + skip, end), start, next);
    } else {
      let at = start;
      while (at < end && space(bytes[at]!)) at++;
      const first = at < end ? bytes[at]! : undefined;
      if (first !== undefined && first !== 0x21 && first !== 0x23 && first !== SEMICOLON) {
        if (needDefline) {
          // The input's first record, without a defline.
          needDefline = false;
          current = newBuilder('', 0, 0);
        }
        dataLine(bytes, at, end, lineNumber, current!, parser);
      }
    }
    start = next;
  }
  close(bytes.length);
  return records;
}

function newBuilder(id: string, headerOffset: number, sequenceOffset: number): Builder {
  return { id, headerOffset, sequenceOffset, length: 0, counts: new Array<number>(256).fill(0), checkpoints: [] };
}

/** A data line from its first byte other than white space: its stored residues until a `;`. */
function dataLine(bytes: Uint8Array, at: number, end: number, line: number, builder: Builder, parser: FastaParserKind): void {
  let checked: number[] = [...bytes.subarray(at, end)];
  if (bytes[at] === GT) {
    const skip = angle(bytes, at, end);
    if (skip === undefined) throw gapError(line);
    checked = [GT, ...checked.slice(skip)];
  }
  // NCBI checks the data lines of a record until it has a residue.
  if (builder.length === 0 && !plausible(checked)) {
    throw new Error(`CFastaReader: Near line ${line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`);
  }
  for (let i = at; i < end && bytes[i] !== SEMICOLON; i++) {
    const residue = stored(bytes[i]!, parser);
    if (residue === undefined) continue;
    if (builder.length % CHECKPOINT_EVERY === 0) builder.checkpoints.push(i);
    builder.counts[residue]!++;
    builder.length++;
  }
}

function finish(builder: Builder, index: number, end: number): IndexedRecord {
  const residueCounts: Record<string, number> = {};
  builder.counts.forEach((count, byte) => {
    if (count > 0) residueCounts[String.fromCharCode(byte)] = count;
  });
  return {
    index,
    id: builder.id,
    header_offset: builder.headerOffset,
    sequence_offset: builder.sequenceOffset,
    end_offset: end,
    length: builder.length,
    line_layout: { kind: 'checkpoints', every: CHECKPOINT_EVERY, offsets: builder.checkpoints },
    residue_counts: residueCounts,
  };
}

/**
 * The input check of the development build (ports/input-check.ts). It is not the engine's:
 * it accepts what the FakeScanner reads with the program's kind, except a record with `!`
 * in its sequence lines, which it refuses with a message shaped as the engine's (the
 * record's number and ID), so that the screen's handling of a refused record can be tried
 * without the engine.
 */
export class FakeInputChecker implements InputChecker {
  async check(program: ProgramId, role: InputRole, bytes: Uint8Array): Promise<InputCheck> {
    let records: IndexedRecord[];
    try {
      records = fakeScan(bytes, fastaParser(programById(program), role));
    } catch (error) {
      return { ok: false, message: `failed to read ${role} FASTA: ${error instanceof Error ? error.message : String(error)}` };
    }
    const refused = records.find((record) => bytes.subarray(record.sequence_offset, record.end_offset).includes(0x21));
    if (refused !== undefined) {
      return {
        ok: false,
        message: `${role} record ${refused.index + 1} (${refused.id}) has '!' in its sequence (FAKE ENGINE check, not LOSAT's)`,
      };
    }
    return { ok: true, records: records.map(recordKey) };
  }
}
