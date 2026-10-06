// A fake of the ABI v2 index scan (docs/web/abi_v2.md §9) for tests and for the
// development build until S09 connects the adapter's scan. It is not the engine's FASTA
// reader and must not be used for anything else. It follows the parser kind 0 rules of
// the ABI document (the `bio::io::fasta` reader) closely enough to pass
// tests/contract/record-scanner.contract.ts, with these simplifications: it validates
// UTF-8 line by line (only the lines that bio reads), and it always reports the
// checkpoints layout.
import { recordKey, type IndexedRecord, type RecordKey } from '../../domain/dataset';
import type { InputRole, ProgramId } from '../../domain/programs';
import type { InputCheck, InputChecker } from '../../ports/input-check';
import type { RecordScanner, ScanResponse } from '../../ports/scan';
import { concatBytes } from '../bytes';

const GT = 0x3e;
const LF = 0x0a;
const CHECKPOINT_EVERY = 65_536;
// Rust's char::is_whitespace (the White_Space property), which bio's trim_end uses.
const WHITESPACE = /[\t\n\v\f\r \u0085\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]/u;
const TRAILING_WHITESPACE = /[\t\n\v\f\r \u0085\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]+$/u;

export class FakeScanner implements RecordScanner {
  async scan(parser: number, chunks: AsyncIterable<Uint8Array>): Promise<ScanResponse> {
    if (parser !== 0) throw new Error(`unknown or unavailable FASTA parser kind ${parser}`);
    const parts: Uint8Array[] = [];
    for await (const chunk of chunks) parts.push(chunk.slice());
    return { records: fakeScan(concatBytes(parts)) };
  }
}

/** The records that the fake reads, as the FakeEngine's `register` reports them. */
export function fakeRecordKeys(bytes: Uint8Array): RecordKey[] {
  return fakeScan(bytes).map(recordKey);
}

interface Builder {
  readonly index: number;
  readonly id: string;
  readonly hasDescription: boolean;
  readonly headerOffset: number;
  readonly sequenceOffset: number;
  length: number;
  readonly counts: number[];
  readonly checkpoints: number[];
}

/** Decodes one line as bio reads it: UTF-8, a leading U+FEFF kept. */
function decodeLine(line: Uint8Array): string {
  try {
    return new TextDecoder('utf-8', { fatal: true, ignoreBOM: true }).decode(line);
  } catch {
    throw new Error('stream did not contain valid UTF-8');
  }
}

export function fakeScan(bytes: Uint8Array): IndexedRecord[] {
  const records: IndexedRecord[] = [];
  let current: Builder | undefined;
  let offset = 0;
  while (offset < bytes.length) {
    const newline = bytes.indexOf(LF, offset);
    const lineEnd = newline < 0 ? bytes.length : newline;
    const next = newline < 0 ? bytes.length : newline + 1;
    decodeLine(bytes.subarray(offset, lineEnd));
    if (offset === 0 && bytes[0] !== GT) throw new Error('Expected > at record start.');
    if (bytes[offset] === GT) {
      if (current !== undefined) {
        const record = finish(current, offset);
        // bio's records() ends at the first empty record.
        if (record === undefined) return records;
        records.push(record);
      }
      current = startRecord(bytes, offset, lineEnd, next, records.length);
    } else if (current !== undefined) {
      addLine(current, bytes.subarray(offset, lineEnd), offset);
    }
    offset = next;
  }
  const last = current === undefined ? undefined : finish(current, bytes.length);
  if (last !== undefined) records.push(last);
  return records;
}

function startRecord(bytes: Uint8Array, start: number, lineEnd: number, next: number, index: number): Builder {
  const header = decodeLine(bytes.subarray(start + 1, lineEnd)).replace(TRAILING_WHITESPACE, '');
  const split = header.search(WHITESPACE);
  return {
    index,
    id: split < 0 ? header : header.slice(0, split),
    hasDescription: split >= 0,
    headerOffset: start,
    sequenceOffset: next,
    length: 0,
    counts: new Array<number>(256).fill(0),
    checkpoints: [],
  };
}

/** A sequence line: every byte except the trailing white space is a residue. */
function addLine(builder: Builder, line: Uint8Array, lineOffset: number): void {
  let keep = line.length;
  if (line.every((byte) => byte < 0x80)) {
    while (keep > 0 && isAsciiWhitespace(line[keep - 1]!)) keep--;
  } else {
    keep = new TextEncoder().encode(decodeLine(line).replace(TRAILING_WHITESPACE, '')).length;
  }
  for (let i = 0; i < keep; i++) {
    if ((builder.length + i) % CHECKPOINT_EVERY === 0) builder.checkpoints.push(lineOffset + i);
    builder.counts[line[i]!]!++;
  }
  builder.length += keep;
}

function finish(builder: Builder, end: number): IndexedRecord | undefined {
  if (builder.id === '' && !builder.hasDescription && builder.length === 0) return undefined;
  const residueCounts: Record<string, number> = {};
  builder.counts.forEach((count, byte) => {
    if (count === 0) return;
    const key = byte >= 0x21 && byte <= 0x7e ? String.fromCharCode(byte) : `0x${byte.toString(16).padStart(2, '0')}`;
    residueCounts[key] = count;
  });
  return {
    index: builder.index,
    id: builder.id,
    header_offset: builder.headerOffset,
    sequence_offset: builder.sequenceOffset,
    end_offset: end,
    length: builder.length,
    line_layout: { kind: 'checkpoints', every: CHECKPOINT_EVERY, offsets: builder.checkpoints },
    residue_counts: residueCounts,
  };
}

function isAsciiWhitespace(byte: number): boolean {
  return byte === 0x20 || (byte >= 0x09 && byte <= 0x0d);
}


/**
 * The input check of the development build (ports/input-check.ts). It is not the engine's:
 * it accepts what the FakeScanner reads, except a record with `!` in its sequence, which
 * it refuses with a message shaped as the engine's (the record's number and ID), so that
 * the screen's handling of a refused record can be tried without the engine.
 */
export class FakeInputChecker implements InputChecker {
  async check(program: ProgramId, role: InputRole, bytes: Uint8Array): Promise<InputCheck> {
    if (program === 'blastx') return { ok: false, message: 'blastx is not available in LOSAT Web ABI v2 yet' };
    let records: IndexedRecord[];
    try {
      records = fakeScan(bytes);
    } catch (error) {
      return { ok: false, message: `failed to read ${role} FASTA: ${error instanceof Error ? error.message : String(error)}` };
    }
    const refused = records.find((record) => (record.residue_counts['!'] ?? 0) > 0);
    if (refused !== undefined) {
      return {
        ok: false,
        message: `${role} record ${refused.index + 1} (${refused.id}) has '!' in its sequence (FAKE ENGINE check, not LOSAT's)`,
      };
    }
    return { ok: true, records: records.map(recordKey) };
  }
}
