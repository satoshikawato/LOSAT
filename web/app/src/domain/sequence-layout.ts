// Where the letters of a record are in its source file (docs/web/abi_v2.md §9 `line_layout`):
// the byte of residue `i` of a uniform record from the formula, the checkpoint to start from in
// any other record, and the residues read forward from there under the rules of the reader kind
// that scanned the source (1: NCBI's reader with nucleotide flags, 2: with protein flags). The
// Data worker reads the File's bytes with these (infra/data/data-service.ts), so an extracted
// sequence is the file's own letters (design §6.1, §11.4 "大文字小文字を保つ"): the scan's stored
// value is upper-cased and reads `U` as `T`, but nothing here changes a byte.
import type { LineLayout } from './dataset';

/** The reader kinds whose records extraction reads (abi_v2.md §9; kind 0 is not used by the app). */
export type ReaderKind = 1 | 2;

/** The reader kind of a revision's parser; any other kind is refused. */
export function readerKind(parser: number): ReaderKind {
  if (parser !== 1 && parser !== 2) {
    throw new RangeError(`records indexed with FASTA reader kind ${parser} cannot be extracted: only kinds 1 and 2 (NCBI's reader) can`);
  }
  return parser;
}

/** The byte offset of residue `i` (0-based) of a uniform record (any positive `eol`). */
export function uniformOffset(layout: { readonly width: number; readonly eol: number }, sequenceOffset: number, i: number): number {
  if (layout.width <= 0) throw new RangeError('a record without residues has no residue offset');
  return sequenceOffset + Math.floor(i / layout.width) * (layout.width + layout.eol) + (i % layout.width);
}

/**
 * Where to start reading residue `first` (0-based) of a record: the byte of a residue, and how
 * many residues to pass before `first`. A uniform record starts at `first` itself; any other at
 * the last checkpoint at or before it.
 */
export function readStart(layout: LineLayout, sequenceOffset: number, first: number): { readonly offset: number; readonly skip: number } {
  if (layout.kind === 'uniform') return { offset: uniformOffset(layout, sequenceOffset, first), skip: 0 };
  const k = Math.floor(first / layout.every);
  const offset = layout.offsets[k];
  if (offset === undefined) throw new RangeError(`the record has no checkpoint for residue ${first + 1}`);
  return { offset, skip: first - k * layout.every };
}

const LF = 0x0a;
const CR = 0x0d;
const SEMICOLON = 0x3b;

/** Kind 1 stores these letters in either case (abi_v2.md §9). */
const NUCLEOTIDE_LETTERS = 'ABCDGHKMNRSTUVWY';
const STORED: Readonly<Record<ReaderKind, Uint8Array>> = {
  1: table((byte) => byte < 0x80 && NUCLEOTIDE_LETTERS.includes(String.fromCharCode(byte).toUpperCase())),
  2: table((byte) => (byte >= 0x41 && byte <= 0x5a) || (byte >= 0x61 && byte <= 0x7a) || byte === 0x2a),
};
/** Space, tab, VT and FF: skipped at the start of a line. */
const LINE_START_SKIP = table((byte) => byte === 0x20 || byte === 0x09 || byte === 0x0b || byte === 0x0c);
/** `!`, `#` and `;` as the first other byte of a line make it a comment line. */
const COMMENT_START = table((byte) => byte === 0x21 || byte === 0x23 || byte === SEMICOLON);

function table(test: (byte: number) => boolean): Uint8Array {
  return Uint8Array.from({ length: 256 }, (_, byte) => (test(byte) ? 1 : 0));
}

/** Whether reader `kind` stores `byte` as a residue (in either case). */
export const stores = (kind: ReaderKind, byte: number): boolean => STORED[kind][byte] === 1;

/**
 * Reads the residues of a record forward from a checkpoint (abi_v2.md §9 "Reading forward from a
 * checkpoint"), fed chunk by chunk, and keeps `count` of them after passing `skip`. The chunks are
 * the source's bytes from the checkpoint's offset on, in order, up to the record's `end_offset`;
 * any chunk sizes give the same residues. CR and LF end a line; at the start of a line, space,
 * tab, VT and FF are skipped, and a line whose first other byte is `!`, `#` or `;` is skipped
 * whole; elsewhere `;` skips the rest of the line; a byte that the kind stores is the next
 * residue, kept as the file has it; every other byte is skipped.
 */
export class ForwardReader {
  private readonly stored: Uint8Array;
  private readonly kept: Uint8Array;
  private keptCount = 0;
  /** Residues passed before `skip`. */
  private passed = 0;
  /** In the bytes that start a line (before its first other byte). A checkpoint's byte is a residue. */
  private lineStart = false;
  /** The rest of this line is skipped. */
  private skippingLine = false;

  constructor(
    kind: ReaderKind,
    private readonly skip: number,
    private readonly count: number,
  ) {
    if (!Number.isSafeInteger(skip) || skip < 0 || !Number.isSafeInteger(count) || count < 0) {
      throw new RangeError(`cannot read ${count} residues after ${skip}`);
    }
    this.stored = STORED[kind];
    this.kept = new Uint8Array(count);
  }

  /** Whether `count` residues have been read. */
  get done(): boolean {
    return this.keptCount >= this.count;
  }

  /** Reads the next bytes; returns true while more residues are wanted. */
  feed(bytes: Uint8Array): boolean {
    const { stored, kept, count } = this;
    let { keptCount, passed, lineStart, skippingLine } = this;
    for (let i = 0; i < bytes.length && keptCount < count; i++) {
      const byte = bytes[i]!;
      if (byte === LF || byte === CR) {
        lineStart = true;
        skippingLine = false;
        continue;
      }
      if (skippingLine) continue;
      if (lineStart) {
        if (LINE_START_SKIP[byte] === 1) continue;
        lineStart = false;
        if (COMMENT_START[byte] === 1) {
          skippingLine = true;
          continue;
        }
      }
      if (byte === SEMICOLON) {
        skippingLine = true;
      } else if (stored[byte] === 1) {
        if (passed < this.skip) passed++;
        else kept[keptCount++] = byte;
      }
    }
    this.keptCount = keptCount;
    this.passed = passed;
    this.lineStart = lineStart;
    this.skippingLine = skippingLine;
    return keptCount < count;
  }

  /** The residues kept: `count` of them once the reader is done, fewer if the bytes ran out first. */
  residues(): Uint8Array {
    return this.keptCount === this.kept.length ? this.kept : this.kept.subarray(0, this.keptCount);
  }
}

/**
 * The counts of residues as the scan's `residue_counts` keeps them (abi_v2.md §9): upper-cased,
 * and `U` counted as `T` in kind 1. Every byte given must be a residue of the kind.
 */
export function residueCounts(kind: ReaderKind, residues: Uint8Array): Record<string, number> {
  const counts = new Float64Array(256);
  for (const byte of residues) counts[byte]!++;
  const out: Record<string, number> = {};
  for (let byte = 0; byte < 256; byte++) {
    const n = counts[byte]!;
    if (n === 0) continue;
    let key = String.fromCharCode(byte).toUpperCase();
    if (kind === 1 && key === 'U') key = 'T';
    out[key] = (out[key] ?? 0) + n;
  }
  return out;
}
