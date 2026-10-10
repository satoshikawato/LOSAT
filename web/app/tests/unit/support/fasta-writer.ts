// Writes FASTA text for the extraction tests and keeps, for each record it writes, where each of
// its residues is: the test's own reference for docs/web/abi_v2.md §9 (reader kinds 1 and 2),
// computed from the bytes as they are written, not by the app's code. Every byte is ASCII, so a
// string's length is its length in bytes.
import type { IndexedRecord, LineLayout } from '../../../src/domain/dataset';

export type ReaderKind = 1 | 2;

/** A seeded random source (mulberry32), so that a failing case can be run again. */
export function seeded(seed: number): () => number {
  let state = seed >>> 0;
  return () => {
    state = (state + 0x6d2b79f5) >>> 0;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export const pick = <T>(random: () => number, items: readonly T[]): T => items[Math.floor(random() * items.length)]!;
export const between = (random: () => number, low: number, high: number): number => low + Math.floor(random() * (high - low + 1));

/** Random residues from `alphabet`, each in lower case with probability `lower`. */
export function residues(random: () => number, length: number, alphabet: string, lower = 0): string {
  let out = '';
  for (let i = 0; i < length; i++) {
    const letter = alphabet[Math.floor(random() * alphabet.length)]!;
    out += random() < lower ? letter.toLowerCase() : letter;
  }
  return out;
}

/** How a record's residues are laid out in lines. */
export type Lines =
  /** Every line `width` residues and then the bytes `eol` (spaces or tabs, then LF or CR LF); the last line may be shorter. */
  | { readonly kind: 'uniform'; readonly width: number; readonly eol: string }
  /**
   * Lines of `minWidth` to `maxWidth` residues ended by `eol`. With `noise`, white space at line
   * starts, skipped bytes between residues (digits, `-`, `.`, `?`, `_`, `!`, `#`, and letters
   * the kind does not store), `;` comments at line ends, and comment and blank lines between them.
   */
  | { readonly kind: 'ragged'; readonly minWidth: number; readonly maxWidth: number; readonly eol: '\n' | '\r\n'; readonly noise?: boolean; readonly comments?: boolean };

export interface WrittenRecord {
  readonly title: string;
  readonly id: string;
  /** The residues as written. */
  readonly letters: string;
  readonly headerOffset: number;
  readonly sequenceOffset: number;
  /** Set when the next record starts or the text ends. */
  endOffset: number;
  /** The byte offset of each residue. */
  readonly offsets: number[];
  readonly uniform?: { readonly width: number; readonly eol: number };
}

const NOT_STORED: Readonly<Record<ReaderKind, string>> = { 1: 'EFIJLOPQXZefijlopqxz*', 2: '' };
const SKIPPED = '0123456789-.?_!#';

export class FastaWriter {
  private text = '';
  readonly records: WrittenRecord[] = [];

  constructor(
    private readonly random: () => number,
    private readonly kind: ReaderKind,
  ) {}

  /** Text that is not part of a record's residues (before the first record: comments, blank lines). */
  raw(text: string): this {
    this.text += text;
    return this;
  }

  record(title: string, letters: string, lines: Lines, headerEnd: '\n' | '\r\n' = lines.eol.endsWith('\r\n') ? '\r\n' : '\n'): WrittenRecord {
    this.close();
    const headerOffset = this.text.length;
    this.text += `>${title}${headerEnd}`;
    const sequenceOffset = this.text.length;
    const offsets: number[] = [];
    const put = (letter: string) => {
      offsets.push(this.text.length);
      this.text += letter;
    };
    if (lines.kind === 'uniform') {
      for (let i = 0; i < letters.length; i += lines.width) {
        for (const letter of letters.slice(i, i + lines.width)) put(letter);
        this.text += lines.eol;
      }
    } else {
      let i = 0;
      while (i < letters.length) {
        if (lines.comments === true || lines.noise === true) this.betweenLines(lines.eol);
        const width = between(this.random, lines.minWidth, lines.maxWidth);
        if (lines.noise === true && this.random() < 0.3) this.text += pick(this.random, [' ', '\t', '\v', '\f', '  ']);
        for (const letter of letters.slice(i, i + width)) {
          put(letter);
          if (lines.noise === true && this.random() < 0.05) this.text += this.skippedByte();
        }
        i += width;
        if (lines.noise === true && this.random() < 0.2) this.text += pick(this.random, [' ', '  ', '\t']);
        if (lines.noise === true && this.random() < 0.1) this.text += `;${residues(this.random, 5, 'ACGTMK')} note`;
        this.text += lines.eol;
      }
    }
    const record: WrittenRecord = {
      title,
      id: title.split(' ')[0]!,
      letters,
      headerOffset,
      sequenceOffset,
      endOffset: -1,
      offsets,
      ...(lines.kind === 'uniform' ? { uniform: uniformLayout(letters.length, lines.width, lines.eol.length) } : {}),
    };
    this.records.push(record);
    return record;
  }

  bytes(): Uint8Array {
    this.close();
    return new TextEncoder().encode(this.text);
  }

  toString(): string {
    this.close();
    return this.text;
  }

  /** A comment line, a blank line or nothing, before a line of residues. */
  private betweenLines(eol: string): void {
    const roll = this.random();
    if (roll < 0.04) this.text += eol;
    else if (roll < 0.06) this.text += ` \t${eol}`;
    else if (roll < 0.1) this.text += `${pick(this.random, ['', ' ', '\t'])}${pick(this.random, [';', '!', '#'])} ACGT ${residues(this.random, 4, 'ACDEFGHIKLMNPQRSTVWY')}${eol}`;
  }

  private skippedByte(): string {
    return pick(this.random, [...SKIPPED, ...NOT_STORED[this.kind]]);
  }

  private close(): void {
    const last = this.records[this.records.length - 1];
    if (last !== undefined && last.endOffset < 0) last.endOffset = this.text.length;
  }
}

/** The uniform layout of `length` residues in lines of `width` and `eol` bytes, as the scan states it (abi_v2.md §9). */
function uniformLayout(length: number, width: number, eol: number): { width: number; eol: number } {
  if (length === 0) return { width: 0, eol: 1 };
  return length <= width ? { width: length, eol: 1 } : { width, eol };
}

/** What the scan of kind `kind` reports for a written record, with checkpoints `every` residues unless it is uniform. */
export function indexedRecord(record: WrittenRecord, index: number, kind: ReaderKind, every = 65_536): IndexedRecord {
  if (record.endOffset < 0) throw new Error('the record is not closed: call bytes() first');
  const line_layout: LineLayout =
    record.uniform !== undefined
      ? { kind: 'uniform', ...record.uniform }
      : { kind: 'checkpoints', every, offsets: record.offsets.filter((_, i) => i % every === 0) };
  return {
    index,
    id: record.id,
    header_offset: record.headerOffset,
    sequence_offset: record.sequenceOffset,
    end_offset: record.endOffset,
    length: record.letters.length,
    line_layout,
    residue_counts: counts(record.letters, kind),
  };
}

/** Upper-cased counts, `U` counted as `T` in kind 1 (abi_v2.md §9). */
export function counts(letters: string, kind: ReaderKind): Record<string, number> {
  const out: Record<string, number> = {};
  for (const letter of letters) {
    let key = letter.toUpperCase();
    if (kind === 1 && key === 'U') key = 'T';
    out[key] = (out[key] ?? 0) + 1;
  }
  return out;
}
