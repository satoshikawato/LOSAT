// A source file named `.gz`, `.bz2` or `.zst` is inflated when it is opened (DW-45, line GZ). The
// rule is the one of NCBI's makeblastdb: the name ends in the suffix, in any case. The inflated
// bytes are kept as a new File under the name without the suffix, so the engine never sees a plain
// text under a compressed name.
// gzip: fflate's streaming Gunzip reads every member of a concatenated file (bgzip).
// `DecompressionStream('gzip')` rejects concatenated members in Chromium ("Failed to fetch"),
// Firefox ("The operation was aborted") and WebKit ("Extra bytes past the end"), checked 2026-11-03.
// zstd: fzstd's streaming Decompress reads every frame and passes over skippable frames (pzstd).
// bzip2: the `bz2` package inflates one stream per call and has no streaming API, so this path
// holds the whole compressed file in memory and splits it at its stream heads (pbzip2 and
// `cat a.bz2 b.bz2` write several streams).
import { Gunzip } from 'fflate';
import { Decompress } from 'fzstd';
// `bz2` is a CommonJS module that assigns `window.bz2` instead of `module.exports` when a global
// `window` exists, so it can be imported only in a Worker (or Node). Never import this module from
// code that runs on the page thread.
import bz2 from 'bz2';

export type Compression = 'gzip' | 'bzip2' | 'zstd';

/** Bytes read from the file at once, and inflated bytes gathered into one Blob part. */
const SLICE_BYTES = 8 * 1024 * 1024;
const PART_BYTES = 8 * 1024 * 1024;
const SUFFIXES: readonly (readonly [RegExp, Compression, string])[] = [
  [/\.gz$/i, 'gzip', '.gz'],
  [/\.bz2$/i, 'bzip2', '.bz2'],
  [/\.zst$/i, 'zstd', '.zst'],
];

/** The compression that the name ends in, in any case; undefined for any other name. */
export function compressionOfName(name: string): Compression | undefined {
  return SUFFIXES.find(([suffix]) => suffix.test(name))?.[1];
}

/** The name without `.gz`, `.bz2` or `.zst`; `input` when nothing is left. */
export function plainName(name: string): string {
  const found = SUFFIXES.find(([suffix]) => suffix.test(name));
  const plain = found === undefined ? name : name.replace(found[0], '');
  return plain === '' ? 'input' : plain;
}

/** Inflated chunks gathered into Blob parts of about PART_BYTES. */
class PartCollector {
  private readonly parts: Blob[] = [];
  private pending: Uint8Array[] = [];
  private pendingBytes = 0;

  /** The chunk is copied, as the decompressors may reuse their buffer later. */
  add(chunk: Uint8Array): void {
    this.pending.push(chunk.slice());
    this.pendingBytes += chunk.length;
    if (this.pendingBytes >= PART_BYTES) this.flush();
  }

  private flush(): void {
    if (this.pendingBytes === 0) return;
    this.parts.push(new Blob(this.pending as BlobPart[]));
    this.pending = [];
    this.pendingBytes = 0;
  }

  finish(): Blob[] {
    this.flush();
    return this.parts;
  }
}

async function pushSlices(file: File, push: (bytes: Uint8Array, last: boolean) => void): Promise<void> {
  for (let start = 0; start < file.size; start += SLICE_BYTES) {
    const end = Math.min(start + SLICE_BYTES, file.size);
    push(new Uint8Array(await file.slice(start, end).arrayBuffer()), end === file.size);
  }
}

const BLOCK_MAGIC = [0x31, 0x41, 0x59, 0x26, 0x53, 0x59];
const END_MAGIC = [0x17, 0x72, 0x45, 0x38, 0x50, 0x90];

/** Whether a bzip2 stream head (`BZh1`-`BZh9` and a block or end-of-stream magic) starts at `at`. */
function isBzip2Head(bytes: Uint8Array, at: number): boolean {
  if (bytes[at] !== 0x42 || bytes[at + 1] !== 0x5a || bytes[at + 2] !== 0x68) return false;
  const level = bytes[at + 3];
  if (level === undefined || level < 0x31 || level > 0x39) return false;
  const matches = (magic: number[]) => magic.every((b, i) => bytes[at + 4 + i] === b);
  return matches(BLOCK_MAGIC) || matches(END_MAGIC);
}

/** Start offsets of the bzip2 streams in `bytes`; the first stream must start at 0. */
function bzip2Streams(bytes: Uint8Array): number[] {
  if (!isBzip2Head(bytes, 0)) throw new Error('no bzip2 stream head at the start');
  const starts = [0];
  // A head inside a stream's bit-packed body is possible in principle (the blocks are not byte
  // aligned); 10 matching bytes by chance are not, and a false split fails the CRC of its piece.
  for (let at = bytes.indexOf(0x42, 1); at !== -1; at = bytes.indexOf(0x42, at + 1)) {
    if (isBzip2Head(bytes, at)) starts.push(at);
  }
  return starts;
}

async function bunzip2(file: File): Promise<Blob[]> {
  const bytes = new Uint8Array(await file.arrayBuffer());
  const starts = bzip2Streams(bytes);
  const collector = new PartCollector();
  starts.forEach((start, i) => {
    collector.add(bz2.decompress(bytes.subarray(start, starts[i + 1] ?? bytes.length), true));
  });
  return collector.finish();
}

async function inflateParts(file: File, compression: Compression): Promise<Blob[]> {
  if (compression === 'bzip2') return bunzip2(file);
  const collector = new PartCollector();
  const decoder =
    compression === 'gzip'
      ? new Gunzip((chunk) => collector.add(chunk))
      : new Decompress((chunk) => collector.add(chunk));
  await pushSlices(file, (bytes, last) => decoder.push(bytes, last));
  return collector.finish();
}

const SUFFIX_OF: Record<Compression, string> = { gzip: '.gz', bzip2: '.bz2', zstd: '.zst' };

/** The file inflated under its plain name; throws when the data is not of the named compression. */
export async function inflateFile(file: File): Promise<File> {
  const compression = compressionOfName(file.name);
  if (compression === undefined) return file;
  let parts: Blob[];
  try {
    if (file.size === 0) throw new Error('the file is empty');
    parts = await inflateParts(file, compression);
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    throw new Error(`${file.name} is named ${SUFFIX_OF[compression]} but is not ${compression} data (${message})`);
  }
  return new File(parts, plainName(file.name), { type: 'text/plain' });
}
