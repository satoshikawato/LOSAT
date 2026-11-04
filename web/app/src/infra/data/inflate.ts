// A source file named `.gz` or `.zst` is inflated when it is opened (DW-45, line GZ). The rule is
// the one of NCBI's makeblastdb: the name ends in the suffix, in any case. The inflated bytes are
// kept as a new File under the name without the suffix, so the engine never sees a plain text under
// a compressed name. A `.bz2` name is not read by the web app (decision 11, 2026-11-04): the
// pure-JS bzip2 decoder ran at 1.1 MB/s and held about 3x the inflated size; the command line reads it.
// gzip: fflate's streaming Gunzip reads every member of a concatenated file (bgzip).
// `DecompressionStream('gzip')` rejects concatenated members in Chromium ("Failed to fetch"),
// Firefox ("The operation was aborted") and WebKit ("Extra bytes past the end"), checked 2026-11-03.
// zstd: fzstd's streaming Decompress reads every frame and passes over skippable frames (pzstd).
import { Gunzip } from 'fflate';
import { Decompress } from 'fzstd';

export type Compression = 'gzip' | 'zstd';

/** Bytes read from the file at once, and inflated bytes gathered into one Blob part. */
const SLICE_BYTES = 8 * 1024 * 1024;
const PART_BYTES = 8 * 1024 * 1024;
const SUFFIXES: readonly (readonly [RegExp, Compression, string])[] = [
  [/\.gz$/i, 'gzip', '.gz'],
  [/\.zst$/i, 'zstd', '.zst'],
];

/** The compression that the name ends in, in any case; undefined for any other name. */
export function compressionOfName(name: string): Compression | undefined {
  return SUFFIXES.find(([suffix]) => suffix.test(name))?.[1];
}

/** `.bz2` when the name ends in it, in any case: the web app refuses such a file; undefined otherwise. */
export function refusedCompressionOfName(name: string): string | undefined {
  return /\.bz2$/i.test(name) ? '.bz2' : undefined;
}

/** The name without `.gz` or `.zst`; `input` when nothing is left. */
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

async function inflateParts(file: File, compression: Compression): Promise<Blob[]> {
  const collector = new PartCollector();
  const decoder =
    compression === 'gzip'
      ? new Gunzip((chunk) => collector.add(chunk))
      : new Decompress((chunk) => collector.add(chunk));
  await pushSlices(file, (bytes, last) => decoder.push(bytes, last));
  return collector.finish();
}

const SUFFIX_OF: Record<Compression, string> = { gzip: '.gz', zstd: '.zst' };

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
