// A source file named `.gz` or `.zst` is inflated when it is opened (DW-45, line GZ), by the
// engine's own decoders through the FileInflater port (swebok-improve SW8 D2): the command line and
// the web app read such a file with the same code, so it inflates to the same bytes, or fails with
// the same sentence (`'x.fa.zst' is named .zst but is not zstd data (…); …`). The engine's name rule
// (makeblastdb's: the name ends in the suffix, in any case) decides; the inflated bytes are kept as
// a new File under the engine's plain name, so the engine never sees a plain text under a compressed
// name. Every member, stream or frame of a concatenated file is read (bgzip, pzstd, `cat`).
// A `.bz2` name is not read by the web app (decision 11, 2026-11-04): the pure-JS bzip2 decoder of
// that time ran at 1.1 MB/s and held about 3x the inflated size; the command line reads it.
import type { FileInflater } from '../../ports/inflate';

/** Bytes read from the file at once, and inflated bytes gathered into one Blob part. */
const SLICE_BYTES = 8 * 1024 * 1024;
const PART_BYTES = 8 * 1024 * 1024;

/** Inflated chunks gathered into Blob parts of about PART_BYTES. */
class PartCollector {
  private readonly parts: Blob[] = [];
  private pending: Uint8Array[] = [];
  private pendingBytes = 0;

  /** The chunk is copied: the engine's bytes are valid only during its call (ABI v2 §5). */
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

/**
 * The file inflated under its plain name, or the file itself when its name is plain. Throws the
 * engine's sentence when the data is not of the named compression (or is cut short, or empty), and
 * the web app's refusal of a `.bz2` name.
 */
export async function inflateFile(file: File, inflater: FileInflater): Promise<File> {
  const collector = new PartCollector();
  const inflation = await inflater.open(file.name, (chunk) => collector.add(chunk));
  if (inflation === undefined) return file;
  try {
    if (inflation.compression === 'bzip2') {
      throw new Error(`'${file.name}' is named .bz2; the web app reads .gz and .zst, not bzip2 (the command line does)`);
    }
    for (let start = 0; start < file.size; start += SLICE_BYTES) {
      const end = Math.min(start + SLICE_BYTES, file.size);
      await inflation.push(new Uint8Array(await file.slice(start, end).arrayBuffer()));
    }
    await inflation.end();
  } finally {
    // An inflation that stopped early (a refused name, a file that could not be read) is freed.
    await inflation.cancel();
  }
  return new File(collector.finish(), inflation.plainName, { type: 'text/plain' });
}
