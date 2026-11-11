// A source file named `.gz`, `.bz2` or `.zst` is inflated when it is opened (DW-45, line GZ), by the
// engine's own decoders through the FileInflater port (swebok-improve SW8 D2): the command line and
// the web app read such a file with the same code, so it inflates to the same bytes, or fails with
// the same sentence (`'x.fa.zst' is named .zst but is not zstd data (…); …`). The engine's name rule
// (makeblastdb's: the name ends in the suffix, in any case) decides; the inflated bytes are kept as
// a new File under the engine's plain name, so the engine never sees a plain text under a compressed
// name. Every member, stream or frame of a concatenated file is read (bgzip, pzstd, `cat`). A
// downloaded Example set (a gzip body) takes the same path in the Data worker, with its byte limit
// (WholeFileInflater, ports/inflate.ts; infra/browser/example-downloader.ts).
// A `.bz2` name is read as well (Owner 2026-11-11, swebok-improve SW8 DECISIONS 51). It was refused
// from 2026-11-04 (line GZ decision 11) because the pure-JS bzip2 decoder of that time ran at
// 1.1 MB/s and held about 3x the inflated size; the engine's decoder has neither cost.
import { INFLATED_TOO_LARGE, type FileInflater } from '../../ports/inflate';

/** Bytes read from the file at once, and inflated bytes gathered into one Blob part. */
const SLICE_BYTES = 8 * 1024 * 1024;
const PART_BYTES = 8 * 1024 * 1024;

/** The inflation gave more than its limit of bytes. */
export class InflatedTooLargeError extends Error {
  constructor(maxBytes: number) {
    super(`the file inflates to more than ${maxBytes} bytes`);
    this.name = INFLATED_TOO_LARGE;
  }
}

/** Inflated chunks gathered into Blob parts of about PART_BYTES, up to a limit of bytes. */
class PartCollector {
  private readonly parts: Blob[] = [];
  private pending: Uint8Array[] = [];
  private pendingBytes = 0;
  private bytes = 0;

  constructor(private readonly maxBytes: number) {}

  /** More than the limit came; the bytes past it are not kept. */
  get over(): boolean {
    return this.bytes > this.maxBytes;
  }

  /**
   * The chunk is copied: the engine's bytes are valid only during its call (ABI v2 §5). It is
   * called from within the engine's call, so it never throws (that would stop the instance): past
   * the limit it drops the chunk, and the caller ends the inflation after the call.
   */
  add(chunk: Uint8Array): void {
    this.bytes += chunk.length;
    if (this.over) return;
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
 * InflatedTooLargeError once more than `maxBytes` bytes have come out (a downloaded Example set; a
 * source has no limit).
 */
export async function inflateFile(file: File, inflater: FileInflater, maxBytes = Number.POSITIVE_INFINITY): Promise<File> {
  const collector = new PartCollector(maxBytes);
  const inflation = await inflater.open(file.name, (chunk) => collector.add(chunk));
  if (inflation === undefined) return file;
  try {
    for (let start = 0; start < file.size; start += SLICE_BYTES) {
      const end = Math.min(start + SLICE_BYTES, file.size);
      await inflation.push(new Uint8Array(await file.slice(start, end).arrayBuffer()));
      if (collector.over) throw new InflatedTooLargeError(maxBytes);
    }
    await inflation.end();
    if (collector.over) throw new InflatedTooLargeError(maxBytes);
  } finally {
    // An inflation that stopped early (a file that could not be read, or over its limit) is freed.
    await inflation.cancel();
  }
  return new File(collector.finish(), inflation.plainName, { type: 'text/plain' });
}
