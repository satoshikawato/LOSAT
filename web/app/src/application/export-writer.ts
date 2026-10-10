// The Writer of the exports (design §12.1): formatters hand it text and bytes in the order of
// the file, and it writes them to the download's sink in blocks of about `blockChars` bytes, so
// that no export holds its whole file in memory. A failure at any point aborts the sink, so a
// partial file is never saved.
//
// Text is encoded into the block as it comes (fix round 2), not collected and joined into one
// string of the block first: the joined strings of a large export were most of what the page's
// garbage collector had to deal with while the file was written.
import type { Downloader, ExportSink } from '../ports/download';

/** About 1 MiB per block: few blocks, and little held at once. */
export const BLOCK_CHARS = 1 << 20;

const encoder = new TextEncoder();

export class ExportWriter {
  /** The block being filled; the sink may keep none of it once its write settles. */
  private readonly block: Uint8Array;
  private filled = 0;
  private bytesWritten = 0;

  constructor(
    private readonly sink: ExportSink,
    private readonly blockChars = BLOCK_CHARS,
  ) {
    // A UTF-16 code unit takes at most 3 bytes in UTF-8, so text of blockChars units always fits.
    this.block = new Uint8Array(3 * blockChars + 4);
  }

  /** Bytes written to the sink so far. */
  get written(): number {
    return this.bytesWritten;
  }

  /** Appends text (UTF-8 in the file). */
  async text(text: string): Promise<void> {
    let rest = this.put(text);
    while (rest !== '') {
      await this.flush();
      rest = this.put(rest);
    }
    if (this.filled >= this.blockChars) await this.flush();
  }

  /** Appends texts, in order, as `text` would one after another (without an await for each). */
  async texts(texts: readonly string[]): Promise<void> {
    for (const text of texts) {
      let rest = this.put(text);
      while (rest !== '') {
        await this.flush();
        rest = this.put(rest);
      }
      if (this.filled >= this.blockChars) await this.flush();
    }
  }

  /** Appends bytes as they are, after the text before them. */
  async bytes(bytes: Uint8Array): Promise<void> {
    await this.flush();
    if (bytes.length === 0) return;
    await this.sink.write(bytes);
    this.bytesWritten += bytes.length;
  }

  /** Writes the text collected so far as one block. */
  async flush(): Promise<void> {
    if (this.filled === 0) return;
    await this.sink.write(this.block.subarray(0, this.filled));
    this.bytesWritten += this.filled;
    this.filled = 0;
  }

  /** Encodes as much of `text` as the block has room for; returns the rest. */
  private put(text: string): string {
    const { read, written } = encoder.encodeInto(text, this.block.subarray(this.filled));
    this.filled += written;
    return read < text.length ? text.slice(read) : '';
  }
}

/**
 * Writes one file through `write` and saves it. A failure of `write` (or of the sink) aborts
 * the file, so nothing is saved, and is thrown again. Resolves with the file's length in bytes.
 */
export async function writeFile(
  downloader: Pick<Downloader, 'open'>,
  fileName: string,
  mimeType: string,
  write: (writer: ExportWriter) => Promise<void>,
): Promise<number> {
  const sink = downloader.open(fileName, mimeType);
  try {
    const writer = new ExportWriter(sink);
    await write(writer);
    await writer.flush();
    await sink.close();
    return writer.written;
  } catch (error) {
    sink.abort();
    throw error;
  }
}
