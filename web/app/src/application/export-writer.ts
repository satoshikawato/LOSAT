// The Writer of the exports (design §12.1): formatters hand it text and bytes in the order of
// the file, and it writes them to the download's sink in blocks of about `blockChars`
// characters, so that no export holds its whole file in memory. A failure at any point aborts
// the sink, so a partial file is never saved.
import type { Downloader, ExportSink } from '../ports/download';

/** About 1 MiB of text per block: few blocks, and little held at once. */
export const BLOCK_CHARS = 1 << 20;

const encoder = new TextEncoder();

export class ExportWriter {
  private pending: string[] = [];
  private pendingChars = 0;
  private bytesWritten = 0;

  constructor(
    private readonly sink: ExportSink,
    private readonly blockChars = BLOCK_CHARS,
  ) {}

  /** Bytes written to the sink so far. */
  get written(): number {
    return this.bytesWritten;
  }

  /** Appends text (UTF-8 in the file). */
  async text(text: string): Promise<void> {
    if (text === '') return;
    this.pending.push(text);
    this.pendingChars += text.length;
    if (this.pendingChars >= this.blockChars) await this.flush();
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
    if (this.pending.length === 0) return;
    const block = encoder.encode(this.pending.join(''));
    this.pending = [];
    this.pendingChars = 0;
    await this.sink.write(block);
    this.bytesWritten += block.length;
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
