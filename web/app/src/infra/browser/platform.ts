// Browser implementations of small platform services used by the composition root.
import type { Downloader, ExportSink } from '../../ports/download';

/**
 * Blocks are joined into one Blob after this many bytes: the browser then holds them in its
 * own Blob storage (which can page large Blobs out of memory), and the copies on the script's
 * heap are freed. A file is never assembled as one buffer (design §12.1).
 */
const BLOB_JOIN_BYTES = 16 * 1024 * 1024;

/** Hands a Blob to the browser's download. */
function download(fileName: string, blob: Blob): void {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  link.click();
  // Revoking immediately can cancel the download in some browsers.
  setTimeout(() => URL.revokeObjectURL(url), 60_000);
}

/** Lets the page paint between blocks of a long export. */
const nextTask = () => new Promise<void>((resolve) => setTimeout(resolve, 0));

class BlobSink implements ExportSink {
  private parts: BlobPart[] = [];
  private unjoined = 0;
  private state: 'open' | 'closed' | 'aborted' = 'open';

  constructor(
    private readonly fileName: string,
    private readonly mimeType: string,
  ) {}

  async write(bytes: Uint8Array): Promise<void> {
    if (this.state !== 'open') throw new Error(`the file ${this.fileName} is no longer open`);
    this.parts.push(bytes.slice());
    this.unjoined += bytes.length;
    if (this.unjoined >= BLOB_JOIN_BYTES) {
      this.parts = [new Blob(this.parts)];
      this.unjoined = 0;
    }
    await nextTask();
  }

  async close(): Promise<void> {
    if (this.state !== 'open') throw new Error(`the file ${this.fileName} is no longer open`);
    this.state = 'closed';
    const blob = new Blob(this.parts, { type: this.mimeType });
    this.parts = [];
    download(this.fileName, blob);
  }

  abort(): void {
    if (this.state === 'open') this.state = 'aborted';
    this.parts = [];
  }
}

export const browserDownloader: Downloader = {
  open(fileName, mimeType) {
    return new BlobSink(fileName, mimeType);
  },
  save(fileName, bytes, mimeType) {
    download(fileName, new Blob([bytes as BlobPart], { type: mimeType }));
  },
};

export async function sha256Hex(bytes: Uint8Array): Promise<string> {
  const digest = await crypto.subtle.digest('SHA-256', bytes as BufferSource);
  return Array.from(new Uint8Array(digest), (b) => b.toString(16).padStart(2, '0')).join('');
}
