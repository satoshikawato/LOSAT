// The browser's download of an Example set (ports/example.ts): the file is asked straight from
// the page (NCBI's FTP and UniProt send CORS headers; the CSP's connect-src names them,
// public/_headers), with no cookies and no referrer. NCBI sends its .gz files as
// application/x-gzip and UniProt its compressed=true answers without Content-Encoding
// (2026-11-03, losat-baselines/web-fetch-20261103/NOTES.md), so the browser leaves them
// compressed: a body that starts with the gzip magic is decompressed with DecompressionStream,
// and a body that something on the way already decoded is kept as it is.
import type { DownloadProgress, DownloadRequest, DownloadResult, FileDownloader } from '../../ports/example';

/** Milliseconds to the answer's headers; the body may take longer. */
const TIMEOUT_MS = 60_000;
/** Decompressed bytes gathered into one Blob, so that the parts held in memory stay small. */
const PART_BYTES = 8 * 1024 * 1024;
const GZIP_MAGIC = [0x1f, 0x8b] as const;

export interface DownloaderOptions {
  readonly fetch?: typeof fetch;
  readonly timeoutMs?: number;
}

export class BrowserFileDownloader implements FileDownloader {
  constructor(private readonly options: DownloaderOptions = {}) {}

  async download(
    request: DownloadRequest,
    onProgress: (progress: DownloadProgress) => void,
    signal: AbortSignal,
  ): Promise<DownloadResult> {
    signal.throwIfAborted();
    const controller = new AbortController();
    const onAbort = () => controller.abort(signal.reason);
    signal.addEventListener('abort', onAbort, { once: true });
    let timedOut = false;
    const timer = setTimeout(() => {
      timedOut = true;
      controller.abort(new DOMException('no answer in time', 'TimeoutError'));
    }, this.options.timeoutMs ?? TIMEOUT_MS);
    try {
      let response: Response;
      try {
        response = await (this.options.fetch ?? fetch)(request.url, {
          signal: controller.signal,
          credentials: 'omit',
          referrerPolicy: 'no-referrer',
        });
      } catch (error) {
        if (signal.aborted) throw error;
        return { ok: false, reason: 'unreachable', detail: timedOut ? 'no answer in time' : messageOf(error) };
      } finally {
        clearTimeout(timer);
      }
      if (!response.ok) {
        await response.body?.cancel().catch(() => undefined);
        return { ok: false, reason: 'refused', status: response.status };
      }
      const result = await readBody(response, request.maxBytes, onProgress, signal);
      if (!result.ok) return result;
      const version = versionOf(response);
      return { ...result, ...(version === undefined ? {} : { version }) };
    } finally {
      clearTimeout(timer);
      signal.removeEventListener('abort', onAbort);
    }
  }
}

type Bytes = Uint8Array<ArrayBuffer>;
type BodyResult = Exclude<DownloadResult, { readonly ok: true }> | { readonly ok: true; readonly data: Blob; readonly downloadedBytes: number };

async function readBody(
  response: Response,
  maxBytes: number,
  onProgress: (progress: DownloadProgress) => void,
  signal: AbortSignal,
): Promise<BodyResult> {
  const total = contentLength(response);
  const reader = (response.body ?? new Blob([]).stream()).getReader();
  let downloaded = 0;
  let networkError: unknown;
  const next = async (): Promise<ReadableStreamReadResult<Bytes>> => {
    let result: ReadableStreamReadResult<Bytes>;
    try {
      result = await reader.read();
    } catch (error) {
      networkError = error;
      throw error;
    }
    if (!result.done) {
      downloaded += result.value.byteLength;
      onProgress(total === undefined ? { bytes: downloaded } : { bytes: downloaded, total });
    }
    return result;
  };
  // The fetch's signal ends a browser's body too; cancelling the reader ends any body.
  const onAbort = () => void reader.cancel(signal.reason).catch(() => undefined);
  signal.addEventListener('abort', onAbort, { once: true });
  const parts: Blob[] = [];
  let chunks: Bytes[] = [];
  let chunkBytes = 0;
  let kept = 0;
  try {
    // The first bytes say whether the body is gzip.
    const head: Bytes[] = [];
    let headBytes = 0;
    let ended = false;
    while (headBytes < GZIP_MAGIC.length) {
      const { done, value } = await next();
      if (done) {
        ended = true;
        break;
      }
      head.push(value);
      headBytes += value.byteLength;
    }
    const gzip = startsWith(head, GZIP_MAGIC);
    const raw = new ReadableStream<Bytes>({
      start(stream) {
        for (const chunk of head) stream.enqueue(chunk);
        if (ended) stream.close();
      },
      async pull(stream) {
        const { done, value } = await next();
        if (done) stream.close();
        else stream.enqueue(value);
      },
      cancel: (reason) => reader.cancel(reason),
    });
    const decoded = (gzip ? raw.pipeThrough(new DecompressionStream('gzip')) : raw).getReader();
    for (;;) {
      const { done, value } = await decoded.read();
      if (done) break;
      kept += value.byteLength;
      if (kept > maxBytes) {
        await decoded.cancel().catch(() => undefined);
        return { ok: false, reason: 'too-large' };
      }
      chunks.push(value);
      chunkBytes += value.byteLength;
      if (chunkBytes >= PART_BYTES) {
        parts.push(new Blob(chunks));
        chunks = [];
        chunkBytes = 0;
      }
    }
    signal.throwIfAborted();
  } catch (error) {
    if (signal.aborted) throw signal.reason;
    if (networkError !== undefined) return { ok: false, reason: 'unreachable', detail: messageOf(networkError) };
    return { ok: false, reason: 'corrupt', detail: messageOf(error) };
  } finally {
    signal.removeEventListener('abort', onAbort);
  }
  parts.push(new Blob(chunks));
  return { ok: true, data: new Blob(parts, { type: 'text/plain' }), downloadedBytes: downloaded };
}

function startsWith(chunks: readonly Uint8Array[], magic: readonly number[]): boolean {
  const bytes: number[] = [];
  for (const chunk of chunks) {
    for (const byte of chunk) {
      if (bytes.length === magic.length) break;
      bytes.push(byte);
    }
  }
  return bytes.length === magic.length && bytes.every((byte, i) => byte === magic[i]);
}

function contentLength(response: Response): number | undefined {
  const value = Number(response.headers.get('Content-Length') ?? Number.NaN);
  return Number.isSafeInteger(value) && value >= 0 ? value : undefined;
}

/** UniProt's release (an exposed header), else the date of Last-Modified (a CORS-safelisted header). */
function versionOf(response: Response): string | undefined {
  const release = response.headers.get('X-UniProt-Release');
  if (release !== null && release.trim() !== '') return `UniProt release ${release.trim()}`;
  const modified = Date.parse(response.headers.get('Last-Modified') ?? '');
  return Number.isNaN(modified) ? undefined : new Date(modified).toISOString().slice(0, 10);
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
