// The browser's download of an Example set (ports/example.ts): the file is asked straight from
// the page (NCBI's FTP and UniProt send CORS headers; the CSP's connect-src names them,
// public/_headers), with no cookies and no referrer. NCBI sends its .gz files as
// application/x-gzip and UniProt its compressed=true answers without Content-Encoding
// (2026-11-03, losat-baselines/web-fetch-20261103/NOTES.md), so the browser leaves them
// compressed: a body that starts with the gzip magic is inflated by the engine's decoders in the
// Data worker (WholeFileInflater, the path of a source named .gz: one decoder in the app, the
// command line's; swebok-improve SW8 D2), and a body that something on the way already decoded is
// kept as it is. Nothing is sent but the request for the file's URL.
import type { DownloadProgress, DownloadRequest, DownloadResult, FileDownloader } from '../../ports/example';
import { INFLATED_TOO_LARGE, type WholeFileInflater } from '../../ports/inflate';
import { dateOf } from '../../domain/download-sentences';
import { ask, contentLength, messageOf, PART_BYTES } from './http';

const GZIP_MAGIC = [0x1f, 0x8b] as const;

export interface DownloaderOptions {
  /** The Data worker's inflation of a whole file (the engine's decoders). */
  readonly inflater: WholeFileInflater;
  readonly fetch?: typeof fetch;
  readonly timeoutMs?: number;
}

export class BrowserFileDownloader implements FileDownloader {
  constructor(private readonly options: DownloaderOptions) {}

  async download(
    request: DownloadRequest,
    onProgress: (progress: DownloadProgress) => void,
    signal: AbortSignal,
  ): Promise<DownloadResult> {
    const answer = await ask(request.url, {}, signal, { fetch: this.options.fetch, timeoutMs: this.options.timeoutMs });
    if (!answer.ok) return { ok: false, reason: 'unreachable', detail: answer.detail };
    const { response } = answer;
    if (!response.ok) {
      await response.body?.cancel().catch(() => undefined);
      return { ok: false, reason: 'refused', status: response.status };
    }
    const result = await readBody(response, request, this.options.inflater, onProgress, signal);
    if (!result.ok) return result;
    const version = versionOf(response);
    return { ...result, ...(version === undefined ? {} : { version }) };
  }
}

type Bytes = Uint8Array<ArrayBuffer>;
type BodyResult = Exclude<DownloadResult, { readonly ok: true }> | { readonly ok: true; readonly data: Blob; readonly downloadedBytes: number };

async function readBody(
  response: Response,
  request: DownloadRequest,
  inflater: WholeFileInflater,
  onProgress: (progress: DownloadProgress) => void,
  signal: AbortSignal,
): Promise<BodyResult> {
  const total = contentLength(response);
  const reader = (response.body ?? new Blob([]).stream()).getReader();
  // The fetch's signal ends a browser's body too; cancelling the reader ends any body.
  const onAbort = () => void reader.cancel(signal.reason).catch(() => undefined);
  signal.addEventListener('abort', onAbort, { once: true });
  const parts: Blob[] = [];
  let chunks: Bytes[] = [];
  let chunkBytes = 0;
  let downloaded = 0;
  /** The first bytes, which say whether the body is gzip. */
  const head: number[] = [];
  try {
    for (;;) {
      let result: ReadableStreamReadResult<Bytes>;
      try {
        result = await reader.read();
      } catch (error) {
        if (signal.aborted) throw signal.reason;
        return { ok: false, reason: 'unreachable', detail: messageOf(error) };
      }
      if (result.done) break;
      const chunk = result.value;
      downloaded += chunk.byteLength;
      onProgress(total === undefined ? { bytes: downloaded } : { bytes: downloaded, total });
      // A gzip body is no larger than about what it inflates to (deflate adds at most 0.03 %), so
      // a body over the limit is over it once inflated too.
      if (downloaded > request.maxBytes) {
        await reader.cancel().catch(() => undefined);
        return { ok: false, reason: 'too-large' };
      }
      for (let i = 0; head.length < GZIP_MAGIC.length && i < chunk.byteLength; i++) head.push(chunk[i]!);
      chunks.push(chunk);
      chunkBytes += chunk.byteLength;
      if (chunkBytes >= PART_BYTES) {
        parts.push(new Blob(chunks));
        chunks = [];
        chunkBytes = 0;
      }
    }
    signal.throwIfAborted();
  } finally {
    signal.removeEventListener('abort', onAbort);
  }
  parts.push(new Blob(chunks));
  const body = new Blob(parts, { type: 'text/plain' });
  if (!GZIP_MAGIC.every((byte, i) => head[i] === byte)) return { ok: true, data: body, downloadedBytes: downloaded };
  // The engine's decoders inflate the body in the Data worker, as they inflate the source
  // `<fileName>.gz` (the set's file on its site), with the limit on the inflated bytes.
  let inflated: File;
  try {
    inflated = await inflater.inflateFile(new File([body], `${request.fileName}.gz`), request.maxBytes);
  } catch (error) {
    if (signal.aborted) throw signal.reason;
    if (error instanceof Error && error.name === INFLATED_TOO_LARGE) return { ok: false, reason: 'too-large' };
    return { ok: false, reason: 'corrupt', detail: messageOf(error) };
  }
  signal.throwIfAborted();
  return { ok: true, data: new Blob([inflated], { type: 'text/plain' }), downloadedBytes: downloaded };
}

/** UniProt's release (an exposed header), else the date of Last-Modified (a CORS-safelisted header). */
function versionOf(response: Response): string | undefined {
  const release = response.headers.get('X-UniProt-Release');
  if (release !== null && release.trim() !== '') return `UniProt release ${release.trim()}`;
  return dateOf(response.headers.get('Last-Modified'));
}
