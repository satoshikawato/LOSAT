// The browser's download of WindowMasker counts (ports/counts.ts, docs/web/fetch_design.md §10): the
// file is asked straight from the page (NCBI's FTP sends CORS headers and accepts a Range; the CSP's
// connect-src names it), with no cookies and no referrer, and only when the user asks. The body is
// handed to the store in parts of about 8 MiB as it arrives, so the 724 MB of a human is never held.
// Continue asks for `Range: bytes=<saved>-`; the saved part is dropped, and the file asked from its
// start, when the answer is not a 206, its Content-Range starts elsewhere, or its Last-Modified is
// not the first answer's (the file changed on the site). There is no automatic retry.
import { COUNTS_BUSY, COUNTS_CHANGED, type CountsCopy, type CountsDownloader, type CountsRequest, type CountsResult, type CountsStore } from '../../ports/counts';
import type { DownloadProgress } from '../../ports/example';
import { ask, contentLength, messageOf, PART_BYTES } from './http';

export interface CountsDownloaderOptions {
  readonly store: Pick<CountsStore, 'beginCounts' | 'appendCounts' | 'finishCounts' | 'closeCounts'>;
  readonly fetch?: typeof fetch;
  readonly timeoutMs?: number;
}

type Bytes = Uint8Array<ArrayBuffer>;

/** The answer that the body is read from, and the byte of the file it starts at. */
interface Body {
  readonly response: Response;
  readonly from: number;
  readonly size: number;
  readonly lastModified: string;
}

export class BrowserCountsDownloader implements CountsDownloader {
  constructor(private readonly options: CountsDownloaderOptions) {}

  async download(request: CountsRequest, onProgress: (progress: DownloadProgress) => void, signal: AbortSignal): Promise<CountsResult> {
    const body = await this.answer(request, signal);
    if ('ok' in body) return body;
    let token: string | undefined;
    try {
      token = await this.options.store.beginCounts(request.taxid, {
        url: request.url,
        lastModified: body.lastModified,
        size: body.size,
        from: body.from,
      });
    } catch (error) {
      await body.response.body?.cancel().catch(() => undefined);
      const name = (error as { name?: unknown }).name;
      if (name === COUNTS_BUSY) return { ok: false, reason: 'busy' };
      if (name === COUNTS_CHANGED) return { ok: false, reason: 'changed' };
      return { ok: false, reason: 'not-kept', detail: messageOf(error) };
    }
    return this.read(request, body, token, onProgress, signal);
  }

  /**
   * Asks for the file, from the saved byte when there is a partial copy. An answer that does not
   * continue that copy is not read: the file is asked again from its start.
   */
  private async answer(request: CountsRequest, signal: AbortSignal): Promise<Body | Extract<CountsResult, { ok: false }>> {
    const resume = request.resume !== undefined && request.resume.savedBytes > 0 ? request.resume : undefined;
    for (const range of resume === undefined ? [undefined] : [resume, undefined]) {
      // A copy whose bytes all came but that was not finished asks for its last byte again: the
      // site answers 416 for a range that starts at the end.
      const from = range === undefined ? 0 : Math.min(range.savedBytes, range.size - 1);
      const answer = await ask(request.url, from === 0 ? {} : { Range: `bytes=${from}-` }, signal, {
        fetch: this.options.fetch,
        timeoutMs: this.options.timeoutMs,
      });
      if (!answer.ok) return { ok: false, reason: 'unreachable', detail: answer.detail };
      const { response } = answer;
      const lastModified = response.headers.get('Last-Modified') ?? '';
      if (from > 0 && response.status !== 200) {
        // A 206 that starts at the saved byte of the same file (its size and Last-Modified) continues
        // the saved part. Anything
        // else (416 when the saved part is longer than the file now, another start, another
        // file) is not the file's rest. A 200 is the whole file: the site ignored the Range.
        const start = contentRangeStart(response);
        const size = contentRangeSize(response);
        if (response.status === 206 && start === from && size === range!.size && lastModified === range!.lastModified) {
          return { response, from, size, lastModified };
        }
        await response.body?.cancel().catch(() => undefined);
        // Another file (a 206 of another start, size or date; 416 for a range past its end) is asked
        // from its start. Any other status (404, 503) fails as a first answer would, and the saved
        // part stays for another Continue: there is no automatic retry.
        if (response.status !== 206 && response.status !== 416) return { ok: false, reason: 'refused', status: response.status };
        continue;
      }
      return this.fromStart(response, lastModified);
    }
    // Unreachable: the last attempt starts from 0 and returns above.
    return { ok: false, reason: 'unreachable', detail: 'no answer' };
  }

  /** The answer that holds the file from its first byte. */
  private async fromStart(response: Response, lastModified: string): Promise<Body | Extract<CountsResult, { ok: false }>> {
    if (response.status !== 200 && response.status !== 206) {
      await response.body?.cancel().catch(() => undefined);
      return { ok: false, reason: 'refused', status: response.status };
    }
    const size = response.status === 206 ? contentRangeSize(response) : contentLength(response);
    if (size === undefined || (response.status === 206 && contentRangeStart(response) !== 0)) {
      await response.body?.cancel().catch(() => undefined);
      return { ok: false, reason: 'unreachable', detail: 'the answer does not say how large the file is' };
    }
    return { response, from: 0, size, lastModified };
  }

  private async read(
    request: CountsRequest,
    body: Body,
    token: string | undefined,
    onProgress: (progress: DownloadProgress) => void,
    signal: AbortSignal,
  ): Promise<CountsResult> {
    const store = this.options.store;
    const reader = (body.response.body ?? new Blob([]).stream()).getReader();
    // The fetch's signal ended with the answer's headers; cancelling the reader ends the body.
    const onAbort = () => void reader.cancel(signal.reason).catch(() => undefined);
    signal.addEventListener('abort', onAbort, { once: true });
    // A Cancel while the store opened the copy (after the answer's headers, before this listener).
    if (signal.aborted) onAbort();
    /** The parts of a browser that keeps nothing: this tab's copy. */
    const kept: Blob[] = [];
    let chunks: Bytes[] = [];
    let chunkBytes = 0;
    let received = 0;
    let failure: Extract<CountsResult, { ok: false }> | undefined;
    const flush = async () => {
      if (chunkBytes === 0) return;
      const part = joined(chunks, chunkBytes);
      chunks = [];
      chunkBytes = 0;
      if (token === undefined) kept.push(new Blob([part]));
      else await store.appendCounts(token, part);
    };
    const total = body.size;
    try {
      onProgress({ bytes: body.from, total });
      for (;;) {
        let result: ReadableStreamReadResult<Bytes>;
        try {
          result = await reader.read();
        } catch (error) {
          if (signal.aborted) break;
          failure = { ok: false, reason: 'unreachable', detail: messageOf(error) };
          break;
        }
        if (result.done) break;
        chunks.push(result.value);
        chunkBytes += result.value.byteLength;
        received += result.value.byteLength;
        onProgress({ bytes: body.from + received, total });
        if (chunkBytes >= PART_BYTES) await flush();
      }
      // What came is kept, also when the download stops: Continue asks for the rest.
      await flush();
    } catch (error) {
      // The storage refused a part (it is full), or could not be written.
      await reader.cancel().catch(() => undefined);
      if (token !== undefined) await store.closeCounts(token).catch(() => undefined);
      return { ok: false, reason: 'not-kept', detail: messageOf(error) };
    } finally {
      signal.removeEventListener('abort', onAbort);
    }
    if (signal.aborted) {
      if (token !== undefined) await store.closeCounts(token).catch(() => undefined);
      throw signal.reason;
    }
    if (failure === undefined && body.from + received !== body.size) {
      failure = { ok: false, reason: 'unreachable', detail: 'the connection ended early' };
    }
    if (failure !== undefined) {
      if (token !== undefined) await store.closeCounts(token).catch(() => undefined);
      return failure;
    }
    const copyOf = (savedAt: number): CountsCopy => ({
      taxid: request.taxid,
      url: request.url,
      lastModified: body.lastModified,
      size: body.size,
      savedBytes: body.size,
      complete: true,
      savedAt,
    });
    if (token === undefined) {
      return { ok: true, copy: copyOf(Date.now()), file: new File(kept, 'wmasker.obinary', { type: 'application/octet-stream' }), kept: false };
    }
    try {
      const saved = await store.finishCounts(token);
      if (saved.file === undefined) throw new Error('the saved file cannot be read');
      return { ok: true, copy: saved.copy, file: saved.file, kept: true };
    } catch (error) {
      await store.closeCounts(token).catch(() => undefined);
      return { ok: false, reason: 'not-kept', detail: messageOf(error) };
    }
  }
}

function joined(chunks: readonly Bytes[], bytes: number): Bytes {
  if (chunks.length === 1) return chunks[0]!;
  const part = new Uint8Array(bytes);
  let at = 0;
  for (const chunk of chunks) {
    part.set(chunk, at);
    at += chunk.byteLength;
  }
  return part;
}

/** `bytes a-b/total` of a 206's Content-Range. */
function contentRange(response: Response): { readonly start: number; readonly size: number | undefined } | undefined {
  const match = /^bytes (\d+)-(\d+)\/(\d+|\*)$/.exec(response.headers.get('Content-Range') ?? '');
  if (match === null) return undefined;
  return { start: Number(match[1]), size: match[3] === '*' ? undefined : Number(match[3]) };
}

function contentRangeStart(response: Response): number | undefined {
  return contentRange(response)?.start;
}

function contentRangeSize(response: Response): number | undefined {
  return contentRange(response)?.size;
}
