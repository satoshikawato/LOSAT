// What the page's downloads share (example-downloader.ts for the Example sets, counts-downloader.ts
// for the WindowMasker counts): the request (the file's URL only: no cookies, no referrer, a time
// limit on the answer's headers), and the reading of the headers that the sites send.

/** Milliseconds to the answer's headers; the body may take longer. */
export const TIMEOUT_MS = 60_000;
/** Downloaded bytes gathered into one part, so that what is held in memory stays small. */
export const PART_BYTES = 8 * 1024 * 1024;

export type Headers = Readonly<Record<string, string>>;

export type Answer =
  | { readonly ok: true; readonly response: Response }
  /** No answer: offline, blocked, or timed out. */
  | { readonly ok: false; readonly detail: string };

/**
 * Asks for the URL and waits for the answer's headers. Rejects only when `signal` aborts. The
 * request is a CORS request without credentials or referrer (the sites send CORS headers; the CSP's
 * connect-src names them, public/_headers).
 */
export async function ask(
  url: string,
  headers: Headers,
  signal: AbortSignal,
  options: { readonly fetch?: typeof fetch | undefined; readonly timeoutMs?: number | undefined } = {},
): Promise<Answer> {
  signal.throwIfAborted();
  const controller = new AbortController();
  const onAbort = () => controller.abort(signal.reason);
  signal.addEventListener('abort', onAbort, { once: true });
  let timedOut = false;
  const timer = setTimeout(() => {
    timedOut = true;
    controller.abort(new DOMException('no answer in time', 'TimeoutError'));
  }, options.timeoutMs ?? TIMEOUT_MS);
  try {
    const response = await (options.fetch ?? fetch)(url, {
      mode: 'cors',
      signal: controller.signal,
      credentials: 'omit',
      referrerPolicy: 'no-referrer',
      ...(Object.keys(headers).length === 0 ? {} : { headers }),
    });
    return { ok: true, response };
  } catch (error) {
    if (signal.aborted) throw error;
    return { ok: false, detail: timedOut ? 'no answer in time' : messageOf(error) };
  } finally {
    clearTimeout(timer);
    signal.removeEventListener('abort', onAbort);
  }
}

/** The Content-Length of an answer, if it sends one. */
export function contentLength(response: Response): number | undefined {
  const value = Number(response.headers.get('Content-Length') ?? Number.NaN);
  return Number.isSafeInteger(value) && value >= 0 ? value : undefined;
}

export function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
