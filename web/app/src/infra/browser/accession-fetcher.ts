// The browser's accession fetch (ports/fetch.ts): NCBI E-utilities, ENA and UniProt (EBI), and
// CNCB GenBase, each asked straight from the page (they send CORS headers; the CSP's
// connect-src names them, public/_headers). A request carries only the accession: no
// cookies, no referrer. The answers to an unknown accession, measured on 2026-11-03
// (losat-baselines/web-fetch-20261103): NCBI efetch 200 with "Error: …" (nuccore) or 400
// (protein), ENA 404 (400 for a RefSeq accession), UniProt 400, GenBase 200 with no body.
import type { FetchSourceId } from '../../domain/accession';
import type { SequenceKind } from '../../domain/programs';
import type { AccessionFetcher, FetchRequest, FetchResult } from '../../ports/fetch';

/** NCBI takes 3 requests a second without an API key (E-utilities usage guidelines). */
const NCBI_INTERVAL_MS = 340;
/** Milliseconds to the answer's headers; the body (up to FETCH_BYTE_LIMIT) may take longer. */
const TIMEOUT_MS = 60_000;

export function fetchUrl(source: FetchSourceId, kind: SequenceKind, accession: string): string {
  const id = encodeURIComponent(accession);
  switch (source) {
    case 'ncbi':
      return `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${kind === 'nucleotide' ? 'nuccore' : 'protein'}&id=${id}&rettype=fasta&retmode=text`;
    case 'ebi':
      return kind === 'nucleotide'
        ? `https://www.ebi.ac.uk/ena/browser/api/fasta/${id}`
        : `https://rest.uniprot.org/uniprotkb/${id}.fasta`;
    case 'cncb':
      return `https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=${id}`;
  }
}

export interface FetcherOptions {
  readonly fetch?: typeof fetch;
  readonly now?: () => number;
  readonly timeoutMs?: number;
}

export class BrowserAccessionFetcher implements AccessionFetcher {
  /** When the next NCBI request may start. */
  private ncbiNext = 0;

  constructor(private readonly options: FetcherOptions = {}) {}

  async fetch(request: FetchRequest, signal: AbortSignal): Promise<FetchResult> {
    if (request.source === 'ncbi') await this.ncbiTurn(signal);
    signal.throwIfAborted();
    const timeout = new AbortController();
    const onAbort = () => timeout.abort(signal.reason);
    signal.addEventListener('abort', onAbort, { once: true });
    let timedOut = false;
    const timer = setTimeout(() => {
      timedOut = true;
      timeout.abort(new DOMException('no answer in time', 'TimeoutError'));
    }, this.options.timeoutMs ?? TIMEOUT_MS);
    try {
      const response = await (this.options.fetch ?? fetch)(fetchUrl(request.source, request.kind, request.accession), {
        signal: timeout.signal,
        credentials: 'omit',
        referrerPolicy: 'no-referrer',
        ...(request.fresh === true ? { cache: 'no-store' as const } : {}),
      });
      clearTimeout(timer);
      if (response.status === 400 || response.status === 404) return { ok: false, reason: 'not-found' };
      if (!response.ok) return { ok: false, reason: 'refused', status: response.status };
      const text = await readText(response, request.maxBytes);
      if (text === undefined) return { ok: false, reason: 'too-large' };
      // NCBI's "Error: …" and GenBase's empty body: no such record.
      if (!text.trimStart().startsWith('>')) return { ok: false, reason: 'not-found' };
      return { ok: true, fasta: text.endsWith('\n') ? text : `${text}\n` };
    } catch (error) {
      if (signal.aborted) throw error;
      return { ok: false, reason: 'unreachable', detail: timedOut ? 'no answer in time' : messageOf(error) };
    } finally {
      clearTimeout(timer);
      signal.removeEventListener('abort', onAbort);
    }
  }

  private async ncbiTurn(signal: AbortSignal): Promise<void> {
    const now = (this.options.now ?? (() => performance.now()))();
    const wait = Math.max(0, this.ncbiNext - now);
    this.ncbiNext = now + wait + NCBI_INTERVAL_MS;
    if (wait > 0) await delay(wait, signal);
  }
}

/** The body as text, or undefined once it passes `maxBytes` (the rest is not read). */
async function readText(response: Response, maxBytes: number): Promise<string | undefined> {
  if (response.body === null) return '';
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  let bytes = 0;
  let text = '';
  for (;;) {
    const { done, value } = await reader.read();
    if (done) return text + decoder.decode();
    bytes += value.byteLength;
    if (bytes > maxBytes) {
      await reader.cancel();
      return undefined;
    }
    text += decoder.decode(value, { stream: true });
  }
}

function delay(ms: number, signal: AbortSignal): Promise<void> {
  return new Promise((resolve, reject) => {
    if (signal.aborted) {
      reject(signal.reason as Error);
      return;
    }
    const timer = setTimeout(() => {
      signal.removeEventListener('abort', abort);
      resolve();
    }, ms);
    const abort = () => {
      clearTimeout(timer);
      reject(signal.reason as Error);
    };
    signal.addEventListener('abort', abort, { once: true });
  });
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
