// An accession fetch double for the tests (ports/fetch.ts): it answers from a table of FASTA
// texts by accession, keeps the requests, and can hold its answers until `release`.
import type { AccessionFetcher, FetchRequest, FetchResult } from '../../ports/fetch';

export class FakeAccessionFetcher implements AccessionFetcher {
  readonly requests: FetchRequest[] = [];
  private held: Array<() => void> | undefined;

  constructor(private readonly records: Readonly<Record<string, string | FetchResult>> = {}) {}

  /** Holds the answers from now on, until `release`. */
  hold(): void {
    this.held ??= [];
  }

  release(): void {
    const held = this.held ?? [];
    this.held = undefined;
    for (const resume of held) resume();
  }

  async fetch(request: FetchRequest, signal: AbortSignal): Promise<FetchResult> {
    this.requests.push(request);
    const held = this.held;
    await (held === undefined
      ? Promise.resolve()
      : new Promise<void>((resolve) => {
          held.push(resolve);
          signal.addEventListener('abort', () => resolve(), { once: true });
        }));
    signal.throwIfAborted();
    const record = this.records[request.accession];
    if (record === undefined) return { ok: false, reason: 'not-found' };
    if (typeof record !== 'string') return record;
    return record.length > request.maxBytes ? { ok: false, reason: 'too-large' } : { ok: true, fasta: record };
  }
}
