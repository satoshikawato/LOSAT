// An accession fetch double for the tests (ports/fetch.ts): it answers from a table of FASTA
// texts by accession (or by "source:accession", which comes first), keeps the requests, and can
// hold the answers of some or all sources until `release`.
import { FETCH_SOURCES, type FetchSourceId } from '../../domain/accession';
import type { AccessionFetcher, FetchRequest, FetchResult } from '../../ports/fetch';

const ALL_SOURCES = FETCH_SOURCES.map((source) => source.id);

export class FakeAccessionFetcher implements AccessionFetcher {
  readonly requests: FetchRequest[] = [];
  private readonly holding = new Set<FetchSourceId>();
  private readonly held = new Map<FetchSourceId, Array<() => void>>();

  constructor(private readonly records: Readonly<Record<string, string | FetchResult>> = {}) {}

  /** Holds the answers of the sources (all when none is named) from now on, until `release`. */
  hold(...sources: FetchSourceId[]): void {
    for (const source of sources.length === 0 ? ALL_SOURCES : sources) this.holding.add(source);
  }

  /** Lets the held answers of the sources (all when none is named) come, in the order they were asked. */
  release(...sources: FetchSourceId[]): void {
    for (const source of sources.length === 0 ? ALL_SOURCES : sources) {
      this.holding.delete(source);
      const resumes = this.held.get(source) ?? [];
      this.held.delete(source);
      for (const resume of resumes) resume();
    }
  }

  async fetch(request: FetchRequest, signal: AbortSignal): Promise<FetchResult> {
    this.requests.push(request);
    if (this.holding.has(request.source)) {
      await new Promise<void>((resolve) => {
        const held = this.held.get(request.source) ?? [];
        held.push(resolve);
        this.held.set(request.source, held);
        signal.addEventListener('abort', () => resolve(), { once: true });
      });
    }
    signal.throwIfAborted();
    const record = this.records[`${request.source}:${request.accession}`] ?? this.records[request.accession];
    if (record === undefined) return { ok: false, reason: 'not-found' };
    if (typeof record !== 'string') return record;
    return record.length > request.maxBytes ? { ok: false, reason: 'too-large' } : { ok: true, fasta: record };
  }
}
