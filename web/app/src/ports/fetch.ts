// Accession fetch port (line WF-1, docs/web/fetch_design.md): the FASTA of one accession from
// a public database. It runs only when the user asks ("Fetch"), never on load or on Run, and
// its request carries only the accession (web/AGENTS.md rule 4). The browser implementation
// is src/infra/browser/accession-fetcher.ts; tests use src/infra/fake/fake-fetcher.ts.
import type { FetchSourceId } from '../domain/accession';
import type { SequenceKind } from '../domain/programs';

export interface FetchRequest {
  readonly source: FetchSourceId;
  readonly kind: SequenceKind;
  readonly accession: string;
  /** The answer is dropped ('too-large') once it passes this many bytes (FETCH_BYTE_LIMIT in all). */
  readonly maxBytes: number;
  /** Ask the site, not the browser's HTTP cache: the race times the site (application/source-race.ts). */
  readonly fresh?: boolean;
}

export type FetchResult =
  | { readonly ok: true; readonly fasta: string }
  /** The source has no such record. */
  | { readonly ok: false; readonly reason: 'not-found' }
  | { readonly ok: false; readonly reason: 'too-large' }
  /** No answer: offline, blocked, or timed out. */
  | { readonly ok: false; readonly reason: 'unreachable'; readonly detail: string }
  /** An HTTP status other than "no such record" (such as 429, too many requests). */
  | { readonly ok: false; readonly reason: 'refused'; readonly status: number };

export interface AccessionFetcher {
  /** Rejects only when `signal` aborts. */
  fetch(request: FetchRequest, signal: AbortSignal): Promise<FetchResult>;
}
