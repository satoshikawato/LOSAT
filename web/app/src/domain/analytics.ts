// Usage statistics (plan §5.10 (b), docs/web/analytics.md): the one event that LOSAT Web
// sends after consent, and the choice that the browser keeps. The page posts the event to
// the site's own endpoint (a Cloudflare Pages Function, functions/api/visit.ts), which
// forwards it to the GA4 Measurement Protocol; no document loads a Google script. The event
// says only that LOSAT Web was opened and which release; it never holds research data.
import { RELEASE_ID_PATTERN } from './site-index';

/** A random client ID: `crypto.randomUUID()`, a lower-case UUID version 4. */
export const CLIENT_ID_PATTERN = /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/;

/** The largest body that the endpoint accepts, in bytes; the event is under 130. */
export const VISIT_BODY_MAX_BYTES = 512;

/** The event: `{"v":1,"release":<release ID>,"client":<client ID>}`. */
export interface VisitEvent {
  readonly v: 1;
  readonly release: string;
  readonly client: string;
}

/** The body of the event, exactly as the page sends it; throws on an ID that does not match its pattern. */
export function visitBody(release: string, client: string): string {
  if (!RELEASE_ID_PATTERN.test(release)) throw new Error('visit: the release ID is not valid');
  if (!CLIENT_ID_PATTERN.test(client)) throw new Error('visit: the client ID is not valid');
  return JSON.stringify({ v: 1, release, client });
}

/**
 * Reads a body that the endpoint received. Only the exact text that `visitBody` writes is an
 * event: other keys, other values, white space, duplicate keys or a different key order are
 * not, so nothing but the two IDs can pass through.
 */
export function parseVisitBody(text: string): VisitEvent | undefined {
  let value: unknown;
  try {
    value = JSON.parse(text);
  } catch {
    return undefined;
  }
  if (typeof value !== 'object' || value === null || Array.isArray(value)) return undefined;
  const { release, client } = value as Record<string, unknown>;
  if (typeof release !== 'string' || typeof client !== 'string') return undefined;
  let canonical: string;
  try {
    canonical = visitBody(release, client);
  } catch {
    return undefined;
  }
  return canonical === text ? Object.freeze({ v: 1, release, client }) : undefined;
}

/** The choice that the browser keeps: consent with the client ID made at consent, or none. */
export type AnalyticsChoice = { readonly consent: 'granted'; readonly client: string } | { readonly consent: 'denied' };

/** The stored form of a choice. */
export function serializeAnalyticsChoice(choice: AnalyticsChoice): string {
  return JSON.stringify(choice.consent === 'granted' ? { v: 1, consent: 'granted', client: choice.client } : { v: 1, consent: 'denied' });
}

/** Reads a stored choice; undefined for none or anything else (the page then asks again). */
export function parseAnalyticsChoice(text: string | null): AnalyticsChoice | undefined {
  if (text === null) return undefined;
  let value: unknown;
  try {
    value = JSON.parse(text);
  } catch {
    return undefined;
  }
  if (typeof value !== 'object' || value === null || Array.isArray(value)) return undefined;
  const stored = value as Record<string, unknown>;
  if (stored['v'] !== 1) return undefined;
  if (stored['consent'] === 'denied') return Object.freeze({ consent: 'denied' });
  if (stored['consent'] === 'granted' && typeof stored['client'] === 'string' && CLIENT_ID_PATTERN.test(stored['client'])) {
    return Object.freeze({ consent: 'granted', client: stored['client'] });
  }
  return undefined;
}
