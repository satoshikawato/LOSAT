// POST /api/visit: the Cloudflare Pages Function that receives LOSAT Web's one usage event
// and forwards it to the GA4 Measurement Protocol (plan §5.10 (b), docs/web/analytics.md).
// Pages runs it from `functions/` beside the build output when the project's root
// directory is web/app. It accepts only the exact event that the page sends, from the
// site's own origin, and forwards nothing else: not the request's headers, address,
// cookies or referrer. Without the two secrets (set by the maintainer in Cloudflare, never
// in the repository) it forwards nothing and still answers 204, so a preview behaves as
// production does without measuring.
import { parseVisitBody, VISIT_BODY_MAX_BYTES, type VisitEvent } from '../../src/domain/analytics';

/** The environment variables that the Function reads; both are needed to forward. */
export interface VisitEnv {
  readonly GA4_MEASUREMENT_ID?: string;
  readonly GA4_API_SECRET?: string;
}

/** The part of a Pages Function's context that this Function uses. */
export interface VisitContext {
  readonly request: Request;
  readonly env: VisitEnv;
  waitUntil(promise: Promise<unknown>): void;
}

export const MEASUREMENT_PROTOCOL_URL = 'https://www.google-analytics.com/mp/collect';

/**
 * Pages' `_headers` file does not apply to a Function's responses, so the answers carry their
 * own: never cached, readable only by the site itself.
 */
const HEADERS = {
  'Cache-Control': 'no-store',
  'Cross-Origin-Resource-Policy': 'same-origin',
  'X-Content-Type-Options': 'nosniff',
} as const;

const answer = (status: number, extra: Record<string, string> = {}) => new Response(null, { status, headers: { ...HEADERS, ...extra } });

/** `application/json`, with at most a UTF-8 charset parameter. */
function isJson(contentType: string | null): boolean {
  if (contentType === null) return false;
  const [type, ...parameters] = contentType.split(';').map((part) => part.trim().toLowerCase());
  return type === 'application/json' && parameters.every((parameter) => parameter === 'charset=utf-8');
}

/** The body as UTF-8 text, or undefined when it is longer than the limit or not UTF-8. */
async function readBody(request: Request): Promise<string | undefined> {
  const declared = request.headers.get('Content-Length');
  if (declared !== null && !(Number(declared) <= VISIT_BODY_MAX_BYTES)) return undefined;
  if (request.body === null) return '';
  const reader = request.body.getReader();
  const chunks: Uint8Array[] = [];
  let length = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    length += value.byteLength;
    if (length > VISIT_BODY_MAX_BYTES) {
      await reader.cancel();
      return undefined;
    }
    chunks.push(value);
  }
  const bytes = new Uint8Array(length);
  let offset = 0;
  for (const chunk of chunks) {
    bytes.set(chunk, offset);
    offset += chunk.byteLength;
  }
  try {
    return new TextDecoder('utf-8', { fatal: true, ignoreBOM: true }).decode(bytes);
  } catch {
    return undefined;
  }
}

/** The Measurement Protocol request of an event: its IDs, and nothing of the visitor's request. */
export function measurementRequest(event: VisitEvent, measurementId: string, apiSecret: string): Request {
  const url = `${MEASUREMENT_PROTOCOL_URL}?measurement_id=${encodeURIComponent(measurementId)}&api_secret=${encodeURIComponent(apiSecret)}`;
  return new Request(url, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      client_id: event.client,
      non_personalized_ads: true,
      events: [{ name: 'app_visit', params: { release_id: event.release } }],
    }),
  });
}

/** Checks a request and forwards its event; `send` is `fetch` (tests give a fake). */
export async function handleVisit(context: VisitContext, send: (request: Request) => Promise<unknown>): Promise<Response> {
  const { request, env } = context;
  if (request.method !== 'POST') return answer(405, { Allow: 'POST' });
  if (request.headers.get('Origin') !== new URL(request.url).origin) return answer(400);
  if (!isJson(request.headers.get('Content-Type'))) return answer(400);
  const text = await readBody(request);
  const event = text === undefined ? undefined : parseVisitBody(text);
  if (event === undefined) return answer(400);
  const measurementId = env.GA4_MEASUREMENT_ID ?? '';
  const apiSecret = env.GA4_API_SECRET ?? '';
  if (measurementId !== '' && apiSecret !== '') {
    // After the answer; a failure is dropped (no retry, no report).
    context.waitUntil(send(measurementRequest(event, measurementId, apiSecret)).then(
      () => undefined,
      () => undefined,
    ));
  }
  return answer(204);
}

/** Pages calls this for every method of /api/visit. */
export const onRequest = (context: VisitContext): Promise<Response> => handleVisit(context, (request) => fetch(request));
