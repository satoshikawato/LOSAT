// The Pages Function of /api/visit (functions/api/visit.ts) with a fake `fetch`.
import { describe, expect, it } from 'vitest';
import { handleVisit, MEASUREMENT_PROTOCOL_URL, onRequest, type VisitEnv } from '../../functions/api/visit';
import { visitBody } from '../../src/domain/analytics';

const SITE = 'https://losat.example';
const RELEASE = '0123456789ab-1a2b3c4d';
const CLIENT = '3f2b8c1e-5d4a-4b7c-9e1f-0a2b3c4d5e6f';
const SECRETS: VisitEnv = { GA4_MEASUREMENT_ID: 'G-TEST000000', GA4_API_SECRET: 'test-secret' };

interface Call {
  readonly response: Response;
  readonly forwarded: Request[];
}

async function call(init: RequestInit, env: VisitEnv = SECRETS, path = '/api/visit'): Promise<Call> {
  const forwarded: Request[] = [];
  const pending: Promise<unknown>[] = [];
  const request = new Request(`${SITE}${path}`, init);
  const response = await handleVisit({ request, env, waitUntil: (promise) => pending.push(promise) }, async (outgoing) => {
    forwarded.push(outgoing);
    return new Response(null, { status: 204 });
  });
  await Promise.all(pending);
  return { response, forwarded };
}

const headers = (extra: Record<string, string> = {}) => ({ Origin: SITE, 'Content-Type': 'application/json', ...extra });
const post = (body: BodyInit | null, extra: Record<string, string> = {}) => ({ method: 'POST', headers: headers(extra), body });

describe('/api/visit', () => {
  it('answers 204 to the exact event and forwards only its two IDs to the Measurement Protocol', async () => {
    const { response, forwarded } = await call(
      post(visitBody(RELEASE, CLIENT), { Cookie: 'a=b', Referer: `${SITE}/r/x/`, 'User-Agent': 'probe', 'CF-Connecting-IP': '192.0.2.1' }),
    );
    expect(response.status).toBe(204);
    expect(response.headers.get('Cache-Control')).toBe('no-store');
    expect(forwarded).toHaveLength(1);
    const outgoing = forwarded[0]!;
    expect(outgoing.method).toBe('POST');
    const url = new URL(outgoing.url);
    expect(`${url.origin}${url.pathname}`).toBe(MEASUREMENT_PROTOCOL_URL);
    expect(Object.fromEntries(url.searchParams)).toEqual({ measurement_id: 'G-TEST000000', api_secret: 'test-secret' });
    // Nothing of the visitor's request: only the content type.
    expect([...outgoing.headers.keys()]).toEqual(['content-type']);
    expect(await outgoing.json()).toEqual({
      client_id: CLIENT,
      non_personalized_ads: true,
      events: [{ name: 'app_visit', params: { release_id: RELEASE } }],
    });
  });

  it('accepts a UTF-8 charset on the content type', async () => {
    const { response, forwarded } = await call(post(visitBody(RELEASE, CLIENT), { 'Content-Type': 'application/json; charset=UTF-8' }));
    expect(response.status).toBe(204);
    expect(forwarded).toHaveLength(1);
  });

  it('forwards nothing without both secrets, and still answers 204', async () => {
    for (const env of [{}, { GA4_MEASUREMENT_ID: 'G-TEST000000' }, { GA4_API_SECRET: 'test-secret' }, { GA4_MEASUREMENT_ID: '', GA4_API_SECRET: '' }]) {
      const { response, forwarded } = await call(post(visitBody(RELEASE, CLIENT)), env);
      expect(response.status, JSON.stringify(env)).toBe(204);
      expect(forwarded).toEqual([]);
    }
  });

  it('answers 405 to other methods and forwards nothing', async () => {
    for (const method of ['GET', 'HEAD', 'PUT', 'DELETE', 'OPTIONS']) {
      const { response, forwarded } = await call({ method, headers: headers() });
      expect(response.status, method).toBe(405);
      expect(response.headers.get('Allow')).toBe('POST');
      expect(forwarded).toEqual([]);
    }
  });

  it('refuses another origin, no origin, another content type, and every body but the exact event, forwarding nothing', async () => {
    const exact = visitBody(RELEASE, CLIENT);
    const refused: [string, RequestInit][] = [
      ['another origin', post(exact, { Origin: 'https://evil.example' })],
      ['an origin of another scheme', post(exact, { Origin: 'http://losat.example' })],
      ['the null origin', post(exact, { Origin: 'null' })],
      ['no origin', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: exact }],
      ['text/plain', post(exact, { 'Content-Type': 'text/plain' })],
      ['a form', post(exact, { 'Content-Type': 'application/x-www-form-urlencoded' })],
      ['another charset', post(exact, { 'Content-Type': 'application/json; charset=latin1' })],
      ['no content type', { method: 'POST', headers: { Origin: SITE }, body: new Blob([exact]) }],
      ['no body', post(null)],
      ['an empty body', post('')],
      ['another key', post(`{"v":1,"release":"${RELEASE}","client":"${CLIENT}","name":"seq1"}`)],
      ['another version', post(`{"v":2,"release":"${RELEASE}","client":"${CLIENT}"}`)],
      ['a bad release ID', post(`{"v":1,"release":"../x","client":"${CLIENT}"}`)],
      ['an upper-case client ID', post(`{"v":1,"release":"${RELEASE}","client":"${CLIENT.toUpperCase()}"}`)],
      ['white space', post(`{ "v":1,"release":"${RELEASE}","client":"${CLIENT}"}`)],
      ['an upper-case release ID', post(`{"v":1,"release":"ACGT","client":"${CLIENT}"}`)],
      ['more than 512 bytes', post(`${exact}${' '.repeat(600)}`)],
      ['a declared length over 512', post(exact, { 'Content-Length': '513' })],
      ['bytes that are not UTF-8', post(new Uint8Array([0x7b, 0xff, 0x7d]))],
    ];
    for (const [what, init] of refused) {
      const { response, forwarded } = await call(init);
      expect(response.status, what).toBe(400);
      expect(forwarded, what).toEqual([]);
    }
  });

  it('a failing Measurement Protocol request changes nothing for the page', async () => {
    const request = new Request(`${SITE}/api/visit`, post(visitBody(RELEASE, CLIENT)));
    const pending: Promise<unknown>[] = [];
    const response = await handleVisit({ request, env: SECRETS, waitUntil: (promise) => pending.push(promise) }, () => Promise.reject(new Error('down')));
    expect(response.status).toBe(204);
    await expect(Promise.all(pending)).resolves.toEqual([undefined]);
  });

  it('is exported as the handler of every method', () => {
    expect(typeof onRequest).toBe('function');
  });
});
