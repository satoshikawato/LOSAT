// The local site server (build/site-server.ts) answers as Cloudflare Pages does: the site's
// `_headers` by path on every response, ETag and 304, `/` -> index.html, the .html redirects,
// 404 with 404.html and no single-page fallback, and a site that changes while it runs.
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { DEFAULT_CACHE_CONTROL, startSiteServer, type SiteRequest, type SiteServer } from '../../build/site-server';

const HEADERS = `/*
  Cross-Origin-Opener-Policy: same-origin
  Cross-Origin-Embedder-Policy: require-corp

/r/*
  Cache-Control: public, max-age=31536000, immutable

/
  Cache-Control: no-cache
`;

let site: string;
let server: SiteServer;
const requests: SiteRequest[] = [];

function write(path: string, text: string | Uint8Array): void {
  mkdirSync(join(site, path, '..'), { recursive: true });
  writeFileSync(join(site, path), text);
}

beforeAll(async () => {
  site = mkdtempSync(join(tmpdir(), 'losat-site-server-'));
  write('_headers', HEADERS);
  write('index.html', '<!doctype html><title>boot</title>\n');
  write('404.html', '<!doctype html><title>Not found</title>\n');
  write('releases.json', '{"format":1}\n');
  write('boot.css', '#boot{}\n');
  write('r/rel-1/assets/index-1.js', 'console.log(1);\n');
  write('r/rel-1/engine/x.wasm', new Uint8Array([0, 0x61, 0x73, 0x6d]));
  server = await startSiteServer({ root: site, onRequest: (request) => requests.push(request) });
});

afterAll(async () => {
  await server?.close();
  rmSync(site, { recursive: true, force: true });
});

const get = (path: string, init: RequestInit = {}) => fetch(`${server.origin}${path}`, { redirect: 'manual', ...init });

describe('the local site server', () => {
  it('serves / from index.html with the rules of / and /*, and an ETag', async () => {
    const response = await get('/');
    expect(response.status).toBe(200);
    expect(await response.text()).toContain('<title>boot</title>');
    expect(response.headers.get('content-type')).toBe('text/html; charset=utf-8');
    expect(response.headers.get('cross-origin-opener-policy')).toBe('same-origin');
    expect(response.headers.get('cross-origin-embedder-policy')).toBe('require-corp');
    expect(response.headers.get('cache-control')).toBe('no-cache');
    expect(response.headers.get('etag')).toMatch(/^"[0-9a-f]{32}"$/);
  });

  it('answers a matching If-None-Match with 304 and the same rule headers', async () => {
    const etag = (await get('/r/rel-1/assets/index-1.js')).headers.get('etag')!;
    const response = await get('/r/rel-1/assets/index-1.js', { headers: { 'If-None-Match': `W/"other", ${etag}` } });
    expect(response.status).toBe(304);
    expect(await response.text()).toBe('');
    expect(response.headers.get('etag')).toBe(etag);
    expect(response.headers.get('cross-origin-embedder-policy')).toBe('require-corp');
    expect(response.headers.get('cache-control')).toBe('public, max-age=31536000, immutable');
    expect((await get('/r/rel-1/assets/index-1.js', { headers: { 'If-None-Match': '"other"' } })).status).toBe(200);
  });

  it('gives each file its content type, and Pages’ default Cache-Control where no rule sets one', async () => {
    const wasm = await get('/r/rel-1/engine/x.wasm');
    expect(wasm.headers.get('content-type')).toBe('application/wasm');
    expect(new Uint8Array(await wasm.arrayBuffer())).toEqual(new Uint8Array([0, 0x61, 0x73, 0x6d]));
    expect((await get('/r/rel-1/assets/index-1.js')).headers.get('content-type')).toBe('text/javascript; charset=utf-8');
    const json = await get('/releases.json');
    expect(json.headers.get('content-type')).toBe('application/json');
    expect(json.headers.get('cache-control')).toBe(DEFAULT_CACHE_CONTROL);
    expect((await get('/boot.css')).headers.get('content-type')).toBe('text/css; charset=utf-8');
  });

  it('answers a missing file with 404 and 404.html, with the rules of the path asked for, never with index.html', async () => {
    for (const path of ['/r/rel-1/assets/missing.js', '/search', '/r/rel-2/release-manifest.json']) {
      const response = await get(path);
      expect(response.status).toBe(404);
      expect(await response.text()).toContain('<title>Not found</title>');
      expect(response.headers.get('cross-origin-embedder-policy')).toBe('require-corp');
    }
    expect((await get('/r/rel-2/release-manifest.json')).headers.get('cache-control')).toBe('public, max-age=31536000, immutable');
  });

  it('redirects .html paths as Pages does, and does not serve its own configuration files', async () => {
    const index = await get('/index.html?x=1');
    expect(index.status).toBe(308);
    expect(index.headers.get('location')).toBe('/?x=1');
    expect(index.headers.get('cross-origin-opener-policy')).toBe('same-origin');
    expect((await get('/404.html')).headers.get('location')).toBe('/404');
    expect((await get('/404')).status).toBe(200);
    expect((await get('/_headers')).status).toBe(404);
  });

  it('serves nothing outside the site directory', async () => {
    for (const path of ['/%2e%2e/%2e%2e/etc/passwd', '/r/..%2f..%2fetc/passwd', '/r/%00']) {
      expect((await get(path)).status).toBe(404);
    }
  });

  it('answers HEAD without a body and other methods with 405, and reports every request with its body', async () => {
    const head = await get('/', { method: 'HEAD' });
    expect(head.status).toBe(200);
    expect(head.headers.get('content-length')).not.toBe('0');
    const post = await get('/api/visit', { method: 'POST', body: 'release_id=rel-1' });
    expect(post.status).toBe(405);
    expect(post.headers.get('cross-origin-opener-policy')).toBe('same-origin');
    const logged = requests.find((request) => request.method === 'POST');
    expect(logged).toMatchObject({ path: '/api/visit', status: 405 });
    expect(new TextDecoder().decode(logged!.body)).toBe('release_id=rel-1');
  });

  it('reads the site on every request, so a test can change it while the server runs', async () => {
    const before = await get('/releases.json');
    write('releases.json', '{"format":1,"latest":"rel-2"}\n');
    write('r/rel-2/release-manifest.json', '{}\n');
    const after = await get('/releases.json', { headers: { 'If-None-Match': before.headers.get('etag')! } });
    expect(after.status).toBe(200);
    expect(await after.text()).toContain('rel-2');
    expect((await get('/r/rel-2/release-manifest.json')).status).toBe(200);
    write('_headers', `${HEADERS}\n/releases.json\n  Cache-Control: no-cache\n`);
    expect((await get('/releases.json')).headers.get('cache-control')).toBe('no-cache');
  });

  it('answers 404 in plain text without a 404.html: never a single-page fallback', async () => {
    rmSync(join(site, '404.html'));
    const response = await get('/search');
    expect(response.status).toBe(404);
    expect(await response.text()).toBe('Not found\n');
  });
});
