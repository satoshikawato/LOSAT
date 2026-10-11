// The service worker's preparation of a release (src/sw/prepare.ts): every listed file is
// checked before anything is stored, and a failure leaves the caches as they were.
import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import { isComplete, prepareRelease, type CacheLike, type CacheStorageLike, type ListedFile, type PrepareDeps } from '../../src/sw/prepare';
import { SHELL_CACHE, type PrepareProgress } from '../../src/sw/protocol';

const ORIGIN = 'http://site.test';
const SITE_HEADERS = {
  'Cross-Origin-Opener-Policy': 'same-origin',
  'Cross-Origin-Embedder-Policy': 'require-corp',
  'Cross-Origin-Resource-Policy': 'same-origin',
  'Content-Security-Policy': "default-src 'self'",
};

class FakeCache implements CacheLike {
  readonly entries = new Map<string, { body: Uint8Array; status: number; headers: [string, string][] }>();
  failPut: ((key: string) => boolean) | undefined;
  async match(key: string): Promise<Response | undefined> {
    const entry = this.entries.get(key);
    return entry === undefined ? undefined : new Response(entry.body.slice(), { status: entry.status, headers: entry.headers });
  }
  async put(key: string, response: Response): Promise<void> {
    if (this.failPut?.(key) === true) throw new Error('quota exceeded');
    const body = new Uint8Array(await response.arrayBuffer());
    this.entries.set(key, { body, status: response.status, headers: [...response.headers] });
  }
  async delete(key: string): Promise<boolean> {
    return this.entries.delete(key);
  }
  async keys(): Promise<Array<{ url: string }>> {
    return [...this.entries.keys()].map((url) => ({ url }));
  }
}

class FakeCaches implements CacheStorageLike {
  readonly stores = new Map<string, FakeCache>();
  failPut: ((key: string) => boolean) | undefined;
  async has(name: string): Promise<boolean> {
    return this.stores.has(name);
  }
  async open(name: string): Promise<FakeCache> {
    let cache = this.stores.get(name);
    if (cache === undefined) {
      cache = new FakeCache();
      this.stores.set(name, cache);
    }
    cache.failPut = (key) => this.failPut?.(key) === true;
    return cache;
  }
  async delete(name: string): Promise<boolean> {
    return this.stores.delete(name);
  }
  async keys(): Promise<string[]> {
    return [...this.stores.keys()];
  }
  /** Every cache, key and body, to compare before and after. */
  snapshot(): Record<string, Record<string, string>> {
    return Object.fromEntries(
      [...this.stores].map(([name, cache]) => [
        name,
        Object.fromEntries([...cache.entries].map(([key, entry]) => [key, Buffer.from(entry.body).toString('hex')])),
      ]),
    );
  }
}

const sha256 = (bytes: Uint8Array | ArrayBuffer): string =>
  createHash('sha256').update(new Uint8Array(bytes)).digest('hex');
const bytesOf = (text: string): Uint8Array => new TextEncoder().encode(text);
const listed = (path: string, body: Uint8Array): ListedFile => ({ path, size: body.length, sha256: sha256(body) });

interface Served {
  body: Uint8Array;
  status?: number;
  redirected?: boolean;
  /** The body breaks after this many bytes (the connection is cut). */
  cutAfter?: number;
  /** The request never ends. */
  hang?: boolean;
}

/** A site with release `a` (two assets) and its shell; `serve` changes what a path answers. */
function site(release = 'a') {
  const shellFiles: Record<string, Uint8Array> = { '/': bytesOf('<!doctype html>boot'), '/boot.js': bytesOf('boot();') };
  const assets: Record<string, Uint8Array> = {
    'assets/index-1.js': bytesOf('console.log("app");'),
    'engine/serial.wasm': new Uint8Array(4096).map((_, n) => n % 251),
  };
  const index = {
    format: 1,
    latest: release,
    releases: [{ id: release, version: '0.0.0', build: 'abc' }],
    shell: Object.entries(shellFiles).map(([path, body]) => listed(path, body)),
  };
  const manifest = {
    format: 1,
    release,
    assets: Object.entries(assets).map(([path, body]) => listed(path, body)),
  };
  const served = new Map<string, Served>();
  for (const [path, body] of Object.entries(shellFiles)) served.set(path, { body });
  for (const [path, body] of Object.entries(assets)) served.set(`/r/${release}/${path}`, { body });
  served.set('/releases.json', { body: bytesOf(JSON.stringify(index)) });
  served.set(`/r/${release}/release-manifest.json`, { body: bytesOf(JSON.stringify(manifest)) });
  const requests: Array<{ path: string; cache: RequestCache | undefined }> = [];
  const fetch = async (url: string, init: RequestInit): Promise<Response> => {
    const path = new URL(url, ORIGIN).pathname;
    requests.push({ path, cache: init.cache });
    const entry = served.get(path);
    if (entry === undefined) return new Response('not found', { status: 404, headers: SITE_HEADERS });
    if (entry.hang === true) {
      return new Promise((_, reject) => init.signal?.addEventListener('abort', () => reject(init.signal?.reason)));
    }
    let body: BodyInit = entry.body.slice();
    if (entry.cutAfter !== undefined) {
      const head = entry.body.slice(0, entry.cutAfter);
      body = new ReadableStream<Uint8Array>({
        start(controller) {
          controller.enqueue(head);
          controller.error(new TypeError('network error: connection reset'));
        },
      });
    }
    const response = new Response(body, { status: entry.status ?? 200, headers: SITE_HEADERS });
    if (entry.redirected === true) Object.defineProperty(response, 'redirected', { value: true });
    return response;
  };
  return { shellFiles, assets, index, manifest, served, requests, fetch };
}

function deps(s: ReturnType<typeof site>, caches: FakeCaches, extra: Partial<PrepareDeps> = {}): PrepareDeps {
  return {
    origin: ORIGIN,
    fetch: s.fetch,
    caches,
    sha256: async (bytes) => sha256(bytes),
    readShell: (index) => (index as { shell: ListedFile[] }).shell,
    readManifest: (manifest, release) => {
      const m = manifest as { release: string; assets: ListedFile[] };
      if (m.release !== release) throw new Error('the release manifest names another release');
      return m.assets;
    },
    ...extra,
  };
}

async function prepare(s: ReturnType<typeof site>, caches: FakeCaches, release = 'a', extra: Partial<PrepareDeps> = {}) {
  const progress: PrepareProgress[] = [];
  const outcome = await prepareRelease(deps(s, caches, extra), release, (p) => progress.push(p));
  return { outcome, progress };
}

describe('preparing a release', () => {
  it('stores every listed file and the shell after checking them, the manifest last, with status and headers', async () => {
    const s = site();
    const caches = new FakeCaches();
    const { outcome, progress } = await prepare(s, caches);
    expect(outcome).toEqual({ ok: true, progress: { files: 4, filesTotal: 4, bytes: expect.any(Number), bytesTotal: expect.any(Number) } });
    const release = caches.stores.get('losat-web-release-a')!;
    expect([...release.entries.keys()].sort()).toEqual([
      `${ORIGIN}/r/a/assets/index-1.js`,
      `${ORIGIN}/r/a/engine/serial.wasm`,
      `${ORIGIN}/r/a/release-manifest.json`,
    ]);
    // The manifest is stored last: it marks the cache complete.
    expect([...release.entries.keys()].at(-1)).toBe(`${ORIGIN}/r/a/release-manifest.json`);
    expect([...caches.stores.get(SHELL_CACHE)!.entries.keys()].sort()).toEqual([`${ORIGIN}/`, `${ORIGIN}/boot.js`]);
    expect([...caches.stores.keys()].sort()).toEqual(['losat-web-release-a', SHELL_CACHE]);
    const wasm = await release.match(`${ORIGIN}/r/a/engine/serial.wasm`);
    expect(wasm!.status).toBe(200);
    for (const [name, value] of Object.entries(SITE_HEADERS)) expect(wasm!.headers.get(name)).toBe(value);
    expect(sha256(await wasm!.arrayBuffer())).toBe(sha256(s.assets['engine/serial.wasm']!));
    expect(await isComplete(caches, ORIGIN, 'a')).toBe(true);
    // Progress: files and bytes done of total, from none to all.
    const total = Object.values(s.shellFiles).concat(Object.values(s.assets)).reduce((sum, body) => sum + body.length, 0);
    expect(progress[0]).toEqual({ files: 0, filesTotal: 4, bytes: 0, bytesTotal: total });
    expect(progress.at(-1)).toEqual({ files: 4, filesTotal: 4, bytes: total, bytesTotal: total });
    expect(progress.map((p) => p.files)).toEqual([0, 1, 2, 3, 4]);
    // The site index is read from the network only; the shell is revalidated.
    expect(s.requests.find((r) => r.path === '/releases.json')?.cache).toBe('no-store');
    expect(s.requests.find((r) => r.path === '/boot.js')?.cache).toBe('no-cache');
  });

  const failures: Array<[string, (s: ReturnType<typeof site>) => void, RegExp]> = [
    [
      'a file with other bytes of the same size (wrong hash)',
      (s) => s.served.set('/r/a/engine/serial.wasm', { body: new Uint8Array(4096) }),
      /engine\/serial\.wasm: its SHA-256 is not the one that the list gives/,
    ],
    ['a file of another size', (s) => s.served.set('/r/a/assets/index-1.js', { body: bytesOf('short') }), /assets\/index-1\.js: 5 bytes, not the 19/],
    [
      'a redirected response',
      (s) => s.served.set('/r/a/assets/index-1.js', { body: s.assets['assets/index-1.js']!, redirected: true }),
      /redirected/,
    ],
    ['a missing file (404)', (s) => s.served.delete('/r/a/engine/serial.wasm'), /engine\/serial\.wasm: the server answered 404/],
    [
      'a download cut in the middle',
      (s) => s.served.set('/r/a/engine/serial.wasm', { body: s.assets['engine/serial.wasm']!, cutAfter: 100 }),
      /engine\/serial\.wasm: the download failed \(network error: connection reset\)/,
    ],
    ['a shell file that does not match', (s) => s.served.set('/boot.js', { body: bytesOf('other();') }), /\/boot\.js: 8 bytes, not the 7/],
    ['a missing release manifest', (s) => s.served.delete('/r/a/release-manifest.json'), /the release manifest: the server answered 404/],
    [
      'a manifest of another release',
      (s) => s.served.set('/r/a/release-manifest.json', { body: bytesOf(JSON.stringify({ ...s.manifest, release: 'b' })) }),
      /names another release/,
    ],
    ['an unreadable site index', (s) => s.served.set('/releases.json', { body: bytesOf('<html>') }), /the site index is not valid/],
  ];

  it.each(failures)('stores nothing after %s, and says why', async (_, change, reason) => {
    const s = site();
    change(s);
    const caches = new FakeCaches();
    // A release and a shell stored before stay as they were.
    const before = site('old');
    await prepare(before, caches, 'old');
    const snapshot = caches.snapshot();
    const { outcome } = await prepare(s, caches);
    expect(outcome.ok).toBe(false);
    expect(outcome.ok ? '' : outcome.reason).toMatch(reason);
    expect(caches.snapshot()).toEqual(snapshot);
    expect(await isComplete(caches, ORIGIN, 'a')).toBe(false);
  });

  it('tries a mismatching file once more past the HTTP cache before it fails', async () => {
    const s = site();
    s.served.set('/r/a/engine/serial.wasm', { body: new Uint8Array(4096) });
    await prepare(s, new FakeCaches());
    expect(s.requests.filter((r) => r.path === '/r/a/engine/serial.wasm').map((r) => r.cache)).toEqual(['default', 'reload']);
  });

  it('fails when no file arrives for a while', async () => {
    const s = site();
    s.served.set('/r/a/engine/serial.wasm', { body: new Uint8Array(1), hang: true });
    const caches = new FakeCaches();
    const { outcome } = await prepare(s, caches, 'a', { stallMs: 50 });
    expect(outcome).toEqual({ ok: false, reason: expect.stringMatching(/engine\/serial\.wasm: the download failed \(no file arrived for a while\)/) });
    expect(caches.snapshot()).toEqual({});
  });

  it('puts the caches back when storing fails', async () => {
    const s = site();
    const caches = new FakeCaches();
    await prepare(site('old'), caches, 'old');
    const snapshot = caches.snapshot();
    caches.failPut = (key) => key.endsWith('/release-manifest.json');
    const { outcome } = await prepare(s, caches);
    expect(outcome).toEqual({ ok: false, reason: 'the files could not be stored (quota exceeded)' });
    expect(caches.snapshot()).toEqual(snapshot);
    // A failure while the shell is stored puts the old shell back and removes the new release.
    caches.failPut = (key) => key === `${ORIGIN}/boot.js`;
    s.served.set('/', { body: bytesOf('<!doctype html>new boot page') });
    s.index.shell[0] = listed('/', bytesOf('<!doctype html>new boot page'));
    s.served.set('/releases.json', { body: bytesOf(JSON.stringify(s.index)) });
    const second = await prepare(s, caches);
    expect(second.outcome).toEqual({ ok: false, reason: 'the start-up files could not be stored (quota exceeded)' });
    expect(caches.snapshot()).toEqual(snapshot);
  });

  it('starts again over a cache that an interrupted preparation left without its manifest', async () => {
    const s = site();
    const caches = new FakeCaches();
    const partial = await caches.open('losat-web-release-a');
    await partial.put(`${ORIGIN}/r/a/assets/stale.js`, new Response('stale'));
    expect(await isComplete(caches, ORIGIN, 'a')).toBe(false);
    const { outcome } = await prepare(s, caches);
    expect(outcome.ok).toBe(true);
    expect([...caches.stores.get('losat-web-release-a')!.entries.keys()]).not.toContain(`${ORIGIN}/r/a/assets/stale.js`);
  });

  it('does not download a complete release again, but reads the shell again', async () => {
    const s = site();
    const caches = new FakeCaches();
    await prepare(s, caches);
    s.requests.length = 0;
    const newBoot = bytesOf('boot(2);');
    s.served.set('/boot.js', { body: newBoot });
    s.index.shell[1] = listed('/boot.js', newBoot);
    s.served.set('/releases.json', { body: bytesOf(JSON.stringify(s.index)) });
    const { outcome, progress } = await prepare(s, caches);
    expect(outcome.ok).toBe(true);
    expect(progress.at(-1)).toMatchObject({ files: 2, filesTotal: 2 });
    expect(s.requests.map((r) => r.path).sort()).toEqual(['/', '/boot.js', '/releases.json']);
    const stored = await caches.stores.get(SHELL_CACHE)!.match(`${ORIGIN}/boot.js`);
    expect(await stored!.text()).toBe('boot(2);');
  });

  it('keeps a complete release and its stored shell when the network is gone', async () => {
    const s = site();
    const caches = new FakeCaches();
    await prepare(s, caches);
    const snapshot = caches.snapshot();
    const offline = await prepare(s, caches, 'a', { fetch: () => Promise.reject(new TypeError('Failed to fetch')) });
    expect(offline.outcome.ok).toBe(true);
    expect(caches.snapshot()).toEqual(snapshot);
    // Without a stored release, the same failure is a failure.
    const fresh = new FakeCaches();
    const first = await prepare(s, fresh, 'a', { fetch: () => Promise.reject(new TypeError('Failed to fetch')) });
    expect(first.outcome).toEqual({ ok: false, reason: 'the site index could not be read (Failed to fetch)' });
    expect(fresh.snapshot()).toEqual({});
  });

  it('refuses a site index without a start page', async () => {
    const s = site();
    s.index.shell = s.index.shell.filter((file) => file.path !== '/');
    s.served.set('/releases.json', { body: bytesOf(JSON.stringify(s.index)) });
    const caches = new FakeCaches();
    expect((await prepare(s, caches)).outcome).toEqual({ ok: false, reason: 'the site index lists no start page' });
    expect(caches.snapshot()).toEqual({});
  });

  it('removes shell entries that the site no longer lists', async () => {
    const s = site();
    const caches = new FakeCaches();
    await prepare(s, caches);
    s.index.shell = s.index.shell.filter((file) => file.path !== '/boot.js');
    s.served.set('/releases.json', { body: bytesOf(JSON.stringify(s.index)) });
    await prepare(s, caches);
    expect([...caches.stores.get(SHELL_CACHE)!.entries.keys()]).toEqual([`${ORIGIN}/`]);
  });
});
