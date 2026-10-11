// What the service worker answers (src/sw/serve.ts) and which release caches it removes
// (src/sw/cleanup.ts).
import { describe, expect, it } from 'vitest';
import { releaseCachesToRemove } from '../../src/sw/cleanup';
import type { CacheLike, CacheStorageLike } from '../../src/sw/prepare';
import { routeOf, serveRelease, serveShell, type ServeDeps } from '../../src/sw/serve';

const ORIGIN = 'http://site.test';

class MapCache implements CacheLike {
  readonly entries = new Map<string, string>();
  async match(key: string): Promise<Response | undefined> {
    const body = this.entries.get(key);
    return body === undefined ? undefined : new Response(body, { headers: { 'X-From': 'cache' } });
  }
  async put(key: string, response: Response): Promise<void> {
    this.entries.set(key, await response.text());
  }
  async delete(key: string): Promise<boolean> {
    return this.entries.delete(key);
  }
  async keys(): Promise<Array<{ url: string }>> {
    return [...this.entries.keys()].map((url) => ({ url }));
  }
}

class MapCaches implements CacheStorageLike {
  readonly stores = new Map<string, MapCache>();
  async has(name: string): Promise<boolean> {
    return this.stores.has(name);
  }
  async open(name: string): Promise<MapCache> {
    if (!this.stores.has(name)) this.stores.set(name, new MapCache());
    return this.stores.get(name)!;
  }
  async delete(name: string): Promise<boolean> {
    return this.stores.delete(name);
  }
  async keys(): Promise<string[]> {
    return [...this.stores.keys()];
  }
}

describe('routes', () => {
  const route = (url: string, method = 'GET') => routeOf({ method, url }, ORIGIN);

  it('answers the files of a release, the start page and the shell files', () => {
    expect(route(`${ORIGIN}/r/a1-fake/assets/index-AbC.js`)).toEqual({ kind: 'release', release: 'a1-fake' });
    expect(route(`${ORIGIN}/r/a1-fake/release-manifest.json`)).toEqual({ kind: 'release', release: 'a1-fake' });
    expect(route(`${ORIGIN}/`)).toEqual({ kind: 'shell', path: '/' });
    expect(route(`${ORIGIN}/?from=bookmark`)).toEqual({ kind: 'shell', path: '/' });
    expect(route(`${ORIGIN}/boot.js`)).toEqual({ kind: 'shell', path: '/boot.js' });
    expect(route(`${ORIGIN}/boot.css`)).toEqual({ kind: 'shell', path: '/boot.css' });
    expect(route(`${ORIGIN}/404.html`)).toEqual({ kind: 'shell', path: '/404.html' });
  });

  it.each([
    [`${ORIGIN}/releases.json`, 'GET'],
    [`${ORIGIN}/api/visit`, 'GET'],
    [`${ORIGIN}/api/visit`, 'POST'],
    [`${ORIGIN}/sw.js`, 'GET'],
    [`${ORIGIN}/r/a1-fake/assets/index.js`, 'POST'],
    [`${ORIGIN}/r/Not_An_Id/assets/index.js`, 'GET'],
    [`${ORIGIN}/r/a1-fake/`, 'GET'],
    [`${ORIGIN}/somewhere/else.js`, 'GET'],
    ['http://other.test/boot.js', 'GET'],
    ['https://www.google-analytics.com/g/collect', 'GET'],
    ['blob:http://site.test/5f0c-4e1d', 'GET'],
    ['blob:http://site.test/r/a1-fake/x', 'GET'],
    ['data:text/plain,hello', 'GET'],
  ])('does not answer %s (%s)', (url, method) => {
    expect(route(url, method)).toEqual({ kind: 'none' });
  });
});

describe('serving', () => {
  function setup(network: (request: Request) => Promise<Response>) {
    const caches = new MapCaches();
    const requests: string[] = [];
    const deps: ServeDeps = {
      origin: ORIGIN,
      caches,
      fetch: (request) => {
        requests.push(request.url);
        return network(request);
      },
      shellTimeoutMs: 30,
    };
    return { caches, requests, deps };
  }
  const fromNetwork = (status = 200) => () => Promise.resolve(new Response('network', { status }));

  it('answers a file of a complete release from its cache, without the network', async () => {
    const { caches, requests, deps } = setup(fromNetwork());
    const cache = await caches.open('losat-web-release-a');
    await cache.put(`${ORIGIN}/r/a/assets/x.js`, new Response('stored'));
    await cache.put(`${ORIGIN}/r/a/release-manifest.json`, new Response('{}'));
    const response = await serveRelease(deps, new Request(`${ORIGIN}/r/a/assets/x.js`), 'a');
    expect(await response.text()).toBe('stored');
    expect(requests).toEqual([]);
  });

  it('answers from the network, storing nothing, while the release is incomplete or lacks the file', async () => {
    const { caches, requests, deps } = setup(fromNetwork());
    const cache = await caches.open('losat-web-release-a');
    await cache.put(`${ORIGIN}/r/a/assets/x.js`, new Response('stored'));
    expect(await (await serveRelease(deps, new Request(`${ORIGIN}/r/a/assets/x.js`), 'a')).text()).toBe('network');
    await cache.put(`${ORIGIN}/r/a/release-manifest.json`, new Response('{}'));
    expect(await (await serveRelease(deps, new Request(`${ORIGIN}/r/a/assets/y.js`), 'a')).text()).toBe('network');
    expect(await (await serveRelease(deps, new Request(`${ORIGIN}/r/b/assets/x.js`), 'b')).text()).toBe('network');
    expect(requests).toHaveLength(3);
    expect([...cache.entries.keys()].sort()).toEqual([`${ORIGIN}/r/a/assets/x.js`, `${ORIGIN}/r/a/release-manifest.json`]);
    expect(await caches.has('losat-web-release-b')).toBe(false);
  });

  async function withShell(network: (request: Request) => Promise<Response>) {
    const s = setup(network);
    await (await s.caches.open('losat-web-shell')).put(`${ORIGIN}/`, new Response('stored boot page'));
    return s;
  }

  it('answers the start page from the network while it answers', async () => {
    const { deps } = await withShell(fromNetwork());
    expect(await (await serveShell(deps, new Request(`${ORIGIN}/`), '/')).text()).toBe('network');
  });

  it.each([
    ['fails', () => Promise.reject(new TypeError('Failed to fetch'))],
    ['answers a server error', fromNetwork(503)],
    ['takes too long', () => new Promise<Response>(() => {})],
  ])('answers the stored start page when the network %s', async (_, network) => {
    const { deps } = await withShell(network);
    expect(await (await serveShell(deps, new Request(`${ORIGIN}/`), '/')).text()).toBe('stored boot page');
  });

  it('passes a 404 of the network through: only a failure falls back', async () => {
    const { deps } = await withShell(fromNetwork(404));
    expect((await serveShell(deps, new Request(`${ORIGIN}/`), '/')).status).toBe(404);
  });

  it('waits for the network when nothing is stored', async () => {
    const { deps } = setup(() => new Promise((resolve) => setTimeout(() => resolve(new Response('slow network')), 80)));
    expect(await (await serveShell(deps, new Request(`${ORIGIN}/boot.js`), '/boot.js')).text()).toBe('slow network');
    const failing = setup(() => Promise.reject(new TypeError('Failed to fetch')));
    await expect(serveShell(failing.deps, new Request(`${ORIGIN}/`), '/')).rejects.toThrow('Failed to fetch');
  });
});

describe('cleanup', () => {
  const names = ['losat-web-release-a', 'losat-web-release-b', 'losat-web-release-c', 'losat-web-shell', 'other-cache'];
  const held = (...releases: string[]) => releases.map((release) => ({ name: `losat-web:release:${release}` }));

  it('removes only release caches that no tab holds or waits for, and that the page does not keep', () => {
    expect(releaseCachesToRemove(names, { held: held('a'), pending: [] }, ['b'], [])).toEqual(['losat-web-release-c']);
    expect(releaseCachesToRemove(names, { held: [], pending: held('c') }, ['b'], [])).toEqual(['losat-web-release-a']);
    expect(releaseCachesToRemove(names, { held: held('a', 'b', 'c') }, [], [])).toEqual([]);
    expect(releaseCachesToRemove(names, {}, ['a'], ['b'])).toEqual(['losat-web-release-c']);
  });

  it('ignores the tmp/ locks and other locks', () => {
    expect(releaseCachesToRemove(names, { held: [{ name: 'losat-web:tmp:c' }, { name: 'c' }, {}] }, ['a', 'b'], [])).toEqual([
      'losat-web-release-c',
    ]);
  });

  it('removes nothing without the Web Locks API', () => {
    expect(releaseCachesToRemove(names, undefined, [], [])).toEqual([]);
  });
});
