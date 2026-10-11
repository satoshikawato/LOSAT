import { describe, expect, it } from 'vitest';
import { parseSiteIndex } from '../../src/domain/site-index';
import { sharedSiteReads } from '../../src/infra/browser/site';

const index = (extra: Record<string, unknown> = {}) => ({
  format: 1,
  latest: 'b',
  releases: [
    { id: 'a', version: '0.0.0', build: '08e4557b' },
    { id: 'b', version: '0.0.1', build: '55c924d8' },
  ],
  ...extra,
});

describe('site index', () => {
  it('reads the releases, the latest and the measurement endpoint', () => {
    const parsed = parseSiteIndex(index({ analytics: { endpoint: '/api/visit' } }));
    expect(parsed.latest).toBe('b');
    expect(parsed.releases.map((release) => release.id)).toEqual(['a', 'b']);
    expect(parsed.analytics).toEqual({ endpoint: '/api/visit' });
  });

  it('has no measurement unless the site names an endpoint', () => {
    expect(parseSiteIndex(index()).analytics).toBeUndefined();
  });

  it('ignores fields it does not know', () => {
    expect(parseSiteIndex(index({ notices: [] })).latest).toBe('b');
  });

  it('reads the shell: the site files a browser needs to start offline', () => {
    const sha256 = 'c'.repeat(64);
    const shell = [
      { path: '/', size: 451, sha256 },
      { path: '/boot.css', size: 544, sha256 },
      { path: '/boot.js', size: 7202, sha256 },
    ];
    expect(parseSiteIndex(index({ shell })).shell).toEqual(shell);
    expect(parseSiteIndex(index()).shell).toBeUndefined();
    for (const path of ['boot.js', '//boot.js', '/../boot.js', '/r/a/../b.js', 'https://example.org/boot.js']) {
      expect(() => parseSiteIndex(index({ shell: [{ path, size: 1, sha256 }] }))).toThrow(/path/);
    }
    expect(() => parseSiteIndex(index({ shell: [shell[1], shell[0]] }))).toThrow(/order/);
    expect(() => parseSiteIndex(index({ shell: [{ path: '/', size: -1, sha256 }] }))).toThrow(/size/);
  });

  it('refuses an index whose latest is not served, an unknown format and bad IDs', () => {
    expect(() => parseSiteIndex(index({ latest: 'c' }))).toThrow(/latest/);
    expect(() => parseSiteIndex(index({ format: 2 }))).toThrow(/format/);
    expect(() => parseSiteIndex(index({ releases: [{ id: '../x', version: '0', build: 'x' }], latest: '../x' }))).toThrow(/id/);
  });

  it('accepts only a same-origin path under /api/ as the endpoint', () => {
    for (const endpoint of ['https://example.org/api/visit', '//example.org/api/visit', '/collect', '/api/../x', '/api/Visit?x=1']) {
      expect(() => parseSiteIndex(index({ analytics: { endpoint } }))).toThrow(/endpoint/);
    }
  });
});

describe('shared site index reads', () => {
  it('shares a read with the calls made within the age, then reads again', async () => {
    let reads = 0;
    let clock = 0;
    const port = { index: async () => (reads++, parseSiteIndex(index())) };
    const shared = sharedSiteReads(port, 60_000, () => clock);
    await Promise.all([shared.index(), shared.index()]);
    clock = 59_000;
    await shared.index();
    expect(reads).toBe(1);
    clock = 120_000;
    await shared.index();
    expect(reads).toBe(2);
  });
});
