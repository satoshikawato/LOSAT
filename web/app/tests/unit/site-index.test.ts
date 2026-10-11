import { describe, expect, it } from 'vitest';
import { parseSiteIndex } from '../../src/domain/site-index';

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
    expect(parseSiteIndex(index({ shell: [] })).latest).toBe('b');
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
