import { describe, expect, it } from 'vitest';
import { parseReleaseManifest, releaseDirectory, releaseFileUrl, releaseManifestUrl } from '../../src/domain/release-manifest';

const SHA = 'a'.repeat(64);
const manifest = (extra: Record<string, unknown> = {}) => ({
  format: 1,
  release: 'abc123-fake',
  app: { version: '0.0.0', build: 'abc1234' },
  engine: null,
  entry: { script: 'assets/index-1.js', styles: ['assets/index-2.css'] },
  assets: [
    { path: 'assets/data-worker-3.js', size: 10, sha256: SHA },
    { path: 'assets/index-1.js', size: 20, sha256: SHA },
    { path: 'assets/index-2.css', size: 0, sha256: SHA },
    { path: 'engine/losat-web-serial-0123.wasm', size: 30, sha256: SHA },
  ],
  ...extra,
});

describe('release manifest', () => {
  it('reads the release, the entry and the assets', () => {
    const parsed = parseReleaseManifest(manifest(), 'abc123-fake');
    expect(parsed.release).toBe('abc123-fake');
    expect(parsed.engine).toBeNull();
    expect(parsed.entry).toEqual({ script: 'assets/index-1.js', styles: ['assets/index-2.css'] });
    expect(parsed.assets.map((file) => file.path)).toHaveLength(4);
  });

  it('reads the engine build names of a build with the engine', () => {
    const engine = { serial: 'losat-web-serial.wasm sha256:0123', threads: 'losat-web-threads.wasm sha256:4567 (shared-memory guard)' };
    expect(parseReleaseManifest(manifest({ engine })).engine).toEqual(engine);
  });

  it("ignores fields it does not know, and refuses another release's manifest or format", () => {
    expect(parseReleaseManifest(manifest({ signature: 'x' })).release).toBe('abc123-fake');
    expect(() => parseReleaseManifest(manifest(), 'other')).toThrow(/not other's/);
    expect(() => parseReleaseManifest(manifest({ format: 2 }))).toThrow(/format/);
    expect(() => parseReleaseManifest(manifest({ release: '../x' }))).toThrow(/release/);
  });

  it('refuses an entry that is not one of the assets', () => {
    expect(() => parseReleaseManifest(manifest({ entry: { script: 'assets/other.js', styles: [] } }))).toThrow(/not one of the assets/);
    expect(() => parseReleaseManifest(manifest({ entry: { script: 'assets/index-1.js', styles: ['x.css'] } }))).toThrow(/not one of the assets/);
  });

  it('refuses paths that leave the release directory, are out of order, repeated or the manifest itself', () => {
    const withAssets = (paths: string[]) =>
      manifest({ entry: { script: paths[0], styles: [] }, assets: paths.map((path) => ({ path, size: 1, sha256: SHA })) });
    for (const path of ['../x.js', '/assets/x.js', 'assets//x.js', 'assets/./x.js', 'assets\\x.js', 'assets/x.js?v=1', '.hidden', '']) {
      expect(() => parseReleaseManifest(withAssets([path]))).toThrow(/path is not valid/);
    }
    expect(() => parseReleaseManifest(withAssets(['b.js', 'a.js']))).toThrow(/order/);
    expect(() => parseReleaseManifest(withAssets(['a.js', 'a.js']))).toThrow(/order/);
    expect(() => parseReleaseManifest(withAssets(['a.js', 'release-manifest.json']))).toThrow(/lists itself/);
  });

  it('refuses a size or digest that is not one', () => {
    const asset = (file: Record<string, unknown>) => manifest({ entry: { script: 'a.js', styles: [] }, assets: [{ path: 'a.js', size: 1, sha256: SHA, ...file }] });
    for (const size of [-1, 1.5, '1', Number.MAX_SAFE_INTEGER + 1]) expect(() => parseReleaseManifest(asset({ size }))).toThrow(/size/);
    for (const sha256 of ['A'.repeat(64), 'a'.repeat(63), 1]) expect(() => parseReleaseManifest(asset({ sha256 }))).toThrow(/sha256/);
  });

  it('names the URLs of a release', () => {
    expect(releaseDirectory('r1')).toBe('/r/r1/');
    expect(releaseFileUrl('r1', 'assets/a.js')).toBe('/r/r1/assets/a.js');
    expect(releaseManifestUrl('r1')).toBe('/r/r1/release-manifest.json');
  });
});
