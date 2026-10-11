// The site that the build writes (build/site.ts, build/site-files.ts): the entry taken from
// Vite's index.html, the release manifest that lists every file of a release, the check of a
// release directory against it, and the site index with the shell.
import { createHash } from 'node:crypto';
import { mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { entryOf } from '../../build/site';
import { analyticsEndpoint, listFiles, verifyRelease, writeReleaseManifest, writeSiteIndex } from '../../build/site-files';

const BASE = '/r/rel-1/';
const VITE_HTML = `<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <title>LOSAT Web</title>
    <script type="module" crossorigin src="/r/rel-1/assets/index-A.js"></script>
    <link rel="modulepreload" crossorigin href="/r/rel-1/assets/vendor-C.js">
    <link rel="stylesheet" crossorigin href="/r/rel-1/assets/index-B.css">
  </head>
  <body><div id="app"></div></body>
</html>`;

let dirs: string[] = [];
function tempDir(): string {
  const dir = mkdtempSync(join(tmpdir(), 'losat-site-'));
  dirs.push(dir);
  return dir;
}
function write(dir: string, path: string, text: string): void {
  mkdirSync(join(dir, path, '..'), { recursive: true });
  writeFileSync(join(dir, path), text);
}
afterEach(() => {
  for (const dir of dirs) rmSync(dir, { recursive: true, force: true });
  dirs = [];
});

describe('the entry of a release', () => {
  it('is the module script and the style sheets of the index.html that Vite built', () => {
    expect(entryOf(VITE_HTML, BASE)).toEqual({ script: 'assets/index-A.js', styles: ['assets/index-B.css'] });
  });

  it('stops the build for anything the boot page would not load', () => {
    expect(() => entryOf(VITE_HTML.replace('</head>', '<script>alert(1)</script></head>'), BASE)).toThrow(/not a module script file/);
    expect(() => entryOf(VITE_HTML.replace('</head>', '<style>p{}</style></head>'), BASE)).toThrow(/inline style/);
    expect(() => entryOf(VITE_HTML.replace('</head>', '<link rel="icon" href="/r/rel-1/x.png"></head>'), BASE)).toThrow(/rel="icon"/);
    expect(() => entryOf(VITE_HTML.replace('/r/rel-1/assets/index-A.js', '/assets/index-A.js'), BASE)).toThrow(/not under/);
    expect(() => entryOf(VITE_HTML.replace(/<script[^>]*><\/script>/, ''), BASE)).toThrow(/0 module scripts/);
  });
});

describe('the release manifest', () => {
  const release = {
    release: 'rel-1',
    app: { version: '0.0.0', build: 'abc1234' },
    engine: null,
    entry: { script: 'assets/index-A.js', styles: ['assets/index-B.css'] },
  };

  function releaseDir(): string {
    const dir = tempDir();
    write(dir, 'assets/index-A.js', 'console.log(1);\n');
    write(dir, 'assets/index-B.css', 'p{}\n');
    write(dir, 'engine/losat-web-serial-0.wasm', '\0asm');
    return dir;
  }

  it('lists every file of the release directory but itself, sorted, with size and SHA-256', () => {
    const dir = releaseDir();
    const manifest = writeReleaseManifest(dir, release);
    expect(manifest.assets.map((file) => file.path)).toEqual(['assets/index-A.js', 'assets/index-B.css', 'engine/losat-web-serial-0.wasm']);
    expect(manifest.assets[1]).toEqual({ path: 'assets/index-B.css', size: 4, sha256: createHash('sha256').update('p{}\n').digest('hex') });
    expect(JSON.parse(readFileSync(join(dir, 'release-manifest.json'), 'utf8'))).toEqual(manifest);
    // The same files give the same bytes.
    const again = readFileSync(join(dir, 'release-manifest.json'));
    writeReleaseManifest(dir, release);
    expect(readFileSync(join(dir, 'release-manifest.json'))).toEqual(again);
  });

  it('refuses an entry that the release does not have', () => {
    expect(() => writeReleaseManifest(releaseDir(), { ...release, entry: { script: 'assets/missing.js', styles: [] } })).toThrow(/not one of the assets/);
  });

  it('checks a release directory against its manifest', () => {
    const dir = releaseDir();
    writeReleaseManifest(dir, release);
    expect(verifyRelease(dir, 'rel-1').release).toBe('rel-1');
    expect(() => verifyRelease(dir, 'rel-2')).toThrow(/not rel-2's/);

    write(dir, 'assets/index-B.css', 'p{color:red}\n');
    expect(() => verifyRelease(dir, 'rel-1')).toThrow(/assets\/index-B.css is not the file its manifest lists/);
    write(dir, 'assets/index-B.css', 'p{}\n');
    write(dir, 'assets/extra.js', '');
    expect(() => verifyRelease(dir, 'rel-1')).toThrow(/assets\/extra.js is not in its manifest/);
    rmSync(join(dir, 'assets/extra.js'));
    rmSync(join(dir, 'engine/losat-web-serial-0.wasm'));
    expect(() => verifyRelease(dir, 'rel-1')).toThrow(/engine\/losat-web-serial-0.wasm is missing/);
  });
});

describe('the site index', () => {
  it('lists the releases, the latest and the shell of the site directory', () => {
    const site = tempDir();
    write(site, 'index.html', '<!doctype html>\n');
    write(site, 'boot.js', 'boot();\n');
    write(site, 'boot.css', '#boot{}\n');
    write(site, '404.html', 'Not found\n');
    const index = writeSiteIndex(site, { releases: [{ id: 'rel-1', version: '0.0.0', build: 'abc1234' }], latest: 'rel-1' });
    const files = new Map(listFiles(site).map((file) => [file.path, file]));
    expect(index.shell).toEqual([
      { ...files.get('index.html')!, path: '/' },
      { ...files.get('boot.css')!, path: '/boot.css' },
      { ...files.get('boot.js')!, path: '/boot.js' },
    ]);
    expect(index.analytics).toBeUndefined();
    expect(JSON.parse(readFileSync(join(site, 'releases.json'), 'utf8'))).toEqual(index);
  });

  it('names the measurement endpoint only when LOSAT_WEB_ANALYTICS_ENDPOINT is a path under /api/', () => {
    const site = tempDir();
    const index = writeSiteIndex(site, { releases: [{ id: 'rel-1', version: '0', build: 'x' }], latest: 'rel-1', analyticsEndpoint: '/api/visit' });
    expect(index.analytics).toEqual({ endpoint: '/api/visit' });
    expect(analyticsEndpoint({})).toBeUndefined();
    expect(analyticsEndpoint({ LOSAT_WEB_ANALYTICS_ENDPOINT: '/api/visit' })).toBe('/api/visit');
    for (const value of ['https://example.org/api/visit', '/collect', '/api/../x']) {
      expect(() => analyticsEndpoint({ LOSAT_WEB_ANALYTICS_ENDPOINT: value })).toThrow(/LOSAT_WEB_ANALYTICS_ENDPOINT/);
    }
  });
});
