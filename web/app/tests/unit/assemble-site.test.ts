// Assembling a site with older releases (build/assemble-site.ts) on small fake sites: the
// releases of earlier sites are copied after every file is checked against its manifest, and
// releases.json lists every release with the new site's latest and shell.
import { existsSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { assembleSite } from '../../build/assemble-site';
import { listFiles, readSiteIndex, siteShell, writeReleaseManifest, writeSiteIndex } from '../../build/site-files';

let dirs: string[] = [];
afterEach(() => {
  for (const dir of dirs) rmSync(dir, { recursive: true, force: true });
  dirs = [];
});

function write(dir: string, path: string, text: string): void {
  mkdirSync(join(dir, path, '..'), { recursive: true });
  writeFileSync(join(dir, path), text);
}

/** A site as the build writes it, with one release `id` (and its boot page saying which build made it). */
function fakeSite(id: string, options: { analytics?: string } = {}): string {
  const site = mkdtempSync(join(tmpdir(), `losat-assemble-${id}-`));
  dirs.push(site);
  write(site, 'index.html', `<!doctype html><!-- ${id} -->\n`);
  write(site, 'boot.js', `/* boot of ${id} */\n`);
  write(site, 'boot.css', '#boot{}\n');
  write(site, '404.html', 'Not found\n');
  write(site, '_headers', '/*\n  X: 1\n');
  const release = join(site, 'r', id);
  write(release, 'assets/index-1.js', `console.log(${JSON.stringify(id)});\n`);
  write(release, 'assets/index-2.css', 'p{}\n');
  const manifest = writeReleaseManifest(release, {
    release: id,
    app: { version: '0.0.0', build: `build-${id}` },
    engine: null,
    entry: { script: 'assets/index-1.js', styles: ['assets/index-2.css'] },
  });
  writeSiteIndex(site, {
    releases: [{ id, version: manifest.app.version, build: manifest.app.build }],
    latest: id,
    analyticsEndpoint: options.analytics,
  });
  return site;
}

const snapshot = (site: string) => listFiles(site).map((file) => `${file.path} ${file.sha256}`);

describe('assembling a site with older releases', () => {
  it("copies the earlier site's release and lists both, the new one latest, with the new site's shell", () => {
    const a = fakeSite('rel-a');
    const b = fakeSite('rel-b');
    const index = assembleSite({ site: b, previous: [a] });
    expect(index.latest).toBe('rel-b');
    expect(index.releases).toEqual([
      { id: 'rel-b', version: '0.0.0', build: 'build-rel-b' },
      { id: 'rel-a', version: '0.0.0', build: 'build-rel-a' },
    ]);
    expect(index.shell).toEqual(siteShell(b));
    expect(readFileSync(join(b, 'boot.js'), 'utf8')).toContain('rel-b');
    expect(listFiles(join(b, 'r', 'rel-a'))).toEqual(listFiles(join(a, 'r', 'rel-a')));
    expect(readSiteIndex(b)).toEqual(index);
    expect(readdirSync(b).filter((name) => name.startsWith('.'))).toEqual([]);
  });

  it('keeps every release of a chain of sites, and takes --latest', () => {
    const a = fakeSite('rel-a');
    const b = fakeSite('rel-b');
    assembleSite({ site: b, previous: [a] });
    const c = fakeSite('rel-c');
    const index = assembleSite({ site: c, previous: [b], latest: 'rel-b' });
    expect(index.releases.map((release) => release.id)).toEqual(['rel-c', 'rel-b', 'rel-a']);
    expect(index.latest).toBe('rel-b');
    expect(() => assembleSite({ site: c, previous: [], latest: 'rel-x' })).toThrow(/no such release/);
  });

  it('refuses a release whose files do not match its manifest, and changes nothing', () => {
    const a = fakeSite('rel-a');
    const b = fakeSite('rel-b');
    const before = snapshot(b);
    write(join(a, 'r', 'rel-a'), 'assets/index-2.css', 'p{color:red}\n');
    expect(() => assembleSite({ site: b, previous: [a] })).toThrow(/assets\/index-2.css is not the file its manifest lists/);
    expect(snapshot(b)).toEqual(before);

    write(join(a, 'r', 'rel-a'), 'assets/index-2.css', 'p{}\n');
    write(join(a, 'r', 'rel-a'), 'assets/extra.js', '');
    expect(() => assembleSite({ site: b, previous: [a] })).toThrow(/extra.js is not in its manifest/);
    expect(existsSync(join(b, 'r', 'rel-a'))).toBe(false);
  });

  it('refuses two different releases of the same ID, and keeps one copy of the same release', () => {
    const a = fakeSite('rel-a');
    const b = fakeSite('rel-b');
    assembleSite({ site: b, previous: [a] });
    expect(assembleSite({ site: b, previous: [a] }).releases).toHaveLength(2);
    const other = fakeSite('rel-a');
    write(join(other, 'r', 'rel-a'), 'assets/index-1.js', 'console.log("another rel-a");\n');
    writeReleaseManifest(join(other, 'r', 'rel-a'), {
      release: 'rel-a',
      app: { version: '0.0.0', build: 'other' },
      engine: null,
      entry: { script: 'assets/index-1.js', styles: [] },
    });
    expect(() => assembleSite({ site: b, previous: [other] })).toThrow(/not the release rel-a that the site already has/);
  });

  it("keeps the new site's measurement setting unless an endpoint is given", () => {
    const a = fakeSite('rel-a');
    const b = fakeSite('rel-b', { analytics: '/api/visit' });
    expect(assembleSite({ site: b, previous: [a] }).analytics).toEqual({ endpoint: '/api/visit' });
    expect(assembleSite({ site: b, previous: [], analyticsEndpoint: '/api/other' }).analytics).toEqual({ endpoint: '/api/other' });
  });
});
