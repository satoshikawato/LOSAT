// Sites for V-OFF (tests/e2e/offline.spec.ts): releases built with `vite build` into temporary
// directories (the same build as `npm run build`, with LOSAT_WEB_RELEASE_ID and, in an engine
// build, LOSAT_WEB_REACTORS from the environment), a site that serves an older release next to
// a newer one, and what the browser should store of them.
import { execFileSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { parseReleaseManifest, type ReleaseManifest } from '../../../src/domain/release-manifest';
import { assembleSite } from '../../../build/assemble-site';
import { readSiteIndex, RELEASES_DIR, sha256Hex } from '../../../build/site-files';

const APP_DIR = fileURLToPath(new URL('../../../', import.meta.url));
const VITE = join(APP_DIR, 'node_modules/vite/bin/vite.js');

/** Builds a site with one release, `release`, into `dir` (cleared first by the build). */
export function buildSite(dir: string, release: string): void {
  execFileSync(process.execPath, [VITE, 'build', '--outDir', dir, '--logLevel', 'warn'], {
    cwd: APP_DIR,
    env: { ...process.env, LOSAT_WEB_RELEASE_ID: release },
    stdio: ['ignore', 'ignore', 'pipe'],
  });
}

/** Adds the releases of `previous` to `site` (build/assemble-site.ts); the site's own latest stays latest. */
export function assemble(site: string, previous: readonly string[]): void {
  assembleSite({ site, previous });
}

/** The manifest of a release of a site directory. */
export function manifestOf(site: string, release: string): ReleaseManifest {
  return parseReleaseManifest(JSON.parse(readFileSync(join(site, RELEASES_DIR, release, 'release-manifest.json'), 'utf8')), release);
}

/** The largest file of a release: the one whose damage the update must notice. */
export function largestAsset(manifest: ReleaseManifest): string {
  return [...manifest.assets].sort((a, b) => b.size - a.size || (a.path < b.path ? -1 : 1))[0]!.path;
}

/** Replaces a file of a release with other bytes of the same size (its manifest is not changed). */
export function damage(site: string, release: string, path: string): void {
  const file = join(site, RELEASES_DIR, release, path);
  const bytes = readFileSync(file);
  for (let i = 0; i < bytes.length; i += 97) bytes[i] = bytes[i]! ^ 0xff;
  writeFileSync(file, bytes);
}

/** Makes the site's service worker another script (a new deployment's worker), byte for byte different. */
export function changeServiceWorker(site: string, mark: string): void {
  const file = join(site, 'sw.js');
  writeFileSync(file, `${readFileSync(file, 'utf8')}\n// ${mark}\n`);
}

/**
 * What a browser should store for offline use of `release` of `site`: every file of its
 * manifest and the manifest itself in the release's cache, the shell in the shell cache; by
 * URL path, with the SHA-256 of the bytes.
 */
export function expectedCaches(site: string, release: string): Record<string, Record<string, string>> {
  const manifest = manifestOf(site, release);
  const files: Record<string, string> = Object.fromEntries(manifest.assets.map((file) => [`/r/${release}/${file.path}`, file.sha256]));
  files[`/r/${release}/release-manifest.json`] = sha256Hex(readFileSync(join(site, RELEASES_DIR, release, 'release-manifest.json')));
  const shell = Object.fromEntries((readSiteIndex(site).shell ?? []).map((file) => [file.path, file.sha256]));
  return { [`losat-web-release-${release}`]: files, 'losat-web-shell': shell };
}
