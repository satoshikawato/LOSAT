// The files of a LOSAT Web site on disk (docs/web/release_manifest.md): the release manifest
// of each `r/<release-id>/` directory, the shell (the site-level files a browser needs to start
// offline) and the site index `releases.json`. The site build (build/site.ts) writes them and
// build/assemble-site.ts checks and rewrites them. Every file is written from sorted lists and
// holds no time, so the same inputs give the same bytes.
import { createHash } from 'node:crypto';
import { existsSync, readdirSync, readFileSync, statSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { parseReleaseManifest, RELEASE_MANIFEST_NAME, type ReleaseManifest } from '../src/domain/release-manifest.ts';
import type { SiteFile } from '../src/domain/site-file.ts';
import { ENDPOINT_PATTERN, parseSiteIndex, type SiteIndex, type SiteRelease } from '../src/domain/site-index.ts';

export const ANALYTICS_ENV = 'LOSAT_WEB_ANALYTICS_ENDPOINT';
export const SITE_INDEX_FILE = 'releases.json';
export const RELEASES_DIR = 'r';

/** The shell: each site path and the file of the site directory that the site serves for it. */
export const SHELL_FILES: readonly (readonly [path: string, file: string])[] = [
  ['/', 'index.html'],
  ['/boot.css', 'boot.css'],
  ['/boot.js', 'boot.js'],
];

export function sha256Hex(bytes: Uint8Array): string {
  return createHash('sha256').update(bytes).digest('hex');
}

/** Every file under `dir` with its size and SHA-256; paths relative to `dir` with `/`, sorted by path. */
export function listFiles(dir: string): SiteFile[] {
  const files: SiteFile[] = [];
  const walk = (relative: string): void => {
    for (const entry of readdirSync(join(dir, relative), { withFileTypes: true })) {
      const path = relative === '' ? entry.name : `${relative}/${entry.name}`;
      if (entry.isDirectory()) walk(path);
      else if (entry.isFile()) {
        const bytes = readFileSync(join(dir, path));
        files.push({ path, size: bytes.length, sha256: sha256Hex(bytes) });
      } else throw new Error(`${join(dir, path)}: not a file or a directory`);
    }
  };
  walk('');
  return files.sort((a, b) => (a.path < b.path ? -1 : a.path > b.path ? 1 : 0));
}

/** JSON as the site writes it: two-space indentation and a final line feed. */
function writeJson(path: string, value: unknown): void {
  writeFileSync(path, `${JSON.stringify(value, null, 2)}\n`);
}

/**
 * Writes `release-manifest.json` into a release directory that holds every file of the release
 * and nothing else; the manifest lists them all. Returns the manifest (checked as a reader checks it).
 */
export function writeReleaseManifest(releaseDir: string, release: Omit<ReleaseManifest, 'format' | 'assets'>): ReleaseManifest {
  const assets = listFiles(releaseDir).filter((file) => file.path !== RELEASE_MANIFEST_NAME);
  const manifest = parseReleaseManifest(
    { format: 1, release: release.release, app: release.app, engine: release.engine, entry: release.entry, assets },
    release.release,
  );
  writeJson(join(releaseDir, RELEASE_MANIFEST_NAME), manifest);
  return manifest;
}

/**
 * Reads the manifest of a release directory and checks that the directory holds exactly the
 * files it lists, each with its size and SHA-256; throws an Error that names the first difference.
 */
export function verifyRelease(releaseDir: string, release: string): ReleaseManifest {
  const manifest = parseReleaseManifest(JSON.parse(readFileSync(join(releaseDir, RELEASE_MANIFEST_NAME), 'utf8')), release);
  const listed = new Map(manifest.assets.map((file) => [file.path, file]));
  const found = listFiles(releaseDir).filter((file) => file.path !== RELEASE_MANIFEST_NAME);
  for (const file of found) {
    const expected = listed.get(file.path);
    if (expected === undefined) throw new Error(`release ${release}: ${file.path} is not in its manifest`);
    if (expected.size !== file.size || expected.sha256 !== file.sha256) throw new Error(`release ${release}: ${file.path} is not the file its manifest lists`);
    listed.delete(file.path);
  }
  const [missing] = listed.keys();
  if (missing !== undefined) throw new Error(`release ${release}: ${missing} is missing`);
  return manifest;
}

/** The release IDs of the `r/` directory of a site, sorted. */
export function siteReleaseIds(siteDir: string): string[] {
  const dir = join(siteDir, RELEASES_DIR);
  if (!existsSync(dir)) return [];
  return readdirSync(dir, { withFileTypes: true })
    .filter((entry) => entry.isDirectory())
    .map((entry) => entry.name)
    .sort();
}

/** The shell entries of a site directory, for the files of SHELL_FILES that it has. */
export function siteShell(siteDir: string): SiteFile[] {
  return SHELL_FILES.filter(([, file]) => existsSync(join(siteDir, file))).map(([path, file]) => {
    const bytes = readFileSync(join(siteDir, file));
    return { path, size: bytes.length, sha256: sha256Hex(bytes) };
  });
}

/** The site index entry of a release, from its manifest. */
export function siteRelease(manifest: ReleaseManifest): SiteRelease {
  return { id: manifest.release, version: manifest.app.version, build: manifest.app.build };
}

/**
 * The measurement endpoint named by LOSAT_WEB_ANALYTICS_ENDPOINT, or undefined when it is not set;
 * throws for anything but a same-origin path under /api/ (plan §5.10 (b)).
 */
export function analyticsEndpoint(env: NodeJS.ProcessEnv = process.env): string | undefined {
  const value = env[ANALYTICS_ENV];
  if (value === undefined || value === '') return undefined;
  if (!ENDPOINT_PATTERN.test(value)) throw new Error(`${ANALYTICS_ENV}=${value}: not a same-origin path under /api/ (${String(ENDPOINT_PATTERN)})`);
  return value;
}

/** Writes `releases.json` with the shell of the site directory; returns the index (checked as a reader checks it). */
export function writeSiteIndex(
  siteDir: string,
  site: { readonly releases: readonly SiteRelease[]; readonly latest: string; readonly analyticsEndpoint?: string | undefined },
): SiteIndex {
  const index = parseSiteIndex({
    format: 1,
    latest: site.latest,
    releases: site.releases,
    ...(site.analyticsEndpoint === undefined ? {} : { analytics: { endpoint: site.analyticsEndpoint } }),
    shell: siteShell(siteDir),
  });
  writeJson(join(siteDir, SITE_INDEX_FILE), index);
  return index;
}

/** Reads the `releases.json` of a site directory. */
export function readSiteIndex(siteDir: string): SiteIndex {
  return parseSiteIndex(JSON.parse(readFileSync(join(siteDir, SITE_INDEX_FILE), 'utf8')));
}

/** Whether a path is a directory. */
export function isDirectory(path: string): boolean {
  return existsSync(path) && statSync(path).isDirectory();
}
