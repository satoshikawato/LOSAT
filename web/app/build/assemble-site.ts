// Assembles a site with older releases (S16 item 3, docs/web/release_manifest.md): copies the
// release directories of earlier sites into a new site, so that tabs still running an older
// release keep its files, and rewrites the new site's releases.json with every release and the
// new site's shell. Every release directory is checked against its manifest first (every file
// listed, no file unlisted, each with its size and SHA-256); one that does not match stops the
// command before anything is copied. How many releases a deployment keeps is the Owner's
// decision (pending): every release given is kept.
//
//   node build/assemble-site.ts --site dist --previous <earlier site> [--previous ...] [--latest <id>]
//
// The latest release stays the new site's unless --latest names another. The measurement
// endpoint is LOSAT_WEB_ANALYTICS_ENDPOINT when it is set, otherwise the new site's own setting.
import { cpSync, existsSync, readFileSync, renameSync, rmSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { RELEASE_MANIFEST_NAME } from '../src/domain/release-manifest.ts';
import type { SiteIndex, SiteRelease } from '../src/domain/site-index.ts';
import { analyticsEndpoint, readSiteIndex, RELEASES_DIR, siteRelease, siteReleaseIds, verifyRelease, writeSiteIndex } from './site-files.ts';

export interface AssembleOptions {
  /** The new site (a site build's output directory); changed in place. */
  readonly site: string;
  /** Earlier sites whose releases the new site keeps, in this order. */
  readonly previous: readonly string[];
  /** The release that first visits pin and updates offer; the new site's latest by default. */
  readonly latest?: string | undefined;
  /** The measurement endpoint; undefined keeps the new site's own setting. */
  readonly analyticsEndpoint?: string | undefined;
}

/** The release IDs of a site in the order of its releases.json (when it has one), then the others sorted. */
function orderedReleaseIds(siteDir: string): string[] {
  const found = siteReleaseIds(siteDir);
  let listed: string[] = [];
  try {
    listed = readSiteIndex(siteDir).releases.map((release) => release.id);
  } catch {
    // A site without a readable index: the directories alone.
  }
  return [...listed.filter((id) => found.includes(id)), ...found.filter((id) => !listed.includes(id))];
}

export function assembleSite(options: AssembleOptions): SiteIndex {
  const site = resolve(options.site);
  const index = readSiteIndex(site);
  const releases = new Map<string, SiteRelease>();
  for (const id of orderedReleaseIds(site)) releases.set(id, siteRelease(verifyRelease(join(site, RELEASES_DIR, id), id)));
  for (const release of index.releases) {
    if (!releases.has(release.id)) throw new Error(`${site}: releases.json lists ${release.id}, which the site does not have`);
  }

  // Check every release of the earlier sites before anything is copied.
  const copies: { id: string; source: string }[] = [];
  for (const previous of options.previous.map((dir) => resolve(dir))) {
    for (const id of orderedReleaseIds(previous)) {
      const source = join(previous, RELEASES_DIR, id);
      const manifest = verifyRelease(source, id);
      if (releases.has(id) || copies.some((copy) => copy.id === id)) {
        const existing = copies.find((copy) => copy.id === id)?.source ?? join(site, RELEASES_DIR, id);
        if (!readFileSync(join(existing, RELEASE_MANIFEST_NAME)).equals(readFileSync(join(source, RELEASE_MANIFEST_NAME)))) {
          throw new Error(`release ${id} of ${previous} is not the release ${id} that the site already has`);
        }
        continue;
      }
      copies.push({ id, source });
      releases.set(id, siteRelease(manifest));
    }
  }

  const latest = options.latest ?? index.latest;
  if (!releases.has(latest)) throw new Error(`--latest ${latest}: the site has no such release`);

  for (const { id, source } of copies) {
    // Copied beside the releases, checked, then moved in: a release directory is complete or absent.
    const staging = join(site, `.assemble-${id}`);
    rmSync(staging, { recursive: true, force: true });
    try {
      cpSync(source, staging, { recursive: true });
      verifyRelease(staging, id);
      renameSync(staging, join(site, RELEASES_DIR, id));
    } finally {
      rmSync(staging, { recursive: true, force: true });
    }
  }

  return writeSiteIndex(site, {
    releases: [...releases.values()],
    latest,
    analyticsEndpoint: options.analyticsEndpoint ?? index.analytics?.endpoint,
  });
}

function main(args: readonly string[]): void {
  let site: string | undefined;
  let latest: string | undefined;
  const previous: string[] = [];
  for (let i = 0; i < args.length; i++) {
    const value = args[i + 1];
    if (value === undefined) throw new Error(`assemble-site: ${args[i]} needs a value`);
    if (args[i] === '--site') site = value;
    else if (args[i] === '--previous') previous.push(value);
    else if (args[i] === '--latest') latest = value;
    else throw new Error(`assemble-site: unknown argument ${args[i]}`);
    i++;
  }
  if (site === undefined || !existsSync(site)) {
    throw new Error('usage: node build/assemble-site.ts --site <site> --previous <earlier site> [--previous ...] [--latest <id>]');
  }
  const index = assembleSite({ site, previous, latest, analyticsEndpoint: analyticsEndpoint() });
  console.log(`${resolve(site)}: ${index.releases.length} releases (${index.releases.map((release) => release.id).join(', ')}); latest ${index.latest}`);
}

if (process.argv[1] !== undefined && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  try {
    main(process.argv.slice(2));
  } catch (error) {
    console.error(error instanceof Error ? error.message : error);
    process.exit(1);
  }
}
