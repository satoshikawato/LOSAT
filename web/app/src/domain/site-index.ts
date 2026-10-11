// The site index `/releases.json` (plan §5.9, §5.10): the releases that the site serves,
// the one a first visit pins and an update offers, and the site's settings. It is written
// when the site is assembled, not by a release, so the same release bytes serve a site with
// or without measurement. It holds no research data. Fields not listed here are ignored, so
// that an older release can read the index of a newer site.
// The imports name the `.ts` files: Node runs build/assemble-site.ts, which imports this module.
import { parseSiteFiles, SITE_PATH_PATTERN, type SiteFile } from './site-file.ts';

/** Where the site serves its index. */
export const SITE_INDEX_PATH = '/releases.json';

/** Lower-case letters, digits, `.` and `-`, starting with a letter or digit, at most 64 characters (build/release.ts). */
export const RELEASE_ID_PATTERN = /^[0-9a-z][0-9a-z.-]{0,63}$/;

/** A same-origin path under /api/ (plan §5.10 (b)). */
export const ENDPOINT_PATTERN = /^\/api\/[a-z0-9][a-z0-9/-]{0,63}$/;

export interface SiteRelease {
  readonly id: string;
  /** package.json's version of the release. */
  readonly version: string;
  /** The commit the release was built from (short SHA). */
  readonly build: string;
}

export interface SiteIndex {
  readonly format: 1;
  /** The release a first visit pins and that an update offers; one of `releases`. */
  readonly latest: string;
  readonly releases: readonly SiteRelease[];
  /** The measurement endpoint; absent when the site measures nothing. */
  readonly analytics?: { readonly endpoint: string };
  /**
   * The site-level files a browser needs to start offline (`/` as served from index.html, the
   * boot script and its style sheet), sorted by path; absent in an index without the list.
   */
  readonly shell?: readonly SiteFile[];
}

/** Checks a parsed `/releases.json`; throws an Error that says what is wrong. */
export function parseSiteIndex(value: unknown): SiteIndex {
  const object = (v: unknown, what: string): Record<string, unknown> => {
    if (typeof v !== 'object' || v === null || Array.isArray(v)) throw new Error(`site index: ${what} is not an object`);
    return v as Record<string, unknown>;
  };
  const text = (v: unknown, what: string, pattern?: RegExp): string => {
    if (typeof v !== 'string' || v.length > 200 || (pattern !== undefined && !pattern.test(v))) {
      throw new Error(`site index: ${what} is not valid`);
    }
    return v;
  };
  const index = object(value, 'the index');
  if (index['format'] !== 1) throw new Error('site index: unknown format');
  const list = index['releases'];
  if (!Array.isArray(list) || list.length === 0 || list.length > 1000) throw new Error('site index: releases is not a list of 1 to 1000');
  const releases = list.map((entry, n) => {
    const release = object(entry, `release ${n + 1}`);
    return Object.freeze({
      id: text(release['id'], `release ${n + 1}'s id`, RELEASE_ID_PATTERN),
      version: text(release['version'], `release ${n + 1}'s version`),
      build: text(release['build'], `release ${n + 1}'s build`),
    });
  });
  const latest = text(index['latest'], 'latest', RELEASE_ID_PATTERN);
  if (!releases.some((release) => release.id === latest)) throw new Error('site index: latest is not one of the releases');
  const analytics =
    index['analytics'] === undefined
      ? undefined
      : Object.freeze({ endpoint: text(object(index['analytics'], 'analytics')['endpoint'], 'analytics endpoint', ENDPOINT_PATTERN) });
  const shell = index['shell'] === undefined ? undefined : parseSiteFiles(index['shell'], 'site index shell', SITE_PATH_PATTERN, { min: 0, max: 100 });
  return Object.freeze({
    format: 1,
    latest,
    releases: Object.freeze(releases),
    ...(analytics === undefined ? {} : { analytics }),
    ...(shell === undefined ? {} : { shell }),
  });
}
