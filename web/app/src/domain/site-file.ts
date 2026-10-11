// A file of the site with its size and SHA-256, as the release manifest
// (domain/release-manifest.ts) lists a release's files and the site index
// (domain/site-index.ts) lists the site-level files a browser needs to start offline. The
// service worker fetches each listed file and stores it only when both match.

export interface SiteFile {
  /** A path relative to the release directory (manifest) or an absolute site path (site index). */
  readonly path: string;
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes the site serves. */
  readonly sha256: string;
}

/**
 * A relative path inside a release directory: segments of letters, digits, `_`, `.`, `~` and
 * `-`, separated by `/`, none empty and none starting with `.` (so no `.` or `..`).
 */
export const RELATIVE_PATH_PATTERN = /^[A-Za-z0-9_~-][A-Za-z0-9._~-]*(?:\/[A-Za-z0-9_~-][A-Za-z0-9._~-]*)*$/;

/** An absolute path of the site: `/`, or `/` followed by a relative path (RELATIVE_PATH_PATTERN). */
export const SITE_PATH_PATTERN = /^\/(?:[A-Za-z0-9_~-][A-Za-z0-9._~-]*(?:\/[A-Za-z0-9_~-][A-Za-z0-9._~-]*)*)?$/;

export const SHA256_PATTERN = /^[0-9a-f]{64}$/;

const MAX_PATH = 512;

/**
 * Checks a list of files: each `{ path, size, sha256 }` with a path that matches `pattern`, in
 * strictly ascending order of path (code units, as Array.prototype.sort orders strings), so
 * that no path is listed twice. Fields not listed here are ignored. `what` names the list in
 * the error.
 */
export function parseSiteFiles(value: unknown, what: string, pattern: RegExp, limits: { min: number; max: number }): readonly SiteFile[] {
  if (!Array.isArray(value) || value.length < limits.min || value.length > limits.max) {
    throw new Error(`${what} is not a list of ${limits.min} to ${limits.max} files`);
  }
  let previous: string | undefined;
  const files = value.map((entry: unknown, n) => {
    if (typeof entry !== 'object' || entry === null || Array.isArray(entry)) throw new Error(`${what}: file ${n + 1} is not an object`);
    const { path, size, sha256 } = entry as Record<string, unknown>;
    if (typeof path !== 'string' || path.length > MAX_PATH || !pattern.test(path)) throw new Error(`${what}: file ${n + 1}'s path is not valid`);
    if (previous !== undefined && !(previous < path)) throw new Error(`${what}: ${path} is not in order of path, or listed twice`);
    previous = path;
    if (typeof size !== 'number' || !Number.isSafeInteger(size) || size < 0) throw new Error(`${what}: ${path}'s size is not valid`);
    if (typeof sha256 !== 'string' || !SHA256_PATTERN.test(sha256)) throw new Error(`${what}: ${path}'s sha256 is not valid`);
    return Object.freeze({ path, size, sha256 });
  });
  return Object.freeze(files);
}
