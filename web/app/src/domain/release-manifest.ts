// The release manifest `/r/<release-id>/release-manifest.json` (plan §5.9, design §14,
// docs/web/release_manifest.md): every file of a release with its size and SHA-256, and the
// entry that the boot page loads. The build writes it (build/site.ts); the boot page
// (src/boot) reads it to start the pinned release, and the service worker checks every file
// against it before it calls a release stored. Fields not listed here are ignored, so that
// an older boot page or service worker can read the manifest of a newer release of format 1.
// The imports name the `.ts` files: Node runs build/assemble-site.ts, which imports this module.
import { parseSiteFiles, RELATIVE_PATH_PATTERN, type SiteFile } from './site-file.ts';
import { RELEASE_ID_PATTERN } from './site-index.ts';

export const RELEASE_MANIFEST_NAME = 'release-manifest.json';

export interface ReleaseManifest {
  readonly format: 1;
  readonly release: string;
  /** package.json's version and the commit of the build (short SHA, or `unknown`). */
  readonly app: { readonly version: string; readonly build: string };
  /** The build names of the two engine modules; null in a FakeEngine build. */
  readonly engine: { readonly serial: string; readonly threads: string } | null;
  /** What the boot page loads, relative to the release directory; every path is one of `assets`. */
  readonly entry: { readonly script: string; readonly styles: readonly string[] };
  /** Every file of the release directory except the manifest, sorted by path. */
  readonly assets: readonly SiteFile[];
}

/** The directory of a release on the site: `/r/<id>/`. */
export function releaseDirectory(id: string): string {
  return `/r/${id}/`;
}

/** The URL path of a file of a release: `/r/<id>/<path>`. */
export function releaseFileUrl(id: string, path: string): string {
  return `${releaseDirectory(id)}${path}`;
}

/** The URL path of a release's manifest. */
export function releaseManifestUrl(id: string): string {
  return releaseFileUrl(id, RELEASE_MANIFEST_NAME);
}

/**
 * Checks a parsed release manifest; throws an Error that says what is wrong. With `release`,
 * the manifest must be that release's.
 */
export function parseReleaseManifest(value: unknown, release?: string): ReleaseManifest {
  const object = (v: unknown, what: string): Record<string, unknown> => {
    if (typeof v !== 'object' || v === null || Array.isArray(v)) throw new Error(`release manifest: ${what} is not an object`);
    return v as Record<string, unknown>;
  };
  const text = (v: unknown, what: string, pattern?: RegExp): string => {
    if (typeof v !== 'string' || v === '' || v.length > 200 || (pattern !== undefined && !pattern.test(v))) {
      throw new Error(`release manifest: ${what} is not valid`);
    }
    return v;
  };
  const manifest = object(value, 'the manifest');
  if (manifest['format'] !== 1) throw new Error('release manifest: unknown format');
  const id = text(manifest['release'], 'release', RELEASE_ID_PATTERN);
  if (release !== undefined && id !== release) throw new Error(`release manifest: it is release ${id}'s, not ${release}'s`);
  const app = object(manifest['app'], 'app');
  const engine = manifest['engine'] === null ? null : object(manifest['engine'], 'engine');
  const entry = object(manifest['entry'], 'entry');
  const styles = entry['styles'];
  if (!Array.isArray(styles) || styles.length > 100) throw new Error('release manifest: entry styles is not a list');
  const assets = parseSiteFiles(manifest['assets'], 'release manifest assets', RELATIVE_PATH_PATTERN, { min: 1, max: 100_000 });
  if (assets.some((file) => file.path === RELEASE_MANIFEST_NAME)) throw new Error('release manifest: it lists itself');
  const listed = new Set(assets.map((file) => file.path));
  const entryPath = (v: unknown, what: string): string => {
    const path = text(v, what, RELATIVE_PATH_PATTERN);
    if (!listed.has(path)) throw new Error(`release manifest: ${what} ${path} is not one of the assets`);
    return path;
  };
  return Object.freeze({
    format: 1,
    release: id,
    app: Object.freeze({ version: text(app['version'], 'app version'), build: text(app['build'], 'app build') }),
    engine:
      engine === null ? null : Object.freeze({ serial: text(engine['serial'], 'engine serial'), threads: text(engine['threads'], 'engine threads') }),
    entry: Object.freeze({
      script: entryPath(entry['script'], 'entry script'),
      styles: Object.freeze(styles.map((style: unknown, n) => entryPath(style, `entry style ${n + 1}`))),
    }),
    assets,
  });
}
