// `vite build` writes a site (plan §5.9, docs/web/release_manifest.md):
//
//   <site>/_headers, ...         the files of public/ (the site root for Cloudflare Pages)
//   <site>/index.html            the boot page (src/boot/index.html)
//   <site>/boot.js, boot.css     the boot script (src/boot/boot.ts, built here) and its style sheet
//   <site>/404.html              the page of a missing file (src/boot/not-found.html)
//   <site>/releases.json         the site index: this release, the shell, the measurement setting
//   <site>/r/<release-id>/       the Vite build with base /r/<release-id>/ and its release-manifest.json
//
// <site> is build.outDir (`dist`, or `vite build --outDir <dir>`); the build clears it first, so
// a site never carries an older release (build/assemble-site.ts adds older releases). The dev
// server is not affected: it serves index.html -> src/main.ts at `/`, without the boot page.
import { copyFileSync, cpSync, existsSync, mkdirSync, readdirSync, rmSync, writeFileSync } from 'node:fs';
import { join, parse, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import { build, type Plugin } from 'vite';
import { releaseDirectory, type ReleaseManifest } from '../src/domain/release-manifest.ts';
import { engineBuildName } from '../src/infra/reactor/assets.ts';
import type { Reactors, ReactorFile } from './reactors.ts';
import { analyticsEndpoint, siteRelease, writeReleaseManifest, writeSiteIndex } from './site-files.ts';

const APP_DIR = fileURLToPath(new URL('..', import.meta.url));
const BOOT_DIR = join(APP_DIR, 'src/boot');
/** The site-level files copied from src/boot: source name -> name in the site. */
const BOOT_PAGES: readonly (readonly [string, string])[] = [
  ['index.html', 'index.html'],
  ['boot.css', 'boot.css'],
  ['not-found.html', '404.html'],
];

export interface SiteOptions {
  /** The release ID (build/release.ts). */
  readonly release: string;
  /** package.json's version and the build's commit. */
  readonly app: { readonly version: string; readonly build: string };
  /** The engine modules of the build; undefined for a FakeEngine build. */
  readonly reactors: Reactors | undefined;
}

/** The engine build names of the release manifest (src/infra/reactor/assets.ts), or null without the engine. */
export function engineNames(reactors: Reactors | undefined): ReleaseManifest['engine'] {
  if (reactors === undefined) return null;
  const name = (kind: 'serial' | 'threads', file: ReactorFile) =>
    engineBuildName(kind, {
      url: file.name,
      sha256: file.sha256,
      size: file.bytes.length,
      artifactSha256: file.artifactSha256,
      ...(file.transform === undefined ? {} : { transform: file.transform }),
    });
  return { serial: name('serial', reactors.serial), threads: name('threads', reactors.threads) };
}

/**
 * The entry of a release from the index.html that Vite built for it: the one module script and
 * the style sheets, as paths relative to the release directory. Anything else that the boot
 * page would not load (an inline script or style, another kind of link) stops the build.
 */
export function entryOf(html: string, base: string): ReleaseManifest['entry'] {
  const attributes = (text: string): Map<string, string> => {
    const map = new Map<string, string>();
    for (const [, name, double, single, bare] of text.matchAll(/([^\s=/]+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s"'>]+)))?/g)) {
      map.set(name!.toLowerCase(), double ?? single ?? bare ?? '');
    }
    return map;
  };
  const relative = (url: string | undefined, what: string): string => {
    if (url === undefined || !url.startsWith(base) || url.length === base.length) throw new Error(`index.html: ${what} ${String(url)} is not under ${base}`);
    return url.slice(base.length);
  };
  if (/<style\b/i.test(html)) throw new Error('index.html has an inline style, which the boot page does not load');
  const scripts = [...html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script\s*>/gi)].map(([, attrs, body]) => {
    const map = attributes(attrs!);
    if (map.get('type') !== 'module' || !map.has('src') || body!.trim() !== '') throw new Error('index.html has a script that is not a module script file');
    return relative(map.get('src'), 'script');
  });
  if (scripts.length !== 1) throw new Error(`index.html has ${scripts.length} module scripts, not 1`);
  const styles: string[] = [];
  for (const [, attrs] of html.matchAll(/<link\b([^>]*)>/gi)) {
    const map = attributes(attrs!);
    const rel = map.get('rel');
    if (rel === 'stylesheet') styles.push(relative(map.get('href'), 'style sheet'));
    // A module preload only hints at a static import; the entry script imports it anyway.
    else if (rel !== 'modulepreload') throw new Error(`index.html has <link rel="${String(rel)}">, which the boot page does not load (add it to src/boot/index.html)`);
  }
  return { script: scripts[0]!, styles };
}

/** Removes a site directory before a build; refuses a directory that does not look like a site. */
function clearSite(siteDir: string): void {
  if (!existsSync(siteDir)) return;
  if (siteDir === parse(siteDir).root || APP_DIR.startsWith(siteDir + sep) || `${siteDir}${sep}` === APP_DIR) {
    throw new Error(`${siteDir}: the site directory may not hold the application`);
  }
  const names = readdirSync(siteDir);
  if (names.length > 0 && !names.some((name) => name === 'releases.json' || name === 'index.html' || name === '_headers')) {
    throw new Error(`${siteDir} is not empty and is not a LOSAT Web site: remove it, or build into another directory`);
  }
  rmSync(siteDir, { recursive: true, force: true });
}

/** Builds src/boot/boot.ts into one module, `<site>/boot.js`, and copies the boot pages. */
export async function writeBoot(siteDir: string): Promise<void> {
  const result = await build({
    configFile: false,
    root: BOOT_DIR,
    publicDir: false,
    logLevel: 'warn',
    build: {
      write: false,
      target: 'es2022',
      modulePreload: false,
      copyPublicDir: false,
      rolldownOptions: { input: join(BOOT_DIR, 'boot.ts'), output: { format: 'es', entryFileNames: 'boot.js' } },
    },
  });
  const outputs = (Array.isArray(result) ? result : [result]).flatMap((output) => {
    if (!('output' in output)) throw new Error('the boot build did not return its output');
    return output.output;
  });
  const [chunk, ...rest] = outputs;
  if (chunk === undefined || chunk.type !== 'chunk' || chunk.fileName !== 'boot.js' || rest.length > 0) {
    throw new Error(`the boot build wrote ${outputs.map((item) => item.fileName).join(', ')}, not boot.js alone`);
  }
  mkdirSync(siteDir, { recursive: true });
  writeFileSync(join(siteDir, 'boot.js'), chunk.code);
  for (const [source, name] of BOOT_PAGES) copyFileSync(join(BOOT_DIR, source), join(siteDir, name));
}

/** The plugin; it applies to `vite build` only. */
export function losatSite(options: SiteOptions): Plugin {
  const base = releaseDirectory(options.release);
  let site = 'dist';
  let publicSource: string | false = 'public';
  let siteDir = '';
  let releaseDir = '';
  let publicDir: string | undefined;
  let endpoint: string | undefined;
  let entry: ReleaseManifest['entry'] | undefined;
  return {
    name: 'losat-site',
    apply: 'build',
    // After Vite's HTML plugin, so that generateBundle sees the built index.html.
    enforce: 'post',
    config(user) {
      site = user.build?.outDir ?? 'dist';
      publicSource = user.publicDir ?? 'public';
      // public/ is the site root: its files go to <site>, never into a release directory.
      return { base, publicDir: false, build: { outDir: join(site, 'r', options.release), emptyOutDir: false } };
    },
    configResolved(config) {
      siteDir = resolve(config.root, site);
      releaseDir = resolve(config.root, config.build.outDir);
      publicDir = publicSource === false ? undefined : resolve(config.root, publicSource);
      // Checked before the build starts: a wrong endpoint stops it.
      endpoint = analyticsEndpoint();
    },
    buildStart() {
      clearSite(siteDir);
      entry = undefined;
    },
    generateBundle(_, bundle) {
      const html = bundle['index.html'];
      if (html === undefined || html.type !== 'asset') this.error('the build has no index.html');
      entry = entryOf(String(html.source), base);
      // The boot page loads the entry; a release has no page of its own.
      delete bundle['index.html'];
    },
    async writeBundle() {
      if (entry === undefined) this.error('the build has no entry');
      const manifest = writeReleaseManifest(releaseDir, { release: options.release, app: options.app, engine: engineNames(options.reactors), entry });
      await writeBoot(siteDir);
      if (publicDir !== undefined && existsSync(publicDir)) cpSync(publicDir, siteDir, { recursive: true });
      writeSiteIndex(siteDir, { releases: [siteRelease(manifest)], latest: options.release, analyticsEndpoint: endpoint });
    },
  };
}
