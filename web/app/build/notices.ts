// The license and the notices that the site carries, at the root of the build, with `.txt`
// names so that a browser shows them as text:
// - LICENSE.txt, THIRD_PARTY_NOTICES.txt and THIRD_PARTY_LICENSES.txt: copied from the
//   repository's root at every build, as the package's build does
//   (web/packages/losat/build/build.ts), so that the site's copies are never older than the
//   root's. The last two cover the engine modules: the code from other projects in LOSAT's
//   source, and the crates, the Rust standard library and wasi-libc that the modules link.
// - THIRD_PARTY_LICENSES_NPM.txt: the npm packages whose code is in the site's JavaScript,
//   with their license files, from the modules of every bundle: the main bundle and the
//   worker bundles. Vite's own `build.license` reads the main bundle only and skips the helper
//   code that Vite and Rolldown add as virtual modules, so it is not used.
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import type { Plugin } from 'vite';

const APP = fileURLToPath(new URL('..', import.meta.url));
const REPOSITORY = join(APP, '..', '..');

/** The files copied from the repository's root: the name there, and the name in the build. */
export const NOTICE_FILES = [
  { source: 'LICENSE', name: 'LICENSE.txt' },
  { source: 'THIRD_PARTY_NOTICES', name: 'THIRD_PARTY_NOTICES.txt' },
  { source: 'THIRD_PARTY_LICENSES', name: 'THIRD_PARTY_LICENSES.txt' },
] as const;

/** The file in the build that lists the npm packages of the bundles and their licenses. */
export const NPM_LICENSES_FILE = 'THIRD_PARTY_LICENSES_NPM.txt';

/**
 * Virtual modules (ids that start with "\0") whose code comes from a project that no installed
 * package's license files name: the build stops at them. Vite's module preload polyfill was
 * adapted from es-module-shims (MIT, Guy Bedford), as Vite 5's source of it says; Rolldown's
 * copy of it in Vite 8 says nothing. vite.config.ts turns it off: the site has one chunk and no
 * modulepreload link, so it would do nothing.
 */
const VIRTUAL_WITHOUT_NOTICE: readonly (readonly [prefix: string, why: string])[] = [
  ['\0vite/modulepreload-polyfill', "Vite's module preload polyfill, adapted from es-module-shims (MIT, Guy Bedford), whose notice no package carries; keep build.modulePreload.polyfill off or give that notice"],
];

/**
 * The packages whose code a virtual module holds. In Vite 8, Rolldown's built-in plugins write
 * Vite's helpers, so those name both. LOSAT's own virtual modules (build/reactors.ts,
 * build/verification.ts) need no notice. A virtual module of any other owner stops the build.
 */
const VIRTUAL_OWNERS: readonly (readonly [prefix: string, packages: readonly string[]])[] = [
  ['\0virtual:losat-', []],
  ['\0vite/', ['vite', 'rolldown']],
  ['\0rolldown/', ['rolldown']],
  ['\0plugin-vue:', ['@vitejs/plugin-vue']],
];

/**
 * Where a package's license file appends the licenses of what is bundled into the package's
 * own build, none of which reaches the site: only the text before the marker is printed.
 */
const LICENSE_ENDS: Readonly<Record<string, string>> = {
  vite: '\n# Licenses of bundled dependencies',
};

/** License, NOTICE (Apache-2.0 section 4(d)) and third-party license files, such as Rolldown's
 * THIRD-PARTY-LICENSE (the notices of Rollup and esbuild, whose helpers it writes). */
const LICENSE_FILE = /^(?:licen[cs]e|copying|notice|third[-_]party[-_]licen[cs]es?)(?:[-_.].*)?$/i;
const RULE = '='.repeat(80);

interface NpmPackage {
  readonly name: string;
  readonly version: string;
  readonly license: string;
  /** The author and the repository that package.json names, if it does. */
  readonly author: string | undefined;
  readonly repository: string | undefined;
  readonly dir: string;
}

type Named = string | { readonly name?: string; readonly url?: string };

function readPackage(dir: string): NpmPackage | undefined {
  const file = join(dir, 'package.json');
  if (!existsSync(file)) return undefined;
  const data = JSON.parse(readFileSync(file, 'utf8')) as { name?: string; version?: string; license?: string; author?: Named; repository?: Named };
  if (data.name === undefined) return undefined;
  const text = (value: Named | undefined, key: 'name' | 'url') => (typeof value === 'string' ? value : value?.[key]);
  return {
    name: data.name,
    version: data.version ?? '',
    license: data.license?.trim() ?? 'no license field',
    author: text(data.author, 'name'),
    repository: text(data.repository, 'url')?.replace(/^git\+/, ''),
    dir,
  };
}

/** The package of a module in node_modules: the nearest package.json with a name. */
function packageOfModule(file: string): NpmPackage {
  for (let dir = dirname(file); !dir.endsWith('node_modules') && dir !== dirname(dir); dir = dirname(dir)) {
    const found = readPackage(dir);
    if (found !== undefined) return found;
  }
  throw new Error(`${file}: no package.json with a name above this module`);
}

/** The installed package of a name, from the node_modules directories above the application. */
function packageNamed(name: string, from: string = APP): NpmPackage {
  for (let dir = from; ; dir = dirname(dir)) {
    const found = readPackage(join(dir, 'node_modules', name));
    if (found !== undefined) return found;
    if (dir === dirname(dir)) throw new Error(`the package ${name} is not installed above ${from}`);
  }
}

/** The packages whose code the modules (the `moduleIds` of the bundles' chunks) hold. */
export function npmPackages(moduleIds: Iterable<string>, from: string = APP): NpmPackage[] {
  const found = new Map<string, NpmPackage>();
  const add = (pkg: NpmPackage) => found.set(`${pkg.name}@${pkg.version}`, pkg);
  for (const id of moduleIds) {
    if (id.startsWith('\0')) {
      const unnoticed = VIRTUAL_WITHOUT_NOTICE.find(([prefix]) => id.startsWith(prefix));
      if (unnoticed !== undefined) throw new Error(`the virtual module ${JSON.stringify(id)} is in a bundle: ${unnoticed[1]} (web/app/build/notices.ts)`);
      const owner = VIRTUAL_OWNERS.find(([prefix]) => id.startsWith(prefix));
      if (owner === undefined) {
        throw new Error(`the virtual module ${JSON.stringify(id)} is in a bundle; name the packages of its code in VIRTUAL_OWNERS (web/app/build/notices.ts)`);
      }
      for (const name of owner[1]) add(packageNamed(name, from));
      continue;
    }
    const file = id.split('?')[0]!;
    if (file.split(/[\\/]/).includes('node_modules')) add(packageOfModule(file));
  }
  // Ordered by code unit, so that the file does not depend on the build machine's locale.
  const key = (pkg: NpmPackage) => `${pkg.name}\0${pkg.version}`;
  return [...found.values()].sort((a, b) => (key(a) < key(b) ? -1 : key(a) > key(b) ? 1 : 0));
}

function licenseTexts(pkg: NpmPackage): string[] {
  const files = readdirSync(pkg.dir)
    .filter((name) => LICENSE_FILE.test(name) && statSync(join(pkg.dir, name)).isFile())
    .sort();
  if (files.length === 0) throw new Error(`${pkg.name} ${pkg.version}: no license file in ${pkg.dir}`);
  return files.map((name) => {
    let text = readFileSync(join(pkg.dir, name), 'utf8');
    const end = LICENSE_ENDS[pkg.name];
    if (end !== undefined) {
      const at = text.indexOf(end);
      if (at < 0) throw new Error(`${pkg.name}: ${name} no longer has ${JSON.stringify(end.trim())}; check what it holds`);
      text = text.slice(0, at);
    }
    return `${name}:\n\n${text.trim()}\n`;
  });
}

/** THIRD_PARTY_LICENSES_NPM.txt for the packages. */
export function npmLicensesText(packages: readonly NpmPackage[]): string {
  const title = (pkg: NpmPackage) => `${pkg.name} ${pkg.version} (${pkg.license})`;
  const out = [
    'THIRD_PARTY_LICENSES_NPM',
    '',
    'The npm packages whose code is in the JavaScript of LOSAT Web, with their license files as',
    'each package ships them. This file is generated at every build of the site',
    '(web/app/build/notices.ts) from the modules of the main bundle and of the worker bundles.',
    "LOSAT's own work is under the MIT License (LICENSE.txt). The code from other projects in",
    "LOSAT's source is in THIRD_PARTY_NOTICES.txt, and the code that the engine modules (the",
    'two .wasm files) link is in THIRD_PARTY_LICENSES.txt.',
    '',
    'Contents:',
    ...packages.map((pkg) => `  - ${title(pkg)}`),
  ];
  for (const pkg of packages) {
    const about = [pkg.author === undefined ? [] : [`Author: ${pkg.author}`], pkg.repository === undefined ? [] : [`Source: ${pkg.repository}`]].flat();
    out.push('', RULE, title(pkg), RULE, '', ...about, ...(about.length > 0 ? [''] : []), licenseTexts(pkg).join('\n'));
  }
  return `${out.join('\n').trimEnd()}\n`;
}

/**
 * The plugins: `site` for the main bundle (it emits the files) and `worker()` for each worker
 * bundle (it hands the worker's modules to `site`). Worker bundles are built while the main
 * bundle is, before its files are written.
 */
export function losatNotices(repository: string = REPOSITORY): { readonly site: Plugin; readonly worker: () => Plugin } {
  const workerModules = new Set<string>();
  let building = false;
  const site: Plugin = {
    name: 'losat-notices',
    configResolved(config) {
      building = config.command === 'build';
    },
    buildStart() {
      if (!building) return;
      for (const file of NOTICE_FILES) {
        this.emitFile({ type: 'asset', fileName: file.name, source: readFileSync(join(repository, file.source)) });
      }
    },
    generateBundle(_, bundle) {
      const modules = new Set(workerModules);
      for (const output of Object.values(bundle)) {
        if (output.type === 'chunk') for (const id of output.moduleIds) modules.add(id);
      }
      this.emitFile({ type: 'asset', fileName: NPM_LICENSES_FILE, source: npmLicensesText(npmPackages(modules)) });
    },
  };
  const worker = (): Plugin => ({
    name: 'losat-notices-worker',
    generateBundle(_, bundle) {
      for (const output of Object.values(bundle)) {
        if (output.type === 'chunk') for (const id of output.moduleIds) workerModules.add(id);
      }
    },
  });
  return { site, worker };
}
