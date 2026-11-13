// Builds the npm package `losat` into dist/ (plan LIB-4, docs/losat_library_surface_plan.md §4.2):
//
// - the JavaScript and the type declarations of src/ (tsc, tsconfig.build.json), with the test
//   hooks of the Engine worker defined as false, as LOSAT Web's build defines them
//   (src/test-hooks.d.ts; web/app/vite.config.ts);
// - the two reactors of LOSAT_WEB_REACTORS (web/adapter/tools/build_reactors.py's output), read
//   and guarded by build/reactors.ts as LOSAT Web's build reads them, with their identity record
//   (reactors/artifacts.json) and reactors/index.js, which gives their URLs and digests as
//   ENGINE_ASSETS;
// - abi.json: the package's version, the engine's version (LOSAT/Cargo.toml; the package's
//   version must be the same) and the ABI version that the serial reactor reports
//   (`losat_web2_abi_version`, read here through the package's own host);
// - package.json for publishing (its exports name dist's files), README.md and LICENSE.
//
// `npm pack` and `npm publish` run in dist/. Publishing is the maintainer's act (plan §4.2).
import { execFileSync } from 'node:child_process';
import { copyFileSync, mkdirSync, readdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { REACTOR_KINDS, REACTORS_ENV, readReactors, type ReactorBytes } from './reactors.ts';

const PACKAGE = fileURLToPath(new URL('..', import.meta.url));
const DIST = join(PACKAGE, 'dist');
const REPOSITORY = join(PACKAGE, '../../..');
const TEST_HOOKS = '__LOSAT_TEST_HOOKS__';

interface Manifest {
  readonly name: string;
  readonly version: string;
  readonly description: string;
  readonly license: string;
  readonly repository: unknown;
  readonly dependencies: Readonly<Record<string, string>>;
}

/** The version of the engine crate (LOSAT/Cargo.toml [package]). */
export function engineVersion(): string {
  const toml = readFileSync(join(REPOSITORY, 'LOSAT/Cargo.toml'), 'utf8');
  const version = /^\[package\][^[]*?^version\s*=\s*"([^"]+)"/m.exec(toml)?.[1];
  if (version === undefined) throw new Error('LOSAT/Cargo.toml has no [package] version');
  return version;
}

/** The files under `dir`, relative to it, with forward slashes. */
function files(dir: string, prefix = ''): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((entry) =>
    entry.isDirectory() ? files(join(dir, entry.name), `${prefix}${entry.name}/`) : [`${prefix}${entry.name}`],
  );
}

/** The offsets of the identifiers `name` in the JavaScript `text` (comments and strings excluded). */
export function identifierOffsets(text: string, name: string): number[] {
  const ts = createRequire(import.meta.url)('typescript') as typeof import('typescript');
  const source = ts.createSourceFile('module.js', text, ts.ScriptTarget.Latest, true, ts.ScriptKind.JS);
  const offsets: number[] = [];
  const visit = (node: import('typescript').Node): void => {
    if (ts.isIdentifier(node) && node.text === name) offsets.push(node.getStart(source));
    ts.forEachChild(node, visit);
  };
  visit(source);
  return offsets;
}

/** Replaces the test hooks' flag by `false` in the emitted JavaScript; returns how many it replaced. */
function defineTestHooks(): number {
  let replaced = 0;
  for (const file of files(DIST)) {
    if (!file.endsWith('.js')) continue;
    const path = join(DIST, file);
    let text = readFileSync(path, 'utf8');
    const offsets = identifierOffsets(text, TEST_HOOKS);
    for (const offset of offsets.reverse()) text = `${text.slice(0, offset)}false${text.slice(offset + TEST_HOOKS.length)}`;
    if (offsets.length > 0) writeFileSync(path, text);
    replaced += offsets.length;
  }
  return replaced;
}

function assetSource(kind: string, reactor: ReactorBytes): string {
  const fields = [
    `url: new URL('./losat-web-${kind}.wasm', import.meta.url).href`,
    `sha256: '${reactor.sha256}'`,
    `size: ${reactor.bytes.length}`,
    `artifactSha256: '${reactor.artifactSha256}'`,
    ...(reactor.transform === undefined ? [] : [`transform: ${JSON.stringify(reactor.transform)}`]),
  ];
  return `  ${kind}: Object.freeze({\n${fields.map((field) => `    ${field},\n`).join('')}  }),\n`;
}

export async function build(env: NodeJS.ProcessEnv = process.env): Promise<void> {
  const manifest = JSON.parse(readFileSync(join(PACKAGE, 'package.json'), 'utf8')) as Manifest;
  const engine = engineVersion();
  if (manifest.version !== engine) {
    throw new Error(`package.json has version ${manifest.version}, but the engine (LOSAT/Cargo.toml) has ${engine}: the package takes the engine's version`);
  }
  const directory = env[REACTORS_ENV];
  if (directory === undefined || directory === '') {
    throw new Error(`set ${REACTORS_ENV} to the output directory of web/adapter/tools/build_reactors.py: the package carries the two reactors`);
  }
  const reactors = readReactors(directory);

  rmSync(DIST, { recursive: true, force: true });
  const tsc = createRequire(import.meta.url).resolve('typescript/bin/tsc');
  execFileSync(process.execPath, [tsc, '-p', join(PACKAGE, 'tsconfig.build.json')], { stdio: 'inherit' });
  if (defineTestHooks() === 0) throw new Error(`no ${TEST_HOOKS} in the emitted JavaScript: the Engine worker's test hooks moved`);

  const at = join(DIST, 'reactors');
  mkdirSync(at, { recursive: true });
  for (const kind of REACTOR_KINDS) writeFileSync(join(at, `losat-web-${kind}.wasm`), reactors[kind].bytes);
  writeFileSync(join(at, 'artifacts.json'), `${JSON.stringify(reactors.record, null, 2)}\n`);
  writeFileSync(
    join(at, 'index.js'),
    '// The two reactors of this package (written by build/build.ts). A bundler that follows\n' +
      '// `new URL(..., import.meta.url)` copies the modules into its build.\n' +
      `export const ENGINE_ASSETS = Object.freeze({\n${REACTOR_KINDS.map((kind) => assetSource(kind, reactors[kind])).join('')}});\n`,
  );
  writeFileSync(
    join(at, 'index.d.ts'),
    "import type { EngineAssets } from '../index.js';\n\n/** The URLs, sizes and SHA-256 of this package's two reactors. */\nexport declare const ENGINE_ASSETS: EngineAssets;\n",
  );

  // The ABI version, from the serial reactor through the package's own host.
  const host = (await import(pathToFileURL(join(DIST, 'index.js')).href)) as typeof import('../src/index.ts');
  const compiled = await host.compileReactor(reactors.serial.bytes);
  const { abi } = await host.instantiateSerial(compiled.module);
  const abiVersion = abi.abiVersion();
  writeFileSync(join(DIST, 'abi.json'), `${JSON.stringify({ package: manifest.version, engine, abi_version: abiVersion }, null, 2)}\n`);

  const entry = (path: string) => ({ types: `./${path}.d.ts`, default: `./${path}.js` });
  const published = {
    name: manifest.name,
    version: manifest.version,
    description: manifest.description,
    license: manifest.license,
    repository: manifest.repository,
    type: 'module',
    exports: {
      '.': entry('index'),
      './domain': entry('domain/index'),
      './ports': entry('ports/index'),
      './reactors': entry('reactors/index'),
      './reactors/*.wasm': './reactors/*.wasm',
      './reactors/artifacts.json': './reactors/artifacts.json',
      './abi.json': './abi.json',
    },
    dependencies: manifest.dependencies,
  };
  writeFileSync(join(DIST, 'package.json'), `${JSON.stringify(published, null, 2)}\n`);
  copyFileSync(join(PACKAGE, 'README.md'), join(DIST, 'README.md'));
  copyFileSync(join(REPOSITORY, 'LICENSE'), join(DIST, 'LICENSE'));
}

if (process.argv[1] !== undefined && import.meta.url === pathToFileURL(process.argv[1]).href) await build();
