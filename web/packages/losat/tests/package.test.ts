// The package's README and versions, and the built package (npm run build, dist/): what `npm pack`
// puts in it, its manifest and ABI table, the reactors' digests, how a consumer resolves its entries,
// and the README's own example on the serial reactor, whose outfmt 6 must be NCBI BLAST+ 2.17.0's
// frozen bytes for the same search (LOSAT/tests/outfmt0_manifest.tsv). The tests of dist/ skip when
// it is not built; web.yml's v-abi job builds it from the reactors of its run.
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, symlinkSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import ts from 'typescript';
import { describe, expect, it } from 'vitest';
import { engineVersion, identifierOffsets } from '../build/build.ts';

const PACKAGE = fileURLToPath(new URL('..', import.meta.url));
const DIST = join(PACKAGE, 'dist');
const REPOSITORY = join(PACKAGE, '../../..');
const built = existsSync(join(DIST, 'package.json'));
const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');
const README = readFileSync(join(PACKAGE, 'README.md'), 'utf8');
const VERSION = (JSON.parse(readFileSync(join(PACKAGE, 'package.json'), 'utf8')) as { version: string }).version;

/** The README's TypeScript examples, in order: the page's, then the worker's or Node's. */
function readmeExamples(): string[] {
  return [...README.matchAll(/^```ts\n([\s\S]*?)^```$/gm)].map((match) => match[1]!);
}

/** The source modules (src/**.ts without declarations), relative to src, without the extension. */
function sourceModules(dir = join(PACKAGE, 'src'), prefix = ''): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((entry) => {
    if (entry.isDirectory()) return sourceModules(join(dir, entry.name), `${prefix}${entry.name}/`);
    return entry.name.endsWith('.ts') && !entry.name.endsWith('.d.ts') ? [`${prefix}${entry.name.slice(0, -3)}`] : [];
  });
}

describe('the README and the versions', () => {
  it("takes the engine's version, which the app's dependency and the README's table name", () => {
    expect(VERSION).toBe(engineVersion());
    const app = JSON.parse(readFileSync(join(PACKAGE, '../../app/package.json'), 'utf8')) as { dependencies: Record<string, string> };
    expect(app.dependencies['losat']).toBe(VERSION);
    expect(README).toContain(`| ${VERSION} | ${engineVersion()} | `);
  });

  it('has a page example that typechecks against the package', () => {
    const [page] = readmeExamples();
    const options = ts.parseJsonConfigFileContent(ts.readConfigFile(join(PACKAGE, 'tsconfig.json'), ts.sys.readFile).config, ts.sys, PACKAGE).options;
    const virtual = new Map([
      [join(PACKAGE, 'readme-page.ts'), page!],
      [join(PACKAGE, 'readme-reactors.d.ts'), "import type { EngineAssets } from 'losat';\nexport declare const ENGINE_ASSETS: EngineAssets;\n"],
    ]);
    const host = ts.createCompilerHost(options);
    const { getSourceFile, fileExists, readFile } = host;
    host.fileExists = (file) => virtual.has(file) || fileExists(file);
    host.readFile = (file) => virtual.get(file) ?? readFile(file);
    host.getSourceFile = (file, language, ...rest) =>
      virtual.has(file) ? ts.createSourceFile(file, virtual.get(file)!, language) : getSourceFile(file, language, ...rest);
    const program = ts.createProgram(
      [join(PACKAGE, 'readme-page.ts'), join(PACKAGE, 'src/test-hooks.d.ts')],
      { ...options, noEmit: true, paths: { losat: ['./src/index.ts'], 'losat/reactors': ['./readme-reactors.d.ts'] } },
      host,
    );
    const problems = ts.getPreEmitDiagnostics(program).map((d) => ts.flattenDiagnosticMessageText(d.messageText, '\n'));
    expect(problems).toEqual([]);
  });
});

describe.skipIf(!built)('the built package', () => {
  it('packs its modules, their types, the reactors and nothing else', () => {
    type Pack = { name: string; version: string; files: { path: string }[] };
    // npm 10 prints a list of packs, npm 12 an object keyed by the package's name.
    const printed = JSON.parse(execFileSync('npm', ['pack', '--dry-run', '--json'], { cwd: DIST, encoding: 'utf8' })) as Pack[] | Record<string, Pack>;
    const [pack] = Array.isArray(printed) ? printed : Object.values(printed);
    const expected = [
      ...sourceModules().flatMap((module) => [`${module}.js`, `${module}.d.ts`]),
      'reactors/index.js',
      'reactors/index.d.ts',
      'reactors/losat-web-serial.wasm',
      'reactors/losat-web-threads.wasm',
      'reactors/artifacts.json',
      'abi.json',
      'package.json',
      'README.md',
      'LICENSE',
    ].sort();
    expect(pack!.name).toBe('losat');
    expect(pack!.version).toBe(engineVersion());
    expect(pack!.files.map((file) => file.path).sort()).toEqual(expected);
  });

  it('has the manifest for publishing: exact dependencies, exports that exist', () => {
    const manifest = JSON.parse(readFileSync(join(DIST, 'package.json'), 'utf8')) as {
      version: string;
      dependencies: Record<string, string>;
      exports: Record<string, string | Record<string, string>>;
      scripts?: unknown;
      devDependencies?: unknown;
    };
    expect(manifest.version).toBe(engineVersion());
    expect(manifest.scripts).toBeUndefined();
    expect((manifest as { private?: unknown }).private).toBeUndefined();
    expect(manifest.devDependencies).toBeUndefined();
    const source = JSON.parse(readFileSync(join(PACKAGE, 'package.json'), 'utf8')) as { dependencies: Record<string, string> };
    expect(manifest.dependencies).toEqual(source.dependencies);
    for (const version of Object.values(manifest.dependencies)) expect(version).toMatch(/^\d+\.\d+\.\d+$/);
    for (const target of Object.values(manifest.exports).flatMap((value) => (typeof value === 'string' ? [value] : Object.values(value)))) {
      if (!target.includes('*')) expect(existsSync(join(DIST, target)), target).toBe(true);
    }
  });

  it('defines the test hooks as false', () => {
    const hooks = readdirSync(DIST, { recursive: true, encoding: 'utf8' }).filter(
      (file) => file.endsWith('.js') && identifierOffsets(readFileSync(join(DIST, file), 'utf8'), '__LOSAT_TEST_HOOKS__').length > 0,
    );
    expect(hooks).toEqual([]);
  });

  it("names the engine's version and the reactors' ABI version, as the README's table does", () => {
    const abi = JSON.parse(readFileSync(join(DIST, 'abi.json'), 'utf8')) as { package: string; engine: string; abi_version: number };
    expect(abi).toEqual({ package: engineVersion(), engine: engineVersion(), abi_version: 2 });
    expect(README).toContain(`| ${abi.package} | ${abi.engine} | ${abi.abi_version} |`);
  });

  it('resolves for a consumer: every entry of its exports', () => {
    const consumer = mkdtempSync(join(tmpdir(), 'losat-consumer-'));
    try {
      mkdirSync(join(consumer, 'node_modules'));
      symlinkSync(DIST, join(consumer, 'node_modules', 'losat'), 'dir');
      const script = [
        "import { WasmEngine, PROGRAM_INPUTS } from 'losat';",
        "import { ENGINE_ASSETS } from 'losat/reactors';",
        "import { OUTPUT_FORMATS } from 'losat/domain';",
        "import { InputMismatchError } from 'losat/ports';",
        "import abi from 'losat/abi.json' with { type: 'json' };",
        "import record from 'losat/reactors/artifacts.json' with { type: 'json' };",
        "const wasm = import.meta.resolve('losat/reactors/losat-web-serial.wasm');",
        'console.log(JSON.stringify([typeof WasmEngine, Object.keys(PROGRAM_INPUTS).length, OUTPUT_FORMATS, typeof InputMismatchError,',
        '  abi.abi_version, typeof record.artifacts, ENGINE_ASSETS.serial.url === wasm]));',
      ].join('\n');
      const printed = execFileSync(process.execPath, ['--input-type=module', '-e', script], { cwd: consumer, encoding: 'utf8' });
      expect(JSON.parse(printed)).toEqual(['function', 5, [0, 6, 7], 'function', 2, 'object', true]);
    } finally {
      rmSync(consumer, { recursive: true, force: true });
    }
  });

  it("runs the README's worker or Node example as written: outfmt 6 is NCBI's frozen bytes", async () => {
    const fixture = readFileSync(join(REPOSITORY, 'LOSAT/tests/outfmt0_manifest.tsv'), 'utf8')
      .split('\n')
      .map((line) => line.split('\t'))
      .find((row) => row[0] === 'e2e.blastp.o_identity');
    // fixture_id program task query subject extra_args outfmt stdout_sha256
    const [, program, task, queryPath, subjectPath, extra, outfmt, frozen] = fixture!;
    expect([program, task, extra, outfmt]).toEqual(['blastp', '', '', '6']);
    // The example as the README writes it: its import from 'losat' becomes the arguments below.
    const example = readmeExamples()[1]!;
    const imported = /^import \{ ([^}]+) \} from 'losat';\n/m.exec(example);
    const names = imported![1]!.split(',').map((name) => name.trim());
    const body = ts.transpileModule(example.replace(imported![0], ''), { compilerOptions: { target: ts.ScriptTarget.ES2022, module: ts.ModuleKind.ESNext } }).outputText;
    const host = (await import(pathToFileURL(join(DIST, 'index.js')).href)) as Record<string, unknown>;
    const AsyncFunction = (async () => undefined).constructor as new (...args: string[]) => (...values: unknown[]) => Promise<Uint8Array[]>;
    const run = new AsyncFunction(...names, 'serialBytes', 'queryFasta', 'subjectFasta', `${body}\nreturn outfmt6;`);
    const outfmt6 = await run(
      ...names.map((name) => host[name]),
      readFileSync(join(DIST, 'reactors/losat-web-serial.wasm')),
      readFileSync(join(REPOSITORY, 'LOSAT', queryPath!)),
      readFileSync(join(REPOSITORY, 'LOSAT', subjectPath!)),
    );
    expect(sha256(Buffer.concat(outfmt6))).toBe(frozen);
  });

  it('describes the reactors it carries', async () => {
    const { ENGINE_ASSETS } = (await import(pathToFileURL(join(DIST, 'reactors/index.js')).href)) as typeof import('../src/index.ts') & {
      ENGINE_ASSETS: import('../src/index.ts').EngineAssets;
    };
    const record = JSON.parse(readFileSync(join(DIST, 'reactors/artifacts.json'), 'utf8')) as { artifacts: Record<string, { sha256: string }> };
    for (const kind of ['serial', 'threads'] as const) {
      const asset = ENGINE_ASSETS[kind];
      const bytes = readFileSync(fileURLToPath(asset.url));
      expect(asset.sha256).toBe(sha256(bytes));
      expect(asset.size).toBe(bytes.length);
      expect(asset.artifactSha256).toBe(record.artifacts[kind]!.sha256);
    }
    // The served threaded module is the guarded one; the serial module is the artifact itself.
    expect(ENGINE_ASSETS.serial.sha256).toBe(ENGINE_ASSETS.serial.artifactSha256);
    expect(ENGINE_ASSETS.threads.transform).toMatch(/^shared-memory guard/);
  });

});

