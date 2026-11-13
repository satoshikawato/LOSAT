// The built package (npm run build, dist/): what `npm pack` puts in it, its manifest and ABI table,
// the reactors' digests, and the README's example on the serial reactor, whose outfmt 6 must be
// NCBI BLAST+ 2.17.0's frozen bytes for the same search (LOSAT/tests/outfmt0_manifest.tsv). The
// tests skip when dist/ is not built; web.yml's package job builds it from the reactors of its run.
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { describe, expect, it } from 'vitest';
import { engineVersion, identifierOffsets } from '../build/build.ts';

const PACKAGE = fileURLToPath(new URL('..', import.meta.url));
const DIST = join(PACKAGE, 'dist');
const REPOSITORY = join(PACKAGE, '../../..');
const built = existsSync(join(DIST, 'package.json'));
const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');

/** The source modules (src/**.ts without declarations), relative to src, without the extension. */
function sourceModules(dir = join(PACKAGE, 'src'), prefix = ''): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((entry) => {
    if (entry.isDirectory()) return sourceModules(join(dir, entry.name), `${prefix}${entry.name}/`);
    return entry.name.endsWith('.ts') && !entry.name.endsWith('.d.ts') ? [`${prefix}${entry.name.slice(0, -3)}`] : [];
  });
}

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

  it("names the engine's version and the reactors' ABI version", () => {
    expect(JSON.parse(readFileSync(join(DIST, 'abi.json'), 'utf8'))).toEqual({ package: engineVersion(), engine: engineVersion(), abi_version: 2 });
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

  it("runs the README's example: outfmt 6 is NCBI's frozen bytes", async () => {
    const fixture = readFileSync(join(REPOSITORY, 'LOSAT/tests/outfmt0_manifest.tsv'), 'utf8')
      .split('\n')
      .map((line) => line.split('\t'))
      .find((row) => row[0] === 'e2e.blastp.o_identity');
    // fixture_id program task query subject extra_args outfmt stdout_sha256
    const [, program, task, queryPath, subjectPath, extra, outfmt, frozen] = fixture!;
    expect([program, task, extra, outfmt]).toEqual(['blastp', '', '', '6']);
    const host = (await import(pathToFileURL(join(DIST, 'index.js')).href)) as typeof import('../src/index.ts');
    const { compileReactor, instantiateSerial, ROLE_QUERY, ROLE_SUBJECT } = host;
    const serialBytes = readFileSync(join(DIST, 'reactors/losat-web-serial.wasm'));
    const queryFasta = readFileSync(join(REPOSITORY, 'LOSAT', queryPath!));
    const subjectFasta = readFileSync(join(REPOSITORY, 'LOSAT', subjectPath!));

    const { module } = await compileReactor(serialBytes);
    const { abi } = await instantiateSerial(module);
    const query = abi.register('blastp', ROLE_QUERY, queryFasta);
    const subject = abi.register('blastp', ROLE_SUBJECT, subjectFasta);
    const outfmt6: Uint8Array[] = [];
    abi.run(['blastp', '-query', queryPath!, '-subject', subjectPath!, '-num_threads', '1'], query.handle, subject.handle, (stream, bytes) => {
      if (stream === 6) outfmt6.push(bytes.slice());
    });
    abi.release(query.handle);
    abi.release(subject.handle);

    expect(sha256(Buffer.concat(outfmt6))).toBe(frozen);
  });
});

/** The README's page example, checked by the typecheck only: WasmEngine runs its searches in workers. */
export async function pageExample(assets: import('../src/index.ts').EngineAssets): Promise<import('../src/index.ts').EngineGateway> {
  const { WasmEngine, instantiateSerial, loadReactor, reactorControl } = await import('../src/index.ts');
  const ENGINE_ASSETS = assets;
  const control = reactorControl(async () => (await instantiateSerial((await loadReactor(ENGINE_ASSETS.serial)).module)).abi);
  const engine = new WasmEngine({ assets: ENGINE_ASSETS, control });
  return engine;
}
