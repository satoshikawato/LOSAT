// The engine reactors (docs/web/abi_v2.md §2) enter the build from the directory that
// web/adapter/tools/build_reactors.py writes: `LOSAT_WEB_REACTORS=<output dir>`. The
// plugin checks the two modules against the SHA-256 of that build's artifacts.json, copies
// them into the build as engine/losat-web-<serial|threads>-<sha256 prefix>.wasm, and
// provides their URLs and digests to the application as the module `virtual:losat-engine`.
// Without the variable, the build has no engine (ENGINE_ASSETS is null) and the
// application uses the FakeEngine, with its banner (web/AGENTS.md rule 7).
//
// The threaded module is served with the shared-memory guard of the Node host
// (LOSAT/tests/wasi_shared_memory.js, applied by LOSAT/tests/wasi_thread_host.js before it
// compiles the module): V8 (Node, and Chromium in S09) can run a bulk memory operation
// against stale bounds after another thread grew the shared memory. The guard is the same
// function, so the browsers run the bytes that V-ABI runs; it changes no search code. The
// build stops rather than serve an unguarded threaded module: the guard throws on an encoding
// that it declines, and a module in which it finds nothing to guard is refused here. An
// artifact that is guarded already (it carries the guard's record) is served as it is, once
// that guard verifies.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import type { IncomingMessage, ServerResponse } from 'node:http';
import { join, resolve } from 'node:path';
import type { Plugin } from 'vite';

interface Guarded {
  readonly bytes: Uint8Array;
  /** The operations that this call guarded (0 for an artifact that is guarded already). */
  readonly fill: number;
  readonly copy: number;
  readonly atomic: number;
  /** The guard's record in the module, or null when it found nothing to guard. */
  readonly record: { readonly version: number; readonly fill: number; readonly copy: number; readonly atomic: number } | null;
}

/** The guard of LOSAT/tests/wasi_shared_memory.js, loaded only for a build with the engine. */
function guardSharedMemory(bytes: Uint8Array): Guarded {
  const guard = createRequire(import.meta.url)('../../../LOSAT/tests/wasi_shared_memory.js') as {
    guardSharedMemory: typeof guardSharedMemory;
    verifyGuarded(bytes: Uint8Array): unknown;
  };
  const guarded = guard.guardSharedMemory(bytes);
  // The guard returns an artifact that carries its record as it is, so the record is checked here.
  if (guarded.record !== null) guard.verifyGuarded(guarded.bytes);
  return guarded;
}

export const REACTORS_ENV = 'LOSAT_WEB_REACTORS';
const VIRTUAL_ID = 'virtual:losat-engine';
const RESOLVED_ID = `\0${VIRTUAL_ID}`;
const NAMES = ['serial', 'threads'] as const;

export interface ReactorFile {
  /** The file name in the build. */
  readonly name: string;
  /** The bytes that the build serves. */
  readonly bytes: Uint8Array;
  /** SHA-256 of the served bytes. */
  readonly sha256: string;
  /** SHA-256 of the artifact of build_reactors.py (artifacts.json). */
  readonly artifactSha256: string;
  /** The transformation of the artifact, if any. */
  readonly transform?: string;
}

export interface Reactors {
  readonly dir: string;
  readonly serial: ReactorFile;
  readonly threads: ReactorFile;
}

/** The reactors named by LOSAT_WEB_REACTORS, or undefined when it is not set. */
export function findReactors(env: NodeJS.ProcessEnv = process.env): Reactors | undefined {
  const value = env[REACTORS_ENV];
  if (value === undefined || value === '') return undefined;
  const dir = resolve(value);
  const artifacts = JSON.parse(readFileSync(join(dir, 'artifacts.json'), 'utf8')) as {
    artifacts?: Record<string, { sha256?: string }>;
  };
  const file = (kind: (typeof NAMES)[number]): ReactorFile => {
    const path = join(dir, `losat-web-${kind}.wasm`);
    const artifact = readFileSync(path);
    const artifactSha256 = sha256Hex(artifact);
    const recorded = artifacts.artifacts?.[kind]?.sha256;
    if (artifactSha256 !== recorded) {
      throw new Error(`${path}: its SHA-256 ${artifactSha256} is not the one in artifacts.json (${String(recorded)})`);
    }
    if (kind === 'serial') {
      return { name: `engine/losat-web-serial-${artifactSha256.slice(0, 16)}.wasm`, bytes: artifact, sha256: artifactSha256, artifactSha256 };
    }
    const guarded = guardSharedMemory(artifact);
    if (guarded.record === null) throw new Error(`${path}: the shared-memory guard found no memory.fill, memory.copy or atomic instruction to guard`);
    const where = guarded.fill + guarded.copy + guarded.atomic === 0 ? 'in the artifact' : 'LOSAT/tests/wasi_shared_memory.js';
    const sha256 = sha256Hex(guarded.bytes);
    return {
      name: `engine/losat-web-threads-${sha256.slice(0, 16)}.wasm`,
      bytes: guarded.bytes,
      sha256,
      artifactSha256,
      transform: `shared-memory guard (${where}: ${guarded.record.fill} fills, ${guarded.record.copy} copies, ${guarded.record.grow} grows, ${guarded.record.atomic} atomics)`,
    };
  };
  return { dir, serial: file('serial'), threads: file('threads') };
}

/**
 * The plugin. `emit` copies the modules into the build; worker bundles use the same
 * plugin with `emit: false`, because the main bundle already copies them. Pass the same
 * `reactors` (findReactors()) to both, so that the modules are read and guarded once.
 */
export function losatEngine(options: { readonly reactors?: Reactors | undefined; readonly emit?: boolean } = {}): Plugin {
  const reactors = 'reactors' in options ? options.reactors : findReactors();
  const emit = options.emit ?? true;
  let base = '/';
  let building = false;
  return {
    name: 'losat-engine',
    configResolved(config) {
      base = config.base;
      building = config.command === 'build';
    },
    resolveId(id) {
      return id === VIRTUAL_ID ? RESOLVED_ID : undefined;
    },
    load(id) {
      if (id !== RESOLVED_ID) return undefined;
      if (reactors === undefined) return 'export const ENGINE_ASSETS = null;\n';
      const asset = (file: ReactorFile) => ({
        url: `${base}${file.name}`,
        sha256: file.sha256,
        size: file.bytes.length,
        artifactSha256: file.artifactSha256,
        ...(file.transform === undefined ? {} : { transform: file.transform }),
      });
      const assets = { serial: asset(reactors.serial), threads: asset(reactors.threads) };
      return `export const ENGINE_ASSETS = Object.freeze(${JSON.stringify(assets)});\n`;
    },
    buildStart() {
      if (!emit || !building || reactors === undefined) return;
      for (const kind of NAMES) {
        const file = reactors[kind];
        this.emitFile({ type: 'asset', fileName: file.name, source: file.bytes });
      }
    },
    configureServer(server) {
      // The development server serves the modules from the reactor directory.
      server.middlewares.use((request: IncomingMessage, response: ServerResponse, next: () => void) => {
        const file = reactors === undefined ? undefined : NAMES.map((kind) => reactors[kind]).find((f) => request.url === `${base}${f.name}`);
        if (file === undefined) {
          next();
          return;
        }
        response.setHeader('Content-Type', 'application/wasm');
        response.end(file.bytes);
      });
    },
  };
}

function sha256Hex(bytes: Uint8Array): string {
  return createHash('sha256').update(bytes).digest('hex');
}
