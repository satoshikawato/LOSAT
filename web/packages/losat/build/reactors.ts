// The engine reactors (docs/web/abi_v2.md §2) as web/adapter/tools/build_reactors.py writes them:
// the directory holds losat-web-serial.wasm, losat-web-threads.wasm and artifacts.json, the
// record of their SHA-256 and of the build's identity (plan TD-6). Both LOSAT Web's build
// (web/app/build/reactors.ts, `LOSAT_WEB_REACTORS`) and this package's build (build/build.ts) read
// them here, so that both serve the same bytes.
//
// The threaded module is served with the shared-memory guard of the Node host
// (LOSAT/tests/wasi_shared_memory.js, applied by LOSAT/tests/wasi_thread_host.js before it
// compiles the module): V8 (Node, and Chromium in S09) can run a bulk memory operation
// against stale bounds after another thread grew the shared memory. The guard is the same
// function, so the browsers run the bytes that V-ABI runs; it changes no search code. A build
// stops rather than serve an unguarded threaded module: the guard throws on an encoding that it
// declines, and a module in which it finds nothing to guard is refused here. An artifact that is
// guarded already (it carries the guard's record) is served as it is, once that guard verifies.
// The guard also lets the waits of the module end with their worker: WebKit never ends a worker
// whose thread waits without a timeout, and keeps its shared memory, so each wait records its
// address, and once the keeper ends the waits (src/infra/reactor/shared-memory.ts endWaits) they
// wait in slices with a call of the module's import sched_yield between them. It can do so only
// in a module that imports sched_yield, so a threaded module without that import is refused here.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { join, resolve } from 'node:path';

interface Guarded {
  readonly bytes: Uint8Array;
  /** The operations that this call guarded (0 for an artifact that is guarded already). */
  readonly fill: number;
  readonly copy: number;
  readonly grow: number;
  readonly atomic: number;
  /** The guard's record in the module, or null when it found nothing to guard. */
  readonly record: { readonly version: number; readonly fill: number; readonly copy: number; readonly grow: number; readonly atomic: number } | null;
}

/** The guard of LOSAT/tests/wasi_shared_memory.js, loaded only for a build with the engine. */
function guardSharedMemory(bytes: Uint8Array): Guarded & { readonly slicesWaits: boolean } {
  const guard = createRequire(import.meta.url)('../../../../LOSAT/tests/wasi_shared_memory.js') as {
    guardSharedMemory: (bytes: Uint8Array) => Guarded;
    verifyGuarded(bytes: Uint8Array): unknown;
    slicesWaits(bytes: Uint8Array): boolean;
  };
  const guarded = guard.guardSharedMemory(bytes);
  // The guard returns an artifact that carries its record as it is, so the record is checked here.
  if (guarded.record !== null) guard.verifyGuarded(guarded.bytes);
  return { ...guarded, slicesWaits: guard.slicesWaits(guarded.bytes) };
}

export const REACTORS_ENV = 'LOSAT_WEB_REACTORS';
export const REACTOR_KINDS = ['serial', 'threads'] as const;
export type ReactorKind = (typeof REACTOR_KINDS)[number];

export interface ReactorBytes {
  /** The bytes that a build serves. */
  readonly bytes: Uint8Array;
  /** SHA-256 of the served bytes. */
  readonly sha256: string;
  /** SHA-256 of the artifact of build_reactors.py (artifacts.json). */
  readonly artifactSha256: string;
  /** The transformation of the artifact, if any. */
  readonly transform?: string;
}

export interface ReactorBuild {
  readonly dir: string;
  /** artifacts.json of build_reactors.py, as it is. */
  readonly record: unknown;
  readonly serial: ReactorBytes;
  readonly threads: ReactorBytes;
}

/** The reactors of build_reactors.py's output directory `dir`, the threaded one guarded. */
export function readReactors(dir: string): ReactorBuild {
  const at = resolve(dir);
  const record = JSON.parse(readFileSync(join(at, 'artifacts.json'), 'utf8')) as {
    artifacts?: Record<string, { sha256?: string }>;
  };
  const file = (kind: ReactorKind): ReactorBytes => {
    const path = join(at, `losat-web-${kind}.wasm`);
    const artifact = readFileSync(path);
    const artifactSha256 = sha256Hex(artifact);
    const recorded = record.artifacts?.[kind]?.sha256;
    if (artifactSha256 !== recorded) {
      throw new Error(`${path}: its SHA-256 ${artifactSha256} is not the one in artifacts.json (${String(recorded)})`);
    }
    if (kind === 'serial') return { bytes: artifact, sha256: artifactSha256, artifactSha256 };
    const guarded = guardSharedMemory(artifact);
    if (guarded.record === null) throw new Error(`${path}: the shared-memory guard found no memory.fill, memory.copy, memory.grow or atomic instruction to guard`);
    if (!guarded.slicesWaits) {
      throw new Error(`${path}: the threaded module does not import wasi_snapshot_preview1.sched_yield, so its waits would not end when the page ends their worker`);
    }
    const where = guarded.fill + guarded.copy + guarded.grow + guarded.atomic === 0 ? 'in the artifact' : 'LOSAT/tests/wasi_shared_memory.js';
    return {
      bytes: guarded.bytes,
      sha256: sha256Hex(guarded.bytes),
      artifactSha256,
      transform: `shared-memory guard (${where}: ${guarded.record.fill} fills, ${guarded.record.copy} copies, ${guarded.record.grow} grows, ${guarded.record.atomic} atomics)`,
    };
  };
  return { dir: at, record, serial: file('serial'), threads: file('threads') };
}

/** The reactors named by LOSAT_WEB_REACTORS, or undefined when it is not set. */
export function reactorsFromEnv(env: NodeJS.ProcessEnv = process.env): ReactorBuild | undefined {
  const value = env[REACTORS_ENV];
  return value === undefined || value === '' ? undefined : readReactors(value);
}

export function sha256Hex(bytes: Uint8Array): string {
  return createHash('sha256').update(bytes).digest('hex');
}
