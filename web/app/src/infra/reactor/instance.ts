// Loading and instantiating the engine modules (docs/web/abi_v2.md §2-§3).
import { sha256Hex } from '../browser/platform';
import { ReactorAbi, type AbiExports } from './abi';
import { inspectReactor, type ReactorArtifact } from './artifact';
import type { ReactorAsset } from './assets';
import { createWasi, type WasiContext } from './wasi';

export interface CompiledReactor {
  readonly module: WebAssembly.Module;
  readonly artifact: ReactorArtifact;
}

export interface ReactorInstance {
  readonly abi: ReactorAbi;
  readonly wasi: WasiContext;
}

/** Fetches a module of this build and checks its SHA-256 before anything uses it. */
export async function fetchReactor(asset: ReactorAsset): Promise<Uint8Array> {
  const response = await fetch(asset.url);
  if (!response.ok) throw new Error(`the engine file ${asset.url} could not be loaded (HTTP ${response.status})`);
  const bytes = new Uint8Array(await response.arrayBuffer());
  const digest = await sha256Hex(bytes);
  if (digest !== asset.sha256) {
    throw new Error(`the engine file ${asset.url} is not the one of this release (SHA-256 ${digest})`);
  }
  return bytes;
}

export async function compileReactor(bytes: Uint8Array): Promise<CompiledReactor> {
  const artifact = inspectReactor(bytes);
  return { module: await WebAssembly.compile(bytes as BufferSource), artifact };
}

/** Imports that are neither WASI nor the shared memory: `losat_host.emit` (ABI v2 §3). */
function hostImports(emit: (stream: number, ptr: number, len: number) => void): WebAssembly.Imports {
  return { losat_host: { emit } };
}

function checkVersion(abi: ReactorAbi): ReactorAbi {
  const version = abi.abiVersion();
  if (version !== 2) throw new Error(`the engine module has ABI version ${version}; LOSAT Web needs version 2`);
  return abi;
}

/** Instantiates the serial module and runs its `_initialize`. */
export async function instantiateSerial(module: WebAssembly.Module): Promise<ReactorInstance> {
  const wasi = createWasi();
  const bound: { abi?: ReactorAbi } = {};
  const instance = await WebAssembly.instantiate(module, {
    ...hostImports((stream, ptr, len) => bound.abi!.emit(stream, ptr, len)),
    wasi_snapshot_preview1: wasi.imports,
  });
  const exports = instance.exports as unknown as AbiExports;
  const abi = new ReactorAbi(exports, () => wasi.output());
  bound.abi = abi;
  wasi.initialize(instance);
  return { abi: checkVersion(abi), wasi };
}

/**
 * Instantiates the threaded module on the thread that calls the exports (the Engine
 * worker), with the shared memory and the `wasi.thread-spawn` of a ThreadHost, and runs
 * its `_initialize`.
 */
export async function instantiateThreadedMain(
  module: WebAssembly.Module,
  memory: WebAssembly.Memory,
  spawn: (startArg: number) => number,
): Promise<ReactorInstance> {
  const wasi = createWasi();
  const bound: { abi?: ReactorAbi } = {};
  const instance = await WebAssembly.instantiate(module, {
    ...hostImports((stream, ptr, len) => bound.abi!.emit(stream, ptr, len)),
    env: { memory },
    wasi: { 'thread-spawn': spawn },
    wasi_snapshot_preview1: wasi.imports,
  });
  const abi = new ReactorAbi(instance.exports as unknown as AbiExports, () => wasi.output());
  bound.abi = abi;
  wasi.initialize(instance);
  return { abi: checkVersion(abi), wasi };
}

/**
 * Instantiates the threaded module for one thread (a thread worker). The instance shares
 * the memory; it must not run `_initialize`, and its `losat_host.emit` is never called
 * (the engine formats on the calling thread, ABI v2 §5), so calling it is an error.
 */
export async function instantiateThread(
  module: WebAssembly.Module,
  memory: WebAssembly.Memory,
  wasi: WasiContext,
): Promise<{ readonly start: (tid: number, startArg: number) => void }> {
  const instance = await WebAssembly.instantiate(module, {
    ...hostImports(() => {
      throw new Error('losat_host.emit is not available on a thread of the search pool');
    }),
    env: { memory },
    wasi: { 'thread-spawn': () => -1 },
    wasi_snapshot_preview1: wasi.imports,
  });
  wasi.bind(instance);
  const start = instance.exports['wasi_thread_start'];
  if (typeof start !== 'function') throw new Error('the threaded engine module lacks wasi_thread_start');
  return { start: (tid, startArg) => void (start as (tid: number, arg: number) => void)(tid, startArg) };
}
