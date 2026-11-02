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

export async function compileReactor(bytes: Uint8Array): Promise<CompiledReactor> {
  const artifact = inspectReactor(bytes);
  return { module: await WebAssembly.compile(bytes as BufferSource), artifact };
}

/**
 * Fetches and compiles a module of this build. Where the browser can, the module compiles
 * while it arrives (`WebAssembly.compileStreaming`) and a second reader of the same
 * response checks the SHA-256: the module is returned only after the check, so nothing
 * instantiates other bytes. A browser may keep the code that it compiled from a streamed
 * response with the cached response, and a later visit then starts from that code (Chromium
 * does once V8 hands it the module's optimized code). The engine URLs carry their digest and
 * are served `immutable` (public/_headers), so the cached response stays valid until a build
 * changes the module.
 */
export async function loadReactor(asset: ReactorAsset): Promise<CompiledReactor> {
  const response = await fetchAsset(asset);
  if (!streamable(response)) return compileReactor(await checkedBytes(asset, response));
  const copy = response.clone();
  const compiling = WebAssembly.compileStreaming(response);
  // Reported below, after the digest check, which comes first.
  compiling.catch(() => undefined);
  const artifact = inspectReactor(await checkedBytes(asset, copy));
  return { module: await compiling, artifact };
}

async function fetchAsset(asset: ReactorAsset): Promise<Response> {
  const response = await fetch(asset.url);
  if (!response.ok) throw new Error(`the engine file ${asset.url} could not be loaded (HTTP ${response.status})`);
  return response;
}

async function checkedBytes(asset: ReactorAsset, response: Response): Promise<Uint8Array> {
  const bytes = new Uint8Array(await response.arrayBuffer());
  const digest = await sha256Hex(bytes);
  if (digest !== asset.sha256) {
    throw new Error(`the engine file ${asset.url} is not the one of this release (SHA-256 ${digest})`);
  }
  return bytes;
}

/** compileStreaming needs the Wasm media type; other responses compile from their bytes. */
function streamable(response: Response): boolean {
  const type = response.headers.get('Content-Type')?.split(';')[0]?.trim().toLowerCase();
  return typeof WebAssembly.compileStreaming === 'function' && type === 'application/wasm';
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
  // The instance is prepared before its thread exists, so it may be older than the memory
  // the engine grew for that thread's stack and TLS. In the Node host of the same ABI
  // (LOSAT/tests/wasi_thread_host.js), such a worker's view of the memory still ended below
  // them and the thread's first atomic store trapped "memory access out of bounds" in about
  // 0.3 % of runs of one 16-thread check; memory.grow(0) updates the worker's view. Not
  // measured in a browser (docs/evidence/wasm_oob_20261025/). Later atomic accesses are
  // covered by the build's guard (web/app/build/reactors.ts, docs/evidence/wasm_oob_20261101/).
  return {
    start: (tid, startArg) => {
      memory.grow(0);
      (start as (tid: number, arg: number) => void)(tid, startArg);
    },
  };
}
