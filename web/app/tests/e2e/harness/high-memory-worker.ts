// The threaded reactor with its heap above 2 GiB (plan TD-7; docs/web/abi_v2.md §3, §6), in
// a worker of its own, as the Engine worker is. It runs each search twice with the
// application's modules (createSharedMemory, ThreadHost and its thread workers, the WASI
// imports, ReactorAbi): once on the shared memory as the Engine worker creates it, and
// once on a memory that starts above 2 GiB and whose free memory below that is used up,
// so that the search's allocations, its thread stacks and every pointer that crosses the
// ABI and WASI are high. The page compares the two.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import { sha256Hex } from '../../../src/infra/browser/platform';
import { ThreadHost } from '../../../src/infra/engine-worker/thread-host';
import { ReactorAbi, ROLE_QUERY, ROLE_SUBJECT, type AbiExports } from '../../../src/infra/reactor/abi';
import type { MemoryLimits } from '../../../src/infra/reactor/artifact';
import { compileReactor, fetchReactor } from '../../../src/infra/reactor/instance';
import { createSharedMemory } from '../../../src/infra/reactor/shared-memory';
import { createWasi } from '../../../src/infra/reactor/wasi';
import type { HighMemoryCase, HighMemoryMode, HighMemoryRequest, HighMemoryResult } from './protocol';

const PAGE = 65_536;
const TWO_GIB = 2 ** 31;
/** The initial size of the high memory: 2 GiB and 64 MiB. */
const HIGH_INITIAL_PAGES = TWO_GIB / PAGE + 1024;

const scope = self as unknown as {
  postMessage(message: unknown): void;
  onmessage: ((event: MessageEvent<HighMemoryRequest>) => void) | null;
};

/** Uses up the free memory below 2 GiB; returns the first pointer above it. */
function fillLow(exports: AbiExports): number {
  let high = 0;
  for (const size of [PAGE, 4096, 256, 16]) {
    for (let i = 0; i < 100_000; i++) {
      const ptr = exports.losat_web2_alloc(size) >>> 0;
      if (ptr === 0) throw new Error(`the engine could not allocate ${size} bytes`);
      if (ptr >= TWO_GIB) {
        exports.losat_web2_dealloc(ptr, size);
        high = ptr;
        break;
      }
    }
  }
  return high;
}

async function runMode(
  module: WebAssembly.Module,
  limits: MemoryLimits,
  request: HighMemoryRequest,
): Promise<HighMemoryMode> {
  const { memory, maximumPages } = createSharedMemory(limits);
  const host = new ThreadHost({ module, memory, faultChannel: request.faultChannel });
  try {
    // As instantiateThreadedMain, keeping the exports for fillLow.
    const wasi = createWasi();
    const bound: { abi?: ReactorAbi } = {};
    const instance = await WebAssembly.instantiate(module, {
      losat_host: { emit: (stream: number, ptr: number, len: number) => bound.abi!.emit(stream, ptr, len) },
      env: { memory },
      wasi: { 'thread-spawn': host.spawn },
      wasi_snapshot_preview1: wasi.imports,
    });
    const exports = instance.exports as unknown as AbiExports;
    const abi = new ReactorAbi(exports, () => wasi.output());
    bound.abi = abi;
    wasi.initialize(instance);
    await host.prepare(request.threads - 1);
    const initialBytes = abi.memoryBytes();
    const firstHighPointer = initialBytes > TWO_GIB ? fillLow(exports) : 0;
    const outputs: HighMemoryCase[] = [];
    for (const search of request.cases) {
      const query = abi.register(search.program, ROLE_QUERY, search.query);
      const subject = abi.register(search.program, ROLE_SUBJECT, search.subject);
      const streams = new Map<number, Uint8Array[]>();
      abi.run([...search.argv, '-num_threads', String(request.threads)], query.handle, subject.handle, (stream, bytes) =>
        streams.set(stream, [...(streams.get(stream) ?? []), bytes.slice()]),
      );
      abi.release(query.handle);
      abi.release(subject.handle);
      const sha256: Record<number, string> = {};
      const bytes: Record<number, number> = {};
      for (const [stream, chunks] of streams) {
        const joined = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
        let offset = 0;
        for (const chunk of chunks) {
          joined.set(chunk, offset);
          offset += chunk.length;
        }
        sha256[stream] = await sha256Hex(joined);
        bytes[stream] = joined.length;
      }
      outputs.push({ program: search.program, sha256, bytes });
    }
    return { maximumPages, initialBytes, firstHighPointer, linearBytes: abi.memoryBytes(), outputs };
  } finally {
    host.terminate();
  }
}

async function run(request: HighMemoryRequest): Promise<HighMemoryResult> {
  if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
  const { module, artifact } = await compileReactor(await fetchReactor(ENGINE_ASSETS.threads));
  const declaredPages = artifact.memory.maximum!;
  const normal = await runMode(module, artifact.memory, request);
  const high = await runMode(module, { ...artifact.memory, initial: HIGH_INITIAL_PAGES }, request);
  return { declaredPages, normal, high };
}

scope.onmessage = (event) => {
  void run(event.data).then(
    (result) => scope.postMessage({ type: 'result', result }),
    (error: unknown) => scope.postMessage({ type: 'error', message: error instanceof Error ? error.message : String(error) }),
  );
};
