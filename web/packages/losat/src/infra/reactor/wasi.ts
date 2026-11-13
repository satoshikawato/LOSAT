// WASI preview1 for the engine modules, from @bjorn3/browser_wasi_shim (the shim that
// gbdraw uses; the version is pinned in package.json). The reactors get no arguments, an
// empty environment (docs/web/abi_v2.md §3), no preopened directory (so the filesystem
// functions that the engine links but the ABI does not use fail instead of touching
// files, docs/web/abi_v2.md §3), an empty standard input, and standard output and error
// kept in a small buffer. The engine's warnings do not go there (they are ABI stream 3);
// the buffer holds what the engine writes when it stops, such as a panic message, so that
// the host can report it.
//
// The shim uses its pointer arguments as they arrive. JavaScript receives a Wasm i32 as a
// signed number, so a pointer of 2 GiB or more is negative there: the shim's DataView
// accesses then throw, and its `random_get` (`subarray` with a negative start) writes
// somewhere else. Every i32 parameter of WASI preview1 is unsigned (pointer, size,
// descriptor, flags), so the imports pass each one on as its unsigned value; the i64
// parameters are BigInt and pass unchanged.
import { ConsoleStdout, File as WasiFile, OpenFile, WASI } from '@bjorn3/browser_wasi_shim';

/** How much of the engine's standard output and error is kept. */
const LOG_BYTES = 64 * 1024;

export interface WasiContext {
  /** The `wasi_snapshot_preview1` imports. */
  readonly imports: WebAssembly.ModuleImports;
  /** Binds a thread's instance of the threaded module, which must not run `_initialize`. */
  bind(instance: WebAssembly.Instance): void;
  /** Binds a reactor instance and runs its `_initialize` once (ABI v2 §2). */
  initialize(instance: WebAssembly.Instance): void;
  /** What the engine wrote to its standard output and error, last bytes first kept. */
  output(): string;
  /** Forgets what the engine wrote so far. */
  clear(): void;
}

type WasiInstance = WASI['inst'] & { exports: { _initialize?: () => unknown } };

type WasiFunction = (...args: Array<number | bigint>) => unknown;

/** The imports with every i32 argument passed on as an unsigned value. */
export function unsignedImports(imports: Record<string, unknown>): WebAssembly.ModuleImports {
  const wrapped: Record<string, WasiFunction> = {};
  for (const [name, value] of Object.entries(imports)) {
    if (typeof value !== 'function') continue;
    const call = value as WasiFunction;
    wrapped[name] = (...args) => {
      for (let i = 0; i < args.length; i++) {
        const arg = args[i];
        if (typeof arg === 'number') args[i] = arg >>> 0;
      }
      return call(...args);
    };
  }
  return wrapped as WebAssembly.ModuleImports;
}

export function createWasi(): WasiContext {
  let log = new Uint8Array(0);
  const keep = (bytes: Uint8Array) => {
    const joined = new Uint8Array(log.length + bytes.length);
    joined.set(log);
    joined.set(bytes, log.length);
    log = joined.length > LOG_BYTES ? joined.slice(joined.length - LOG_BYTES) : joined;
  };
  const output = new ConsoleStdout(keep);
  const stdin = new OpenFile(new WasiFile(new Uint8Array(0), { readonly: true }));
  // The shim logs its calls to the console unless debug is off.
  const wasi = new WASI([], [], [stdin, output, output], { debug: false });
  return {
    imports: unsignedImports(wasi.wasiImport),
    bind(instance) {
      wasi.inst = instance as unknown as WasiInstance;
    },
    initialize(instance) {
      wasi.initialize(instance as unknown as WasiInstance);
    },
    output: () => new TextDecoder().decode(log),
    clear: () => {
      log = new Uint8Array(0);
    },
  };
}
