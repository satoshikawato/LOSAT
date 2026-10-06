// WASI preview1 for the engine modules, from @bjorn3/browser_wasi_shim (the shim that
// gbdraw uses; the version is pinned in package.json). The reactors get no arguments, an
// empty environment (docs/web/abi_v2.md §3), no preopened directory (so the filesystem
// functions that the engine links but the ABI does not use fail instead of touching
// files, docs/web/abi_v2.md §3), an empty standard input, and standard output and error
// kept in a small buffer. The engine's warnings do not go there (they are ABI stream 3);
// the buffer holds what the engine writes when it stops, such as a panic message, so that
// the host can report it.
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
    imports: wasi.wasiImport as WebAssembly.ModuleImports,
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
