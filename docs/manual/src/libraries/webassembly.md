# WebAssembly and Browser

LOSAT builds to two WebAssembly families from the same Rust codebase:

1. **Command-line WASI modules (`losat.wasm`)** — run `losat` subcommands directly in WASI runtimes such as Node.js (`wasm32-wasip1` serial and `wasm32-wasip1-threads` multi-threaded). Released directly in [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases).
2. **Browser / npm host package (`web/packages/losat`)** — TypeScript `WasmEngine` host wrapping the ABI v2 serial and threaded reactor builds for local in-browser BLAST searches (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`) inside Web Workers. Used by LOSAT Web and third-party web tools such as [gbdraw](https://github.com/satoshikawato/gbdraw).

## WASI Command Modules (`losat.wasm`)

Both WASI archives (`LOSAT-0.2.0-wasm32-wasip1.tar.gz` and `LOSAT-0.2.0-wasm32-wasip1-threads.tar.gz`) accept programs directly as the first argument (`blastn ...` rather than `losat blastn ...`).

Use the bundled Node.js 24+ launchers in `LOSAT/tests/` to pass working directory (`PWD`), time zones (`TZ`), file handles, and >2 GiB unsigned linear-memory pointers:

```bash
# Serial WASI build
node LOSAT/tests/run_losat_wasi.js losat.wasm blastn -query query.fna -subject subject.fna -outfmt 6

# Threaded WASI build
node LOSAT/tests/run_losat_wasi_threads.js losat-threads.wasm blastn -query query.fna -subject subject.fna -num_threads 4 -outfmt 6
```

See [Platforms and WebAssembly](../platforms.md#wasi-losatwasm) for required WASI host capabilities, compression limits, and throughput benchmarks vs native execution.

## Embedding the Browser Reactors (`WasmEngine`)

The TypeScript package in `web/packages/losat` exports a worker-backed `WasmEngine` alongside `ENGINE_ASSETS`:

```ts
import { WasmEngine, instantiateSerial, loadReactor, reactorControl } from 'losat';
import { ENGINE_ASSETS } from 'losat/reactors';

const open = async () =>
  (await instantiateSerial((await loadReactor(ENGINE_ASSETS.serial)).module)).abi;

const control = reactorControl(open, open);
const engine = new WasmEngine({ assets: ENGINE_ASSETS, control });
```

### Key requirements for browsers

- **Cross-origin isolation for threads:** Browsers require `SharedArrayBuffer` for `losat-web-threads.wasm`. Serve your app with:
  ```text
  Cross-Origin-Opener-Policy: same-origin
  Cross-Origin-Embedder-Policy: require-corp
  ```
  Non-isolated pages automatically fall back to `losat-web-serial.wasm` and record the fallback in `RuntimeInfo.fallbackReason`.
- **MIME type and digests:** Serve `.wasm` assets with `application/wasm`. `WasmEngine` checks each reactor's SHA-256 digest against `ENGINE_ASSETS` before startup.
- **Memory probing:** Shared linear memory defaults up to 65,535 pages (4 GiB minus 64 KiB) and steps down (`49,152`, `32,768`, `24,576`, `16,384` pages) if the browser restricts large shared allocations.
- **Stream discards (`DISCARD_STREAM`):** When using lower-level direct ABI calls, listen for stream `9` (`DISCARD_STREAM`) via `RunOutputReceiver` so aborted partial format chunks are discarded cleanly, preserving byte parity on surviving outputs.

Browser reactors accept FASTA payloads in memory (`-subject`) and do not open `-db` database volumes. For full API signatures, see `web/packages/losat/README.md` and the [ABI v2 specification](https://github.com/satoshikawato/LOSAT/blob/main/docs/web/abi_v2.md).
