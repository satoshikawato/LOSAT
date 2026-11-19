# losat (npm / WebAssembly)

Local NCBI BLAST+ compatible sequence search (`blastn`, `blastp`, `blastx`, `tblastn`, and `tblastx`) compiled to WebAssembly for browsers and Node.js.

Searches execute locally on the user's machine so input sequences never leave the client. This package bundles two WebAssembly builds of the [LOSAT](https://github.com/satoshikawato/LOSAT) engine along with the TypeScript host (`WasmEngine`) that powers [LOSAT Web](https://satoshikawato.github.io/LOSAT/platforms.html#losat-web):

- `losat-web-serial.wasm`: Single-threaded build that works in any modern browser or Node.js.
- `losat-web-threads.wasm`: Multi-threaded build using `SharedArrayBuffer` on cross-origin isolated pages.

Both reactors speak ABI v2 (version 3, documented in the [ABI specification](https://github.com/satoshikawato/LOSAT/blob/main/docs/web/abi_v2.md)). Searches accept CLI `argv` tokens and FASTA bytes, emitting `-outfmt 0`, `6`, and `7` byte-identically to the native CLI (when discarding any aborted stream chunks signaled on `DISCARD_STREAM`), alongside structured HSP records.

## Versions

| package | engine | ABI |
| :--- | :--- | :--- |
| `0.2.0` | `0.2.0` | `3` |

`losat/abi.json` records package and ABI versions, and `losat/reactors/artifacts.json` records the SHA-256 digests and build provenance for both `.wasm` reactors.

## Serving the Reactors and Cross-Origin Isolation

`import { ENGINE_ASSETS } from 'losat/reactors'` returns `new URL('./losat-web-<serial|threads>.wasm', import.meta.url)` URLs, sizes, and SHA-256 digests (`WasmEngine` checks each reactor's digest before instantiation). Modern bundlers (such as Vite) copy both `.wasm` modules and worker scripts automatically when served as `application/wasm`.

To enable multi-threaded execution in browsers via `SharedArrayBuffer`, serve the HTML document and scripts with:

```text
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

If a page is not cross-origin isolated (or cannot allocate shared linear memory of at least 16,384 pages / 1 GiB), multi-thread requests fall back automatically to `losat-web-serial.wasm` and report the reason in `RuntimeInfo.fallbackReason`.

## Memory Limits

The threaded reactor imports shared memory configured up to 65,535 pages (4 GiB minus 64 KiB). Because browser address-space limits vary, the host probes `65,535`, `49,152`, `32,768`, `24,576`, and `16,384` pages and selects the largest shared memory it can allocate and grow. Both serial and threaded builds support up to 4 GiB total linear memory.

## Usage

### Browser (`WasmEngine` in Dedicated Workers)

```ts
import { WasmEngine, instantiateSerial, loadReactor, reactorControl } from 'losat';
import { ENGINE_ASSETS } from 'losat/reactors';

const open = async () =>
  (await instantiateSerial((await loadReactor(ENGINE_ASSETS.serial)).module)).abi;

const control = reactorControl(open, open);
const engine = new WasmEngine({ assets: ENGINE_ASSETS, control });
```

- `engine.describe(program)`: inspects option schemas and defaults.
- `engine.validate(argv, registry)`: checks options using the engine's CLI parser.
- `engine.run(request, port, onPhase)`: executes searches in a worker and streams outputs over `MessagePort` to a `RunOutputReceiver` (pass `onDiscard` so `DISCARD_STREAM` stream 9 drops partially emitted output if a secondary `-outfmt` stream fails).
- `ReactorScanner`: parses FASTA records and sizes (`EngineInput.records`) prior to execution.
- `engine.start(argv)`: resolves `-conffile`, `-logfile`, and `-dryrun`; `-conffile` contents pass via `EngineRunRequest.conffile`, and `-logfile` outputs return on `RuntimeInfo.runFiles` (or `runFilesOf(error)` on failure).

### Direct Serial Call (Node.js >= 22.18 or Worker)

```ts
import {
  compileReactor,
  discardedStream,
  DISCARD_STREAM,
  instantiateSerial,
  ROLE_QUERY,
  ROLE_SUBJECT,
} from 'losat';

const { module } = await compileReactor(serialBytes);
const { abi } = await instantiateSerial(module);
const query = abi.register('blastp', ROLE_QUERY, queryFasta);
const subject = abi.register('blastp', ROLE_SUBJECT, subjectFasta);

let outfmt6: Uint8Array[] = [];
abi.run(
  ['blastp', '-query', 'query.faa', '-subject', 'subject.faa', '-num_threads', '1'],
  query.handle,
  subject.handle,
  (stream, bytes) => {
    if (stream === 6) outfmt6.push(bytes.slice());
    else if (stream === DISCARD_STREAM && discardedStream(bytes) === 6) outfmt6 = [];
  },
);

abi.release(query.handle);
abi.release(subject.handle);
```

See the [WebAssembly documentation](https://satoshikawato.github.io/LOSAT/platforms.html) in the LOSAT user manual for performance notes and browser limits.

## License

Released under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/main/LICENSE).
