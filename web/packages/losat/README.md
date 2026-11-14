# losat

LOSAT's local sequence search (BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX) for the browser.
Searches run on the user's machine, so sequences never leave it. The package holds two
WebAssembly builds of the [LOSAT](https://github.com/satoshikawato/LOSAT) engine, the same
engine as the `losat` command line, together with the TypeScript host that runs them. LOSAT Web
runs its searches through this package.

- `losat-web-serial.wasm` runs one thread and works in any browser and in Node.
- `losat-web-threads.wasm` runs a search on several threads. It needs memory that is shared
  between threads, which a browser allows only on a cross-origin isolated page (below).

The two builds speak the same ABI (version 2, `docs/web/abi_v2.md` in the repository). A search
gets an argv, as on the command line, and FASTA bytes. It writes outfmt 0, 6 and 7 byte for byte
as the command line does, and it also writes an HSP record for each hit.

## Versions

The package's version is the engine's version. `losat/abi.json` names both, together with the
ABI version that the reactors report (`losat_web2_abi_version`):

| package | engine | ABI |
|---|---|---|
| 0.2.0 | 0.2.0 | 2 |

`losat/reactors/artifacts.json` is the record of the build of the two modules: how they were built
and the SHA-256 of each as it was built. The threaded module in the package is that build with a guard
added around its shared-memory operations, so its own SHA-256 is the one in `ENGINE_ASSETS` (below).

## Serving the reactors

`import { ENGINE_ASSETS } from 'losat/reactors'` gives the URL, size and SHA-256 of each module
(`artifactSha256` is the record's, and `transform` names the guard of the threaded one). The URLs are
`new URL('./losat-web-<serial|threads>.wasm', import.meta.url)`, a form that bundlers follow: Vite's
production build, for one, copies both modules and the three worker modules of the host into its
output. Serve the modules as `application/wasm`. The host checks their SHA-256 before it runs them.

The threaded module needs `SharedArrayBuffer`, and the browser provides it only on a page that is
cross-origin isolated. Serve the page and its scripts with:

```
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

On a page that is not isolated, a search that asks for more than one thread runs on the serial
module with one thread, and `RuntimeInfo.fallbackReason` says why.

## Memory

The threaded module imports one shared memory. The module declares a maximum of 65,535 pages
(4 GiB less 64 KiB). Some browsers cannot reserve that much address space, so the host tries
65,535, 49,152, 32,768, 24,576 and 16,384 pages and takes the first one that it can create and
grow. Below 16,384 pages (1 GiB), searches run on the serial module. The serial module defines
its own memory, which can also grow to 4 GiB.

## Use

In a page, `WasmEngine` runs each search in a worker of its own. It uses the threaded module
when the page is isolated and the search asks for more than one thread:

```ts
import { WasmEngine, instantiateSerial, loadReactor, reactorControl } from 'losat';
import { ENGINE_ASSETS } from 'losat/reactors';

// describe and validate run on a serial instance in the page (here a new one for each call).
const open = async () => (await instantiateSerial((await loadReactor(ENGINE_ASSETS.serial)).module)).abi;
const control = reactorControl(open, open);
const engine = new WasmEngine({ assets: ENGINE_ASSETS, control });
```

`engine.describe(program)` lists a program's options with their defaults, and
`engine.validate(argv)` checks an argv the way the command line checks it. `engine.run(request,
port, onPhase)` runs one search. The outputs arrive on the `MessagePort` as chunks of each
stream; `RunOutputReceiver` reads them on the other end of a `MessageChannel`. `EngineInput.records`
holds the ID and length of every record of the input. `ReactorScanner` produces them, and the
engine checks them before it searches.

`engine.start(argv)` reads what the command line reads before its arguments: the file that
`-conffile` names, the file that `-logfile` names, and whether the argv has `-dryrun`. Give the
file of `-conffile` under that name to `engine.validate(argv, registry)` and as
`EngineRunRequest.conffile`; `reactorControl` validates it on an instance of its own, opened with its
second argument. A run returns the `-logfile` file in `RuntimeInfo.runFiles`. When a
run fails, `runFilesOf(error)` returns it.

In a worker or in Node, the serial module can also be called directly:

```ts
import { compileReactor, instantiateSerial, ROLE_QUERY, ROLE_SUBJECT } from 'losat';

const { module } = await compileReactor(serialBytes);
const { abi } = await instantiateSerial(module);
const query = abi.register('blastp', ROLE_QUERY, queryFasta);
const subject = abi.register('blastp', ROLE_SUBJECT, subjectFasta);
const outfmt6: Uint8Array[] = [];
// -query and -subject name the inputs in the outputs; the sequences come from the handles.
abi.run(['blastp', '-query', 'query.faa', '-subject', 'subject.faa', '-num_threads', '1'], query.handle, subject.handle, (stream, bytes) => {
  if (stream === 6) outfmt6.push(bytes.slice());
});
abi.release(query.handle);
abi.release(subject.handle);
```

The package's tests run this example on one of the repository's fixtures (`tests/package.test.ts`):
its outfmt 6 is byte for byte what NCBI BLAST+ 2.17.0 writes for the same search.

## License

MIT
