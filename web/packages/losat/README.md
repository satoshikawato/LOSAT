# losat

LOSAT's local sequence search (BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX) for the browser.
Searches run on the user's machine, so sequences never leave it. The package holds two
WebAssembly builds of the [LOSAT](https://github.com/satoshikawato/LOSAT) engine, the same
engine as the `losat` command line, together with the TypeScript host that runs them. LOSAT Web
runs its searches through this package.

- `losat-web-serial.wasm` runs one thread and works in any browser and in Node.
- `losat-web-threads.wasm` runs a search on several threads. It needs memory that is shared
  between threads, which a browser allows only on a cross-origin isolated page (below).

The two builds speak the same ABI (ABI v2 at version 11, `docs/web/abi_v2.md` in the repository). A search
gets an argv, as on the command line, and FASTA bytes. It writes outfmt 0, 6 and 7 byte for byte
as the command line does, and it also writes an HSP record for each hit. This holds for a host that
drops the output that the engine tells it to discard (stream 9, in "Use").

The package is not published on npm. Build it from the [LOSAT repository](https://github.com/satoshikawato/LOSAT):
make the two modules with `web/adapter/tools/build_reactors.py`, then run
`LOSAT_WEB_REACTORS=<modules directory> npm run build` in `web/packages/losat`; `npm pack` there makes a
file to install ([web/AGENTS.md](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/web/AGENTS.md)).

## Versions

The package's version is the engine's version. `losat/abi.json` names both, together with the
ABI version that the reactors report (`losat_web2_abi_version`):

| package | engine | ABI |
|---|---|---|
| 0.2.0 | 0.2.0 | 11 |

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
engine checks them before it searches. When the reader fails on a query record, the scan ends there:
its records are those before it, and `reader_error` gives the command line's message, the offset
of that record, and whether it starts with a defline (`defline`). Pass the query's bytes from that
offset on after the records you keep, and the run searches the query batches before the record and
ends with that message, as the command line does. Bytes from an offset whose record has no defline
(an input's first line, such as a byte order mark) would be read as letters of the record before
them: do not pass them after another input's records.
A subject with such a record is refused when it is registered, as the command line reads every
subject before its first query.

The engine tells the host which outputs a run could not write, and where their errors are: for
each, it sends stream 9 (`DISCARD_STREAM`), whose 8 bytes are that output's stream number and then
the place of its error, each a little-endian unsigned integer (`discardSignal` reads them), and the
host drops everything that arrived of that stream, which then stays empty. The place is where the
run wrote the error, as the command line writes the error of that output. A BLAST+ error goes to
stream 3, or nowhere (`PLACE_NOWHERE`) where the run's messages go nowhere (TBLASTX's `-logfile` of
a name whose messages go nowhere: `""`, `/dev/null`, `/dev/fd/2`, `/dev/full`...). In TBLASTX with `-logfile NAME` it goes to stream 3 too: the file holds
what the command line's file holds for the run's own `-outfmt`, which does not write that output.
LOSAT's own refusal goes to stream 3 whatever `-logfile` is, as the command line writes it to its
standard error. The signal comes for each output that fails on its own,
whether or not some of its chunks were sent, and for each output of a search that fails in a run of
two searches, and the run goes on with the other outputs. These signals are the only sign of an
output that could not be written: an empty stream is not one (outfmt 6 of a search without hits is
empty). `RunOutputReceiver` reports each signal through `onDiscard`, its required third argument
(`new RunOutputReceiver(port, onChunk, onDiscard)`; `onDiscard(stream, error)` drops what arrived
of that output, whose error is at `error`: `'diagnostics'`, `'logfile'` or `'nowhere'`), and
`RunOutputMessage` has a `discard` variant for it.

A host that keeps every chunk would keep output that the command line does not write.

`engine.start(argv)` reads what the command line reads before its arguments: the file that
`-conffile` names, the file that `-logfile` names, the search strategy file that
`-export_search_strategy` names (`StartInfo.strategyFile`), and whether the argv has `-dryrun`. Give the
file of `-conffile` under that name to `engine.validate(argv, registry)` and as
`EngineRunRequest.conffile`; `reactorControl` validates it on an instance of its own, opened with its
second argument. A run returns the `-logfile` file and the search strategy file in
`RuntimeInfo.runFiles`, each `RunFile` with the option that names it (`runFileOption(file)`:
`'logfile'` or `'export_search_strategy'`). When a run fails, `runFilesOf(error)` returns them. A run that the engine failed after it began to write ends
its output with `fail` in place of `end` (`runOutputFailedOf(error)` is then true): the stream 3 that
arrived holds the warnings that the command line writes to its standard error before the error. The engine reads a `-logfile` name as the command line on
Linux reads it where the name's directories exist: the file's name is the whole name, which may
contain `/` (`logs/run.log`, `/tmp/run.log`, `/dev/stdout`, `logs/`, which the command line makes as
the file `logs`); show and save it as a name, not as a path. The standard error (`/dev/stderr`), a
directory (`.`, `logs/.`) and a name where a user's command line cannot make or open a file give
no file: their messages are on stream 3. A user makes files in the working directory, in a home
(`/home/<user>/...`) and where every user does (`/tmp`, `/var/tmp`, `/dev/shm`, `/run/lock`,
`/var/lock`, `/var/crash`); every other directory is root's (`/etc/run.log`, `/data/run.log`,
`/dev/run.log`, `//run.log`), and the system refuses a component of more than 255 bytes. A name
may climb out of a directory of root's (`/etc/default/../../tmp/run.log` is `/tmp/run.log`), but not out of one
that a user cannot enter (`/root/../tmp/run.log`). `/dev/fd/2` and
`/proc/self/fd/2` give no file and send the messages that go to the file nowhere, as `/dev/null`
does: the command line writes nothing to a log on a file system without free space, which these
are.

The search strategy file holds the command line's bytes for the argv. Its name is read as the
`-logfile` name is, as the command line opens it for writing: a name the command line cannot open
(`""`, a directory, `/etc/s.asn`) is the command line's `File is not accessible` error at once;
`/dev/null` and the devices whose writes are lost give no file; `-` and the standard output are the
command line's refusals, and the standard error, the other devices and files of the system, and the
file of `-query`, `-subject`, `-logfile` or `-window_masker_db` are LOSAT Web's, raised where the
command line writes its first request. A run that fails after the file is made keeps it, empty or cut
back as the command line leaves it.

`RuntimeInfo.messages` says where the run's messages went (the warnings that the command line writes
to its standard error without `-logfile`): `'diagnostics'` (stream 3), `'logfile'` (the file of
`-logfile`, where TBLASTX with `-logfile NAME` writes them) or `'nowhere'` (TBLASTX's `-logfile` of
a name whose messages go nowhere: `""`, `/dev/null`, `/dev/fd/2`, `/dev/full`...). The engine sends it once per run, before the outputs, on stream 10
(`MESSAGES_STREAM`): 4 bytes, the place as a little-endian unsigned integer (`messagesSignal` reads
it). An empty stream 3 means no warnings only where it is `'diagnostics'`; where it is `'logfile'`,
an empty file means none. These are the warnings that the command line writes for the run's own
`-outfmt`. What the run's other outputs add to them (their option warnings, and in a run of two
searches the other search's warnings and errors) comes on stream 3 where it is `'diagnostics'` or
`'logfile'`, and is not kept where it is `'nowhere'`.

In a worker or in Node, the serial module can also be called directly:

```ts
import { compileReactor, discardedStream, DISCARD_STREAM, instantiateSerial, ROLE_QUERY, ROLE_SUBJECT } from 'losat';

const { module } = await compileReactor(serialBytes);
const { abi } = await instantiateSerial(module);
const argv = ['blastp', '-query', 'query.faa', '-subject', 'subject.faa', '-num_threads', '1'];
// The verdict on the argv says whether its run reads the inputs with NCBI's data loaders on.
const validation = abi.validate(argv);
if (!validation.ok) throw new Error(validation.message);
const query = abi.register('blastp', ROLE_QUERY, validation.dataLoaders, queryFasta);
const subject = abi.register('blastp', ROLE_SUBJECT, validation.dataLoaders, subjectFasta);
let outfmt6: Uint8Array[] = [];
// -query and -subject name the inputs in the outputs; the sequences come from the handles.
abi.run(argv, query.handle, subject.handle, (stream, bytes) => {
  if (stream === 6) outfmt6.push(bytes.slice());
  // Stream 9: the engine could not write the output that it names; if it is outfmt 6, drop what arrived.
  else if (stream === DISCARD_STREAM && discardedStream(bytes) === 6) outfmt6 = [];
});
abi.release(query.handle);
abi.release(subject.handle);
```

The package's tests run this example on one of the repository's fixtures (`tests/package.test.ts`):
its outfmt 6 is byte for byte what NCBI BLAST+ 2.17.0 writes for the same search.

## License

LOSAT is a Rust translation of the search programs and database tools of NCBI BLAST+ 2.17.0. The code written by NCBI is in the public domain in the United States. LOSAT's own work (the changes and additions of the translation, threading, the WebAssembly builds, the tests, the Web application and the Python package) is under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LICENSE). Parts that come from other projects (LMDB, Boost, Arm Optimized Routines, Sun's fdlibm, Microsoft's C runtime, BitMagic, the Rust crates bzip2 and ruzstd) keep their own notices: see [THIRD_PARTY_NOTICES](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/THIRD_PARTY_NOTICES).
