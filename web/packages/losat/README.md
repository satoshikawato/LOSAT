# losat

LOSAT's local sequence search (BLASTN, BLASTP, BLASTX, PSIBLAST, TBLASTN and TBLASTX) for the browser.
Searches run on the user's machine, so sequences never leave it. The package holds two
WebAssembly builds of the [LOSAT](https://github.com/satoshikawato/LOSAT) engine, the same
engine as the `losat` command line, together with the TypeScript host that runs them. LOSAT Web
runs its searches through this package.

- `losat-web-serial.wasm` runs one thread and works in any browser and in Node.
- `losat-web-threads.wasm` runs a search on several threads. It needs memory that is shared
  between threads, which a browser allows only on a cross-origin isolated page (below).

The two builds speak the same ABI (ABI v2 at version 17, `docs/web/abi_v2.md` in the repository). A search
gets an argv, as on the command line, and FASTA bytes (or `-in_pssm`/`-in_msa` files in the run directory). It writes outfmt 0, 6 and 7 byte for byte
as the command line does, and it also writes an HSP record for each hit. This holds for a host that
drops the output that the engine tells it to discard (stream 9, in "Use").

The package is not published on npm. Build it from the [LOSAT repository](https://github.com/satoshikawato/LOSAT):
make the two modules with `web/adapter/tools/build_reactors.py`, then run
`LOSAT_WEB_REACTORS=<modules directory> npm run build` in `web/packages/losat`; `npm pack` there makes a
file to install (see the [WebAssembly manual page](https://satoshikawato.github.io/LOSAT/libraries/webassembly.html) and [CONTRIBUTING.md](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/CONTRIBUTING.md)).

## Versions

The package's version is the engine's version. `losat/abi.json` names both, together with the
ABI version that the reactors report (`losat_web2_abi_version`):

| package | engine | ABI |
|---|---|---|
| 0.2.0 | 0.2.0 | 17 |

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
`-conffile` names, the file that `-logfile` names (`StartInfo.logFile` and its normalized host
name `StartInfo.logFileName`), the search strategy file that
`-export_search_strategy` names (`StartInfo.strategyFile`), and whether the argv has `-dryrun`. Give the
file of `-conffile` under that name to `engine.validate(argv, registry)` and as
`EngineRunRequest.conffile`; `reactorControl` validates it on an instance of its own, opened with its
second argument. A run returns the `-logfile` file and the search strategy file in
`RuntimeInfo.runFiles`, each `RunFile` with the option that names it (`runFileOption(file)`:
`'logfile'` or `'export_search_strategy'`). When a run fails, `runFilesOf(error)` returns them. A run that the engine failed after it began to write ends
its output with `fail` in place of `end` (`runOutputFailedOf(error)` is then true): the stream 3 that
arrived holds the warnings that the command line writes to its standard error before the error. The engine reads a `-logfile` name as the command line on
Linux reads it in a working directory holding the run's files (the run directory's files, plus a
search's `-query`, `-subject`, `-conffile` and `-window_masker_db` names): a name in a folder that
is missing (`logs/run.log` without `logs`) or naming a folder gives no file and keeps its messages
on stream 3; a name that is an existing file of the run or a file the run opens or checks for reading even when that file is missing (`q.fa`, `.ncbirc`, `/tmp/q.fa`, a `-db` list file, a database lookup file of `-db`, `-filtering_db` or `windowmasker -infmt blastdb`, `taxdb.bti` when `-outfmt` consults the taxonomy database, an input of `windowmasker` including entries listed by `-fa_list true`, inside or outside the working directory) is refused, because LOSAT Web does not write to the run's files. When a file is
made, its name is the command line's file name (`run.log` for ` run.log` and `./run.log`, `logs` for
`logs/`, `sub/run.log` when `sub` is a folder of the run, `/tmp/run.log` outside the working
directory, and `/dev/stdout` for a search's `/dev/stdout`); show and save it as a name, not as a
path. The standard error (`/dev/stderr`), a
directory (`.`, `logs/.`) and a name where a user's command line cannot make or open a file give
no file: their messages are on stream 3. A user makes files in the working directory, in a home
(`/home/<user>/...`) and where every user does (`/tmp`, `/var/tmp`, `/dev/shm`, `/run/lock`,
`/var/lock`, `/var/crash`); every other directory is root's (`/etc/run.log`, `/data/run.log`,
`/dev/run.log`, `//run.log`), and the system refuses a component of more than 255 bytes. A name
may climb out of a directory of root's (`/etc/default/../../tmp/run.log` is `/tmp/run.log`), but not out of one
that a user cannot enter (`/root/../tmp/run.log`). `/dev/fd/2`, `/dev/fd/1` and
`/proc/self/fd/2` give no file and send the messages that go to the file nowhere, as `/dev/null`
does: the command line writes nothing to a log on a file system without free space, which these
are.

The search strategy file holds the command line's bytes for the argv. Its name is read as the
`-logfile` name is, as the command line opens it for writing: a name the command line cannot open
(`""`, a directory, `/etc/s.asn`) is the command line's `File is not accessible` error at once;
`/dev/null` and the devices whose writes are lost give no file; `-` and the standard output (where
the command line writes the requests; a run's standard output is its report), the standard error,
the other devices and files of the system, and the file of `-query`, `-subject`, `-logfile` or
`-window_masker_db` are LOSAT Web's refusals, raised where the command line writes its first
request. A run that fails after the file is made keeps it, empty or cut back as the command line
leaves it.

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

## Databases

A search of `-db` reads a BLAST database from the user's device. Give the run the database's files, and
the lists of its filters (`-seqidlist`, `-taxidlist`, `-gilist` and their negatives), as
`EngineRunRequest.files`: each a `DirectoryFile`, a `File` named by its path in the run's working
directory, as the command line run from that directory names it (`-db db/nt` reads `db/nt.nal`,
`db/nt.00.nin`, ...). The engine looks for the taxonomy files (`taxdb.bti`, `taxdb.btd`,
`taxonomy4blast.sqlite3`) as the command line does: at the directory's root, or in the directory that
`[BLAST] BLASTDB` of the run's `.ncbirc` or `-conffile` names (a registry file), not beside the
database; files named `db/taxdb.*` give the command line's `N/A` names. The directory is the engine's
current directory, `/work`: a relative name is a file of it, an absolute one (`/db/nt`) names none, and
a path is looked up as the OS looks it up (`nofolder/../db/nt` fails), as on the command line run from
that directory. Such a run has no `subject`, and `engine.start(argv)` says so (`database`). The engine opens the files by those names and reads each in place with `FileReaderSync` in the
workers: the bytes that it needs and, for a file read in order, up to 1 MiB ahead (a read of a file costs
about a millisecond whatever its size, so a file read in small pieces is read in fewer, larger ones). A
database is never copied whole into memory or into the browser's storage. Every thread of a threaded search reads the same files.

BLASTN's `-filtering_db` (a database of repeats that masks the query) runs the same way, with `-db` or
with a `-subject`: give the run the repeat database's files (every file of that database, as for `-db`)
in `EngineRunRequest.files` under the name of the option's value, and the engine searches it for repeats
as the command line does. LOSAT Web's form has no input for these files yet and refuses the option.

The limits are the browser's: a search must fit in WebAssembly's 4 GiB of memory, and a browser keeps
no access to the user's files across a reload, so the user gives them again. `-remote` and
`-entrez_query` are refused, as on the command line.

In Node, `instantiateSerial(module, directory)` takes a `RunDirectory`; give it a function that reads a
`Blob`'s bytes synchronously, since Node has no `FileReaderSync` (the package's tests read each Blob's
file with `fs.readSync`, `tests/database.test.ts`).

## Calling a module directly

In a worker or in Node, the serial module can also be called directly:

```ts
import { compileReactor, discardedStream, DISCARD_STREAM, instantiateSerial, ROLE_QUERY, ROLE_SUBJECT } from 'losat';

const { module } = await compileReactor(serialBytes);
const { abi } = await instantiateSerial(module);
const argv = ['blastp', '-query', 'query.faa', '-subject', 'subject.faa', '-num_threads', '1'];
// The verdict on the argv says whether its run reads the inputs with NCBI's data loaders on.
const validation = abi.validate(argv);
if (!validation.ok) throw new Error(validation.message);
const query = abi.register('blastp', ROLE_QUERY, validation.dataLoaders, validation.parseDeflines, queryFasta);
const subject = abi.register('blastp', ROLE_SUBJECT, validation.dataLoaders, validation.parseDeflines, subjectFasta);
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

## Running CLI tools (`windowmasker`, `blast_formatter`, `dustmasker`, `segmasker`, `convert2blastmask`)

`engine.runTool(argv, options)` runs the command line's browser-supported tools (`windowmasker`, `blast_formatter`,
`dustmasker`, `segmasker` and `convert2blastmask` here; `makeblastdb`, which writes its database into the run's directory, in "Running makeblastdb" below) on the user's
files, as the command line runs them in a folder that holds them:

- `windowmasker` (`losat windowmasker`: `-mk_counts` in each `-sformat`, `-convert`, and the masking of `-ustat` in each
  `-outfmt`, with their options, help and version), reading `-in`, `-ustat`, `-fa_list` and the files it lists, `-ids`,
  `-exclude_ids`, a BLAST database of `-infmt blastdb`, and `-conffile` from `options.directory`.
- `blast_formatter` (`losat blast_formatter`: re-formatting `-outfmt 11` BLAST archives via `-archive` into every output
  format that LOSAT's `blast_formatter` writes, with `-out`, `-max_target_seqs`, `-num_alignments`, `-parse_deflines`,
  help and version; `-rid` is refused), reading `-archive`, any local BLAST database referenced by a `-db` archive, and
  `-conffile` from `options.directory`.
- `dustmasker` and `segmasker` (`losat dustmasker`, `losat segmasker`: the masking of low-complexity regions of nucleotide
  and protein sequences in each `-outfmt`, with their options, help and version), reading `-in` (FASTA, or a BLAST database
  of `-infmt blastdb`) and `-conffile` from `options.directory`.
- `convert2blastmask` (`losat convert2blastmask`: the lower-case masks of FASTA as masking information in each `-outfmt`,
  with its options, help and version), reading `-in` and `-conffile` from `options.directory`.

The outputs of `-outfmt maskinfo_asn1_bin` and `maskinfo_asn1_text` are the files that `makeblastdb -mask_data` reads.

`argv[0]` is the tool name (`windowmasker`, `blast_formatter`, `dustmasker`, `segmasker` or `convert2blastmask`), then the command line's words. `options.directory` holds
the files that the words name, each a `DirectoryFile` that the engine reads in place, as for a search of `-db` (above). A
folder of the directory opens as it does on Linux, and a read of it fails: `-in <folder>` ends as the command line does,
with "Is a directory". The promise resolves with the command line's
exit status (`status`, whatever it is) and what the tool wrote: `stdout`, `stderr` and `files` (`-out FILE`, `-logfile FILE`,
in the order the run opened them, under the command line's file name: `x.counts` for `./x.counts`, `logs` for `-logfile logs/`,
`sub/x` when `sub` is a folder of the directory; `-logfile /dev/stdout` writes diagnostic messages into `stdout` in the command
line's order, and a `-logfile` naming an existing file or an input the run opens or checks for reading (including `windowmasker`'s `-fa_list true` entries and `-infmt blastdb` database lookup files, or `blast_formatter`'s `-archive`, or the `-in` of a masking tool), or a `-out` naming an existing file or a missing input that the command line reads after opening `-out` (`-in` of `-mk_counts`, including `-fa_list true` entries and `-infmt blastdb` lookup files, `-ustat` of masking, and `dustmasker`'s `-in`), is refused with status 1). The status is the command line's: an argument error and a refusal by LOSAT are status 1, NCBI's exceptions
(a missing file, a bad counts file) are status 3, a `-conffile` that is not there is status 255, and a registry file that is
malformed (a `-conffile`, or a `windowmasker.ini` / `blast_formatter.ini` in the directory, with a bad section or entry line) is status 2; their text
is on `stderr`, but with `-logfile FILE` NCBI's report is in that file, as it is on the command line. `stderr` arrives as
the tool writes it, not when the run ends, so `onStderr` can show a warning while the run goes on. A name that is not a tool
rejects. The run is one at a time with the searches; `options.signal` cancels it.

```ts
// Re-formatting a BLAST archive (-outfmt 11) with blast_formatter:
const formatted = await engine.runTool(
  ['blast_formatter', '-archive', 'results.asn', '-outfmt', '12'],
  { directory: [{ name: 'results.asn', file: archiveFile }] },
);
if (formatted.status !== 0) throw new Error(await formatted.stderr.text());
const jsonText = await formatted.stdout.text();
```

```ts
import { WasmEngine } from 'losat';

declare const engine: WasmEngine; // as in "Use"
declare const fasta: File; // the genome, from a file picker

const argv = ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-unit', '15', '-out', 'genome.counts', '-logfile', 'genome.log'];
const directory = [{ name: 'genome.fa', file: fasta }];
// The command line's check of the words, before any file is opened. It reads the files that the run reads
// (`-conffile`, `.ncbirc`, `windowmasker.ini`), so give it the run's directory.
const verdict = await engine.validateTool(argv, { directory });
if (!verdict.ok) throw new Error(verdict.message);

// Without sinks, each output is collected as a Blob (a list of the engine's chunks, never one joined array).
const { status, stderr, files } = await engine.runTool(argv, { directory });
const counts = files.find((file) => file.name === 'genome.counts')?.file; // a Blob
if (status !== 0) throw new Error(await stderr.text());
```

A tool's output may be larger than the browser's memory (human unit counts are 724 MB, a genome masked as FASTA is
3.1 GB), so the engine sends it in pieces as it writes it (a file's bytes in frames of up to 1 MiB), and a sink takes each
chunk instead of the host keeping it. A file that grows slowly, such as a `-logfile`'s messages, reaches the sink in a frame
that is sent when it fills, before the run writes elsewhere in that file or cuts it, and when the run ends: a small log may
arrive only then. The standard error is the exception: it arrives at each write.

```ts
declare const folder: FileSystemDirectoryHandle; // the user's choice, from showDirectoryPicker()

// One TextDecoder for the whole stream, with `stream: true`, so that a character that two chunks split is not garbled.
const decoder = new TextDecoder();
const result = await engine.runTool(
  ['windowmasker', '-ustat', 'genome.counts', '-in', 'genome.fa', '-outfmt', 'fasta', '-out', 'masked.fa'],
  {
    directory: [
      { name: 'genome.counts', file: counts! },
      { name: 'genome.fa', file: fasta },
    ],
    sinks: {
      onStderr: (chunk) => console.warn(decoder.decode(chunk, { stream: true })),
      // One sink for each file that the run opens (one call for a file that it opens twice, as `-out X -logfile X` does).
      // The host calls `write` for each chunk, in order, then `close`.
      onFile: async (name) => {
        // A name may hold `/`: `sub/x` is the file `x` of the folder `sub`.
        let parent = folder;
        const parts = name.split('/');
        for (const part of parts.slice(0, -1)) parent = await parent.getDirectoryHandle(part, { create: true });
        const writable = await (await parent.getFileHandle(parts.at(-1)!, { create: true })).createWritable();
        return {
          // The bytes go at `offset` of the file: an append in the usual case, but not always (see below).
          write: (chunk, offset) => writable.write({ type: 'write', position: offset, data: chunk }),
          // The file's length becomes `length`: the run empties a file when it opens it again for writing.
          truncate: (length) => writable.truncate(length),
          close: () => writable.close(),
          abort: (reason) => writable.abort(reason),
        };
      },
    },
  },
);
decoder.decode(); // the end of the standard error: the decoder flushes a character that was left incomplete
// The files and the standard error went to the sinks, so `result.files` is empty and `result.stderr` is an empty Blob;
// `result.status` is the exit status. A sink that returns undefined for a name has that file collected in `files`.
```

The host awaits each call of a sink before the next one (`write`, `truncate` and `close` may return a promise), and the chunks
are the sink's own. A run that opens one file twice (`-out X -logfile X`, or `-out ./X -logfile X`: the same file) writes
one file, as the command line does: each handle has its own position, so a write may come at an offset before the end of
the file, and one past the end leaves zero bytes between. A sink that writes at `offset` and cuts at `truncate` is exact for
that; the default collection, a Blob, is too. `onFile` is called once for each file that the run opens, also then. A cancel, a
failure of the engine or an error of a sink ends the run and calls `abort` of each sink whose file is still open, and the
promise rejects.

A sink slower than the engine does not make the page hold the output. Where the page can share memory with the Engine worker
(a cross-origin isolated page, which the threaded module needs too; Node), the worker stops the engine while more than 64 MiB
of the bytes that it posted wait for the sinks, and goes on as they take them. Where the page is not cross-origin isolated
(no `SharedArrayBuffer`), nothing can wait: the chunks queue on the page until the sink takes them, so a sink that is slower
than the engine holds what the engine has written and the sink has not taken. Serve the page isolated to bound it, or take the
output with a sink that keeps up. A file that no sink takes is collected as a Blob, and does not queue.

`engine.validateTool(argv, { directory })` checks the words as the command line does before the tool opens a file, and reads
the run's registry files (`-conffile`, `.ncbirc`, `windowmasker.ini`) as the command line does, so it takes the same
`directory` as `runTool`: it opens the files on the instance that the Engine worker has loaded (a threaded runtime that the
searches use keeps its threads, its retained subject and counts file; the serial one is loaded only when none is). It takes
its turn there like a run: it waits behind a run in progress. It does not make the engine busy: a run may start while a
check is pending, and a check may be asked during a run. When the worker that has it ends (a cancel, a renewal, a stop), the
check goes to the next worker; a check that three workers ended through no run's doing is given up. A check whose instance
stops ends the worker as a run's does: that check rejects, a run posted behind it on that worker fails with "Run the search
again.", and the next run starts a new worker.

`engine.describeTool(tool)` reads a tool's options, defaults, value lists, ranges and exclusions by running
`[tool, '-xmlhelp']` (on the control reactor when `reactorControl` provides it, else through `runTool`) and parsing the
XML with `parseToolDescription`. `parseToolDescription(xml)` is also exported from `losat` and `losat/domain` and needs
no `DOMParser`, so it works in workers and in Node.

## Running makeblastdb

`engine.runTool(['makeblastdb', ...words], options)` runs the command line's `makeblastdb` (`losat makeblastdb`, with every
option: nucleotide and protein, BLAST database versions 4 and 5, `-parse_seqids`, `-hash_index`, volumes of `-max_file_sz`,
`-mask_data`, `-gi_mask`, `-taxid`, `-taxid_map`, `-input_type` fasta, blastdb, asn1_txt and asn1_bin, gzip input, `-logfile`,
help and version) on the files of `options.directory`, and writes the same bytes as the command line. The run's directory
is the command line's working directory, and it is writable for a tool's run: the database's files go into it as the
command line writes them, the run reads them back, renames and removes files (a run that replaces a database of
`directory` removes the old volumes), and makes the folders of `-out sub/dir/db`. What it writes is kept in the Engine
worker's memory, and the user's files (the Blobs of `directory`) are never changed. When the run has ended, the files that
it made or changed and that are still there come to `sinks.onFile`, or to `files` as Blobs, in the byte order of their
names: a file that the run removed never arrives, nor does a file of `directory` that it left as it was. The standard
input is empty (`-in -` reads no sequence, as the command line with `< /dev/null`).

The dates that the database holds (the index files' "Nov 27, 2026  12:23 AM", an alias file's "Alias file created", the
mask column's and the report's dates) are the user's local time, as on the user's machine: the host gives the engine the
browser's time zone (`Intl`) as a POSIX `TZ` rule of the current year, with its summer time, made at the start of each
run, so a zone that the user changed since an earlier run is the new run's.

`-logfile` writes where the command line writes it: a file of `directory` (made or appended to, then a file of the run),
`/dev/stdout` and `/dev/stderr` among the run's standard output and error, `/dev/null` nowhere, and a file outside
`directory` (`/tmp/mk.log`) as a file of the run under that path; a name that does not open (a missing folder, a folder)
leaves the messages on the standard error, as on the command line.

```ts
import { WasmEngine } from 'losat';

declare const engine: WasmEngine; // as in "Use"
declare const proteins: File; // a FASTA file, from a file picker

const argv = ['makeblastdb', '-in', 'proteins.faa', '-dbtype', 'prot', '-parse_seqids', '-out', 'proteins'];
const directory = [{ name: 'proteins.faa', file: proteins }];
const { status, stdout, stderr, files } = await engine.runTool(argv, { directory });
if (status !== 0) throw new Error(await stderr.text());
console.log(await stdout.text()); // "Building a new DB, current time: ..." and the counts, as the command line prints them
// files: proteins.pdb, proteins.phr, proteins.pin, ... (a Blob each, by name). Each is a `DirectoryFile`, so a search
// of `-db proteins` takes them as they are (`EngineRunRequest.files`, in "Databases"), as does another tool's `directory`.
```

`engine.validateTool(argv, { directory })` makes the command line's checks of the words that come before its run starts:
an option that is not known, a missing `-dbtype`, a `-dbtype` that is not `nucl` or `prot`, and the registry files of the
directory. The errors that the command line reports once it has started (a missing `-in`, a mask file that cannot be read)
are the run's: `runTool` returns them as the command line does.

## License

LOSAT is a Rust translation of the search programs and database tools of NCBI BLAST+ 2.17.0. The code written by NCBI is in the public domain in the United States. LOSAT's own work (the changes and additions of the translation, threading, the WebAssembly builds, the tests, the Web application and the Python package) is under the [MIT License](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/LICENSE). Parts that come from other projects (LMDB, Boost, Arm Optimized Routines, Sun's fdlibm, Microsoft's C runtime, BitMagic, the Rust crates bzip2 and ruzstd) keep their own notices: see [THIRD_PARTY_NOTICES](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/THIRD_PARTY_NOTICES).
