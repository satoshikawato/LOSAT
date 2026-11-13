# LOSAT Web ABI v2 (engine adapter ↔ application)

Status: **version 2, implemented** in `web/adapter/` (session S05). Drafted in S01
([`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md) §4, §7). The application
depends on the TypeScript ports under `web/app/src/ports/` (the searches and the HSP record in
`engine.ts`, the index scan in `scan.ts`, the inflation in `inflate.ts`), not on the exported
functions. The HSP record (§8) and the JSON responses (§9) cross into the
application, so any change to them updates the port that types them in the same
commit, and an incompatible change increments the version (§10).

Programs: BLASTP, TBLASTN, BLASTN, TBLASTX and BLASTX (since session SX).

## 1. Purpose

The application (`web/app`) runs LOSAT searches through a Wasm module built from the
adapter crate `web/adapter`, which links the engine crate `LOSAT/` as a library.

- Search conditions are one CLI argv, parsed by the same clap parser as the CLI.
- One run writes every output format the program supports (listed by `describe`),
  the HSP records and the diagnostics from the same search result (§7), and the output of
  the argv's own `-outfmt` when it is another format that the CLI writes, or 6 or 7 with fields
  (§5, §7).
- All BLAST-defined values and all compatibility text come from the engine.

ABI v1 (`losat_web_*` in `LOSAT/src/web_api.rs`, used by gbdraw) stays in the same
module unchanged, apart from fail-fast defect fixes (plan TD-1). The two ABIs share no
state.

## 2. Artifacts and module kind

| Artifact | Target | Cargo features of `LOSAT` | Notes |
|---|---|---|---|
| `losat-web-serial.wasm` | `wasm32-wasip1` | none | Serial. The engine rejects `-num_threads` above 1 |
| `losat-web-threads.wasm` | `wasm32-wasip1-threads` | `parallel`, `wasm-threads` | Needs a cross-origin isolated page and shared memory |

Both are WASI preview1 **reactors**: the host calls `_initialize` once after
instantiation, then calls the exports below any number of times. The build must match
the certified command-WASI builds in profile, rustflags, link arguments and shared
dependency versions (plan TD-6). `web/adapter/tools/check_build_identity.py` checks
this, and `web/adapter/tools/build_reactors.py` builds both artifacts and records their
identity. The reactor start-up (`crt1-reactor.o`, `--entry=_initialize`) comes from
`LOSAT/build.rs` through the `LOSAT` dependency: Cargo passes that build script's
`rustc-cdylib-link-arg` to the adapter's cdylib as well, so the adapter has no build
script of its own (a second copy would define `_initialize` twice).

The reactors' bytes depend on where they are built. `LOSAT` is a path dependency
outside the adapter's workspace, so rustc gets its sources by absolute path, and every
crate from the registry is compiled from `CARGO_HOME`. `build_reactors.py` replaces both
prefixes in the source paths that the modules embed (panic locations) with `/losat` and
`/cargo` (`--remap-path-prefix`, plan TD-11), so a module carries no path of the build
machine. The checkout path still changes the bytes: Cargo's metadata hash of the path
dependency, which appears in the symbol names of the name section, includes it. Two
builds of the same commit from the same checkout path are identical; the identity record
(`losat-web-*.json`) records the checkout path, `CARGO_HOME` and the remapping.

## 3. Imports

| Module | Name | Used by | Purpose |
|---|---|---|---|
| `wasi_snapshot_preview1` | `fd_write`, `fd_close`, `fd_fdstat_get`, `clock_time_get`, `random_get`, `sched_yield`, `environ_get`, `environ_sizes_get`, `proc_exit`, and the filesystem functions `path_open`, `path_create_directory`, `path_filestat_get`, `fd_prestat_get`, `fd_prestat_dir_name` | both | stderr, clocks, random, environment (the host gives an empty environment: variables such as `BATCH_SIZE` change the CLI's results, whose parity is defined without them). The filesystem functions are linked by engine code that the exports do not use (for example the CLI's `-out` file, which `validate` rejects); the host gives no preopened directory, so such a call fails instead of touching files. In the threaded module, the shared-memory guard also calls `sched_yield` between the slices of a wait once the host ended the memory's waits (or while its table of waits is full) (`LOSAT/tests/wasi_shared_memory.js`, `endWaits` in `web/packages/losat/src/infra/reactor/shared-memory.ts`): WebKit ends a worker that the page ends only when its thread is in a host function |
| `wasi` | `thread-spawn` | threads | starts a rayon worker thread |
| `env` | `memory` | threads | shared linear memory, created by the host. The module declares a maximum of 65,535 pages (4 GiB less one 64 KiB page; plan TD-7, DW-28); the host uses that maximum or, where it cannot be reserved, a smaller one (below) |
| `losat_host` | `emit(stream: u32, ptr: u32, len: u32)` | both | receives output bytes (§5). Worker instances of the threaded module import it too, but never call it |

The shared memory's maximum. The threaded module is linked with `--max-memory=4294901760`
(`LOSAT/.cargo/config.toml` and its copy `web/adapter/.cargo/config.toml`). One page less
than wasm32's 4 GiB keeps the byte size of the memory within 32 bits, in the module's C
library and in the host. A host must not provide a memory with a larger maximum than the
declared one; it may provide a smaller one. An engine reserves the maximum as address
space when the memory is created. If it cannot, it either throws `RangeError`
(JavaScriptCore without a fast memory) or keeps a smaller reservation without an error,
down to the initial pages alone, so that `memory.grow` fails (V8 without guard regions,
as on Android or under an address-space limit). The hosts of this repository therefore
try 65,535, 49,152, 32,768, 24,576 and 16,384 pages: they create each candidate one page
short of the module's initial size, grow that page, and use the first candidate for
which both steps succeed. Below 16,384 pages (1 GiB, the maximum declared until plan
DW-28) a host does what it did before: LOSAT Web runs the search on the serial module,
and the Node host fails. The serial module defines its own memory without a maximum; that
memory is not reserved in advance and can grow to wasm32's 4 GiB as well.

## 4. Exports

All functions return `i32`: `0` (or a non-negative handle) on success, `-1` on failure.
On failure, `losat_web2_last_error_ptr/len` hold a UTF-8 message. Engine errors use the
CLI's wording: argv errors are the message of the engine's parser
(`LOSAT::web_adapter::render_message`), search errors are the error and its causes joined by
`: ` (as ABI v1 reports them). The adapter parses the argv exactly as the command line
(every program accepts the CLI's default `-outfmt 0`), so an argv error of `validate`, an
unknown program's included, is the CLI's message (session S08). Only the adapter's own
rules give other messages: `-out` is not accepted, and an `-outfmt` of a format other than 0, 6
and 7 is refused (§7). `run` finds the
handles of the argv's program before it parses the argv, so it reports an unknown program
as `unknown program '<name>'`; a host validates the argv first.

| Export | Arguments | Result |
|---|---|---|
| `losat_web2_abi_version()` | — | `2` |
| `losat_web2_alloc(len)` / `losat_web2_dealloc(ptr, len)` | — | memory for inputs (§6) |
| `losat_web2_describe(program_ptr, program_len)` | program name | emits a *describe* JSON on stream 2 |
| `losat_web2_validate(argv_ptr, argv_len)` | argv (§7) | `0` if the argv is valid, else `-1` with the CLI error. For BLASTN it also checks the scoring options as NCBI does before a search (the checks that do not depend on the query), and returns NCBI's message (`BLAST query/options error: …`), or LOSAT's message with `not supported by LOSAT's BLASTN` for options that NCBI runs and LOSAT does not (`docs/evidence/losat_web_e2c/AUTHORITY.md` §G). For TBLASTX it checks the options alone, with the messages of a run: NCBI's `-evalue` check (`BLAST query/options error: expect value or cutoff score must be greater than zero`). For BLASTX it checks the options of each search's formats (§7) as the run checks them before it searches (`blastx::check_options`), with the messages of a run: LOSAT Web's refusals, LOSAT's BLASTX limits (for example `-matrix IDENTITY` with the compressed lookup table of `-ungapped -word_size 5`, with which NCBI crashes) and NCBI's errors of the options (`BLAST query/options error: …`); as for BLASTN, the argv is valid when one search is. For every program it reads `-subject_loc` and `-query_loc` as NCBI's `ParseSequenceRange` does, in the command line's order (for BLASTN, BLASTP, TBLASTN and TBLASTX the subject range before the program's option checks, the query range among them; for BLASTX both after LOSAT's limits of its options, the subject range first, then the query range, then NCBI's errors of the options, as its run reads them), with NCBI's messages (`BLAST engine error: Invalid specification of query location (…)`, `subject location` for the subject) or LOSAT's for a part that NCBI cannot convert to an int (session S11). Whether a range fits the records is checked by `run`: a start past a subject record's end (NCBI's `Invalid from coordinate (greater than sequence length)`), a start just past a record's end (LOSAT's rejection; BLASTX searches that record as a sequence without data, as NCBI does, with NCBI's `Sequence contains no data` or `Subject sequence contains no data` warning), every query record of a batch skipped (NCBI's `Empty CBlastQueryVector`). The CLI reads the subjects before it checks the query options, so for an argv and inputs with both kinds of error the CLI reports the record's error and `validate` the argv's (`docs/evidence/losat_web_e2d/AUTHORITY.md` §A, §B). The checks of the inputs come with `register` and `run` |
| `losat_web2_register(program_ptr, program_len, role, bytes_ptr, bytes_len)` | role `0` query, `1` subject; original FASTA bytes | handle ≥ 1; emits a *register* JSON on stream 2. A handle belongs to the program that registered it. The records are those of the program's own reader, the engine's port of NCBI BLAST+'s `CBlastFastaInputSource` (session SF): nucleotide flags for BLASTN, TBLASTX, the TBLASTN subject and the BLASTX query, protein flags for BLASTP, the TBLASTN query and the BLASTX subject, the bytes read as a file, the data loaders on (as the CLI without `.ncbirc` or `DATA_LOADERS`). Every record that the reader returns is registered, records without residues and a first record without a defline included; the reader's warnings stay with the records and come on stream 3 of a run, where the CLI writes them. `register` first scans the bytes with the program's scan kind (1 or 2, §9) and fails with the scan's error: LOSAT Web's rejections (`not supported in the browser`: a `>?` gap line, which the CLI reads; a first line that NCBI's data loaders would fetch as a Seq-id; a record whose residues two joined lines hide from the index; a record over 2147483647 letters) and the reader's own error with the CLI's text (`BLAST query error: CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`). A query with blank and comment lines only (`!`, `#`, `;`) fails with NCBI's `BLAST engine error: Empty CBlastQueryVector`, which every run of it ends with. A query of white space only has no record: a run gives NCBI's `Query is Empty!` warning. A subject without records is registered: a run gives NCBI's `BLAST engine error: Empty CBlastQueryVector`. `register` then checks that the scan found the reader's records (§9). `register` judges one input alone, in the order the host calls it, and never looks at the other role's input. The CLI reads the subjects while it processes its arguments, before it reads any query: a subject that fails to read, and a subject without records (`BLAST engine error: Empty CBlastQueryVector`, exit 3), end the CLI before the query is read, in BLASTN, BLASTP, TBLASTN and TBLASTX (checked with NCBI BLAST+ 2.17.0; audit finding D-1 of session SFc), and in LOSAT's BLASTX (its command line and `run_local`). So for a query and a subject that both fail, the CLI reports the subject's error, while a host that registers the query first sees the query's. A host that wants the CLI's first error registers the subject first, stops at its error, and takes a subject registered without records as the run's `BLAST engine error: Empty CBlastQueryVector` before it registers the query |
| `losat_web2_release(handle)` | handle | `0` |
| `losat_web2_file_begin(kind, name_ptr, name_len, size)` / `losat_web2_file_chunk(file, ptr, len)` / `losat_web2_file_end(file)` / `losat_web2_file_release(file)` | *named input files* (line NP WW1, WW2): kind `1` = a WindowMasker unit counts file, the file that BLASTN's `-window_masker_db` names, or the file of `-window_masker_taxid <taxid>`; the name (UTF-8) is that option's value, or `windowmasker/<taxid>/wmasker.obinary` (the taxid in decimal: LOSAT Web's directory of `-window_masker_taxid` is `windowmasker`, where the CLI reads NCBI's registry entry `[WINDOW_MASKER] WINDOW_MASKER_PATH`, and the engine joins the names as NCBI does), compared byte for byte with the argv's or the joined one (no path is normalized: `./x` and `x` are two names); `size` is the file's length in bytes; then the file's bytes in order, in chunks of any size | `file_begin` returns a file handle ≥ 1, or fails for another kind (`unknown input file kind N`), an empty name, a name that is registered or loading (`an input file named N is registered; release it before registering another`: the host releases the old file first, nothing is replaced silently), or `size` bytes that cannot be had (`not enough memory for the WindowMasker counts N (S bytes)`; the bytes are reserved here, before any chunk). `file_chunk` fails, and drops the file, when the chunks pass `size`. `file_end` fails, and drops the file, when fewer than `size` bytes came; otherwise it registers the file and emits a *file* JSON on stream 2 (§9). It checks the bytes as NCBI's `DiscoverStatType` does and, for the two binary formats that BLASTN reads, reads the counts with BLASTN's settings (`s_BuildSeqMasker`), keeps both results and drops the bytes. An error of the contents (a format other than the two binary ones, a file too short for its tables) does not fail `file_end`: it is in the JSON's `error` and `stat_type`, and the run reports it as the CLI reports that file, with the CLI's message, exit status and place among the other errors. A registered file is seen by every later `run` of the instance, as a file of that name in the working directory of the CLI, until `file_release`, which frees its memory. A run whose `-window_masker_db` names no registered file fails as the CLI fails for a file it cannot open (`Error: … could not open N`, exit 255); a run of a taxid without a registered file is the CLI's run of a taxid without a file (each query batch reports `could not open ` with the empty name, exit 3). `validate` reads the two options where the CLI reads them (after `-dust`, before the query range and the other options: the check of the file of `-window_masker_db`, "at most one of" the options, the file of the taxid), so a host that registers the files before it validates gets the CLI's first error; for a `-window_masker_db` name that is not registered `validate` leaves the file to `run` (a host may validate without the file). A run without the options is unchanged |
| `losat_web2_scan_begin(parser)` / `losat_web2_scan_chunk(scanner, ptr, len)` / `losat_web2_scan_end(scanner)` | parser kind: `1` NCBI BLAST+'s reader with nucleotide flags (BLASTN, TBLASTX, the TBLASTN subject, the BLASTX query), `2` with protein flags (BLASTP, the TBLASTN query, the BLASTX subject) (kind `0`, the `bio::io::fasta` 1.6.0 reader of the four programs before session SF, was removed after the application switched to kinds 1 and 2, plan DW-38); FASTA bytes in chunks of any size | `scan_begin` returns a scanner handle (`unknown or unavailable FASTA parser kind N` for another kind); `scan_end` emits a *scan* JSON on stream 2 (§9), or fails with the reader's error or LOSAT Web's rejection |
| `losat_web2_inflate_begin(name_ptr, name_len)` / `losat_web2_inflate_chunk(inflation, ptr, len)` / `losat_web2_inflate_end(inflation)` | the file's name (UTF-8); the file's bytes in chunks of any size | `inflate_begin` reads the name with the engine's rule, the one by which the command line inflates `-query`, `-subject` and makeblastdb's `-in` (`LOSAT::web_adapter::compression_of_name`: the name ends in `.gz`, `.bz2` or `.zst`, in any case). A name that ends in none of them returns `0`, no handle and no response: the file is read as it is. Otherwise it returns an inflation handle and emits an *inflate* JSON on stream 2 (§9). `inflate_chunk` inflates the next bytes, and `inflate_end` ends the file; both emit the inflated bytes on stream 5 (§5) with the command line's decoder (`PushInflater`, the push form of the reader `Compression::decoder`; swebok-improve SW8 DECISIONS 44), so a file gives the bytes that the command line reads, however it is cut into chunks. A file that is not what its name says, is cut short or is empty fails `inflate_chunk` or `inflate_end` with the command line's sentence without its `-query: ` (`Compression::not_data`: `'x.fa.zst' is named .zst but is not zstd data (a zstd frame's checksum does not match its data); LOSAT reads a file named .zst through zstd`), where the decoder meets the fault, and the host discards what stream 5 gave for the file. The handle is dropped after an error and by `inflate_end`; a host that stops before the end of the file calls `inflate_end` and ignores its result. The application refuses a `.bz2` name itself (line GZ decision 11), after `inflate_begin`. A Wasm linear memory never shrinks, and the decoders grow it (a zstd frame's window is reserved at its first frame, up to 128 MiB), so the application runs each file on an instance of its own, opened for that file from the compiled module and dropped when the inflation ends, fails or is cancelled; the session's shared instance only reads the name rule (`inflate_begin`, then `inflate_end` with nothing pushed) |
| `losat_web2_run(argv_ptr, argv_len, query_handle, subject_handle)` | argv, handles registered for the argv's program | emits the program's supported format streams (0, 6, 7), stream 4 when the argv's `-outfmt` is another format or 6 or 7 with fields (§7), stream 1 (§8) and stream 3; returns after the run ends. It searches once for each distinct hit list size of the formats (§7): one search, or two when `-num_descriptions` / `-num_alignments` give format 0 another size than formats 6 and 7. The registered records enter the engine's `run_local` as the reader's records, so a run reads them as the CLI reads the same file: the subjects' warnings first, then the query batches with their warnings; a first batch that fails (for example a batch of queries without residues only, NCBI's `BLAST engine error: Warning: Sequence contains no data ...`) fails the run after the CLI would have written the outfmt 0 prolog to its file, and the host discards what it received (§5). Subject titles are written by the report layer ported in E2h (`NStr::HtmlDecode`, `GuessEncoding`, `x_CleanAndCompress`), so a title that NCBI decodes no longer fails a run (BLASTN, TBLASTX and TBLASTN outfmt 0/6 byte-identical to NCBI 2.17.0, SFc 2026-10-10); a run fails as the CLI does only for the CLI's explicit rejections that happen after the search starts (`docs/evidence/losat_web_e2h/AUTHORITY.md` §J: RP-20 in the tabular formats, BLASTP's punctuation-only subject titles in outfmt 0) and for NCBI's own errors (for example `The average subject length is too short` when every subject is without residues) |
| `losat_web2_last_error_ptr()` / `losat_web2_last_error_len()` | — | last error message |

`validate` has no query, so it does not check the Karlin-Altschul tables of BLASTN (scores and gap costs without a row,
LOSAT's limit of greedy gap costs). A CLI run meets them only when some query context has an ungapped
Karlin-Altschul block, which depends on the query's letters; `run` makes them and stops with the CLI's message and
exit code (`BLAST engine error: Error: …`), or writes NCBI's warning `Could not calculate ungapped Karlin-Altschul
parameters …` and no hit. What `validate` rejects, the CLI rejects.

## After a trap

The release profile sets `panic = "abort"` (`LOSAT/Cargo.toml:141-145`). The adapter has
a copy of it (`web/adapter/Cargo.toml:27-32`, checked by
`web/adapter/tools/check_build_identity.py`), and the adapter builds the shipped
reactors. So a Rust panic inside any export of this ABI or of ABI v1
([`abi_v1.md`](abi_v1.md)) does not unwind and does not return `-1`. It ends as a Wasm
trap, and the host sees an exception (in JavaScript, a `WebAssembly.RuntimeError`). The
same holds for an allocation failure and for a pointer outside the memory (§6).

After a trap:

- Do not call the instance again, not even `losat_web2_dealloc` or
  `losat_web2_last_error_ptr`. Its linear memory is in an unknown state.
- Every handle (`register`, `scan_begin`, `inflate_begin`), the stored records and the last error are
  lost. There is no error message to read for the call that trapped. The chunks of
  streams that the host received before the trap belong to a run that did not finish:
  the host discards them (§5).
- The host discards the instance and instantiates the module again, calls
  `_initialize` (§2), and registers its inputs again. The new handles differ from the
  old ones.
- In the threaded module the main instance and the worker instances share one memory
  that the host created (§3). The host discards the instance, the shared memory and the
  workers together, and creates new ones. The old shared memory is not reused: its
  contents are in an unknown state.

A failure that is reported (`-1` and a message) is not a trap. The instance stays usable
after it.

## 5. Output streams

`emit` may be called many times per stream. Chunks of one stream arrive in order. The
bytes are valid only during the call, so the host must copy them. `emit` is always
called on the thread that called the export: the engine formats on the calling thread,
which runs slot zero of the search's thread pool (`LOSAT/src/utils/threading.rs`). `emit` must
not call an export of the instance, and an exception thrown by `emit` leaves the instance as
after a trap (§4, "After a trap").

| Stream | Content |
|---|---|
| `0` | outfmt 0 text (only for programs that support outfmt 0) |
| `6` | outfmt 6 text |
| `7` | outfmt 7 text (only for programs that support outfmt 7) |
| `4` | the output of the argv's `-outfmt` when it is a format besides 0, 6 and 7 that the CLI writes (5, 8 to 16, 18, 20, BLASTN's 17; with or without words after the number), or 6 or 7 with words after the number (fields, `delim=`, §7): the CLI's output for the same argv, text, or binary ASN.1 for outfmt 9. For 13 and 14 it is a ZIP archive of the files that the CLI writes with `-out results`: the side files `results_<n>.json` (13) or `results_<n>.xml` (14), one per query, each added as its report is written, then the master file that lists them, as `results.json` or `results.xml`, then the central directory. The entries are stored (no compression) with their CRC-32 and the date 1980-01-01 00:00, so that the archive of a run is always the same bytes; ZIP64 records come only where the classic fields cannot hold a value (65,535 entries or more, or an offset of 4 GiB - 1 or more). Only such a run writes it. An output without bytes (outfmt 6 of a search without hits) comes as no chunk at all, so a host tells from the argv, as the adapter does, whether the run has this output |
| `1` | HSP records as JSON Lines, one object per HSP (§8); every program |
| `2` | JSON response of `describe`, `register`, `scan_end` or `inflate_begin` |
| `3` | diagnostics that the CLI writes to stderr (warnings), UTF-8, each written once: the CLI's standard error for `-outfmt 6` (that of `-outfmt 7` is the same, and that of `-outfmt 0` has the same lines without the formatting handler's `-num_descriptions` and `Examining 5 or more matches` lines, which `-num_descriptions` / `-num_alignments` can make differ, and without its `-line_length`, `-sorthits` and `-sorthsps` lines, which it writes for formats above 4 or other than 0 only). With stream 4 it is the CLI's standard error for the argv's `-outfmt` instead: for 6, 7, 10, 18 and 20 the same lines, and the warnings of the fields' writer (NCBI's `Taxonomy name lookup from taxid requires installation of taxdb database`); for the structured formats (5, 8, 9, 11 to 17) without the queries' warnings, which their writers do not write; as the CLI writes them. In a run of two searches in which one failed before it wrote anything (§7), it is that standard error followed by the CLI's standard error for `-outfmt 0`, error included |
| `5` | the inflated bytes of `inflate_chunk` and `inflate_end` (§4): each call sends what the decoder released during it, in chunks of at most 1 MiB, in order; the host copies them during the call |

The adapter sends each output stream in chunks of 1 MiB as they fill, and the rest when
the run ends, so an output never has to fit in linear memory as a whole. Stream 2 is
chunked the same way: a large *scan* or *register* response (thousands of records)
arrives in several chunks, which the host joins before it parses the JSON (S09). Stream 3 and
the HSP records are sent after the search. If `run` fails, the host discards what it
received.

## 6. Memory ownership

The host allocates input buffers with `losat_web2_alloc`, writes the bytes, passes the
pointer, and frees the buffer with `losat_web2_dealloc` after the call returns. The
adapter never keeps a pointer to a host buffer after a call: `register` copies what it
keeps, and so does `file_chunk` (one buffer can serve every chunk of a file).

A named input file holds its size in linear memory from `file_begin` to `file_end` (the bytes),
then about the size again while its counts are read (the tables), and keeps the tables until
`file_release` (the bytes are freed at `file_end`). The memory of a Wasm instance does not
shrink: the freed bytes are reused by the searches. A host that keeps counts files across runs
keeps at most what its memory allows; the application keeps one per instance (plan §4.6 R1).

Pointers and lengths are unsigned 32-bit values. JavaScript receives a Wasm `i32` (a
result of an export, an argument of an import) as a signed number, so a pointer of 2 GiB
or more arrives negative. A JavaScript host converts every pointer and length to its
unsigned value (`ptr >>> 0`) before it uses it as an offset into the memory: those of
this ABI (`losat_web2_alloc`, `losat_web2_last_error_ptr`, `emit`) and the arguments of
the WASI imports.

## 7. argv

- UTF-8 words separated by NUL (`\0`). The first word is the program name
  (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`).
- `-query <name>` and `-subject <name>` are required. The names are used only in
  output text (for example the `# Database:` line of outfmt 7); the sequences come from
  the handles.
- `-out` is rejected: the adapter writes every format the program supports (`describe` lists
  them).
- `-outfmt` (once, as on the command line; `-outfmt=value` too) is read as the CLI reads it, with
  the CLI's messages: a value that is not a format (`'abc' is not a valid output format`,
  `Formatting choice is out of range`), a format that LOSAT's program does not write, NCBI's
  formatting handler's errors (another program's format such as SAM, a custom delimiter that a
  field uses to separate its entries). Outfmt 13 and 14 do not need the CLI's `-out`: the run gives
  the engine a receiver of their side files (`OutputSink::Files`), and stream 4 is the archive of
  their files (§5). The value 0, 6 or 7 without
  words after the number adds nothing: the run writes those formats anyway; nor does 0 with words
  (`0 qseqid`), which the CLI reads as 0. Every other format that the CLI writes, with or without
  words, and 6 and 7 with words after the number (after NCBI's white space around the value is
  removed, a space follows the number: fields, `-field`, `std`, `delim=`) add stream 4, written by
  the tabular search (the first search, below) as its first format, so that it is the CLI's output
  for the same argv and adds no search (the hit list size of every format above 4 is that of 6 and
  7). What the CLI's run of that format has besides, the run has too (session W11): SAM's `@PG`
  line shows the command line (`CL:losat ` and the run's argv with the host's `-num_threads`, each
  word followed by a space, as NCBI's `GetCmdlineArgs`; the CLI writes its own argv there); the BLAST archive's Bioseqs hold
  the whole records (with `-query_loc` and `-subject_loc` too) and the deflines as read (`register`
  keeps a defline that `>` and the title do not give). A custom delimiter that a field uses to separate its entries
  (`6 delim=; qseqid staxids`) is NCBI's formatting handler's error, which the CLI makes after the
  ranges and the query options: `validate` reports the error of the CLI's own checks of the argv
  with that `-outfmt` (a `-subject_loc` of bad syntax first), whatever the hit list sizes of the
  searches, and `run` fails with it, never tolerating the failure of the tabular search (below).
  TBLASTN's and TBLASTX's engine reads every format of a run before the ranges and the option
  checks that come before the handler (`-seg`, `-query_loc`), so their `run` of an argv with both
  errors reports the delimiter's, where the CLI and `validate` report the other.
  The formats that the handler rejects (SAM of another program, 19, 21) are rejected with the
  handler's message as the argv is parsed, before the ranges. The rule reads the words, not the fields that they give, so that a
  host reads it from the argv without the table of fields: `6 std` adds a stream 4 with the bytes of
  stream 6. The engine is given the run's formats, and the checks of the options read them as
  without `-outfmt`: streams 0, 6, 7 and 1 are those of the argv without `-outfmt` (sessions W10,
  W11).
- `-num_descriptions` and `-num_alignments` are accepted (every program; `-max_target_seqs` excludes both, as on
  the command line) and described (§9). NCBI derives the hit list size of a format from them
  (`blast_args.cpp:2894-2973`): format 0 keeps the larger of `-num_descriptions` (500) and `-num_alignments` (250) where
  either is given, else `-max_target_seqs`; formats 6 and 7 keep `-max_target_seqs`, else `-num_alignments`, else 500.
  The sizes are equal without the options, with `-max_target_seqs`, and for example with `-num_descriptions 2
  -num_alignments 4`; they differ for `-num_alignments 10` and `-num_descriptions 10`.
- `-num_threads` is appended by the host. When the host falls back to the serial
  module, it passes `-num_threads 1` and records the reason (plan §4.7). The host
  validates the argv that it will run, with its `-num_threads`, so a `-num_threads` in
  the user's words is reported by `validate` (the parser rejects the repeated option).
- Every other word is parsed by the engine's CLI parser, so options and errors behave
  exactly as on the command line. LOSAT's own debug output (environment variables such
  as `LOSAT_TIMING`) goes to the module's WASI standard error, not to stream 3; stream 3
  carries the warnings that the CLI writes.
- Options are resolved separately for each output format, as the CLI resolves its single `-outfmt`. `run` makes
  one search for each distinct hit list size of the formats 0, 6 and 7 (the engine's `run_local` writes the formats of
  one hit list size from one search, and refuses formats whose search options differ otherwise): the tabular
  formats first, then format 0 if its size differs. Without `-num_descriptions` and `-num_alignments` this is the one
  search of today. Each stream is the CLI's output for that format alone with the same options, because the
  hit list size takes part in the search (a smaller list is not a prefix of a larger one: the preliminary list evicts
  subjects before their alignments are extended). With two searches the run costs about twice the time of one. Its
  memory peak is that of the larger search plus the HSP list of the tabular formats, which stays in memory (about the size of
  stream 1) while the second search runs (measured: +21.7 MiB on 142.6 MiB for BLASTN of Sakai against EDL933 with
  `-num_alignments 3`).
- A search of a two-search run that fails before it has written a byte of its formats does not fail the run when the
  other search succeeds, but only when the options' checks reject its hit list size, as the CLI rejects that format
  alone: for BLASTN and BLASTX the checks of `validate` for that search's formats fail (the adapter's `blastn_group_check`:
  `check_scoring_options` and `check_losat_limits`; `blastx::check_options`), for the other programs the search's hit list size is 0. This
  is `-num_alignments 0` for the tabular formats (NCBI's error `BLAST query/options error: No hits are being saved`
  for `-outfmt 6` and `-outfmt 7`, while `-outfmt 0` keeps its size, 500 or `-num_descriptions`), and for BLASTN's format
  0 a `-num_descriptions` of 2,147,483,598 or more (N + 50 overflows the 32-bit preliminary hit list size; LOSAT's
  error, where NCBI crashes; the CLI accepts up to 2,147,483,597), while the tabular formats ignore
  `-num_descriptions`. Any other error of a search, for example a thread pool that cannot start, fails the run, as it
  did when the search was alone. The streams of the failed search's formats stay
  empty; if it is the tabular search there are no HSP records (stream 1), else the records have no `out0`. Stream 3 is
  the CLI's standard error for the tabular formats (the tabular search's warnings, and NCBI's error when it failed)
  followed by the CLI's standard error for `-outfmt 0` (the format 0 search's warnings, and its error when it failed):
  the lines that the CLI user sees running the two formats one after the other, none altered or dropped. When both searches
  fail, `run` fails with the tabular search's error, as the CLI fails every format (so does `-num_descriptions 0
  -num_alignments 0`, whose hit lists are both empty and make one search); so does a search that fails after it has
  written a byte. `validate` accepts the argv when at least one search passes the checks, and reports the tabular
  search's error otherwise.

## 8. HSP record

One JSON object per HSP of the final, sorted result. Field names follow the Rust
`common::Hit` and `report::PairwiseHit` fields. Every program emits them (TBLASTX since
session S08, BLASTX since SX). BLASTN rows are nucleotide rows: the query on its
plus strand, residues of masked regions in lowercase, no frames. TBLASTX rows are the
translated rows of outfmt 0 (both sequences translated again from the nucleotides with
`-query_gencode` and `-db_gencode`, SEG-masked query residues in lowercase), with both
frames. BLASTX rows are the rows of outfmt 0: the query translated in its frame with
`-query_gencode` (`query_frame` 1 to 3 or -1 to -3; residues of the query's masked
regions in lowercase) and the protein subject (`subject_frame` null); with `-query_loc` and
`-subject_loc` the coordinates and `query_frame` are those of the whole records, as outfmt 6
prints them.

| Field | Type | Meaning |
|---|---|---|
| `index` | integer | 0-based position of the HSP in the final, sorted hit list of the formats 6 and 7 (`HspIndex` in `LOSAT/src/api/local_blast.rs`); the same HSP has the same index in every output format of the same search |
| `q_idx`, `s_idx` | integer | 0-based record index of the query and the subject in the registered inputs (the `index` of the *register* response: records without residues count; with `-query_loc`, a query record that NCBI skips because its interval starts past its end keeps its place, so `q_idx` is not the position among the searched queries) |
| `rank` | integer | 0-based position of the HSP among the HSPs of its query (derived from `index`) |
| `raw_score`, `bit_score`, `e_value` | number | engine values, not rounded |
| `q_start`, `q_end`, `s_start`, `s_end` | integer | 1-based coordinates as printed in outfmt 6 (start > end means minus strand; a BLASTN HSP of one letter has start = end on either strand, and only its `Strand=` line in the stream 0 section shows the strand, S07+ AUTHORITY §P) |
| `query_frame`, `subject_frame` | integer or null | translation frames where applicable |
| `subject_length` | integer or null | subject length |
| `query_aligned`, `subject_aligned` | string or null | aligned sequences with `-` for gaps |
| `out6` | [integer, integer] or null | byte range [start, end) of this HSP's row in the stream 6 text; null if outfmt 6 does not show it |
| `out0` | [integer, integer] or null | byte range [start, end) of this HSP's section in the stream 0 text: its score lines and alignment. The subject heading that precedes the first HSP of each subject is not part of any section; null if outfmt 0 does not show the HSP (for example BLASTX shows alignments for the first 250 subjects by default). In a run of two searches (§7) outfmt 0 has its own hit list: the range is set only where outfmt 0 shows exactly the HSPs that the records hold for that query and subject (the same count and every field of each HSP equal, in order; the subject has a heading in stream 0), else null (for example when outfmt 0 shows only some of the pair's HSPs, or other HSPs besides), and stream 0 may hold sections and headings that no record points to |
| `out0_subject` | [integer, integer] or null | byte range [start, end) of the heading of this HSP's subject in the stream 0 text: the defline and the `Length=` line that NCBI's `x_ShowAlnvecInfo` writes before the first HSP of each subject (`showalign.cpp:3613-3632`). Every HSP of the subject in that query has the same range |

The ranges come from formatter observer events: a formatter reports when it starts and
ends the row or section of an HSP, identified by its `index`, and the heading of a
subject, identified by the index of its first HSP; the adapter records the byte
positions at those moments. The formatter's output bytes do not
change (plan TD-3, §4.5). Displayed numbers are taken from the outfmt 6 row, not
formatted from the raw values (plan §4.4).

## 9. JSON responses (stream 2)

- *describe*: `{ "program", "formats": [0, 6, 7], "parameters": [{ "flag", "help", "takes_value", "default"?, "choices"? }], "query_gencodes"?, "subject_gencodes"? }`, generated from the clap definitions (without `-out`, `-outfmt` and `-num_threads`, which the adapter and the host set or refuse, or read from the argv's words as the CLI does: `-outfmt`, §7; `-num_descriptions` and `-num_alignments` are listed with the CLI's help text, which carries the defaults 500 and 250), the program's supported output formats and the genetic codes that the engine's parser accepts for `-query_gencode` and `-db_gencode`, and `"defaults"`.
  - `"defaults": { "tasks": { "<task>": { "<flag>": "<value>" } }, "by": { "<flag>": { "<value>": { "<task>": { "<flag>": "<value>" } } } }, "unless": { "<flag>": ["<flag>"] }, "given": ["<flag>"] }` are the values that the engine resolves for the options a search leaves out (`web/adapter/src/defaults.rs`), from the engine's own resolution of an argv (BLASTP's, TBLASTN's and BLASTX's `resolve`, BLASTN's `coordination` functions, TBLASTX's getters), printed as the CLI's help prints them (`"10"`, `"19.3"`, `"20 64 1"`, `"no"`, `"true"`). `tasks`: each task's values with no other option given (the key `""` for TBLASTX, which has no `-task`). `by`: for `-matrix` and `-word_size`, for each value that the engine resolves (the nine matrices of `ScoringMatrix::ALL`, word sizes 1-32), the values by task of the options that some value changes (a matrix's gap costs, threshold and window; the threshold of a given word size); a program lists only the options that change something (BLASTN neither, BLASTX the word size). `unless`: an option's default holds only while none of these options is given (`-max_target_seqs` excludes `-num_descriptions` and `-num_alignments`, §7). `given`: the options that BLAST+ tells given from not given, with the tasks where it does, so that a host writes one that is set even at its default: the hit list options in every task (§7), and each option that, given at the value shown for it (the task's, or the table's value for a given `-matrix` or `-word_size`), resolves another search (the word size of blastp-fast, tblastn-fast and blastx-fast, whose given 5 sets NCBI's threshold 19.3 where the task's own 5 keeps 20). The response of each program is computed once (some hundred probe argvs) and kept. An option without an entry has the clap `default` of `parameters`, or no default.
- *register*: `{ "handle", "records": [{ "index", "id", "length" }] }`, from the records of the program's own reader (§4). `id` is the record's title up to its first space (0x20), empty when the record has no title (never NCBI's `Query_N` or `Subject_N`); bytes that are not UTF-8 become U+FFFD in the JSON only (the output streams keep the title's bytes). `length` is the number of residues the reader stores. `register` also scans the same bytes with the program's kind (1 or 2) and fails if the records, IDs, lengths or residue counts differ.
- *scan*: `{ "records": [{ "index", "id", "header_offset", "sequence_offset", "end_offset", "length", "line_layout", "residue_counts" }] }`. Offsets are byte offsets in the scanned input and may exceed 2³²; they are encoded as JSON numbers and stay below 2⁵³. `line_layout` is `{ "kind": "uniform", "width", "eol" }` when residue `i` is at `sequence_offset + floor(i / width) * (width + eol) + i mod width`, or `{ "kind": "checkpoints", "every": 65536, "offsets": [...] }`, where `offsets[k]` is the byte offset of residue `k * 65536` and the residues are read forward from there under the kind's rules (below). `residue_counts` maps each byte value of the stored residues to its count: the bytes 0x21-0x7E as the character itself, every other byte (including the space) as `"0xNN"`.
  - Kinds 1 and 2 (session SF) give exactly the records that NCBI BLAST+'s reader returns from the bytes read as a file, with the program's flags and the data loaders on: records without residues and a first record without a defline included. The scan writes no messages. An input of white space, blank and comment lines only has no record and scans to `{"records":[]}`.
    - `id`: the title up to its first space (0x20), as in *register*. The title is the defline after `>` (after `>?_`, which NCBI reads as `>`) and the white space that follows `>`, up to the first byte below 0x20 after its first byte, without trailing white space; empty without a title; bytes that are not UTF-8 become U+FFFD.
    - `header_offset`: the defline's `>`. A first record without a defline has `header_offset` = `sequence_offset` = 0, so the extracted bytes of the record start the input as before.
    - `sequence_offset`: the first byte that NCBI's line reader reads after the defline's end of line (CR, LF or CR LF; a CR or LF that the line reader drops right after it is skipped), or the end of the input. `end_offset`: the next record's `>`, or the end of the input (a tail that NCBI's line reader loses is inside the last record's range).
    - `length` and `residue_counts`: the stored residues, upper-cased. Kind 1 stores `A B C D G H K M N R S T U V W Y` in either case and counts `U` as `T` (keys `A B C D G H K M N R S T V W Y`); kind 2 stores every ASCII letter and `*` (keys `A`-`Z`, `*`). Nothing else is counted: no `-`, digits, `*` in kind 1, invalid letters, white space or other bytes.
    - Uniform layout: `width` is the number of residues before the first jump, `eol` the number of bytes between the first line's last residue and the next residue (any positive count, not only 1 or 2); a record on one line has `width` = `length` and `eol` 1, a record without residues `width` 0 and `eol` 1. The layout is uniform only when the formula holds for every residue.
    - Reading forward from a checkpoint (its byte is a residue): read the bytes in order; CR and LF end a line (CR LF is two ends, harmless); at the start of a line skip space, tab, VT (0x0B) and FF (0x0C), and skip the whole line when its first other byte is `!`, `#` or `;`; elsewhere `;` skips the rest of the line; a byte that the kind stores (kind 1: the 16 letters above in either case; kind 2: `A`-`Z`, `a`-`z`, `*`) is the next residue (its value is the byte upper-cased, `U` read as `T` in kind 1); every other byte is skipped (white space, `-`, `>`, `?`, `_`, digits, other ASCII, bytes 0x80-0xFF, a byte order mark). A record ends at `end_offset`; read at most `length` residues.
    - `scan_end` fails, in the reader's order, with the reader's error `CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` (NCBI's text without the CLI's `BLAST query error: `; N is NCBI's line number: every line counts, CR LF is one line end), or with LOSAT Web's explicit rejection (`not supported in the browser`) of a first line that NCBI's data loaders would fetch as a Seq-id (start the input with a `>` defline), of a `>?` gap line (its residues have no bytes to locate; maintainer decision 3 of session SF: the CLI reads gap lines, LOSAT Web does not), of a record over 2147483647 letters, and of a record that is not uniform and whose residues NCBI's line reader joins from two lines of the file (a CR line end in a file of LF line ends, or an LF in a file of CR line ends), which the read-forward rules cannot follow (use one kind of line end). The rejection names the line (N 1-based) or the record (`record N ("id")`).
  - The scan exists only to index original records for extraction; the search always reads the input with the program's own reader (plan TD-8). `web/adapter/tests/scan_ncbi_properties.rs` checks kinds 1 and 2 against the engine's reader and an independent reference walk, in chunks of any size.
- *file*: `{ "handle", "name", "kind": "windowmasker_counts", "size", "sha256", "stat_type", "unit_size", "error" }` (`file_end`, §4). `sha256` is the lowercase hex SHA-256 of the bytes given, counted by the adapter over the chunks. `stat_type` is what NCBI's format check found: `"obinary"`, `"binary"`, `"ascii"`, `"oascii"` or `"unknown"`, or null when the check failed. `unit_size` is the counts' unit size when they were read, else null. `error` is null, or the message of the check's or the reader's error without NCBI's framing (for example `not enough data to fill the hash table`); a `stat_type` other than `"obinary"` and `"binary"` with a null `error` is a format that BLASTN refuses (`Only optimized binary windowmasker stat files are supported`). The run gives these errors with the CLI's framing (§4).
- *inflate*: `{ "compression", "plain_name" }` of `inflate_begin`: `compression` is `"gzip"`, `"bzip2"` or `"zstd"`, the format that the name's suffix announces; `plain_name` is the name without the suffix (`a.fa` for `a.fa.GZ`), or `input` when nothing is left (`.gz`), the name under which the application keeps the inflated file. `web/adapter/tests/inflate.rs` checks the three formats whole, joined, damaged, cut and empty, in chunks of one byte and whole, against the command line's reader and sentences.

## 10. Versioning

`losat_web2_abi_version` returns `2`. Any incompatible change to a signature, a stream,
or a JSON field increments the version, and the application refuses a module whose
version it does not know.

Stream 4 and the `-outfmt` of §7 (session W10, 2026-11-05) kept version 2: a host that does not
know them never sends `-outfmt` (an earlier adapter refused it) and drops a stream it does not
know, and a host that sends `-outfmt` to an earlier adapter gets its refusal
(`-outfmt is not accepted`), not other outputs.

The named input files of §4 (line NP WW1, 2026-11-12) kept version 2: they are four exports more,
with no change to a signature, a stream or a JSON field of the others. The application ships with
its module and checks that the exports exist when it instantiates the module; a host that does
not know them never registers a file, and a run of `-window_masker_db` then fails as the CLI's
does for a file it cannot open.

`-window_masker_taxid` (line NP WW2, 2026-11-12) kept version 2: it is a name of the kind 1 files of
WW1 and adds no export, stream or field. A host that does not know it never registers a file of that
name, and a run of a taxid then reports the CLI's taxid without a file.

The inflation exports, the *inflate* response and stream 5 (swebok-improve SW8 D2, 2026-11-11) kept
version 2 as stream 4 did: they add calls and a stream that a host which does not know them never
calls and never receives, and change no other export or stream. The application requires them in the
modules it loads (`ABI_EXPORTS`, `web/packages/losat/src/infra/reactor/artifact.ts`).
