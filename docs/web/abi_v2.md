# LOSAT Web ABI v2 (engine adapter ↔ application)

Status: **version 2, implemented** in `web/adapter/` (session S05). Drafted in S01
([`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md) §4, §7). The application
depends on the TypeScript port `web/app/src/ports/engine.ts`, not on the exported
functions. The HSP record (§8) and the JSON responses (§9) cross into the
application, so any change to them updates `web/app/src/ports/engine.ts` in the same
commit, and an incompatible change increments the version (§10).

Programs: BLASTP, TBLASTN, BLASTN and TBLASTX. BLASTX joins in session SX; until
then `describe`, `validate`, `register` and `run` reject it.

## 1. Purpose

The application (`web/app`) runs LOSAT searches through a Wasm module built from the
adapter crate `web/adapter`, which links the engine crate `LOSAT/` as a library.

- Search conditions are one CLI argv, parsed by the same clap parser as the CLI.
- One run writes every output format the program supports (listed by `describe`),
  the HSP records and the diagnostics from the same search result (§7).
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
outside the adapter's workspace, so rustc gets its sources by absolute path: the panic
locations of the engine code (data section) contain the checkout path, which changes the
code and data sections as well as the symbol hashes of the name section. Every crate
from the registry is compiled from `CARGO_HOME` (`<CARGO_HOME>/registry/src/...` in its
panic locations), as in the engine's own Wasm builds. Two builds of the same commit with
the same checkout path and `CARGO_HOME` are identical; the identity record
(`losat-web-*.json`) records both. These paths are local file names of the machine that
built the module; they are removed from published modules before the first publication
(plan TD-11).

## 3. Imports

| Module | Name | Used by | Purpose |
|---|---|---|---|
| `wasi_snapshot_preview1` | `fd_write`, `fd_close`, `fd_fdstat_get`, `clock_time_get`, `random_get`, `sched_yield`, `environ_get`, `environ_sizes_get`, `proc_exit`, and the filesystem functions `path_open`, `path_create_directory`, `path_filestat_get`, `fd_prestat_get`, `fd_prestat_dir_name` | both | stderr, clocks, random, environment. The filesystem functions are linked by engine code that the exports do not use (for example the CLI's `-out` file, which `validate` rejects); the host gives no preopened directory, so such a call fails instead of touching files |
| `wasi` | `thread-spawn` | threads | starts a rayon worker thread |
| `env` | `memory` | threads | shared linear memory, created by the host with the maximum the module declares (16384 pages, 1 GiB, the same as the certified threaded builds; plan TD-7) |
| `losat_host` | `emit(stream: u32, ptr: u32, len: u32)` | both | receives output bytes (§5). Worker instances of the threaded module import it too, but never call it |

## 4. Exports

All functions return `i32`: `0` (or a non-negative handle) on success, `-1` on failure.
On failure, `losat_web2_last_error_ptr/len` hold a UTF-8 message. Engine errors use the
CLI's wording: argv errors are the message of the engine's parser
(`LOSAT::cli::render_message`), search errors are the error and its causes joined by
`: ` (as ABI v1 reports them). Two details of argv errors differ from the command line:
the adapter parses the argv with `losat` as the binary name and `-outfmt 6` inserted
after the program name (§7), so a usage line in a message names `losat` and can list
`-outfmt`; and an unknown program gives `unknown program '<name>'`.

| Export | Arguments | Result |
|---|---|---|
| `losat_web2_abi_version()` | — | `2` |
| `losat_web2_alloc(len)` / `losat_web2_dealloc(ptr, len)` | — | memory for inputs (§6) |
| `losat_web2_describe(program_ptr, program_len)` | program name | emits a *describe* JSON on stream 2 |
| `losat_web2_validate(argv_ptr, argv_len)` | argv (§7) | `0` if the argv is valid, else `-1` with the CLI error |
| `losat_web2_register(program_ptr, program_len, role, bytes_ptr, bytes_len)` | role `0` query, `1` subject; original FASTA bytes | handle ≥ 1; emits a *register* JSON on stream 2 |
| `losat_web2_release(handle)` | handle | `0` |
| `losat_web2_scan_begin(parser)` / `losat_web2_scan_chunk(scanner, ptr, len)` / `losat_web2_scan_end(scanner)` | parser kind (`0` the `bio::io::fasta` reader of BLASTP, TBLASTN, BLASTN and TBLASTX; `1`, the NCBI-style reader of BLASTX, joins in SX); FASTA bytes in chunks of any size | `scan_begin` returns a scanner handle; `scan_end` emits a *scan* JSON on stream 2 (§9), or fails with the parser's error |
| `losat_web2_run(argv_ptr, argv_len, query_handle, subject_handle)` | argv, handles | emits the program's supported format streams (0, 6, 7), stream 1 (BLASTP and TBLASTN; §8) and stream 3; returns after the run ends |
| `losat_web2_last_error_ptr()` / `losat_web2_last_error_len()` | — | last error message |

## 5. Output streams

`emit` may be called many times per stream. Chunks of one stream arrive in order. The
bytes are valid only during the call, so the host must copy them. `emit` is always
called on the thread that called the export: the engine formats on the calling thread,
which runs slot zero of the search's thread pool (`LOSAT/src/utils/threading.rs`).

| Stream | Content |
|---|---|
| `0` | outfmt 0 text (only for programs that support outfmt 0) |
| `6` | outfmt 6 text |
| `7` | outfmt 7 text (only for programs that support outfmt 7) |
| `1` | HSP records as JSON Lines, one object per HSP (§8); BLASTP and TBLASTN only |
| `2` | JSON response of `describe`, `register` or `scan_end` |
| `3` | diagnostics that the CLI writes to stderr (warnings), UTF-8, each written once |

The adapter sends each output stream in chunks of 1 MiB as they fill, and the rest when
the run ends, so an output never has to fit in linear memory as a whole. Stream 3 and
the HSP records are sent after the search. If `run` fails, the host discards what it
received.

## 6. Memory ownership

The host allocates input buffers with `losat_web2_alloc`, writes the bytes, passes the
pointer, and frees the buffer with `losat_web2_dealloc` after the call returns. The
adapter never keeps a pointer to a host buffer after a call: `register` copies what it
keeps.

## 7. argv

- UTF-8 words separated by NUL (`\0`). The first word is the program name
  (`blastn`, `blastp`, `blastx`, `tblastn`, `tblastx`).
- `-query <name>` and `-subject <name>` are required. The names are used only in
  output text (for example the `# Database:` line of outfmt 7); the sequences come from
  the handles.
- `-out` and `-outfmt` are rejected: the adapter writes every format the program supports
  (`describe` lists them).
- `-num_threads` is appended by the host. When the host falls back to the serial
  module, it passes `-num_threads 1` and records the reason (plan §4.7).
- Every other word is parsed by the engine's CLI parser, so options and errors behave
  exactly as on the command line. LOSAT's own progress output (`-verbose`) and its
  debug output (environment variables such as `LOSAT_TIMING`) go to the module's WASI
  standard error, not to stream 3; stream 3 carries the warnings that the CLI writes.
- Options are resolved separately for each output format. If the resolved *search*
  options differ between formats (for example the hitlist size), `run` fails with an
  explicit error instead of searching more than once (plan TD-4). No current program
  can reach this case.

## 8. HSP record

One JSON object per HSP of the final, sorted result. Field names follow the Rust
`common::Hit` and `report::PairwiseHit` fields. BLASTP and TBLASTN emit them. BLASTN and
TBLASTX emit none until their final HSP lists become `PairwiseHit` lists (plan S07 and
S08); their `run` writes the format streams and stream 3 only.

| Field | Type | Meaning |
|---|---|---|
| `index` | integer | 0-based position of the HSP in the final, sorted hit list of the run (`HspIndex` in `LOSAT/src/api/local_blast.rs`); the same HSP has the same index in every output format |
| `q_idx`, `s_idx` | integer | 0-based record index of the query and the subject in the registered inputs |
| `rank` | integer | 0-based position of the HSP among the HSPs of its query (derived from `index`) |
| `raw_score`, `bit_score`, `e_value` | number | engine values, not rounded |
| `q_start`, `q_end`, `s_start`, `s_end` | integer | 1-based coordinates as printed in outfmt 6 (start > end means minus strand) |
| `query_frame`, `subject_frame` | integer or null | translation frames where applicable |
| `subject_length` | integer or null | subject length |
| `query_aligned`, `subject_aligned` | string or null | aligned sequences with `-` for gaps |
| `out6` | [integer, integer] or null | byte range [start, end) of this HSP's row in the stream 6 text; null if outfmt 6 does not show it |
| `out0` | [integer, integer] or null | byte range [start, end) of this HSP's section in the stream 0 text: its score lines and alignment. The subject heading that precedes the first HSP of each subject is not part of any section; null if outfmt 0 does not show the HSP (for example BLASTX shows alignments for the first 250 subjects by default) |
| `out0_subject` | [integer, integer] or null | byte range [start, end) of the heading of this HSP's subject in the stream 0 text: the defline and the `Length=` line that NCBI's `x_ShowAlnvecInfo` writes before the first HSP of each subject (`showalign.cpp:3613-3632`). Every HSP of the subject in that query has the same range |

The ranges come from formatter observer events: a formatter reports when it starts and
ends the row or section of an HSP, identified by its `index`, and the heading of a
subject, identified by the index of its first HSP; the adapter records the byte
positions at those moments. The formatter's output bytes do not
change (plan TD-3, §4.5). Displayed numbers are taken from the outfmt 6 row, not
formatted from the raw values (plan §4.4).

## 9. JSON responses (stream 2)

- *describe*: `{ "program", "formats": [0, 6, 7], "parameters": [{ "flag", "help", "takes_value", "default"?, "choices"? }], "query_gencodes"?, "subject_gencodes"? }`, generated from the clap definitions (without `-out`, `-outfmt` and `-num_threads`, which the adapter and the host set), the program's supported output formats and the genetic codes that the engine's parser accepts for `-query_gencode` and `-db_gencode`.
- *register*: `{ "handle", "records": [{ "index", "id", "length" }] }`, produced by the same FASTA parser that the program uses in a search. `register` also scans the same bytes (parser kind of the program) and fails if the record IDs or lengths differ.
- *scan*: `{ "records": [{ "index", "id", "header_offset", "sequence_offset", "end_offset", "length", "line_layout", "residue_counts" }] }`. Offsets are byte offsets in the scanned input and may exceed 2³²; they are encoded as JSON numbers and stay below 2⁵³. `header_offset` is the `>` of the header line, `sequence_offset` the first byte after it, `end_offset` the next header line or the end of the input; `length` is the sequence length in bytes as the parser reports it.
  - `line_layout` is `{ "kind": "uniform", "width", "eol" }` when residue `i` is at `sequence_offset + floor(i / width) * (width + eol) + i mod width`, or `{ "kind": "checkpoints", "every": 65536, "offsets": [...] }`, where `offsets[k]` is the byte offset of residue `k * 65536`; from there the residues are read forward under the parser's rules (for kind 0: every byte of a line except its trailing whitespace, lines starting with `>` end the record).
  - `residue_counts` maps each byte value of the sequence to its count: the bytes 0x21-0x7E as the character itself, every other byte (including the space) as `"0xNN"`.
  - The scan exists only to index original records for extraction; the search always parses the input with the program's own reader (plan TD-8). Kind 0 reproduces `bio::io::fasta` 1.6.0; `web/adapter/tests/scan_properties.rs` checks it against the parser.

## 10. Versioning

`losat_web2_abi_version` returns `2`. Any incompatible change to a signature, a stream,
or a JSON field increments the version, and the application refuses a module whose
version it does not know.
