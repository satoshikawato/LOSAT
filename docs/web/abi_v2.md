# LOSAT Web ABI v2 (engine adapter ↔ application)

Status: **draft.** Written in session S01 as the input for session S05, which
implements it in `web/adapter/` and finalizes the names and types
([`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md) §4, §7). The application
depends on the TypeScript port `web/app/src/ports/engine.ts`, not on the exported
functions, so S05 may rename or reshape the functions (§4) without changing the
application. The HSP record (§8) and the JSON responses (§9) cross into the
application, so any change to them updates `web/app/src/ports/engine.ts` in the same
commit.

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
dependency versions (plan TD-6); S05 adds a check for this.

## 3. Imports

| Module | Name | Used by | Purpose |
|---|---|---|---|
| `wasi_snapshot_preview1` | standard functions | both | stderr, clocks, random, environment. No filesystem access is needed |
| `wasi` | `thread-spawn` | threads | starts a rayon worker thread |
| `env` | `memory` | threads | shared linear memory, created by the host with the maximum the module declares (16384 pages, 1 GiB, the same as the certified threaded builds; plan TD-7) |
| `losat_host` | `emit(stream: u32, ptr: u32, len: u32)` | both | receives output bytes (§5) |

## 4. Exports

All functions return `i32`: `0` (or a non-negative handle) on success, `-1` on failure.
On failure, `losat_web2_last_error_ptr/len` hold a UTF-8 message. Engine errors use the
same wording as the CLI.

| Export | Arguments | Result |
|---|---|---|
| `losat_web2_abi_version()` | — | `2` |
| `losat_web2_alloc(len)` / `losat_web2_dealloc(ptr, len)` | — | memory for inputs (§6) |
| `losat_web2_describe(program_ptr, program_len)` | program name | emits a *describe* JSON on stream 2 |
| `losat_web2_validate(argv_ptr, argv_len)` | argv (§7) | `0` if the argv is valid, else `-1` with the CLI error |
| `losat_web2_register(program_ptr, program_len, role, bytes_ptr, bytes_len)` | role `0` query, `1` subject; original FASTA bytes | handle ≥ 1; emits a *register* JSON on stream 2 |
| `losat_web2_release(handle)` | handle | `0` |
| `losat_web2_scan_begin(parser)` / `losat_web2_scan_chunk(scanner, ptr, len)` / `losat_web2_scan_end(scanner)` | parser kind (`0` the `bio::io::fasta` reader, `1` the NCBI-style reader of BLASTX); FASTA bytes in chunks | `scan_end` emits a *scan* JSON on stream 2 (§9) |
| `losat_web2_run(argv_ptr, argv_len, query_handle, subject_handle)` | argv, handles | emits the program's supported format streams (0, 6, 7), and streams 1 and 3; returns after the run ends |
| `losat_web2_last_error_ptr()` / `losat_web2_last_error_len()` | — | last error message |

## 5. Output streams

`emit` may be called many times per stream. Chunks of one stream arrive in order. The
bytes are valid only during the call, so the host must copy them.

| Stream | Content |
|---|---|
| `0` | outfmt 0 text (only for programs that support outfmt 0) |
| `6` | outfmt 6 text |
| `7` | outfmt 7 text (only for programs that support outfmt 7) |
| `1` | HSP records as JSON Lines, one object per HSP (§8) |
| `2` | JSON response of `describe`, `register` or `scan_end` |
| `3` | diagnostics that the CLI writes to stderr (warnings), UTF-8, each written once |

The adapter flushes each stream in chunks of at most 1 MiB, so an output never has to
fit in linear memory as a whole.

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
- `-out` and `-outfmt` are rejected: the adapter writes every format the program supports.
- `-num_threads` is appended by the host. When the host falls back to the serial
  module, it passes `-num_threads 1` and records the reason (plan §4.7).
- Every other word is parsed by the engine's CLI parser, so options and errors behave
  exactly as on the command line.
- Options are resolved separately for each output format. If the resolved *search*
  options differ between formats (for example the hitlist size), `run` fails with an
  explicit error instead of searching more than once (plan TD-4). No current program
  can reach this case.

## 8. HSP record

One JSON object per HSP of the final, sorted result. Field names follow the Rust
`common::Hit` and `report::PairwiseHit` fields.

| Field | Type | Meaning |
|---|---|---|
| `q_idx`, `s_idx` | integer | 0-based record index of the query and the subject in the registered inputs |
| `rank` | integer | 0-based position of the HSP in the final result of its query; with `q_idx` it identifies the HSP within the run |
| `raw_score`, `bit_score`, `e_value` | number | engine values, not rounded |
| `q_start`, `q_end`, `s_start`, `s_end` | integer | 1-based coordinates as printed in outfmt 6 (start > end means minus strand) |
| `query_frame`, `subject_frame` | integer or null | translation frames where applicable |
| `subject_length` | integer or null | subject length |
| `query_aligned`, `subject_aligned` | string or null | aligned sequences with `-` for gaps |
| `out6` | [integer, integer] or null | byte range [start, end) of this HSP's row in the stream 6 text; null if outfmt 6 does not show it |
| `out0` | [integer, integer] or null | byte range [start, end) of this HSP's section in the stream 0 text; null if outfmt 0 does not show it (for example BLASTX shows alignments for the first 250 subjects by default) |

The ranges come from formatter observer events: a formatter reports when it starts and
ends the row or section of an HSP, identified by (`q_idx`, `rank`), and the adapter
records the byte positions at those moments. The formatter's output bytes do not
change (plan TD-3, §4.5). Displayed numbers are taken from the outfmt 6 row, not
formatted from the raw values (plan §4.4).

## 9. JSON responses (stream 2)

- *describe*: `{ "program", "formats": [0, 6, 7], "parameters": [{ "flag", "help", "takes_value", "default"?, "choices"? }], "query_gencodes"?, "subject_gencodes"? }`, generated from the clap definitions, the program's supported output formats and the engine's genetic-code allow lists.
- *register*: `{ "handle", "records": [{ "index", "id", "length" }] }`, produced by the same FASTA parser that the program uses in a search.
- *scan*: `{ "records": [{ "index", "id", "header_offset", "sequence_offset", "end_offset", "length", "line_layout", "residue_counts" }] }`. Offsets are byte offsets in the scanned input and may exceed 2³²; they are encoded as JSON numbers and stay below 2⁵³. The scan exists only to index original records for extraction; the search always parses the input with the program's own reader, and `register` checks that both agree (plan TD-8).

## 10. Versioning

`losat_web2_abi_version` returns `2`. Any incompatible change to a signature, a stream,
or a JSON field increments the version, and the application refuses a module whose
version it does not know.
