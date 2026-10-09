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
| `wasi_snapshot_preview1` | `fd_write`, `fd_close`, `fd_fdstat_get`, `clock_time_get`, `random_get`, `sched_yield`, `environ_get`, `environ_sizes_get`, `proc_exit`, and the filesystem functions `path_open`, `path_create_directory`, `path_filestat_get`, `fd_prestat_get`, `fd_prestat_dir_name` | both | stderr, clocks, random, environment (the host gives an empty environment: variables such as `BATCH_SIZE` change the CLI's results, whose parity is defined without them). The filesystem functions are linked by engine code that the exports do not use (for example the CLI's `-out` file, which `validate` rejects); the host gives no preopened directory, so such a call fails instead of touching files |
| `wasi` | `thread-spawn` | threads | starts a rayon worker thread |
| `env` | `memory` | threads | shared linear memory, created by the host with the maximum the module declares (16384 pages, 1 GiB, the same as the certified threaded builds; plan TD-7) |
| `losat_host` | `emit(stream: u32, ptr: u32, len: u32)` | both | receives output bytes (§5). Worker instances of the threaded module import it too, but never call it |

## 4. Exports

All functions return `i32`: `0` (or a non-negative handle) on success, `-1` on failure.
On failure, `losat_web2_last_error_ptr/len` hold a UTF-8 message. Engine errors use the
CLI's wording: argv errors are the message of the engine's parser
(`LOSAT::cli::render_message`), search errors are the error and its causes joined by
`: ` (as ABI v1 reports them). The adapter parses the argv exactly as the command line
(every program accepts the CLI's default `-outfmt 0`), so an argv error of `validate`, an
unknown program's included, is the CLI's message (session S08). Only the adapter's own
rules give other messages: `-out` and `-outfmt` are not accepted (§7), and `blastx` is not
available until SX. `run` finds the handles of the argv's program before it parses the
argv, so it reports an unknown program as `unknown program '<name>'` and `blastx` as not
available; a host validates the argv first.

| Export | Arguments | Result |
|---|---|---|
| `losat_web2_abi_version()` | — | `2` |
| `losat_web2_alloc(len)` / `losat_web2_dealloc(ptr, len)` | — | memory for inputs (§6) |
| `losat_web2_describe(program_ptr, program_len)` | program name | emits a *describe* JSON on stream 2 |
| `losat_web2_validate(argv_ptr, argv_len)` | argv (§7) | `0` if the argv is valid, else `-1` with the CLI error. For BLASTN it also checks the scoring options as NCBI does before a search and against NCBI's Karlin-Altschul tables, and returns NCBI's message (`BLAST query/options error: …`, or `BLAST engine error: Error: …` as for one query), or LOSAT's message with `not supported by LOSAT's BLASTN` for options that NCBI runs and LOSAT does not (`docs/evidence/losat_web_e2c/AUTHORITY.md` §G). For TBLASTX it checks the options alone, with the messages of a run: NCBI's `-evalue` check (`BLAST query/options error: expect value or cutoff score must be greater than zero`) and the options that LOSAT's TBLASTX rejects (`-window_size 0`). A run whose query has no record ends with NCBI's `Query is Empty!` before LOSAT's limits, so `validate` is stricter there. For BLASTN, BLASTP, TBLASTN and TBLASTX it reads `-subject_loc` and `-query_loc` as NCBI's `ParseSequenceRange` does, in the command line's order (the subject range before the program's option checks, the query range among them), with NCBI's messages (`BLAST engine error: Invalid specification of query location (…)`, `subject location` for the subject) or LOSAT's for a part that NCBI cannot convert to an int (session S11). Whether a range fits the records is checked by `run`: a start past a subject record's end (NCBI's `Invalid from coordinate (greater than sequence length)`), a start just past a record's end (LOSAT's rejection), every query record of a batch skipped (NCBI's `Empty CBlastQueryVector`). The CLI reads the subjects before it checks the query options, so for an argv and inputs with both kinds of error the CLI reports the record's error and `validate` the argv's (`docs/evidence/losat_web_e2d/AUTHORITY.md` §A, §B). The checks of the inputs come with `register` and `run` |
| `losat_web2_register(program_ptr, program_len, role, bytes_ptr, bytes_len)` | role `0` query, `1` subject; original FASTA bytes | handle ≥ 1; emits a *register* JSON on stream 2. A handle belongs to the program that registered it. The records are those of the program's own reader, the engine's port of NCBI BLAST+'s `CBlastFastaInputSource` (session SF): nucleotide flags for BLASTN, TBLASTX and the TBLASTN subject, protein flags for BLASTP and the TBLASTN query, the bytes read as a file, the data loaders on (as the CLI without `.ncbirc` or `DATA_LOADERS`). Every record that the reader returns is registered, records without residues and a first record without a defline included; the reader's warnings stay with the records and come on stream 3 of a run, where the CLI writes them. `register` first scans the bytes with the program's scan kind (1 or 2, §9) and fails with the scan's error: LOSAT Web's rejections (`not supported by LOSAT Web`: a `>?` gap line, which the CLI reads; a first line that NCBI's data loaders would fetch as a Seq-id; a record whose residues two joined lines hide from the index; a record over 2147483647 letters) and the reader's own error with the CLI's text (`BLAST query error: CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`). A query with blank and comment lines only (`!`, `#`, `;`) fails with NCBI's `BLAST engine error: Empty CBlastQueryVector`, which every run of it ends with. A query of white space only has no record: a run gives NCBI's `Query is Empty!` warning. A subject without records is registered: a run gives NCBI's `BLAST engine error: Empty CBlastQueryVector`. `register` then checks that the scan found the reader's records (§9) |
| `losat_web2_release(handle)` | handle | `0` |
| `losat_web2_scan_begin(parser)` / `losat_web2_scan_chunk(scanner, ptr, len)` / `losat_web2_scan_end(scanner)` | parser kind: `1` NCBI BLAST+'s reader with nucleotide flags (BLASTN, TBLASTX, the TBLASTN subject), `2` with protein flags (BLASTP, the TBLASTN query), `0` the `bio::io::fasta` 1.6.0 reader that the four programs used before session SF (kept until the application switches); FASTA bytes in chunks of any size | `scan_begin` returns a scanner handle (`unknown or unavailable FASTA parser kind N` for another kind); `scan_end` emits a *scan* JSON on stream 2 (§9), or fails with the reader's error or LOSAT Web's rejection |
| `losat_web2_run(argv_ptr, argv_len, query_handle, subject_handle)` | argv, handles registered for the argv's program | emits the program's supported format streams (0, 6, 7), stream 1 (§8) and stream 3; returns after the run ends. The registered records enter the engine's `run_local` as the reader's records, so a run reads them as the CLI reads the same file: the subjects' warnings first, then the query batches with their warnings; a first batch that fails (for example a batch of queries without residues only, NCBI's `BLAST engine error: Warning: Sequence contains no data ...`) fails the run after the CLI would have written the outfmt 0 prolog to its file, and the host discards what it received (§5). Because a run writes outfmt 0, it fails as the CLI's outfmt 0 does for a subject title that LOSAT does not write as NCBI: for TBLASTX and TBLASTN, a subject in a query's final hit list (which the outfmt 0 description table lists) whose title NCBI's `NStr::HtmlDecode` changes or its `x_CleanAndCompress` reads past (`not supported by LOSAT's TBLASTX` or `TBLASTN`); for BLASTN, any subject whose title NCBI decodes (`docs/evidence/losat_web_e2b/AUTHORITY.md` §F). Outfmt 6 and 7 alone would accept such subjects |
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
| `1` | HSP records as JSON Lines, one object per HSP (§8); BLASTP, TBLASTN, BLASTN and TBLASTX |
| `2` | JSON response of `describe`, `register` or `scan_end` |
| `3` | diagnostics that the CLI writes to stderr (warnings), UTF-8, each written once |

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
  module, it passes `-num_threads 1` and records the reason (plan §4.7). The host
  validates the argv that it will run, with its `-num_threads`, so a `-num_threads` in
  the user's words is reported by `validate` (the parser rejects the repeated option).
- Every other word is parsed by the engine's CLI parser, so options and errors behave
  exactly as on the command line. LOSAT's own debug output (environment variables such
  as `LOSAT_TIMING`) goes to the module's WASI standard error, not to stream 3; stream 3
  carries the warnings that the CLI writes.
- Options are resolved separately for each output format. If the resolved *search*
  options differ between formats (for example the hitlist size), `run` fails with an
  explicit error instead of searching more than once (plan TD-4). No current program
  can reach this case.

## 8. HSP record

One JSON object per HSP of the final, sorted result. Field names follow the Rust
`common::Hit` and `report::PairwiseHit` fields. BLASTP, TBLASTN, BLASTN and TBLASTX
emit them (TBLASTX since session S08). BLASTN rows are nucleotide rows: the query on its
plus strand, residues of masked regions in lowercase, no frames. TBLASTX rows are the
translated rows of outfmt 0 (both sequences translated again from the nucleotides with
`-query_gencode` and `-db_gencode`, SEG-masked query residues in lowercase), with both
frames.

| Field | Type | Meaning |
|---|---|---|
| `index` | integer | 0-based position of the HSP in the final, sorted hit list of the run (`HspIndex` in `LOSAT/src/api/local_blast.rs`); the same HSP has the same index in every output format |
| `q_idx`, `s_idx` | integer | 0-based record index of the query and the subject in the registered inputs (the `index` of the *register* response: records without residues count; with `-query_loc`, a query record that NCBI skips because its interval starts past its end keeps its place, so `q_idx` is not the position among the searched queries) |
| `rank` | integer | 0-based position of the HSP among the HSPs of its query (derived from `index`) |
| `raw_score`, `bit_score`, `e_value` | number | engine values, not rounded |
| `q_start`, `q_end`, `s_start`, `s_end` | integer | 1-based coordinates as printed in outfmt 6 (start > end means minus strand; a BLASTN HSP of one letter has start = end on either strand, and only its `Strand=` line in the stream 0 section shows the strand, S07+ AUTHORITY §P) |
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
- *register*: `{ "handle", "records": [{ "index", "id", "length" }] }`, from the records of the program's own reader (§4). `id` is the record's title up to its first space (0x20), empty when the record has no title (never NCBI's `Query_N` or `Subject_N`); bytes that are not UTF-8 become U+FFFD in the JSON only (the output streams keep the title's bytes). `length` is the number of residues the reader stores. `register` also scans the same bytes with the program's kind (1 or 2) and fails if the records, IDs, lengths or residue counts differ.
- *scan*: `{ "records": [{ "index", "id", "header_offset", "sequence_offset", "end_offset", "length", "line_layout", "residue_counts" }] }`. Offsets are byte offsets in the scanned input and may exceed 2³²; they are encoded as JSON numbers and stay below 2⁵³. `line_layout` is `{ "kind": "uniform", "width", "eol" }` when residue `i` is at `sequence_offset + floor(i / width) * (width + eol) + i mod width`, or `{ "kind": "checkpoints", "every": 65536, "offsets": [...] }`, where `offsets[k]` is the byte offset of residue `k * 65536` and the residues are read forward from there under the kind's rules (below). `residue_counts` maps each byte value of the stored residues to its count: the bytes 0x21-0x7E as the character itself, every other byte (including the space) as `"0xNN"`.
  - Kinds 1 and 2 (session SF) give exactly the records that NCBI BLAST+'s reader returns from the bytes read as a file, with the program's flags and the data loaders on: records without residues and a first record without a defline included. The scan writes no messages. An input of white space, blank and comment lines only has no record and scans to `{"records":[]}`.
    - `id`: the title up to its first space (0x20), as in *register*. The title is the defline after `>` (after `>?_`, which NCBI reads as `>`) and the white space that follows `>`, up to the first byte below 0x20 after its first byte, without trailing white space; empty without a title; bytes that are not UTF-8 become U+FFFD.
    - `header_offset`: the defline's `>`. A first record without a defline has `header_offset` = `sequence_offset` = 0, so the extracted bytes of the record start the input as before.
    - `sequence_offset`: the first byte that NCBI's line reader reads after the defline's end of line (CR, LF or CR LF; a CR or LF that the line reader drops right after it is skipped), or the end of the input. `end_offset`: the next record's `>`, or the end of the input (a tail that NCBI's line reader loses is inside the last record's range).
    - `length` and `residue_counts`: the stored residues, upper-cased. Kind 1 stores `A B C D G H K M N R S T U V W Y` in either case and counts `U` as `T` (keys `A B C D G H K M N R S T V W Y`); kind 2 stores every ASCII letter and `*` (keys `A`-`Z`, `*`). Nothing else is counted: no `-`, digits, `*` in kind 1, invalid letters, white space or other bytes.
    - Uniform layout: `width` is the number of residues before the first jump, `eol` the number of bytes between the first line's last residue and the next residue (any positive count, not only 1 or 2); a record on one line has `width` = `length` and `eol` 1, a record without residues `width` 0 and `eol` 1. The layout is uniform only when the formula holds for every residue.
    - Reading forward from a checkpoint (its byte is a residue): read the bytes in order; CR and LF end a line (CR LF is two ends, harmless); at the start of a line skip space, tab, VT (0x0B) and FF (0x0C), and skip the whole line when its first other byte is `!`, `#` or `;`; elsewhere `;` skips the rest of the line; a byte that the kind stores (kind 1: the 16 letters above in either case; kind 2: `A`-`Z`, `a`-`z`, `*`) is the next residue (its value is the byte upper-cased, `U` read as `T` in kind 1); every other byte is skipped (white space, `-`, `>`, `?`, `_`, digits, other ASCII, bytes 0x80-0xFF, a byte order mark). A record ends at `end_offset`; read at most `length` residues.
    - `scan_end` fails, in the reader's order, with the reader's error `CFastaReader: Near line N, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` (NCBI's text without the CLI's `BLAST query error: `; N is NCBI's line number: every line counts, CR LF is one line end), or with LOSAT Web's explicit rejection (`not supported by LOSAT Web`) of a first line that NCBI's data loaders would fetch as a Seq-id (start the input with a `>` defline), of a `>?` gap line (its residues have no bytes to locate; maintainer decision 3 of session SF: the CLI reads gap lines, LOSAT Web does not), of a record over 2147483647 letters, and of a record that is not uniform and whose residues NCBI's line reader joins from two lines of the file (a CR line end in a file of LF line ends, or an LF in a file of CR line ends), which the read-forward rules cannot follow (use one kind of line end). The rejection names the line (N 1-based) or the record (`record N ("id")`).
  - Kind 0 (`bio::io::fasta` 1.6.0, the reader before session SF, kept for the application until it switches to kinds 1 and 2): `header_offset` is the `>` of the header line, `sequence_offset` the first byte after it, `end_offset` the next header line or the end of the input; `length` is the sequence length in bytes as the parser reports it; residues are read forward as every byte of a line except its trailing whitespace, and lines starting with `>` end the record. An input of white space only has no record: the kind-0 scan fails with `Expected > at record start.`, although `register` accepts it without a record (§4). LOSAT Web refuses such an input with the scan's message before it is queued (S09). `web/adapter/tests/scan_properties.rs` checks kind 0 against the parser.
  - The scan exists only to index original records for extraction; the search always reads the input with the program's own reader (plan TD-8). `web/adapter/tests/scan_ncbi_properties.rs` checks kinds 1 and 2 against the engine's reader and an independent reference walk, in chunks of any size.

## 10. Versioning

`losat_web2_abi_version` returns `2`. Any incompatible change to a signature, a stream,
or a JSON field increments the version, and the application refuses a module whose
version it does not know.
