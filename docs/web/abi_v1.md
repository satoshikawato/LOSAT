# LOSAT Web ABI v1 (`losat_web_*`)

Status: **version 1, frozen** (plan TD-1), apart from fail-fast defect fixes. The
exports are in `LOSAT/src/web_api.rs`, which `LOSAT/src/lib.rs:28` builds for
`wasm32` only. The program modules are `LOSAT/src/web_api/v1_blastn.rs`,
`v1_blastp.rs`, `v1_tblastx.rs` and `v1_bio.rs`; BLASTX has its own
`LOSAT/src/algorithm/blastx/web.rs`. gbdraw calls this ABI. The adapter
(`web/adapter/`) links it in unchanged next to ABI v2 ([`abi_v2.md`](abi_v2.md)), and
the two ABIs share no state. The v1 exports are in both the serial and the threaded
module.

This file describes what the code does. It does not add a rule.

## 1. Purpose

A host gives two FASTA inputs, a program name, an output format and extra options to
the module. The module runs one search and keeps the report text, or an error message.
The host then reads the text or the message from linear memory.

- One search compares one query input with one subject input. There is no database
  search.
- The result is the report text of the search. It is not JSON.
- The module keeps one result buffer, one error buffer and one FASTA store for the whole
  instance (`LOSAT/src/web_api.rs:25-27`).

## 2. Exports

Every export is `extern "C"`. Wasm sees `usize` and pointers as `i32` (wasm32).

| Export | Parameters | Returns |
|---|---|---|
| `losat_web_alloc` | `len: i32` | pointer `i32` to `len` bytes of uninitialised memory (§4). It has no failure value; a failed allocation is a trap (§9) |
| `losat_web_dealloc` | `ptr: i32`, `len: i32` | nothing. A null `ptr` does nothing (§4) |
| `losat_web_store_fasta` | `fasta_ptr`, `fasta_len`, `label_ptr`, `label_len` (all `i32`) | handle `≥ 1` on success; `-1` on failure (§6). Success clears the error buffer |
| `losat_web_release_fasta` | `handle: i32` | `0` on success; `-1` if the handle is unknown (§6) |
| `losat_web_clear_fasta_store` | — | `0`; `-1` only if the store's mutex is poisoned |
| `losat_web_run_pair` | `program`, `query`, `subject`, `outfmt`, `extra_args`, each a `ptr, len` pair (10 × `i32`) | `0` on success, with the report in the result buffer; `-1` on failure, with the message in the error buffer |
| `losat_web_run_pair_handles` | `program_ptr`, `program_len`, `query_handle`, `subject_handle`, `outfmt_ptr`, `outfmt_len`, `extra_args_ptr`, `extra_args_len` (all `i32`) | as `losat_web_run_pair`. Only BLASTP and BLASTX (§7) |
| `losat_web_result_ptr` | — | pointer `i32` to the result buffer (§3) |
| `losat_web_result_len` | — | length `i32` of the result buffer in bytes |
| `losat_web_error_ptr` | — | pointer `i32` to the error buffer (§3) |
| `losat_web_error_len` | — | length `i32` of the error buffer in bytes |
| `losat_web_clear` | — | nothing. Empties the result and the error buffer (§3) |
| `losat_web_has_direct_api` | — | `1`. A probe that the handle exports exist; it has no other state |

Locations: `losat_web_alloc` `LOSAT/src/web_api.rs:897`, `losat_web_dealloc` `:905`,
`losat_web_store_fasta` `:915`, `losat_web_release_fasta` `:955`,
`losat_web_clear_fasta_store` `:978`, `losat_web_run_pair` `:995`,
`losat_web_run_pair_handles` `:1038`, `losat_web_result_ptr` `:1071`,
`losat_web_result_len` `:1079`, `losat_web_error_ptr` `:1084`,
`losat_web_error_len` `:1089`, `losat_web_clear` `:1094`,
`losat_web_has_direct_api` `:1100`.

There is no version export. The version is not a value that a host can read.

The only status values are `0`, a handle `≥ 1` and `-1`. `-1` carries no reason: the
reason is the error buffer (§3). A host reads the status before it reads either buffer.

## 3. Result and error buffers

`set_result` and `set_error` (`LOSAT/src/web_api.rs:132-141`) are the only writers.

| Event | Result buffer | Error buffer |
|---|---|---|
| `run_pair` / `run_pair_handles` succeeds | replaced by the report | emptied |
| any export that returns `-1` | emptied | replaced by the message |
| `store_fasta` succeeds | kept | emptied |
| `release_fasta` / `clear_fasta_store` returns `0` | kept | emptied |
| `losat_web_clear` | emptied | emptied |
| `alloc`, `dealloc`, the `*_ptr` / `*_len` exports | unchanged | unchanged |

So a failed `store_fasta` or `release_fasta` also empties the last result. A host that
needs a result later copies it out first.

- After `run_pair` or `run_pair_handles` returns `0`, `losat_web_result_ptr()` and
  `losat_web_result_len()` describe the report. After `-1`, `losat_web_error_ptr()` and
  `losat_web_error_len()` describe the message, and the result length is `0`.
- The report is the bytes the program's report layer wrote: text in the format given
  by `outfmt` (§7). The module adds no header and no terminator. An empty report has
  length `0`.
- The error message is UTF-8 text. Engine errors are the error and its causes joined by
  `: ` (`engine_error`, `LOSAT/src/web_api.rs:128-130`).
- A pointer is valid until the next call into the module that writes the same buffer
  (the table above). A host's view of the memory (`memory.buffer` in JavaScript) is
  replaced when the memory grows, and any call that allocates (`losat_web_alloc`
  included) can grow it. So the host takes a new view after such a call before it reads,
  or copies the bytes out first.
- For an empty buffer the pointer is still a valid non-null address (a Rust `Vec`
  that holds no bytes); it must not be read.
- `losat_web_clear` empties the two buffers. It sets their length to `0`. It does not
  free their capacity, and it does not touch the FASTA store.
- Diagnostics (warnings that the CLI writes to standard error) go to the module's
  standard error, which is WASI file descriptor 2. They are not in either buffer.

## 4. Memory ownership

- `losat_web_alloc(len)` makes a `Vec<u8>` of capacity `len`, forgets it, and returns
  its pointer (`LOSAT/src/web_api.rs:897-902`). The caller owns the block. The block has
  alignment 1 and its bytes are not initialised.
- `losat_web_dealloc(ptr, len)` rebuilds the `Vec` with capacity `len` and drops it
  (`LOSAT/src/web_api.rs:905-911`). `len` must be the value given to `losat_web_alloc`
  for that pointer. A pointer from anywhere else, or a second `dealloc` of the same
  pointer, is undefined behaviour. A null pointer does nothing.
- `losat_web_alloc(0)` returns a non-null dangling pointer; `losat_web_dealloc` of it
  with `0` frees nothing.
- If the allocation fails, Rust aborts, and the export ends in a trap (§9).
- The module never keeps a host pointer after a call returns. `store_fasta` copies the
  FASTA bytes (and parses them); `run_pair` copies what it needs.
- The module owns the result, error and FASTA store memory. The host never frees it.
  The host frees only what it allocated with `losat_web_alloc`, normally after the
  call that used it returns.
- Pointers and lengths are unsigned 32-bit values. A JavaScript host converts a pointer
  that it receives as a signed `i32` with `ptr >>> 0` (as for ABI v2, `abi_v2.md` §6).

## 5. Input reading (`read_str`, `read_bytes`)

All `ptr, len` inputs go through two helpers (`LOSAT/src/web_api.rs:143-166`). `name` is
the label in the error message.

| Case | `read_str` | `read_bytes` |
|---|---|---|
| `ptr` null, `len` 0 | `""` | empty slice |
| `ptr` null, `len` > 0 | error `<name> pointer is null` | error `<name> pointer is null` |
| `ptr` not null, `len` 0 | `""` | empty slice |
| bytes are not UTF-8 | error `<name> is not valid UTF-8: <detail>` | not checked |
| otherwise | the text | the bytes |

- There is no upper limit on `len` and no check that `ptr + len` lies inside the
  memory. The host must pass a block it allocated with `losat_web_alloc` and written.
  A range outside the memory ends in a trap (§9).
- `read_bytes` reads only `fasta` of `losat_web_store_fasta`. `read_str` reads the label
  of `store_fasta`, and `program`, `query`, `subject`, `outfmt` and `extra args` of the
  run exports. So `run_pair` needs UTF-8 FASTA; `store_fasta` takes any bytes and the
  FASTA parser then judges them.
- Names in messages: `program`, `query FASTA`, `subject FASTA`, `outfmt`, `extra args`
  (`run_pair`, `LOSAT/src/web_api.rs:995`), `program`, `outfmt`, `extra args`
  (`run_pair_handles`), and `FASTA`, `FASTA label` (`store_fasta`).
- The run exports read their inputs in the order of the parameters and stop at the first
  error. `run_pair_handles` does not read FASTA text; it reads handles.

## 6. FASTA handles

The store (`FastaStore`, `LOSAT/src/web_api.rs:72-123`) holds parsed FASTA inputs, so a
host can search one input many times without sending it again.

- `losat_web_store_fasta(fasta, label)` copies the bytes, parses them with the `bio`
  FASTA reader, and returns a handle. Handles are `i32` and start at `1`. The next
  handle is the previous one plus one.
    - Empty input, or input with no record, fails with
      `FASTA store input did not contain any records`.
    - A parse failure gives `failed to parse stored FASTA: <detail>`
      (`LOSAT/src/web_api.rs:168-181`).
    - When the counter would pass `i32::MAX`, the call fails with
      `FASTA store handle space exhausted`.
    - The store keeps the original bytes and the parsed records.
    - `label` is any UTF-8 text, empty allowed. BLASTP and BLASTX use the subject's
      label as the subject name in the report; the query label is used by BLASTP only.
- `losat_web_release_fasta(handle)` removes one entry and returns `0`. An unknown
  handle (never issued, already released, `0`, negative) returns `-1` with
  `unknown FASTA handle: <handle>`. Releasing twice is an error the second time. A
  release does not lower the counter, so a released number is not issued again.
- `losat_web_clear_fasta_store()` removes every entry and sets the counter back to `1`.
  The next `store_fasta` returns `1` again. A handle held from before the clear is
  therefore not rejected after a new store: it names the new entry that has its number.
  A host that clears the store drops all its old handles.
- `losat_web_run_pair_handles` looks up the query handle and then the subject handle,
  and fails with `unknown query FASTA handle: <n>` or `unknown subject FASTA handle: <n>`.
  The lookup happens after `outfmt` and `extra args` are read, and for BLASTP and BLASTX
  only. The same handle may be used as both the query and the subject.
- The store lives as long as the instance. Nothing releases an entry but the two
  release exports. Handles are not valid in another instance (§9).
- The search holds the store's lock while it runs (BLASTP and BLASTX).

## 7. Request: program, `outfmt`, extra arguments

`run_pair` (`LOSAT/src/web_api.rs:725-803`) and `run_pair_handles` (`:806-894`) build the
same request.

- `program` is one of `blastp`, `blastn`, `tblastx`, `blastx`, in lower case. Anything
  else (`tblastn` included) fails with `unsupported LOSAT program: <name>`.
- `extra_args` is a list of words separated by NUL (`\0`). Empty words are dropped
  (`parse_extra_args`, `LOSAT/src/web_api.rs:183`). A value is the next word
  (`-evalue` `1e-5`) or joined with `=` (`-evalue=1e-5`).
- `outfmt` is a string. If it is empty or white space only, it is `6`
  (`web_outfmt_or_default`, `:201`). The module then appends `-outfmt` and this value
  after the extra words (`:734-735`, `:815-816`).
- The report is for one output format only.

| Program | Module | `outfmt` accepted | Extra words accepted |
|---|---|---|---|
| `blastn` | `LOSAT/src/web_api/v1_blastn.rs`; arguments `parse_blastn_args` (`LOSAT/src/web_api.rs:235`) | `6`, `7`. `0` fails with `unsupported BLASTN output format: 0`; a value with more than one word fails with `unsupported BLASTN custom outfmt specification` | `-task` (`megablast` or `blastn`), `-outfmt`, `-word_size`, `-num_threads`, `-evalue` |
| `blastp` | `LOSAT/src/web_api/v1_blastp.rs`; arguments `parse_blastp_args` (`:527`) | the formats that BLASTP accepts, with field lists for the tabular ones (`check_web_v1_request`, `v1_blastp.rs:104`) | `-task`, `-outfmt`, `-evalue`, `-threshold`, `-word_size`, `-num_threads`, `-max_target_seqs`, `-max_hsps`, `-window_size`, `-matrix`, `-gapopen`, `-gapextend`, `-comp_based_stats`, `-seg`, `-ungapped`, `-use_sw_tback` |
| `tblastx` | `LOSAT/src/web_api/v1_tblastx.rs`; arguments `parse_tblastx_args` (`:361`) | `6` only (`tblastx_v1_outfmt`, `LOSAT/src/blastinput/value_parsers.rs:1136`) | `-query_gencode`, `-db_gencode`, `-num_threads`, `-outfmt`, `-evalue` |
| `blastx` | `LOSAT/src/algorithm/blastx/web.rs` | `0`, `6`, `7`, with field lists for `6` and `7` (`v1_outfmt`, `web.rs:68`) | the BLASTX options of the CLI's parser, except `-db`, `-query`, `-subject`, `-out`, `-help`, `--help`, `--version` (`web.rs:14-30`) |

- For `blastn`, `blastp` and `tblastx` any other word fails with
  `unsupported <program> argument for web API: <word>` (`:339`, `:463`, `:717`). A flag
  at the end of the list without a value fails with `<flag> requires a value`. A value
  that does not parse fails with `<flag> parse error: <detail>`. An `-outfmt` word in
  the extra words is replaced by the `outfmt` parameter, because the parameter is
  appended last.
- `blastn`: `-task` other than `megablast` or `blastn` fails with
  `unsupported blastn task for web API: <task>` (`:347`). The default task is
  `megablast`. `-word_size 28` is read as no word size, so `-task blastn -word_size 28`
  has the blastn default of 11 (`v1_blastn_word_size`, `:357`). The `E`-value default is
  10. The strands are both.
- `blastp`: the default task is `blastp` and `-comp_based_stats` is `2`.
- The other options of each program take the defaults of the args structure in the
  parser. Ranges (`-query_loc`, `-subject_loc`), masking, and database options are not
  available in v1.
- `-num_threads` defaults to `1`. The serial module rejects a value above 1
  (`abi_v2.md` §2); `validate_threads` is in `LOSAT/src/utils/threading.rs:104`.
- `run_pair_handles` supports only `blastp` and `blastx`. For `blastn` and `tblastx` it
  fails with `FASTA handle API is not yet ported for <program>; use losat_web_run_pair`
  (`LOSAT/src/web_api.rs:889`).
- Besides the parser, the v1 modules check the inputs with helpers in
  `LOSAT/src/web_api/v1_bio.rs` (its module header lists them): deflines that are
  empty, start with white space, hold a control or non-ASCII byte, or start with `?`;
  non-ASCII bytes in sequence lines; records without residues; residues that are not
  IUPAC letters. Which checks run depends on the program: `v1_blastn.rs`, `v1_blastp.rs`
  and `v1_tblastx.rs` each call their own set. These messages say that NCBI BLAST+
  reads the input differently and that LOSAT's program does not support it.

## 8. Order of one call

1. `losat_web_alloc` for every input, write the bytes.
2. `losat_web_run_pair` (or `store_fasta`, then `run_pair_handles`).
3. `losat_web_dealloc` for every input block.
4. Read the status. On `0`, copy `losat_web_result_len()` bytes from
   `losat_web_result_ptr()`. On `-1`, copy the message from `losat_web_error_ptr()` and
   `losat_web_error_len()`.
5. `losat_web_clear()` when the host no longer needs the buffers.

## 9. After a trap

The release profile sets `panic = "abort"` (`LOSAT/Cargo.toml:141-145`). The adapter has
a copy of it (`web/adapter/Cargo.toml:27-32`, checked by
`web/adapter/tools/check_build_identity.py`), and the adapter builds the shipped
modules. So a Rust panic inside any export does not unwind and does not return an error
code. It ends as a Wasm trap, and the host sees an exception (in JavaScript, a
`WebAssembly.RuntimeError`). The same holds for an allocation failure and for a read
outside the memory (§5).

After a trap:

- Do not call the instance again, not even `losat_web_clear` or `losat_web_dealloc`.
  Its linear memory is in an unknown state.
- The FASTA store, every handle, the result buffer and the error buffer are lost. A
  handle from the old instance means nothing in a new one. The status codes and buffers
  of the call that trapped do not exist; there is no error message to read.
- The host discards the instance and instantiates the module again (including
  `_initialize` for the reactor, `abi_v2.md` §2). It stores the FASTA inputs again and
  gets new handles. It runs the request again only if it wants to.
- In the threaded module (`wasm32-wasip1-threads`) the main instance and the worker
  instances share one memory that the host created (`abi_v2.md` §3). The host
  discards the instance, the shared memory and the workers together, and creates new
  ones. The old shared memory is not reused: its contents are in an unknown state.

A failure that the code reports (a status of `-1`) is not a trap. The instance stays
usable after it.
