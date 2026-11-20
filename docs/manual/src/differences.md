# Differences from NCBI BLAST+

LOSAT v0.2.0 reproduces the byte-for-byte output of NCBI BLAST+ 2.17.0 across the release test set. This page documents every intentional difference and improvement, plus rare cases where NCBI exhibits undefined behaviour, overflow, or platform drift.

Cases that LOSAT refuses with an error are listed in [What LOSAT refuses](refusals.md).

## Search results

### Summary of practical differences

| Area | LOSAT behaviour | NCBI BLAST+ behaviour | Record |
| :--- | :--- | :--- | :--- |
| **`tblastx` / `tblastn -subject -db_gencode`** | Translates `-subject` FASTA sequences using the `-db_gencode` table you requested (also adds code `32` for `tblastn`). | Silently translates local `-subject` with code `1` regardless of `-db_gencode` (often producing very different alignments/E-values for non-standard codes). | [Local `-db_gencode`](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TBLASTX-LOCAL-SUBJECT-GENCODE.md) |
| **`-db` thread determinism** | Output is identical to NCBI `-num_threads 1` at every thread count (`-num_threads` and `-mt_mode` affect speed only). | Multi-threaded `-db` runs can alter composition-based redos (e.g. `blastp`) or batch XML headers (`-mt_mode 1`). | [Database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md) |
| **`tblastx` word cutoff overflow** | Retains the finite word cutoff on large genome pairs (e.g., 1.2 Mb x 5.5 Mb finishes in ~2 min). | Overflows a double-to-int conversion to cutoff `1`, consuming tens/hundreds of GiB of RAM or crashing. | [NCBI defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| **`-taxids 1` & high-copy accessions** | Root `-taxids 1` and accessions shared by 500+ sequences run cleanly to completion. | `-taxids 1` hangs endlessly; 500+ matching OIDs fail with `Accessions to Oids lookup error`. | [NCBI defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) |
| **SAM (`-outfmt 17`) `@PG CL:` tag** | Records the `losat` command line actually invoked. | Records NCBI's binary invocation. | [CLI differences](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md) |

### A local `-subject` search of `tblastx` and `tblastn` uses `-db_gencode`

- **What you see.** With `-subject`, `tblastx` and `tblastn` translate the subject sequences with `-db_gencode`. NCBI BLAST+ translates local `-subject` sequences with genetic code `1` regardless of `-db_gencode`, passing your `-db_gencode` only to the display metadata. With any code other than `1`, LOSAT's alignments, scores, and E-values reflect the requested genetic code rather than code `1`. (`-db` searches, code `1`, and `-query_gencode` match NCBI identically.)
- **Genetic code 32.** `tblastn` also supports genetic code `32` for local `-subject` searches, which NCBI's `tblastn` does not expose.
- **Records.** [PD-TBLASTX-LOCAL-SUBJECT-GENCODE](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TBLASTX-LOCAL-SUBJECT-GENCODE.md) and [PD-TLOSAN-LOCAL-GENCODE-32](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TLOSAN-LOCAL-GENCODE-32.md).

### `-db` output does not depend on `-num_threads`

- **What you see.** LOSAT's `-db` output matches NCBI `-num_threads 1` identically regardless of `-num_threads` (`-mt_mode` `0`, `1`, and `2` are accepted without changing output). In NCBI BLAST+, multi-threaded `blastp`, `psiblast`, `rpsblast`, or `rpstblastn` runs can change composition-based cutoff sets or batch boundaries compared with single-thread runs.
- **Record.** [PD-LOSAT-DATABASE-SEARCH](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).

<details>
<summary>Rare NCBI BLAST+ edge cases (overflows, out-of-bounds reads, and damaged inputs)</summary>

All cases below stem from memory-safety errors, integer overflows, or uninitialized reads in NCBI BLAST+ 2.17.0; see [PD-LOSAT-NCBI-DEFECTS](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

- **`tblastx` word cutoff on large genome pairs:** Where NCBI overflows a double-to-int cutoff to `1` and exhausts memory (crashing at 44 GiB or hanging at 300 GiB), LOSAT keeps the exact finite cutoff (under extreme synthetic `-searchsp` + `-evalue` on tiny inputs, NCBI reports extra weak hits).
- **One-hit gapped starts past the end of a sequence (`-window_size 0`):** `blastp`, `blastx`, and `tblastn` stop at the end marker instead of reading past the sequence buffer; output matches NCBI run under Valgrind.
- **`blastn -db_hard_mask` when a sequence's first mask runs to its end (`word_size >= 11`):** LOSAT scans the unmasked prefix cleanly without scanning past-end masked letters into subsequent sequences (for `-task blastn-short`, LOSAT refuses such databases cleanly in [What LOSAT refuses](refusals.md)).
- **Titles made only of punctuation:** In `-outfmt 0` (`blastn`, `tblastn`, `tblastx`) or `lcl|Subject_N , ,` IDs, LOSAT cleans titles without reading past the string end where NCBI crashes.
- **Repeated query chunking (`CHUNK_SIZE` / `OVERLAP_CHUNK_SIZE`):** In `blastn` where a chunk would trigger recursive splitting in NCBI (raising a null-pointer crash), LOSAT searches the chunk once.
- **`-culling_limit` near 2^31 - 1 (`2147483645`–`2147483647`):** LOSAT applies no cull just like `2147483644`, avoiding NCBI's signed 32-bit overflow.
- **`gb`, `emb`, and `dbj` ID collisions across `-db` or `-parse_deflines`:** LOSAT preserves each stored prefix where NCBI's shared ID table can overwrite prefixes (matches NCBI with `OBJECTS_PACK_TEXTID=0`).
- **Extreme ASN.1/XML/PSSM integers and deep PSSM nesting:** LOSAT parses values accurately without 32-bit wrapping or 8 MiB stack crashes.
- **`rpsbproc` coordinates and PSSM IDs past 2^31:** Avoids NCBI signed 32-bit printing overflow (`-2147483648`) and out-of-bounds motif reads.
- **Damaged `taxdb.btd` records extending past NCBI's file map:** Writes NCBI's output when safe and refuses outputs that would expose unmapped bytes (see [What LOSAT refuses](refusals.md)).

</details>

## Platforms

### Last digit of a full-precision value on macOS

- **What you see.** The last (15th) digit of a value printed with all its digits can differ by 1 from the Linux and Windows value.
- **When.** macOS: the E value in `-outfmt 11` and `-outfmt 8` (tabular output is not affected). Values rounded to 6 digits (XML, JSON) or printed with 2 or 3 digits (tabular) are not affected in practice.
- **Why.** The platform's maths library differs by 1 ULP. NCBI BLAST+ on macOS is expected to show the same (not tested).
- **Record.** [The decision record on the last digit of full-precision values](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PLATFORM-LIBM.md).

### NCBI's own builds differ across platforms

- **What you see.** LOSAT produces identical Linux x86_64 bytes on Windows, macOS, and Linux aarch64. NCBI's own non-x86_64/non-Linux binaries can diverge slightly from NCBI Linux x86_64 within documented ULP/tie-break tolerances.
- **Record.** [PD-NCBI-PLATFORM-VARIANCE](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-NCBI-PLATFORM-VARIANCE.md).

## Command line and files

Full provenance: [PD-LOSAT-CLI-NONSEARCH-DIFFERENCES](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

### Threads with `-subject`

With `-subject`, NCBI BLAST+ ignores `-num_threads` and prints `'num_threads' is currently ignored when 'subject' is specified.`. LOSAT parallelizes `-subject` searches across `-num_threads` cleanly without that warning (and does not emit thread-reduction chatter when `-num_threads` exceeds core count), keeping outputs byte-identical to 1 thread.

### Compressed `-query`, `-subject` and `-in`

Native LOSAT automatically detects and reads `.gz`, `.bz2`, and `.zst` inputs by file extension for `-query`, `-subject`, and `makeblastdb -in` (matching NCBI `makeblastdb` for search commands too). WASI modules have no temporary directory and require uncompressed files or `-query -`. Invalid compressed payloads fail cleanly with exit status `1`.

### Argument errors and help

Syntax errors (unknown flags, missing values, out-of-range numeric arguments) print LOSAT's parser error to `stderr` and exit with status `2` (NCBI prints `USAGE` and exits `1`). Post-parse engine/input/validation messages and exit codes match NCBI identically.

### Version text

`-version` and `-version-full*` match NCBI (`blastn: 2.17.0+`) while naming LOSAT on the `Package:` line (`Package: LOSAT 0.2.0, output compatible with blast 2.17.0`).

### Write failures, closed pipes, and out-of-memory

When output writes succeed, results match NCBI byte for byte. If disk writes fail (`/dev/full`, permission errors) or a downstream pipe closes (`... | head`), LOSAT reports clean I/O errors or exits cleanly rather than aborting (`makeblastdb` / `rpsbproc` silent drops vs clean exit `4`–`6`). Out-of-memory conditions abort cleanly via Rust's allocator rather than C++ exception catch (`exit 4`).

## Databases and tools

### `makeblastdb`

- Reads multi-stream `.bz2` and multi-frame `.zst` inputs completely to EOF (NCBI stops after the first stream silently).
- Properly treats `-out` paths containing spaces or quotes (or overwritten DB paths) as a single database.
- Linux x86_64 `.pin`/`.nin`/`.psq`/`.nsq` files match NCBI `makeblastdb` sequence/header bytes; LMDB/metadata timestamps reflect creation time.

Records: [PD-LOSAT-NCBI-DEFECTS](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) and [PD-LOSAT-DATABASE-SEARCH](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).

### `blastdbcmd`

Fixes six NCBI extraction/masking defects (`-long_seqids` minus-strand `-range`, multi-volume v4 `gnl|BL_ORD_ID|N`, `-mask_sequence_with` on minus strand, repeated `%m` formatters, protein `-mask_sequence_with` without range, and mixed masked/unmasked batch rows); see [PD-LOSAT-BLASTDBCMD-DEVIATIONS](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

## The library

### `-outfmt 11` and preserved definition lines

The CLI keeps raw FASTA definition lines verbatim for `-outfmt 11`. Calling `losat::api` with pre-parsed records that do not preserve the raw definition line reconstructs `>` + title cleanly; see [PD-LOSAT-ARCHIVE-DEFLINE](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-ARCHIVE-DEFLINE.md).
