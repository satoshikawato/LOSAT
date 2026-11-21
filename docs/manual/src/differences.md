# Differences from NCBI BLAST+

LOSAT v0.2.0 writes the output of NCBI BLAST+ 2.17.0. This page lists every case in which it does not, with what you see, when it happens, and why the difference is acceptable or better. Each difference has a decision record on GitHub, linked from its item.

Cases that LOSAT refuses with an error are on [What LOSAT refuses](refusals.md), not here.

## Search results

### A local `-subject` search of `tblastx` and `tblastn` uses `-db_gencode`

- **What you see.** With `-subject`, `tblastx` and `tblastn` translate the subject sequences with the genetic code of `-db_gencode`. NCBI BLAST+ translates a local subject with code 1 whatever `-db_gencode` says, and only hands your code to the report. With a code other than 1, LOSAT's alignments, scores and E values differ from NCBI's.
- **When.** `tblastx` or `tblastn`, with `-subject` (not `-db`), and `-db_gencode` other than 1. `-db` searches, code 1 and `-query_gencode` match NCBI. In the 20 `tblastx` cases of one test profile, the 6 that use another subject code differ and 14 are identical. For example, one case gives 12,672 rows with NCBI and 13,644 with LOSAT.
- **Why.** NCBI searches the sequence with a code you did not ask for. LOSAT searches with the code you asked for. `tblastn` also accepts genetic code 32 (NCBI's `tblastn` does not list it).
- **Records.** [The decision record on genetic codes for local tblastx subjects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TBLASTX-LOCAL-SUBJECT-GENCODE.md), [the decision record on genetic code 32 for local tblastn subjects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-TLOSAN-LOCAL-GENCODE-32.md).

### `tblastx` word cutoff on large genome pairs

- **What you see.** LOSAT writes hits where NCBI BLAST+ 2.17.0 runs out of memory or time without output.
- **When.** `tblastx` with default options on ordinary genome pairs, when a cutoff in NCBI's statistics overflows a double. In one test (a 1.2 Mb query against a 5.5 Mb genome) NCBI crashed after 44 GB, or wrote nothing in 4 hours with 300 GB; LOSAT wrote 88,562 rows in 2 minutes. Where the value stays finite, both give the same bytes. With an extreme `-searchsp` and `-evalue` on small inputs, NCBI reports more rows than LOSAT.
- **Why.** NCBI turns the overflow into a cutoff of 1 through a C conversion that is undefined. LOSAT keeps the mathematically equal finite cutoff.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### One-hit gapped starts past the end of a sequence

- **What you see.** The same output as NCBI where NCBI does not crash, and a result where NCBI sometimes crashes.
- **When.** `blastp`, `blastx` or `tblastn` with `-window_size 0` (one-hit word finder) and a low `-threshold`, when the 11-letter window of a one-hit alignment passes the end of a sequence.
- **Why.** NCBI reads memory outside its buffer there. LOSAT reads the letters past the end as the end marker. Under valgrind, NCBI's output is LOSAT's output in all the cases of the test.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `blastn -db_hard_mask` when a sequence's first mask runs to its end

- **What you see.** LOSAT searches such a sequence before its mask, as NCBI searches every other piece. NCBI scans the masked letters past that piece and can report an extra alignment, lose the alignments of the next sequence, or crash.
- **When.** `-db_hard_mask` with a word size of 11 or more (megablast, blastn, dc-megablast), on a sequence whose first mask with letters runs to its end with nothing before it but one piece (a sequence masked entirely, or a masked poly-A tail). With `-task blastn-short` LOSAT refuses such a database (see [What LOSAT refuses](refusals.md)).
- **Why.** NCBI uses two lengths for one subject there. LOSAT gives NCBI's result wherever its scan past the piece finds nothing more.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### Titles made only of punctuation

- **What you see.** LOSAT writes the report. NCBI crashes.
- **When.** `blastn`, `tblastn` or `tblastx`, outfmt 0, for a subject whose definition line is commas, semicolons, tildes and spaces ending in a run of them (such as `, ,`) and that has hits. Also the tabular id of a database record `lcl|Subject_N , ,`.
- **Why.** NCBI's title clean-up reads past the end of the string. LOSAT stops the clean-up at the end (`, ,` gives `, `).
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### Query chunks that NCBI would split again

- **What you see.** LOSAT searches each chunk once. NCBI stops with a null-pointer exception (exit status 3).
- **When.** `blastn` with `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` set so that a chunk is long enough to be split again.
- **Why.** NCBI's failure is an unchecked condition. LOSAT's output equals NCBI's output at the largest overlap that does not split a chunk again, within the weak alignments at chunk boundaries that NCBI's own splits also lose.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `-culling_limit` near the largest 32-bit integer

- **What you see.** Limits of 2147483645 to 2147483647 cull nothing, as 2147483644 does. NCBI's 32-bit sum overflows and acts as a limit of 1 for the last two.
- **When.** Any program with `-culling_limit` at those values. Every other limit matches NCBI.
- **Why.** NCBI's behaviour is a signed overflow, which is undefined in C.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `gb`, `emb` and `dbj` ids

- **What you see.** LOSAT prints the id type that is stored. NCBI can print `gb|AF000002.1|` for a sequence stored as `emb|AF000002.1|`.
- **When.** Only in a custom tabular field list (`sseqid` with `-db`, `qseqid` with `-parse_deflines`, and with `-subject -parse_deflines` also `sseqid` and `sallseqid`) and in the ids of outfmt 5, 8, 9, 11, 12 and 13 to 16, and only when definition lines give one accession prefix two of these types. NCBI's core_nt has no such ids. The default tabular fields do not show it.
- **Why.** NCBI keeps the three types in one table for the whole run. LOSAT's output equals NCBI's with `OBJECTS_PACK_TEXTID=0`, which turns this off.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `-taxids 1`, and an accession on very many sequences

- **What you see.** `-taxids 1` (the root) with `-db` expands to every tax id of the taxonomy file and the search ends; NCBI searches without end. A `-seqidlist` or `-negative_seqidlist` naming an accession stored under more than about 500 sequences of a version 5 database works; NCBI stops with `Accessions to Oids lookup error`. The same two cases apply to `blastdbcmd -taxids` and `-entry`.
- **Why.** Both are failures of NCBI's lookup code. LOSAT gives what NCBI gives for a nearby input.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### Numbers in a PSSM, archive or XML file that overflow NCBI's `int`

- **What you see.** LOSAT reads the number with its exact value. NCBI wraps it. `psiblast` and `tblastn` read a PSSM (`-in_pssm`) whose values nest to any depth, where NCBI may crash with its default 8 MiB stack. A real number such as `1.55e-2147483647` is read as the smallest normal double.
- **When.** Only with files that hold such values: they do not occur in files that NCBI writes.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `rpsbproc`: ids and positions past 2^31

- **What you see.** LOSAT prints positions and orders PSSM ids without overflow (NCBI's binary prints `-2147483648` for a base at 2^31). A motif check that would read NCBI's list of motif positions past its end counts as failed, and the site is dropped.
- **When.** Only with data that reaches these limits. No CDD release has a PSSM id of 2^31.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### A taxonomy record past NCBI's mapping of `taxdb.btd`

- **What you see.** LOSAT writes what NCBI writes when NCBI completes. Whether NCBI completes depends on its memory map: it may crash on the same files.
- **When.** A damaged `taxdb.btd` whose record ends past NCBI's mapping, where the output does not write the bytes read there. Output that would write those bytes is refused (see [What LOSAT refuses](refusals.md)).
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `-db` output does not depend on `-num_threads`

- **What you see.** LOSAT's `-db` output is NCBI's output for `-num_threads 1`, at every `-num_threads`. NCBI's own `-db` output can depend on the thread count (with composition-based statistics, `blastp` ends the redo of a query early by the hit lists of its batch, and with more threads NCBI keeps one set per thread: 42 rows differed in one test).
- **`-mt_mode`.** `-mt_mode 0`, `1` and `2` are accepted by `blastn`, `blastp`, `blastx` and `tblastn` and do not change the output, which is always that of `-num_threads 1`. NCBI's `-mt_mode 1` writes the end of an outfmt 0 report, or each XML document, once per batch of queries; LOSAT does not.
- **Record.** [The decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).

## Platforms

### Last digit of a full-precision value on macOS and WASI

- **What you see.** The last (15th) digit of a value printed with all its digits can differ by 1 from the Linux and Windows value.
- **When.** macOS: the E value in `-outfmt 11` and `-outfmt 8` (tabular output is not affected). WASI (`wasm32-wasip1` and `wasm32-wasip1-threads`): the values of the ASN.1 PSSM of `psiblast -out_pssm` written with `-out_ascii_pssm`. Values rounded to 6 digits (XML, JSON) or printed with 2 or 3 digits (tabular) are not affected in practice.
- **Why.** The platform's maths library differs by 1 ULP. NCBI BLAST+ on macOS is expected to show the same (not tested). On WASI, NCBI BLAST+ 2.17.0 writes the same digit on an x86-64 CPU without FMA.
- **Record.** [The decision record on the last digit of full-precision values](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PLATFORM-LIBM.md).

### NCBI's own builds differ across platforms

- **What you see.** If you compare LOSAT with NCBI's own binary on Windows, macOS or Linux aarch64, NCBI's output may differ from its Linux x64 output within bounds recorded in the decision. LOSAT's output is the same bytes on every platform as on Linux x64.
- **Record.** [The decision record on differences between NCBI builds](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-NCBI-PLATFORM-VARIANCE.md).

## Command line and files

These cases are about what the program prints around the search, not the search results. They are described in [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

### Threads with `-subject`

- **What you see.** With `-subject`, NCBI BLAST+ ignores `-num_threads` and writes `'num_threads' is currently ignored when 'subject' is specified.` on the standard error. LOSAT uses the threads you ask for and does not print the line. Above the CPU count NCBI also writes `Number of threads was reduced to N to match the number of available CPUs`; LOSAT does not.
- **Why.** The output is the same at every thread count, so using the threads only makes the search faster. Every other byte of standard output and standard error equals NCBI's.

### Argument errors and help

- **What you see.** An error found while reading the arguments (an unknown option, a missing value, a value outside its range such as `-word_size 3`) prints LOSAT's message and exits with status 2. NCBI prints its USAGE text and the argument exception, and exits with status 1. `-help`, `-help-full` and `-h` print LOSAT's help; `-xmlhelp` prints NCBI's XML for the arguments LOSAT runs.
- **What does not change.** Errors that NCBI raises after it has read the arguments (`BLAST query/options error: ...`, input errors, engine errors) keep NCBI's condition, error prefix (`BLAST query/options error:`, `BLAST engine error:`, `BLAST Database error:`, `BLAST options error:`), and message, with LOSAT's user-manual footer (`Please refer to the LOSAT user manual.`).

### Version text

- **What you see.** `-version` writes `blastn: 2.17.0+` and a `Package:` line that names LOSAT. `-version-full`, `-version-full-xml` and `-version-full-json` give NCBI's forms with LOSAT's package; the lines in which NCBI's build describes itself (configuration, build signature, dates, revisions) are not written.

### A write that fails, a closed pipe, memory

- **Tabular output (6, 7) cannot be written** (`-out /dev/full`, a full disk): NCBI aborts. LOSAT reports the write error and exits with a non-zero status. For outfmt 0 LOSAT gives NCBI's message and exit status (`BLAST failed to write output`, 6).
- **A side file of outfmt 13 or 14 cannot be written**: NCBI reports it on the standard error and goes on, exit status 0. LOSAT stops with `BLAST failed to write output` (exit status 6).
- **Memory runs out**: NCBI catches it and writes `BLAST ran out of memory` (exit status 4). LOSAT aborts.
- **The reader of a pipe closes it** (`... | head`): NCBI is ended by the signal (exit status 141 in a shell, no message). LOSAT stops at once and reports the write failure for outfmt 0, 6 and 7 (exit status 6 for outfmt 0), or stops with exit status 0 for `blastdbcmd`, `dustmasker`, `segmasker`, `convert2blastmask` and for the version and help texts. A short report that LOSAT has written before the reader closes ends with status 0 where NCBI is ended later.
- **Standard output closed at the start** (`>&-`): NCBI fails on its first write. LOSAT writes the report to `/dev/null` and exits as the search ends.
- **`rpsbproc` cannot write its report**: NCBI's `rpsbproc` drops the writes and exits with 0. LOSAT stops with `rpsbproc failed to write output` (exit status 4).

These are error paths: when the output can be written, the output matches NCBI's.

### Compressed `-query`, `-subject` and `-in`

- **What you see.** A `-query`, `-subject` or `makeblastdb -in` file whose name ends in `.gz`, `.bz2` or `.zst` is read through gzip, bzip2 or zstd. NCBI's search programs read such a file as text and stop at its first line (`CFastaReader: Near line 1, there's a line that doesn't look like plausible data`).
- **When.** A search program inflates the whole file into an unnamed temporary file beside the input (under `TMPDIR` if the directory refuses it or has no room) before the search; the WASI command stops instead. `makeblastdb` streams. The name alone decides, in any case. A file that is not what its name says stops the run before the search (`-query: 'x.fa.gz' is named .gz but is not gzip data (...)`, exit status 1).
- **Why.** NCBI's own `makeblastdb` reads `.gz`, `.bz2` and `.zst` names, so LOSAT does the same for the searches. A plain file, or a compressed file read through `-query -`, is unchanged.

## Databases and tools

### `makeblastdb`

- A `-in` file made of several bzip2 streams or zstd frames is read to the end; NCBI's `makeblastdb` stops after the first and drops the rest silently. A single-stream file builds the same files as NCBI.
- An `-out` path with a space or a double quote, or an existing database that has to be replaced, is handled as one database. NCBI treats the path as a list of names and can keep volume files of the old database.
- Database files are compared with NCBI's `makeblastdb` 2.17.0 on Linux x86-64. On Windows and macOS LOSAT follows NCBI's source where the file layout depends on the platform; those files have not been compared with NCBI's builds there.

Records: [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) and [the decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).

### `blastdbcmd`

LOSAT writes correct output where NCBI's `blastdbcmd` writes letters or sequences that do not belong to the header or id you asked for. The six cases are `-long_seqids` with `-range` on the minus strand; `gnl|BL_ORD_ID|N` in a version 4 database of several volumes; `-mask_sequence_with` with `-strand minus`; two or more `%m` in one `-outfmt`; `-mask_sequence_with` on a protein without a range; a batch line without masks after a line with masks for the same sequence. In each, LOSAT's output is NCBI's own output for the neighbouring input that does not reach the defect.

Record: [The decision record on blastdbcmd output that NCBI gets wrong](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

## The library

### `-outfmt 11` and a definition line kept as read

- **What you see.** The command line stores the definition line of each query or subject as read, as NCBI does. A caller of the Rust library that hands over already-parsed records without a kept definition line gets `>` and the title, which differs from NCBI's archive for a definition line with extra spaces or control bytes.
- **Record.** [The decision record on the definition line in the BLAST archive](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-ARCHIVE-DEFLINE.md).
