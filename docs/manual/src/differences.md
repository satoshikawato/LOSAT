# Differences from NCBI BLAST+

LOSAT v0.2.0 writes the output of NCBI BLAST+ 2.17.0. This page lists every case in which it does not, with what you see, when it happens, and why the difference is acceptable or better. Each difference has a decision record on GitHub, linked from its item.

Cases that LOSAT refuses with an error are on [What LOSAT refuses](refusals.md), not here.

The differences you are most likely to meet in everyday use:

- [A local `-subject` search of `tblastx` and `tblastn` translates the subject with `-db_gencode`](#a-local--subject-search-of-tblastx-and-tblastn-uses--db_gencode); NCBI uses genetic code 1 there.
- [`-db` output is the same at every `-num_threads`](#-db-output-does-not-depend-on--num_threads), and [`-subject` searches use the threads you ask for](#threads-with--subject).
- [`tblastx` reports hits on large genome pairs](#tblastx-word-cutoff-on-large-genome-pairs) where NCBI BLAST+ runs out of memory or time.
- [A `-query`, `-subject` or `makeblastdb -in` file named `.gz`, `.bz2` or `.zst` is read through its decompressor](#compressed--query--subject-and--in).
- [On macOS, the last digit of a full-precision E value can differ by 1](#last-digit-of-a-full-precision-value-on-macos).
- [SAM output (`-outfmt 17`) records LOSAT's command line](#sam-output-the-cl-tag).
- [The version line and the error prefixes name LOSAT](#the-version-line-and-the-error-prefixes-name-losat), where NCBI writes its own version and `BLAST`.
- [`rpsbproc` keeps hits that tie in E-value and bit score in their input order](#rpsbproc-hits-that-tie); NCBI's order of them can differ.

The other items are rare cases in which NCBI BLAST+ crashes, reads outside its memory, or overflows an integer.

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

### `tblastx` with very many short queries against one subject

- **What you see.** LOSAT writes linked hits. NCBI BLAST+ 2.17.0 crashes (exit status 139) and writes nothing.
- **When.** `tblastx` against a single subject sequence (`-subject` with one sequence, or a database of one; other subjects of fewer than 3 bases in all count as none), when one query batch holds so many very short queries that the batch's average translated length rounds down to 0: for example a 1,200-nt query with 1,300 queries of 3 nt, or with `-strand plus` or `minus` only 300 of them. When the other subjects add at least one codon, LOSAT and NCBI give the same results (apart from the version line, which names LOSAT).
- **Why.** NCBI's statistics for linking hits then give a cutoff of minus infinity, and a sum with it overflows a 32-bit integer. LOSAT links the hits by NCBI's rules with the sums in 64-bit integers. No run that NCBI completes can check these rows: NCBI's nearest one, with a second subject that has no hits, links the hits differently.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### One-hit gapped starts past the end of a sequence

- **What you see.** The same output as NCBI (apart from the version line, which names LOSAT) where NCBI does not crash, and a result where NCBI sometimes crashes.
- **When.** `blastp`, `blastx` or `tblastn` with `-window_size 0` (one-hit word finder) and a low `-threshold`, when the 11-letter window of a one-hit alignment passes the end of a sequence.
- **Why.** NCBI reads memory outside its buffer there. LOSAT reads the letters past the end as the end marker. Under valgrind, NCBI's output is LOSAT's output in all the cases of the test.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `-db_hard_mask` when a sequence's first mask runs to its end

- **What you see.** LOSAT searches such a sequence before its mask, as NCBI searches every other piece. NCBI scans the masked letters past that piece and can report an extra alignment or a different one, lose the alignments of the next sequence, write one alignment many times (`tblastx`), or crash.
- **When.** `-db_hard_mask` on a sequence whose first mask with letters runs to its end with nothing before it but one piece (a sequence masked entirely, or a masked tail such as a poly-A tail): `blastn` with a word size of 11 or more (megablast, blastn, dc-megablast), and `blastp` and `blastx` at every word size. `tblastn` and `tblastx` read this in each reading frame: a frame whose first piece, from its first letter, ends at a mask with no other piece of the frame after it (masks at or near either end of a sequence). With `blastn -task blastn-short` LOSAT refuses such a database (see [What LOSAT refuses](refusals.md)).
- **Why.** NCBI uses two lengths for one subject there. In `blastn` LOSAT gives NCBI's result wherever its scan past the piece finds nothing more. In the other programs LOSAT searches the piece as NCBI searches it in the same database with the mask shortened until a letter follows the piece (one letter in `blastp` and `blastx`, a few bases in `tblastn` and `tblastx`). The masked letters stay masked in the alignments and the statistics, so LOSAT's result is NCBI's result for that database with the freed letters masked (`X`, or `N` in nucleotides). Where no such shorter mask exists (a mask of one or two bases at the end of a very short sequence), LOSAT's result is NCBI's own.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `-db_hard_mask` with a short mask at an end of a sequence (`tblastn`, `tblastx`)

- **What you see.** LOSAT searches the sequence and writes its alignments. NCBI crashes.
- **When.** `tblastn` or `tblastx` with `-db_hard_mask` on a sequence with a mask within its first two bases, or within the one or two bases after its last whole codon (for example [1200, 1202) of 1,202 bases), or with an empty mask at its start. Such a mask leaves a reading frame an empty part at its start.
- **Why.** NCBI then searches the rest of that frame as the whole frame and frees memory that it did not allocate, in every run. LOSAT searches the frame as a frame without a mask, and the masked bases stay masked in every frame, as for every other mask. In the cases checked, LOSAT's result is NCBI's result for the same database with that short mask dropped and its bases stored as `N`, which is also what NCBI gives where another mask of the sequence avoids the crash. NCBI's preliminary search reads such a stored `N` as a random base, where LOSAT reads the real base, so a result can differ from that run.
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

- **What you see.** Limits of 2147483645 to 2147483647 cull nothing, as 2147483644 does. NCBI adds 3 to the limit in a signed 32-bit integer; the sum overflows, and the last two limits act as a limit of 1.
- **When.** Any program with `-culling_limit` at those values. Every other limit matches NCBI.
- **Why.** NCBI's behaviour is a signed overflow, which is undefined in C. LOSAT computes the sum without overflow.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### Hit list sizes from 2^30

- **What you see.** A hit list size (`-max_target_seqs`, `-num_alignments`, or `-num_descriptions` in `-outfmt` 0 to 4) of 2^30 or more, or 2^30 - 25 or more with composition-based statistics, keeps every hit, as 2^30 - 26 does. NCBI computes the preliminary hit list size in a signed 32-bit integer, which overflows there: NCBI keeps 10 hits per query (without composition-based statistics, up to 2^31 - 51) or 2 to 48 (with them, from 2147483624), and crashes at the other sizes. With `-best_hit_overhang` or `-best_hit_score_edge` NCBI's best-hit filter counts twice that size in the same integer, which overflows from 2^30 - 50 (2^29 - 25 with composition-based statistics, 2^30 in ungapped searches): there NCBI keeps some alignments that the filter removes at smaller sizes, and LOSAT gives the result of the size just below.
- **When.** Gapped `blastn`, `blastp`, `psiblast`, `rpsblast`, `blastx`, `rpstblastn` and `tblastn` searches with such a size, and searches with `-best_hit_overhang` or `-best_hit_score_edge`. Ungapped searches and `tblastx` otherwise keep the size itself, as NCBI does. Every smaller size matches NCBI.
- **Why.** NCBI's behaviour is a signed overflow, which is undefined in C. LOSAT computes the sizes without overflow.
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

### Numbers and nesting in a PSSM, archive or XML file

- **What you see.**
  - A real number whose exponent does not fit in a 32-bit integer, such as `1.55e-2147483647`, is read with its exact exponent: the smallest normal double. NCBI's reader overflows the exponent and reads the largest double.
  - `psiblast` and `tblastn` read a PSSM (`-in_pssm`) whose values nest to any depth, and `psiblast -out_pssm` writes it back. NCBI reads and writes such values by recursion and can crash (exit status 139, no output) with its default 8 MiB stack; with a stack large enough it writes LOSAT's output.
  - When `psiblast` writes a PSSM (`-out_pssm`) that holds -2147483648 or -9223372036854775808, LOSAT writes the number. NCBI's writer overflows on these two numbers and prints `-8` or `-0`.
  - A database length of -9223372036854775808 is written as that number in BLAST XML and XML2 (`-outfmt` 5, 14 and 16), in the search strategy (`-export_search_strategy`) and in the BLAST archive. NCBI's writer prints `-0` there.
- **When.** Only with a file that holds such a number or nest (a PSSM given to `-in_pssm`, an archive or XML given to `rpsbproc`). A PSSM that `psiblast` makes from a FASTA query does not nest this way, and NCBI's writers do not write such exponents. A database has that length only with an alias file's `LENGTH` of 9223372036854775808 or with `-dbsize -9223372036854775808`.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `rpsbproc`: ids and positions past 2^31

- **What you see.** LOSAT prints positions and orders PSSM ids without overflow (NCBI's binary prints `-2147483648` for a base at 2^31). A motif check that would read NCBI's list of motif positions past its end counts as failed, and the site is dropped.
- **When.** Only with data that reaches these limits. No CDD release has a PSSM id of 2^31.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `rpsbproc`: hits that tie

- **What you see.** Hits of one query with equal E-value and bit score keep their order in the archive or XML. NCBI's `rpsbproc` orders them with its C++ library's unstable sort, so the order of such hits can differ, and so can which of two tied hits is `Specific`: `-m std` then prints the other hit's line, and `-m full` swaps their hit types. A hit whose E-value is not a number (`-nan`) comes after every hit whose E-value is a number; among hits of one E-value, one whose bit score is not a number comes last. NCBI's place for such a hit depends on the steps of its library's sort, and in such a query the representative and `Specific` hits, and their sites and motifs, can differ too.
- **When.** In the DOMAINS lines of `-m std` and `-m full`. In 74 reports of E. coli K-12 proteins and genome pieces against CDD v3.22, 18 differed, by 1 to 6 lines each; in 7 of them one query's `Specific` hit changed. None of them has a `-nan` value.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### A taxonomy record past NCBI's mapping of `taxdb.btd`

- **What you see.** LOSAT writes what NCBI writes (apart from the version line, which names LOSAT) when NCBI completes. Whether NCBI completes depends on its memory map: it may crash on the same files.
- **When.** A damaged `taxdb.btd` whose record ends past NCBI's mapping, where the output does not write the bytes read there. Output that would write those bytes is refused (see [What LOSAT refuses](refusals.md)).
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `blastx -use_sw_tback` in e-value mode

- **What you see.** LOSAT writes a report where NCBI BLAST+ 2.17.0 crashes, and the same report in every run and at every `-num_threads`. Where NCBI ends, LOSAT's report is NCBI's in the cases we ran, except that NCBI's rows can change from one run to the next.
- **When.** `blastx -use_sw_tback` with `-sum_stats false` (or a `-max_intron_length` of 1 to 4), on a query batch with a hit in a frame whose number is the batch's number of queries or more. The number is 6 × the query's place in the batch (counted from 0) plus 0 to 5 for frames +1, +2, +3, -1, -2, -3. With one query in a batch (for example a search of a single query, or a query that fills a batch of 10,002 letters, BLASTX's batch size, by itself), every hit outside frame +1 is such a case.
- **Why.** NCBI checks each new Smith-Waterman alignment against a list of kept alignments that it picks by that number, but it keeps one list per query. Past its lists it reads other memory: it crashes (with more than one thread, in every run we made), fails the batch with `Traceback failed`, or goes on with what it read. LOSAT checks such an alignment against an empty list. Its report equals NCBI's report of the same queries followed by enough queries without hits to keep NCBI inside its lists, in all 9 cases of our test in which NCBI ends.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `blastn -use_index true` with `-db_hard_mask`

- **What you see.** LOSAT reports the alignments in the unmasked parts of each masked sequence, the same in every run and at every `-num_threads`. NCBI BLAST+ 2.17.0 reports fewer or other alignments for such sequences, and with more than one thread its rows can change from run to run.
- **When.** `blastn -use_index true` with a MegaBLAST index of the database and `-db_hard_mask`, for sequences with masked letters. A sequence without masked letters gets NCBI's result, except for a match that runs on past its end into masked letters of the sequence after it in the index.
- **Why.** NCBI's search takes the index's seeds of a masked sequence from another part of it, or from the sequence after it, because it numbers the index's pieces of a sequence over the whole sequence but the pieces of a masked search from each unmasked stretch. LOSAT reads the index as if the masked letters matched nothing and searches each unmasked stretch with its own seeds, so no masked letter changes the result. With the default `-window_size 0` of megablast, LOSAT's rows were those of NCBI's search without the index with the same `-db_hard_mask` in every case of our tests (the test suite's cases and 308 random databases and queries). With a `-window_size` above 0 the index and the scan pick different seeds even without a mask, and LOSAT keeps the index's choice.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

### `-db` output does not depend on `-num_threads`

- **What you see.** LOSAT's `-db` output is NCBI's output for `-num_threads 1` (apart from the version line, which names LOSAT), at every `-num_threads`. NCBI's own `-db` output can depend on the thread count (with composition-based statistics, `blastp` ends the redo of a query early by the hit lists of its batch, and with more threads NCBI keeps one set per thread: 42 rows differed in one test).
- **`-mt_mode`.** `-mt_mode 0`, `1` and `2` are accepted by `blastn`, `blastp`, `blastx` and `tblastn` and do not change the output, which is always that of `-num_threads 1`. NCBI's `-mt_mode 1` writes the end of an outfmt 0 report, or each XML document, once per batch of queries; LOSAT does not.
- **Record.** [The decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).

## Platforms

### Last digit of a full-precision value on macOS

- **What you see.** The last (15th) digit of a value printed with all its digits can differ by 1 from the Linux and Windows value.
- **When.** macOS: the E value in `-outfmt 11` and `-outfmt 8` (tabular output is not affected). Values rounded to 6 digits (XML, JSON) or printed with 2 or 3 digits (tabular) are not affected in practice. The WASI modules and the Web engine give the Linux value.
- **Why.** The platform's maths library differs by 1 ULP. NCBI BLAST+ on macOS is expected to show the same (not tested).
- **Record.** [The decision record on the last digit of full-precision values](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-PLATFORM-LIBM.md).

### NCBI's own builds differ across platforms

- **What you see.** If you compare LOSAT with NCBI's own binary on Windows, macOS or Linux aarch64, NCBI's output may differ from its Linux x86_64 output within bounds recorded in the decision. LOSAT's tests compare its output on every platform with the Linux x86_64 bytes: on the release test set it is the same, apart from the macOS last digit above.
- **Record.** [The decision record on differences between NCBI builds](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-NCBI-PLATFORM-VARIANCE.md).

## Command line and files

These cases are about what the program prints around the search, not the search results. They are described in [the decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md).

### Threads with `-subject`

- **What you see.** With `-subject`, NCBI BLAST+ ignores `-num_threads` and writes `'num_threads' is currently ignored when 'subject' is specified.` on the standard error. LOSAT uses the threads you ask for and does not print the line. Above the CPU count NCBI also writes `Number of threads was reduced to N to match the number of available CPUs`; LOSAT does not.
- **Why.** The output is the same at every thread count, so using the threads makes the search faster without changing its output. Every other byte of standard output and standard error equals NCBI's, apart from the version line and the error prefixes, which name LOSAT, and `blastn`'s progress bar on a terminal (next section).
- **Memory with more than 64 threads.** LOSAT's memory grows with each thread, while NCBI BLAST+'s does not grow with `-num_threads` in a `-subject` search, which it runs on one thread. Above 64 threads LOSAT can then use more memory than NCBI BLAST+: in the v0.2.0 tests (development build `ff07c6a7d`, 64 logical CPUs), up to 1.21 times NCBI's peak memory at 96 threads and up to 1.44 times at 128 threads (`tblastn -subject`; also `blastx`, `blastn` and `blastp` with `-subject`, and `tblastn` with a database of one sequence), while LOSAT finished 8 to 69 times sooner.

### `blastn`'s progress bar on a terminal

- **What you see.** While `blastn` searches with `-num_threads` 2 or more against two subjects or more (the records of `-subject` or the sequences of `-db`), it draws a progress bar on the standard error: the elapsed time, a bar and the number of subjects done (`19/33`). When the search of each batch of queries ends, the bar is erased, so nothing of it stays on the screen: the report and the shell's prompt appear where they would without it. NCBI BLAST+ writes nothing there during a search.
- **When.** Only when the standard error is a terminal and, outside Windows, `TERM` is set and is not `dumb`. A standard error sent to a file or a pipe (`2> log.txt`, `2>&1 | less`) gets no byte of it. The library, pylosat and LOSAT Web never draw it.
- **What does not change.** The bytes of the standard output, the `-out` file and the exit status.

### Argument errors and help

- **What you see.** An error found while reading the arguments (an unknown option, a missing value, a value outside its range such as `-word_size 3`) prints LOSAT's message and exits with status 2. NCBI prints its USAGE text and the argument exception, and exits with status 1. `-help`, `-help-full` and `-h` print LOSAT's help; `-xmlhelp` prints NCBI's XML for the arguments LOSAT runs, with LOSAT's name and version in its description line.
- **What does not change.** Errors that NCBI raises after it has read the arguments (where NCBI writes `BLAST query/options error: ...`, input errors and engine errors) keep NCBI's condition, message and exit status. Their prefix names LOSAT (`LOSAT query/options error:`, `LOSAT engine error:`, `LOSAT database error:`, `LOSAT options error:`); see [the version line and the error prefixes](#the-version-line-and-the-error-prefixes-name-losat).

### Wording of messages

- **What you see.** Messages refer to LOSAT, not to NCBI's product: the last line of a query or options error reads `Please refer to the LOSAT user manual.`, and an error that NCBI reports as a C++ exception shows its message only, without NCBI's `NCBI C++ Exception:` header, source file, line, class name and function (`Error: Cannot convert string 'abc' to Int8 (position 0)`).
- **What does not change.** The message texts, NCBI's diagnostic codes such as `(106.15)`, the exit statuses and every search result. The error prefixes name LOSAT, as [the version line and the error prefixes](#the-version-line-and-the-error-prefixes-name-losat) describes.
- **Record.** [The decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md), exception 12.

### Version text

- **What you see.** `-version` writes `blastn: LOSAT 0.2.0` and a `Package:` line that names LOSAT (`Package: LOSAT 0.2.0, output compatible with blast 2.17.0`). `-version-full`, `-version-full-xml` and `-version-full-json` give NCBI's forms with LOSAT's version and package; the lines in which NCBI's build describes itself (configuration, build signature, dates, revisions) are not written.

### The version line and the error prefixes name LOSAT

- **What you see.** Where NCBI writes its own version, LOSAT writes its own name and version, and where NCBI starts an error with `BLAST`, LOSAT starts it with `LOSAT`.
  - After a program's name, `BLASTP (LOSAT 0.2.0)` where NCBI writes `BLASTP 2.17.0+`: in the first line of `-outfmt 0`, in the `# BLASTP (LOSAT 0.2.0)` line of `-outfmt 7`, in `<BlastOutput_version>` of `-outfmt 5`, in the `"version"` of `-outfmt 13` and `15`, in the `<version>` of `-outfmt 14` and `16`, and in the description line of the `-help` texts of `blastdbcmd` and `rpstblastn` and of `-xmlhelp`.
  - On its own, `LOSAT 0.2.0` where NCBI writes `2.17.0+`: in the `ident` of an `-outfmt 11` archive and of `-export_search_strategy`, in the `VN:` of the `@PG` line of `-outfmt 17`, in the first line of `-version` (`blastn: LOSAT 0.2.0`; also `dustmasker`, `segmasker` and `rpsbproc`, where NCBI writes its own versions 1.0.0 and 0.5.1), and in `BLAST database client, version LOSAT 0.2.0` of `blastdbcmd -help`. The banner that `rpsbproc` writes on the standard error reads `Post-RPSBLAST Processing Utility (LOSAT 0.2.0)`.
  - As a number, `0.2.0` in the program's `<version>` of `-xmlhelp`, where NCBI writes `2.17.0` (`1.0.0` for `dustmasker` and `segmasker`, `0.5.1` for `rpsbproc`).
  - The error prefixes are `LOSAT query/options error:`, `LOSAT engine error:`, `LOSAT database error:`, `LOSAT database creation error:`, `LOSAT options error:` and `LOSAT query error:`, where NCBI writes `BLAST query/options error:`, `BLAST engine error:`, `BLAST Database error:`, `BLAST Database creation error:`, `BLAST options error:` and `BLAST query error:`.
- **When.** In every run that writes one of these lines or errors. Every other byte is NCBI's: `-outfmt 6` and `10` (they have no such line), `8`, `9` and `12`, the output of `dustmasker`, `segmasker` and `convert2blastmask`, what `blastdbcmd` retrieves, the files of `makeblastdb`, the `Reference:` lines, `# BLAST processed N queries`, and the structure of the XML and JSON documents (DOCTYPE, namespaces, schema locations, element names). The messages after the prefix, NCBI's diagnostic codes and the exit statuses are unchanged.
- **Why.** The outputs of LOSAT name LOSAT, not NCBI's program. The program's name stays first in each line, so a tool that reads the first word as the program (Biopython's reader of `-outfmt 7`, for example) keeps working.
- **Comparing with NCBI.** To compare LOSAT's output with an NCBI file, map the NCBI file to LOSAT's wording with `python3 -I LOSAT/tests/losat_wording.py map-output <file>` from a LOSAT checkout, or leave these lines out of the comparison.

### SAM output: the `CL:` tag

- **What you see.** In `blastn -outfmt 17`, the `CL:` tag of the `@PG` line holds LOSAT's own command line, as it was run. `VN:` reads `LOSAT 0.2.0` where NCBI writes `2.17.0+` (see [the version line](#the-version-line-and-the-error-prefixes-name-losat)). The rest of the document is NCBI's, byte for byte. A search run through the library has no command line, so its `@PG` line has no `CL:` tag.
- **Record.** [The decision record on command-line behavior outside the search results](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md), exception 13. See also [Output formats](output-formats.md).

### A write that fails, a closed pipe, memory

- **Tabular output (6, 7) cannot be written** (`-out /dev/full`, a full disk): NCBI aborts. LOSAT reports the write error and exits with a non-zero status. For outfmt 0 LOSAT gives NCBI's message and exit status (`BLAST failed to write output`, 6).
- **A side file of outfmt 13 or 14 cannot be written**: NCBI reports it on the standard error and goes on, exit status 0. LOSAT stops with `BLAST failed to write output` (exit status 6).
- **A search program runs out of memory**: NCBI catches it and writes `BLAST ran out of memory` (exit status 4). LOSAT aborts.
- **The reader of a pipe closes it** (`... | head`): NCBI is ended by the signal (exit status 141 in a shell, no message). LOSAT stops at once and reports the write failure for outfmt 0, 6 and 7 (exit status 6 for outfmt 0), or stops with exit status 0 for `blastdbcmd`, `dustmasker`, `segmasker`, `windowmasker`, `convert2blastmask`, `makembindex -show_filters` and for the version and help texts. A short report that LOSAT has written before the reader closes ends with status 0 where NCBI is ended later.
- **Standard output closed at the start** (`>&-`): NCBI fails on its first write. LOSAT writes the report to `/dev/null` and exits as the search ends.
- **`rpsbproc` cannot write its report**: NCBI's `rpsbproc` drops the writes and exits with 0. LOSAT stops with `rpsbproc failed to write output` (exit status 4).

These are error paths: when the output can be written, the output matches NCBI's.

### Compressed `-query`, `-subject` and `-in`

- **What you see.** A `-query`, `-subject` or `makeblastdb -in` file whose name ends in `.gz`, `.bz2` or `.zst` is read through gzip, bzip2 or zstd. NCBI's search programs read such a file as text and stop at its first line (`CFastaReader: Near line 1, there's a line that doesn't look like plausible data`).
- **When.** A search program inflates the whole file into an unnamed temporary file beside the input before the search. If that directory refuses the file or has no room, the native program uses the temporary directory (`TMPDIR`; see [Environment variables](environment.md#temporary-files)), and the WASI command, which has no temporary directory, stops with an error. `makeblastdb` streams. The name alone decides, in any case. A file that is not what its name says stops the run before the search (`-query: 'x.fa.gz' is named .gz but is not gzip data (...)`, exit status 1).
- **Why.** NCBI's own `makeblastdb` reads `.gz`, `.bz2` and `.zst` names, so LOSAT does the same for the searches. A plain file, or a compressed file read through `-query -`, is unchanged.

## Databases and tools

### `makeblastdb`

- A `-in` file made of several bzip2 streams or zstd frames is read to the end; NCBI's `makeblastdb` stops after the first and drops the rest silently. A single-stream file builds the same files as NCBI.
- `-gi_mask` on a database in one volume writes an alias file that names the volume by its file name, so the database opens and the build ends with exit status 0 and its metadata file. NCBI's `makeblastdb` names the database by its full path there, which no program can open, and ends with `DB alias files are mutually recursive` (exit status 2). Every other file is NCBI's, and NCBI's programs read the database and its masks.
- An `-out` path with a space or a double quote, or an existing database that has to be replaced, is handled as one database. NCBI treats the path as a list of names and can keep volume files of the old database.
- Database files are compared with NCBI's `makeblastdb` 2.17.0 on Linux x86-64. On Windows and macOS LOSAT follows NCBI's source where the file layout depends on the platform; those files have not been compared with NCBI's builds there. The order of `-randomize` is the Linux builds' on every platform: NCBI's Windows and macOS builds shuffle with their own C++ libraries.

Records: [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md) and [the decision record on database search](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).

### `blastdbcmd`

LOSAT writes correct output where NCBI's `blastdbcmd` writes letters or sequences that do not belong to the header or id you asked for. The six cases are `-long_seqids` with `-range` on the minus strand; `gnl|BL_ORD_ID|N` in a version 4 database of several volumes; `-mask_sequence_with` with `-strand minus`; two or more `%m` in one `-outfmt`; `-mask_sequence_with` on a protein without a range; a batch line without masks after a line with masks for the same sequence. In each, LOSAT's output is NCBI's own output for the neighbouring input that does not reach the defect.

Record: [The decision record on blastdbcmd output that NCBI gets wrong](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-BLASTDBCMD-DEVIATIONS.md).

### `windowmasker`: an empty id in `-ids` or `-exclude_ids`

- **What you see.** LOSAT writes NCBI's error `bad id: ...: ignoring` for the line, skips it, and masks (or counts) as NCBI does for the file without that line. NCBI BLAST+ 2.17.0 writes the same error and then crashes (exit status 139) without output.
- **When.** `windowmasker` with an `-ids` or `-exclude_ids` file matched as text (without `-text_match F`), in masking or with `-mk_counts`, when a line gives an empty id: a line that starts with a space or a tab, or a line `>`, `|`, `>|`, `>>` or `>>|` (alone or before a space or a tab).
- **Why.** NCBI's error says that it ignores the id, but NCBI then stores the id far past the end of its list of ids, which is undefined in C++. LOSAT does what the error says.
- **Record.** [The decision record on NCBI BLAST+ defects](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md).

## The library

### `-outfmt 11` and a definition line kept as read

- **What you see.** The command line stores the definition line of each query or subject as read, as NCBI does. A caller of the Rust library that hands over already-parsed records without a kept definition line gets `>` and the title, which differs from NCBI's archive for a definition line with extra spaces or control bytes.
- **Record.** [The decision record on the definition line in the BLAST archive](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/product_decisions/PD-LOSAT-ARCHIVE-DEFLINE.md).
