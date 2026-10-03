# S08 result notes, range D (tblastx `-outfmt 7` and the shared parts of `-outfmt 6`)

Result table: `result_D.tsv` (60 rows, same order as `D.tsv`). Counts: ported 31, n/a 11, faithful 8, reused 4, GAP 3 (rows 2, 17, 37), exception 1 (row 8), rejected 1 (row 38), deferred 1 (row 51). No UNSURE rows (the three UNSURE items of `D_notes.md` section 5 are resolved, see below).

## Provenance of the binary (read this first)

`/home/kawato/.cache/losat-web-gui-target/s08/native/release/LOSAT` was rebuilt at 18:53:31 (while I was working) from a working tree that has uncommitted edits on top of `0d533ba76`: `algorithm/blastn/input.rs` and `algorithm/tblastx/blast_engine/run_impl.rs` (`search()`: `-culling_limit` rejection, rejection of non-IUPAC residues and of records without residues, `U` read as `T`). About 445 of my ~490 comparisons ran before the rebuild with the commit build; 45 ran after it (`allval_*`, `q1_allval_*`, `e1_*`, `pq2_*`, `em_ws*`, `blank_in`, `lead_nl`, `rnd_*`, `rndb_*`, `se_*`, `th4_*`), none of which has a U, a non-IUPAC letter, a record without residues or a culling limit, so the two builds behave alike on them. The ENGINE-1 evidence was also reproduced before the rebuild (`pre_*`, `cod_*`, `qn3_*`). `rscratch_D/out/*.losat.*` keep the bytes. All row classifications refer to the committed code (`git show 0d533ba76:...`); line numbers are those of `0d533ba76` (checked against the saved copy `rscratch_D/run_impl_0d533.rs`; `report.rs` is identical to the working tree).

## Method

`rscratch_D/cmp.sh NAME args...` runs NCBI 2.17.0 `tblastx` and `LOSAT tblastx` with the same arguments and compares stdout, stderr, exit status and the `2>&1` interleaving (`out/NAME.{ncbi,losat}.{out,err,rc,both}`). Inputs are under `rscratch_D/in/` (copies of the inventory's inputs plus new ones). About 490 comparisons were run; all rows of the inventory that depend on bytes were exercised at least through the frozen fixtures (`tblastx.*.7`, `warnings.*`, `env.*`, `hitlist.*`, `empty.*`) and by live runs.

Batteries with equal stdout, stderr and rc (except where listed below):
- full pair LC738874 x LC738875 outfmt 7 (4319 hits, 322018 bytes); `tblastx_multi`, `tblastx_batch`, `tblastx_ambig`, `tblastx_code4` (query gencode 4) pairs in outfmt 6 and 7; `-num_threads 4` (only NCBI's thread warning differs: approved exception).
- `tblastx_many_query` x `tblastx_many_subject` (260 subjects) with `-max_target_seqs` 1, 2, 3, 4, 5, 6, 10, 100, 259, 260, 261, 1000 and default, outfmt 7.
- 121 tiny/odd queries each alone (`-evalue 1e-100`): all 64 3-nt codons, 1..40 nt random, homopolymers, repeats, stops, IUPAC-only, N-padded, lowercase; and 119 of them concatenated into one file, alone / before q1 / after q1.
- all-N / 2-nt / poly-A / ACGT x75 / IUPAC-only alone, between valid queries, in a second batch (mb1..mb5, `BATCH_SIZE` 900/3000/12000 on a 10-query random file), long-label warnings (25 bytes + `.. `), empty / blank / whitespace-only query files (`Query is Empty!`), `-out file`, subject path spellings, `-max_target_seqs 1` with an empty query and with an unsearched query (warning order).
- `-outfmt 6` of the same inputs; `-outfmt 7` rows equal the non-comment lines of `-outfmt 6` in NCBI and LOSAT.

## GAP rows

### Row 2 (IsIStreamEmpty, pipe branch) - low
The seekable branch is ported (row 1). `-query -` is rejected by the argument parser (`blastinput/value_parsers.rs:539`, "a file path is required; stdin is not implemented", clap exit 2). A pipe passed as a path is not recognised:
```
printf '' | tblastx -query /dev/stdin -subject in/s1.fa -outfmt 7
  NCBI : stdout "# BLAST processed 0 queries", stderr empty, rc 0
  LOSAT: stdout empty, stderr "Warning: [tblastx] Query is Empty!", rc 0
printf '  \n' | tblastx -query /dev/stdin -subject in/s1.fa -outfmt 7
  NCBI : stdout empty, stderr "BLAST engine error: Empty CBlastQueryVector", rc 3
  LOSAT: stdout empty, stderr "Warning: [tblastx] Query is Empty!", rc 0
```
Only reachable through `/dev/stdin`, FIFOs or process substitution. If named pipes are out of scope, treat as rejected-equivalent.

### Row 17 (s_ReplaceLocalId: id columns) - low-medium, outfmt 6 and 7
Cause: ids come from `bio`'s `id()` (`run_impl.rs:1298,1302`: `split_whitespace().next().unwrap_or("unknown")`), and `check_report_titles` (`run_impl.rs:889`) returns early when every format is tabular (outfmt 6), and for outfmt 7 it checks bio's id+desc, not the raw defline. Input: 1500 nt of LC738874 as the query (`rscratch_D/in/defl3/*.fa`), `-subject in/s1.fa -outfmt 6 -evalue 1e-20`, first row, first column:

| defline | NCBI | LOSAT (outfmt 6) |
|---|---|---|
| `>` (empty) | `Query_1` | `unknown` |
| `>  leading space title` | `leading` | `unknown` |
| `>\tx title` | `x` | `unknown` |
| `>x<0x01>y title` | `x` | `x<0x01>y` |
| `>x<NBSP>title` (c2 a0) | `x<NBSP>title` | `x` |

Subject side (`-subject` with the same deflines, outfmt 6): `>  lead title`: NCBI `lead`, LOSAT `unknown`; `>x<0x01>y t`: NCBI `x`, LOSAT `x<0x01>y`; empty defline: NCBI `Subject_1`, LOSAT `unknown` (the deferred `Subject_N` item). With outfmt 7 the first four query cases and the first three subject cases stop with LOSAT's rejection ("query record 1 has a defline that is empty, starts with white space or has a control character ..."), but the NBSP/Unicode-space cases are not rejected (see row 37) and the subject NBSP case differs too. The stderr label of the Karlin warning has the same differences with outfmt 6 (row 12): empty defline NCBI `Warning: [tblastx] Query_1: Could not...` vs LOSAT `Query_1 : Could not...`; `>   leading title` NCBI `Query_1 leading title:` vs LOSAT `Query_1    leading title:`; `>x<TAB>TAB title` NCBI `Query_1 x:` vs LOSAT `Query_1 x TAB title:`.
Fix idea: run the same title check for outfmt 6 (BLASTN rejects such deflines in every format) and run it on the raw defline bytes.

### Row 37 (`# Query:` title) - low, outfmt 7
`bio` splits id/description at the first `char::is_whitespace` and `fasta_defline` re-joins them with a plain space, so tab, VT, FF and any Unicode white space become `' '` before `check_report_titles` looks for control characters and non-ASCII bytes. Input: 300 nt as query, `-subject in/s1.fa -outfmt 7 -evalue 1e-100`, exit 0 and empty stderr in both:
```
>x<TAB>TAB title   NCBI "# Query: x"             LOSAT "# Query: x TAB title"
>x<VT>TAB title    NCBI "# Query: x"             LOSAT "# Query: x TAB title"
>x<FF>TAB title    NCBI "# Query: x"             LOSAT "# Query: x TAB title"
>x<NBSP>title      NCBI "# Query: x<c2 a0>title" LOSAT "# Query: x title"   (also NEL c2 85, EM SPACE e2 80 83)
```
A tab after a space (`>x <TAB>title`), 0x01 or DEL in the title, a leading space/tab are rejected (fine); a trailing tab is equal.

## Other differences found while testing (not tied to a D row)

- ENGINE-1 (medium; linking/e-values; outfmt 6 and 7): a query of 3 or 4 nt (so that some of the 6 frames are missing but not all) placed BEFORE another query in the same batch changes the later query's HSPs in LOSAT but not (or differently) in NCBI. Repro (committed build): `printf '>a\nACG\n' ; q2 3000 nt` as `rscratch_D/in/pq2_n3.fa` vs `-subject in/ss.fa -outfmt 6`: NCBI output is identical to q2 alone (first rows `q2 s1 37.931 29 ... 0.032 23.0` and `q2 s1 77.778 9 ... 0.032 20.3`, two linked HSPs); LOSAT prints those two HSPs with e-values 0.57 and 3.6 (not linked) and the `q2 s2 ... 0.11` row first. Same HSP set, other sum-statistics e-values and order. Also with `in/pre_ACG.fa` (ACG + 4500 nt of q1): NCBI 42 hits at `-evalue 1e-100`, LOSAT 43 (extra row `q1 s1 100.000 49 0 0 3 149 3444 3298 0.0 116`; NCBI prints that HSP with e-value 9.38e-30 at the default e-value, LOSAT 0.0). Prefix of 1 or 2 nt, of 5 or 6 nt, or the short query AFTER the long one: equal. NNN and NNNN prefixes also differ, so it is not the Lambda of the 1-aa context (all 16 tested codons differ, also those whose contexts have the ideal Lambda). All frames missing (1-2 nt) or none missing (5+ nt) are equal; partially missing frames differ. I did not find the cause (the compact context list of `prepare_lookup_query`, `lookup/backbone.rs:334-420`, omits zero-length frames; `compute_avg_query_length_ncbi` and the context offsets agree with `blast_setup_cxx.cpp:55-90` on paper). Reaches rows 46/47 only by input, not by code.
- O-2 (outside the IUPAC input scope; being rejected in the working tree): characters that are not nucleotide letters (`X`, digits, `*`, `?`, `-`) in a query: NCBI ignores them with `FASTA-Reader: Ignoring invalid residues at position(s): On line 200: 11` on stderr and searches the shortened sequence; committed LOSAT keeps them silently (different hits, e.g. `bad2 s1 75.000 8 2 0 26 3 ...` vs `# 0 hits found`). A query of only `X` gives NCBI `BLAST engine error: Warning: Sequence contains no data ` rc 3 (LOSAT: rc 0).
- O-3 (U is an IUPAC letter; being fixed in the working tree by `with_u_as_t`): queries or subjects with `U`: `in/q1U.fa` (q1T with T->U): NCBI 70 hits, committed LOSAT 54; `in/s1U.fa` as subject: NCBI 70, LOSAT 40; a 45-nt `ACGUACGU...` query: NCBI 0 hits, LOSAT 1 hit (`U s2 7.692 13 12 0 44 6 1985 2023 0.15 18.0`).
- O-4 (query reading): a query file with blank lines before the first `>` (`printf '\n\n>x\nACGT\n'`) is processed by NCBI (rc 0, normal block); LOSAT stops with `Error: failed to read query FASTA ...: Expected > at record start.` rc 1. A sequence without any `>` line (`nodef.fa`) is a query with an empty title in NCBI (`# Query: ` + `# 0 hits found`, rc 0); LOSAT stops with the same error. Not in the explicit-rejection list; belongs to the query-reading ranges.
- O-5 (low): under extreme options, `-evalue 100000 -threshold 8 -seg no` on `in/q3.fa` x `in/ss.fa` (48604 rows, 2-aa HSPs), 10 pairs of exactly tied HSPs (same e-value, score, subject offsets and frame-relative query offsets, different query frames) come out in the other order (rows `q3 s2 100.000 2 0 0 540 535 20 15 9447 7.9` and `... 538 533 ...`). The sort of row 54 is as NCBI; the pre-sort order of fully tied HSPs is the engine's. Default options and the weaker `-evalue 1000 -threshold 8` run (16242 lines) are equal.
- O-6 (cross-range, in scope option; fixed in the working tree by a rejection): in the committed build `-culling_limit 1|2|5` with outfmt 7 gives `# 0 hits found` for every query and exit 0 while NCBI prints 29 hits for the first query of `tblastx_multi` (`-culling_limit 0` is equal). `-word_size 2` and `-window_size 0` are explicit rejections (equal to the COMMON list / pre-existing argument validator).

## Inventory corrections and resolved UNSURE items

- D_notes section 5 UNSURE 1 (validity of poly-A, ACGT repeat, IUPAC-only, < 3 nt): resolved, LOSAT equals NCBI on all tested compositions (row 47).
- UNSURE 2 (`-max_target_seqs` cut between equal e-values): resolved by the 260-subject fixtures and the live runs (row 25).
- UNSURE 3 (empty-sequence record first in a later batch): not testable in LOSAT (deferred); in the committed code any record without residues gives no error and a header-only block with the Karlin warning (`donly.fa`, `e2.fa`).
- Inventory row 15/23/54/58/59 said the common.rs writer was reusable/faithful; the final TBLASTX code does not call it at all (`write_output_ncbi_order_evalue_hsp_order_to_writer` is only re-exported at `blast_engine/mod.rs:71`); the hit list and sort are `report::final_hit_order` (uses the shared `HitList`).
- Row 18/29: the inventory said identities were already faithful; the final code recomputes identity, mismatches and pident from the displayed rows (`set_displayed_identities`), which matters for IUPAC codons only.
- Row 17 was marked `divergent` for outfmt 6/7 deflines in the inventory; S08 rejected the cases only for outfmt 0/7 and only when bio keeps the offending byte, hence the residual GAPs above.
- Row 3: the output of all batches is written after the last batch is searched, not after each batch. Bytes and 2>&1 order are equal; only the failure-in-a-later-batch case differs (deferred input only).

## "What the port must do" (D_notes section 3) status
1 dispatch on the format: done (row 13). 2 one block per query for all queries: done (row 60). 3 block lines: done (rows 32-36, 41). 4 footer: done (rows 43/44). 5 per-batch validity and unsearched queries: done (rows 46/47). 6 warnings only for unsearched batches: done (rows 12, 49, 50; outfmt 6 odd-defline labels are the row 17 GAP). 7 `# Query:` from the reader title: done for ordinary deflines, GAP row 37 for tab/VT/FF/Unicode white space. 8 `-max_target_seqs`: done (row 25). 9 `Query is Empty!`: done (row 1). 10 write failure: approved exception (row 8). 11 FormatProbe hooks: present in `write_tabular` (`probe.begin/end` around each row, `report.rs:583-610`).
