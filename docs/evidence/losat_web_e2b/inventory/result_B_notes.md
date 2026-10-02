# Result notes, range B (tblastx `-outfmt 0` alignment section, showalign.cpp / align_format_util.cpp)

Result TSV: `result_B.tsv` (61 rows, one per row of `B.tsv`). Scratch: `rscratch_B/`.
NCBI source: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++` (598d8ae6). Final LOSAT code: commit `0d533ba76`.

## Counts

ported 39, reused 12, n/a 8, rejected 2. faithful 0, exception 0, deferred 0, **GAP 0**, UNSURE 0.

(There is no `exception` row: the approved `-db_gencode` exception concerns the SEARCH. The DISPLAY of a non-default `-db_gencode` equals NCBI's local `-subject` display, see rows 5 and 9 below.)

## Which binary I used (important)

The shared binary `/home/kawato/.cache/losat-web-gui-target/s08/native/release/LOSAT` was replaced at 18:53 JST during my work by a build that contains UNCOMMITTED working-tree changes of the parent session (`git status` then showed modified `blastn/input.rs`, `tblastx/blast_engine/run_impl.rs`, `tests/tblastx_regression_fixtures.py`, new `input.rna_*` fixtures). The first part of my runs used the 18:22 binary. To have a result that belongs to commit `0d533ba76` alone, I exported that commit with `git archive` into `rscratch_B/final_src`, built it with `RUSTUP_TOOLCHAIN=1.92.0-x86_64-unknown-linux-gnu CARGO_TARGET_DIR=rscratch_B/final_target cargo build --release --offline` (same size as the 18:22 binary, 4604232 bytes), and re-ran the whole differential battery below with `rscratch_B/final_target/release/LOSAT` (scripts `cmpf.sh`, outputs in `rscratch_B/final_runs/`, summaries `final_rnd_results.txt`, `final_misc_results.txt`). Every statement below that names "identical" is from that binary unless it says otherwise.

## Method

For each of the 61 rows I read the NCBI lines (showalign.cpp, alnvec.cpp/.hpp, alnmap.cpp, Seq_align.cpp, align_format_util.cpp, blast_format.cpp, blast_aux.cpp, blast_filter.c, blast_query_info.c, setup_factory.cpp) and the final LOSAT code (`algorithm/tblastx/report.rs`, `report/pairwise.rs` TBLASTX section, `run_impl.rs`), then ran NCBI 2.17.0 `tblastx` and LOSAT on small inputs and compared the `-out` files byte for byte (`rscratch_B/cmp.sh`, `cmpf.sh`).

Differential battery (outfmt 0, NCBI 2.17.0 vs LOSAT `0d533ba76`):

| set | what | result |
|---|---|---|
| `rnd/r1..r360` (`gen.py`) | 360 random pairs: a window around a random HSP of the LC738874/LC738875 pair (200-2400 nt), random reverse complement of query and subject, 0-3 low-complexity inserts in the query (some touching), IUPAC substitutions in query/subject, N islands, lowercase input; options none, `-query_gencode 4`, `-query_gencode 2`, `-seg no`, `-evalue 10`, `-threshold 14 -window_size 60`, `-query_gencode 11 -seg no`, `-max_target_seqs 1` | 360 of 360 byte-identical (17182 HSP blocks, all 36 frame pairs, 7490 `Expect(n)` blocks, 411 Query rows with lowercase, 1909 Sbjct rows with B/Z/J/X) |
| `w/pp pm mp mm p10 amb amb2 lc dup two` | the B-range oracle cases (`+/+`, `+/-`, `-/+`, `-/-`, coordinate width at 10000, IUPAC, lowercase input, two subjects with the same title, two queries x two subjects) | identical |
| `ttl/*` | 26 subject deflines: no title, trailing spaces/punctuation, `gi|123|gb|X1.1| t`, `lcl|abc t`, `[organism=..]`, `TPA:`/`MAG:`/`UNVERIFIED:`/`PREDICTED:`, 130-char word, 60-byte words, `gnl|`, `sp|` | identical (the 27th, a non-ASCII byte, is the documented rejection) |
| `big/b1..b7` | 12-140 kb identical or reverse-complemented windows, up to 4000-residue HSPs, 23238 HSP blocks | identical except `b5` (see observation 2: a tie in the sum-statistics linking, also in outfmt 6) |
| `multi/m1..m8`, `m9`, `m10` | 4 queries x 4 subjects built from random pairs, with `-max_target_seqs 2` | identical |
| `-evalue 1000 / 100000 / 1e-50` on 5 pairs; `-seg yes`, `-seg "10 1.8 2.0"`, `-seg "20 2.0 2.5"` on 7 pairs; `lcx/*` low-complexity self pairs | E-value cascade, SEG parameters | identical |

The frozen fixtures that cover the rows are listed in the TSV `evidence` column (mainly `tblastx.lc.0`, `tblastx.multi.0`, `tblastx.ambig.0`, `tblastx.code4q.0`, `tblastx.segno.0`, `tblastx.many.0`).

## Rows that need a remark

* **5, 9 (`-db_gencode`).** NCBI's display translates the subject row with the given `-db_gencode` even for a local `-subject` (oracle: `tblastx -query w/pm_q.fa -subject w/pm_s.fa -db_gencode 4` prints the code-4 subject row). LOSAT's display does the same (`displayed_rows`, `db_code`). LOSAT also applies the code to the SEARCH (approved exception), so scores, HSP boundaries and `Identities` differ from NCBI's local `-subject` output in that case (`rscratch_B/o_gcs4.*`, `o_gcr3/9/21.*`); that is the exception, not a display defect. Fixture `tblastx.code4.0` compares with NCBI's database-mode output.
* **12, 28 (frames).** NCBI recomputes the frame from the row coordinates and the full sequence length (`s_GetFrame`); LOSAT prints the engine frames (`Hit.query_frame`, `TblastxHsp.subject_frame`). They agree for every HSP seen (queries and subjects of all lengths mod 3, all 36 frame pairs).
* **22 (title).** Reuse of `ncbi_nucleotide_title`; subjects whose NCBI title cannot be reproduced are stopped by `check_report_titles` (`run_impl.rs:889`) with: `subject record N has an HTML character reference (such as &amp;) ...`, `subject record N has a defline of punctuation that NCBI BLAST+ reads past its end ...`, and (queries and subjects) `... has a defline that is empty, starts with white space or has a control character or a non-ASCII byte ...`.
* **26 (bit score > 99999).** The `%5.3le` branch did not occur in any run (largest HSP 7014 bits, a 140 kb self comparison); the function is shared and unchanged, and it is exercised by the other programs' fixtures.
* **48, 49, 60 (masks).** Code reading plus oracle:
  * `s_LocalQueryData2Packed_seqint` (setup_factory.cpp:95-119) builds kTarget as `[0, len]` for a whole-sequence query location and `[0, len-1]` for an interval, so the inventory's "[0, L-1]" is only one of the two cases. For SEG masks of a translated query neither is reachable: after `BlastMaskLocProteinToDNA` a plus-frame mask ends at most at L-3 and a minus-frame mask never starts below 3 (frame -1) or ends below L-1 (frames -2, -3), so `range == kTarget` never happens. LOSAT's `from == 0 && to == query_length - 1` test (`report.rs:403`) is therefore dead but harmless.
  * A one-residue run on a minus frame gives from = to + 1: NCBI's `Map()` returns the whole target for an empty range and the mask is dropped; LOSAT skips `from > to`. Same result. SEG itself does not produce one-residue runs (window 12), so no oracle case exists.
  * Touching intervals `[a,b],[b+1,c]` (or sharing the endpoint `[a,b],[b,c]`) are not merged by `BlastSeqLocCombine(link 0)` (blast_filter.c:993 `stop > left`); on a minus frame the residue `b` would then be displayed uppercase between two lowercase runs. The standalone `segmasker` output for LC738874 (`scratch_B/q6.seg`) contains such pairs (e.g. frame -3, 560-571 and 572-584), but I found no output in the NCBI library path that shows the pattern: in the two 140 kb self comparisons (`big6`, `big7`; 11703 lowercase Query rows) there is no lower-UPPER-lower pattern at all, and LOSAT equals NCBI. So this path is verified only by code reading (identical `combine_masks` rule); I do not consider it UNSURE because the rule is the same code and no input produced a difference.
* **50 (subject masks).** 0 of 4837 Sbjct rows in `tblastx.lc.0`, and lowercase subject input (`w/lc`, 15% of the random subjects) prints uppercase subject rows in both programs.

## Observations outside range B (not GAP rows of this range; for the other ranges / the parent)

1. **`-culling_limit` differs from NCBI in outfmt 6 and 0 (HSP selection, engine).** On `rscratch_B/rnd/r5_{q,s}.fa`: NCBI `-culling_limit 1` prints 10 outfmt-6 rows, LOSAT 0; `-culling_limit 2`: NCBI 12, LOSAT 10; `-culling_limit 3`: both 12. Over r1..r12 LOSAT's limit N behaves like NCBI's N-1 (r3: NCBI c1 = 18 rows = LOSAT c2). Command: `LOSAT tblastx -query rnd/r5_q.fa -subject rnd/r5_s.fa -outfmt 6 -culling_limit 1` vs `tblastx` (NCBI). `-culling_limit` is in the option scope of the inventory and is also listed under S08+ (`session_s08p_e2e_protein_options.md`, TBLASTX options), so it may already be planned; in outfmt 0 it shows as `***** No hits found *****` for LOSAT (`o_o5_culling_limit1.*`).
2. **Tie in the sum-statistics linking (engine, outfmt 6 as well).** Query = subject = `LC738874.fasta` residues 20001-100000 (`rscratch_B/big/b5_{q,s}.fa`), `-outfmt 6`: NCBI rows 3612 and 3628 are `23655 23684 1558 1587 3.80e-99` and `23653 23682 1558 1587 1.62e-93`; LOSAT has the two query intervals exchanged (`23653..23682` with `3.80e-99`, `23655..23684` with `1.62e-93`). The two HSPs have equal score 42, the same subject interval (frame +1) and the same frame-relative query offset (7884) in frames +1 and +3; the choice of which one enters which linked set differs. The outfmt 0 text differs only by this swap (`Frame = +3/+1` / `+1/+1` and `Expect(12)` / `Expect(11)`). Not seen in any of the 360 random pairs, the 140 kb self comparisons, or the genome pair.
3. **Input letters (range C/E territory; display rows depend on it).** At commit `0d533ba76`: a FASTA nucleotide `U`/`u` is searched differently from NCBI (`rscratch_B/odd/q_U.fa` vs `odd/s_ok.fa`, outfmt 6 first row `... 1.30e-30 113` in LOSAT vs `... 6.35e-30 111` in NCBI; same for a subject `U`), and non-IUPAC letters (`X`, `x`, `*`, `-`, `J`, `O`, `Z`, digits, `.`, `_`) are kept by LOSAT without a message, while NCBI drops them with `FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 41-42` (or `CFastaReader: Hyphens are invalid ...`) and then reports a shorter sequence (query `Length=312` vs `Length=320` for 8 `X` in a 320-nt query). `B` is identical. The working tree of the parent session (uncommitted) already contains a change in this area: with the later binary `U` gave identical output and `X` was stopped with `Error: query record 1 (oq) has 'X' at residue 41, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not supported by LOSAT's TBLASTX ...`. My rows 31 and 32 only concern the display of the Bioseq letters (`display_base`: U as T, X as N, other as gap), which is as NCBI; the difference above is in the reading of the FASTA.
4. `-word_size 2` is stopped by clap (`unsupported TBLASTX word_size: only 3 is implemented`) and `-num_threads N` prints no thread warning (approved); both were seen in passing and need no action in this range.

## Inventory remarks

* B.tsv is accurate for what I checked. Corrections: row 48's `range != kTarget` text assumes `[0, L-1]` (see above, kTarget can be `[0, L]`; unreachable either way). The status "needs-param" for rows 16, 35, 42 became separate TBLASTX writers instead of parameters (`write_tblastx_alignment`), and `write_tblastn_hsp_info`/`write_sequence_row`/`write_blastx_alignment` stayed unchanged.
* B_notes.md section 4 asked for `Hit.num_positives` to be recomputed: it is (`pairwise_hits`, and `set_displayed_identities` for outfmt 6/7).

## "What the port must do" (B_notes section 5), where each item is in the final code

1. Every HSP of every kept subject in Seq-align order, header only before the first HSP: `write_tblastx_pairwise_report`, `report/pairwise.rs:2813-2866`.
2. Rows from the nucleotides with the CTrans_table rule, minus rows reverse complemented, minus query not flipped: `algorithm/tblastx/report.rs:56,104,137`; `report/pairwise.rs:2514`.
3. Identities, Positives, Gaps from the rows: `report.rs:182,423`.
4. `Frame = qf/sf` with `+`: `report/pairwise.rs:2422`.
5. `Score = ... Expect[(n)] = E`, `(n)` iff `sum_n >= 2`: `report/pairwise.rs:2422`.
6. Chunks of 60 residues = 180 nt, coordinates and padding: `report/pairwise.rs:2514`.
7. Middle line: `report/pairwise.rs:2514`.
8. Blank lines (chunk, HSP, footer): `report/pairwise.rs:2514,2863,2612`.
9. Lowercase of query residues covered by same-frame SEG masks with the minus-frame quirk: `algorithm/tblastx/report.rs:283,374`; `algorithm/blastx/query_setup.rs:517,557`.
10. Counts before lowercasing: `algorithm/tblastx/report.rs:423` (counts from `query_row` before `lowercase_query_row`).
