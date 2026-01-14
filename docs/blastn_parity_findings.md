# BLASTN parity findings (LOSAT vs NCBI BLAST+)

## Scope and sources
- LOSAT repo: `/mnt/c/Users/genom/GitHub/LOSAT`
- NCBI BLAST+ repo: `/mnt/c/Users/genom/GitHub/ncbi-blast/`
- Focus: blastn and megablast tasks, pipeline stages (setup, lookup, scanning, ungapped, gapped, traceback, hitlist pruning).

## Confirmed discrepancies that can affect output

### F1. Re-evaluation uses approximate search space
- LOSAT recalculates E-values after trimming using `eff_searchsp = db_len_total * queries[0].len`, which is only an approximation.
- NCBI uses context-aware effective lengths and `Blast_HSPListGetEvalues` for each HSP list.
- Evidence:
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2375`
  - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234`
- Impact: bit score and E-value mismatches for trimmed HSPs (post-purge and re-evaluation).

### F2. Interval tree bounds use max query length, not concatenated query length
- LOSAT sizes the interval tree using max query length, not the full concatenated query length across contexts.
- NCBI uses `query->length + 1`, which is the concatenated length for both strands.
- Evidence:
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:780`
  - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3811`
- Impact: containment checks can be wrong for multiple queries, affecting which HSPs are extended and kept.

### F3. Subject scanning and stride behavior diverges from NCBI
- LOSAT uses a custom rolling ASCII scan and applies `kmer_start % scan_step` to skip positions.
- NCBI uses compressed-subject scan paths with stride and mask handling in `blast_nascan.c`.
- Evidence:
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:902`
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1565`
  - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:147`
- Impact: different seed sets for megablast and for any non-default stride, especially with subject masking or ambiguous bases.

### F4. Hitlist size handling differs from BLAST+ CLI behavior
- LOSAT uses `min(hitlist_size, max_target_seqs)`.
- NCBI overrides hitlist size with `max_target_seqs` when explicitly set.
- Evidence:
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:356`
  - NCBI: `ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960`
- Impact: subject list pruning and ordering differences when user sets `-max_target_seqs`.

### F5. Subject-best-hit option not wired
- LOSAT includes `subject_best_hit` but does not expose or invoke it.
- NCBI runs `Blast_HSPListSubjectBestHit` when `-subject_besthit` is set.
- Evidence:
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/mod.rs:198`
  - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:586`
- Impact: option-level parity gap for users relying on subject-best-hit filtering.

### F6. Database word-count filtering and subject masking scan paths missing
- NCBI supports db word-count filtering and mask-aware scanning paths.
- No LOSAT equivalent found for `db_filter`, `max_db_word_count`, or subject `mask_type`.
- Evidence:
  - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:991`
  - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:178`
- Impact: option-level parity gap when those NCBI options are enabled.

## Proposed fix plan (priority order)
1. Replace approximate `eff_searchsp` in re-evaluation with NCBI-equivalent effective length search space per context.
   - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2375`
   - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:234`
2. Size interval tree with concatenated query length (`query->length + 1`) rather than max query length.
   - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:780`
   - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3811`
3. Align subject scanning to NCBI compressed-subject scan logic (stride, mask type, scan ranges).
   - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:902`
   - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:147`
4. Match BLAST+ hitlist size behavior when `-max_target_seqs` is set.
   - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:356`
   - NCBI: `ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960`
5. Wire `-subject_besthit` and db-word-count filtering options to match NCBI feature set.
   - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/mod.rs:198`
   - NCBI: `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:586`

## Deep-diff work still to do (gapped/greedy/ungapped)
These are the files to diff next for residual parity gaps:
- Ungapped extension: `LOSAT/src/algorithm/blastn/extension.rs` vs `ncbi-blast/c++/src/algo/blast/core/na_ungapped.c`
- Gapped DP: `LOSAT/src/algorithm/blastn/alignment/gapped.rs` vs `ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c`
- Greedy gapped: `LOSAT/src/algorithm/blastn/alignment/greedy.rs` vs `ncbi-blast/c++/src/algo/blast/core/greedy_align.c`

The next pass should verify boundary conditions, x-drop handling, packed subject offsets, and seed-selection logic are identical to NCBI.
