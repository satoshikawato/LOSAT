# BLASTN Performance Notes (LOSAT vs NCBI)

This note captures the main reasons LOSAT BLASTN (especially megablast,
single-thread) can be slower than NCBI BLASTN, with direct code references.
It is intended as a checklist for parity-safe performance work.

## Observed gap

- Example report: NCBI BLASTN < 1s vs LOSAT BLASTN ~6-7s on the same inputs.

## Likely hotspots (with code references)

1) Subject scan fallback walks every base even when `scan_step` is large.
   - In `scan_subject_kmers_range`, the fallback path runs `pos += 1` and only
     emits on `pos == next_emit`, so it still scans O(subject_len) bases for
     each subject when `lut_word_length > 8`.
   - Evidence: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:328`,
     `LOSAT/src/algorithm/blastn/blast_engine/run.rs:342`.
   - This path is hit for megablast two-stage lookup where `lut_word_length`
     is often 9-12 (see adaptive lookup selection).
   - NCBI uses specialized megablast scan loops that advance by `scan_step`
     using packed windows (no per-base loop):
     `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1482`.

2) Per-subject re-encoding of subject sequences.
   - Each subject is encoded to BLASTNA and packed ncbi2na inside the subject
     loop, adding extra full passes over the subject per scan.
   - Evidence: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1451`.
   - This is extra overhead vs NCBI’s pre-packed DB access (or cached
     sequences) in `BlastSeqSrc`.

3) Lookup table allocation and structure overhead.
   - Two-stage and PV direct lookups allocate a `Vec<Vec<...>>` of size 4^k
     (k up to 12), which is heavy to allocate and pointer-chase.
   - Evidence: `LOSAT/src/algorithm/blastn/lookup.rs:576`,
     `LOSAT/src/algorithm/blastn/lookup.rs:593`.
   - NCBI uses compact chain arrays and PV bitsets (see `BlastMBLookupTableNew`).
     Reference: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1227`.

4) Seed mask checking in the hot loop re-derives kmers and scans hit lists.
   - `is_seed_masked` calls `packed_kmer_at` and then searches the hit list for
     the exact `(q_idx, q_pos)` match.
   - Evidence: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2045`.
   - This is correct for parity, but may be expensive when many seeds are
     tested.

## Quick verification steps

- Run with `--verbose` to confirm adaptive lookup (`lut_word_length`,
  `scan_step`), which determines whether the fallback scan path is used.
  Evidence: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:935`.
- Set `LOSAT_DEBUG_BLASTN=1` and look for `[PERF] Subject scan took ...` to
  confirm scan time dominates.
  Evidence: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2428`.
- Ensure inputs are comparable: NCBI `-db` uses preformatted data; LOSAT reads
  FASTA for each run, which adds IO/encoding overhead.
  Evidence: `LOSAT/src/algorithm/blastn/coordination.rs:405`.

## Parity-safe improvement directions (NCBI-based)

- Port the NCBI megablast scan variants (`s_MBScanSubject_*`) so scan stride
  advances by `scan_step` and uses packed windows for `lut_word_length` 9-12.
  Reference: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nascan.c:1482`.
- Use NCBI-like lookup table layouts (contiguous chains + PV) to reduce memory
  and improve cache locality.
  Reference: `/mnt/c/Users/kawato/Documents/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1227`.
- Cache/retain subject encodings during `prepare_sequence_data` to avoid
  per-subject re-encoding in the hot path.
  Evidence: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1451`.
