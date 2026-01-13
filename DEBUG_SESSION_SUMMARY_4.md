# Debug Session Summary 4

## Scope
- Investigate LOSATN megablast crash in greedy traceback (`reduce_gaps`) while keeping NCBI BLAST parity.

## Current Repro
- Command: `LOSAT/target/release/LOSAT blastn -q ./fasta/EDL933.fna -s ./fasta/Sakai.fna -o ./losat_out/EDL933.Sakai.losatn.megablast.out -n 1`
- Log: `LOSAT/tests/losat_out/EDL933.Sakai.losatn.megablast.debug.log`
- Prior failure: panic at `LOSAT/src/algorithm/blastn/alignment/greedy.rs:1431` (index out of bounds, len 7665, index 7676). Resolved: command now completes after `reduce_gaps` offset fix.

## Latest Result
- User reran the release megablast command from `LOSAT/tests`; it completed successfully (no panic; progress 1/1).

## Code Changes (this session)
- `LOSAT/src/algorithm/blastn/alignment/greedy.rs`
  - `reduce_gaps` uses signed indices (`isize`) to mirror NCBI pointer arithmetic.
  - `reduce_gaps` now indexes into full query/subject buffers using `q_start/q_end` and `s_start/s_end`, matching NCBI `s_ReduceGaps` pointer semantics (`blast_gapalign.c:2669-2758`, call at `blast_gapalign.c:2881`).
  - Megablast repro now completes (panic resolved).

## Existing Relevant Changes in Working Tree (from earlier progress)
- `LOSAT/src/core/blast_encoding.rs`
  - Added BLASTNA encoding table and `encode_iupac_to_blastna` / `encode_iupac_to_ncbi2na`.
  - NCBI refs: `blast_encoding.c:85-93`, `blast_util.c:476-489`.
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
  - Greedy traceback uses BLASTNA subject (uncompressed), while greedy score-only uses NCBI2NA.
  - NCBI refs: `blast_stat.h:866-869`, `blast_traceback.c:503-507`, `blast_util.c:806-833`.
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
  - Removed recursion in HSP insertion to prevent stack overflow.
  - NCBI ref: `blast_itree.c:703-746`.

## Debug Notes
- `reduce_gaps` is entered from `greedy_gapped_alignment_internal` with `gap_open=0`, `gap_extend=0` (non-affine greedy path, megablast defaults).
- GDB breakpoint at `reduce_gaps` fires, but the failing call state is not yet captured.
  - Initial breakpoint inspection showed `esp.size=1` for one call (not necessarily the crashing call).
  - Attempts to dump `esp.num` / `esp.op_type` with a scripted GDB loop failed due to Vec layout access and `max-value-size` limits.
- NCBI reference for `s_ReduceGaps` confirms it uses raw pointer arithmetic with `q`/`s` and end pointers (`qf`/`sf`), no bounds checks. Rust slice bounds are now the crash point, indicating an edit-script/extent mismatch.

## Commands Run
- `cargo build --release` (success, warnings only)
- `cargo build` (debug build, success)
- Megablast repro command above
- GDB breakpoints on `LOSAT/src/algorithm/blastn/alignment/greedy.rs:1412`

## Next Steps Toward NCBI Parity
1. Fix BLASTN coordinate off-by-one in `LOSAT/src/algorithm/blastn/alignment/gapped.rs` by matching NCBI `blast_gapalign.c:735-962` (`Blast_SemiGappedAlign`), especially `extend_gapped_one_direction_ex`.
2. Remove non-NCBI `filter_contained_hsps` and cleanup exports/imports (`LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`, `LOSAT/src/algorithm/blastn/filtering/mod.rs`, `LOSAT/src/algorithm/blastn/blast_engine/mod.rs`).
3. Re-run parity comparisons (EDL933 vs Sakai megablast; NZ_CP006932 self blastn) and use `LOSAT_DEBUG_COORDS`/`LOSAT_DEBUG_BLASTN` to localize remaining missing hits.
4. Investigate TBLASTX 2x hit inflation on long sequences (gencode=4) by matching `Blast_HSPListReevaluateUngapped` and extension boundary/termination logic.
