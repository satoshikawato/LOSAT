# Debug Session Summary 4

## Scope
- Investigate LOSATN megablast crash in greedy traceback (`reduce_gaps`) while keeping NCBI BLAST parity.

## Current Repro
- Command: `LOSAT/target/release/LOSAT blastn -q ./fasta/EDL933.fna -s ./fasta/Sakai.fna -o ./losat_out/EDL933.Sakai.losatn.megablast.out -n 1`
- Log: `LOSAT/tests/losat_out/EDL933.Sakai.losatn.megablast.debug.log`
- Failure: panic at `LOSAT/src/algorithm/blastn/alignment/greedy.rs:1431` (index out of bounds, len 7665, index 7676).

## Code Changes (this session)
- `LOSAT/src/algorithm/blastn/alignment/greedy.rs`
  - `reduce_gaps` uses signed indices (`isize`) to mirror NCBI pointer arithmetic.
  - Added NCBI reference comment: `ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:2669-2758`.
  - Crash still reproduces at the same site.

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

## Open Questions / Next Steps
1. Capture the *failing* `reduce_gaps` call state (edit script + lengths).
   - Run GDB interactively and print:
     - `q.length`, `s.length`, `esp.size`
     - `esp.num` and `esp.op_type` data via Vec pointer:
       - `set $num_ptr = (int*)esp.num.buf.inner.ptr.pointer.pointer`
       - `set $op_ptr = (int*)esp.op_type.buf.inner.ptr.pointer.pointer`
     - Sum lengths implied by the edit script:
       - Query consumed: `Sub + Ins`
       - Subject consumed: `Sub + Del`
   - Confirm whether the script implies more bases than `q.len()` / `s.len()`.
2. Verify greedy gapped extents against NCBI:
   - `q_end = q_off + q_ext_r`, `s_end = s_off + s_ext_r` should be within `q.len()` / `s.len()` for the slices used in `reduce_gaps`.
   - If not, trace `q_ext_*`/`s_ext_*` outputs from `blast_greedy_align` (non-affine) against NCBI `BLAST_GreedyAlign`.
3. Compare non-affine greedy traceback step-by-step with NCBI `greedy_align.c:379-751`.
   - Focus on the edit-script construction (`get_next_non_affine_tback` path) and `last_seq2_off[0][diag_origin]` usage.
4. Only after the failing edit-script mismatch is identified, adjust implementation (with NCBI refs) to restore parity and eliminate the panic.

