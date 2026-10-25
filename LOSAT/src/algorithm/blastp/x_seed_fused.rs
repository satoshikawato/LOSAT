//! EXPERIMENT (LOSAT_X_BPSEED): the blastp two-hit seed stage with the subject scan and the
//! diagonal test in one loop.
//!
//! The reference path copies every word hit of the subject into the offset-pair buffer
//! (`s_blast_aa_scan_subject_one_range`, NCBI's `s_BlastAaScanSubject`) and then reads the buffer
//! back for the two-hit test (NCBI's `s_BlastAaWordFinder_TwoHit`). Here the test reads the hit
//! list of the lookup cell in place, in the same order (subject offset ascending, the cell's list
//! in its stored order), and applies the outcomes that only touch the hit's diagonal entry. A hit
//! that needs the rest of the test (the context check and the ungapped extension) goes to the
//! caller's closure, in scan order, with the values the reference computes at that point.
//!
//! Why the outcome is the reference's: the test of a hit reads and writes only the entry of its
//! diagonal `(q_off - s_off) & diag_mask`, and the closure of a candidate writes only that entry
//! and the HSP list. The hits are taken in the reference's order, so every entry sees the same
//! sequence of reads and writes and the HSP list gets the same appends in the same order.
//! The buffer's capacity (`array_size`) only cuts the reference's scan into rounds; it changes
//! neither the order nor the set of hits.
//!
//! `LOSAT_X_BPSEED` is a switch of the strict set (`utils/xswitch.rs`): on unless it is `0`.
//!
//! EXPERIMENT (LOSAT_X_BPDIAG): with the fused loop on, the hits of one subject offset (one lookup
//! cell) are tested in batches of eight instead of one after the other: for every hit of a batch
//! the entry's new value and whether the hit is a candidate are computed without branches in
//! AVX2 registers (eight lanes), every entry of the batch is stored, and then the candidates go to
//! the closure in the cell's order.
//! Why the outcome is the per-hit loop's:
//! - The query offsets of one cell are distinct: NCBI's `s_AddWordHits` adds a query offset to a
//!   cell once, as a neighbour of the offset's word or (when the word's self score is below the
//!   threshold, or the threshold is 0) as the word itself, which the neighbour enumeration then
//!   does not return. The diagonal table is at least the query's length plus the window
//!   (`BlastDiagTableNew`), so the hits of one subject offset fall on distinct entries.
//! - The test of a hit reads and writes only its entry, and the closure writes only the
//!   candidate's entry and the HSP list. So the hits of one subject offset are independent, and
//!   only the order of the closure's calls matters, which is the cell's order.
//! - A candidate's entry keeps its value in the test (the batch stores its old value back), so
//!   storing the batch before the calls changes nothing the closure reads, and a call's write is
//!   never overwritten by the batch.
//! - The values are the bit-field's (31-bit `last_hit`, sign-extended on read, truncated on
//!   write) and the i32 arithmetic of `x_test_hit`, which wraps (see there).
//!
//! `LOSAT_X_BPDIAG` is a switch of the strict set: on unless it is `0`.
//! It needs `LOSAT_X_BPSEED`. Without AVX2
//! (another architecture, or an x86-64 CPU without it) the switch keeps the per-hit loop: the
//! scalar batch was slower than the per-hit loop on both machines it was timed on (2026-10-14,
//! human block 1 x Species4, one thread, scan_ungapped: workstation 12.56 s against 11.59 s,
//! local 7.53 / 7.08 s against 6.78 / 6.88 s). The scalar batch remains as the AVX2 path's handler
//! of a cell's last batch of fewer than eight hits and as the cell-by-cell version the tests
//! compare with; `LOSAT_X_BPDIAG_SIMD=0` forces it (diagnostic, used by the tests; not in the
//! strict set). With `LOSAT_TIMING=1` each search of a query set prints
//! `[TIMING] x_bpdiag: mode=<avx2|avx2-defer|avx2-defer-masked|scalar|per-hit> cells=<n>
//! batched_hits=<m> candidates=<k>` (the path that ran, the cells tested in batches, their hits,
//! and the hits passed to the closure; zeros with `per-hit`).
//!
//! EXPERIMENT (LOSAT_X_BPDEFER): with the AVX2 batches on, a batch's candidates are not called at
//! the end of the batch: their query offsets are packed into a buffer (the candidate lanes moved to
//! the front of the register by a permutation from a table of the 256 lane masks, the register
//! stored, the buffer's length advanced by the number of candidates), and the closure is called for
//! the buffer in order after the cell's last batch of eight. The batch loop then has no branch on
//! the candidates (on the workstation 60 % of the cell test's mispredicts were the per-batch "no
//! candidate" branch, 2026-10-16). Why the outcome is the per-batch calls': the hits of one cell lie
//! on distinct entries and a call writes only its candidate's entry (above), so the batches after a
//! candidate neither read nor write what its call reads or writes, and the calls keep their order
//! (the batches' order, then the lanes'); a candidate's entry still holds the value its test read.
//! A cell's last batch of fewer than eight hits is the scalar batch, called after the buffer.
//! `LOSAT_X_BPDEFER` is a switch of the strict set: on unless it is `0`.
//! It acts only where the AVX2 batch runs (`LOSAT_X_BPDIAG` with AVX2, not
//! with `LOSAT_X_BPDIAG_SIMD=0`).
//!
//! EXPERIMENT (LOSAT_X_BXDEFER, blastx's switch): blastx's fused loop (`blastx/x_seed_diag.rs`)
//! takes LOSAT_X_BPDEFER's deferred calls with two changes of the batches
//! ([`x_cell_avx2_defer_masked`]): a batch of eight computes its lanes' entry coordinates from
//! scalar loads of its hits instead of taking them out of the hits' register, which starts the
//! batch's chain of entry loads earlier; and a cell's last batch of fewer than eight hits is one
//! masked AVX2 batch whose missing lanes repeat the first hit, its candidates packed into the same
//! buffer, instead of the scalar batch after the calls. A blastx cell holds about 19 hits, a blastp
//! cell about 106 (v4-F s3), so the last batch is a larger share of blastx's cell test. Why the
//! outcome is LOSAT_X_BPDEFER's: the coordinates are the same values; a repeated lane reads the
//! first hit's entry before any lane of the batch stores (every load of a batch comes before its
//! stores), so it stores the value lane 0 stores, and its candidate bit is dropped; the hits of the
//! last batch lie on entries that no other hit of the cell and no other call touches (above; the
//! closure is handed only its candidate's entry, as the table is borrowed by the cell test), so
//! testing them before the calls changes nothing a call reads, and the calls keep their order (the
//! batches', then the last batch's lanes'). The two changes come from line R's probes of blastp's
//! cell test (v4-R R2, 2026-10-18); blastp's fused loop takes them too since LOSAT_X_BPMASKED
//! (below), and LOSAT_X_BPDEFER's loops keep their code (the masked batch has a scan function of
//! its own, never inlined, and a cell test of its own). `LOSAT_X_BXDEFER` is a switch of
//! the strict set: on unless it is `0`; off, blastx's
//! fused loop takes the AVX2 batch ([`XDiagMode::Avx2`]). It acts only where that loop runs the AVX2
//! batch and LOSAT_X_BXCROSS's cell test does not (below: `LOSAT_X_BXCROSS=0`, or a CPU without
//! POPCNT), on tables with at least `LOSAT_X_BXDEFERMIN` hits per used cell
//! (`blastx/x_seed_diag.rs`).
//!
//! EXPERIMENT (LOSAT_X_BPMASKED): blastp's fused loop takes LOSAT_X_BXDEFER's masked batch
//! ([`x_cell_avx2_defer_masked`], [`XDiagMode::Avx2DeferMasked`]) in place of LOSAT_X_BPDEFER's
//! deferred calls: the same buffer and calls, with the batches' entry coordinates from scalar loads
//! of the hits and a cell's last batch of fewer than eight hits as one masked AVX2 batch. Why the
//! outcome is LOSAT_X_BPDEFER's: the argument of LOSAT_X_BXDEFER above, which rests only on what
//! both loops share (the hits of one cell lie on distinct entries; the closure is handed only its
//! candidate's entry, as the table is borrowed by the cell test). A blastp cell holds about 100
//! hits on the dense tables of all-against-all searches, where the masked batch saves 4-6 % of the
//! job's cycles (one core, against `LOSAT_X_BPMASKED=0`; the cell test 7.28 -> 6.32 cycles per
//! hit on an OrthoFinder job, v4-A14 s6). On a table with about one hit per used cell (a short
//! query against a database) the masked last batch and the buffer cost more than the scalar batch
//! they replace (x1.07-1.19 at 1.1-1.5 hits per used cell), so a table with fewer than
//! `LOSAT_X_BPMASKEDMIN` hits per used cell (hundredths, built-in 250; LOSAT_X_BXDIAGMIN's measure,
//! `blastx/x_seed_diag.rs`) keeps the deferred calls ([`x_bp_table_mode`], once per thread and
//! table).
//! `LOSAT_X_BPMASKED` is a switch of the strict set: on unless it is `0`; off, blastp's fused loop
//! takes LOSAT_X_BPDEFER's deferred calls ([`XDiagMode::Avx2Defer`]). It acts only where those
//! run (`LOSAT_X_BPDEFER` with the AVX2 batch, [`x_diag_mode`]; not with LOSAT_X_BPSEED=0,
//! LOSAT_X_BPDIAG=0 or LOSAT_X_BPDIAG_SIMD=0, not in the per-hit loop of a CPU without AVX2), and
//! never in blastp-fast's loops on the compressed table, which have no AVX2 batches.
//!
//! EXPERIMENT (LOSAT_X_BXCROSS, blastx's switch): blastx's fused loop tests
//! the hits in batches of eight cut across the cells' borders ([`x_scan_cells_avx2_cross`]). The
//! scan's call for a cell only stages its hits (per hit the entry coordinate `(q - s) & mask` and
//! the subject offset, in two arrays of one allocation: three masked loads of eight lanes copy up
//! to 24 hits, a loop the rest); when [`X_CROSS_CHUNK`] (1,024) hits are staged, and at the end of
//! the subject, the staged stream is cut into windows of eight from its first hit. A window whose
//! eight hits lie on eight distinct entries is one batch; otherwise (hits of different cells on one
//! diagonal: the hits of one cell lie on distinct entries, above) it is tested in parts, each part
//! the hits before the first one that repeats an earlier one's entry of the part, as a batch whose
//! other lanes repeat the part's first hit. A batch loads its entries, stores their new values (the
//! lanes of [`x_batch_avx2_at`] with a subject offset per lane) and calls the closure for its
//! candidates in lane order, all before the next batch loads. The fewer than eight hits left at the
//! end of the subject go one after the other through [`x_test_hit`], each candidate called at once;
//! those left when a chunk is tested move to the front of the stage, before the hits staged next.
//! The query offset is not staged: the table has a power of two entries, at least the query's
//! length plus the window (`blastx/seed.rs`, `Diagonals::new`, NCBI's `diag_array_length`), so
//! every query offset `q` of the lookup is at most `mask`, and `q = (coord + s) & mask` exactly
//! ([`x_cross_query`]; checked once per subject scan in [`x_scan_cells_cross`]). Why the outcome
//! is the per-hit loop's: the staged stream is the scan order, and a hit's test and a candidate's
//! call read and write only the hit's entry (the closure is handed only that entry; blastx's
//! closure writes nothing else of the table). In a batch of distinct entries no lane reads an entry
//! that an earlier lane's store or call writes, so loading every entry, then storing every entry,
//! then calling the candidates in lane order gives every entry the reads and writes, and every call
//! the arguments, of the hit-by-hit order; a candidate's entry keeps its value in the store, so its
//! call reads the value its test read. The repeated lanes of a part read their entry before any
//! lane stores and store, after the first lane, the value it stores, and their candidate bits are
//! dropped (as in LOSAT_X_BXDEFER's masked batch), so they change nothing. The last hits take the
//! hit-by-hit order itself, and every batch ends (stores and calls) before the next one loads.
//! `LOSAT_X_BXCROSS` is a switch of the strict set (since 2026-10-23): on unless it is `0`; on,
//! it is the cell test of every table blastx's fused loop takes where the CPU has AVX2 and POPCNT, in place of LOSAT_X_BXDEFER and its rule
//! (`blastx/x_seed_diag.rs`); off, or without POPCNT, blastx's fused loop takes the cell test it
//! took before (LOSAT_X_BXDEFER's choice). blastp's loops never take it: it is not a cell test of
//! [`XDiagMode`], the choice that blastp's loops share ([`x_scan_two_hit_in_with`]), but blastx's
//! [`XCellTest::Cross`], which blastx's [`x_scan_two_hit_with`] runs, and every function of it is
//! its own. The chunk is a constant of each instance; the program runs
//! the one of 1,024, and the tests also one of 64 (`X_CROSS_TEST_CHUNK`), which drains the stage
//! every few cells.
//!
//! EXPERIMENT (LOSAT_X_BXCROSSDEFER, blastx's switch): LOSAT_X_BXCROSS's drain defers the calls of
//! each unit (a batch, or a part of a window tested in parts) to the unit that next reads one of
//! their entries ([`x_cross_drain`] with `DEFER`). A unit's candidate lanes store the marker
//! ([`X_CROSS_MARKER`], raw `0xFFFF_FFFF`: flag 1, last_hit -1) in place of the entry they read
//! and store without it; the entries it read, its first staged hit and its candidate lanes go to
//! a record past the stage's two arrays (written by every unit without a branch; the count
//! advances by a unit with candidates). A unit that reads the marker in one of its entries stops
//! before it stores; the deferred calls are made in order ([`x_cross_flush`]: records in unit
//! order, lanes in lane order), each candidate's entry first given back the value its unit read,
//! and the unit runs again from the drain loop's head, so that none of the batch's values is live
//! across the calls (a call inside the loop made it spill them: v4-A14 s8 r1). The drain makes the
//! calls left after its last window, before the hits left are tested one by one or moved. Why the
//! outcome is that of the calls made at once: the calls keep their order; a call reads and writes
//! only its candidate's entry (above), which holds the marker from its unit to its call, so no unit
//! reads or stores it in between (a unit that would read it has the calls made first), and the call
//! sees the value its unit read and, at once, stored back; every other entry a unit reads holds
//! what it would hold with the calls made at once, as no deferred call writes it. A stored value
//! can equal the marker (a flagged last_hit of -1: a call stores `s_last_off - (word - 1) +
//! offset`, -1 in 31 bits near offsets of 2^31; any bits in the tests' tables): a unit that meets
//! it with no call deferred runs as without the deferral (no marker, its calls at once; a unit
//! that meets it with calls deferred makes them first, which changes nothing), so any table
//! content gives the same calls and the same table. On blastx's dense tables a later unit reads a
//! deferred candidate's entry 0.017-0.029 times per batch, and the deferral cuts the drain's
//! branch mispredictions by 0.08-0.15 per batch for 3-5 more instructions (v4-A14 s8).
//! `LOSAT_X_BXCROSSDEFER` is a switch of the strict set (since 2026-10-25): on unless it is `0`;
//! off, each unit's candidates are called before the next unit loads (the instance without
//! `DEFER`, whose functions are those before the switch). It acts only where LOSAT_X_BXCROSS's
//! cell test runs, and is a field of blastx's [`XCellTest::Cross`], so blastp's loops keep their
//! code.
//!
//! BXSCANBF (graduated 2026-10-26: blastx always takes it): `x_scan_two_hit_with(.., branch_free: true, ..)`
//! finds the cells with hits by the candidate loop of tblastn's branch-free scan
//! (`tblastn/x_scan_bf.rs`, `x_block_candidates`): the positions of a block of words are written
//! to a list whose length advances by the presence bit, without a branch on it, and the listed
//! cells then go to the cell test in position order. The cells, their order and their hit lists
//! are the presence-test loop's, so every outcome above is unchanged. blastp's own loop
//! (`x_scan_two_hit`) takes the same cells under its own switch, LOSAT_X_BPSCANBF (a switch of the
//! strict set; [`x_bp_scan_bf`]); blastp-fast's loop on the flat compressed table lists its cells
//! the same way under LOSAT_X_BPCLUTSCANBF (`x_seed_fused_clut.rs`).
//!
//! EXPERIMENT (LOSAT_X_BPHITBF): the per-hit loop (`x_cell_per_hit`: the fused loop without
//! LOSAT_X_BPDIAG's batches, the default where the CPU has no AVX2, i.e. aarch64 and WebAssembly)
//! takes the two outcomes most two-hit tests end in, too far and overlap, without a branch (tblastn's
//! tblastn's branch-free two-hit loop and TBLASTX's LOSAT_X_TBXHITBF do the same in their loops, on tables whose entries
//! need no sign extension):
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-551
//! ```c
//!                 if (diff >= window) {
//!                     /* We are beyond the window for this diagonal; start a
//!                        new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//!                 }
//!
//!                 /* If the difference is less than the wordsize (i.e. last
//!                    hit and this hit overlap), give up */
//!
//!                 if (diff < wordsize) {
//!                     continue;
//!                 }
//! ```
//! With the per-hit loop on every subject (`LOSAT_X_BPDIAG=0`) 72-74 % of the hits of the two blastp
//! cases of v4-A10 probe bp1 are too far and 18.6-18.8 % overlap, which of the two is hard to
//! predict, and 2.5-3.7 % meet a flag-1 entry and 4.7-5.3 % are candidates for the extension.
//!
//! `XHitBf::test` reads the packed entry `raw` (`last_hit` in bits 0-30, the flag in bit 31) and
//! works with twice the values: `d2 = 2 (s + diag_offset) - (raw << 1)`, in wrapping `u32`
//! arithmetic read as `i32`. For a flag-0 entry `raw << 1` is twice the 31-bit field's signed value
//! modulo 2^32 (a cleared entry holds `-window`), so `d2` is twice the reference's `diff = s -
//! (last_hit - diag_offset)` whenever that fits in an `i32`. `XHitBf::new` (once per subject) makes
//! it fit: with `0 <= diag_offset` and `diag_offset + length + 2 window < 2^30` (`length` the
//! subject's letters), every flag-0 `last_hit` the table holds is in `[-window, diag_offset +
//! length)`: the cleared value (`DiagStruct::clear`, when the offset passes `INT4_MAX / 4`), 0 in a
//! new or refilled table, or an `s + diag_offset` of this subject or of an earlier one. Every store
//! with flag 0 is such an `s + diag_offset` (the test's, and the extension's for a hit outside its
//! context or one that does not extend to the right; the extension's `s_last_off - (wordsize - 1) +
//! diag_offset` is stored with flag 1, and the flag goes back to 0 only with a store of `s +
//! diag_offset`). The offset only grows between clears (`BlastpSubjectScratch::finish_subject`; a
//! subject searched on its own starts at the serial run's offset on a refilled table,
//! `prepare_independent_subject`), so an earlier subject's stores are below this subject's offset,
//! and none of them wrapped in the 31-bit field, as they were below 2^30. `s + diag_offset` is in
//! `[diag_offset, diag_offset + length)`, so `diff` is in `(-length, diag_offset + length + window)`
//! and `2 diff` in `(-2^31, 2^31)`. Then too far (`diff >= window`) is `d2 >= 2 window`, and a
//! candidate (`wordsize <= diff < window`) is `(d2 - 2 wordsize) as u32 < 2 (window - wordsize)`,
//! one unsigned compare: below `2 wordsize` the difference is negative (it does not wrap, as `2
//! length + 2 wordsize < 2^31`), so at least 2^31 as a `u32`. The rule also asks `1 <= wordsize <
//! window`, a two-hit search. A flag-1 entry (`raw as i32 < 0`) or a candidate goes on to
//! [`x_test_hit`], unchanged; any other hit stores `s + diag_offset` (too far; the flag stays 0) or
//! what it read (overlap), selected without a branch. So the table holds after each hit what
//! `x_test_hit` leaves there, the extension is called for the same hits with the same values and
//! entries, and the initial HSPs are the same, in the same order. blastx's direct loop takes the same
//! test (`XHitBf` is shared; LOSAT_X_BXHITBF, `blastx/x_seed_scan.rs`, whose header shows that its
//! table holds the same range of flag-0 values), and so does blastp-fast's loop on the flat
//! compressed table (LOSAT_X_BPFSEED, `x_seed_fused_clut.rs`: the same diagonal table and
//! extension as blastp's loop).
//!
//! The switch is read once; it is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_BPHITBF=0`. The loop with the test is its own function, never inlined and with a
//! closure of its own type (`x_scan_two_hit_hitbf`), so the loops without it keep their code; the
//! AVX2 and scalar batches keep their own tests, and so does blastp-fast's loop on the reference
//! compressed table (`x_seed_fused_compressed.rs`). blastp-fast's loop on the flat compressed table
//! takes the test, without the read-ahead, under its own switch (LOSAT_X_BPFSEED; LOSAT_X_BPHITBF
//! does not act there). `x_scan_two_hit_hitbf` also reads ahead without a bound check per hit
//! (`x_cell_per_hit_bf`). With `LOSAT_X_STATS` (and the `xstats` feature) `main` prints
//! `[X_STATS] BPHITBF subjects=N` at exit when N is not zero, N the subjects the loop with the
//! test scanned.

use crate::algorithm::tblastn::x_scan_bf::{x_block_candidates, X_BLOCK};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::{
    BlastAaLookupTable, AA_HITS_PER_CELL, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use crate::utils::matrix::BLASTAA_SIZE;
use std::mem::MaybeUninit;

// No NCBI counterpart: reads the LOSAT_X_BPSEED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSEED): a strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_seed() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPSEED"))
}

// No NCBI counterpart: how many hits ahead the fused loop reads the diagonal entries
// (LOSAT_X_BPSEED_PREFETCH, default 16; 0 = no prefetch). Tuning only.
pub(super) fn x_prefetch_distance() -> usize {
    use std::sync::OnceLock;
    static D: OnceLock<usize> = OnceLock::new();
    *D.get_or_init(|| {
        std::env::var("LOSAT_X_BPSEED_PREFETCH")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(16)
    })
}

// No NCBI counterpart: a cache hint for one diagonal entry; it reads nothing and changes no value.
#[inline(always)]
pub(super) fn x_prefetch(diag_array: &[DiagStruct], index: usize) {
    #[cfg(target_arch = "x86_64")]
    if let Some(entry) = diag_array.get(index) {
        // SAFETY: a prefetch of a valid address has no effect on program state.
        unsafe {
            std::arch::x86_64::_mm_prefetch(
                (entry as *const DiagStruct).cast::<i8>(),
                std::arch::x86_64::_MM_HINT_T0,
            );
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = (diag_array, index);
}

// No NCBI counterpart: reads the LOSAT_X_BPDIAG switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDIAG): a strict-set switch, on unless it is `0`.
/// It acts in the fused loop only.
pub(crate) fn x_bp_diag() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDIAG"))
}

// No NCBI counterpart: reads the LOSAT_X_BPDEFER switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPDEFER): a strict-set switch, on unless it is `0`.
/// It acts on the AVX2 batches only.
pub(crate) fn x_bp_defer() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDEFER"))
}

// No NCBI counterpart: reads the LOSAT_X_BPMASKED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPMASKED): a strict-set switch, on unless it is `0`.
/// It acts on LOSAT_X_BPDEFER's deferred calls only.
pub(crate) fn x_bp_masked() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPMASKED"))
}

// No NCBI counterpart: reads LOSAT_X_BPMASKEDMIN once; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): the fewest hits per used lookup cell, in hundredths, for
/// which blastp's fused loop takes LOSAT_X_BPMASKED's masked batch (built-in
/// [`X_BP_MASKED_MIN`]); a table with fewer takes LOSAT_X_BPDEFER's deferred calls. With about one
/// hit per cell the masked last batch and the buffer cost more than the scalar batch they replace
/// (one core, cycles against `LOSAT_X_BPMASKED=0` with the masked batch on every table: 1.11 hits
/// per used cell x1.148, 2.08 x1.020, 2.60 x0.987, 4.52 x0.924; `docs/agents/diagnostics.md`).
/// `0` takes it for every table; a value that is not a number takes the built-in one.
#[cfg(target_arch = "x86_64")]
fn x_bp_masked_min() -> u64 {
    use std::sync::OnceLock;
    static MIN: OnceLock<u64> = OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BPMASKEDMIN")
            .ok()
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(X_BP_MASKED_MIN)
    })
}

/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): the built-in value of [`x_bp_masked_min`], between the
/// measured loss at 2.08 hits per used cell and the gain at 2.60.
#[cfg(target_arch = "x86_64")]
const X_BP_MASKED_MIN: u64 = 250;

#[cfg(target_arch = "x86_64")]
thread_local! {
    /// EXPERIMENT (LOSAT_X_BPMASKEDMIN): the last lookup table this thread chose a cell test for
    /// (its address, its backbone's address and its overflow length) and its hits per used cell
    /// (see [`x_bp_table_ratio`]).
    static X_BP_CHOICE: std::cell::Cell<(usize, usize, usize, u64)> =
        const { std::cell::Cell::new((0, 0, 0, 0)) };
}

// No NCBI counterpart: counters of LOSAT_X_BPMASKEDMIN for its timing line: the (thread, table)
// choices of the masked batch and of the deferred calls, and the least and most hits per used
// cell (hundredths) of the tables.
#[cfg(target_arch = "x86_64")]
static X_BP_TABLES: [std::sync::atomic::AtomicU64; 2] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];
#[cfg(target_arch = "x86_64")]
static X_BP_RATIO_MIN: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(u64::MAX);
#[cfg(target_arch = "x86_64")]
static X_BP_RATIO_MAX: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

// No NCBI counterpart: a property of the lookup table that chooses LOSAT's cell test; it does
// not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): `table`'s hits per used backbone cell, in hundredths (0 for
/// a table without hits; LOSAT_X_BXDIAGMIN's measure, `blastx/x_seed_diag.rs`), worked out once
/// per table and thread (the table's addresses and overflow length identify it; a later table
/// that matched all three would reuse the value, which changes only the speed), and whether this
/// call worked it out.
#[cfg(target_arch = "x86_64")]
fn x_bp_table_ratio(table: &BlastAaLookupTable) -> (u64, bool) {
    let key = (
        table as *const BlastAaLookupTable as usize,
        table.backbone.as_ptr() as usize,
        table.overflow.len(),
    );
    X_BP_CHOICE.with(|choice| {
        let (a, b, c, ratio) = choice.get();
        if (a, b, c) == key {
            return (ratio, false);
        }
        let (mut hits, mut used) = (0u64, 0u64);
        for cell in &table.backbone {
            if cell.num_used > 0 {
                hits += cell.num_used as u64;
                used += 1;
            }
        }
        let ratio = (hits * 100).checked_div(used).unwrap_or(0);
        choice.set((key.0, key.1, key.2, ratio));
        (ratio, true)
    })
}

// No NCBI counterpart: chooses LOSAT's cell test; every cell test gives the reference's hits and
// table.
/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): the cell test for the subjects of `table`: `mode`, but
/// LOSAT_X_BPDEFER's deferred calls instead of LOSAT_X_BPMASKED's masked batch where the table has
/// fewer hits per used cell than [`x_bp_masked_min`]. Counted for the timing line and marked for
/// `LOSAT_X_REACH` once per table and thread.
fn x_bp_table_mode(mode: XDiagMode, table: &BlastAaLookupTable) -> XDiagMode {
    #[cfg(target_arch = "x86_64")]
    if mode == XDiagMode::Avx2DeferMasked {
        let (ratio, first) = x_bp_table_ratio(table);
        let masked = ratio >= x_bp_masked_min();
        if first {
            use std::sync::atomic::Ordering::Relaxed;
            X_BP_TABLES[usize::from(!masked)].fetch_add(1, Relaxed);
            X_BP_RATIO_MIN.fetch_min(ratio, Relaxed);
            X_BP_RATIO_MAX.fetch_max(ratio, Relaxed);
            crate::utils::xreach::x_reach!("LOSAT_X_BPMASKED", masked);
        }
        if !masked {
            return XDiagMode::Avx2Defer;
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = table;
    mode
}

// No NCBI counterpart: reads the LOSAT_X_BPSCANBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSCANBF): blastp's fused loop finds the cells with [`x_scan_cells_bf`]
/// (module header). A strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_scan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPSCANBF"))
}

// No NCBI counterpart: reads the LOSAT_X_BPHITBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPHITBF): whether blastp's per-hit loop takes the too-far and overlap
/// outcomes without a branch where [`XHitBf::new`] finds the subject in range (module header); a
/// switch of the strict set, on unless it is `0`.
pub(crate) fn x_bp_hit_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPHITBF"))
}

/// EXPERIMENT (LOSAT_X_BPHITBF, LOSAT_X_BXHITBF): the constants of the branch-free two-hit test
/// ([`XHitBf::test`]) for one subject: the diagonal offset, `2 wordsize`, `2 (window - wordsize)`
/// and `2 window`. blastp's per-hit loop, blastx's direct loop (`blastx/x_seed_scan.rs`) and
/// blastp-fast's loop on the flat compressed table (LOSAT_X_BPFSEED, `x_seed_fused_clut.rs`) use it.
#[derive(Clone, Copy, Debug)]
pub(crate) struct XHitBf {
    diag_offset: i32,
    word2: i32,
    span2: u32,
    window2: i32,
}

impl XHitBf {
    // No NCBI counterpart: the range in which LOSAT's branch-free test gives the reference's
    // outcomes; it does not change any value NCBI computes.
    /// Whether the test gives the reference's outcomes for a subject of `length` letters scanned
    /// at `diag_offset` (module header): `1 <= word_size < window`, `0 <= diag_offset` and
    /// `diag_offset + length + 2 window < 2^30`.
    pub(crate) fn fits(diag_offset: i32, length: usize, window: i32, word_size: i32) -> bool {
        let room = (1i64 << 30) - i64::from(diag_offset) - 2 * i64::from(window);
        1 <= word_size
            && word_size < window
            && 0 <= diag_offset
            && i64::try_from(length).is_ok_and(|n| n < room)
    }

    /// The test for a subject of `length` letters scanned at `diag_offset` where
    /// [`XHitBf::fits`] holds; `None` elsewhere.
    pub(crate) fn new(
        diag_offset: i32,
        length: usize,
        window: i32,
        word_size: i32,
    ) -> Option<Self> {
        Self::fits(diag_offset, length, window, word_size)
            .then(|| Self::in_range(diag_offset, window, word_size))
    }

    /// The constants for a caller that checked [`XHitBf::fits`] for the subject: there `2 window
    /// < 2^30`, so they do not overflow. The test is exact only there.
    #[inline(always)]
    pub(crate) fn in_range(diag_offset: i32, window: i32, word_size: i32) -> Self {
        Self {
            diag_offset,
            word2: word_size.wrapping_mul(2),
            span2: window.wrapping_sub(word_size).wrapping_mul(2) as u32,
            window2: window.wrapping_mul(2),
        }
    }

    /// The subject offset plus the diagonal offset, the `here` of [`XHitBf::test`].
    #[inline(always)]
    pub(crate) fn here(self, subject_offset: i32) -> u32 {
        subject_offset.wrapping_add(self.diag_offset) as u32
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:519-551 (quoted at `x_test_hit`)
    /// The two-hit test of one hit on `entry`, with `here` = [`XHitBf::here`]: when the entry's
    /// flag is 0 and the hit is too far or overlaps, stores what [`x_test_hit`] stores and returns
    /// `true`; otherwise (a flag-1 entry or a candidate) leaves the entry and returns `false`.
    /// Exact where [`XHitBf::fits`] holds for the subject (module header).
    #[inline(always)]
    pub(crate) fn test(self, entry: &mut DiagStruct, here: u32) -> bool {
        let raw = entry.raw_bits();
        let d2 = here.wrapping_mul(2).wrapping_sub(raw << 1) as i32;
        let candidate = (d2.wrapping_sub(self.word2) as u32) < self.span2;
        // `|`, not `||`: both tests are taken without a branch.
        if ((raw as i32) < 0) | candidate {
            return false;
        }
        entry.set_raw_bits(if d2 >= self.window2 { here } else { raw });
        true
    }
}

/// How the fused loop tests the hits of one lookup cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XDiagMode {
    /// One hit after the other ([`x_test_hit`]): the fused loop without LOSAT_X_BPDIAG.
    PerHit,
    /// LOSAT_X_BPDIAG with LOSAT_X_BPDIAG_SIMD=0 (diagnostic): batches of [`X_DIAG_LANES`] hits in
    /// scalar code ([`x_test_batch_scalar`]).
    Scalar,
    /// LOSAT_X_BPDIAG with AVX2 and LOSAT_X_BPDEFER off (and blastx's fused loop with
    /// LOSAT_X_BXDEFER off and LOSAT_X_BXCROSS off or no POPCNT): batches of eight hits in AVX2 registers ([`x_test_batch_avx2`]); a
    /// cell's last batch of fewer than eight hits is a scalar batch.
    #[cfg(target_arch = "x86_64")]
    Avx2,
    /// LOSAT_X_BPDIAG with AVX2 and LOSAT_X_BPDEFER, with LOSAT_X_BPMASKED off or on a table below
    /// LOSAT_X_BPMASKEDMIN ([`x_bp_table_mode`]): the AVX2 batches with the candidates of a cell's
    /// batches of eight called after its last one ([`x_cell_avx2_defer`]).
    #[cfg(target_arch = "x86_64")]
    Avx2Defer,
    /// LOSAT_X_BPDIAG with AVX2, LOSAT_X_BPDEFER and LOSAT_X_BPMASKED, on tables with at least
    /// LOSAT_X_BPMASKEDMIN hits per used cell (and blastx's fused loop with LOSAT_X_BXDEFER,
    /// [`x_cell_test`], and LOSAT_X_BXCROSS off or no POPCNT): [`XDiagMode::Avx2Defer`] with
    /// the lanes' entry coordinates from scalar loads of the hits and a cell's last batch of fewer
    /// than eight hits as one masked AVX2 batch ([`x_cell_avx2_defer_masked`]).
    #[cfg(target_arch = "x86_64")]
    Avx2DeferMasked,
}

impl XDiagMode {
    /// The name in the `[TIMING] x_bpdiag` line (and blastx's `[TIMING] blastx_diag`).
    pub(crate) fn name(self) -> &'static str {
        match self {
            XDiagMode::PerHit => "per-hit",
            XDiagMode::Scalar => "scalar",
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2 => "avx2",
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2Defer => "avx2-defer",
            #[cfg(target_arch = "x86_64")]
            XDiagMode::Avx2DeferMasked => "avx2-defer-masked",
        }
    }
}

/// EXPERIMENT (LOSAT_X_BXDIAG, LOSAT_X_BXCROSS): the cell test of blastx's fused seed loop
/// ([`x_cell_test`]): one of the cell tests of [`XDiagMode`], which blastx's loop runs through the
/// choice that blastp's loops share ([`x_scan_two_hit_in_with`]), or LOSAT_X_BXCROSS's, which
/// only blastx's [`x_scan_two_hit_with`] runs. A type of its own, so that the shared choice has no
/// arm for LOSAT_X_BXCROSS and blastp's loops keep their code.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XCellTest {
    /// A cell test of the shared choice.
    Diag(XDiagMode),
    /// LOSAT_X_BXCROSS: the subject's hits staged and tested in AVX2 batches of eight cut across
    /// the cells' borders ([`x_scan_cells_avx2_cross`]); with `defer` (LOSAT_X_BXCROSSDEFER) each
    /// batch's calls deferred to a later batch that meets their entries, or the drain's end.
    #[cfg(target_arch = "x86_64")]
    Cross { defer: bool },
}

impl XCellTest {
    /// The name in blastx's `[TIMING] blastx_diag` line.
    pub(crate) fn name(self) -> &'static str {
        match self {
            XCellTest::Diag(mode) => mode.name(),
            #[cfg(target_arch = "x86_64")]
            XCellTest::Cross { defer: true } => "avx2-cross-defer",
            #[cfg(target_arch = "x86_64")]
            XCellTest::Cross { defer: false } => "avx2-cross",
        }
    }
}

// No NCBI counterpart: whether the CPU has AVX2; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
fn x_avx2_available() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

// No NCBI counterpart: whether the CPU has what LOSAT_X_BXCROSS's functions are compiled for; it
// does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSS): AVX2 and POPCNT (the cross functions count a batch's candidates
/// with `popcnt` instead of the multiply sequence of a build for x86-64 without it).
#[cfg(target_arch = "x86_64")]
fn x_cross_available() -> bool {
    x_avx2_available() && std::arch::is_x86_feature_detected!("popcnt")
}

// No NCBI counterpart: chooses LOSAT's code path once (LOSAT_X_BPDIAG; the AVX2 batch where the
// CPU has AVX2, with LOSAT_X_BPDEFER the deferred calls, with LOSAT_X_BPMASKED also the masked
// batch, otherwise the per-hit loop; LOSAT_X_BPDIAG_SIMD=0 forces the scalar batch); it does not
// change any value NCBI computes.
fn x_diag_mode() -> XDiagMode {
    use std::sync::OnceLock;
    static MODE: OnceLock<XDiagMode> = OnceLock::new();
    *MODE.get_or_init(|| {
        if !x_bp_diag() {
            return XDiagMode::PerHit;
        }
        if std::env::var_os("LOSAT_X_BPDIAG_SIMD").is_some_and(|v| v.as_os_str() == "0") {
            return XDiagMode::Scalar;
        }
        #[cfg(target_arch = "x86_64")]
        if x_avx2_available() {
            return match (x_bp_defer(), x_bp_masked()) {
                (true, true) => XDiagMode::Avx2DeferMasked,
                (true, false) => XDiagMode::Avx2Defer,
                (false, _) => XDiagMode::Avx2,
            };
        }
        // Without AVX2 the scalar batch would be slower than the per-hit loop (module header).
        XDiagMode::PerHit
    })
}

// No NCBI counterpart: chooses LOSAT's code path once per call site (LOSAT_X_BXDIAG,
// LOSAT_X_BXDEFER, LOSAT_X_BXCROSS, LOSAT_X_BXCROSSDEFER); it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): the cell test of blastx's fused seed loop
/// (`blastx/x_seed_diag.rs`): `"per-hit"` or `"scalar"` when asked by name (diagnostics), otherwise
/// the AVX2 batch where the CPU has AVX2, else the per-hit loop, as LOSAT_X_BPDIAG chooses for
/// blastp (the scalar batch was slower than the per-hit loop; module header). The AVX2 batch is
/// [`XCellTest::Cross`] with `cross` (LOSAT_X_BXCROSS) where the CPU also has POPCNT (its calls
/// deferred with `defer`, LOSAT_X_BXCROSSDEFER), otherwise [`XDiagMode::Avx2DeferMasked`] with
/// `masked` (LOSAT_X_BXDEFER), [`XDiagMode::Avx2`] without.
pub(crate) fn x_cell_test(requested: &str, masked: bool, cross: bool, defer: bool) -> XCellTest {
    XCellTest::Diag(match requested {
        "per-hit" => XDiagMode::PerHit,
        "scalar" => XDiagMode::Scalar,
        _ => {
            #[cfg(target_arch = "x86_64")]
            if x_avx2_available() {
                if cross && x_cross_available() {
                    return XCellTest::Cross { defer };
                }
                return XCellTest::Diag(if masked {
                    XDiagMode::Avx2DeferMasked
                } else {
                    XDiagMode::Avx2
                });
            }
            #[cfg(not(target_arch = "x86_64"))]
            let _ = (masked, cross, defer);
            XDiagMode::PerHit
        }
    })
}

// No NCBI counterpart: counters of LOSAT_X_BPDIAG for its timing line.
static X_DIAG_CELLS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_DIAG_HITS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static X_DIAG_CANDIDATES: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BPDIAG): the counts of one subject's scan (per cell, not per hit), added to
/// the totals once.
#[derive(Default)]
struct XDiagCounts {
    cells: u64,
    hits: u64,
    candidates: u64,
}

impl XDiagCounts {
    #[inline(always)]
    fn add_cell(&mut self, hits: usize, candidates: u32) {
        self.cells += 1;
        self.hits += hits as u64;
        self.candidates += u64::from(candidates);
    }

    /// EXPERIMENT (LOSAT_X_BXCROSS): a cell and its hits as they are staged (their candidates are
    /// counted when the batches run, [`XDiagCounts::add_candidates`]).
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    fn add_staged(&mut self, hits: usize) {
        self.cells += 1;
        self.hits += hits as u64;
    }

    /// EXPERIMENT (LOSAT_X_BXCROSS): the candidates of the staged hits tested in one drain.
    #[cfg(target_arch = "x86_64")]
    #[inline(always)]
    fn add_candidates(&mut self, candidates: u32) {
        self.candidates += u64::from(candidates);
    }

    fn add_to_totals(&self) {
        if self.cells != 0 {
            use std::sync::atomic::Ordering::Relaxed;
            X_DIAG_CELLS.fetch_add(self.cells, Relaxed);
            X_DIAG_HITS.fetch_add(self.hits, Relaxed);
            X_DIAG_CANDIDATES.fetch_add(self.candidates, Relaxed);
        }
    }
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the counts of one subject's scan with the batches across the
/// cells' borders: [`XDiagCounts`]'s (the cells and hits as they are staged, the candidates per
/// drain) and the windows of eight staged hits with a repeated entry, tested in parts
/// ([`x_cross_window`], added once per drain). blastx adds them to the totals of its timing line
/// once per subject (`blastx/x_seed_diag.rs`), as it does the other cell tests' counts.
#[cfg(target_arch = "x86_64")]
#[derive(Default)]
struct XCrossCounts {
    diag: XDiagCounts,
    splits: u64,
}

// No NCBI counterpart: a LOSAT timing line.
/// EXPERIMENT (LOSAT_X_BPDIAG): with the switch on, the timing line of the cells tested in batches
/// since the last one.
pub(crate) fn x_bp_diag_print_timing() {
    if x_bp_diag() {
        use std::sync::atomic::Ordering::Relaxed;
        eprintln!(
            "[TIMING] x_bpdiag: mode={} cells={} batched_hits={} candidates={}",
            x_diag_mode().name(),
            X_DIAG_CELLS.swap(0, Relaxed),
            X_DIAG_HITS.swap(0, Relaxed),
            X_DIAG_CANDIDATES.swap(0, Relaxed)
        );
        // EXPERIMENT (LOSAT_X_BPMASKEDMIN): the tables that took the masked batch and those its
        // rule gave the deferred calls (once per thread and table).
        #[cfg(target_arch = "x86_64")]
        if x_diag_mode() == XDiagMode::Avx2DeferMasked {
            let least = X_BP_RATIO_MIN.swap(u64::MAX, Relaxed);
            eprintln!(
                "[TIMING] x_bpmasked: min_hits_per_cell={} tables_masked={} tables_deferred={} hits_per_cell={}-{}",
                x_bp_masked_min(),
                X_BP_TABLES[0].swap(0, Relaxed),
                X_BP_TABLES[1].swap(0, Relaxed),
                if least == u64::MAX { 0 } else { least },
                X_BP_RATIO_MAX.swap(0, Relaxed)
            );
        }
    }
}

// No NCBI counterpart: NCBI indexes `diag_array` directly; this keeps the table with its mask so
// that the loops below read entries without a bounds check on a range checked once.
/// The diagonal table of the fused loop and its mask; `new` checks that every
/// `(q_off - s_off) & mask` names an entry.
struct XDiagTable<'a> {
    entries: &'a mut [DiagStruct],
    mask: i32,
}

impl<'a> XDiagTable<'a> {
    fn new(entries: &'a mut [DiagStruct], mask: i32) -> Self {
        assert!(mask >= 0 && entries.len() > mask as usize);
        Self { entries, mask }
    }

    /// The entry index of a hit.
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:516
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    // ```
    #[inline(always)]
    fn coord(&self, query_offset: i32, subject_offset: i32) -> usize {
        (query_offset.wrapping_sub(subject_offset) & self.mask) as usize
    }

    #[inline(always)]
    fn entry(&mut self, query_offset: i32, subject_offset: i32) -> &mut DiagStruct {
        let coord = self.coord(query_offset, subject_offset);
        // SAFETY: coord <= mask < entries.len() (checked in `new`).
        unsafe { self.entries.get_unchecked_mut(coord) }
    }

    #[inline(always)]
    fn prefetch(&self, query_offset: i32, subject_offset: i32) {
        x_prefetch(self.entries, self.coord(query_offset, subject_offset));
    }

    /// EXPERIMENT (LOSAT_X_BPHITBF): [`XDiagTable::entry`] with the index taken as a `u32` (it is
    /// at most `mask`, so not negative), for the loop of [`x_cell_per_hit_bf`].
    #[inline(always)]
    fn entry_bf(&mut self, query_offset: i32, subject_offset: i32) -> &mut DiagStruct {
        let coord = (query_offset.wrapping_sub(subject_offset) & self.mask) as u32 as usize;
        // SAFETY: coord <= mask < entries.len() (checked in `new`).
        unsafe { self.entries.get_unchecked_mut(coord) }
    }

    // No NCBI counterpart: a cache hint for one diagonal entry; it reads nothing and changes no value.
    /// EXPERIMENT (LOSAT_X_BPHITBF): [`XDiagTable::prefetch`] without the bound check of
    /// [`x_prefetch`] (the index is at most `mask`), for the loop of [`x_cell_per_hit_bf`].
    #[inline(always)]
    fn prefetch_bf(&self, query_offset: i32, subject_offset: i32) {
        #[cfg(target_arch = "x86_64")]
        {
            let coord = (query_offset.wrapping_sub(subject_offset) & self.mask) as u32 as usize;
            // SAFETY: coord <= mask < entries.len() (checked in `new`), so the address is inside
            // the table; a prefetch has no effect on program state.
            unsafe {
                std::arch::x86_64::_mm_prefetch(
                    self.entries.as_ptr().add(coord).cast::<i8>(),
                    std::arch::x86_64::_MM_HINT_T0,
                );
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        let _ = (query_offset, subject_offset);
    }
}

/// The values of the diagonal test that stay the same for a subject.
#[derive(Clone, Copy, Debug)]
struct XDiagTest {
    diag_offset: i32,
    window: i32,
    word_size: i32,
    /// LOSAT_X_BPSEED_PREFETCH (the per-hit loop only).
    prefetch: usize,
}

/// EXPERIMENT (LOSAT_X_BPHITBF): of the subjects, those the per-hit loop scanned with the
/// branch-free test (counted with the `xstats` feature, for the line at exit).
// No NCBI counterpart: a counter; it does not change any value NCBI computes.
static X_HITBF_SUBJECTS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BPHITBF): with LOSAT_X_STATS, the `[X_STATS] BPHITBF subjects=N` line at
/// exit (only when N is not zero).
// No NCBI counterpart: prints a counter; it does not change any value NCBI computes.
pub fn x_bp_hit_bf_print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_none() {
        return;
    }
    let subjects = X_HITBF_SUBJECTS.load(std::sync::atomic::Ordering::Relaxed);
    if subjects > 0 {
        eprintln!("[X_STATS] BPHITBF subjects={subjects}");
    }
}

/// The diagonal test of one word hit as NCBI applies it before an extension, on the hit's own
/// entry. Returns `Some((last_hit, diff))` when the hit needs the rest of the test.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
//             /* If the reset bit is set, an extension just happened. */
//             if (diag_array[diag_coord].flag) {
//                 /* If we've already extended past this hit, skip it. */
//                 if ((Int4) (subject_offset + diag_offset) <
//                     diag_array[diag_coord].last_hit) {
//                     continue;
//                 }
//                 /* Otherwise, start a new hit. */
//                 else {
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     diag_array[diag_coord].flag = 0;
//                 }
//             }
//             /* If the reset bit is cleared, try to start an extension. */
//             else {
//                 /* find the distance to the last hit on this diagonal */
//                 last_hit = diag_array[diag_coord].last_hit - diag_offset;
//                 diff = subject_offset - last_hit;
//
//                 if (diff >= window) {
//                     /* We are beyond the window for this diagonal; start a
//                        new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
//
//                 /* If the difference is less than the wordsize (i.e. last
//                    hit and this hit overlap), give up */
//
//                 if (diff < wordsize) {
//                     continue;
//                 }
// ```
//
// The sums and differences wrap (two's complement i32): that is what `+` and `-` do in release
// builds (no overflow checks), and what NCBI's arithmetic on its `Uint4 subject_offset` gives.
// LOSAT_X_BPDIAG's batches compute the same wrapping values per lane, and their tests compare
// with this function on offsets whose sum passes 2^31 - 1.
#[inline(always)]
pub(super) fn x_test_hit(
    entry: &mut DiagStruct,
    subject_offset: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
) -> Option<(i32, i32)> {
    if entry.flag() != 0 {
        let subject_plus_offset = subject_offset.wrapping_add(diag_offset);
        if subject_plus_offset >= entry.last_hit() {
            entry.set_last_hit(subject_plus_offset);
            entry.set_flag(0);
        }
        return None;
    }
    let last_hit = entry.last_hit().wrapping_sub(diag_offset);
    let diff = subject_offset.wrapping_sub(last_hit);
    if diff >= window {
        entry.set_last_hit(subject_offset.wrapping_add(diag_offset));
        return None;
    }
    if diff < word_size {
        return None;
    }
    Some((last_hit, diff))
}

/// Scans `subject` with `lookup` (two-hit mode) and applies the diagonal test to every word hit
/// in scan order; `candidate(query_offset, subject_offset, last_hit, diff, entry)` receives each
/// hit that needs the rest of the test, with its diagonal entry.
///
/// `subject` starts at the first residue and keeps the sentinel after the last one, as the
/// reference scan reads it; `subject_length` is the number of residues.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:495-512
// ```c
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
//     scan_range[2] = subject->seq_ranges[0].right - wordsize;
//
//     if (scan_range[2] < scan_range[1])
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
//
//         totalhits += hits;
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ```
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit<F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let test = XDiagTest {
        diag_offset,
        window,
        word_size,
        prefetch: x_prefetch_distance(),
    };
    let table = &mut XDiagTable::new(diag_array, diag_mask);
    // EXPERIMENT (LOSAT_X_BPHITBF): the per-hit loop with the branch-free two-hit test, where it
    // gives the reference's outcomes (module header).
    if let Some(bf) = x_hit_bf_for(
        x_diag_mode() == XDiagMode::PerHit,
        diag_offset,
        subject_length,
        window,
        word_size,
    ) {
        x_scan_two_hit_hitbf(
            x_bp_scan_bf(),
            lookup,
            subject,
            subject_length,
            table,
            test,
            bf,
            &mut candidate,
        );
        return;
    }
    // EXPERIMENT (LOSAT_X_BPMASKEDMIN): the masked batch only on tables with enough hits per cell.
    let mode = x_bp_table_mode(x_diag_mode(), lookup);
    // EXPERIMENT (LOSAT_X_BPSCANBF): the same cells found without a branch on the presence bit.
    let counts = if x_bp_scan_bf() {
        x_scan_two_hit_bf(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    } else {
        x_scan_two_hit_in(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    };
    counts.add_to_totals();
}

/// EXPERIMENT (LOSAT_X_BPSCANBF): [`x_scan_two_hit_in`] with the cells found by [`x_scan_cells_bf`].
/// Out of line, and with a closure of its own type, so that the default instance in
/// [`x_scan_two_hit`] and the cell tests it calls are compiled as without the switch: with the
/// same closure type the cell test had two callers, and the default instance ran 2-6 % more
/// instructions.
#[inline(never)]
fn x_scan_two_hit_bf<F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut own = |query_offset, subject_offset, last_hit, diff, entry: &mut DiagStruct| {
        candidate(query_offset, subject_offset, last_hit, diff, entry)
    };
    x_scan_two_hit_in_with::<true, _>(mode, lookup, subject, subject_length, table, test, &mut own)
}

// No NCBI counterpart: chooses LOSAT's loop; both loops give the reference's outcomes.
/// EXPERIMENT (LOSAT_X_BPHITBF): the branch-free test for a subject when the switch is on, the
/// loop is a per-hit loop (`per_hit`) and [`XHitBf::new`] finds the subject in range; the subject
/// is counted for the line of LOSAT_X_STATS.
fn x_hit_bf_for(
    per_hit: bool,
    diag_offset: i32,
    subject_length: usize,
    window: i32,
    word_size: i32,
) -> Option<XHitBf> {
    if !(per_hit && x_bp_hit_bf()) {
        return None;
    }
    let bf = XHitBf::new(diag_offset, subject_length, window, word_size)?;
    crate::utils::xstats::add(&X_HITBF_SUBJECTS, 1);
    Some(bf)
}

/// EXPERIMENT (LOSAT_X_BPHITBF): the per-hit loop of [`x_scan_two_hit`] with the branch-free
/// two-hit test ([`x_cell_per_hit_bf`]), the cells found by [`x_scan_cells_bf`] when
/// `branch_free` (LOSAT_X_BPSCANBF). Out of line, and with a closure of its own type, so that the
/// loops without the test keep their code (see [`x_scan_two_hit_bf`]).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_two_hit_hitbf<F>(
    branch_free: bool,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut own = |query_offset, subject_offset, last_hit, diff, entry: &mut DiagStruct| {
        candidate(query_offset, subject_offset, last_hit, diff, entry)
    };
    if branch_free {
        x_scan_cells_with::<true, _>(lookup, subject, subject_length, |s, hits| {
            x_cell_per_hit_bf(hits, s, table, test, bf, &mut own)
        })
    } else {
        x_scan_cells_with::<false, _>(lookup, subject, subject_length, |s, hits| {
            x_cell_per_hit_bf(hits, s, table, test, bf, &mut own)
        })
    }
}

/// EXPERIMENT (LOSAT_X_BXDIAG): [`x_scan_two_hit`] with the cell test `mode` (blastx's choice,
/// [`x_cell_test`], not LOSAT_X_BPDIAG's); returns the cells tested in batches, their hits, the
/// hits passed to `candidate` (zeros for [`XDiagMode::PerHit`]) and the windows of eight staged
/// hits that LOSAT_X_BXCROSS's batches tested in parts (zero for the other cell tests). Nothing is
/// added to the totals of the `[TIMING] x_bpdiag` line. Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx):
/// `branch_free` finds the cells with [`x_scan_cells_bf`] (module header). EXPERIMENT
/// (LOSAT_X_BXCROSS): [`XCellTest::Cross`] runs here ([`x_scan_cells_cross`]); the other cell
/// tests run through the choice that blastp's loops share ([`x_scan_two_hit_in_with`]).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:495-512 (quoted at `x_scan_two_hit`)
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit_with<F>(
    mode: XCellTest,
    branch_free: bool,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    diag_array: &mut [DiagStruct],
    diag_mask: i32,
    diag_offset: i32,
    window: i32,
    word_size: i32,
    mut candidate: F,
) -> [u64; 4]
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let test = XDiagTest {
        diag_offset,
        window,
        word_size,
        prefetch: x_prefetch_distance(),
    };
    let table = &mut XDiagTable::new(diag_array, diag_mask);
    let mode = match mode {
        XCellTest::Diag(mode) => mode,
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { defer } => {
            let counts = x_scan_cells_cross::<X_CROSS_CHUNK, F>(
                branch_free,
                defer,
                lookup,
                subject,
                subject_length,
                table,
                test,
                &mut candidate,
            );
            let XDiagCounts {
                cells,
                hits,
                candidates,
            } = counts.diag;
            return [cells, hits, candidates, counts.splits];
        }
    };
    // A constant per instantiation, so that each scan loop has one call of the cell test and the
    // compiler inlines it, as in blastp's loop.
    let counts = if branch_free {
        x_scan_two_hit_in_with::<true, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    } else {
        x_scan_two_hit_in_with::<false, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    };
    [counts.cells, counts.hits, counts.candidates, 0]
}

/// [`x_scan_two_hit`] with the cell test of `mode`; returns the counts of LOSAT_X_BPDIAG (zero
/// for [`XDiagMode::PerHit`]).
fn x_scan_two_hit_in<F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_scan_two_hit_in_with::<false, F>(
        mode,
        lookup,
        subject,
        subject_length,
        table,
        test,
        candidate,
    )
}

/// [`x_scan_two_hit_in`], with the cells found by [`x_scan_cells_bf`] when `BRANCH_FREE`
/// (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx).
fn x_scan_two_hit_in_with<const BRANCH_FREE: bool, F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XDiagCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut counts = XDiagCounts::default();
    match mode {
        XDiagMode::PerHit => {
            x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
                x_cell_per_hit(hits, s, table, test, candidate)
            })
        }
        XDiagMode::Scalar => {
            x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
                let candidates = x_cell_scalar(hits, s, table, test, candidate);
                counts.add_cell(hits.len(), candidates);
            })
        }
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2 => {
            assert!(x_avx2_available(), "the AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                x_scan_cells_avx2::<BRANCH_FREE, F>(
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                    &mut counts,
                )
            }
        }
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2Defer => {
            assert!(x_avx2_available(), "the deferred AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                x_scan_cells_avx2_defer::<BRANCH_FREE, F>(
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                    &mut counts,
                )
            }
        }
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2DeferMasked => {
            assert!(x_avx2_available(), "the masked AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                x_scan_cells_avx2_defer_masked::<BRANCH_FREE, F>(
                    lookup,
                    subject,
                    subject_length,
                    table,
                    test,
                    candidate,
                    &mut counts,
                )
            }
        }
    }
    counts
}

/// EXPERIMENT (LOSAT_X_BXCROSS): [`x_scan_cells_avx2_cross`] with the chunk `CHUNK` (the program
/// runs [`X_CROSS_CHUNK`]; the tests also `X_CROSS_TEST_CHUNK`) and the scan of `branch_free`
/// (BXSCANBF, graduated for blastx), after the checks that the CPU has AVX2 and POPCNT and of the bound the
/// stage rests on (below); returns the subject's counts. Out of line, and called by blastx's
/// [`x_scan_two_hit_with`] before the choice of the cell test that blastp's loops share
/// ([`x_scan_two_hit_in_with`]), so that blastp's loops keep their code.
///
/// The stage keeps no query offset ([`x_cross_query`] recovers it as `(coord + s) & mask`), which
/// is exact for the query offsets `0..=mask`. blastx's fused loop takes `mask` from its diagonal
/// table (`blastx/x_seed_diag.rs`, `x_word_finder_fused`: `diagonals.entries.len() - 1`), made by
/// `Diagonals::new(last.offset + last.length, window)` (`blastx/search.rs`; `blastx/seed.rs`:
/// the next power of two at least the query's length plus the window, NCBI's
/// `diag_array_length`), so `mask >= last.offset + last.length`. The lookup is built on the same
/// concatenated query (`blastx/seed.rs`, `Lookup::new`: `batch.sequence_start[1..]`) from the
/// segments of `complement_masks` (`blastx/query_setup.rs`), which lie inside the contexts, so
/// its offsets are word starts below `lookup.query_length = last.offset + last.length`
/// (`tblastx/lookup/backbone.rs`, `build_lookup_from_prepared`). The assert below checks that
/// bound once per subject scan (one compare; none per hit).
// No NCBI counterpart: chooses LOSAT's cell test; it does not change any value NCBI computes.
#[cfg(target_arch = "x86_64")]
#[inline(never)]
fn x_scan_cells_cross<const CHUNK: usize, F>(
    branch_free: bool,
    defer: bool,
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> XCrossCounts
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    assert!(
        x_cross_available(),
        "the cross-cell AVX2 batch needs AVX2 and POPCNT"
    );
    assert!(
        (0..=i64::from(table.mask) + 1).contains(&i64::from(lookup.query_length)),
        "the cross-cell stage needs every query offset at most the diagonal mask"
    );
    // `mask + 1` a power of two (written in i64, so that `mask + 1` cannot overflow): the
    // recovery `(coord + s) & mask` of `x_cross_query` equals `q` only then (mask 5, q 2, s 3
    // stages coord (2 - 3) & 5 = 5 and recovers (5 + 3) & 5 = 0).
    assert!(
        i64::from(table.mask) >= 0 && i64::from(table.mask) & (i64::from(table.mask) + 1) == 0,
        "the cross-cell stage needs a diagonal mask of the form 2^k - 1"
    );
    let mut counts = XCrossCounts::default();
    // SAFETY: the CPU has AVX2 and POPCNT, and every query offset of the lookup is below
    // `query_length`, so at most the mask (asserted above).
    unsafe {
        match (branch_free, defer) {
            (true, true) => x_scan_cells_avx2_cross::<true, true, CHUNK, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
            (true, false) => x_scan_cells_avx2_cross::<true, false, CHUNK, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
            (false, true) => x_scan_cells_avx2_cross::<false, true, CHUNK, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
            (false, false) => x_scan_cells_avx2_cross::<false, false, CHUNK, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
        }
    }
    counts
}

/// [`x_scan_cells`] with [`x_cell_avx2`], compiled for AVX2 as a whole so that the cell test is
/// inlined into the scan.
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_scan_cells_avx2<const BRANCH_FREE: bool, F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates = unsafe { x_cell_avx2(hits, s, table, test, candidate) };
        counts.add_cell(hits.len(), candidates);
    });
}

thread_local! {
    /// EXPERIMENT (LOSAT_X_BPDEFER): the candidate buffer of a thread's scans, kept between them.
    static X_DEFER_BUFFER: std::cell::Cell<Vec<i32>> = const { std::cell::Cell::new(Vec::new()) };
}

/// [`x_scan_cells_avx2`] with [`x_cell_avx2_defer`] (LOSAT_X_BPDEFER).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_scan_cells_avx2_defer<const BRANCH_FREE: bool, F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut buffer = X_DEFER_BUFFER.take();
    x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates = unsafe { x_cell_avx2_defer(hits, s, table, test, candidate, &mut buffer) };
        counts.add_cell(hits.len(), candidates);
    });
    X_DEFER_BUFFER.set(buffer);
}

/// [`x_scan_cells_avx2`] with [`x_cell_avx2_defer_masked`] (LOSAT_X_BPMASKED in blastp's fused
/// loop, LOSAT_X_BXDEFER in blastx's). Never inlined, and its cell test is a function of its own,
/// so the loops of LOSAT_X_BPDEFER (`LOSAT_X_BPMASKED=0`) keep their code.
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn x_scan_cells_avx2_defer_masked<const BRANCH_FREE: bool, F>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut buffer = X_DEFER_BUFFER.take();
    x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates =
            unsafe { x_cell_avx2_defer_masked(hits, s, table, test, candidate, &mut buffer) };
        counts.add_cell(hits.len(), candidates);
    });
    X_DEFER_BUFFER.set(buffer);
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the number of staged hits at which the batches run (a constant
/// of each instance, so a cell's stores and the chunk test take no register); the stage's two
/// arrays then take 8 KiB and stay in the first cache level.
#[cfg(target_arch = "x86_64")]
const X_CROSS_CHUNK: usize = 1024;

/// EXPERIMENT (LOSAT_X_BXCROSS, tests only): the chunk of the second instance the tests run,
/// which drains the stage every few cells (and inside long cells), so that the hits left by a
/// drain and the windows with a repeated entry occur across drains on small inputs.
#[cfg(all(test, target_arch = "x86_64"))]
const X_CROSS_TEST_CHUNK: usize = 64;

/// EXPERIMENT (LOSAT_X_BXCROSS): the hits a piece of a cell's copy stages with fixed vector work
/// (three registers of eight lanes).
#[cfg(target_arch = "x86_64")]
const X_CROSS_PIECE: usize = 3 * X_DIAG_LANES;

/// EXPERIMENT (LOSAT_X_BXCROSS): the room in each array past the chunk: a piece starts with fewer
/// than the chunk staged and stores [`X_CROSS_PIECE`] lanes (those past its hits are never read as hits),
/// so a drain starts with fewer than the chunk plus 24 staged, and its move of the fewer than eight
/// left reads eight lanes from below that count ([`x_cross_drain`]).
#[cfg(target_arch = "x86_64")]
const X_CROSS_SLACK: usize = 32;

#[cfg(target_arch = "x86_64")]
const _: () = assert!(X_CROSS_SLACK >= X_CROSS_PIECE + X_DIAG_LANES - 1);

#[cfg(target_arch = "x86_64")]
thread_local! {
    /// EXPERIMENT (LOSAT_X_BXCROSS): the stage of a thread's scans (the two arrays of [`XCross`]
    /// in one allocation), kept between them.
    static X_CROSS_STAGE: std::cell::Cell<Vec<i32>> = const { std::cell::Cell::new(Vec::new()) };
}

// No NCBI counterpart: NCBI keeps the offset pair of every hit in `offset_pairs`; the stage keeps
// the entry coordinate and the subject offset of the hits not yet tested.
/// EXPERIMENT (LOSAT_X_BXCROSS): one subject's view of the stage: two arrays of
/// [`XCross::STRIDE`] values in one allocation, so that a cell's stores take one base register
/// and constant displacements; per staged hit its entry coordinate `(q - s) & mask` (the array at
/// `base`) and its subject offset (at `base + STRIDE`); its query offset is recovered from the two
/// ([`x_cross_query`]). The number of hits staged is the scan's (a local, passed to and returned by
/// [`x_cross_cell`] and [`x_cross_drain`]); between two cells it is below `CHUNK`.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct XCross<'a, const CHUNK: usize> {
    base: *mut i32,
    _stage: std::marker::PhantomData<&'a mut Vec<i32>>,
}

#[cfg(target_arch = "x86_64")]
impl<'a, const CHUNK: usize> XCross<'a, CHUNK> {
    /// The values of each array: the chunk and the slack.
    const STRIDE: usize = CHUNK + X_CROSS_SLACK;

    /// EXPERIMENT (LOSAT_X_BXCROSSDEFER): the records of deferred calls a drain can hold: one per
    /// unit (a batch or a part of a window) with a candidate, and a drain tests fewer than
    /// `CHUNK + X_CROSS_PIECE` hits, each unit at least one.
    const PENDING: usize = CHUNK + X_CROSS_PIECE;

    /// The values of the stage with `DEFER` (LOSAT_X_BXCROSSDEFER): the two arrays, then per record
    /// the eight entries its unit read, then per record its first staged hit and candidate lanes;
    /// without, the two arrays.
    const fn len(defer: bool) -> usize {
        if defer {
            2 * Self::STRIDE + (X_DIAG_LANES + 2) * Self::PENDING
        } else {
            2 * Self::STRIDE
        }
    }

    fn new<const DEFER: bool>(stage: &'a mut Vec<i32>) -> Self {
        if stage.len() < Self::len(DEFER) {
            stage.resize(Self::len(DEFER), 0);
        }
        Self {
            base: stage.as_mut_ptr(),
            _stage: std::marker::PhantomData,
        }
    }

    /// The entry coordinates.
    #[inline(always)]
    fn coord(self) -> *mut i32 {
        self.base
    }

    /// The subject offsets.
    #[inline(always)]
    fn subject(self) -> *mut i32 {
        // SAFETY: STRIDE < 2 STRIDE, the allocation's length (`new`).
        unsafe { self.base.add(Self::STRIDE) }
    }

    /// EXPERIMENT (LOSAT_X_BXCROSSDEFER): the eight entries a deferred record's unit read (record
    /// `r` at `8 r`). Only in a stage made with `DEFER`.
    #[inline(always)]
    fn pending_old(self) -> *mut i32 {
        // SAFETY: 2 STRIDE < len(true), the allocation's length (`new::<true>`).
        unsafe { self.base.add(2 * Self::STRIDE) }
    }

    /// EXPERIMENT (LOSAT_X_BXCROSSDEFER): a deferred record's first staged hit and candidate lanes
    /// (record `r` at `2 r`). Only in a stage made with `DEFER`.
    #[inline(always)]
    fn pending_at(self) -> *mut i32 {
        // SAFETY: 2 STRIDE + 8 PENDING < len(true), the allocation's length (`new::<true>`).
        unsafe {
            self.base
                .add(2 * Self::STRIDE + X_DIAG_LANES * Self::PENDING)
        }
    }
}

// No NCBI counterpart: NCBI keeps the query offset of every hit (`offset_pairs`); the stage keeps
// the entry coordinate and the subject offset and recovers the query offset from them.
/// EXPERIMENT (LOSAT_X_BXCROSS): the query offset of a staged hit from its entry coordinate
/// `coord = (q - s) & mask` and its subject offset `s`: `(coord + s) & mask` is `q & mask` (the
/// sums wrap modulo 2^32, a multiple of the table's power-of-two length), which is `q` for the
/// query offsets `0..=mask` of the lookup ([`x_scan_cells_cross`] checks the bound and that the
/// table's length is a power of two once per subject; [`x_cross_cell`] checks each cell's offsets in
/// debug builds).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn x_cross_query(coord: i32, s: i32, mask: i32) -> i32 {
    coord.wrapping_add(s) & mask
}

/// [`x_scan_cells_avx2`] with LOSAT_X_BXCROSS's batches across the cells' borders (blastx's fused
/// loop only; module header): `on_cell` stages the cell's hits ([`x_cross_cell`], which runs the
/// batches when `CHUNK` hits are staged), and the hits staged at the end of the subject are tested
/// by [`x_cross_drain`]. The number of hits staged is kept in a local during the scan. Never
/// inlined, and every function of its cell test is its own, so blastp's loops keep their code.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, and every query offset of `lookup` must be at most
/// `table.mask` ([`x_scan_cells_cross`]).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline(never)]
unsafe fn x_scan_cells_avx2_cross<
    const BRANCH_FREE: bool,
    const DEFER: bool,
    const CHUNK: usize,
    F,
>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XCrossCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut stage = X_CROSS_STAGE.take();
    {
        let cross = XCross::<CHUNK>::new::<DEFER>(&mut stage);
        let mut staged = 0usize;
        x_scan_cells_with::<BRANCH_FREE, _>(lookup, subject, subject_length, |s, hits| {
            // SAFETY: the CPU has AVX2 and POPCNT (this function's contract); `staged` is 0 at the
            // start and below the chunk after every call of `x_cross_cell`.
            staged = unsafe {
                x_cross_cell::<DEFER, CHUNK, F>(
                    cross, staged, hits, s, table, test, candidate, counts,
                )
            };
        });
        x_cross_drain::<DEFER, CHUNK, F>(cross, staged, true, table, test, candidate, counts);
    }
    X_CROSS_STAGE.set(stage);
}

/// Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): [`x_scan_cells`], or [`x_scan_cells_bf`] when `BRANCH_FREE` (the
/// same calls). One loop per instantiation, so `on_cell` has one call site and is inlined.
#[inline(always)]
fn x_scan_cells_with<const BRANCH_FREE: bool, C>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    on_cell: C,
) where
    C: FnMut(i32, &[i32]),
{
    if BRANCH_FREE {
        x_scan_cells_bf(lookup, subject, subject_length, on_cell);
    } else {
        x_scan_cells(lookup, subject, subject_length, on_cell);
    }
}

/// The subject scan of [`x_scan_two_hit`]: `on_cell(subject_offset, hits)` for every scan position
/// whose lookup cell has hits, in scan order, with the cell's list in its stored order.
// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99
// ```c
//    for (s = s_first; s <= s_last; s++) {
//       /* compute the index value */
//       index = ComputeTableIndexIncremental(word_length, lut->charsize,
//                                            lut->mask, s, index);
//       /* if there are hits... */
//       if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
//          Int4 numhits = bbcell[index].num_used;
//          ...
//          if (numhits <= AA_HITS_PER_CELL)
//             /* hits live in thick_backbone */
//             src = bbcell[index].payload.entries;
//          else
//             /* hits live in overflow array */
//             src = & (lut->overflow[bbcell[index].payload.overflow_cursor]);
//
//          /* copy the hits. */
//          for (i = 0; i < numhits; i++) {
//             offset_pairs[i + totalhits].qs_offsets.q_off = src[i];
//             offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
//          }
// ```
#[inline(always)]
fn x_scan_cells<C>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    mut on_cell: C,
) where
    C: FnMut(i32, &[i32]),
{
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    let s_first = 0i32;
    let mut s_last =
        i32::try_from(subject_length).expect("subject length must fit Int4") - lookup.word_length;
    if s_last < s_first {
        s_last = s_first;
    }
    // A subject shorter than a word has one scan position whose word runs past the sentinel; the
    // reference reads that byte unchecked (one byte past the buffer for a one-residue subject).
    // It is read as a sentinel here. With one scan position no two hits share a diagonal, and the
    // next subject starts window + length further on the diagonals, so the byte reaches no output.
    let byte_at = |i: usize| subject.get(i).copied().unwrap_or(0) as usize;
    let pv = &lookup.pv;
    let backbone = &lookup.backbone;
    let overflow = &lookup.overflow;

    let mut index = 0usize;
    for i in 0..word_length - 1 {
        index = (index << charsize) | byte_at(i);
    }
    for s in s_first..=s_last {
        let new_char = byte_at(s as usize + word_length - 1);
        index = ((index << charsize) | new_char) & mask;
        if pv[index >> PV_ARRAY_BTS] & (1u32 << (index & PV_ARRAY_MASK)) == 0 {
            continue;
        }
        let cell = &backbone[index];
        let numhits = cell.num_used as usize;
        let hits: &[i32] = if numhits <= AA_HITS_PER_CELL {
            &cell.entries[..numhits]
        } else {
            let cursor = cell.entries[0] as usize;
            &overflow[cursor..cursor + numhits]
        };
        on_cell(s, hits);
    }
}

/// Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): [`x_scan_cells`] with the cells that have hits found without a
/// branch on the presence bit: for each block of [`X_BLOCK`] positions [`x_block_candidates`] lists
/// `position << 32 | index`, the list's length advanced by the bit, and the listed cells then go to
/// `on_cell` in position order. The positions, indices and presence bits are those of
/// [`x_scan_cells`]'s loop (the same rolling index over the same bytes), so the calls are the same.
// NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99 (quoted at `x_scan_cells`)
#[inline(always)]
fn x_scan_cells_bf<C>(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    subject_length: usize,
    mut on_cell: C,
) where
    C: FnMut(i32, &[i32]),
{
    let word_length = lookup.word_length as usize;
    let charsize = lookup.charsize as usize;
    let mask = lookup.mask as usize;
    i32::try_from(subject_length).expect("subject length must fit Int4");
    // The bytes x_scan_cells reads: positions 0..=s_last, s_last = max(subject_length -
    // word_length, 0), each with its word. A subject shorter than a word has one position whose
    // word runs past the subject, and x_scan_cells reads the missing bytes as 0: here they come
    // from a zero-padded copy.
    let mut padded = [0u8; 8];
    let (bytes, s_last) = if subject_length >= word_length {
        (&subject[..subject_length], subject_length - word_length)
    } else {
        assert!(word_length <= padded.len(), "word size {word_length}");
        let n = subject.len().min(word_length);
        padded[..n].copy_from_slice(&subject[..n]);
        (&padded[..word_length], 0)
    };
    // The backbone has NCBI's `backbone_size` cells (the index of a word of letters below
    // BLASTAA_SIZE, plus one) and the presence vector a bit for each, so the unchecked presence
    // reads of `x_block_candidates` stay inside it when every byte read is below BLASTAA_SIZE
    // (x_scan_cells would panic on an index past the vector instead).
    let largest = bytes.iter().copied().max().unwrap_or(0);
    assert!(
        usize::from(largest) < BLASTAA_SIZE
            && lookup.pv.len() << PV_ARRAY_BTS >= lookup.backbone.len(),
        "letter {largest} outside the lookup table's alphabet"
    );
    let backbone = &lookup.backbone;
    let overflow = &lookup.overflow;

    let mut index = 0usize;
    for &byte in &bytes[..word_length - 1] {
        index = (index << charsize) | byte as usize;
    }
    let mut cand = [MaybeUninit::<u64>::uninit(); X_BLOCK];
    let mut block_first = 0usize;
    while block_first <= s_last {
        let block_last = s_last.min(block_first + X_BLOCK - 1);
        // Reads bytes[block_last + word_length - 1] at most, the last byte of `bytes`.
        let n = x_block_candidates(
            bytes,
            lookup.pv.as_ptr(),
            (word_length, charsize, mask),
            (block_first, block_last),
            &mut index,
            &mut cand,
        );
        for c in &cand[..n] {
            // SAFETY: the first `n` candidates were written by `x_block_candidates`.
            let c = unsafe { c.assume_init() };
            let cell = &backbone[(c & 0xffff_ffff) as usize];
            let numhits = cell.num_used as usize;
            let hits: &[i32] = if numhits <= AA_HITS_PER_CELL {
                &cell.entries[..numhits]
            } else {
                let cursor = cell.entries[0] as usize;
                &overflow[cursor..cursor + numhits]
            };
            on_cell((c >> 32) as i32, hits);
        }
        block_first = block_last + 1;
    }
}

/// The diagonal test of a cell's hits one after the other (the fused loop without
/// LOSAT_X_BPDIAG).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[inline(always)]
fn x_cell_per_hit<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    for (i, &query_offset) in hits.iter().enumerate() {
        // Read ahead the entry of a later hit of this cell: the hits of one subject offset
        // fall on distinct diagonals spread over the table, so each is a cache miss that the
        // hardware cannot predict. A prefetch reads nothing the test depends on.
        if let Some(&ahead) = hits.get(i + test.prefetch) {
            table.prefetch(ahead, s);
        }
        let entry = table.entry(query_offset, s);
        if let Some((last_hit, diff)) =
            x_test_hit(entry, s, test.diag_offset, test.window, test.word_size)
        {
            candidate(query_offset, s, last_hit, diff, entry);
        }
    }
}

/// EXPERIMENT (LOSAT_X_BPHITBF): [`x_cell_per_hit`] with each hit first through the branch-free
/// test ([`XHitBf::test`]): a flag-1 entry or a candidate goes on to [`x_test_hit`], unchanged
/// (module header).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[inline(always)]
fn x_cell_per_hit_bf<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let here = bf.here(s);
    // The read-ahead of `x_cell_per_hit` (a cache hint, x86-64 only): the hits with a later hit
    // `prefetch` places ahead read its entry ahead, in a loop of their own so that no hit checks
    // whether there is one; the last `prefetch` hits read none. (With `prefetch` 0 no hit reads
    // ahead; `x_cell_per_hit` then hints the hit's own entry, which it reads at once.)
    let p = test.prefetch;
    let ahead = if p == 0 {
        0
    } else {
        hits.len().saturating_sub(p)
    };
    let (first, rest) = hits.split_at(ahead);
    for (&query_offset, &later) in first.iter().zip(&hits[p.min(hits.len())..]) {
        table.prefetch_bf(later, s);
        x_hit_bf(query_offset, s, here, table, test, bf, candidate);
    }
    for &query_offset in rest {
        x_hit_bf(query_offset, s, here, table, test, bf, candidate);
    }
}

/// EXPERIMENT (LOSAT_X_BPHITBF): one hit of [`x_cell_per_hit_bf`]: [`XHitBf::test`], then
/// [`x_test_hit`] and the candidate's call for a flag-1 entry or a candidate.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[inline(always)]
fn x_hit_bf<F>(
    query_offset: i32,
    s: i32,
    here: u32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    bf: XHitBf,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let entry = table.entry_bf(query_offset, s);
    if bf.test(entry, here) {
        return;
    }
    if let Some((last_hit, diff)) =
        x_test_hit(entry, s, test.diag_offset, test.window, test.word_size)
    {
        candidate(query_offset, s, last_hit, diff, entry);
    }
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the hits of one batch, the lanes of an AVX2 register.
const X_DIAG_LANES: usize = 8;

/// EXPERIMENT (LOSAT_X_BPDIAG): the 31-bit `last_hit` field of a packed entry.
const X_LAST_HIT_MASK: u32 = 0x7fff_ffff;

// No NCBI counterpart: a check of the property the batches rest on (module header).
/// EXPERIMENT (LOSAT_X_BPDIAG): whether the entries of a batch's hits are distinct.
fn x_batch_coords_distinct(coords: &[usize]) -> bool {
    coords
        .iter()
        .enumerate()
        .all(|(i, coord)| !coords[i + 1..].contains(coord))
}

/// EXPERIMENT (LOSAT_X_BPDIAG): calls `candidate` for the lanes set in `lanes` (bit `i` = hit `i`
/// of `hits`) in lane order, which is the cell's order, with the values [`x_test_hit`] returns.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[inline(always)]
fn x_call_candidates<F>(
    mut lanes: u32,
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    while lanes != 0 {
        let lane = lanes.trailing_zeros() as usize;
        lanes &= lanes - 1;
        let query_offset = hits[lane];
        let entry = table.entry(query_offset, s);
        // The batch stored a candidate's old value back, so this is the entry the test read.
        let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
        let diff = s.wrapping_sub(last_hit);
        candidate(query_offset, s, last_hit, diff, entry);
    }
}

/// EXPERIMENT (LOSAT_X_BPDIAG): [`x_test_hit`] for up to [`X_DIAG_LANES`] hits of one cell: for
/// each hit the new entry and whether it is a candidate, without branches on the entry; every
/// entry is stored, then the candidates are called in order. Returns the number of candidates.
///
/// With `flag` the entry's bit 31, `last` its sign-extended 31-bit `last_hit`,
/// `sum = s + diag_offset` and `diff = s - (last - diag_offset)` (wrapping, as in `x_test_hit`):
/// - write = flag ? sum >= last : diff >= window; the new entry is then `sum` truncated to 31
///   bits with the flag cleared (`set_last_hit` + `set_flag(0)`, or `set_last_hit` on an entry
///   without the flag), otherwise the old entry;
/// - candidate = !flag && word_size <= diff < window (its entry is not written).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[inline(always)]
fn x_test_batch_scalar<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    debug_assert!(hits.len() <= X_DIAG_LANES);
    let subject_plus_offset = s.wrapping_add(test.diag_offset);
    let written = (subject_plus_offset as u32) & X_LAST_HIT_MASK;
    let mut lanes = 0u32;
    for (lane, &query_offset) in hits.iter().enumerate() {
        let entry = table.entry(query_offset, s);
        let old = entry.raw_bits();
        let flag = old >> 31 != 0;
        let last = entry.last_hit();
        let diff = s.wrapping_sub(last.wrapping_sub(test.diag_offset));
        let write = if flag {
            subject_plus_offset >= last
        } else {
            diff >= test.window
        };
        entry.set_raw_bits(if write { written } else { old });
        let is_candidate = !flag & (diff >= test.word_size) & (diff < test.window);
        lanes |= u32::from(is_candidate) << lane;
    }
    debug_assert!(x_batch_coords_distinct(
        &hits
            .iter()
            .map(|&query_offset| table.coord(query_offset, s))
            .collect::<Vec<_>>()
    ));
    x_call_candidates(lanes, hits, s, table, test, candidate);
    lanes.count_ones()
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the cell's hits in batches of [`X_DIAG_LANES`], each by
/// [`x_test_batch_scalar`]; returns the number of candidates. The batches read no entry ahead
/// (LOSAT_X_BPSEED_PREFETCH is the per-hit loop's): prefetching the batch 8 or 16 hits ahead was
/// not faster than none on 2026-10-14 (local, human block 1 x Species4 and the OrthoFinder
/// sample job, one thread).
#[inline(always)]
fn x_cell_scalar<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut candidates = 0;
    for batch in hits.chunks(X_DIAG_LANES) {
        candidates += x_test_batch_scalar(batch, s, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPDIAG): [`x_test_batch_scalar`] for exactly eight hits in AVX2 registers.
/// The lane arithmetic is the scalar code's: `_mm256_sub_epi32` wraps like `wrapping_sub`, and
/// `_mm256_cmpgt_epi32` is the signed i32 comparison. All eight entries are stored (the old value
/// where the test writes nothing): storing only the written lanes, from the movemask, was slower
/// on 2026-10-14 (local, human block 1 x Species4, one thread).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_test_batch_avx2<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let (_, lanes) = x_batch_avx2(hits, s, table, test);
    x_call_candidates(lanes, hits, s, table, test, candidate);
    lanes.count_ones()
}

/// The test of [`x_test_batch_avx2`] without the calls: stores the batch's entries and returns the
/// hits' register and the candidate lanes (bit `i` = hit `i`).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_batch_avx2(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> (std::arch::x86_64::__m256i, u32) {
    use std::arch::x86_64::*;
    assert!(hits.len() == X_DIAG_LANES);
    let subject_plus_offset = s.wrapping_add(test.diag_offset);
    let written = (subject_plus_offset as u32) & X_LAST_HIT_MASK;
    let query = _mm256_loadu_si256(hits.as_ptr().cast());
    // (q - s) & diag_mask per lane
    let coords = _mm256_and_si256(
        _mm256_sub_epi32(query, _mm256_set1_epi32(s)),
        _mm256_set1_epi32(table.mask),
    );
    let mut lane_coords = [0u32; X_DIAG_LANES];
    _mm256_storeu_si256(lane_coords.as_mut_ptr().cast(), coords);
    // The eight entries with scalar loads, not `_mm256_i32gather_epi32`: on the workstation (Xeon
    // 6338, Gather Data Sampling microcode mitigation on) the gather was slower, 2026-10-14, one
    // pinned core, scan_ungapped: human block 1 per-hit 11.59 s, gather 10.17 s, scalar loads
    // 9.13 s; C. elegans block 5 32.68, 28.55, 26.45 s. Locally (no mitigation) both were equal.
    let mut lane_old = [0u32; X_DIAG_LANES];
    for (old, &coord) in lane_old.iter_mut().zip(&lane_coords) {
        // coord <= mask < entries.len() (XDiagTable::new)
        *old = table.entries.get_unchecked(coord as usize).raw_bits();
    }
    let old = _mm256_loadu_si256(lane_old.as_ptr().cast());
    // flag: all ones where bit 31 is set; last: the 31-bit field sign-extended (last_hit()).
    let flag = _mm256_srai_epi32::<31>(old);
    let last = _mm256_srai_epi32::<1>(_mm256_slli_epi32::<1>(old));
    let diff = _mm256_sub_epi32(
        _mm256_set1_epi32(s),
        _mm256_sub_epi32(last, _mm256_set1_epi32(test.diag_offset)),
    );
    // keep (no write) = flag ? sum < last : diff < window
    let sum_below_last = _mm256_cmpgt_epi32(last, _mm256_set1_epi32(subject_plus_offset));
    let in_window = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.window), diff);
    let keep = _mm256_blendv_epi8(in_window, sum_below_last, flag);
    // candidate = !flag && diff >= word_size && diff < window
    let below_word = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.word_size), diff);
    let candidates = _mm256_andnot_si256(_mm256_or_si256(flag, below_word), in_window);
    let new = _mm256_blendv_epi8(_mm256_set1_epi32(written as i32), old, keep);
    let mut lane_new = [0u32; X_DIAG_LANES];
    _mm256_storeu_si256(lane_new.as_mut_ptr().cast(), new);
    for (&coord, &raw) in lane_coords.iter().zip(&lane_new) {
        // coord <= mask < entries.len()
        table
            .entries
            .get_unchecked_mut(coord as usize)
            .set_raw_bits(raw);
    }
    debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
    let lanes = _mm256_movemask_ps(_mm256_castsi256_ps(candidates)) as u32;
    (query, lanes)
}

/// EXPERIMENT (LOSAT_X_BPDEFER): for each mask of eight lanes, the set lanes in order, then zeros:
/// the permutation that moves a batch's candidates to the front of its register.
const X_PACK_LANES: [[u8; X_DIAG_LANES]; 256] = {
    let mut table = [[0u8; X_DIAG_LANES]; 256];
    let mut mask = 0;
    while mask < 256 {
        let mut n = 0;
        let mut lane = 0;
        while lane < X_DIAG_LANES {
            if mask >> lane & 1 != 0 {
                table[mask][n] = lane as u8;
                n += 1;
            }
            lane += 1;
        }
        mask += 1;
    }
    table
};

/// EXPERIMENT (LOSAT_X_BPDEFER): [`x_cell_avx2`] with the calls of the batches of eight after the
/// last of them (module header): each batch packs its candidates' query offsets into `buffer`, then
/// the closure is called for the buffer in order, as [`x_call_candidates`] calls it, and the last
/// batch of fewer than eight hits is a scalar batch. Returns the number of candidates.
///
/// # Safety
/// The CPU must have AVX2.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2_defer<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    buffer: &mut Vec<i32>,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let mut batches = hits.chunks_exact(X_DIAG_LANES);
    // Room for a full register past the last candidate (each batch stores eight lanes).
    buffer.clear();
    buffer.reserve(hits.len() + X_DIAG_LANES);
    let out = buffer.as_mut_ptr();
    let mut packed = 0usize;
    for batch in &mut batches {
        let (query, lanes) = x_batch_avx2(batch, s, table, test);
        let order = _mm256_cvtepu8_epi32(_mm_loadl_epi64(
            X_PACK_LANES[lanes as usize].as_ptr().cast(),
        ));
        // packed <= the hits of the batches before this one, so the eight stores stay inside
        // the reserved capacity.
        _mm256_storeu_si256(
            out.add(packed).cast(),
            _mm256_permutevar8x32_epi32(query, order),
        );
        packed += lanes.count_ones() as usize;
    }
    // The first `packed` values were stored above.
    buffer.set_len(packed);
    for &query_offset in buffer.iter() {
        let entry = table.entry(query_offset, s);
        // The batch stored a candidate's old value back, and no later batch of the cell touches
        // its entry, so this is the entry the test read.
        let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
        let diff = s.wrapping_sub(last_hit);
        candidate(query_offset, s, last_hit, diff, entry);
    }
    let mut candidates = packed as u32;
    let rest = batches.remainder();
    if !rest.is_empty() {
        candidates += x_test_batch_scalar(rest, s, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPMASKED, LOSAT_X_BXDEFER): the test of [`x_batch_avx2`] for eight lanes
/// whose entry coordinates the caller computed: reads the eight entries, stores their new values
/// and returns the candidate lanes (bit `i` = lane `i`). Every load comes before every store, so
/// lanes with the same coordinate (the masked batch of [`x_cell_avx2_defer_masked`]) read the same
/// value and store the same new one.
///
/// # Safety
/// The CPU must have AVX2, and every coordinate must be at most `table.mask`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_batch_avx2_at(
    lane_coords: [u32; X_DIAG_LANES],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> u32 {
    use std::arch::x86_64::*;
    let subject_plus_offset = s.wrapping_add(test.diag_offset);
    let written = (subject_plus_offset as u32) & X_LAST_HIT_MASK;
    // coord <= mask < entries.len() (the caller's contract; XDiagTable::new)
    let e = |i: usize| {
        table
            .entries
            .get_unchecked(lane_coords[i] as usize)
            .raw_bits() as i32
    };
    let old = _mm256_setr_epi32(e(0), e(1), e(2), e(3), e(4), e(5), e(6), e(7));
    // From here to the stores: x_batch_avx2's test.
    let flag = _mm256_srai_epi32::<31>(old);
    let last = _mm256_srai_epi32::<1>(_mm256_slli_epi32::<1>(old));
    let diff = _mm256_sub_epi32(
        _mm256_set1_epi32(s),
        _mm256_sub_epi32(last, _mm256_set1_epi32(test.diag_offset)),
    );
    let sum_below_last = _mm256_cmpgt_epi32(last, _mm256_set1_epi32(subject_plus_offset));
    let in_window = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.window), diff);
    let keep = _mm256_blendv_epi8(in_window, sum_below_last, flag);
    let below_word = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.word_size), diff);
    let candidates = _mm256_andnot_si256(_mm256_or_si256(flag, below_word), in_window);
    let new = _mm256_blendv_epi8(_mm256_set1_epi32(written as i32), old, keep);
    let mut lane_new = [0u32; X_DIAG_LANES];
    _mm256_storeu_si256(lane_new.as_mut_ptr().cast(), new);
    for (&coord, &raw) in lane_coords.iter().zip(&lane_new) {
        // coord <= mask < entries.len()
        table
            .entries
            .get_unchecked_mut(coord as usize)
            .set_raw_bits(raw);
    }
    _mm256_movemask_ps(_mm256_castsi256_ps(candidates)) as u32
}

/// EXPERIMENT (LOSAT_X_BPMASKED, LOSAT_X_BXDEFER): [`x_cell_avx2_defer`] with two changes of the
/// batches (module header), for blastp's and blastx's fused loops: a batch of eight takes its
/// lanes' entry coordinates from
/// scalar loads of its hits, and the last batch of fewer than eight hits is one masked AVX2 batch
/// (the missing lanes repeat the first hit, their candidate bits are dropped) whose candidates go
/// into the buffer after the others'. The closure is then called for the buffer in order, so the
/// calls are [`x_cell_avx2_defer`]'s, in the same order. Returns the number of candidates.
///
/// # Safety
/// The CPU must have AVX2.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2_defer_masked<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    buffer: &mut Vec<i32>,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let mut batches = hits.chunks_exact(X_DIAG_LANES);
    // Room for a full register past the last candidate (each batch stores eight lanes).
    buffer.clear();
    buffer.reserve(hits.len() + X_DIAG_LANES);
    let out = buffer.as_mut_ptr();
    let mut packed = 0usize;
    let mask = table.mask;
    for batch in &mut batches {
        let p = batch.as_ptr();
        let query = _mm256_loadu_si256(p.cast());
        // The coordinates from scalar loads of the hits, not taken out of `query`: the volatile
        // reads keep the compiler from extracting them from the register, which would start the
        // chain of the entry loads 3-6 cycles later (line R's probe of blastp's cell test,
        // 2026-10-18).
        // i < X_DIAG_LANES = batch.len()
        let lane_coords: [u32; X_DIAG_LANES] = std::array::from_fn(|i| {
            (std::ptr::read_volatile(p.add(i)).wrapping_sub(s) & mask) as u32
        });
        debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
        let lanes = x_batch_avx2_at(lane_coords, s, table, test);
        let order = _mm256_cvtepu8_epi32(_mm_loadl_epi64(
            X_PACK_LANES[lanes as usize].as_ptr().cast(),
        ));
        // packed <= the hits of the batches before this one, so the eight stores stay inside
        // the reserved capacity.
        _mm256_storeu_si256(
            out.add(packed).cast(),
            _mm256_permutevar8x32_epi32(query, order),
        );
        packed += lanes.count_ones() as usize;
    }
    let rest = batches.remainder();
    if !rest.is_empty() {
        let n = rest.len();
        let in_rest = _mm256_cmpgt_epi32(
            _mm256_set1_epi32(n as i32),
            _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7),
        );
        // Lanes 0..n hold the hits (the masked load reads nothing past them), the others the
        // first hit again: they read and store its entry as lane 0 does (x_batch_avx2_at).
        let query = _mm256_blendv_epi8(
            _mm256_set1_epi32(rest[0]),
            _mm256_maskload_epi32(rest.as_ptr(), in_rest),
            in_rest,
        );
        let mut lane_hits = [0i32; X_DIAG_LANES];
        _mm256_storeu_si256(lane_hits.as_mut_ptr().cast(), query);
        let lane_coords: [u32; X_DIAG_LANES] =
            std::array::from_fn(|i| (lane_hits[i].wrapping_sub(s) & mask) as u32);
        debug_assert!(x_batch_coords_distinct(
            &lane_coords.map(|c| c as usize)[..n]
        ));
        // The repeated lanes' candidate bits are lane 0's; only lanes 0..n are hits.
        let lanes = x_batch_avx2_at(lane_coords, s, table, test) & ((1u32 << n) - 1);
        let order = _mm256_cvtepu8_epi32(_mm_loadl_epi64(
            X_PACK_LANES[lanes as usize].as_ptr().cast(),
        ));
        // packed <= hits.len() - n, so the eight stores end inside the reserved capacity.
        _mm256_storeu_si256(
            out.add(packed).cast(),
            _mm256_permutevar8x32_epi32(query, order),
        );
        packed += lanes.count_ones() as usize;
    }
    // The first `packed` values were stored above.
    buffer.set_len(packed);
    for &query_offset in buffer.iter() {
        let entry = table.entry(query_offset, s);
        // The batch stored a candidate's old value back, and no other batch of the cell touches
        // its entry, so this is the entry the test read.
        let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
        let diff = s.wrapping_sub(last_hit);
        candidate(query_offset, s, last_hit, diff, entry);
    }
    packed as u32
}

/// EXPERIMENT (LOSAT_X_BXCROSS): stages one cell's hits (module header) in pieces of
/// [`X_CROSS_PIECE`]: per piece three masked loads of eight lanes of the list (the lanes past the
/// cell's hits read nothing), the coordinates `(q - s) & mask` by lane, and three stores into each
/// of the two arrays at the staged count, all from one base register at constant displacements;
/// then the count advances by the piece's hits, and the batches run ([`x_cross_drain`]) when
/// `CHUNK` or more are staged. A cell of up to 24 hits is one piece. The cell and its hits are
/// counted here, its candidates by the drain. Returns the number of hits staged.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, `staged` must be below `CHUNK` (as it is at the start of a
/// subject and after every call of this function), and every query offset of `hits` must be at
/// most `table.mask` ([`x_cross_query`]).
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-516 (the coordinate, quoted at
// `XDiagTable::coord`) and blast_aascan.c:60-99 (the copy of a cell's hits, quoted at
// `x_scan_cells`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
unsafe fn x_cross_cell<const DEFER: bool, const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    mut staged: usize,
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XCrossCounts,
) -> usize
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    // The stage keeps no query offset ([`x_cross_query`]): each must be at most the mask (the
    // bound `x_scan_cells_cross` asserts from the lookup's length), or the recovery gives another.
    debug_assert!(
        hits.iter().all(|q| (0..=table.mask).contains(q)),
        "a query offset of the cell at subject offset {s} is above the mask {}",
        table.mask
    );
    counts.diag.add_staged(hits.len());
    let subject = _mm256_set1_epi32(s);
    let mask = _mm256_set1_epi32(table.mask);
    // The positions in a piece of the lanes of its three registers.
    let positions = [
        _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7),
        _mm256_setr_epi32(8, 9, 10, 11, 12, 13, 14, 15),
        _mm256_setr_epi32(16, 17, 18, 19, 20, 21, 22, 23),
    ];
    let mut p = hits.as_ptr();
    let mut left = hits.len();
    loop {
        let piece = left.min(X_CROSS_PIECE);
        let n = _mm256_set1_epi32(piece as i32);
        // staged < CHUNK, so the piece's lanes end at staged + 24 < STRIDE in each array
        // (`XCross::new`); the two stores of a register are STRIDE values apart.
        let at = cross.coord().add(staged);
        for (k, &position) in positions.iter().enumerate() {
            // The lanes at positions below `piece`; the pointer may lie past the list
            // (`wrapping_add`), where the masked load reads nothing.
            let in_piece = _mm256_cmpgt_epi32(n, position);
            let query = _mm256_maskload_epi32(p.wrapping_add(X_DIAG_LANES * k), in_piece);
            let coords = _mm256_and_si256(_mm256_sub_epi32(query, subject), mask);
            let lane = X_DIAG_LANES * k;
            _mm256_storeu_si256(at.add(lane).cast(), coords);
            _mm256_storeu_si256(at.add(XCross::<CHUNK>::STRIDE + lane).cast(), subject);
        }
        staged += piece;
        if staged >= CHUNK {
            staged = x_cross_drain::<DEFER, CHUNK, F>(
                cross, staged, false, table, test, candidate, counts,
            );
        }
        if left <= X_CROSS_PIECE {
            break;
        }
        left -= X_CROSS_PIECE;
        p = p.add(X_CROSS_PIECE);
    }
    staged
}

/// EXPERIMENT (LOSAT_X_BXCROSS): tests the `staged` hits of the stage in windows of eight from the
/// first one: a window whose eight entries are distinct ([`x_cross_conflict`]) is one batch
/// ([`x_cross_batch`]); otherwise it is tested in parts ([`x_cross_window`], cold and out of line,
/// so that the loop carries nothing of it). Each batch's candidates are called before the next
/// batch loads. EXPERIMENT (LOSAT_X_BXCROSSDEFER, with `DEFER`): each unit's candidates are
/// deferred instead (the batch [`x_cross_batch_marked`] and [`x_cross_defer`], the window
/// [`x_cross_window_defer`]) and called in order ([`x_cross_flush`]) when a later unit meets a
/// deferred candidate's entry (the unit then runs again from the loop's head, so that nothing of
/// it is live across the calls), and after the last window. Then with `end` (the subject's scan is
/// over) the fewer than eight left go one after the other ([`x_cross_hit`], the reference's
/// order); otherwise they move to the front of the stage by one eight-lane load and store per
/// array (the lanes past them are never read), where the hits staged next follow them. The
/// candidates and the windows tested in parts are added to `counts` (once per drain). Returns the
/// hits left staged. Out of line: it runs once per chunk of staged hits and once per subject.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, `staged` must be below `CHUNK + X_CROSS_PIECE`, every
/// staged hit's query offset must be at most `table.mask`, and with `DEFER` the stage made by
/// `XCross::new::<true>`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline(never)]
unsafe fn x_cross_drain<const DEFER: bool, const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    staged: usize,
    end: bool,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XCrossCounts,
) -> usize
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let full = staged & !(X_DIAG_LANES - 1);
    let mut candidates = 0u32;
    let mut windows = 0u32;
    // The records of deferred calls (LOSAT_X_BXCROSSDEFER; always 0 without `DEFER`).
    let mut pending = 0usize;
    let mut b = 0;
    while b < full {
        let coords = _mm256_loadu_si256(cross.coord().add(b).cast());
        if x_cross_conflict(coords) {
            if DEFER {
                let (window, left) =
                    x_cross_window_defer(cross, b, pending, table, test, candidate);
                candidates += window;
                pending = left;
            } else {
                candidates += x_cross_window(cross, b, table, test, candidate);
            }
            windows += 1;
        } else if !DEFER {
            candidates += x_cross_batch(cross, b, table, test, candidate);
        } else {
            match x_cross_batch_marked(cross, b, table, test) {
                Some((lanes, old)) => {
                    x_cross_defer(cross, &mut pending, b, lanes, old);
                    candidates += lanes.count_ones();
                }
                None if pending != 0 => {
                    // A deferred candidate's entry (or a stored value equal to the marker): the
                    // calls, then the batch again from the loop's head.
                    x_cross_flush_cold(cross, pending, table, test, candidate);
                    pending = 0;
                    continue;
                }
                // A stored value equal to the marker, no call deferred: the batch with its calls
                // at once.
                None => candidates += x_cross_batch_cold(cross, b, table, test, candidate),
            }
        }
        b += X_DIAG_LANES;
    }
    if DEFER {
        // The deferred calls, before the hits left are tested (they read the entries) or moved
        // (over the staged hits the records name).
        x_cross_flush(cross, pending, table, test, candidate);
    }
    let left = if end {
        for i in full..staged {
            candidates += x_cross_hit(cross, i, table, test, candidate);
        }
        0
    } else {
        // Fewer than eight left from `full`: the loads end at full + 8 <= staged + 8 < CHUNK + 32
        // = STRIDE, inside each array. Each register is loaded before it is stored.
        for row in [cross.coord(), cross.subject()] {
            let lanes = _mm256_loadu_si256(row.add(full).cast());
            _mm256_storeu_si256(row.cast(), lanes);
        }
        staged - full
    };
    counts.diag.add_candidates(candidates);
    counts.splits += u64::from(windows);
    left
}

/// EXPERIMENT (LOSAT_X_BXCROSS): [`x_test_hit`] for staged hit `i`, its candidate called at once
/// (the reference's order) with the query offset of [`x_cross_query`]; returns 1 for a candidate,
/// else 0.
///
/// # Safety
/// `i` must be below the staged count, and the hit's query offset at most `table.mask`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn x_cross_hit<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    i: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let coord = *cross.coord().add(i);
    let s = *cross.subject().add(i);
    // coord <= mask < entries.len() (staged masked; XDiagTable::new)
    let entry = table.entries.get_unchecked_mut(coord as u32 as usize);
    match x_test_hit(entry, s, test.diag_offset, test.window, test.word_size) {
        Some((last_hit, diff)) => {
            let query_offset = x_cross_query(coord, s, table.mask);
            candidate(query_offset, s, last_hit, diff, entry);
            1
        }
        None => 0,
    }
}

// No NCBI counterpart: a check of the property the batches rest on (module header).
/// EXPERIMENT (LOSAT_X_BXCROSS): whether two of the eight coordinates are equal: each lane against
/// the lanes one to four places further round the register (a pair `d` places apart, `d` in 5-7,
/// is the pair `8 - d` places apart taken from its other lane), and one test of the four compares'
/// OR.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
unsafe fn x_cross_conflict(coords: std::arch::x86_64::__m256i) -> bool {
    use std::arch::x86_64::*;
    let equal_at = |r: i32| {
        let rotate = _mm256_and_si256(
            _mm256_add_epi32(
                _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7),
                _mm256_set1_epi32(r),
            ),
            _mm256_set1_epi32(7),
        );
        _mm256_cmpeq_epi32(coords, _mm256_permutevar8x32_epi32(coords, rotate))
    };
    let equal = _mm256_or_si256(
        _mm256_or_si256(equal_at(1), equal_at(2)),
        _mm256_or_si256(equal_at(3), equal_at(4)),
    );
    _mm256_testz_si256(equal, equal) == 0
}

// No NCBI counterpart: a check of the property the batches rest on (module header).
/// EXPERIMENT (LOSAT_X_BXCROSS): the end of the part of a window's eight coordinates that starts at
/// lane `start` (below 8): the first lane after it whose coordinate a lane of the part before it
/// holds, or 8. The lanes `start..end` hold distinct coordinates.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn x_cross_part_end(coords: &[i32; X_DIAG_LANES], start: usize) -> usize {
    (start + 1..X_DIAG_LANES)
        .find(|&lane| coords[start..lane].contains(&coords[lane]))
        .unwrap_or(X_DIAG_LANES)
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the test of [`x_batch_avx2_at`] for eight lanes with a subject
/// offset each (`s`, lane `i` the subject offset of the hit at `lane_coords[i]`): reads the eight
/// entries, stores their new values and returns the candidate lanes (bit `i` = lane `i`). The lane
/// formulas are `x_batch_avx2_at`'s with the broadcast subject offset replaced by `s`:
/// `subject_plus_offset = s + diag_offset`, `written = (s + diag_offset) & X_LAST_HIT_MASK`,
/// `diff = s - (last - diag_offset)`. Every load comes before every store, so lanes with the same
/// coordinate and subject offset (the repeated lanes of [`x_cross_window`]) read the same value
/// and store the same new one.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, and every coordinate must be at most `table.mask`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
unsafe fn x_batch_avx2_lanes(
    lane_coords: [u32; X_DIAG_LANES],
    s: std::arch::x86_64::__m256i,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> u32 {
    let Some((lanes, _)) = x_batch_avx2_lanes_mark::<false>(lane_coords, s, table, test) else {
        unreachable!("without the marker the batch always stores")
    };
    lanes
}

/// EXPERIMENT (LOSAT_X_BXCROSS): [`x_batch_avx2_lanes`]'s test, which returns with the candidate
/// lanes the eight entries it read. EXPERIMENT (LOSAT_X_BXCROSSDEFER, with `MARK`): it returns
/// `None` before any store where an entry holds [`X_CROSS_MARKER`], and its candidate lanes store
/// the marker instead of the entry they read (which they store without `MARK`).
///
/// # Safety
/// As [`x_batch_avx2_lanes`].
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
unsafe fn x_batch_avx2_lanes_mark<const MARK: bool>(
    lane_coords: [u32; X_DIAG_LANES],
    s: std::arch::x86_64::__m256i,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> Option<(u32, std::arch::x86_64::__m256i)> {
    use std::arch::x86_64::*;
    let diag_offset = _mm256_set1_epi32(test.diag_offset);
    let subject_plus_offset = _mm256_add_epi32(s, diag_offset);
    let written = _mm256_and_si256(
        subject_plus_offset,
        _mm256_set1_epi32(X_LAST_HIT_MASK as i32),
    );
    // coord <= mask < entries.len() (the caller's contract; XDiagTable::new)
    let e = |i: usize| {
        table
            .entries
            .get_unchecked(lane_coords[i] as usize)
            .raw_bits() as i32
    };
    let old = _mm256_setr_epi32(e(0), e(1), e(2), e(3), e(4), e(5), e(6), e(7));
    // LOSAT_X_BXCROSSDEFER: an entry of a deferred candidate (the marker) ends the batch before
    // any store; the caller makes the deferred calls and tests the eight again.
    if MARK {
        let marked = _mm256_cmpeq_epi32(old, _mm256_set1_epi32(X_CROSS_MARKER));
        if _mm256_movemask_ps(_mm256_castsi256_ps(marked)) != 0 {
            return None;
        }
    }
    // From here to the stores: x_batch_avx2's test, with `s` by lane.
    let flag = _mm256_srai_epi32::<31>(old);
    let last = _mm256_srai_epi32::<1>(_mm256_slli_epi32::<1>(old));
    let diff = _mm256_sub_epi32(s, _mm256_sub_epi32(last, diag_offset));
    let sum_below_last = _mm256_cmpgt_epi32(last, subject_plus_offset);
    let in_window = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.window), diff);
    let keep = _mm256_blendv_epi8(in_window, sum_below_last, flag);
    let below_word = _mm256_cmpgt_epi32(_mm256_set1_epi32(test.word_size), diff);
    let candidates = _mm256_andnot_si256(_mm256_or_si256(flag, below_word), in_window);
    let new = _mm256_blendv_epi8(written, old, keep);
    // With MARK a candidate lane (all ones) stores the marker, the others x_batch_avx2's new value.
    let new = if MARK {
        _mm256_or_si256(new, candidates)
    } else {
        new
    };
    let mut lane_new = [0u32; X_DIAG_LANES];
    _mm256_storeu_si256(lane_new.as_mut_ptr().cast(), new);
    for (&coord, &raw) in lane_coords.iter().zip(&lane_new) {
        // coord <= mask < entries.len()
        table
            .entries
            .get_unchecked_mut(coord as usize)
            .set_raw_bits(raw);
    }
    Some((
        _mm256_movemask_ps(_mm256_castsi256_ps(candidates)) as u32,
        old,
    ))
}

/// EXPERIMENT (LOSAT_X_BXCROSS): calls `candidate` for the candidate lanes of the batch at staged
/// hit `b` (bit `i` = staged hit `b + i`) in lane order, with the values [`x_test_hit`] returns (as
/// [`x_call_candidates`] calls them) and the query offset of [`x_cross_query`].
///
/// # Safety
/// `b + i` must be below the staged count for every lane `i` set in `lanes`, and those hits' query
/// offsets at most `table.mask`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn x_cross_calls<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    mut lanes: u32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    while lanes != 0 {
        let i = b + lanes.trailing_zeros() as usize;
        lanes &= lanes - 1;
        let coord = *cross.coord().add(i);
        let s = *cross.subject().add(i);
        let query_offset = x_cross_query(coord, s, table.mask);
        // coord <= mask < entries.len() (staged masked; XDiagTable::new). The batch stored a
        // candidate's old value back, and no other lane of the batch or call touches its entry,
        // so this is the entry the test read.
        let entry = table.entries.get_unchecked_mut(coord as u32 as usize);
        let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
        let diff = s.wrapping_sub(last_hit);
        candidate(query_offset, s, last_hit, diff, entry);
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
thread_local! {
    /// EXPERIMENT (LOSAT_X_BXCROSSDEFER, tests only): the units of this thread that met the marker
    /// with calls deferred (the calls made, the unit run again) and with none (a stored value, the
    /// unit run with its calls at once).
    static X_CROSS_MARKS_MET: std::cell::Cell<[u64; 2]> = const { std::cell::Cell::new([0; 2]) };
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER, tests only): counts a unit that met the marker with calls
/// deferred (`stored` false) or none (true).
#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn x_cross_mark_met(stored: bool) {
    #[cfg(test)]
    X_CROSS_MARKS_MET.with(|met| {
        let mut counts = met.get();
        counts[usize::from(stored)] += 1;
        met.set(counts);
    });
    #[cfg(not(test))]
    let _ = stored;
}

// No NCBI counterpart: the marker of an entry whose candidate's call is deferred.
/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): the raw value a candidate lane stores while its call is
/// deferred: flag 1 with last_hit -1. A call stores it only near offsets of 2^31 (a flagged
/// last_hit is `s_last_off - (word - 1) + offset`); a unit that meets it with no call deferred
/// takes it as a stored value and runs with its calls at once, so any table content gives the same
/// calls (module header).
#[cfg(target_arch = "x86_64")]
const X_CROSS_MARKER: i32 = -1;

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): records a unit's candidates for deferred calls: its first
/// staged hit `b`, its candidate lanes and the eight entries it read go to record `pending`
/// (written for every unit, without a branch), and the count advances when the unit has
/// candidates.
///
/// # Safety
/// The CPU must have AVX2, and `pending` must be below [`XCross::PENDING`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
unsafe fn x_cross_defer<const CHUNK: usize>(
    cross: XCross<'_, CHUNK>,
    pending: &mut usize,
    b: usize,
    lanes: u32,
    old: std::arch::x86_64::__m256i,
) {
    use std::arch::x86_64::*;
    let r = *pending;
    debug_assert!(r < XCross::<CHUNK>::PENDING);
    _mm256_storeu_si256(cross.pending_old().add(X_DIAG_LANES * r).cast(), old);
    *cross.pending_at().add(2 * r) = b as i32;
    *cross.pending_at().add(2 * r + 1) = lanes as i32;
    *pending = r + usize::from(lanes != 0);
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): calls `candidate` for the deferred records in order, each
/// record's candidate lanes in lane order (bit `i` = staged hit `b + i`), with the values
/// [`x_test_hit`] returns (as [`x_call_candidates`] calls them) and the query offset of
/// [`x_cross_query`]. A candidate's entry holds the marker since its unit; it gets back the value
/// the unit read, which the call then sees (as when the call followed the unit at once: no unit
/// between read the entry, or it would have met the marker). Leaves no record.
///
/// # Safety
/// The records must be those of [`x_cross_defer`] in this drain, and the staged hits they name
/// not yet moved.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[inline(never)]
unsafe fn x_cross_flush<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    pending: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    for r in 0..pending {
        let b = *cross.pending_at().add(2 * r) as usize;
        let mut lanes = *cross.pending_at().add(2 * r + 1) as u32;
        let old = cross.pending_old().add(X_DIAG_LANES * r);
        while lanes != 0 {
            let lane = lanes.trailing_zeros() as usize;
            lanes &= lanes - 1;
            let i = b + lane;
            let coord = *cross.coord().add(i);
            let s = *cross.subject().add(i);
            let query_offset = x_cross_query(coord, s, table.mask);
            // coord <= mask < entries.len() (staged masked; XDiagTable::new)
            let entry = table.entries.get_unchecked_mut(coord as u32 as usize);
            debug_assert_eq!(entry.raw_bits(), X_CROSS_MARKER as u32);
            entry.set_raw_bits(*old.add(lane) as u32);
            let last_hit = entry.last_hit().wrapping_sub(test.diag_offset);
            let diff = s.wrapping_sub(last_hit);
            candidate(query_offset, s, last_hit, diff, entry);
        }
    }
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_flush`] where a unit met the marker (rare), out of
/// the loop.
#[cfg(target_arch = "x86_64")]
#[cold]
#[inline(never)]
unsafe fn x_cross_flush_cold<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    pending: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_cross_mark_met(false);
    x_cross_flush(cross, pending, table, test, candidate);
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_batch_avx2_lanes_mark`] until it meets no marker (the
/// deferred calls made in between), then its candidates deferred ([`x_cross_defer`]); a marker met
/// with no call deferred is a stored value, and the unit then runs with its calls at once. Returns
/// the candidate lanes (within `keep`). For the parts of [`x_cross_window_defer`].
///
/// # Safety
/// As [`x_batch_avx2_lanes`], and `pending` below [`XCross::PENDING`] with the records of this
/// drain.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
#[allow(clippy::too_many_arguments)]
unsafe fn x_cross_unit<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    pending: &mut usize,
    b: usize,
    lane_coords: [u32; X_DIAG_LANES],
    s: std::arch::x86_64::__m256i,
    keep: u32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    loop {
        if let Some((lanes, old)) = x_batch_avx2_lanes_mark::<true>(lane_coords, s, table, test) {
            let lanes = lanes & keep;
            x_cross_defer(cross, pending, b, lanes, old);
            return lanes;
        }
        if *pending == 0 {
            // A stored value equal to the marker, no call deferred: the unit with its calls at
            // once (as without the deferral).
            x_cross_mark_met(true);
            let lanes = x_batch_avx2_lanes(lane_coords, s, table, test) & keep;
            x_cross_calls(cross, b, lanes, table, test, candidate);
            return lanes;
        }
        x_cross_flush_cold(cross, *pending, table, test, candidate);
        *pending = 0;
    }
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the staged hits `b..b + 8`, whose entries are distinct
/// ([`x_cross_conflict`]), as one batch (module header): the lanes' coordinates from scalar
/// loads of the stage, [`x_batch_avx2_lanes`] with the staged subject offsets, and the candidates
/// called in lane order ([`x_cross_calls`]). Returns the number of candidates.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, `b + 8` must be at most the staged count, and the hits'
/// query offsets at most `table.mask`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
unsafe fn x_cross_batch<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let p = cross.coord().add(b);
    // The coordinates from scalar loads, not taken out of the register of the conflict test (the
    // volatile reads keep the compiler from extracting them; x_cell_avx2_defer_masked).
    let lane_coords: [u32; X_DIAG_LANES] =
        std::array::from_fn(|i| std::ptr::read_volatile(p.add(i)) as u32);
    debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
    let s = _mm256_loadu_si256(cross.subject().add(b).cast());
    let lanes = x_batch_avx2_lanes(lane_coords, s, table, test);
    x_cross_calls(cross, b, lanes, table, test, candidate);
    lanes.count_ones()
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_batch`] with the marker
/// ([`x_batch_avx2_lanes_mark`]) and without the calls: returns its candidate lanes and the entries
/// it read, or `None` (nothing stored) where an entry holds the marker ([`x_cross_drain`] then makes
/// the deferred calls and runs the batch again).
///
/// # Safety
/// As [`x_cross_batch`].
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[inline]
unsafe fn x_cross_batch_marked<const CHUNK: usize>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> Option<(u32, std::arch::x86_64::__m256i)> {
    use std::arch::x86_64::*;
    let p = cross.coord().add(b);
    // The coordinates from scalar loads, as in x_cross_batch.
    let lane_coords: [u32; X_DIAG_LANES] =
        std::array::from_fn(|i| std::ptr::read_volatile(p.add(i)) as u32);
    debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
    let s = _mm256_loadu_si256(cross.subject().add(b).cast());
    x_batch_avx2_lanes_mark::<true>(lane_coords, s, table, test)
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_batch`] (no marker, its candidates called at
/// once) out of the drain's loop, for a batch that met a stored value equal to the marker while
/// no call was deferred. Returns the number of candidates.
///
/// # Safety
/// As [`x_cross_batch`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[cold]
#[inline(never)]
unsafe fn x_cross_batch_cold<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_cross_mark_met(true);
    x_cross_batch(cross, b, table, test, candidate)
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the staged hits `b..b + 8`, two of which name one entry (hits of
/// different cells on one diagonal), in parts (module header). From lane `start` (0 first), the
/// part is the lanes before the first one that repeats an entry of the part
/// ([`x_cross_part_end`]), so its entries are distinct, and the part as one batch gives what its
/// hits give one after the other. The batch's lanes past the part repeat the part's first hit (its
/// coordinate and subject offset): they read its entry and store, after it, the value it stores
/// (every load of [`x_batch_avx2_lanes`] comes before its stores), as the masked batch of
/// [`x_cell_avx2_defer_masked`] does; only the part's candidate bits are kept, and their
/// candidates are called in lane order ([`x_cross_calls`]) before the next part, which starts at
/// the repeating lane, loads. Nothing past `b + 8` is read. Returns the number of candidates. Cold
/// and never inlined, so that the drain's loop keeps none of its values.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, `b + 8` must be at most the staged count, and the hits'
/// query offsets at most `table.mask`.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[cold]
#[inline(never)]
unsafe fn x_cross_window<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let coords: [i32; X_DIAG_LANES] = std::array::from_fn(|i| *cross.coord().add(b + i));
    let subjects: [i32; X_DIAG_LANES] = std::array::from_fn(|i| *cross.subject().add(b + i));
    let mut candidates = 0;
    let mut start = 0;
    while start < X_DIAG_LANES {
        let end = x_cross_part_end(&coords, start);
        // Batch lane i: window lane start + i inside the part, else the part's first lane.
        let lane = |i: usize| if start + i < end { start + i } else { start };
        let lane_coords: [u32; X_DIAG_LANES] = std::array::from_fn(|i| coords[lane(i)] as u32);
        debug_assert!(x_batch_coords_distinct(
            &lane_coords.map(|c| c as usize)[..end - start]
        ));
        let lane_subjects: [i32; X_DIAG_LANES] = std::array::from_fn(|i| subjects[lane(i)]);
        let s = _mm256_loadu_si256(lane_subjects.as_ptr().cast());
        let lanes = x_batch_avx2_lanes(lane_coords, s, table, test) & ((1u32 << (end - start)) - 1);
        x_cross_calls(cross, b + start, lanes, table, test, candidate);
        candidates += lanes.count_ones();
        start = end;
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_window`] with the parts' candidates deferred
/// ([`x_cross_unit`]): a part that repeats a deferred candidate's entry (the next part starts at
/// the lane that repeats an entry of the one before) meets the marker, and the deferred calls are
/// made before it is tested. Takes and returns the drain's count of records; returns first the
/// number of candidates. Cold and never inlined, as [`x_cross_window`].
///
/// # Safety
/// As [`x_cross_window`], and `pending` below [`XCross::PENDING`] with the records of this drain.
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,popcnt")]
#[cold]
#[inline(never)]
unsafe fn x_cross_window_defer<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    mut pending: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> (u32, usize)
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use std::arch::x86_64::*;
    let coords: [i32; X_DIAG_LANES] = std::array::from_fn(|i| *cross.coord().add(b + i));
    let subjects: [i32; X_DIAG_LANES] = std::array::from_fn(|i| *cross.subject().add(b + i));
    let mut candidates = 0;
    let mut start = 0;
    while start < X_DIAG_LANES {
        let end = x_cross_part_end(&coords, start);
        // Batch lane i: window lane start + i inside the part, else the part's first lane.
        let lane = |i: usize| if start + i < end { start + i } else { start };
        let lane_coords: [u32; X_DIAG_LANES] = std::array::from_fn(|i| coords[lane(i)] as u32);
        debug_assert!(x_batch_coords_distinct(
            &lane_coords.map(|c| c as usize)[..end - start]
        ));
        let lane_subjects: [i32; X_DIAG_LANES] = std::array::from_fn(|i| subjects[lane(i)]);
        let s = _mm256_loadu_si256(lane_subjects.as_ptr().cast());
        let part = (1u32 << (end - start)) - 1;
        let lanes = x_cross_unit(
            cross,
            &mut pending,
            b + start,
            lane_coords,
            s,
            part,
            table,
            test,
            candidate,
        );
        candidates += lanes.count_ones();
        start = end;
    }
    (candidates, pending)
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the cell's hits in batches of eight by [`x_test_batch_avx2`], the
/// last batch of fewer than eight by [`x_test_batch_scalar`]; returns the number of candidates.
/// No entry is read ahead (see [`x_cell_scalar`]).
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2<F>(
    hits: &[i32],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut candidates = 0;
    let mut batches = hits.chunks_exact(X_DIAG_LANES);
    for batch in &mut batches {
        candidates += x_test_batch_avx2(batch, s, table, test, candidate);
    }
    let rest = batches.remainder();
    if !rest.is_empty() {
        candidates += x_test_batch_scalar(rest, s, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPSEED): the context of an offset of the concatenated query without the
/// binary search over the contexts: a coarse table gives the context at the start of the offset's
/// block, and a walk over the (small, contiguous) frame bases moves on to the last context starting
/// at or before the offset. The result is the largest context index whose `frame_base` is at most
/// the offset (0 below the first), which is what the binary search of
/// `blastp_context_idx_for_query_offset_with_bounds` returns. The coarse table has at most
/// `COARSE_ENTRIES` entries so that both arrays stay in the first cache levels: a table with one
/// entry per offset (as large as the diagonal table) cost a cache miss for every candidate hit.
// NCBI reference: c++/src/algo/blast/core/blast_query_info.c:211-232
// ```c
// Int4 BSearchContextInfo(Int4 n, const BlastQueryInfo * A)
// {
//     Int4 m=0, b=0, e=0, size=0;
//
//     size = A->last_context+1;
//     ...
//     while (b < e - 1) {
//         m = (b + e) / 2;
//         if (A->contexts[m].query_offset > n)
//             e = m;
//         else
//             b = m;
//     }
//     return b;
// }
// ```
pub(crate) struct XContextTable {
    /// Frame bases in context order, then `i32::MAX`.
    starts: Vec<i32>,
    /// Context index at offset `b << shift`.
    coarse: Vec<u32>,
    shift: u32,
    span: usize,
}

impl XContextTable {
    const COARSE_ENTRIES: usize = 1 << 14;

    /// `frame_bases` in context order (non-decreasing); `span` = the offsets to cover.
    pub(crate) fn new(frame_bases: &[i32], span: usize) -> Self {
        let mut shift = 0u32;
        while (span >> shift) >= Self::COARSE_ENTRIES {
            shift += 1;
        }
        Self::with_shift(frame_bases, span, shift)
    }

    pub(crate) fn with_shift(frame_bases: &[i32], span: usize, shift: u32) -> Self {
        let mut starts = Vec::with_capacity(frame_bases.len() + 2);
        starts.extend_from_slice(frame_bases);
        starts.push(i32::MAX);
        if frame_bases.is_empty() {
            // no context: every offset gives index 0, as the binary search does
            starts.push(i32::MAX);
        }
        let blocks = (span >> shift) + 1;
        let mut coarse = Vec::with_capacity(blocks);
        let mut idx = 0usize;
        for b in 0..blocks {
            let offset = i64::try_from(b << shift).expect("query offset must fit i64");
            while idx + 1 < frame_bases.len() && i64::from(starts[idx + 1]) <= offset {
                idx += 1;
            }
            coarse.push(u32::try_from(idx).expect("context index must fit u32"));
        }
        Self {
            starts,
            coarse,
            shift,
            span,
        }
    }

    /// The context of `query_offset`, or `None` outside the table.
    #[inline(always)]
    pub(crate) fn get(&self, query_offset: i32) -> Option<usize> {
        let q = usize::try_from(query_offset)
            .ok()
            .filter(|&q| q < self.span)?;
        let mut idx = self.coarse[q >> self.shift] as usize;
        // starts[last] = i32::MAX ends the walk; the coarse entry is a context index, so idx + 1
        // is always inside `starts`.
        while self.starts[idx + 1] <= query_offset {
            idx += 1;
        }
        Some(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- EXPERIMENT (LOSAT_X_BPDIAG): the batched cell test against the per-hit test ----

    /// One call of the candidate closure: query offset, subject offset, last_hit, diff, and the
    /// entry before the call.
    type Call = (i32, i32, i32, i32, u32);

    fn next(seed: &mut u64) -> u32 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (*seed >> 32) as u32
    }

    fn entry_of(raw: u32) -> DiagStruct {
        let mut entry = DiagStruct::default();
        entry.set_raw_bits(raw);
        entry
    }

    /// The modes this machine runs: the per-hit test (the reference), the scalar batch and, with
    /// AVX2, the AVX2 batches (blastp's and blastx's, LOSAT_X_BXCROSS's among them where the CPU
    /// also has POPCNT, with its calls deferred, LOSAT_X_BXCROSSDEFER, and not).
    fn cell_modes() -> Vec<XCellTest> {
        #[allow(unused_mut)]
        let mut modes = vec![XDiagMode::PerHit, XDiagMode::Scalar];
        #[cfg(target_arch = "x86_64")]
        if x_avx2_available() {
            modes.push(XDiagMode::Avx2);
            modes.push(XDiagMode::Avx2Defer);
            modes.push(XDiagMode::Avx2DeferMasked);
        } else {
            eprintln!("no AVX2 on this CPU: the AVX2 batch is not compared");
        }
        #[allow(unused_mut)]
        let mut modes: Vec<XCellTest> = modes.into_iter().map(XCellTest::Diag).collect();
        #[cfg(target_arch = "x86_64")]
        if x_cross_available() {
            modes.push(XCellTest::Cross { defer: false });
            modes.push(XCellTest::Cross { defer: true });
        }
        modes
    }

    /// The per-hit test, the reference of the other cell tests.
    const PER_HIT: XCellTest = XCellTest::Diag(XDiagMode::PerHit);

    /// The runs of [`cell_modes`] with the LOSAT_X_BXCROSS chunk each takes: LOSAT_X_BXCROSS's
    /// modes with the chunk of the test hook and with the built-in one, the other modes once (they
    /// take no chunk; 0).
    fn cell_runs() -> Vec<(XCellTest, usize)> {
        let mut runs = Vec::new();
        for mode in cell_modes() {
            #[cfg(target_arch = "x86_64")]
            if let XCellTest::Cross { .. } = mode {
                runs.push((mode, X_CROSS_TEST_CHUNK));
                runs.push((mode, X_CROSS_CHUNK));
                continue;
            }
            runs.push((mode, 0));
        }
        runs
    }

    /// The candidate closure of the tests: writes the candidate's entry as the extension does
    /// (with the flag for some), so a batch that stored its entries after the call would leave
    /// another table.
    fn extend(query_offset: i32, subject_offset: i32, entry: &mut DiagStruct, diag_offset: i32) {
        if (query_offset ^ subject_offset) & 1 == 0 {
            entry.set_flag(1);
            entry.set_last_hit(subject_offset.wrapping_add(diag_offset).wrapping_add(5));
        } else {
            entry.set_last_hit(subject_offset.wrapping_add(diag_offset));
        }
    }

    /// LOSAT_X_BXCROSS's cell test over `cells` as for one subject's scan, its calls deferred with
    /// `DEFER` (LOSAT_X_BXCROSSDEFER): every cell staged with the batches run at `CHUNK` staged
    /// hits, and the stage drained after the last; checks the cells and hits counted and returns
    /// the candidates and the windows tested in parts counted.
    #[cfg(target_arch = "x86_64")]
    fn run_cells_cross<const DEFER: bool, const CHUNK: usize, F>(
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
        candidate: &mut F,
    ) -> (u32, u64)
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        assert!(x_cross_available());
        assert!(
            cells
                .iter()
                .all(|(_, hits)| hits.iter().all(|q| (0..=table.mask).contains(q))),
            "the stage recovers query offsets up to the mask only"
        );
        let mut stage = Vec::new();
        let cross = XCross::<CHUNK>::new::<DEFER>(&mut stage);
        let mut counts = XCrossCounts::default();
        let mut staged = 0usize;
        for (s, hits) in cells {
            // SAFETY: the CPU has AVX2 and POPCNT, the query offsets are at most the mask
            // (asserted above), and `staged` is the count the last call returned (0 at first).
            staged = unsafe {
                x_cross_cell::<DEFER, CHUNK, F>(
                    cross,
                    staged,
                    hits,
                    *s,
                    table,
                    test,
                    candidate,
                    &mut counts,
                )
            };
        }
        // SAFETY: as above; `staged` is below the chunk, and the stage made for `DEFER`.
        let left = unsafe {
            x_cross_drain::<DEFER, CHUNK, F>(
                cross,
                staged,
                true,
                table,
                test,
                candidate,
                &mut counts,
            )
        };
        assert_eq!(left, 0, "the end of the subject drains the stage");
        assert_eq!(counts.diag.cells, cells.len() as u64, "cells staged");
        let hits: usize = cells.iter().map(|(_, hits)| hits.len()).sum();
        assert_eq!(counts.diag.hits, hits as u64, "hits staged");
        (counts.diag.candidates as u32, counts.splits)
    }

    /// [`run_cells_cross`] with the instances of `chunk` and `defer` (LOSAT_X_BXCROSSDEFER);
    /// returns the candidates and the windows tested in parts counted.
    #[cfg(target_arch = "x86_64")]
    fn run_cells_cross_with<F>(
        (chunk, defer): (usize, bool),
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
        candidate: &mut F,
    ) -> (u32, u64)
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        match (chunk, defer) {
            (X_CROSS_TEST_CHUNK, false) => {
                run_cells_cross::<false, X_CROSS_TEST_CHUNK, F>(table, test, cells, candidate)
            }
            (X_CROSS_TEST_CHUNK, true) => {
                run_cells_cross::<true, X_CROSS_TEST_CHUNK, F>(table, test, cells, candidate)
            }
            (X_CROSS_CHUNK, false) => {
                run_cells_cross::<false, X_CROSS_CHUNK, F>(table, test, cells, candidate)
            }
            (X_CROSS_CHUNK, true) => {
                run_cells_cross::<true, X_CROSS_CHUNK, F>(table, test, cells, candidate)
            }
            _ => panic!("no instance of chunk {chunk}"),
        }
    }

    /// LOSAT_X_BXCROSS's windows tested in parts for one subject's `cells` with the instance of
    /// `chunk` ([`run_cells_cross`] on a copy of `entries`), the same with the calls deferred
    /// (LOSAT_X_BXCROSSDEFER) and not. They depend on the entry coordinates only, so the test can
    /// work them out from the hits ([`CrossCases`]).
    #[cfg(target_arch = "x86_64")]
    fn cross_splits(
        entries: &[DiagStruct],
        mask: i32,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
        chunk: usize,
    ) -> u64 {
        let splits = [false, true].map(|defer| {
            let mut entries = entries.to_vec();
            let mut table = XDiagTable::new(&mut entries, mask);
            let mut candidate = |q: i32, s: i32, _: i32, _: i32, entry: &mut DiagStruct| {
                extend(q, s, entry, test.diag_offset);
            };
            run_cells_cross_with((chunk, defer), &mut table, test, cells, &mut candidate).1
        });
        assert_eq!(
            splits[0], splits[1],
            "chunk {chunk}: windows with deferred calls"
        );
        splits[0]
    }

    /// Tests the cells `(subject_offset, hits)` in order with `mode`; returns the closure's calls.
    /// A batch mode's candidate count must be its number of calls. LOSAT_X_BXCROSS's mode runs
    /// [`run_cells_cross`] with the instance of `chunk` (the other modes ignore it).
    fn run_cells(
        mode: XCellTest,
        chunk: usize,
        entries: &mut [DiagStruct],
        mask: i32,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
    ) -> Vec<Call> {
        let mut calls = Vec::new();
        let mut counted = 0u32;
        #[cfg(target_arch = "x86_64")]
        let mut buffer = Vec::new();
        #[cfg(not(target_arch = "x86_64"))]
        let _ = chunk;
        let mut table = XDiagTable::new(entries, mask);
        let mut candidate = |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
            calls.push((q, s, last_hit, diff, entry.raw_bits()));
            extend(q, s, entry, test.diag_offset);
        };
        let diag_mode = match mode {
            XCellTest::Diag(mode) => Some(mode),
            #[cfg(target_arch = "x86_64")]
            XCellTest::Cross { .. } => None,
        };
        #[cfg(target_arch = "x86_64")]
        if let XCellTest::Cross { defer } = mode {
            counted =
                run_cells_cross_with((chunk, defer), &mut table, test, cells, &mut candidate).0;
        }
        for (s, hits) in cells {
            // LOSAT_X_BXCROSS's cell test ran above.
            let Some(mode) = diag_mode else { break };
            match mode {
                XDiagMode::PerHit => x_cell_per_hit(hits, *s, &mut table, test, &mut candidate),
                XDiagMode::Scalar => {
                    counted += x_cell_scalar(hits, *s, &mut table, test, &mut candidate);
                }
                #[cfg(target_arch = "x86_64")]
                XDiagMode::Avx2 => {
                    assert!(x_avx2_available());
                    // SAFETY: the CPU has AVX2 (asserted above).
                    counted += unsafe { x_cell_avx2(hits, *s, &mut table, test, &mut candidate) };
                }
                #[cfg(target_arch = "x86_64")]
                XDiagMode::Avx2Defer => {
                    assert!(x_avx2_available());
                    // SAFETY: the CPU has AVX2 (asserted above).
                    counted += unsafe {
                        x_cell_avx2_defer(hits, *s, &mut table, test, &mut candidate, &mut buffer)
                    };
                }
                #[cfg(target_arch = "x86_64")]
                XDiagMode::Avx2DeferMasked => {
                    assert!(x_avx2_available());
                    // SAFETY: the CPU has AVX2 (asserted above).
                    counted += unsafe {
                        x_cell_avx2_defer_masked(
                            hits,
                            *s,
                            &mut table,
                            test,
                            &mut candidate,
                            &mut buffer,
                        )
                    };
                }
            }
        }
        if mode != PER_HIT {
            assert_eq!(counted as usize, calls.len(), "{mode:?}: candidate count");
        }
        calls
    }

    /// Every mode against the per-hit test from the same table: the same calls in the same order
    /// with the same arguments, and the same table bits. Returns the number of calls.
    fn same_in_every_mode(
        entries: &[DiagStruct],
        mask: i32,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
    ) -> usize {
        let mut reference = entries.to_vec();
        let reference_calls = run_cells(PER_HIT, 0, &mut reference, mask, test, cells);
        for (mode, chunk) in cell_runs() {
            let mut got = entries.to_vec();
            let calls = run_cells(mode, chunk, &mut got, mask, test, cells);
            assert_eq!(
                calls, reference_calls,
                "{mode:?} chunk {chunk} {test:?}: candidate calls"
            );
            let at = got
                .iter()
                .zip(&reference)
                .position(|(a, b)| a.raw_bits() != b.raw_bits());
            assert!(
                at.is_none(),
                "{mode:?} chunk {chunk} {test:?}: entry {at:?}: {:#x} reference {:#x}",
                got[at.unwrap_or(0)].raw_bits(),
                reference[at.unwrap_or(0)].raw_bits()
            );
        }
        reference_calls.len()
    }

    /// `n` distinct query offsets from `base .. base + size` taken modulo `size` (a power of two),
    /// so on distinct entries of a table of `size` entries and at most its mask, as the lookup's
    /// offsets are (LOSAT_X_BXCROSS recovers them from the entry coordinate), in random order.
    fn distinct_offsets(seed: &mut u64, size: usize, n: usize, base: i32) -> Vec<i32> {
        assert!(size.is_power_of_two());
        let mut all: Vec<i32> = (0..size as i32).collect();
        for i in 0..n {
            let j = i + next(seed) as usize % (size - i);
            all.swap(i, j);
        }
        all.truncate(n);
        all.into_iter()
            .map(|q| (q + base) & (size as i32 - 1))
            .collect()
    }

    /// Diagonal offsets with sums past 2^30 (the 31-bit field truncates) and past 2^31 - 1 (the
    /// i32 sum wraps) for subject offsets up to a few thousand.
    const DIAG_OFFSETS: [i32; 6] = [
        40,
        1_000,
        (1 << 30) - 600,
        (1 << 30) + 77,
        i32::MAX / 4,
        i32::MAX - 700,
    ];

    // Every last-batch size alone (1-7) and after full batches (9-15, 17).
    const HIT_COUNTS: [usize; 22] = [
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 24, 64, 333, 1024,
    ];

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[test]
    fn batched_cells_equal_the_per_hit_test_on_random_tables() {
        let mut seed = 0x5eed_u64;
        let size = 1024usize;
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for (window, word_size) in [(40, 3), (16, 2), (3, 5), (1, 0)] {
            for diag_offset in DIAG_OFFSETS {
                for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: [0, 3, 16][k % 3],
                    };
                    let s0 = (next(&mut seed) % 4000) as i32;
                    // Half the entries any bit pattern, half near the cells' subject offsets
                    // (flag set or not), so every outcome of the test occurs.
                    let entries: Vec<DiagStruct> = (0..size)
                        .map(|_| {
                            let r = next(&mut seed);
                            if r & 1 == 0 {
                                entry_of(next(&mut seed))
                            } else {
                                let back = (next(&mut seed) % (window as u32 + 20)) as i32 - 8;
                                let last = s0.wrapping_add(diag_offset).wrapping_sub(back);
                                entry_of((last as u32 & X_LAST_HIT_MASK) | (r & 2) << 30)
                            }
                        })
                        .collect();
                    // Later cells read what earlier cells (and their candidates) wrote.
                    let cells: Vec<(i32, Vec<i32>)> = [0, 1, 1, 3, 10, 41]
                        .iter()
                        .map(|&step| {
                            let base = (next(&mut seed) % 3) as i32 * 7;
                            (s0 + step, distinct_offsets(&mut seed, size, n, base))
                        })
                        .collect();
                    calls += same_in_every_mode(&entries, mask, test, &cells);
                }
            }
        }
        assert!(calls > 1000, "the tables must give candidates: {calls}");
    }

    /// The entry that gives the hit `(query_offset, s)` the outcome `case` of the test.
    fn adversarial_entry(case: usize, s: i32, test: XDiagTest) -> u32 {
        let sum = s.wrapping_add(test.diag_offset);
        // last_hit for a hit without the flag at distance `diff`: s - diff + diag_offset
        let at_distance = |diff: i32| s.wrapping_sub(diff).wrapping_add(test.diag_offset);
        let flagged = |last: i32| 0x8000_0000 | (last as u32 & X_LAST_HIT_MASK);
        let clear = |last: i32| last as u32 & X_LAST_HIT_MASK;
        match case % 11 {
            // flag set: s + diag_offset above, equal to, below last_hit
            0 => flagged(sum.wrapping_add(1)),
            1 => flagged(sum),
            2 => flagged(sum.wrapping_sub(1)),
            // flag clear: diff at the window and word size edges, 0, negative, far
            3 => clear(at_distance(test.window)),
            4 => clear(at_distance(test.window - 1)),
            5 => clear(at_distance(test.word_size)),
            6 => clear(at_distance(test.word_size - 1)),
            7 => clear(at_distance(0)),
            8 => clear(at_distance(-3)),
            9 => clear(at_distance(1 << 29)),
            _ => flagged(at_distance(test.window)),
        }
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[test]
    fn batched_cells_equal_the_per_hit_test_at_every_edge_of_the_test() {
        let mut seed = 0xed9e_u64;
        let size = 2048usize;
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for (window, word_size) in [(40, 3), (4, 3), (2, 2)] {
            for diag_offset in DIAG_OFFSETS {
                for s in [0, 5, 599, 600, 601, 699, 700, 701, 3000] {
                    for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                        let test = XDiagTest {
                            diag_offset,
                            window,
                            word_size,
                            prefetch: [16, 0, 1][k % 3],
                        };
                        let hits = distinct_offsets(&mut seed, size, n, 0);
                        let mut entries: Vec<DiagStruct> =
                            (0..size).map(|_| entry_of(next(&mut seed))).collect();
                        let rotate = next(&mut seed) as usize;
                        for (j, &q) in hits.iter().enumerate() {
                            let coord = ((q - s) & mask) as usize;
                            entries[coord] = entry_of(adversarial_entry(j + rotate, s, test));
                        }
                        // The same hits again one and window offsets later: the entries the first
                        // cell and its candidates wrote are tested again.
                        let cells =
                            vec![(s, hits.clone()), (s + 1, hits.clone()), (s + window, hits)];
                        calls += same_in_every_mode(&entries, mask, test, &cells);
                    }
                }
            }
        }
        assert!(calls > 1000, "the edges must give candidates: {calls}");
    }

    // ---- EXPERIMENT (LOSAT_X_BPMASKED): the masked batch against the deferred calls ----

    /// blastp's default cell test (the masked batch) against the one `LOSAT_X_BPMASKED=0` runs
    /// (LOSAT_X_BPDEFER's deferred calls), from the same table, on cells of every size, with
    /// entries at every edge of the test and random bits, and later cells on the entries earlier
    /// cells and their calls wrote: the same calls in the same order with the same arguments, and
    /// the same bits in every entry of the table.
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn the_masked_batch_equals_the_deferred_calls() {
        if !x_avx2_available() {
            return;
        }
        let deferred = XCellTest::Diag(XDiagMode::Avx2Defer);
        let masked = XCellTest::Diag(XDiagMode::Avx2DeferMasked);
        let mut seed = 0xa5ced_u64;
        let size = 2048usize;
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for (window, word_size) in [(40, 3), (4, 3), (2, 2)] {
            for diag_offset in DIAG_OFFSETS {
                for n in HIT_COUNTS {
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: 0,
                    };
                    let s0 = (next(&mut seed) % 4000) as i32;
                    let entries: Vec<DiagStruct> = (0..size)
                        .map(|j| {
                            let r = next(&mut seed);
                            if r & 3 == 0 {
                                entry_of(next(&mut seed))
                            } else {
                                entry_of(adversarial_entry(j + r as usize, s0, test))
                            }
                        })
                        .collect();
                    let cells: Vec<(i32, Vec<i32>)> = [0, 1, 2, window]
                        .iter()
                        .map(|&step| (s0 + step, distinct_offsets(&mut seed, size, n, 0)))
                        .collect();
                    let mut want_table = entries.clone();
                    let want = run_cells(deferred, 0, &mut want_table, mask, test, &cells);
                    let mut got_table = entries;
                    let got = run_cells(masked, 0, &mut got_table, mask, test, &cells);
                    assert_eq!(got, want, "{n} hits {test:?}: candidate calls");
                    let at = got_table
                        .iter()
                        .zip(&want_table)
                        .position(|(a, b)| a.raw_bits() != b.raw_bits());
                    assert!(at.is_none(), "{n} hits {test:?}: entry {at:?}");
                    calls += want.len();
                }
            }
        }
        assert!(calls > 1000, "the tables must give candidates: {calls}");
    }

    // ---- EXPERIMENT (LOSAT_X_BPHITBF): the branch-free test against the per-hit test ----

    /// One hit by the reference (`x_test_hit`) and by `bf` then `x_test_hit` where `bf` leaves the
    /// hit: the entry's bits and the result of each, and whether `bf` took the hit.
    #[allow(clippy::type_complexity)]
    fn hit_outcomes(
        raw: u32,
        s: i32,
        test: XDiagTest,
        bf: XHitBf,
    ) -> ((u32, Option<(i32, i32)>), (u32, Option<(i32, i32)>), bool) {
        let (offset, window, word) = (test.diag_offset, test.window, test.word_size);
        let mut reference = entry_of(raw);
        let expected = x_test_hit(&mut reference, s, offset, window, word);
        let mut got = entry_of(raw);
        let taken = bf.test(&mut got, bf.here(s));
        let result = if taken {
            None
        } else {
            x_test_hit(&mut got, s, offset, window, word)
        };
        (
            (reference.raw_bits(), expected),
            (got.raw_bits(), result),
            taken,
        )
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[test]
    fn the_branch_free_test_gives_the_reference_outcomes_in_its_range() {
        // too far, overlap, candidate, flag 1 kept, flag 1 reset
        let mut reached = [0usize; 5];
        let pairs = [
            (2, 1),
            (3, 1),
            (3, 2),
            (40, 3),
            (41, 40),
            (1 << 20, 6),
            ((1 << 28) - 1, 1),
            ((1 << 28) - 1, (1 << 28) - 2),
        ];
        // A one-word window, a word size of 0 or of the window or more: no test.
        for (window, word_size) in [(3, 3), (40, 0), (2, 3), (1, 0)] {
            assert!(XHitBf::new(40, 300, window, word_size).is_none());
        }
        for (window, word_size) in pairs {
            // The longest subject the rule allows (at the offset 0).
            let longest = ((1i64 << 30) - 2 * i64::from(window) - 1) as usize;
            for length in [1usize, 2, 300, 1 << 20, longest] {
                // The most diagonal offset the rule allows for this subject, and no more.
                let most = ((1i64 << 30) - 2 * i64::from(window) - length as i64 - 1) as i32;
                assert!(XHitBf::new(most + 1, length, window, word_size).is_none());
                assert!(XHitBf::new(-1, length, window, word_size).is_none());
                let offsets = [0, 1, window, most / 2, most];
                for diag_offset in offsets.into_iter().filter(|&offset| offset <= most) {
                    let bf = XHitBf::new(diag_offset, length, window, word_size).expect("in range");
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: 0,
                    };
                    // The flag-0 values a table holds (module header): [-window, diag_offset +
                    // length).
                    let (lowest, highest) = (-window, diag_offset + length as i32 - 1);
                    for s in [0, length as i32 - 1] {
                        let here = s + diag_offset;
                        let flag0 = [
                            here - (word_size - 1),
                            here - word_size,
                            here - (window - 1),
                            here - window,
                            here - window - 1,
                            here,
                            here + 1,
                            0,
                            lowest,
                            highest,
                        ];
                        let flag1 = [here + 1, here, here - 1, lowest, highest];
                        let entries = flag0
                            .into_iter()
                            .filter(|last| (lowest..=highest).contains(last))
                            .map(|last| last as u32 & X_LAST_HIT_MASK)
                            .chain(flag1.map(|last| 0x8000_0000 | (last as u32 & X_LAST_HIT_MASK)));
                        for raw in entries {
                            let (expected, got, taken) = hit_outcomes(raw, s, test, bf);
                            assert_eq!(
                                got, expected,
                                "window {window} word {word_size} length {length} offset {diag_offset} s {s} raw {raw:#x}"
                            );
                            let outcome = match (raw >> 31, expected) {
                                (0, (_, Some(_))) => 2,
                                (0, (after, None)) if after != raw => 0,
                                (0, _) => 1,
                                (_, (after, _)) if after == raw => 3,
                                _ => 4,
                            };
                            assert_eq!(taken, outcome < 2, "raw {raw:#x}: taken");
                            reached[outcome] += 1;
                        }
                    }
                }
            }
        }
        assert!(reached.iter().all(|&n| n > 0), "every outcome: {reached:?}");
        // Just past the range the doubled difference wraps: a cleared entry at the offset 2^30 is
        // too far for the reference and an overlap for the test.
        let (window, word_size, diag_offset) = (40, 3, 1 << 30);
        let bf = XHitBf {
            diag_offset,
            word2: 2 * word_size,
            span2: 2 * (window - word_size) as u32,
            window2: 2 * window,
        };
        let test = XDiagTest {
            diag_offset,
            window,
            word_size,
            prefetch: 0,
        };
        let (expected, got, _) = hit_outcomes(-window as u32 & X_LAST_HIT_MASK, 0, test, bf);
        assert_ne!(got, expected, "outside the range");
    }

    /// [`run_cells`] with the per-hit loop of the branch-free test.
    fn run_cells_bf(
        entries: &mut [DiagStruct],
        mask: i32,
        test: XDiagTest,
        bf: XHitBf,
        cells: &[(i32, Vec<i32>)],
    ) -> Vec<Call> {
        let mut calls = Vec::new();
        let mut table = XDiagTable::new(entries, mask);
        let mut candidate = |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
            calls.push((q, s, last_hit, diff, entry.raw_bits()));
            extend(q, s, entry, test.diag_offset);
        };
        for (s, hits) in cells {
            x_cell_per_hit_bf(hits, *s, &mut table, test, bf, &mut candidate);
        }
        calls
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    #[test]
    fn the_branch_free_per_hit_loop_equals_the_per_hit_loop() {
        let mut seed = 0xb7_u64;
        let size = 1024usize;
        let mask = (size - 1) as i32;
        let length = 5_000usize;
        let mut calls = 0;
        for (window, word_size) in [(40, 3), (16, 2), (4, 3), (2, 1)] {
            // From the first offset to the most the rule allows for the subject.
            let most = ((1i64 << 30) - 2 * i64::from(window) - length as i64 - 1) as i32;
            for diag_offset in [window, 1_000, i32::MAX / 4, most] {
                let bf = XHitBf::new(diag_offset, length, window, word_size).expect("in range");
                let (lowest, highest) = (-window, diag_offset + length as i32 - 1);
                for (k, n) in HIT_COUNTS.into_iter().enumerate() {
                    // The read-ahead distances of the other tests, and one equal to the cell's hits
                    // and one past them (no hit reads ahead).
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: [0, 3, 16, n, n + 1][k % 5],
                    };
                    let s0 = (next(&mut seed) % 4000) as i32;
                    // The entries a table holds (module header): flag 0 with last_hit in [-window,
                    // diag_offset + length), most of them near the cells' subject offsets, or
                    // flag 1.
                    let entries: Vec<DiagStruct> = (0..size)
                        .map(|_| {
                            let r = next(&mut seed);
                            let last = if r % 4 == 0 {
                                lowest + (next(&mut seed) % (highest - lowest + 1) as u32) as i32
                            } else {
                                let back = (next(&mut seed) % (window as u32 + 20)) as i32 - 8;
                                (s0 + diag_offset - back).clamp(lowest, highest)
                            };
                            let flag = if r % 4 == 1 { 0x8000_0000 } else { 0 };
                            entry_of(flag | (last as u32 & X_LAST_HIT_MASK))
                        })
                        .collect();
                    // Later cells read what earlier cells (and their candidates) wrote.
                    let cells: Vec<(i32, Vec<i32>)> = [0, 1, 1, 3, 10, 41]
                        .iter()
                        .map(|&step| {
                            let base = (next(&mut seed) % 3) as i32 * 7;
                            (s0 + step, distinct_offsets(&mut seed, size, n, base))
                        })
                        .collect();
                    let mut reference = entries.clone();
                    let expected = run_cells(PER_HIT, 0, &mut reference, mask, test, &cells);
                    let mut got = entries.clone();
                    let got_calls = run_cells_bf(&mut got, mask, test, bf, &cells);
                    assert_eq!(got_calls, expected, "{test:?}: candidate calls");
                    let at = got
                        .iter()
                        .zip(&reference)
                        .position(|(a, b)| a.raw_bits() != b.raw_bits());
                    assert!(at.is_none(), "{test:?}: entry {at:?}");
                    calls += expected.len();
                }
            }
        }
        assert!(calls > 1000, "the tables must give candidates: {calls}");
    }

    /// One subject's scan with the cell test `mode` and the scan of `branch_free`
    /// ([`x_scan_two_hit_in_with`]; BXSCANBF, graduated for blastx), LOSAT_X_BXCROSS's cell test as blastx's
    /// entry runs it ([`x_scan_cells_cross`]) with the instance of `chunk`; returns the counts of
    /// the cells tested in batches.
    #[allow(clippy::too_many_arguments)]
    fn scan_subject<F>(
        (mode, chunk, branch_free): (XCellTest, usize, bool),
        lookup: &BlastAaLookupTable,
        subject: &[u8],
        subject_length: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> XDiagCounts
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        let mode = match mode {
            XCellTest::Diag(mode) => mode,
            #[cfg(target_arch = "x86_64")]
            XCellTest::Cross { defer } => {
                return match chunk {
                    X_CROSS_TEST_CHUNK => x_scan_cells_cross::<X_CROSS_TEST_CHUNK, F>(
                        branch_free,
                        defer,
                        lookup,
                        subject,
                        subject_length,
                        table,
                        test,
                        candidate,
                    ),
                    X_CROSS_CHUNK => x_scan_cells_cross::<X_CROSS_CHUNK, F>(
                        branch_free,
                        defer,
                        lookup,
                        subject,
                        subject_length,
                        table,
                        test,
                        candidate,
                    ),
                    _ => panic!("no instance of chunk {chunk}"),
                }
                .diag;
            }
        };
        #[cfg(not(target_arch = "x86_64"))]
        let _ = chunk;
        if branch_free {
            x_scan_two_hit_in_with::<true, F>(
                mode,
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
            )
        } else {
            x_scan_two_hit_in(
                mode,
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
            )
        }
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:495-512 and blast_aascan.c:60-99
    // (quoted at `x_scan_two_hit` and `x_scan_cells`)
    #[test]
    fn batched_scans_of_real_lookup_tables_equal_the_per_hit_scan() {
        use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
        use crate::algorithm::tblastx::translation::QueryFrame;
        use crate::config::ScoringMatrix;
        use crate::stats::lookup_protein_params_ungapped;

        // NCBISTDAA codes of the 20 standard amino acids
        const STANDARD: [u8; 20] = [
            1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
        ];
        let mut seed = 0x1007_u64;
        let mut protein = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| STANDARD[next(&mut seed) as usize % STANDARD.len()])
                .collect()
        };
        let mut queries: Vec<Vec<u8>> = (0..6).map(|i| protein(80 + 90 * i)).collect();
        // Low complexity: cells far above AA_HITS_PER_CELL (the overflow array).
        let mut low = vec![11u8; 150];
        low.extend(protein(30));
        low.extend(std::iter::repeat_n([17u8, 18], 60).flatten());
        queries.push(low);
        let frames: Vec<Vec<QueryFrame>> = queries
            .iter()
            .map(|seq| {
                let mut aa_seq = vec![0u8];
                aa_seq.extend_from_slice(seq);
                aa_seq.push(0);
                vec![QueryFrame {
                    frame: 1,
                    aa_seq,
                    aa_seq_nomask: None,
                    aa_len: seq.len(),
                    orig_len: seq.len(),
                    seg_masks: Vec::new(),
                }]
            })
            .collect();
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (lookup, _) = build_ncbi_lookup(&frames, 11, &params, false);
        assert!(
            lookup
                .backbone
                .iter()
                .any(|cell| cell.num_used as usize > 3 * X_DIAG_LANES),
            "the query set must give long cells"
        );
        // Subjects: unrelated proteins, copies of the queries with every 5th letter changed (two
        // hits on a diagonal, extensions), the low-complexity query, and short ones.
        let mut subjects: Vec<Vec<u8>> = queries
            .iter()
            .map(|q| {
                q.iter()
                    .enumerate()
                    .map(|(i, &a)| if i % 5 == 4 { STANDARD[i % 20] } else { a })
                    .collect()
            })
            .collect();
        subjects.extend((0..4).map(|i| protein(50 + 200 * i)));
        subjects.push(queries[6].clone());
        subjects.push(protein(1));
        subjects.push(protein(3));
        let window = 40;
        let mut size = 1usize;
        while size < lookup.query_length as usize + window as usize {
            size <<= 1;
        }
        let mask = (size - 1) as i32;
        let mut calls = 0;
        for start in [window, (1 << 30) - 2_000, i32::MAX - 3_000] {
            for prefetch in [16, 0, 5] {
                let mut results = Vec::new();
                // Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): every cell test also with the branch-free scan
                // (LOSAT_X_BXCROSS's with each chunk).
                let runs = cell_runs()
                    .into_iter()
                    .flat_map(|(mode, chunk)| [(mode, chunk, false), (mode, chunk, true)]);
                for (mode, chunk, branch_free) in runs {
                    let mut entries = vec![DiagStruct::clear(window); size];
                    let mut log: Vec<(usize, Call)> = Vec::new();
                    let mut diag_offset: i32 = start;
                    let mut counts = (0u64, 0u64, 0u64);
                    for (k, subject) in subjects.iter().enumerate() {
                        let mut buffer = subject.clone();
                        buffer.push(0);
                        let test = XDiagTest {
                            diag_offset,
                            window,
                            word_size: 3,
                            prefetch,
                        };
                        let mut candidate =
                            |q: i32, s: i32, last_hit: i32, diff: i32, entry: &mut DiagStruct| {
                                log.push((k, (q, s, last_hit, diff, entry.raw_bits())));
                                extend(q, s, entry, diag_offset);
                            };
                        let table = &mut XDiagTable::new(&mut entries, mask);
                        let c = scan_subject(
                            (mode, chunk, branch_free),
                            &lookup,
                            &buffer,
                            subject.len(),
                            table,
                            test,
                            &mut candidate,
                        );
                        counts = (
                            counts.0 + c.cells,
                            counts.1 + c.hits,
                            counts.2 + c.candidates,
                        );
                        diag_offset = diag_offset.wrapping_add(subject.len() as i32 + window);
                    }
                    if mode == PER_HIT {
                        assert_eq!(counts, (0, 0, 0), "the per-hit loop counts nothing");
                    } else {
                        assert!(counts.0 > 0 && counts.1 >= counts.0, "{mode:?}: {counts:?}");
                        assert_eq!(counts.2 as usize, log.len(), "{mode:?}: candidates");
                    }
                    let bits: Vec<u32> = entries.iter().map(|e| e.raw_bits()).collect();
                    results.push(((mode, chunk, branch_free), log, bits, counts));
                }
                let (_, reference_log, reference_bits, _) = &results[0];
                // Every batch mode counts the same cells, hits and candidates (LOSAT_X_BXCROSS's
                // as it stages the cells and drains the stage, the others per cell).
                let batched = results
                    .iter()
                    .find(|(run, ..)| run.0 != PER_HIT)
                    .map(|(.., counts)| *counts);
                for (mode, log, bits, counts) in &results[1..] {
                    assert!(
                        log == reference_log,
                        "{mode:?} start {start} prefetch {prefetch}: candidate calls differ"
                    );
                    assert!(
                        bits == reference_bits,
                        "{mode:?} start {start} prefetch {prefetch}: tables differ"
                    );
                    if mode.0 != PER_HIT {
                        assert_eq!(Some(*counts), batched, "{mode:?}: counts");
                    }
                }
                calls += reference_log.len();
            }
        }
        assert!(calls > 100, "the subjects must give candidates: {calls}");
    }

    // NCBI reference: c++/src/algo/blast/core/blast_aascan.c:60-99 (quoted at `x_scan_cells`)
    // Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): the branch-free scan calls `on_cell` with the positions and
    // hit lists of the presence-test loop, on subjects across the block boundaries, with the
    // letters of translated frames (X, stop; all of 1-27) and with the word sizes 2-4 (the general
    // loop).
    #[test]
    fn the_branch_free_scan_finds_the_cells_of_the_presence_test() {
        use crate::algorithm::tblastx::lookup::build_ncbi_lookup_for_profile;
        use crate::algorithm::tblastx::translation::QueryFrame;
        use crate::config::ScoringMatrix;
        use crate::stats::lookup_protein_params_ungapped;

        let mut seed = 0x5ca7_u64;
        let mut letters = |n: usize, top: u64| -> Vec<u8> {
            (0..n)
                .map(|_| 1 + (u64::from(next(&mut seed)) % top) as u8)
                .collect()
        };
        let mut compared = 0usize;
        for (case, word_length) in [3, 3, 2, 4, 3].into_iter().enumerate() {
            let frames: Vec<Vec<QueryFrame>> = (0..1 + 3 * case)
                .map(|i| {
                    let seq = letters(60 + 70 * i, 22);
                    let mut aa_seq = vec![0u8];
                    aa_seq.extend_from_slice(&seq);
                    aa_seq.push(0);
                    vec![QueryFrame {
                        frame: 1,
                        aa_seq,
                        aa_seq_nomask: None,
                        aa_len: seq.len(),
                        orig_len: seq.len(),
                        seg_masks: Vec::new(),
                    }]
                })
                .collect();
            let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
            let (lookup, _) = build_ncbi_lookup_for_profile(
                &frames,
                [11, 13][case % 2],
                &params,
                false,
                ScoringMatrix::Blosum62,
                word_length,
            );
            for length in [0, 1, 2, 3, 4, 255, 256, 257, 258, 511, 512, 513, 2_000] {
                let mut subject = letters(length, 27);
                subject.push(0);
                let mut want: Vec<(i32, Vec<i32>)> = Vec::new();
                let mut got: Vec<(i32, Vec<i32>)> = Vec::new();
                x_scan_cells(&lookup, &subject, length, |s, hits| {
                    want.push((s, hits.to_vec()))
                });
                x_scan_cells_with::<true, _>(&lookup, &subject, length, |s, hits| {
                    got.push((s, hits.to_vec()))
                });
                assert_eq!(got, want, "case {case} length {length}");
                compared += want.len();
            }
        }
        assert!(compared > 1_000, "{compared} cells compared");
    }

    #[test]
    fn the_pack_table_lists_the_set_lanes_in_order() {
        for (mask, order) in X_PACK_LANES.iter().enumerate() {
            let lanes: Vec<u8> = (0..X_DIAG_LANES as u8)
                .filter(|&lane| mask >> lane & 1 != 0)
                .collect();
            assert_eq!(&order[..lanes.len()], &lanes[..], "mask {mask:#x}");
            assert!(order[lanes.len()..].iter().all(|&lane| lane == 0));
        }
    }

    #[test]
    fn the_batch_check_sees_a_repeated_entry() {
        assert!(x_batch_coords_distinct(&[]));
        assert!(x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 4]));
        assert!(!x_batch_coords_distinct(&[3, 1, 2, 0, 7, 6, 5, 3]));
    }

    // ---- EXPERIMENT (LOSAT_X_BXCROSS): the batches across the cells' borders ----

    /// What the cells of subjects give LOSAT_X_BXCROSS's stage, worked out from the hits and the
    /// reference's calls by following the pieces and drains of `x_cross_cell` and `x_cross_drain`
    /// and the parts of `x_cross_window`, so that a test can require each case to occur.
    #[cfg(target_arch = "x86_64")]
    #[derive(Debug, Default)]
    struct CrossCases {
        /// Windows of eight with a repeated entry (tested in parts).
        windows: usize,
        /// Parts of such windows cut short before a repeated entry, by the hits they take (1-7).
        cuts: [usize; X_DIAG_LANES],
        /// Parts that start past their window's first lane (after a cut).
        late_parts: usize,
        /// Windows of eight distinct entries with a candidate in a lane before a lane of a later
        /// cell.
        candidate_before_later_cell: usize,
        /// Cells of more than 24 hits (more than one piece), and of more than 48.
        long_cells: usize,
        longer_cells: usize,
        /// Drains while staging (the chunk reached), and those in the middle of a cell.
        drains: usize,
        drains_in_a_cell: usize,
        /// The hits a drain while staging moves to the front of the stage, by count (0-7).
        moved: [usize; X_DIAG_LANES],
        /// The hits left for the one-by-one test at the end of a subject, by count (0-7).
        left: [usize; X_DIAG_LANES],
    }

    #[cfg(target_arch = "x86_64")]
    impl CrossCases {
        /// Adds the cases of one subject's `cells`, whose reference calls are `calls`, with the
        /// batches run at `chunk` staged hits.
        fn add(&mut self, cells: &[(i32, Vec<i32>)], mask: i32, calls: &[Call], chunk: usize) {
            // The stream of staged hits (q, s) and their entries.
            let stream: Vec<(i32, i32)> = cells
                .iter()
                .flat_map(|(s, hits)| hits.iter().map(move |&q| (q, *s)))
                .collect();
            let coords: Vec<usize> = stream
                .iter()
                .map(|&(q, s)| (q.wrapping_sub(s) & mask) as usize)
                .collect();
            let called: std::collections::HashSet<(i32, i32)> =
                calls.iter().map(|&(q, s, ..)| (q, s)).collect();
            // x_cross_drain over the stream's hits `front..front + staged` (`front` at the stage's
            // first position): windows of eight from the front; returns the hits left.
            let drain = |cases: &mut Self, front: usize, staged: usize| {
                let full = staged - staged % X_DIAG_LANES;
                for b in (front..front + full).step_by(X_DIAG_LANES) {
                    let lanes = &coords[b..b + X_DIAG_LANES];
                    if (1..X_DIAG_LANES).any(|j| lanes[..j].contains(&lanes[j])) {
                        cases.windows += 1;
                        let mut start = 0;
                        while start < X_DIAG_LANES {
                            let mut seen = std::collections::HashSet::new();
                            let end = (start..X_DIAG_LANES)
                                .find(|&j| !seen.insert(lanes[j]))
                                .unwrap_or(X_DIAG_LANES);
                            if end < X_DIAG_LANES {
                                cases.cuts[end - start] += 1;
                            }
                            cases.late_parts += usize::from(start > 0);
                            start = end;
                        }
                        continue;
                    }
                    let batch = &stream[b..b + X_DIAG_LANES];
                    if batch.iter().enumerate().any(|(i, hit)| {
                        called.contains(hit) && batch[i + 1..].iter().any(|later| later.1 != hit.1)
                    }) {
                        cases.candidate_before_later_cell += 1;
                    }
                }
                staged - full
            };
            let (mut front, mut staged) = (0usize, 0usize);
            for (_, hits) in cells {
                self.long_cells += usize::from(hits.len() > X_CROSS_PIECE);
                self.longer_cells += usize::from(hits.len() > 2 * X_CROSS_PIECE);
                let mut left = hits.len();
                loop {
                    let piece = left.min(X_CROSS_PIECE);
                    staged += piece;
                    if staged >= chunk {
                        self.drains += 1;
                        self.drains_in_a_cell += usize::from(left > X_CROSS_PIECE);
                        let rest = drain(self, front, staged);
                        self.moved[rest] += 1;
                        front += staged - rest;
                        staged = rest;
                    }
                    if left <= X_CROSS_PIECE {
                        break;
                    }
                    left -= X_CROSS_PIECE;
                }
            }
            let rest = drain(self, front, staged);
            self.left[rest] += 1;
        }
    }

    /// The reference's calls for `cells` from `entries` (the per-hit test).
    #[cfg(target_arch = "x86_64")]
    fn reference_calls(
        entries: &[DiagStruct],
        mask: i32,
        test: XDiagTest,
        cells: &[(i32, Vec<i32>)],
    ) -> Vec<Call> {
        run_cells(PER_HIT, 0, &mut entries.to_vec(), mask, test, cells)
    }

    /// Entries near `s0`'s diagonal values (flag set or not) for half the table and any bits for
    /// the other half, so that every outcome of the test occurs.
    fn entries_near(seed: &mut u64, size: usize, s0: i32, test: XDiagTest) -> Vec<DiagStruct> {
        (0..size)
            .map(|_| {
                let r = next(seed);
                if r & 1 == 0 {
                    entry_of(next(seed))
                } else {
                    let back = (next(seed) % (test.window as u32 + 20)) as i32 - 8;
                    let last = s0.wrapping_add(test.diag_offset).wrapping_sub(back);
                    entry_of((last as u32 & X_LAST_HIT_MASK) | (r & 2) << 30)
                }
            })
            .collect()
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    // Subjects of 7 to 200 cells at consecutive subject offsets (some gaps), each cell of 0 to 40
    // hits or of 25 to 120 (more than one piece of the copy), one of 333 in some: the batches cross
    // the cells' borders, the stage is drained between cells and inside long ones (with the chunk
    // of the test hook and the built-in one), every count of hits is moved to the front by a drain
    // and left for the end. Later cells hit the entries earlier cells and their candidates wrote.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn cross_batches_equal_the_per_hit_test_on_cells_of_every_size() {
        if !x_cross_available() {
            eprintln!("no AVX2 and POPCNT on this CPU: the cross-cell batches are not compared");
            return;
        }
        let mut seed = 0xc055_u64;
        let size = 4096usize;
        let mask = (size - 1) as i32;
        let mut cases = CrossCases::default();
        let mut calls = 0;
        let mut sizes = std::collections::BTreeSet::new();
        for (window, word_size) in [(40, 3), (16, 2), (40, 1), (3, 5)] {
            for diag_offset in DIAG_OFFSETS {
                for (round, cell_count) in [7usize, 40, 90, 200, 13, 120].into_iter().enumerate() {
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: 0,
                    };
                    let s0 = (next(&mut seed) % 4000) as i32;
                    let entries = entries_near(&mut seed, size, s0, test);
                    let mut s = s0;
                    let cells: Vec<(i32, Vec<i32>)> = (0..cell_count)
                        .map(|k| {
                            let n = match next(&mut seed) % 8 {
                                0 => 25 + next(&mut seed) as usize % 96,
                                1 if round % 3 == 2 && k == cell_count / 2 => 333,
                                _ => next(&mut seed) as usize % 41,
                            };
                            sizes.insert(n);
                            let base = (next(&mut seed) % 3) as i32 * 7;
                            let cell = (s, distinct_offsets(&mut seed, size, n, base));
                            s += 1 + (next(&mut seed) % 4 / 3 * (next(&mut seed) % 6)) as i32;
                            cell
                        })
                        .collect();
                    let reference = reference_calls(&entries, mask, test, &cells);
                    for chunk in [X_CROSS_TEST_CHUNK, X_CROSS_CHUNK] {
                        let windows = cases.windows;
                        cases.add(&cells, mask, &reference, chunk);
                        assert_eq!(
                            cross_splits(&entries, mask, test, &cells, chunk),
                            (cases.windows - windows) as u64,
                            "chunk {chunk} {test:?}: windows tested in parts"
                        );
                    }
                    calls += same_in_every_mode(&entries, mask, test, &cells);
                }
            }
        }
        assert!(calls > 1000, "the cells must give candidates: {calls}");
        assert!((0..=40).all(|n| sizes.contains(&n)), "{sizes:?}");
        assert!(cases.windows > 0, "{cases:?}");
        assert!(cases.candidate_before_later_cell > 100, "{cases:?}");
        assert!(
            cases.long_cells > 100 && cases.longer_cells > 20,
            "{cases:?}"
        );
        assert!(cases.drains > 50 && cases.drains_in_a_cell > 0, "{cases:?}");
        assert!(cases.moved.iter().all(|&n| n > 0), "{cases:?}");
        assert!(cases.left.iter().all(|&n| n > 0), "{cases:?}");
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    // Cells of 1 to 5 hits at consecutive subject offsets, where about half the hits of a cell
    // continue a hit of the cell before on its diagonal (q + 1 at s + 1): eight staged hits then
    // often hold hits of several cells on one entry, and the window is tested in parts, cut before
    // each lane that repeats an entry of its part (parts of every length 1-7, parts that start past
    // the window's first lane; with word size 1 the second hit of a diagonal is a candidate, and its
    // call's write is what the third reads). Both chunks run (64 drains every few cells).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn cross_batches_with_hits_of_one_diagonal_in_consecutive_cells() {
        if !x_cross_available() {
            eprintln!("no AVX2 and POPCNT on this CPU: the cross-cell batches are not compared");
            return;
        }
        let mut seed = 0xd1a9_u64;
        let size = 2048usize;
        let mask = (size - 1) as i32;
        let mut cases = CrossCases::default();
        let mut calls = 0;
        for (window, word_size) in [(40, 1), (40, 3), (4, 1), (2, 1)] {
            for diag_offset in DIAG_OFFSETS {
                for _ in 0..6 {
                    let test = XDiagTest {
                        diag_offset,
                        window,
                        word_size,
                        prefetch: 0,
                    };
                    let s0 = (next(&mut seed) % 4000) as i32;
                    let entries = entries_near(&mut seed, size, s0, test);
                    let mut cells: Vec<(i32, Vec<i32>)> = Vec::new();
                    for k in 0..300 {
                        let n = 1 + next(&mut seed) as usize % 5;
                        let mut hits: Vec<i32> = Vec::new();
                        if let Some((_, before)) = cells.last() {
                            for &q in before {
                                if hits.len() < n && next(&mut seed) % 2 == 0 {
                                    // The next offset of the diagonal, at most the mask (as the
                                    // lookup's offsets are).
                                    hits.push((q + 1) & mask);
                                }
                            }
                        }
                        while hits.len() < n {
                            let q = (next(&mut seed) % size as u32) as i32;
                            if !hits.contains(&q) {
                                hits.push(q);
                            }
                        }
                        // The cell's order: the continued hits are not always first.
                        let rotate = next(&mut seed) as usize % n;
                        hits.rotate_left(rotate);
                        cells.push((s0 + k, hits));
                    }
                    let reference = reference_calls(&entries, mask, test, &cells);
                    for chunk in [X_CROSS_TEST_CHUNK, X_CROSS_CHUNK] {
                        let windows = cases.windows;
                        cases.add(&cells, mask, &reference, chunk);
                        assert_eq!(
                            cross_splits(&entries, mask, test, &cells, chunk),
                            (cases.windows - windows) as u64,
                            "chunk {chunk} {test:?}: windows tested in parts"
                        );
                    }
                    calls += same_in_every_mode(&entries, mask, test, &cells);
                }
            }
        }
        assert!(calls > 1000, "the cells must give candidates: {calls}");
        assert!(cases.windows > 1000, "{cases:?}");
        assert!(cases.cuts[1..].iter().all(|&n| n > 0), "{cases:?}");
        assert!(cases.late_parts > 100, "{cases:?}");
        assert!(cases.candidate_before_later_cell > 100, "{cases:?}");
        assert!(cases.drains > 100, "{cases:?}");
    }

    // No NCBI counterpart: LOSAT's stage; the checks of its precondition have no NCBI form.
    /// Runs [`x_scan_cells_cross`] on a one-protein lookup whose `query_length` is set to
    /// `query_length` over a table of `mask + 1` entries (the scan stops at its first assert).
    #[cfg(target_arch = "x86_64")]
    fn scan_cross_with_bounds(query_length: i32, entries: usize, mask: i32) {
        use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
        use crate::algorithm::tblastx::translation::QueryFrame;
        use crate::config::ScoringMatrix;
        use crate::stats::lookup_protein_params_ungapped;

        let seq = [1u8, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];
        let mut aa_seq = vec![0u8];
        aa_seq.extend_from_slice(&seq);
        aa_seq.push(0);
        let frames = vec![vec![QueryFrame {
            frame: 1,
            aa_seq,
            aa_seq_nomask: None,
            aa_len: seq.len(),
            orig_len: seq.len(),
            seg_masks: Vec::new(),
        }]];
        let params = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
        let (mut lookup, _) = build_ncbi_lookup(&frames, 11, &params, false);
        lookup.query_length = query_length;
        let mut entries = vec![DiagStruct::clear(40); entries];
        let mut table = XDiagTable::new(&mut entries, mask);
        let test = XDiagTest {
            diag_offset: 40,
            window: 40,
            word_size: 3,
            prefetch: 0,
        };
        let subject = [1u8, 3, 4, 5, 6, 7, 0];
        let mut candidate = |_: i32, _: i32, _: i32, _: i32, _: &mut DiagStruct| {};
        x_scan_cells_cross::<X_CROSS_TEST_CHUNK, _>(
            false,
            false,
            &lookup,
            &subject,
            subject.len() - 1,
            &mut table,
            test,
            &mut candidate,
        );
    }

    // The scan asserts, once per subject and before it stages a hit, that every query offset of
    // the lookup is at most the mask: here the lookup is longer than the table.
    #[cfg(target_arch = "x86_64")]
    #[test]
    #[should_panic(expected = "needs every query offset at most the diagonal mask")]
    fn cross_scan_asserts_that_the_lookup_fits_the_table() {
        if !x_cross_available() {
            // The scan would panic on the CPU check with another message: fail as expected.
            panic!(
                "needs every query offset at most the diagonal mask (no AVX2 and POPCNT: skipped)"
            );
        }
        scan_cross_with_bounds(64, 16, 15);
    }

    // The scan also asserts that `mask + 1` is a power of two: the recovery of the query offset
    // needs it (here the table has 21 entries).
    #[cfg(target_arch = "x86_64")]
    #[test]
    #[should_panic(expected = "needs a diagonal mask of the form 2^k - 1")]
    fn cross_scan_asserts_that_the_mask_is_a_power_of_two_less_one() {
        if !x_cross_available() {
            panic!("needs a diagonal mask of the form 2^k - 1 (no AVX2 and POPCNT: skipped)");
        }
        scan_cross_with_bounds(5, 21, 20);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:513-552 (quoted at `x_test_hit`)
    // LOSAT_X_BXCROSSDEFER's two rare paths, on tables made for them, against the per-hit test
    // (every mode, both chunks, [`same_in_every_mode`]), counted by the test probe of the deferred
    // arm: cells at subject offsets 5 apart on the same diagonals (q at s, q + 5 at s + 5), each
    // diagonal's entry first 5 behind, so a cell's hits are candidates and the next cell's batch
    // (cells of eight hits) or the next part of the window (cells of four: two cells per window,
    // on the same entries) meets the marker of their deferred calls; and entries that hold the
    // marker's value (flag 1, last_hit -1, as a call near offsets of 2^31 stores) under a drain's
    // first batch or window, nothing deferred, and under a batch after a cell of candidates (the
    // calls made, then the batch again, nothing deferred).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn deferred_calls_meet_their_marker_and_a_stored_marker_is_a_value() {
        if !x_cross_available() {
            eprintln!("no AVX2 and POPCNT on this CPU: the deferred calls are not compared");
            return;
        }
        let size = 1024usize;
        let mask = (size - 1) as i32;
        let test = XDiagTest {
            diag_offset: 40,
            window: 40,
            word_size: 3,
            prefetch: 0,
        };
        let s0 = 100;
        // `count` cells of `hits` hits on the diagonals 7 + 3 j (j below `hits`) from `s0`,
        // subject offsets `step` apart, and a table whose diagonals last hit 5 before `s0`.
        let cells_on = |count: i32, hits: i32, step: i32, s0: i32| -> Vec<(i32, Vec<i32>)> {
            (0..count)
                .map(|k| {
                    let s = s0 + step * k;
                    (s, (0..hits).map(|j| (7 + 3 * j + s) & mask).collect())
                })
                .collect()
        };
        let behind = |hits: i32, s0: i32| {
            let mut entries = vec![DiagStruct::default(); size];
            for j in 0..hits {
                let last = s0 - 5 + test.diag_offset;
                entries[(7 + 3 * j) as usize] = entry_of(last as u32 & X_LAST_HIT_MASK);
            }
            entries
        };
        let marker = entry_of(X_CROSS_MARKER as u32);
        // The candidate calls and the probe's [met with calls deferred, met with none].
        let run = |entries: &[DiagStruct], cells: &[(i32, Vec<i32>)]| {
            X_CROSS_MARKS_MET.set([0; 2]);
            let calls = same_in_every_mode(entries, mask, test, cells);
            (calls, X_CROSS_MARKS_MET.get())
        };
        for hits in [8, 4] {
            let (calls, met) = run(&behind(hits, s0), &cells_on(40, hits, 5, s0));
            assert!(calls > 40, "cells of {hits}: {calls} calls");
            assert!(met[0] > 10, "cells of {hits}: {met:?}");
            assert_eq!(met[1], 0, "cells of {hits}: {met:?}");
        }
        // A stored marker under the first batch or window of a drain, then under a batch after a
        // cell of candidates (on other diagonals: 8 + 3 j).
        for hits in [8, 4] {
            let mut entries = vec![DiagStruct::default(); size];
            entries[7] = marker;
            let (calls, met) = run(&entries, &cells_on(2, hits, 5, s0));
            assert!(
                met[1] >= 2 && met[0] == 0,
                "cells of {hits}: {calls} calls, {met:?}"
            );
        }
        let mut entries = behind(8, s0);
        entries[8] = marker;
        let mut cells = cells_on(1, 8, 5, s0);
        cells.push((
            s0 + 1,
            (0..8).map(|j| (8 + 3 * j + s0 + 1) & mask).collect(),
        ));
        let (calls, met) = run(&entries, &cells);
        assert!(calls >= 8, "{calls} calls");
        assert!(met[0] >= 2 && met[1] >= 2, "{met:?}");
    }

    // Every pattern of equal lanes of eight coordinates (the 4,140 partitions of the lanes, as
    // restricted growth strings): the conflict test sees whether two lanes are equal, and the part
    // from each lane ends at the first later lane whose coordinate a lane of the part holds.
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn the_conflict_test_and_the_parts_follow_the_equal_lanes() {
        if !x_cross_available() {
            eprintln!("no AVX2 and POPCNT on this CPU: the conflict test is not run");
            return;
        }
        let conflict = |coords: [i32; X_DIAG_LANES]| {
            // SAFETY: the CPU has AVX2 and POPCNT (checked above); the load reads the eight lanes
            // of `coords`.
            unsafe {
                x_cross_conflict(std::arch::x86_64::_mm256_loadu_si256(
                    coords.as_ptr().cast(),
                ))
            }
        };
        let mut classes = [0i32; X_DIAG_LANES];
        let mut patterns = 0;
        loop {
            let coords = classes.map(|class| class * 4_099 + 17);
            let repeated = (1..X_DIAG_LANES).any(|j| coords[..j].contains(&coords[j]));
            assert_eq!(conflict(coords), repeated, "{classes:?}");
            for start in 0..X_DIAG_LANES {
                // The part's classes so far; it ends at the first class it already holds.
                let mut held = [false; X_DIAG_LANES];
                let mut end = start;
                while end < X_DIAG_LANES && !held[classes[end] as usize] {
                    held[classes[end] as usize] = true;
                    end += 1;
                }
                assert_eq!(x_cross_part_end(&coords, start), end, "{classes:?} {start}");
            }
            patterns += 1;
            // The next restricted growth string: lane i takes a class up to one above the
            // largest of the lanes before it.
            let Some(i) = (1..X_DIAG_LANES)
                .rev()
                .find(|&i| classes[i] <= *classes[..i].iter().max().unwrap())
            else {
                break;
            };
            classes[i] += 1;
            classes[i + 1..].fill(0);
        }
        assert_eq!(patterns, 4_140);
        assert!(!conflict([3, 1, 2, 0, 7, 6, 5, 4096]));
        assert!(conflict([9; X_DIAG_LANES]));
    }

    // EXPERIMENT (LOSAT_X_BXCROSS): the query offset recovered from the entry coordinate and the
    // subject offset is the hit's for every query offset up to the mask, with subject offsets
    // below, at and past it, up to i32::MAX (the sum wraps).
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn the_query_offset_is_recovered_from_the_coordinate() {
        for mask in [0i32, 1, 7, 1_023, 65_535, (1 << 30) - 1] {
            let queries = [0, 1, mask / 3, mask / 2, mask - 1, mask];
            let past = mask.saturating_mul(2).saturating_add(3);
            let subjects = [0, 1, 5, mask, mask + 1, past, 5_000, i32::MAX - 1, i32::MAX];
            for q in queries.into_iter().filter(|q| (0..=mask).contains(q)) {
                for s in subjects {
                    let coord = q.wrapping_sub(s) & mask;
                    assert_eq!(x_cross_query(coord, s, mask), q, "q {q} s {s} mask {mask}");
                }
            }
        }
    }
}
