//! The blastp two-hit seed stage with the subject scan and the diagonal test in one loop.
//!
//! NCBI copies every word hit of the subject into the offset-pair buffer
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
//! The loop is the two-hit seed stage of every blastp search with a window (`-window_size` other
//! than 0); the one-hit search (`-window_size 0`) keeps the buffer of offset pairs
//! (`blastp_for_each_offset_pair`).
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
//! Without AVX2
//! (another architecture, or an x86-64 CPU without it) the switch keeps the per-hit loop, except
//! where LOSAT_X_BPCELL128's 128-bit batches run (below; WebAssembly, x86-64 with SSSE3, SSE4.1 and
//! POPCNT): the
//! scalar batch was slower than the per-hit loop on both machines it was timed on (2026-10-14,
//! human block 1 x Species4, one thread, scan_ungapped: workstation 12.56 s against 11.59 s,
//! local 7.53 / 7.08 s against 6.78 / 6.88 s). The scalar batch remains as the AVX2 path's handler
//! of a cell's last batch of fewer than eight hits and as the cell-by-cell version the tests
//! compare with; `LOSAT_X_BPDIAG_SIMD=0` forces it (diagnostic, used by the tests; not in the
//! strict set). With `LOSAT_TIMING=1` each search of a query set prints
//! `[TIMING] x_bpdiag: mode=<avx2|avx2-defer|avx2-defer-masked|w128|scalar|per-hit> cells=<n>
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
//! The masked batch: blastx's fused loop (`blastx/x_seed_diag.rs`), with a value of LOSAT_X_BXDIAG or
//! LOSAT_X_BXDIAGSHADOW where LOSAT_X_BXCROSS's cell test does not run, takes LOSAT_X_BPDEFER's deferred calls with two
//! changes of the batches
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
//! its own, never inlined, and a cell test of its own).
//!
//! EXPERIMENT (LOSAT_X_BPMASKED): blastp's fused loop takes the masked batch
//! ([`x_cell_avx2_defer_masked`], [`XDiagMode::Avx2DeferMasked`]) in place of LOSAT_X_BPDEFER's
//! deferred calls: the same buffer and calls, with the batches' entry coordinates from scalar loads
//! of the hits and a cell's last batch of fewer than eight hits as one masked AVX2 batch. Why the
//! outcome is LOSAT_X_BPDEFER's: the argument of the masked batch above, which rests only on what
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
//! run (`LOSAT_X_BPDEFER` with the AVX2 batch, [`x_diag_mode`]; not with
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
//! dropped (as in the masked batch), so they change nothing. The last hits take the
//! hit-by-hit order itself, and every batch ends (stores and calls) before the next one loads.
//! `LOSAT_X_BXCROSS` is a switch of the strict set (since 2026-10-23): on unless it is `0`; on,
//! it is the cell test of every table blastx's fused loop takes where the CPU has AVX2 and POPCNT
//! (`blastx/x_seed_diag.rs`), and elsewhere LOSAT_X_BXCROSS128's (below) where that runs; off, or
//! without either, blastx takes the direct loop, or with a value of LOSAT_X_BXDIAG or
//! LOSAT_X_BXDIAGSHADOW the masked batch (the per-hit loop without AVX2). blastp's loops never take it: it is not a cell test of
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
//! EXPERIMENT (LOSAT_X_BXCROSS128, blastx's switch): LOSAT_X_BXCROSS's batches in 128-bit vectors,
//! for the targets without its AVX2 kernel: the WebAssembly build (simd128) and x86-64 with SSSE3,
//! SSE4.1 and POPCNT but without AVX2 ([`XCellTest::Cross128`], [`x_scan_cells_w128_cross`]). Both
//! kernels run one driver, generic over [`XCrossKernel`]: the scan and its stage (`x_cross_scan!`,
//! [`x_cross_cell`]), the drain ([`x_cross_drain_with`]) and its batches ([`x_cross_batch`],
//! [`x_batch_lanes`]); the deferred records ([`x_cross_defer`], [`x_cross_batch_marked`]) are the
//! AVX2 kernel's only. A kernel gives a piece's stage, the conflict test, the test of eight lanes
//! and its window with a repeated entry, and each kernel's functions are compiled for its
//! instructions (the trait's notes say why the AVX2 functions keep their instructions). The
//! 128-bit kernel ([`XCrossW128`], over `utils/xdrop_simd.rs` `arch`, one source for x86-64 SSE
//! and wasm32 simd128) holds the eight
//! lanes of a window in two vectors of four. Its conflict test compares each half with itself one
//! and two lanes round (the six pairs within it) and the low half with the high one zero to three
//! lanes round (the sixteen pairs across), so it sees all 28 pairs of lanes, as the AVX2 test does.
//! Its batch reads and stores the entries by lane, as the AVX2 batch does, and computes the AVX2
//! batch's lane formulas (the same wrapping 32-bit operations) on each half. A piece of a cell is
//! staged as six vectors of four lanes (`XHit::x_load24_w128`: `u16` hits by three 16-byte loads,
//! the reads of the AVX2 stage; `i32` hits by a 16-byte load for each four hits inside the list and
//! scalar loads for the rest, as SSE and simd128 have no masked load). A window with a repeated entry
//! is tested hit by hit ([`x_cross_window_w128`]: [`x_cross_hit`], the reference's order) instead of
//! in parts: at the fewest hits per used cell the fused loop takes a fifth of the windows repeat an
//! entry, and the parts' batches cost the 128-bit kernel more than the eight scalar tests. So the
//! staged stream, the batches of distinct entries and their calls are the AVX2 kernel's, the
//! windows with a repeated entry take the per-hit order itself, and the argument above holds as it
//! stands. The 128-bit kernel calls each batch's candidates at once (the drain without `DEFER`,
//! which its `drain` refuses when it is compiled; LOSAT_X_BXCROSSDEFER does not act on it): in the
//! WebAssembly build the deferral was slower (the prototype build bx1 of v4-W W2 s3, whose windows
//! with a repeated entry were still tested in parts, on a 50-kb piece of the fly genome x the fly
//! proteome at one thread: x0.780 of `dev` with it, x0.749 without), while on x86-64 with `=2` it
//! was 3 % faster, which a CPU without AVX2 forgoes so that the kernel has one form.
//! `LOSAT_X_BXCROSS128` is a switch of the strict set (class platform): on unless it is `0`; on,
//! where LOSAT_X_BXCROSS is on and its AVX2 kernel does not run, the 128-bit kernel is the cell test
//! of every table blastx's fused loop takes (`blastx/x_seed_diag.rs`; LOSAT_X_BXDIAGMIN chooses the
//! tables, as for the AVX2 kernel); `=2` (diagnostic) runs it in place of the AVX2 kernel on any
//! x86-64 CPU with SSSE3, SSE4.1 and POPCNT (the tests); `=0`, or LOSAT_X_BXCROSS=0, leaves the
//! choice before the switch (the direct loop, or with a value of LOSAT_X_BXDIAG the per-hit loop).
//! The timing line names `mode=w128-cross` and counts `cross_splits=` (the windows with a repeated
//! entry, which this kernel tests hit by hit) as for the AVX2 kernel.
//!
//! EXPERIMENT (LOSAT_X_BPCELL128): LOSAT_X_BPDIAG's batches of eight in 128-bit vectors, for the
//! targets without AVX2: the WebAssembly build (simd128) and x86-64 with SSSE3, SSE4.1 and POPCNT
//! but without AVX2, which ran the per-hit loop on every table before ([`XDiagMode::W128`],
//! [`x_scan_cells_w128`], [`x_cell_w128`]). A cell's hits go in batches of eight, each the test of
//! eight lanes of LOSAT_X_BXCROSS128's kernel ([`x_batch_lanes`] with [`XCrossW128`], every lane at
//! the cell's subject offset, so its lane formulas are the AVX2 batch's, [`x_batch_avx2_at`]); a
//! batch's candidates are called at once, in lane order ([`x_call_candidates`], as [`x_cell_avx2`]
//! calls them); the fewer than eight hits after the last batch are tested one after the other,
//! with LOSAT_X_BPHITBF's branch-free test where the per-hit loop takes it ([`XHitBf::new`] once
//! per subject, with LOSAT_X_BPHITBF on) and with [`x_test_hit`] elsewhere. Why the outcome is the
//! per-hit loop's: LOSAT_X_BPDIAG's argument above for the batches (the hits of one cell lie on
//! distinct entries, a candidate's entry keeps its value in the test, the calls keep the cell's
//! order), and the last hits are the per-hit loop's own tests (LOSAT_X_BPHITBF's argument below)
//! on the same entries in the same order. LOSAT_X_BPDEFER's deferral and the masked last batch are
//! not taken: in the WebAssembly build both cost more than they save (v4-W W2 s4, instructions of
//! the compiled module under callgrind, V8 TurboFan, one thread, 1-100 queries of Species4 x
//! Species6: at 118 hits per used cell the batches with deferred calls 0.887 of the per-hit loop's,
//! with calls at once 0.831, at 3.10 the masked last batch 1.059 against 1.004 for the last hits
//! one after the other; prototype builds). A table with fewer hits per used cell than LOSAT_X_BPMASKEDMIN (250
//! hundredths, the masked batch's value) keeps the per-hit loop ([`x_bp_table_mode`]): most of its
//! cells hold one or two hits, which the batches' loop only adds to (with every table batched
//! 1.012-1.024 of the per-hit loop's instructions at 1.25-2.42 hits per used cell; the final build
//! at 2.52 and 3.10 hits 1.002 and 1.013, 3.72 0.968, 6.37 0.940, 14.7 0.889, 118 0.842, 122
//! 0.821).
//! `LOSAT_X_BPCELL128` is a switch of the strict set (class platform): on unless it is `0`; on,
//! where LOSAT_X_BPDIAG is on and the AVX2 batches do not run, the 128-bit batches are the cell
//! test of the tables LOSAT_X_BPMASKEDMIN gives them; `=2` (a diagnostic) runs them in place of the
//! AVX2 batches on any x86-64 CPU with SSSE3, SSE4.1 and POPCNT (the tests); `=0` leaves the
//! per-hit loop. The timing line names `mode=w128` and prints LOSAT_X_BPMASKEDMIN's line of tables
//! (`tables_masked=` the tables batched, `tables_deferred=` those left to the per-hit loop).
//!
//! BXSCANBF (graduated 2026-10-26: blastx always takes it): `x_scan_two_hit_with(.., branch_free: true, ..)`
//! finds the cells with hits by the candidate loop of tblastn's branch-free scan
//! (`tblastx/x_aascan_bf.rs`, `x_block_candidates`): the positions of a block of words are written
//! to a list whose length advances by the presence bit, without a branch on it, and the listed
//! cells then go to the cell test in position order. The cells, their order and their hit lists
//! are the presence-test loop's, so every outcome above is unchanged. blastp's own loop
//! (`x_scan_two_hit`) takes the same cells ([`x_scan_two_hit_bf`]); blastp-fast's loop on the flat compressed table lists its cells
//! the same way under LOSAT_X_BPCLUTSCANBF (`x_seed_fused_clut.rs`).
//!
//! EXPERIMENT (LOSAT_X_BPHITBF): the per-hit loop (`x_cell_per_hit`: the fused loop without
//! LOSAT_X_BPDIAG's batches, the default where the CPU has no AVX2, i.e. aarch64, and where
//! LOSAT_X_BPCELL128's 128-bit batches do not take the table, i.e. WebAssembly's tables with few
//! hits per used cell; LOSAT_X_BPCELL128 takes the same test for a cell's last hits)
//! takes the two outcomes most two-hit tests end in, too far and overlap, without a branch (tblastn's
//! tblastn's branch-free two-hit loop and TBLASTX's LOSAT_X_TBXHITBF do the same in their loops, on tables whose entries
//! need no sign extension):
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:538-551
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
//! length)`: the cleared value (`DiagStruct::clear`, when the offset passes `INT4_MAX` divided by 4), 0 in a
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
//! compressed table (`x_seed_fused_clut.rs`, `x_scan_flat_fseed_loop`: the same diagonal table and
//! extension as blastp's loop).
//!
//! The switch is read once; it is in the strict set (`utils/xswitch.rs`): on unless
//! `LOSAT_X_BPHITBF=0`. The loop with the test is its own function, never inlined and with a
//! closure of its own type (`x_scan_two_hit_hitbf`), so the loops without it keep their code; the
//! AVX2 and scalar batches keep their own tests, and so does blastp-fast's loop on the reference
//! compressed table (`x_seed_fused_compressed.rs`). blastp-fast's loop on the flat compressed table
//! takes the test, without the read-ahead, wherever the target is not wasm32 and the table takes
//! the branch-free list (LOSAT_X_BPCLUTSCANBF; LOSAT_X_BPHITBF does not act there).
//! `x_scan_two_hit_hitbf` also reads ahead without a bound check per hit (`x_cell_per_hit_bf`).
//! With `LOSAT_X_STATS` (and the `xstats` feature) `main` prints `[X_STATS] BPHITBF subjects=N` at
//! exit when N is not zero, N the subjects the loop with the test scanned.

use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::{
    BlastAaLookupTable, XCell, XHit, PV_ARRAY_BTS, PV_ARRAY_MASK,
};
use crate::algorithm::tblastx::x_aascan_bf::{x_block_candidates, X_BLOCK};
use crate::utils::matrix::BLASTAA_SIZE;
use std::mem::MaybeUninit;

// No NCBI counterpart: how many hits ahead the fused loop reads the diagonal entries (0 = no
// prefetch).
pub(super) const X_PREFETCH_DISTANCE: usize = 16;

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
#[cfg(target_arch = "x86_64")]
pub(crate) fn x_bp_defer() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPDEFER"))
}

// No NCBI counterpart: reads the LOSAT_X_BPMASKED switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPMASKED): a strict-set switch, on unless it is `0`.
/// It acts on LOSAT_X_BPDEFER's deferred calls only.
#[cfg(target_arch = "x86_64")]
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
/// The same value chooses the tables of LOSAT_X_BPCELL128's 128-bit batches; a table with fewer
/// keeps the per-hit loop, which the batches' cells of one or two hits only add to (module
/// header). `0` takes the batches for every table; a value that is not a number takes the
/// built-in one.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const X_BP_MASKED_MIN: u64 = 250;

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
static X_BP_TABLES: [std::sync::atomic::AtomicU64; 2] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
static X_BP_RATIO_MIN: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(u64::MAX);
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
static X_BP_RATIO_MAX: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

// No NCBI counterpart: a property of the lookup table that chooses LOSAT's cell test; it does
// not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): `table`'s hits per used backbone cell, in hundredths (0 for
/// a table without hits; LOSAT_X_BXDIAGMIN's measure, `blastx/x_seed_diag.rs`), worked out once
/// per table and thread (the table's addresses and overflow length identify it; a later table
/// that matched all three would reuse the value, which changes only the speed), and whether this
/// call worked it out.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_bp_table_ratio<K: XCell>(table: &BlastAaLookupTable<K>) -> (u64, bool) {
    let key = (
        table as *const BlastAaLookupTable<K> as usize,
        table.backbone.as_ptr() as usize,
        table.overflow.len(),
    );
    X_BP_CHOICE.with(|choice| {
        let (a, b, c, ratio) = choice.get();
        if (a, b, c) == key {
            return (ratio, false);
        }
        let (hits, used) = K::x_hits_and_used(&table.backbone);
        let ratio = (hits * 100).checked_div(used).unwrap_or(0);
        choice.set((key.0, key.1, key.2, ratio));
        (ratio, true)
    })
}

// No NCBI counterpart: chooses LOSAT's cell test; every cell test gives the reference's hits and
// table.
/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): the cell test for the subjects of `table`: `mode`, but
/// where the table has fewer hits per used cell than [`x_bp_masked_min`] LOSAT_X_BPDEFER's
/// deferred calls instead of LOSAT_X_BPMASKED's masked batch, and the per-hit loop instead of
/// LOSAT_X_BPCELL128's 128-bit batches. Counted for the timing line once per table and thread,
/// and the masked batch's choice marked for `LOSAT_X_REACH`.
fn x_bp_table_mode<K: XCell>(mode: XDiagMode, table: &BlastAaLookupTable<K>) -> XDiagMode {
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
    // EXPERIMENT (LOSAT_X_BPCELL128): the 128-bit batches on the tables the masked batch would
    // take; a table with fewer hits per used cell keeps the per-hit loop (module header).
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    if mode == XDiagMode::W128 {
        let (ratio, first) = x_bp_table_ratio(table);
        let batched = ratio >= x_bp_masked_min();
        if first {
            use std::sync::atomic::Ordering::Relaxed;
            X_BP_TABLES[usize::from(!batched)].fetch_add(1, Relaxed);
            X_BP_RATIO_MIN.fetch_min(ratio, Relaxed);
            X_BP_RATIO_MAX.fetch_max(ratio, Relaxed);
        }
        if !batched {
            return XDiagMode::PerHit;
        }
    }
    #[cfg(not(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )))]
    let _ = table;
    mode
}

// No NCBI counterpart: whether LOSAT's cell test is chosen per table; it computes no value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPMASKEDMIN): whether [`x_bp_table_mode`] chooses the cell test of `mode`
/// per table: LOSAT_X_BPMASKED's masked batch and LOSAT_X_BPCELL128's 128-bit batches.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_bp_ruled(mode: XDiagMode) -> bool {
    #[cfg(target_arch = "x86_64")]
    if mode == XDiagMode::Avx2DeferMasked {
        return true;
    }
    mode == XDiagMode::W128
}

// No NCBI counterpart: reads the LOSAT_X_BPCELL128 switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPCELL128): 0 = off, 1 = the 128-bit batches where the AVX2 batches do not
/// run, 2 = those in place of the AVX2 batches too (a diagnostic, for the tests); a switch of the
/// strict set (class platform), on unless it is `0`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_bp_cell128() -> usize {
    use std::sync::OnceLock;
    static ON: OnceLock<usize> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BPCELL128", 1..=2, 1).unwrap_or(0))
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
/// blastp-fast's loop on the flat compressed table (`x_seed_fused_clut.rs`) use it.
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551 (cited at `x_test_hit`)
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
    /// LOSAT_X_BPDIAG with AVX2 and LOSAT_X_BPDEFER off: batches of eight hits in AVX2 registers ([`x_test_batch_avx2`]); a
    /// cell's last batch of fewer than eight hits is a scalar batch.
    #[cfg(target_arch = "x86_64")]
    Avx2,
    /// LOSAT_X_BPDIAG with AVX2 and LOSAT_X_BPDEFER, with LOSAT_X_BPMASKED off or on a table below
    /// LOSAT_X_BPMASKEDMIN ([`x_bp_table_mode`]): the AVX2 batches with the candidates of a cell's
    /// batches of eight called after its last one ([`x_cell_avx2_defer`]).
    #[cfg(target_arch = "x86_64")]
    Avx2Defer,
    /// LOSAT_X_BPDIAG with AVX2, LOSAT_X_BPDEFER and LOSAT_X_BPMASKED, on tables with at least
    /// LOSAT_X_BPMASKEDMIN hits per used cell (and blastx's fused loop, [`x_cell_test`], with a
    /// value of LOSAT_X_BXDIAG or LOSAT_X_BXDIAGSHADOW and LOSAT_X_BXCROSS off or no POPCNT):
    /// [`XDiagMode::Avx2Defer`] with
    /// the lanes' entry coordinates from scalar loads of the hits and a cell's last batch of fewer
    /// than eight hits as one masked AVX2 batch ([`x_cell_avx2_defer_masked`]).
    #[cfg(target_arch = "x86_64")]
    Avx2DeferMasked,
    /// LOSAT_X_BPCELL128: batches of eight in two 128-bit vectors of four lanes, their calls at
    /// once, and a cell's last hits one after the other ([`x_cell_w128`]), where the AVX2 batches
    /// do not run; LOSAT_X_BPMASKEDMIN's tables (fewer hits per used cell keep the per-hit loop).
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    W128,
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
            #[cfg(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            ))]
            XDiagMode::W128 => "w128",
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
    /// LOSAT_X_BXCROSS128: [`XCellTest::Cross`]'s batches in two 128-bit vectors of four lanes
    /// ([`x_scan_cells_w128_cross`]; wasm32 simd128, x86-64 with SSSE3, SSE4.1 and POPCNT), where the
    /// AVX2 batches do not run, with each batch's calls at once (LOSAT_X_BXCROSSDEFER does not act
    /// here: its deferral made the WebAssembly build slower, module header) and a window with a
    /// repeated entry tested hit by hit ([`x_cross_window_w128`]).
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    Cross128,
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
            #[cfg(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            ))]
            XCellTest::Cross128 => "w128-cross",
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

// No NCBI counterpart: whether the CPU has what LOSAT_X_BXCROSS128's (and LOSAT_X_BPCELL128's)
// functions are compiled for; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSS128, LOSAT_X_BPCELL128): SSSE3, SSE4.1 and POPCNT on x86-64 (the
/// 128-bit layer of `utils/xdrop_simd.rs`, and `popcnt` as for LOSAT_X_BXCROSS); always on wasm32
/// simd128.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_cross128_available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("ssse3")
            && std::arch::is_x86_feature_detected!("sse4.1")
            && std::arch::is_x86_feature_detected!("popcnt")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        true
    }
}

// No NCBI counterpart: chooses LOSAT's code path once (LOSAT_X_BPDIAG; the AVX2 batch where the
// CPU has AVX2, with LOSAT_X_BPDEFER the deferred calls, with LOSAT_X_BPMASKED also the masked
// batch; without AVX2 LOSAT_X_BPCELL128's 128-bit batches where the target has them; otherwise
// the per-hit loop; LOSAT_X_BPDIAG_SIMD=0 forces the scalar batch); it does not change any value
// NCBI computes.
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
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        // EXPERIMENT (LOSAT_X_BPCELL128, `=2`, a diagnostic): the 128-bit batches in place of the
        // AVX2 ones.
        if x_bp_cell128() == 2 && x_cross128_available() {
            return XDiagMode::W128;
        }
        #[cfg(target_arch = "x86_64")]
        if x_avx2_available() {
            return match (x_bp_defer(), x_bp_masked()) {
                (true, true) => XDiagMode::Avx2DeferMasked,
                (true, false) => XDiagMode::Avx2Defer,
                (false, _) => XDiagMode::Avx2,
            };
        }
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        // EXPERIMENT (LOSAT_X_BPCELL128): without AVX2, the 128-bit batches where the target has
        // them (wasm32 simd128; x86-64 with SSSE3, SSE4.1 and POPCNT).
        if x_bp_cell128() != 0 && x_cross128_available() {
            return XDiagMode::W128;
        }
        // Elsewhere the scalar batch would be slower than the per-hit loop (module header).
        XDiagMode::PerHit
    })
}

// No NCBI counterpart: chooses LOSAT's code path once per call site (LOSAT_X_BXDIAG,
// LOSAT_X_BXCROSS, LOSAT_X_BXCROSSDEFER); it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAG): the cell test of blastx's fused seed loop
/// (`blastx/x_seed_diag.rs`): `"per-hit"` or `"scalar"` when asked by name (diagnostics), otherwise
/// the AVX2 batch where the CPU has AVX2, else the per-hit loop, as LOSAT_X_BPDIAG chooses for
/// blastp (the scalar batch was slower than the per-hit loop; module header). The AVX2 batch is
/// [`XCellTest::Cross`] with `cross` (LOSAT_X_BXCROSS) where the CPU also has POPCNT (its calls
/// deferred with `defer`, LOSAT_X_BXCROSSDEFER), otherwise [`XDiagMode::Avx2DeferMasked`].
///
/// EXPERIMENT (LOSAT_X_BXCROSS128): with `cross` and `cross128` 1 (the switch's value), where
/// [`XCellTest::Cross`] does not run (no AVX2 or no POPCNT, or another architecture) and the CPU has
/// what [`x_cross128_available`] asks, the cell test is [`XCellTest::Cross128`] (which does not
/// defer its calls: `defer` acts on [`XCellTest::Cross`] only); with `cross128` 2 (diagnostic) it
/// is that in place of [`XCellTest::Cross`] too; with 0 the choice is the one before the switch.
/// The choice for each CPU is [`x_cell_test_on`]'s, with this CPU's instructions.
pub(crate) fn x_cell_test(requested: &str, cross: bool, defer: bool, cross128: usize) -> XCellTest {
    x_cell_test_on(requested, cross, defer, cross128, XCellCpu::this())
}

// No NCBI counterpart: the instructions LOSAT's cell tests are compiled for; it does not change any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSS128): what [`x_cell_test_on`] asks of the CPU, so that its choice
/// for any CPU can be tested on any machine.
#[derive(Clone, Copy, Debug)]
struct XCellCpu {
    /// AVX2 (the masked AVX2 batch; x86-64 only).
    avx2: bool,
    /// AVX2 and POPCNT ([`x_cross_available`]: LOSAT_X_BXCROSS's AVX2 kernel; x86-64 only).
    cross: bool,
    /// [`x_cross128_available`]: LOSAT_X_BXCROSS128's kernel (x86-64 SSSE3, SSE4.1 and POPCNT;
    /// wasm32 simd128).
    cross128: bool,
}

impl XCellCpu {
    /// This CPU's.
    fn this() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            XCellCpu {
                avx2: x_avx2_available(),
                cross: x_cross_available(),
                cross128: x_cross128_available(),
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            XCellCpu {
                avx2: false,
                cross: false,
                #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
                cross128: x_cross128_available(),
                #[cfg(not(all(target_arch = "wasm32", target_feature = "simd128")))]
                cross128: false,
            }
        }
    }
}

// No NCBI counterpart: chooses LOSAT's code path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXDIAG, LOSAT_X_BXCROSS128): [`x_cell_test`] on a CPU with the instructions
/// `cpu`. The AVX2 cell tests exist on x86-64 only and the 128-bit one on x86-64 and wasm32
/// simd128, so elsewhere `cpu` asks for none.
fn x_cell_test_on(
    requested: &str,
    cross: bool,
    defer: bool,
    cross128: usize,
    cpu: XCellCpu,
) -> XCellTest {
    XCellTest::Diag(match requested {
        "per-hit" => XDiagMode::PerHit,
        "scalar" => XDiagMode::Scalar,
        _ => {
            #[cfg(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            ))]
            if cross && cross128 == 2 && cpu.cross128 {
                return XCellTest::Cross128;
            }
            #[cfg(target_arch = "x86_64")]
            if cross && cpu.cross {
                return XCellTest::Cross { defer };
            }
            #[cfg(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            ))]
            if cross && cross128 != 0 && cpu.cross128 {
                return XCellTest::Cross128;
            }
            #[cfg(target_arch = "x86_64")]
            if cpu.avx2 {
                return XCellTest::Diag(XDiagMode::Avx2DeferMasked);
            }
            #[cfg(not(target_arch = "x86_64"))]
            let _ = (defer, cpu.avx2, cpu.cross);
            #[cfg(not(any(
                target_arch = "x86_64",
                all(target_arch = "wasm32", target_feature = "simd128")
            )))]
            let _ = (cross, cross128, cpu.cross128);
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
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
    #[inline(always)]
    fn add_staged(&mut self, hits: usize) {
        self.cells += 1;
        self.hits += hits as u64;
    }

    /// EXPERIMENT (LOSAT_X_BXCROSS): the candidates of the staged hits tested in one drain.
    #[cfg(any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    ))]
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
/// drain) and the windows of eight staged hits with a repeated entry (tested in parts by the AVX2
/// kernel, [`x_cross_window`], and hit by hit by the 128-bit one, [`x_cross_window_w128`]; added
/// once per drain). blastx adds them to the totals of its timing line
/// once per subject (`blastx/x_seed_diag.rs`), as it does the other cell tests' counts.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
        // EXPERIMENT (LOSAT_X_BPMASKEDMIN): the tables that took the masked batch (or
        // LOSAT_X_BPCELL128's 128-bit batches) and those its rule gave the deferred calls (or the
        // per-hit loop), once per thread and table.
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        if x_bp_ruled(x_diag_mode()) {
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516
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
    /// [`X_PREFETCH_DISTANCE`] (the per-hit loop only).
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

// The two-hit diagonal test of one hit (`common::diag_entry`, shared with the one word of a
// hard chunk shorter than the word).
pub(super) use crate::algorithm::common::diag_entry::x_test_hit;

/// Scans `subject` with `lookup` (two-hit mode) and applies the diagonal test to every word hit
/// in scan order; `candidate(query_offset, subject_offset, last_hit, diff, entry)` receives each
/// hit that needs the rest of the test, with its diagonal entry.
///
/// `subject` starts at the first residue and keeps the sentinel after the last one, as the
/// reference scan reads it; `subject_length` is the number of residues.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-512
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit<K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
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
        prefetch: X_PREFETCH_DISTANCE,
    };
    let table = &mut XDiagTable::new(diag_array, diag_mask);
    // EXPERIMENT (LOSAT_X_BPMASKEDMIN): the masked batch (or LOSAT_X_BPCELL128's 128-bit batches)
    // only on tables with enough hits per cell; chosen before the per-hit loop's test, which the
    // tables it gives the per-hit loop take too.
    let mode = x_bp_table_mode(x_diag_mode(), lookup);
    // EXPERIMENT (LOSAT_X_BPHITBF): the per-hit loop with the branch-free two-hit test, where it
    // gives the reference's outcomes (module header).
    if let Some(bf) = x_hit_bf_for(
        mode == XDiagMode::PerHit,
        diag_offset,
        subject_length,
        window,
        word_size,
    ) {
        x_scan_two_hit_hitbf(
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
    // The cells are found without a branch on the presence bit ([`x_scan_cells_bf`]).
    let counts = x_scan_two_hit_bf(
        mode,
        lookup,
        subject,
        subject_length,
        table,
        test,
        &mut candidate,
    );
    counts.add_to_totals();
}

/// [`x_scan_two_hit_in_with`] with the cells found by [`x_scan_cells_bf`]. Out of line, and with
/// a closure of its own type, so that its cell tests keep a single caller (when the closure type
/// was shared with another loop, the cell test had two callers and that instance ran 2-6 % more
/// instructions).
#[inline(never)]
fn x_scan_two_hit_bf<K: XCell, F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable<K>,
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
    x_scan_two_hit_in_with::<true, _, _>(
        mode,
        lookup,
        subject,
        subject_length,
        table,
        test,
        &mut own,
    )
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
/// two-hit test ([`x_cell_per_hit_bf`]), the cells found by [`x_scan_cells_bf`]. Out of line, and
/// with a closure of its own type, so that the loops without the test keep their code (see
/// [`x_scan_two_hit_bf`]).
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_two_hit_hitbf<K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
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
    x_scan_cells_with::<true, _, _>(lookup, subject, subject_length, |s, hits| {
        x_cell_per_hit_bf(hits, s, table, test, bf, &mut own)
    })
}

/// EXPERIMENT (LOSAT_X_BXDIAG): [`x_scan_two_hit`] with the cell test `mode` (blastx's choice,
/// [`x_cell_test`], not LOSAT_X_BPDIAG's); returns the cells tested in batches, their hits, the
/// hits passed to `candidate` (zeros for [`XDiagMode::PerHit`]) and the windows of eight staged
/// hits with a repeated entry of LOSAT_X_BXCROSS's batches (in parts) and LOSAT_X_BXCROSS128's
/// (hit by hit; zero for the other cell tests). Nothing is
/// added to the totals of the `[TIMING] x_bpdiag` line. Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx):
/// `branch_free` finds the cells with [`x_scan_cells_bf`] (module header). EXPERIMENT
/// (LOSAT_X_BXCROSS): [`XCellTest::Cross`] runs here ([`x_scan_cells_cross`]); the other cell
/// tests run through the choice that blastp's loops share ([`x_scan_two_hit_in_with`]).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-512 (cited at `x_scan_two_hit`)
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_scan_two_hit_with<K: XCell, F>(
    mode: XCellTest,
    branch_free: bool,
    lookup: &BlastAaLookupTable<K>,
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
        prefetch: X_PREFETCH_DISTANCE,
    };
    let table = &mut XDiagTable::new(diag_array, diag_mask);
    #[cfg_attr(
        not(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        )),
        expect(
            clippy::infallible_destructuring_match,
            reason = "the `Cross` and `Cross128` variants exist on x86_64 and wasm32 simd128 only, so elsewhere `Diag` is the only pattern"
        )
    )]
    let mode = match mode {
        XCellTest::Diag(mode) => mode,
        #[cfg(target_arch = "x86_64")]
        XCellTest::Cross { defer } => {
            let counts = x_scan_cells_cross::<X_CROSS_CHUNK, K, F>(
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
        // EXPERIMENT (LOSAT_X_BXCROSS128): the same batches in 128-bit vectors.
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        XCellTest::Cross128 => {
            let counts = x_scan_cells_cross128::<X_CROSS_CHUNK, K, F>(
                branch_free,
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
        x_scan_two_hit_in_with::<true, _, F>(
            mode,
            lookup,
            subject,
            subject_length,
            table,
            test,
            &mut candidate,
        )
    } else {
        x_scan_two_hit_in_with::<false, _, F>(
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

/// [`x_scan_two_hit`] with the cell test of `mode`; returns the counts of LOSAT_X_BPDIAG (zero for
/// [`XDiagMode::PerHit`]). The cells are found by [`x_scan_cells_bf`] when `BRANCH_FREE`
/// (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx).
fn x_scan_two_hit_in_with<const BRANCH_FREE: bool, K: XCell, F>(
    mode: XDiagMode,
    lookup: &BlastAaLookupTable<K>,
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
            x_scan_cells_with::<BRANCH_FREE, _, _>(lookup, subject, subject_length, |s, hits| {
                x_cell_per_hit(hits, s, table, test, candidate)
            })
        }
        XDiagMode::Scalar => {
            x_scan_cells_with::<BRANCH_FREE, _, _>(lookup, subject, subject_length, |s, hits| {
                let candidates = x_cell_scalar(hits, s, table, test, candidate);
                counts.add_cell(hits.len(), candidates);
            })
        }
        #[cfg(target_arch = "x86_64")]
        XDiagMode::Avx2 => {
            assert!(x_avx2_available(), "the AVX2 batch needs AVX2");
            // SAFETY: the CPU has AVX2 (asserted above).
            unsafe {
                x_scan_cells_avx2::<BRANCH_FREE, K, F>(
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
                x_scan_cells_avx2_defer::<BRANCH_FREE, K, F>(
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
                x_scan_cells_avx2_defer_masked::<BRANCH_FREE, K, F>(
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
        #[cfg(any(
            target_arch = "x86_64",
            all(target_arch = "wasm32", target_feature = "simd128")
        ))]
        XDiagMode::W128 => {
            assert!(
                x_cross128_available(),
                "the 128-bit batches need SSSE3, SSE4.1 and POPCNT"
            );
            // SAFETY: on x86-64 the CPU has SSSE3, SSE4.1 and POPCNT (asserted above).
            unsafe {
                x_scan_cells_w128::<BRANCH_FREE, K, F>(
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
fn x_scan_cells_cross<const CHUNK: usize, K: XCell, F>(
    branch_free: bool,
    defer: bool,
    lookup: &BlastAaLookupTable<K>,
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
    x_cross_check_bounds(lookup.query_length, table.mask);
    let mut counts = XCrossCounts::default();
    // SAFETY: the CPU has AVX2 and POPCNT, and every query offset of the lookup is below
    // `query_length`, so at most the mask (asserted above).
    unsafe {
        match (branch_free, defer) {
            (true, true) => x_scan_cells_avx2_cross::<true, true, CHUNK, K, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
            (true, false) => x_scan_cells_avx2_cross::<true, false, CHUNK, K, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
            (false, true) => x_scan_cells_avx2_cross::<false, true, CHUNK, K, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            ),
            (false, false) => x_scan_cells_avx2_cross::<false, false, CHUNK, K, F>(
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

// No NCBI counterpart: a check of the bound LOSAT_X_BXCROSS's stage rests on.
/// EXPERIMENT (LOSAT_X_BXCROSS): the checks of [`x_scan_cells_cross`] (and of
/// [`x_scan_cells_cross128`]) that every query offset of a lookup of `query_length` letters is at
/// most `mask`, and that `mask + 1` is a power of two, once per subject scan.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
fn x_cross_check_bounds(query_length: i32, mask: i32) {
    assert!(
        (0..=i64::from(mask) + 1).contains(&i64::from(query_length)),
        "the cross-cell stage needs every query offset at most the diagonal mask"
    );
    // `mask + 1` a power of two (written in i64, so that `mask + 1` cannot overflow): the
    // recovery `(coord + s) & mask` of `x_cross_query` equals `q` only then (mask 5, q 2, s 3
    // stages coord (2 - 3) & 5 = 5 and recovers (5 + 3) & 5 = 0).
    assert!(
        i64::from(mask) >= 0 && i64::from(mask) & (i64::from(mask) + 1) == 0,
        "the cross-cell stage needs a diagonal mask of the form 2^k - 1"
    );
}

/// EXPERIMENT (LOSAT_X_BXCROSS128): [`x_scan_cells_cross`] with the 128-bit kernel
/// ([`x_scan_cells_w128_cross`], its calls at once): the same checks, the instances of
/// `branch_free`. Out of line and called by blastx's [`x_scan_two_hit_with`] only, as
/// [`x_scan_cells_cross`].
// No NCBI counterpart: chooses LOSAT's cell test; it does not change any value NCBI computes.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn x_scan_cells_cross128<const CHUNK: usize, K: XCell, F>(
    branch_free: bool,
    lookup: &BlastAaLookupTable<K>,
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
        x_cross128_available(),
        "the cross-cell 128-bit batch needs SSSE3, SSE4.1 and POPCNT"
    );
    x_cross_check_bounds(lookup.query_length, table.mask);
    let mut counts = XCrossCounts::default();
    // SAFETY: on x86-64 the CPU has SSSE3, SSE4.1 and POPCNT, and every query offset of the
    // lookup is below `query_length`, so at most the mask (asserted above).
    unsafe {
        if branch_free {
            x_scan_cells_w128_cross::<true, CHUNK, K, F>(
                lookup,
                subject,
                subject_length,
                table,
                test,
                candidate,
                &mut counts,
            )
        } else {
            x_scan_cells_w128_cross::<false, CHUNK, K, F>(
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
    counts
}

/// [`x_scan_cells`] with [`x_cell_avx2`], compiled for AVX2 as a whole so that the cell test is
/// inlined into the scan.
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn x_scan_cells_avx2<const BRANCH_FREE: bool, K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_scan_cells_with::<BRANCH_FREE, _, _>(lookup, subject, subject_length, |s, hits| {
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
unsafe fn x_scan_cells_avx2_defer<const BRANCH_FREE: bool, K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
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
    x_scan_cells_with::<BRANCH_FREE, _, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: the CPU has AVX2 (this function's contract).
        let candidates = unsafe { x_cell_avx2_defer(hits, s, table, test, candidate, &mut buffer) };
        counts.add_cell(hits.len(), candidates);
    });
    X_DEFER_BUFFER.set(buffer);
}

/// [`x_scan_cells_avx2`] with [`x_cell_avx2_defer_masked`] (LOSAT_X_BPMASKED in blastp's fused
/// loop, and blastx's fused loop with a value of LOSAT_X_BXDIAG or LOSAT_X_BXDIAGSHADOW). Never inlined, and its cell test is a function of its own,
/// so the loops of LOSAT_X_BPDEFER (`LOSAT_X_BPMASKED=0`) keep their code.
///
/// # Safety
/// The CPU must have AVX2.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline(never)]
unsafe fn x_scan_cells_avx2_defer_masked<const BRANCH_FREE: bool, K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
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
    x_scan_cells_with::<BRANCH_FREE, _, _>(lookup, subject, subject_length, |s, hits| {
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const X_CROSS_CHUNK: usize = 1024;

/// EXPERIMENT (LOSAT_X_BXCROSS, tests only): the chunk of the second instance the tests run,
/// which drains the stage every few cells (and inside long cells), so that the hits left by a
/// drain and the windows with a repeated entry occur across drains on small inputs.
#[cfg(all(
    test,
    any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )
))]
const X_CROSS_TEST_CHUNK: usize = 64;

/// EXPERIMENT (LOSAT_X_BXCROSS): the hits a piece of a cell's copy stages with fixed vector work
/// (three registers of eight lanes).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const X_CROSS_PIECE: usize = 3 * X_DIAG_LANES;

/// EXPERIMENT (LOSAT_X_BXCROSS): the room in each array past the chunk: a piece starts with fewer
/// than the chunk staged and stores [`X_CROSS_PIECE`] lanes (those past its hits are never read as hits),
/// so a drain starts with fewer than the chunk plus 24 staged, and its move of the fewer than eight
/// left reads eight lanes from below that count ([`x_cross_drain`]).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const X_CROSS_SLACK: usize = 32;

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const _: () = assert!(X_CROSS_SLACK >= X_CROSS_PIECE + X_DIAG_LANES - 1);
// The `u16` loads of `x_cross_cell` (three of `X_DIAG_LANES` hits, 16 bytes each, from the start of
// a piece of a list) end `X_CROSS_PIECE` hits from it: the padding of a table of 16-bit hits holds
// them, past a cell's own entries (2 bytes into the 8-byte cell, one piece) and past an overflow
// list (its last piece has at least one hit left).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const _: () = {
    use crate::algorithm::tblastx::lookup::BackboneCell16;
    assert!(
        2 + 2 * X_CROSS_PIECE <= (1 + <BackboneCell16 as XCell>::PAD_CELLS) * 8
            && X_CROSS_PIECE - 1 <= <BackboneCell16 as XCell>::PAD_HITS
    )
};

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[derive(Clone, Copy)]
struct XCross<'a, const CHUNK: usize> {
    base: *mut i32,
    _stage: std::marker::PhantomData<&'a mut Vec<i32>>,
}

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
fn x_cross_query(coord: i32, s: i32, mask: i32) -> i32 {
    coord.wrapping_add(s) & mask
}

// No NCBI counterpart: the body of LOSAT's scan with the cross-cell batches, written once.
/// EXPERIMENT (LOSAT_X_BXCROSS, LOSAT_X_BXCROSS128): the scan of a kernel `$kernel` of
/// LOSAT_X_BXCROSS's batches ([`XCrossKernel`]): `on_cell` stages the cell's hits ([`x_cross_cell`],
/// which drains the stage when `CHUNK` hits are staged), and the hits staged at the end of the
/// subject are tested by the kernel's drain. A macro, expanded in each kernel's scan function
/// ([`x_scan_cells_avx2_cross`], [`x_scan_cells_w128_cross`]), so that the closure of `on_cell`
/// is written in that function and has its target features (a closure takes those of the
/// function it is written in; written in a generic function without them, the AVX2 scan's
/// closure made LLVM inline the drain without `DEFER` into the scan, whose machine code then
/// differed from the build before LOSAT_X_BXCROSS128). `$branch_free`, `$defer`, `$chunk`, `$k`
/// and `$f` are the scan function's parameters `BRANCH_FREE`, `DEFER`, `CHUNK`, `K` and `F`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-509 (the scan and its hits in
// order; the rest of the loop at `x_scan_two_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
macro_rules! x_cross_scan {
    (
        $kernel:ty, $branch_free:ident, $defer:ident, $chunk:ident, $k:ident, $f:ident,
        ($lookup:ident, $subject:ident, $subject_length:ident, $table:ident, $test:ident,
         $candidate:ident, $counts:ident)
    ) => {{
        let mut stage = X_CROSS_STAGE.take();
        {
            let cross = XCross::<$chunk>::new::<$defer>(&mut stage);
            let mut staged = 0usize;
            x_scan_cells_with::<$branch_free, _, _>(
                $lookup,
                $subject,
                $subject_length,
                |s, hits| {
                    // SAFETY: as the scan function's contract; `staged` is 0 at the start and below
                    // the chunk after every call of `x_cross_cell`.
                    staged = unsafe {
                        x_cross_cell::<$kernel, $k::Hit, $defer, $chunk, $f>(
                            cross, staged, hits, s, $table, $test, $candidate, $counts,
                        )
                    };
                },
            );
            <$kernel as XCrossKernel>::drain::<$defer, $chunk, $f>(
                cross, staged, true, $table, $test, $candidate, $counts,
            );
        }
        X_CROSS_STAGE.set(stage);
    }};
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
    K: XCell,
    F,
>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XCrossCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    x_cross_scan!(
        XCrossAvx2,
        BRANCH_FREE,
        DEFER,
        CHUNK,
        K,
        F,
        (
            lookup,
            subject,
            subject_length,
            table,
            test,
            candidate,
            counts
        )
    );
}

/// Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): [`x_scan_cells`], or [`x_scan_cells_bf`] when `BRANCH_FREE` (the
/// same calls). One loop per instantiation, so `on_cell` has one call site and is inlined.
#[inline(always)]
fn x_scan_cells_with<const BRANCH_FREE: bool, K: XCell, C>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    on_cell: C,
) where
    C: FnMut(i32, &[K::Hit]),
{
    if BRANCH_FREE {
        x_scan_cells_bf(lookup, subject, subject_length, on_cell);
    } else {
        x_scan_cells(lookup, subject, subject_length, on_cell);
    }
}

/// The subject scan of [`x_scan_two_hit`]: `on_cell(subject_offset, hits)` for every scan position
/// whose lookup cell has hits, in scan order, with the cell's list in its stored order.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:60-113
#[inline(always)]
fn x_scan_cells<K: XCell, C>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    mut on_cell: C,
) where
    C: FnMut(i32, &[K::Hit]),
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
        let hits = backbone[index].hits(overflow);
        on_cell(s, hits);
    }
}

/// Branch-free scan (LOSAT_X_BXSCANBF, graduated 2026-10-26 for blastx): [`x_scan_cells`] with the cells that have hits found without a
/// branch on the presence bit: for each block of [`X_BLOCK`] positions [`x_block_candidates`] lists
/// `position << 32 | index`, the list's length advanced by the bit, and the listed cells then go to
/// `on_cell` in position order. The positions, indices and presence bits are those of
/// [`x_scan_cells`]'s loop (the same rolling index over the same bytes), so the calls are the same.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:60-99 (cited at `x_scan_cells`)
#[inline(always)]
fn x_scan_cells_bf<K: XCell, C>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    mut on_cell: C,
) where
    C: FnMut(i32, &[K::Hit]),
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
    // BLASTAA_SIZE, plus one; a table of 16-bit hits also has `XCell::PAD_CELLS` zero cells after
    // them) and the presence vector a bit for each, so the unchecked presence
    // reads of `x_block_candidates` stay inside it when every byte read is below BLASTAA_SIZE
    // (x_scan_cells would panic on an index past the vector instead).
    let largest = bytes.iter().copied().max().unwrap_or(0);
    assert!(
        usize::from(largest) < BLASTAA_SIZE
            && lookup.pv.len() << PV_ARRAY_BTS >= lookup.backbone.len() - K::PAD_CELLS,
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
            let hits = backbone[(c & 0xffff_ffff) as usize].hits(overflow);
            on_cell((c >> 32) as i32, hits);
        }
        block_first = block_last + 1;
    }
}

/// The diagonal test of a cell's hits one after the other (the fused loop without
/// LOSAT_X_BPDIAG).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[inline(always)]
fn x_cell_per_hit<H: XHit, F>(
    hits: &[H],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    for (i, &hit) in hits.iter().enumerate() {
        let query_offset = hit.q();
        // Read ahead the entry of a later hit of this cell: the hits of one subject offset
        // fall on distinct diagonals spread over the table, so each is a cache miss that the
        // hardware cannot predict. A prefetch reads nothing the test depends on.
        if let Some(&ahead) = hits.get(i + test.prefetch) {
            table.prefetch(ahead.q(), s);
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[inline(always)]
fn x_cell_per_hit_bf<H: XHit, F>(
    hits: &[H],
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
        table.prefetch_bf(later.q(), s);
        x_hit_bf(query_offset.q(), s, here, table, test, bf, candidate);
    }
    for &query_offset in rest {
        x_hit_bf(query_offset.q(), s, here, table, test, bf, candidate);
    }
}

/// EXPERIMENT (LOSAT_X_BPHITBF): one hit of [`x_cell_per_hit_bf`]: [`XHitBf::test`], then
/// [`x_test_hit`] and the candidate's call for a flag-1 entry or a candidate.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[inline(always)]
fn x_call_candidates<H: XHit, F>(
    mut lanes: u32,
    hits: &[H],
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
        let query_offset = hits[lane].q();
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
/// - write when `flag` is set and `sum >= last`, or `flag` is clear and `diff >= window`; the new entry is then `sum` truncated to 31
///   bits with the flag cleared (`set_last_hit` + `set_flag(0)`, or `set_last_hit` on an entry
///   without the flag), otherwise the old entry;
/// - candidate = !flag && word_size <= diff < window (its entry is not written).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[inline(always)]
fn x_test_batch_scalar<H: XHit, F>(
    hits: &[H],
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
        let entry = table.entry(query_offset.q(), s);
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
            .map(|&query_offset| table.coord(query_offset.q(), s))
            .collect::<Vec<_>>()
    ));
    x_call_candidates(lanes, hits, s, table, test, candidate);
    lanes.count_ones()
}

/// EXPERIMENT (LOSAT_X_BPDIAG): the cell's hits in batches of [`X_DIAG_LANES`], each by
/// [`x_test_batch_scalar`]; returns the number of candidates. The batches read no entry ahead
/// ([`X_PREFETCH_DISTANCE`] is the per-hit loop's): prefetching the batch 8 or 16 hits ahead was
/// not faster than none on 2026-10-14 (local, human block 1 x Species4 and the OrthoFinder
/// sample job, one thread).
#[inline(always)]
fn x_cell_scalar<H: XHit, F>(
    hits: &[H],
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
unsafe fn x_test_batch_avx2<H: XHit, F>(
    hits: &[H],
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
unsafe fn x_batch_avx2<H: XHit>(
    hits: &[H],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> (std::arch::x86_64::__m256i, u32) {
    use std::arch::x86_64::*;
    assert!(hits.len() == X_DIAG_LANES);
    let subject_plus_offset = s.wrapping_add(test.diag_offset);
    let written = (subject_plus_offset as u32) & X_LAST_HIT_MASK;
    // SAFETY: (caller) the CPU has AVX2; `hits` has eight hits (asserted above).
    let query = unsafe { H::x_load8(hits.as_ptr()) };
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
    // keep (no write): `sum < last` when `flag` is set, else `diff < window`
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
#[cfg(any(target_arch = "x86_64", test))]
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2_defer<H: XHit, F>(
    hits: &[H],
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

/// EXPERIMENT (LOSAT_X_BPMASKED): the test of [`x_batch_avx2`] for eight lanes
/// whose entry coordinates the caller computed: reads the eight entries, stores their new values
/// and returns the candidate lanes (bit `i` = lane `i`). Every load comes before every store, so
/// lanes with the same coordinate (the masked batch of [`x_cell_avx2_defer_masked`]) read the same
/// value and store the same new one.
///
/// # Safety
/// The CPU must have AVX2, and every coordinate must be at most `table.mask`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
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

/// EXPERIMENT (LOSAT_X_BPMASKED): [`x_cell_avx2_defer`] with two changes of the
/// batches (module header), for blastp's and blastx's fused loops: a batch of eight takes its
/// lanes' entry coordinates from
/// scalar loads of its hits, and the last batch of fewer than eight hits is one masked AVX2 batch
/// (the missing lanes repeat the first hit, their candidate bits are dropped) whose candidates go
/// into the buffer after the others'. The closure is then called for the buffer in order, so the
/// calls are [`x_cell_avx2_defer`]'s, in the same order. Returns the number of candidates.
///
/// # Safety
/// The CPU must have AVX2.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn x_cell_avx2_defer_masked<H: XHit, F>(
    hits: &[H],
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
        // SAFETY: the batch has eight hits.
        let query = H::x_load8(p);
        // The coordinates from scalar loads of the hits, not taken out of `query`: the volatile
        // reads keep the compiler from extracting them from the register, which would start the
        // chain of the entry loads 3-6 cycles later (line R's probe of blastp's cell test,
        // 2026-10-18).
        // i < X_DIAG_LANES = batch.len()
        let lane_coords: [u32; X_DIAG_LANES] =
            std::array::from_fn(|i| (H::x_read_volatile(p.add(i)).wrapping_sub(s) & mask) as u32);
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
        // SAFETY (`x_load_tail`): the lanes below `n` are set in `in_rest`; for `u16` hits the 16
        // bytes at `rest` lie inside the table's allocation (the cells and the overflow array
        // are padded, `XCell::PAD_CELLS` and `XCell::PAD_HITS`).
        let query = _mm256_blendv_epi8(
            _mm256_set1_epi32(rest[0].q()),
            H::x_load_tail(rest.as_ptr(), in_rest),
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
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
    x_cross_drain_with::<XCrossAvx2, DEFER, CHUNK, F>(
        cross, staged, end, table, test, candidate, counts,
    )
}

/// EXPERIMENT (LOSAT_X_BXCROSS): [`x_test_hit`] for staged hit `i`, its candidate called at once
/// (the reference's order) with the query offset of [`x_cross_query`]; returns 1 for a candidate,
/// else 0.
///
/// # Safety
/// `i` must be below the staged count, and the hit's query offset at most `table.mask`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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

/// EXPERIMENT (LOSAT_X_BXCROSS): [`x_batch_lanes`]'s test, which returns with the candidate
/// lanes the eight entries it read. EXPERIMENT (LOSAT_X_BXCROSSDEFER, with `MARK`): it returns
/// `None` before any store where an entry holds [`X_CROSS_MARKER`], and its candidate lanes store
/// the marker instead of the entry they read (which they store without `MARK`).
///
/// # Safety
/// As [`x_batch_lanes`].
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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

#[cfg(all(
    test,
    any(
        target_arch = "x86_64",
        all(target_arch = "wasm32", target_feature = "simd128")
    )
))]
thread_local! {
    /// EXPERIMENT (LOSAT_X_BXCROSSDEFER, tests only): the units of this thread that met the marker
    /// with calls deferred (the calls made, the unit run again) and with none (a stored value, the
    /// unit run with its calls at once).
    static X_CROSS_MARKS_MET: std::cell::Cell<[u64; 2]> = const { std::cell::Cell::new([0; 2]) };
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER, tests only): counts a unit that met the marker with calls
/// deferred (`stored` false) or none (true).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const X_CROSS_MARKER: i32 = -1;

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
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`): a hit
// that passes the test reaches the context check with `last_hit` and `diff`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
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

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_batch`] (no marker, its candidates called at
/// once) out of the drain's loop, for a batch that met a stored value equal to the marker while
/// no call was deferred. Returns the number of candidates.
///
/// # Safety
/// As [`x_cross_batch`].
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
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
    <XCrossAvx2 as XCrossKernel>::batch::<CHUNK, F>(cross, b, table, test, candidate)
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the staged hits `b..b + 8`, two of which name one entry (hits of
/// different cells on one diagonal), in parts (module header). From lane `start` (0 first), the
/// part is the lanes before the first one that repeats an entry of the part
/// ([`x_cross_part_end`]), so its entries are distinct, and the part as one batch gives what its
/// hits give one after the other. The batch's lanes past the part repeat the part's first hit (its
/// coordinate and subject offset): they read its entry and store, after it, the value it stores
/// (every load of [`x_batch_lanes`] comes before its stores), as the masked batch of
/// [`x_cell_avx2_defer_masked`] does; only the part's candidate bits are kept, and their
/// candidates are called in lane order ([`x_cross_calls`]) before the next part, which starts at
/// the repeating lane, loads. Nothing past `b + 8` is read. Returns the number of candidates. Cold
/// and never inlined, so that the drain's loop keeps none of its values.
///
/// # Safety
/// The CPU must have AVX2 and POPCNT, `b + 8` must be at most the staged count, and the hits'
/// query offsets at most `table.mask`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
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
        let lanes = <XCrossAvx2 as XCrossKernel>::lanes(lane_coords, s, table, test)
            & ((1u32 << (end - start)) - 1);
        x_cross_calls(cross, b + start, lanes, table, test, candidate);
        candidates += lanes.count_ones();
        start = end;
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_batch_avx2_lanes_mark`] until it meets no marker (the
/// deferred calls made in between), then its candidates deferred ([`x_cross_defer`]); a marker met
/// with no call deferred is a stored value, and the unit then runs with its calls at once. Returns
/// the candidate lanes (within `keep`). For the parts of [`x_cross_window_defer`].
///
/// # Safety
/// As [`x_batch_lanes`], and `pending` below [`XCross::PENDING`] with the records of this
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
            <XCrossAvx2 as XCrossKernel>::defer::<CHUNK>(cross, pending, b, lanes, old);
            return lanes;
        }
        if *pending == 0 {
            // A stored value equal to the marker, no call deferred: the unit with its calls at
            // once (as without the deferral).
            x_cross_mark_met(true);
            let lanes = <XCrossAvx2 as XCrossKernel>::lanes(lane_coords, s, table, test) & keep;
            x_cross_calls(cross, b, lanes, table, test, candidate);
            return lanes;
        }
        x_cross_flush_cold(cross, *pending, table, test, candidate);
        *pending = 0;
    }
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_window`] with the parts' candidates deferred
/// ([`x_cross_unit`]): a part that repeats a deferred candidate's entry (the next part starts at
/// the lane that repeats an entry of the one before) meets the marker, and the deferred calls are
/// made before it is tested. Takes and returns the drain's count of records; returns first the
/// number of candidates. Cold and never inlined, as [`x_cross_window`].
///
/// # Safety
/// As [`x_cross_window`], and `pending` below [`XCross::PENDING`] with the records of this drain.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`)
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

// ---- EXPERIMENT (LOSAT_X_BXCROSS128): LOSAT_X_BXCROSS's batches in 128-bit vectors ----

// No NCBI counterpart: the vector kernel of LOSAT's cell test; it holds no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSS, LOSAT_X_BXCROSS128): a kernel of LOSAT_X_BXCROSS's batches of
/// eight staged hits: what the scan's stage (`x_cross_scan!`, [`x_cross_cell`]), the drain
/// ([`x_cross_drain_with`]) and its batches run in the kernel's vectors, and the kernel's window
/// with a repeated entry. [`XCrossAvx2`]: one AVX2 register of eight lanes; [`XCrossW128`]: two
/// vectors of four lanes of the 128-bit layer (`utils/xdrop_simd.rs` `arch`, x86-64 SSE and wasm32
/// simd128). The kernel's functions are compiled for its instructions and inline the generic
/// bodies with their kernel.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
trait XCrossKernel: Sized {
    /// Eight `i32` lanes: eight staged values, the subject offsets of a batch, the entries a unit
    /// read (LOSAT_X_BXCROSSDEFER's records).
    type Lanes: Copy;
    /// The vectors a cell's stage uses (its subject offset and the mask, made once per cell).
    type Cell: Copy;

    /// The vectors of a cell at subject offset `s` with the table's `mask`.
    ///
    /// # Safety
    /// Inlined into code compiled for the kernel's instructions.
    unsafe fn cell(s: i32, mask: i32) -> Self::Cell;

    /// Stages the first `piece` (at most [`X_CROSS_PIECE`]) hits from `p`: their coordinates
    /// `(q - s) & mask` at `at` and the subject offset at `at + stride`, [`X_CROSS_PIECE`] lanes
    /// each (those past `piece` are never read as hits).
    ///
    /// # Safety
    /// As [`x_cross_cell`]; the [`X_CROSS_PIECE`] lanes from `at` and from `at + stride`
    /// writable.
    unsafe fn stage<H: XHit>(
        cell: Self::Cell,
        p: *const H,
        piece: usize,
        at: *mut i32,
        stride: usize,
    );

    /// The eight values at `p`.
    ///
    /// # Safety
    /// Inlined into code compiled for the kernel's instructions; eight values readable at `p`.
    unsafe fn load8(p: *const i32) -> Self::Lanes;

    /// `v` stored at `p`.
    ///
    /// # Safety
    /// Inlined into code compiled for the kernel's instructions; eight values writable at `p`.
    unsafe fn store8(p: *mut i32, v: Self::Lanes);

    /// Whether two of the eight coordinates `coords` are equal.
    ///
    /// # Safety
    /// Inlined into code compiled for the kernel's instructions.
    unsafe fn conflict(coords: Self::Lanes) -> bool;

    /// The test of eight lanes (lane `i` the hit at `lane_coords[i]` with the subject offset of
    /// lane `i` of `s`) of [`x_batch_avx2_lanes_mark`]: reads the eight entries, stores their new
    /// values and returns the candidate lanes and the entries read; with `MARK` it returns `None`
    /// before any store where an entry holds [`X_CROSS_MARKER`], and its candidate lanes store the
    /// marker.
    ///
    /// # Safety
    /// Inlined into code compiled for the kernel's instructions; every coordinate at most
    /// `table.mask`.
    unsafe fn lanes_mark<const MARK: bool>(
        lane_coords: [u32; X_DIAG_LANES],
        s: Self::Lanes,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
    ) -> Option<(u32, Self::Lanes)>;

    /// The kernel's [`x_batch_lanes`]. This method and [`XCrossKernel::batch`] (and the AVX2
    /// kernel's `batch_marked` and `defer`) are the kernel's functions of the batches (compiled for
    /// its instructions, `#[target_feature]` with `#[inline]`, each the generic body), so that the
    /// drain and the windows call them as LOSAT_X_BXCROSS's AVX2 functions called theirs before
    /// LOSAT_X_BXCROSS128 (with the bodies inlined into the drain as generic functions, rustc
    /// 1.92.0 and its LLVM laid out the drain and the deferred window otherwise).
    ///
    /// # Safety
    /// As [`x_batch_lanes`].
    unsafe fn lanes(
        lane_coords: [u32; X_DIAG_LANES],
        s: Self::Lanes,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
    ) -> u32;

    /// The kernel's [`x_cross_batch`].
    ///
    /// # Safety
    /// As [`x_cross_batch`].
    unsafe fn batch<const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> u32
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct);

    /// The kernel's [`x_cross_batch_marked`]. This method, [`XCrossKernel::defer`],
    /// [`XCrossKernel::window_defer`] and [`XCrossKernel::batch_cold`] are the deferral's
    /// (LOSAT_X_BXCROSSDEFER), which only a drain with `DEFER` runs: the AVX2 kernel gives them;
    /// their defaults are for a kernel that never defers (LOSAT_X_BXCROSS128's, whose `drain`
    /// refuses `DEFER` when it is compiled), and no program path reaches them.
    ///
    /// # Safety
    /// As [`x_cross_batch_marked`].
    unsafe fn batch_marked<const CHUNK: usize>(
        _: XCross<'_, CHUNK>,
        _: usize,
        _: &mut XDiagTable<'_>,
        _: XDiagTest,
    ) -> Option<(u32, Self::Lanes)> {
        unreachable!("{X_CROSS_NO_DEFER}")
    }

    /// The kernel's [`x_cross_defer`] (the deferral's, as [`XCrossKernel::batch_marked`]).
    ///
    /// # Safety
    /// As [`x_cross_defer`].
    unsafe fn defer<const CHUNK: usize>(
        _: XCross<'_, CHUNK>,
        _: &mut usize,
        _: usize,
        _: u32,
        _: Self::Lanes,
    ) {
        unreachable!("{X_CROSS_NO_DEFER}")
    }

    /// The kernel's [`x_cross_drain_with`], out of line. A kernel's method that calls its drain is
    /// compiled for its instructions (`#[target_feature]` with `#[inline]`, not
    /// `#[inline(always)]`): rustc 1.92.0 (the toolchain of `rust-toolchain.toml`, which this
    /// layout was checked with; nothing in CI guards it) marks a call of a `#[target_feature]`
    /// function's `#[inline(never)]` only at call sites in functions with those features, so a call
    /// written in a function without them (the generic driver) loses it once inlined into the
    /// kernel's function, and LLVM then inlined the drain without `DEFER` into the AVX2 scan. The
    /// methods that call the kernel's cold functions (`window`, and the AVX2 kernel's
    /// `window_defer` and `batch_cold`) are `#[inline(always)]`: LLVM does not inline a cold
    /// function, and with them inlined early the AVX2 drain without `DEFER` keeps the machine code
    /// of the build before LOSAT_X_BXCROSS128 (with `#[inline]` it did not); the drain with `DEFER`
    /// keeps its instructions, with its blocks laid out otherwise and a frame 16 bytes smaller.
    ///
    /// # Safety
    /// As [`x_cross_drain_with`].
    #[allow(clippy::too_many_arguments)]
    unsafe fn drain<const DEFER: bool, const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        staged: usize,
        end: bool,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
        counts: &mut XCrossCounts,
    ) -> usize
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct);

    /// The kernel's window of eight staged hits with a repeated entry, cold and out of line (the
    /// AVX2 kernel's parts, [`x_cross_window`]; the 128-bit kernel's hits one after the other,
    /// [`x_cross_window_w128`]).
    ///
    /// # Safety
    /// As [`x_cross_window`].
    unsafe fn window<const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> u32
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct);

    /// The kernel's window with a repeated entry in a drain that defers its calls (the deferral's,
    /// as [`XCrossKernel::batch_marked`]: the AVX2 kernel's [`x_cross_window_defer`], cold and out
    /// of line).
    ///
    /// # Safety
    /// As [`x_cross_window_defer`].
    unsafe fn window_defer<const CHUNK: usize, F>(
        _: XCross<'_, CHUNK>,
        _: usize,
        _: usize,
        _: &mut XDiagTable<'_>,
        _: XDiagTest,
        _: &mut F,
    ) -> (u32, usize)
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        unreachable!("{X_CROSS_NO_DEFER}")
    }

    /// The kernel's [`x_cross_batch`] for a batch that met a stored value equal to the marker
    /// while no call was deferred (the deferral's, as [`XCrossKernel::batch_marked`]: the AVX2
    /// kernel's [`x_cross_batch_cold`], cold and out of line).
    ///
    /// # Safety
    /// As [`x_cross_batch`].
    unsafe fn batch_cold<const CHUNK: usize, F>(
        _: XCross<'_, CHUNK>,
        _: usize,
        _: &mut XDiagTable<'_>,
        _: XDiagTest,
        _: &mut F,
    ) -> u32
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        unreachable!("{X_CROSS_NO_DEFER}")
    }
}

// No NCBI counterpart: a message of LOSAT's kernel trait.
/// EXPERIMENT (LOSAT_X_BXCROSS128): the message of the deferral's defaults of [`XCrossKernel`],
/// which a kernel that never defers its calls leaves in place and no drain reaches.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
const X_CROSS_NO_DEFER: &str = "the kernel does not defer its calls (LOSAT_X_BXCROSSDEFER)";

/// EXPERIMENT (LOSAT_X_BXCROSS): stages one cell's hits (module header) in pieces of
/// [`X_CROSS_PIECE`] with the kernel `Kn` ([`XCrossKernel::stage`]): per piece the coordinates
/// `(q - s) & mask` by lane and the subject offset, stored into each of the two arrays at the
/// staged count, all from one base register at constant displacements (AVX2: three loads of eight
/// lanes of the list, `XHit::x_load_tail`: for `i32` hits the lanes past the cell's hits read
/// nothing, for `u16` hits each is a 16-byte widening load whose lanes past the hits are not used;
/// LOSAT_X_BXCROSS128: the six vectors of four lanes of `XHit::x_load24_w128`); then the count
/// advances by the piece's hits, and the batches run (`Kn`'s drain, [`x_cross_drain`] or
/// [`x_cross_drain_w128`]) when `CHUNK` or more are staged. A cell of up to 24 hits is one piece.
/// The cell and its hits are counted here, its candidates by the drain. Returns the number of hits
/// staged.
///
/// # Safety
/// Inlined into code compiled for `Kn`'s instructions ([`x_scan_cells_avx2_cross`]: AVX2 and
/// POPCNT; [`x_scan_cells_w128_cross`]); `staged` must be below `CHUNK` (as it is at the start of a
/// subject and after every call of this function), and every query offset of `hits` must be at
/// most `table.mask` ([`x_cross_query`]). For `u16` hits the 48 bytes from the list's start must
/// lie inside the table's allocations (`XCell::PAD_CELLS`, `XCell::PAD_HITS`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:107-114 (the copy of a cell's
// hits) and aa_ungapped.c:516 (the coordinate)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn x_cross_cell<Kn: XCrossKernel, H: XHit, const DEFER: bool, const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    mut staged: usize,
    hits: &[H],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XCrossCounts,
) -> usize
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    debug_assert!(
        hits.iter().all(|q| (0..=table.mask).contains(&q.q())),
        "a query offset of the cell at subject offset {s} is above the mask {}",
        table.mask
    );
    counts.diag.add_staged(hits.len());
    let cell = Kn::cell(s, table.mask);
    let mut p = hits.as_ptr();
    let mut left = hits.len();
    loop {
        let piece = left.min(X_CROSS_PIECE);
        // staged < CHUNK, so the piece's lanes end at staged + 24 < STRIDE in each array
        // (`XCross::new`).
        Kn::stage(
            cell,
            p,
            piece,
            cross.coord().add(staged),
            XCross::<CHUNK>::STRIDE,
        );
        staged += piece;
        if staged >= CHUNK {
            staged =
                Kn::drain::<DEFER, CHUNK, F>(cross, staged, false, table, test, candidate, counts);
        }
        if left <= X_CROSS_PIECE {
            break;
        }
        left -= X_CROSS_PIECE;
        p = p.add(X_CROSS_PIECE);
    }
    staged
}

/// EXPERIMENT (LOSAT_X_BXCROSS128): [`x_cross_drain`] with the kernel `Kn`: the staged hits in
/// windows of eight from the first one, a window with distinct entries as one batch, otherwise in
/// parts (`Kn`'s cold window functions); with `DEFER` (LOSAT_X_BXCROSSDEFER) each unit's calls
/// deferred to a later unit that meets their entries, or the end of the windows. Then with `end`
/// the fewer than eight left go one after the other, otherwise they move to the front of the
/// stage. Returns the hits left staged.
///
/// # Safety
/// Inlined into code compiled for `Kn`'s instructions; as [`x_cross_drain`].
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
unsafe fn x_cross_drain_with<Kn: XCrossKernel, const DEFER: bool, const CHUNK: usize, F>(
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
    let full = staged & !(X_DIAG_LANES - 1);
    let mut candidates = 0u32;
    let mut windows = 0u32;
    // The records of deferred calls (LOSAT_X_BXCROSSDEFER; always 0 without `DEFER`).
    let mut pending = 0usize;
    let mut b = 0;
    while b < full {
        if Kn::conflict(Kn::load8(cross.coord().add(b))) {
            if DEFER {
                let (window, left) = Kn::window_defer(cross, b, pending, table, test, candidate);
                candidates += window;
                pending = left;
            } else {
                candidates += Kn::window(cross, b, table, test, candidate);
            }
            windows += 1;
        } else if !DEFER {
            candidates += Kn::batch::<CHUNK, F>(cross, b, table, test, candidate);
        } else {
            match Kn::batch_marked::<CHUNK>(cross, b, table, test) {
                Some((lanes, old)) => {
                    Kn::defer::<CHUNK>(cross, &mut pending, b, lanes, old);
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
                None => candidates += Kn::batch_cold(cross, b, table, test, candidate),
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
        // = STRIDE, inside each array. Each row's lanes are loaded before they are stored.
        for row in [cross.coord(), cross.subject()] {
            let lanes = Kn::load8(row.add(full));
            Kn::store8(row, lanes);
        }
        staged - full
    };
    counts.diag.add_candidates(candidates);
    counts.splits += u64::from(windows);
    left
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the test of [`x_batch_avx2_at`] for eight lanes with a subject
/// offset each (`s`, lane `i` the subject offset of the hit at `lane_coords[i]`), with the kernel
/// `Kn` ([`XCrossKernel::lanes_mark`] without the marker): reads the eight entries, stores their
/// new values and returns the candidate lanes (bit `i` = lane `i`). The lane formulas are
/// `x_batch_avx2_at`'s with the broadcast subject offset replaced by `s`:
/// `subject_plus_offset = s + diag_offset`, `written = (s + diag_offset) & X_LAST_HIT_MASK`,
/// `diff = s - (last - diag_offset)`. Every load comes before every store, so lanes with the same
/// coordinate and subject offset (the repeated lanes of [`x_cross_window`]) read the same value
/// and store the same new one.
///
/// # Safety
/// Inlined into code compiled for `Kn`'s instructions; every coordinate must be at most
/// `table.mask`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
unsafe fn x_batch_lanes<Kn: XCrossKernel>(
    lane_coords: [u32; X_DIAG_LANES],
    s: Kn::Lanes,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> u32 {
    let Some((lanes, _)) = Kn::lanes_mark::<false>(lane_coords, s, table, test) else {
        unreachable!("without the marker the batch always stores")
    };
    lanes
}

/// EXPERIMENT (LOSAT_X_BXCROSS): the staged hits `b..b + 8`, whose entries are distinct
/// (`Kn`'s conflict test, [`x_cross_conflict`]), as one batch (module header): the lanes'
/// coordinates from scalar loads of the stage, [`x_batch_lanes`] with the staged subject offsets,
/// and the candidates called in lane order ([`x_cross_calls`]). Returns the number of candidates.
///
/// # Safety
/// Inlined into code compiled for `Kn`'s instructions; `b + 8` must be at most the staged count,
/// and the hits' query offsets at most `table.mask`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
unsafe fn x_cross_batch<Kn: XCrossKernel, const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let p = cross.coord().add(b);
    // The coordinates from scalar loads, not taken out of the conflict test's vectors (the
    // volatile reads keep the compiler from extracting them; x_cell_avx2_defer_masked).
    let lane_coords: [u32; X_DIAG_LANES] =
        std::array::from_fn(|i| std::ptr::read_volatile(p.add(i)) as u32);
    debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
    let s = Kn::load8(cross.subject().add(b));
    let lanes = Kn::lanes(lane_coords, s, table, test);
    x_cross_calls(cross, b, lanes, table, test, candidate);
    lanes.count_ones()
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): [`x_cross_batch`] with the marker
/// ([`XCrossKernel::lanes_mark`], [`x_batch_avx2_lanes_mark`]) and without the calls: returns its
/// candidate lanes and the entries it read, or `None` (nothing stored) where an entry holds the
/// marker ([`x_cross_drain`] then makes the deferred calls and runs the batch again).
///
/// # Safety
/// As [`x_cross_batch`].
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn x_cross_batch_marked<Kn: XCrossKernel, const CHUNK: usize>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
) -> Option<(u32, Kn::Lanes)> {
    let p = cross.coord().add(b);
    // The coordinates from scalar loads, as in x_cross_batch.
    let lane_coords: [u32; X_DIAG_LANES] =
        std::array::from_fn(|i| std::ptr::read_volatile(p.add(i)) as u32);
    debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
    let s = Kn::load8(cross.subject().add(b));
    Kn::lanes_mark::<true>(lane_coords, s, table, test)
}

/// EXPERIMENT (LOSAT_X_BXCROSSDEFER): records a unit's candidates for deferred calls: its first
/// staged hit `b`, its candidate lanes and the eight entries it read go to record `pending`
/// (written for every unit, without a branch), and the count advances when the unit has
/// candidates.
///
/// # Safety
/// Inlined into code compiled for `Kn`'s instructions; `pending` must be below
/// [`XCross::PENDING`].
// No NCBI counterpart: NCBI calls each candidate at once (aa_ungapped.c:513-552, cited at
// `x_test_hit`); the record of a deferred call is LOSAT's and holds no value NCBI computes.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
unsafe fn x_cross_defer<Kn: XCrossKernel, const CHUNK: usize>(
    cross: XCross<'_, CHUNK>,
    pending: &mut usize,
    b: usize,
    lanes: u32,
    old: Kn::Lanes,
) {
    let r = *pending;
    debug_assert!(r < XCross::<CHUNK>::PENDING);
    Kn::store8(cross.pending_old().add(X_DIAG_LANES * r), old);
    *cross.pending_at().add(2 * r) = b as i32;
    *cross.pending_at().add(2 * r + 1) = lanes as i32;
    *pending = r + usize::from(lanes != 0);
}

// No NCBI counterpart: writes two methods of LOSAT's kernel trait; they hold no value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXCROSS128): a kernel's [`XCrossKernel::lanes`] and [`XCrossKernel::batch`],
/// written once for both kernels: each the generic body ([`x_batch_lanes`], [`x_cross_batch`])
/// with the kernel's attributes (its `target_feature`) and `#[inline]` (the trait's notes say why).
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
macro_rules! x_cross_batch_methods {
    ($(#[$feature:meta])*) => {
        $(#[$feature])*
        #[inline]
        unsafe fn lanes(
            lane_coords: [u32; X_DIAG_LANES],
            s: Self::Lanes,
            table: &mut XDiagTable<'_>,
            test: XDiagTest,
        ) -> u32 {
            x_batch_lanes::<Self>(lane_coords, s, table, test)
        }

        $(#[$feature])*
        #[inline]
        unsafe fn batch<const CHUNK: usize, F>(
            cross: XCross<'_, CHUNK>,
            b: usize,
            table: &mut XDiagTable<'_>,
            test: XDiagTest,
            candidate: &mut F,
        ) -> u32
        where
            F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
        {
            x_cross_batch::<Self, CHUNK, F>(cross, b, table, test, candidate)
        }
    };
}

// No NCBI counterpart: the vector kernel of LOSAT's cell test; it holds no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSS): the kernel of one AVX2 register of eight lanes per batch: the
/// conflict test [`x_cross_conflict`], the batch [`x_batch_avx2_lanes_mark`], and three registers
/// of eight lanes per piece of a cell's stage. Its out-of-line functions ([`x_scan_cells_avx2_cross`],
/// [`x_cross_drain`], [`x_cross_window`], [`x_cross_window_defer`], [`x_cross_batch_cold`]) are
/// compiled for AVX2 and POPCNT.
#[cfg(target_arch = "x86_64")]
#[derive(Clone, Copy)]
struct XCrossAvx2;

#[cfg(target_arch = "x86_64")]
impl XCrossKernel for XCrossAvx2 {
    type Lanes = std::arch::x86_64::__m256i;
    type Cell = [std::arch::x86_64::__m256i; 2];

    // No NCBI counterpart: the vectors of LOSAT's stage of a cell (its subject offset, the mask).
    #[inline(always)]
    unsafe fn cell(s: i32, mask: i32) -> Self::Cell {
        use std::arch::x86_64::*;
        [_mm256_set1_epi32(s), _mm256_set1_epi32(mask)]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:107-114 (the copy of a cell's
    // hits) and aa_ungapped.c:516 (the coordinate)
    #[inline(always)]
    unsafe fn stage<H: XHit>(
        cell: Self::Cell,
        p: *const H,
        piece: usize,
        at: *mut i32,
        stride: usize,
    ) {
        use std::arch::x86_64::*;
        let [subject, mask] = cell;
        // The positions in a piece of the lanes of its three registers.
        let positions = [
            _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7),
            _mm256_setr_epi32(8, 9, 10, 11, 12, 13, 14, 15),
            _mm256_setr_epi32(16, 17, 18, 19, 20, 21, 22, 23),
        ];
        let n = _mm256_set1_epi32(piece as i32);
        for (k, &position) in positions.iter().enumerate() {
            // The lanes at positions below `piece`; the pointer may lie past the list
            // (`wrapping_add`), where the masked load reads nothing.
            let in_piece = _mm256_cmpgt_epi32(n, position);
            // `i32` hits: `vpmaskmovd`, which reads nothing past the list. `u16` hits: one
            // `vpmovzxwd ymm, m128` of the 16 bytes at the lane, whose lanes past the list are
            // never read (the stage's count advances by the piece only) and which lie inside the
            // table (`BackboneCell16::PAD_CELLS`, `PAD_HITS`: the three loads end 48 bytes from
            // the list's start).
            let query = H::x_load_tail(p.wrapping_add(X_DIAG_LANES * k), in_piece);
            let coords = _mm256_and_si256(_mm256_sub_epi32(query, subject), mask);
            let lane = X_DIAG_LANES * k;
            // The two stores of a register are `stride` values apart.
            _mm256_storeu_si256(at.add(lane).cast(), coords);
            _mm256_storeu_si256(at.add(stride + lane).cast(), subject);
        }
    }

    // No NCBI counterpart: a vector load of LOSAT's stage.
    #[inline(always)]
    unsafe fn load8(p: *const i32) -> Self::Lanes {
        std::arch::x86_64::_mm256_loadu_si256(p.cast())
    }

    // No NCBI counterpart: a vector store of LOSAT's stage.
    #[inline(always)]
    unsafe fn store8(p: *mut i32, v: Self::Lanes) {
        std::arch::x86_64::_mm256_storeu_si256(p.cast(), v)
    }

    // No NCBI counterpart: a check of the property the batches rest on (module header).
    #[inline(always)]
    unsafe fn conflict(coords: Self::Lanes) -> bool {
        x_cross_conflict(coords)
    }

    // No NCBI counterpart: the kernel's `x_batch_avx2_lanes_mark`, which holds the NCBI reference.
    #[inline(always)]
    unsafe fn lanes_mark<const MARK: bool>(
        lane_coords: [u32; X_DIAG_LANES],
        s: Self::Lanes,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
    ) -> Option<(u32, Self::Lanes)> {
        x_batch_avx2_lanes_mark::<MARK>(lane_coords, s, table, test)
    }

    x_cross_batch_methods!(#[target_feature(enable = "avx2,popcnt")]);

    // No NCBI counterpart: the kernel's `x_cross_batch_marked`, which holds the NCBI reference.
    #[target_feature(enable = "avx2,popcnt")]
    #[inline]
    unsafe fn batch_marked<const CHUNK: usize>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
    ) -> Option<(u32, Self::Lanes)> {
        x_cross_batch_marked::<Self, CHUNK>(cross, b, table, test)
    }

    // No NCBI counterpart: the kernel's `x_cross_defer` (LOSAT's record of a deferred call).
    #[target_feature(enable = "avx2,popcnt")]
    #[inline]
    unsafe fn defer<const CHUNK: usize>(
        cross: XCross<'_, CHUNK>,
        pending: &mut usize,
        b: usize,
        lanes: u32,
        old: Self::Lanes,
    ) {
        x_cross_defer::<Self, CHUNK>(cross, pending, b, lanes, old)
    }

    // No NCBI counterpart: the kernel's `x_cross_drain`, which holds the NCBI reference.
    #[target_feature(enable = "avx2,popcnt")]
    #[inline]
    unsafe fn drain<const DEFER: bool, const CHUNK: usize, F>(
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
        x_cross_drain::<DEFER, CHUNK, F>(cross, staged, end, table, test, candidate, counts)
    }

    // No NCBI counterpart: the kernel's `x_cross_window`, which holds the NCBI reference.
    #[inline(always)]
    unsafe fn window<const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> u32
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        x_cross_window(cross, b, table, test, candidate)
    }

    // No NCBI counterpart: the kernel's `x_cross_window_defer`, which holds the NCBI reference.
    #[inline(always)]
    unsafe fn window_defer<const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        pending: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> (u32, usize)
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        x_cross_window_defer(cross, b, pending, table, test, candidate)
    }

    // No NCBI counterpart: the kernel's `x_cross_batch_cold`, which holds the NCBI reference.
    #[inline(always)]
    unsafe fn batch_cold<const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> u32
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        x_cross_batch_cold(cross, b, table, test, candidate)
    }
}

// No NCBI counterpart: the vector kernel of LOSAT's cell test; it holds no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BXCROSS128): the kernel of two 128-bit vectors of four lanes per batch of
/// eight (`utils/xdrop_simd.rs` `arch`): the conflict test compares each half with itself one and
/// two places round and the low half with the high one zero to three places round; the batch's
/// entries are read and stored by lane, as the AVX2 kernel's; the candidate lanes are the two
/// halves' masks. It calls each batch's candidates at once: it keeps the trait's defaults of the
/// deferral, which no drain reaches (its `drain` and `lanes_mark` refuse `DEFER` and the marker
/// when they are compiled). On x86-64 its functions are compiled for SSSE3, SSE4.1 and POPCNT:
/// SSE4.1 for `ptest`, `pblendvb` and `pmovzxwd`, POPCNT for the count of a batch's candidates, and
/// SSSE3 because it is part of the 128-bit layer's contract (the layer's other users, the X-drop
/// kernel and LOSAT_X_BPGAP128, ask for SSSE3 and SSE4.1), though no operation of this kernel needs
/// it.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[derive(Clone, Copy)]
struct XCrossW128;

#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
impl XCrossKernel for XCrossW128 {
    type Lanes = [crate::utils::xdrop_simd::arch::V; 2];
    type Cell = [crate::utils::xdrop_simd::arch::V; 2];

    // No NCBI counterpart: a vector load of LOSAT's stage.
    #[inline(always)]
    unsafe fn load8(p: *const i32) -> Self::Lanes {
        use crate::utils::xdrop_simd::arch::*;
        [load(p), load(p.add(4))]
    }

    // No NCBI counterpart: a vector store of LOSAT's stage.
    #[inline(always)]
    unsafe fn store8(p: *mut i32, v: Self::Lanes) {
        use crate::utils::xdrop_simd::arch::*;
        store(p, v[0]);
        store(p.add(4), v[1]);
    }

    // No NCBI counterpart: the vectors of LOSAT's stage of a cell (its subject offset, the mask).
    #[inline(always)]
    unsafe fn cell(s: i32, mask: i32) -> Self::Cell {
        use crate::utils::xdrop_simd::arch::*;
        [splat_i32(s), splat_i32(mask)]
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:107-114 (the copy of a cell's
    // hits) and aa_ungapped.c:516 (the coordinate)
    #[inline(always)]
    unsafe fn stage<H: XHit>(
        cell: Self::Cell,
        p: *const H,
        piece: usize,
        at: *mut i32,
        stride: usize,
    ) {
        use crate::utils::xdrop_simd::arch::*;
        let [subject, mask] = cell;
        let queries = H::x_load24_w128(p, piece);
        for (k, &query) in queries.iter().enumerate() {
            store(at.add(4 * k), and(sub32(query, subject), mask));
            store(at.add(stride + 4 * k), subject);
        }
    }

    // No NCBI counterpart: a check of the property the batches rest on (module header).
    #[inline(always)]
    unsafe fn conflict(coords: Self::Lanes) -> bool {
        use crate::utils::xdrop_simd::arch::*;
        let [lo, hi] = coords;
        // The six pairs within a half: one and two places round (three places round is one
        // place round taken from the other lane).
        let within = or(
            or(eq32(lo, rot32_1(lo)), eq32(lo, rot32_2(lo))),
            or(eq32(hi, rot32_1(hi)), eq32(hi, rot32_2(hi))),
        );
        // The sixteen pairs across the halves.
        let across = or(
            or(eq32(lo, hi), eq32(lo, rot32_1(hi))),
            or(eq32(lo, rot32_2(hi)), eq32(lo, rot32_3(hi))),
        );
        any_set(or(within, across))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (each lane's test)
    #[inline(always)]
    unsafe fn lanes_mark<const MARK: bool>(
        lane_coords: [u32; X_DIAG_LANES],
        s: Self::Lanes,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
    ) -> Option<(u32, Self::Lanes)> {
        // Only the deferral's batch asks for the marker, which this kernel never runs.
        const { assert!(!MARK, "the 128-bit kernel does not defer its calls") };
        use crate::utils::xdrop_simd::arch::*;
        let diag_offset = splat_i32(test.diag_offset);
        let subject_plus_offset = s.map(|s| add32(s, diag_offset));
        // coord <= mask < entries.len() (the caller's contract; XDiagTable::new)
        let e = |i: usize| {
            table
                .entries
                .get_unchecked(lane_coords[i] as usize)
                .raw_bits() as i32
        };
        let old = [set32(e(0), e(1), e(2), e(3)), set32(e(4), e(5), e(6), e(7))];
        // From here to the stores: x_batch_avx2_lanes_mark's test without the marker, by half.
        let mut lanes = 0u32;
        let mut lane_new = [0u32; X_DIAG_LANES];
        for half in 0..2 {
            let old = old[half];
            let flag = sra32::<31>(old);
            let last = sra32::<1>(shl32::<1>(old));
            let diff = sub32(s[half], sub32(last, diag_offset));
            let sum_below_last = gt_i32(last, subject_plus_offset[half]);
            let in_window = gt_i32(splat_i32(test.window), diff);
            let keep = blend8(in_window, sum_below_last, flag);
            let below_word = gt_i32(splat_i32(test.word_size), diff);
            let candidates = andnot(or(flag, below_word), in_window);
            let written = and(subject_plus_offset[half], splat_i32(X_LAST_HIT_MASK as i32));
            let new = blend8(written, old, keep);
            store(lane_new.as_mut_ptr().add(4 * half), new);
            lanes |= mask32(candidates) << (4 * half);
        }
        for (&coord, &raw) in lane_coords.iter().zip(&lane_new) {
            // coord <= mask < entries.len()
            table
                .entries
                .get_unchecked_mut(coord as usize)
                .set_raw_bits(raw);
        }
        Some((lanes, old))
    }

    x_cross_batch_methods!(
        #[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
    );

    // No NCBI counterpart: the kernel's `x_cross_drain_w128`, which holds the NCBI reference.
    #[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
    #[inline]
    unsafe fn drain<const DEFER: bool, const CHUNK: usize, F>(
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
        // The kernel calls each batch's candidates at once (module header): no drain defers them.
        const { assert!(!DEFER, "the 128-bit kernel does not defer its calls") };
        x_cross_drain_w128::<DEFER, CHUNK, F>(cross, staged, end, table, test, candidate, counts)
    }

    // No NCBI counterpart: the kernel's `x_cross_window_w128`, which holds the NCBI reference.
    #[inline(always)]
    unsafe fn window<const CHUNK: usize, F>(
        cross: XCross<'_, CHUNK>,
        b: usize,
        table: &mut XDiagTable<'_>,
        test: XDiagTest,
        candidate: &mut F,
    ) -> u32
    where
        F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
    {
        x_cross_window_w128(cross, b, table, test, candidate)
    }
}

/// EXPERIMENT (LOSAT_X_BXCROSS128): [`x_scan_cells_avx2_cross`] with the 128-bit kernel
/// ([`XCrossW128`]) and each unit's calls at once (the drain without `DEFER`), compiled for its
/// instructions as a whole (x86-64: SSSE3, SSE4.1 and POPCNT; wasm32: simd128, the target's).
/// Never inlined.
///
/// # Safety
/// On x86-64 the CPU must have SSSE3, SSE4.1 and POPCNT; every query offset of `lookup` must be at
/// most `table.mask` ([`x_scan_cells_cross128`]).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-509 (the scan and its hits in
// order; the rest of the loop at `x_scan_two_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
#[inline(never)]
unsafe fn x_scan_cells_w128_cross<const BRANCH_FREE: bool, const CHUNK: usize, K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XCrossCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    const DEFER: bool = false;
    x_cross_scan!(
        XCrossW128,
        BRANCH_FREE,
        DEFER,
        CHUNK,
        K,
        F,
        (
            lookup,
            subject,
            subject_length,
            table,
            test,
            candidate,
            counts
        )
    );
}

/// EXPERIMENT (LOSAT_X_BXCROSS128): [`x_cross_drain`] with the 128-bit kernel, out of line (it
/// runs once per chunk of staged hits and once per subject).
///
/// # Safety
/// As [`x_scan_cells_w128_cross`] and [`x_cross_drain`].
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
#[inline(never)]
unsafe fn x_cross_drain_w128<const DEFER: bool, const CHUNK: usize, F>(
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
    x_cross_drain_with::<XCrossW128, DEFER, CHUNK, F>(
        cross, staged, end, table, test, candidate, counts,
    )
}

/// EXPERIMENT (LOSAT_X_BXCROSS128): the 128-bit kernel's window of eight staged hits two of which
/// name one entry: its hits one after the other ([`x_cross_hit`]: the reference's order, each
/// candidate called at once), in place of the AVX2 kernel's parts ([`x_cross_window`]). A short
/// query has a small diagonal table, on which the hits of consecutive cells meet on one diagonal
/// often: against the fly proteome a 2.1-kb cDNA (6.77 hits per used cell, the fewest the fused
/// loop takes) repeats an entry in a fifth of its windows (2.1 M of 10.2 M), a 50-kb piece of the
/// genome in 2 %; there the parts' batches of two vectors each cost the 128-bit kernel more than
/// the eight scalar tests (in the WebAssembly build the parts took 0.235 G of the cDNA's 2.06 G
/// cycles, its drain 0.247 G). Returns the number of candidates. Cold and never inlined, so that
/// the drain's loop keeps none of its values.
///
/// # Safety
/// As [`x_scan_cells_w128_cross`] and [`x_cross_hit`]: `b + 8` at most the staged count.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:508-516 (the hits in order, each
// on its entry; the test of each, :518-552, is cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
#[cold]
#[inline(never)]
unsafe fn x_cross_window_w128<const CHUNK: usize, F>(
    cross: XCross<'_, CHUNK>,
    b: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let mut candidates = 0;
    for i in b..b + X_DIAG_LANES {
        candidates += x_cross_hit(cross, i, table, test, candidate);
    }
    candidates
}

/// EXPERIMENT (LOSAT_X_BPCELL128): the scan of [`XDiagMode::W128`]: each cell's hits by
/// [`x_cell_w128`], with LOSAT_X_BPHITBF's branch-free test for the hits after the cell's last
/// batch of eight where it is on and [`XHitBf::new`] finds the subject in range (once per
/// subject, as the per-hit loop does, [`x_hit_bf_for`]). Adds the cells, their hits and the
/// candidates to `counts`. Out of line, compiled for the 128-bit layer's instructions on x86-64.
///
/// # Safety
/// On x86-64 the CPU must have SSSE3, SSE4.1 and POPCNT ([`x_cross128_available`]).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-552 (the scan's cells, cited at
// `x_scan_two_hit`; the test of each hit, cited at `x_test_hit`)
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
#[inline(never)]
unsafe fn x_scan_cells_w128<const BRANCH_FREE: bool, K: XCell, F>(
    lookup: &BlastAaLookupTable<K>,
    subject: &[u8],
    subject_length: usize,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    candidate: &mut F,
    counts: &mut XDiagCounts,
) where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    let bf = if x_bp_hit_bf() {
        XHitBf::new(
            test.diag_offset,
            subject_length,
            test.window,
            test.word_size,
        )
    } else {
        None
    };
    x_scan_cells_with::<BRANCH_FREE, _, _>(lookup, subject, subject_length, |s, hits| {
        // SAFETY: this function's contract; the lookup's query offsets are at most the table's
        // mask (`XDiagTable::new` and the lookup's build).
        let candidates = unsafe { x_cell_w128(hits, s, table, test, bf, candidate) };
        counts.add_cell(hits.len(), candidates);
    });
}

/// EXPERIMENT (LOSAT_X_BPCELL128): the hits of one cell (subject offset `s`) in batches of eight,
/// each the test of eight lanes of the 128-bit kernel ([`x_batch_lanes`] with [`XCrossW128`],
/// LOSAT_X_BXCROSS128's batch, every lane at `s`: the lane formulas of [`x_batch_avx2_at`]), its
/// candidates called in lane order at once ([`x_call_candidates`], as [`x_cell_avx2`] calls them);
/// the fewer than eight hits after the last batch one after the other, with `bf`
/// ([`XHitBf::test`], then [`x_test_hit`] for a flag-1 entry or a candidate: LOSAT_X_BPHITBF's
/// per-hit loop, [`x_hit_bf`]) or without it ([`x_test_hit`]). Returns the number of candidates.
///
/// # Safety
/// Inlined into code compiled for the 128-bit layer (on x86-64 SSSE3, SSE4.1 and POPCNT); every
/// hit's query offset must be at most `table.mask`; with `bf`, [`XHitBf::fits`] must hold for the
/// subject the cell is in.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552 (cited at `x_test_hit`): the
// hits of the cell in order, each tested on its entry; a hit that passes reaches the context check
// with `last_hit` and `diff`.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[inline(always)]
unsafe fn x_cell_w128<H: XHit, F>(
    hits: &[H],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    bf: Option<XHitBf>,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    use crate::utils::xdrop_simd::arch::splat_i32;
    let mask = table.mask;
    let lanes_s = [splat_i32(s), splat_i32(s)];
    let mut batches = hits.chunks_exact(X_DIAG_LANES);
    let mut candidates = 0u32;
    for batch in &mut batches {
        let p = batch.as_ptr();
        // SAFETY: i < X_DIAG_LANES = batch.len().
        let lane_coords: [u32; X_DIAG_LANES] =
            std::array::from_fn(|i| (unsafe { (*p.add(i)).q() }.wrapping_sub(s) & mask) as u32);
        debug_assert!(x_batch_coords_distinct(&lane_coords.map(|c| c as usize)));
        // SAFETY: this function's contract (the coordinates are masked).
        let lanes = unsafe { x_batch_lanes::<XCrossW128>(lane_coords, lanes_s, table, test) };
        x_call_candidates(lanes, batch, s, table, test, candidate);
        candidates += lanes.count_ones();
    }
    let rest = batches.remainder();
    match bf {
        Some(bf) => {
            let here = bf.here(s);
            for &hit in rest {
                let query_offset = hit.q();
                let entry = table.entry_bf(query_offset, s);
                if bf.test(entry, here) {
                    continue;
                }
                if let Some((last_hit, diff)) =
                    x_test_hit(entry, s, test.diag_offset, test.window, test.word_size)
                {
                    candidate(query_offset, s, last_hit, diff, entry);
                    candidates += 1;
                }
            }
        }
        None => {
            for &hit in rest {
                let query_offset = hit.q();
                let entry = table.entry(query_offset, s);
                if let Some((last_hit, diff)) =
                    x_test_hit(entry, s, test.diag_offset, test.window, test.word_size)
                {
                    candidate(query_offset, s, last_hit, diff, entry);
                    candidates += 1;
                }
            }
        }
    }
    candidates
}

// No NCBI counterpart: a test hook; it computes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPCELL128): [`x_cell_w128`] compiled for the 128-bit layer, for the tests.
///
/// # Safety
/// [`x_cell_w128`]'s contract; on x86-64 the CPU must have SSSE3, SSE4.1 and POPCNT.
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
#[cfg(test)]
#[cfg_attr(target_arch = "x86_64", target_feature(enable = "ssse3,sse4.1,popcnt"))]
unsafe fn x_cell_w128_test<H: XHit, F>(
    hits: &[H],
    s: i32,
    table: &mut XDiagTable<'_>,
    test: XDiagTest,
    bf: Option<XHitBf>,
    candidate: &mut F,
) -> u32
where
    F: FnMut(i32, i32, i32, i32, &mut DiagStruct),
{
    // SAFETY: this function's contract.
    unsafe { x_cell_w128(hits, s, table, test, bf, candidate) }
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
unsafe fn x_cell_avx2<H: XHit, F>(
    hits: &[H],
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

// The context table (`XContextTable`) is in algorithm/common (`x_context_table.rs`), where TBLASTX,
// below blastp in the program layering, takes it too; blastp and blastx name it here.
pub(crate) use crate::algorithm::common::x_context_table::XContextTable;

#[cfg(test)]
#[path = "x_seed_fused_tests.rs"]
mod tests;
