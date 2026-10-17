//! EXPERIMENT (LOSAT_X_BXSCAN / LOSAT_X_BXSCANSHADOW): the BLASTX seed loop reads the scan
//! buffer directly, and the buffer is kept per thread.
//!
//! NCBI allocates the offset-pair array once per search thread and its word finder reads the
//! pairs straight from that array:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
//! ```c
//!     aux_struct->offset_pairs =
//!       (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:502-511
//! ```c
//!     while (scan_range[1] <= scan_range[2]) {
//!         /* scan the subject sequence for hits */
//!         hits = scansub(lookup_wrap, subject,
//!                                   offset_pairs, array_size, scan_range);
//!
//!         totalhits += hits;
//!         /* for each hit, */
//!         for (i = 0; i < hits; ++i) {
//!             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//!             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//! ```
//!
//! The reference (`Lookup::scan` + `word_finder_plain`) allocates and zero-fills a buffer of
//! `capacity()` pairs per (chunk, subject) and copies every batch into a new `(i32, i32)` vector
//! before the diagonal loop reads it. `x_word_finder_scan` is `word_finder_plain` with the
//! scan loop of `Lookup::scan` inlined: the same `scan_range` state machine, the same
//! `array_size` (`capacity()`), the scanner writing into a per-thread buffer (grown when a table
//! needs more, never shrunk), and the diagonal statements reading the pairs in place.
//!
//! Identity argument: the scanners (`s_blast_aa_scan_subject_one_range`, the compressed
//! `scan_subject`) only write the buffer and return the number of pairs written, and only
//! `pairs[0..hits)` is read afterwards, so stale contents beyond `hits` are never read. The
//! pairs are visited in the same order with the same values (`u32 as i32`, as the copy did),
//! and the diagonal statements below are those of `word_finder_plain`. The caller uses this
//! path only when nothing traces the seeds (not in the diagnostic stages) and with the
//! reference diagonal order (LOSAT_X_SEEDBUCKET mode 0).
//!
//! EXPERIMENT (LOSAT_X_BXSCANBF): with this switch (any value; off unless set, not in the strict
//! set) the scan of the standard lookup table is tblastn's branch-free scan
//! (`tblastn/x_scan_bf.rs`, `x_scan_one_range_bf`, LOSAT_X_TBNSCANBF's loop): the presence bits of
//! a block of words advance a candidate list instead of a branch, and the listed words' hits are
//! copied in position order with NCBI's stop and resume when the offset array is full. Its pairs
//! and `scan_range` values are those of `s_blast_aa_scan_subject_one_range`, so the loop above
//! reads the same pairs; it writes up to `X_SCAN_BF_SLACK` pairs past `array_size`, which the
//! per-thread buffer then holds and nothing reads. With LOSAT_X_BXDIAG the fused loop takes the
//! same candidate loop (`x_seed_diag.rs`). The compressed lookup table keeps its scan.

use super::*;
use crate::algorithm::tblastn::x_scan_bf::{x_scan_one_range_bf, X_SCAN_BF_SLACK};
use std::cell::RefCell;

/// EXPERIMENT (LOSAT_X_BXSCAN): 0 = off, 1 = direct scan loop, 2 = shadow
/// (`LOSAT_X_BXSCANSHADOW`: also the reference on a copy of the diagonal table, compared).
pub(crate) fn x_bxscan_mode() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BXSCANSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_BXSCAN") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: reads the LOSAT_X_BXSCANBF switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BXSCANBF): whether the seed scans (BXSCAN's loop and BXDIAG's fused loop)
/// find the words with hits by the branch-free candidate loop (module header). Off unless set.
pub(crate) fn x_bxscan_bf() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXSCANBF").is_some())
}

/// EXPERIMENT (LOSAT_X_BXSCANSHADOW): subjects compared so far.
pub(crate) static X_SHADOW_SUBJECTS: std::sync::atomic::AtomicU64 =
    std::sync::atomic::AtomicU64::new(0);

/// EXPERIMENT (LOSAT_X_BXLUTSHADOW / LOSAT_X_BXSCANSHADOW / LOSAT_X_BXSEGMEMOSHADOW): the
/// shadow-mode summary of the S-C switches (printed by `main` at exit, one line per switch in
/// shadow mode).
// No NCBI counterpart: prints the shadow-mode counters; it does not change any value NCBI computes.
pub fn print_shadow_stats() {
    use std::sync::atomic::Ordering;
    super::x_seed_diag::x_print_shadow_summary();
    #[cfg(feature = "parallel")]
    super::super::search::x_subject_flow::x_print_shadow_summary();
    if let Some(tables) = crate::algorithm::tblastx::lookup::x_bxlut_shadow_tables() {
        eprintln!("[X_SHADOW] BXLUT tables_compared={tables} (lookup tables identical)");
    }
    if x_bxscan_mode() == 2 {
        eprintln!(
            "[X_SHADOW] BXSCAN subjects_compared={} (hit lists and diagonal tables identical)",
            X_SHADOW_SUBJECTS.load(Ordering::Relaxed)
        );
    }
    if super::super::x_seg_memo::x_bxsegmemo_mode() == 2 {
        eprintln!(
            "[X_SHADOW] BXSEGMEMO memo_hits_compared={} (SEG intervals identical)",
            super::super::x_seg_memo::X_SHADOW_CALLS.load(Ordering::Relaxed)
        );
    }
}

thread_local! {
    /// The per-thread offset-pair array (blast_engine.c:1040-1041).
    static X_PAIRS: RefCell<Vec<OffsetPair>> = const { RefCell::new(Vec::new()) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
// ...
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
// ...
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
// ```
/// EXPERIMENT (LOSAT_X_BXSCAN): `word_finder` (mode 0) without the per-call buffer and the
/// batch copy; with the shadow switch the reference also runs on a copy of the diagonal table
/// and the hit lists and tables are compared.
pub(crate) fn x_word_finder_scan(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
) -> Result<Vec<InitHsp>> {
    // EXPERIMENT (LOSAT_X_BXDIAG): the fused loop where it acts (x_seed_diag.rs); where it does
    // not, the loop LOSAT_X_BXSCAN chooses (the reference word finder when BXSCAN is off).
    if super::x_seed_diag::x_bxdiag_mode() != 0 {
        if let Some(hits) = super::x_seed_diag::x_word_finder_diag(
            batch, parameters, options, subject, lookup, diagonals,
        )? {
            return Ok(hits);
        }
        if x_bxscan_mode() == 0 {
            return word_finder(
                batch,
                parameters,
                options,
                subject,
                lookup,
                diagonals,
                |_| {},
            );
        }
    }
    if x_bxscan_mode() == 2 {
        let mut copy = Diagonals {
            entries: diagonals.entries.clone(),
            offset: diagonals.offset,
            window: diagonals.window,
        };
        let reference = word_finder_plain(
            batch,
            parameters,
            options,
            subject,
            lookup,
            &mut copy,
            |_| {},
        )?;
        let direct = x_word_finder_direct(
            batch,
            parameters,
            options,
            subject,
            lookup,
            diagonals,
            x_bxscan_bf(),
        )?;
        assert!(
            copy.offset == diagonals.offset && copy.entries == diagonals.entries,
            "LOSAT_X_BXSCANSHADOW: diagonal table differs"
        );
        assert!(
            direct == reference,
            "LOSAT_X_BXSCANSHADOW: hit list differs ({} vs {} hits)",
            direct.len(),
            reference.len()
        );
        X_SHADOW_SUBJECTS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        return Ok(direct);
    }
    x_word_finder_direct(
        batch,
        parameters,
        options,
        subject,
        lookup,
        diagonals,
        x_bxscan_bf(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-511
// ```c
//         wordsize = lookup->word_length;
//     }
//
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
// ...
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
// ```
// `Lookup::scan` (its range set-up and loop) and `word_finder_plain` (the statements per pair)
// in one function; the pairs are read from the scanner's buffer. EXPERIMENT (LOSAT_X_BXSCANBF):
// `branch_free` scans the standard table with `x_scan_one_range_bf` (module header).
fn x_word_finder_direct(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    branch_free: bool,
) -> Result<Vec<InitHsp>> {
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let mask = diagonals.entries.len() as i32 - 1;
    let mut hits = Vec::new();
    let lazy_context = {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXLAZYCTX"))
    };
    // Lookup::scan's set-up, with length = subject.len() - 1 as word_finder_plain passes it.
    let length = subject.len().saturating_sub(1);
    let window = options.window_size;
    let lookup_word = match lookup {
        Lookup::Standard(l) => l.word_length,
        Lookup::Compressed(l) => l.word_length,
    };
    let ranges = [(0, length as i32)];
    let mut range = [0, 0, length as i32 - lookup_word];
    if window > 0 && range[2] < range[1] {
        range[2] = range[1];
    }
    let capacity = lookup.capacity();
    // EXPERIMENT (LOSAT_X_BXSCANBF): room for the pairs the branch-free scan writes past
    // `array_size` (never read).
    let room = capacity + if branch_free { X_SCAN_BF_SLACK } else { 0 };
    X_PAIRS.with(|buffer| {
        let mut buffer = buffer.borrow_mut();
        if buffer.len() < room {
            buffer.resize(room, OffsetPair::default());
        }
        let pairs = &mut buffer[..room];
        let array_size = capacity as i32;
        while range[1] <= range[2] {
            let n = match lookup {
                Lookup::Standard(l) if branch_free => {
                    x_scan_one_range_bf(l, subject, pairs, array_size, &mut range)
                }
                Lookup::Standard(l) => {
                    s_blast_aa_scan_subject_one_range(l, subject, pairs, array_size, &mut range)
                }
                Lookup::Compressed(l) => {
                    l.scan_subject(subject, &ranges, pairs, array_size, &mut range)
                }
            };
            for pair in &pairs[..n as usize] {
                let (q, s) = (pair.q_off as i32, pair.s_off as i32);
                // From here to the end of the loop body: the statements of word_finder_plain.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,783-784
                // ```c
                //                 curr_context = BSearchContextInfo(query_offset, query_info);
                // ...
                //                 Int4 curr_context = BSearchContextInfo(query_offset,
                //                                                        query_info);
                // ```
                let context_of = |q: i32| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                };
                let eager_context = if lazy_context {
                    usize::MAX
                } else {
                    context_of(q)
                };
                let context;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516,776
                // ```c
                //             diag_coord = (query_offset - subject_offset) & diag_mask;
                // ...
                //             diag_coord = (subject_offset - query_offset) & diag_mask;
                // ```
                let index = if options.window_size == 0 {
                    (s - q) & mask
                } else {
                    (q - s) & mask
                } as usize;
                let entry = &mut diagonals.entries[index];
                let (u, end, extended) = if options.window_size == 0 {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:777-781
                    // ```c
                    //             diff = subject_offset -
                    //                 (diag_array[diag_coord].last_hit - diag_offset);
                    //
                    //             /* do an extension, but only if we have not already extended this
                    //                far */
                    // ```
                    if s.wrapping_sub(entry.last_hit().wrapping_sub(diagonals.offset)) < 0 {
                        continue;
                    }
                    context = if lazy_context {
                        context_of(q)
                    } else {
                        eager_context
                    };
                    let p = &parameters[context];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:787-790
                    // ```c
                    //                 score = s_BlastAaExtendOneHit(matrix, subject, query,
                    //                                               subject_offset, query_offset,
                    //                                               cutoffs->x_dropoff,
                    //                                               &hsp_q, &hsp_s, &hsp_len,
                    // ```
                    let Some(result) = extend_one_hit_blosum62(
                        query,
                        subject,
                        q as usize,
                        s as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, true)
                } else {
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551
                    // ```c
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
                    if entry.flag() != 0 {
                        if s.wrapping_add(diagonals.offset) < entry.last_hit() {
                            continue;
                        }
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        entry.set_flag(0);
                        continue;
                    }
                    let previous = entry.last_hit().wrapping_sub(diagonals.offset);
                    let diff = s.wrapping_sub(previous);
                    if diff >= options.window_size {
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        continue;
                    }
                    if diff < word {
                        continue;
                    }
                    context = if lazy_context {
                        context_of(q)
                    } else {
                        eager_context
                    };
                    let p = &parameters[context];
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
                    // ```c
                    //                 if (query_offset - diff <
                    //                     query_info->contexts[curr_context].query_offset) {
                    //
                    //                     /* there was no last hit for this diagnol; start a new hit */
                    //                     diag_array[diag_coord].last_hit =
                    //                         subject_offset + diag_offset;
                    //                     continue;
                    //                 }
                    // ```
                    if q - diff < batch.contexts[context].offset as i32 {
                        entry.set_last_hit(s.wrapping_add(diagonals.offset));
                        continue;
                    }
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:576-582
                    // ```c
                    //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
                    //                                               last_hit + wordsize,
                    //                                               subject_offset, query_offset,
                    //                                               cutoffs->x_dropoff,
                    //                                               &hsp_q, &hsp_s,
                    //                                               &hsp_len, use_pssm,
                    //                                               wordsize, &right_extend,
                    // ```
                    let Some(result) = extend_two_hit_blosum62(
                        query,
                        subject,
                        (previous + word) as usize,
                        s as usize,
                        q as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, result.right_extend)
                };
                let p = &parameters[context];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-606
                // ```c
                //                 if (score >= cutoffs->cutoff_score)
                //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
                //                                      query_offset, subject_offset, hsp_len,
                //                                      score);
                //
                //                 /* If an extension to the right happened, reset the last hit
                //                    so that future hits to this diagonal must start over. */
                //
                //                 if (right_extend) {
                //                     diag_array[diag_coord].flag = 1;
                //                     diag_array[diag_coord].last_hit =
                //                         s_last_off - (wordsize - 1) + diag_offset;
                //                 }
                //                 /* Otherwise, make the present hit into the previous hit for
                //                    this diagonal */
                //                 else {
                //                     diag_array[diag_coord].last_hit =
                //                         subject_offset + diag_offset;
                //                 }
                // ```
                if u.score >= p.word_cutoff {
                    hits.push(InitHsp {
                        q_seed: q,
                        s_seed: s,
                        q_start: u.q_start,
                        s_start: u.s_start,
                        length: u.length,
                        score: u.score,
                    });
                }
                if options.window_size == 0 {
                    entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
                } else if extended {
                    entry.set_flag(1);
                    entry.set_last_hit((end - (word - 1)).wrapping_add(diagonals.offset));
                } else {
                    entry.set_last_hit(s.wrapping_add(diagonals.offset));
                }
            }
        }
    });
    // The rest of word_finder_plain, unchanged.
    diagonals.finish_subject(subject.len().saturating_sub(1))?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // The frozen pinned-NCBI INIT rows of the session F fixture (one-hit and two-hit windows,
    // fresh, retained and reset diagonal tables): the direct loop gives the NCBI hits and the
    // same diagonal table as the reference loop.
    /// An entry with `last_hit` and the flag set.
    fn entry(last_hit: i32) -> DiagStruct {
        let mut entry = DiagStruct::default();
        entry.set_last_hit(last_hit);
        entry.set_flag(1);
        entry
    }

    #[test]
    fn direct_scan_loop_matches_ncbi_hits_and_reference_diagonal_table() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let q = parse_fasta(&std::fs::read(dir.join("query.fna")).unwrap(), false, false).unwrap();
        let s = parse_fasta(
            &std::fs::read(dir.join("subject.faa")).unwrap(),
            true,
            false,
        )
        .unwrap();
        let mut subject = s[0]
            .sequence
            .iter()
            .copied()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect::<Vec<_>>();
        subject.push(0);
        for window in [0, 1, 39, 40, 41] {
            let args = [
                "losat",
                "blastx",
                "-query",
                "query",
                "-subject",
                "subject",
                "-seg",
                "no",
                "-comp_based_stats",
                "0",
                "-window_size",
                &window.to_string(),
            ];
            let crate::cli::Commands::Blastx(args) =
                crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                    .unwrap()
                    .command
            else {
                unreachable!()
            };
            let options = args.resolve().unwrap();
            let batch = prepare_queries(&q, &options).unwrap();
            let mut params = parameters::score_block(&batch);
            parameters::effective_lengths(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                None,
            )
            .unwrap();
            let lookup = Lookup::new(&batch, &params, &options).unwrap();
            parameters::subject_parameters(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                s[0].sequence.len(),
                None,
            )
            .unwrap();
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            let expected =
                std::fs::read_to_string(dir.join(format!("window{window}.stage"))).unwrap();
            let initial = expected
                .lines()
                .filter(|x| x.starts_with("INIT\t"))
                .map(|x| {
                    let f = x.split('\t').collect::<Vec<_>>();
                    InitHsp {
                        q_seed: f[3].parse().unwrap(),
                        s_seed: f[4].parse().unwrap(),
                        q_start: f[5].parse().unwrap(),
                        s_start: f[6].parse().unwrap(),
                        length: f[7].parse().unwrap(),
                        score: f[8].parse().unwrap(),
                    }
                })
                .collect::<Vec<_>>();
            let mut direct = Diagonals::new(length, window).unwrap();
            let mut reference = Diagonals::new(length, window).unwrap();
            // EXPERIMENT (LOSAT_X_BXSCANBF): the direct loop with the branch-free scan, on its own
            // table, compared alike.
            let mut branch_free = Diagonals::new(length, window).unwrap();
            for state in 0..4 {
                for d in [&mut direct, &mut reference, &mut branch_free] {
                    if state == 2 {
                        d.offset = 536_870_910;
                        d.entries.fill(entry(123));
                    }
                    if state == 3 {
                        d.offset = 536_870_911;
                        d.entries.fill(entry(536_870_900));
                        d.finish_subject(83).unwrap();
                    }
                }
                let hits = x_word_finder_direct(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut direct,
                    false,
                )
                .unwrap();
                let hits_bf = x_word_finder_direct(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut branch_free,
                    true,
                )
                .unwrap();
                let plain = word_finder_plain(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap();
                assert_eq!(hits, initial, "window={window} state={state}");
                assert_eq!(hits, plain, "window={window} state={state}");
                assert_eq!(hits_bf, plain, "window={window} state={state} branch-free");
                for d in [&direct, &branch_free] {
                    assert!(
                        d.offset == reference.offset && d.entries == reference.entries,
                        "window={window} state={state}"
                    );
                }
            }
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
            // ```c
            // typedef struct DiagStruct {
            //    signed int last_hit   : 31; /**< Offset of the last hit */
            //    unsigned int flag      : 1 ; /**< Reset the next extension? */
            // } DiagStruct;
            // ```
            // An offset just below 2^30 (NCBI reaches such offsets at a -window_size near 2^29,
            // from the second subject on): the stores past 2^30 wrap in the 31-bit field, and both
            // loops read the wrapped values back alike. tests/blastx_diag_bitfield.rs compares the
            // rows with NCBI 2.17.0 at such a window. (The offset is past INT4_MAX / 4, so the
            // tables are cleared at the end of the subject.)
            for d in [&mut direct, &mut reference, &mut branch_free] {
                d.offset = (1 << 30) - 64;
                d.entries.fill(DiagStruct::clear(window));
            }
            let hits = x_word_finder_direct(
                &batch,
                &params,
                &options,
                &subject,
                &lookup,
                &mut direct,
                false,
            )
            .unwrap();
            let hits_bf = x_word_finder_direct(
                &batch,
                &params,
                &options,
                &subject,
                &lookup,
                &mut branch_free,
                true,
            )
            .unwrap();
            let plain = word_finder_plain(
                &batch,
                &params,
                &options,
                &subject,
                &lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            assert_eq!(hits, plain, "window={window} wrapped");
            assert_eq!(hits_bf, plain, "window={window} wrapped branch-free");
            for d in [&direct, &branch_free] {
                assert!(
                    d.offset == reference.offset && d.entries == reference.entries,
                    "window={window} wrapped"
                );
            }
            if !initial.is_empty() {
                assert_ne!(
                    hits, initial,
                    "window={window}: the wrapped stores change the seeds"
                );
            }
        }
    }
}
