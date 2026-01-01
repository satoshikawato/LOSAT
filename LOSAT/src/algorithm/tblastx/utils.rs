//! TBLASTX - LINE-BY-LINE NCBI BLAST PORT
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:439-619
//!            ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-131

use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::sync::atomic::{AtomicI32, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::mpsc::channel;

use crate::common::{write_output, Hit};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::{lookup_protein_params, search_space::SearchSpace};
use crate::utils::dust::{DustMasker, MaskedInterval};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegMasker;

use super::args::TblastxArgs;
use super::chaining::UngappedHit;
use super::constants::{CUTOFF_E_TBLASTX, GAP_TRIGGER_BIT_SCORE, X_DROP_UNGAPPED};
use super::diagnostics::{
    diagnostics_enabled, print_summary as print_diagnostics_summary, DiagnosticCounters,
};
use super::extension::{convert_coords, extend_hit_two_hit};
use super::lookup::{
    build_ncbi_lookup, decode_kmer, encode_kmer, get_charsize, get_mask, pv_test,
    BlastAaLookupTable, NeighborLookup, AA_HITS_PER_CELL, LOOKUP_ALPHABET_SIZE,
};
use crate::utils::matrix::blosum62_score;
use super::sum_stats_linking::apply_sum_stats_even_gap_linking;
use super::translation::{generate_frames, QueryFrame};
use crate::stats::karlin::{
    bit_score as calc_bit_score, raw_score_from_bit_score, raw_score_from_evalue_with_decay,
};
use crate::stats::sum_statistics::defaults::GAP_DECAY_RATE_UNGAPPED;

// [C] typedef struct DiagStruct { Int4 last_hit; Uint1 flag; } DiagStruct;
#[derive(Clone, Copy)]
struct DiagStruct {
    last_hit: i32,
    flag: u8,
}
impl Default for DiagStruct {
    fn default() -> Self {
        Self { last_hit: 0, flag: 0 }
    }
}

/// Query/subject offset pair (NCBI `BlastOffsetPair` equivalent for the hot path).
#[derive(Clone, Copy, Default)]
struct OffsetPair {
    q_off: i32,
    s_off: i32,
}

struct WorkerState {
    tx: std::sync::mpsc::Sender<Vec<Hit>>,
    offset_pairs: Vec<OffsetPair>,
    diag_array: Vec<DiagStruct>,
    diag_offset: i32,
}

#[inline]
fn atomic_max_usize(dst: &AtomicUsize, val: usize) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
fn atomic_min_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val < cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

#[inline]
fn atomic_max_i32(dst: &AtomicI32, val: i32) {
    let mut cur = dst.load(AtomicOrdering::Relaxed);
    while val > cur {
        match dst.compare_exchange(cur, val, AtomicOrdering::Relaxed, AtomicOrdering::Relaxed) {
            Ok(_) => return,
            Err(next) => cur = next,
        }
    }
}

/// s_BlastAaScanSubject - NCBI BLAST style scan
/// Reference: blast_aascan.c:48-131
///
/// Optimized version: stop codons and invalid residues are handled by the
/// presence vector (pv_test) - buckets for k-mers containing invalid residues
/// are empty, so pv_test returns false and they are naturally skipped.
#[inline(never)]
/// NCBI-style subject sequence scan for k-mer hits.
///
/// This function scans the subject sequence for k-mer matches in the lookup table,
/// using a rolling index computation (NCBI `ComputeTableIndexIncremental`).
///
/// Optimizations applied:
/// - Unsafe array access to eliminate bounds checking in hot loop
/// - Pointer-style iteration matching NCBI's C implementation
/// - Inline presence vector test
#[inline]
fn s_blast_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3], // [C] Int4 *s_range
) -> i32 {
    let mut totalhits: i32 = 0;

    // NCBI subject masking support (masksubj.inl). In NCBI this walks
    // `subject->seq_ranges[]`; in LOSAT each translated frame is one contiguous
    // range already baked into `s_range[1..=2]`.
    while s_range[1] <= s_range[2] {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        // Fast fail: nothing to scan.
        if s_first > s_last || subject.len() < 5 {
            return totalhits;
        }

        let charsize = lookup.charsize as usize;
        let mask = lookup.mask as usize;
        
        // Cache pointers for hot loop - matches NCBI pointer-based iteration
        let pv = lookup.pv.as_ptr();
        let backbone = lookup.backbone.as_ptr();
        let overflow = lookup.overflow.as_ptr();

        // [C] index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
        // Prime the rolling index with first (word_length - 1) residues
        // SAFETY: s_first + 1 < subject.len() is guaranteed by s_first <= s_last and subject.len() >= 5
        let mut index: usize = unsafe {
            (*subject.get_unchecked(s_first) as usize) << charsize
                | (*subject.get_unchecked(s_first + 1) as usize)
        };

        // [C] for (s = s_first; s <= s_last; s++)
        let mut s = s_first;
        while s <= s_last {
            // [C] index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
            // Rolling index computation: shift in the new character
            // SAFETY: s + 2 <= s_last + 2 < subject.len() is guaranteed by loop bounds
            let new_char = unsafe { *subject.get_unchecked(s + 2) } as usize;
            index = ((index << charsize) | new_char) & mask;

            // [C] if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
            // Inline presence vector test for performance
            // SAFETY: index is masked, pv array is sized to cover all possible indices
            let pv_word = unsafe { *pv.add(index >> 6) };
            if (pv_word & (1u64 << (index & 63))) != 0 {
                // SAFETY: index is within backbone bounds (masked)
                let cell = unsafe { &*backbone.add(index) };
                let numhits = cell.num_used;

                // [C] if (numhits <= (array_size - totalhits)) { ... }
                if numhits <= array_size - totalhits {
                    let s_off = s as i32;
                    let dest_base = totalhits as usize;

                    if numhits as usize <= AA_HITS_PER_CELL {
                        // Hits in backbone cell - unroll for common cases
                        // SAFETY: dest_base + numhits <= offset_pairs.len() guaranteed by array_size check
                        unsafe {
                            let dest = offset_pairs.as_mut_ptr().add(dest_base);
                            match numhits {
                                1 => {
                                    (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                }
                                2 => {
                                    (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                }
                                3 => {
                                    (*dest) = OffsetPair { q_off: cell.entries[0], s_off };
                                    (*dest.add(1)) = OffsetPair { q_off: cell.entries[1], s_off };
                                    (*dest.add(2)) = OffsetPair { q_off: cell.entries[2], s_off };
                                }
                                _ => {}
                            }
                        }
        } else {
                        // Hits in overflow array
                        let cursor = cell.entries[0] as usize;
                        // SAFETY: overflow bounds checked during lookup table construction
                        unsafe {
                            let src = overflow.add(cursor);
                            let dest = offset_pairs.as_mut_ptr().add(dest_base);
                            for i in 0..numhits as usize {
                                (*dest.add(i)) = OffsetPair {
                                    q_off: *src.add(i),
                                    s_off,
                                };
                            }
                        }
                    }
                    totalhits += numhits;
        } else {
                    // Not enough space in the destination array; return early
                    s_range[1] = s as i32;
                    return totalhits;
                }
            }

            s += 1;
        }

        s_range[1] = s as i32;
    }

    totalhits
}

/// Lazy neighbor scan - compute neighbors dynamically during scan.
///
/// For each subject k-mer, we compute all neighboring query k-mers that would
/// produce a score >= threshold, then check if those neighbors exist in the
/// lookup table (which contains only exact query k-mers).
///
/// This trades scan-time computation for drastically reduced lookup table size.
#[inline(never)]
fn s_blast_aa_scan_subject_lazy(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    let mut totalhits: i32 = 0;
    let threshold = lookup.threshold;
    let row_max = &lookup.row_max;
    
    while s_range[1] <= s_range[2] {
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        if s_first > s_last || subject.len() < 5 {
            return totalhits;
        }

        let charsize = get_charsize();
        let mask = get_mask();
        let alphabet_size = LOOKUP_ALPHABET_SIZE;
        
        let pv = lookup.pv.as_ptr();
        let backbone = lookup.backbone.as_ptr();
        let overflow = lookup.overflow.as_ptr();

        // Rolling index for subject k-mer
        let mut index: usize = unsafe {
            (*subject.get_unchecked(s_first) as usize) << charsize
                | (*subject.get_unchecked(s_first + 1) as usize)
        };

        let mut s = s_first;
        while s <= s_last {
            let new_char = unsafe { *subject.get_unchecked(s + 2) } as usize;
            index = ((index << charsize) | new_char) & mask;
            
            // Decode subject k-mer
            let (s0, s1, s2) = decode_kmer(index);
            
            // Skip invalid residues
            if s0 >= alphabet_size || s1 >= alphabet_size || s2 >= alphabet_size {
                s += 1;
                continue;
            }

            // For lazy mode, we need to find all query k-mers Q such that
            // score(Q, S) >= threshold, where S is the subject k-mer.
            // 
            // This is the REVERSE of the standard neighbor generation:
            // - Standard: for each query k-mer Q, find all S such that score(Q, S) >= threshold
            // - Lazy: for each subject k-mer S, find all Q such that score(Q, S) >= threshold
            //
            // We enumerate all possible Q and check if:
            // 1. score(Q, S) >= threshold
            // 2. Q exists in the lookup table (pv_test)
            //
            // Pruning with row_max:
            // score(Q, S) = blosum62(q0, s0) + blosum62(q1, s1) + blosum62(q2, s2)
            // max_score_for_s = row_max[s0] + row_max[s1] + row_max[s2]
            // If max_score_for_s < threshold, no Q can match.
            
            let max_possible = row_max[s0] + row_max[s1] + row_max[s2];
            if max_possible < threshold {
                    s += 1;
                    continue;
                }
            
            // Enumerate potential query k-mers with pruning
            let rm12 = row_max[s1] + row_max[s2];
            let rm2 = row_max[s2];
            
            for q0 in 0..alphabet_size {
                let sc0 = blosum62_score(q0 as u8, s0 as u8);
                if sc0 + rm12 < threshold {
                    continue;
                }
                for q1 in 0..alphabet_size {
                    let sc1 = sc0 + blosum62_score(q1 as u8, s1 as u8);
                    if sc1 + rm2 < threshold {
                        continue;
                    }
                    for q2 in 0..alphabet_size {
                        let total_score = sc1 + blosum62_score(q2 as u8, s2 as u8);
                        if total_score < threshold {
                            continue;
                        }
                        
                        // This query k-mer Q = (q0, q1, q2) is a neighbor of subject k-mer S
                        let q_idx = encode_kmer(q0, q1, q2);
                        
                        // Check if Q exists in lookup table
                        let pv_word = unsafe { *pv.add(q_idx >> 6) };
                        if (pv_word & (1u64 << (q_idx & 63))) == 0 {
                            continue;
                        }
                        
                        // Q exists - get its hits
                        let cell = unsafe { &*backbone.add(q_idx) };
                let numhits = cell.num_used;
                        
                        if numhits <= 0 {
                            continue;
                        }
                        
                        if numhits <= array_size - totalhits {
                        let s_off = s as i32;
                            let dest_base = totalhits as usize;

                            if numhits as usize <= AA_HITS_PER_CELL {
                                unsafe {
                                    let dest = offset_pairs.as_mut_ptr().add(dest_base);
                        for i in 0..numhits as usize {
                                        (*dest.add(i)) = OffsetPair {
                                            q_off: cell.entries[i],
                                            s_off,
                                        };
                                    }
                                }
                            } else {
                                let cursor = cell.entries[0] as usize;
                                unsafe {
                                    let src = overflow.add(cursor);
                                    let dest = offset_pairs.as_mut_ptr().add(dest_base);
                                    for i in 0..numhits as usize {
                                        (*dest.add(i)) = OffsetPair {
                                            q_off: *src.add(i),
                                            s_off,
                                        };
                                    }
                                }
                        }
                        totalhits += numhits;
                    } else {
                        s_range[1] = s as i32;
                        return totalhits;
                        }
                    }
                }
            }

            s += 1;
        }

        s_range[1] = s as i32;
    }

    totalhits
}

pub fn run(args: TblastxArgs) -> Result<()> {
    // Use neighbor_map mode for faster scanning with pre-computed neighbor relationships
    if args.neighbor_map {
        return run_with_neighbor_map(args);
    }
    
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let only_qframe = args.only_qframe;
    let only_sframe = args.only_sframe;

    let valid_frame = |f: i8| matches!(f, 1 | 2 | 3 | -1 | -2 | -3);
    if let Some(f) = only_qframe {
        if !valid_frame(f) {
            anyhow::bail!("Invalid --only-qframe");
        }
    }
    if let Some(f) = only_sframe {
        if !valid_frame(f) {
            anyhow::bail!("Invalid --only-sframe");
        }
    }

    // [C] window = diag->window;
    let window = args.window_size as i32;
    // [C] wordsize = lookup->word_length;
    let wordsize: i32 = 3;
    let dropoff = X_DROP_UNGAPPED;
    let diag_enabled = diagnostics_enabled();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    eprintln!("Reading queries...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();

    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        queries_raw.iter().map(|r| masker.mask_sequence(r.seq())).collect()
    } else {
        queries_raw.iter().map(|_| Vec::new()).collect()
    };

    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| {
            let mut frames = generate_frames(r.seq(), &query_code);
            if let Some(f) = only_qframe {
                frames.retain(|x| x.frame == f);
            }
            frames
        })
        .collect();

    if args.seg {
        let seg = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        for frames in &mut query_frames {
            for frame in frames {
                if frame.aa_seq.len() >= 3 {
                    for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }

                    // NCBI BLAST query masking semantics (SEG, etc.): keep an unmasked copy
                    // (`sequence_nomask`), then overwrite masked residues in the working
                    // query sequence buffer with X.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1 kProtMask = 21;     /* X in NCBISTDAA */
                    //   query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
                    //   query_blk->sequence_nomask = query_blk->sequence_start_nomask + 1;
                    //   buffer[index] = kMaskingLetter;
                    //
                    // LOSAT uses NCBISTDAA encoding where X = 21.
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        const X_MASK_NCBISTDAA: u8 = 21; // NCBI: kProtMask = 21
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1); // keep last sentinel untouched
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Building lookup table...");
    let (lookup, contexts) = build_ncbi_lookup(
        &query_frames,
        &query_masks,
        args.threshold,
        args.include_stop_seeds,
        args.ncbi_stop_stop_score,
        args.max_hits_per_kmer,
        false, // lazy_neighbors disabled - use neighbor_map mode instead
    );

    eprintln!("Reading subjects...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    let scoring = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: i32::MAX,
        gap_extend: i32::MAX,
    };
    let params = lookup_protein_params(&scoring);

    eprintln!(
        "Searching {} queries vs {} subjects...",
        queries_raw.len(),
        subjects_raw.len()
    );

    let bar = ProgressBar::new(subjects_raw.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    let (tx, rx) = channel::<Vec<Hit>>();
    let out_path = args.out.clone();
    let evalue_threshold = args.evalue;

    let writer = std::thread::spawn(move || -> Result<()> {
        let mut all: Vec<Hit> = Vec::new();
        while let Ok(h) = rx.recv() {
            all.extend(h);
        }
        all.retain(|h| h.e_value <= evalue_threshold);
        all.sort_by(|a, b| {
            b.bit_score
                .partial_cmp(&a.bit_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        write_output(&all, out_path.as_ref())?;
        Ok(())
    });

    // [C] diag_mask = diag->diag_mask; (power of 2 - 1)
    let max_q: i32 = *lookup.frame_bases.last().unwrap_or(&0)
        + contexts
            .last()
            .map(|c| c.aa_seq.len() as i32)
            .unwrap_or(0);
    let max_s: i32 = subjects_raw
        .iter()
        .map(|r| r.seq().len() as i32 / 3 + 10)
        .max()
        .unwrap_or(0);
    let mut diag_array_size = 1i32;
    while diag_array_size < (max_q + max_s + window + 1) {
        diag_array_size <<= 1;
    }
    let diag_mask = diag_array_size - 1;

    // [C] array_size for offset_pairs
    // NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain
    // Reference: ncbi-blast/c++/include/algo/blast/core/lookup_wrap.h + lookup_wrap.c
    const OFFSET_ARRAY_SIZE: i32 = 4096;
    let offset_array_size: i32 = OFFSET_ARRAY_SIZE + lookup.longest_chain.max(0);

    let lookup_ref = &lookup;
    let contexts_ref = &contexts;
    let query_ids_ref = &query_ids;

    subjects_raw.par_iter().enumerate().for_each_init(
        || WorkerState {
            tx: tx.clone(),
            offset_pairs: vec![OffsetPair::default(); offset_array_size as usize],
            diag_array: vec![DiagStruct::default(); diag_array_size as usize],
            diag_offset: max_q + max_s,
        },
        |st, (s_idx, s_rec)| {
            let mut s_frames = generate_frames(s_rec.seq(), &db_code);
            if let Some(f) = only_sframe {
                s_frames.retain(|x| x.frame == f);
            }

            let s_id = s_rec
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            let s_len = s_rec.seq().len();

            // [C] BlastOffsetPair *offset_pairs
            let offset_pairs = &mut st.offset_pairs;

            // [C] DiagStruct *diag_array = diag->hit_level_array;
            let diag_array = &mut st.diag_array;

            // [C] diag_offset = diag->offset;  (per-thread, reused across subjects)
            let mut diag_offset: i32 = st.diag_offset;

            // Precompute per-(context, subject-frame) cutoff scores.
            // This is a hot-path replacement for NCBI's `cutoffs = word_params->cutoffs + curr_context`.
            //
            // We do this once per subject record to avoid per-extension hashing/branching.
            let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; contexts_ref.len()];
            for (ctx_idx, ctx) in contexts_ref.iter().enumerate() {
                for (sf_idx, sf) in s_frames.iter().enumerate() {
                    let s_aa_len = sf.aa_len;
                    let ss = SearchSpace {
                        effective_query_len: ctx.aa_len.min(s_aa_len) as f64,
                        effective_db_len: s_aa_len as f64,
                        effective_space: (ctx.aa_len.min(s_aa_len) * s_aa_len) as f64,
                        length_adjustment: 0,
                    };
                    let c = raw_score_from_evalue_with_decay(
                        CUTOFF_E_TBLASTX,
                        &params,
                        &ss,
                        true,
                        GAP_DECAY_RATE_UNGAPPED,
                    );
                    let g = raw_score_from_bit_score(GAP_TRIGGER_BIT_SCORE, &params);
                    cutoff_scores[ctx_idx][sf_idx] = c.min(g as i32);
                }
            }
            let mut ungapped_hits: Vec<UngappedHit> = Vec::new();

            for (s_f_idx, s_frame) in s_frames.iter().enumerate() {
                let subject = &s_frame.aa_seq;
                let s_aa_len = s_frame.aa_len;
                if subject.len() < 5 {
                    continue;
                }

                // [C] scan_range[1] = subject->seq_ranges[0].left;
                // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                let mut scan_range: [i32; 3] = [0, 1, (subject.len() - 4) as i32];

                // [C] while (scan_range[1] <= scan_range[2])
                while scan_range[1] <= scan_range[2] {
                    let prev_scan_left = scan_range[1];
                    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                    let hits = s_blast_aa_scan_subject(
                        lookup_ref,
                        subject,
                        offset_pairs,
                        offset_array_size,
                        &mut scan_range,
                    );

                    if diag_enabled && hits > 0 {
                        diagnostics
                            .base
                            .kmer_matches
                            .fetch_add(hits as usize, AtomicOrdering::Relaxed);
                    }

                    if hits == 0 && scan_range[1] == prev_scan_left {
                        // Safety guard: with correct NCBI-sized offset arrays, this should not happen.
                        // If it does, breaking avoids an infinite loop.
                        break;
                    }

                    // [C] for (i = 0; i < hits; ++i)
                    // OPTIMIZATION: Use raw pointers to eliminate bounds checking in hot loop
                    let diag_ptr = diag_array.as_mut_ptr();
                    let offset_pairs_ptr = offset_pairs.as_ptr();
                    
                    for i in 0..hits as usize {
                        // SAFETY: i < hits, and hits <= offset_array_size (checked by scan)
                        let pair = unsafe { &*offset_pairs_ptr.add(i) };
                        let query_offset = pair.q_off;
                        let subject_offset = pair.s_off;

                        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                        let diag_coord = ((query_offset - subject_offset) & diag_mask) as usize;
                        
                        // SAFETY: diag_coord is masked by diag_mask, which is < diag_array.len()
                        let diag_entry = unsafe { &mut *diag_ptr.add(diag_coord) };

                        // [C] if (diag_array[diag_coord].flag)
                        if diag_entry.flag != 0 {
                            // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                            if subject_offset + diag_offset < diag_entry.last_hit {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_masked
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                continue;
                            }
                            // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            // [C] diag_array[diag_coord].flag = 0;
                            diag_entry.last_hit = subject_offset + diag_offset;
                            diag_entry.flag = 0;
                        }
                        // [C] else
                        else {
                            // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                            let last_hit = diag_entry.last_hit - diag_offset;
                            // [C] diff = subject_offset - last_hit;
                            let diff = subject_offset - last_hit;

                            // [C] if (diff >= window)
                            if diff >= window {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_too_far
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }

                            // [C] if (diff < wordsize)
                            if diff < wordsize {
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .seeds_second_hit_overlap
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                continue;
                            }

                            // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                            let ctx_idx = lookup_ref.get_context_idx(query_offset);
                            let ctx = unsafe { contexts_ref.get_unchecked(ctx_idx) };
                            let q_raw = (query_offset - ctx.frame_base) as usize;
                            // NCBI uses masked sequence for extension (same as neighbor-map mode)
                            let query = &ctx.aa_seq;

                            // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                            if query_offset - diff < ctx.frame_base {
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }

                            if diag_enabled {
                                diagnostics
                                    .base
                                    .seeds_second_hit_window
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                diagnostics
                                    .base
                                    .seeds_passed
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                            }

                            // [C] cutoffs = word_params->cutoffs + curr_context;
                            let cutoff = unsafe { *cutoff_scores.get_unchecked(ctx_idx).get_unchecked(s_f_idx) };

                            // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
                            //                                   last_hit + wordsize, subject_offset, query_offset, ...)
                            // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                            let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
                                extend_hit_two_hit(
                                    query,
                                    subject,
                                    (last_hit + wordsize) as usize, // s_left_off (end of first hit word)
                                    subject_offset as usize, // s_right_off (second hit start)
                                    q_raw as usize,          // q_right_off (second hit start, local)
                                    dropoff,                 // x_dropoff (raw)
                                );

                            let hsp_q: i32 = hsp_q_u as i32;
                            let hsp_s: i32 = hsp_s_u as i32;
                            let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
                            let s_last_off: i32 = s_last_off_u as i32;

                            if diag_enabled {
                                diagnostics
                                    .base
                                    .ungapped_extensions
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                if right_extend {
                                    diagnostics
                                        .base
                                        .ungapped_two_hit_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                } else {
                                    diagnostics
                                        .base
                                        .ungapped_one_hit_extensions
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                if hsp_len > 0 {
                                    diagnostics
                                        .base
                                        .extension_total_length
                                        .fetch_add(hsp_len as usize, AtomicOrdering::Relaxed);
                                    atomic_max_usize(
                                        &diagnostics.base.extension_max_length,
                                        hsp_len as usize,
                                    );
                                }
                            }

                            // [C] if (right_extend)
                            if right_extend {
                                // [C] diag_array[diag_coord].flag = 1;
                                // [C] diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                                diag_entry.flag = 1;
                                diag_entry.last_hit =
                                    s_last_off - (wordsize - 1) + diag_offset;
                                if diag_enabled {
                                    diagnostics
                                        .base
                                        .mask_updates
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                            } else {
                                // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                diag_entry.last_hit = subject_offset + diag_offset;
                            }

                            // [C] if (score >= cutoffs->cutoff_score)
                            if score >= cutoff {
                                if diag_enabled {
                                    diagnostics
                                        .ungapped_only_hits
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                let qs = hsp_q as usize;
                                let qe = (hsp_q + hsp_len) as usize;
                                let ss = hsp_s as usize;
                                let se = (hsp_s + hsp_len) as usize;

                                // Convert to logical coords (subtract sentinel)
                                let (qs_l, qe_l) = (qs.saturating_sub(1), qe.saturating_sub(1));
                                let (ss_l, se_l) = (ss.saturating_sub(1), se.saturating_sub(1));

                                // Hot-path record: postpone string clones, identity, and bit/evalue until after
                                // sum-statistics linking + final filtering.
                                ungapped_hits.push(UngappedHit {
                                    q_idx: ctx.q_idx,
                                    s_idx: s_idx as u32,
                                    ctx_idx,
                                    s_f_idx,
                                    q_frame: ctx.frame,
                                    s_frame: s_frame.frame,
                                    q_aa_start: qs_l,
                                    q_aa_end: qe_l,
                                    s_aa_start: ss_l,
                                    s_aa_end: se_l,
                                    q_orig_len: ctx.orig_len,
                                    s_orig_len: s_len,
                                    raw_score: score,
                                    e_value: f64::INFINITY,
                                });
                            } else if diag_enabled {
                                diagnostics
                                    .ungapped_cutoff_failed
                                    .fetch_add(1, AtomicOrdering::Relaxed);
                                atomic_min_i32(&diagnostics.ungapped_cutoff_failed_min_score, score);
                                atomic_max_i32(&diagnostics.ungapped_cutoff_failed_max_score, score);
                            }
                        }
                    }
                }

                // [C] Blast_ExtendWordExit(ewp, subject->length);
                diag_offset += s_aa_len as i32;
            }

            if !ungapped_hits.is_empty() {
                if diag_enabled {
                    diagnostics
                        .base
                        .hsps_before_chain
                        .fetch_add(ungapped_hits.len(), AtomicOrdering::Relaxed);
                }
                let linked = apply_sum_stats_even_gap_linking(ungapped_hits, &params);
                if diag_enabled {
                    diagnostics
                        .base
                        .hsps_after_chain
                        .fetch_add(linked.len(), AtomicOrdering::Relaxed);
                }
                let mut final_hits: Vec<Hit> = Vec::new();
                for h in linked {
                    if h.e_value > evalue_threshold {
                        if diag_enabled {
                            diagnostics
                                .ungapped_evalue_failed
                                .fetch_add(1, AtomicOrdering::Relaxed);
                        }
                        continue;
                    }
                    if diag_enabled {
                        diagnostics
                            .ungapped_evalue_passed
                            .fetch_add(1, AtomicOrdering::Relaxed);
                    }

                    let ctx = &contexts_ref[h.ctx_idx];
                    let s_frame = &s_frames[h.s_f_idx];

                    // Compute identity/mismatch on demand (final hits only)
                    // NCBI computes identities using the *unmasked* query sequence buffer
                    // (`sequence_nomask`), even when the working query sequence was masked
                    // (e.g., SEG) to influence seeding/extension.
                    //
                    // NCBI reference (verbatim):
                    //   const Uint1* query_nomask = query_blk->sequence_nomask +
                    //       query_info->contexts[context].query_offset;
                    let len = h.q_aa_end.saturating_sub(h.q_aa_start);
                    let q0 = h.q_aa_start + 1; // +1 for sentinel
                    let s0 = h.s_aa_start + 1; // +1 for sentinel
                    let q_seq_nomask: &[u8] = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
                    let mut matches = 0usize;
                    for k in 0..len {
                        if q_seq_nomask[q0 + k] == s_frame.aa_seq[s0 + k] {
                            matches += 1;
                        }
                    }
                    let mismatch = len.saturating_sub(matches);
                    let identity = if len > 0 {
                        (matches as f64 / len as f64) * 100.0
                    } else {
                        0.0
                    };

                    let bit = calc_bit_score(h.raw_score, &params);
                    let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, ctx.frame, ctx.orig_len);
                    let (s_start, s_end) =
                        convert_coords(h.s_aa_start, h.s_aa_end, s_frame.frame, s_len);

                    final_hits.push(Hit {
                        query_id: query_ids_ref[ctx.q_idx as usize].clone(),
                        subject_id: s_id.clone(),
                        identity,
                        length: len,
                        mismatch,
                        gapopen: 0,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        e_value: h.e_value,
                        bit_score: bit,
                    });
                }

                if !final_hits.is_empty() {
                    if diag_enabled {
                        diagnostics
                            .base
                            .hsps_after_overlap_filter
                            .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                        diagnostics
                            .output_from_ungapped
                            .fetch_add(final_hits.len(), AtomicOrdering::Relaxed);
                    }
                    st.tx.send(final_hits).unwrap();
                }
            }
            bar.inc(1);
            st.diag_offset = diag_offset;
        },
    );

    // Close channel so the writer can exit.
    drop(tx);

    bar.finish();
    writer.join().unwrap()?;
    if diag_enabled {
        print_diagnostics_summary(&diagnostics);
    }
    Ok(())
}

/// Run TBLASTX with pre-computed neighbor map using subject-side indexing.
/// This approach:
/// 1. Index query k-mers: query_lookup[kmer] = [(q_idx, f_idx, pos), ...]
/// 2. Pre-compute neighbor relationships once
/// 3. For each subject k-mer, find all matching query positions via expanded lookup
/// 4. Apply NCBI-style sum_stats_even_gap_linking for HSP merging
///
/// Two-hit/diag logic is ported verbatim from NCBI aa_ungapped.c:s_BlastAaWordFinder_TwoHit
fn run_with_neighbor_map(args: TblastxArgs) -> Result<()> {
    use super::lookup::LOOKUP_ALPHABET_SIZE;

    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    // NCBI: window = diag->window
    let window: i32 = args.window_size as i32;
    // NCBI: wordsize = lookup->word_length (always 3 for protein/tblastx)
    let wordsize: i32 = 3;
    let threshold = args.threshold;
    let x_drop = X_DROP_UNGAPPED;
    let diag_enabled = diagnostics_enabled();

    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("Failed to build thread pool")?;

    eprintln!("Reading queries (neighbor map mode)...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw
        .iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();

    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        queries_raw.iter().map(|r| masker.mask_sequence(r.seq())).collect()
    } else {
        queries_raw.iter().map(|_| Vec::new()).collect()
    };

    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();

    // Apply SEG masking if enabled
    if args.seg {
        let seg_masker = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        const X_MASK_NCBISTDAA: u8 = 21;
        for frames in query_frames.iter_mut() {
            for frame in frames.iter_mut() {
                if frame.aa_seq.len() >= 3 {
                    for m in seg_masker.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len() - 1]) {
                        frame.seg_masks.push((m.start, m.end));
                    }
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1);
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBISTDAA;
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Building neighbor lookup...");
    let neighbor_lookup = NeighborLookup::build(&query_frames, &query_masks, threshold);
    
    // Use compressed neighbor index: no expanded_lookup pre-computation
    // Instead, resolve neighbors on-the-fly during scan
    eprintln!("Using compressed neighbor index (no expanded_lookup)...");

    eprintln!("Reading subjects...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    let scoring = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: i32::MAX,
        gap_extend: i32::MAX,
    };
    let params = lookup_protein_params(&scoring);

    eprintln!(
        "Searching {} queries vs {} subjects (neighbor map mode)...",
        queries_raw.len(),
        subjects_raw.len()
    );

    let bar = ProgressBar::new(subjects_raw.len() as u64 * 6);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    let evalue_threshold = args.evalue;
    let query_frames_ref = &query_frames;
    let query_ids_ref = &query_ids;
    let neighbor_map_ref = &neighbor_lookup.neighbor_map.map;
    let query_lookup_ref = &neighbor_lookup.query_lookup;
    let params_ref = &params;

    // Collect UngappedHit for sum_stats_linking
    // Key: (q_idx, s_idx)
    let all_ungapped: std::sync::Mutex<Vec<UngappedHit>> = std::sync::Mutex::new(Vec::new());

    // Process each subject
    for (s_idx, s_rec) in subjects_raw.iter().enumerate() {
        let s_id = s_rec
            .id()
            .split_whitespace()
            .next()
            .unwrap_or("unknown")
            .to_string();
        let s_len = s_rec.seq().len();
        let s_frames = generate_frames(s_rec.seq(), &db_code);
        
        // Precompute per-(query_frame, subject_frame) cutoff scores
        // This matches NCBI's `cutoffs = word_params->cutoffs + curr_context`
        let mut cutoff_scores: Vec<Vec<i32>> = vec![vec![0; s_frames.len()]; query_frames.iter().map(|f| f.len()).sum()];
        let mut ctx_idx = 0;
        for q_frames in query_frames.iter() {
            for q_frame in q_frames.iter() {
                for (sf_idx, sf) in s_frames.iter().enumerate() {
                    let q_aa_len = q_frame.aa_len;
                    let s_aa_len = sf.aa_len;
                    let ss = SearchSpace {
                        effective_query_len: q_aa_len.min(s_aa_len) as f64,
                        effective_db_len: s_aa_len as f64,
                        effective_space: (q_aa_len.min(s_aa_len) * s_aa_len) as f64,
                        length_adjustment: 0,
                    };
                    let c = raw_score_from_evalue_with_decay(
                        CUTOFF_E_TBLASTX,
                        &params,
                        &ss,
                        true,
                        GAP_DECAY_RATE_UNGAPPED,
                    );
                    let g = raw_score_from_bit_score(GAP_TRIGGER_BIT_SCORE, &params);
                    cutoff_scores[ctx_idx][sf_idx] = c.min(g as i32);
                }
                ctx_idx += 1;
            }
        }
        let cutoff_scores_ref = &cutoff_scores;

        // Process each subject frame in parallel
        let frame_hits: Vec<UngappedHit> = s_frames
            .par_iter()
            .enumerate()
            .flat_map(|(s_f_idx, s_frame)| {
                let s_aa = &s_frame.aa_seq;
                if s_aa.len() < 5 {
                    bar.inc(1);
                    return Vec::new();
                }

                let mut local_hits: Vec<UngappedHit> = Vec::new();
                let s_aa_len = s_aa.len();
                
                // Count total query contexts for Vec-based diagonal tracking
                let total_q_contexts: usize = query_frames_ref.iter().map(|f| f.len()).sum();
                let max_q_aa_len = query_frames_ref.iter()
                    .flat_map(|frames| frames.iter().map(|f| f.aa_seq.len()))
                    .max()
                    .unwrap_or(1);
                
                // Compute diag_array dimensions
                let diag_array_len = s_aa_len + max_q_aa_len + 1;
                let diag_offset = max_q_aa_len as i32;
                
                // Build context index mapping: (q_idx, q_f_idx) -> flat index
                let mut ctx_base: Vec<usize> = Vec::with_capacity(query_frames_ref.len() + 1);
                ctx_base.push(0);
                for frames in query_frames_ref.iter() {
                    ctx_base.push(ctx_base.last().unwrap() + frames.len());
                }
                
                // NCBI-style DiagStruct arrays: [ctx_flat_idx][diag] = DiagStruct { last_hit, flag }
                // NCBI: last_hit stores (subject_offset + diag_offset) to avoid subtraction in hot path
                // NCBI: flag=1 means we just extended, flag=0 means ready for two-hit check
                let mut diag_array: Vec<Vec<DiagStruct>> = vec![vec![DiagStruct::default(); diag_array_len]; total_q_contexts];
                
                // NCBI aa_ungapped.c:502-610 - main scan loop
                for s_pos in 1..=(s_aa_len.saturating_sub(4)) {
                    let c0 = unsafe { *s_aa.get_unchecked(s_pos) } as usize;
                    let c1 = unsafe { *s_aa.get_unchecked(s_pos + 1) } as usize;
                    let c2 = unsafe { *s_aa.get_unchecked(s_pos + 2) } as usize;

                    if c0 >= LOOKUP_ALPHABET_SIZE || c1 >= LOOKUP_ALPHABET_SIZE || c2 >= LOOKUP_ALPHABET_SIZE {
                        continue;
                    }

                    let s_kmer = encode_kmer(c0, c1, c2);
                    
                    // Compressed neighbor index: resolve neighbors on-the-fly
                    let neighbor_kmers = &neighbor_map_ref[s_kmer];
                    if neighbor_kmers.is_empty() {
                        continue;
                    }
                    
                    // subject_offset in NCBI terms (0-based position in subject AA sequence)
                    let subject_offset: i32 = (s_pos - 1) as i32;
                    
                    // For each neighbor query k-mer, get query positions
                    for &neighbor_kmer in neighbor_kmers {
                        let query_hits = &query_lookup_ref[neighbor_kmer as usize];
                        for &(q_idx, q_f_idx, q_pos) in query_hits {
                            let ctx_flat = ctx_base[q_idx as usize] + q_f_idx as usize;
                            
                            // query_offset in NCBI terms (0-based position in query AA sequence)
                            let query_offset: i32 = q_pos as i32;
                            
                            // NCBI: diag_coord = (query_offset - subject_offset) & diag_mask
                            // We use a Vec so we compute a direct index instead of masking
                            let diag_coord = (query_offset - subject_offset + diag_offset) as usize;
                            
                            if diag_coord >= diag_array_len {
                                continue;
                            }
                            
                            let diag_entry = &mut diag_array[ctx_flat][diag_coord];
                            
                            // NCBI aa_ungapped.c:519-530 - flag==1 block
                            // "If the reset bit is set, an extension just happened."
                            if diag_entry.flag != 0 {
                                // "If we've already extended past this hit, skip it."
                                if subject_offset + diag_offset < diag_entry.last_hit {
                                    continue;
                                }
                                // "Otherwise, start a new hit." - reset and CONTINUE (don't extend this hit)
                                diag_entry.last_hit = subject_offset + diag_offset;
                                diag_entry.flag = 0;
                                continue; // NCBI: after reset, move to next hit, don't extend
                            }
                            
                            // NCBI aa_ungapped.c:533-606 - flag==0 block
                            // "If the reset bit is cleared, try to start an extension."
                            
                            // NCBI: last_hit = diag_array[diag_coord].last_hit - diag_offset
                            let last_hit = diag_entry.last_hit - diag_offset;
                            // NCBI: diff = subject_offset - last_hit
                            let diff = subject_offset - last_hit;
                            
                            // NCBI aa_ungapped.c:538-544: "if (diff >= window)"
                            if diff >= window {
                                // "We are beyond the window for this diagonal; start a new hit"
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:549-551: "if (diff < wordsize)"
                            if diff < wordsize {
                                // "If the difference is less than the wordsize (i.e. last hit and this hit overlap), give up"
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:560-573: Check if last hit is in current query context
                            // "if (query_offset - diff < query_info->contexts[curr_context].query_offset)"
                            // In neighbor-map mode, each (q_idx, q_f_idx) is its own context with frame_base=0
                            // So we check: query_offset - diff < 0
                            if query_offset - diff < 0 {
                                // "there was no last hit for this diagonal; start a new hit"
                                diag_entry.last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // PASSED all two-hit checks - now extend
                            let q_frame = &query_frames_ref[q_idx as usize][q_f_idx as usize];
                            let q_aa = &q_frame.aa_seq;
                            
                            // Convert to raw positions (+1 for sentinel)
                            let q_raw = (query_offset + 1) as usize;
                            let s_raw = (subject_offset + 1) as usize;
                            
                            if q_raw + 2 >= q_aa.len() || s_raw + 2 >= s_aa.len() {
                                continue;
                            }
                            
                            // NCBI aa_ungapped.c:576-583: s_BlastAaExtendTwoHit
                            // s_left_off = last_hit + wordsize (end of first hit word)
                            // s_right_off = subject_offset (second hit start)
                            // q_right_off = query_offset (second hit start)
                            let s_left_off = (last_hit + wordsize) as usize;
                            let (hsp_q, hsp_qe, hsp_s, _hsp_se, score, right_extend, s_last_off) =
                                extend_hit_two_hit(
                                    q_aa,
                                    s_aa,
                                    s_left_off,      // s_left_off (end of first hit word, raw coords)
                                    s_raw,           // s_right_off (second hit start, raw coords)
                                    q_raw,           // q_right_off (second hit start, raw coords)
                                    x_drop,
                                );
                            
                            let hsp_len = hsp_qe.saturating_sub(hsp_q);
                            
                            // NCBI aa_ungapped.c:596-606: Update diag_array after extension
                            if right_extend {
                                // "If an extension to the right happened, reset the last hit so that
                                //  future hits to this diagonal must start over."
                                diag_entry.flag = 1;
                                diag_entry.last_hit = (s_last_off as i32) - (wordsize - 1) + diag_offset;
                            } else {
                                // "Otherwise, make the present hit into the previous hit for this diagonal"
                                diag_entry.last_hit = subject_offset + diag_offset;
                            }
                            
                            // NCBI aa_ungapped.c:588-591: Check score threshold
                            // cutoffs = word_params->cutoffs + curr_context
                            // if (score >= cutoffs->cutoff_score) ...
                            let cutoff = cutoff_scores_ref[ctx_flat][s_f_idx];
                            if score < cutoff {
                                continue;
                            }
                            
                            // Convert to logical coordinates (subtract sentinel)
                            let qs_logical = hsp_q.saturating_sub(1);
                            let qe_logical = hsp_qe.saturating_sub(1);
                            let ss_logical = hsp_s.saturating_sub(1);
                            let se_logical = ss_logical + hsp_len;
                            
                            if hsp_len == 0 {
                                continue;
                            }
                            
                            // Debug: print details for top-scoring hits (diagnostics only)
                            if diag_enabled && score > 9000 && hsp_len > 1900 {
                                eprintln!("[DEBUG TOP HIT] score={} len={} q_frame={} s_frame={}", 
                                    score, hsp_len, q_frame.frame, s_frame.frame);
                                eprintln!("  coords: q_start={} q_end={} s_start={} s_end={}", 
                                    qs_logical, qe_logical, ss_logical, se_logical);
                                eprintln!("  seed positions: s_left_off={} s_raw={}", s_left_off, s_raw);
                                
                                // Show what's at the boundaries (use raw coords since q_aa has sentinels)
                                // q_aa[hsp_q] = logical position qs_logical = hsp_q - 1
                                if hsp_q > 0 && hsp_q < q_aa.len() {
                                    let q_at_start = q_aa[hsp_q];  // Raw position = logical + 1
                                    let q_before_start = if hsp_q > 1 { q_aa[hsp_q - 1] } else { 0 };
                                    eprintln!("  q_aa[raw {}]={}  q_aa[raw-1 {}]={} (X=21)", 
                                        hsp_q, q_at_start, hsp_q - 1, q_before_start);
                                }
                                if hsp_qe > 0 && hsp_qe < q_aa.len() {
                                    let q_at_end = q_aa[hsp_qe - 1];  // Last position in alignment
                                    let q_after_end = if hsp_qe < q_aa.len() { q_aa[hsp_qe] } else { 0 };
                                    eprintln!("  q_aa[raw end-1 {}]={}  q_aa[raw end {}]={} (X=21)", 
                                        hsp_qe - 1, q_at_end, hsp_qe, q_after_end);
                                }
                                
                                // Show masked (X=21) positions near the boundaries
                                let left_boundary = qs_logical.saturating_sub(30);
                                let right_boundary = (qe_logical + 30).min(q_aa.len().saturating_sub(1));
                                let left_xs: Vec<usize> = (left_boundary..=qs_logical)
                                    .filter(|&i| i < q_aa.len() && q_aa[i] == 21)
                                    .collect();
                                let right_xs: Vec<usize> = (qe_logical..=right_boundary)
                                    .filter(|&i| i < q_aa.len() && q_aa[i] == 21)
                                    .collect();
                                eprintln!("  X's near left boundary [{}-{}]: {:?}", left_boundary, qs_logical, left_xs);
                                eprintln!("  X's near right boundary [{}-{}]: {:?}", qe_logical, right_boundary, right_xs);
                            }
                            
                            local_hits.push(UngappedHit {
                                q_idx,
                                s_idx: s_idx as u32,
                                ctx_idx: q_f_idx as usize,
                                s_f_idx,
                                q_frame: q_frame.frame,
                                s_frame: s_frame.frame,
                                q_aa_start: qs_logical,
                                q_aa_end: qe_logical,
                                s_aa_start: ss_logical,
                                s_aa_end: se_logical,
                                q_orig_len: q_frame.orig_len,
                                s_orig_len: s_len,
                                raw_score: score,
                                e_value: f64::INFINITY, // Will be computed by linking
                            });
                        } // end for query_hits
                    } // end for neighbor_kmers
                } // end for s_pos

                bar.inc(1);
                local_hits
            })
            .collect();

        // Add to global collection
        all_ungapped.lock().unwrap().extend(frame_hits);
    }

    bar.finish();

    let all_ungapped = all_ungapped.into_inner().unwrap();
    if diag_enabled {
        eprintln!("=== Stage Counters ===");
        eprintln!("[1] Raw ungapped hits (after extension): {}", all_ungapped.len());
    }

    // Apply sum-stats linking (assigns E-values to chains)
    let linked_hits = apply_sum_stats_even_gap_linking(all_ungapped, &params);
    if diag_enabled {
        eprintln!("[2] After sum_stats_linking: {} hits", linked_hits.len());

        // Count E-value distribution before final filter
        let mut e_0 = 0usize;      // E <= 0 (0.0)
        let mut e_tiny = 0usize;   // 0 < E <= 1e-50
        let mut e_small = 0usize;  // 1e-50 < E <= 1e-10
        let mut e_med = 0usize;    // 1e-10 < E <= 10
        let mut e_large = 0usize;  // E > 10
        let mut e_inf = 0usize;    // E == INF
        for h in &linked_hits {
            if h.e_value == f64::INFINITY {
                e_inf += 1;
            } else if h.e_value <= 0.0 {
                e_0 += 1;
            } else if h.e_value <= 1e-50 {
                e_tiny += 1;
            } else if h.e_value <= 1e-10 {
                e_small += 1;
            } else if h.e_value <= 10.0 {
                e_med += 1;
            } else {
                e_large += 1;
            }
        }
        eprintln!("[3] E-value distribution (before final filter):");
        eprintln!(
            "    E=0: {}, E<=1e-50: {}, E<=1e-10: {}, E<=10: {}, E>10: {}, E=INF: {}",
            e_0, e_tiny, e_small, e_med, e_large, e_inf
        );
    }

    // Pre-generate subject frames to avoid redundant translation
    if diag_enabled {
        eprintln!("Pre-generating subject frames for identity calculation...");
    }
    let subject_frames_cache: Vec<Vec<QueryFrame>> = subjects_raw
        .iter()
        .map(|s| generate_frames(s.seq(), &db_code))
        .collect();
    
    // Convert UngappedHit to Hit
    let mut final_hits: Vec<Hit> = Vec::with_capacity(linked_hits.len());
    
    for h in linked_hits {
        if h.e_value > evalue_threshold {
            continue;
        }
        
        let q_id = &query_ids_ref[h.q_idx as usize];
        let s_id = subjects_raw[h.s_idx as usize]
            .id()
            .split_whitespace()
            .next()
            .unwrap_or("unknown");
        
        let q_frame_obj = &query_frames_ref[h.q_idx as usize][h.ctx_idx];
        let q_aa = q_frame_obj.aa_seq_nomask.as_ref().unwrap_or(&q_frame_obj.aa_seq);
        
        // Get subject frame from cache (no redundant translation)
        let s_frame_obj = &subject_frames_cache[h.s_idx as usize][h.s_f_idx];
        let s_aa = &s_frame_obj.aa_seq;
        
        let len = h.q_aa_end.saturating_sub(h.q_aa_start);
        
        // Calculate identity using raw positions (+1 for sentinel)
        let q0 = h.q_aa_start + 1;
        let s0 = h.s_aa_start + 1;
        let mut matches = 0;
        for k in 0..len {
            if q0 + k < q_aa.len() && s0 + k < s_aa.len() && q_aa[q0 + k] == s_aa[s0 + k] {
                matches += 1;
            }
        }
        let identity = if len > 0 {
            (matches as f64 / len as f64) * 100.0
        } else {
            0.0
        };

        let bit_score = calc_bit_score(h.raw_score, &params);
        let (q_start, q_end) = convert_coords(h.q_aa_start, h.q_aa_end, h.q_frame, h.q_orig_len);
        let (s_start, s_end) = convert_coords(h.s_aa_start, h.s_aa_end, h.s_frame, h.s_orig_len);

        final_hits.push(Hit {
            query_id: q_id.clone(),
            subject_id: s_id.to_string(),
            identity,
            length: len,
            mismatch: len.saturating_sub(matches),
            gapopen: 0,
            q_start,
            q_end,
            s_start,
            s_end,
            e_value: h.e_value,
            bit_score,
        });
    }

    if diag_enabled {
        eprintln!(
            "[4] Final hits after E-value filter (threshold={}): {}",
            evalue_threshold,
            final_hits.len()
        );

        // Report top hit alignment length
        if !final_hits.is_empty() {
            let max_len = final_hits.iter().map(|h| h.length).max().unwrap_or(0);
            eprintln!("[5] Top alignment length: {} AA", max_len);
        }
        eprintln!("=== End Stage Counters ===");
    }
    
    // Sort by bit score descending
    final_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    write_output(&final_hits, args.out.as_ref())?;
    Ok(())
}
