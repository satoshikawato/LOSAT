//! TBLASTX - LINE-BY-LINE NCBI BLAST PORT
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:439-619
//!            ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:48-131

use anyhow::{Context, Result};
use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::mpsc::channel;

use crate::common::{write_output, Hit};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::{lookup_protein_params, search_space::SearchSpace};
use crate::utils::dust::{DustMasker, MaskedInterval};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegMasker;

use super::args::TblastxArgs;
use super::constants::{CUTOFF_E_TBLASTX, GAP_TRIGGER_BIT_SCORE, X_DROP_UNGAPPED};
use super::chaining::UngappedHit;
use super::diagnostics::{diagnostics_enabled, DiagnosticCounters, print_summary as print_diagnostics_summary};
use super::sum_stats_linking::apply_sum_stats_even_gap_linking;
use super::extension::{convert_coords, extend_hit_two_hit};
use crate::stats::karlin::{bit_score as calc_bit_score, raw_score_from_evalue_with_decay, raw_score_from_bit_score};
use crate::stats::sum_statistics::defaults::GAP_DECAY_RATE_UNGAPPED;
use super::lookup::{build_ncbi_lookup, pv_test, BlastAaLookupTable, AA_HITS_PER_CELL};
use super::translation::{generate_frames, QueryFrame};

// [C] typedef struct DiagStruct { Int4 last_hit; Uint1 flag; } DiagStruct;
#[derive(Clone, Copy)]
struct DiagStruct {
    last_hit: i32,
    flag: u8,
}
impl Default for DiagStruct {
    fn default() -> Self { Self { last_hit: 0, flag: 0 } }
}

// [C] BlastOffsetPair with qs_offsets.q_off and qs_offsets.s_off
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

/// s_BlastAaScanSubject - LINE BY LINE PORT
/// Reference: blast_aascan.c:48-131
#[inline(never)]
fn s_blast_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [OffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],  // [C] Int4 *s_range
) -> i32 {
    // [C] Int4 totalhits = 0;
    let mut totalhits: i32 = 0;
    
    // NCBI subject masking support (masksubj.inl). In NCBI this walks
    // `subject->seq_ranges[]`; in LOSAT each translated frame is one contiguous
    // range already baked into `s_range[1..=2]`.
    //
    // [C] while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) { ... }
    while s_range[1] <= s_range[2] {
        // [C] s_first = subject->sequence + s_range[1];
        // [C] s_last  = subject->sequence + s_range[2];
        let s_first = s_range[1] as usize;
        let s_last = s_range[2] as usize;

        // Fast fail: nothing to scan.
        if s_first > s_last || subject.len() < 5 {
            return totalhits;
        }

        let word_length = lookup.word_length as usize; // [C] word_length = lookup->word_length;
        let charsize = lookup.charsize as usize;       // [C] lookup->charsize
        let mask = lookup.mask as usize;               // [C] lookup->mask
        let alphabet_size = lookup.alphabet_size as usize;

        // [C] index = ComputeTableIndex(word_length - 1, lookup->charsize, s_first);
        // LOSAT note: because we exclude stop-codon (24) from the lookup alphabet,
        // we must treat residues >= alphabet_size as invalid and rebuild the rolling
        // index just like the old base-24 code path did.
        let c0 = subject[s_first] as usize;
        let c1 = subject[s_first + 1] as usize;
        let mut index: usize = if c0 < alphabet_size && c1 < alphabet_size {
            (c0 << charsize) | c1
        } else {
            0
        };
        let mut skip: i32 = if c0 >= alphabet_size || c1 >= alphabet_size {
            word_length as i32
        } else {
            1
        };

        // [C] for (s = s_first; s <= s_last; s++)
        let mut s = s_first;
        while s <= s_last {
            // [C] index = ComputeTableIndexIncremental(word_length, lookup->charsize, lookup->mask, s, index);
            // [C] return ((index << charsize) | word[wordsize - 1]) & mask;
            let new_char = subject[s + (word_length - 1)] as usize;
            if new_char >= alphabet_size {
                // reset rolling index after an invalid residue (stop codon)
                skip = word_length as i32;
                index = 0;
                s += 1;
                continue;
            }

            index = ((index << charsize) | new_char) & mask;
            if skip > 0 {
                skip -= 1;
                if skip > 0 {
                    s += 1;
                    continue;
                }
            }

            let kmer = index;

            // [C] if (PV_TEST(pv, index, PV_ARRAY_BTS)) {
            if pv_test(&lookup.pv, kmer) {
                // [C] numhits = bbc[index].num_used;
                let cell = unsafe { lookup.backbone.get_unchecked(kmer) };
                let numhits = cell.num_used;
                if numhits != 0 {
                    // [C] if (numhits <= (array_size - totalhits)) { ... } else { s_range[1]=...; return totalhits; }
                    if numhits <= array_size - totalhits {
                        // [C] if (numhits <= AA_HITS_PER_CELL) src = bbc[index].payload.entries;
                        // [C] else src = &(ovfl[bbc[index].payload.overflow_cursor]);
                        let src: &[i32] = if numhits as usize <= AA_HITS_PER_CELL {
                            &cell.entries[..numhits as usize]
                        } else {
                            let cursor = cell.entries[0] as usize;
                            &lookup.overflow[cursor..cursor + numhits as usize]
                        };

                        // [C] Int4 s_off = (Int4)(s - subject->sequence);
                        let s_off = s as i32;

                        // [C] for (i = 0; i < numhits; i++) { offset_pairs[i+totalhits].q_off=src[i]; offset_pairs[i+totalhits].s_off=s_off; }
                        for i in 0..numhits as usize {
                            offset_pairs[(totalhits as usize) + i] = OffsetPair { q_off: src[i], s_off };
                        }
                        totalhits += numhits;
                    } else {
                        s_range[1] = s as i32;
                        return totalhits;
                    }
                }
            }

            s += 1;
        }

        // [C] s_range[1] = (Int4)(s - subject->sequence);
        s_range[1] = s as i32;
    }

    totalhits
}

pub fn run(args: TblastxArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 { num_cpus::get() } else { args.num_threads };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let only_qframe = args.only_qframe;
    let only_sframe = args.only_sframe;
    
    let valid_frame = |f: i8| matches!(f, 1 | 2 | 3 | -1 | -2 | -3);
    if let Some(f) = only_qframe { if !valid_frame(f) { anyhow::bail!("Invalid --only-qframe"); } }
    if let Some(f) = only_sframe { if !valid_frame(f) { anyhow::bail!("Invalid --only-sframe"); } }
    
    // [C] window = diag->window;
    let window = args.window_size as i32;
    // [C] wordsize = lookup->word_length;
    let wordsize: i32 = 3;
    let dropoff = X_DROP_UNGAPPED;
    let diag_enabled = diagnostics_enabled();
    let diagnostics = std::sync::Arc::new(DiagnosticCounters::default());

    rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global()
        .context("Failed to build thread pool")?;

    eprintln!("Reading queries...");
    let query_reader = fasta::Reader::from_file(&args.query)?;
    let queries_raw: Vec<fasta::Record> = query_reader.records().filter_map(|r| r.ok()).collect();
    let query_ids: Vec<String> = queries_raw.iter()
        .map(|r| r.id().split_whitespace().next().unwrap_or("unknown").to_string())
        .collect();

    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        queries_raw.iter().map(|r| masker.mask_sequence(r.seq())).collect()
    } else {
        queries_raw.iter().map(|_| Vec::new()).collect()
    };

    let mut query_frames: Vec<Vec<QueryFrame>> = queries_raw.iter()
        .map(|r| {
            let mut frames = generate_frames(r.seq(), &query_code);
            if let Some(f) = only_qframe { frames.retain(|x| x.frame == f); }
            frames
        })
        .collect();

    if args.seg {
        let seg = SegMasker::new(args.seg_window, args.seg_locut, args.seg_hicut);
        for frames in &mut query_frames {
            for frame in frames {
                if frame.aa_seq.len() >= 3 {
                    for m in seg.mask_sequence(&frame.aa_seq[1..frame.aa_seq.len()-1]) {
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
                    // LOSAT note: we encode amino acids in NCBI packed matrix order
                    // (ARNDCQEGHILKMFPSTWYVBJZX*) where X == 23, and we include
                    // sentinel bytes at aa_seq[0] and aa_seq[len-1].
                    if !frame.seg_masks.is_empty() {
                        if frame.aa_seq_nomask.is_none() {
                            frame.aa_seq_nomask = Some(frame.aa_seq.clone());
                        }
                        const X_MASK_NCBI_MATRIX: u8 = 23;
                        let raw_end_exclusive = frame.aa_seq.len().saturating_sub(1); // keep last sentinel untouched
                        for &(s, e) in &frame.seg_masks {
                            let raw_s = 1usize.saturating_add(s);
                            let raw_e = 1usize.saturating_add(e).min(raw_end_exclusive);
                            for pos in raw_s..raw_e {
                                frame.aa_seq[pos] = X_MASK_NCBI_MATRIX;
                            }
                        }
                    }
                }
            }
        }
    }

    eprintln!("Building lookup table...");
    let (lookup, contexts) = build_ncbi_lookup(&query_frames, &query_masks, args.threshold);

    eprintln!("Reading subjects...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> = subject_reader.records().filter_map(|r| r.ok()).collect();
    if queries_raw.is_empty() || subjects_raw.is_empty() { return Ok(()); }

    let scoring = ProteinScoringSpec { matrix: ScoringMatrix::Blosum62, gap_open: i32::MAX, gap_extend: i32::MAX };
    let params = lookup_protein_params(&scoring);

    eprintln!("Searching {} queries vs {} subjects...", queries_raw.len(), subjects_raw.len());

    let bar = ProgressBar::new(subjects_raw.len() as u64);
    bar.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len}")
        .unwrap());

    let (tx, rx) = channel::<Vec<Hit>>();
    let out_path = args.out.clone();
    let evalue_threshold = args.evalue;
    
    let writer = std::thread::spawn(move || -> Result<()> {
        let mut all: Vec<Hit> = Vec::new();
        while let Ok(h) = rx.recv() { all.extend(h); }
        all.retain(|h| h.e_value <= evalue_threshold);
        all.sort_by(|a, b| b.bit_score.partial_cmp(&a.bit_score).unwrap_or(std::cmp::Ordering::Equal));
        write_output(&all, out_path.as_ref())?;
        Ok(())
    });

    // [C] diag_mask = diag->diag_mask; (power of 2 - 1)
    let max_q: i32 = *lookup.frame_bases.last().unwrap_or(&0) + 
        contexts.last().map(|c| c.aa_seq.len() as i32).unwrap_or(0);
    let max_s: i32 = subjects_raw.iter().map(|r| r.seq().len() as i32 / 3 + 10).max().unwrap_or(0);
    let mut diag_array_size = 1i32;
    while diag_array_size < (max_q + max_s + window + 1) { diag_array_size <<= 1; }
    let diag_mask = diag_array_size - 1;

    // [C] array_size for offset_pairs
    const ARRAY_SIZE: i32 = 32768;

    let lookup_ref = &lookup;
    let contexts_ref = &contexts;
    let query_ids_ref = &query_ids;

    subjects_raw
        .par_iter()
        .enumerate()
        .for_each_init(
            || WorkerState {
                tx: tx.clone(),
                offset_pairs: vec![OffsetPair::default(); ARRAY_SIZE as usize],
                diag_array: vec![DiagStruct::default(); diag_array_size as usize],
                diag_offset: max_q + max_s,
            },
            |st, (s_idx, s_rec)| {
            let mut s_frames = generate_frames(s_rec.seq(), &db_code);
            if let Some(f) = only_sframe { s_frames.retain(|x| x.frame == f); }
            
            let s_id = s_rec.id().split_whitespace().next().unwrap_or("unknown").to_string();
            let s_len = s_rec.seq().len();
            
            // [C] BlastOffsetPair *offset_pairs
            let offset_pairs = &mut st.offset_pairs;
            
            // [C] DiagStruct *diag_array = diag->hit_level_array;
            let diag_array = &mut st.diag_array;
            
            // [C] diag_offset = diag->offset;  (per-thread, reused across subjects)
            let mut diag_offset: i32 = st.diag_offset;
            
            let mut cutoff_cache: FxHashMap<(usize, u8), i32> = FxHashMap::default();
            let mut ungapped_hits: Vec<UngappedHit> = Vec::new();
            
            for (s_f_idx, s_frame) in s_frames.iter().enumerate() {
                let subject = &s_frame.aa_seq;
                let s_aa_len = s_frame.aa_len;
                if subject.len() < 5 { continue; }
                
                // [C] scan_range[1] = subject->seq_ranges[0].left;
                // [C] scan_range[2] = subject->seq_ranges[0].right - wordsize;
                let mut scan_range: [i32; 3] = [0, 1, (subject.len() - 4) as i32];
                
                // [C] while (scan_range[1] <= scan_range[2])
                while scan_range[1] <= scan_range[2] {
                    // [C] hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
                    let hits = s_blast_aa_scan_subject(
                        lookup_ref,
                        subject,
                        offset_pairs,
                        ARRAY_SIZE,
                        &mut scan_range,
                    );
                    
                    if hits == 0 { break; }
                    
                    // [C] for (i = 0; i < hits; ++i)
                    for i in 0..hits as usize {
                        // [C] Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
                        // [C] Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
                        let query_offset = offset_pairs[i].q_off;
                        let subject_offset = offset_pairs[i].s_off;
                        
                        // [C] diag_coord = (query_offset - subject_offset) & diag_mask;
                        let diag_coord = ((query_offset - subject_offset) & diag_mask) as usize;
                        
                        // [C] if (diag_array[diag_coord].flag)
                        if diag_array[diag_coord].flag != 0 {
                            // [C] if ((Int4)(subject_offset + diag_offset) < diag_array[diag_coord].last_hit)
                            if subject_offset + diag_offset < diag_array[diag_coord].last_hit {
                                continue;
                            }
                            // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            // [C] diag_array[diag_coord].flag = 0;
                            diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            diag_array[diag_coord].flag = 0;
                        }
                        // [C] else
                        else {
                            // [C] last_hit = diag_array[diag_coord].last_hit - diag_offset;
                            let last_hit = diag_array[diag_coord].last_hit - diag_offset;
                            // [C] diff = subject_offset - last_hit;
                            let diff = subject_offset - last_hit;
                            
                            // [C] if (diff >= window)
                            if diff >= window {
                                diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // [C] if (diff < wordsize)
                            if diff < wordsize {
                                continue;
                            }
                            
                            // [C] curr_context = BSearchContextInfo(query_offset, query_info);
                            let ctx_idx = lookup_ref.get_context_idx(query_offset);
                            let ctx = &contexts_ref[ctx_idx];
                            let q_raw = (query_offset - ctx.frame_base) as usize;
                            let query = &ctx.aa_seq;
                            
                            // [C] if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                            if query_offset - diff < ctx.frame_base {
                                diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                continue;
                            }
                            
                            // [C] cutoffs = word_params->cutoffs + curr_context;
                            let key = (ctx_idx, s_f_idx as u8);
                            let cutoff = *cutoff_cache.entry(key).or_insert_with(|| {
                                let ss = SearchSpace {
                                    effective_query_len: ctx.aa_len.min(s_aa_len) as f64,
                                    effective_db_len: s_aa_len as f64,
                                    effective_space: (ctx.aa_len.min(s_aa_len) * s_aa_len) as f64,
                                    length_adjustment: 0,
                                };
                                let c = raw_score_from_evalue_with_decay(CUTOFF_E_TBLASTX, &params, &ss, true, GAP_DECAY_RATE_UNGAPPED);
                                let g = raw_score_from_bit_score(GAP_TRIGGER_BIT_SCORE, &params);
                                c.min(g as i32)
                            });
                            
                            // [C] score = s_BlastAaExtendTwoHit(matrix, subject, query,
                            //                                   last_hit + wordsize, subject_offset, query_offset, ...)
                            // Two-hit ungapped extension (NCBI `s_BlastAaExtendTwoHit`)
                            // Reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:1089-1158
                            let (hsp_q_u, hsp_qe_u, hsp_s_u, _hsp_se_u, score, right_extend, s_last_off_u) =
                                extend_hit_two_hit(
                                    query,
                                    subject,
                                    (last_hit + wordsize) as usize, // s_left_off (end of first hit word)
                                    subject_offset as usize,        // s_right_off (second hit start)
                                    q_raw as usize,                 // q_right_off (second hit start, local)
                                    dropoff,                        // x_dropoff (raw)
                                );

                            let hsp_q: i32 = hsp_q_u as i32;
                            let hsp_s: i32 = hsp_s_u as i32;
                            let hsp_len: i32 = (hsp_qe_u - hsp_q_u) as i32;
                            let s_last_off: i32 = s_last_off_u as i32;
                            
                            // [C] if (right_extend)
                            if right_extend {
                                // [C] diag_array[diag_coord].flag = 1;
                                // [C] diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                                diag_array[diag_coord].flag = 1;
                                diag_array[diag_coord].last_hit = s_last_off - (wordsize - 1) + diag_offset;
                            } else {
                                // [C] diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                                diag_array[diag_coord].last_hit = subject_offset + diag_offset;
                            }
                            
                            // [C] if (score >= cutoffs->cutoff_score)
                            if score >= cutoff {
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
                            }
                        }
                    }
                }
                
                // [C] Blast_ExtendWordExit(ewp, subject->length);
                diag_offset += s_aa_len as i32;
            }
            
            if !ungapped_hits.is_empty() {
                let linked = apply_sum_stats_even_gap_linking(ungapped_hits, &params);
                let mut final_hits: Vec<Hit> = Vec::new();
                for h in linked {
                    if h.e_value > evalue_threshold {
                        continue;
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
                    let (s_start, s_end) = convert_coords(h.s_aa_start, h.s_aa_end, s_frame.frame, s_len);

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
                    st.tx.send(final_hits).unwrap();
                }
            }
            bar.inc(1);
            st.diag_offset = diag_offset;
        });

    // Close channel so the writer can exit.
    drop(tx);

    bar.finish();
    writer.join().unwrap()?;
    if diag_enabled { print_diagnostics_summary(&diagnostics); }
    Ok(())
}
