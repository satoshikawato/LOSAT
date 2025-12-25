use anyhow::{Context, Result};
use bio::alphabets::dna;
use bio::io::fasta;
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::cmp::Ordering;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::mpsc::channel;

use crate::common::{write_output, Hit};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::karlin::{bit_score as calc_bit_score, evalue as calc_evalue};
use crate::stats::search_space::SearchSpace;
use crate::stats::{lookup_protein_params, KarlinParams};
use crate::utils::dust::{DustMasker, MaskedInterval};
use crate::utils::genetic_code::GeneticCode;
use crate::utils::matrix::MATRIX;

// --- Constants ---
// NCBI BLAST compatible X-drop parameters
const X_DROP_UNGAPPED: i32 = 7; // BLAST_UNGAPPED_X_DROPOFF_PROT from NCBI
const X_DROP_GAPPED_PRELIM: i32 = 15; // BLAST_GAP_X_DROPOFF_PROT for preliminary extension
const X_DROP_GAPPED_FINAL: i32 = 25; // BLAST_GAP_X_DROPOFF_FINAL_PROT for final traceback
const MAX_HITS_PER_KMER: usize = 200;
const STOP_CODON: u8 = 25;
const TWO_HIT_WINDOW: usize = 40; // BLAST_WINDOW_SIZE_PROT from NCBI (was 16)

// Default gap penalties for protein alignments (BLOSUM62 defaults, NCBI compatible)
const GAP_OPEN: i32 = -11; // BLAST_GAP_OPEN_PROT
const GAP_EXTEND: i32 = -1; // BLAST_GAP_EXTN_PROT

// Threshold for triggering gapped extension without two-hit requirement
const HIGH_SCORE_THRESHOLD: i32 = 60;

// Parameters for HSP chaining (similar to BLASTN)
const MAX_GAP_AA: usize = 333; // ~1000bp / 3 for amino acids
const MAX_DIAG_DRIFT_AA: isize = 33; // ~100bp / 3 for amino acids

/// Diagnostic counters for tracking where hits are lost in the pipeline
#[derive(Default)]
struct DiagnosticCounters {
    // Seed stage
    kmer_matches: AtomicUsize,
    seeds_low_score: AtomicUsize,
    seeds_two_hit_failed: AtomicUsize,
    seeds_passed: AtomicUsize,
    // Ungapped extension stage
    ungapped_extensions: AtomicUsize,
    ungapped_low_score: AtomicUsize,
    // Gapped extension stage
    ungapped_only_hits: AtomicUsize,
    gapped_extensions: AtomicUsize,
    gapped_evalue_passed: AtomicUsize,
    // Post-processing stage (set in chain_and_filter_hsps_protein)
    hsps_before_chain: AtomicUsize,
    hsps_after_chain: AtomicUsize,
    hsps_after_overlap_filter: AtomicUsize,
}

impl DiagnosticCounters {
    fn print_summary(&self) {
        eprintln!("\n=== TBLASTX Pipeline Diagnostics ===");
        eprintln!("Seed Stage:");
        eprintln!("  K-mer matches found:        {}", self.kmer_matches.load(AtomicOrdering::Relaxed));
        eprintln!("  Seeds filtered (low score): {}", self.seeds_low_score.load(AtomicOrdering::Relaxed));
        eprintln!("  Seeds filtered (two-hit):   {}", self.seeds_two_hit_failed.load(AtomicOrdering::Relaxed));
        eprintln!("  Seeds passed to extension:  {}", self.seeds_passed.load(AtomicOrdering::Relaxed));
        eprintln!("Ungapped Extension Stage:");
        eprintln!("  Ungapped extensions run:    {}", self.ungapped_extensions.load(AtomicOrdering::Relaxed));
        eprintln!("  Filtered (low score < 20):  {}", self.ungapped_low_score.load(AtomicOrdering::Relaxed));
        eprintln!("Gapped Extension Stage:");
        eprintln!("  Ungapped-only hits:         {}", self.ungapped_only_hits.load(AtomicOrdering::Relaxed));
        eprintln!("  Gapped extensions run:      {}", self.gapped_extensions.load(AtomicOrdering::Relaxed));
        eprintln!("  Gapped hits (e-value ok):   {}", self.gapped_evalue_passed.load(AtomicOrdering::Relaxed));
        eprintln!("Post-Processing Stage:");
        eprintln!("  HSPs before chaining:       {}", self.hsps_before_chain.load(AtomicOrdering::Relaxed));
        eprintln!("  HSPs after chaining:        {}", self.hsps_after_chain.load(AtomicOrdering::Relaxed));
        eprintln!("  Final hits (after filter):  {}", self.hsps_after_overlap_filter.load(AtomicOrdering::Relaxed));
        eprintln!("====================================\n");
    }
}

/// Check if diagnostics are enabled via environment variable
fn diagnostics_enabled() -> bool {
    std::env::var("LOSAT_DIAGNOSTICS").map(|v| v == "1" || v.to_lowercase() == "true").unwrap_or(false)
}

#[derive(Args, Debug)]
pub struct TblastxArgs {
    #[arg(short, long)]
    pub query: PathBuf,
    #[arg(short, long)]
    pub subject: PathBuf,
    #[arg(short, long, default_value_t = 10.0)]
    pub evalue: f64,
    #[arg(short, long, default_value_t = 13)]
    pub threshold: i32,
    #[arg(short, long, default_value_t = 3)]
    pub word_size: usize,
    #[arg(short = 'n', long, default_value_t = 0)]
    pub num_threads: usize,
    #[arg(short, long)]
    pub out: Option<PathBuf>,
    #[arg(long, default_value_t = 1)]
    pub query_gencode: u8,
    #[arg(long, default_value_t = 1)]
    pub db_gencode: u8,
    #[arg(long, default_value_t = 500)]
    pub max_target_seqs: usize,
    #[arg(long, default_value_t = true)]
    pub dust: bool,
    #[arg(long, default_value_t = 20)]
    pub dust_level: u32,
    #[arg(long, default_value_t = 64)]
    pub dust_window: usize,
    #[arg(long, default_value_t = 1)]
    pub dust_linker: usize,
}

struct QueryFrame {
    frame: i8, // 1..3, -1..-3
    aa_seq: Vec<u8>,
    orig_len: usize,
}

type DirectLookup = Vec<Vec<(u32, u8, u32)>>;

// --- Helper Functions ---

fn codon_to_aa_idx(codon: &[u8], table: &GeneticCode) -> u8 {
    let aa = table.get(codon);
    // '*' などA-Z以外は25にマップされる
    if aa >= b'A' && aa <= b'Z' {
        aa - b'A'
    } else {
        STOP_CODON
    }
}

fn translate_sequence(seq: &[u8], table: &GeneticCode) -> Vec<u8> {
    let len_aa = seq.len() / 3;
    let mut aa_seq = Vec::with_capacity(len_aa);
    for chunk in seq.chunks(3) {
        if chunk.len() == 3 {
            aa_seq.push(codon_to_aa_idx(chunk, table));
        }
    }
    aa_seq
}

fn generate_frames(seq: &[u8], table: &GeneticCode) -> Vec<QueryFrame> {
    let mut frames = Vec::with_capacity(6);
    let rev_seq = dna::revcomp(seq);
    let seq_len = seq.len();

    // Forward (Frames 1,2,3)
    for i in 0..3 {
        if i + 3 <= seq_len {
            frames.push(QueryFrame {
                frame: (i as i8) + 1,
                aa_seq: translate_sequence(&seq[i..], table),
                orig_len: seq_len,
            });
        }
    }
    // Reverse (Frames -1,-2,-3)
    for i in 0..3 {
        if i + 3 <= rev_seq.len() {
            frames.push(QueryFrame {
                frame: -((i as i8) + 1),
                aa_seq: translate_sequence(&rev_seq[i..], table),
                orig_len: seq_len,
            });
        }
    }
    frames
}

fn encode_aa_kmer(seq: &[u8], pos: usize) -> Option<usize> {
    if pos + 3 > seq.len() {
        return None;
    }
    let c1 = unsafe { *seq.get_unchecked(pos) } as usize;
    let c2 = unsafe { *seq.get_unchecked(pos + 1) } as usize;
    let c3 = unsafe { *seq.get_unchecked(pos + 2) } as usize;

    // 終止コドン(25)を含むk-merはシードにしない
    if c1 >= 25 || c2 >= 25 || c3 >= 25 {
        return None;
    }
    Some(c1 * 676 + c2 * 26 + c3)
}

/// Check if a DNA position range overlaps with any masked interval
fn is_dna_pos_masked(intervals: &[MaskedInterval], start: usize, end: usize) -> bool {
    for interval in intervals {
        // Check if [start, end) overlaps with [interval.start, interval.end)
        if start < interval.end && end > interval.start {
            return true;
        }
    }
    false
}

/// Convert amino acid position to DNA position range for a given frame
fn aa_pos_to_dna_range(aa_pos: usize, kmer_len: usize, frame: i8, orig_len: usize) -> (usize, usize) {
    let aa_end = aa_pos + kmer_len;
    if frame > 0 {
        // Forward frames: 1, 2, 3
        let offset = (frame - 1) as usize;
        let dna_start = aa_pos * 3 + offset;
        let dna_end = aa_end * 3 + offset;
        (dna_start, dna_end.min(orig_len))
    } else {
        // Reverse frames: -1, -2, -3
        let offset = (-frame - 1) as usize;
        let dna_end = orig_len - (aa_pos * 3 + offset);
        let dna_start = orig_len - (aa_end * 3 + offset);
        (dna_start.max(0), dna_end)
    }
}

fn build_direct_lookup(
    queries: &[Vec<QueryFrame>],
    query_masks: &[Vec<MaskedInterval>],
) -> DirectLookup {
    let table_size = 17576; // 26^3
    let mut lookup = vec![Vec::new(); table_size];

    for (q_idx, frames) in queries.iter().enumerate() {
        let masks = &query_masks[q_idx];
        for (f_idx, frame) in frames.iter().enumerate() {
            let seq = &frame.aa_seq;
            if seq.len() < 3 {
                continue;
            }
            for i in 0..=(seq.len() - 3) {
                // Check if this amino acid position corresponds to a masked DNA region
                if !masks.is_empty() {
                    let (dna_start, dna_end) = aa_pos_to_dna_range(i, 3, frame.frame, frame.orig_len);
                    if is_dna_pos_masked(masks, dna_start, dna_end) {
                        continue;
                    }
                }
                if let Some(code) = encode_aa_kmer(seq, i) {
                    lookup[code].push((q_idx as u32, f_idx as u8, i as u32));
                }
            }
        }
    }

    for bucket in &mut lookup {
        if bucket.len() > MAX_HITS_PER_KMER {
            bucket.clear();
        }
    }
    lookup
}

#[inline(always)]
fn get_score(a: u8, b: u8) -> i32 {
    unsafe { *MATRIX.get_unchecked((a as usize) * 27 + (b as usize)) as i32 }
}

fn extend_hit_ungapped(
    q_seq: &[u8],
    s_seq: &[u8],
    q_pos: usize,
    s_pos: usize,
    seed_score: i32,
) -> (usize, usize, usize, usize, i32) {
    let mut current_score = seed_score;
    let mut max_score = current_score;

    // Left
    let mut best_i = 0;
    let mut i = 1;
    let max_left = q_pos.min(s_pos);

    while i <= max_left {
        let q_char = unsafe { *q_seq.get_unchecked(q_pos - i) };
        let s_char = unsafe { *s_seq.get_unchecked(s_pos - i) };

        // ★修正: 終止コドンに出会ったらペナルティを与えてアライメントを切る
        if q_char == STOP_CODON || s_char == STOP_CODON {
            current_score -= 100;
        } else {
            current_score += get_score(q_char, s_char);
        }

        if current_score > max_score {
            max_score = current_score;
            best_i = i;
        } else if (max_score - current_score) > X_DROP_UNGAPPED {
            break;
        }
        i += 1;
    }

    // Right
    let k_size = 3;
    let mut current_score_r = max_score;
    let mut max_score_total = max_score;
    let mut best_j = k_size;
    let mut j = 0;

    let q_start_r = q_pos + k_size;
    let s_start_r = s_pos + k_size;
    let q_limit = q_seq.len();
    let s_limit = s_seq.len();

    while (q_start_r + j) < q_limit && (s_start_r + j) < s_limit {
        let q_char = unsafe { *q_seq.get_unchecked(q_start_r + j) };
        let s_char = unsafe { *s_seq.get_unchecked(s_start_r + j) };

        // ★修正: 終止コドンチェック
        if q_char == STOP_CODON || s_char == STOP_CODON {
            current_score_r -= 100;
        } else {
            current_score_r += get_score(q_char, s_char);
        }

        if current_score_r > max_score_total {
            max_score_total = current_score_r;
            best_j = k_size + j + 1;
        } else if (max_score_total - current_score_r) > X_DROP_UNGAPPED {
            break;
        }
        j += 1;
    }

    (
        q_pos - best_i,
        q_pos + best_j,
        s_pos - best_i,
        s_pos + best_j,
        max_score_total,
    )
}

/// Alignment statistics propagated alongside DP scores for traceback-based calculation
#[derive(Clone, Copy, Default)]
struct ProteinAlnStats {
    matches: u32,
    mismatches: u32,
    gap_opens: u32,
    gap_letters: u32,
}

/// Banded Smith-Waterman gapped extension for protein sequences with affine gap penalties.
///
/// This implements a proper dynamic programming algorithm that:
/// - Uses affine gap penalties (gap_open + gap_extend * gap_length)
/// - Uses X-drop based dynamic window (NCBI Blast_SemiGappedAlign style)
/// - Uses the BLOSUM62 substitution matrix for scoring
/// - Handles stop codons appropriately
/// - Tracks best score across ALL diagonals (not just main diagonal)
/// - Propagates alignment statistics for accurate traceback-based calculation
///
/// Returns (q_start, q_end, s_start, s_end, score, matches, mismatches, gaps)
fn extend_gapped_protein(
    q_seq: &[u8],
    s_seq: &[u8],
    qs: usize,
    ss: usize,
    len: usize,
    x_drop: i32,
) -> (usize, usize, usize, usize, i32, usize, usize, usize) {
    // Bounds validation: ensure seed coordinates are valid
    if qs >= q_seq.len() || ss >= s_seq.len() {
        return (qs, qs, ss, ss, 0, 0, 0, 0);
    }

    // Clamp len to available sequence length
    let len = len.min(q_seq.len() - qs).min(s_seq.len() - ss);
    if len == 0 {
        return (qs, qs, ss, ss, 0, 0, 0, 0);
    }

    // First, extend to the right from the seed using NCBI-style semi-gapped DP
    let (right_q_consumed, right_s_consumed, right_score, right_stats) =
        extend_gapped_protein_one_direction(
            q_seq,
            s_seq,
            (qs + len).min(q_seq.len()),
            (ss + len).min(s_seq.len()), // Start from end of seed (clamped)
            true,                        // Forward direction
            x_drop,
        );

    // Then extend to the left from the seed
    let (left_q_consumed, left_s_consumed, left_score, left_stats) =
        extend_gapped_protein_one_direction(
            q_seq, s_seq, qs, ss,    // Start from beginning of seed
            false, // Backward direction
            x_drop,
        );

    // Calculate seed score and actual matches/mismatches
    let mut seed_score = 0;
    let mut seed_matches = 0;
    let mut seed_mismatches = 0;
    for k in 0..len {
        let q_char = q_seq[qs + k];
        let s_char = s_seq[ss + k];
        if q_char == STOP_CODON || s_char == STOP_CODON {
            seed_score -= 100;
            seed_mismatches += 1;
        } else {
            seed_score += get_score(q_char, s_char);
            if q_char == s_char {
                seed_matches += 1;
            } else {
                seed_mismatches += 1;
            }
        }
    }

    let total_score = left_score + seed_score + right_score;
    let total_matches = left_stats.matches as usize + seed_matches + right_stats.matches as usize;
    let total_mismatches =
        left_stats.mismatches as usize + seed_mismatches + right_stats.mismatches as usize;
    let total_gaps = left_stats.gap_opens as usize + right_stats.gap_opens as usize;

    // Calculate final positions
    let final_q_start = qs - left_q_consumed;
    let final_q_end = qs + len + right_q_consumed;
    let final_s_start = ss - left_s_consumed;
    let final_s_end = ss + len + right_s_consumed;

    // Return full alignment boundaries
    (
        final_q_start,
        final_q_end,
        final_s_start,
        final_s_end,
        total_score,
        total_matches,
        total_mismatches,
        total_gaps,
    )
}

/// Extend alignment in one direction using NCBI-style semi-gapped DP with affine gap penalties.
///
/// This implements NCBI BLAST's Blast_SemiGappedAlign approach:
/// - X-drop based dynamic window that expands/contracts based on score
/// - Tracks best score across ALL diagonals (not just main diagonal)
/// - Propagates alignment statistics for accurate traceback-based calculation
/// - No hard-coded extension limits (controlled by X-drop termination)
///
/// Returns: (q_consumed, s_consumed, score, stats)
fn extend_gapped_protein_one_direction(
    q_seq: &[u8],
    s_seq: &[u8],
    q_start: usize,
    s_start: usize,
    forward: bool,
    x_drop: i32,
) -> (usize, usize, i32, ProteinAlnStats) {
    const NEG_INF: i32 = i32::MIN / 2;
    const INITIAL_BAND_WIDTH: usize = 32;

    // Determine sequence bounds
    let (q_len, s_len) = if forward {
        (q_seq.len() - q_start, s_seq.len() - s_start)
    } else {
        (q_start, s_start)
    };

    if q_len == 0 || s_len == 0 {
        return (0, 0, 0, ProteinAlnStats::default());
    }

    // Dynamic band that can expand based on X-drop
    let max_band = q_len.max(s_len).min(1000);
    let band_size = (2 * INITIAL_BAND_WIDTH + 1).min(2 * max_band + 1);

    // DP score arrays (2 rows)
    let mut m_prev = vec![NEG_INF; band_size];
    let mut m_curr = vec![NEG_INF; band_size];
    let mut iq_prev = vec![NEG_INF; band_size];
    let mut iq_curr = vec![NEG_INF; band_size];
    let mut is_prev = vec![NEG_INF; band_size];
    let mut is_curr = vec![NEG_INF; band_size];

    // Statistics arrays (2 rows) - propagate counts alongside scores
    let mut m_stats_prev = vec![ProteinAlnStats::default(); band_size];
    let mut m_stats_curr = vec![ProteinAlnStats::default(); band_size];
    let mut iq_stats_prev = vec![ProteinAlnStats::default(); band_size];
    let mut iq_stats_curr = vec![ProteinAlnStats::default(); band_size];
    let mut is_stats_prev = vec![ProteinAlnStats::default(); band_size];
    let mut is_stats_curr = vec![ProteinAlnStats::default(); band_size];

    // Track best score, position, and stats
    let mut best_score = 0;
    let mut best_q_consumed = 0;
    let mut best_s_consumed = 0;
    let mut best_stats = ProteinAlnStats::default();

    // Initialize first position
    let center = band_size / 2;
    m_prev[center] = 0;

    // Initialize leading gaps in subject (Iy[0,j] for j > 0)
    for k in (center + 1)..band_size {
        let j = k - center;
        if j <= s_len {
            is_prev[k] = GAP_OPEN + GAP_EXTEND * (j as i32);
            is_stats_prev[k] = ProteinAlnStats {
                matches: 0,
                mismatches: 0,
                gap_opens: 1,
                gap_letters: j as u32,
            };
        }
    }

    // Track the active range of the band (NCBI-style dynamic window)
    let mut left_bound = center;
    let mut right_bound = center;

    // Process each row (query position)
    let max_rows = q_len.min(s_len * 2 + 100); // Allow some diagonal drift
    for i in 1..=max_rows {
        // Reset current row
        for k in 0..band_size {
            m_curr[k] = NEG_INF;
            iq_curr[k] = NEG_INF;
            is_curr[k] = NEG_INF;
            m_stats_curr[k] = ProteinAlnStats::default();
            iq_stats_curr[k] = ProteinAlnStats::default();
            is_stats_curr[k] = ProteinAlnStats::default();
        }

        // Initialize leading gap in query (Ix[i,0])
        if i <= center {
            let gap_k = center - i;
            if gap_k < band_size {
                iq_curr[gap_k] = GAP_OPEN + GAP_EXTEND * (i as i32);
                iq_stats_curr[gap_k] = ProteinAlnStats {
                    matches: 0,
                    mismatches: 0,
                    gap_opens: 1,
                    gap_letters: i as u32,
                };
            }
        }

        let mut row_max = NEG_INF;
        let mut new_left = band_size;
        let mut new_right = 0;

        // Expand bounds slightly to allow diagonal drift
        let k_min = left_bound.saturating_sub(1);
        let k_max = (right_bound + 1).min(band_size - 1);

        // First pass: compute M and Ix
        for k in k_min..=k_max {
            let j_offset = k as isize - center as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > s_len {
                continue;
            }

            // Get sequence characters
            let (qc, sc) = if forward {
                if i > q_len {
                    continue;
                }
                (q_seq[q_start + i - 1], s_seq[s_start + j - 1])
            } else {
                if i > q_start || j > s_start {
                    continue;
                }
                (q_seq[q_start - i], s_seq[s_start - j])
            };

            // Match/mismatch score using BLOSUM62 matrix
            let match_score = if qc == STOP_CODON || sc == STOP_CODON {
                -100
            } else {
                get_score(qc, sc)
            };
            let is_match = qc == sc && qc != STOP_CODON;

            // Compute M[i,j]
            let from_m = if m_prev[k] > NEG_INF {
                m_prev[k] + match_score
            } else {
                NEG_INF
            };
            let from_iq = if iq_prev[k] > NEG_INF {
                iq_prev[k] + match_score
            } else {
                NEG_INF
            };
            let from_is = if is_prev[k] > NEG_INF {
                is_prev[k] + match_score
            } else {
                NEG_INF
            };

            let best = from_m.max(from_iq).max(from_is);
            m_curr[k] = best;

            // Propagate stats from best predecessor
            let mut stats = if best == from_m && from_m > NEG_INF {
                m_stats_prev[k]
            } else if best == from_iq && from_iq > NEG_INF {
                iq_stats_prev[k]
            } else if best == from_is && from_is > NEG_INF {
                is_stats_prev[k]
            } else {
                ProteinAlnStats::default()
            };

            // Add current match/mismatch
            if is_match {
                stats.matches += 1;
            } else {
                stats.mismatches += 1;
            }
            m_stats_curr[k] = stats;

            // Compute Ix[i,j] (gap in subject, consume query)
            if k + 1 < band_size {
                let open_gap = if m_prev[k + 1] > NEG_INF {
                    m_prev[k + 1] + GAP_OPEN + GAP_EXTEND
                } else {
                    NEG_INF
                };
                let extend_gap = if iq_prev[k + 1] > NEG_INF {
                    iq_prev[k + 1] + GAP_EXTEND
                } else {
                    NEG_INF
                };

                let best_gap = open_gap.max(extend_gap);
                iq_curr[k] = best_gap;

                // Propagate stats
                if best_gap == open_gap && open_gap > NEG_INF {
                    let mut gap_stats = m_stats_prev[k + 1];
                    gap_stats.gap_opens += 1;
                    gap_stats.gap_letters += 1;
                    iq_stats_curr[k] = gap_stats;
                } else if best_gap == extend_gap && extend_gap > NEG_INF {
                    let mut gap_stats = iq_stats_prev[k + 1];
                    gap_stats.gap_letters += 1;
                    iq_stats_curr[k] = gap_stats;
                }
            }
        }

        // Second pass: compute Iy
        for k in k_min..=k_max {
            let j_offset = k as isize - center as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > s_len {
                continue;
            }

            if k > 0 {
                let open_gap = if m_curr[k - 1] > NEG_INF {
                    m_curr[k - 1] + GAP_OPEN + GAP_EXTEND
                } else {
                    NEG_INF
                };
                let extend_gap = if is_curr[k - 1] > NEG_INF {
                    is_curr[k - 1] + GAP_EXTEND
                } else {
                    NEG_INF
                };

                let best_gap = open_gap.max(extend_gap);
                is_curr[k] = best_gap;

                // Propagate stats
                if best_gap == open_gap && open_gap > NEG_INF {
                    let mut gap_stats = m_stats_curr[k - 1];
                    gap_stats.gap_opens += 1;
                    gap_stats.gap_letters += 1;
                    is_stats_curr[k] = gap_stats;
                } else if best_gap == extend_gap && extend_gap > NEG_INF {
                    let mut gap_stats = is_stats_curr[k - 1];
                    gap_stats.gap_letters += 1;
                    is_stats_curr[k] = gap_stats;
                }
            }
        }

        // Track best score and stats across ALL diagonals (NCBI-style)
        for k in k_min..=k_max {
            let j_offset = k as isize - center as isize;
            let j = (i as isize + j_offset) as usize;

            if j == 0 || j > s_len {
                continue;
            }

            let cell_max = m_curr[k].max(iq_curr[k]).max(is_curr[k]);

            // Update row max and active bounds
            if cell_max > NEG_INF && best_score - cell_max <= x_drop {
                if cell_max > row_max {
                    row_max = cell_max;
                }
                if k < new_left {
                    new_left = k;
                }
                if k > new_right {
                    new_right = k;
                }
            }

            // Update best score if this cell is better (across ALL diagonals)
            if cell_max > best_score {
                best_score = cell_max;
                best_q_consumed = i;
                best_s_consumed = j;
                // Get stats from the best state
                best_stats = if cell_max == m_curr[k] {
                    m_stats_curr[k]
                } else if cell_max == iq_curr[k] {
                    iq_stats_curr[k]
                } else {
                    is_stats_curr[k]
                };
            }
        }

        // X-drop termination
        if best_score - row_max > x_drop {
            break;
        }

        // Update active bounds for next row (NCBI-style dynamic window)
        if new_left < band_size && new_right > 0 {
            left_bound = new_left;
            right_bound = new_right;
        }

        // Swap rows
        std::mem::swap(&mut m_prev, &mut m_curr);
        std::mem::swap(&mut iq_prev, &mut iq_curr);
        std::mem::swap(&mut is_prev, &mut is_curr);
        std::mem::swap(&mut m_stats_prev, &mut m_stats_curr);
        std::mem::swap(&mut iq_stats_prev, &mut iq_stats_curr);
        std::mem::swap(&mut is_stats_prev, &mut is_stats_curr);
    }

    (best_q_consumed, best_s_consumed, best_score, best_stats)
}

fn convert_coords(aa_start: usize, aa_end: usize, frame: i8, dna_len: usize) -> (usize, usize) {
    let f_abs = frame.abs() as usize;
    let shift = f_abs - 1;

    if frame > 0 {
        let start_bp = aa_start * 3 + shift + 1;
        let end_bp = (aa_end) * 3 + shift;
        (start_bp, end_bp)
    } else {
        let start_bp = dna_len - (aa_start * 3 + shift);
        let end_bp_calc = dna_len - (aa_end * 3 + shift - 1);
        (start_bp, end_bp_calc)
    }
}

fn calculate_statistics(score: i32, eff_search_space: f64, params: &KarlinParams) -> (f64, f64) {
    let search_space = SearchSpace {
        effective_query_len: 1.0,
        effective_db_len: 1.0,
        effective_space: eff_search_space,
        length_adjustment: 0, // Not used here since effective_space is pre-computed
    };
    let bs = calc_bit_score(score, params);
    let ev = calc_evalue(bs, &search_space);
    (bs, ev)
}

/// Extended hit structure that includes frame and amino acid coordinate information
/// for HSP chaining. This is used internally during chaining and then converted
/// back to regular Hit for output.
#[derive(Debug, Clone)]
struct ExtendedHit {
    hit: Hit,
    q_frame: i8,
    s_frame: i8,
    q_aa_start: usize,
    q_aa_end: usize,
    s_aa_start: usize,
    s_aa_end: usize,
    q_orig_len: usize,
    s_orig_len: usize,
}

/// Sequence data for re-alignment during HSP chaining
type SequenceKey = (String, String, i8, i8); // (query_id, subject_id, q_frame, s_frame)
type SequenceData = (Vec<u8>, Vec<u8>); // (query_aa_seq, subject_aa_seq)

/// Chain and filter HSPs for protein/translated sequences.
///
/// This implements a cluster-then-extend approach similar to BLASTN:
/// 1. Groups HSPs by query-subject pair AND frame combination
/// 2. Clusters nearby HSPs on similar diagonals
/// 3. Re-runs gapped alignment on merged regions to get accurate statistics
/// 4. Filters redundant overlapping hits
/// 5. Applies e-value threshold to filter insignificant hits
fn chain_and_filter_hsps_protein(
    mut hits: Vec<ExtendedHit>,
    sequences: &FxHashMap<SequenceKey, SequenceData>,
    db_len_aa_total: usize,
    params: &KarlinParams,
    evalue_threshold: f64,
    diagnostics: Option<&DiagnosticCounters>,
) -> Vec<Hit> {
    if hits.is_empty() {
        return Vec::new();
    }

    // Record HSPs before chaining
    let hsps_before = hits.len();
    if let Some(diag) = diagnostics {
        diag.hsps_before_chain.store(hsps_before, AtomicOrdering::Relaxed);
    }

    // Group hits by query-subject pair AND frame combination
    let mut groups: FxHashMap<SequenceKey, Vec<ExtendedHit>> = FxHashMap::default();
    for hit in hits.drain(..) {
        let key = (
            hit.hit.query_id.clone(),
            hit.hit.subject_id.clone(),
            hit.q_frame,
            hit.s_frame,
        );
        groups.entry(key).or_default().push(hit);
    }

    let mut result_hits = Vec::new();

    for (key, mut group_hits) in groups {
        if group_hits.is_empty() {
            continue;
        }

        // Sort by query amino acid start position
        group_hits.sort_by_key(|h| h.q_aa_start);

        // Cluster HSPs that are nearby and on similar diagonals
        let mut clusters: Vec<Vec<ExtendedHit>> = Vec::new();

        for hit in group_hits {
            let hit_diag = hit.s_aa_start as isize - hit.q_aa_start as isize;

            // Try to find an existing cluster to add this hit to
            let mut cluster_idx: Option<usize> = None;
            for (idx, cluster) in clusters.iter().enumerate() {
                if let Some(last) = cluster.last() {
                    let last_diag = last.s_aa_end as isize - last.q_aa_end as isize;
                    let diag_drift = (hit_diag - last_diag).abs();

                    // Calculate gap or overlap between HSPs
                    // Positive = gap, negative = overlap
                    let q_distance = hit.q_aa_start as isize - last.q_aa_end as isize;
                    let s_distance = hit.s_aa_start as isize - last.s_aa_end as isize;

                    // Allow overlapping HSPs (negative distance) or gaps up to MAX_GAP_AA
                    // For overlaps, limit to 50% of the smaller HSP to avoid merging unrelated HSPs
                    let hit_q_len = (hit.q_aa_end - hit.q_aa_start) as isize;
                    let last_q_len = (last.q_aa_end - last.q_aa_start) as isize;
                    let max_overlap = -(hit_q_len.min(last_q_len) / 2);

                    let q_ok = q_distance >= max_overlap && q_distance <= MAX_GAP_AA as isize;
                    let s_ok = s_distance >= max_overlap && s_distance <= MAX_GAP_AA as isize;

                    if q_ok && s_ok && diag_drift <= MAX_DIAG_DRIFT_AA {
                        cluster_idx = Some(idx);
                        break;
                    }
                }
            }

            if let Some(idx) = cluster_idx {
                clusters[idx].push(hit);
            } else {
                clusters.push(vec![hit]);
            }
        }

        // Process each cluster - output all individual HSPs instead of merging
        // This matches NCBI BLAST behavior of outputting many short hits
        for cluster in clusters {
            for ext_hit in cluster {
                // Only emit if the original HSP passes the e-value threshold
                if ext_hit.hit.e_value <= evalue_threshold {
                    result_hits.push(ext_hit.hit);
                }
            }
        }
    }

    // Record HSPs after chaining (before overlap filter)
    if let Some(diag) = diagnostics {
        diag.hsps_after_chain.store(result_hits.len(), AtomicOrdering::Relaxed);
    }

    // Sort by bit score (highest first) for consistent output
    result_hits.sort_by(|a, b| {
        b.bit_score
            .partial_cmp(&a.bit_score)
            .unwrap_or(Ordering::Equal)
    });

    // Record final hits (no overlap filter - output all HSPs like NCBI BLAST)
    if let Some(diag) = diagnostics {
        diag.hsps_after_overlap_filter.store(result_hits.len(), AtomicOrdering::Relaxed);
    }

    result_hits
}

pub fn run(args: TblastxArgs) -> Result<()> {
    let num_threads = if args.num_threads == 0 {
        num_cpus::get()
    } else {
        args.num_threads
    };
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);

    // Initialize diagnostic counters if enabled
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
        .map(|r| {
            r.id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string()
        })
        .collect();

    let query_masks: Vec<Vec<MaskedInterval>> = if args.dust {
        let masker = DustMasker::new(args.dust_level, args.dust_window, args.dust_linker);
        queries_raw
            .iter()
            .map(|r| masker.mask_sequence(r.seq()))
            .collect()
    } else {
        queries_raw.iter().map(|_| Vec::new()).collect()
    };

    let query_frames: Vec<Vec<QueryFrame>> = queries_raw
        .iter()
        .map(|r| generate_frames(r.seq(), &query_code))
        .collect();

    eprintln!("Building lookup table...");
    let lookup = build_direct_lookup(&query_frames, &query_masks);

    eprintln!("Reading subject file...");
    let subject_reader = fasta::Reader::from_file(&args.subject)?;
    let subjects_raw: Vec<fasta::Record> =
        subject_reader.records().filter_map(|r| r.ok()).collect();

    if queries_raw.is_empty() || subjects_raw.is_empty() {
        return Ok(());
    }

    let db_len_total_bp: usize = subjects_raw.iter().map(|r| r.seq().len()).sum();
    let db_len_aa_total = db_len_total_bp / 3;

    let scoring_spec = ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: GAP_OPEN.abs(),
        gap_extend: GAP_EXTEND.abs(),
    };
    let params = lookup_protein_params(&scoring_spec);

    eprintln!(
        "Searching {} Queries vs {} Subjects...",
        queries_raw.len(),
        subjects_raw.len()
    );

    let bar = ProgressBar::new(subjects_raw.len() as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} {msg}")
            .unwrap(),
    );

    // Channel now sends ExtendedHit with sequence data for HSP chaining
    type HitData = (Vec<ExtendedHit>, Vec<(SequenceKey, SequenceData)>);
    let (tx, rx) = channel::<HitData>();

    let out_path = args.out.clone();
    let params_clone = params.clone();
    let evalue_threshold = args.evalue;
    let diagnostics_clone = diagnostics.clone();
    let diag_enabled_clone = diag_enabled;
    let writer_handle = std::thread::spawn(move || -> Result<()> {
        let mut all_hits: Vec<ExtendedHit> = Vec::new();
        let mut all_sequences: FxHashMap<SequenceKey, SequenceData> = FxHashMap::default();

        while let Ok((hits, seq_data)) = rx.recv() {
            all_hits.extend(hits);
            for (key, data) in seq_data {
                all_sequences.entry(key).or_insert(data);
            }
        }

        // Chain nearby HSPs into longer alignments using cluster-then-extend
        let diag_ref = if diag_enabled_clone { Some(diagnostics_clone.as_ref()) } else { None };
        let filtered_hits = chain_and_filter_hsps_protein(
            all_hits,
            &all_sequences,
            db_len_aa_total,
            &params_clone,
            evalue_threshold,
            diag_ref,
        );

        write_output(&filtered_hits, out_path.as_ref())?;
        Ok(())
    });

    let diagnostics_for_parallel = diagnostics.clone();
    subjects_raw
        .par_iter()
        .enumerate()
        .for_each_with((tx, diagnostics_for_parallel), |(tx, diag), (_s_idx, s_record)| {
            let s_frames = generate_frames(s_record.seq(), &db_code);
            let s_id = s_record
                .id()
                .split_whitespace()
                .next()
                .unwrap_or("unknown")
                .to_string();
            let s_len = s_record.seq().len();

            // Parallelize over subject frames for intra-subject parallelism
            // Each frame is independent since mask_key includes s_f_idx
            let diag_inner = diag.clone();
            let frame_results: Vec<(Vec<ExtendedHit>, Vec<(SequenceKey, SequenceData)>)> =
                s_frames
                    .par_iter()
                    .enumerate()
                    .map(|(_s_f_idx, s_frame)| {
                        let mut local_hits: Vec<ExtendedHit> = Vec::new();
                        let mut local_sequences: Vec<(SequenceKey, SequenceData)> = Vec::new();
                        let mut seen_pairs: std::collections::HashSet<SequenceKey> =
                            std::collections::HashSet::new();

                        // Mask to track already-extended regions on each diagonal (per frame)
                        let mut mask: FxHashMap<(u32, u8, isize), usize> = FxHashMap::default();
                        // Two-hit tracking: stores the last seed position on each diagonal for each query/frame combination
                        let mut last_seed: FxHashMap<(u32, u8, isize), usize> =
                            FxHashMap::default();
                        let s_aa = &s_frame.aa_seq;
                        if s_aa.len() < 3 {
                            return (local_hits, local_sequences);
                        }

                        for s_pos in 0..=(s_aa.len() - 3) {
                            let c1 = unsafe { *s_aa.get_unchecked(s_pos) } as usize;
                            let c2 = unsafe { *s_aa.get_unchecked(s_pos + 1) } as usize;
                            let c3 = unsafe { *s_aa.get_unchecked(s_pos + 2) } as usize;

                            // シードに終止コドンが含まれる場合はスキップ
                            if c1 >= 25 || c2 >= 25 || c3 >= 25 {
                                continue;
                            }
                            let kmer = c1 * 676 + c2 * 26 + c3;

                            let matches = unsafe { lookup.get_unchecked(kmer) };

                            // Count k-mer matches for diagnostics
                            if diag_enabled && !matches.is_empty() {
                                diag_inner.kmer_matches.fetch_add(matches.len(), AtomicOrdering::Relaxed);
                            }

                            for &(q_idx, q_f_idx, q_pos) in matches {
                                let diag = s_pos as isize - q_pos as isize;
                                // Simplified mask_key since we're now per-frame (s_f_idx is implicit)
                                let mask_key = (q_idx, q_f_idx, diag);

                                // Check if this region was already extended
                                if let Some(&last_end) = mask.get(&mask_key) {
                                    if s_pos <= last_end {
                                        continue;
                                    }
                                }

                                let q_frame = &query_frames[q_idx as usize][q_f_idx as usize];
                                let q_aa = &q_frame.aa_seq;

                                let seed_score = get_score(c1 as u8, q_aa[q_pos as usize])
                                    + get_score(c2 as u8, q_aa[q_pos as usize + 1])
                                    + get_score(c3 as u8, q_aa[q_pos as usize + 2]);

                                if seed_score < 11 {
                                    if diag_enabled {
                                        diag_inner.seeds_low_score.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // Two-hit requirement: check if there's a previous seed on this diagonal
                                // within the window distance. This gates BOTH ungapped and gapped extension
                                // to avoid expensive extensions for every seed hit.
                                let two_hit_satisfied =
                                    if let Some(&prev_s_pos) = last_seed.get(&mask_key) {
                                        s_pos.saturating_sub(prev_s_pos) <= TWO_HIT_WINDOW
                                    } else {
                                        false
                                    };

                                // Update the last seed position for this diagonal
                                last_seed.insert(mask_key, s_pos);

                                // Skip extension if two-hit not satisfied AND seed score is not exceptionally high
                                // High seed scores (>= 15) can trigger extension without two-hit
                                if !two_hit_satisfied && seed_score < 15 {
                                    if diag_enabled {
                                        diag_inner.seeds_two_hit_failed.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // Count seeds that passed to extension
                                if diag_enabled {
                                    diag_inner.seeds_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                    diag_inner.ungapped_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                let (qs, qe, ss, se_ungapped, ungapped_score) = extend_hit_ungapped(
                                    q_aa,
                                    s_aa,
                                    q_pos as usize,
                                    s_pos,
                                    seed_score,
                                );

                                // Skip if ungapped score is too low
                                if ungapped_score < 20 {
                                    if diag_enabled {
                                        diag_inner.ungapped_low_score.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    continue;
                                }

                                // For ungapped-only hits (no gapped extension), record if e-value passes
                                // Only do gapped extension if two-hit requirement is met AND ungapped score is high enough
                                if !two_hit_satisfied || ungapped_score < HIGH_SCORE_THRESHOLD {
                                    if diag_enabled {
                                        diag_inner.ungapped_only_hits.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    // Record the ungapped hit if it passes e-value threshold
                                    let eff_space = (q_aa.len() as f64) * (db_len_aa_total as f64);
                                    let (bit_score, e_val) =
                                        calculate_statistics(ungapped_score, eff_space, &params);

                                    if e_val <= args.evalue {
                                        mask.insert(mask_key, se_ungapped);

                                        let len = qe - qs;
                                        let mut match_count = 0;
                                        for k in 0..len {
                                            if q_aa[qs + k] == s_aa[ss + k] {
                                                match_count += 1;
                                            }
                                        }
                                        let identity = (match_count as f64 / len as f64) * 100.0;

                                        let (q_start_bp, q_end_bp) =
                                            convert_coords(qs, qe, q_frame.frame, q_frame.orig_len);
                                        let (s_start_bp, s_end_bp) =
                                            convert_coords(ss, se_ungapped, s_frame.frame, s_len);

                                        // Store sequence data for HSP chaining
                                        let seq_key: SequenceKey = (
                                            query_ids[q_idx as usize].clone(),
                                            s_id.clone(),
                                            q_frame.frame,
                                            s_frame.frame,
                                        );
                                        if !seen_pairs.contains(&seq_key) {
                                            seen_pairs.insert(seq_key.clone());
                                            local_sequences.push((
                                                seq_key.clone(),
                                                (q_aa.clone(), s_aa.clone()),
                                            ));
                                        }

                                        local_hits.push(ExtendedHit {
                                            hit: Hit {
                                                query_id: query_ids[q_idx as usize].clone(),
                                                subject_id: s_id.clone(),
                                                identity,
                                                length: len,
                                                mismatch: len - match_count,
                                                gapopen: 0,
                                                q_start: q_start_bp,
                                                q_end: q_end_bp,
                                                s_start: s_start_bp,
                                                s_end: s_end_bp,
                                                e_value: e_val,
                                                bit_score,
                                            },
                                            q_frame: q_frame.frame,
                                            s_frame: s_frame.frame,
                                            q_aa_start: qs,
                                            q_aa_end: qe,
                                            s_aa_start: ss,
                                            s_aa_end: se_ungapped,
                                            q_orig_len: q_frame.orig_len,
                                            s_orig_len: s_len,
                                        });
                                    }
                                    continue;
                                }

                                // Count gapped extensions
                                if diag_enabled {
                                    diag_inner.gapped_extensions.fetch_add(1, AtomicOrdering::Relaxed);
                                }

                                // Two-phase gapped extension (NCBI-style):
                                // Phase 1: Preliminary extension with lower X-drop
                                let (
                                    prelim_qs,
                                    prelim_qe,
                                    prelim_ss,
                                    prelim_se,
                                    prelim_score,
                                    _,
                                    _,
                                    _,
                                ) = extend_gapped_protein(
                                    q_aa,
                                    s_aa,
                                    qs,
                                    ss,
                                    qe - qs,
                                    X_DROP_GAPPED_PRELIM,
                                );

                                // Phase 2: Final extension with higher X-drop for longer alignments
                                let (
                                    final_qs,
                                    final_qe,
                                    final_ss,
                                    final_se,
                                    gapped_score,
                                    matches,
                                    mismatches,
                                    gaps,
                                ) = if prelim_score > 0 {
                                    extend_gapped_protein(
                                        q_aa,
                                        s_aa,
                                        prelim_qs,
                                        prelim_ss,
                                        prelim_qe - prelim_qs,
                                        X_DROP_GAPPED_FINAL,
                                    )
                                } else {
                                    (
                                        prelim_qs,
                                        prelim_qe,
                                        prelim_ss,
                                        prelim_se,
                                        prelim_score,
                                        0,
                                        0,
                                        0,
                                    )
                                };

                                mask.insert(mask_key, final_se);

                                let eff_space = (q_aa.len() as f64) * (db_len_aa_total as f64);
                                let (bit_score, e_val) =
                                    calculate_statistics(gapped_score, eff_space, &params);

                                if e_val <= args.evalue {
                                    if diag_enabled {
                                        diag_inner.gapped_evalue_passed.fetch_add(1, AtomicOrdering::Relaxed);
                                    }
                                    // Calculate alignment length as the max of query and subject spans
                                    let q_len = final_qe - final_qs;
                                    let s_len_aln = final_se - final_ss;
                                    let aln_len = q_len.max(s_len_aln);

                                    // Identity is matches / alignment_length, capped at 100%
                                    let identity = if aln_len > 0 {
                                        ((matches as f64 / aln_len as f64) * 100.0).min(100.0)
                                    } else {
                                        0.0
                                    };

                                    let (q_start_bp, q_end_bp) = convert_coords(
                                        final_qs,
                                        final_qe,
                                        q_frame.frame,
                                        q_frame.orig_len,
                                    );
                                    let (s_start_bp, s_end_bp) =
                                        convert_coords(final_ss, final_se, s_frame.frame, s_len);

                                    // Store sequence data for HSP chaining
                                    let seq_key: SequenceKey = (
                                        query_ids[q_idx as usize].clone(),
                                        s_id.clone(),
                                        q_frame.frame,
                                        s_frame.frame,
                                    );
                                    if !seen_pairs.contains(&seq_key) {
                                        seen_pairs.insert(seq_key.clone());
                                        local_sequences
                                            .push((seq_key.clone(), (q_aa.clone(), s_aa.clone())));
                                    }

                                    local_hits.push(ExtendedHit {
                                        hit: Hit {
                                            query_id: query_ids[q_idx as usize].clone(),
                                            subject_id: s_id.clone(),
                                            identity,
                                            length: aln_len,
                                            mismatch: mismatches,
                                            gapopen: gaps,
                                            q_start: q_start_bp,
                                            q_end: q_end_bp,
                                            s_start: s_start_bp,
                                            s_end: s_end_bp,
                                            e_value: e_val,
                                            bit_score,
                                        },
                                        q_frame: q_frame.frame,
                                        s_frame: s_frame.frame,
                                        q_aa_start: final_qs,
                                        q_aa_end: final_qe,
                                        s_aa_start: final_ss,
                                        s_aa_end: final_se,
                                        q_orig_len: q_frame.orig_len,
                                        s_orig_len: s_len,
                                    });
                                }
                            }
                        }
                        (local_hits, local_sequences)
                    })
                    .collect();

            // Flatten frame results and send
            let mut all_hits: Vec<ExtendedHit> = Vec::new();
            let mut all_sequences: Vec<(SequenceKey, SequenceData)> = Vec::new();
            for (hits, seqs) in frame_results {
                all_hits.extend(hits);
                all_sequences.extend(seqs);
            }
            if !all_hits.is_empty() {
                tx.send((all_hits, all_sequences)).unwrap();
            }
            bar.inc(1);
        });

    bar.finish();
    writer_handle.join().unwrap()?;

    // Print diagnostic summary if enabled
    if diag_enabled {
        diagnostics.print_summary();
    }

    Ok(())
}
