//! The word finder of a window below 0 (`WindowSize` of an imported search strategy; the command
//! line takes 0 or more): NCBI's one-hit steps (`two_hits` is `window > 0`) over a diagonal table
//! that starts at the window and keeps its state from one word finder call to the next.
//!
//! With a window of 0 or more every entry of an earlier call is older than the window, so the
//! other word finders of this folder need no state from the call before. With a window below 0
//! the table's offset starts at the window and its cells start at 0 (the array's and the hash's
//! missing record), so a word with `s_off + window < 0` on an untouched diagonal is never
//! extended, and what an earlier chunk or subject stored (`last_hit`) blocks words of the next:
//! the table is one per search thread and its state carries (`Blast_ExtendWordExit`, which also
//! clears it once the offset passes a quarter of `INT4_MAX`). The search runs on one thread for such a
//! window (`run.rs` `search`), the subjects in NCBI's one-thread order on one table.
//!
//! Arithmetic: NCBI's `Int4` values are held in `i64` here, so a sum that NCBI's `Int4` would
//! overflow (undefined behaviour there) is its value without the wrap and no operation panics;
//! the stored `last_hit` and `level` are NCBI's 31-bit signed fields (`s31`).
//!
//! The steps are those of `pair_extend::extend_pairs_in` (two-stage tables) and
//! `na_scan::scan_na_lookup` (the direct and the na tables) for `two_hits` false, in the same
//! order; the scans that collect the offset pairs are the reference scans (`two_stage_scan`,
//! `na_scan`'s scan loop). No diagnostics are written here.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,90-106,109-185
//! (the table, its clear, `BlastExtendWordNew`, `Blast_ExtendWordExit`),
//! c++/src/algo/blast/core/na_ungapped.c:631-775 (`s_BlastnDiagTableExtendInitialHit`),
//! :799-944 (`s_BlastnDiagHashExtendInitialHit`), :342-449 (the hash's retrieve and insert) and
//! :1597-1695 (`BlastNaWordFinder`), and c++/include/algo/blast/core/blast_extend.h:57-60,62-71
//! (the 31-bit fields).
//!
//! This module is a child of `run` for the items it keeps private.

use super::pair_extend::{word_extend_pair, LutHit, PairExtendInputs, WordExtend};
use super::prelim_subject::{chunk_hsp_list, PrelimSubject, ScanBatch};
use super::two_stage_scan::{scan_two_stage_lookup, TwoStageScanInputs, TwoStageScanState};
use super::{
    extend_hit_ungapped_approx_ncbi, extend_hit_ungapped_exact_ncbi, mask_lookup_index,
    packed_kmer_at_seed_mask, scan_subject_kmers_with_ranges, select_mb_scan_kind, type_of_word,
    x_push_ungapped, GapAlignScratch, OffsetPair, PrelimHit, QueryContext, SubjectChunk,
    SubjectScratch, UngappedHit, COMPRESSION_RATIO,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:164,173
const OFFSET_RESET: i64 = (i32::MAX / 4) as i64;

/// NCBI's 31-bit signed bit-field (`DiagStruct::last_hit`, `DiagHashCell::level`): the low 31
/// bits of `value`, sign-extended.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60,62-71
#[inline(always)]
fn s31(value: i64) -> i32 {
    (((value as u32) << 1) as i32) >> 1
}

/// One cell of the diagonal hash's chain (NCBI's `DiagHashCell`; the hit length and the saved
/// flag are read by the two-hit steps only).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:62-71
#[derive(Clone, Copy, Default)]
struct HashCell {
    diag: i32,
    level: i32,
    next: u32,
}

/// The diagonal hash (NCBI's `BLAST_DiagHash`): the buckets' heads and the chain, whose first
/// cell is the empty one (`occupancy` starts at 1).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:97-106
struct NegHash {
    backbone: Vec<u32>,
    chain: Vec<HashCell>,
    /// The window of the stale-cell rule of the insert: `window + Delta + 1`.
    stale_window: i64,
}

enum Body {
    /// The diagonal array: one 31-bit `last_hit` per diagonal; the table's length and mask.
    Array {
        levels: Vec<i32>,
        length: i64,
        mask: i64,
    },
    Hash(NegHash),
}

/// The diagonal table of a window below 0 (NCBI's `BLAST_DiagTable` or `BLAST_DiagHash`, with the
/// offset that starts at the window).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,109-157
pub(super) struct NegTable {
    window: i64,
    offset: i64,
    body: Body,
}

impl NegTable {
    /// A new table (`BlastExtendWordNew`): the array of `array_length` cells (a power of 2) when
    /// it is given, else the hash with the insert's stale window `stale_window`. The offset is
    /// the window and every cell is 0 (calloc'd).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,122-157
    pub(super) fn new(window: i32, array_length: Option<usize>, stale_window: i64) -> Self {
        let body = match array_length {
            Some(length) => Body::Array {
                levels: vec![0; length],
                length: length as i64,
                mask: length as i64 - 1,
            },
            None => Body::Hash(NegHash {
                backbone: vec![0; DIAG_HASH_BUCKETS],
                chain: vec![HashCell::default(); 1],
                stale_window,
            }),
        };
        NegTable {
            window: i64::from(window),
            offset: i64::from(window),
            body,
        }
    }

    /// Whether this table is the one that a search with these values makes (a scratch that
    /// outlived its batch gets a new table).
    fn fits(&self, window: i32, array_length: Option<usize>, stale_window: i64) -> bool {
        self.window == i64::from(window)
            && match (&self.body, array_length) {
                (Body::Array { length, .. }, Some(wanted)) => *length == wanted as i64,
                (Body::Hash(hash), None) => hash.stale_window == stale_window,
                _ => false,
            }
    }

    /// The table's offset (NCBI's `offset`).
    #[inline(always)]
    fn offset(&self) -> i64 {
        self.offset
    }

    /// The last hit of diagonal `diag` (`s_off - q_off`): the array's cell of the diagonal masked
    /// to the table's length, or the hash's record (0 when it has none).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:655-659,820-831
    #[inline(always)]
    fn last_hit(&self, diag: i64) -> i32 {
        match &self.body {
            Body::Array {
                levels,
                length,
                mask,
            } => levels[((diag + *length) & *mask) as usize],
            Body::Hash(hash) => {
                let mut index = hash.backbone[bucket(diag)];
                while index != 0 {
                    let cell = &hash.chain[index as usize];
                    if cell.diag == diag as i32 {
                        return cell.level;
                    }
                    index = cell.next;
                }
                0
            }
        }
    }

    /// The seed at subject offset `s_off` on diagonal `diag` is within the explored area when
    /// `s_off + offset < last_hit`: `None` then, else the seed's `s_off + offset`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:670-671,830-836
    #[inline(always)]
    fn seed_pos(&self, s_off: usize, diag: i64) -> Option<i64> {
        let s_off_pos = s_off as i64 + self.offset;
        (s_off_pos >= i64::from(self.last_hit(diag))).then_some(s_off_pos)
    }

    /// Stores the end position of the seed's hit on `diag` (the table's `last_hit`, 31 bits), the
    /// seed having started at `s_off_pos`; the hash replaces a cell of the bucket whose level
    /// lies more than the stale window before `s_off_pos`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:768-772 (array),
    /// :396-449,939-941 (hash)
    #[inline(always)]
    fn store(&mut self, diag: i64, end_pos: i64, s_off_pos: i64) {
        let level = s31(end_pos);
        match &mut self.body {
            Body::Array {
                levels,
                length,
                mask,
            } => levels[((diag + *length) & *mask) as usize] = level,
            Body::Hash(hash) => {
                let diag = diag as i32;
                let bucket = bucket(i64::from(diag));
                let mut index = hash.backbone[bucket];
                while index != 0 {
                    let cell = &mut hash.chain[index as usize];
                    if cell.diag == diag {
                        cell.level = level;
                        return;
                    }
                    if s_off_pos - i64::from(cell.level) > hash.stale_window {
                        cell.diag = diag;
                        cell.level = level;
                        return;
                    }
                    index = cell.next;
                }
                let cell_index = hash.chain.len() as u32;
                hash.chain.push(HashCell {
                    diag,
                    level,
                    next: hash.backbone[bucket],
                });
                hash.backbone[bucket] = cell_index;
            }
        }
    }

    /// The end of a word finder call on a subject (chunk) of `subject_length` letters: the offset
    /// passes on by the length and the window, or, once it is a quarter of `INT4_MAX` or more, the table
    /// is cleared and the offset starts at the window again (the array's cells at
    /// `-window`, the hash empty).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:90-106,162-185
    pub(super) fn exit(&mut self, subject_length: usize) {
        if self.offset >= OFFSET_RESET {
            self.offset = self.window;
            match &mut self.body {
                Body::Array { levels, .. } => levels.fill(s31(-self.window)),
                Body::Hash(hash) => {
                    hash.chain.truncate(1);
                    hash.backbone.fill(0);
                }
            }
        } else {
            self.offset += subject_length as i64 + self.window;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:123-124
const DIAG_HASH_BUCKETS: usize = 512;

/// The hash's bucket of a diagonal.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-363
#[inline(always)]
fn bucket(diag: i64) -> usize {
    ((diag as i32 as u32).wrapping_mul(0x9E37_0001) % DIAG_HASH_BUCKETS as u32) as usize
}

/// The ungapped extension of a seed that the table let through, the cutoff test, the save of the
/// initial hit and the table's update: `last_hit` is the end of the saved hit's extension, else
/// the end of the seed's word.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:728-773 (array),
/// :896-943 (hash)
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn extend_and_store(
    scan: &ScanBatch<'_>,
    table: &mut NegTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    search_seq_packed: &[u8],
    s_len: usize,
    context_idx: usize,
    ctx: &QueryContext,
    (q_off, s_off, s_end): (usize, usize, usize),
    diag: i64,
    (s_off_pos, s_end_pos): (i64, i64),
) {
    let x_dropoff = scan.x_dropoff_scores[context_idx];
    let ungapped = if scan.exact_ungapped {
        extend_hit_ungapped_exact_ncbi(
            scan.encoded_query_concat_blastna,
            search_seq_packed,
            q_off,
            s_off,
            s_len,
            x_dropoff,
            scan.score_matrix,
        )
    } else {
        extend_hit_ungapped_approx_ncbi(
            scan.encoded_query_concat_blastna,
            search_seq_packed,
            q_off,
            s_off,
            s_end,
            s_len,
            x_dropoff,
            scan.nucl_score_table,
            scan.reduced_cutoff_scores[context_idx],
            scan.score_matrix,
        )
    };
    let end_pos = if ungapped.score >= scan.cutoff_scores[context_idx] {
        let saved_end = ungapped.s_start + ungapped.length;
        x_push_ungapped(
            ungapped_hits,
            UngappedHit::saved(context_idx as u32, ctx, q_off, s_off, &ungapped),
        );
        saved_end as i64 + table.offset()
    } else {
        s_end_pos
    };
    table.store(diag, end_pos, s_off_pos);
}

/// The steps of a batch of offset pairs of a two-stage table (`pair_extend::extend_offset_pairs`
/// for a window below 0): for each pair in order its context and bounds, the word extension, the
/// diagonal's last hit and the test against it, the word's type (masks, lookup word), the
/// ungapped extension, the cutoff test, the save and the table's update. Then the batch is empty.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1093-1173 (`s_BlastNaExtend`),
/// :631-775 and :799-944
pub(super) fn extend_offset_pairs(
    scan: &ScanBatch<'_>,
    x: &PairExtendInputs<'_>,
    table: &mut NegTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    offset_pairs: &[OffsetPair],
    offset_pairs_len: &mut usize,
) {
    let &PairExtendInputs {
        two_stage,
        small_na_word,
        subject_masked,
        search_seq_packed,
        s_len,
        ..
    } = x;
    let word_length = scan.word_length;
    let lut_word_length = scan.lut_word_length;
    let word = WordExtend {
        small_na_word,
        small_na_aligned_one_byte: scan.small_na_aligned_one_byte,
        subject_masked,
        word_length,
        lut_word_length,
        search_seq_packed,
        encoded_query_concat_blastna_with_sentinels: scan
            .encoded_query_concat_blastna_with_sentinels,
        debug_enabled: false,
    };
    let (mut left_iters, mut right_iters) = (0usize, 0usize);
    for pair in &offset_pairs[..*offset_pairs_len] {
        let q_off0 = pair.q_off;
        let kmer_start = pair.s_off;
        let s_range = pair.s_range;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:730-733
        let context_idx = scan.query_context_index.context_for_offset(q_off0);
        let ctx = &scan.query_contexts[context_idx];
        let q_context_start = ctx.query_offset as usize;
        let q_context_end = q_context_start + ctx.len;
        if q_off0 + lut_word_length > q_context_end || kmer_start + lut_word_length > s_len {
            continue;
        }
        let Some((q_ext_start, s_ext_start, _, _)) = word_extend_pair(
            &word,
            LutHit {
                q_off0,
                kmer_start,
                q_context_start,
                q_context_end,
                s_range,
            },
            &mut left_iters,
            &mut right_iters,
        ) else {
            continue;
        };
        // The diagonal and its last hit after the word extension, as the two-stage table's
        // steps take them.
        let diag = s_ext_start as i64 - q_ext_start as i64;
        let Some(s_off_pos) = table.seed_pos(s_ext_start, diag) else {
            continue;
        };
        let mut is_seed_masked = |s_pos: usize, q_pos: usize| -> bool {
            if s_pos + lut_word_length > s_len {
                return true;
            }
            let kmer = mask_lookup_index(
                packed_kmer_at_seed_mask(search_seq_packed, s_pos, lut_word_length),
                lut_word_length,
            );
            if !two_stage.has_hits(kmer) {
                return true;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1027-1034
            !two_stage.contains_hit(kmer, (q_pos + 1) as u32)
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:718-726
        let (word_type, extended, q_off, s_off) = type_of_word(
            q_ext_start,
            s_ext_start,
            small_na_word.is_some() || !ctx.masks.is_empty(),
            q_context_end,
            s_range,
            word_length,
            lut_word_length,
            false,
            &mut is_seed_masked,
        );
        if word_type == 0 {
            continue;
        }
        let s_end = s_ext_start + word_length + extended;
        let s_end_pos = s_end as i64 + table.offset();
        extend_and_store(
            scan,
            table,
            ungapped_hits,
            search_seq_packed,
            s_len,
            context_idx,
            ctx,
            (q_off, s_off, s_end),
            diag,
            (s_off_pos, s_end_pos),
        );
    }
    *offset_pairs_len = 0;
}

/// The scan of one chunk with a table that is not two-stage (the direct table or the na table)
/// and the steps of each seed it finds, for a window below 0: the lookup width is the word size,
/// so a seed has no word extension and its word is the lookup word.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:963-1030
/// (`s_BlastNaExtendDirect`), :631-775 and :799-944
fn scan_na_lookup(
    scan: &ScanBatch<'_>,
    table: &mut NegTable,
    ungapped_hits: &mut Vec<UngappedHit>,
    search_seq_packed: &[u8],
    s_len: usize,
    subject_seq_ranges: &[(i32, i32)],
    subject_masked: bool,
) {
    let safe_k = scan.safe_k;
    let mb_scan_kind = select_mb_scan_kind(safe_k, scan.scan_step, subject_masked);
    scan_subject_kmers_with_ranges(
        search_seq_packed,
        s_len,
        safe_k,
        safe_k,
        scan.scan_step,
        subject_seq_ranges,
        subject_masked,
        mb_scan_kind,
        |kmer_start, _s_range, current_kmer| {
            let matches_slice: &[u32] = if scan.use_direct_lookup {
                scan.pv_direct_lookup_ref
                    .map(|pv_dl| pv_dl.get_hits_checked(current_kmer))
                    .unwrap_or(&[])
            } else {
                scan.na_lookup_ref
                    .map(|na_lt| na_lt.get_hits_checked(current_kmer))
                    .unwrap_or(&[])
            };
            for &q_off_1 in matches_slice {
                if q_off_1 == 0 {
                    continue;
                }
                let q_off0 = q_off_1 as usize - 1;
                let context_idx = scan.query_context_index.context_for_offset(q_off0);
                let ctx = &scan.query_contexts[context_idx];
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:828-836
                let diag = kmer_start as i64 - q_off0 as i64;
                let Some(s_off_pos) = table.seed_pos(kmer_start, diag) else {
                    continue;
                };
                let s_end = kmer_start + safe_k;
                let s_end_pos = s_end as i64 + table.offset();
                extend_and_store(
                    scan,
                    table,
                    ungapped_hits,
                    search_seq_packed,
                    s_len,
                    context_idx,
                    ctx,
                    (q_off0, kmer_start, s_end),
                    diag,
                    (s_off_pos, s_end_pos),
                );
            }
        },
    );
}

/// The search of one chunk of a subject for a window below 0: `prelim_subject::prelim_chunk_hits`
/// with the reference scans and the steps of this module over the search's one diagonal table
/// (made by the first chunk, kept in the scratch), whose offset passes on after the chunk (also
/// a chunk shorter than the word, which has no scan). The chunk's HSP list is the other word
/// finders' (`prelim_subject::chunk_hsp_list`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:487-500,522-527,555,560-563
/// (one turn of the chunk loop of `s_BlastSearchEngineOneContext`) and
/// c++/src/algo/blast/core/na_ungapped.c:1597-1695
pub(super) fn neg_chunk_hits<const UNGAPPED: bool>(
    scan: &ScanBatch<'_>,
    subject: &PrelimSubject<'_>,
    chunk: &SubjectChunk,
    soft_ranges: &[(i32, i32)],
    gap_scratch: &mut GapAlignScratch,
    subject_scratch: &mut SubjectScratch,
    reuse_prelim_hits: bool,
) -> Vec<PrelimHit> {
    let batch = scan.batch;
    subject_scratch.hits.clear();
    subject_scratch.prelim_hits.clear();
    subject_scratch.ungapped_hits.clear();
    let s_len = chunk.length;
    let subject_seq_ranges = chunk.seq_ranges.as_slice(
        soft_ranges,
        &mut subject_scratch.seq_ranges_scratch,
        chunk.offset as i32,
        s_len as i32,
    );
    let packed_offset = chunk.offset / COMPRESSION_RATIO;
    let s_seq_packed = &subject.s_seq_packed_full[packed_offset..];

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,109-157
    // The table: the array of the query block's length plus the window (a power of 2, which
    // can be shorter than the query), or the hash for a query block over 8000 letters whose
    // stale window is `window + Delta + 1`, `Delta` not clamped at 0.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:822-823,939-941
    let window = scan.neg_window;
    let array_length = scan.use_array_indexing.then_some(scan.diag_table_length);
    let delta = (scan.scan_range as i64).min(i64::from(window) - scan.word_length as i64);
    let stale_window = i64::from(window) + delta + 1;
    let fits = matches!(
        &subject_scratch.neg,
        Some(table) if table.fits(window, array_length, stale_window)
    );
    if !fits {
        subject_scratch.neg = Some(NegTable::new(window, array_length, stale_window));
    }
    let Some(table) = subject_scratch.neg.as_mut() else {
        return Vec::new();
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-176
    if s_len < batch.effective_word_size {
        table.exit(s_len);
        return Vec::new();
    }
    if let Some(two_stage) = scan.two_stage_lookup_ref {
        let (mut positions, mut seeds, mut window_seeds, mut low, mut gapped) =
            (0usize, 0usize, 0usize, 0usize, 0usize);
        scan_two_stage_lookup(
            &TwoStageScanInputs {
                scan,
                two_stage,
                search_seq_packed: s_seq_packed,
                s_len,
                subject_seq_ranges,
                subject_masked: chunk.masked,
                chunk,
                diag_offset: 0,
                x_lv_ref: false,
                s_idx: subject.s_idx,
                s_id: subject.s_id,
            },
            &mut TwoStageScanState {
                ungapped_hits: &mut subject_scratch.ungapped_hits,
                hit_level_array: &mut subject_scratch.hit_level_array,
                hit_len_array: &mut subject_scratch.hit_len_array,
                diag_hash: &mut subject_scratch.diag_hash,
                offset_pairs: &mut subject_scratch.offset_pairs,
                neg: Some(&mut *table),
                dbg_total_s_positions: &mut positions,
                dbg_seeds_found: &mut seeds,
                dbg_window_seeds: &mut window_seeds,
                dbg_ungapped_low: &mut low,
                dbg_gapped_attempted: &mut gapped,
            },
        );
    } else {
        scan_na_lookup(
            scan,
            table,
            &mut subject_scratch.ungapped_hits,
            s_seq_packed,
            s_len,
            subject_seq_ranges,
            chunk.masked,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-500,560-563
    if let Some(index) = batch.mb_index {
        index.note_scanned_hits(subject.s_idx as u32, subject_scratch.ungapped_hits.len());
    }
    let chunk_prelim_hits = chunk_hsp_list::<UNGAPPED>(
        batch,
        subject,
        chunk,
        s_seq_packed,
        gap_scratch,
        subject_scratch,
        reuse_prelim_hits,
        batch.good_init_extends,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-176
    if let Some(table) = subject_scratch.neg.as_mut() {
        table.exit(s_len);
    }
    chunk_prelim_hits
}

#[cfg(test)]
#[path = "neg_window_tests.rs"]
mod tests;
