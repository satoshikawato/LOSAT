//! The word finders of RPS-BLAST (`s_BlastRPSWordFinder_TwoHit`, `s_BlastRPSWordFinder_OneHit`):
//! the scan's hits, bucket by bucket, on NCBI's diagonal table of the concatenation, extended
//! with PS4a's PSSM ungapped kernels on the whole concatenation.
//!
//! The diagonal table has `L` entries, `L` the smallest power of two at least the concatenated
//! length `Lc` plus the window (`s_BlastDiagTableNew` of `concat_db.length`,
//! `RpsRun::concat_length`: the volume's rows, or an alias file's `LENGTH` and `NSEQ` over one
//! volume, so that `L` may be below the rows or above them), and is indexed with `mask = L - 1`:
//! the two-hit finder reads `(q - s) & mask`, the one-hit finder `(s - q) & mask`. Where a query
//! is longer than `L` less the concatenation, two diagonals share an entry, as in NCBI.
//!
//! NCBI keeps one table for the queries of a batch, moving `diag_offset` on by the query's
//! length and the window after each (`Blast_ExtendWordExit`); an entry that an earlier query
//! wrote then reads as "no hit on this diagonal yet" in both finders, as an entry never written
//! does. A worker's table here follows the same rule for the queries it searches, so every query
//! gets NCBI's hits whatever queries came before it.

use super::scan::{
    q_place, s_place, scan_group0, scan_subject, RpsBlockPool, RpsBuckets, RpsHits, RpsReadCursor,
    RpsTable, MAX_SCAN_HITS,
};
use super::RpsRows;
use crate::algorithm::blastp::extension::BlastpUngappedData;
use crate::algorithm::pssm_search::rows::PssmRows;
use crate::algorithm::pssm_search::ungapped::{pssm_extend_one_hit, pssm_extend_two_hit};
use crate::seqdb::rps::{RPS_BUCKET_BITS, RPS_BUCKET_SIZE, RPS_WORDSIZE};

/// One entry of the diagonal table (`DiagStruct`: a 31-bit signed `last_hit` and a 1-bit
/// `flag`), packed in a `u32` as NCBI's bit fields are. The table is a `Vec<u32>` so that a
/// zeroed one comes from the allocator's zeroed pages (`vec![0u32; n]`, `calloc`), as NCBI's.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
trait RpsDiag {
    fn last_hit(&self) -> i32;
    fn flag(&self) -> bool;
    fn set(&mut self, last_hit: i32, flag: bool);
}

impl RpsDiag for u32 {
    #[inline(always)]
    fn last_hit(&self) -> i32 {
        ((*self << 1) as i32) >> 1
    }

    #[inline(always)]
    fn flag(&self) -> bool {
        *self >> 31 != 0
    }

    #[inline(always)]
    fn set(&mut self, last_hit: i32, flag: bool) {
        *self = (last_hit as u32 & 0x7fff_ffff) | (u32::from(flag) << 31);
    }
}

/// The diagonal table of one worker (`BLAST_DiagTable` of `ewp`).
pub(crate) struct RpsDiagTable {
    array: Vec<u32>,
    mask: i32,
    offset: i32,
    window: i32,
    /// NCBI's full table length (`diag_array_length(concat_length, window)`), or `0` when
    /// initialized directly via `reset_for` (tests that fix a table length).
    ncbi_len: i32,
    /// Concatenated rows of the current volume (`prepare_for_volume`).
    num_rows: usize,
}

impl RpsDiagTable {
    /// An empty table (sized by `reset_for` or `prepare_for_volume`).
    pub(crate) fn new() -> Self {
        Self {
            array: Vec::new(),
            mask: 0,
            offset: 0,
            window: 0,
            ncbi_len: 0,
            num_rows: 0,
        }
    }

    /// `BlastExtendWordNew` for `concat_length`: a zeroed table of the unit's length. The
    /// table is allocated zeroed and not written (NCBI's `calloc`): only the pages of the
    /// diagonals that the queries hit become resident, whatever the window.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-62
    /// The table has `length` entries (`diag_table_length`).
    #[cfg(test)]
    pub(crate) fn reset_for(&mut self, length: i32, window: i32) {
        self.array = vec![0u32; length as usize];
        self.mask = length - 1;
        self.offset = window;
        self.window = window;
        self.ncbi_len = 0;
        self.num_rows = 0;
    }

    /// The table for a volume of `num_rows` concatenated rows, concatenated length
    /// `concat_length`, word window `window`, and queries of at most `longest` letters.
    ///
    /// When `num_rows + longest <= diag_array_length(concat_length, window)`, no two distinct
    /// diagonals collide modulo NCBI's table length, so the word finder can process the volume
    /// bucket by bucket on a sliding ring of `(RPS_BUCKET_SIZE + subject.len()).next_power_of_two()`
    /// entries (16 KiB for typical proteins, vs. 32 MiB on `cdd_delta`), keeping the exact same
    /// per-diagonal hit sequence and sorted `init_hsps` as NCBI. Where diagonals wrap modulo
    /// NCBI's length, the full NCBI length is kept.
    pub(crate) fn prepare_for_volume(
        &mut self,
        num_rows: usize,
        concat_length: i32,
        window: i32,
        longest: usize,
    ) {
        let ncbi = diag_array_length(concat_length, window);
        let no_wrap = (num_rows as u64).saturating_add(longest as u64) <= ncbi as u64;
        let length = if no_wrap {
            bucket_ring_length(longest).min(ncbi as u64) as i32
        } else {
            diag_table_length(concat_length, window, num_rows, longest)
        };
        if self.array.len() != length as usize || self.window != window {
            self.array = vec![0u32; length as usize];
            self.mask = length - 1;
            self.offset = window;
            self.window = window;
        }
        self.ncbi_len = ncbi;
        self.num_rows = num_rows;
    }

    /// `Blast_ExtendWordExit`: the offset after a query of `subject_length` letters.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-178
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:102-104 (`s_BlastDiagClear`)
    fn exit(&mut self, subject_length: i32) {
        if self.offset >= i32::MAX / 4 {
            self.offset = self.window;
            for entry in &mut self.array {
                entry.set(-self.window, false);
            }
        } else {
            self.offset += subject_length + self.window;
        }
    }
}

/// Power-of-two ring length for sliding bucket-by-bucket over a query of `subject_len` letters:
/// spans one bucket (`RPS_BUCKET_SIZE` rows) plus the `subject_len - 1` diagonals shared with
/// earlier buckets, so every diagonal in the current bucket maps to a unique slot and every new
/// bucket's `RPS_BUCKET_SIZE` slots are disjoint from the carried-over overlap.
#[inline]
pub(crate) fn bucket_ring_length(subject_len: usize) -> u64 {
    (RPS_BUCKET_SIZE as u64)
        .saturating_add(subject_len as u64)
        .next_power_of_two()
}

/// NCBI's length of the diagonal table of a concatenated length `concat_length` (`Lc`, also 0
/// or negative) and the word window `window` (`s_BlastDiagTableNew`, cited at
/// `RpsDiagTable::reset_for`): the power of 2 from 1 up that is at least their sum (1 entry where
/// the sum is 1 or less; the sum is at most 2^30, `RpsUnit::check_engine_length`).
pub(crate) fn diag_array_length(concat_length: i32, window: i32) -> i32 {
    let mut length = 1i32;
    while length < concat_length + window {
        length <<= 1;
    }
    length
}

/// The entries of the diagonal table of a search of `rows` concatenated rows and queries (frames)
/// of at most `longest` letters with NCBI's length `diag_array_length(concat_length, window)`:
/// that length, or the power of 2 from 1 up that is at least `rows + longest` where that is
/// smaller (review B of D1a-4, SF2). The word finders' diagonal of a hit is `(q - s) & mask` with
/// `q` within the rows and `s` within the query: the `rows + longest - 1` diagonals of a query are
/// apart in a table of either length, so no two share an entry in NCBI's table where it is the
/// longer one, nor in this one, and each entry reads and writes as NCBI's (the offset of `exit`
/// and its clear are the same for every entry). Where NCBI's table is the shorter one (an `Lc`
/// below the rows) its diagonals share entries as NCBI's do. One worker's table of NCBI's length
/// is 4 GiB at `Lc + window` of 2^30, where NCBI searches the volume with one table.
pub(crate) fn diag_table_length(
    concat_length: i32,
    window: i32,
    rows: usize,
    longest: usize,
) -> i32 {
    let ncbi = diag_array_length(concat_length, window);
    let apart = (rows as u64)
        .saturating_add(longest as u64)
        .next_power_of_two();
    if apart < ncbi as u64 {
        apart as i32
    } else {
        ncbi
    }
}

/// NCBI's `BlastInitHSP` with its `BlastUngappedData` (`BlastSaveInitHsp`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RpsInitHsp {
    pub(crate) q_off: i32,
    pub(crate) s_off: i32,
    pub(crate) ungapped: BlastpUngappedData,
}

/// What the word finder of one (query, volume) reads.
pub(crate) struct RpsWordParams<'a> {
    pub(crate) rows: RpsRows<'a>,
    pub(crate) table: RpsTable<'a>,
    /// `-window_size`; 0 is the one-hit finder.
    pub(crate) window: i32,
    /// `cutoffs->cutoff_score` of the query.
    pub(crate) cutoff: i32,
    /// `cutoffs->x_dropoff` of the query.
    pub(crate) dropoff: i32,
}

/// The scratch of one worker's word finder: its diagonal table and buckets.
pub(crate) struct RpsWordScratch {
    pub(crate) diag: RpsDiagTable,
    pub(crate) buckets: RpsBuckets,
}

impl RpsWordScratch {
    /// A scratch without a table, its buckets' blocks from the run's `pool` (sized by
    /// `prepare_for`).
    pub(crate) fn new(pool: std::sync::Arc<RpsBlockPool>) -> Self {
        Self {
            diag: RpsDiagTable::new(),
            buckets: RpsBuckets::new(0, pool),
        }
    }

    /// The table of a search of the concatenated length `concat_length` (`Lc`,
    /// `RpsRun::concat_length`: `BlastExtendWordNew` for `concat_db.length`) whose queries (frames)
    /// have at most `longest` letters (`diag_table_length`) and the buckets of a volume of
    /// `num_rows` concatenated rows (`RpsDiagTable::prepare_for`, `RpsBuckets::reset_for`): what
    /// the scratch kept stays.
    pub(crate) fn prepare_for(
        &mut self,
        num_rows: usize,
        concat_length: i32,
        window: i32,
        longest: usize,
    ) {
        self.diag
            .prepare_for_volume(num_rows, concat_length, window, longest);
        self.buckets.reset_for(num_rows);
    }
}

/// `BlastRPSWordFinder` on `subject` (the query's letters): the initial HSPs into `init_hsps`
/// (not sorted; the caller sorts them, `Blast_InitHitListSortByScore`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:265-283
///
/// The extensions read the rows as `RpsRowsWhole` where `Lc` is the rows' count (every search
/// but one over an alias file whose `LENGTH` and `NSEQ` give another `Lc`), else as `RpsRows`
/// with `Lc` and NCBI's step past the rows, out of line (`word_finder_other_length`): the same
/// scan, whose common instantiation carries neither the bound nor the step (fix rounds 6 and 8
/// of line RG D1a-4).
pub(crate) fn word_finder(
    params: &RpsWordParams<'_>,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    match params.rows.whole() {
        Some(rows) => scan(params, &rows, scratch, subject, init_hsps),
        None => word_finder_other_length(params, scratch, subject, init_hsps),
    }
    scratch.diag.exit(subject.len() as i32);
}

/// `word_finder`'s scan with an `Lc` that is not the rows' count: the extensions' bound below
/// the rows (`RpsRows::length`) or NCBI's step past them (`RpsRows::past`, an `Lc` above them).
#[cold]
#[inline(never)]
fn word_finder_other_length(
    params: &RpsWordParams<'_>,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    scan(params, &params.rows, scratch, subject, init_hsps);
}

/// The two-hit or the one-hit finder (`-window_size` 0), whose extensions read `rows`.
///
/// Uses the bucket-sliding diagonal ring whenever NCBI's full table has no diagonal collisions
/// (`num_rows + subject.len() <= ncbi_len`) and the ring is strictly smaller than `ncbi_len`;
/// otherwise uses NCBI's full table so any diagonal wrap-around matches NCBI exactly.
#[inline(always)]
fn scan<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    if subject.len() < RPS_WORDSIZE {
        return;
    }
    let ncbi = scratch.diag.ncbi_len;
    let ring = bucket_ring_length(subject.len());
    let use_sliding = ncbi > 0
        && (scratch.diag.num_rows as u64).saturating_add(subject.len() as u64) <= ncbi as u64
        && ring < ncbi as u64;
    if use_sliding {
        let ring_len = ring as usize;
        if scratch.diag.array.len() < ring_len {
            scratch.diag.array.resize(ring_len, 0);
        }
        if params.window > 0 {
            two_hit_sliding(params, rows, scratch, subject, ring_len, init_hsps);
        } else {
            one_hit_sliding(params, rows, scratch, subject, ring_len, init_hsps);
        }
    } else if params.window > 0 {
        two_hit(params, rows, scratch, subject, init_hsps);
    } else {
        one_hit(params, rows, scratch, subject, init_hsps);
    }
}

/// Advances the bucket-sliding ring from `*prev` to `bucket`, clearing the `RPS_BUCKET_SIZE` new
/// diagonal slots for every bucket in `(*prev + 1)..=bucket`. `*prev == usize::MAX` marks the
/// first non-empty bucket of a query (whose ring entries are already stale via `diag_offset`).
#[inline(never)]
fn advance_bucket_cold(table: &mut [u32], mask: usize, prev: &mut usize, bucket: usize) {
    let prev_b = *prev;
    *prev = bucket;
    if prev_b == usize::MAX {
        return;
    }
    let gap = bucket - prev_b;
    let slots = table.len() >> RPS_BUCKET_BITS;
    if gap >= slots {
        table.fill(0);
    } else {
        for b in (prev_b + 1)..=bucket {
            let start = (b << RPS_BUCKET_BITS) & mask;
            table[start..start + RPS_BUCKET_SIZE].fill(0);
        }
    }
}

/// Processes one block of hits for `two_hit` (both full-table and bucket-sliding modes).
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn two_hit_block<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    subject: &[u8],
    table: &mut [u32],
    diag_mask: i32,
    diag_offset: i32,
    hits: RpsHits<'_>,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let window = params.window;
    // Slice to `mask + 1` so `(coord & diag_mask) as usize` is statically in bounds.
    let table = &mut table[..=diag_mask as usize];
    let RpsHits {
        entries,
        q_base,
        s_base,
    } = hits;
    // Each hit as its places in the block (`RpsHits`): `query_offset - subject_offset` is
    // `q - s + diag_base` and `subject_offset + diag_offset` is `s + s_diag_offset`, the
    // same 32-bit values; the extension reads the offsets whole.
    let (q_base, s_base) = (q_base as i32, s_base as i32);
    let diag_base = q_base.wrapping_sub(s_base);
    let s_diag_offset = diag_offset.wrapping_add(s_base);
    // An index loop (the slice iterator's step took six instructions per hit here).
    let mut next = 0;
    while next < entries.len() {
        let hit_entry = entries[next];
        next += 1;
        let q = q_place(hit_entry) as i32;
        let s = s_place(hit_entry) as i32;
        let diag_coord = (q.wrapping_sub(s).wrapping_add(diag_base) & diag_mask) as usize;
        let entry = &mut table[diag_coord];
        if entry.flag() {
            if s + s_diag_offset < entry.last_hit() {
                continue;
            }
            entry.set(s + s_diag_offset, false);
        } else {
            // `last_hit - diag_offset` less the segment's first offset.
            let last_hit = entry.last_hit() - s_diag_offset;
            let diff = s - last_hit;
            if diff >= window {
                entry.set(s + s_diag_offset, entry.flag());
                continue;
            }
            if diff < wordsize {
                continue;
            }
            let query_offset = q_base + q;
            let subject_offset = s_base + s;
            let hit = pssm_extend_two_hit(
                rows,
                subject,
                s_base + last_hit + wordsize,
                subject_offset,
                query_offset,
                params.dropoff,
                wordsize,
            );
            if hit.ungapped_data.score >= params.cutoff {
                init_hsps.push(RpsInitHsp {
                    q_off: query_offset,
                    s_off: subject_offset,
                    ungapped: hit.ungapped_data,
                });
            }
            if hit.right_extend {
                entry.set(hit.s_last_off - (wordsize - 1) + diag_offset, true);
            } else {
                entry.set(s + s_diag_offset, entry.flag());
            }
        }
    }
}

/// Processes one block of hits for `one_hit`: when `Q_MINUS_S` is false, indexes `(s - q) & mask`
/// (NCBI's full-table sign); when true (bucket-sliding mode, where no diagonals collide), indexes
/// `(q - s) & mask` so bucket transitions clear the same contiguous `[b * B, (b + 1) * B)` slots
/// as `two_hit`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn one_hit_block<R: PssmRows, const Q_MINUS_S: bool>(
    params: &RpsWordParams<'_>,
    rows: &R,
    subject: &[u8],
    table: &mut [u32],
    diag_mask: i32,
    diag_offset: i32,
    hits: RpsHits<'_>,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    // Slice to `mask + 1` so `(coord & diag_mask) as usize` is statically in bounds.
    let table = &mut table[..=diag_mask as usize];
    let RpsHits {
        entries,
        q_base,
        s_base,
    } = hits;
    let (q_base, s_base) = (q_base as i32, s_base as i32);
    let diag_base = q_base.wrapping_sub(s_base);
    let s_diag_offset = diag_offset.wrapping_add(s_base);
    // An index loop (the slice iterator's step took six instructions per hit here).
    let mut next = 0;
    while next < entries.len() {
        let hit_entry = entries[next];
        next += 1;
        let q = q_place(hit_entry) as i32;
        let s = s_place(hit_entry) as i32;
        let diag_coord = if Q_MINUS_S {
            (q.wrapping_sub(s).wrapping_add(diag_base) & diag_mask) as usize
        } else {
            (s.wrapping_sub(q).wrapping_sub(diag_base) & diag_mask) as usize
        };
        let entry = &mut table[diag_coord];
        let diff = s - (entry.last_hit() - s_diag_offset);
        if diff >= 0 {
            let query_offset = q_base + q;
            let subject_offset = s_base + s;
            let hit = pssm_extend_one_hit(
                rows,
                subject,
                subject_offset,
                query_offset,
                params.dropoff,
                wordsize,
            );
            if hit.ungapped_data.score >= params.cutoff {
                init_hsps.push(RpsInitHsp {
                    q_off: query_offset,
                    s_off: subject_offset,
                    ungapped: hit.ungapped_data,
                });
            }
            entry.set(hit.s_last_off - (wordsize - 1) + diag_offset, entry.flag());
        }
    }
}

/// Group-major, bucket-sliding word finder (both two-hit and one-hit, parameterized by `TWO_HIT`):
/// scans each group of buckets over the entire subject and slides the `ring_len` diagonal ring
/// bucket by bucket. Because no two distinct diagonals collide in NCBI's table when this runs,
/// every diagonal sees the exact same hit sequence in ascending `s` order as NCBI's scan, and
/// `sort_init_hsps_by_score` produces the exact same sorted `init_hsps`.
#[inline(always)]
fn run_sliding<R: PssmRows, const TWO_HIT: bool>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    ring_len: usize,
    max_hits: i32,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let last_offset = subject.len() as i32 - wordsize;
    let RpsWordScratch { diag, buckets } = scratch;
    let diag_offset = diag.offset;
    let diag_mask = (ring_len - 1) as i32;
    let table = &mut diag.array[..=diag_mask as usize];
    let mut prev_bucket = usize::MAX;

    // Group 0: scan `subject` for group 0's hits, recording overflow words and inline groups for
    // groups 1..num_groups.
    let mut s = 0i32;
    scan_group0(&params.table, subject, &mut s, buckets, max_hits);
    let cursor = buckets.start_group(&params.table, subject, 0);
    if s > last_offset {
        drain_sliding_group::<R, TWO_HIT>(
            params,
            rows,
            subject,
            table,
            diag_mask,
            diag_offset,
            buckets,
            cursor,
            &mut prev_bucket,
            init_hsps,
        );
    } else {
        finish_multipass_group::<R, TWO_HIT>(
            params,
            rows,
            subject,
            table,
            diag_mask,
            diag_offset,
            buckets,
            0,
            cursor,
            s,
            0,
            last_offset,
            max_hits,
            &mut prev_bucket,
            init_hsps,
        );
    }

    // Remaining groups in ascending order.
    let num_groups = buckets.groups();
    for group in 1..num_groups {
        if !buckets.needs_pass(group) {
            continue;
        }
        let mut s = 0i32;
        let mut word_pos = 0usize;
        let cursor = buckets.fill_group_bounded(
            &params.table,
            subject,
            group,
            &mut s,
            &mut word_pos,
            max_hits,
        );
        if s > last_offset {
            drain_sliding_group::<R, TWO_HIT>(
                params,
                rows,
                subject,
                table,
                diag_mask,
                diag_offset,
                buckets,
                cursor,
                &mut prev_bucket,
                init_hsps,
            );
        } else {
            finish_multipass_group::<R, TWO_HIT>(
                params,
                rows,
                subject,
                table,
                diag_mask,
                diag_offset,
                buckets,
                group,
                cursor,
                s,
                word_pos,
                last_offset,
                max_hits,
                &mut prev_bucket,
                init_hsps,
            );
        }
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn drain_sliding_group<R: PssmRows, const TWO_HIT: bool>(
    params: &RpsWordParams<'_>,
    rows: &R,
    subject: &[u8],
    table: &mut [u32],
    diag_mask: i32,
    diag_offset: i32,
    buckets: &mut RpsBuckets,
    mut cursor: RpsReadCursor,
    prev_bucket: &mut usize,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    while let Some(hits) = buckets.next_block(&mut cursor) {
        let bucket = (hits.q_base >> RPS_BUCKET_BITS) as usize;
        if bucket != *prev_bucket {
            advance_bucket_cold(table, diag_mask as usize, prev_bucket, bucket);
        }
        if TWO_HIT {
            two_hit_block(
                params,
                rows,
                subject,
                table,
                diag_mask,
                diag_offset,
                hits,
                init_hsps,
            );
        } else {
            one_hit_block::<R, true>(
                params,
                rows,
                subject,
                table,
                diag_mask,
                diag_offset,
                hits,
                init_hsps,
            );
        }
    }
}

/// Handles a single group whose hits alone exceed `MAX_SCAN_HITS`: runs all `s`-passes of that
/// group on a group-wide ring `(group_span + subject.len()).next_power_of_two()` (1 MiB for
/// typical proteins), carrying the `subject.len() - 1` overlap diagonals across the group
/// boundaries.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn finish_multipass_group<R: PssmRows, const TWO_HIT: bool>(
    params: &RpsWordParams<'_>,
    rows: &R,
    subject: &[u8],
    table: &mut [u32],
    diag_mask: i32,
    diag_offset: i32,
    buckets: &mut RpsBuckets,
    group: usize,
    mut cursor: RpsReadCursor,
    mut s: i32,
    mut word_pos: usize,
    last_offset: i32,
    max_hits: i32,
    prev_bucket: &mut usize,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let s_len = subject.len();
    let (g0, g1) = buckets.group_bucket_range(group);
    let g0_q = (g0 << RPS_BUCKET_BITS) as i32;
    let g1_q = (g1 << RPS_BUCKET_BITS) as i32;
    let group_span = (g1 - g0) << RPS_BUCKET_BITS;
    let g_len = group_span.saturating_add(s_len).next_power_of_two();
    let g_mask = (g_len - 1) as i32;
    let mut g_table = vec![0u32; g_len];
    if *prev_bucket != usize::MAX {
        let prev_b = *prev_bucket;
        let carry_start = g0_q - (s_len as i32 - 1);
        let carry_end = ((prev_b + 1) << RPS_BUCKET_BITS) as i32;
        for d in carry_start..carry_end {
            g_table[(d & g_mask) as usize] = table[(d & diag_mask) as usize];
        }
    }
    loop {
        while let Some(hits) = buckets.next_block(&mut cursor) {
            if TWO_HIT {
                two_hit_block(
                    params,
                    rows,
                    subject,
                    &mut g_table,
                    g_mask,
                    diag_offset,
                    hits,
                    init_hsps,
                );
            } else {
                one_hit_block::<R, true>(
                    params,
                    rows,
                    subject,
                    &mut g_table,
                    g_mask,
                    diag_offset,
                    hits,
                    init_hsps,
                );
            }
        }
        if s > last_offset {
            break;
        }
        if group == 0 {
            scan_group0(&params.table, subject, &mut s, buckets, max_hits);
            cursor = buckets.start_group(&params.table, subject, 0);
        } else {
            cursor = buckets.fill_group_bounded(
                &params.table,
                subject,
                group,
                &mut s,
                &mut word_pos,
                max_hits,
            );
        }
    }
    table.fill(0);
    let tail_start = g1_q - (s_len as i32 - 1);
    for d in tail_start..g1_q {
        table[(d & diag_mask) as usize] = g_table[(d & g_mask) as usize];
    }
    *prev_bucket = g1 - 1;
}

#[inline(never)]
fn two_hit_sliding<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    ring_len: usize,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    run_sliding::<R, true>(
        params,
        rows,
        scratch,
        subject,
        ring_len,
        MAX_SCAN_HITS,
        init_hsps,
    );
}

#[inline(never)]
fn one_hit_sliding<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    ring_len: usize,
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    run_sliding::<R, false>(
        params,
        rows,
        scratch,
        subject,
        ring_len,
        MAX_SCAN_HITS,
        init_hsps,
    );
}

/// `s_BlastRPSWordFinder_TwoHit`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:331-426
#[cold]
#[inline(never)]
fn two_hit<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let RpsWordScratch { diag, buckets } = scratch;
    let diag_offset = diag.offset;
    let diag_mask = diag.mask;
    // The table's entries: `mask + 1` of them (`RpsDiagTable::reset_for`), so that a diagonal
    // masked with `mask` is inside.
    let table = &mut diag.array[..=diag_mask as usize];
    let last_offset = subject.len() as i32 - wordsize;
    let mut first_offset = 0i32;
    while first_offset <= last_offset {
        scan_subject(
            &params.table,
            subject,
            &mut first_offset,
            buckets,
            MAX_SCAN_HITS,
        );
        // The hits in passes of groups of buckets (`RpsBuckets::start_group`), each read as NCBI reads its
        // buckets: the order of the hits is NCBI's whatever the grouping.
        for group in 0..buckets.groups() {
            if !buckets.needs_pass(group) {
                continue;
            }
            let mut cursor = buckets.start_group(&params.table, subject, group);
            while let Some(hits) = buckets.next_block(&mut cursor) {
                two_hit_block(
                    params,
                    rows,
                    subject,
                    table,
                    diag_mask,
                    diag_offset,
                    hits,
                    init_hsps,
                );
            }
        }
    }
}

/// `s_BlastRPSWordFinder_OneHit`: the diagonal is `(s - q) & mask`, the opposite sign of the
/// two-hit finder's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:657-703
#[cold]
#[inline(never)]
fn one_hit<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let RpsWordScratch { diag, buckets } = scratch;
    let diag_offset = diag.offset;
    let diag_mask = diag.mask;
    // (`two_hit`.)
    let table = &mut diag.array[..=diag_mask as usize];
    let last_offset = subject.len() as i32 - wordsize;
    let mut first_offset = 0i32;
    while first_offset <= last_offset {
        scan_subject(
            &params.table,
            subject,
            &mut first_offset,
            buckets,
            MAX_SCAN_HITS,
        );
        // The hits in passes of groups of buckets (`RpsBuckets::start_group`), each read as NCBI reads its
        // buckets: the order of the hits is NCBI's whatever the grouping.
        for group in 0..buckets.groups() {
            if !buckets.needs_pass(group) {
                continue;
            }
            let mut cursor = buckets.start_group(&params.table, subject, group);
            while let Some(hits) = buckets.next_block(&mut cursor) {
                one_hit_block::<R, false>(
                    params,
                    rows,
                    subject,
                    table,
                    diag_mask,
                    diag_offset,
                    hits,
                    init_hsps,
                );
            }
        }
    }
}

#[cfg(test)]
#[path = "wordfinder_tests.rs"]
mod tests;
