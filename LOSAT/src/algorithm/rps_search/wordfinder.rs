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
    q_place, s_place, scan_subject, RpsBlockPool, RpsBuckets, RpsHits, RpsTable, MAX_SCAN_HITS,
};
use super::RpsRows;
use crate::algorithm::blastp::extension::BlastpUngappedData;
use crate::algorithm::pssm_search::rows::PssmRows;
use crate::algorithm::pssm_search::ungapped::{pssm_extend_one_hit, pssm_extend_two_hit};
use crate::seqdb::rps::RPS_WORDSIZE;

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
}

impl RpsDiagTable {
    /// An empty table (sized by `reset_for`).
    pub(crate) fn new() -> Self {
        Self {
            array: Vec::new(),
            mask: 0,
            offset: 0,
            window: 0,
        }
    }

    /// `BlastExtendWordNew(concat_length, ...)`: a zeroed table of the unit's length. The
    /// table is allocated zeroed and not written (NCBI's `calloc`): only the pages of the
    /// diagonals that the queries hit become resident, whatever the window.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-62
    /// The table has `length` entries (`diag_table_length`).
    pub(crate) fn reset_for(&mut self, length: i32, window: i32) {
        self.array = vec![0u32; length as usize];
        self.mask = length - 1;
        self.offset = window;
        self.window = window;
    }

    /// The table of a search of `length` entries (`diag_table_length`): the table as it is when
    /// it has that length and the window already (series PR, item PS12s, D3: a worker keeps its
    /// table from volume to volume), else a zeroed one (`reset_for`). The entries an earlier
    /// volume's queries wrote read as "no hit on this diagonal yet", as those of an earlier query
    /// of the same volume do (the offset rule of `exit`, which keeps NCBI's clear at
    /// `INT4_MAX / 4`): a query's hits do not depend on what the table held.
    pub(crate) fn prepare_for(&mut self, length: i32, window: i32) {
        if self.array.len() != length as usize || self.window != window {
            self.reset_for(length, window);
        }
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
    /// `RpsRun::concat_length`: `BlastExtendWordNew(concat_db.length, ...)`) whose queries (frames)
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
        self.diag.prepare_for(
            diag_table_length(concat_length, window, num_rows, longest),
            window,
        );
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
#[inline(always)]
fn scan<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    if params.window > 0 {
        two_hit(params, rows, scratch, subject, init_hsps);
    } else {
        one_hit(params, rows, scratch, subject, init_hsps);
    }
}

/// `s_BlastRPSWordFinder_TwoHit`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:331-426
fn two_hit<R: PssmRows>(
    params: &RpsWordParams<'_>,
    rows: &R,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let window = params.window;
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
            while let Some(RpsHits {
                entries,
                q_base,
                s_base,
            }) = buckets.next_block(&mut cursor)
            {
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
                    let diag_coord =
                        (q.wrapping_sub(s).wrapping_add(diag_base) & diag_mask) as usize;
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
        }
    }
}

/// `s_BlastRPSWordFinder_OneHit`: the diagonal is `(s - q) & mask`, the opposite sign of the
/// two-hit finder's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:657-703
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
            while let Some(RpsHits {
                entries,
                q_base,
                s_base,
            }) = buckets.next_block(&mut cursor)
            {
                // (`two_hit`: `subject_offset - query_offset` is `s - q - diag_base`.)
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
                    let diag_coord =
                        (s.wrapping_sub(q).wrapping_sub(diag_base) & diag_mask) as usize;
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
        }
    }
}

#[cfg(test)]
#[path = "wordfinder_tests.rs"]
mod tests;
