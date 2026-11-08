//! The word finders of RPS-BLAST (`s_BlastRPSWordFinder_TwoHit`, `s_BlastRPSWordFinder_OneHit`):
//! the scan's hits, bucket by bucket, on NCBI's diagonal table of the concatenation, extended
//! with PS4a's PSSM ungapped kernels on the whole concatenation.
//!
//! The diagonal table has `L` entries, `L` the smallest power of two at least the
//! concatenation's length plus the window (`s_BlastDiagTableNew`), and is indexed with `mask =
//! L - 1`: the two-hit finder reads `(q - s) & mask`, the one-hit finder `(s - q) & mask`. Where a
//! query is longer than `L` less the concatenation, two diagonals share an entry, as in NCBI.
//!
//! NCBI keeps one table for the queries of a batch, moving `diag_offset` on by the query's
//! length and the window after each (`Blast_ExtendWordExit`); an entry that an earlier query
//! wrote then reads as "no hit on this diagonal yet" in both finders, as an entry never written
//! does. A worker's table here follows the same rule for the queries it searches, so every query
//! gets NCBI's hits whatever queries came before it.

use super::scan::{scan_subject, RpsBuckets, RpsTable, MAX_SCAN_HITS};
use super::RpsRows;
use crate::algorithm::blastp::extension::BlastpUngappedData;
use crate::algorithm::pssm_search::ungapped::{pssm_extend_one_hit, pssm_extend_two_hit};
use crate::seqdb::rps::RPS_WORDSIZE;

/// One entry of the diagonal table (`DiagStruct`: a 31-bit signed `last_hit` and a 1-bit
/// `flag`), packed in a `u32` as NCBI's bit fields are. The table is a `Vec<u32>` so that a
/// zeroed one comes from the allocator's zeroed pages (`vec![0u32; n]`, `calloc`), as NCBI's.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
/// ```c
/// typedef struct DiagStruct {
///    signed int last_hit   : 31; /**< Offset of the last hit */
///    unsigned int flag      : 1 ; /**< Reset the next extension? */
/// } DiagStruct;
/// ```
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

    /// `BlastExtendWordNew(concat_length, ...)`: a zeroed table of the volume's length. The
    /// table is allocated zeroed and not written (NCBI's `calloc`): only the pages of the
    /// diagonals that the queries hit become resident, whatever the window.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:145-146
    /// ```c
    ///         diag_table->hit_level_array = (DiagStruct *)
    ///             calloc(diag_table->diag_array_length, sizeof(DiagStruct));
    /// ```
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-62
    /// ```c
    ///                 diag_array_length = 1;
    ///                 /* What power of 2 is just longer than the query? */
    ///                 while (diag_array_length < (qlen+window_size))
    ///                 {
    ///                         diag_array_length = diag_array_length << 1;
    ///                 }
    ///                 /* These are used in the word finders to shift and mask
    ///                 rather than dividing and taking the remainder. */
    ///                 diag_table->diag_array_length = diag_array_length;
    ///                 diag_table->diag_mask = diag_array_length-1;
    ///                 diag_table->multiple_hits = multiple_hits;
    ///                 diag_table->offset = window_size;
    ///                 diag_table->window = window_size;
    /// ```
    pub(crate) fn reset_for(&mut self, concat_length: i32, window: i32) {
        let length = diag_array_length(concat_length, window);
        self.array = vec![0u32; length as usize];
        self.mask = length - 1;
        self.offset = window;
        self.window = window;
    }

    /// `Blast_ExtendWordExit`: the offset after a query of `subject_length` letters.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:171-178
    /// ```c
    ///     if (ewp->diag_table) {
    ///         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    ///             ewp->diag_table->offset = ewp->diag_table->window;
    ///             s_BlastDiagClear(ewp->diag_table);
    ///         } else {
    ///             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    ///         }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:72-89 (`s_BlastDiagClear`)
    /// ```c
    ///         for (i=0; i < n; i++) {
    ///             diag_struct_array[i].flag = 0;
    ///             diag_struct_array[i].last_hit = -diag->window;
    ///             if (diag->hit_len_array) diag->hit_len_array[i] = 0;
    ///         }
    /// ```
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

/// The length of the diagonal table of a concatenation of `concat_length` rows.
pub(crate) fn diag_array_length(concat_length: i32, window: i32) -> i32 {
    let mut length = 1i32;
    while length < concat_length + window {
        length <<= 1;
    }
    length
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

/// `BlastRPSWordFinder` on `subject` (the query's letters): the initial HSPs into `init_hsps`
/// (not sorted; the caller sorts them, `Blast_InitHitListSortByScore`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:265-283
/// ```c
///     if (ewp->diag_table->multiple_hits) {
///         status = s_BlastRPSWordFinder_TwoHit(subject, query,
///                                              lut_wrap, ewp,
///                                              matrix,
///                                              cutoffs->cutoff_score,
///                                              cutoffs->x_dropoff,
///                                              init_hitlist, ungapped_stats);
///     } else {
///         status = s_BlastRPSWordFinder_OneHit(subject, query,
///                                              lut_wrap, ewp,
///                                              matrix,
///                                              cutoffs->cutoff_score,
///                                              cutoffs->x_dropoff,
///                                              init_hitlist, ungapped_stats);
///     }
/// ```
pub(crate) fn word_finder(
    params: &RpsWordParams<'_>,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    if params.window > 0 {
        two_hit(params, scratch, subject, init_hsps);
    } else {
        one_hit(params, scratch, subject, init_hsps);
    }
    scratch.diag.exit(subject.len() as i32);
}

/// `s_BlastRPSWordFinder_TwoHit`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:331-426
/// ```c
///     while (first_offset <= last_offset) {
///         /* scan the subject sequence for hits */
///         hits = BlastRPSScanSubject(lookup_wrap, subject, &first_offset);
///
///         totalhits += hits;
///         /* for each region of the concatenated database */
///         for (i = 0; i < lookup->num_buckets; ++i) {
///             ...
///             for (j = 0; j < hits; ++j) {
///                 Uint4 query_offset = offset_pairs[j].qs_offsets.q_off;
///                 Uint4 subject_offset = offset_pairs[j].qs_offsets.s_off;
///
///                 /* calculate the diagonal associated with this query-subject
///                    pair */
///
///                 diag_coord = (query_offset - subject_offset) & diag_mask;
///
///                 /* If the reset bit is set, an extension just happened. */
///                 if (diag_array[diag_coord].flag) {
///                     /* If we've already extended past this hit, skip it. */
///                     if ((Int4) (subject_offset + diag_offset) <
///                         diag_array[diag_coord].last_hit) {
///                         continue;
///                     }
///                     /* Otherwise, start a new hit. */
///                     else {
///                         diag_array[diag_coord].last_hit =
///                             subject_offset + diag_offset;
///                         diag_array[diag_coord].flag = 0;
///                     }
///                 }
///                 /* If the reset bit is cleared, try to start an extension. */
///                 else {
///                     /* find the distance to the last hit on this diagonal */
///                     last_hit = diag_array[diag_coord].last_hit - diag_offset;
///                     diff = subject_offset - last_hit;
///
///                     if (diff >= window) {
///                         /* We are beyond the window for this diagonal; start
///                            a new hit */
///                         diag_array[diag_coord].last_hit =
///                             subject_offset + diag_offset;
///                         continue;
///                     }
///
///                     /* If the difference is less than the wordsize (i.e.
///                        last hit and this hit overlap), give up */
///
///                     if (diff < wordsize) {
///                         continue;
///                     }
///
///                     /* Extend this pair of hits. The extension to the left
///                        must reach the end of the first word in order for
///                        extension to the right to proceed. */
///
///                     score = s_BlastAaExtendTwoHit(matrix, subject, query,
///                                                   last_hit + wordsize,
///                                                   subject_offset, query_offset,
///                                                   dropoff, &hsp_q, &hsp_s,
///                                                   &hsp_len, TRUE,
///                                                   wordsize, &right_extend,
///                                                   &s_last_off);
///
///                     ++hits_extended;
///
///                     /* if the hsp meets the score threshold, report it */
///                     if (score >= cutoff)
///                         BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
///                                          query_offset, subject_offset,
///                                          hsp_len, score);
///
///                     /* If an extension to the right happened, reset the last
///                        hit so that future hits to this diagonal must start
///                        over. */
///
///                     if (right_extend) {
///                         diag_array[diag_coord].flag = 1;
///                         diag_array[diag_coord].last_hit =
///                             s_last_off - (wordsize - 1) + diag_offset;
///                     }
///                     /* Otherwise, make the present hit into the previous hit
///                        for this diagonal */
///                     else {
///                         diag_array[diag_coord].last_hit =
///                             subject_offset + diag_offset;
///                     }
///                 }               /* end else */
///
///             }                   /* end for - done with this batch of hits */
///         }                       /* end for - done with this bucket, call
///                                    scansubject again */
///     }                           /* end while - done with the entire sequence.
///                                  */
/// ```
fn two_hit(
    params: &RpsWordParams<'_>,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let window = params.window;
    let RpsWordScratch { diag, buckets } = scratch;
    let diag_offset = diag.offset;
    let diag_mask = diag.mask;
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
        for bucket in &buckets.buckets {
            for pair in bucket {
                let query_offset = pair.q_off as i32;
                let subject_offset = pair.s_off as i32;
                let diag_coord = (query_offset.wrapping_sub(subject_offset) & diag_mask) as usize;
                let entry = &mut diag.array[diag_coord];
                if entry.flag() {
                    if subject_offset + diag_offset < entry.last_hit() {
                        continue;
                    }
                    entry.set(subject_offset + diag_offset, false);
                } else {
                    let last_hit = entry.last_hit() - diag_offset;
                    let diff = subject_offset - last_hit;
                    if diff >= window {
                        entry.set(subject_offset + diag_offset, entry.flag());
                        continue;
                    }
                    if diff < wordsize {
                        continue;
                    }
                    let hit = pssm_extend_two_hit(
                        &params.rows,
                        subject,
                        last_hit + wordsize,
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
                        entry.set(subject_offset + diag_offset, entry.flag());
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
/// ```c
///     while (first_offset <= last_offset) {
///         /* scan the subject sequence for hits */
///         hits = BlastRPSScanSubject(lookup_wrap, subject, &first_offset);
///
///         totalhits += hits;
///         for (i = 0; i < lookup->num_buckets; i++) {
///             ...
///             for (j = 0; j < hits; ++j) {
///                 Uint4 query_offset = offset_pairs[j].qs_offsets.q_off;
///                 Uint4 subject_offset = offset_pairs[j].qs_offsets.s_off;
///                 diag_coord = (subject_offset - query_offset) & diag_mask;
///                 diff = subject_offset -
///                     (diag_array[diag_coord].last_hit - diag_offset);
///
///                 /* do an extension, but only if we have not already extended
///                    this far */
///                 if (diff >= 0) {
///                     ++hits_extended;
///                     score = s_BlastAaExtendOneHit(matrix, subject, query,
///                                                   subject_offset, query_offset,
///                                                   dropoff, &hsp_q, &hsp_s,
///                                                   &hsp_len, wordsize, TRUE,
///                                                   &s_last_off);
///
///                     /* if the hsp meets the score threshold, report it */
///                     if (score >= cutoff) {
///                         BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
///                                          query_offset, subject_offset,
///                                          hsp_len, score);
///                     }
///                     diag_array[diag_coord].last_hit =
///                         s_last_off - (wordsize - 1) + diag_offset;
///                 }
///             }                   /* end for (one bucket) */
///         }                       /* end for (all buckets) */
///     }                           /* end while */
/// ```
fn one_hit(
    params: &RpsWordParams<'_>,
    scratch: &mut RpsWordScratch,
    subject: &[u8],
    init_hsps: &mut Vec<RpsInitHsp>,
) {
    let wordsize = RPS_WORDSIZE as i32;
    let RpsWordScratch { diag, buckets } = scratch;
    let diag_offset = diag.offset;
    let diag_mask = diag.mask;
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
        for bucket in &buckets.buckets {
            for pair in bucket {
                let query_offset = pair.q_off as i32;
                let subject_offset = pair.s_off as i32;
                let diag_coord = (subject_offset.wrapping_sub(query_offset) & diag_mask) as usize;
                let entry = &mut diag.array[diag_coord];
                let diff = subject_offset - (entry.last_hit() - diag_offset);
                if diff >= 0 {
                    let hit = pssm_extend_one_hit(
                        &params.rows,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_is_the_power_of_two_past_the_concatenation_and_the_window() {
        assert_eq!(diag_array_length(5035, 40), 8192);
        assert_eq!(diag_array_length(1130, 40), 2048);
        assert_eq!(diag_array_length(984, 40), 1024);
        assert_eq!(diag_array_length(985, 40), 2048);
        assert_eq!(diag_array_length(1024, 0), 1024);
    }

    #[test]
    fn the_entry_keeps_ncbis_31_bit_signed_field_and_flag() {
        let mut entry = 0u32;
        entry.set(-40, true);
        assert_eq!((entry.last_hit(), entry.flag()), (-40, true));
        entry.set(i32::MAX / 4 + 10, false);
        assert_eq!((entry.last_hit(), entry.flag()), (i32::MAX / 4 + 10, false));
    }

    // Diagonals `L` apart share an entry (`(q - s) & mask`): across buckets, a hit of a later
    // bucket at a smaller `s` reads the entry that an earlier bucket's hit at a larger `s` wrote
    // (exactness point 1 of the design). A hit of bucket 0 on diagonal -4180 makes the two hits of
    // bucket 1 on diagonal -84 (4,096 apart) read as overlapping it; with a table of 8,192 entries
    // the second of them extends.
    #[test]
    fn diagonals_a_table_length_apart_share_an_entry_across_buckets() {
        use crate::seqdb::rps::RPS_BACKBONE_SIZE;
        let concat = 2100usize;
        let mut rows = vec![[3i32; 28]; concat];
        rows[concat - 1] = [-32767; 28];
        let word = |l: usize| (l << 10) | (l << 5) | l;
        let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
        let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
        // The last letter's rows of the words A A A (q 10), C C C (q 2060), D D D (q 2065).
        for (letter, entry) in [(1usize, 12), (3, 2062), (4, 2067)] {
            cells[word(letter)] = [1, entry, 0, 0];
            pv[word(letter) / 64] |= 1 << (word(letter) % 64);
        }
        let mut query = vec![20u8; 4200];
        query[4190..4193].fill(1);
        query[2144..2147].fill(3);
        query[2149..2152].fill(4);
        let params = RpsWordParams {
            rows: RpsRows { rows: &rows },
            table: RpsTable {
                cells: &cells,
                overflow: &[],
                pv: &pv,
            },
            window: 40,
            cutoff: 1,
            dropoff: 20,
        };
        let hits_with_table = |length: i32| {
            let mut scratch = RpsWordScratch {
                diag: RpsDiagTable::new(),
                buckets: RpsBuckets::new(concat),
            };
            scratch.diag.reset_for(length, 40);
            let mut hits = Vec::new();
            word_finder(&params, &mut scratch, &query, &mut hits);
            (scratch.diag.mask, hits)
        };
        let (mask, ncbi) = hits_with_table(concat as i32);
        assert_eq!(mask, 4095);
        assert!(ncbi.is_empty(), "{ncbi:?}");
        let (wide_mask, wide) = hits_with_table(8000);
        assert_eq!(wide_mask, 8191);
        assert_eq!(wide.len(), 1);
        assert_eq!((wide[0].q_off, wide[0].s_off), (2065, 2149));
    }

    // Stale entries of earlier queries read as no hit: a table that searched another query
    // gives the hits of a fresh one (exactness point 2 of the design).
    #[test]
    fn a_used_table_gives_the_hits_of_a_fresh_one() {
        use crate::seqdb::rps::RPS_BACKBONE_SIZE;
        // One profile of 8 positions that scores +5 on letter 1 (A) and -4 elsewhere, then its
        // sentinel row; every word A A A hits the rows of the profile.
        let mut rows = vec![[-4i32; 28]; 9];
        for row in rows.iter_mut().take(8) {
            row[1] = 5;
        }
        rows[8] = [-32767; 28];
        let mut cells = vec![[0i32; 4]; RPS_BACKBONE_SIZE + 1];
        let mut pv = vec![0u64; RPS_BACKBONE_SIZE / 64];
        let w = (1 << 10) | (1 << 5) | 1;
        let overflow: Vec<i32> = vec![3, 4, 5, 6, 7];
        cells[w] = [6, 2, 0, 0];
        pv[w / 64] |= 1 << (w % 64);
        let params = RpsWordParams {
            rows: RpsRows { rows: &rows },
            table: RpsTable {
                cells: &cells,
                overflow: &overflow,
                pv: &pv,
            },
            window: 40,
            cutoff: 1,
            dropoff: 20,
        };
        let query = [1u8, 1, 1, 1, 1, 1, 20, 1, 1, 1, 1];
        let other = [1u8; 30];
        let run = |before: &[&[u8]]| {
            let mut scratch = RpsWordScratch {
                diag: RpsDiagTable::new(),
                buckets: RpsBuckets::new(rows.len()),
            };
            scratch.diag.reset_for(rows.len() as i32, 40);
            let mut ignored = Vec::new();
            for q in before {
                word_finder(&params, &mut scratch, q, &mut ignored);
            }
            let mut hits = Vec::new();
            word_finder(&params, &mut scratch, &query, &mut hits);
            hits
        };
        let fresh = run(&[]);
        assert!(!fresh.is_empty());
        assert_eq!(run(&[&other]), fresh);
        assert_eq!(run(&[&other, &query, &other]), fresh);
    }
}
