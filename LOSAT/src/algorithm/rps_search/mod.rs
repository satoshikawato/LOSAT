//! The search of RPS-BLAST (`losat rpsblast`, series PR, item PS11a): the user's queries
//! against a database of profiles (`makeprofiledb -dbtype rps`), whose concatenated profiles are
//! the search's "query" (one context per profile) and each user query its "subject". NCBI's
//! `s_RPSPreliminarySearchEngine` scans each user query with the prebuilt lookup table
//! (`scan.rs`), extends the word hits on NCBI's diagonal table (`wordfinder.rs`, with PS4a's PSSM
//! ungapped kernels), aligns them on each profile's rows (`prelim.rs`, with PS4a's score-only
//! gapped kernels), and its RPS traceback rescales each profile to the query's composition
//! (`rescale.rs`) before aligning again and swapping the roles back (`traceback.rs`). The
//! engine (`engine.rs`) runs every (query, volume) pair alone, the queries in parallel, and
//! computes NCBI's per-batch Karlin-Altschul K of the outfmt 0 footers after.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1131-1158
//! ```c
//!     num_db_seqs = BlastSeqSrcGetNumSeqs(seq_src);
//!     dbsize = BlastSeqSrcGetTotLen(seq_src) + num_db_seqs;
//!     if (dbsize > INT4_MAX)
//!       return -3;
//!
//!     /* Concatenate all of the DB sequences together, and pretend
//!       this is a large multiplexed sequence. Note that because the
//!       scoring is position-specific, the actual sequence data is
//!       not needed */
//!
//!     memset(&concat_db, 0, sizeof(concat_db)); /* fill in SequenceBlk */
//!     concat_db.length = (Int4) dbsize;
//!
//!     /* Change the table of diagonals that will be used for the
//!       search; we need a diag table that can fit the entire
//!       concatenated DB */
//!     BlastExtendWordFree(aux_struct->ewp);
//!     BlastExtendWordNew(concat_db.length, word_params, &aux_struct->ewp);
//!
//!     /* Run the search; the input query is what gets scanned
//!       and the concatenated DB is the sequence associated with
//!       the score matrix. This essentially means that 'query'
//!       and 'subject' have opposite conventions for the search.
//! ```

pub(crate) mod engine;
pub(crate) mod prelim;
pub(crate) mod redo;
pub(crate) mod rescale;
pub(crate) mod scan;
pub(crate) mod traceback;
pub(crate) mod volumes;
pub(crate) mod wordfinder;

use crate::algorithm::pssm_search::rows::PssmRows;
use crate::core::blast_util::BLASTAA_SIZE;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:948, 1066-1069 and
// c++/src/algo/blast/core/blast_aalookup.c:184-190
// ```c
//     n = MIN(subject->length - s_off, query_size - q_off);
// ...
//         total_score = s_BlastPSSMExtendRight(matrix, subject, query->length,
//                                              s_right_off + 1, q_right_off + 1,
// ...
//     num_pssm_rows = lookup->rps_seq_offsets[lookup->num_profiles];
//     lookup->rps_pssm = (Int4 **) malloc((num_pssm_rows + 1) * sizeof(Int4 *));
//     pssm_start = profile_header->start_offsets + lookup->num_profiles + 1;
//
//     for (i = 0; i < num_pssm_rows + 1; i++) {
//         lookup->rps_pssm[i] = pssm_start;
// ```
/// The concatenated profile rows of a volume as the ungapped extensions read them: one row per
/// consensus position and one sentinel row after each profile (`lookup->rps_pssm`, the query of
/// the extensions), and their `query->length`, the unit's concatenated length `Lc`
/// (`RpsRun::concat_length`), which bounds the extensions to the right.
///
/// `length` is `Lc` within 0 and the rows: an `Lc` below the rows stops the extensions there (and
/// at 0 or below none runs to the right of a position, as NCBI's negative bound), one above them
/// reads NCBI's rows up to the last profile's sentinel row, which ends an extension unless the
/// X-drop is above its scores (`makeprofiledb` writes `-BLAST_SCORE_MAX`). An extension that
/// scores that row and goes on reads in NCBI the row after the rows (`rps_pssm[rows]`, the bytes
/// after the last row in the `.rps` file's mapping), and then past the pointer array: `past`
/// (where `Lc` is above the rows) holds that row as NCBI reads it, which the extension scores
/// (`PssmRows::past_end`; DECISIONS D1a4-9 of line RG), and the flag of LOSAT's refusal of the
/// query batch (`prelim::RpsPrelim::search_frame`) where NCBI's read faults: a letter whose bytes
/// lie past the file's last page, or a step after that row (`Lc` of the rows plus 2 or more).
#[derive(Clone, Copy)]
pub(crate) struct RpsRows<'a> {
    pub(crate) rows: &'a [[i32; BLASTAA_SIZE]],
    pub(crate) length: usize,
    pub(crate) past: Option<&'a RpsPastRows<'a>>,
}

/// What NCBI's extensions read past the rows where `Lc` is above them (`RpsRows::past`), one for
/// each context of a volume's search.
pub(crate) struct RpsPastRows<'a> {
    /// NCBI's row after the rows for each letter (`seqdb::rps::RpsVolumeFiles::row_after_rows`;
    /// `None` where its read faults).
    pub(crate) row: &'a [Option<i32>; BLASTAA_SIZE],
    /// `Lc`, NCBI's bound of the extensions to the right.
    pub(crate) concat_length: i32,
    /// Set where an extension of the context reads what NCBI cannot (LOSAT's refusal of the
    /// batch).
    pub(crate) overran: std::sync::atomic::AtomicBool,
}

impl<'a> RpsRows<'a> {
    /// The rows as `RpsRowsWhole` where `Lc` is their count (`length` is the rows' count and
    /// `past` is `None`), else `None`.
    #[inline(always)]
    pub(crate) fn whole(&self) -> Option<RpsRowsWhole<'a>> {
        (self.past.is_none() && self.length == self.rows.len())
            .then_some(RpsRowsWhole { rows: self.rows })
    }
}

/// `RpsRows` whose `Lc` is the rows' count (`RpsRows::whole`), as a type: the rows alone bound
/// the extensions (`len` is the slice's length, no second value) and there is no step past them
/// (`PssmRows::past_end` is `None`), so the word finders' loops instantiated with it are those
/// of a search without `Lc`: no code of the step, and one value fewer to keep in the loops.
/// `wordfinder::word_finder` scans with it in every search whose `Lc` is the rows' count (every
/// database but one whose alias file's `LENGTH` and `NSEQ` give another `Lc`), and with `RpsRows`
/// (the bound and the step) out of line otherwise (`wordfinder::word_finder_other_length`): one
/// generic scan, the same function wherever both are valid.
#[derive(Clone, Copy)]
pub(crate) struct RpsRowsWhole<'a> {
    rows: &'a [[i32; BLASTAA_SIZE]],
}

impl PssmRows for RpsRowsWhole<'_> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        &self.rows[pos]
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.rows.len()
    }
}

impl PssmRows for RpsRows<'_> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        &self.rows[pos]
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.length
    }

    #[inline(always)]
    fn past_end(&self, end: i32, letter: Option<u8>) -> Option<i32> {
        let past = self.past?;
        let letter = letter?;
        if end as usize != self.length {
            return None;
        }
        let score = past.row[usize::from(letter)];
        if score.is_none() {
            past.overran
                .store(true, std::sync::atomic::Ordering::Relaxed);
        }
        score
    }

    #[inline(always)]
    fn past_end_continues(&self, end: i32, subject_left: bool) {
        if let Some(past) = self.past {
            if subject_left && end < past.concat_length {
                past.overran
                    .store(true, std::sync::atomic::Ordering::Relaxed);
            }
        }
    }
}

/// The rows of one profile (one context of the concatenation): the PSSM pointer that
/// `BLAST_GetGappedScore` moves to the profile's start, and the context's length.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3891-3894
/// ```c
///       if (rpsblast_pssms)
///          gap_align->sbp->psi_matrix->pssm->data = rpsblast_pssms +
///              query_info->contexts[context].query_offset;
/// ```
#[derive(Clone, Copy)]
pub(crate) struct ProfileSlice<'a> {
    pub(crate) rows: &'a [[i32; BLASTAA_SIZE]],
    pub(crate) start: usize,
    pub(crate) len: usize,
}

impl PssmRows for ProfileSlice<'_> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        &self.rows[self.start + pos]
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.len
    }
}

/// `RPS_K_MULT`: the factor of the profile's `.aux` K in the gapped Karlin block of its
/// traceback.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1011
/// ```c
/// #define RPS_K_MULT 1.2
/// ```
pub(crate) const RPS_K_MULT: f64 = 1.2;

/// `BLAST_Nint`: the nearest integer, halves away from zero.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
/// ```c
/// long BLAST_Nint(double x)
/// {
///    x += (x >= 0. ? 0.5 : -0.5);
///    return (long)x;
/// }
/// ```
#[inline]
pub(crate) fn blast_nint(x: f64) -> i64 {
    let x = x + if x >= 0.0 { 0.5 } else { -0.5 };
    x as i64
}

#[cfg(test)]
mod tests {
    use super::{RpsPastRows, RpsRows};
    use crate::algorithm::pssm_search::rows::PssmRows;
    use crate::algorithm::pssm_search::ungapped::pssm_extend_one_hit;
    use crate::core::blast_util::BLASTAA_SIZE;
    use std::sync::atomic::{AtomicBool, Ordering};

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:938-975, 1036-1083 (quoted
    // at `pssm_search::ungapped`): the extension to the right reads up to `query->length`, the
    // concatenated length `Lc`, which `RpsRows::length` holds within the rows, and past them the
    // row after the rows (`RpsRows::past`). Nine rows that score letter 1 (A) 5,000 and one
    // sentinel row of -32,767 (`-BLAST_SCORE_MAX`); a word of three letters at 0 against a query
    // of A.
    #[test]
    fn the_extension_to_the_right_stops_at_the_concatenated_length() {
        let mut row = [-4i32; BLASTAA_SIZE];
        row[1] = 5_000;
        let mut rows = vec![row; 9];
        rows.push([-32_767; BLASTAA_SIZE]);
        let query = [1u8; 20];
        let zeros = [Some(0i32); BLASTAA_SIZE];
        let faults = [None; BLASTAA_SIZE];
        let extend = |length: usize,
                      past: Option<(&[Option<i32>; BLASTAA_SIZE], i32)>,
                      subject: &[u8],
                      dropoff: i32| {
            let past = past.map(|(row, concat_length)| RpsPastRows {
                row,
                concat_length,
                overran: AtomicBool::new(false),
            });
            let rows = RpsRows {
                rows: &rows,
                length,
                past: past.as_ref(),
            };
            let hit = pssm_extend_one_hit(&rows, subject, 0, 0, dropoff, 3);
            let flagged = past.is_some_and(|past| past.overran.load(Ordering::Relaxed));
            (hit, flagged)
        };
        // `Lc` 6 below the 10 rows: the walk from 3 stops at 6.
        let (short, _) = extend(6, None, &query, 40_000);
        assert_eq!(
            (short.ungapped_data.length, short.ungapped_data.score),
            (6, 30_000)
        );
        assert_eq!(short.s_last_off, 6);
        // `Lc` 11: the sentinel row stops an X-drop of 20.
        let (stopped, flagged) = extend(10, Some((&zeros, 11)), &query, 20);
        assert_eq!(
            (stopped.ungapped_data.length, stopped.ungapped_data.score),
            (9, 45_000)
        );
        assert_eq!(stopped.s_last_off, 9);
        assert!(!flagged);
        // An X-drop above the sentinel's scores: the walk scores the sentinel row, then NCBI's row
        // after the rows (zeros of the `.rps` file's last page) and ends at `Lc` 11.
        let (through, flagged) = extend(10, Some((&zeros, 11)), &query, 40_000);
        assert_eq!(through.s_last_off, 11);
        assert_eq!(through.ungapped_data.length, 9);
        assert!(!flagged);
        // `Lc` 10 (the rows): NCBI's bound is the rows' end.
        let (at_rows, _) = extend(10, None, &query, 40_000);
        assert_eq!(at_rows.s_last_off, 10);
        // The row after the rows past the file's last page: NCBI's read faults.
        assert!(extend(10, Some((&faults, 11)), &query, 40_000).1);
        // `Lc` 12: NCBI reads past the pointer array after that row.
        assert!(extend(10, Some((&zeros, 12)), &query, 40_000).1);
        // The query ends with the rows: NCBI's bound is the query's end too.
        assert!(!extend(10, Some((&faults, 12)), &query[..10], 40_000).1);
        // A row after the rows of other bytes (a `.rps` file with bytes after its rows) scores
        // them: 40,000 for A, above the sentinel's -32,767, takes the best score to that row.
        let mut high = zeros;
        high[1] = Some(40_000);
        let (scored, flagged) = extend(10, Some((&high, 11)), &query, 1_000_000);
        assert_eq!(scored.ungapped_data.length, 11);
        assert!(!flagged);
    }

    // `RpsRows::whole` gives the rows as one slice only where the slice's own length is every
    // bound of the scan: `Lc` is the rows' count and there is no row after the rows. An `Lc`
    // below the rows (an alias file's `LENGTH`) and an `Lc` past them keep `RpsRows`.
    #[test]
    fn the_rows_are_one_slice_only_where_the_concatenated_length_is_their_count() {
        let rows = vec![[0i32; BLASTAA_SIZE]; 10];
        let zeros = [Some(0i32); BLASTAA_SIZE];
        let past = RpsPastRows {
            row: &zeros,
            concat_length: 11,
            overran: AtomicBool::new(false),
        };
        let whole = |length: usize, past: Option<&RpsPastRows<'_>>| {
            let rows = RpsRows {
                rows: &rows,
                length,
                past,
            };
            rows.whole().map(|whole| PssmRows::len(&whole))
        };
        assert_eq!(whole(10, None), Some(10));
        assert_eq!(whole(9, None), None);
        assert_eq!(whole(0, None), None);
        assert_eq!(whole(10, Some(&past)), None);
    }
}
