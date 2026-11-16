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
/// scores that row and goes on reads in NCBI past its rows (the pointer after the last row
/// points past them in the `.rps` file's mapping, the next ones past `rps_pssm`): with
/// `overran` (where `Lc` is above the rows) such an extension, one that runs to `length` with
/// letters of the query left, sets the flag, LOSAT's refusal of its query batch
/// (`prelim::RpsPrelim::search_frame`).
#[derive(Clone, Copy)]
pub(crate) struct RpsRows<'a> {
    pub(crate) rows: &'a [[i32; BLASTAA_SIZE]],
    pub(crate) length: usize,
    pub(crate) overran: Option<&'a std::sync::atomic::AtomicBool>,
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
    fn ran_to_end(&self, end: i32, subject_left: bool) {
        if let Some(overran) = self.overran {
            if subject_left && end as usize == self.length {
                overran.store(true, std::sync::atomic::Ordering::Relaxed);
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
    use super::RpsRows;
    use crate::algorithm::pssm_search::ungapped::pssm_extend_one_hit;
    use crate::core::blast_util::BLASTAA_SIZE;
    use std::sync::atomic::{AtomicBool, Ordering};

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:938-975, 1036-1083 (quoted
    // at `pssm_search::ungapped`): the extension to the right reads up to `query->length`, the
    // concatenated length `Lc`, which `RpsRows::length` holds within the rows. Nine rows that score
    // letter 1 (A) 5,000 and one sentinel row of -32,767 (`-BLAST_SCORE_MAX`); a word of three
    // letters at 0 against a query of A.
    #[test]
    fn the_extension_to_the_right_stops_at_the_concatenated_length() {
        let mut row = [-4i32; BLASTAA_SIZE];
        row[1] = 5_000;
        let mut rows = vec![row; 9];
        rows.push([-32_767; BLASTAA_SIZE]);
        let query = [1u8; 20];
        let flag = AtomicBool::new(false);
        let extend = |length: usize, overran: Option<&AtomicBool>, subject: &[u8], dropoff: i32| {
            let rows = RpsRows {
                rows: &rows,
                length,
                overran,
            };
            pssm_extend_one_hit(&rows, subject, 0, 0, dropoff, 3)
        };
        // `Lc` 6 below the 10 rows: the walk from 3 stops at 6.
        let short = extend(6, None, &query, 40_000);
        assert_eq!(
            (short.ungapped_data.length, short.ungapped_data.score),
            (6, 30_000)
        );
        assert_eq!(short.s_last_off, 6);
        // The rows: the sentinel row stops an X-drop of 20.
        let stopped = extend(10, Some(&flag), &query, 20);
        assert_eq!(
            (stopped.ungapped_data.length, stopped.ungapped_data.score),
            (9, 45_000)
        );
        assert_eq!(stopped.s_last_off, 9);
        assert!(!flag.load(Ordering::Relaxed));
        // An X-drop above the sentinel's scores: the walk scores the sentinel row and runs to the
        // rows' end with letters of the query left, where NCBI's bound `Lc` lies past them and
        // NCBI reads past its rows: the flag of an `Lc` above the rows.
        let through = extend(10, Some(&flag), &query, 40_000);
        assert_eq!(through.s_last_off, 10);
        assert!(flag.swap(false, Ordering::Relaxed));
        let _ = extend(10, None, &query, 40_000);
        assert!(!flag.load(Ordering::Relaxed));
        // The query ends with the rows: NCBI's bound is the query's end too.
        let _ = extend(10, Some(&flag), &query[..10], 40_000);
        assert!(!flag.load(Ordering::Relaxed));
    }
}
