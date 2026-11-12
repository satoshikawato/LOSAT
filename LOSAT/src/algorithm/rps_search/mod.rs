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

/// The concatenated profile rows of a volume as the preliminary stage reads them: one row per
/// consensus position and one sentinel row after each profile (`lookup->rps_pssm`, the query of
/// the ungapped extensions, whose `query->length` is the concatenation's).
#[derive(Clone, Copy)]
pub(crate) struct RpsRows<'a> {
    pub(crate) rows: &'a [[i32; BLASTAA_SIZE]],
}

impl PssmRows for RpsRows<'_> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        &self.rows[pos]
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.rows.len()
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
