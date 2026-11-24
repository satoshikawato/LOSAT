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

pub(crate) mod cdd_input;
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
pub(crate) const RPS_K_MULT: f64 = 1.2;

/// `BLAST_Nint`: the nearest integer, halves away from zero.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
#[inline]
pub(crate) fn blast_nint(x: f64) -> i64 {
    let x = x + if x >= 0.0 { 0.5 } else { -0.5 };
    x as i64
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
