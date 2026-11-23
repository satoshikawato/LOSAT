//! One profile as NCBI's RPS traceback sees it (`s_RPSComputeTraceback`, where the profile is
//! the traceback's "query"): the one owner of that view (series PR, item PS10e). The preliminary
//! stage reads a profile's `.rps` rows and offsets; the traceback reads its `.psq` letters and
//! takes its `.pin` length `L_pin` as the profile's length: the rows it reads from the profile's
//! offset (`L_pin` of them), the composition, the words and the copy of the near-identical test,
//! the redo's matrix and query, the e-value length and the reports' letters. A profile whose
//! `L_pin` is the `.rps` column count `L_rps` (every profile `makeprofiledb` writes from a protein
//! `seq-data` of the profile's length) reads nothing of another; one whose `.psq` letters are
//! nucleotide-coded or more or fewer than its columns makes NCBI's traceback read past its data
//! at one of five reads (`ProfileRead`), asked where NCBI makes each:
//!
//! - A, the rows from the profile's offset (`.freq` with composition-based statistics, `.rps` for
//!   rpsblast's rescale without them): NCBI's row-pointer arrays have `R + 1` entries (`R` the
//!   volume's rows; entry `R` the row after the rows, in the file's mapping), so the profile's
//!   offset plus `L_pin` past `R + 1`, or a column of row `R` past the file's last page, faults
//!   (`rows`, `freq_rows`).
//! - B, the composition of the letters (composition-based statistics): a letter of 28 or more
//!   indexes past NCBI's 28-entry tables (`composition`).
//! - C, the near-identical walk on NCBI's copy of the letters (two bytes more than `L_pin`) and
//!   its words (seven fewer than `L_pin`; `redo_fault`, the walk's reads are
//!   `RedoRun::with_profile_reads`').
//! - D, the start of an alignment at or past the `L_pin` rows of NCBI's heap matrix (the redo;
//!   rpsblast's rescaled matrix without composition-based statistics, `trace_reads`: an
//!   alignment that starts below `L_pin` reads only rows below `L_pin`, so the traceback reads
//!   that matrix as it is).
//! - E, the identities and positives of a traced alignment without composition-based statistics:
//!   a profile letter of 28 or more in a substituted column reads past NCBI's 28 score-matrix
//!   rows (`trace_reads`).
//!
//! rpstblastn without composition-based statistics reads the mapped `.rps` rows from the offset,
//! not a heap matrix: there an alignment that starts at or past `L_pin` stays in NCBI's mapping
//! but gives HSPs whose coordinates, scores and identities disagree (PD-LOSAT-NCBI-DEFECTS, the
//! rule-2 item pending the maintainer's decision). LOSAT refuses it for now, as its own arm
//! (`ProfileRead::StartInMapping`), the one place to change.
//!
//! The first read that leaves NCBI's data refuses the query batch (`TracebackProfile::check`, a
//! traceback error as `RpsPastRows::overran`'s): NCBI ends there with a fault or reads what it does
//! not own. A consistent profile never faults: its preliminary HSPs stay in its `L_pin` columns.

use std::borrow::Cow;
use std::sync::OnceLock;

use anyhow::Result;

use crate::algorithm::pssm_search::redo1::RedoFault;
use crate::algorithm::pssm_search::rows::PssmRows;
use crate::algorithm::pssm_search::traceback0::TraceReads;
use crate::common::GapEditOp;
use crate::core::blast_util::BLASTAA_SIZE;
use crate::seqdb::rps::RpsVolumeFiles;

/// The row a view gives for a read that faults (the caller refuses before it uses the result).
static ZERO_ROW: [i32; BLASTAA_SIZE] = [0; BLASTAA_SIZE];

/// Where NCBI's traceback of a profile leaves its data (module documentation: reads A to E).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ProfileRead {
    /// A: a row past NCBI's row-pointer array, or row `R`'s bytes past the file's last page.
    Rows,
    /// B: this letter (28 or more) of the composition.
    Composition(u8),
    /// C: the near-identical walk past its copy of the letters or its words.
    Copy,
    /// D: an alignment that starts at or past the `L_pin` rows of NCBI's heap matrix.
    Start,
    /// D of rpstblastn without composition-based statistics: the start at or past `L_pin` in
    /// NCBI's mapped rows (rule 2, pending the maintainer's decision; module documentation).
    StartInMapping,
    /// E: this profile letter (28 or more) in a substituted column of a traced alignment.
    Substitution(u8),
}

impl ProfileRead {
    /// The read as the refusal names it.
    fn read(self) -> String {
        match self {
            Self::Rows => "reads past its rows".to_string(),
            Self::Composition(letter) => {
                format!("reads letter {letter} past its composition tables")
            }
            Self::Copy => "reads past its copy of the consensus".to_string(),
            Self::Start | Self::StartInMapping => {
                "starts an alignment past the consensus".to_string()
            }
            Self::Substitution(letter) => format!("reads letter {letter} past its score matrix"),
        }
    }
}

/// One profile of a volume in the traceback of one HSP list (module documentation).
pub(crate) struct TracebackProfile<'a> {
    /// The profile's OID in its volume.
    oid: u32,
    /// The profile's offset in the volume's rows.
    start: usize,
    /// `L_rps`, its `.rps` columns (its rows without the sentinel row).
    columns: usize,
    /// `R`, the volume's rows (its last offset).
    num_rows: usize,
    /// Its `.psq` letters (`L_pin` of them).
    letters: &'a [u8],
    /// The first letter of 28 or more.
    first_wide: Option<u8>,
    /// The program of LOSAT's refusals (`RPSBLAST`, `RPSTBLASTN`).
    program: &'static str,
    /// The first read that left NCBI's data (shared by the row views the alignments read).
    fault: OnceLock<ProfileRead>,
}

impl<'a> TracebackProfile<'a> {
    /// The profile `oid` of a volume of offsets `offsets` (`num_profiles + 1`) whose `.psq`
    /// letters are `letters`.
    pub(crate) fn new(oid: u32, offsets: &[i32], letters: &'a [u8], program: &'static str) -> Self {
        let oid_index = oid as usize;
        Self {
            oid,
            start: offsets[oid_index] as usize,
            columns: (offsets[oid_index + 1] - offsets[oid_index] - 1) as usize,
            num_rows: offsets[offsets.len() - 1] as usize,
            letters,
            first_wide: first_wide(letters),
            program,
            fault: OnceLock::new(),
        }
    }

    /// `L_pin`: the traceback's length of the profile (its alignments' bound, the end of the
    /// tree's root one past it, the redo's query length, the e-value length).
    #[inline]
    pub(crate) fn length(&self) -> usize {
        self.letters.len()
    }

    /// The `L_pin` letters (the traceback's and the redo's query, the reports' identities).
    #[inline]
    pub(crate) fn letters(&self) -> &'a [u8] {
        self.letters
    }

    fn record(&self, read: ProfileRead) {
        // The first read stays (a later `set` fails).
        let _ = self.fault.set(read);
    }

    /// The first read that left NCBI's data.
    #[inline]
    pub(crate) fn fault(&self) -> Option<ProfileRead> {
        self.fault.get().copied()
    }

    /// LOSAT's refusal of the read `read`.
    pub(crate) fn refusal(&self, read: ProfileRead) -> String {
        format!(
            "a hit on profile {} of an RPS-BLAST volume, whose consensus has {} letters in .psq for {} columns of scores in .rps (makeprofiledb writes it for a scored profile whose seq-data is nucleotide-coded or has more or fewer letters than its length), where NCBI BLAST+ 2.17.0's traceback {}, is not supported by LOSAT's {}",
            self.oid,
            self.length(),
            self.columns,
            read.read(),
            self.program
        )
    }

    /// The refusal of the first read that left NCBI's data, if any.
    pub(crate) fn check(&self) -> Result<()> {
        match self.fault() {
            Some(read) => Err(anyhow::anyhow!(self.refusal(read))),
            None => Ok(()),
        }
    }

    /// The refusal of the first read that left NCBI's data, after recording `read`.
    fn refuse(&self, read: ProfileRead) -> Result<()> {
        self.record(read);
        Err(anyhow::anyhow!(self.refusal(self.fault().unwrap_or(read))))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:968-982
    /// Whether NCBI's `R + 1` row pointers (of a volume of `R` rows) hold the profile's `L_pin`
    /// from its offset (read A).
    fn rows_in_pointers(&self) -> bool {
        self.start + self.length() <= self.num_rows + 1
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4719-4730
    /// The `L_pin` `.rps` rows from the profile's offset that rpsblast's rescale reads without
    /// composition-based statistics (`RPSRescalePssm`, read A): the volume's `rows` (`R` of
    /// them), then NCBI's row after the rows (`row_after`, `RpsVolumeFiles::row_after_rows`) as
    /// row `R`, each of `alphabet` columns. For `L_pin` up to one more than `L_rps` they are the
    /// profile's own rows (with its sentinel row); longer, the next profile's rows.
    pub(crate) fn rows<'r>(
        &self,
        rows: &'r [[i32; BLASTAA_SIZE]],
        row_after: &[Option<i32>; BLASTAA_SIZE],
        alphabet: usize,
    ) -> Result<Cow<'r, [[i32; BLASTAA_SIZE]]>> {
        let end = self.start + self.length();
        if end <= rows.len() {
            return Ok(Cow::Borrowed(&rows[self.start..end]));
        }
        match row_after_letters(row_after, alphabet) {
            Some(after) if self.rows_in_pointers() => {
                let mut owned = rows[self.start..].to_vec();
                owned.push(after);
                Ok(Cow::Owned(owned))
            }
            _ => self
                .refuse(ProfileRead::Rows)
                .map(|()| Cow::Borrowed(&[][..])),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:985-1006
    /// The `L_pin` frequency-ratio rows from the profile's offset that the redo reads
    /// (`s_RPSFillFreqRatiosInPsiMatrix`, read A) into `out` (`RpsVolumeFiles::freq_ratios`, row
    /// `R` from the file's last page).
    pub(crate) fn freq_rows(
        &self,
        files: &RpsVolumeFiles,
        bytes: &mut Vec<u8>,
        out: &mut Vec<[f64; BLASTAA_SIZE]>,
    ) -> Result<()> {
        if !self.rows_in_pointers() {
            return self.refuse(ProfileRead::Rows);
        }
        let mapped = files
            .freq_ratios(self.start, self.length(), bytes, out)
            .map_err(|error| {
                anyhow::anyhow!(
                    "a profile database whose .freq file cannot be read ({error}) is not supported by LOSAT's {}",
                    self.program
                )
            })?;
        if !mapped {
            return self.refuse(ProfileRead::Rows);
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2309-2340
    // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/composition_adjustment.c:1162-1199
    /// The composition of the letters that the redo reads first (`s_GetQueryInfo`,
    /// `Blast_ReadAaComposition`, read B): every letter below 28.
    pub(crate) fn composition(&self) -> Result<()> {
        match self.first_wide {
            Some(letter) => self.refuse(ProfileRead::Composition(letter)),
            None => Ok(()),
        }
    }

    /// The redo's reads of the profile (`RedoRun::profile_fault`: reads C and D).
    pub(crate) fn redo_fault(&self, fault: Option<RedoFault>) -> Result<()> {
        match fault {
            Some(RedoFault::Copy) => self.refuse(ProfileRead::Copy),
            Some(RedoFault::Start) => self.refuse(ProfileRead::Start),
            None => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1194-1198
    /// rpstblastn's rows without composition-based statistics: NCBI's `.rps` row pointers from
    /// the profile's offset, `L_pin` the traceback's length; a row past the volume's `rows` is
    /// NCBI's row after the rows (`row_after`, of `alphabet` columns), and one past it read A.
    pub(crate) fn mapped_rows<'m>(
        &'m self,
        rows: &'m [[i32; BLASTAA_SIZE]],
        row_after: &'m [Option<i32>; BLASTAA_SIZE],
        alphabet: usize,
    ) -> MappedRows<'m> {
        MappedRows {
            rows: rows.get(self.start..).unwrap_or(&[]),
            row_after,
            alphabet,
            after: OnceLock::new(),
            profile: self,
        }
    }

    /// What the traceback without composition-based statistics asks of the profile
    /// (`TraceReads`): a start at or past `L_pin` is `past_start` (D, or rpstblastn's
    /// `StartInMapping`), a traced alignment's substituted columns E.
    pub(crate) fn trace_reads(&self, past_start: ProfileRead) -> ProfileTraceReads<'_> {
        ProfileTraceReads {
            profile: self,
            past_start,
        }
    }
}

/// The first of `letters` of 28 or more: their largest letter first (a loop that compiles to
/// vector maxima, once per HSP list), the first such letter only when there is one.
fn first_wide(letters: &[u8]) -> Option<u8> {
    let largest = letters
        .iter()
        .fold(0u8, |largest, &letter| largest.max(letter));
    if usize::from(largest) < BLASTAA_SIZE {
        return None;
    }
    letters
        .iter()
        .copied()
        .find(|&letter| usize::from(letter) >= BLASTAA_SIZE)
}

/// NCBI's row after the rows as `alphabet` columns are read (`None` where one lies past the file's
/// last page); the other columns 0.
fn row_after_letters(
    row_after: &[Option<i32>; BLASTAA_SIZE],
    alphabet: usize,
) -> Option<[i32; BLASTAA_SIZE]> {
    if row_after[..alphabet].iter().any(Option::is_none) {
        return None;
    }
    Some(std::array::from_fn(|letter| row_after[letter].unwrap_or(0)))
}

/// `TracebackProfile::mapped_rows`.
pub(crate) struct MappedRows<'m> {
    /// The volume's rows from the profile's offset.
    rows: &'m [[i32; BLASTAA_SIZE]],
    row_after: &'m [Option<i32>; BLASTAA_SIZE],
    alphabet: usize,
    /// NCBI's row after the rows as `alphabet` columns are read, made when first read.
    after: OnceLock<Option<[i32; BLASTAA_SIZE]>>,
    profile: &'m TracebackProfile<'m>,
}

impl MappedRows<'_> {
    /// Row `pos` from the offset at or past the volume's rows: NCBI's row after the rows, or
    /// read A. Out of the gapped alignment's loop, which reads a profile's own rows.
    #[cold]
    #[inline(never)]
    fn past_rows(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        let after = self
            .after
            .get_or_init(|| row_after_letters(self.row_after, self.alphabet));
        match after {
            Some(row) if pos == self.rows.len() => row,
            _ => {
                self.profile.record(ProfileRead::Rows);
                &ZERO_ROW
            }
        }
    }
}

impl PssmRows for MappedRows<'_> {
    #[inline(always)]
    fn row(&self, pos: usize) -> &[i32; BLASTAA_SIZE] {
        match self.rows.get(pos) {
            Some(row) => row,
            None => self.past_rows(pos),
        }
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.profile.length()
    }
}

/// `TracebackProfile::trace_reads`.
pub(crate) struct ProfileTraceReads<'p> {
    profile: &'p TracebackProfile<'p>,
    past_start: ProfileRead,
}

impl TraceReads for ProfileTraceReads<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512
    fn starts(&self, q_start: usize) -> bool {
        if self.profile.fault.get().is_some() {
            return false;
        }
        if q_start >= self.profile.length() {
            self.profile.record(self.past_start);
            return false;
        }
        true
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:586-594
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:796-812
    // Each substituted column whose letters differ reads the matrix row of the profile's letter;
    // the user query's letters are below 28, so a profile letter of 28 or more always differs.
    fn traced(&self, query_start: i32, edit_script: &[GapEditOp]) -> bool {
        if self.profile.fault.get().is_some() {
            return false;
        }
        if self.profile.first_wide.is_none() {
            return true;
        }
        let mut q = query_start.max(0) as usize;
        for op in edit_script {
            match *op {
                GapEditOp::Sub(n) => {
                    let end = q + n as usize;
                    if let Some(&letter) = self
                        .profile
                        .letters
                        .get(q..end.min(self.profile.length()))
                        .and_then(|letters| {
                            letters
                                .iter()
                                .find(|&&letter| usize::from(letter) >= BLASTAA_SIZE)
                        })
                    {
                        self.profile.record(ProfileRead::Substitution(letter));
                        return false;
                    }
                    q = end;
                }
                GapEditOp::Ins(n) => q += n as usize,
                GapEditOp::Del(_) => {}
            }
        }
        true
    }
}

#[cfg(test)]
#[path = "profile_tests.rs"]
mod tests;
