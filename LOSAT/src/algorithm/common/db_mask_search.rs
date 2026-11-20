//! The database masking of a search of a protein database (`-db_soft_mask`, `-db_hard_mask`;
//! series NP, M4b): what BLASTP and BLASTX share, and what TBLASTN and TBLASTX take for their
//! own chunks and reports. The options, their resolution and their refusals are
//! `blastinput::db_mask`; the archive's and the search strategy's program options are
//! `algorithm::common::database::note_archive_database_masking`; the `Masked using:` lines are
//! `algorithm::common::database::pairwise_database_masked`. Here are the parts of the search:
//!
//! - [`SearchMask`]: the resolved algorithm and kind that the search reads;
//! - [`SubjectMasks`]: a reader of the stored masks of database sequences, one block of them at a
//!   time, with the preliminary stage's read errors;
//! - [`SubjectPlan`]: the scan ranges (`-db_soft_mask`) or the hard chunks (`-db_hard_mask`) of a
//!   sequence with masks, from the complement of its masks (`core::subject_split`);
//! - [`mask_protein_letters`]: the letters of NCBI's later fetches with `-db_hard_mask`;
//! - [`shown_subject_masks`]: the masks that the pairwise report shows in lower case.
//!
//! Nothing here runs in a search without the options: each program chooses its masked kernel once
//! per call.

use anyhow::Result;

use crate::blastinput::db_mask::{MaskKind, ResolvedMask};
use crate::core::subject_split::{complement_masks, HardChunk, HardSplitState};
use crate::seqdb::{AmbiguityRun, BlockMasks, BlockRange, MaskBlock, SeqDb, SeqDbError};
use crate::sequence::packed_subject::PackedNucleotide;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1631-1636
// A protein's copied fetch writes the letter 21 over every masked residue of the slice it reads
// (`s_SeqDBMaskSequence`, seqdbvol.cpp:1518-1539).
/// NCBIstdaa `X` (21): the letter of every masked residue in NCBI's fetches of a protein database
/// sequence after the preliminary stage of a search with `-db_hard_mask`.
pub(crate) const HARD_MASK_PROTEIN_LETTER: u8 = 21;

/// The database masking that a search reads ([`ResolvedMask`]): the algorithm's database-wide id
/// and the kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct SearchMask {
    /// The filtering algorithm (`MaskAlgorithms::ids`).
    pub(crate) algo_id: i32,
    /// `-db_soft_mask` or `-db_hard_mask`.
    pub(crate) kind: MaskKind,
}

impl SearchMask {
    /// The masking of a search resolved to `resolved` (`None`: a search without it, also after
    /// the formatter's warning).
    pub(crate) fn of(resolved: Option<&ResolvedMask>) -> Option<Self> {
        resolved.map(|mask| SearchMask {
            algo_id: mask.id,
            kind: mask.kind,
        })
    }

    /// Whether the search hard-masks its subjects (`-db_hard_mask`).
    pub(crate) fn hard(self) -> bool {
        self.kind == MaskKind::Hard
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321
// Every fetch of a sequence of a masked search reads the sequence's masks of the search's
// algorithm (`GetMaskData`); LOSAT reads them for a block of sequences at a time.
/// The stored masks of a database's sequences of one algorithm, read for a window of up to
/// `window` sequences at a time (`SeqDb::read_mask_block`): a preliminary stage reads its
/// sequences in order and takes a large window, a later stage reads them in the order of its
/// matches and takes a window of one sequence. One per worker.
#[derive(Debug)]
pub(crate) struct SubjectMasks {
    block: MaskBlock,
    /// The sequences that `block` holds (the ordinal id of the first, the end), once read.
    held: Option<BlockRange>,
    window: u32,
}

impl Default for SubjectMasks {
    fn default() -> Self {
        Self::with_window(PRELIMINARY_WINDOW)
    }
}

/// The window of [`SubjectMasks`] in a preliminary stage, which reads the sequences in order.
pub(crate) const PRELIMINARY_WINDOW: u32 = 4096;

impl SubjectMasks {
    /// A reader that reads the masks of `window` sequences at a time (at least one).
    pub(crate) fn with_window(window: u32) -> Self {
        SubjectMasks {
            block: MaskBlock::default(),
            held: None,
            window: window.max(1),
        }
    }

    /// The stored masks of database sequence `ordinal` for the algorithm `algo_id`
    /// (`[start, end)`, in the order stored; none for a sequence without mask data), with the
    /// reader's error as `error` makes it for `program` (`shared::search_read_error_of` in a
    /// preliminary stage, `shared::database_error_of` after it). A volume without the algorithm
    /// was refused before the search (`SeqDb::check_search_masks`).
    pub(crate) fn of(
        &mut self,
        db: &SeqDb,
        ordinal: u32,
        algo_id: i32,
        error: fn(&SeqDbError, &str) -> anyhow::Error,
        program: &str,
    ) -> Result<&[(u32, u32)]> {
        let held = self
            .held
            .filter(|range| range.first_ordinal <= ordinal && ordinal < range.end_ordinal());
        let range = match held {
            Some(range) => range,
            None => {
                let end = ordinal
                    .saturating_add(self.window)
                    .min(db.num_seqs())
                    .max(ordinal.saturating_add(1));
                let Some(range) = db.x_block_range(ordinal, end) else {
                    anyhow::bail!(
                        "LOSAT internal error: database sequence {ordinal} is not in the database"
                    );
                };
                self.held = Some(range);
                range
            }
        };
        db.read_mask_block(range, algo_id, &mut self.block);
        match self
            .block
            .ranges(range.first + (ordinal - range.first_ordinal))
        {
            Ok(BlockMasks::Ranges(stored)) => Ok(stored),
            Ok(BlockMasks::AlgorithmNotInVolume) => Err(anyhow::anyhow!(
                "LOSAT internal error: the volume of sequence {ordinal} has no mask algorithm \
                 {algo_id}"
            )),
            Err(read) => Err(error(read.error(), program)),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:162-176
// A soft-masked sequence's complement ranges are its scan ranges; a hard-masked one's are its
// hard ranges, each searched as its own chunk (`s_GetNextSubjectChunk`, blast_engine.c:221-311),
// whose HSPs move by the chunk's offset and join the sequence's list (blast_engine.c:575-584).
/// How a database sequence with masks is searched: its scan ranges or its hard chunks. A sequence
/// without masks is searched as in a search without the options (its one range is the whole
/// sequence, which every scan reads without ranges).
#[derive(Debug, Default)]
pub(crate) struct SubjectPlan {
    /// `-db_soft_mask`: the scan ranges (`seq_ranges`), the complement of the masks in the order
    /// stored (`complement_masks`).
    pub(crate) ranges: Vec<(i32, i32)>,
    /// `-db_hard_mask`: the chunker over the hard ranges (residual 0: the protein programs').
    pub(crate) split: HardSplitState,
}

impl SubjectPlan {
    /// The scan ranges of a sequence of `length` letters with the stored masks `masks` (never
    /// empty), for `-db_soft_mask`.
    pub(crate) fn soft(&mut self, masks: &[(u32, u32)], length: usize) -> &[(i32, i32)] {
        complement_masks(masks, length, &mut self.ranges);
        &self.ranges
    }

    /// Starts the hard chunks of a sequence of `length` letters with the stored masks `masks`
    /// (never empty), for `-db_hard_mask`; [`SubjectPlan::next_chunk`] gives them in order.
    pub(crate) fn start_hard(&mut self, masks: &[(u32, u32)], length: usize) {
        let mut ranges = std::mem::take(&mut self.split.ranges);
        complement_masks(masks, length, &mut ranges);
        self.split.start(length, ranges, false);
    }

    /// The next hard chunk (`None`: none left), with the chunk overlap of the protein programs
    /// (`DBSEQ_CHUNK_OVERLAP`; no merge across a mask reads it, `merge_chunk_lists`).
    pub(crate) fn next_chunk(&mut self) -> Option<HardChunk> {
        self.split
            .next_chunk(crate::core::subject_split::DBSEQ_CHUNK_OVERLAP)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1539
// The mask letter goes over each mask's letters inside the fetched range, the masks taken in
// the order stored from the first that ends after the range's start until the first that starts
// at or after its end.
/// Writes [`HARD_MASK_PROTEIN_LETTER`] over the residues of `residues` (a whole protein sequence,
/// from its first residue) that the stored masks `masks` cover, as NCBI's fetch of the sequence
/// after the preliminary stage does with `-db_hard_mask` (masks past the sequence are cut at its
/// end; `HardMaskShape::Unsorted` refused the masks out of order before the search).
pub(crate) fn mask_protein_letters(residues: &mut [u8], masks: &[(u32, u32)]) {
    let end = residues.len();
    let mut i = 0;
    while i < masks.len() && masks[i].1 == 0 {
        i += 1;
    }
    while i < masks.len() && (masks[i].0 as usize) < end {
        let from = masks[i].0 as usize;
        let to = (masks[i].1 as usize).min(end);
        if from < to {
            residues[from..to].fill(HARD_MASK_PROTEIN_LETTER);
        }
        i += 1;
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1539, 1672
// A nucleotide's copied fetch writes ncbi4na B (C, G or T) over every masked base after its
// ambiguities.
/// ncbi4na `B` (14): the code of every masked base in NCBI's fetches of a nucleotide database
/// sequence after the preliminary stage of a search with `-db_hard_mask` (TBLASTN, TBLASTX).
pub(crate) const HARD_MASK_NUCLEOTIDE_CODE: u8 = 14;

/// The store `store` of a nucleotide database sequence as NCBI's fetches after the preliminary
/// stage read it with `-db_hard_mask`: its ambiguity runs, then [`HARD_MASK_NUCLEOTIDE_CODE`] over
/// the bases that the stored masks `masks` cover, taken in the order of
/// [`mask_protein_letters`] (masks past the sequence are cut at its end). The packed bases are
/// the store's: the preliminary stage reads them unmasked.
pub(crate) fn mask_nucleotide_store(
    store: &PackedNucleotide,
    masks: &[(u32, u32)],
) -> PackedNucleotide {
    let end = store.len();
    let mut runs = store.ambiguity_runs().to_vec();
    let mut i = 0;
    while i < masks.len() && masks[i].1 == 0 {
        i += 1;
    }
    while i < masks.len() && (masks[i].0 as usize) < end {
        let from = masks[i].0;
        let to = (masks[i].1 as usize).min(end) as u32;
        if from < to {
            runs.push(AmbiguityRun {
                offset: from,
                length: to - from,
                ncbi4na: HARD_MASK_NUCLEOTIDE_CODE,
            });
        }
        i += 1;
    }
    PackedNucleotide::x_from_stored_buffer(store.packed().to_vec(), end, &runs)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1586-1601 and
// c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:109-139
// The report's masks of a database sequence are those of its stored masks that hold a letter
// and meet the closed range from an HSP's subject offset to its end (one past its last letter).
/// The stored masks `masks` of a database sequence that the pairwise report shows in lower case
/// with the HSPs whose subject ranges are `hsp_ranges` (`[offset, end)` each, end one past the
/// last letter): each mask with a letter that meets the closed range `[offset, end]` of one of
/// them, in the order stored.
pub(crate) fn shown_subject_masks(
    masks: &[(u32, u32)],
    hsp_ranges: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    masks
        .iter()
        .map(|&(start, end)| (start as usize, end as usize))
        .filter(|&(start, end)| {
            start < end
                && hsp_ranges
                    .iter()
                    .any(|&(offset, hsp_end)| start <= hsp_end && offset < end)
        })
        .collect()
}

#[cfg(test)]
#[path = "db_mask_search_tests.rs"]
mod db_mask_search_tests;
