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
//! - [`shown_subject_masks`], [`read_shown_stored`] and [`shown_db_subject_masks`]: the masks
//!   that the pairwise report shows in lower case, by NCBI's Seq-id ([`subject_id_keys`]).
//!
//! Nothing here runs in a search without the options: each program chooses its masked kernel once
//! per call.

use std::collections::HashMap;

use anyhow::Result;

use crate::blastinput::db_mask::{MaskKind, ResolvedMask};
use crate::core::subject_split::{complement_masks, HardChunk, HardSplitState};
use crate::report::pairwise::PairwiseHit;
use crate::seqdb::{AmbiguityRun, BlockMasks, BlockRange, MaskBlock, SeqDb, SeqDbError};
use crate::seqid::defline_set::ordered_key;
use crate::seqid::model::{Dbtag, ObjectId};
use crate::seqid::ParsedSeqId;
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
            Ok(BlockMasks::NoMaskData) => Ok(&[]),
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

    /// The next hard chunk (`None`: none left), with NCBI's chunk overlap (`DBSEQ_CHUNK_OVERLAP`,
    /// which NCBI's merge takes after every chunk that ends a range, as it compares the chunk's
    /// offset with the next range's left). BLASTP and BLASTX search every chunk, whatever its
    /// length (`algorithm::common::hard_chunk_walk::DiagonalWalk::chunk`), and merge each chunk's
    /// HSPs across `chunk.overlap` via `core::blast_hits::merge_subject_chunk_hsps`
    /// (`Blast_HSPListsMerge`, series NP, M4b-3d).
    pub(crate) fn next_chunk(&mut self) -> Option<HardChunk> {
        self.split
            .next_chunk(crate::core::subject_split::DBSEQ_CHUNK_OVERLAP)
    }

    /// Every hard chunk of the sequence started by [`SubjectPlan::start_hard`], in order, into
    /// `chunks` (which it replaces).
    pub(crate) fn hard_chunks(&mut self, chunks: &mut Vec<HardChunk>) {
        chunks.clear();
        while let Some(chunk) = self.next_chunk() {
            chunks.push(chunk);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1539
// The mask letter goes over each mask's letters inside the fetched range, the masks taken in
// the order stored from the first that ends after the range's start until the first that starts
// at or after its end.
/// Writes [`HARD_MASK_PROTEIN_LETTER`] over the residues of `residues` (a whole protein sequence,
/// from its first residue) that the stored masks `masks` cover, as NCBI's fetch of the sequence
/// after the preliminary stage does with `-db_hard_mask` (masks past the sequence are cut at its
/// end; masks that overlap or nest, or are out of the order of their start, are each written, and
/// NCBI's scan stops at the first mask that starts at or past the end).
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

/// The stored masks of the algorithm `algo_id` of each subject rank of `hits` (the pairwise rows),
/// read from `db` at the rank's ordinal id (`ordinal`; a rank without one has none); `program`
/// names the read errors (`database_error_of`).
pub(crate) fn read_shown_stored(
    hits: &[PairwiseHit],
    db: &SeqDb,
    algo_id: i32,
    program: &str,
    ordinal: impl Fn(u32) -> Option<u32>,
) -> Result<HashMap<u32, Vec<(u32, u32)>>> {
    let mut reader = SubjectMasks::with_window(1);
    let mut stored = HashMap::new();
    for hit in hits {
        let rank = hit.hit.s_idx;
        if stored.contains_key(&rank) {
            continue;
        }
        let Some(ordinal) = ordinal(rank) else {
            continue;
        };
        let masks = reader.of(
            db,
            ordinal,
            algo_id,
            crate::algorithm::common::database::database_error_of,
            program,
        )?;
        stored.insert(rank, masks.to_vec());
    }
    Ok(stored)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:115-116
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:51-69
// The masks of a sequence carry its first Seq-id; its Seq-aligns carry its first GI, else that
// first Seq-id. A database without parsed Seq-ids names each sequence `gnl|BL_ORD_ID|oid`.
/// The [`ordered_key`]s of the first Seq-id of a database sequence (`CSeqDB::GetSeqIDs` of `oid`,
/// its first: NCBI's subject masks are kept under it) and of its Seq-align id
/// (`DeflineSet::alignment_id`: a row looks the masks up by it): those of its parsed Seq-ids
/// `parsed`, or `gnl|BL_ORD_ID|oid` for both.
pub(crate) fn subject_id_keys(
    parsed: Option<&crate::report::db_subject::ParsedSubject>,
    oid: crate::seqdb::Oid,
) -> (Vec<u8>, Vec<u8>) {
    match parsed {
        Some(parsed) => {
            let alignment = ordered_key(parsed.alignment_id());
            let first = parsed
                .deflines()
                .iter()
                .flat_map(|defline| defline.ids.iter())
                .next()
                .map_or_else(|| alignment.clone(), ordered_key);
            (first, alignment)
        }
        None => {
            let id = ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(oid.get() as i32),
            });
            let key = ordered_key(&id);
            (key.clone(), key)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1568-1601
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1046-1052 and 1247-1250
// A query's subject masks go, in the order of its hit list, into a map by the first Seq-id of
// each sequence (a sequence with none that meets its HSPs adds nothing; a later sequence's
// replace an earlier one's of an equal Seq-id), and each row of the pairwise report shows the
// entry of its Seq-align id. Two database sequences share a Seq-id only across databases with
// parsed Seq-ids (`-db "a b"`, an alias of both): each row then shows the masks of the last.
/// The subject masks that the pairwise rows `hits` (in the order of the hit lists) show in lower
/// case, by (query, subject rank), inclusive: each subject's stored masks (`stored`, by rank) that
/// meet its HSPs with the query ([`shown_subject_masks`]), under the first Seq-id of its
/// [`subject_id_keys`] (`id_keys`), read back by its Seq-align id.
pub(crate) fn shown_db_subject_masks<'a>(
    hits: &[PairwiseHit],
    stored: impl Fn(u32) -> Option<&'a [(u32, u32)]>,
    mut id_keys: impl FnMut(u32) -> (Vec<u8>, Vec<u8>),
) -> HashMap<(u32, u32), Vec<(i32, i32)>> {
    // Each query's subject ranks in the order of the rows, and the HSP ranges of each pair.
    let mut queries: Vec<(u32, Vec<u32>)> = Vec::new();
    let mut query_at: HashMap<u32, usize> = HashMap::new();
    let mut hsp_ranges: HashMap<(u32, u32), Vec<(usize, usize)>> = HashMap::new();
    for hit in hits {
        let (query, rank) = (hit.hit.q_idx, hit.hit.s_idx);
        let ranges = hsp_ranges.entry((query, rank)).or_insert_with(|| {
            let at = *query_at.entry(query).or_insert_with(|| {
                queries.push((query, Vec::new()));
                queries.len() - 1
            });
            queries[at].1.push(rank);
            Vec::new()
        });
        ranges.push((hit.hit.s_start - 1, hit.hit.s_end));
    }
    let mut keys: HashMap<u32, (Vec<u8>, Vec<u8>)> = HashMap::new();
    let mut shown = HashMap::new();
    for (query, ranks) in &queries {
        let mut by_id: HashMap<Vec<u8>, Vec<(i32, i32)>> = HashMap::new();
        for &rank in ranks {
            let Some(masks) = stored(rank) else {
                continue;
            };
            let runs: Vec<(i32, i32)> = shown_subject_masks(masks, &hsp_ranges[&(*query, rank)])
                .into_iter()
                .map(|(start, end)| (start as i32, end as i32 - 1))
                .collect();
            if !runs.is_empty() {
                let first = keys.entry(rank).or_insert_with(|| id_keys(rank)).0.clone();
                by_id.insert(first, runs);
            }
        }
        if by_id.is_empty() {
            continue;
        }
        for &rank in ranks {
            let alignment = &keys.entry(rank).or_insert_with(|| id_keys(rank)).1;
            if let Some(runs) = by_id.get(alignment) {
                shown.insert((*query, rank), runs.clone());
            }
        }
    }
    shown
}

#[cfg(test)]
#[path = "db_mask_search_tests.rs"]
mod db_mask_search_tests;
