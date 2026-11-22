//! The input of `makeblastdb -input_type blastdb`: the sequences of existing databases as they
//! are stored (`CRawSeqDBSource`, `c++/src/app/blastdb/makeblastdb.cpp`), with their definition
//! lines, the ranges of their masking algorithms and the blobs of their columns, which
//! `BuildDatabase::add_raw_sequences` writes; and the scan of the input's leaf taxids that
//! makeblastdb makes before the copy.
//!
//! The reader is LOSAT's (`crate::seqdb`), the one owner of the input's files: the sequences are
//! read a block at a time (`SeqDb::read_block`), the masks and the blobs of each block with one
//! read of each column (`SeqDb::read_mask_block`, `SeqDb::read_column_block`), and every error
//! is given at the sequence where NCBI's reads of one sequence at a time meet it.

use std::collections::{BTreeMap, BTreeSet};

use super::build_db::BuildDatabase;
use super::mask_data::registry::FILTER_PROGRAM_OTHER;
use super::mask_data::MaskedRanges;
use super::WriteDbError;
use crate::ncbi_str::string_to_int;
use crate::seqdb::column::MASK_DATA_TITLE;
use crate::seqdb::mask::BlockMaskError;
use crate::seqdb::{
    BlockMasks, BlockRange, ColumnBlock, MaskBlock, SeqDb, SeqDbError, SeqType, SequenceBlock,
};
use crate::seqid::DeflineSet;
use crate::utils::std_shuffle::{shuffle, Mt19937};

/// The bytes of the sequence file read at a time in the input's order: a block of one volume,
/// at least one sequence (`SeqDb::next_block`). In the order of `-randomize` the next sequence is
/// anywhere, and a block is the one sequence.
const BLOCK_BYTES: u64 = 1 << 20;

/// Why the copy of a database stopped.
#[derive(Debug)]
pub(crate) enum RawCopyError {
    /// An error of the reader of the input databases.
    Read(SeqDbError),
    /// An error of the build: of the writer, of makeblastdb's checks, or LOSAT's refusal.
    Build(WriteDbError),
}

impl From<SeqDbError> for RawCopyError {
    fn from(error: SeqDbError) -> Self {
        RawCopyError::Read(error)
    }
}

impl From<WriteDbError> for RawCopyError {
    fn from(error: WriteDbError) -> Self {
        RawCopyError::Build(error)
    }
}

/// makeblastdb's map of leaf taxids (`TIdToLeafs`): the `AsFastaString` of the first Seq-id of
/// a definition line to the leaf taxids of the definition lines that start with it.
pub(crate) type LeafTaxIds = BTreeMap<Vec<u8>, BTreeSet<i32>>;

/// LOSAT's refusal of the definition lines `set` of the input's ordinal id `ordinal`, before
/// makeblastdb reads them: lines with bytes NCBI reads as '#' (its reader writes an error for
/// them each time it reads them, twice here), a set without lines (an empty header is a null
/// object that makeblastdb dereferences; an empty set makes it read the first line of none),
/// and a line with leaf taxids and no Seq-id (makeblastdb reads its first Seq-id).
fn refuse_deflines(db: &SeqDb, ordinal: u32, set: &DeflineSet) -> Result<(), WriteDbError> {
    let refuse = |what: &str| {
        Err(WriteDbError::Unsupported(format!(
            "OID {} of the input database: {what}",
            db.oid(ordinal).get()
        )))
    };
    if set.deflines.is_empty() {
        return refuse(
            "the sequence has no definition line, which NCBI makeblastdb cannot copy (it reads \
             the first line of none)",
        );
    }
    if set
        .deflines
        .iter()
        .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
    {
        return refuse(
            "a definition line with bytes NCBI reads as '#' (it writes an error each time it \
             reads them)",
        );
    }
    if set.deflines.iter().any(|defline| {
        defline.ids.is_empty()
            && defline
                .links
                .as_ref()
                .is_some_and(|links| !links.is_empty())
    }) {
        return refuse(
            "a definition line with leaf taxids and no Seq-id, whose first Seq-id NCBI \
             makeblastdb reads",
        );
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1101-1113
// NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:149-163
// makeblastdb reads every header of the input once before the copy (with a reader of its own,
// opened on the same names) and keys the leaf taxids (`links`) of each definition line by its
// first Seq-id's `AsFastaString`. The copy looks the keys of each line up by `GetSeqIdKey`
// (`BuildDatabase::add_raw_sequences`), a label of another form, so the map changes a line only
// where the two forms meet.
/// The scan of makeblastdb's leaf taxids over the input `db`, in the order of its sequences. Its
/// errors are those of the first read of each header; LOSAT's refusals of a set of definition
/// lines come here, before anything is written.
pub(crate) fn scan_leaf_tax_ids(db: &SeqDb) -> Result<LeafTaxIds, RawCopyError> {
    let mut leafs = LeafTaxIds::new();
    for ordinal in 0..db.num_seqs() {
        let set = db.defline_set(ordinal)?;
        refuse_deflines(db, ordinal, &set)?;
        for defline in &set.deflines {
            let Some(links) = defline.links.as_ref().filter(|links| !links.is_empty()) else {
                continue;
            };
            let mut key = Vec::new();
            crate::seqid::label::fasta(&defline.ids[0], &mut key);
            leafs.entry(key).or_default().extend(links.iter().copied());
        }
    }
    Ok(leafs)
}

/// One sequence of the input (`GetNext`'s outputs).
pub(crate) struct RawSequence<'s> {
    /// The sequence as stored: a protein's NCBIstdaa bytes without the separating zero bytes, a
    /// nucleotide's packed bases with the count of the last byte.
    pub(crate) sequence: &'s [u8],
    /// A nucleotide's ambiguity data as stored; empty for a protein.
    pub(crate) ambiguities: &'s [u8],
    /// The definition lines (`GetHdr`: a `gnl|BL_ORD_ID|n` holds the OID in the input).
    pub(crate) deflines: DeflineSet,
    /// The ranges of each algorithm of the input, in the order of the input's ids, under the id
    /// the output gave the algorithm; empty where the sequence has none.
    pub(crate) masks: Vec<MaskedRanges>,
    /// The blob of each column of [`RawSeqDbSource::column_titles`], as stored; empty where the
    /// sequence has none (and always for the mask column of an input with algorithms, whose
    /// ranges are `masks`).
    pub(crate) column_blobs: Vec<Vec<u8>>,
}

/// `CRawSeqDBSource` over the input database `db` (all the names of `-in`).
pub(crate) struct RawSeqDbSource<'a> {
    db: &'a SeqDb,
    /// `m_MaskIds` and `m_MaskIdMap`: each algorithm id of the input
    /// (`GetAvailableMaskAlgorithms`) and the id the output registered it under.
    masks: Vec<(i32, i32)>,
    /// `m_ColumnNames`: the titles of the input's columns (`ListColumns`).
    columns: Vec<Vec<u8>>,
    /// `m_RandomOids` of `-randomize`: the ordinal ids in the order of the copy (`CheckOrFindOID`
    /// walks the OIDs the input's alias files include, which are LOSAT's ordinal ids); `None` for
    /// the input's order.
    order: Option<Vec<u32>>,
    /// `m_Oid`: the place of the next sequence in the walk, its ordinal id in the input's order.
    next: u32,
    /// The block read and its sequences.
    block: Option<(BlockRange, SequenceBlock)>,
    /// The masks of the block, one per algorithm of `masks`.
    mask_blocks: Vec<MaskBlock>,
    /// The blobs of the block, one per column of `columns`.
    column_blocks: Vec<ColumnBlock>,
    /// NCBI's memory of its last lookup of a volume's id of an algorithm.
    last_lookup: VolAlgoMemory,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2234-2242
/// `CSeqDB_IdRemapper::GetVolAlgo`'s memory: the volume and the algorithm (its place in the
/// input's list) asked last, and whether that volume knows that algorithm. The question is
/// kept before the lookup, which may fail; the id found is kept only when it does not.
#[derive(Default)]
struct VolAlgoMemory(Option<(usize, usize, bool)>);

impl VolAlgoMemory {
    /// The lookup of the algorithm `algorithm` in the volume `volume`, which knows it or not
    /// (`known`): `true` when NCBI answers from its memory with the id of another lookup (the
    /// same question as the last one, which failed).
    fn reads_a_kept_id(&mut self, volume: usize, algorithm: usize, known: bool) -> bool {
        let kept = !known && self.0 == Some((volume, algorithm, false));
        self.0 = Some((volume, algorithm, known));
        kept
    }
}

impl<'a> RawSeqDbSource<'a> {
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:640-675
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1263-1274
    // The program of an algorithm is its description's first item read as a number
    // (`StringToNumeric` without throwing: 0 for one that is not). The output's name of an
    // algorithm other than `other` is its program's name followed by the input's id ("dust11"),
    // so that an algorithm of the input registers apart from the same one of `-mask_data`, which
    // registered before; `other`'s name is its id, and its description is not passed on.
    //
    // With `-randomize` the copy walks the OIDs the input includes in the order `std::shuffle`
    // leaves them with a `std::mt19937` seeded with the seconds of the clock (`time` with a null
    // pointer; `crate::utils::std_shuffle`; the number of OIDs that ends NCBI's list is the end of
    // the walk here).
    /// The constructor: registers the masking algorithms of `db` on the output `outdb`, in the
    /// order of their ids, lists the columns of `db` and, with `-randomize`, shuffles the order
    /// of the copy with the seed `randomize` (the seconds of the clock, modulo 2^32).
    pub(crate) fn new(
        db: &'a SeqDb,
        outdb: &mut BuildDatabase<'_>,
        randomize: Option<u32>,
    ) -> Result<Self, RawCopyError> {
        let algorithms = db.mask_algorithms()?;
        let mut masks = Vec::new();
        for id in algorithms.ids() {
            let details = algorithms.details(id)?;
            let program = string_to_int(details.program.as_bytes(), 0).unwrap_or(0);
            let mut name = details.program_name.into_bytes();
            if program != FILTER_PROGRAM_OTHER {
                name.extend_from_slice(id.to_string().as_bytes());
            }
            let output_id =
                outdb.register_masking_algorithm(program, details.options.as_bytes(), &name)?;
            masks.push((id, output_id));
        }
        let columns = db.column_titles()?;
        let order = randomize.map(|seed| {
            let mut order: Vec<u32> = (0..db.num_seqs()).collect();
            shuffle(&mut order, &mut Mt19937::new(seed));
            order
        });
        Ok(Self {
            db,
            mask_blocks: masks.iter().map(|_| MaskBlock::default()).collect(),
            column_blocks: columns.iter().map(|_| ColumnBlock::default()).collect(),
            masks,
            columns,
            order,
            next: 0,
            block: None,
            last_lookup: VolAlgoMemory::default(),
        })
    }

    /// `GetColumnNames`: the titles of the input's columns, in their byte order.
    pub(crate) fn column_titles(&self) -> &[Vec<u8>] {
        &self.columns
    }

    /// `GetColumnMetaData` of the `GetColumnId` of a title, for the column `index` of
    /// [`Self::column_titles`].
    pub(crate) fn column_meta_data(
        &self,
        index: usize,
    ) -> Result<BTreeMap<Vec<u8>, Vec<u8>>, SeqDbError> {
        self.db.column_meta_data(&self.columns[index])
    }

    /// The block that holds the ordinal id `ordinal`, read with its masks and blobs when it is
    /// not the one held. A block that cannot be read is read again as its first sequence alone,
    /// so that an error is the one of the sequence NCBI would meet it at.
    fn block_of(&mut self, ordinal: u32) -> Result<BlockRange, SeqDbError> {
        if let Some((range, _)) = &self.block {
            if range.first_ordinal <= ordinal && ordinal < range.end_ordinal() {
                return Ok(*range);
            }
        }
        self.block = None;
        let bytes = if self.order.is_some() { 0 } else { BLOCK_BYTES };
        let mut range = self.db.next_block(ordinal, bytes).ok_or_else(|| {
            SeqDbError::Range(format!(
                "ordinal id {ordinal} is outside the database ({} sequences)",
                self.db.num_seqs()
            ))
        })?;
        let block = match self.db.read_block(range) {
            Ok(block) => block,
            Err(error) if range.len() == 1 => return Err(error),
            Err(_) => {
                range.end = range.first + 1;
                self.db.read_block(range)?
            }
        };
        for (block, &(algo_id, _)) in self.mask_blocks.iter_mut().zip(&self.masks) {
            self.db.read_mask_block(range, algo_id, block);
        }
        for (block, title) in self.column_blocks.iter_mut().zip(&self.columns) {
            if title == MASK_DATA_TITLE && !self.masks.is_empty() {
                continue;
            }
            self.db.read_column_block(title, range, block);
        }
        self.block = Some((range, block));
        Ok(range)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:678-738
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2997-3072
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2234-2242,2531-2555
    // Each sequence is read in NCBI's order: its data (a protein without residues is NCBI's
    // file error), its header, its masks for each algorithm, then its blob of each column.
    //
    // The volume of a sequence with mask data may not know an algorithm of the database: NCBI's
    // lookup of the volume's id then fails and the sequence has no ranges of it. The lookup keeps
    // the last volume and algorithm it was asked for, and the volume's id it found for them, and
    // answers a repeated question from that memory without asking again; a failed lookup keeps
    // the question and the id found before it. With two algorithms or more the questions
    // alternate, so every failed lookup is asked again; with one, a sequence with mask data whose
    // volume was the last one asked, and failed, reads the ranges of the id remembered from
    // another lookup, which LOSAT refuses. (In the input's order that is the second sequence
    // with mask data in such a volume; with `-randomize`, unless a sequence with mask data of
    // another volume came between them.)
    /// `GetNext`: the next sequence of the input, `None` after the last.
    pub(crate) fn next(&mut self) -> Result<Option<RawSequence<'_>>, RawCopyError> {
        let ordinal = match &self.order {
            Some(order) => match order.get(self.next as usize) {
                Some(&ordinal) => ordinal,
                None => return Ok(None),
            },
            None if self.next < self.db.num_seqs() => self.next,
            None => return Ok(None),
        };
        self.next += 1;
        let range = self.block_of(ordinal)?;
        let i = (ordinal - range.first_ordinal) as usize;
        let vol_oid = range.first + i as u32;
        let protein = self.db.seqtype() == SeqType::Protein;
        if let Some((_, block)) = &self.block {
            if protein && block.protein(i).is_empty() {
                return Err(WriteDbError::Database(
                    "File error: could not get sequence data.".to_string(),
                )
                .into());
            }
        }
        let deflines = self.db.defline_set(ordinal)?;

        let mut masks = Vec::with_capacity(self.masks.len());
        for k in 0..self.masks.len() {
            let offsets = match self.mask_blocks[k].ranges(vol_oid) {
                Ok(BlockMasks::NoMaskData) => Vec::new(),
                Ok(BlockMasks::Ranges(ranges)) => {
                    self.last_lookup.reads_a_kept_id(range.volume, k, true);
                    ranges.to_vec()
                }
                Ok(BlockMasks::AlgorithmNotInVolume) => {
                    if self.last_lookup.reads_a_kept_id(range.volume, k, false) {
                        return Err(WriteDbError::Unsupported(format!(
                            "OID {} of the input database: its volume does not know the \
                             database's only masking algorithm, and NCBI then reads the masks \
                             of an algorithm id it kept from another sequence",
                            self.db.oid(ordinal).get()
                        ))
                        .into());
                    }
                    Vec::new()
                }
                Err(BlockMaskError::Ranges(error)) => return Err(error.into()),
                Err(BlockMaskError::Read(_)) => {
                    let block = std::mem::take(&mut self.mask_blocks[k]);
                    return Err(block.into_read_error(vol_oid).into());
                }
            };
            masks.push(MaskedRanges {
                algorithm_id: self.masks[k].1,
                offsets,
            });
        }

        let mut column_blobs = Vec::with_capacity(self.columns.len());
        for (block, title) in self.column_blocks.iter_mut().zip(&self.columns) {
            if title == MASK_DATA_TITLE && !self.masks.is_empty() {
                column_blobs.push(Vec::new());
                continue;
            }
            column_blobs.push(block.blob(vol_oid)?.to_vec());
        }

        let Some((_, block)) = &self.block else {
            unreachable!("block_of holds the block");
        };
        let (sequence, ambiguities): (&[u8], &[u8]) = if protein {
            (block.protein(i), &[])
        } else {
            (block.nucleotide_packed(i), block.ambiguity_bytes(i))
        };
        Ok(Some(RawSequence {
            sequence,
            ambiguities,
            deflines,
            masks,
            column_blobs,
        }))
    }
}

#[cfg(test)]
#[path = "raw_source_tests.rs"]
mod tests;
