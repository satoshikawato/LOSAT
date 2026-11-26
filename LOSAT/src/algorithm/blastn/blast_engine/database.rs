//! The database path of BLASTN (`-db NAME`; tasks megablast, blastn, dc-megablast and
//! blastn-short): the search of an NCBI BLAST nucleotide database read in blocks through
//! `crate::seqdb`, as NCBI BLAST+ 2.17.0 `blastn -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, stage DB-4).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the set-level numbers come from the index files (`BlastnDatabase::metadata`): the
//!   number of sequences and the total length (an `Int8`, not the `int` of the subject
//!   set's report);
//! - the preliminary stage reads the stored 2-bit data in blocks and searches them in the
//!   order of their ordinal ids (`search_blocks`); a sequence's preliminary HSPs go into
//!   the preliminary hit lists as soon as its block is handed over (LOSAT_X_BNDBSTREAM; without
//!   it, when its block's round is searched), so nothing but the hit lists is kept for the
//!   whole database;
//! - the traceback reads again, by ordinal id, the sequences of the preliminary hit lists
//!   (`traceback_listed`), with their ambiguity codes (BLASTNA);
//! - the definition lines of the final hit lists are read after each query batch's
//!   traceback (`read_titles`), and the reports show the database
//!   (`BlastnDatabase::pairwise_database`, `name`) and the stored titles and letters of the
//!   reported sequences only (`reported_subjects`).
//!
//! The per-subject kernel, the hit lists and the report writers are those of the
//! `-subject` path (`subject_passes.rs`, a child of `run.rs`); the opening of the database, the reader's errors and the
//! rules for the stored titles are shared with BLASTP (`algorithm::common::database`).

use std::cell::RefCell;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

use anyhow::{Context, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::algorithm::blastn::coordination::SubjectMetadata;
use crate::algorithm::blastn::hsp::{BlastnHitList, BlastnHspList, BlastnOutputFormat};
pub(super) use crate::algorithm::common::database::PreliminaryReadError;
use crate::algorithm::common::database::{self as shared, ParsedShown, ReportedHeader};
use crate::blastinput::db_mask::{DbMaskOption, ResolvedMask};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::ncbi_exception::NcbiException;
use crate::report::db_subject::{HitLoads, ParsedSubject};
use crate::report::pairwise::PairwiseDatabase;
use crate::report::tabular::{PrintedIds, SubjectFetch};
use crate::seqdb::{
    BlockMasks, BlockRange, MaskBlock, NucleotideAmbiguity, SeqDb, SeqDbError, SeqType,
    SequenceBlock,
};
use crate::sequence::packed_subject::PackedNucleotide;
use crate::utils::threading::SearchPool;

use super::BlastnSubjectRef;

/// The program that LOSAT's refusals name.
const PROGRAM: &str = "BLASTN";

/// Bytes of the sequence file that one block of the preliminary stage holds at most
/// (`SeqDb::next_block`; a longer sequence is a block of its own). No NCBI counterpart:
/// NCBI maps the volume files and hands its threads chunks of ordinal ids
/// (seqdbimpl.cpp:328-350). The size does not change which sequences are searched or
/// their order.
const BLOCK_BYTES: u64 = 1 << 20;

/// Blocks that one parallel round of the preliminary stage reads and searches per
/// thread before their results are merged in ordinal order. No NCBI counterpart (as
/// `BLOCK_BYTES`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const BLOCKS_PER_THREAD: usize = 8;

// EXPERIMENT (LOSAT_X_BNDBCUT): the blocks of the preliminary stage cut
// by the database size and the thread count (a child module, so that it uses the private
// driver of this file).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_db_cut.rs"]
mod x_db_cut;

// EXPERIMENT (LOSAT_X_BNDBSTREAM): the parallel preliminary stage without rounds, the blocks'
// results handed to the hit lists in ordinal order while the threads search the next blocks (a
// child module, so that it uses the private block reader of this file).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_db_stream.rs"]
mod x_db_stream;

/// No NCBI counterpart (scheduling only): the heap bytes of a result of the preliminary stage
/// while it waits for the owner of `x_db_stream.rs`'s window.
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    allow(dead_code)
)]
pub(super) trait XWaitingBytes {
    fn x_waiting_bytes(&self) -> usize;
}

impl<T> XWaitingBytes for Vec<T> {
    fn x_waiting_bytes(&self) -> usize {
        self.capacity().saturating_mul(std::mem::size_of::<T>())
    }
}

#[cfg(test)]
impl XWaitingBytes for u32 {
    fn x_waiting_bytes(&self) -> usize {
        0
    }
}

#[cfg(test)]
impl XWaitingBytes for u64 {
    fn x_waiting_bytes(&self) -> usize {
        0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:75-78
/// The IUPAC letter of a BLASTNA code.
#[cfg(test)]
const BLASTNA_TO_IUPACNA: [u8; 16] = *b"ACGTRYMKWSBDHVN-";

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:375,388-391
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1376-1378
// NCBI's preliminary stage reads each sequence in place in the mapped volume (`GetFileDataPtr`,
// seqdbvol.cpp:1817, cited above `block_subject`), and each search thread allocates its search
// structures once. LOSAT reads a block with one positioned read (`seqdb::file`, no mapping, so
// that native builds and WebAssembly share the reader): each worker keeps one buffer for its
// blocks, read again for each block. It is allocated once per worker and grows only for a larger
// block, and only its bytes past the last block's length are zeroed.
/// `SeqDb::read_block` (the preliminary stage's blocks). With LOSAT_X_BNMAPPED, into the worker's
/// block buffer `kept` (`SeqDb::x_read_block_into`: below 1 MiB the worker's own buffer, kept from
/// one block to the next; a block of 1 MiB or more, a long sequence or at one thread a block of
/// exactly 1 MiB, in a buffer of `utils::x_mapped::X_BYTES`); `x_done_kept_block` gives the
/// buffer back. Out of line, as the read of a block was before: inlined into the one-thread loop
/// (`search_blocks_of`), its call of `Volume::read_block_in` with arguments on the stack made the
/// loop pass the stack argument of every subject's `search_subject` call with pushes (+3
/// instructions per subject on syn megablast).
#[inline(never)]
fn x_read_block(
    db: &SeqDb,
    range: BlockRange,
    kept: &mut Vec<u8>,
) -> Result<SequenceBlock, SeqDbError> {
    if super::x_bn_mapped() {
        db.x_read_block_into(range, kept)
    } else {
        db.read_block(range)
    }
}

/// With LOSAT_X_BNMAPPED, gives the buffer of a block of `x_read_block` back for the worker's next
/// block (`SequenceBlock::x_into_kept`: into `kept` below 1 MiB, else to
/// `utils::x_mapped::X_BYTES`).
fn x_done_kept_block(block: SequenceBlock, kept: &mut Vec<u8>) {
    if super::x_bn_mapped() {
        block.x_into_kept(kept);
    }
}

/// No NCBI counterpart (where LOSAT keeps its bytes): a worker of the parallel preliminary stage
/// (the rounds' `map_init`, `x_db_stream.rs`'s window), its search scratch and its block buffer
/// (`x_read_block`), both kept from one block to the next.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XWorker<T> {
    scratch: T,
    bytes: Vec<u8>,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<T> XWorker<T> {
    fn new(scratch: T) -> Self {
        Self {
            scratch,
            bytes: Vec::new(),
        }
    }
}

/// `SeqDb::read_sequence`, mapped with LOSAT_X_BNMAPPED.
fn x_read_sequence(db: &SeqDb, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
    if super::x_bn_mapped() {
        db.x_read_sequence_mapped(ordinal)
    } else {
        db.read_sequence(ordinal)
    }
}

/// With LOSAT_X_BNMAPPED, gives a sequence read by `x_read_sequence` back to
/// `utils::x_mapped::X_BYTES`.
fn x_done_block(block: SequenceBlock) {
    if super::x_bn_mapped() {
        block.x_recycle();
    }
}

/// With LOSAT_X_BNMAPPED, gives letters of `x_blastna` back to `utils::x_mapped::X_BYTES`.
fn x_done_bytes(mut bytes: Vec<u8>) {
    if super::x_bn_mapped() {
        super::x_bn_mapped_reach(bytes.len());
        crate::utils::x_mapped::X_BYTES.give(&mut bytes);
    }
}

/// `SequenceBlock::nucleotide_blastna(0)`, mapped with LOSAT_X_BNMAPPED.
fn x_blastna(block: &SequenceBlock) -> Result<Vec<u8>, SeqDbError> {
    if super::x_bn_mapped() {
        block.x_nucleotide_blastna_mapped(0)
    } else {
        block.nucleotide_blastna(0)
    }
}

/// A reader error outside the preliminary search as the command line reports it
/// (`algorithm::common::database::database_error`).
fn database_error(error: SeqDbError) -> anyhow::Error {
    shared::database_error(error, PROGRAM)
}

/// A reader error while the preliminary search reads the sequences
/// (`algorithm::common::database::search_read_error`).
fn search_read_error(error: SeqDbError) -> anyhow::Error {
    shared::search_read_error(error, PROGRAM)
}

/// The database of a `-db` search: the opened database, its name as NCBI keeps it
/// (`BlastnArgs::database_name`: the `-db` value without white space at its ends) and its
/// numbers as the search uses them.
pub(super) struct BlastnDatabase {
    pub(crate) db: SeqDb,
    name: String,
    metadata: SubjectMetadata,
    /// The taxonomy database of this run, opened at the first name asked for (outfmt 15 and 16,
    /// and the batches' prefetch of outfmt 0, 5, 6, 7 and 17).
    taxdb: shared::TaxdbRun,
    /// The run's record of its loads (`shared::TaxdbRun::loaded`) as it was before the last batch
    /// that `note_loaded` recorded: the loads before the reports that the loop writes next
    /// (LOSAT_X_RPTSTREAM, which records each batch before it writes the batch's reports; empty
    /// without it, where the loop records nothing), for `reported_subjects`.
    loaded_before_batch: std::sync::Mutex<crate::report::db_subject::LoadsBefore>,
    /// Whether the database has parsed Seq-ids (`shared::check_seqids`), whose subjects NCBI's
    /// prefetch can fail at (`batch_prefetch`).
    parsed_seqids: bool,
}

/// The length and the ambiguity runs of the database sequences that one report call wrote
/// letters of (`BlastnDatabase::sam_letters`), by ordinal id.
#[derive(Default)]
pub(super) struct SamLetterCache(RefCell<HashMap<u32, NucleotideAmbiguity>>);

impl BlastnDatabase {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2645-2668
    // (`Blast_HSPListReevaluateUngapped`)
    /// Calls `use_letters` with the letters of database sequence `ordinal` in BLASTNA, without
    /// the sentinel bytes: the subject that an ungapped (`-ungapped`) search re-evaluates its
    /// HSPs on in the preliminary stage (NCBI gets the sequence again, unpacked, for a sequence
    /// whose HSP list is not empty), with `-db_hard_mask` (the algorithm `hard_mask`) its masked
    /// letters 'B' (`mask_hard_letters`; NCBI's fetch in BLASTNA applies the masks). The letters
    /// stay in the kernel's buffer with its sentinel bytes (no copy). A read error is the
    /// preliminary search's.
    pub(super) fn with_blastna<T>(
        &self,
        ordinal: u32,
        hard_mask: Option<i32>,
        use_letters: impl FnOnce(&[u8]) -> T,
    ) -> Result<T> {
        let block = x_read_sequence(&self.db, ordinal).map_err(search_read_error)?;
        let length = block.nucleotide_len(0).map_err(search_read_error)?;
        let mut blastna = x_blastna(&block).map_err(search_read_error)?;
        if let Some(algo_id) = hard_mask {
            let mut masks = MaskBlock::default();
            let stored =
                self.sequence_masks(ordinal, algo_id, &mut masks, shared::search_read_error_of)?;
            mask_hard_letters(&mut blastna, stored, &[(0, length)]);
        }
        // The kernel's BLASTNA has a sentinel byte at each end.
        let result = use_letters(&blastna[1..=length]);
        x_done_bytes(blastna);
        x_done_block(block);
        Ok(result)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:185
    /// The OIDs of the database's ordinal ids, which name the reported sequences without
    /// parsed Seq-ids (`gnl|BL_ORD_ID|n`).
    pub(super) fn oid_map(&self) -> crate::seqdb::OidMap<'_> {
        self.db.oid_map()
    }

    /// The length in bases of each of the database sequences `ordinals` (outfmt 17: the
    /// Bioseq's length for `@SQ` and the clips of a record), from the index and the last byte
    /// of each sequence's packed bases; no letters are read.
    pub(super) fn sam_lengths(&self, ordinals: &[u32]) -> Result<Vec<usize>> {
        ordinals
            .iter()
            .map(|&ordinal| self.db.nucleotide_length(ordinal).map_err(database_error))
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-851
    /// The taxonomy of the tabular fields (`staxid`, `staxids`;
    /// `report::tabular::TabularTaxonomy::of_deflines`) of the reported sequences `ordinals`, by
    /// rank: from the definition lines of a sequence with parsed Seq-ids (`parsed[rank]`), else
    /// from its stored ones (`SeqDb::defline_set`, which keeps whether a line has a `taxid`), for
    /// the ranks that `hit_lists` name and that a report loads (`loads`); the others are left
    /// empty: NCBI reads no definition line of a subject whose load throws (its rows keep the
    /// taxonomy of the rows before) nor of one that no report reaches, and looks none of their
    /// taxids up (review 1 of line DBID pfoff-2, H1 and M1). A header that does not decode is
    /// refused before the report.
    pub(super) fn tabular_taxonomy(
        &self,
        ordinals: &[u32],
        parsed: &[Option<Arc<ParsedSubject>>],
        hit_lists: &[Option<BlastnHitList>],
        loads: &HitLoads,
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        use crate::report::tabular::TabularTaxonomy;
        let mut taxonomy = vec![TabularTaxonomy::default(); ordinals.len()];
        let mut read = vec![false; ordinals.len()];
        for hsp_list in hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| &hit_list.hsplist_array)
        {
            let Some(rank) = hsp_list.hsps.first().map(|hsp| hsp.s_idx as usize) else {
                continue;
            };
            if rank >= ordinals.len()
                || !loads.loaded(rank)
                || std::mem::replace(&mut read[rank], true)
            {
                continue;
            }
            taxonomy[rank] = match parsed.get(rank).and_then(Option::as_deref) {
                Some(subject) => TabularTaxonomy::of_deflines(subject.deflines()),
                None => {
                    let set = self
                        .db
                        .defline_set(ordinals[rank])
                        .map_err(database_error)?;
                    TabularTaxonomy::of_deflines(&set.deflines)
                }
            };
        }
        Ok(taxonomy)
    }

    /// The names of the taxonomy database that the rows of `spec` look up for the subjects
    /// `taxonomy` (`shared::tabular_names`).
    pub(super) fn tabular_names(
        &self,
        spec: &crate::report::tabular::TabularSpec,
        taxonomy: &[crate::report::tabular::TabularTaxonomy],
    ) -> Result<Option<crate::report::tabular::TabularNames>> {
        shared::tabular_names(&self.db, &self.taxdb, spec, taxonomy, PROGRAM)
    }

    /// The open of the taxonomy database by the rows of a query batch (`shared::TaxdbRun::open_for_rows`).
    pub(super) fn open_for_tabular_rows(&self) -> Result<()> {
        self.taxdb.open_for_rows(&self.db, PROGRAM)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:222-234
    /// The IUPAC letters `[from, to)` of database sequence `ordinal` (outfmt 17 with `SQ`):
    /// only the bytes of the range are read, and the sequence's length and ambiguity runs are
    /// read once per call and kept in `cache` (NCBI's `CSeqVector` also keeps them).
    pub(super) fn sam_letters(
        &self,
        cache: &SamLetterCache,
        ordinal: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>> {
        let mut entries = cache.0.borrow_mut();
        let ambiguity = match entries.entry(ordinal) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => entry.insert(
                self.db
                    .nucleotide_ambiguity(ordinal)
                    .map_err(database_error)?,
            ),
        };
        self.db
            .nucleotide_iupac_range_with(ordinal, from, to, ambiguity)
            .map_err(database_error)
    }

    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(super) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    /// The database of a call of the run after other calls whose loads its record does not hold
    /// where `unrecorded` (`shared::TaxdbRun::note_calls_before`: the call's first query is not the
    /// run's first, `first_query_ordinal` above 0).
    pub(super) fn after_calls(mut self, unrecorded: bool) -> Self {
        self.taxdb.note_calls_before(unrecorded);
        self
    }

    /// Whether NCBI's formatter keeps its scope across query batches (`shared::TaxdbRun::
    /// set_keeps_scope`: BLAST XML, outfmt 5), and whether XML2/JSON taxonomy reads belong to a
    /// secondary format (`shared::TaxdbRun::set_separate_xml2`).
    pub(super) fn keeping_scope(mut self, keeps: bool, separate_xml2: bool) -> Self {
        self.taxdb.set_keeps_scope(keeps);
        self.taxdb.set_separate_xml2(separate_xml2);
        self
    }

    /// Mirrors `CBlastFormat::ResetScopeHistory` at the end of a query batch
    /// (`shared::TaxdbRun::reset_scope`).
    pub(super) fn reset_scope(&self) {
        self.taxdb.reset_scope();
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:182-187
    /// Opens the database of `-db name` (`algorithm::common::database::open_database`)
    /// with the command line's errors. The total length is NCBI's `Int8`, which the engine holds
    /// as such (`SubjectMetadata::db_len_total`) on every target.
    pub(super) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        Self::open_restricted(name, settings, None, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line (`shared::restrict_database`, whose warnings go to `diagnostics`), and the id of the
    /// database mask `mask` checked as NCBI checks it there (`DbMaskOption::validate_at_open`).
    pub(super) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        list: Option<&crate::seqdb::DbList>,
        mask: Option<&DbMaskOption>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        // A database with parsed Seq-ids is read (line DBID, ID2b).
        let parsed_seqids =
            shared::restrict_database(&mut db, name, list, mask, diagnostics, settings, PROGRAM)?;
        // BLASTN has no Gumbel block: the setup of a batch with a valid context divides by the
        // statistics' number of sequences, a count of 0 LOSAT's refusal there
        // (`context_statistics`, `shared::check_setup_count`).
        // The setup's total is not negative (`DbStatistics::setup_total`).
        let (db_len_total, num_seqs) = db.statistics_numbers();
        let metadata = SubjectMetadata {
            db_len_total,
            db_num_seqs: num_seqs,
        };
        Ok(Self {
            db,
            name: name.to_string(),
            metadata,
            taxdb: shared::TaxdbRun::default(),
            loaded_before_batch: Default::default(),
            parsed_seqids,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:643-651
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_dbindex.cpp:664-668
    /// The volumes that NCBI's loader of a new-style MegaBLAST index finds for the `-db` value
    /// `indexname` (`EnumerateDbVolumes`: `FindVolumePaths` of `vols` with `true`, on the database opened
    /// again by that name, as given, on the search path of this one), each with its number of
    /// sequences (`GetNumOIDs` of the `CSeqDB` of the volume, the volume's in this database). Where that
    /// name does not open (a `-db` value with white space other than spaces at its ends, which
    /// the search trims and the loader keeps), the database reader's exception, or `None` where
    /// LOSAT does not know its place or the name opens other volumes than this database's.
    pub(super) fn index_volumes(
        &self,
        indexname: &str,
    ) -> Result<Vec<(Vec<u8>, u32)>, Option<NcbiException>> {
        let path = self.db.search_path().ok_or(None)?;
        let (paths, _) = crate::seqdb::SeqDb::volume_and_alias_paths(
            indexname.as_bytes(),
            SeqType::Nucleotide,
            path,
        )
        .map_err(|error| {
            if let Some(exception) = error.ncbi_exception() {
                return Some(exception.clone());
            }
            let text = error.ncbi_message()?;
            let site = error.throw_site()?;
            Some(NcbiException::new(site, text))
        })?;
        paths
            .into_iter()
            .map(|vol_name| {
                let volume = self
                    .db
                    .volumes()
                    .iter()
                    .find(|volume| {
                        volume.base() == &*crate::utils::os_bytes::path_of_bytes(&vol_name)
                    })
                    .ok_or(None)?;
                Ok((vol_name, volume.num_oids()))
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:453-462
    /// `SeqDB_ResolveDbPath` of `file` on this database's search path: the path of the first
    /// directory that holds the regular file `file` (an absolute `file` as it is), `None` for
    /// NCBI's empty string. The path is the bytes of NCBI's `std::string`.
    pub(super) fn resolve_db_path(&self, file: &[u8]) -> Option<Vec<u8>> {
        match self.db.search_path() {
            Some(search) => search.resolve_file(file),
            None => Some(file.to_vec())
                .filter(|_| crate::utils::os_bytes::path_of_bytes(file).is_file()),
        }
    }

    /// Whether NCBI's prefetch of a batch can fail at a subject (a database with parsed
    /// Seq-ids, `batch_prefetch`).
    pub(super) fn parsed_seqids(&self) -> bool {
        self.parsed_seqids
    }

    /// The Blast-def-line taxid of each reported sequence by rank (`ordinals`: the ordinal id of
    /// each rank, `ReportedSubjects::ordinals`; `parsed`: the sequences with parsed Seq-ids,
    /// `ReportedSubjects::parsed`) and what the taxonomy database says about it, for the XML2
    /// and JSON formats (`algorithm::common::database::read_reported_taxonomy`); the sequences
    /// of `hit_lists` with HSPs are the reported ones, of which those that the reports load
    /// (`loads`, `ReportedSubjects::loads`: the reads of a subject whose `Hit` throws, and of one
    /// shown only after it, are never made, so the texts of opening the taxonomy database and of a
    /// volume's first lookup of a taxid go to the loads that make them; review 1 of line DBID
    /// pfoff-1, H1). `(spans, gapped)`: the queries of each batch of the call (the whole call as
    /// one where empty) and whether the search is gapped, for NCBI's loads by batch
    /// (`shared::XmlBatchLoads`, `batches_prefetch`; asked only where they decide the texts).
    #[allow(clippy::type_complexity)]
    pub(super) fn read_taxids(
        &self,
        ordinals: &[u32],
        parsed: &[Option<Arc<ParsedSubject>>],
        hit_lists: &[Option<BlastnHitList>],
        loads: &HitLoads,
        (spans, gapped, format): (&[std::ops::Range<usize>], bool, shared::PrefetchFormat),
    ) -> Result<(
        Vec<i32>,
        Vec<crate::report::blast_xml::HitTaxonomy>,
        Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
        shared::TaxdbNotes,
    )> {
        let reported: std::collections::BTreeSet<usize> = hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| &hit_list.hsplist_array)
            .filter_map(|hsp_list| hsp_list.hsps.first())
            .map(|hsp| hsp.s_idx as usize)
            .filter(|&rank| loads.loaded(rank))
            .collect();
        // Each query's `Hit`s, the subjects of its HSP lists in order.
        let batches = || {
            shared::XmlBatchLoads::of(
                hit_lists.len(),
                spans,
                &batches_prefetch(format, hit_lists, 0, spans, gapped),
                loads,
                |query| {
                    hit_lists[query]
                        .iter()
                        .flat_map(|hit_list| &hit_list.hsplist_array)
                        .filter_map(|hsp_list| hsp_list.hsps.first())
                        .map(|hsp| hsp.s_idx as usize)
                        .collect()
                },
            )
        };
        shared::read_reported_taxonomy(
            &self.db,
            &self.taxdb,
            ordinals,
            &reported,
            parsed,
            &batches,
            PROGRAM,
        )
    }

    /// Whether `taxdb.bti` is on the search path, asked once for the run
    /// (`shared::TaxdbRun::found`): the batches' prefetch reads the taxonomy database only then
    /// (`prefetch_messages`).
    pub(super) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
    }

    /// NCBI's messages about the taxonomy database when the query batch whose final hit lists
    /// are `hit_lists` (after their HSPs of score 0 were dropped, still naming the sequences by
    /// ordinal id) prefetches the subjects it reports (`shared::prefetch_messages`): the
    /// sequences with HSPs in them that the reports load (`loads`, by ordinal id), the ones
    /// `read_taxids` reads, with the definition lines of those with parsed Seq-ids from `titles`
    /// (`read_titles`). Empty, and nothing read, without `taxdb.bti` on the search path. A
    /// taxonomy database that LOSAT cannot read as NCBI does, and a header whose taxid does not
    /// decode, are refused.
    pub(super) fn prefetch_messages(
        &self,
        titles: &DatabaseTitles,
        hit_lists: &[Option<BlastnHitList>],
        loads: &HitLoads,
    ) -> Result<Vec<u8>> {
        let subjects = hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| {
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                &hit_list.hsplist_array[..count]
            })
            .filter(|hsp_list| !hsp_list.hsps.is_empty())
            .filter(|hsp_list| loads.reached(hsp_list.oid as usize))
            .map(|hsp_list| {
                let parsed = match titles.get(&hsp_list.oid) {
                    Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
                    _ => None,
                };
                (hsp_list.oid, parsed)
            });
        shared::prefetch_messages(&self.db, &self.taxdb, subjects, PROGRAM)
    }

    /// [`Self::prefetch_messages`] by subject (`shared::subject_read_texts`; `loads`: as there): the text of the
    /// first read of each subject that writes one, with its ordinal id, where a query's report
    /// reads its subjects one by one: in the order of the hit lists, or, for the pairwise report
    /// (`pairwise`: how it reads them), those that it shows in the order that it reads them
    /// (`pairwise_read_order`, `PairwiseReads::subjects`).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn subject_read_texts(
        &self,
        titles: &DatabaseTitles,
        hit_lists: &[Option<BlastnHitList>],
        rows_open: bool,
        pairwise: Option<crate::report::pairwise::PairwiseReading>,
        placements: &crate::blastinput::seq_range::Placements,
        first_query: usize,
        ungapped: bool,
        loads: &HitLoads,
    ) -> Result<crate::report::query_warnings::QueryReads> {
        let parsed = |oid: u32| match titles.get(&oid) {
            Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
            _ => None,
        };
        let ordinals: Vec<u32> = match pairwise {
            Some(reading) => pairwise_read_order(
                hit_lists,
                |query| placements.offset(first_query + query),
                ungapped,
                reading,
                |oid| titled_parsed(titles, oid as usize),
            )
            .subjects()
            .into_iter()
            .map(|oid| oid as u32)
            .collect(),
            None => hit_lists
                .iter()
                .flat_map(shown_lists)
                .filter(|hsp_list| !hsp_list.hsps.is_empty())
                .map(|hsp_list| hsp_list.oid)
                .collect(),
        };
        let subjects = ordinals
            .into_iter()
            .filter(|&oid| loads.reached(oid as usize))
            .map(|oid| (oid, oid, parsed(oid)));
        shared::subject_read_texts(&self.db, &self.taxdb, subjects, rows_open, PROGRAM)
    }

    /// The XML2 and JSON formats, a batch of `hit_lists` (by query) whose prefetch is off: each
    /// query's report's `Hit` loads (`shared::scope_utf8_hit_reads`, the subjects in the hit
    /// list's order), into `reads` by query.
    fn scope_utf8_hit_reads(
        &self,
        titles: &DatabaseTitles,
        hit_lists: &[Option<BlastnHitList>],
        after: Option<crate::report::db_subject::RunThrow>,
        reads: &mut [crate::report::query_warnings::QueryReads],
    ) {
        if !self.db.scope_utf8_seen() {
            return;
        }
        for (index, hit_list) in hit_lists.iter().enumerate() {
            let subjects = hit_list.iter().flat_map(|hit_list| {
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                hit_list.hsplist_array[..count]
                    .iter()
                    .filter(|hsp_list| !hsp_list.hsps.is_empty())
                    .map(|hsp_list| {
                        let parsed = match titles.get(&hsp_list.oid) {
                            Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
                            _ => None,
                        };
                        (hsp_list.oid, hsp_list.oid, parsed)
                    })
            });
            let query =
                shared::scope_utf8_hit_reads(&self.db, &self.taxdb, after, subjects, PROGRAM);
            if !query.is_empty() {
                reads[index] = query;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:307-308
    // A gapped search makes a Seq-align of each HSP; an ungapped one makes one Seq-align of a
    // subject's HSPs, a Dense-diag each (`shared::subject_seq_aligns`), whose strand is its
    // first one's: the subject's first HSP (probe of line DBID pfwarn, `R/pfwarn` W4). BLASTN's
    // Seq-aligns are on the subject's plus strand.
    /// NCBI's prefetch of the query batch whose final hit lists are `hit_lists` on standard error
    /// (`shared::batch_prefetch`): the taxonomy texts of [`BlastnDatabase::prefetch_messages`],
    /// then the warning where the prefetch fails, at a subject whose lookup fails or one without
    /// a Bioseq that a Seq-align on the query's minus strand flips (`gapped`: each HSP is a
    /// Seq-align, else a subject's first HSP decides).
    /// `before`: the batches of the call before it that the loop still holds (none where it wrote
    /// and dropped them, LOSAT_X_RPTSTREAM, which records each batch instead: `note_loaded`), as
    /// NCBI's data loader loaded them (`LoadedBefore`), read only where the batch stops; the outfmt
    /// 5 reports before the batch are the queries' (`shared::PrefetchQueries::after_xml_reports`).
    pub(super) fn batch_prefetch(
        &self,
        titles: &DatabaseTitles,
        before: LoadedBefore<'_>,
        hit_lists: &[Option<BlastnHitList>],
        gapped: bool,
        queries: &shared::PrefetchQueries<'_>,
    ) -> Result<crate::report::diag_log::Posts> {
        shared::batch_prefetch(
            &self.db,
            &self.taxdb,
            Self::loaded_subjects(titles, before, gapped),
            Self::prefetch_subjects(titles, hit_lists, gapped, SubjectFetch::Every),
            PROGRAM,
            queries,
        )
    }

    /// NCBI's prefetch of the query batch of `queries` whose final hit lists are `hit_lists`, for
    /// the XML2 and JSON formats, which write the taxonomy texts at each `Hit`: where
    /// `PRE_FETCH_SEQS_LIMIT` turns it off (`shared::prefetch_off`), the `Hit` loads of each
    /// query's report (`scope_utf8_hit_reads`, into `reads`, by query of the batch) and `None`;
    /// where it is on and loads a sequence (a batch with Seq-aligns), its stop
    /// (`shared::held_batch_prefetch_stop`): after the lookups of the reports of the batches
    /// before, whose final hit lists are `before` (by ordinal id, as `hit_lists`; `spans`: the
    /// queries of each of those batches, numbered from 0, as BLASTN holds them all; the reports
    /// read every reported subject that they load: `hit_loads` of `batches_prefetch`, the reports'
    /// own answer, which leaves out a subject whose `Hit` throws and one shown only after it, with
    /// `after`; review 1 of line DBID pfoff-1, H1; review 2, N-M2), which BLASTN writes
    /// after the stop (the XML formats are not streamed, `x_rpt_stream`). `Some(Ok)`: the warning
    /// of a UTF8String header that a load of a batch that does not stop posts
    /// (`shared::batch_prefetch_stop`), for the batch's prefetch posts.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn xml2_prefetch(
        &self,
        titles: &DatabaseTitles,
        before: &[Option<BlastnHitList>],
        spans: &[std::ops::Range<usize>],
        (after, format): (
            crate::report::db_subject::FormatThrow,
            shared::PrefetchFormat,
        ),
        hit_lists: &[Option<BlastnHitList>],
        gapped: bool,
        queries: &shared::PrefetchQueries<'_>,
        reads: &mut [crate::report::query_warnings::QueryReads],
    ) -> Option<Result<crate::report::diag_log::Posts>> {
        let seq_aligns = seq_aligns_of(hit_lists, gapped);
        if shared::prefetch_off(format, seq_aligns) {
            self.scope_utf8_hit_reads(titles, hit_lists, after.run(), reads);
            return None;
        }
        (seq_aligns > 0).then(|| {
            let prefetch = batches_prefetch(format, before, 0, spans, gapped);
            let no_placements = crate::blastinput::seq_range::Placements::default();
            let loads = hit_loads(
                before,
                &prefetch,
                (after, None),
                &no_placements,
                !gapped,
                self.taxdb.calls_before(),
                |ordinal| titled_parsed(titles, ordinal),
            );
            shared::held_batch_prefetch_stop(
                &self.db,
                &self.taxdb,
                Self::prefetch_subjects(titles, before, gapped, SubjectFetch::Every)
                    .filter(|subject| loads.loaded(subject.ordinal as usize)),
                std::iter::empty(),
                Self::prefetch_subjects(titles, hit_lists, gapped, SubjectFetch::Every),
                PROGRAM,
                queries,
            )
            .map(crate::report::diag_log::Posts::of_lines)
        })
    }

    /// The run's record of the sequences that NCBI's BLAST database data loader loaded for a query
    /// batch whose final hit lists are `hit_lists` (`shared::TaxdbRun::note_loaded`;
    /// `loaded_subjects` of the batch alone, `rows`: the first format's), for a later batch's
    /// stop: where the loop writes each batch's reports and drops its lists (LOSAT_X_RPTSTREAM),
    /// so that `batch_prefetch`'s `before` no longer holds them. The stream never writes outfmt 5
    /// (`x_rpt_stream`), so no report of the batch is cut (`shared::XmlCuts::NONE`). `start`: the
    /// batch's first query in the call, its number in the record (`shared::TaxdbRun::note_loaded`,
    /// as `reported_subjects` and `QueryFetches::batch` number the batches).
    pub(super) fn note_loaded(
        &self,
        titles: &DatabaseTitles,
        hit_lists: &[Option<BlastnHitList>],
        start: usize,
        gapped: bool,
        (format, rows): (shared::PrefetchFormat, SubjectFetch),
    ) {
        if self.parsed_seqids() {
            *self
                .loaded_before_batch
                .lock()
                .expect("the loads before the batch") = self.taxdb.loads_before(Some(start));
            let batch = [0..hit_lists.len()];
            self.taxdb.note_loaded(&self.db, start, || {
                Self::loaded_subjects(
                    titles,
                    LoadedBefore {
                        hits: hit_lists,
                        spans: &batch,
                        format,
                        rows,
                        cuts: shared::XmlCuts::NONE,
                    },
                    gapped,
                )
            });
        }
    }

    /// The sequences that NCBI's BLAST database data loader loaded for the batches `before`
    /// (`shared::LoadedBefore::subjects`: the subjects of each batch's Seq-aligns where its
    /// prefetch is on, else those that the rows of the reports NCBI writes read,
    /// `shared::BatchLoads::reported_rows`; review 3 of line DBID loadord, H2; review 6, H2).
    fn loaded_subjects<'s>(
        titles: &'s DatabaseTitles,
        before: LoadedBefore<'s>,
        gapped: bool,
    ) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
        before.subjects(
            |hit_lists, span| hit_lists.get(span),
            move |_, hit_lists| seq_aligns_of(hit_lists, gapped),
            // BLASTN's prefetched set is its reported set.
            move |_, hit_lists, loads| {
                Self::prefetch_subjects(titles, hit_lists, gapped, loads.reported_rows())
            },
        )
    }

    /// The sequences that the prefetch of the query batch whose final hit lists are `hit_lists`
    /// reads (`batch_prefetch`), or, as its reports read them where the prefetch is off, those
    /// whose rows read them (`rows`, by each HSP's `num_ident`; BLASTN's formatter does not fetch
    /// for the identities, so its rows read every subject or none).
    fn prefetch_subjects<'s>(
        titles: &'s DatabaseTitles,
        hit_lists: &'s [Option<BlastnHitList>],
        gapped: bool,
        rows: SubjectFetch,
    ) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
        hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| {
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                &hit_list.hsplist_array[..count]
            })
            .filter(move |hsp_list| hsp_list.hsps.iter().any(|hsp| rows.reads(hsp.num_ident)))
            .map(move |hsp_list| {
                let parsed = match titles.get(&hsp_list.oid) {
                    Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
                    _ => None,
                };
                let flipped = if gapped {
                    hsp_list.hsps.iter().any(|hsp| hsp.query_frame < 0)
                } else {
                    hsp_list.hsps[0].query_frame < 0
                };
                shared::PrefetchSubject {
                    ordinal: hsp_list.oid,
                    parsed,
                    flipped,
                }
            })
    }

    /// The `-db` value as NCBI keeps it (outfmt 7's `# Database:` line).
    pub(super) fn name(&self) -> &str {
        &self.name
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(super) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
    }

    /// Whether the first query batch with a valid context fails with NCBI's `Database memory map
    /// file error` (`SeqDb::first_oid_read_fails`: no OID, and statistics that read OID 0).
    pub(super) fn first_oid_read_fails(&self) -> bool {
        self.db.first_oid_read_fails()
    }

    /// The error of the preliminary search's first read of the sequences after a batch's set-up
    /// checks, where a sequence file cannot be opened (`shared::sequence_file_read_error`).
    pub(super) fn sequence_file_read_error(&self) -> Option<anyhow::Error> {
        shared::sequence_file_read_error(&self.db, PROGRAM)
    }

    /// `BlastSeqSrcGetAvgSeqLen` (`SeqDb::average_length`).
    pub(super) fn average_length(&self) -> i32 {
        self.db.average_length()
    }

    /// The numbers the statistics read, with NCBI's types (`SeqDb::statistics`); the BLAST
    /// archive's database length reads them too.
    pub(super) fn statistics(&self) -> crate::seqdb::DbStatistics {
        self.db.statistics()
    }

    /// The taxids of the search's list (`SeqDb::list_taxids`).
    pub(super) fn list_taxids(&self) -> &[i32] {
        self.db.list_taxids()
    }

    /// The GIs of the search's list (`SeqDb::list_gis`).
    pub(super) fn list_gis(&self) -> &[i64] {
        self.db.list_gis()
    }

    /// The number of sequences the search reads: the included ones (`SeqDb::num_seqs`), which
    /// the schedules count; the statistics' number (`metadata`) can be another.
    pub(super) fn num_subjects(&self) -> usize {
        self.db.num_seqs() as usize
    }

    /// Whether NCBI's traceback sorts the results by e-value at its end
    /// (`DbStatistics::sorts_results`).
    pub(super) fn sorts_results(&self) -> bool {
        self.db.statistics().sorts_results()
    }

    /// The number of sequences and the total length that the statistics read
    /// (`SeqDb::statistics_numbers`): the real numbers, also with `-dbsize`.
    pub(super) fn metadata(&self) -> &SubjectMetadata {
        &self.metadata
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`).
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// `pairwise_database` of a search masked with `mask` (`None`: not masked): each name's
    /// algorithm, which the `Masked using:` lines show (`shared::pairwise_database_masked`).
    pub(super) fn pairwise_database_masked(&self, mask: Option<&ResolvedMask>) -> PairwiseDatabase {
        shared::pairwise_database_masked(&self.db, mask.map(|mask| mask.shown.as_slice()))
    }

    /// The formatter's resolution of `option` (`DbMaskOption::resolve`; its warning goes to
    /// `warning`), and LOSAT's refusals of the masked search (`ResolvedMask::check_losat_limits`).
    pub(super) fn resolve_mask_option(
        &self,
        option: &DbMaskOption,
        warning: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
    ) -> Result<Option<ResolvedMask>> {
        option.resolve(&self.db, &self.name, PROGRAM, warning)
    }

    /// LOSAT's refusals of the search masked with `mask` (`ResolvedMask::check_losat_limits`),
    /// whose scan reads words of `scanned_word` letters.
    pub(super) fn check_mask_limits(&self, mask: &ResolvedMask, scanned_word: usize) -> Result<()> {
        mask.check_losat_limits(
            &self.db,
            PROGRAM,
            crate::blastinput::db_mask::HardScan::Nucleotide { scanned_word },
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:197-218 and
    // c++/include/objtools/blast/seqdb_reader/seqdb.hpp:286-306 (cited at
    // `core::subject_split::complement_masks`): the complement of the masks in the order they
    // are stored, nothing sorted or merged (decision D3 of M4a-2).
    /// The ranges (`seq_ranges`) of database sequence `ordinal` of length `length`, read in the
    /// block `block` of the preliminary stage (`block_subject`), in a search with `-db_soft_mask`
    /// or `-db_hard_mask` (the algorithm `algo_id`): the complement of its masks, which `masks`
    /// reads with the masks of the whole block (`SeqDb::read_mask_block`, once for all its
    /// sequences), into `ranges`, which it replaces. They are the subject's soft ranges or, with
    /// `-db_hard_mask`, its hard ranges (blast_engine.c:165-176).
    pub(super) fn seq_ranges_into(
        &self,
        block: BlockRange,
        ordinal: u32,
        algo_id: i32,
        length: usize,
        masks: &mut MaskBlock,
        ranges: &mut Vec<(i32, i32)>,
    ) -> Result<()> {
        let stored =
            self.stored_masks(block, ordinal, algo_id, masks, shared::search_read_error_of)?;
        crate::core::subject_split::complement_masks(stored, length, ranges);
        Ok(())
    }

    /// The stored masks of the database sequences `ordinals` for the algorithm `algo_id`, read
    /// in blocks of at most `MASK_BLOCK_SEQS` sequences (`stored_masks`), each sequence's given to
    /// `each` with its ordinal id: the masked letters of a MegaBLAST index volume's search with
    /// `-db_hard_mask` (`mb_index_seeds::IndexSeeds::volume_masks`). A sequence whose masks do not
    /// read is left out; its own search reads them again and reports the error.
    pub(super) fn masks_of(
        &self,
        ordinals: std::ops::Range<u32>,
        algo_id: i32,
        mut each: impl FnMut(u32, &[(u32, u32)]),
    ) {
        let mut masks = MaskBlock::default();
        let mut ordinal = ordinals.start;
        while ordinal < ordinals.end {
            let end = ordinals.end.min(ordinal.saturating_add(MASK_BLOCK_SEQS));
            let Some(block) = self.db.x_block_range(ordinal, end) else {
                break;
            };
            if block.is_empty() {
                break;
            }
            for at in ordinal..ordinal + block.len() {
                if let Ok(stored) =
                    self.stored_masks(block, at, algo_id, &mut masks, shared::search_read_error_of)
                {
                    each(at, stored);
                }
            }
            ordinal += block.len();
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqinfosrc_seqdb.cpp:121-133
    /// The stored masks of database sequence `ordinal` for the algorithm `algo_id`, in the order
    /// they are stored (`[start, end)`), which the reports of a masked search show in lowercase
    /// where they meet an HSP (`pairwise::Rows`); `masks` is the reports' block buffer.
    pub(super) fn display_masks(
        &self,
        ordinal: u32,
        algo_id: i32,
        masks: &mut MaskBlock,
    ) -> Result<Vec<crate::utils::dust::MaskedInterval>> {
        let Some(block) = self.db.x_block_range(ordinal, ordinal.saturating_add(1)) else {
            anyhow::bail!(
                "LOSAT internal error: the reported sequence {ordinal} is not in the database"
            );
        };
        Ok(self
            .stored_masks(block, ordinal, algo_id, masks, shared::database_error_of)?
            .iter()
            .map(|&(start, end)| {
                crate::utils::dust::MaskedInterval::new(start as usize, end as usize)
            })
            .collect())
    }

    /// The stored masks of database sequence `ordinal` for the algorithm `algo_id`, read alone
    /// (a block of one sequence in `masks`), with the reader's errors as `error` reports them
    /// (`stored_masks`): those of the letters of a hard-masked search's fetch.
    fn sequence_masks<'m>(
        &self,
        ordinal: u32,
        algo_id: i32,
        masks: &'m mut MaskBlock,
        error: fn(&SeqDbError, &str) -> anyhow::Error,
    ) -> Result<&'m [(u32, u32)]> {
        let Some(block) = self.db.x_block_range(ordinal, ordinal.saturating_add(1)) else {
            anyhow::bail!("LOSAT internal error: sequence {ordinal} is not in the database");
        };
        self.stored_masks(block, ordinal, algo_id, masks, error)
    }

    /// The stored masks of database sequence `ordinal` of the block `block` for the algorithm
    /// `algo_id` (`SeqDb::read_mask_block`, `MaskBlock::ranges`), with the reader's errors as
    /// `error` reports them: the preliminary search's (`shared::search_read_error_of`, as its
    /// reads of the sequences: NCBI reads the masks in the search's `GetSequence`), the reports'
    /// (`shared::database_error_of`). A volume without the algorithm is a refusal that the search
    /// made before it started (`SeqDb::check_search_masks`).
    fn stored_masks<'m>(
        &self,
        block: BlockRange,
        ordinal: u32,
        algo_id: i32,
        masks: &'m mut MaskBlock,
        error: fn(&SeqDbError, &str) -> anyhow::Error,
    ) -> Result<&'m [(u32, u32)]> {
        self.db.read_mask_block(block, algo_id, masks);
        match masks.ranges(block.first + (ordinal - block.first_ordinal)) {
            Ok(BlockMasks::Ranges(stored)) => Ok(stored),
            Ok(BlockMasks::NoMaskData) => Ok(&[]),
            Ok(BlockMasks::AlgorithmNotInVolume) => Err(anyhow::anyhow!(
                "LOSAT internal error: the volume of sequence {ordinal} has no mask algorithm {algo_id}"
            )),
            Err(read) => Err(error(read.error(), PROGRAM)),
        }
    }
}

/// The sequences of one read of masks in `BlastnDatabase::masks_of` (its blobs are held at once).
const MASK_BLOCK_SEQS: u32 = 1 << 16;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1518-1538
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1665-1673
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:590-606
/// BLASTNA 'B' (NCBI4na 14 mapped to BLASTNA).
const HARD_MASK_BLASTNA: u8 = 10;

/// With `-db_hard_mask`, NCBI's fetch of a sequence in BLASTNA (the traceback's, the ungapped
/// re-evaluation's) writes the mask letter NCBI4na 14 over every masked position after the
/// ambiguities, before the BLASTNA mapping: the letters are 'B' (BLASTNA 10) there (series NP,
/// M4a-2b decision H2). `blastna` has a sentinel byte before the first letter; `masks` are the
/// sequence's stored masks (`[start, end)`, sorted and disjoint: `SeqDb::check_hard_masks`
/// refuses the others, whose result NCBI's scan of sorted masks decides); `ranges` are the decoded
/// parts (`[from, to)`: the whole sequence, or the slices of the partial decode), and only their
/// letters are masked, as NCBI masks each slice.
fn mask_hard_letters(blastna: &mut [u8], masks: &[(u32, u32)], ranges: &[(usize, usize)]) {
    for &(begin, end) in ranges {
        let mut i = 0;
        while i < masks.len() && masks[i].1 as usize <= begin {
            i += 1;
        }
        while i < masks.len() && (masks[i].0 as usize) < end {
            let from = (masks[i].0 as usize).max(begin);
            let to = (masks[i].1 as usize).min(end);
            if from < to {
                blastna[1 + from..1 + to].fill(HARD_MASK_BLASTNA);
            }
            i += 1;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1405
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1805-1838
// The preliminary stage searches the stored 2-bit data as makeblastdb wrote it (its
// ambiguous bases already replaced by makeblastdb's random bases), in the layout of the
// subject set's packed letters (four bases per byte from the high bits down, the count of
// the last byte's bases in its two low bits; `encode_subject_ncbi2na_packed`).
/// Sequence `i` of a block as the preliminary stage reads it: its length and its stored
/// 2-bit data, and the block (`range`), whose masks a search with `-db_soft_mask` reads once for
/// all its sequences (`BlastnDatabase::seq_ranges_into`). The block is lent, not copied: every
/// sequence of every database search carries it, and a copy of its 24 bytes into each sequence's
/// reference cost megablast's default path 0.3 % of its instructions (ws, SILVA, ad-M4a-5).
fn block_subject<'a>(
    block: &'a SequenceBlock,
    range: &'a BlockRange,
    i: usize,
) -> Result<BlastnSubjectRef<'a>> {
    let (packed, length) = block.nucleotide_packed_len(i).map_err(search_read_error)?;
    Ok(BlastnSubjectRef {
        length,
        packed,
        letters: None,
        blastna: None,
        store: None,
        x_fenced: None,
        block: Some(range),
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
/// The preliminary stage of a database search: every sequence, in the order of its
/// ordinal id, through `search` (the per-subject kernel), and what it keeps of a sequence
/// (`Some`) through `merge` in that order (NCBI's loop with one thread; `merge` puts the
/// sequence's preliminary HSPs into the preliminary hit lists).
///
/// Without `parallel`: the blocks in order, one scratch (`new_scratch`) for the whole
/// database, as the one-thread loop of the `-subject` path.
///
/// With `parallel`: rounds of blocks. The blocks of a round are read and searched in
/// parallel with a scratch per worker (`map_init`), as the parallel `-subject` path
/// searches its subjects (the result of a subject does not depend on what the worker's
/// scratch searched before it, which that path relies on). The round's results are merged
/// in ordinal order and its blocks dropped before the next round is read: a worker holds
/// one block at a time, so the sequence memory is bounded by the threads and the block
/// size, not by the database. A block that cannot be read does not stop the other blocks
/// of its round; the round's first error in ordinal order is returned, the error of the
/// one-thread order (the blocks are cut alike).
///
/// EXPERIMENT (LOSAT_X_BNDBCUT): with more than one thread the blocks are cut by the
/// database size and the thread count, and a failed search returns the one-thread order's
/// error (`x_db_cut.rs`, exactness argument in its header).
///
/// EXPERIMENT (LOSAT_X_BNDBSTREAM): with more than one thread there are no rounds: the
/// threads claim the blocks in ordinal order and one of them hands the results over in that
/// order while the others search the next blocks, at most a window of blocks ahead
/// (`x_db_stream.rs`, exactness argument in its header).
#[allow(clippy::too_many_arguments)]
pub(super) fn search_blocks<T, R, I, K, M>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send + XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let blocks_per_thread = BLOCKS_PER_THREAD;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let blocks_per_thread = 1;
    // EXPERIMENT (LOSAT_X_BNDBCUT): `x_db_cut.rs`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_db_cut::x_bn_db_cut() {
        return x_db_cut::x_search_blocks_cut(
            database,
            pool,
            parallel,
            blocks_per_thread,
            new_scratch,
            search,
            merge,
        );
    }
    search_blocks_of(
        database,
        pool,
        parallel,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        search,
        merge,
    )
}

// No NCBI counterpart (scheduling only): the rounds of LOSAT's parallel stage; the merges in
// ordinal order are NCBI's one-thread order (blast_engine.c:1553-1554, cited above
// `search_blocks`).
/// The rounds of `search_blocks_of` with `parallel` (without LOSAT_X_BNDBSTREAM):
/// `blocks_per_thread` blocks of at most `block_bytes` per thread
/// searched with a scratch per worker (`map_init`), and the round's results merged in ordinal
/// order when the whole round is searched.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn search_rounds<T, R, I, K, M>(
    db: &SeqDb,
    pool: &SearchPool<'_>,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send + XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R),
{
    let round_blocks = pool.threads().saturating_mul(blocks_per_thread).max(1);
    let mut next_ordinal = 0u32;
    loop {
        let mut round: Vec<BlockRange> = Vec::with_capacity(round_blocks);
        while round.len() < round_blocks {
            let Some(range) = db.next_block(next_ordinal, block_bytes) else {
                break;
            };
            round.push(range);
            next_ordinal = range.end_ordinal();
        }
        if round.is_empty() {
            break;
        }
        let found: Vec<Result<Vec<(u32, R)>>> = pool.install(|| {
            round
                .par_iter()
                .map_init(
                    || XWorker::new(new_scratch()),
                    |worker, &range| -> Result<Vec<(u32, R)>> {
                        let block = x_read_block(db, range, &mut worker.bytes)
                            .map_err(search_read_error)?;
                        let mut kept = Vec::new();
                        for i in 0..block.len() {
                            let ordinal = range.first_ordinal + i as u32;
                            if let Some(result) = search(
                                &mut worker.scratch,
                                ordinal,
                                block_subject(&block, &range, i)?,
                            ) {
                                kept.push((ordinal, result));
                            }
                        }
                        x_done_kept_block(block, &mut worker.bytes);
                        Ok(kept)
                    },
                )
                .collect()
        });
        // The sequences that keep nothing change no hit list, so the order of the
        // merges is the one-thread order.
        for block_results in found {
            for (ordinal, result) in block_results? {
                merge(ordinal, result);
            }
        }
    }
    Ok(())
}

/// `search_blocks` with blocks of at most `block_bytes` and rounds of `blocks_per_thread`
/// blocks per thread (with LOSAT_X_BNDBSTREAM, the window of `x_db_stream.rs` instead of
/// rounds).
#[allow(clippy::too_many_arguments)]
fn search_blocks_of<T, R, I, K, M>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send + XWaitingBytes,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R) + Send,
{
    let db = &database.db;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parallel {
        // EXPERIMENT (LOSAT_X_BNDBSTREAM): `x_db_stream.rs`.
        if x_db_stream::x_bn_db_stream() {
            crate::utils::xreach::x_reach!("LOSAT_X_BNDBSTREAM", true);
            return x_db_stream::x_search_blocks_stream(
                db,
                pool,
                block_bytes,
                new_scratch,
                search,
                merge,
            );
        }
        return search_rounds(
            db,
            pool,
            block_bytes,
            blocks_per_thread,
            new_scratch,
            search,
            merge,
        );
    }
    // One thread: the blocks are searched one after another (the reference path).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_db_stream::x_bn_db_stream() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNDBSTREAM", false);
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, parallel, blocks_per_thread);

    let mut scratch = new_scratch();
    // The block buffer, kept from one block to the next (`x_read_block`).
    let mut bytes = Vec::new();
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block = x_read_block(db, range, &mut bytes).map_err(search_read_error)?;
        for i in 0..block.len() {
            let ordinal = range.first_ordinal + i as u32;
            if let Some(result) = search(&mut scratch, ordinal, block_subject(&block, &range, i)?) {
                merge(ordinal, result);
            }
        }
        x_done_kept_block(block, &mut bytes);
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:196-203
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:819-828
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1567-1573
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
// The sequence comes in BLASTNA with its ambiguity codes, as the whole sequence (NCBI's
// partial fetch of the HSP ranges with fences gives the same alignments; decision 87 of
// the campaign: the whole sequence first).
// The traceback is outside the handler of the reader's errors of the preliminary search
// (local_blast.cpp:240-251), so a sequence that cannot be read there is the application's
// `BLAST Database error` (`database_error`).
/// The traceback of a database search: the sequences of `listed` (their ordinal ids,
/// ascending, with what the preliminary stage kept of each) read again one at a time in
/// the worker that searches them (the stored 2-bit data and the BLASTNA letters with the
/// ambiguity codes, `SequenceBlock::nucleotide_blastna`, dropped after the sequence), and
/// what `traceback` returns of each, in the order of `listed`. With `parallel` the
/// sequences are searched in parallel, and each takes a scratch that no worker holds (or makes
/// one when none is idle) and gives it back when its traceback ends: the stage makes about one
/// scratch per worker (a few more when a worker waiting in a sequence's nested parallel stage
/// traces another sequence). Before, `map_init` made one per rayon task, and rayon cuts a list
/// into about twice as many tasks as threads and cuts a task again each time a thread steals
/// from it (megablast, 100 queries against 16S rRNA `-db`, 16 threads: 22,032 scratches for
/// 49,137 listed sequences, each with its offset pairs, about 100 KiB filled with zeros, and
/// the speculative and look-ahead alignment scratches its first traceback makes, one per
/// thread). A sequence's result does not depend on what its scratch traced before, as with the
/// tasks of `map_init`, which also traced runs of sequences that change from one execution to
/// the next.
///
/// With `-db_hard_mask` (the algorithm `hard_mask`) the masked letters of each decoded part are
/// 'B' (`mask_hard_letters`; NCBI's fetch with `check_oid_exclusion` applies the masks,
/// seqsrc_seqdb.cpp:324-353), the masks read with the sequence (`BlastnDatabase::sequence_masks`).
///
/// EXPERIMENT (LOSAT_X_BNTBRANGE): a sequence of at least `x_tb_range_min` bases is decoded in
/// the ranges of its kept HSPs only (`spans`: `[from, to)` on the sequence; `x_partial_ranges`),
/// with fences; when a traceback alignment meets one, the sequence is decoded whole and traced
/// again from `kept`.
#[allow(clippy::too_many_arguments)]
pub(super) fn traceback_listed<T, P, R, I, K, S>(
    database: &BlastnDatabase,
    hard_mask: Option<i32>,
    pool: &SearchPool<'_>,
    parallel: bool,
    listed: &[(u32, P)],
    new_scratch: &I,
    traceback: &K,
    spans: &S,
) -> Result<Vec<R>>
where
    T: Send,
    P: Sync,
    R: Send + std::fmt::Debug,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>, &P) -> Option<R> + Sync + Send,
    S: Fn(&P) -> Vec<(usize, usize)> + Sync + Send,
{
    let db = &database.db;
    let x_ranged = super::x_bn_tb_range();
    let x_min = super::x_tb_range_min();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-353
    // The masks are read before the letters, and outside the handler of the preliminary
    // search's read errors, as the letters are.
    let hard_masks = |ordinal: u32, algo_id: i32| -> Result<Vec<(u32, u32)>> {
        let mut masks = MaskBlock::default();
        Ok(database
            .sequence_masks(ordinal, algo_id, &mut masks, shared::database_error_of)?
            .to_vec())
    };
    let one = |scratch: &mut T, &(ordinal, ref kept): &(u32, P)| -> Result<Option<R>> {
        let masks = match hard_mask {
            Some(algo_id) => Some(hard_masks(ordinal, algo_id)?),
            None => None,
        };
        let stored = x_read_sequence(db, ordinal).map_err(database_error)?;
        let length = stored.nucleotide_len(0).map_err(database_error)?;
        let x_in_ranges = x_ranged && length >= x_min;
        if x_ranged {
            crate::utils::xreach::x_reach!("LOSAT_X_BNTBRANGE", x_in_ranges);
        }
        if x_in_ranges {
            let ranges = x_partial_ranges(spans(kept), length);
            let mut blastna = stored
                .x_nucleotide_blastna_ranges(0, &ranges)
                .map_err(database_error)?;
            if let Some(masks) = &masks {
                mask_hard_letters(&mut blastna, masks, &ranges);
            }
            let fenced = std::sync::atomic::AtomicBool::new(false);
            let subject = BlastnSubjectRef {
                length,
                packed: stored.nucleotide_packed(0),
                letters: None,
                blastna: Some(&blastna[1..=length]),
                store: None,
                x_fenced: Some(&fenced),
                block: None,
            };
            let result = traceback(scratch, ordinal, subject, kept);
            drop(blastna);
            let fenced = fenced.into_inner();
            x_tb_range_count(length, &ranges, fenced);
            if !fenced {
                x_done_block(stored);
                return Ok(result);
            }
        }
        let result = x_traceback_whole(
            &stored,
            length,
            masks.as_deref(),
            scratch,
            ordinal,
            kept,
            traceback,
        )?;
        x_done_block(stored);
        Ok(result)
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parallel {
        // No NCBI counterpart (scheduling only): the scratches that no worker holds.
        let scratches: std::sync::Mutex<Vec<T>> = std::sync::Mutex::new(Vec::new());
        let idle = || {
            scratches
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
        };
        let found: Vec<Result<Option<R>>> = pool.install(|| {
            listed
                .par_iter()
                .map(|entry| {
                    // The lock is released before a new scratch is made.
                    let held = idle().pop();
                    let mut scratch = held.unwrap_or_else(new_scratch);
                    let result = one(&mut scratch, entry);
                    idle().push(scratch);
                    result
                })
                .collect()
        });
        let mut results = Vec::with_capacity(found.len());
        for result in found {
            results.extend(result?);
        }
        return Ok(results);
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, parallel);

    let mut scratch = new_scratch();
    let mut results = Vec::new();
    for entry in listed {
        results.extend(one(&mut scratch, entry)?);
    }
    Ok(results)
}

/// The traceback of `traceback_listed` on the whole sequence of `stored`.
/// `hard_masks`: the masks of a search with `-db_hard_mask` (`mask_hard_letters`).
fn x_traceback_whole<T, P, R, K>(
    stored: &SequenceBlock,
    length: usize,
    hard_masks: Option<&[(u32, u32)]>,
    scratch: &mut T,
    ordinal: u32,
    kept: &P,
    traceback: &K,
) -> Result<Option<R>>
where
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>, &P) -> Option<R>,
{
    let mut blastna = x_blastna(stored).map_err(database_error)?;
    if let Some(masks) = hard_masks {
        mask_hard_letters(&mut blastna, masks, &[(0, length)]);
    }
    let subject = BlastnSubjectRef {
        length,
        packed: stored.nucleotide_packed(0),
        letters: None,
        // The kernel's BLASTNA has no sentinel bytes.
        blastna: Some(&blastna[1..=length]),
        store: None,
        x_fenced: None,
        block: None,
    };
    let result = traceback(scratch, ordinal, subject, kept);
    x_done_bytes(blastna);
    Ok(result)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seqsrc.c:467-519
// (blast_seqsrc.h:203-204: BLAST_SEQSRC_MINGAP 1024, BLAST_SEQSRC_OVERHANG 1024.) The ranges
// come from every HSP that the batch keeps of the sequence (`BLAST_SetupPartialFetching`,
// blast_traceback.c:1320-1370); the merged ranges do not depend on the order of equal starts.
/// EXPERIMENT (LOSAT_X_BNTBRANGE): NCBI's partial-fetch ranges of a sequence of `length` bases
/// whose HSPs cover `spans` (`[from, to)`): each widened by 1024 bases on both sides within the
/// sequence, sorted by start, and merged into the range before unless it starts more than 1024
/// bases after that one's end. Sorted and disjoint.
pub(super) fn x_partial_ranges(
    spans: impl IntoIterator<Item = (usize, usize)>,
    length: usize,
) -> Vec<(usize, usize)> {
    x_partial_ranges_padded(spans, length, x_tb_range_pad())
}

/// EXPERIMENT (`LOSAT_X_BNTBRANGEPAD=<bases>`, a test knob): NCBI's 1024 of `x_partial_ranges`
/// (both the widening and the gap), 16 to 1,048,576; a small one makes the traceback meet
/// fences, so that the tests reach the whole sequence's second traceback. Not below 16: the
/// start of a gapped alignment is found by reading up to 12 bases on each side of the seed
/// (`blast_get_start_for_gapped_alignment_nucl`), outside the kernels' fence tests.
fn x_tb_range_pad() -> usize {
    static PAD: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *PAD.get_or_init(|| {
        std::env::var("LOSAT_X_BNTBRANGEPAD")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|pad| (16..=1 << 20).contains(pad))
            .unwrap_or(1024)
    })
}

fn x_partial_ranges_padded(
    spans: impl IntoIterator<Item = (usize, usize)>,
    length: usize,
    pad: usize,
) -> Vec<(usize, usize)> {
    #[allow(non_snake_case)]
    let (OVERHANG, MINGAP) = (pad, pad);
    let mut ranges: Vec<(usize, usize)> = spans
        .into_iter()
        .map(|(from, to)| {
            (
                from.saturating_sub(OVERHANG).min(length),
                to.saturating_add(OVERHANG).min(length),
            )
        })
        .collect();
    ranges.sort_unstable_by_key(|range| range.0);
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
    for (from, to) in ranges {
        match merged.last_mut() {
            Some(last) if from <= last.1.saturating_add(MINGAP) => last.1 = last.1.max(to),
            _ => merged.push((from, to)),
        }
    }
    merged
}

/// EXPERIMENT (LOSAT_X_BNTBRANGE): sequences decoded in ranges, those traced again whole, their
/// bases and the bases of their ranges (`x_tb_range_print_summary`).
static X_TB_RANGE_COUNTS: [std::sync::atomic::AtomicU64; 4] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 4];

fn x_tb_range_count(length: usize, ranges: &[(usize, usize)], fenced: bool) {
    use std::sync::atomic::Ordering::Relaxed;
    let in_ranges: usize = ranges.iter().map(|range| range.1 - range.0).sum();
    X_TB_RANGE_COUNTS[0].fetch_add(1, Relaxed);
    X_TB_RANGE_COUNTS[1].fetch_add(u64::from(fenced), Relaxed);
    X_TB_RANGE_COUNTS[2].fetch_add(length as u64, Relaxed);
    X_TB_RANGE_COUNTS[3].fetch_add(in_ranges as u64, Relaxed);
}

/// EXPERIMENT (LOSAT_X_BNTBRANGELOG=1): one stderr line with the counts of LOSAT_X_BNTBRANGE.
pub(super) fn x_tb_range_print_summary() {
    use std::sync::atomic::Ordering::Relaxed;
    if std::env::var_os("LOSAT_X_BNTBRANGELOG").is_none_or(|v| v != "1") {
        return;
    }
    let [ranged, fenced, bases, in_ranges] = X_TB_RANGE_COUNTS.each_ref().map(|c| c.load(Relaxed));
    eprintln!(
        "[X_BNTBRANGE] sequences in ranges {ranged}, traced again whole {fenced}, bases {bases}, in ranges {in_ranges}"
    );
}

/// The headers of the database sequences that the final hit lists of the query batches searched
/// so far name, by ordinal id (`read_titles`): the stored title, or the definition lines of a
/// sequence with parsed Seq-ids.
pub(super) type DatabaseTitles = BTreeMap<u32, ReportedHeader>;

/// NCBI's fetches with `check_oid_exclusion` of a query batch's search of `database`
/// (`shared::FetchExcluded`; `None` where they exclude nothing): the traceback's of the
/// subjects of the preliminary hit lists (`SubjectPasses::traceback_pass`), or in an
/// `-ungapped` search the re-evaluation's of each subject with HSPs, before the preliminary
/// hit lists (`SubjectPasses::prelim_pass`).
pub(super) fn fetch_excluded(database: &BlastnDatabase) -> Option<shared::FetchExcluded<'_>> {
    shared::FetchExcluded::of(&database.db)
}

/// LOSAT's refusal of a search whose fetch with `check_oid_exclusion` cannot read the header
/// of a subject (`shared::FetchExcluded::refusal`), where NCBI's search of the batch fails: after the
/// reports of the batches before, which NCBI's batch loop has printed (`super::after_batch_reports`;
/// for the first batch, the prologs), as TBLASTN's and TBLASTX's (line RG item X-short-bat-2).
pub(super) fn fetch_refusal(failed: (u32, crate::seqdb::SeqDbError)) -> anyhow::Error {
    super::after_batch_reports(shared::FetchExcluded::refusal(failed, PROGRAM))
}

impl BlastnDatabase {
    /// The titles of the sequences that the final hit lists `hit_lists` of a query batch name,
    /// read where NCBI makes the batch's Seq-aligns ([`read_titles`], into the titles of `reads`),
    /// after those of no Seq-align leave the lists ([`drop_unaligned`]) and after the batch's
    /// fetches (`fetched`: the first that read a UTF8String header). NCBI's warning of the batch's
    /// first read of such a header goes to the posts of `reads` (the batch's reader messages and
    /// request posts, before its prefetch), also where the read fails after it
    /// (`shared::Utf8Read::of_read`). The error is `read_titles`'.
    pub(super) fn batch_titles(
        &self,
        (fetched, one_by_one): (Option<crate::seqdb::SeqLocation>, bool),
        hit_lists: &mut [Option<BlastnHitList>],
        (titles, posts): (&mut DatabaseTitles, &mut crate::report::diag_log::Posts),
    ) -> Result<()> {
        let aligned = drop_unaligned(self, hit_lists);
        let earlier = shared::EarlierUtf8::of_batch(0, fetched, aligned);
        let read = read_titles(self, hit_lists, (earlier, one_by_one), titles);
        let first = shared::Utf8Read::of_read(&read, Option::as_ref);
        if let Some(warning) = shared::Utf8Read::post_of(first, &self.db, "blastn") {
            posts.append(&warning);
        }
        read.map(|_| ())
    }
}

/// The final hit lists `hit_lists` of a query batch without the sequences that NCBI makes no
/// Seq-align of (`shared::SeqAligned`): after their cut to the hit list size, before
/// `read_titles`. The first of its reads that took a UTF8String for a VisibleString is returned
/// (`shared::SeqAligned::first_utf8`).
pub(super) fn drop_unaligned(
    database: &BlastnDatabase,
    hit_lists: &mut [Option<BlastnHitList>],
) -> Option<crate::seqdb::SeqLocation> {
    // Every blastn search fetches its subjects with `check_oid_exclusion`: the gapped traceback,
    // or the ungapped re-evaluation of a nucleotide subject (`fetch_excluded`).
    let mut aligned = shared::SeqAligned::of(&database.db, true)?;
    for hit_list in hit_lists.iter_mut().flatten() {
        hit_list.hsplist_array.truncate(hit_list.hsplist_count);
        aligned.retain(&mut hit_list.hsplist_array, |list| Some(list.oid));
        hit_list.hsplist_count = hit_list.hsplist_array.len();
    }
    aligned.first_utf8()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1582
// NCBI first reads a header file (and opens it, c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:127-133) when the traceback's
// results become Seq-aligns: after the batch's traceback, for the final hit lists only
// (after their sort by e-value and their cut to the hit list size, before the HSPs of
// score 0 are dropped), query by query and in the order of each hit list. A header file
// that is missing or cannot be read fails the batch there (`BLAST Database error`, exit
// status 2, blast_app_util.hpp:185-203, outside the search stage's handler of
// local_blast.cpp:240-251), after the outfmt 0 prolog and the reports of the batches
// before; a batch without hits reads no header. A title read in an earlier batch is not
// read again: the header file stays open, and the title is the same.
/// Reads the definition lines of the sequences that the final hit lists `hit_lists` of a
/// query batch name, in NCBI's order, into `titles` (`common::database::read_header`). A reader
/// error of NCBI's own is an `ErrorAfterBatchReports` (it fails the batch); LOSAT's refusals are errors
/// of their own. The batch's first read of a record with a UTF8String where a VisibleString is
/// expected, its fetches' and Seq-aligns' (`earlier`) or these, is returned
/// (`shared::HeaderReader::first_utf8`): NCBI posts its warning there.
pub(super) fn read_titles(
    database: &BlastnDatabase,
    hit_lists: &[Option<BlastnHitList>],
    (earlier, one_by_one): (shared::EarlierUtf8, bool),
    titles: &mut DatabaseTitles,
) -> Result<Option<shared::Utf8Read>> {
    let headers = shared::HeaderReader::new(&database.db, one_by_one);
    headers.note_earlier(earlier);
    for (query, hit_list) in hit_lists.iter().enumerate() {
        let Some(hit_list) = hit_list else {
            continue;
        };
        let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
        for hsp_list in &hit_list.hsplist_array[..count] {
            if titles.contains_key(&hsp_list.oid) {
                continue;
            }
            let title = shared::read_header(&headers, hsp_list.oid, query, PROGRAM)?;
            titles.insert(hsp_list.oid, title);
        }
    }
    Ok(headers.first_utf8())
}

/// What the reports `output_formats` show of a reported sequence with parsed Seq-ids
/// (`reported_subjects`); outfmt 13 to 16 and the tabular reports of taxonomy fields write its
/// taxids; the first format's batches prefetch the subjects.
pub(super) fn parsed_shown(output_formats: &[BlastnOutputFormat]) -> ParsedShown {
    let pairwise = output_formats.iter().any(|format| format.is_pairwise());
    ParsedShown {
        tabular: output_formats
            .iter()
            .any(|format| format.tabular_reads_subject()),
        stops_at_failed_lookup: output_formats.iter().any(|format| {
            format.is_tabular()
                && crate::report::tabular::may_stop_at_failed_lookup(format.tabular_spec())
        }),
        pairwise,
        structured: output_formats
            .iter()
            .any(|format| format.is_structured_xml()),
        joined: pairwise
            || output_formats.iter().any(|format| {
                format.is_structured_xml()
                    || format
                        .tabular_spec()
                        .is_some_and(|spec| spec.needs().joins_subjects())
            }),
        sam: output_formats
            .iter()
            .any(|format| matches!(format, BlastnOutputFormat::Sam(_))),
        taxonomy: output_formats.iter().any(|format| {
            format.xml2_kind().is_some()
                || format
                    .tabular_spec()
                    .is_some_and(|spec| spec.needs().taxonomy)
        }),
        words: output_formats
            .iter()
            .filter(|format| format.is_tabular())
            .fold(PrintedIds::default(), |words, format| {
                words.or(PrintedIds::of(format.tabular_spec()))
            }),
        prefetched: output_formats
            .first()
            .is_some_and(|format| format.prefetches_subjects())
            && !shared::prefetch_always_off(prefetch_format(output_formats)),
        one_by_one: shared::prefetch_may_be_off(prefetch_format(output_formats)),
    }
}

/// The first of the formats `output_formats` of a run as NCBI's prefetch of each query batch asks
/// it (`shared::PrefetchFormat`, the one owner of whether a batch prefetches); a run writes at
/// least one format, outfmt 0 stands for none.
pub(super) fn prefetch_format(output_formats: &[BlastnOutputFormat]) -> shared::PrefetchFormat {
    output_formats.first().map_or(
        shared::PrefetchFormat::of(
            crate::blastinput::app::output_format::PAIRWISE,
            crate::report::archive::options::RequestKind::Blastn,
        ),
        |format| format.prefetch_format(),
    )
}

/// The Seq-aligns of the results `hit_lists` that NCBI's prefetch counts against
/// `PRE_FETCH_SEQS_LIMIT` (`shared::prefetch_off`): one for each HSP of a `gapped` search, one
/// for each subject of an ungapped one (`shared::subject_seq_aligns`).
/// The query batches of a call before a batch's prefetch that the loop holds
/// (`shared::LoadedBefore`, read by `BlastnDatabase::loaded_subjects`): their final hit lists by
/// query (`hits`), the queries of each batch (`spans`; those past `hits` are not held), and which
/// rows of the first format read a subject where a batch's prefetch is off (`rows`).
pub(super) type LoadedBefore<'a> = shared::LoadedBefore<'a, Option<BlastnHitList>>;

pub(super) fn seq_aligns_of(hit_lists: &[Option<BlastnHitList>], gapped: bool) -> usize {
    hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| {
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            &hit_list.hsplist_array[..count]
        })
        .map(|hsp_list| shared::subject_seq_aligns(hsp_list.hsps.len(), gapped))
        .sum()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:728-753 (cited at
// `common::database::prefetch_may_be_off`: `s_PreFetchSeqs` counts the Seq-aligns of one batch's
// results against `PRE_FETCH_SEQS_LIMIT`)
/// Whether the batch of each query of `hit_lists` (final hit lists by query) prefetches
/// (`shared::BatchPrefetch::of`): `batches` are the queries of each batch, numbered from `first`
/// (the query of `hit_lists[0]`), and `seq_aligns_of` counts a batch's Seq-aligns (`gapped`). The
/// one producer of the prefetch that the reports' loads read (`hit_loads`), for both of its
/// readers: the stop of a later batch, of the batches before it
/// (`BlastnDatabase::batch_prefetch_stop`), and the reports (`reported_subjects`), so the two
/// cannot disagree on what the batches before loaded (review 2 of line DBID pfoff-1, N-M2).
fn batches_prefetch(
    format: shared::PrefetchFormat,
    hit_lists: &[Option<BlastnHitList>],
    first: usize,
    batches: &[std::ops::Range<usize>],
    gapped: bool,
) -> shared::BatchPrefetch {
    shared::BatchPrefetch::of(format, hit_lists.len(), first, batches, |span| {
        seq_aligns_of(&hit_lists[span], gapped)
    })
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(super) struct ReportedSubjects {
    /// The ordinal id of each rank, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The record of each rank: the id `gnl|BL_ORD_ID|n` and the stored title (a sequence with
    /// parsed Seq-ids: the title of its first definition line), without letters.
    pub(super) records: Vec<FastaRecord>,
    /// When the reports show letters, the stored data of each rank (`reported_store`): its
    /// 2-bit data and ambiguity runs, from which a report makes the letters of the ranges it
    /// shows (`SubjectLetters`), not the letters of the whole sequence.
    pub(super) stores: Option<Vec<PackedNucleotide>>,
    /// The sequences with parsed Seq-ids by rank (empty when none is reported): their ids and
    /// titles come from their definition lines (`report::db_subject`).
    pub(super) parsed: Vec<Option<Arc<ParsedSubject>>>,
    /// The names of the sources from which NCBI makes the empty stored titles of the others, by
    /// rank (`shared::StoredTitle::source`; empty when none has one).
    pub(super) sources: Vec<crate::report::defline::SourceTaxname>,
    /// What the XML2 and JSON reports do at each subject's `Hit`, by rank (`hit_loads`, the one
    /// answer: a `Hit` whose lookup throws, a subject that no report loads), which the refusals,
    /// the names and the taxonomy reads take.
    pub(super) loads: HitLoads,
}

/// The hit lists of `hit_list` that the reports show (its first `hsplist_count`).
fn shown_lists(hit_list: &Option<BlastnHitList>) -> &[BlastnHspList] {
    hit_list.as_ref().map_or(&[][..], |hit_list| {
        &hit_list.hsplist_array[..hit_list.hsplist_count.min(hit_list.hsplist_array.len())]
    })
}

/// The most identities of an HSP of `hsp_list` (the Seq-align's `num_ident`), which tells whether
/// a row of `nident` alone reads its subject (`SubjectFetch::WithIdentities`).
fn most_identities(hsp_list: &BlastnHspList) -> usize {
    hsp_list
        .hsps
        .iter()
        .map(|hsp| hsp.num_ident)
        .max()
        .unwrap_or(0)
}

/// What a pairwise report (outfmt 0 to 4) of the hit lists `hit_lists` reads, by the key that
/// `hit_lists` names the subjects by (`QueryGroup::read_order` of the HSPs as the writer reads
/// them, `reading`; `parsed` gives the subject with parsed Seq-ids of a key, whose copies join a
/// hit). The rows that the writer reads are these HSPs (`pairwise::Rows::list`: shifted to the
/// records of a ranged query by `query_offset`, with the run's `ungapped` flag deciding whether a
/// negative Seq-align start is lifted by `-query_loc` or stays class N with an unsigned `TSeqPos`
/// range; the rows name the subjects by rank, in the order of the ordinal ids).
fn pairwise_read_order<'p>(
    hit_lists: &[Option<BlastnHitList>],
    query_offset: impl Fn(usize) -> usize,
    ungapped: bool,
    reading: crate::report::pairwise::PairwiseReading,
    parsed: impl Fn(u32) -> Option<&'p ParsedSubject>,
) -> crate::report::db_subject::PairwiseReads {
    crate::report::x_pairwise_par::QueryGroup::read_order(
        hit_lists.iter().enumerate().flat_map(|(query, hit_list)| {
            let offset = query_offset(query);
            shown_lists(hit_list)
                .iter()
                .filter(|hsp_list| !hsp_list.hsps.is_empty())
                .flat_map(move |hsp_list| {
                    hsp_list
                        .hsps
                        .iter()
                        .map(move |hsp| (hsp_list.oid, (hsp, offset)))
                })
        }),
        |&(hsp, offset)| {
            crate::algorithm::blastn::left_of_query::blastn_hsp_scores(hsp, offset, ungapped)
        },
        |previous, oid| match (parsed(previous), parsed(oid)) {
            (Some(earlier), Some(later)) => {
                crate::report::db_subject::joins_as_copy(earlier, later)
            }
            _ => false,
        },
        reading,
        false,
        None,
    )
}

/// The subjects of the query report of `hit_list` for `HitLoads::of` where the run's formats
/// are pairwise (`pairwise`: how its report reads them, with `after` a pairwise one): the
/// report's reads (`pairwise_read_order`, with `query_offset` and `ungapped` matching the
/// writer). `None` otherwise (the HSP lists' order).
fn pairwise_keys<'p>(
    hit_list: &Option<BlastnHitList>,
    query_offset: usize,
    ungapped: bool,
    after: crate::report::db_subject::FormatThrow,
    pairwise: Option<crate::report::pairwise::PairwiseReading>,
    parsed: impl Fn(u32) -> Option<&'p ParsedSubject>,
) -> Option<crate::report::db_subject::ReportKeys> {
    let reading = pairwise.filter(|_| after.run().is_some_and(|run| run.after.is_pairwise()))?;
    Some(crate::report::db_subject::ReportKeys::Pairwise(
        pairwise_read_order(
            std::slice::from_ref(hit_list),
            |_| query_offset,
            ungapped,
            reading,
            parsed,
        ),
    ))
}

/// What the reports of `hit_lists` load where their batch's prefetch is off (`RunLoads::loads` of
/// those reports alone, by ordinal id, after the reports before them in the run, `run`): the
/// subjects of its reads one by one (`DatabaseTitles`' texts of a batch whose prefetch is off),
/// not one whose load throws, nor one that no report reaches after it. `pairwise`: how a
/// pairwise report reads them (`pairwise_keys`, with `placements` from `first_query` and
/// `ungapped`); `before`: the loads of the run's record and of the calls before
/// (`RunLoads::loads`), whose last report is of the batch that `before` names
/// (`TaxdbRun::loads_before`), in whose scope these reports read.
fn report_loads(
    run: &crate::report::db_subject::RunLoads,
    titles: &DatabaseTitles,
    hit_lists: &[Option<BlastnHitList>],
    (pairwise, before): (
        Option<crate::report::pairwise::PairwiseReading>,
        crate::report::db_subject::LoadsBefore,
    ),
    placements: &crate::blastinput::seq_range::Placements,
    first_query: usize,
    ungapped: bool,
) -> HitLoads {
    let read = crate::report::db_subject::ReportRead {
        one_by_one: true,
        batch: before.batch(),
    };
    run.loads(
        before,
        hit_lists.iter().enumerate().map(|(query, hit_list)| {
            let parsed = |oid: u32| titled_parsed(titles, oid as usize);
            let keys = pairwise_keys(
                hit_list,
                placements.offset(first_query + query),
                ungapped,
                run.answer(),
                pairwise,
                parsed,
            )
            .unwrap_or_else(|| {
                crate::report::db_subject::ReportKeys::Entries(
                    shown_lists(hit_list)
                        .iter()
                        .map(|hsp_list| {
                            (!hsp_list.hsps.is_empty()).then_some(hsp_list.oid as usize)
                        })
                        .collect(),
                )
            });
            (read, keys)
        }),
        |ordinal| titled_parsed(titles, ordinal),
        |report, list| most_identities(&shown_lists(&hit_lists[report])[list]),
    )
}

/// The fetches of the query batches of a run whose first format reads the subjects one by one
/// where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off (`shared::prefetch_off`), and the
/// open of the taxonomy database by the rows of a batch (`rows_open`).
pub(super) struct QueryFetches<'a> {
    pub(super) database: &'a BlastnDatabase,
    pub(super) titles: &'a DatabaseTitles,
    /// The run's record of the loads after a throw (`report_loads`).
    pub(super) run: &'a crate::report::db_subject::RunLoads,
    /// Whether the first format writes the texts where it reads each subject
    /// (`ReportFormat::places_subject_reads`), else with the query's warnings.
    pub(super) places: bool,
    /// Whether the first format is a table with a plural name, which looks every taxid of its
    /// rows up, taxid 0 too (`TabularSpec::looks_up_every_taxid`), with `taxdb.bti` found.
    pub(super) plural_names: bool,
    /// How the run's pairwise reports read the subjects (`pairwise_reading`), for the loads.
    pub(super) pairwise: Option<crate::report::pairwise::PairwiseReading>,
    /// Whether the first format is a pairwise one, which places the texts in the order of its
    /// reads (`places`).
    pub(super) pairwise_first: bool,
    /// Whether each HSP is a Seq-align (`seq_aligns_of`) and the search is gapped (`!ungapped`).
    pub(super) gapped: bool,
    /// The `-query_loc` offsets of the run's queries (for `pairwise_read_order`).
    pub(super) placements: &'a crate::blastinput::seq_range::Placements,
}

impl QueryFetches<'_> {
    /// The open of the taxonomy database by the rows of `hit_lists` that look every taxid up
    /// (`plural_names`, `BlastnDatabase::open_for_tabular_rows`): rows of a subject that the
    /// reports load (`loads`; all of them where the batch's prefetch loads them, `None`). A row
    /// whose load throws keeps the names of the rows before and looks nothing up (review 1 of
    /// line DBID pfoff-2, H1).
    pub(super) fn rows_open(
        &self,
        hit_lists: &[Option<BlastnHitList>],
        loads: Option<&HitLoads>,
    ) -> Result<()> {
        let looks_up = |hsp_list: &BlastnHspList| {
            loads.is_none_or(|loads| {
                !hsp_list.hsps.is_empty() && loads.loaded(hsp_list.oid as usize)
            })
        };
        if self.plural_names
            && hit_lists
                .iter()
                .flatten()
                .any(|hit_list| hit_list.hsplist_array.iter().any(looks_up))
        {
            self.database.open_for_tabular_rows()?;
        }
        Ok(())
    }

    /// The fetches of each query of the batch whose final hit lists are `batch` (the query of
    /// `batch[0]` is the run's `first`) where its prefetch is off: each query's report loads the
    /// subjects that it shows (`report_loads`, after the reports before it in the run), where
    /// the first format reads them (`fetches`), with the taxonomy texts of each first read
    /// placed by subject in `query_reads` (`places`) or returned for the query's warnings; then
    /// its rows may open the taxonomy database (`rows_open`).
    pub(super) fn batch(
        &self,
        batch: &[Option<BlastnHitList>],
        first: usize,
        fetches: bool,
        query_reads: &mut [crate::report::query_warnings::QueryReads],
    ) -> Result<Vec<Vec<u8>>> {
        batch
            .iter()
            .enumerate()
            .map(|(index, hit_list)| {
                let shown = std::slice::from_ref(hit_list);
                // The loop records a batch after its texts here (`note_loaded`).
                let before = self.database.taxdb.loads_before(Some(first));
                let loads = report_loads(
                    self.run,
                    self.titles,
                    shown,
                    (self.pairwise, before),
                    self.placements,
                    first + index,
                    !self.gapped,
                );
                let mut messages = Vec::new();
                if fetches && seq_aligns_of(shown, self.gapped) > 0 {
                    if self.places {
                        query_reads[first + index] = self.database.subject_read_texts(
                            self.titles,
                            shown,
                            self.plural_names,
                            self.pairwise.filter(|_| self.pairwise_first),
                            self.placements,
                            first + index,
                            !self.gapped,
                            &loads,
                        )?;
                    } else {
                        messages = self
                            .database
                            .prefetch_messages(self.titles, shown, &loads)?;
                    }
                }
                self.rows_open(shown, Some(&loads))?;
                Ok(messages)
            })
            .collect()
    }
}

/// What the formats `output_formats` do after a subject whose lookup throws where its batch's
/// prefetch is off (`FormatThrow::of_formats`, for `hit_loads` and the run's `RunLoads`), with the
/// call's refusal over the `-outfmt` values of its other searches, `others`
/// (`FormatThrow::of_search`, `ReportOutputs::other_outfmts`); `html`: the run has `-html` links
/// (`BlastnArgs::html_links`, which the writers' `PairwiseConfig::html` gets too).
pub(super) fn after_throw(
    output_formats: &[BlastnOutputFormat],
    others: &[&str],
    html: bool,
) -> crate::report::db_subject::FormatThrow {
    crate::report::db_subject::FormatThrow::of_formats(
        output_formats
            .iter()
            .map(|format| format.format_throw(html)),
    )
    .of_search(
        others,
        crate::blastinput::app::Program::Blastn,
        crate::report::archive::options::RequestKind::Blastn,
        html,
    )
}

/// How the pairwise reports of the run's formats `output_formats` read each query's subjects
/// (`QueryGroup::read_order`): as the first pairwise format's view reads them
/// (`PairwiseView::reading` of `-sorthits` and `-sorthsps` `sort` and the run's numbers
/// `counts`; the formats of a run whose loads they make share one view's order,
/// `FormatThrow::and`); `None` without one.
pub(super) fn pairwise_reading(
    output_formats: &[BlastnOutputFormat],
    sort: crate::report::hit_sort::ReportSort,
    counts: crate::blastinput::report_counts::ReportCounts,
) -> Option<crate::report::pairwise::PairwiseReading> {
    output_formats.iter().find_map(|format| match format {
        BlastnOutputFormat::Pairwise(view) => Some(view.reading(sort, counts.shown())),
        _ => None,
    })
}

/// Whether every format of `output_formats` writes the same where every load of a subject with a
/// positive taxid throws (`BlastnOutputFormat::unchanged_by_load_throws`, for the batch's
/// prefetch, `shared::PrefetchQueries::where_reports_unchanged_by_load_throws`).
pub(super) fn reports_unchanged_by_load_throws(output_formats: &[BlastnOutputFormat]) -> bool {
    output_formats
        .iter()
        .all(BlastnOutputFormat::unchanged_by_load_throws)
}

/// What the reports of `hit_lists` (final hit lists by query, in the order of the reports; the
/// subjects named by `key`, an ordinal id or a rank) do at each subject's `Hit` (`HitLoads::of`,
/// the one answer): `prefetch` tells by query whether its batch reads one by one
/// (`BatchPrefetch::reads_one_by_one`), `parsed` gives the subject with parsed Seq-ids of a
/// key; the `Hit`s are those that `blast_xml::write_blastn_xml` makes, or the pairwise reports'
/// reads (`pairwise`: how they read them, `pairwise_keys`), after the loads `before`
/// (`HitLoads::of`).
pub(super) fn hit_loads<'p>(
    hit_lists: &[Option<BlastnHitList>],
    prefetch: &shared::BatchPrefetch,
    (after, pairwise): (
        crate::report::db_subject::FormatThrow,
        Option<crate::report::pairwise::PairwiseReading>,
    ),
    placements: &crate::blastinput::seq_range::Placements,
    ungapped: bool,
    before: crate::report::db_subject::LoadsBefore,
    parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
) -> HitLoads {
    let lists = |query: usize| {
        hit_lists[query]
            .as_ref()
            .map_or(&[][..], |hit_list| hit_list.hsplist_array.as_slice())
    };
    HitLoads::of(
        before,
        after,
        (0..hit_lists.len()).map(|query| {
            // A batch that prefetches loads every subject, in any order.
            let read = prefetch.report_read(query, 0);
            let one_by_one = read.one_by_one;
            let keys = pairwise_keys(
                &hit_lists[query],
                placements.offset(query),
                ungapped,
                after,
                pairwise.filter(|_| one_by_one),
                |oid| parsed(oid as usize),
            )
            .unwrap_or_else(|| {
                crate::report::db_subject::ReportKeys::Entries(
                    lists(query)
                        .iter()
                        .map(|hsp_list| hsp_list.hsps.first().map(|hsp| hsp.s_idx as usize))
                        .collect(),
                )
            });
            (read, keys)
        }),
        &parsed,
        |query, list| most_identities(&lists(query)[list]),
    )
}

/// The subject with parsed Seq-ids of the ordinal id `ordinal` in `titles`.
fn titled_parsed(titles: &DatabaseTitles, ordinal: usize) -> Option<&ParsedSubject> {
    match titles.get(&u32::try_from(ordinal).ok()?) {
        Some(ReportedHeader::Parsed(subject)) => Some(subject),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1406-1414
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:80-83
// The formatter reads a database sequence from its Bioseq: the stored 2-bit data, or with
// ambiguity codes their NCBI4na values (`x_GetAmbigSeq`'s letters before their mapping to
// BLASTNA). LOSAT keeps the stored 2-bit data and ambiguity runs (`reported_store`), and a
// report writes `NCBI4NA_TO_IUPACNA` of each letter of the ranges it shows (the letters of
// `SequenceBlock::nucleotide_iupac` there). Upper case: a database has no lowercase letters.
/// The reported database sequences of `hit_lists` (the final hit lists of all query
/// batches, after their HSPs of score 0 were dropped): the hit lists name them by rank
/// from here on (the rank grows with the ordinal id, so every comparison that breaks a tie
/// by the ordinal id gives the same order), and their records are made for the reports
/// from `titles`, with their stored data when `letters` (outfmt 0, the hit records and the
/// formats with letters read them; `reported_store`). A reported sequence whose stored title LOSAT
/// does not report as NCBI does is LOSAT's refusal, and so is one with parsed Seq-ids that the
/// reports `shown` do not write as NCBI (`common::database::check_reported_parsed`, and with
/// `shown.joined` `refuse_joined_hits`; `first_query`: the index of the query of `hit_lists[0]`;
/// `batches`: the queries of each batch, numbered from `first_query`, whose prefetch, where only
/// the prefetch reads a subject, is `batches_prefetch`'s, with the Seq-aligns of a `gapped`
/// search or not; `after`: what the formats do after a subject whose lookup throws
/// (`after_throw`), for the one answer of the reports' loads, `hit_loads`, which the refusals
/// take, with how the pairwise reports read them, `pairwise`, `pairwise_reading`, after the loads
/// of the batches before that the run recorded, `BlastnDatabase::loaded_before_batch`, and of the
/// calls before, `shared::TaxdbRun::unrecorded_calls`). The records are made on `pool`
/// (`x_reported_records`).
#[allow(clippy::too_many_arguments)]
pub(super) fn reported_subjects(
    database: &BlastnDatabase,
    titles: &DatabaseTitles,
    hit_lists: &mut [Option<BlastnHitList>],
    letters: bool,
    (shown, format): (ParsedShown, shared::PrefetchFormat),
    (after, pairwise): (
        crate::report::db_subject::FormatThrow,
        Option<crate::report::pairwise::PairwiseReading>,
    ),
    placements: &crate::blastinput::seq_range::Placements,
    batches: &[std::ops::Range<usize>],
    first_query: usize,
    gapped: bool,
    pool: Option<&SearchPool<'_>>,
) -> Result<ReportedSubjects> {
    let mut ordinals: Vec<u32> = hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| hit_list.hsplist_array.iter())
        .filter(|hsp_list| !hsp_list.hsps.is_empty())
        .map(|hsp_list| hsp_list.oid)
        .collect();
    ordinals.sort_unstable();
    ordinals.dedup();
    // What the batches read of each subject (`prefetch`: a batch whose Seq-aligns outnumber a
    // positive `PRE_FETCH_SEQS_LIMIT` reads nothing of them before its reports, which read
    // them one by one; `shared::BatchReads::shown`); where the prefetch fails, the batch writes
    // NCBI's warning (`BlastnDatabase::batch_prefetch`).
    let prefetch = batches_prefetch(format, hit_lists, first_query, batches, gapped);
    let reads = prefetch.reads(
        hit_lists
            .iter()
            .enumerate()
            .filter_map(|(query, hit_list)| Some((query, hit_list.as_ref()?)))
            .flat_map(|(query, hit_list)| {
                hit_list
                    .hsplist_array
                    .iter()
                    .filter(|hsp_list| !hsp_list.hsps.is_empty())
                    .map(move |hsp_list| (query, hsp_list.oid))
            }),
    );
    // What the reports do at each subject's `Hit`, by ordinal id (renamed to ranks below).
    let recorded = database
        .loaded_before_batch
        .lock()
        .expect("the loads before the batch")
        .clone();
    let loads = hit_loads(
        hit_lists,
        &prefetch,
        (after, pairwise),
        placements,
        !gapped,
        recorded.with_calls_before(database.taxdb.unrecorded_calls()),
        |ordinal| titled_parsed(titles, ordinal),
    );
    let rank = |ordinal: u32| -> u32 {
        ordinals
            .binary_search(&ordinal)
            .expect("every reported subject was listed") as u32
    };
    for hit_list in hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            if hsp_list.hsps.is_empty() {
                continue;
            }
            hsp_list.oid = rank(hsp_list.oid);
            for hsp in &mut hsp_list.hsps {
                hsp.s_idx = rank(hsp.s_idx);
            }
        }
    }
    // Whether `taxdb.bti` is on the search path, looked up at the first sequence with parsed
    // Seq-ids that outfmt 0 shows.
    let taxdb_found = std::sync::OnceLock::new();
    let db = &database.db;
    let record = |&ordinal: &u32| -> Result<Reported> {
        let header = titles
            .get(&ordinal)
            .with_context(|| format!("database sequence {} has no title read", db.oid(ordinal)))?;
        let (title, subject, source) = match header {
            ReportedHeader::Stored(stored) => {
                shared::check_reported_title(stored.replaced_bytes, || db.oid(ordinal), PROGRAM)?;
                if shown.pairwise {
                    stored.check_shown_source(|| db.oid(ordinal), PROGRAM)?;
                }
                (stored.title.clone(), None, stored.source.clone())
            }
            ReportedHeader::Parsed(subject) => {
                let taxdb = shown.pairwise && *taxdb_found.get_or_init(|| db.taxdb_found());
                let shown = reads.shown(shown, ordinal);
                shared::check_reported_parsed(
                    db,
                    subject,
                    shown,
                    loads.load(ordinal as usize),
                    false,
                    taxdb,
                    ordinal,
                    PROGRAM,
                )?;
                (
                    subject.first_title().to_vec(),
                    Some(Arc::clone(subject)),
                    Default::default(),
                )
            }
        };
        let store = if letters {
            Some(reported_store(db, ordinal)?)
        } else {
            None
        };
        let record = FastaRecord {
            // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the OID
            // (line DI-8).
            local_id: db.oid(ordinal).fasta_id(),
            title,
            sequence: Vec::new(),
            warnings: Vec::new(),
            source: Default::default(),
        };
        Ok((record, store, subject, source))
    };
    let made = x_reported_records(&ordinals, pool, &record)?;
    let mut records = Vec::with_capacity(ordinals.len());
    let mut stores = letters.then(|| Vec::with_capacity(ordinals.len()));
    let mut parsed: Vec<Option<Arc<ParsedSubject>>> = Vec::new();
    let mut sources: Vec<crate::report::defline::SourceTaxname> = Vec::new();
    for (rank, (record, store, subject, source)) in made.into_iter().enumerate() {
        records.push(record);
        if let (Some(stores), Some(store)) = (stores.as_mut(), store) {
            stores.push(store);
        }
        if let Some(subject) = subject {
            if parsed.is_empty() {
                parsed = vec![None; ordinals.len()];
            }
            parsed[rank] = Some(subject);
        }
        shared::keep_source(&mut sources, rank, ordinals.len(), source);
    }
    if shown.joined && !parsed.is_empty() {
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            shared::refuse_joined_hits(
                first_query + query,
                hit_list
                    .hsplist_array
                    .iter()
                    .filter(|hsp_list| !hsp_list.hsps.is_empty())
                    .map(|hsp_list| {
                        let rank = hsp_list.oid as usize;
                        Some((ordinals[rank], parsed[rank].as_deref()?))
                    }),
                db.oid_map(),
                PROGRAM,
            )?;
        }
    }
    let loads = loads.renamed(|ordinal| {
        let ordinal = u32::try_from(ordinal).ok()?;
        ordinals.binary_search(&ordinal).ok()
    });
    Ok(ReportedSubjects {
        ordinals,
        records,
        stores,
        parsed,
        sources,
        loads,
    })
}

// No NCBI counterpart: NCBI's formatter reads the Bioseq of its scope (the references above
// `reported_subjects`); this is the form in which LOSAT holds a reported sequence's data.
/// The stored data of database sequence `ordinal` that a report reads its letters from: its
/// 2-bit data and ambiguity runs (`PackedNucleotide::from_stored`, whose letters are those of
/// `SequenceBlock::nucleotide_iupac`), with the errors of `nucleotide_iupac` in the same order
/// (`SequenceBlock::nucleotide_stored`). The store is the read's buffer cut to the 2-bit data, a
/// quarter of a byte per letter, as `common::database::x_read_nucleotide_store` keeps it: with
/// LOSAT_X_BNMAPPED a buffer of 1 MiB or more is mapped and leaves at the store's drop (it is not
/// given back to `x_mapped::X_BYTES`), and no copy is made in the heap.
fn reported_store(db: &SeqDb, ordinal: u32) -> Result<PackedNucleotide> {
    let block = x_read_sequence(db, ordinal).map_err(database_error)?;
    let (length, runs) = block.nucleotide_stored(0).map_err(database_error)?;
    let packed = block.x_into_nucleotide_packed();
    if crate::utils::x_mapped::x_is_mapped(packed.capacity()) {
        crate::utils::x_mapped::X_BYTES.x_disown(&packed);
    }
    Ok(PackedNucleotide::x_from_stored_buffer(
        packed, length, &runs,
    ))
}

/// A reported sequence's record, its stored data (`reported_store`, when the reports show
/// letters) and its parsed Seq-ids (`reported_subjects`).
type Reported = (
    FastaRecord,
    Option<PackedNucleotide>,
    Option<Arc<ParsedSubject>>,
    crate::report::defline::SourceTaxname,
);

// No NCBI counterpart: NCBI's formatters read the sequence of each alignment when they write it.
/// The sequences of a report that one read at a time cost more than a task of the pool.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const X_REPORTED_RECORDS_MIN: usize = 64;

/// `record` of each of `ordinals`, in their order, on `pool` while the search's threads wait
/// between the query batches (each batch's report reads its own reported sequences with
/// LOSAT_X_RPTSTREAM); the error of the first one in that order that has one.
fn x_reported_records(
    ordinals: &[u32],
    pool: Option<&SearchPool<'_>>,
    record: &(dyn Fn(&u32) -> Result<Reported> + Sync),
) -> Result<Vec<Reported>> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if let Some(pool) = pool.filter(|pool| pool.enabled()) {
        use rayon::prelude::*;
        let made: Vec<_> = pool.install(|| {
            ordinals
                .par_iter()
                .with_min_len(X_REPORTED_RECORDS_MIN)
                .map(record)
                .collect()
        });
        return made.into_iter().collect();
    }
    let _ = pool;
    ordinals.iter().map(record).collect()
}

impl ReportedSubjects {
    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every
    /// reported sequence, by rank (`algorithm::common::database::tabular_subject_ids`). A
    /// sequence whose id LOSAT cannot write as NCBI is an error before any format is
    /// written where a tabular report `reads` the subject, else its id is empty; `None` where the
    /// lookup fails (`common::tabular::SubjectId`).
    pub(super) fn tabular_subject_ids(
        &self,
        oids: crate::seqdb::OidMap<'_>,
        reads: bool,
    ) -> Result<Vec<crate::algorithm::common::tabular::SubjectId>> {
        shared::tabular_subject_ids(
            &self.records,
            &self.ordinals,
            oids,
            &self.parsed,
            &self.sources,
            false,
            reads,
            PROGRAM,
        )
    }

    /// The reported sequences as the library's records name them (`shared::library_subjects`).
    pub(super) fn library_subjects(
        &self,
        database: &BlastnDatabase,
    ) -> Result<Vec<crate::api::local_blast::ReportedSubject>> {
        let lengths = database.sam_lengths(&self.ordinals)?;
        Ok(shared::library_subjects(
            &database.db,
            &self.records,
            &self.ordinals,
            &self.parsed,
            &self.sources,
            |rank| lengths[rank],
            PROGRAM,
        ))
    }
}

#[cfg(test)]
#[path = "database_tests.rs"]
mod tests;
