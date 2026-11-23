//! The database path of BLASTP (`-db NAME`): the search of an NCBI BLAST database read
//! in blocks through `crate::seqdb`, as NCBI BLAST+ 2.17.0 `blastp -db NAME
//! -num_threads 1` searches it (PD-LOSAT-DATABASE-SEARCH, P1 and P1.1).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the set-level numbers come from the index files (`BlastpDatabase`): the number of
//!   sequences, the total length, and 10 for the shortest sequence behind the hit-saving
//!   cutoff (`BLAST_SEQSRC_MINLENGTH`);
//! - the preliminary stage reads the sequences in blocks and searches them in the order
//!   of their ordinal ids (`search_blocks`); nothing but the hit lists is kept;
//! - the sequences of the hit lists are read again by ordinal id after the preliminary
//!   stage (`ListedSubjects`) for the redo, and the definition lines of the final hit
//!   lists after the redo (`ListedSubjects::read_titles`) for the reports;
//! - the reports show the database (`BlastpDatabase::pairwise_database`, `name`) and
//!   the subject ids and titles of sequences without parsed Seq-ids
//!   (`ListedTitles::tabular_subject_ids`, `stitles`).
//!
//! The per-subject kernel, the composition redo and the report writers are those of the
//! `-subject` path. The opening of the database, the reader's errors and the rules for the
//! stored titles are shared with BLASTN (`algorithm::common::database`).

use std::sync::Arc;

use anyhow::{bail, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::algorithm::blastp::hsp::BlastCompoHeap;
use crate::algorithm::common::database as shared;
pub(super) use crate::algorithm::common::database::{
    HeaderReadError, PreliminaryReadError, Utf8Read,
};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
#[cfg(all(test, unix))]
use crate::cli::NativeError;
use crate::report::blast_xml::HitTaxonomy;
use crate::report::db_subject::{HitLoads, ParsedSubject};
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{BlockRange, SeqDb, SeqType, SequenceBlock};
use crate::stats::search_space::SearchSpace;
use crate::stats::KarlinParams;
use crate::utils::matrix::BLASTAA_SIZE;

use super::{
    ncbistdaa_to_ascii, BlastpHitList, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch,
    EncodedProtein, PrelimSink, XPrelimHsp,
};

// EXPERIMENT (LOSAT_X_BPDBCUT): the blocks of the preliminary stage cut by the database size and
// the thread count.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_db_cut;
// EXPERIMENT (LOSAT_X_BPDBKEEP): one search scratch and one block buffer per worker for the whole
// preliminary stage.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_db_keep;
mod x_listed_lean;
mod x_redo_db_subjects;
#[cfg(test)]
pub(super) use x_listed_lean::x_assert_same_titles_and_letters;
pub(super) use x_listed_lean::x_bp_listed_lean;
pub(super) use x_redo_db_subjects::XShownResidues;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-119
/// NCBIstdaa code of the gap character (`-`), which is also the sentinel byte.
const NCBISTDAA_GAP: u8 = 0;

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

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`ResolvedBlastpArgs::db`: the `-db` value without white space at its ends).
pub(super) struct BlastpDatabase {
    db: SeqDb,
    name: String,
    /// The taxonomy database of this run, opened at the first name asked for (outfmt 15 and 16,
    /// and the batches' prefetch of outfmt 0, 5, 6 and 7).
    taxdb: shared::TaxdbRun,
    /// The program that LOSAT's refusals name: BLASTP or PSIBLAST (`BlastpApp::name`).
    program: &'static str,
    /// rpsblast: the search units of the database and their statistics, built once per run where
    /// the database opens (`rps_hooks::with_units`); `None` for every other search.
    rps_units: Option<crate::algorithm::rps_search::volumes::RpsUnits>,
}

impl BlastpDatabase {
    /// The database with rpsblast's search units `units` (`rps_hooks::with_units`).
    pub(super) fn with_rps_units(
        mut self,
        units: crate::algorithm::rps_search::volumes::RpsUnits,
    ) -> Self {
        self.rps_units = Some(units);
        self
    }

    /// rpsblast's search units (`with_rps_units`), `None` for every other search.
    pub(super) fn rps_units(&self) -> Option<&crate::algorithm::rps_search::volumes::RpsUnits> {
        self.rps_units.as_ref()
    }

    /// The taxonomy database's lookups back to a run's start (`shared::TaxdbRun::restart`), and the
    /// loads' posts of NCBI's warning of a UTF8String header (`SeqDb::restart_utf8_loads`), for the
    /// search again of the batches before a failed one.
    pub(super) fn restart_taxonomy(&self) {
        self.taxdb.restart();
        self.db.restart_utf8_loads();
    }

    /// The OIDs of the database's ordinal ids, which name the reported sequences without
    /// parsed Seq-ids (`gnl|BL_ORD_ID|n`).
    pub(super) fn oid_map(&self) -> crate::seqdb::OidMap<'_> {
        self.db.oid_map()
    }

    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(super) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    /// Opens the database of `-db name` (`algorithm::common::database::open_database`)
    /// with the command line's errors, which name `program`.
    #[cfg(test)]
    pub(super) fn open(
        name: &str,
        settings: &ApplicationSettings,
        program: &'static str,
    ) -> Result<Self> {
        Self::open_restricted(name, settings, program, None, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line, and with the id of the database mask `mask` checked (`-db_soft_mask`,
    /// `-db_hard_mask`; series NP, M4b) (`shared::restrict_database`, whose warnings go to
    /// `diagnostics`).
    pub(super) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        program: &'static str,
        list: Option<&crate::seqdb::DbList>,
        mask: Option<&crate::blastinput::db_mask::DbMaskOption>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Protein, settings, program)?;
        // A database with parsed Seq-ids (line DBID): BLASTP (ID2a) and psiblast read it.
        shared::restrict_database(&mut db, name, list, mask, diagnostics, settings, program)?;
        Ok(Self {
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
            program,
            rps_units: None,
        })
    }

    /// The opened database: rpsblast reads its volume's path and its profiles' consensus
    /// sequences (`rps_hooks`).
    pub(super) fn seqdb(&self) -> &SeqDb {
        &self.db
    }

    /// The `-db` value as NCBI keeps it.
    pub(super) fn name(&self) -> &str {
        &self.name
    }

    /// Whether `taxdb.bti` is on the search path, asked once for the run
    /// (`shared::TaxdbRun::found`): the batches' prefetch reads the taxonomy database only then
    /// (`ListedTitles::prefetch_messages`).
    pub(super) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
    }

    /// Whether a fetch of the run may be refused for NCBI's taxonomy cache
    /// (`shared::TaxdbRun::cache_may_refuse`), which can come at a later batch.
    pub(super) fn cache_may_refuse(&self) -> bool {
        self.taxdb.cache_may_refuse(&self.db)
    }

    /// The number of included sequences (`SeqDb::num_seqs`: the OIDs that the alias files'
    /// ranges and a list include, not NCBI's `GetNumSeqs`).
    pub(super) fn num_sequences(&self) -> u32 {
        self.db.num_seqs()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:145-150
    /// The shortest sequence length of the preliminary hit cutoff (`BlastSeqSrcGetMinSeqLen`,
    /// `min_subject_length` of `BLAST_GapAlignSetUp`): `CSeqDB::GetMinLength`
    /// (`SeqDb::min_length`: the alias files' `MIN_SEQ_LENGTH` or `BLAST_SEQSRC_MINLENGTH`
    /// without a list, the scanned shortest length with one, line DI-8).
    pub(super) fn cutoff_subject_length(&self) -> i32 {
        self.db.min_length()
    }

    /// The numbers the statistics read, with NCBI's types (`SeqDb::statistics`).
    pub(super) fn statistics(&self) -> crate::seqdb::DbStatistics {
        self.db.statistics()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-906
    /// The total length the statistics' setup reads (`DbStatistics::setup_total`), NCBI's `Int8`
    /// on every target: an alias file's `LENGTH` declares it, so it is not bounded by the
    /// letters LOSAT holds (2^32 or more on a 32-bit target).
    pub(super) fn total_letters(&self) -> i64 {
        self.statistics().setup_total
    }

    /// The length of the included sequences (`SeqDb::included_length`), which the schedules read.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) fn included_letters(&self) -> Result<u64, crate::seqdb::SeqDbError> {
        self.db.included_length()
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(super) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
    }

    /// The taxids of the search's list (`SeqDb::list_taxids`).
    pub(super) fn list_taxids(&self) -> &[i32] {
        self.db.list_taxids()
    }

    /// The number of sequences the statistics' setup reads (`DbStatistics::setup_count`), NCBI's
    /// `Int4`.
    pub(super) fn statistics_num_sequences(&self) -> i32 {
        self.statistics().setup_count
    }

    /// Whether the first query batch with a valid context fails with NCBI's `Database memory map
    /// file error` (`SeqDb::first_oid_read_fails`: no OID, and statistics that read OID 0).
    pub(super) fn first_oid_read_fails(&self) -> bool {
        self.db.first_oid_read_fails()
    }

    /// The error of the preliminary search's first read of the sequences after a batch's set-up
    /// checks, where a sequence file cannot be opened (`shared::sequence_file_read_error`).
    pub(super) fn sequence_file_read_error(&self, program: &str) -> Option<anyhow::Error> {
        shared::sequence_file_read_error(&self.db, program)
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`): the title (the name as given
    /// for an empty title), the date and the real totals, also with `-dbsize`.
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// `pairwise_database` of a search masked with `mask` (`None`: not masked): each name's
    /// algorithm, which the `Masked using:` lines show (`shared::pairwise_database_masked`).
    pub(super) fn pairwise_database_masked(
        &self,
        mask: Option<&crate::blastinput::db_mask::ResolvedMask>,
    ) -> PairwiseDatabase {
        shared::pairwise_database_masked(&self.db, mask.map(|mask| mask.shown.as_slice()))
    }

    /// The names of the taxonomy database that the rows of the custom tabular report `spec`
    /// look up for the subjects `taxonomy` (`shared::tabular_names`, with the run's `TaxdbRun`);
    /// `program` names NCBI's messages and LOSAT's refusals (`BLASTP`, `PSIBLAST`).
    pub(super) fn tabular_names(
        &self,
        spec: &crate::report::tabular::TabularSpec,
        taxonomy: &[crate::report::tabular::TabularTaxonomy],
        program: &str,
    ) -> Result<Option<crate::report::tabular::TabularNames>> {
        shared::tabular_names(&self.db, &self.taxdb, spec, taxonomy, program)
    }

    /// The open of the taxonomy database by the rows of a query batch (`shared::TaxdbRun::open_for_rows`).
    pub(super) fn open_for_tabular_rows(&self, program: &str) -> Result<()> {
        self.taxdb.open_for_rows(&self.db, program)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-848
/// The length adjustment and the effective search space of a query of `query_len`
/// letters against a database of `db_length` letters (the real total, or `-dbsize`) and
/// `db_num_seqs` sequences (NCBI's `Int4`), in NCBI's `Int8` (`stats::setup::calc_eff_lengths`):
/// the search space as the reports print it, and as the statistics read it (a `double`). A
/// product that does not fit `Int8`, where NCBI's arithmetic overflows, is rejected (the refusal
/// names `program`).
pub(super) fn blastp_database_search_space(
    query_len: usize,
    db_length: i64,
    db_num_seqs: i32,
    params: &KarlinParams,
    program: &str,
) -> Result<(SearchSpace, i64)> {
    let query_length = i64::try_from(query_len).unwrap_or(i64::MAX);
    let lengths = crate::stats::setup::calc_eff_lengths(
        query_length,
        db_length,
        db_num_seqs,
        params,
        0,
        program,
    )?;
    Ok((
        SearchSpace {
            effective_space: lengths.eff_searchsp as f64,
            length_adjustment: i64::from(lengths.length_adjustment),
        },
        lengths.eff_searchsp,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:494-508
// NCBI searches the stored bytes in place, O (26) included. Its lookup table indexes all 28 letters
// (blast_aalookup.c:244), and the query never holds O (the query encoder reads it as X,
// blast_setup_cxx.cpp:894-898), so a stored O meets a query word only as a neighbour: a word whose
// score with O's row and column (those of X, blast_stat.c:1583-1593) is at the threshold or more. A
// query word whose self score is under the threshold sits under its own letters only, so `AXA` (7
// with BLOSUM62, threshold 11) meets a stored `AXA` and not a stored `AOA`. The search therefore
// reads a stored O as it is: the score rows of the extension give it X's row and column, a PSSM its
// own column O, and the redo and the reports read the stored letter (`ListedSubjects`). A stored
// gap character (0, which makeblastdb writes for `-` in a protein FASTA) is searched as the residue
// it is in NCBI: row and column 0 stay `BLAST_SCORE_MIN`, so an ungapped extension stops at it and
// a gapped alignment goes round it with a gap, and NCBI's reports show it as `-`
// (`recount_shown_identities`). A byte that is no NCBIstdaa residue (28 or more) indexes NCBI's
// matrix past its rows and is rejected.
/// Checks the block `block` (the sequences from the ordinal id `range.first_ordinal`) for the
/// search, which reads its bytes as they are stored: a sequence that LOSAT cannot search as NCBI
/// does is an error that names its OID (`oids`, the database's) and `program`.
fn prepare_search_block(
    block: &SequenceBlock,
    oids: crate::seqdb::OidMap<'_>,
    program: &str,
) -> Result<()> {
    let sequences = block.len();
    let invalid = holds_non_residue_byte(block.protein_bytes());
    // A well-formed block holds one sentinel byte before its first sequence and one
    // after every sequence.
    let sentinels_in_place = (0..sequences).all(|i| {
        let with_sentinels = block.protein_with_sentinels(i);
        (i > 0 || with_sentinels[0] == NCBISTDAA_GAP)
            && with_sentinels[with_sentinels.len() - 1] == NCBISTDAA_GAP
    });
    if !invalid && sentinels_in_place {
        return Ok(());
    }
    for i in 0..sequences {
        let ordinal = block.first_ordinal() + i as u32;
        shared::check_protein_sequence(
            block.protein_with_sentinels(i),
            || oids.oid(ordinal),
            program,
        )?;
    }
    bail!("a database block is not well-formed, which is not supported by LOSAT's {program}")
}

// No NCBI counterpart: LOSAT's check of a block's bytes (`prepare_search_block`). The flag is made
// with `u8` lanes and no branch per byte, so the compiler works on 16 bytes (32 with AVX2) at a
// time; a loop with a branch per byte and a `usize` count of the gap bytes had taken 38 % of the
// cycles of a one-thread blastp-fast `-db` search (one 357-residue query against UniRef50 volume
// 00).
/// Whether a byte of `bytes` is no NCBIstdaa residue (28 or more).
fn holds_non_residue_byte(bytes: &[u8]) -> bool {
    const _: () = assert!(BLASTAA_SIZE <= u8::MAX as usize);
    const BLASTAA_SIZE_U8: u8 = BLASTAA_SIZE as u8;
    let mut invalid = 0u8;
    for &byte in bytes {
        invalid |= u8::from(byte >= BLASTAA_SIZE_U8);
    }
    invalid != 0
}

/// Reads the block `range` and prepares it for the search (`prepare_search_block`); the errors
/// name `program`.
fn read_search_block(db: &SeqDb, range: BlockRange, program: &str) -> Result<SequenceBlock> {
    let block = db
        .read_block(range)
        .map_err(|error| shared::search_read_error(error, program))?;
    prepare_search_block(&block, db.oid_map(), program)?;
    Ok(block)
}

/// Sequence `i` of a prepared block as the kernel reads it.
fn block_subject(block: &SequenceBlock, i: usize) -> BlastpSubjectRef<'_> {
    let aa_seq = block.protein_with_sentinels(i);
    BlastpSubjectRef {
        aa_seq,
        aa_len: aa_seq.len() - 2,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:160-172
/// The offset of the diagonal table after the sequences of `range` when it is
/// `diag_offset` before them: the recurrence of `BlastpSubjectScratch::finish_subject`
/// over the sequence lengths of the index. A sequence without letters is not searched
/// (blast_engine.c:1429-1432), so the offset does not move for it.
///
/// No search result depends on this value: from a cleared table, every starting offset
/// from `window` up to INT4_MAX / 4 gives the same decisions (a cleared entry reads as
/// `diff >= window`, or `diff >= 0` for one hit, and its flag is 0). A block starts at
/// the one-thread offset so that its table holds the values of the one-thread table and
/// the INT4_MAX / 4 reset (which keeps the 31-bit `last_hit` in range) falls on the
/// same sequences.
///
/// The lengths come from the index entries that the block's read reads (`Volume::read_block`):
/// when they cannot be read here, the offset is returned unchanged, and the block's read meets
/// the same error, which the round returns.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn diag_offset_after(db: &SeqDb, range: BlockRange, mut diag_offset: i32, window: i32) -> i32 {
    let Ok(offsets) = db.volumes()[range.volume]
        .index()
        .sequence_offsets(range.first, range.end)
    else {
        return diag_offset;
    };
    for bounds in offsets.windows(2) {
        // The sequence and the sentinel byte after it.
        let subject_length = (bounds[1] - bounds[0]).saturating_sub(1);
        if subject_length == 0 {
            continue;
        }
        if diag_offset >= i32::MAX / 4 {
            diag_offset = window;
        } else {
            diag_offset += subject_length as i32 + window;
        }
    }
    diag_offset
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
/// The preliminary stage of a database search: every sequence, in the order of its
/// ordinal id, through the per-subject kernel `process_subject`, with its HSP lists put
/// into `x_compact_hit_lists` in that order (NCBI's loop with one thread).
///
/// One thread: the blocks in order, one scratch (one diagonal table) for the whole
/// database, the hit lists updated after every subject, as the one-thread loop of the
/// `-subject` path.
///
/// More threads: rounds of blocks. The blocks of a round are read and searched in
/// parallel; each worker starts its block from a cleared diagonal table at the offset
/// that the one-thread order gives the block's first sequence
/// (`prepare_independent_subject`, `diag_offset_after`) and advances it inside the block
/// as the one-thread path does. A table entry left by an earlier subject is below the
/// offset of every later subject by at least the window, so it reads as a cleared entry
/// does (the property the parallel `-subject` path uses per subject). The round's
/// results are put into the hit lists in ordinal order, and its blocks are dropped
/// before the next round is read: during this stage a worker holds one block at a
/// time, so the sequence memory is bounded by the threads and the block size, not by
/// the database (the results of a round wait for the merge; after the stage,
/// `ListedSubjects` holds the sequences of the listed subjects). A block that cannot be
/// read or searched does not stop the other blocks of its round; the round's first
/// error in ordinal order is returned, the error of the one-thread order (the blocks
/// are cut alike). The hit lists hold the compact element (`x_prelim_compact.rs`).
/// EXPERIMENT (LOSAT_X_BPDBCUT): with more than one thread the blocks are cut by the
/// database size and the thread count, and a failed search returns the one-thread order's
/// error (`database/x_db_cut.rs`, exactness argument in its header).
#[allow(clippy::too_many_arguments)]
pub(super) fn search_blocks<S, K, P>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
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
    // EXPERIMENT (LOSAT_X_BPDBCUT): `database/x_db_cut.rs`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        if x_db_cut::x_bp_db_cut() {
            return x_db_cut::x_search_blocks_cut(
                database,
                pool,
                window,
                blocks_per_thread,
                new_scratch,
                process_subject,
                sink,
            );
        }
    }
    search_blocks_of(
        database,
        pool,
        window,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        process_subject,
        sink,
    )
}

/// `search_blocks` with blocks of at most `block_bytes` and rounds of `blocks_per_thread`
/// blocks per thread.
#[allow(clippy::too_many_arguments)]
fn search_blocks_of<S, K, P>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
{
    // EXPERIMENT (LOSAT_X_BPDBKEEP): `database/x_db_keep.rs`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_db_keep::x_bp_db_keep() {
        return x_db_keep::x_search_blocks_kept(
            database,
            pool,
            window,
            block_bytes,
            blocks_per_thread,
            new_scratch,
            process_subject,
            sink,
        );
    }
    search_blocks_made(
        database,
        pool,
        window,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        process_subject,
        sink,
    )
}

/// `search_blocks_of` without LOSAT_X_BPDBKEEP: each task makes its scratch and reads its block
/// into a new buffer.
#[allow(clippy::too_many_arguments)]
fn search_blocks_made<S, K, P>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink,
{
    let db = &database.db;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if super::blastp_should_parallelize_stage(
        "prelim",
        pool.threads(),
        pool.threads(),
        db.num_seqs() as usize,
        db.num_seqs() as usize,
    ) && pool.enabled()
    {
        let round_blocks = pool.threads().saturating_mul(blocks_per_thread).max(1);
        let mut next_ordinal = 0u32;
        let mut diag_offset = window;
        loop {
            // The blocks of the round and the diagonal offset of each one's first
            // sequence in the one-thread order.
            let mut round: Vec<(BlockRange, i32)> = Vec::with_capacity(round_blocks);
            while round.len() < round_blocks {
                let Some(range) = db.next_block(next_ordinal, block_bytes) else {
                    break;
                };
                round.push((range, diag_offset));
                diag_offset = diag_offset_after(db, range, diag_offset, window);
                next_ordinal = range.end_ordinal();
            }
            if round.is_empty() {
                break;
            }
            let found: Vec<Result<Vec<BlastpSubjectResult>>> = pool.install(|| {
                round
                    .par_iter()
                    .map_init(
                        new_scratch,
                        |scratch,
                         &(range, first_diag_offset)|
                         -> Result<Vec<BlastpSubjectResult>> {
                            let block = read_search_block(db, range, database.program)?;
                            scratch.prepare_independent_subject(first_diag_offset);
                            let mut with_hits = Vec::new();
                            let mut subject_result = BlastpSubjectResult::new(0);
                            for i in 0..block.len() {
                                let s_idx = range.first_ordinal as usize + i;
                                process_subject(
                                    scratch,
                                    s_idx,
                                    block_subject(&block, i),
                                    &mut subject_result,
                                );
                                if !subject_result.x_compact_hsp_lists.is_empty() {
                                    with_hits.push(std::mem::replace(
                                        &mut subject_result,
                                        BlastpSubjectResult::new(0),
                                    ));
                                }
                            }
                            Ok(with_hits)
                        },
                    )
                    .collect()
            });
            // The subjects without HSP lists change no hit list, so the order of the
            // updates is the one-thread order.
            for block_results in found {
                for mut subject_result in block_results? {
                    sink.add(&mut subject_result);
                }
            }
        }
        return Ok(());
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, window, blocks_per_thread);

    let mut scratch = new_scratch();
    let mut subject_result = BlastpSubjectResult::new(0);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block = read_search_block(db, range, database.program)?;
        for i in 0..block.len() {
            process_subject(
                &mut scratch,
                range.first_ordinal as usize + i,
                block_subject(&block, i),
                &mut subject_result,
            );
            sink.add(&mut subject_result);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

/// The database sequences that the preliminary hit lists name, read by ordinal id after
/// the preliminary stage: they stand for the subject set in the composition redo. Their
/// definition lines are read after the redo, for the final hit lists only
/// (`read_titles`).
///
/// From here on a hit list names a sequence by its rank among the listed ordinal ids
/// (`ordinals[rank]`), so that the code of the `-subject` path indexes these vectors as
/// it indexes the encoded subject set. The rank grows with the ordinal id, so every
/// comparison that breaks a tie by the ordinal id (`s_EvalueCompareHSPLists`,
/// `ScoreCompareHSPs`; `hsp.rs`) gives the same order.
pub(super) struct ListedSubjects {
    /// The listed ordinal ids, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The sequences as stored (O stays O), between two sentinel bytes, by rank: empty
    /// until the redo has run (it reads each match's subject from the database, `x_list`), and
    /// then what the reports read (`x_read_shown`): a subject whose residues no report reads holds
    /// its length only (`aa_seq` empty).
    pub(super) encoded: Vec<EncodedProtein>,
    /// The number of residues of each subject, by rank (`x_list`).
    x_lengths: Vec<u32>,
    /// Whether a read of the redo or of the first Newton problems found a pyrrolysine (O) or a
    /// gap character, by rank (`x_list`).
    x_recounted: Vec<std::sync::atomic::AtomicBool>,
    /// The subjects of the preliminary hit lists whose header holds a UTF8String where a
    /// VisibleString is expected, which NCBI's fetches with `check_oid_exclusion` read
    /// (`exclude_fetched`; `note_fetched`).
    fetched_utf8: XFetchedUtf8,
    /// The call's reads before those of its reported sequences that took a UTF8String for a
    /// VisibleString: of its fetches (`note_fetched`) and of its Seq-aligns (`drop_unaligned`).
    earlier_utf8: std::sync::Mutex<shared::EarlierUtf8>,
    /// The places of the lists that the redo reaches and its fetch excludes, which left the
    /// preliminary hit lists (`exclude_fetched`), by query and place.
    pub(super) fetch_excluded: Vec<FetchExcludedMatch>,
}

/// A preliminary hit list that NCBI's composition-based redo reaches in its stream and whose
/// subject the redo's fetch then excludes (`exclude_fetched`): the list leaves the query's hit
/// list, and this mark keeps its place, `position` kept lists of the query before it, with what
/// the redo reads of it there (its best e-value, the early termination test's). The
/// Smith-Waterman redo's status reads the marks, in the one-thread loop
/// (`redo_run::one_thread_redo`) and in the replay of the parallel redo's windows
/// (`x_redo_par.rs`, `XSwBatch`), both through `redo_run::query_fetch_excluded_reached`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(in crate::algorithm::blastp) struct FetchExcludedMatch {
    pub(in crate::algorithm::blastp) query: usize,
    pub(in crate::algorithm::blastp) position: usize,
    pub(in crate::algorithm::blastp) best_evalue: f64,
}

/// What `exclude_fetched` keeps of the subjects whose header holds a UTF8String where a
/// VisibleString is expected (`reads`) or cannot be read (`unreadable`, whose lists stay:
/// `shared::FetchedReads::unreadable`): each with its query and its place in the order in which
/// the redo takes its query's lists (`XRedoKey`), and, for each query with one, the ordinal ids
/// of its kept lists in that order (the redo's positions of its heap's records, `note_fetched`).
#[derive(Debug, Default)]
struct XFetchedUtf8 {
    reads: Vec<XFetchedRead>,
    unreadable: Vec<XFetchedRead>,
    orders: Vec<(usize, Vec<u32>)>,
}

/// One subject of `XFetchedUtf8`: the query of its list, its key in the redo's order, and the
/// kept lists of its query that come before it there (`shared::redo_fetches`).
#[derive(Clone, Copy, Debug)]
struct XFetchedRead {
    query: usize,
    key: XRedoKey,
    kept: bool,
    before: usize,
}

/// What `compare_hsp_lists` (`Blast_HitListSortByEvalue`, the order in which the redo takes a
/// query's lists, cited at `BlastpHitList::sort_by_evalue`) reads of a list: lists without HSPs
/// last, then the best E-value (`s_EvalueComp`), the first HSP's score, larger first, and the
/// subject, larger first (the ordinal id: the ranks that the lists hold after `x_list` keep its
/// order).
#[derive(Clone, Copy, Debug, PartialEq)]
struct XRedoKey {
    has_hsps: bool,
    best_evalue: f64,
    score: i32,
    ordinal: u32,
}

impl XRedoKey {
    fn of<H: crate::algorithm::blastp::hsp::HitListHsp>(
        list: &crate::algorithm::blastp::hsp::BlastpHspList<H>,
    ) -> Self {
        XRedoKey {
            has_hsps: !list.hsps.is_empty(),
            best_evalue: list.best_evalue,
            score: list.hsps.first().map_or(0, |hsp| hsp.raw_score()),
            ordinal: list.oid,
        }
    }

    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .has_hsps
            .cmp(&self.has_hsps)
            .then_with(|| {
                crate::algorithm::blastp::hsp::evalue_comp(self.best_evalue, other.best_evalue)
            })
            .then_with(|| other.score.cmp(&self.score))
            .then_with(|| other.ordinal.cmp(&self.ordinal))
    }
}

impl XFetchedUtf8 {
    /// Keeps the subjects `found` of the list of `query` (whose kept lists are `kept`).
    #[cold]
    #[inline(never)]
    fn add<H: crate::algorithm::blastp::hsp::HitListHsp>(
        &mut self,
        query: usize,
        kept: &[crate::algorithm::blastp::hsp::BlastpHspList<H>],
        found: shared::FetchedReads<XRedoKey>,
    ) {
        let mut keys: Vec<XRedoKey> = kept.iter().map(XRedoKey::of).collect();
        keys.sort_unstable_by(XRedoKey::cmp);
        let place = |(read, key): (shared::FetchedUtf8, XRedoKey)| XFetchedRead {
            query,
            key,
            kept: read.kept,
            before: keys.partition_point(|other| other.cmp(&key).is_lt()),
        };
        self.reads.extend(found.utf8.into_iter().map(place));
        self.unreadable
            .extend(found.unreadable.into_iter().map(place));
        self.orders
            .push((query, keys.iter().map(|key| key.ordinal).collect()));
    }

    /// The first, in NCBI's order of the fetches, of the reads `reads` (`self.reads` or
    /// `self.unreadable`) that the ordinary traceback (`shared::FetchedUtf8::traceback_first`:
    /// every subject, in increasing ordinal) or, with `heaps` (as the redo left them), the
    /// composition-based redo makes (`shared::redo_fetches`: query by query, each in its lists'
    /// order), the batches of the queries (`batch_of`) one after the other: its query and ordinal
    /// id.
    #[cold]
    #[inline(never)]
    fn first(
        &self,
        reads: &[XFetchedRead],
        heaps: Option<&[BlastCompoHeap]>,
        batch_of: &dyn Fn(usize) -> std::ops::Range<usize>,
        ordinals: &[u32],
    ) -> Option<(usize, u32)> {
        let Some(heaps) = heaps else {
            // A batch's traceback, then the next batch's.
            return reads
                .iter()
                .min_by_key(|read| (batch_of(read.query).start, read.key.ordinal))
                .map(|read| (read.query, read.key.ordinal));
        };
        self.first_redone(reads, heaps, batch_of, ordinals, &|_| true)
    }

    /// `first` of `reads` in the composition-based redo with `heaps`, among the reads that
    /// `reached` keeps (the fetches before the failure of a batch's redo,
    /// `ListedSubjects::fetched_before_failure`).
    fn first_redone(
        &self,
        reads: &[XFetchedRead],
        heaps: &[BlastCompoHeap],
        batch_of: &dyn Fn(usize) -> std::ops::Range<usize>,
        ordinals: &[u32],
        reached: &dyn Fn(&XFetchedRead) -> bool,
    ) -> Option<(usize, u32)> {
        // The heap of `query` at the turn of a match after `before` of its lists (`order`), as
        // `shared::redo_fetches` reads it.
        let heap_at = |query: usize, before: Option<(&[u32], usize)>| {
            let heap = &heaps[query];
            let inserted = |record: &&crate::algorithm::blastp::hsp::BlastCompoHeapRecord| {
                before.is_none_or(|(order, before)| {
                    let ordinal = ordinals[record.subject_index as usize];
                    order[..before.min(order.len())].contains(&ordinal)
                })
            };
            shared::RedoHeap {
                ecutoff: heap.ecutoff,
                heap_threshold: heap.heap_threshold,
                inserted_before: heap
                    .records()
                    .iter()
                    .filter(|record| record.best_evalue <= heap.ecutoff)
                    .filter(inserted)
                    .count(),
                ends_early: heap.x_ends_early,
            }
        };
        reads
            .iter()
            .filter(|read| read.key.has_hsps && reached(read))
            .filter(|read| {
                let order = &self
                    .orders
                    .iter()
                    .find(|(of, _)| *of == read.query)
                    .expect("a query with a read has its order")
                    .1;
                // The heaps of the queries before the match's are as the redo left them, those
                // after it untouched (the redo takes the lists query by query).
                let redo_heaps: Vec<shared::RedoHeap> = batch_of(read.query)
                    .filter(|&query| query <= read.query)
                    .map(|query| {
                        heap_at(
                            query,
                            (query == read.query).then_some((order.as_slice(), read.before)),
                        )
                    })
                    .chain(
                        batch_of(read.query)
                            .filter(|&query| query > read.query)
                            .map(|query| shared::RedoHeap {
                                inserted_before: 0,
                                ..heap_at(query, None)
                            }),
                    )
                    .collect();
                shared::redo_fetches(read.key.best_evalue, &redo_heaps)
            })
            .min_by(|a, b| {
                a.query
                    .cmp(&b.query)
                    .then(a.before.cmp(&b.before))
                    .then(a.kept.cmp(&b.kept))
                    .then_with(|| a.key.cmp(&b.key))
            })
            .map(|read| (read.query, read.key.ordinal))
    }
}

/// The notes of the subjects that a call's fetches with `check_oid_exclusion` read (a seqid
/// list's, `ListedSubjects::note_fetched`) around its traceback: the PSSM's position-based redo
/// notes those of its own heaps (`pssm_redo`, `pssm_hooks::pssm_traceback`); after the traceback
/// (`after`) every subject of the ordinary traceback, or those that the composition-based redo's
/// early termination test leaves (its heaps as it left them). The replay of the parallel redo's
/// windows (`x_redo_par.rs` `XSwBatch`) reads `refused_before_failure` and `fetched_before_failure`
/// through the one-thread loop's errors (`redo_run::sw_batch_failed_after`,
/// `redo_run::sw_endless_after`).
pub(in crate::algorithm::blastp) struct FetchNotes<'a> {
    listed: Option<&'a ListedSubjects>,
    database: Option<&'a BlastpDatabase>,
    /// Whether the PSSM's position-based redo ran and noted its own.
    pssm_redone: std::cell::Cell<bool>,
}

impl<'a> FetchNotes<'a> {
    pub(super) fn new(
        listed: Option<&'a ListedSubjects>,
        database: Option<&'a BlastpDatabase>,
    ) -> Self {
        Self {
            listed,
            database,
            pssm_redone: std::cell::Cell::new(false),
        }
    }

    /// The notes of a search with neither a seqid list nor a database: none to read.
    #[cfg(all(
        test,
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(in crate::algorithm::blastp) fn none() -> FetchNotes<'static> {
        FetchNotes::new(None, None)
    }

    /// The PSSM's position-based redo of modes 1 to 3 ran, whose heaps `heaps` tell which
    /// subjects its fetches read; LOSAT's refusal where one of them has a header that cannot be
    /// read (`ListedSubjects::note_fetched`).
    pub(super) fn pssm_redo(&self, heaps: &[BlastCompoHeap]) -> Result<()> {
        self.pssm_redone.set(true);
        if let (Some(listed), Some(database)) = (self.listed, self.database) {
            let all = |_: usize| 0..heaps.len();
            listed.note_fetched(database, Some(heaps), &all)?;
        }
        Ok(())
    }

    /// LOSAT's refusal where the composition-based redo of the batch `batch`, before its
    /// Smith-Waterman redo failed or did not end at the turn of the list after `before` of the
    /// kept lists of query `query` (as `fetched_before_failure`), with `heaps` as the redo left
    /// them, fetched a subject whose header cannot be read: NCBI's search ended there
    /// (`shared::FetchExcluded::refusal_at`).
    pub(super) fn refused_before_failure(
        &self,
        heaps: &[BlastCompoHeap],
        batch: std::ops::Range<usize>,
        query: usize,
        before: usize,
    ) -> Option<anyhow::Error> {
        let (listed, database) = (self.listed?, self.database?);
        if listed.fetched_utf8.unreadable.is_empty() {
            return None;
        }
        let batch_of = |_: usize| batch.clone();
        let reached = |read: &XFetchedRead| {
            batch.contains(&read.query)
                && (read.query < query || (read.query == query && read.before <= before))
        };
        let (query, ordinal) = listed.fetched_utf8.first_redone(
            &listed.fetched_utf8.unreadable,
            heaps,
            &batch_of,
            &listed.ordinals,
            &reached,
        )?;
        Some(fetch_refused(database, query, ordinal))
    }

    /// The first fetch of the composition-based redo of the batch `batch` that read a header with a
    /// UTF8String where a VisibleString is expected, before the batch's Smith-Waterman redo
    /// failed (`redo_run::sw_batch_failed`) or did not end (`redo_run::sw_endless`) at the turn of
    /// the list after `before` of the kept lists of query `query` (`usize::MAX`: after the batch's
    /// last list), with `heaps` as the redo left them: NCBI posted its warning before its error,
    /// its crash or its endless redo (`HeaderReadError::utf8`).
    pub(in crate::algorithm::blastp) fn fetched_before_failure(
        &self,
        heaps: &[BlastCompoHeap],
        batch: std::ops::Range<usize>,
        query: usize,
        before: usize,
    ) -> Option<Utf8Read> {
        let (listed, database) = (self.listed?, self.database?);
        if listed.fetched_utf8.reads.is_empty() {
            return None;
        }
        let batch_of = |_: usize| batch.clone();
        let reached = |read: &XFetchedRead| {
            batch.contains(&read.query)
                && (read.query < query || (read.query == query && read.before <= before))
        };
        let (_, ordinal) = listed.fetched_utf8.first_redone(
            &listed.fetched_utf8.reads,
            heaps,
            &batch_of,
            &listed.ordinals,
            &reached,
        )?;
        let at = database.db.location_of(ordinal)?;
        Some(Utf8Read {
            query: batch.start,
            path: database.db.utf8_string_path_at(at)?,
        })
    }

    /// After the traceback: `heaps`, the composition-based redo's, which made the hit lists where
    /// the ordinary traceback made none (`ordinary`), the batch of each query being in `batches`
    /// (`redo_run::kappa_query_batches`); LOSAT's refusal where a fetch read a header that
    /// cannot be read (`ListedSubjects::note_fetched`).
    pub(super) fn after<T>(
        &self,
        heaps: &[BlastCompoHeap],
        ordinary: &Option<T>,
        batches: &[std::ops::Range<usize>],
    ) -> Result<()> {
        let (Some(listed), Some(database)) = (self.listed, self.database) else {
            return Ok(());
        };
        let batch_of = |q_idx: usize| batches[q_idx].clone();
        if ordinary.is_none() {
            listed.note_fetched(database, Some(heaps), &batch_of)?;
        } else if !self.pssm_redone.get() {
            listed.note_fetched(database, None, &batch_of)?;
        }
        Ok(())
    }
}

/// LOSAT's refusal of a search whose redo or traceback fetches the sequence `ordinal`, whose
/// header cannot be read (`shared::FetchExcluded::refusal_at`), from the list of query `query`:
/// NCBI's search ends in that query's batch, after the reports of the batches before, so the
/// refusal fails that batch (`HeaderReadError`, whose batches before are searched again and
/// reported), as `-use_sw_tback`'s failed batch does (the Owner's form of 2026-11-16 for
/// line RG's Q-X2c).
#[cold]
#[inline(never)]
fn fetch_refused(database: &BlastpDatabase, query: usize, ordinal: u32) -> anyhow::Error {
    HeaderReadError {
        query,
        error: shared::FetchExcluded::refusal_at(&database.db, ordinal, database.program),
        utf8: None,
    }
    .into()
}

/// The reads of a call's reported subjects' titles (`titles`), one batch after the other, and
/// their first read of a header with a UTF8String where a VisibleString is expected, the call's
/// (`utf8`, `shared::CallUtf8`).
pub(super) struct CallTitles<'a> {
    database: Option<&'a BlastpDatabase>,
    /// The call's query batches, each with its searched queries.
    batches: &'a [crate::blastinput::seq_range::RangedBatch],
    utf8: shared::CallUtf8,
}

impl<'a> CallTitles<'a> {
    pub(super) fn new(
        database: Option<&'a BlastpDatabase>,
        batches: &'a [crate::blastinput::seq_range::RangedBatch],
    ) -> Self {
        Self {
            database,
            batches,
            utf8: shared::CallUtf8::default(),
        }
    }

    /// The titles of the call's reported database subjects, the `listed` ones
    /// (`ListedSubjects::read_titles`, or with LOSAT_X_BPLISTEDLEAN `x_lean`,
    /// `x_read_titles_lean`), read one batch after the other: the first queries of the call's
    /// batches (`shared::HeaderReader::batch_starts`). Their `utf8` is NCBI's first read, in the
    /// Seq-aligns of the call's batches, of a header with a UTF8String where a VisibleString is
    /// expected, the call's (`utf8`): its warning comes before anything the batch writes after its
    /// search (`header_read`), and before an error of the batch's prefetch
    /// (`HeaderReadError::after_reads`).
    pub(super) fn titles<'l>(
        &self,
        listed: Option<&'l ListedSubjects>,
        hit_lists: &[Option<BlastpHitList>],
        x_lean: bool,
    ) -> Result<Option<ListedTitles<'l>>>
    where
        'a: 'l,
    {
        let (Some(listed), Some(database)) = (listed, self.database) else {
            return Ok(None);
        };
        let starts: Vec<usize> = self
            .batches
            .iter()
            .map(|batch| batch.searched.start)
            .collect();
        let titles = if x_lean {
            listed.x_read_titles_lean(database, hit_lists, &starts)?
        } else {
            listed.read_titles(database, hit_lists, &starts)?
        };
        if let Some(read) = &titles.utf8 {
            let batch = self
                .batches
                .iter()
                .find(|batch| batch.searched.contains(&read.query))
                .map_or(read.query, |batch| batch.searched.start);
            self.utf8.note(batch, Some(read));
        }
        Ok(Some(titles))
    }

    /// The call's first read of such a header in its batches' Seq-aligns (`titles`).
    pub(super) fn utf8(&self) -> Option<&Utf8Read> {
        self.utf8.read()
    }

    /// NCBI's warning of that read in the first query of the batch whose Seq-aligns read it, with
    /// the posts that follow the reader's messages of the batch (`query_warnings::query_posts`),
    /// where the run has not posted it (a load of the model of the loads may have,
    /// `shared::load_read_texts`), with the post prefix `program`.
    pub(super) fn header_read(
        &self,
        program: &str,
    ) -> Option<(usize, crate::report::diag_log::Posts)> {
        self.utf8.post(
            self.database.map(BlastpDatabase::seqdb),
            program,
            usize::MAX,
        )
    }
}

impl ListedSubjects {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1381-1385
    // NCBI's redo fetches each sequence of the hit lists again from the database, in the
    // encoding of the search: the stored bytes.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1620-1636
    // Here the stored letter matters: a stored O (26) is not the X (21) that the query
    // holds for its own O (blast_setup_cxx.cpp:894-898), so an O aligned to an O counts
    // against "near identical" and the subject is then SEG-filtered, where the same
    // records given as `-subject` (O read as X on both sides) are near identical. Shown
    // with the oracle (campaign dbsearch-20261012, `$BUILD_ROOT/dbsearch/db3/exp_o_near`;
    // fixture rows `db.o_stored.*`): NCBI `-db` gives 266 bits where `-subject` gives
    // 307. So the redo reads the stored bytes; its matrices score O as X
    // (`adjust_scores.rs`, composition_adjustment.c:681-703). The preliminary stage reads the
    // stored bytes too (`prepare_search_block`).
    /// The listed subjects of a search with every residue read, as the reports of a search
    /// that makes records read them after the redo (`x_list`, then `x_read_shown` with
    /// `XShownResidues::Listed`); the unit tests of the titles start here. The searches
    /// themselves read each match's subject from the database in the redo (`x_list`).
    #[cfg(test)]
    pub(super) fn read(
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
        x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Self> {
        let mut listed = Self::x_list(database, hit_lists, x_compact_hit_lists)?;
        listed.x_read_shown(database, &[], XShownResidues::Listed, &[])?;
        Ok(listed)
    }

    /// The final hit lists of the call (subjects named by rank) without the subjects that
    /// NCBI makes no Seq-align of (`shared::SeqAligned`): right after the cut to the hit list
    /// size, before the reads of the reported subjects.
    pub(super) fn drop_unaligned(
        &self,
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
    ) {
        // blastp and psiblast fetch every subject of the redo or the traceback with
        // `check_oid_exclusion` (`exclude_fetched`).
        let Some(mut aligned) = shared::SeqAligned::of(&database.db, true) else {
            return;
        };
        let mut first = None;
        for (query, hit_list) in hit_lists.iter_mut().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            hit_list.hsplist_array.truncate(hit_list.hsplist_count);
            aligned.retain(&mut hit_list.hsplist_array, |list| {
                Some(self.ordinals[list.oid as usize])
            });
            hit_list.hsplist_count = hit_list.hsplist_array.len();
            if first.is_none() {
                first = aligned.first_utf8().map(|at| (query, at));
            }
        }
        self.earlier_utf8
            .lock()
            .expect("the UTF8String reads are not poisoned")
            .aligned = first;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3575 (cited at
    // `exclude_fetched`) and blast_traceback.c:1567-1577
    /// Notes the call's first fetch with `check_oid_exclusion` that read a header with a
    /// UTF8String where a VisibleString is expected (`XFetchedUtf8::first`): of the ordinary
    /// traceback, or with `heaps` of the composition-based redo (its heaps as it left them, by
    /// query), the batch of each query being `batch_of` (`redo_run::kappa_query_batches`). NCBI
    /// posts its warning there, before the reads of the batch's Seq-aligns (`read_titles`).
    /// Where these fetches read a header that cannot be read (`XFetchedUtf8::unreadable`), NCBI's
    /// search ends there: LOSAT's refusal (`shared::FetchExcluded::refusal_at`), as the failure of
    /// the batch of the query whose list it is (`fetch_refused`).
    pub(super) fn note_fetched(
        &self,
        database: &BlastpDatabase,
        heaps: Option<&[BlastCompoHeap]>,
        batch_of: &dyn Fn(usize) -> std::ops::Range<usize>,
    ) -> Result<()> {
        if !self.fetched_utf8.unreadable.is_empty() {
            if let Some((query, ordinal)) = self.fetched_utf8.first(
                &self.fetched_utf8.unreadable,
                heaps,
                batch_of,
                &self.ordinals,
            ) {
                return Err(fetch_refused(database, query, ordinal));
            }
        }
        if self.fetched_utf8.reads.is_empty() {
            return Ok(());
        }
        if let Some((query, ordinal)) =
            self.fetched_utf8
                .first(&self.fetched_utf8.reads, heaps, batch_of, &self.ordinals)
        {
            self.earlier_utf8
                .lock()
                .expect("the UTF8String reads are not poisoned")
                .fetched = database.db.location_of(ordinal).map(|at| (query, at));
        }
        Ok(())
    }

    /// The reader of the reported sequences' headers of the call (`read_titles`), whose batches
    /// start at the queries `batch_starts`, with the call's reads before them that took a
    /// UTF8String for a VisibleString (`note_fetched`, `drop_unaligned`).
    fn header_reader<'a>(
        &self,
        database: &'a BlastpDatabase,
        batch_starts: &[usize],
    ) -> shared::HeaderReader<'a> {
        let headers = shared::HeaderReader::new(&database.db);
        headers.batch_starts(batch_starts.to_vec());
        headers.note_earlier(
            *self
                .earlier_utf8
                .lock()
                .expect("the UTF8String reads are not poisoned"),
        );
        headers
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1582
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310
    // NCBI first reads a header file (and opens it, c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:127-133) when the
    // traceback's results become Seq-aligns: after the batch's traceback, for the final
    // hit lists only, query by query and in the order of each hit list. A header file
    // that is missing or cannot be read fails the batch there (`BLAST Database error`,
    // exit status 2, blast_app_util.hpp:185-203, outside the search stage's handler of
    // local_blast.cpp:240-251), after the outfmt 0 prolog and the reports of the batches
    // before; a batch without hits reads no header. Shown with the oracle (campaign
    // dbsearch-20261012, `scripts/db3_nophr_probe.py`; comparison cases `cli.path.no_phr*`).
    //
    // NCBI reference for the titles:
    // c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:970-1017 (`s_GetBioseqTitle`: the
    // title of the first definition line, the others appended as ` >id title`).
    /// Reads the definition lines of the subjects of the final hit lists (`hit_lists`,
    /// after the redo and before anything is dropped from them), in NCBI's order. A reader
    /// error of NCBI's own is a `HeaderReadError` that names the query whose hit list
    /// reached it; LOSAT's refusals are errors of their own.
    pub(super) fn read_titles<'a>(
        &'a self,
        database: &'a BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        batch_starts: &[usize],
    ) -> Result<ListedTitles<'a>> {
        let mut titles: Vec<Option<shared::StoredTitle>> = vec![None; self.ordinals.len()];
        let mut parsed = ParsedSubjects::default();
        let headers = self.header_reader(database, batch_starts);
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            for hsp_list in &hit_list.hsplist_array[..count] {
                let rank = hsp_list.oid as usize;
                if titles[rank].is_some() {
                    continue;
                }
                let ordinal = self.ordinals[rank];
                titles[rank] = Some(parsed.read(
                    database,
                    &headers,
                    rank,
                    ordinal,
                    query,
                    self.ordinals.len(),
                )?);
            }
        }
        let mut records = Vec::with_capacity(self.ordinals.len());
        let mut replaced_title_bytes = Vec::with_capacity(self.ordinals.len());
        let mut sources = Vec::new();
        for (rank, ((&ordinal, encoded), title)) in self
            .ordinals
            .iter()
            .zip(&self.encoded)
            .zip(titles)
            .enumerate()
        {
            let letters: Vec<u8> = encoded.aa_seq[1..=encoded.aa_len]
                .iter()
                .map(|&residue| ncbistdaa_to_ascii(residue) as u8)
                .collect();
            // A subject that is in no final hit list keeps an empty title: no report
            // shows it.
            let stored = title.unwrap_or_default();
            replaced_title_bytes.push(stored.replaced_bytes);
            shared::keep_source(&mut sources, rank, self.ordinals.len(), stored.source);
            records.push(FastaRecord::new(
                // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the
                // OID (line DI-8).
                database.db.oid(ordinal).fasta_id(),
                &stored.title,
                &letters,
            ));
        }
        Ok(ListedTitles {
            listed: self,
            oids: database.oid_map(),
            shown: ListedShown::Records(records),
            replaced_title_bytes,
            sources,
            parsed: parsed.by_rank,
            program: database.program,
            utf8: headers.first_utf8(),
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3560-3575
// `s_MatchingSequenceInitialize` fetches with `check_oid_exclusion` (blast_kappa.c:1378), as
// the ordinary traceback does (blast_traceback.c:1567-1577): a subject that the negative seqid
// list keeps without a line (`shared::FetchExcluded`) leaves the redo or the traceback before
// it is aligned, after the early termination test of its match (which reads only the heaps,
// where it never enters), so it takes no place in the final hit lists. LOSAT removes it from
// the preliminary hit lists before they are listed (after the trim of `ADAPTIVE_CBS` and the
// sort into the redo's order); blastp and psiblast are gapped here (`-ungapped` is refused).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3540
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:306-314
// The redo's loop reaches such a list in its stream: when the early termination test holds it
// frees the list and goes on, and otherwise the excluded fetch sets the status back to 0
// (blast_kappa.c:3560-3575). With `-use_sw_tback` the status of the batch is the last reached
// match's (`redo_run::one_thread_redo`), so the list's place counts: each list that the fetch
// excludes (`shared::FetchMarks`) is removed and leaves a mark (`FetchExcludedMatch`) with its
// place in its query's list, in the order the redo reads them (the lists are sorted before they
// are listed). Shown on ws with NCBI 2.17.0 (rows
// `db.x4excl.*`: an excluded match after a failed Smith-Waterman redo ends the batch with
// status 0).
/// The preliminary hit lists of a call (subjects named by ordinal id, each list in the redo's
/// order) without the subjects that NCBI's fetches of the redo or the traceback exclude, in
/// their order, and the marks of the removed lists, by query and place.
/// The subjects whose header holds a UTF8String where a VisibleString is expected or cannot be
/// read come with them (`XFetchedUtf8`): NCBI's fetches read them, and the lists of the latter
/// stay; the call refuses the search where a fetch reaches one (`ListedSubjects::note_fetched`).
fn exclude_fetched(
    database: &BlastpDatabase,
    hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) -> (Vec<FetchExcludedMatch>, XFetchedUtf8) {
    let mut fetched = XFetchedUtf8::default();
    let excluded = shared::FetchExcluded::of(&database.db);
    // Without a negative seqid list the fetches exclude nothing: no list is read here.
    if excluded.is_none() {
        return (Vec::new(), fetched);
    }
    let mut marks = Vec::new();
    for (query, hit_list) in hit_lists.iter_mut().enumerate() {
        let Some(hit_list) = hit_list else {
            continue;
        };
        // The redo reads a query's compact list where it has one (`x_expand_hit_list`).
        let read = x_compact_hit_lists
            .get(query)
            .is_none_or(|compact| compact.is_none());
        leave_fetched(
            excluded.as_ref(),
            query,
            hit_list,
            read.then_some(&mut marks),
            &mut fetched,
        );
    }
    for (query, hit_list) in x_compact_hit_lists.iter_mut().enumerate() {
        if let Some(hit_list) = hit_list {
            leave_fetched(
                excluded.as_ref(),
                query,
                hit_list,
                Some(&mut marks),
                &mut fetched,
            );
        }
    }
    // By query; the marks of a query are in place order already (a stable sort).
    marks.sort_by_key(|mark| mark.query);
    (marks, fetched)
}

/// `hit_list` (query `query`'s) without the lists whose subject the fetch excludes
/// (`shared::FetchMarks`), in their order, each removed list leaving its mark in `marks` (where
/// the redo reads this list) with the number of kept lists before it; nothing changes where
/// none is excluded. The lists whose header holds a UTF8String where a VisibleString is expected
/// or cannot be read go to `fetched` (read in the same pass, `shared::FetchExcluded::marks_noting`;
/// the latter stay in `hit_list`).
///
/// No NCBI counterpart: NCBI keeps these lists in the redo's stream, and its fetch excludes each
/// subject when the redo reaches it (the reference lines are at `exclude_fetched`).
fn leave_fetched<H: crate::algorithm::blastp::hsp::HitListHsp>(
    excluded: Option<&shared::FetchExcluded<'_>>,
    query: usize,
    hit_list: &mut BlastpHitList<H>,
    mut marks: Option<&mut Vec<FetchExcludedMatch>>,
    fetched: &mut XFetchedUtf8,
) {
    let lists = &hit_list.hsplist_array[..hit_list.hsplist_count];
    let (fetch, found) =
        shared::FetchMarks::of_noting(excluded, lists, |list| Some(list.oid), XRedoKey::of);
    let mut kept = 0;
    let mut removed = false;
    for (index, list) in lists.iter().enumerate() {
        if !fetch.excluded(index) {
            kept += 1;
            continue;
        }
        removed = true;
        if let Some(marks) = marks.as_deref_mut() {
            marks.push(FetchExcludedMatch {
                query,
                position: kept,
                best_evalue: list.best_evalue,
            });
        }
    }
    if removed {
        hit_list.hsplist_array.truncate(hit_list.hsplist_count);
        let mut index = 0;
        hit_list.hsplist_array.retain(|_| {
            index += 1;
            !fetch.excluded(index - 1)
        });
        hit_list.hsplist_count = hit_list.hsplist_array.len();
    }
    if !found.utf8.is_empty() || !found.unreadable.is_empty() {
        fetched.add(
            query,
            &hit_list.hsplist_array[..hit_list.hsplist_count],
            found,
        );
    }
}

/// The ordinal ids of the subjects of `hit_lists` and `x_compact_hit_lists`, ascending, each once.
fn listed_ordinals(
    hit_lists: &[Option<BlastpHitList>],
    x_compact_hit_lists: &[Option<BlastpHitList<XPrelimHsp>>],
) -> Vec<u32> {
    let mut ordinals: Vec<u32> = hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| hit_list.hsplist_array.iter().map(|hsp_list| hsp_list.oid))
        .chain(
            x_compact_hit_lists
                .iter()
                .flatten()
                .flat_map(|hit_list| hit_list.hsplist_array.iter().map(|hsp_list| hsp_list.oid)),
        )
        .collect();
    ordinals.sort_unstable();
    ordinals.dedup();
    ordinals
}

/// Stored residues between two sentinel bytes.
fn encoded_stored(residues: &[u8]) -> EncodedProtein {
    // (A sequence with a gap residue or a foreign byte never reaches a hit
    // list: `prepare_search_block` rejects its block.)
    let mut aa_seq = Vec::with_capacity(residues.len() + 2);
    aa_seq.push(NCBISTDAA_GAP);
    aa_seq.extend_from_slice(residues);
    aa_seq.push(NCBISTDAA_GAP);
    EncodedProtein {
        aa_seq,
        aa_len: residues.len(),
    }
}

/// Renames the subjects of `hit_lists` and `x_compact_hit_lists` to their ranks among `ordinals`.
fn rename_to_ranks(
    ordinals: &[u32],
    hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) {
    let rank = |ordinal: u32| -> u32 {
        ordinals
            .binary_search(&ordinal)
            .expect("every subject of the hit lists was listed") as u32
    };
    for hit_list in hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            hsp_list.oid = rank(hsp_list.oid);
            for hsp in &mut hsp_list.hsps {
                hsp.s_idx = rank(hsp.s_idx);
            }
        }
    }
    for hit_list in x_compact_hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            hsp_list.oid = rank(hsp_list.oid);
            for hsp in &mut hsp_list.hsps {
                hsp.hsp.s_idx = rank(hsp.hsp.s_idx);
            }
        }
    }
}

/// The stored titles and letters of the listed subjects as the reports read them, by rank
/// (`ListedSubjects::read_titles`).
pub(super) struct ListedTitles<'a> {
    listed: &'a ListedSubjects,
    /// The OIDs of the database (`SeqDb::oid_map`), by which LOSAT's refusals name the subjects.
    oids: crate::seqdb::OidMap<'a>,
    /// The titles (and letters) as the reports read them, by rank.
    pub(super) shown: ListedShown,
    /// Whether the stored title of a subject holds a byte outside ' ' ..= '~' (read as
    /// '#'), by rank (`check_reported_titles`).
    replaced_title_bytes: Vec<bool>,
    /// The sources from which NCBI makes the empty stored titles of the subjects without parsed
    /// Seq-ids, by rank (`shared::StoredTitle::source`; empty when none has one).
    pub(super) sources: Vec<crate::report::defline::SourceTaxname>,
    /// The subjects with parsed Seq-ids, by rank (empty when no reported subject has them;
    /// their entry of `shown` is the title of their first definition line).
    pub(super) parsed: Vec<Option<Arc<ParsedSubject>>>,
    /// The program that LOSAT's refusals name (`BlastpDatabase::program`).
    program: &'static str,
    /// The first read of the reads above, in NCBI's order, of a record that holds a UTF8String
    /// where a VisibleString is expected (`shared::HeaderReader::first_utf8`): NCBI posts its
    /// warning there, once per process (`shared::Utf8Read::post`).
    pub(super) utf8: Option<shared::Utf8Read>,
}

/// The reported subjects with parsed Seq-ids of `ListedSubjects::read_titles`, by rank.
#[derive(Default)]
pub(super) struct ParsedSubjects {
    pub(super) by_rank: Vec<Option<Arc<ParsedSubject>>>,
}

impl ParsedSubjects {
    /// Reads the header of the subject `ordinal` of rank `rank` (of `ranks`) of `database` with
    /// `headers` (its reader), where the hit list of query `query` reaches it
    /// (`shared::read_header`): its stored title (`shared::StoredTitle`), or for a sequence with
    /// parsed Seq-ids, kept here, the title of its first definition line.
    pub(super) fn read(
        &mut self,
        database: &BlastpDatabase,
        headers: &shared::HeaderReader<'_>,
        rank: usize,
        ordinal: u32,
        query: usize,
        ranks: usize,
    ) -> Result<shared::StoredTitle> {
        match shared::read_header(headers, ordinal, query, database.program)? {
            shared::ReportedHeader::Stored(stored) => Ok(stored),
            shared::ReportedHeader::Parsed(subject) => {
                let title = subject.first_title().to_vec();
                if self.by_rank.is_empty() {
                    self.by_rank = vec![None; ranks];
                }
                self.by_rank[rank] = Some(subject);
                Ok(shared::StoredTitle {
                    title,
                    ..Default::default()
                })
            }
        }
    }
}

/// The titles (and letters) of the listed subjects as the reports read them, by rank.
pub(super) enum ListedShown {
    /// The stored titles (empty for a subject in no final hit list) and the stored letters.
    Records(Vec<FastaRecord>),
    /// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the stored titles, where the subjects in no final hit
    /// list share one empty title; the letters are the stored residues
    /// (`ListedSubjects::encoded`, `x_listed_lean.rs`).
    Lean(Vec<Arc<[u8]>>),
}

impl ListedTitles<'_> {
    /// With `-parse_deflines`, `common::database::refuse_id_of_a_query` of each subject with
    /// parsed Seq-ids of the final hit lists `hit_lists`, with the scope of the query whose list
    /// reaches it (`scopes`): the reports of each query batch ask for the batch's own, also of a
    /// subject that a list of another batch reached first.
    pub(super) fn refuse_ids_of_queries(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        scopes: crate::blastinput::query_batch::QueryScopesFrom<'_>,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            for hsp_list in &hit_list.hsplist_array[..count] {
                let rank = hsp_list.oid as usize;
                if let Some(subject) = &self.parsed[rank] {
                    shared::refuse_id_of_a_query(
                        subject,
                        || self.oids.oid(self.listed.ordinals[rank]),
                        scopes.of_query(query),
                        self.program,
                    )?;
                }
            }
        }
        Ok(())
    }

    /// Outfmt 1 to 4 of queries read without `-parse_deflines` (line RG T3a, review 2 D2): LOSAT's
    /// refusal of a subject with parsed Seq-ids of the final hit lists `hit_lists` whose Seq-align
    /// id is the generated id (`lcl|Query_N`) of a query of the batch that reports it, as NCBI's
    /// scope finds an id (`scope_of`: that batch's scope, by query, `QueryScope::of_generated_ids`;
    /// the keys of `refuse_id_of_a_query`). NCBI's formatter adds the batch's queries to its scope
    /// before their reports and empties it after them (`query_batch::QueryScopeSpan::Batch`; no
    /// outfmt 5 here) and finds such a subject there: its reports show the query's title and
    /// letters in the subject's place, and the query-anchored merge of the query's own Bioseq can
    /// abort (`double free or corruption`, as `pairwise::blastp_section_aborts`). A weak match
    /// (`QueryScopeKeys::weak`, of text Seq-ids) never meets a local id. Then the subjects whose
    /// description lines the reports show (`described`) and whose display id is such an id
    /// (`refuse_display_ids`; review 3 of T3a, D1: `gi|777|lcl|Query_1`): NCBI shows the query's
    /// title in that line.
    pub(super) fn refuse_generated_query_ids<'s>(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        scope_of: impl Fn(usize) -> &'s crate::blastinput::query_batch::QueryScope,
        described: impl FnMut(usize) -> Result<Vec<u32>>,
    ) -> Result<()> {
        self.refuse_generated_ids(
            hit_lists,
            &scope_of,
            |_, _, _| true,
            shared::QUERY_ANCHORED_EFFECT,
        )?;
        self.refuse_display_ids(hit_lists, scope_of, true, described)
    }

    /// The Seq-align ids of the final hit lists `hit_lists` that are the generated id of a query
    /// of the batch that reports them (`scope_of`, `refuse_generated_query_ids`), refused with
    /// `effect` (what NCBI does with them in the formats asked for, `shared::XML_EFFECT` for
    /// rpsblast's BLAST XML) where the formats differ for them: `differs(query, list index,
    /// subject)` (the index among the query's reported lists, `EngineCounts`' key).
    pub(super) fn refuse_generated_ids<'s>(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        scope_of: &impl Fn(usize) -> &'s crate::blastinput::query_batch::QueryScope,
        differs: impl Fn(usize, usize, &crate::report::db_subject::ParsedSubject) -> bool,
        effect: &str,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            for (list_index, hsp_list) in hit_list.hsplist_array[..count].iter().enumerate() {
                let rank = hsp_list.oid as usize;
                let Some(subject) = &self.parsed[rank] else {
                    continue;
                };
                if !differs(query, list_index, subject) {
                    continue;
                }
                shared::refuse_generated_query_id(
                    subject,
                    self.oids.oid(self.listed.ordinals[rank]),
                    scope_of(query),
                    self.program,
                    effect,
                )?;
            }
        }
        Ok(())
    }

    /// With `-parse_deflines` (or a PSSM query shown with its ids) and a description table
    /// (outfmt 0 to 4): `common::database::refuse_display_id_of_a_query` of the subjects whose
    /// description lines the reports show (`described`), with the scope of the query whose list
    /// reaches each (`scopes`, as `refuse_ids_of_queries`, which refuses their Seq-align ids in
    /// every format; review 3 of line RG T3a, D1: `gi|779|lcl|abc` and a query `lcl|abc`).
    pub(super) fn refuse_display_ids_of_queries(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        scopes: crate::blastinput::query_batch::QueryScopesFrom<'_>,
        described: impl FnMut(usize) -> Result<Vec<u32>>,
    ) -> Result<()> {
        self.refuse_display_ids(hit_lists, |query| scopes.of_query(query), false, described)
    }

    /// `common::database::refuse_display_id_of_a_query` (`generated` as there) of each subject
    /// with parsed Seq-ids of the final hit lists `hit_lists` whose description line a report
    /// shows, with the scope `scope_of` of the query whose list reaches it: `described(query)`
    /// gives the subjects (by rank) of the query's description table
    /// (`pairwise::blastp_described_subjects`), asked only for a query whose list holds a subject
    /// whose display id the scope has (`QueryScopeKeys::display`).
    fn refuse_display_ids<'s>(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        scope_of: impl Fn(usize) -> &'s crate::blastinput::query_batch::QueryScope,
        generated: bool,
        mut described: impl FnMut(usize) -> Result<Vec<u32>>,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let scope = scope_of(query);
            let mut shown: Option<Vec<u32>> = None;
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            for hsp_list in &hit_list.hsplist_array[..count] {
                let rank = hsp_list.oid as usize;
                let Some(subject) = &self.parsed[rank] else {
                    continue;
                };
                if !scope.has_key(&subject.query_scope_keys().display) {
                    continue;
                }
                let shown = match &mut shown {
                    Some(shown) => shown,
                    None => shown.insert(described(query)?),
                };
                if shown.contains(&(rank as u32)) {
                    shared::refuse_display_id_of_a_query(
                        subject,
                        self.oids.oid(self.listed.ordinals[rank]),
                        scope,
                        generated,
                        self.program,
                    )?;
                }
            }
        }
        Ok(())
    }

    /// The stored title of the subject of rank `rank` (empty for a subject in no final hit
    /// list).
    fn title(&self, rank: usize) -> &[u8] {
        match &self.shown {
            ListedShown::Records(records) => &records[rank].title,
            ListedShown::Lean(titles) => &titles[rank],
        }
    }

    /// The records of `read_titles`.
    #[cfg(test)]
    fn records(&self) -> &[FastaRecord] {
        match &self.shown {
            ListedShown::Records(records) => records,
            ListedShown::Lean(_) => panic!("lean titles have no records"),
        }
    }

    /// LOSAT's rejection of a reported subject (one with rows in the final `hit_lists`)
    /// whose stored title holds a byte outside ' ' ..= '~'
    /// (`algorithm::common::database::check_reported_title`): NCBI writes an error message
    /// for the byte each time its formatter reads the definition line.
    pub(super) fn check_reported_titles(&self, hit_lists: &[Option<BlastpHitList>]) -> Result<()> {
        for rank in Self::reported(hit_lists) {
            shared::check_reported_title(
                self.replaced_title_bytes[rank],
                || self.oids.oid(self.listed.ordinals[rank]),
                self.program,
            )?;
        }
        Ok(())
    }

    /// The Blast-def-line taxid of each reported subject by rank (0 for none and for a
    /// subject in no final hit list), which the XML2 and JSON formats write in the `Hit`
    /// description, and what the taxonomy database says about it (for a subject with parsed
    /// Seq-ids, about the taxid of each definition line: `shared::read_reported_taxonomy`).
    /// They are read apart from the titles, for those formats only.
    ///
    /// NCBI writes the scientific name of a taxid when it finds `taxdb.bti` on the search path
    /// (`data4xml2format.cpp x_InitCommon`, `blastxml2_format.cpp:380-386`; validity is not
    /// checked): `taxonomy[rank]` is then `Named`, or `NotFound` when the taxonomy database does
    /// not know the taxid (the description is the id alone). Without `taxdb.bti`, or without a
    /// taxid, `taxonomy` is empty. A taxonomy database that LOSAT cannot read as NCBI does is a
    /// refusal made here, before any byte is written (a scientific name of any bytes is written
    /// as NCBI's streams write it).
    ///
    /// NCBI's messages about the taxonomy database go to standard error: those of opening and of
    /// the lookup that the volume makes for each taxid of a fetched subject (`CSeqDBVol::
    /// x_GetTaxonomy`, once per volume and taxid) are returned as `TaxdbNotes`, which
    /// `write_blastp_xml` writes at the first `Hit` that needs them, and those of the lookup of
    /// each `Hit` are written there too.
    ///
    /// Of the reported subjects, those that the reports load (`loads`, `hit_loads`: the reads of
    /// a subject whose `Hit` throws, and of one shown only after it, are never made, so the texts
    /// of opening the taxonomy database and of a volume's first lookup of a taxid go to the loads
    /// that make them; review 1 of line DBID pfoff-1, H1).
    pub(super) fn read_taxids(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        loads: &HitLoads,
    ) -> Result<(
        Vec<i32>,
        Vec<HitTaxonomy>,
        Vec<Vec<HitTaxonomy>>,
        shared::TaxdbNotes,
    )> {
        let mut reported = Self::reported(hit_lists);
        reported.retain(|&rank| loads.loaded(rank));
        shared::read_reported_taxonomy(
            &database.db,
            &database.taxdb,
            &self.listed.ordinals,
            &reported,
            &self.parsed,
            self.program,
        )
    }

    /// NCBI's messages about the taxonomy database when a query batch prefetches the subjects of
    /// its Seq-aligns, of ranks `ranks` (`shared::prefetch_messages`): BLASTP's are the subjects
    /// with rows in the batch's final hit lists (`reported_ranks`), the ones `read_taxids` reads;
    /// psiblast's those of its traceback hit lists, before its formatter's cut. Empty, and
    /// nothing read, without `taxdb.bti` on the search path. A taxonomy database that LOSAT
    /// cannot read as NCBI does, and a header whose taxid does not decode, are refused.
    pub(super) fn prefetch_messages(
        &self,
        database: &BlastpDatabase,
        ranks: impl IntoIterator<Item = usize>,
    ) -> Result<Vec<u8>> {
        shared::prefetch_messages(
            &database.db,
            &database.taxdb,
            ranks
                .into_iter()
                .map(|rank| (self.listed.ordinals[rank], self.parsed_subject(rank))),
            self.program,
        )
    }

    /// [`Self::prefetch_messages`] by subject (`shared::subject_read_texts`): the text of the
    /// first read of each subject of rank `ranks` that writes one, with its rank.
    pub(super) fn subject_read_texts(
        &self,
        database: &BlastpDatabase,
        ranks: impl IntoIterator<Item = usize>,
        rows_open: bool,
    ) -> Result<crate::report::query_warnings::QueryReads> {
        shared::subject_read_texts(
            &database.db,
            &database.taxdb,
            ranks.into_iter().map(|rank| {
                (
                    rank as u32,
                    self.listed.ordinals[rank],
                    self.parsed_subject(rank),
                )
            }),
            rows_open,
            self.program,
        )
    }

    /// `shared::scope_utf8_hit_reads` of the subjects of ranks `ranks` (an XML2 or JSON report's
    /// `Hit`s, in order, of a batch whose prefetch is off).
    pub(super) fn scope_utf8_hit_reads(
        &self,
        database: &BlastpDatabase,
        after: Option<crate::report::db_subject::RunThrow>,
        ranks: impl IntoIterator<Item = usize>,
    ) -> crate::report::query_warnings::QueryReads {
        shared::scope_utf8_hit_reads(
            &database.db,
            &database.taxdb,
            after,
            ranks.into_iter().map(|rank| {
                (
                    rank as u32,
                    self.listed.ordinals[rank],
                    self.parsed_subject(rank),
                )
            }),
            self.program,
        )
    }

    /// The searched queries of the first batch of `batches` (each batch's searched queries, in
    /// order) whose prefetch is on (`prefetch`) and reads a subject whose lookup fails on a file
    /// with NCBI's exception (`ParsedSubject::prefetch_stops`): NCBI's prefetch stops the run
    /// there (`shared::batch_prefetch`), before that batch's reports and those of the batches
    /// after it. `ranks` gives the ranks that a batch's prefetch reads, by its searched queries.
    pub(super) fn prefetch_stop<I: IntoIterator<Item = usize>>(
        &self,
        batches: &[std::ops::Range<usize>],
        prefetch: &shared::BatchPrefetch,
        ranks: impl Fn(std::ops::Range<usize>) -> I,
    ) -> Option<std::ops::Range<usize>> {
        if self.parsed.is_empty() {
            return None;
        }
        batches
            .iter()
            .find(|span| {
                !span.is_empty()
                    && !prefetch.reads_one_by_one(span.start)
                    && ranks((*span).clone()).into_iter().any(|rank| {
                        self.parsed_subject(rank)
                            .is_some_and(|subject| subject.prefetch_stops().is_some())
                    })
            })
            .cloned()
    }

    /// The subject of rank `rank` as a batch's prefetch reads it (a protein has no strand).
    fn prefetch_subject(&self, rank: usize) -> shared::PrefetchSubject<'_> {
        shared::PrefetchSubject::unflipped(self.listed.ordinals[rank], self.parsed_subject(rank))
    }

    /// NCBI's prefetch of a query batch's subjects of ranks `ranks` on standard error
    /// (`shared::batch_prefetch`): the taxonomy texts of `prefetch_messages`, then the warning
    /// where the prefetch fails at a subject whose lookup fails (a protein has no strand); the
    /// stop of the batch of `queries` (`shared::batch_prefetch`), after the batches of the call
    /// before it, for which NCBI's BLAST database data loader loaded the subjects of ranks
    /// `before` (read only where the batch stops; `shared::batch_loads`; the outfmt 5 reports
    /// before it, `shared::PrefetchQueries::after_xml_reports`).
    pub(super) fn batch_prefetch(
        &self,
        database: &BlastpDatabase,
        before: impl IntoIterator<Item = usize>,
        ranks: impl IntoIterator<Item = usize>,
        queries: &shared::PrefetchQueries<'_>,
    ) -> Result<crate::report::diag_log::Posts> {
        shared::batch_prefetch(
            &database.db,
            &database.taxdb,
            before.into_iter().map(|rank| self.prefetch_subject(rank)),
            ranks.into_iter().map(|rank| self.prefetch_subject(rank)),
            self.program,
            queries,
        )
    }

    /// The stop of NCBI's prefetch of a query batch's subjects of ranks `ranks` (the batch of
    /// `queries`), for the XML2 and JSON formats, which write the taxonomy texts at each `Hit`
    /// (`shared::held_batch_prefetch_stop`): after the lookups of the reports of the batches
    /// before, whose final hit lists are `before` (their subjects are those `read_taxids` reads),
    /// which BLASTP writes after the stop (the search of those batches again, `HeaderReadError`),
    /// and the loads of those batches of ranks `loaded` (`shared::batch_loads`, with psiblast's
    /// traceback beyond the cut, which no report reads), recorded where the batch stops. Of both,
    /// the subjects that those reports load (`loads`, `hit_loads` of `before`: not one whose
    /// `Hit` throws, nor one shown only after it; review 1 of line DBID pfoff-1, H1). `Ok`: the
    /// warning of a UTF8String header that a load of a batch that does not stop posts
    /// (`shared::batch_prefetch_stop`), for the batch's prefetch posts.
    pub(super) fn batch_prefetch_stop(
        &self,
        database: &BlastpDatabase,
        before: &[Option<BlastpHitList>],
        loaded: impl IntoIterator<Item = usize>,
        loads: &HitLoads,
        ranks: impl IntoIterator<Item = usize>,
        queries: &shared::PrefetchQueries<'_>,
    ) -> Result<Vec<u8>> {
        shared::held_batch_prefetch_stop(
            &database.db,
            &database.taxdb,
            Self::reported_ranks(before)
                .filter(|&rank| loads.loaded(rank))
                .map(|rank| self.prefetch_subject(rank)),
            loaded
                .into_iter()
                .filter(|&rank| loads.loaded(rank))
                .map(|rank| self.prefetch_subject(rank)),
            ranks.into_iter().map(|rank| self.prefetch_subject(rank)),
            self.program,
            queries,
        )
    }

    /// NCBI's warning of a UTF8String header that a load of the prefetch of the subjects of ranks
    /// `ranks` (a batch of the XML2 and JSON formats whose prefetch is on and does not stop) posts
    /// (`shared::scope_utf8_prefetch_post`); empty where none.
    pub(super) fn scope_utf8_prefetch(
        &self,
        database: &BlastpDatabase,
        ranks: impl IntoIterator<Item = usize>,
    ) -> Result<Vec<u8>> {
        let subjects: Vec<shared::PrefetchSubject<'_>> = ranks
            .into_iter()
            .map(|rank| self.prefetch_subject(rank))
            .collect();
        shared::scope_utf8_prefetch_post(&database.db, &subjects, self.program)
    }

    /// The run's record of the subjects of ranks `ranks` that NCBI's BLAST database data loader
    /// loaded for a call's batches (`shared::TaxdbRun::note_loaded`): psiblast's round of a query,
    /// a call of its own, for the stop of a later round's or query's prefetch, and the batches
    /// before a stop of the XML2 and JSON formats.
    pub(super) fn note_loaded<I: IntoIterator<Item = usize>>(
        &self,
        database: &BlastpDatabase,
        ranks: impl FnOnce() -> I,
    ) {
        if self.parsed.is_empty() {
            return;
        }
        database.taxdb.note_loaded(&database.db, || {
            ranks().into_iter().map(|rank| self.prefetch_subject(rank))
        });
    }

    /// The subjects of ranks `ranks` of psiblast's traceback hit lists as NCBI's prefetch reads
    /// them by their Seq-align ids (`shared::prefetch_reads`), for `batch_prefetch`:
    /// those it finds as themselves, and one whose lookup throws on its id
    /// (`ParsedSubject::prefetch_fails`: the batch's prefetch fails there); one that it finds
    /// nowhere is left out (the prefetch reads nothing for it). Such a failure is refused where
    /// the reports read the subjects they show (`reads`): the prefetch then leaves those after
    /// it unloaded, which NCBI's formatters read one by one (its scope can match an id of a
    /// sequence read before weakly, as with the prefetch off). One whose lookup fails on a file
    /// with NCBI's exception (`ParsedSubject::prefetch_stops`) is kept: the batch's prefetch stops
    /// there (`shared::batch_prefetch`). A lookup that fails on a file with LOSAT's error and
    /// another sequence that NCBI can read for one and whose definition lines do not hold its
    /// Seq-align id (NCBI stops) are refused. With `taxids`
    /// (`taxdb.bti` on the search path), one for which the prefetch reads another sequence that
    /// the lookup gives first is kept (`batch_prefetch` looks that sequence's taxids up), and one
    /// for which it can read either of two (`Prefetched::Another`) is refused where the taxonomy
    /// database's texts can differ between them (`shared::lookup_texts_differ`), and else kept as
    /// itself or left out. The subjects that the reports show were refused before for all of
    /// these (`check_parsed_reports`).
    pub(super) fn prefetched_ranks(
        &self,
        database: &BlastpDatabase,
        ranks: &[usize],
        taxids: bool,
        reads: bool,
    ) -> Result<Vec<usize>> {
        let mut kept = Vec::with_capacity(ranks.len());
        for &rank in ranks {
            let Some(subject) = self.parsed_subject(rank) else {
                kept.push(rank);
                continue;
            };
            if (subject.prefetch_fails() && !reads) || subject.prefetch_stops().is_some() {
                kept.push(rank);
                continue;
            }
            match shared::prefetch_reads(
                &database.db,
                self.listed.ordinals[rank],
                subject,
                taxids,
                self.program,
            )? {
                Ok(Some(_)) => kept.push(rank),
                Ok(None) => {}
                Err(reason) => {
                    return Err(shared::parsed_refusal(
                        self.oids.oid(self.listed.ordinals[rank]),
                        &format!(
                            "{reason}, and is in psiblast's traceback hit list, which NCBI BLAST+ 2.17.0 prefetches (also beyond -max_target_seqs, where no report shows it)"
                        ),
                        self.program,
                    ))
                }
            }
        }
        Ok(kept)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450
    /// What psiblast's PSSM reads of its subjects by their Seq-align ids, made from the
    /// Seq-aligns of `hit_lists` before the formatter's cut to `-max_target_seqs`
    /// (`PssmReads`): NCBI reads the letters of each Seq-align's range by its Seq-align id, also
    /// for subjects beyond the cut that no report shows, and only of the subjects with a
    /// Seq-align below `inclusion_ethresh` (`psi_round::included`). Following subjects whose
    /// Seq-align ids `Match` share a row of the PSSM's multiple alignment, each with the letters
    /// read by its own id (`psiblast::pssm::input_msa`). A subject's read is what NCBI's scope
    /// gets by its id (`shared::scope_bioseq`): its own letters, another sequence's, or none;
    /// LOSAT refuses where it does not read as NCBI (`ParsedSubject::pssm_read`, and another
    /// holder of the id that the scope can get instead) where a PSSM is made
    /// (`PssmReads::check`). A subject without parsed Seq-ids is read by its ordinal id
    /// (`gnl|BL_ORD_ID|n`): its own letters. Where every load of a sequence with a positive taxid
    /// throws (`shared::TaxdbRun::load_throws`), a read whose load looks such a taxid up
    /// (`shared::blob_taxids`) gets none, with parsed Seq-ids or without. `prefetch`: whether
    /// each query's batch prefetches (the round's), which reads them all by their own ids before
    /// the PSSM; one whose prefetch is off reads them one by one (`ParsedSubject::pssm_read`'s
    /// `one_by_one`). `taxids`: whether the batches look the subjects' taxids up (`taxdb.bti` on
    /// the search path); only then, and only for the queries whose batch reads one by one, does
    /// `pssm_load_texts` make the PSSM's loads of a query, from its hit list after the
    /// formatter's cut (`cut`: the length of each query's list after it, `prune_by_size`;
    /// `usize::MAX` for none), with what that cut drops and the blobs of the reads that do not
    /// load the subject's own, kept here (`QueryLoads`, `PssmReads::other_blobs`).
    pub(super) fn pssm_reads(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        prefetch: &shared::BatchPrefetch,
        inclusion_ethresh: f64,
        taxids: bool,
        cut: &dyn Fn(usize) -> usize,
    ) -> Result<PssmReads> {
        let mut reads = PssmReads::default();
        // The queries whose PSSM's reads can load a sequence that the batch did not read before
        // (`pssm_load_texts`).
        let keeps_loads = |query: usize| taxids && prefetch.reads_one_by_one(query);
        let any_loads = (0..hit_lists.len()).any(keeps_loads);
        // Every load of a sequence with a positive taxid throws (`shared::TaxdbRun::load_throws`):
        // the read of such a subject gets no letters, with parsed Seq-ids or without.
        let throws = database.taxdb.load_throws(&database.db);
        if self.parsed.is_empty() && !any_loads && !throws {
            return Ok(reads);
        }
        // The included subjects with parsed Seq-ids, and of each query whose loads are made the
        // included subjects that the formatter's cut drops from its hit list, in its order (each
        // subject once: its HSPs follow each other).
        let mut read: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        // The included subjects without them, where such a load throws.
        let mut by_ordinal: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        // The PSSM's reads one by one are kept for the run's record of loads only where a
        // prefetch of the run can stop (`shared::TaxdbRun::prefetch_can_stop`, asked at the first
        // such read; review 4 of line DBID loadord, M1).
        let mut records_loads: Option<bool> = None;
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let keep = keeps_loads(query);
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            let kept = cut(query).min(count);
            let mut beyond = Vec::new();
            for (index, hsp_list) in hit_list.hsplist_array[..count]
                .iter()
                .enumerate()
                .filter(|(_, hsp_list)| pssm_includes(hsp_list, inclusion_ethresh))
            {
                if self.parsed_subject(hsp_list.oid as usize).is_none() {
                    if throws {
                        by_ordinal.insert(hsp_list.oid as usize);
                    }
                } else {
                    read.insert(hsp_list.oid as usize);
                    if prefetch.reads_one_by_one(query)
                        && *records_loads
                            .get_or_insert_with(|| database.taxdb.prefetch_can_stop(&database.db))
                    {
                        if reads.one_by_one.len() <= query {
                            reads.one_by_one.resize_with(query + 1, Vec::new);
                        }
                        reads.one_by_one[query].push(hsp_list.oid as usize);
                    }
                }
                if keep && index >= kept {
                    beyond.push(hsp_list.oid);
                }
            }
            if keep {
                if reads.loads.len() <= query {
                    reads.loads.resize_with(query + 1, || None);
                }
                reads.loads[query] = Some(QueryLoads { kept, beyond });
            }
        }
        // The blob that each read of a subject with parsed Seq-ids loads, whose taxids NCBI looks
        // up there (`pssm_load_texts`), kept where it is not the subject's own.
        if !read.is_empty() {
            let batch_reads = Self::batch_reads(hit_lists, prefetch);
            for rank in read {
                let Some(subject) = self.parsed_subject(rank) else {
                    continue;
                };
                let ordinal = self.listed.ordinals[rank];
                // A lookup that fails on a file with NCBI's exception (a list of version 5
                // databases with and without the accession table): the read catches it, so the
                // PSSM reads no letters (`Failed to retrieve sequence`) and loads nothing.
                if subject.prefetch_stops().is_some() {
                    reads.read.insert(rank, None);
                    if any_loads {
                        reads.other_blobs.insert(rank, None);
                    }
                    continue;
                }
                // Also another holder of the id whose Bioseq the scope gets where it read it
                // before (`shared::ScopeBioseq::instead`).
                if let Err(reason) = subject.pssm_read(batch_reads.one_by_one(rank as u32)) {
                    reads
                        .refused
                        .get_or_insert((self.oids.oid(ordinal), reason, self.program));
                    if any_loads {
                        reads.other_blobs.insert(rank, None);
                    }
                    continue;
                }
                let scope = shared::scope_bioseq(&database.db, ordinal, subject, self.program)?;
                if let Some(blob) = scope.taxonomy().filter(|_| throws) {
                    if !shared::blob_taxids(
                        &database.db,
                        ordinal,
                        Some(subject),
                        blob,
                        self.program,
                    )?
                    .is_empty()
                    {
                        reads.read.insert(rank, None);
                        if any_loads {
                            reads.other_blobs.insert(rank, None);
                        }
                        continue;
                    }
                }
                if any_loads && scope.taxonomy() != Some(ordinal) {
                    reads.other_blobs.insert(rank, scope.taxonomy());
                }
                match scope.letters() {
                    Some(blob) if blob == ordinal => {}
                    Some(other) => {
                        let letters = database
                            .db
                            .read_sequence(other)
                            .map_err(|error| shared::database_error(error, self.program))?
                            .protein(0)
                            .to_vec();
                        reads.read.insert(rank, Some(letters));
                    }
                    None => {
                        reads.read.insert(rank, None);
                    }
                }
            }
        }
        // A subject without parsed Seq-ids is read by its ordinal id, which loads its own blob.
        for rank in by_ordinal {
            let ordinal = self.listed.ordinals[rank];
            if !shared::blob_taxids(&database.db, ordinal, None, ordinal, self.program)?.is_empty()
            {
                reads.read.insert(rank, None);
            }
        }
        Ok(reads)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1142-1182
    // (cited at `shared::SlotHolds::walk`)
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:513 (cited at
    // `shared::pssm_load_texts`)
    // Where the batch's prefetch is off, the PSSM's read of a Seq-align id loads its blob where
    // no fetch of the run loaded it before, and `x_GetTaxonomy` looks each taxid of its
    // definition lines up in the volume's 256-slot cache, which asks the taxonomy database at the
    // volume's first lookup of a taxid and where the slot holds another taxid; the texts of opening
    // come with the run's first lookup. Such a lookup writes its text between the query's report
    // and the PSSM's warnings of the later reads (review 1 of line DBID pssmread, H1; review 2, H1:
    // a report that read the subject before, or any lookup of the run before for the texts of
    // opening, writes them there; review 3, H1: the reads by an ordinal id of the subjects without
    // parsed Seq-ids load the same way; review 4, H1 and M1: the offset error again at each lookup
    // that misses its slot, in the order of the reads, also inside one load). It is the volume's
    // cache that keeps a taxid for the run's later queries and rounds, not NCBI's scope, which may
    // load a sequence read before again (`shared::pssm_load_texts`).
    /// NCBI's taxonomy texts of the PSSM's reads of query `query` (of the call) whose batch reads
    /// one by one, after that query's report read its subjects (`shared::prefetch_messages`):
    /// every read of the query's PSSM, in the order of the reads (the included subjects of its hit
    /// list `hit_list`, after the formatter's cut, then those that the cut dropped,
    /// `QueryLoads::beyond`), with the texts of the lookups of its load (`shared::blob_taxids`,
    /// `shared::pssm_load_texts`, with the run's `TaxdbRun`; empty for none, and for a read that
    /// loads no blob), for `psiblast::pssm::PssmReadMessages::add_load_texts`; empty for a query
    /// whose loads `pssm_reads` does not make. Only this query's reads are made, here: none are
    /// kept for the batch. The lookups are the run's from then on.
    pub(super) fn pssm_load_texts(
        &self,
        database: &BlastpDatabase,
        reads: &PssmReads,
        query: usize,
        hit_list: Option<&BlastpHitList>,
        inclusion_ethresh: f64,
    ) -> Result<Vec<(u32, Vec<u8>)>> {
        let Some(Some(loads)) = reads.loads.get(query) else {
            return Ok(Vec::new());
        };
        let shown: &[crate::algorithm::blastp::hsp::BlastpHspList] = match hit_list {
            Some(hit_list) => {
                &hit_list.hsplist_array[..hit_list.hsplist_count.min(hit_list.hsplist_array.len())]
            }
            None => &[],
        };
        if shown.len() != loads.kept {
            anyhow::bail!(
                "internal error: query {query}'s hit list holds {} subjects after the cut, {} when psiblast's PSSM reads were made",
                shown.len(),
                loads.kept
            );
        }
        let ranks = shown
            .iter()
            .filter(|hsp_list| pssm_includes(hsp_list, inclusion_ethresh))
            .map(|hsp_list| hsp_list.oid)
            .chain(loads.beyond.iter().copied());
        let mut texts = Vec::new();
        for rank in ranks {
            let ordinal = self.listed.ordinals[rank as usize];
            // A subject without parsed Seq-ids loads its own blob.
            let blob = match self.parsed_subject(rank as usize) {
                Some(_) => reads
                    .other_blobs
                    .get(&(rank as usize))
                    .copied()
                    .unwrap_or(Some(ordinal)),
                None => Some(ordinal),
            };
            let text = match blob {
                Some(blob) => {
                    let lookups = shared::blob_taxids(
                        &database.db,
                        ordinal,
                        self.parsed_subject(rank as usize),
                        blob,
                        self.program,
                    )?;
                    shared::pssm_load_texts(
                        &database.db,
                        &database.taxdb,
                        self.parsed_subject(rank as usize),
                        blob,
                        &lookups,
                        self.program,
                    )?
                }
                // The lookup read the header without a blob (`x_GetOid`'s check of the ids).
                None => shared::pssm_load_texts(
                    &database.db,
                    &database.taxdb,
                    self.parsed_subject(rank as usize),
                    ordinal,
                    &[],
                    self.program,
                )?,
            };
            texts.push((rank, text));
        }
        Ok(texts)
    }

    /// LOSAT's refusals of the subjects of ranks `ranks` that a batch's prefetch reads though its
    /// reports read nothing of them (the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`, a
    /// table of fields that read nothing of the subject), without `taxdb.bti` on the search path
    /// (`shared::check_prefetched_parsed`; with it, `prefetch_messages` refuses them too). A
    /// lookup that fails is the batch's warning (`batch_prefetch`).
    pub(super) fn check_prefetched(
        &self,
        database: &BlastpDatabase,
        ranks: impl IntoIterator<Item = usize>,
    ) -> Result<()> {
        for rank in ranks {
            if let Some(subject) = self.parsed_subject(rank) {
                shared::check_prefetched_parsed(
                    &database.db,
                    subject,
                    self.listed.ordinals[rank],
                    self.program,
                )?;
            }
        }
        Ok(())
    }

    /// The ranks of the subjects with rows in `hit_lists` (the final hit lists).
    fn reported(hit_lists: &[Option<BlastpHitList>]) -> std::collections::BTreeSet<usize> {
        // Inserted one by one: collecting would first gather the rank of every HSP of every query
        // in a vector and sort it (with a scratch of the same size, as `BTreeSet::from_iter` does),
        // 8 MB for 1,000 queries of 500 hits, to keep a few hundred ranks.
        let mut reported = std::collections::BTreeSet::new();
        for rank in Self::reported_ranks(hit_lists) {
            reported.insert(rank);
        }
        reported
    }

    /// What the batches read of the subjects with rows in `hit_lists`, by rank: `prefetch`
    /// tells by query whether the query's batch prefetches (`shared::BatchPrefetch::reads`).
    fn batch_reads(
        hit_lists: &[Option<BlastpHitList>],
        prefetch: &shared::BatchPrefetch,
    ) -> shared::BatchReads {
        prefetch.reads(hit_lists.iter().enumerate().flat_map(|(query, hit_list)| {
            Self::reported_ranks(std::slice::from_ref(hit_list))
                .map(move |rank| (query, rank as u32))
        }))
    }

    /// The rank of the subject of each row of `hit_lists`, in the order of the lists.
    pub(super) fn reported_ranks(
        hit_lists: &[Option<BlastpHitList>],
    ) -> impl Iterator<Item = usize> + '_ {
        Self::fetched_ranks(hit_lists, crate::report::tabular::SubjectFetch::Every)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-904
    /// The ranks of the subjects of `hit_lists`, one for each HSP whose row reads the subject
    /// (`rows`, by the HSP's identities, the Seq-align's `num_ident`; psiblast's Seq-align holds
    /// -1 for none, blast_seqalign.cpp:1492-1496), where a query's report fetches them
    /// (`PRE_FETCH_SEQS_LIMIT`).
    pub(super) fn fetched_ranks(
        hit_lists: &[Option<BlastpHitList>],
        rows: crate::report::tabular::SubjectFetch,
    ) -> impl Iterator<Item = usize> + '_ {
        hit_lists.iter().flatten().flat_map(move |hit_list| {
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            hit_list.hsplist_array[..count]
                .iter()
                .flat_map(move |hsp_list| {
                    hsp_list
                        .hsps
                        .iter()
                        .filter(move |hsp| rows.reads(hsp.num_ident))
                        .map(|hsp| hsp.s_idx as usize)
                })
        })
    }

    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every
    /// listed subject, by rank (`shared::tabular_subject_id`); a
    /// subject without rows gets an empty id, one whose lookup fails none
    /// (`common::tabular::SubjectId`: `TabularSubjectId::Unknown`, a parsed subject's
    /// `ParsedSubject::tabular_ids`; `report::tabular::check_failed_lookups`). A subject with rows
    /// whose id LOSAT cannot write as NCBI is an error before any format is written where a
    /// tabular report `reads` the subject (`report::tabular::reads_subject`), else its ids are
    /// empty: a table of fields that read nothing of the subject prints none of them.
    pub(super) fn tabular_subject_ids(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        reads: bool,
    ) -> Result<(
        Vec<crate::algorithm::common::tabular::SubjectId>,
        Vec<Option<[Arc<[u8]>; 2]>>,
    )> {
        let reported = Self::reported(hit_lists);
        let ranks = self.listed.ordinals.len();
        let mut ids: Vec<crate::algorithm::common::tabular::SubjectId> = Vec::with_capacity(ranks);
        // `sacc` and `saccver` of the subjects with parsed Seq-ids and of those with the general
        // id, by rank (empty without them; the others' are their `sseqid`).
        let mut accessions: Vec<Option<[Arc<[u8]>; 2]>> = Vec::new();
        // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): one empty id for every subject without rows.
        let empty: Arc<[u8]> = Arc::from(&b""[..]);
        for rank in 0..ranks {
            if !reported.contains(&rank) {
                ids.push(Some(match self.shown {
                    ListedShown::Lean(_) => Arc::clone(&empty),
                    ListedShown::Records(_) => Arc::from(&b""[..]),
                }));
                continue;
            }
            let ordinal = self.listed.ordinals[rank];
            if let Some(subject) = self.parsed_subject(rank) {
                // No ids where the lookup throws (`ParsedSubject::tabular_ids`).
                let found = match subject.tabular_ids(true) {
                    Ok(ids) => ids,
                    Err(_) if !reads => Some(Default::default()),
                    Err(problem) => {
                        return Err(shared::parsed_refusal(
                            self.oids.oid(ordinal),
                            &format!("has {problem}"),
                            self.program,
                        ))
                    }
                };
                ids.push(found.as_ref().map(|[id, _, _]| Arc::from(id.as_slice())));
                if accessions.is_empty() {
                    accessions = vec![None; ranks];
                }
                accessions[rank] =
                    found.map(|[_, accession, version]| [Arc::from(accession), Arc::from(version)]);
                continue;
            }
            let id = match shared::tabular_subject_id(
                self.title(rank),
                self.source(rank),
                self.oids.oid(ordinal),
                true,
                self.program,
            ) {
                Err(_) if !reads => {
                    ids.push(Some(Arc::from(&b""[..])));
                    continue;
                }
                id => id?,
            };
            ids.push(match id {
                shared::TabularSubjectId::Unknown(_) => None,
                ref id => Some(Arc::from(id.id())),
            });
            if let shared::TabularSubjectId::General(_) = id {
                if accessions.is_empty() {
                    accessions = vec![None; ranks];
                }
                let accession: Arc<[u8]> = Arc::from(id.accession());
                accessions[rank] = Some([Arc::clone(&accession), accession]);
            }
        }
        Ok((ids, accessions))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-851
    /// The taxonomy of the custom tabular fields (`staxid`, `staxids`;
    /// `report::tabular::TabularTaxonomy::of_deflines`) of the reported subjects (those with
    /// rows in `hit_lists`), by rank: from the definition lines of a subject with parsed Seq-ids,
    /// else from its stored ones (`SeqDb::defline_set`, which keeps whether a line has a
    /// `taxid`); the other ranks are left empty, and so are those of the subjects that no report
    /// loads (`loads`: one whose load throws, whose rows keep the taxonomy of the rows before, and
    /// one that no report reaches; NCBI reads none of their definition lines and looks none of
    /// their taxids up, `report::db_subject::HitLoads`). A header that does not decode is refused
    /// before the report, as BLASTN's (`common::tabular::tabular_taxonomy`).
    pub(super) fn tabular_taxonomy(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
        loads: &HitLoads,
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        use crate::report::tabular::TabularTaxonomy;
        let ranks = self.listed.ordinals.len();
        let mut taxonomy = vec![TabularTaxonomy::default(); ranks];
        for rank in Self::reported(hit_lists) {
            if rank >= ranks || !loads.loaded(rank) {
                continue;
            }
            taxonomy[rank] = match self.parsed_subject(rank) {
                Some(subject) => TabularTaxonomy::of_deflines(subject.deflines()),
                None => {
                    let set = database
                        .db
                        .defline_set(self.listed.ordinals[rank])
                        .map_err(|error| shared::database_error(error, self.program))?;
                    TabularTaxonomy::of_deflines(&set.deflines)
                }
            };
        }
        Ok(taxonomy)
    }

    /// The listed subjects as the library's records name them (`shared::library_subject`); a
    /// subject without rows, which no record names, has an empty title.
    pub(super) fn library_subjects(
        &self,
        database: &BlastpDatabase,
    ) -> Vec<crate::api::local_blast::ReportedSubject> {
        (0..self.listed.ordinals.len())
            .map(|rank| {
                shared::library_subject(
                    &database.db,
                    self.title(rank),
                    self.source(rank),
                    self.listed.ordinals[rank],
                    self.parsed_subject(rank),
                    self.listed.x_lengths[rank] as usize,
                    true,
                    self.program,
                )
            })
            .collect()
    }

    /// The subject of rank `rank` when it has parsed Seq-ids.
    pub(super) fn parsed_subject(&self, rank: usize) -> Option<&ParsedSubject> {
        self.parsed.get(rank).and_then(Option::as_deref)
    }

    /// The OID of the subject of rank `rank`.
    pub(super) fn oid(&self, rank: usize) -> crate::seqdb::Oid {
        self.oids.oid(self.listed.ordinals[rank])
    }

    /// The source from which NCBI makes the empty stored title of the subject of rank `rank`
    /// (`report::pairwise::stored_source`).
    pub(super) fn source(&self, rank: usize) -> &crate::report::defline::SourceTaxname {
        crate::report::pairwise::stored_source(&self.sources, rank as u32)
    }

    /// LOSAT's refusals of the subjects with parsed Seq-ids that the reports of the kinds asked
    /// for do not write as NCBI, before any format is written: a subject that NCBI's formatters
    /// cannot find by its Seq-align id where a format reads its Bioseq (`fetched`: outfmt 0, 5,
    /// 6, 7, 13 to 16; `ParsedSubject::fetched`), two of them that follow each
    /// other in a hit list with matching Seq-align ids where outfmt 0, 5 or 13 to 16 is written,
    /// but for copies of one record, which those reports show in one hit (`joined`,
    /// `shared::refuse_joined_hits`), and a word of a title made from the taxids that a table
    /// prints (`words`, `ParsedSubject::check_printed_word`). `taxonomy`: whether a format
    /// writes the subjects' taxids or their names (13 to 16, a table of taxonomy fields).
    /// `hit_lists` are in the order of the reports; `prefetch` tells by query whether its
    /// batch prefetches, else its reports read the subjects one by one
    /// (`ParsedSubject::fetched`'s `one_by_one`). `loads`: what the reports do at each subject's
    /// `Hit` (`hit_loads`, the one answer: a `Hit` that throws, a subject that no report loads).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn check_parsed_reports(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        prefetch: &shared::BatchPrefetch,
        fetched: bool,
        joined: bool,
        taxonomy: bool,
        loads: &HitLoads,
        words: crate::report::tabular::PrintedIds,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        if fetched {
            let reads = Self::batch_reads(hit_lists, prefetch);
            for rank in Self::reported(hit_lists) {
                if let Some(subject) = self.parsed_subject(rank) {
                    subject
                        .fetched(taxonomy, reads.one_by_one(rank as u32), loads.load(rank))
                        .map_err(|reason| {
                            shared::parsed_refusal(
                                self.oids.oid(self.listed.ordinals[rank]),
                                &reason,
                                self.program,
                            )
                        })?;
                    subject.check_printed_word(true, words).map_err(|problem| {
                        shared::parsed_refusal(
                            self.oids.oid(self.listed.ordinals[rank]),
                            &format!("has {problem}"),
                            self.program,
                        )
                    })?;
                }
            }
        }
        if joined {
            for (query, hit_list) in hit_lists.iter().enumerate() {
                let Some(hit_list) = hit_list else {
                    continue;
                };
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                shared::refuse_joined_hits(
                    query,
                    hit_list.hsplist_array[..count].iter().map(|hsp_list| {
                        let rank = hsp_list.oid as usize;
                        self.parsed_subject(rank)
                            .map(|subject| (self.listed.ordinals[rank], subject))
                    }),
                    self.oids,
                    self.program,
                )?;
            }
        }
        Ok(())
    }

    /// What the XML2 and JSON reports of the final hit lists `hit_lists` (in the order of the
    /// reports) do at each subject's `Hit`, by rank (`HitLoads::of`, the one answer of the
    /// call): `prefetch` tells by query whether its batch reads one by one
    /// (`BatchPrefetch::reads_one_by_one`); the `Hit`s are those that `blast_xml::write_blastp_xml`
    /// makes, of the lists with HSPs. `run`: what the run's formats do after a subject whose
    /// lookup throws (the call's answer, `FormatThrow::of_outfmts` and `FormatThrow::of_search`;
    /// every subject is loaded where it writes nothing, `FormatThrow::run`), at the rows that load
    /// it (by the HSPs' identities, the Seq-align's `num_ident`, as `fetched_ranks`).
    /// A pairwise run (outfmt 0 to 4) reads of query `query` what its writer reads,
    /// `pairwise_reads(query)` (by rank), asked where its batch reads one by one.
    pub(super) fn hit_loads(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        prefetch: &shared::BatchPrefetch,
        run: crate::report::db_subject::FormatThrow,
        pairwise_reads: impl Fn(usize) -> crate::report::db_subject::PairwiseReads,
    ) -> HitLoads {
        let lists = |query: usize| {
            hit_lists[query].as_ref().map_or(&[][..], |hit_list| {
                &hit_list.hsplist_array[..hit_list.hsplist_count.min(hit_list.hsplist_array.len())]
            })
        };
        // Without a subject of parsed Seq-ids no lookup throws.
        let run = if self.parsed.is_empty() {
            crate::report::db_subject::FormatThrow::Unread
        } else {
            run
        };
        let pairwise = run.run().is_some_and(|run| run.after.is_pairwise());
        HitLoads::of(
            run,
            (0..hit_lists.len()).map(|query| {
                // A batch that prefetches loads every subject, in any order.
                let one_by_one = prefetch.reads_one_by_one(query);
                let subjects = if pairwise && one_by_one {
                    crate::report::db_subject::ReportKeys::Pairwise(pairwise_reads(query))
                } else {
                    crate::report::db_subject::ReportKeys::Entries(
                        lists(query)
                            .iter()
                            .map(|hsp_list| hsp_list.hsps.first().map(|hsp| hsp.s_idx as usize))
                            .collect(),
                    )
                };
                (one_by_one, subjects)
            }),
            |rank| self.parsed_subject(rank),
            |query, list| {
                lists(query)[list]
                    .hsps
                    .iter()
                    .map(|hsp| hsp.num_ident)
                    .max()
                    .unwrap_or(0)
            },
        )
    }

    /// The number of HSP lists of each of psiblast's traceback hit lists `hit_lists` (in the
    /// order of the reports) that its cut to `max_target_seqs` subjects keeps for every format,
    /// as NCBI's formatter counts them (`shared::joined_cut`): a copy of one record that follows
    /// another counts with it; two other parsed subjects that follow each other with matching
    /// Seq-align ids where the count reaches them are LOSAT's refusal.
    pub(super) fn joined_cut(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        max_target_seqs: usize,
    ) -> Result<Vec<usize>> {
        hit_lists
            .iter()
            .enumerate()
            .map(|(query, hit_list)| {
                let Some(hit_list) = hit_list else {
                    return Ok(max_target_seqs);
                };
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                if self.parsed.is_empty() || count <= max_target_seqs {
                    return Ok(max_target_seqs);
                }
                shared::joined_cut(
                    query,
                    max_target_seqs,
                    hit_list.hsplist_array[..count].iter().map(|hsp_list| {
                        let rank = hsp_list.oid as usize;
                        self.parsed_subject(rank)
                            .map(|subject| (self.listed.ordinals[rank], subject))
                    }),
                    self.oids,
                    self.program,
                )
            })
            .collect()
    }

    /// The `stitle` text of every listed subject, by rank: the stored title (an empty
    /// one is written `N/A` by the shared tabular writer, `report::tabular`).
    pub(super) fn stitles(&self) -> Vec<Arc<[u8]>> {
        match &self.shown {
            ListedShown::Records(records) => records
                .iter()
                .map(|record| Arc::from(record.title()))
                .collect(),
            ListedShown::Lean(titles) => titles.clone(),
        }
    }

    /// LOSAT's rejection of the outfmt 0 titles that it does not write as NCBI, of the
    /// subjects that the report shows (`shown`, by rank): a title of punctuation that
    /// NCBI's `x_CleanAndCompress` reads past the end of when it makes the title of the
    /// description table (NCBI BLAST+ 2.17.0 crashes there), as for the records of a
    /// subject set (`check_shown_subject_titles`), and a title that NCBI makes from a source
    /// that LOSAT refuses (`shared::StoredTitle::check_shown_source`). A subject whose load
    /// throws (`loads`) shows only its Seq-align id (`shared::check_pairwise_parsed`).
    pub(super) fn check_shown_titles(
        &self,
        database: &BlastpDatabase,
        shown: &std::collections::BTreeSet<usize>,
        loads: &HitLoads,
    ) -> Result<()> {
        for &rank in shown {
            let ordinal = self.listed.ordinals[rank];
            match self.parsed_subject(rank) {
                Some(subject) => shared::check_pairwise_parsed(
                    subject,
                    loads.load(rank),
                    true,
                    database.db.taxdb_found(),
                    || self.oids.oid(ordinal),
                    self.program,
                )?,
                None => {
                    shared::check_shown_protein_title(
                        self.title(rank),
                        || self.oids.oid(ordinal),
                        self.program,
                    )?;
                    shared::source_name(
                        self.source(rank),
                        || self.oids.oid(ordinal),
                        self.program,
                    )?;
                }
            }
        }
        Ok(())
    }
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search writes a block's rows before the next block is
// searched (`x_query_blocks.rs`), so with `-db` it may run only where the search of all queries
// would not fail after some rows. Past the preliminary stage (whose reader errors come in the
// first block, which reads the whole database as every batch does), the failures of a database
// search that depend on its hits are those of the definition lines of the reported subjects:
// NCBI's reader error (`HeaderReadError`, the batch of the query that reaches it fails after the
// reports of the batches before) and LOSAT's refusals (`ListedSubjects::read_titles` and
// `ListedSubjects::x_read_titles_lean`,
// `ListedTitles::check_reported_titles`, `ListedTitles::tabular_subject_ids`, and for outfmt 0
// `ListedTitles::check_shown_titles`). When every definition line of the database reads and none
// of them is refused, none of these failures can happen in any batch. With `taxdb.bti` on the
// search path, the prefetch of each batch (`ListedTitles::prefetch_messages`) also reads the
// reported subjects' taxids and looks them up in the taxonomy database, which can refuse as
// XML2's lookups can (DBID item T): `taxids` then also asks that every sequence's taxid decodes
// and that no lookup can be refused (`shared::taxdb_may_refuse` without the names, which the
// prefetch does not ask for), as `x_titles_taxids_allow_blocks` asks for outfmt 15 and 16.
// The prefetch's refusal of NCBI's taxonomy cache (`shared::TaxdbRun::cache_collision`,
// `shared::prefetch_slot_plan`) can come at any batch that prefetches: a run that may make it
// takes the search of all queries (`titles_allowed` in `blast_engine.rs`).

// No NCBI counterpart: the order of LOSAT's own checks; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): whether `check` holds for every range of
/// `db.header_ranges(HEADER_RANGE_SEQS)`, read in `pool` as the title checks read them, with each
/// range made from its index (`SeqDb::x_header_range`): no list of the ranges is allocated, whose
/// free (1.5 MiB for UniRef50) would raise glibc's mmap threshold before the search.
fn x_all_header_ranges(
    db: &crate::seqdb::SeqDb,
    pool: &crate::utils::threading::SearchPool<'_>,
    check: impl Fn(&BlockRange) -> bool + Sync,
) -> bool {
    let seqs = crate::seqdb::HEADER_RANGE_SEQS;
    let count = db.x_header_range_count(seqs);
    let check_index = |index: usize| check(&db.x_header_range(seqs, index));
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() {
        return pool.install(|| (0..count).into_par_iter().all(check_index));
    }
    let _ = pool;
    (0..count).all(check_index)
}

// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): whether every sequence of the database has one definition line
/// that reads, holds no byte outside ' ' ..= '~' and gives a tabular subject id as NCBI writes
/// it, or with `pairwise` (an outfmt 0 report, which writes no tabular id) a
/// title that the outfmt 0 report writes as NCBI, and with `taxids` (the batches' prefetch reads
/// the taxonomy database) a taxid that decodes, where no lookup of the taxonomy database can be
/// refused. Read in `pool` over the ordinal ids; a reader error is `false` (the search of all
/// queries then reports it where NCBI does). The tabular id is the reports' own question
/// (`shared::tabular_subject_id` of the title and the record's source,
/// `shared::HeaderReader::stored_source` of the taxid read with it, as `shared::read_header`
/// gives them to `ListedTitles::tabular_subject_ids`), so a title is ruled out here wherever the
/// search of all queries can refuse it (a failed lookup, `TabularSubjectId::Unknown`, also where
/// the field list does not stop at it).
/// EXPERIMENT (LOSAT_X_DBIDXREAD): the definition lines are read in ranges of
/// `HEADER_RANGE_SEQS` sequences (`SeqDb::visit_deflines`, the same results as
/// `SeqDb::deflines`, the titles and the taxids in one pass) instead of one sequence at a time.
pub(super) fn x_titles_allow_blocks(
    database: &BlastpDatabase,
    pairwise: bool,
    taxids: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    let allowed_deflines = |headers: &shared::HeaderReader<'_>,
                            ordinal: u32,
                            deflines: std::result::Result<
        Vec<crate::seqdb::OrdinalDefline>,
        crate::seqdb::SeqDbError,
    >|
     -> bool {
        match deflines {
            // A record with a UTF8String where a VisibleString is expected reads with NCBI's
            // warning, which the block whose call reads it first posts where the search of all
            // queries posts it (`shared::Utf8Read`, in the first query of the batch that reads it;
            // the blocks write their reports one after the other, and the warning once per run).
            Ok(deflines) => {
                deflines.len() == 1
                    && deflines[0].replaced_title_bytes == 0
                    && if pairwise {
                        // The check of the search of all queries (`ListedTitles::check_shown_titles`),
                        // which the name of an empty title's source does not change
                        // (`check_shown_protein_title`); a lookup that LOSAT refuses is ruled out
                        // with `taxids` (`taxdb_may_refuse`), and without `taxdb.bti` nothing is
                        // looked up. A sequence with parsed Seq-ids is not `deflines` (`Err`): the
                        // search of all queries checks its title from its own source
                        // (`check_shown_parsed_subject`).
                        shared::check_shown_protein_title(
                            &deflines[0].title,
                            || database.db.oid(ordinal),
                            database.program,
                        )
                        .is_ok()
                    } else {
                        // The tabular id of the title with its source (an empty title made from
                        // it, the end of a title ending in `]` rewritten from it). A lookup that
                        // LOSAT refuses is ruled out with `taxids` too. A failed lookup
                        // (`TabularSubjectId::Unknown`) is refused by some field lists
                        // (`report::tabular::failed_lookup_stop`): not allowed.
                        let defline = &deflines[0];
                        matches!(
                            shared::tabular_subject_id(
                                &defline.title,
                                &headers.stored_source(&defline.title, defline.taxid, ordinal),
                                db.oid(ordinal),
                                true,
                                database.program,
                            ),
                            Ok(shared::TabularSubjectId::Word(_)
                                | shared::TabularSubjectId::General(_))
                        )
                    }
            }
            Err(_) => false,
        }
    };
    if crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        // A record whose taxid does not decode is an `Err` (`visit_deflines`). Each range has its
        // reader of the sources' names (`HeaderReader`, which opens the taxonomy database at its
        // first lookup).
        let allowed_range = |range: &BlockRange| -> bool {
            let headers = shared::HeaderReader::new(db);
            db.visit_deflines(*range, |ordinal, deflines| {
                allowed_deflines(&headers, ordinal, deflines)
            })
        };
        if crate::seqdb::x_bp_db_keep() {
            return x_all_header_ranges(db, pool, allowed_range)
                && !(taxids && shared::taxdb_may_refuse(db));
        }
        let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let titles = if pool.enabled() {
            pool.install(|| ranges.par_iter().all(allowed_range))
        } else {
            ranges.iter().all(allowed_range)
        };
        #[cfg(not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )))]
        let titles = {
            let _ = pool;
            ranges.iter().all(allowed_range)
        };
        return titles && !(taxids && shared::taxdb_may_refuse(db));
    }
    let allowed = |headers: &mut shared::HeaderReader<'_>, ordinal: u32| -> bool {
        allowed_deflines(headers, ordinal, db.deflines(ordinal))
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let titles = if pool.enabled() {
        pool.install(|| {
            (0..db.num_seqs())
                .into_par_iter()
                .with_min_len(4096)
                .map_init(|| shared::HeaderReader::new(db), allowed)
                .all(|allowed| allowed)
        })
    } else {
        let mut headers = shared::HeaderReader::new(db);
        (0..db.num_seqs()).all(|ordinal| allowed(&mut headers, ordinal))
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let titles = {
        let mut headers = shared::HeaderReader::new(db);
        (0..db.num_seqs()).all(|ordinal| allowed(&mut headers, ordinal))
    };
    // One sequence at a time (`LOSAT_X_DBIDXREAD=0`): the taxids after the titles.
    titles && (!taxids || x_taxids_allow_blocks(database, pool))
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search of outfmt 8, 9, 12 and 5.
// These reports show no tabular id and no outfmt 0 title, so of the failures that
// `x_titles_allow_blocks` rules out only those of the definition lines' reading remain: NCBI's
// reader error (`HeaderReadError`), more than one definition line, and a stored title with a byte
// outside ' ' ..= '~' (`ListedSubjects::read_titles`, `ListedTitles::check_reported_titles`;
// outfmt 5's stream writes such a title as NCBI's does, so there the check only keeps such a
// database off the block search). A new function, so that
// `x_titles_allow_blocks` and the per-definition-line code of the outfmt 0, 6 and 7 path stay as
// they are.
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): whether every sequence of the database has
/// one definition line that reads and holds no byte outside ' ' ..= '~'. Read in `pool` as
/// `x_titles_allow_blocks` reads it; a reader error is `false`.
pub(super) fn x_titles_read_allow_blocks(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    let readable = |deflines: std::result::Result<
        Vec<crate::seqdb::OrdinalDefline>,
        crate::seqdb::SeqDbError,
    >|
     -> bool {
        match deflines {
            // A UTF8String record as in `x_titles_allow_blocks`.
            Ok(deflines) => deflines.len() == 1 && deflines[0].replaced_title_bytes == 0,
            Err(_) => false,
        }
    };
    if crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        let readable_range = |range: &BlockRange| -> bool {
            db.visit_deflines(*range, |_ordinal, deflines| readable(deflines))
        };
        if crate::seqdb::x_bp_db_keep() {
            return x_all_header_ranges(db, pool, readable_range);
        }
        let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            return pool.install(|| ranges.par_iter().all(readable_range));
        }
        let _ = pool;
        return ranges.iter().all(readable_range);
    }
    let allowed = |ordinal: u32| -> bool { readable(db.deflines(ordinal)) };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() {
        return pool.install(|| {
            (0..db.num_seqs())
                .into_par_iter()
                .with_min_len(4096)
                .all(allowed)
        });
    }
    let _ = pool;
    (0..db.num_seqs()).all(allowed)
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search of outfmt 15 and 16. Besides what outfmt 5
// reads (`x_titles_read_allow_blocks`) they read, for each subject with hits, the Blast-def-line
// taxid (`SeqDb::defline_taxids`, a decode apart from the titles') and, when `taxdb.bti` is on the
// search path, the taxonomy database (`shared::read_reported_taxonomy`): its opening and the lookup
// of the taxid (a scientific name of any bytes is written as NCBI's stream writes it). The search
// of all queries makes those refusals before its first byte; a block would make them after the
// earlier blocks' bytes. So the block search is taken only where none can fail for any subject
// of the database: every sequence's taxid reads, and `shared::taxdb_may_refuse` (OF5b, the one
// implementation of "can a lookup of this database refuse?", DECISIONS 34 and 38) says no (the
// refusal of NCBI's taxonomy cache, `shared::TaxdbRun::cache_may_refuse`, is asked apart, by
// `titles_allowed` in `blast_engine.rs`). That
// check answers from the smaller side: one pass over the taxonomy database when it is the smaller
// one (NCBI's real taxdb has no refusable record, so no taxid lookup and no pass over the
// database's taxids), else the database's taxids, sorted and deduplicated, looked up. It has its
// own taxonomy handle, so the run's open messages stay for the report. This is conservative,
// since the report reads the taxids of the subjects with hits only.
// The titles and the taxid decodes are read in one pass over the headers (`SeqDb::visit_deflines`,
// one read per `HEADER_RANGE_SEQS` sequences, as `x_titles_read_allow_blocks` reads them); no set of
// taxids is built here.
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): `x_titles_read_allow_blocks` and `x_taxids_allow_blocks`
/// in one pass over the headers: whether every sequence of the database has one definition line
/// that reads and holds no byte outside ' ' ..= '~', its taxid reads, and no lookup of the
/// taxonomy database can be refused (`shared::taxdb_may_refuse`). Read in `pool`; a failure is
/// `false`. It writes nothing and changes no state of
/// the run. With `LOSAT_X_DBIDXREAD=0` the two checks run one after the other, one sequence at a
/// time.
pub(super) fn x_titles_taxids_allow_blocks(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    if crate::seqdb::x_dbidxread() == crate::seqdb::XDbIdxRead::Off {
        return x_titles_read_allow_blocks(database, pool) && x_taxids_allow_blocks(database, pool);
    }
    let readable_range = |range: &BlockRange| -> bool {
        db.visit_deflines(*range, |_ordinal, deflines| match deflines {
            // The taxid `read_reported_taxonomy` takes of a sequence: its first definition
            // line's.
            // A UTF8String record as in `x_titles_allow_blocks`.
            Ok(deflines) => deflines.len() == 1 && deflines[0].replaced_title_bytes == 0,
            Err(_) => false,
        })
    };
    if crate::seqdb::x_bp_db_keep() {
        return x_all_header_ranges(db, pool, readable_range) && !shared::taxdb_may_refuse(db);
    }
    let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let readable = if pool.enabled() {
        pool.install(|| ranges.par_iter().all(readable_range))
    } else {
        ranges.iter().all(readable_range)
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let readable = {
        let _ = pool;
        ranges.iter().all(readable_range)
    };
    readable && !shared::taxdb_may_refuse(db)
}

/// EXPERIMENT (LOSAT_X_BPQBLOCK): whether every sequence's taxid reads and no lookup of the
/// taxonomy database can be refused (`shared::taxdb_may_refuse`, which is `false` where
/// `taxdb.bti` is not found). One sequence at a time (`LOSAT_X_DBIDXREAD=0`); see
/// `x_titles_taxids_allow_blocks`.
fn x_taxids_allow_blocks(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    let reads = |ordinal: u32| db.defline_taxids(ordinal).is_ok();
    let count = db.num_seqs();
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let all_read = if pool.enabled() {
        pool.install(|| (0..count).into_par_iter().with_min_len(4096).all(reads))
    } else {
        (0..count).all(reads)
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let all_read = {
        let _ = pool;
        (0..count).all(reads)
    };
    all_read && !shared::taxdb_may_refuse(db)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450 (cited at
// `algorithm::psiblast::pssm::input_msa`)
/// What psiblast's PSSM reads of the subjects that it includes, by rank, where that is not the
/// subject's own letters (`ListedTitles::pssm_reads`): for one with parsed Seq-ids another
/// sequence's letters in NCBIstdaa, or `None` for no sequence (NCBI's warning `Failed to retrieve
/// sequence`, and the row becomes the query's), which a subject without them also reads where its
/// load throws (`shared::TaxdbRun::load_throws`). `refused`: LOSAT's refusal of the first subject
/// whose read it does not make as NCBI, which stops a run where a PSSM is made (`check`). `loads`: by query
/// (of the call), what `ListedTitles::pssm_load_texts` needs to make the reads of a query whose
/// batch reads one by one with `taxdb.bti` on the search path (`QueryLoads`), `None` elsewhere;
/// `other_blobs`: then, the subjects with parsed Seq-ids whose read loads the blob of another
/// sequence (`shared::ScopeBioseq::taxonomy`), or none (`None`, also for a refused read), by rank.
/// The reads themselves are made one query at a time, so the batch keeps no list of them.
#[derive(Default)]
pub(super) struct PssmReads {
    read: std::collections::BTreeMap<usize, Option<Vec<u8>>>,
    /// The OID, the reason and the program of the refusal.
    refused: Option<(crate::seqdb::Oid, String, &'static str)>,
    loads: Vec<Option<QueryLoads>>,
    other_blobs: std::collections::BTreeMap<usize, Option<u32>>,
    /// By query (of the call), the ranks of the included subjects with parsed Seq-ids that the
    /// PSSM of a query whose batch reads one by one reads by their Seq-align ids, before the cut
    /// (`one_by_one`).
    one_by_one: Vec<Vec<usize>>,
}

/// What the PSSM's reads of a query need beyond its hit list after the formatter's cut
/// (`PssmReads::loads`): `kept`, that list's length, which `pssm_load_texts` checks; `beyond`, the
/// ranks of the included subjects that the cut drops, in their order (NCBI's PSSM reads them too;
/// empty unless the traceback list is longer than `-max_target_seqs`).
struct QueryLoads {
    kept: usize,
    beyond: Vec<u32>,
}

/// Whether psiblast's PSSM includes the subject of `hsp_list`, so reads its letters: one of its
/// HSPs has an e-value below `inclusion_ethresh` (`psi_round::included`).
fn pssm_includes(
    hsp_list: &crate::algorithm::blastp::hsp::BlastpHspList,
    inclusion_ethresh: f64,
) -> bool {
    hsp_list
        .hsps
        .iter()
        .any(|hsp| super::psi_round::included(hsp.e_value, inclusion_ethresh))
}

impl PssmReads {
    /// LOSAT's refusal of the reads, where a PSSM is made (`-out_pssm` and `-out_ascii_pssm`, and
    /// a round that makes the next round's PSSM).
    /// The included subjects with parsed Seq-ids, by rank, that the PSSM of the call's query
    /// `query` reads by their Seq-align ids where its batch reads one by one (prefetch off:
    /// `shared::BatchPrefetch::reads_one_by_one`), before the cut: loads of NCBI's BLAST
    /// database data loader that no prefetch made (`shared::TaxdbRun::note_loaded`; where the
    /// prefetch is on, it loaded them first). Empty for the other queries, where no PSSM reads
    /// are made, and where no prefetch of the run can stop (`shared::TaxdbRun::prefetch_can_stop`:
    /// nothing reads them then).
    pub(super) fn one_by_one(&self, query: usize) -> &[usize] {
        self.one_by_one.get(query).map_or(&[], Vec::as_slice)
    }

    pub(super) fn check(&self) -> Result<()> {
        match &self.refused {
            Some((oid, reason, program)) => Err(shared::parsed_refusal(
                *oid,
                &format!("{reason}, and is in psiblast's PSSM, whose letters NCBI BLAST+ 2.17.0 reads by that id (also beyond -max_target_seqs)"),
                program,
            )),
            None => Ok(()),
        }
    }

    /// The letters that the PSSM reads by the Seq-align id of the subject of rank `s_idx`:
    /// `own()` (its own letters), another sequence's, or `None` for no sequence.
    pub(super) fn letters(
        &self,
        s_idx: u32,
        own: impl FnOnce() -> Result<Vec<u8>>,
    ) -> Result<Option<Vec<u8>>> {
        match self.read.get(&(s_idx as usize)) {
            Some(read) => Ok(read.clone()),
            None => own().map(Some),
        }
    }
}

#[cfg(test)]
#[path = "database_tests.rs"]
mod tests;
