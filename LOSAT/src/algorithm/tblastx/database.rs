//! The database path of TBLASTX (`-db NAME`): the search of an NCBI BLAST nucleotide
//! database, as NCBI BLAST+ 2.17.0 `tblastx -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, stage DB-5).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the subjects are the database's sequences, checked once as a copy of them would be read
//!   (`TblastxDatabase::x_stream`) and read where each query batch searches them, in blocks of
//!   the sequence file shared by the threads, a long sequence alone
//!   (`subjects::Subjects::Stream`, `common::x_db_stream`): the preliminary stage translates the
//!   stored 2-bit data (makeblastdb's bases for the ambiguous ones, NCBI's preliminary input) and
//!   the ungapped re-evaluation reads the stored ambiguity codes, as the packed store of the
//!   sequence (`PackedNucleotide::from_stored`) gives them; the search's subject index is the
//!   ordinal id. A read that fails during the search fails its query batch with the error that
//!   the check gives for it;
//! - the set-level numbers come from the index files (`TblastxDatabase::numbers`): the
//!   number of sequences and the total length (an `Int8`); TBLASTX has no Gumbel block, so
//!   the shortest length a database answers does not enter;
//! - the definition lines of the reported sequences are read after each query batch's search
//!   (`read_titles`), and so are the sequences when a report reads their letters
//!   (`read_reported`); the reports show the database (`TblastxDatabase::name`,
//!   `pairwise_database`) and the stored titles of the reported sequences only
//!   (`reported_subjects`), or the definition lines of those with parsed Seq-ids
//!   (`report::db_subject`, line DBID ID2b).
//!
//! The search, the hit order and the report writers are those of the `-subject` path; the
//! opening of the database, the reader's errors and the rules for the stored titles are
//! shared with BLASTP, BLASTN, TBLASTN and BLASTX (`algorithm::common::database`).

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::algorithm::common::database::{self as shared, ParsedShown, ReportedHeader};
use crate::algorithm::common::x_db_stream::XDbStream;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::db_subject::{HitLoads, ParsedSubject};
use crate::report::pairwise::PairwiseDatabase;
use crate::report::tabular::SubjectFetch;
use crate::seqdb::{SeqDb, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

use super::blast_engine::TblastxHsp;
use super::subjects::SubjectSetNumbers;

/// The program that LOSAT's refusals name.
const PROGRAM: &str = "TBLASTX";

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`TblastxArgs::database_name`: the `-db` value without white space at its ends).
pub(crate) struct TblastxDatabase {
    db: SeqDb,
    name: String,
    /// Whether the first header of some volume holds parsed Seq-ids (`shared::check_seqids`).
    parsed: bool,
    /// The taxonomy database of the run and the taxids its volumes have looked up (XML2 and
    /// JSON, and the batches' prefetch of outfmt 0, 5, 6 and 7), kept across the query batches
    /// of the run.
    taxdb: shared::TaxdbRun,
}

impl TblastxDatabase {
    /// The OIDs of the database's ordinal ids, which name the reported sequences without
    /// parsed Seq-ids (`gnl|BL_ORD_ID|n`).
    pub(crate) fn oid_map(&self) -> crate::seqdb::OidMap<'_> {
        self.db.oid_map()
    }

    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(crate) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    /// Opens the nucleotide database `name` as NCBI's `CSearchDatabase` opens it
    /// (`common::database::open_database`). The total length is NCBI's `Int8`, which the engine
    /// holds as such (`SubjectSetNumbers::total`) on every target.
    #[cfg(test)]
    pub(crate) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        Self::open_restricted(name, settings, None, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line, and with the id of the database mask `mask` checked (`-db_soft_mask`,
    /// `-db_hard_mask`; series NP, M4b PR 2) (`shared::restrict_database`, whose warnings go to
    /// `diagnostics`).
    pub(crate) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        list: Option<&crate::seqdb::DbList>,
        mask: Option<&crate::blastinput::db_mask::DbMaskOption>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        // A database with parsed Seq-ids is read (line DBID, ID2b).
        let parsed =
            shared::restrict_database(&mut db, name, list, mask, diagnostics, settings, PROGRAM)?;
        // TBLASTX has no Gumbel block: the setup of a batch with a valid context divides by the
        // statistics' number of sequences, a count of 0 LOSAT's refusal there
        // (`context_cutoffs`, `shared::check_setup_count`).
        Ok(Self {
            db,
            name: name.to_string(),
            parsed,
            taxdb: shared::TaxdbRun::default(),
        })
    }

    /// The `-db` value (outfmt 7's `# Database:` line).
    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    /// Whether `taxdb.bti` is on the search path, asked once for the run
    /// (`shared::TaxdbRun::found`): the batches' prefetch reads the taxonomy database only then
    /// (`prefetch_messages`).
    pub(crate) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
    }

    /// The names of the taxonomy database that the rows of `spec` look up for the subjects
    /// `taxonomy` (`shared::tabular_names`).
    pub(crate) fn tabular_names(
        &self,
        spec: &crate::report::tabular::TabularSpec,
        taxonomy: &[crate::report::tabular::TabularTaxonomy],
    ) -> Result<Option<crate::report::tabular::TabularNames>> {
        shared::tabular_names(&self.db, &self.taxdb, spec, taxonomy, PROGRAM)
    }

    /// The open of the taxonomy database by the rows of a query batch (`shared::TaxdbRun::open_for_rows`).
    pub(crate) fn open_for_tabular_rows(&self) -> Result<()> {
        self.taxdb.open_for_rows(&self.db, PROGRAM)
    }

    /// Whether the database has parsed Seq-ids, whose subjects NCBI's prefetch can fail at
    /// (`batch_prefetch`).
    pub(crate) fn parsed_seqids(&self) -> bool {
        self.parsed
    }

    /// Whether the database has parsed Seq-ids and more than `max_sequences` sequences in its
    /// volumes (`SeqDb::num_oids`, also under a taxid list or an alias file's range: the range
    /// walks of `reports_may_refuse` decode the excluded sequences of every range they read):
    /// the check before the search does not read its headers
    /// (`shared::PARSED_CHECK_SEQUENCES_MAX`).
    pub(crate) fn parsed_beyond(&self, max_sequences: u32) -> bool {
        self.parsed && self.db.num_oids() > max_sequences
    }

    /// Whether a prefetch of the run is refused because of NCBI's taxonomy cache
    /// (`shared::TaxdbRun::prefetch_refused`), which turns the report stream off.
    pub(crate) fn prefetch_refused(&self) -> bool {
        self.taxdb.prefetch_refused(&self.db)
    }

    /// Whether a report that names any sequence of the database could be refused by LOSAT
    /// (`common::database::reports_may_refuse`, with the reports `shown` and the batches'
    /// prefetch `prefetch`), read in `pool`.
    pub(crate) fn reports_may_refuse(
        &self,
        shown: ParsedShown,
        prefetch: bool,
        queries: Option<&crate::blastinput::query_batch::QueryScope>,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Option<shared::ReportRefusal> {
        shared::reports_may_refuse(
            &self.db,
            shown,
            false,
            self.taxdb_found(),
            prefetch,
            queries,
            pool,
            PROGRAM,
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:164-168,183-187
    /// The numbers of the database that the statistics read (`SeqDb::statistics`: the setup's
    /// total and NCBI's `Int4` count, `BlastSeqSrcGetTotLen` for the linking, the average subject
    /// length), with `-dbsize`.
    pub(crate) fn numbers(&self, dbsize: Option<i64>) -> SubjectSetNumbers {
        let statistics = self.db.statistics();
        SubjectSetNumbers {
            total: statistics.setup_total,
            count: i64::from(statistics.setup_count),
            // NCBI's `Int8` of `GetTotalLength`: negative for a `LENGTH` of 2^63 or more, which
            // takes the subject's own branch of the link cutoffs (blast_engine.c:1451-1456).
            link_total: statistics.engine_db_length,
            average: Some(statistics.average),
            dbsize,
            searchsp: None,
            statistics: Some(statistics),
        }
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(crate) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
    }

    /// The numbers the statistics read, with NCBI's types (`SeqDb::statistics`); the BLAST
    /// archive's database length reads them too.
    pub(crate) fn statistics(&self) -> crate::seqdb::DbStatistics {
        self.db.statistics()
    }

    /// The taxids of the search's list (`SeqDb::list_taxids`).
    pub(crate) fn list_taxids(&self) -> &[i32] {
        self.db.list_taxids()
    }

    /// Whether the first query batch with a valid context fails with NCBI's `Database memory map
    /// file error` (`SeqDb::first_oid_read_fails`: no OID, and statistics that read OID 0).
    pub(crate) fn first_oid_read_fails(&self) -> bool {
        self.db.first_oid_read_fails()
    }

    /// The OIDs of the volumes (`SeqDb::num_oids`).
    #[cfg(test)]
    pub(crate) fn num_oids(&self) -> u32 {
        self.db.num_oids()
    }

    /// The `Database:` lines of the pairwise report.
    pub(crate) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// `pairwise_database` of a search masked with `mask` (`None`: not masked): each name's
    /// algorithm, which the `Masked using:` lines show (`shared::pairwise_database_masked`).
    pub(crate) fn pairwise_database_masked(
        &self,
        mask: Option<&crate::blastinput::db_mask::ResolvedMask>,
    ) -> PairwiseDatabase {
        shared::pairwise_database_masked(&self.db, mask.map(|mask| mask.shown.as_slice()))
    }

    /// The formatter's step of the database masking `option` (`shared::resolve_database_masking`;
    /// series NP, M4b PR 2): its warning goes to `diagnostics` (as one message of the archive
    /// `archive`), and the program options to `archive` and `strategy`.
    pub(crate) fn resolve_mask_option(
        &self,
        option: &crate::blastinput::db_mask::DbMaskOption,
        archive: Option<&crate::report::archive::ArchiveHandle>,
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Option<crate::blastinput::db_mask::ResolvedMask>> {
        shared::resolve_database_masking(
            &self.db,
            &self.name,
            option,
            PROGRAM,
            archive,
            strategy,
            &mut |text| match archive {
                Some(handle) => handle.write_warning(diagnostics, text),
                None => diagnostics.write_all(text),
            },
        )
    }

    /// LOSAT's refusals of the search masked with `mask` (`ResolvedMask::check_losat_limits`)
    /// whose words have `word_size` letters (the shortest subject searched has `word_size * 3 + 2`
    /// bases).
    pub(crate) fn check_mask_limits(
        &self,
        mask: &crate::blastinput::db_mask::ResolvedMask,
        word_size: usize,
    ) -> Result<()> {
        mask.check_losat_limits(
            &self.db,
            PROGRAM,
            crate::blastinput::db_mask::HardScan::Translated {
                min_length: word_size * 3 + 2,
            },
        )
    }

    /// The Blast-def-line taxid of each reported sequence by rank (`reported`: the ordinal id of
    /// each rank and the sequences with parsed Seq-ids; every rank is reported, a rank is made
    /// from a hit) and what the taxonomy database says about it, for the XML2 and JSON formats
    /// (`common::database::read_reported_taxonomy`: the taxid, `sciname`, and the taxonomy
    /// database's refusals before any byte). Only the ranks that the reports load are read
    /// (`ReportedSubjects::loads`: not the subject of a `Hit` that throws, nor those after it in a
    /// query's report, which NCBI never loads, so that the taxonomy database's one-time texts go
    /// with the next `Hit` that it loads; review 1 of DBID pfoff-1, H1).
    pub(crate) fn read_taxids(&self, reported: &ReportedSubjects) -> Result<XmlTaxa> {
        let ranks: std::collections::BTreeSet<usize> = (0..reported.ordinals.len())
            .filter(|rank| reported.loads.loaded(*rank))
            .collect();
        let (taxids, taxonomy, parsed_taxonomy, notes) = shared::read_reported_taxonomy(
            &self.db,
            &self.taxdb,
            &reported.ordinals,
            &ranks,
            &reported.parsed,
            PROGRAM,
        )?;
        Ok(XmlTaxa {
            taxids,
            taxonomy,
            parsed_taxonomy,
            notes,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-851
    /// The taxonomy of the tabular fields (`staxid`, `staxids`) of every reported sequence that a
    /// report loads (`ReportedSubjects::loads`), by rank (every rank is made from a hit,
    /// `reported_subjects`; `algorithm::common::tabular::tabular_taxonomy`). A header that does
    /// not decode is refused before the report.
    pub(crate) fn tabular_taxonomy(
        &self,
        reported: &ReportedSubjects,
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        crate::algorithm::common::tabular::tabular_taxonomy(
            &self.db,
            &reported.ordinals,
            &reported.parsed,
            0..reported.ordinals.len(),
            &reported.loads,
            PROGRAM,
        )
    }

    /// The database's stream (`common::x_db_stream::XDbStream`): the database checked as a copy
    /// of its sequences would be read, with the same errors, and its sequences read where each
    /// query batch searches them (`subjects::Subjects::Stream`) and reports them
    /// (`read_reported`).
    pub(crate) fn x_stream(
        &self,
        threads: usize,
        mask: Option<&crate::blastinput::db_mask::ResolvedMask>,
    ) -> Result<XDbStream<'_>> {
        Ok(XDbStream::new(&self.db, threads, PROGRAM)?.with_mask(
            crate::algorithm::common::db_mask_search::SearchMask::of(mask),
        ))
    }
}

/// What the XML2 and JSON formats show of a database sequence beyond its title: its
/// Blast-def-line taxid and the taxonomy database's answer, by rank (empty for outfmt 5 and
/// for `-subject`: every taxid is then 0), and NCBI's messages about the taxonomy database,
/// which are written at the first `Hit` that needs them.
#[derive(Default)]
pub(crate) struct XmlTaxa {
    pub(crate) taxids: Vec<i32>,
    pub(crate) taxonomy: Vec<crate::report::blast_xml::HitTaxonomy>,
    /// What the taxonomy database says about the taxid of each definition line of a sequence
    /// with parsed Seq-ids, by rank (empty as `taxonomy`).
    pub(crate) parsed_taxonomy: Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
    pub(crate) notes: shared::TaxdbNotes,
}

/// The headers of the database sequences that the reported hits of the query batches searched
/// so far name, by ordinal id (`read_titles`): the stored title, or the definition lines of a
/// sequence with parsed Seq-ids.
pub(crate) type DatabaseTitles = BTreeMap<u32, ReportedHeader>;

/// The final HSPs `ordered` of a query batch (report order) without those of the sequences that
/// NCBI makes no Seq-align of (`shared::SeqAligned`): after the cut of each query's hit list to
/// the hit list size, before `read_titles`.
pub(crate) fn drop_unaligned(database: &TblastxDatabase, ordered: &mut Vec<TblastxHsp>) {
    // tblastx fetches every subject with HSPs with `check_oid_exclusion` (the ungapped
    // re-evaluation of a nucleotide subject: `run_impl::search_query_batch`).
    let Some(mut aligned) = shared::SeqAligned::of(&database.db, true) else {
        return;
    };
    aligned.retain(ordered, |hsp| Some(hsp.hit.s_idx));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
/// Reads the definition lines of the sequences that the reported hits `ordered` of a query
/// batch name (in report order: query by query, in hit-list order; the subject index is the
/// ordinal id, the query index counts from the batch's first query `first_query`) into
/// `titles` (`common::database::read_header`). A reader error of NCBI's own is a
/// `HeaderReadError` (it fails the batch); LOSAT's refusals are errors of their own.
pub(crate) fn read_titles(
    database: &TblastxDatabase,
    ordered: &[TblastxHsp],
    first_query: usize,
    titles: &mut DatabaseTitles,
) -> Result<()> {
    let headers = shared::HeaderReader::new(&database.db);
    // A subject's HSPs follow each other: one look-up per run of them.
    let mut previous = None;
    for hsp in ordered {
        let ordinal = hsp.hit.s_idx;
        if previous.replace(ordinal) == Some(ordinal) || titles.contains_key(&ordinal) {
            continue;
        }
        let title = shared::read_header(
            &headers,
            ordinal,
            first_query + hsp.hit.q_idx as usize,
            PROGRAM,
        )?;
        titles.insert(ordinal, title);
    }
    Ok(())
}

/// The stores of the sequences that the reported hits of the query batches read so far name, by
/// ordinal id (`read_reported`).
pub(crate) type ReportedStores = BTreeMap<u32, ReportedStore>;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2717
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:322-370
// NCBI counts an HSP's identities and positives (the Seq-align's `num_ident`, `num_positives`)
// in its ungapped re-evaluation, on the translation of the subject that the re-evaluation
// fetched: with `-db_hard_mask`, a copy with `B` over the masked bases. The reports show the
// stored letters.
/// The letters of a reported database sequence (`read_reported`): those the reports show and
/// those the search's re-evaluation read, on which the engine's counts are made again
/// (`report::engine_counts`, `report::stored_engine_identities`).
#[derive(Clone)]
pub(crate) struct ReportedStore {
    /// The stored letters (`XDbStream::read_shared`): the rows of the reports.
    pub(crate) shown: Arc<PackedNucleotide>,
    /// The letters of the re-evaluation (`XDbStream::read_stage`: with `-db_hard_mask`, `B` over
    /// the masked bases; otherwise the same store as `shown`): the engine's identities and
    /// positives.
    pub(crate) engine: Arc<PackedNucleotide>,
}

impl ReportedStore {
    /// The letters of sequence `ordinal` of `stream`.
    fn read(stream: &XDbStream<'_>, ordinal: usize) -> Result<Self> {
        let shown = stream.read_shared(ordinal)?;
        // Without a hard mask this is `shown` again (the stream shares a store while it is held).
        let engine = stream.read_stage(ordinal)?;
        Ok(Self { shown, engine })
    }
}

/// Reads the sequences that the reported hits `ordered` of a query batch name into `stores`, once
/// each (`XDbStream::read_shared`: a long sequence that a search holds is not read again), for
/// the letters of the reports. A read that fails is an `XReadError` (it fails the batch).
pub(crate) fn read_reported(
    stream: &XDbStream<'_>,
    ordered: &[TblastxHsp],
    stores: &mut ReportedStores,
) -> Result<()> {
    // A subject's HSPs follow each other: one look-up per run of them.
    let mut previous = None;
    for hsp in ordered {
        let ordinal = hsp.hit.s_idx;
        if previous.replace(ordinal) == Some(ordinal) {
            continue;
        }
        if let std::collections::btree_map::Entry::Vacant(entry) = stores.entry(ordinal) {
            entry.insert(ReportedStore::read(stream, ordinal as usize)?);
        }
    }
    Ok(())
}

/// The lengths of the sequences that the reported hits of the query batches read so far name, by
/// ordinal id, when no report reads their letters (`read_reported_lengths`).
pub(crate) type ReportedLengths = BTreeMap<u32, u32>;

/// Reads the lengths of the sequences that the reported hits `ordered` of a query batch name into
/// `lengths`, once each (`XDbStream::read_length`: the index and the sequence's last byte, as
/// NCBI reads the length of a sequence it does not fetch), for the reports that read no letters
/// (`read_reported` reads the stores, which hold them, for the others). A read that fails is an
/// `XReadError` (it fails the batch).
pub(crate) fn read_reported_lengths(
    stream: &XDbStream<'_>,
    ordered: &[TblastxHsp],
    lengths: &mut ReportedLengths,
) -> Result<()> {
    // A subject's HSPs follow each other: one look-up per run of them.
    let mut previous = None;
    for hsp in ordered {
        let ordinal = hsp.hit.s_idx;
        if previous.replace(ordinal) == Some(ordinal) {
            continue;
        }
        if let std::collections::btree_map::Entry::Vacant(entry) = lengths.entry(ordinal) {
            let length = stream.read_length(ordinal as usize)?;
            entry.insert(u32::try_from(length).expect("a database sequence's length fits Int4"));
        }
    }
    Ok(())
}

/// NCBI's messages about the taxonomy database when the query batch whose reported hits are
/// `ordered` (subjects named by ordinal id, as `read_titles` reads them) prefetches the subjects
/// it reports (`shared::prefetch_messages`): the sequences the hits name, with the definition
/// lines of those with parsed Seq-ids from `titles` (`read_titles`). Empty, and nothing read,
/// without `taxdb.bti` on the search path. A taxonomy database that LOSAT cannot read as NCBI
/// does, and a header whose taxid does not decode, are refused.
pub(crate) fn prefetch_messages<'h>(
    database: &TblastxDatabase,
    ordered: impl IntoIterator<Item = &'h TblastxHsp>,
    titles: &DatabaseTitles,
) -> Result<Vec<u8>> {
    let subjects = ordered.into_iter().map(|hsp| {
        let ordinal = hsp.hit.s_idx;
        let parsed = match titles.get(&ordinal) {
            Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
            _ => None,
        };
        (ordinal, parsed)
    });
    shared::prefetch_messages(&database.db, &database.taxdb, subjects, PROGRAM)
}

/// [`prefetch_messages`] by subject (`shared::subject_read_texts`): the text of the first read
/// of each subject that writes one, with its ordinal id, where a query's report reads its
/// subjects one by one, the subjects `ordinals` in the order of its reads ([`read_ordinals`]).
pub(crate) fn subject_read_texts(
    database: &TblastxDatabase,
    ordinals: Vec<u32>,
    titles: &DatabaseTitles,
    rows_open: bool,
) -> Result<crate::report::query_warnings::QueryReads> {
    let parsed = |ordinal: u32| match titles.get(&ordinal) {
        Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
        _ => None,
    };
    let subjects = ordinals
        .into_iter()
        .map(|ordinal| (ordinal, ordinal, parsed(ordinal)));
    shared::subject_read_texts(&database.db, &database.taxdb, subjects, rows_open, PROGRAM)
}

/// The ordinal ids of the subjects that one query's report reads one by one, in the order of its
/// reads: those of its rows that read them, `shown` (loaded ones, in report order), or for the
/// pairwise report (`pairwise`: how it reads them) those that it shows of the query's hits
/// `ordered` in the order that it reads them (`pairwise_read_order`, `PairwiseReads::subjects`),
/// but those that it does not load (`loads`).
pub(crate) fn read_ordinals(
    shown: &[&TblastxHsp],
    ordered: &[TblastxHsp],
    titles: &DatabaseTitles,
    pairwise: Option<crate::report::pairwise::PairwiseReading>,
    loads: &HitLoads,
) -> Vec<u32> {
    match pairwise {
        Some(reading) => pairwise_read_order(ordered, titles, reading)
            .subjects()
            .into_iter()
            .filter(|&ordinal| loads.loaded(ordinal))
            .map(|ordinal| ordinal as u32)
            .collect(),
        None => shown.iter().map(|hsp| hsp.hit.s_idx).collect(),
    }
}

/// What one query's pairwise report (outfmt 0 to 4) of the hits `ordered` (in report order, by
/// ordinal id) reads (`QueryGroup::read_order` of the rows' hits as the writer reads them,
/// `reading`; a copy of a subject with parsed Seq-ids in `titles` joins the hit before; the rows'
/// hits are `report::PairwiseRows::hit`: the HSP's hit, whose `num_ident` is the engine's that the
/// rows' Seq-align identities copy, `report::seqalign_identities`; the rows' shift to the records
/// of a ranged query moves every HSP of the query alike, which no sort reads, and the rows name
/// the subjects by rank, in the order of the ordinal ids), for `read_ordinals` and `hit_loads`.
fn pairwise_read_order<'h>(
    ordered: impl IntoIterator<Item = &'h TblastxHsp>,
    titles: &DatabaseTitles,
    reading: crate::report::pairwise::PairwiseReading,
) -> crate::report::db_subject::PairwiseReads {
    let parsed = |ordinal: u32| titled_parsed(titles, ordinal as usize);
    crate::report::x_pairwise_par::QueryGroup::read_order(
        ordered.into_iter().map(|hsp| (hsp.hit.s_idx, hsp)),
        |hsp| crate::report::hit_sort::HspScores::of_hit(&hsp.hit, hsp.hit.num_ident),
        |previous, ordinal| match (parsed(previous), parsed(ordinal)) {
            (Some(earlier), Some(later)) => {
                crate::report::db_subject::joins_as_copy(earlier, later)
            }
            _ => false,
        },
        reading,
        true,
        None,
    )
}

/// The subjects of one query's report of the hits `hits` for `HitLoads::of` where the run's
/// formats are pairwise (`pairwise`: how the report reads them, with `after` a pairwise one): the
/// report's reads (`pairwise_read_order`). `None` otherwise.
fn pairwise_keys(
    hits: &[TblastxHsp],
    titles: &DatabaseTitles,
    after: Option<crate::report::db_subject::RunThrow>,
    pairwise: Option<crate::report::pairwise::PairwiseReading>,
) -> Option<crate::report::db_subject::ReportKeys> {
    let reading = pairwise.filter(|_| after.is_some_and(|run| run.after.is_pairwise()))?;
    Some(crate::report::db_subject::ReportKeys::Pairwise(
        pairwise_read_order(hits, titles, reading),
    ))
}

/// NCBI's prefetch of the query batch whose reported hits are `ordered` (in report order, by
/// ordinal id) on standard error (`shared::batch_prefetch`): the taxonomy texts of
/// [`prefetch_messages`], then the warning where the prefetch fails, at a subject whose lookup
/// fails or one without a Bioseq with a Seq-align whose first Std-seg has a minus query frame
/// and a plus subject frame (`CSeq_align::GetSeqStrand` reads the first segment): `ordered`
/// holds each query's HSPs of a subject in one run, in the order of the Std-segs (`seqalign`).
///
/// `before`: the batches of the run before it that the loop still holds (none where it wrote and
/// dropped them, `report_stream`, which records each batch instead: `note_loaded`), as NCBI's
/// data loader loaded them ([`LoadedBefore`]), read only where the batch stops; the outfmt 5
/// reports before the batch are the queries' (`shared::PrefetchQueries::after_xml_reports`).
pub(crate) fn batch_prefetch(
    database: &TblastxDatabase,
    before: LoadedBefore<'_>,
    ordered: &[TblastxHsp],
    titles: &DatabaseTitles,
    queries: &shared::PrefetchQueries<'_>,
) -> Result<crate::report::diag_log::Posts> {
    shared::batch_prefetch(
        &database.db,
        &database.taxdb,
        loaded_subjects(before, titles),
        prefetch_subjects(ordered, titles),
        PROGRAM,
        queries,
    )
}

/// The stop of NCBI's prefetch of the query batch of `queries` whose reported hits are
/// `ordered`, for the XML2 and JSON formats, which write the taxonomy texts at each `Hit`
/// (`shared::held_batch_prefetch_stop`): after the lookups of the reports of the batches before,
/// whose reported hits are `before` (by ordinal id; `before_queries` queries, those of each batch
/// `before_spans`), of the subjects that those reports load (`hit_loads` of the batches before
/// alone, of `batches_prefetch`, the reports' own answer, with `after` and the hits' engine
/// `identities`, `hit_identities`; review 2 of line DBID pfoff-1, N-M2), which TBLASTX writes
/// after the stop (the XML formats are not streamed, `report_stream`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn batch_prefetch_stop(
    database: &TblastxDatabase,
    before: &[TblastxHsp],
    before_queries: usize,
    before_spans: &[std::ops::Range<usize>],
    after: Option<crate::report::db_subject::RunThrow>,
    identities: &dyn Fn(&TblastxHsp) -> usize,
    ordered: &[TblastxHsp],
    titles: &DatabaseTitles,
    queries: &shared::PrefetchQueries<'_>,
) -> Result<()> {
    let prefetch = batches_prefetch(before, before_queries, before_spans);
    let loads = hit_loads(before, &prefetch, (after, None), titles, identities);
    shared::held_batch_prefetch_stop(
        &database.db,
        &database.taxdb,
        prefetch_subjects(before, titles).filter(|subject| loads.loaded(subject.ordinal as usize)),
        std::iter::empty(),
        prefetch_subjects(ordered, titles),
        PROGRAM,
        queries,
    )
}

/// The run's record of the sequences that NCBI's BLAST database data loader loaded for a query
/// batch whose reported hits are `ordered` (`shared::TaxdbRun::note_loaded`; [`loaded_subjects`]
/// of the batch alone, `rows`: the first format's), for a later batch's stop: where the loop
/// writes each batch's reports and drops its hits (`report_stream`), so that `batch_prefetch`'s
/// `before` no longer holds them. The stream never writes outfmt 5 (`report_stream::Off`), so no
/// report of the batch is cut (`shared::XmlCuts::NONE`); `queries` are the batch's queries as its
/// hits number them (`q_idx`, from 0 in the batch).
pub(crate) fn note_loaded(
    database: &TblastxDatabase,
    ordered: &[TblastxHsp],
    queries: std::ops::Range<usize>,
    titles: &DatabaseTitles,
    rows: SubjectFetch,
) {
    if database.parsed_seqids() {
        let batch = [queries];
        database.taxdb.note_loaded(&database.db, || {
            loaded_subjects(
                LoadedBefore {
                    hits: ordered,
                    spans: &batch,
                    rows,
                    cuts: shared::XmlCuts::NONE,
                },
                titles,
            )
        });
    }
}

/// The query batches of the run before a batch's prefetch that the loop holds
/// (`shared::LoadedBefore`, read by [`loaded_subjects`]): their reported hits in report order
/// (`hits`, by query of the run), the queries of each batch (`spans`), which rows of the first
/// format read a subject where a batch's prefetch is off (`rows`), and the outfmt 5 reports cut
/// before their hits (`cuts`).
pub(crate) type LoadedBefore<'a> = shared::LoadedBefore<'a, TblastxHsp>;

/// The sequences that NCBI's BLAST database data loader loaded for reported hits `ordered` whose
/// loads are `loads` (`shared::batch_loads`, review 3 of line DBID loadord, H2): every subject
/// where its batch's prefetch ran (TBLASTX's prefetched set is its reported set,
/// `shared::BatchLoads::reported_rows`); where it is off, those that the rows of a report NCBI
/// writes read (`ordered` is then one query's hits, `shared::LoadedBefore::subjects`), which
/// LOSAT takes as all of them unless no row reads one (`SubjectFetch::Never`): rows that read a
/// subject for its identities only (`nident`) read those of HSPs with identities, which the
/// stores of the batch count (`run_impl`'s fetches of a batch whose prefetch is off), so every
/// subject is a superset (TBLASTX's own: its hits carry no identities before those fetches).
fn batch_loaded<'s>(
    ordered: &'s [TblastxHsp],
    loads: shared::BatchLoads,
    titles: &'s DatabaseTitles,
) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
    let loaded = if loads.reported_rows() == SubjectFetch::Never {
        &[]
    } else {
        ordered
    };
    prefetch_subjects(loaded, titles)
}

/// The sequences that NCBI's BLAST database data loader loaded for the batches `before`
/// (`shared::LoadedBefore::subjects`, [`batch_loaded`] of each; none of an outfmt 5 report that
/// NCBI's XML stream cut, review 6 of line DBID loadord, H2).
fn loaded_subjects<'s>(
    before: LoadedBefore<'s>,
    titles: &'s DatabaseTitles,
) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
    before.subjects(
        |hits, span| Some(query_hits(hits, span)),
        |_, hits| seq_aligns_of(hits),
        move |_, hits, loads| batch_loaded(hits, loads, titles),
    )
}

/// The sequences that the prefetch of the query batch whose reported hits are `ordered` reads
/// ([`batch_prefetch`]).
fn prefetch_subjects<'s>(
    ordered: &'s [TblastxHsp],
    titles: &'s DatabaseTitles,
) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
    let mut previous = None;
    ordered.iter().map(move |hsp| {
        let ordinal = hsp.hit.s_idx;
        let parsed = match titles.get(&ordinal) {
            Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
            _ => None,
        };
        let run = (hsp.hit.q_idx, ordinal);
        let first = previous.replace(run) != Some(run);
        let flipped = first && hsp.hit.query_frame < 0 && hsp.subject_frame > 0;
        shared::PrefetchSubject {
            ordinal,
            parsed,
            flipped,
        }
    })
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(crate) struct ReportedSubjects {
    /// The ordinal id of each rank.
    pub(crate) ordinals: Vec<u32>,
    /// The records of the ranks: `gnl|BL_ORD_ID|n` and the stored title (a sequence with parsed
    /// Seq-ids: the title of its first definition line), without letters (the reports read them
    /// from the stores, `report::SubjectLetters::Database`).
    pub(crate) records: Vec<FastaRecord>,
    /// The number of bases of each rank.
    pub(crate) lengths: Vec<u32>,
    /// The stores of each rank (`read_reported`), or none when no report reads letters.
    pub(crate) stores: Vec<ReportedStore>,
    /// The sequences with parsed Seq-ids by rank (empty when none is reported): their ids and
    /// titles come from their definition lines (`report::db_subject`).
    pub(crate) parsed: Vec<Option<Arc<ParsedSubject>>>,
    /// The sources from which NCBI makes the empty stored titles of the others, by rank
    /// (`shared::StoredTitle::source`; empty when none has one).
    pub(crate) sources: Vec<crate::report::defline::SourceTaxname>,
    /// What the XML2 and JSON reports do at each rank's `Hit` (`HitLoads`, the one answer, by
    /// rank): the refusals, the names (`blast_xml`'s `XmlSubjects`) and the taxonomy reads
    /// (`read_taxids`) take it.
    pub(crate) loads: HitLoads,
}

/// What the XML2 and JSON reports of the reported hits `ordered` (in report order, subjects named
/// by ordinal id) do at each subject's `Hit` (`HitLoads::of`, the one answer, by ordinal id):
/// `prefetch` tells by the hits' query index whether its batch reads one by one
/// (`BatchPrefetch::reads_one_by_one`), `after` what the formats do after a subject whose lookup
/// throws (`FormatThrow::of_outfmts`), `identities` the engine's identities of a hit
/// (`hit_identities`), which tell whether a row of `nident` alone reads its subject; the `Hit`s
/// are those that `blast_xml` makes (the HSPs of each query, a subject's following each other),
/// or the pairwise reports' reads (`pairwise`: how they read them, `pairwise_keys`).
fn hit_loads(
    ordered: &[TblastxHsp],
    prefetch: &shared::BatchPrefetch,
    (after, pairwise): (
        Option<crate::report::db_subject::RunThrow>,
        Option<crate::report::pairwise::PairwiseReading>,
    ),
    titles: &DatabaseTitles,
    identities: &dyn Fn(&TblastxHsp) -> usize,
) -> HitLoads {
    let reports: Vec<&[TblastxHsp]> = ordered
        .chunk_by(|a, b| a.hit.q_idx == b.hit.q_idx)
        .collect();
    HitLoads::of(
        after,
        reports.iter().map(|hits| {
            // A batch that prefetches loads every subject, in any order.
            let one_by_one = prefetch.reads_one_by_one(hits[0].hit.q_idx as usize);
            let keys = pairwise_keys(hits, titles, after, pairwise.filter(|_| one_by_one))
                .unwrap_or_else(|| {
                    crate::report::db_subject::ReportKeys::Entries(
                        hits.iter()
                            .map(|hsp| Some(hsp.hit.s_idx as usize))
                            .collect(),
                    )
                });
            (one_by_one, keys)
        }),
        |ordinal| titled_parsed(titles, ordinal),
        |report, hit| identities(&reports[report][hit]),
    )
}

/// The subject with parsed Seq-ids of the ordinal id `ordinal` in `titles`.
fn titled_parsed(titles: &DatabaseTitles, ordinal: usize) -> Option<&ParsedSubject> {
    match titles.get(&u32::try_from(ordinal).ok()?) {
        Some(ReportedHeader::Parsed(subject)) => Some(subject),
        _ => None,
    }
}

/// The engine's identities of the hit `hsp` (its subject named by ordinal id) of the query whose
/// letters are `query` (`report::stored_engine_identities` against its subject's store in
/// `stores`, which a report with `nident` reads): whether a row of `nident` alone reads its
/// subject (`SubjectFetch::WithIdentities`). A report whose fields include `nident` reads the
/// letters of its subjects (`tabular::reads_letters`: `nident` is one of its `ROW_FIELDS`), so
/// the subjects that this is asked for have a store; one without a store would count as a subject
/// whose rows read it (1), which ends the run where NCBI's `nident` reads it (review 3 of line
/// DBID pfoff-2, L3).
pub(crate) fn hit_identities(
    stores: &ReportedStores,
    hsp: &TblastxHsp,
    query: &[u8],
    query_code: &crate::utils::genetic_code::GeneticCode,
    db_code: &crate::utils::genetic_code::GeneticCode,
) -> usize {
    stores.get(&hsp.hit.s_idx).map_or(1, |store| {
        super::report::stored_engine_identities(
            hsp,
            &store.engine,
            query,
            query_code,
            db_code,
            &mut Vec::new(),
        )
    })
}

/// What the formats of the run, `outfmts` as `-outfmt` gives them, do after a subject whose lookup
/// throws where its batch's prefetch is off (`FormatThrow::of_outfmts`, for `hit_loads` and the
/// run's `RunLoads`): the stop of a batch's prefetch (`batch_prefetch_stop`) and the reports
/// (`reported_subjects`) ask it of the run's formats alike (review 2 of line DBID pfoff-1, N-M2).
pub(crate) fn after_throw<'a>(
    outfmts: impl IntoIterator<Item = &'a str>,
) -> Option<crate::report::db_subject::RunThrow> {
    crate::report::db_subject::FormatThrow::of_outfmts(
        outfmts,
        crate::blastinput::app::Program::Tblastx,
        crate::report::archive::options::RequestKind::Tblastx,
    )
}

/// How the pairwise reports of the formats `formats` read each query's subjects, as the first
/// pairwise format's writer reads them (`PairwiseView::reading` of `-sorthits` and `-sorthsps`
/// `sort` and the run's numbers `counts`; the formats of a run whose loads they make share one
/// view's order, `FormatThrow::and`); `None` without one.
pub(crate) fn pairwise_reading(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
    sort: crate::report::hit_sort::ReportSort,
    counts: crate::blastinput::report_counts::ReportCounts,
) -> Option<crate::report::pairwise::PairwiseReading> {
    let view = formats.iter().find_map(|format| {
        format
            .report_format()
            .and_then(crate::report::pairwise::PairwiseView::of)
    })?;
    Some(view.reading(sort, counts.shown()))
}

/// What one query's report of the hits `ordered` (in report order) loads where its batch's
/// prefetch is off (`RunLoads::loads` of that report alone, by ordinal id, after the reports
/// before it in the run, `run`; `identities`: [`hit_identities`]): the subjects of its reads one
/// by one (`run_impl`'s texts of a batch whose prefetch is off), not one whose load throws, nor
/// one that no report reaches after it. `pairwise`: how a pairwise report reads them
/// (`pairwise_keys`).
pub(crate) fn report_loads(
    run: &crate::report::db_subject::RunLoads,
    ordered: &[TblastxHsp],
    titles: &DatabaseTitles,
    pairwise: Option<crate::report::pairwise::PairwiseReading>,
    identities: impl Fn(&TblastxHsp) -> usize,
) -> HitLoads {
    // A subject's HSPs follow each other: one `Hit`.
    let subjects: Vec<&[TblastxHsp]> = ordered
        .chunk_by(|a, b| a.hit.s_idx == b.hit.s_idx)
        .collect();
    let keys = pairwise_keys(ordered, titles, run.run(), pairwise).unwrap_or_else(|| {
        crate::report::db_subject::ReportKeys::Entries(
            subjects
                .iter()
                .map(|hsps| Some(hsps[0].hit.s_idx as usize))
                .collect(),
        )
    });
    run.loads(
        [(true, keys)],
        |ordinal| titled_parsed(titles, ordinal),
        |_, subject| subjects[subject].iter().map(&identities).max().unwrap_or(0),
    )
}

/// The Seq-aligns of the hits `ordered` (in report order) that NCBI's prefetch counts against
/// `PRE_FETCH_SEQS_LIMIT` (`shared::prefetch_off`): one for each subject of each query (TBLASTX's
/// search is ungapped), whose HSPs follow each other.
pub(crate) fn seq_aligns_of(ordered: &[TblastxHsp]) -> usize {
    ordered
        .chunk_by(|a, b| a.hit.q_idx == b.hit.q_idx && a.hit.s_idx == b.hit.s_idx)
        .map(|hsps| shared::subject_seq_aligns(hsps.len(), false))
        .sum()
}

/// The hits of `ordered` (in report order, by query) of the queries `queries`.
pub(crate) fn query_hits(ordered: &[TblastxHsp], queries: std::ops::Range<usize>) -> &[TblastxHsp] {
    let from = ordered.partition_point(|hsp| (hsp.hit.q_idx as usize) < queries.start);
    let to = ordered.partition_point(|hsp| (hsp.hit.q_idx as usize) < queries.end);
    &ordered[from..to.max(from)]
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:728-753 (cited at
// `common::database::prefetch_may_be_off`: `s_PreFetchSeqs` counts the Seq-aligns of one batch's
// results against `PRE_FETCH_SEQS_LIMIT`)
/// Whether the batch of each query of the reported hits `ordered` (in report order, by query of
/// `queries` queries, numbered from 0) prefetches (`shared::BatchPrefetch::of`): `spans` are the
/// queries of each batch, whose Seq-aligns `seq_aligns_of` counts. The one producer of the
/// prefetch that the reports' loads read (`hit_loads`), for both of its readers: the stop of a
/// later batch, of the batches before it (`batch_prefetch_stop`), and the reports
/// (`reported_subjects`), so the two cannot disagree on what the batches before loaded (review 2
/// of line DBID pfoff-1, N-M2).
fn batches_prefetch(
    ordered: &[TblastxHsp],
    queries: usize,
    spans: &[std::ops::Range<usize>],
) -> shared::BatchPrefetch {
    shared::BatchPrefetch::of(queries, 0, spans, |span| {
        seq_aligns_of(query_hits(ordered, span))
    })
}

/// The reported database sequences of `ordered` (the reported hits of all query batches in
/// report order, subjects named by ordinal id): the hits name them by rank from here on (the
/// rank grows with the ordinal id, so every comparison that breaks a tie by the subject index
/// gives the same order), and their records (`titles`, by ordinal id), lengths and,
/// when a report reads letters, stores (`stores`, `read_reported`, which give the lengths; else
/// `lengths`, `read_reported_lengths`) are kept for the reports. A reported sequence whose stored title LOSAT does not
/// report as NCBI does is LOSAT's refusal, and so is one with parsed Seq-ids that the reports
/// `shown` do not write as NCBI (`common::database::check_reported_parsed`, and with
/// `shown.joined` `refuse_joined_hits`; where only the prefetch reads a subject, a batch whose
/// prefetch is on reports it, `batches_prefetch` of the hits' `queries` queries and the queries of
/// each batch `spans`; `after`: what the formats do after a subject whose lookup throws and how
/// the pairwise reports read the subjects, for the one answer of what the reports load, `hit_loads`, with the
/// hits' engine `identities`, `hit_identities`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn reported_subjects(
    database: &TblastxDatabase,
    titles: &DatabaseTitles,
    ordered: &mut [TblastxHsp],
    stores: Option<&ReportedStores>,
    read_lengths: &ReportedLengths,
    shown: ParsedShown,
    after: (
        Option<crate::report::db_subject::RunThrow>,
        Option<crate::report::pairwise::PairwiseReading>,
    ),
    identities: &dyn Fn(&TblastxHsp) -> usize,
    queries: usize,
    spans: &[std::ops::Range<usize>],
) -> Result<ReportedSubjects> {
    let mut ordinals: Vec<u32> = ordered.iter().map(|hsp| hsp.hit.s_idx).collect();
    ordinals.sort_unstable();
    ordinals.dedup();
    let prefetch = batches_prefetch(ordered, queries, spans);
    // What the batches read of each subject (`prefetch`: a batch whose Seq-aligns outnumber a
    // positive `PRE_FETCH_SEQS_LIMIT` reads nothing of them before its reports, which read them
    // one by one; `shared::BatchReads::shown`); where the prefetch fails, the batch writes
    // NCBI's warning (`batch_prefetch`).
    let reads = prefetch.reads(
        ordered
            .iter()
            .map(|hsp| (hsp.hit.q_idx as usize, hsp.hit.s_idx)),
    );
    // What the reports do at each subject's `Hit`, by ordinal id (the one answer, `hit_loads`).
    let loads = hit_loads(ordered, &prefetch, after, titles, identities);
    for hsp in ordered.iter_mut() {
        let rank = ordinals
            .binary_search(&hsp.hit.s_idx)
            .expect("every reported subject was listed");
        hsp.hit.s_idx = u32::try_from(rank).context("TBLASTX subject rank out of range")?;
    }
    let mut records = Vec::with_capacity(ordinals.len());
    let mut lengths = Vec::with_capacity(ordinals.len());
    let mut kept = Vec::new();
    let mut parsed: Vec<Option<Arc<ParsedSubject>>> = Vec::new();
    let mut sources = Vec::new();
    for (rank, &ordinal) in ordinals.iter().enumerate() {
        let header = titles.get(&ordinal).with_context(|| {
            format!(
                "database sequence {} has no title read",
                database.db.oid(ordinal)
            )
        })?;
        let title = match header {
            ReportedHeader::Stored(stored) => {
                shared::check_reported_title(
                    stored.replaced_bytes,
                    || database.db.oid(ordinal),
                    PROGRAM,
                )?;
                if shown.pairwise {
                    stored.check_shown_source(|| database.db.oid(ordinal), PROGRAM)?;
                }
                shared::keep_source(&mut sources, rank, ordinals.len(), stored.source.clone());
                stored.title.clone()
            }
            ReportedHeader::Parsed(subject) => {
                // Whether `taxdb.bti` is on the search path, looked up at the first sequence
                // with parsed Seq-ids that outfmt 0 shows (once for the run, `taxdb_found`).
                let taxdb = shown.pairwise && database.taxdb_found();
                shared::check_reported_parsed(
                    &database.db,
                    subject,
                    reads.shown(shown, ordinal),
                    loads.load(ordinal as usize),
                    false,
                    taxdb,
                    ordinal,
                    PROGRAM,
                )?;
                if parsed.is_empty() {
                    parsed = vec![None; ordinals.len()];
                }
                parsed[rank] = Some(Arc::clone(subject));
                subject.first_title().to_vec()
            }
        };
        lengths.push(match stores {
            Some(stores) => {
                let store = stores.get(&ordinal).with_context(|| {
                    format!(
                        "database sequence {} was not read",
                        database.db.oid(ordinal)
                    )
                })?;
                kept.push(store.clone());
                u32::try_from(store.shown.len()).expect("a database sequence's length fits Int4")
            }
            None => *read_lengths.get(&ordinal).with_context(|| {
                format!(
                    "database sequence {} was not read",
                    database.db.oid(ordinal)
                )
            })?,
        });
        records.push(FastaRecord {
            // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the OID
            // (line DI-8).
            local_id: database.db.oid(ordinal).fasta_id(),
            title,
            sequence: Vec::new(),
            warnings: Vec::new(),
            source: Default::default(),
        });
    }
    if shown.joined && !parsed.is_empty() {
        // The subjects of each query's hit list in its order: the HSPs of a subject follow
        // each other in report order.
        let mut start = 0;
        while start < ordered.len() {
            let q_idx = ordered[start].hit.q_idx;
            let end = start
                + ordered[start..]
                    .iter()
                    .take_while(|hsp| hsp.hit.q_idx == q_idx)
                    .count();
            let mut ranks: Vec<u32> = ordered[start..end]
                .iter()
                .map(|hsp| hsp.hit.s_idx)
                .collect();
            ranks.dedup();
            shared::refuse_joined_hits(
                q_idx as usize,
                ranks.iter().map(|&rank| {
                    let subject = parsed[rank as usize].as_deref()?;
                    Some((ordinals[rank as usize], subject))
                }),
                database.db.oid_map(),
                PROGRAM,
            )?;
            start = end;
        }
    }
    let loads = loads.renamed(|ordinal| ordinals.binary_search(&u32::try_from(ordinal).ok()?).ok());
    Ok(ReportedSubjects {
        ordinals,
        records,
        lengths,
        stores: kept,
        parsed,
        sources,
        loads,
    })
}

impl ReportedSubjects {
    /// The subject ids of the tabular reports (`saccver`), by rank: the first word of the stored
    /// title or the general id's `BL_ORD_ID:n` (`oids`), or the id of a sequence with parsed
    /// Seq-ids (`common::database::tabular_subject_ids`; an id LOSAT cannot write is empty where no
    /// tabular report `reads` the subject).
    pub(crate) fn tabular_subject_ids(
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
    pub(crate) fn library_subjects(
        &self,
        database: &TblastxDatabase,
    ) -> Vec<crate::api::local_blast::ReportedSubject> {
        shared::library_subjects(
            &database.db,
            &self.records,
            &self.ordinals,
            &self.parsed,
            &self.sources,
            |rank| self.lengths[rank] as usize,
            PROGRAM,
        )
    }
}

#[cfg(test)]
#[path = "database_tests.rs"]
mod tests;
