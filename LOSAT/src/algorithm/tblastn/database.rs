//! The database path of TBLASTN (`-db NAME`): the search of an NCBI BLAST nucleotide
//! database, as NCBI BLAST+ 2.17.0 `tblastn -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, stage DB-5).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the subjects are the database's sequences, each read where it is searched, as packed
//!   subjects (`TblastnDatabase::x_stream`, `common/x_db_stream.rs`): the preliminary stage translates the stored
//!   2-bit data (makeblastdb's bases for the ambiguous ones, NCBI's preliminary input) and
//!   the traceback and the reports read the stored ambiguity codes; the search's subject
//!   index is the ordinal id;
//! - the set-level numbers come from the index files (`TblastnDatabase::numbers`): the
//!   number of sequences, the total length (an `Int8`) and the shortest length a database
//!   answers (`BLAST_SEQSRC_MINLENGTH`; `SubjectSetNumbers`);
//! - the definition lines of the final hit lists are read after each query batch's search
//!   (`read_titles`), and the reports show the database (`TblastnDatabase::name`,
//!   `pairwise_database`) and the stored titles of the reported sequences only
//!   (`reported_subjects`).
//!
//! The search, the hit lists and the report writers are those of the `-subject` path; the
//! opening of the database, the reader's errors and the rules for the stored titles are
//! shared with BLASTP and BLASTN (`algorithm::common::database`).

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::algorithm::common::database::{self as shared, ReportedHeader};
use crate::algorithm::common::x_db_stream::XDbStream;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::db_subject::{HitLoads, ParsedSubject};
use crate::report::pairwise::PairwiseDatabase;
use crate::report::tabular::SubjectFetch;
use crate::seqdb::{SeqDb, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

use super::stage_d_pipeline::DatabaseNumbers;
use super::stage_d_results::{KappaResultHitList, KappaResultList};

/// The program that LOSAT's refusals name.
const PROGRAM: &str = "TBLASTN";

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`TblastnArgs::database_name`: the `-db` value without white space at its ends).
pub(super) struct TblastnDatabase {
    pub(crate) db: SeqDb,
    name: Vec<u8>,
    /// The taxonomy database of the run, opened by the XML2 and JSON formats and by the
    /// batches' prefetch of outfmt 0, 5, 6 and 7 (`prefetch_messages`).
    taxdb: shared::TaxdbRun,
    /// Whether the database has parsed Seq-ids (`shared::check_seqids`), whose subjects NCBI's
    /// prefetch can fail at (`batch_prefetch`).
    parsed_seqids: bool,
}

impl TblastnDatabase {
    /// The OIDs of the database's ordinal ids, which name the reported sequences without
    /// parsed Seq-ids (`gnl|BL_ORD_ID|n`).
    pub(super) fn oid_map(&self) -> crate::seqdb::OidMap<'_> {
        self.db.oid_map()
    }

    /// Whether a load of a sequence with a positive taxid can write a text or throw in this run
    /// (`shared::TaxdbRun::loads_write_texts`).
    pub(super) fn loads_write_texts(&self) -> bool {
        self.taxdb.loads_write_texts(&self.db)
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
    pub(super) fn keeping_scope(
        mut self,
        keeps: bool,
        separate_xml2: bool,
        record_keeps: Option<bool>,
    ) -> Self {
        self.taxdb.set_keeps_scope(keeps);
        self.taxdb.set_separate_xml2(separate_xml2);
        self.taxdb.set_record_keeps_scope(record_keeps);
        self
    }

    /// Mirrors `CBlastFormat::ResetScopeHistory` at the end of a query batch
    /// (`shared::TaxdbRun::reset_scope`).
    pub(super) fn reset_scope(&self) {
        self.taxdb.reset_scope();
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:188-194
    /// Opens the database of `-db name` (`algorithm::common::database::open_database`) with
    /// the command line's errors. The total length is NCBI's `Int8`, which the engine holds as
    /// such (`DatabaseNumbers::total_length`) on every target.
    #[cfg(test)]
    pub(super) fn open(name: impl AsRef<[u8]>, settings: &ApplicationSettings) -> Result<Self> {
        Self::open_restricted(name.as_ref(), settings, None, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line, and with the id of the database mask `mask` checked (`-db_soft_mask`,
    /// `-db_hard_mask`; series NP, M4b PR 2) (`shared::restrict_database`, whose warnings go to
    /// `diagnostics`).
    pub(super) fn open_restricted(
        name: &[u8],
        settings: &ApplicationSettings,
        list: Option<&crate::seqdb::DbList>,
        mask: Option<&crate::blastinput::db_mask::DbMaskOption>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        // A database with parsed Seq-ids is read (line DBID, ID2b).
        let parsed_seqids = shared::restrict_database(
            &mut db,
            name,
            list,
            mask,
            diagnostics,
            settings,
            PROGRAM,
            PROGRAM,
        )?;
        Ok(Self {
            db,
            name: name.to_vec(),
            taxdb: shared::TaxdbRun::default(),
            parsed_seqids,
        })
    }

    /// Whether NCBI's prefetch of a batch can fail at a subject (a database with parsed
    /// Seq-ids, `batch_prefetch`).
    pub(super) fn parsed_seqids(&self) -> bool {
        self.parsed_seqids
    }

    /// The `-db` value as NCBI keeps it (outfmt 7's `# Database:` line).
    pub(super) fn name(&self) -> &[u8] {
        &self.name
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:182-187
    /// The numbers of the statistics (`SeqDb::statistics`): the setup's total length and number
    /// of sequences, the Gumbel block's length and `CSeqDB::GetMinLength`; the real numbers, also
    /// with `-dbsize`.
    pub(super) fn numbers(&self) -> DatabaseNumbers {
        let statistics = self.db.statistics();
        DatabaseNumbers {
            // The setup's total is not negative (a length of OID 0 where the total is not
            // positive).
            total_length: statistics.setup_total,
            num_seqs: statistics.setup_count,
            gumbel_db_length: statistics.gumbel_db_length,
            engine_db_length: statistics.engine_db_length,
            min_length: statistics.min_length,
            request_count: statistics.request_count,
            average: statistics.average,
        }
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(super) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
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

    /// Whether the first query batch with a valid context fails with NCBI's `Database memory map
    /// file error` (`SeqDb::first_oid_read_fails`: no OID, and statistics that read OID 0).
    pub(super) fn first_oid_read_fails(&self) -> bool {
        self.db.first_oid_read_fails()
    }

    /// What a PSSM search's first read of the database meets (`shared::PsiSearchReadFailure`:
    /// NCBI's `CPsiBlast` gives SeqDB's own error).
    pub(super) fn pssm_search_read_failure(&self) -> Option<shared::PsiSearchReadFailure> {
        shared::PsiSearchReadFailure::of(&self.db)
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`).
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

    /// The formatter's step of the database masking `option` (`shared::resolve_database_masking`;
    /// series NP, M4b PR 2): its warning goes to `diagnostics` (as one message of the archive
    /// `archive`), and the program options to `archive` and `strategy`.
    pub(super) fn resolve_mask_option(
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
            PROGRAM.to_ascii_lowercase().as_str(),
            archive,
            strategy,
            &mut |text| match archive {
                Some(handle) => handle.write_warning(diagnostics, text),
                None => diagnostics.write_all(text),
            },
        )
    }

    /// LOSAT's refusals of the search masked with `mask` (`ResolvedMask::check_losat_limits`)
    /// for word size `word`, two-hit window `window` (0 for one-hit), and whether the lookup table
    /// can hold a word with residue 0 (`zero_words`).
    pub(super) fn check_mask_limits(
        &self,
        mask: &crate::blastinput::db_mask::ResolvedMask,
        word: i32,
        window: i32,
        zero_words: bool,
    ) -> Result<()> {
        mask.check_losat_limits(
            &self.db,
            PROGRAM,
            crate::blastinput::db_mask::HardScan::Translated {
                word,
                window,
                zero_words,
            },
        )
    }

    /// The Blast-def-line taxid of each reported sequence by rank and what the taxonomy
    /// database says about it, for the XML2 and JSON formats
    /// (`algorithm::common::database::read_reported_taxonomy`, with the definition lines of the
    /// sequences with parsed Seq-ids); the sequences of `hitlists` with HSPs are the reported
    /// ones (their subject index is the rank, `reported_subjects`), and those that the reports
    /// load are read (`ReportedSubjects::loads`: not the subject of a `Hit` that throws, nor those
    /// after it in a query's report, which NCBI never loads; review 1 of DBID pfoff-1, H1). The
    /// texts of the loads are the batches' prefetches' or reads', modelled before this call.
    pub(super) fn read_taxids(
        &self,
        reported: &ReportedSubjects,
        hitlists: &[KappaResultHitList],
    ) -> Result<XmlTaxa> {
        let ranks: std::collections::BTreeSet<usize> = reported_ranks(hitlists)
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

    /// The taxonomy of the tabular fields (`staxid`, `staxids`) of the reported sequences, by
    /// rank, for the sequences of `hitlists` with HSPs (their subject index is the rank,
    /// `reported_subjects`; `algorithm::common::tabular::tabular_taxonomy`) that a report loads
    /// (`ReportedSubjects::loads`); the others are left empty. A header that does not decode is
    /// refused before the report.
    pub(super) fn tabular_taxonomy(
        &self,
        reported: &ReportedSubjects,
        hitlists: &[KappaResultHitList],
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        crate::algorithm::common::tabular::tabular_taxonomy(
            &self.db,
            &reported.ordinals,
            &reported.parsed,
            hitlists
                .iter()
                .flat_map(KappaResultHitList::lists)
                .filter(|list| !list.hsps.hsps.is_empty())
                .map(|list| list.oid as usize),
            &reported.loads,
            PROGRAM,
        )
    }

    /// Whether `taxdb.bti` is on the search path, asked once for the run
    /// (`shared::TaxdbRun::found`): the batches' prefetch reads the taxonomy database only then
    /// (`prefetch_messages`).
    pub(super) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:750-763
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:877-890
    /// The open of the taxonomy database by the rows `rows` of the final hit lists `hitlists`
    /// (subjects named by rank, `reported_subjects`) after their batch's prefetch
    /// (`shared::RowLookups::open_after_prefetch`): the rows of a subject that the reports load
    /// (`HitLoads::loaded`; a row whose load throws looks nothing up). With `ends_on_throw`
    /// (`AfterThrow::EndsReport` or `EndsRun`), a query's subjects stop at the first whose load
    /// throws: its report loads nothing after it (`shared::one_by_one_query_reads`).
    pub(super) fn open_for_rows(
        &self,
        reported: &ReportedSubjects,
        hitlists: &[KappaResultHitList],
        rows: shared::RowLookups,
        ends_on_throw: bool,
    ) -> Result<()> {
        let mut subjects = Vec::new();
        for hitlist in hitlists {
            for list in hitlist.lists() {
                if list.hsps.hsps.is_empty() {
                    continue;
                }
                let rank = list.oid as usize;
                if ends_on_throw && reported.loads.stops_report(rank) {
                    break;
                }
                if reported.loads.loaded(rank) {
                    subjects.push((reported.ordinals[rank], reported.parsed_subject(rank)));
                }
            }
        }
        rows.open_after_prefetch(&self.db, &self.taxdb, subjects, PROGRAM)
    }

    /// NCBI's messages about the taxonomy database when the query batch whose final hit lists
    /// are `hitlists` (subjects named by rank, `reported_subjects`) prefetches the subjects it
    /// reports (`shared::prefetch_messages`): the sequences with HSPs in them, the ones
    /// `read_taxids` reads, with the definition lines of those with parsed Seq-ids; where a
    /// query's report fetches them (`PRE_FETCH_SEQS_LIMIT`), those with an HSP whose row reads
    /// them (`rows`, by the engine's identities). Empty, and nothing read, without `taxdb.bti` on
    /// the search path. A taxonomy database that LOSAT cannot read as NCBI does, and a header
    /// whose taxid does not decode, are refused.
    pub(super) fn prefetch_messages(
        &self,
        reported: &ReportedSubjects,
        hitlists: &[KappaResultHitList],
        rows: SubjectFetch,
    ) -> Result<Vec<u8>> {
        let subjects = hitlists
            .iter()
            .flat_map(KappaResultHitList::lists)
            .filter(|list| reads_subject(list, rows))
            // Not a subject that no report reaches after one whose load throws
            // (`HitLoads::reached`).
            .filter(|list| reported.loads.reached(list.oid as usize))
            .map(|list| {
                let rank = list.oid as usize;
                (reported.ordinals[rank], reported.parsed_subject(rank))
            });
        shared::prefetch_messages(&self.db, &self.taxdb, subjects, PROGRAM)
    }

    /// [`Self::prefetch_messages`] by subject (`shared::subject_read_texts`): the text of the
    /// first read of each subject that writes one, with its rank, where a query's report reads
    /// its subjects one by one: in the order of the hit lists, or, for the pairwise report
    /// (`pairwise`: how it reads them), those that it shows in the order that it reads them
    /// (`pairwise_read_order`, `PairwiseReads::subjects`).
    pub(super) fn subject_read_texts(
        &self,
        reported: &ReportedSubjects,
        hitlists: &[KappaResultHitList],
        rows: SubjectFetch,
        row_lookups: shared::RowLookups,
        pairwise: Option<crate::report::pairwise::PairwiseReading>,
        (joins_copies, ends_on_throw): (bool, bool),
    ) -> Result<crate::report::query_warnings::QueryReads> {
        let lists = hitlists
            .iter()
            .flat_map(KappaResultHitList::lists)
            .filter(|list| reads_subject(list, rows));
        let ranks: Vec<usize> = match pairwise {
            Some(reading) => pairwise_read_order(
                lists,
                |rank| Ok(reported.lengths[rank]),
                |previous, rank| {
                    crate::report::db_subject::rank_joins_as_copy(
                        &reported.parsed,
                        previous as usize,
                        rank as usize,
                    )
                },
                reading,
            )?
            .subjects()
            .into_iter()
            .filter(|&rank| reported.loads.reached(rank))
            .collect(),
            None => shared::one_by_one_query_reads(
                lists.map(|list| list.oid as usize),
                |rank| reported.parsed_subject(rank),
                &reported.loads,
                joins_copies,
                ends_on_throw,
            ),
        };
        let subjects = ranks.into_iter().map(|rank| {
            (
                rank as u32,
                reported.ordinals[rank],
                reported.parsed_subject(rank),
            )
        });
        shared::subject_read_texts(&self.db, &self.taxdb, subjects, row_lookups, PROGRAM)
    }

    /// The stop of the prefetch of the batch of `queries` whose reported
    /// subjects are `subjects` (`prefetch_stop`): NCBI's error (a `ErrorAfterBatchReports`), or
    /// LOSAT's refusal (`shared::batch_prefetch`, after the batches before, `before`: the run's
    /// record of what NCBI's BLAST database data loader loaded for them, `LoadedBefore`).
    /// The outfmt 5 reports before the batch are the queries'
    /// (`shared::PrefetchQueries::after_xml_reports`).
    pub(super) fn prefetch_stop_error(
        &self,
        reported: Option<&ReportedSubjects>,
        before: LoadedBefore<'_>,
        subjects: &[(u32, Option<Arc<ParsedSubject>>)],
        queries: &shared::PrefetchQueries<'_>,
    ) -> anyhow::Error {
        let subjects = subjects.iter().map(|(ordinal, parsed)| {
            shared::PrefetchSubject::unflipped(*ordinal, parsed.as_deref())
        });
        let before = reported
            .into_iter()
            .flat_map(|reported| Self::loaded_subjects(reported, before));
        let stopped =
            shared::batch_prefetch(&self.db, &self.taxdb, before, subjects, PROGRAM, queries).err();
        stopped.expect("`prefetch_stop` found a subject whose lookup stops the batch's prefetch")
    }

    /// The sequences that the prefetch of a query batch whose final hit lists are `hitlists`
    /// (subjects named by rank, `reported_subjects`) reads: those with HSPs in them (the protein
    /// query has no strand), or, as its reports read them where the prefetch is off, those whose
    /// rows read them (`rows`, `reads_subject`).
    fn prefetch_subjects<'s>(
        reported: &'s ReportedSubjects,
        hitlists: &'s [KappaResultHitList],
        rows: SubjectFetch,
    ) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
        hitlists
            .iter()
            .flat_map(KappaResultHitList::lists)
            .filter(move |list| reads_subject(list, rows))
            .map(|list| {
                let rank = list.oid as usize;
                shared::PrefetchSubject::unflipped(
                    reported.ordinals[rank],
                    reported.parsed_subject(rank),
                )
            })
    }

    /// The sequences that NCBI's BLAST database data loader loaded for the query batches `before`
    /// (`shared::LoadedBefore::subjects`: the subjects of each batch's Seq-aligns where its
    /// prefetch is on, else those that the rows of the reports NCBI writes read,
    /// `shared::BatchLoads::reported_rows`; review 3 of line DBID loadord, H2; review 6, H2).
    fn loaded_subjects<'s>(
        reported: &'s ReportedSubjects,
        before: LoadedBefore<'s>,
    ) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
        let gapped = reported.gapped;
        before.subjects(
            |hitlists, span| hitlists.get(span),
            move |_, hitlists| seq_aligns_of(hitlists, gapped),
            // TBLASTN's prefetched set is its reported set.
            move |_, hitlists, loads| {
                Self::prefetch_subjects(reported, hitlists, loads.reported_rows())
            },
        )
    }

    /// NCBI's prefetch of the query batch whose final hit lists are `hitlists` (subjects named by
    /// rank, `reported_subjects`) on standard error (`shared::batch_prefetch`): the taxonomy
    /// texts of the sequences with HSPs in them, then the warning where the prefetch fails at a
    /// subject whose lookup fails (the protein query has no strand); the stop of the batch of
    /// `queries`, after the batches before, `before` (read only where the batch stops; the outfmt
    /// 5 reports before it, `shared::PrefetchQueries::after_xml_reports`).
    pub(super) fn batch_prefetch(
        &self,
        reported: &ReportedSubjects,
        before: LoadedBefore<'_>,
        hitlists: &[KappaResultHitList],
        queries: &shared::PrefetchQueries<'_>,
    ) -> Result<crate::report::diag_log::Posts> {
        shared::batch_prefetch(
            &self.db,
            &self.taxdb,
            Self::loaded_subjects(reported, before),
            Self::prefetch_subjects(reported, hitlists, SubjectFetch::Every),
            PROGRAM,
            queries,
        )
    }

    /// Whether a reader met a reported sequence whose Seq-align id's lookup reads a header with a
    /// UTF8String where a VisibleString is expected (`SeqDb::scope_utf8_seen`): the loads of the
    /// reports' reads are modelled for its warning also without `taxdb.bti`.
    pub(super) fn scope_utf8_seen(&self) -> bool {
        self.db.scope_utf8_seen()
    }

    /// The database read where it is searched by `threads` threads, after the checks of
    /// `common::database::x_check_nucleotide_stores` (`common/x_db_stream.rs`).
    pub(super) fn x_stream(
        &self,
        threads: usize,
        mask: Option<&crate::blastinput::db_mask::ResolvedMask>,
    ) -> Result<XDbStream<'_>> {
        Ok(XDbStream::new(&self.db, threads, PROGRAM)?
            .with_mask(crate::algorithm::common::db_mask_search::SearchMask::of(
                mask,
            ))
            .with_hard_starts(mask.and_then(|m| m.translated_starts().cloned())))
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-904
/// Whether the rows of the hit list `list` read its subject (`rows`, by each HSP's engine
/// identities, the Seq-align's `num_ident`): for `SubjectFetch::Every` a list with HSPs.
pub(super) fn reads_subject(list: &KappaResultList, rows: SubjectFetch) -> bool {
    match rows {
        SubjectFetch::Every => !list.hsps.hsps.is_empty(),
        SubjectFetch::Never | SubjectFetch::WithIdentities => list
            .payloads
            .iter()
            .any(|payload| rows.reads(usize::try_from(payload.num_ident).unwrap_or(0))),
    }
}

/// The query batches of a call before a batch's prefetch (`shared::LoadedBefore`, read by
/// `TblastnDatabase::loaded_subjects`): their final hit lists by searched query (`hits`, named by
/// rank), the searched queries of each batch (`spans`; those past `hits` are not of them),
/// which rows of the first format read a subject where a batch's prefetch is off (`rows`), and
/// the outfmt 5 reports cut before their hits (`cuts`).
pub(super) type LoadedBefore<'a> = shared::LoadedBefore<'a, KappaResultHitList>;

/// The headers of the database sequences that the final hit lists of the query batches searched
/// so far name, by ordinal id (`read_titles`): the stored title, or the definition lines of a
/// sequence with parsed Seq-ids.
pub(super) type DatabaseTitles = BTreeMap<u32, ReportedHeader>;

/// NCBI's fetches with `check_oid_exclusion` of the redo and the traceback of a search of
/// `database` (`shared::FetchExcluded`; `None` where they exclude nothing): the subjects they
/// exclude leave the preliminary hit lists (`LocalStageDProfile::fetch_excluded`).
pub(super) fn fetch_excluded(database: &TblastnDatabase) -> Option<shared::FetchExcluded<'_>> {
    shared::FetchExcluded::of(&database.db)
}

/// A TBLASTN call's reads of header records: NCBI's fetches with `check_oid_exclusion` of the
/// redo and the traceback of its batches (`excluded`, [`fetch_excluded`]), with the note of each
/// batch's first fetch that read a header with a UTF8String where a VisibleString is expected, by
/// the batch's first query (`notes`: the batches may be searched side by side), then the reads of
/// each batch's Seq-aligns and reported sequences (`titles`), whose first such read is the call's
/// (`utf8`, `shared::CallUtf8`: its warning, the model of the loads after it).
pub(super) struct CallReads<'a> {
    database: Option<&'a TblastnDatabase>,
    excluded: Option<shared::FetchExcluded<'a>>,
    notes: BTreeMap<usize, crate::seqdb::Utf8Reads>,
    utf8: shared::CallUtf8,
    /// NCBI's posts of the titles that each batch's Seq-aligns read, by the batch's first query
    /// (`shared::seq_align_title_posts`, `take_title_posts`).
    titled: std::sync::Mutex<Vec<(usize, crate::report::diag_log::Posts)>>,
    /// Those of the failed batch, by its first query, which its Seq-aligns read before the
    /// failure (`ErrorAfterBatchReports::title_posts`, `write_failed`).
    failed_titled: std::sync::Mutex<Option<(usize, crate::report::diag_log::Posts)>>,
}

impl<'a> CallReads<'a> {
    /// The reads of a call of `database` (`None`: the subjects are records) whose batches
    /// `batches` makes (each with its searched queries).
    pub(super) fn of(
        database: Option<&'a TblastnDatabase>,
        batches: impl FnOnce() -> Vec<crate::blastinput::seq_range::RangedBatch>,
    ) -> Self {
        let excluded = database.and_then(fetch_excluded);
        let notes = match excluded {
            Some(_) => batches()
                .into_iter()
                .map(|batch| (batch.searched.start, Default::default()))
                .collect(),
            None => Default::default(),
        };
        Self {
            database,
            excluded,
            notes,
            utf8: shared::CallUtf8::default(),
            titled: Default::default(),
            failed_titled: Default::default(),
        }
    }

    /// NCBI's posts of the titles that the batches' Seq-aligns read (`titles`), each with the
    /// batch's first query: they follow the batch's request posts and its warning of a UTF8String
    /// header, before its prefetch (`query_warnings::query_posts`, before the prefetches' posts).
    pub(super) fn take_title_posts(&self) -> Vec<(usize, crate::report::diag_log::Posts)> {
        std::mem::take(
            &mut *self
                .titled
                .lock()
                .expect("the title posts are not poisoned"),
        )
    }

    /// The fetches of the batch of the queries `range` (`stage_d_pipeline::BatchFetches`).
    pub(super) fn of_batch(
        &self,
        range: &std::ops::Range<usize>,
    ) -> Option<super::stage_d_pipeline::BatchFetches<'_>> {
        self.excluded
            .as_ref()
            .map(|excluded| super::stage_d_pipeline::BatchFetches {
                excluded,
                utf8: self.notes.get(&range.start),
            })
    }

    /// The batch of the queries `range` failed after its fetches with `check_oid_exclusion`
    /// (a failed Smith-Waterman redo of the batch's last match, `stage_d_pipeline::note_redo_fetched`):
    /// the first of them that read a header with a UTF8String where a VisibleString is expected is
    /// the call's read of this batch, where none came before, so that its warning comes before the
    /// batch's error (`write_failed`).
    pub(super) fn fetched_before_failure(&self, range: &std::ops::Range<usize>) {
        let Some(database) = self.database else {
            return;
        };
        let Some(note) = self
            .notes
            .get(&range.start)
            .and_then(crate::seqdb::Utf8Reads::take)
        else {
            return;
        };
        let read = database
            .db
            .utf8_string_path_at(note.at)
            .map(|path| shared::Utf8Read {
                query: range.start,
                path,
            });
        self.utf8.note(range.start, read.as_ref());
    }

    /// `-db`: the definition lines of the final hit lists `hitlists` of the batch of the queries
    /// `range`, read where NCBI makes the batch's Seq-aligns ([`read_titles`], after the batch's
    /// fetches and Seq-aligns that read a UTF8String header, [`drop_unaligned`]), into `titles`;
    /// the first read of such a header is the call's where none came before (`utf8`).
    /// `Ok(Some(error))`: NCBI's reader error (a header file that NCBI cannot read), which fails
    /// this batch after the reports of the batches before; LOSAT's refusals are errors.
    /// `Ok(None)` without a database.
    pub(super) fn titles(
        &self,
        hitlists: &mut [KappaResultHitList],
        range: &std::ops::Range<usize>,
        titles: &mut DatabaseTitles,
    ) -> Result<Option<anyhow::Error>> {
        let Some(database) = self.database else {
            return Ok(None);
        };
        let fetched = self
            .notes
            .get(&range.start)
            .and_then(crate::seqdb::Utf8Reads::take)
            .map(|note| note.at);
        let aligned = drop_unaligned(database, hitlists);
        let earlier = shared::EarlierUtf8::of_batch(range.start, fetched, aligned);
        let read = read_titles(database, hitlists, range.start, earlier, titles);
        self.utf8.note(
            range.start,
            shared::Utf8Read::of_read(&read, Option::as_ref),
        );
        // Whether this batch posts the call's warning of a UTF8String header (a batch before may
        // have read it).
        let utf8 = self.utf8.of_batch(range.start);
        // NCBI's posts of the titles that the batch's Seq-aligns read, in `read_titles`' order.
        let aligned = hitlists
            .iter()
            .flat_map(KappaResultHitList::lists)
            .filter_map(|list| u32::try_from(list.oid).ok());
        match read {
            Ok(_) => {
                let posts = shared::batch_title_posts(&database.db, aligned, utf8, PROGRAM)?;
                if !posts.is_empty() {
                    self.titled
                        .lock()
                        .expect("the title posts are not poisoned")
                        .push((range.start, posts));
                }
                Ok(None)
            }
            Err(error) => {
                match crate::blastinput::app::ErrorAfterBatchReports::after_seq_aligns(
                    error, aligned,
                )
                .downcast::<crate::blastinput::app::ErrorAfterBatchReports>()
                {
                    Ok(failure) => {
                        let posts = failure.title_posts(&database.db, utf8, PROGRAM)?;
                        *self
                            .failed_titled
                            .lock()
                            .expect("the title posts are not poisoned") =
                            Some((range.start, posts));
                        Ok(Some(failure.error))
                    }
                    Err(error) => Err(error),
                }
            }
        }
    }

    /// `shared::CallUtf8::reach` before the model of the loads of the batch of the queries `span`.
    pub(super) fn reach(&self, span: &std::ops::Range<usize>) {
        if let Some(database) = self.database {
            self.utf8.reach(span, &database.db);
        }
    }

    /// NCBI's warning of the call's first read of a UTF8String header among its batches' own in
    /// the first query of the batch that read it, a reported one (of the `reported` queries), with
    /// the posts that follow the reader's messages of the batch (`query_warnings::query_posts`),
    /// where the run has not posted it (a load of the model of the loads may have,
    /// `shared::load_read_texts`).
    pub(super) fn header_read(
        &self,
        reported: usize,
    ) -> Option<(usize, crate::report::diag_log::Posts)> {
        self.utf8.post(
            self.database.map(|database| &database.db),
            PROGRAM_NAME,
            reported,
        )
    }

    /// Writes NCBI's warning of that read to `diagnostics` (and the archive's messages) where the
    /// failed batch of the queries `range` made it with its Seq-aligns, before the failure: after
    /// the batch's reader messages and request posts; then the posts of the titles that those
    /// Seq-aligns read before the failure (`titles`).
    pub(super) fn write_failed(
        &self,
        range: &std::ops::Range<usize>,
        diagnostics: &mut dyn std::io::Write,
        archive: Option<&crate::report::archive::ArchiveHandle>,
    ) -> Result<()> {
        let db = self.database.map(|database| &database.db);
        if let Some(posts) = self.utf8.failed_post(db, PROGRAM_NAME, range.start) {
            posts.write(diagnostics, archive)?;
        }
        if let Some((_, posts)) = self
            .failed_titled
            .lock()
            .expect("the title posts are not poisoned")
            .as_ref()
            .filter(|(first, _)| *first == range.start)
        {
            posts.write(diagnostics, archive)?;
        }
        Ok(())
    }
}

/// The post prefix of TBLASTN's warnings (`[tblastn]`).
const PROGRAM_NAME: &str = "tblastn";

/// The final hit lists `hitlists` of a query batch without the sequences that NCBI makes no
/// Seq-align of (`shared::SeqAligned`): after their cut to the hit list size, before
/// `read_titles`.
pub(super) fn drop_unaligned(
    database: &TblastnDatabase,
    hitlists: &mut [KappaResultHitList],
) -> Option<crate::seqdb::SeqLocation> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2650-2664
    // The redo, the traceback or (`-ungapped`) the ungapped re-evaluation fetched every subject
    // with `check_oid_exclusion` (`fetch_excluded`).
    let mut aligned = shared::SeqAligned::of(&database.db, true)?;
    for hitlist in hitlists {
        hitlist.with_lists_mut(|lists| aligned.retain(lists, |list| u32::try_from(list.oid).ok()));
    }
    aligned.first_utf8()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
/// Reads the definition lines of the sequences that the final hit lists `hitlists` of a
/// query batch name (their subject index is the ordinal id), in NCBI's order (query by query,
/// in hit-list order) into `titles` (`common::database::read_header`); `first_query` is the
/// index of the batch's first query. A reader error of NCBI's own is an `ErrorAfterBatchReports` (it
/// fails the batch); LOSAT's refusals are errors of their own. The first read of a record with a
/// UTF8String where a VisibleString is expected (`shared::HeaderReader::first_utf8`) is returned:
/// NCBI posts its warning there. `one_by_one`: the reports may read the reported sequences one by
/// one (`shared::HeaderReader::new`).
pub(super) fn read_titles(
    database: &TblastnDatabase,
    hitlists: &[KappaResultHitList],
    first_query: usize,
    earlier: shared::EarlierUtf8,
    titles: &mut DatabaseTitles,
) -> Result<Option<shared::Utf8Read>> {
    let headers = shared::HeaderReader::new(&database.db, database.taxdb.unrecorded_calls());
    headers.note_earlier(earlier);
    for (index, hitlist) in hitlists.iter().enumerate() {
        for list in hitlist.lists() {
            let ordinal = u32::try_from(list.oid).context("negative TBLASTN subject OID")?;
            if titles.contains_key(&ordinal) {
                continue;
            }
            let title = shared::read_header(&headers, ordinal, first_query + index, PROGRAM)?;
            titles.insert(ordinal, title);
        }
    }
    Ok(headers.first_utf8())
}

/// What the XML2 and JSON formats show of a database sequence beyond its title: its
/// Blast-def-line taxid and the taxonomy database's answer, by rank (empty for the other
/// formats and for `-subject`: every taxid is then 0), and NCBI's messages about the taxonomy
/// database, which are written at the first `Hit` that needs them.
#[derive(Default)]
pub(super) struct XmlTaxa {
    pub(super) taxids: Vec<i32>,
    pub(super) taxonomy: Vec<crate::report::blast_xml::HitTaxonomy>,
    /// What the taxonomy database says about the taxid of each definition line of a sequence
    /// with parsed Seq-ids, by rank (empty as `taxonomy`).
    pub(super) parsed_taxonomy: Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
    pub(super) notes: shared::TaxdbNotes,
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(super) struct ReportedSubjects {
    /// The ordinal id of each rank, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The record of each rank: the id `gnl|BL_ORD_ID|n` and the stored title (a sequence with
    /// parsed Seq-ids: the title of its first definition line), no letters.
    pub(super) records: Vec<FastaRecord>,
    /// The number of bases of each rank.
    pub(super) lengths: Vec<usize>,
    /// The packed subject of each rank (its letters for the reports); empty when no report reads
    /// them.
    pub(super) stores: Vec<PackedNucleotide>,
    /// The sequences with parsed Seq-ids by rank (empty when none is reported): their ids and
    /// titles come from their definition lines (`report::db_subject`).
    pub(super) parsed: Vec<Option<Arc<ParsedSubject>>>,
    /// The sources from which NCBI makes the empty stored titles of the others, by rank
    /// (`shared::StoredTitle::source`; empty when none has one).
    pub(super) sources: Vec<crate::report::defline::SourceTaxname>,
    /// What the XML2 and JSON reports do at each rank's `Hit` (`HitLoads`, the one answer, by
    /// rank): the refusals, the names (`blast_xml::prepare`), the taxonomy reads (`read_taxids`)
    /// and the stop's record of the loads (`prefetch_stop_error`) take it.
    pub(super) loads: HitLoads,
    /// Whether the batch of each query prefetches, and whether its prefetch fails
    /// (`shared::BatchPrefetch::of`, the one producer, by the call's query index): the batch loop's
    /// reads of a batch whose reports load one by one (`BatchPrefetch::loads_one_by_one`) take it.
    pub(super) prefetch: shared::BatchPrefetch,
    /// Whether the search is gapped: the Seq-aligns that NCBI's prefetch counts are one for each
    /// HSP, else one for each subject (`seq_aligns_of`).
    pub(super) gapped: bool,
}

impl ReportedSubjects {
    /// The sequence of rank `rank` when it has parsed Seq-ids.
    pub(super) fn parsed_subject(&self, rank: usize) -> Option<&ParsedSubject> {
        self.parsed.get(rank).and_then(Option::as_deref)
    }

    /// The reported sequences as the library's records name them (`shared::library_subjects`).
    pub(super) fn library_subjects(
        &self,
        database: &TblastnDatabase,
    ) -> Vec<crate::api::local_blast::ReportedSubject> {
        shared::library_subjects(
            &database.db,
            &self.records,
            &self.ordinals,
            &self.parsed,
            &self.sources,
            |rank| self.lengths[rank],
            PROGRAM,
        )
    }

    /// The subjects as the reports read them, by rank: their stores, or without stores (no report
    /// reads their letters) their lengths.
    pub(super) fn letters(&self) -> super::x_subject_store::SubjectLetters<'_> {
        use super::x_subject_store::SubjectLetters;
        if self.stores.is_empty() {
            SubjectLetters::Lengths(&self.lengths)
        } else {
            SubjectLetters::Packed(&self.stores)
        }
    }
}

/// The ranks of the reported sequences of the final hit lists `hitlists` (their subject index is
/// the rank, `reported_subjects`): those with HSPs.
fn reported_ranks(hitlists: &[KappaResultHitList]) -> impl Iterator<Item = usize> + '_ {
    hitlists
        .iter()
        .flat_map(KappaResultHitList::lists)
        .filter(|list| !list.hsps.hsps.is_empty())
        .map(|list| list.oid as usize)
}

/// The first of the query batches `batches` (the searched queries of each) of the final hit lists
/// `hitlists` (subjects named by ordinal id, before `reported_subjects`) whose prefetch is on in a
/// run of the first format `format` (its Seq-aligns of a nonzero score do not outnumber a positive
/// `PRE_FETCH_SEQS_LIMIT`, `shared::prefetch_off`) and loads a subject whose lookup fails on a file with NCBI's exception
/// (`ParsedSubject::prefetch_stops`), with the subjects that batch reports (ordinal id, and the
/// definition lines of one with parsed Seq-ids, `titles`): NCBI's prefetch stops the run there
/// (`shared::batch_prefetch`), before that batch's reports. The caller asks it where the first
/// format prefetches.
pub(super) fn prefetch_stop(
    titles: &DatabaseTitles,
    hitlists: &[KappaResultHitList],
    batches: &[std::ops::Range<usize>],
    (format, gapped, has_pssm): (shared::PrefetchFormat, bool, bool),
    batch_scope: &dyn Fn(&std::ops::Range<usize>) -> crate::blastinput::query_batch::QueryScope,
) -> Option<(
    std::ops::Range<usize>,
    Vec<(u32, Option<Arc<ParsedSubject>>)>,
)> {
    // The report drops the HSPs of score 0 first (`reported_subjects`).
    let scored = |list: &KappaResultList| {
        list.hsps
            .hsps
            .iter()
            .filter(|linked| linked.hsp.score != 0)
            .count()
    };
    batches.iter().find_map(|span| {
        let lists = || {
            hitlists
                .get(span.clone())
                .into_iter()
                .flatten()
                .flat_map(KappaResultHitList::lists)
        };
        let seq_aligns = lists()
            .map(|list| shared::subject_seq_aligns(scored(list), gapped))
            .sum();
        if shared::prefetch_off(format, seq_aligns)
            || (has_pssm
                && hitlists
                    .get(span.clone())
                    .is_some_and(super::nullb_report::hitlists_have_empty_subject_row_db))
        {
            return None;
        }
        let subjects: Vec<(u32, Option<Arc<ParsedSubject>>)> = lists()
            .filter(|list| scored(list) > 0)
            .filter_map(|list| u32::try_from(list.oid).ok())
            .map(|ordinal| {
                let parsed = match titles.get(&ordinal) {
                    Some(ReportedHeader::Parsed(subject)) => Some(Arc::clone(subject)),
                    _ => None,
                };
                (ordinal, parsed)
            })
            .collect();
        let mut scope: Option<crate::blastinput::query_batch::QueryScope> = None;
        subjects
            .iter()
            .any(|(_, parsed)| {
                parsed.as_ref().is_some_and(|subject| {
                    shared::prefetch_stops_in_scope(subject, &mut || {
                        scope.get_or_insert_with(|| batch_scope(span)).clone()
                    })
                })
            })
            .then(|| (span.clone(), subjects))
    })
}

/// The Seq-aligns of the results `hitlists` that NCBI's prefetch counts against
/// `PRE_FETCH_SEQS_LIMIT` (`shared::prefetch_off`): one for each HSP of a `gapped` search, one
/// for each subject of an ungapped one (`-ungapped`, `shared::subject_seq_aligns`).
pub(super) fn seq_aligns_of(hitlists: &[KappaResultHitList], gapped: bool) -> usize {
    hitlists
        .iter()
        .flat_map(KappaResultHitList::lists)
        .map(|list| shared::subject_seq_aligns(list.hsps.hsps.len(), gapped))
        .sum()
}

/// What the reports `formats` show of the reported database sequences, for LOSAT's refusals of
/// those with parsed Seq-ids (`common::database::ParsedShown`): the tabular formats (6, 7, 10,
/// 20) whose fields read the subject (`report::tabular::reads_subject`), outfmt 0, the
/// structured documents (5, 13 to 16), the formats that join the alignments of following hits
/// whose Seq-align ids match (0, 5, 13 to 16), those that write the sequences' taxids (13 to 16,
/// the tabular reports of taxonomy fields), and whether the first format's batches prefetch the
/// subjects.
pub(super) fn parsed_shown(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> crate::algorithm::common::database::ParsedShown {
    let pairwise = formats
        .iter()
        .any(|format| (0..=4).contains(&super::stage_e_report::format_number(format.outfmt)));
    crate::algorithm::common::database::ParsedShown {
        stops_at_failed_lookup: formats.iter().any(|format| {
            matches!(
                super::stage_e_report::format_number(format.outfmt),
                6 | 7 | 10 | 20
            ) && crate::report::tabular::stops_at_lookup_failure_for(
                super::tabular::tabular_spec(format.outfmt).as_ref(),
                crate::blastinput::app::Program::Tblastn,
            )
        }),
        stops_at_lineless: formats.iter().any(|format| {
            matches!(
                super::stage_e_report::format_number(format.outfmt),
                6 | 7 | 10 | 20
            ) && crate::report::tabular::stops_at_lineless_accession(
                super::tabular::tabular_spec(format.outfmt).as_ref(),
            )
        }),
        tabular: formats.iter().any(|format| {
            matches!(
                super::stage_e_report::format_number(format.outfmt),
                6 | 7 | 10 | 18 | 20
            ) && crate::report::tabular::writes_rows(
                super::tabular::tabular_spec(format.outfmt).as_ref(),
            ) && crate::report::tabular::reads_subject(
                super::tabular::tabular_spec(format.outfmt).as_ref(),
                crate::blastinput::app::Program::Tblastn,
            )
        }),
        pairwise,
        structured: formats
            .iter()
            .any(|format| super::stage_e_report::is_structured_xml(format.outfmt)),
        joined: pairwise
            || formats.iter().any(|format| {
                super::stage_e_report::is_structured_xml(format.outfmt)
                    || super::tabular::tabular_spec(format.outfmt)
                        .is_some_and(|spec| spec.needs().joins_subjects())
            }),
        sam: false,
        taxonomy: formats.iter().any(|format| {
            super::stage_e_report::xml2_kind(format.outfmt).is_some()
                || super::tabular::tabular_spec(format.outfmt)
                    .is_some_and(|spec| spec.needs().taxonomy)
        }),
        words: formats
            .iter()
            .filter(|format| {
                matches!(
                    super::stage_e_report::format_number(format.outfmt),
                    6 | 7 | 10 | 20
                )
            })
            .fold(
                crate::report::tabular::PrintedIds::default(),
                |words, format| {
                    words.or(crate::report::tabular::PrintedIds::of(
                        super::tabular::tabular_spec(format.outfmt).as_ref(),
                    ))
                },
            ),
        // NCBI's `s_IsPrefetchFormat` prefetches for every tabular format (cited in
        // `tblastn_app::search`).
        prefetched: shared::prefetches_subjects(Some(prefetch_format(formats)))
            && !shared::prefetch_always_off(prefetch_format(formats)),
        titles: shared::TitleReads::of_outfmts(
            crate::report::archive::options::RequestKind::Tblastn,
            formats
                .iter()
                .map(|format| super::stage_e_report::format_number(format.outfmt)),
        ),
    }
}

/// The first of the formats `formats` of a run as NCBI's prefetch of each query batch asks it
/// (`shared::PrefetchFormat`, the one owner of whether a batch prefetches): outfmt 10, 18 and 20
/// as the tabular 6, which NCBI prefetches for alike (`s_IsPrefetchFormat`); a run writes at
/// least one format, outfmt 0 stands for none.
pub(super) fn prefetch_format(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> shared::PrefetchFormat {
    let number =
        formats
            .first()
            .map_or(crate::blastinput::app::output_format::PAIRWISE, |format| {
                match super::stage_e_report::format_number(format.outfmt) {
                    10 | 18 | 20 => 6,
                    number => number,
                }
            });
    shared::PrefetchFormat::of(
        number,
        crate::report::archive::options::RequestKind::Tblastn,
    )
}

/// What the formats `formats` do after a subject whose lookup throws where its batch's prefetch is
/// off (`FormatThrow::of_outfmts`, for `reported_subjects`'s loads), with the call's refusal over
/// the `-outfmt` values of its other searches, `others` (`FormatThrow::of_search`,
/// `ReportOutputs::other_outfmts`); `html`: the run has `-html` links (`ResolvedTblastnArgs::html`,
/// which the writer's `PairwiseConfig::html` gets too).
pub(super) fn after_throw(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
    (others, html): (&[&str], bool),
) -> crate::report::db_subject::FormatThrow {
    crate::report::db_subject::FormatThrow::of_outfmts(
        formats.iter().map(|format| format.outfmt),
        crate::blastinput::app::Program::Tblastn,
        crate::report::archive::options::RequestKind::Tblastn,
        html,
    )
    .of_search(
        others,
        crate::blastinput::app::Program::Tblastn,
        crate::report::archive::options::RequestKind::Tblastn,
        html,
    )
}

/// What one query's pairwise report (outfmt 0 to 4) of the hit lists `lists` reads, by the key
/// that `lists` names the subjects by, a rank or an ordinal id (`QueryGroup::read_order` of the
/// rows' hits as the writer reads them, `reading`: each list's HSPs in the reports' order,
/// `stage_d_results::evalue_order`, with the values of their rows, `stage_e_report::hsp_scores`,
/// of the subject of `length(key)` nucleotides; `joins`: a copy joins the hit before).
fn pairwise_read_order<'l>(
    lists: impl IntoIterator<Item = &'l KappaResultList>,
    length: impl Fn(usize) -> Result<usize>,
    joins: impl Fn(u32, u32) -> bool,
    reading: crate::report::pairwise::PairwiseReading,
) -> Result<crate::report::db_subject::PairwiseReads> {
    let mut hits = Vec::new();
    for list in lists {
        let key = usize::try_from(list.oid).context("negative TBLASTN subject key")?;
        anyhow::ensure!(
            list.payloads.len() == list.hsps.hsps.len(),
            "TBLASTN report payload count mismatch"
        );
        let length = length(key)?;
        for index in super::stage_d_results::evalue_order(list) {
            hits.push((
                key as u32,
                super::stage_e_report::hsp_scores(
                    &list.hsps.hsps[index],
                    &list.payloads[index],
                    length,
                )?,
            ));
        }
    }
    Ok(crate::report::x_pairwise_par::QueryGroup::read_order(
        hits,
        |scores| *scores,
        joins,
        reading,
        false,
        None,
    ))
}

/// The reported database sequences of `hitlists` (the final hit lists of all query batches,
/// subjects named by ordinal id; their HSPs of score 0 are dropped first, as the report
/// drops them): the hit lists name them by rank from here on (the rank grows with the
/// ordinal id, so every comparison that breaks a tie by the subject index gives the same
/// order), and their records (`titles`) and lengths are kept for the reports, with their packed
/// subjects (`stores`, by ordinal id) when `letters` (`stage_e_report::reads_subject_letters`:
/// outfmt 0, the XML family and the tabular rows read them, the other reports the lengths alone).
/// A reported sequence whose stored title LOSAT does not report as
/// NCBI does is LOSAT's refusal, and so is one with parsed Seq-ids that the reports `formats` do
/// not write as NCBI (what they show, `parsed_shown`: `common::database::check_reported_parsed`,
/// and where they join hits `refuse_joined_hits`; where only the prefetch reads a subject, a batch
/// of `batches`, the query ranges of the batches of `hitlists`, whose prefetch is on reports it;
/// what the formats do after a `Hit` that throws, `after_throw`, with the call's refusal over the
/// `-outfmt` values of its other searches and whether it has `-html` links, `call`, for the one
/// answer of what the reports load, `hit_loads`, with how the pairwise reports read each query's
/// subjects, `pairwise`, after the loads of the calls before, whose loads nothing here records,
/// `TaxdbRun::calls_before`). The packed subjects are read from the database
/// (`SubjectLetters::Stream`, one read each) instead of cloned from the copy.
#[allow(clippy::too_many_arguments)]
pub(super) fn reported_subjects(
    database: &TblastnDatabase,
    titles: &DatabaseTitles,
    hitlists: &mut [KappaResultHitList],
    stores: super::x_subject_store::SubjectLetters<'_>,
    letters: bool,
    formats: &[crate::api::local_blast::FormatOutput<'_>],
    call: (&[&str], bool),
    pairwise: Option<crate::report::pairwise::PairwiseReading>,
    batches: &[std::ops::Range<usize>],
    gapped: bool,
    has_nullb: bool,
    seen: crate::report::query_bioseq::SeenFrom<'_>,
) -> Result<ReportedSubjects> {
    use super::x_subject_store::SubjectLetters;
    let shown = parsed_shown(formats);
    for hitlist in hitlists.iter_mut() {
        hitlist.drop_zero_score_hsps();
    }
    let mut ordinals = Vec::new();
    for list in hitlists.iter().flat_map(KappaResultHitList::lists) {
        ordinals.push(u32::try_from(list.oid).context("negative TBLASTN subject OID")?);
    }
    // What the batches read of each subject (a batch whose Seq-aligns outnumber a positive
    // `PRE_FETCH_SEQS_LIMIT` reads nothing of them before its reports, which read them one by
    // one; `shared::BatchPrefetch`, `shared::BatchReads::shown`).
    let format = prefetch_format(formats);
    let titled = |ordinal: usize| match titles.get(&u32::try_from(ordinal).ok()?) {
        Some(ReportedHeader::Parsed(subject)) => Some(subject.as_ref()),
        _ => None,
    };
    let prefetch = shared::BatchPrefetch::of(
        format,
        hitlists.len(),
        0,
        batches,
        |span| seq_aligns_of(&hitlists[span], gapped),
        |span| {
            (has_nullb
                && super::nullb_report::hitlists_have_empty_subject_row_db(&hitlists[span.clone()]))
                || shared::BatchPrefetch::prefetch_fails_among(
                    hitlists[span]
                        .iter()
                        .flat_map(KappaResultHitList::lists)
                        .map(|list| titled(usize::try_from(list.oid).ok()?)),
                )
        },
    );
    let reads = prefetch.reads(hitlists.iter().enumerate().flat_map(|(query, hitlist)| {
        hitlist
            .lists()
            .iter()
            .map(move |list| (query, list.oid as u32))
    }));
    // What the reports do at each subject's `Hit`, by ordinal id (the one answer, `hit_loads`).
    let after = after_throw(formats, call);
    // A pairwise report's reads, of the subjects' lengths (`pairwise_read_order`), where its
    // batch reads them one by one.
    let pairwise_keys = |query: usize| -> Result<Option<crate::report::db_subject::PairwiseReads>> {
        let Some(reading) = pairwise
            .filter(|_| after.run().is_some_and(|run| run.after.is_pairwise()))
            .filter(|_| prefetch.loads_one_by_one(query))
        else {
            return Ok(None);
        };
        let order = pairwise_read_order(
            hitlists[query]
                .lists()
                .iter()
                .filter(|list| !list.hsps.hsps.is_empty()),
            |ordinal| shared::x_nucleotide_length(&database.db, ordinal as u32, PROGRAM),
            |previous, ordinal| match (titled(previous as usize), titled(ordinal as usize)) {
                (Some(earlier), Some(later)) => {
                    crate::report::db_subject::joins_as_copy(earlier, later)
                }
                _ => false,
            },
            reading,
        )?;
        Ok(Some(order))
    };
    let pairwise_order = (0..hitlists.len())
        .map(pairwise_keys)
        .collect::<Result<Vec<_>>>()?;
    let loads = hit_loads(
        hitlists,
        &prefetch,
        (after, database.taxdb.calls_before()),
        pairwise_order,
        titled,
    );
    ordinals.sort_unstable();
    ordinals.dedup();
    let seen_all: std::collections::HashSet<u32> = if seen.any() {
        ordinals
            .iter()
            .copied()
            .filter(|&ordinal| {
                seen.hides_ordinal(
                    hitlists.iter().enumerate().filter_map(|(q_idx, hitlist)| {
                        hitlist
                            .lists()
                            .iter()
                            .any(|list| {
                                !list.hsps.hsps.is_empty()
                                    && u32::try_from(list.oid).ok() == Some(ordinal)
                            })
                            .then_some(q_idx)
                    }),
                    ordinal,
                )
            })
            .collect()
    } else {
        std::collections::HashSet::new()
    };
    for hitlist in hitlists.iter_mut() {
        for list in hitlist.lists_mut() {
            let rank = ordinals
                .binary_search(&(list.oid as u32))
                .expect("every reported subject was listed");
            list.oid = i32::try_from(rank).context("TBLASTN subject rank out of range")?;
        }
    }
    let mut records = Vec::with_capacity(ordinals.len());
    let mut lengths = Vec::with_capacity(ordinals.len());
    let mut kept = Vec::new();
    let mut parsed: Vec<Option<Arc<ParsedSubject>>> = Vec::new();
    let mut sources = Vec::new();
    // Whether `taxdb.bti` is on the search path, looked up at the first sequence with parsed
    // Seq-ids that outfmt 0 shows.
    let mut taxdb_found: Option<bool> = None;
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
                    &database.db,
                    ordinal,
                    stored.replaced_bytes,
                    shown.titles,
                    PROGRAM,
                )?;
                if shown.pairwise {
                    stored.check_shown_source(|| database.db.oid(ordinal), PROGRAM)?;
                }
                shared::keep_source(&mut sources, rank, ordinals.len(), stored.source.clone());
                stored.title.clone()
            }
            ReportedHeader::Parsed(subject) => {
                let taxdb =
                    shown.pairwise && *taxdb_found.get_or_insert_with(|| database.db.taxdb_found());
                shared::check_reported_parsed(
                    &database.db,
                    subject,
                    reads.shown(shown, ordinal),
                    loads.load(ordinal as usize),
                    false,
                    taxdb,
                    ordinal,
                    seen_all.contains(&ordinal),
                    PROGRAM,
                )?;
                if parsed.is_empty() {
                    parsed = vec![None; ordinals.len()];
                }
                parsed[rank] = Some(Arc::clone(subject));
                subject.first_title().to_vec()
            }
        };
        records.push(FastaRecord {
            // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the OID
            // (line DI-8).
            local_id: database.db.oid(ordinal).fasta_id(),
            title,
            sequence: Vec::new(),
            warnings: Vec::new(),
            source: Default::default(),
        });
        // The length of a sequence the reports read is its store's; of another, read from the
        // index and its last byte (`XDbStream::read_length`, NCBI's `GetSeqLength`).
        let stored = || {
            let SubjectLetters::Packed(stores) = stores else {
                unreachable!("a -db search holds its sequences")
            };
            stores.get(ordinal as usize).with_context(|| {
                format!(
                    "database sequence {} was not read",
                    database.db.oid(ordinal)
                )
            })
        };
        lengths.push(if letters {
            let store = match stores {
                SubjectLetters::Stream(stream) => stream.read(ordinal as usize)?,
                _ => stored()?.clone(),
            };
            let length = store.len();
            kept.push(store);
            length
        } else {
            match stores {
                SubjectLetters::Stream(stream) => stream.read_length(ordinal as usize)?,
                _ => stored()?.len(),
            }
        });
    }
    if shown.joined && !parsed.is_empty() {
        for (query, hitlist) in hitlists.iter().enumerate() {
            shared::refuse_joined_hits(
                query,
                hitlist
                    .lists()
                    .iter()
                    .filter(|list| !list.hsps.hsps.is_empty())
                    .map(|list| {
                        let rank = usize::try_from(list.oid).ok()?;
                        Some((ordinals[rank], parsed[rank].as_deref()?))
                    }),
                database.db.oid_map(),
                PROGRAM,
            )?;
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
        prefetch,
        gapped,
    })
}

/// What the XML2 and JSON reports of `hitlists` (the final hit lists by query, in the order of the
/// reports; the subjects named by `key`, an ordinal id or a rank) do at each subject's `Hit`
/// (`HitLoads::of`, the one answer): `prefetch` tells by query whether its batch loads one by one
/// (`BatchPrefetch::report_read`), `parsed` gives the subject with parsed Seq-ids of a key;
/// the `Hit`s are those that `blast_xml::prepare`'s reports make (the lists with HSPs), or a
/// pairwise report's reads (`pairwise_reads[query]`, where given). TBLASTN keeps no record of its
/// loads: the reports are the call's first loads, after `before`, those of the calls before
/// (`shared::TaxdbRun::calls_before`, `HitLoads::of`).
fn hit_loads<'p>(
    hitlists: &[KappaResultHitList],
    prefetch: &shared::BatchPrefetch,
    (after, before): (
        crate::report::db_subject::FormatThrow,
        crate::report::db_subject::LoadsBefore,
    ),
    pairwise_reads: Vec<Option<crate::report::db_subject::PairwiseReads>>,
    parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
) -> HitLoads {
    let mut pairwise_reads = pairwise_reads.into_iter();
    HitLoads::of(
        before,
        after,
        hitlists.iter().enumerate().map(|(query, hitlist)| {
            let keys = match pairwise_reads.next().flatten() {
                Some(reads) => crate::report::db_subject::ReportKeys::Pairwise(reads),
                None => crate::report::db_subject::ReportKeys::Entries(
                    hitlist
                        .lists()
                        .iter()
                        .map(|list| {
                            (!list.hsps.hsps.is_empty())
                                .then(|| usize::try_from(list.oid).ok())
                                .flatten()
                        })
                        .collect(),
                ),
            };
            (prefetch.report_read(query, 0), keys)
        }),
        parsed,
        // The most identities of a list's HSPs (`reads_subject`).
        |query, list| {
            hitlists[query].lists()[list]
                .payloads
                .iter()
                .map(|payload| usize::try_from(payload.num_ident).unwrap_or(0))
                .max()
                .unwrap_or(0)
        },
    )
}

impl ReportedSubjects {
    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every reported
    /// sequence, by rank (`algorithm::common::database::tabular_subject_ids`). A sequence whose
    /// id LOSAT cannot write as NCBI is an error before any format is written where a tabular
    /// report `reads` the subject, else its id is empty.
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
}

#[cfg(test)]
#[path = "database_tests.rs"]
mod tests;
