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
    db: SeqDb,
    name: String,
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

    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(super) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:188-194
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    //         _ASSERT(db_adapter && scope);
    // ```
    /// Opens the database of `-db name` (`algorithm::common::database::open_database`) with
    /// the command line's errors. The total length is NCBI's `Int8`; one that the engine
    /// cannot hold in a `usize` (a 32-bit target) is LOSAT's refusal.
    #[cfg(test)]
    pub(super) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        Self::open_restricted(name, settings, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line (`shared::restrict_database`, whose warnings go to `diagnostics`).
    pub(super) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        list: Option<&crate::seqdb::DbList>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        // A database with parsed Seq-ids is read (line DBID, ID2b).
        let parsed_seqids =
            shared::restrict_database(&mut db, name, list, None, diagnostics, settings, PROGRAM)?;
        // The statistics' setup total (`DbStatistics::setup_total`: not negative, below 2^63; a
        // `LENGTH` of 2^63 or more is NCBI's negative `Int8`, read as such) as the engine holds it.
        let total = db.statistics().setup_total;
        usize::try_from(total)
            .ok()
            .with_context(|| {
                format!(
                    "a database of {total} letters is not supported by LOSAT's TBLASTN on this platform"
                )
            })?;
        Ok(Self {
            db,
            name: name.to_string(),
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
    pub(super) fn name(&self) -> &str {
        &self.name
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:182-187
    // ```c++
    // static Int8
    // s_SeqDbGetTotLen(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetTotalLength();
    // }
    // ```
    /// The numbers of the statistics (`SeqDb::statistics`): the setup's total length and number
    /// of sequences, the Gumbel block's length and `CSeqDB::GetMinLength`; the real numbers, also
    /// with `-dbsize`.
    pub(super) fn numbers(&self) -> DatabaseNumbers {
        let statistics = self.db.statistics();
        DatabaseNumbers {
            // The setup's total is not negative (a length of OID 0 where the total is not
            // positive).
            total_length: statistics.setup_total as u64,
            num_seqs: statistics.setup_count,
            gumbel_db_length: statistics.gumbel_db_length,
            engine_db_length: statistics.engine_db_length,
            min_length: statistics.min_length,
            request_count: statistics.request_count,
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

    /// LOSAT's refusal of statistics of 0 sequences for a search without the Gumbel block
    /// (`shared::check_statistics_count`).
    pub(super) fn check_statistics_count(&self) -> Result<()> {
        shared::check_statistics_count(&self.db, PROGRAM)
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`).
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// The Blast-def-line taxid of each reported sequence by rank and what the taxonomy
    /// database says about it, for the XML2 and JSON formats
    /// (`algorithm::common::database::read_reported_taxonomy`, with the definition lines of the
    /// sequences with parsed Seq-ids); the sequences of `hitlists` with HSPs are the reported
    /// ones (their subject index is the rank, `reported_subjects`), and those that the reports
    /// load are read (`ReportedSubjects::loads`: not the subject of a `Hit` that throws, nor those
    /// after it in a query's report, which NCBI never loads, so that the taxonomy database's
    /// one-time texts go with the next `Hit` that it loads; review 1 of DBID pfoff-1, H1).
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

    /// The open of the taxonomy database by the rows of a query batch (`shared::TaxdbRun::open_for_rows`).
    pub(super) fn open_for_tabular_rows(&self) -> Result<()> {
        self.taxdb.open_for_rows(&self.db, PROGRAM)
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
            // A subject whose load throws loads nothing, nor does one that no report reaches
            // after it (`HitLoads::loaded`).
            .filter(|list| reported.loads.loaded(list.oid as usize))
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
        rows_open: bool,
        pairwise: Option<crate::report::pairwise::PairwiseReading>,
    ) -> Result<crate::report::query_warnings::QueryReads> {
        let lists = hitlists
            .iter()
            .flat_map(KappaResultHitList::lists)
            .filter(|list| reads_subject(list, rows));
        let ranks: Vec<u32> = match pairwise {
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
            .map(|rank| rank as u32)
            .collect(),
            None => lists.map(|list| list.oid as u32).collect(),
        };
        // A subject whose load throws loads nothing, nor does one that no report reaches after
        // it (`HitLoads::loaded`).
        let subjects = ranks
            .into_iter()
            .filter(|&rank| reported.loads.loaded(rank as usize))
            .map(|rank| {
                (
                    rank,
                    reported.ordinals[rank as usize],
                    reported.parsed_subject(rank as usize),
                )
            });
        shared::subject_read_texts(&self.db, &self.taxdb, subjects, rows_open, PROGRAM)
    }

    /// The stop of the prefetch of the batch of `queries` whose reported
    /// subjects are `subjects` (`prefetch_stop`): NCBI's error (a `shared::HeaderReadError`), or
    /// LOSAT's refusal (`shared::batch_prefetch`, after the batches before, `before`: the run's
    /// record of what NCBI's BLAST database data loader loaded for them, `LoadedBefore`).
    /// `at_hits`: the XML2 and JSON formats, which write the taxonomy texts at each `Hit`
    /// (`shared::held_batch_prefetch_stop`), after the lookups of the reports of the batches
    /// before (as `read_taxids` reads them, the reported subjects that those reports load,
    /// `ReportedSubjects::loads`, made of the batches before alone); TBLASTN writes those reports
    /// after the stop. The outfmt 5 reports before the batch are the queries'
    /// (`shared::PrefetchQueries::after_xml_reports`).
    pub(super) fn prefetch_stop_error(
        &self,
        reported: Option<&ReportedSubjects>,
        before: LoadedBefore<'_>,
        subjects: &[(u32, Option<Arc<ParsedSubject>>)],
        queries: &shared::PrefetchQueries<'_>,
        at_hits: bool,
    ) -> anyhow::Error {
        let subjects = subjects.iter().map(|(ordinal, parsed)| {
            shared::PrefetchSubject::unflipped(*ordinal, parsed.as_deref())
        });
        let stopped = if at_hits {
            let before = reported.into_iter().flat_map(|reported| {
                reported_ranks(before.hits)
                    .filter(|rank| reported.loads.loaded(*rank))
                    .map(|rank| {
                        shared::PrefetchSubject::unflipped(
                            reported.ordinals[rank],
                            reported.parsed_subject(rank),
                        )
                    })
            });
            shared::held_batch_prefetch_stop(
                &self.db,
                &self.taxdb,
                before,
                std::iter::empty(),
                subjects,
                PROGRAM,
                queries,
            )
            .err()
        } else {
            let before = reported
                .into_iter()
                .flat_map(|reported| Self::loaded_subjects(reported, before));
            shared::batch_prefetch(&self.db, &self.taxdb, before, subjects, PROGRAM, queries).err()
        };
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
        before.subjects(
            |hitlists, span| hitlists.get(span),
            |_, hitlists| seq_aligns_of(hitlists),
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

    /// The database read where it is searched by `threads` threads, after the checks of
    /// `common::database::x_check_nucleotide_stores` (`common/x_db_stream.rs`).
    pub(super) fn x_stream(&self, threads: usize) -> Result<XDbStream<'_>> {
        XDbStream::new(&self.db, threads, PROGRAM)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:900-904
// ```c++
//     if (x_IsFieldRequested(eQueryStart) || x_IsFieldRequested(eQueryEnd) ||
//     ...
//         x_IsFieldRequested(eSubjectSeq) || (x_IsFieldRequested(eNumIdentical) && num_ident > 0 ) ||
// ```
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

/// The final hit lists `hitlists` of a query batch without the sequences that NCBI makes no
/// Seq-align of (`shared::SeqAligned`): after their cut to the hit list size, before
/// `read_titles`.
pub(super) fn drop_unaligned(database: &TblastnDatabase, hitlists: &mut [KappaResultHitList]) {
    // tblastn is gapped (`-ungapped` is refused): the redo or the traceback fetched every
    // subject with `check_oid_exclusion` (`fetch_excluded`).
    let Some(mut aligned) = shared::SeqAligned::of(&database.db, true) else {
        return;
    };
    for hitlist in hitlists {
        hitlist.with_lists_mut(|lists| aligned.retain(lists, |list| u32::try_from(list.oid).ok()));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
// ```c++
//     m_TbackSearch.Reset(new CBlastTracebackSearch(m_QueryFactory,
//                                                   m_InternalData,
//                                                   m_Opts,
//                                                   seqinfo_src,
//                                                   search_msgs));
//     ...
//     CRef<CSearchResultSet> retval = m_TbackSearch->Run();
// ```
/// Reads the definition lines of the sequences that the final hit lists `hitlists` of a
/// query batch name (their subject index is the ordinal id), in NCBI's order (query by query,
/// in hit-list order) into `titles` (`common::database::read_header`); `first_query` is the
/// index of the batch's first query. A reader error of NCBI's own is a `HeaderReadError` (it
/// fails the batch); LOSAT's refusals are errors of their own.
pub(super) fn read_titles(
    database: &TblastnDatabase,
    hitlists: &[KappaResultHitList],
    first_query: usize,
    titles: &mut DatabaseTitles,
) -> Result<()> {
    let headers = shared::HeaderReader::new(&database.db);
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
    Ok(())
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
}

impl ReportedSubjects {
    /// The sequence of rank `rank` when it has parsed Seq-ids.
    pub(super) fn parsed_subject(&self, rank: usize) -> Option<&ParsedSubject> {
        self.parsed.get(rank).and_then(Option::as_deref)
    }

    /// Whether a row of `hitlists` (subjects by rank) loads its subject in the reports: only such
    /// a row looks a taxonomy name up, and a row whose load throws keeps the names of the rows
    /// before (`loads`; review 1 of line DBID pfoff-2, H1).
    pub(super) fn loads_a_row(
        &self,
        hitlists: &[super::stage_d_results::KappaResultHitList],
    ) -> bool {
        hitlists
            .iter()
            .flat_map(super::stage_d_results::KappaResultHitList::lists)
            .any(|list| !list.hsps.hsps.is_empty() && self.loads.loaded(list.oid as usize))
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
/// `hitlists` (subjects named by ordinal id, before `reported_subjects`) whose prefetch is on (its
/// Seq-aligns of a nonzero score do not outnumber a positive `PRE_FETCH_SEQS_LIMIT`,
/// `shared::prefetch_off`) and loads a subject whose lookup fails on a file with NCBI's exception
/// (`ParsedSubject::prefetch_stops`), with the subjects that batch reports (ordinal id, and the
/// definition lines of one with parsed Seq-ids, `titles`): NCBI's prefetch stops the run there
/// (`shared::batch_prefetch`), before that batch's reports. The caller asks it where the first
/// format prefetches.
pub(super) fn prefetch_stop(
    titles: &DatabaseTitles,
    hitlists: &[KappaResultHitList],
    batches: &[std::ops::Range<usize>],
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
        if shared::prefetch_off(lists().map(scored).sum()) {
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
        subjects
            .iter()
            .any(|(_, parsed)| {
                parsed
                    .as_ref()
                    .is_some_and(|subject| subject.prefetch_stops().is_some())
            })
            .then(|| (span.clone(), subjects))
    })
}

/// The Seq-aligns of the results `hitlists` that NCBI's prefetch counts against
/// `PRE_FETCH_SEQS_LIMIT` (`shared::prefetch_off`): one for each HSP (TBLASTN's search is gapped;
/// LOSAT refuses `-ungapped`, whose Seq-aligns are one per subject).
pub(super) fn seq_aligns_of(hitlists: &[KappaResultHitList]) -> usize {
    hitlists
        .iter()
        .flat_map(KappaResultHitList::lists)
        .map(|list| list.hsps.hsps.len())
        .sum()
}

/// What the reports `formats` show of the reported database sequences, for LOSAT's refusals of
/// those with parsed Seq-ids (`common::database::ParsedShown`): the tabular formats (6, 7, 10,
/// 20) whose fields read the subject (`report::tabular::reads_subject`), outfmt 0, the
/// structured documents (5, 13 to 16), the formats that join the alignments of following hits
/// whose Seq-align ids match (0, 5, 13 to 16), those that write the sequences' taxids (13 to 16,
/// the tabular reports of taxonomy fields), and whether the first format's batches prefetch the
/// subjects.
fn parsed_shown(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> crate::algorithm::common::database::ParsedShown {
    let numbers: Vec<i32> = formats
        .iter()
        .map(|format| super::stage_e_report::format_number(format.outfmt))
        .collect();
    let pairwise = numbers.iter().any(|number| (0..=4).contains(number));
    crate::algorithm::common::database::ParsedShown {
        stops_at_failed_lookup: formats.iter().any(|format| {
            matches!(
                super::stage_e_report::format_number(format.outfmt),
                6 | 7 | 10 | 20
            ) && crate::report::tabular::may_stop_at_failed_lookup(
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
        // NCBI's `s_IsPrefetchFormat` prefetches for every tabular format (quoted in
        // `tblastn_app::search`).
        prefetched: crate::algorithm::common::database::prefetches_subjects(numbers.first().map(
            |&number| match number {
                10 | 18 | 20 => 6,
                number => number,
            },
        )) && !crate::algorithm::common::database::prefetch_always_off(),
        one_by_one: crate::algorithm::common::database::prefetch_may_be_off(),
    }
}

/// What the formats `formats` do after a subject whose lookup throws where its batch's prefetch is
/// off (`FormatThrow::of_outfmts`, for `reported_subjects`'s loads).
fn after_throw(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> Option<crate::report::db_subject::RunThrow> {
    crate::report::db_subject::FormatThrow::of_outfmts(
        formats.iter().map(|format| format.outfmt),
        crate::blastinput::app::Program::Tblastn,
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
/// whether every format writes a `Hit` that throws as NCBI does, `writes_thrown`, for the one
/// answer of what the reports load, `hit_loads`, with how the pairwise reports read each query's
/// subjects, `pairwise`). The packed
/// subjects are read from the database (`SubjectLetters::Stream`, one read each) instead of
/// cloned from the copy.
#[allow(clippy::too_many_arguments)]
pub(super) fn reported_subjects(
    database: &TblastnDatabase,
    titles: &DatabaseTitles,
    hitlists: &mut [KappaResultHitList],
    stores: super::x_subject_store::SubjectLetters<'_>,
    letters: bool,
    formats: &[crate::api::local_blast::FormatOutput<'_>],
    pairwise: Option<crate::report::pairwise::PairwiseReading>,
    batches: &[std::ops::Range<usize>],
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
    let prefetch = shared::BatchPrefetch::of(hitlists.len(), 0, batches, |span| {
        seq_aligns_of(&hitlists[span])
    });
    let reads = prefetch.reads(hitlists.iter().enumerate().flat_map(|(query, hitlist)| {
        hitlist
            .lists()
            .iter()
            .map(move |list| (query, list.oid as u32))
    }));
    // What the reports do at each subject's `Hit`, by ordinal id (the one answer, `hit_loads`).
    let titled = |ordinal: usize| match titles.get(&u32::try_from(ordinal).ok()?) {
        Some(ReportedHeader::Parsed(subject)) => Some(subject.as_ref()),
        _ => None,
    };
    let after = after_throw(formats);
    // A pairwise report's reads, of the subjects' lengths (`pairwise_read_order`), where its
    // batch reads them one by one.
    let pairwise_keys = |query: usize| -> Result<Option<crate::report::db_subject::PairwiseReads>> {
        let Some(reading) = pairwise
            .filter(|_| after.is_some_and(|run| run.after.is_pairwise()))
            .filter(|_| prefetch.reads_one_by_one(query))
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
    let loads = hit_loads(hitlists, &prefetch, after, pairwise_order, titled);
    ordinals.sort_unstable();
    ordinals.dedup();
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
    })
}

/// What the XML2 and JSON reports of `hitlists` (the final hit lists by query, in the order of the
/// reports; the subjects named by `key`, an ordinal id or a rank) do at each subject's `Hit`
/// (`HitLoads::of`, the one answer): `prefetch` tells by query whether its batch reads one by one
/// (`BatchPrefetch::reads_one_by_one`), `parsed` gives the subject with parsed Seq-ids of a key;
/// the `Hit`s are those that `blast_xml::prepare`'s reports make (the lists with HSPs), or a
/// pairwise report's reads (`pairwise_reads[query]`, where given).
fn hit_loads<'p>(
    hitlists: &[KappaResultHitList],
    prefetch: &shared::BatchPrefetch,
    after: Option<crate::report::db_subject::RunThrow>,
    pairwise_reads: Vec<Option<crate::report::db_subject::PairwiseReads>>,
    parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
) -> HitLoads {
    let mut pairwise_reads = pairwise_reads.into_iter();
    HitLoads::of(
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
            (prefetch.reads_one_by_one(query), keys)
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
mod tests {
    use super::*;

    // DECISIONS 3 Q8 of line DBID: NCBI's choice of the Seq-ids of a database with parsed
    // Seq-ids depends on SEQ_ID_PREFER_ACCESSION_OVER_GI and SEQ_ID_AVOID_GI.
    #[test]
    fn a_seq_id_variable_refuses_a_database_with_parsed_seqids_only() {
        let path = |name: &str| {
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(name)
                .to_str()
                .expect("fixture path is UTF-8")
                .to_string()
        };
        let settings = ApplicationSettings {
            seq_id_rank_variable: Some("SEQ_ID_PREFER_ACCESSION_OVER_GI"),
            ..ApplicationSettings::default()
        };
        let parsed = path("tests/db_search/ids/nucl_v5");
        TblastnDatabase::open(&parsed, &ApplicationSettings::default())
            .expect("a database with parsed Seq-ids opens");
        let Err(error) = TblastnDatabase::open(&parsed, &settings) else {
            panic!("a database with parsed Seq-ids opened with SEQ_ID_PREFER_ACCESSION_OVER_GI")
        };
        assert_eq!(
            format!("{error:#}"),
            "the environment variable SEQ_ID_PREFER_ACCESSION_OVER_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's TBLASTN"
        );
        TblastnDatabase::open(&path("tests/seqdb/db/nucl_v4"), &settings)
            .expect("a database without parsed Seq-ids opens with SEQ_ID_PREFER_ACCESSION_OVER_GI");
    }

    /// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
    /// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
    fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> TblastnDatabase {
        use crate::seqdb::alias::SearchPath;
        let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
        let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
        let search =
            SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None).unwrap();
        let path = tests.join("db_search").join(name);
        TblastnDatabase {
            db: {
                let mut db =
                    SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
                db.read_name_records(path.to_str().unwrap(), false).unwrap();
                db
            },
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
            parsed_seqids: true,
        }
    }

    /// The final hit list of one query whose lists name the sequences `lists` by ordinal id, each
    /// with one HSP of the score given (0: no Seq-align, `KappaResultHitList::drop_zero_score_hsps`).
    fn hitlist(lists: &[(i32, i32)]) -> KappaResultHitList {
        use super::super::search_gapped::GappedHsp;
        use super::super::stage_d_linking::{LinkedHsp, LinkedHspList};
        use super::super::stage_d_results::KappaResultList;
        let mut hitlist = KappaResultHitList::new(500).unwrap();
        for &(oid, score) in lists {
            let hsp = LinkedHsp {
                context: 0,
                hsp: GappedHsp {
                    frame: 1,
                    score,
                    q_start: 0,
                    q_end: 1,
                    q_gapped_start: 0,
                    s_start: 0,
                    s_end: 1,
                    s_gapped_start: 0,
                },
                num: 1,
                evalue: 1e-9,
                source_index: 0,
            };
            let list = KappaResultList {
                oid,
                hsps: LinkedHspList {
                    hsps: vec![hsp],
                    best_evalue: 1e-9,
                },
                payloads: Vec::new(),
            };
            assert!(hitlist.update(list).unwrap().is_none());
        }
        hitlist
    }

    // Reports that read no letters get the lengths of the reported sequences and no store; the
    // others get both, by rank.
    #[test]
    fn the_reported_sequences_are_read_only_for_reports_of_letters() {
        use super::super::x_subject_store::SubjectLetters;
        let cwd = crate::utils::test_dir::temp_dir().join(format!(
            "losat-tblastn-db-reported-{}-{:?}",
            crate::utils::test_dir::process_id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&cwd).unwrap();
        let database = open_with_taxdb("x16_ntaxid", None, &cwd);
        let stream = database.x_stream(1).unwrap();
        let of = |letters: bool| {
            let mut hitlists = vec![hitlist(&[(1, 50), (0, 40)])];
            let mut titles = DatabaseTitles::new();
            read_titles(&database, &hitlists, 0, &mut titles).unwrap();
            reported_subjects(
                &database,
                &titles,
                &mut hitlists,
                SubjectLetters::Stream(&stream),
                letters,
                &[],
                None,
                &[],
            )
            .unwrap()
        };
        let whole = [
            stream.read_length(0).unwrap(),
            stream.read_length(1).unwrap(),
        ];
        let reported = of(false);
        assert_eq!(reported.lengths, whole);
        assert!(reported.stores.is_empty());
        let letters = reported.letters();
        assert!(matches!(letters, SubjectLetters::Lengths(_)));
        assert_eq!(letters.count(), 2);
        assert_eq!(letters.length(1), whole[1]);
        let reported = of(true);
        assert_eq!(reported.lengths, whole);
        assert_eq!(reported.stores.len(), 2);
        assert_eq!(reported.stores[0].len(), whole[0]);
        assert_eq!(reported.stores[1].len(), whole[1]);
        let letters = reported.letters();
        assert!(matches!(letters, SubjectLetters::Packed(_)));
        assert_eq!(letters.length(1), whole[1]);
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    // DBID item T: a batch's prefetch looks up the sequences that its reports show (an HSP with a
    // Seq-align), once per volume and taxid over the batches of a run, and opens the taxonomy
    // database at the run's first lookup; nothing without `taxdb.bti`. `x16_ntaxid_99990` (two
    // sequences of taxid 99990) with `offset_last_past_end` (the last entry, 99990, is past the
    // end of `taxdb.btd`), `x16_ntaxid` (taxid 9606) with `bad_magic`.
    #[test]
    fn the_prefetch_looks_up_the_reported_sequences_once_per_run() {
        use super::super::x_subject_store::SubjectLetters;
        let cwd = crate::utils::test_dir::temp_dir().join(format!(
            "losat-tblastn-db-prefetch-{}-{:?}",
            crate::utils::test_dir::process_id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&cwd).unwrap();
        // The reports of one batch of `queries` (the lists of each query), as the run makes them:
        // the titles, the reported sequences by rank, then the batch's prefetch.
        let text = |database: &TblastnDatabase, queries: &[&[(i32, i32)]]| -> String {
            let stream = database.x_stream(1).unwrap();
            let mut hitlists: Vec<KappaResultHitList> =
                queries.iter().map(|&lists| hitlist(lists)).collect();
            let mut titles = DatabaseTitles::new();
            read_titles(database, &hitlists, 0, &mut titles).unwrap();
            let reported = reported_subjects(
                database,
                &titles,
                &mut hitlists,
                SubjectLetters::Stream(&stream),
                true,
                &[],
                None,
                &[],
            )
            .unwrap();
            let messages = database.prefetch_messages(&reported, &hitlists, SubjectFetch::Every);
            String::from_utf8(messages.unwrap()).unwrap()
        };
        let database = open_with_taxdb("x16_ntaxid_99990", Some("offset_last_past_end"), &cwd);
        assert!(database.taxdb_found());
        // A sequence whose HSPs have no Seq-align is not fetched: the taxonomy database is not
        // opened.
        assert_eq!(text(&database, &[&[(0, 0)]]), "");
        let offset = "Error: [tblastn] Error: Offset error at end of taxdb file.\n";
        assert_eq!(text(&database, &[&[(0, 0)], &[(1, 50)]]), offset);
        // The same volume and taxid in a later batch.
        assert_eq!(text(&database, &[&[(0, 50), (1, 50)]]), "");
        let database = open_with_taxdb("x16_ntaxid", Some("bad_magic"), &cwd);
        assert_eq!(
            text(&database, &[&[(1, 50), (0, 40)]]),
            "Error: [tblastn] Error: Tax database file has wrong magic number.\n"
        );
        let database = open_with_taxdb("x16_ntaxid", None, &cwd);
        assert!(!database.taxdb_found());
        assert_eq!(text(&database, &[&[(0, 50), (1, 40)]]), "");
        std::fs::remove_dir_all(&cwd).unwrap();
    }
}
