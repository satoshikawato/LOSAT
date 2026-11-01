//! The database path of TBLASTX (`-db NAME`): the search of an NCBI BLAST nucleotide
//! database, as NCBI BLAST+ 2.17.0 `tblastx -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, stage DB-5).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the subjects are the database's sequences, read once, in ordinal order, as packed
//!   subjects (`TblastxDatabase::read_stores`): the preliminary stage translates the stored
//!   2-bit data (makeblastdb's bases for the ambiguous ones, NCBI's preliminary input) and
//!   the ungapped re-evaluation and the reports read the stored ambiguity codes
//!   (`subjects::Subjects::Stored`); the search's subject index is the ordinal id;
//! - the set-level numbers come from the index files (`TblastxDatabase::numbers`): the
//!   number of sequences and the total length (an `Int8`); TBLASTX has no Gumbel block, so
//!   the shortest length a database answers does not enter;
//! - the definition lines of the reported sequences are read after each query batch's search
//!   (`read_titles`), and the reports show the database (`TblastxDatabase::name`,
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
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::db_subject::ParsedSubject;
use crate::report::pairwise::PairwiseDatabase;
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
    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(crate) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    /// Opens the nucleotide database `name` as NCBI's `CSearchDatabase` opens it
    /// (`common::database::open_database`). A total length that LOSAT cannot hold as an index
    /// (above `usize`) or as NCBI's `Int8` is LOSAT's refusal.
    pub(crate) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        Self::open_restricted(name, settings, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line (`shared::restrict_database`, whose warnings go to `diagnostics`).
    pub(crate) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        list: Option<&crate::seqdb::DbList>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        // A database with parsed Seq-ids is read (line DBID, ID2b), unless a variable changes
        // NCBI's choice of its Seq-ids.
        let parsed = shared::check_seqids(&db, settings, PROGRAM)?;
        shared::restrict_database(&mut db, list, diagnostics, PROGRAM)?;
        let total = db.total_length();
        usize::try_from(total)
            .ok()
            .filter(|&letters| i64::try_from(letters).is_ok())
            .with_context(|| {
                format!(
                    "a database of {total} letters is not supported by LOSAT's TBLASTX on this platform"
                )
            })?;
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

    /// Whether the database has parsed Seq-ids and more than `max_sequences` sequences in its
    /// volumes (`SeqDb::printed_num_seqs`, also under a taxid list: the range walks of
    /// `reports_may_refuse` decode a list's excluded sequences of every range they read): the
    /// check before the search does not read its headers (`shared::PARSED_CHECK_SEQUENCES_MAX`).
    pub(crate) fn parsed_beyond(&self, max_sequences: u32) -> bool {
        self.parsed && self.db.printed_num_seqs() > max_sequences
    }

    /// Whether a report that names any sequence of the database could be refused by LOSAT
    /// (`common::database::reports_may_refuse`, with the reports `shown` and the batches'
    /// prefetch `prefetch`), read in `pool`.
    pub(crate) fn reports_may_refuse(
        &self,
        shown: ParsedShown,
        prefetch: bool,
        pool: &crate::utils::threading::SearchPool<'_>,
    ) -> Option<shared::ReportRefusal> {
        shared::reports_may_refuse(
            &self.db,
            shown,
            false,
            self.taxdb_found(),
            prefetch,
            pool,
            PROGRAM,
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:164-168,183-187
    // ```c++
    // s_SeqDbGetNumSeqs(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetNumSeqs();
    // }
    // ...
    // s_SeqDbGetTotLen(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetTotalLength();
    // }
    // ```
    /// The numbers of the database that the statistics read (the total fits an `Int8`,
    /// `open`), with `-dbsize`.
    pub(crate) fn numbers(&self, dbsize: Option<i64>) -> SubjectSetNumbers {
        let (total, count) = self.db.statistics_numbers();
        SubjectSetNumbers {
            total: i64::try_from(total).expect("checked when opened"),
            count: i64::from(count),
            dbsize,
            searchsp: None,
        }
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(crate) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
    }

    /// The taxids of the search's list (`SeqDb::list_taxids`).
    pub(crate) fn list_taxids(&self) -> &[i32] {
        self.db.list_taxids()
    }

    /// The number of sequences of the volumes (`SeqDb::printed_num_seqs`): none is the database
    /// whose first query batch fails, also with a list.
    pub(crate) fn printed_num_seqs(&self) -> u32 {
        self.db.printed_num_seqs()
    }

    /// The `Database:` lines of the pairwise report.
    pub(crate) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// The Blast-def-line taxid of each reported sequence by rank (`reported`: the ordinal id of
    /// each rank and the sequences with parsed Seq-ids; every rank is reported, a rank is made
    /// from a hit) and what the taxonomy database says about it, for the XML2 and JSON formats
    /// (`common::database::read_reported_taxonomy`: the taxid, `sciname`, the taxonomy
    /// database's refusals and the non-ASCII refusal come before any byte).
    pub(crate) fn read_taxids(&self, reported: &ReportedSubjects) -> Result<XmlTaxa> {
        let ranks: std::collections::BTreeSet<usize> = (0..reported.ordinals.len()).collect();
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
    // ```c++
    //             	CRef<CBlast_def_line_set> bdlRef =
    //             			CSeqDB::ExtractBlastDefline(subject_bh);
    //     ...
    //             	if(setSubjectTaxInfoAll) {
    //             		x_SetTaxInfoAll(subject_bh, bdlRef);
    //             	}
    // ```
    /// The taxonomy of the tabular fields (`staxid`, `staxids`) of every reported sequence, by
    /// rank (every rank is made from a hit, `reported_subjects`;
    /// `algorithm::common::tabular::tabular_taxonomy`). A header that does not decode is
    /// refused before the report.
    pub(crate) fn tabular_taxonomy(
        &self,
        reported: &ReportedSubjects,
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        crate::algorithm::common::tabular::tabular_taxonomy(
            &self.db,
            &reported.ordinals,
            &reported.parsed,
            0..reported.ordinals.len(),
            PROGRAM,
        )
    }

    /// The database's sequences by ordinal id (`common::database::read_nucleotide_stores`).
    pub(crate) fn read_stores(&self) -> Result<Vec<PackedNucleotide>> {
        shared::read_nucleotide_stores(&self.db, PROGRAM)
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
    for hsp in ordered {
        let ordinal = hsp.hit.s_idx;
        if titles.contains_key(&ordinal) {
            continue;
        }
        let title = shared::read_header(
            &database.db,
            ordinal,
            first_query + hsp.hit.q_idx as usize,
            PROGRAM,
        )?;
        titles.insert(ordinal, title);
    }
    Ok(())
}

/// NCBI's messages about the taxonomy database when the query batch whose reported hits are
/// `ordered` (subjects named by ordinal id, as `read_titles` reads them) prefetches the subjects
/// it reports (`shared::prefetch_messages`): the sequences the hits name, with the definition
/// lines of those with parsed Seq-ids from `titles` (`read_titles`). Empty, and nothing read,
/// without `taxdb.bti` on the search path. A taxonomy database that LOSAT cannot read as NCBI
/// does, and a header whose taxid does not decode, are refused.
pub(crate) fn prefetch_messages(
    database: &TblastxDatabase,
    ordered: &[TblastxHsp],
    titles: &DatabaseTitles,
) -> Result<Vec<u8>> {
    let subjects = ordered.iter().map(|hsp| {
        let ordinal = hsp.hit.s_idx;
        let parsed = match titles.get(&ordinal) {
            Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
            _ => None,
        };
        (ordinal, parsed)
    });
    shared::prefetch_messages(&database.db, &database.taxdb, subjects, PROGRAM)
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(crate) struct ReportedSubjects {
    /// The ordinal id of each rank.
    pub(crate) ordinals: Vec<u32>,
    /// The records of the ranks: `gnl|BL_ORD_ID|n` and the stored title (a sequence with parsed
    /// Seq-ids: the title of its first definition line), without letters (the reports read them
    /// from the stored codes, `report::SubjectLetters::Stored`).
    pub(crate) records: Vec<FastaRecord>,
    /// The sequences with parsed Seq-ids by rank (empty when none is reported): their ids and
    /// titles come from their definition lines (`report::db_subject`).
    pub(crate) parsed: Vec<Option<Arc<ParsedSubject>>>,
}

/// The reported database sequences of `ordered` (the reported hits of all query batches in
/// report order, subjects named by ordinal id): the hits name them by rank from here on (the
/// rank grows with the ordinal id, so every comparison that breaks a tie by the subject index
/// gives the same order), and their records (`titles`, by ordinal id; each is a sequence of
/// `stores`) are kept for the reports. A reported sequence whose stored title LOSAT does not
/// report as NCBI does is LOSAT's refusal, and so is one with parsed Seq-ids that the reports
/// `shown` do not write as NCBI (`common::database::check_reported_parsed`, and with
/// `shown.joined` `refuse_joined_hits`).
pub(crate) fn reported_subjects(
    database: &TblastxDatabase,
    titles: &DatabaseTitles,
    ordered: &mut [TblastxHsp],
    stores: &[PackedNucleotide],
    shown: ParsedShown,
) -> Result<ReportedSubjects> {
    let mut ordinals: Vec<u32> = ordered.iter().map(|hsp| hsp.hit.s_idx).collect();
    ordinals.sort_unstable();
    ordinals.dedup();
    for hsp in ordered.iter_mut() {
        let rank = ordinals
            .binary_search(&hsp.hit.s_idx)
            .expect("every reported subject was listed");
        hsp.hit.s_idx = u32::try_from(rank).context("TBLASTX subject rank out of range")?;
    }
    let mut records = Vec::with_capacity(ordinals.len());
    let mut parsed: Vec<Option<Arc<ParsedSubject>>> = Vec::new();
    for (rank, &ordinal) in ordinals.iter().enumerate() {
        let header = titles
            .get(&ordinal)
            .with_context(|| format!("database sequence {ordinal} has no title read"))?;
        let title = match header {
            ReportedHeader::Stored(stored) => {
                shared::check_reported_title(stored.replaced_bytes, ordinal, PROGRAM)?;
                stored.title.clone()
            }
            ReportedHeader::Parsed(subject) => {
                // Whether `taxdb.bti` is on the search path, looked up at the first sequence
                // with parsed Seq-ids that outfmt 0 shows (once for the run, `taxdb_found`).
                let taxdb = shown.pairwise && database.taxdb_found();
                shared::check_reported_parsed(subject, shown, false, taxdb, ordinal, PROGRAM)?;
                if parsed.is_empty() {
                    parsed = vec![None; ordinals.len()];
                }
                parsed[rank] = Some(Arc::clone(subject));
                subject.first_title().to_vec()
            }
        };
        anyhow::ensure!(
            (ordinal as usize) < stores.len(),
            "database sequence {ordinal} was not read"
        );
        records.push(FastaRecord {
            // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the OID
            // (line DI-8).
            local_id: format!("gnl|BL_ORD_ID|{}", database.db.oid(ordinal)),
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
                PROGRAM,
            )?;
            start = end;
        }
    }
    Ok(ReportedSubjects {
        ordinals,
        records,
        parsed,
    })
}

impl ReportedSubjects {
    /// The subject ids of the tabular reports, by rank: the first word of the stored title, or
    /// the id of a sequence with parsed Seq-ids (`common::database::tabular_subject_ids`).
    pub(crate) fn tabular_subject_ids(&self) -> Result<Vec<Arc<[u8]>>> {
        shared::tabular_subject_ids(&self.records, &self.ordinals, &self.parsed, false, PROGRAM)
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
        TblastxDatabase::open(&parsed, &ApplicationSettings::default())
            .expect("a database with parsed Seq-ids opens");
        let Err(error) = TblastxDatabase::open(&parsed, &settings) else {
            panic!("a database with parsed Seq-ids opened with SEQ_ID_PREFER_ACCESSION_OVER_GI")
        };
        assert_eq!(
            format!("{error:#}"),
            "the environment variable SEQ_ID_PREFER_ACCESSION_OVER_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's TBLASTX"
        );
        TblastxDatabase::open(&path("tests/seqdb/db/nucl_v4"), &settings)
            .expect("a database without parsed Seq-ids opens with SEQ_ID_PREFER_ACCESSION_OVER_GI");
    }

    /// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
    /// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
    fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> TblastxDatabase {
        use crate::seqdb::alias::SearchPath;
        let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
        let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
        let search =
            SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None).unwrap();
        let path = tests.join("db_search").join(name);
        let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
        TblastxDatabase {
            parsed: db.has_parsed_seqids().unwrap(),
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
        }
    }

    /// A reported hit of the batch's first query on the database sequence `ordinal`.
    fn hit_on(ordinal: u32) -> TblastxHsp {
        TblastxHsp {
            hit: crate::common::Hit {
                identity: 100.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 30,
                s_start: 1,
                s_end: 30,
                e_value: 1e-9,
                bit_score: 50.0,
                num_ident: 10,
                query_frame: 1,
                query_length: 0,
                q_idx: 0,
                s_idx: ordinal,
                raw_score: 100,
                sort_query_offset: 0,
                sort_query_end: 10,
                sort_subject_offset: 0,
                sort_subject_end: 10,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 10,
            },
            subject_frame: 1,
            num: 1,
        }
    }

    // DBID item T: a batch's prefetch looks up the sequences that its reported hits name, once
    // per volume and taxid over the batches of a run, and opens the taxonomy database at the
    // run's first lookup; nothing without `taxdb.bti`. `x16_ntaxid_99990` (two sequences of
    // taxid 99990) with `offset_last_past_end` (the last entry, 99990, is past the end of
    // `taxdb.btd`), `x16_ntaxid` (taxid 9606) with `bad_magic`.
    #[test]
    fn the_prefetch_looks_up_the_reported_sequences_once_per_run() {
        let cwd = std::env::temp_dir().join(format!(
            "losat-tblastx-db-prefetch-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&cwd).unwrap();
        let text = |database: &TblastxDatabase, ordinals: &[u32]| -> String {
            let ordered: Vec<TblastxHsp> =
                ordinals.iter().map(|&ordinal| hit_on(ordinal)).collect();
            let mut titles = DatabaseTitles::new();
            read_titles(database, &ordered, 0, &mut titles).unwrap();
            String::from_utf8(prefetch_messages(database, &ordered, &titles).unwrap()).unwrap()
        };
        let database = open_with_taxdb("x16_ntaxid_99990", Some("offset_last_past_end"), &cwd);
        assert!(database.taxdb_found());
        // A batch without reported hits fetches nothing: the taxonomy database is not opened.
        assert_eq!(text(&database, &[]), "");
        // Two hits on sequence 1 and one on sequence 0, all of taxid 99990: one lookup.
        let offset = "Error: [tblastx] Error: Offset error at end of taxdb file.\n";
        assert_eq!(text(&database, &[1, 1, 0]), offset);
        // The same volume and taxid in a later batch.
        assert_eq!(text(&database, &[0, 1]), "");
        let database = open_with_taxdb("x16_ntaxid", Some("bad_magic"), &cwd);
        assert_eq!(
            text(&database, &[1, 0]),
            "Error: [tblastx] Error: Tax database file has wrong magic number.\n"
        );
        let database = open_with_taxdb("x16_ntaxid", None, &cwd);
        assert!(!database.taxdb_found());
        assert_eq!(text(&database, &[0, 1]), "");
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    // Line DBID, ID4t fix 2b: the loop writes each query batch's reports when it takes the batch
    // only where no report of any sequence of the database can be refused (`report_stream`): a
    // sequence that NCBI's formatters do not find by its Seq-align id (nrefuse_v5) refuses the
    // tabular formats and outfmt 0, and two sequences whose Seq-align ids match (njoined, in two
    // volumes) refuse outfmt 0 only, read on one thread and on four.
    #[test]
    fn a_report_of_any_sequence_that_could_be_refused_is_found() {
        let open = |name: &str| {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/db_search/ids")
                .join(name);
            TblastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
                .unwrap_or_else(|error| panic!("{name}: {error:#}"))
        };
        let tabular = ParsedShown {
            tabular: true,
            ..ParsedShown::default()
        };
        let pairwise = ParsedShown {
            pairwise: true,
            joined: true,
            ..ParsedShown::default()
        };
        let answers = |pool: &crate::utils::threading::SearchPool<'_>| {
            ["nucl_v5", "nucl_v4", "nrefuse_v5", "njoined"].map(|name| {
                let database = open(name);
                (
                    database.reports_may_refuse(tabular, false, pool),
                    database.reports_may_refuse(pairwise, false, pool),
                )
            })
        };
        use shared::ReportRefusal::{Headers, Joined};
        let serial = answers(&crate::utils::threading::SearchPool::x_serial());
        assert_eq!(
            serial,
            [
                (None, None),
                (None, None),
                (Some(Headers), Some(Headers)),
                (None, Some(Joined)),
            ]
        );
        let parallel =
            crate::utils::threading::with_search_pool(4, "TBLASTX", |pool| Ok(answers(pool)))
                .unwrap();
        assert_eq!(serial, parallel);
    }

    // Line DBID, DECISIONS 19: a database with parsed Seq-ids is beyond a bound of fewer
    // sequences than its volumes hold (`report_stream::Off::LargeParsed`), also where a taxid list
    // keeps fewer (the range walks decode the excluded ones of every range they read); a database
    // without parsed Seq-ids is never beyond it (read whatever its size, as blastn's).
    #[test]
    fn a_database_with_parsed_seqids_is_beyond_a_bound_below_its_volumes_sequences() {
        let open = |name: &str| {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join(name);
            TblastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
                .unwrap_or_else(|error| panic!("{name}: {error:#}"))
        };
        for name in ["db_search/ids/nucl_v5", "db_search/ids/nucl_v4"] {
            let database = open(name);
            let count = database.db.num_seqs();
            assert!(count > 1, "{name}");
            assert!(!database.parsed_beyond(count), "{name}");
            assert!(database.parsed_beyond(count - 1), "{name}");
        }
        let mut restricted = open("db_search/ids/nucl_v5");
        let count = restricted.db.num_seqs();
        restricted
            .db
            .restrict_taxids(&crate::seqdb::TaxidList::Positive(
                [9606].into_iter().collect(),
            ))
            .unwrap();
        let included = restricted.db.num_seqs();
        assert!(included > 0 && included < count);
        assert!(!restricted.parsed_beyond(count));
        assert!(restricted.parsed_beyond(included));
        let stored = open("seqdb/db/nucl_v4");
        assert!(stored.db.num_seqs() > 0);
        assert!(!stored.parsed_beyond(0));
    }
}
