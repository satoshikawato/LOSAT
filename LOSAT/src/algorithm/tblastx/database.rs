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
    /// The taxonomy database of the run, opened by the XML2 and JSON formats only.
    taxdb: shared::TaxdbRun,
}

impl TblastxDatabase {
    /// Opens the nucleotide database `name` as NCBI's `CSearchDatabase` opens it
    /// (`common::database::open_database`). A total length that LOSAT cannot hold as an index
    /// (above `usize`) or as NCBI's `Int8` is LOSAT's refusal.
    pub(crate) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        let db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        // A database with parsed Seq-ids is read (line DBID, ID2b), unless a variable changes
        // NCBI's choice of its Seq-ids.
        shared::check_seqids(&db, settings, PROGRAM)?;
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
            taxdb: shared::TaxdbRun::default(),
        })
    }

    /// The `-db` value (outfmt 7's `# Database:` line).
    pub(crate) fn name(&self) -> &str {
        &self.name
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
        SubjectSetNumbers {
            total: i64::try_from(self.db.total_length()).expect("checked when opened"),
            count: i64::from(self.db.num_seqs()),
            dbsize,
        }
    }

    /// The `Database:` lines of the pairwise report.
    pub(crate) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db, &self.name)
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
    // Whether `taxdb.bti` is on the search path, looked up at the first sequence with parsed
    // Seq-ids that outfmt 0 shows.
    let mut taxdb_found: Option<bool> = None;
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
                let taxdb =
                    shown.pairwise && *taxdb_found.get_or_insert_with(|| database.db.taxdb_found());
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
            local_id: format!("gnl|BL_ORD_ID|{ordinal}"),
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
}
