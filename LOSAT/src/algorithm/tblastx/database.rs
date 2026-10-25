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
//!   (`reported_subjects`).
//!
//! The search, the hit order and the report writers are those of the `-subject` path; the
//! opening of the database, the reader's errors and the rules for the stored titles are
//! shared with BLASTP, BLASTN, TBLASTN and BLASTX (`algorithm::common::database`).

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::{Context, Result};

use crate::algorithm::common::database::{self as shared, StoredTitle};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
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
        // A database with parsed Seq-ids is refused (line DBID: ID2b opens it).
        shared::refuse_parsed_seqids(&db, PROGRAM)?;
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

    /// The Blast-def-line taxid of each reported sequence by rank (`ordinals`: the ordinal id of
    /// each rank, `ReportedSubjects::ordinals`; every rank is reported, a rank is made from a
    /// hit) and what the taxonomy database says about it, for the XML2 and JSON formats
    /// (`common::database::read_reported_taxids`: the taxid, `sciname`, the taxonomy database's
    /// refusals and the non-ASCII refusal come before any byte).
    pub(crate) fn read_taxids(&self, ordinals: &[u32]) -> Result<XmlTaxa> {
        let reported: std::collections::BTreeSet<usize> = (0..ordinals.len()).collect();
        let (taxids, taxonomy, notes) =
            shared::read_reported_taxids(&self.db, &self.taxdb, ordinals, &reported, PROGRAM)?;
        Ok(XmlTaxa {
            taxids,
            taxonomy,
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
    pub(crate) notes: shared::TaxdbNotes,
}

/// The stored titles of the database sequences that the reported hits of the query batches
/// searched so far name, by ordinal id (`read_titles`).
pub(crate) type DatabaseTitles = BTreeMap<u32, StoredTitle>;

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
/// `titles`. A reader error of NCBI's own is a `HeaderReadError` (it fails the batch); LOSAT's
/// refusals are errors of their own.
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
        let title = shared::read_title(
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
    /// The records of the ranks: `gnl|BL_ORD_ID|n` and the stored title, without letters (the
    /// reports read them from the stored codes, `report::SubjectLetters::Stored`).
    pub(crate) records: Vec<FastaRecord>,
}

/// The reported database sequences of `ordered` (the reported hits of all query batches in
/// report order, subjects named by ordinal id): the hits name them by rank from here on (the
/// rank grows with the ordinal id, so every comparison that breaks a tie by the subject index
/// gives the same order), and their records (`titles`, by ordinal id; each is a sequence of
/// `stores`) are kept for the reports. A reported sequence whose stored title LOSAT does not
/// report as NCBI does is LOSAT's refusal.
pub(crate) fn reported_subjects(
    titles: &DatabaseTitles,
    ordered: &mut [TblastxHsp],
    stores: &[PackedNucleotide],
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
    for &ordinal in &ordinals {
        let stored = titles
            .get(&ordinal)
            .with_context(|| format!("database sequence {ordinal} has no title read"))?;
        shared::check_reported_title(stored.replaced_bytes, ordinal, PROGRAM)?;
        anyhow::ensure!(
            (ordinal as usize) < stores.len(),
            "database sequence {ordinal} was not read"
        );
        records.push(FastaRecord {
            local_id: format!("gnl|BL_ORD_ID|{ordinal}"),
            title: stored.title.clone(),
            sequence: Vec::new(),
            warnings: Vec::new(),
            source: Default::default(),
        });
    }
    Ok(ReportedSubjects { ordinals, records })
}

impl ReportedSubjects {
    /// The subject ids of the tabular reports, by rank: the first word of the stored title
    /// (`common::database::tabular_subject_id`).
    pub(crate) fn tabular_subject_ids(&self) -> Result<Vec<Arc<[u8]>>> {
        self.records
            .iter()
            .zip(&self.ordinals)
            .map(|(record, &ordinal)| {
                shared::tabular_subject_id(&record.title, ordinal, false, PROGRAM).map(Arc::from)
            })
            .collect()
    }
}
