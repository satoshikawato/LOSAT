//! The database path of TBLASTN (`-db NAME`): the search of an NCBI BLAST nucleotide
//! database, as NCBI BLAST+ 2.17.0 `tblastn -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, stage DB-5).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the subjects are the database's sequences, read once, in ordinal order, as packed
//!   subjects (`TblastnDatabase::read_stores`): the preliminary stage translates the stored
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

use crate::algorithm::common::database::{self as shared, StoredTitle};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{SeqDb, SeqType};
use crate::sequence::packed_subject::PackedNucleotide;

use super::stage_d_pipeline::DatabaseNumbers;
use super::stage_d_results::KappaResultHitList;

/// The program that LOSAT's refusals name.
pub(super) const PROGRAM: &str = "TBLASTN";

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`TblastnArgs::database_name`: the `-db` value without white space at its ends).
pub(super) struct TblastnDatabase {
    db: SeqDb,
    name: String,
    /// The taxonomy database of the run, opened by the XML2 and JSON formats only.
    taxdb: shared::TaxdbRun,
}

impl TblastnDatabase {
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
    pub(super) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        let db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        let total = db.total_length();
        usize::try_from(total)
            .ok()
            .filter(|&letters| i64::try_from(letters).is_ok())
            .with_context(|| {
                format!(
                    "a database of {total} letters is not supported by LOSAT's TBLASTN on this platform"
                )
            })?;
        Ok(Self {
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
        })
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
    /// The number of sequences (`CSeqDB::GetNumSeqs`) and the total length
    /// (`CSeqDB::GetTotalLength`): the real numbers, also with `-dbsize`.
    pub(super) fn numbers(&self) -> DatabaseNumbers {
        DatabaseNumbers {
            total_length: self.db.total_length(),
            num_seqs: u64::from(self.db.num_seqs()),
        }
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`).
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db, &self.name)
    }

    /// The Blast-def-line taxid of each reported sequence by rank and what the taxonomy
    /// database says about it, for the XML2 and JSON formats
    /// (`algorithm::common::database::read_reported_taxids`); the sequences of `hitlists` with
    /// HSPs are the reported ones (their subject index is the rank, `reported_subjects`).
    pub(super) fn read_taxids(
        &self,
        reported: &ReportedSubjects,
        hitlists: &[KappaResultHitList],
    ) -> Result<XmlTaxa> {
        let ranks: std::collections::BTreeSet<usize> = hitlists
            .iter()
            .flat_map(KappaResultHitList::lists)
            .filter(|list| !list.hsps.hsps.is_empty())
            .map(|list| list.oid as usize)
            .collect();
        let (taxids, taxonomy, notes) = shared::read_reported_taxids(
            &self.db,
            &self.taxdb,
            &reported.ordinals,
            &ranks,
            PROGRAM,
        )?;
        Ok(XmlTaxa {
            taxids,
            taxonomy,
            notes,
        })
    }

    /// Every sequence of the database by ordinal id, packed
    /// (`algorithm::common::database::read_nucleotide_stores`).
    pub(super) fn read_stores(&self) -> Result<Vec<PackedNucleotide>> {
        shared::read_nucleotide_stores(&self.db, PROGRAM)
    }

    /// EXPERIMENT (LOSAT_X_TBNDBSTREAM): the database read where it is searched, after the checks
    /// of `read_stores` (`x_db_stream.rs`).
    pub(super) fn x_stream(&self) -> Result<super::x_db_stream::XDbStream<'_>> {
        super::x_db_stream::XDbStream::new(&self.db, || self.read_stores())
    }
}

/// The stored titles of the database sequences that the final hit lists of the query
/// batches searched so far name, by ordinal id (`read_titles`).
pub(super) type DatabaseTitles = BTreeMap<u32, StoredTitle>;

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
/// in hit-list order) into `titles`; `first_query` is the index of the batch's first query.
/// A reader error of NCBI's own is a `HeaderReadError` (it fails the batch); LOSAT's
/// refusals are errors of their own.
pub(super) fn read_titles(
    database: &TblastnDatabase,
    hitlists: &[KappaResultHitList],
    first_query: usize,
    titles: &mut DatabaseTitles,
) -> Result<()> {
    for (index, hitlist) in hitlists.iter().enumerate() {
        for list in hitlist.lists() {
            let ordinal = u32::try_from(list.oid).context("negative TBLASTN subject OID")?;
            if titles.contains_key(&ordinal) {
                continue;
            }
            let title = shared::read_title(&database.db, ordinal, first_query + index, PROGRAM)?;
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
    pub(super) notes: shared::TaxdbNotes,
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(super) struct ReportedSubjects {
    /// The ordinal id of each rank, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The record of each rank: the id `gnl|BL_ORD_ID|n` and the stored title, no letters.
    pub(super) records: Vec<FastaRecord>,
    /// The packed subject of each rank (its letters for the reports).
    pub(super) stores: Vec<PackedNucleotide>,
}

/// The reported database sequences of `hitlists` (the final hit lists of all query batches,
/// subjects named by ordinal id; their HSPs of score 0 are dropped first, as the report
/// drops them): the hit lists name them by rank from here on (the rank grows with the
/// ordinal id, so every comparison that breaks a tie by the subject index gives the same
/// order), and their records (`titles`) and packed subjects (`stores`, by ordinal id) are
/// kept for the reports. A reported sequence whose stored title LOSAT does not report as
/// NCBI does is LOSAT's refusal. EXPERIMENT (LOSAT_X_TBNDBSTREAM): the packed subjects are read
/// from the database (`SubjectLetters::Stream`, one read each) instead of cloned from the copy.
pub(super) fn reported_subjects(
    titles: &DatabaseTitles,
    hitlists: &mut [KappaResultHitList],
    stores: super::x_subject_store::SubjectLetters<'_>,
) -> Result<ReportedSubjects> {
    use super::x_subject_store::SubjectLetters;
    for hitlist in hitlists.iter_mut() {
        hitlist.drop_zero_score_hsps();
    }
    let mut ordinals = Vec::new();
    for list in hitlists.iter().flat_map(KappaResultHitList::lists) {
        ordinals.push(u32::try_from(list.oid).context("negative TBLASTN subject OID")?);
    }
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
    let mut kept = Vec::with_capacity(ordinals.len());
    for &ordinal in &ordinals {
        let stored = titles
            .get(&ordinal)
            .with_context(|| format!("database sequence {ordinal} has no title read"))?;
        shared::check_reported_title(stored.replaced_bytes, ordinal, PROGRAM)?;
        records.push(FastaRecord {
            local_id: format!("gnl|BL_ORD_ID|{ordinal}"),
            title: stored.title.clone(),
            sequence: Vec::new(),
            warnings: Vec::new(),
            source: Default::default(),
        });
        kept.push(match stores {
            SubjectLetters::Stream(stream) => stream.read(ordinal as usize)?,
            SubjectLetters::Packed(stores) => stores
                .get(ordinal as usize)
                .with_context(|| format!("database sequence {ordinal} was not read"))?
                .clone(),
            SubjectLetters::Records(_) => unreachable!("a -db search holds no records"),
        });
    }
    Ok(ReportedSubjects {
        ordinals,
        records,
        stores: kept,
    })
}

impl ReportedSubjects {
    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every reported
    /// sequence, by rank (`algorithm::common::database::tabular_subject_id`). A sequence whose
    /// id LOSAT cannot write as NCBI is an error before any format is written.
    pub(super) fn tabular_subject_ids(&self) -> Result<Vec<Arc<[u8]>>> {
        self.records
            .iter()
            .zip(&self.ordinals)
            .map(|(record, &ordinal)| {
                shared::tabular_subject_id(&record.title, ordinal, false, PROGRAM).map(Arc::from)
            })
            .collect()
    }
}
