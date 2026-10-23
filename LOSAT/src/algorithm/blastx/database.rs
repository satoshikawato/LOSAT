//! The database path of BLASTX (`-db NAME`): the search of an NCBI BLAST protein database,
//! as NCBI BLAST+ 2.17.0 `blastx -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, P1 and P1.1).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the subjects are the stored bytes of the whole database, read once before the query
//!   batches (`BlastxDatabase::read_stores`, `common::database::read_protein_stores`) and
//!   searched by ordinal id through `subjects::Subjects::Database`; the set-level numbers
//!   are the index's (total, count) and 10 for the shortest sequence (D1);
//! - a database that cannot be read fails the first query batch with a valid context, after
//!   the reports of the batches before it (NCBI does not search a batch without one);
//! - after each batch's search, the titles of its reported sequences are read in NCBI's
//!   order and the batch's hit lists are renamed to ranks over records of those sequences
//!   (`ReportedSubjects::for_batch`), so the report writers keep indexing by subject;
//! - the reports show the database (`BlastxDatabase::pairwise_database`, `name`), the subject
//!   ids made from the stored titles, and `stitle`.
//!
//! The search, the composition redo and the report writers are those of the `-subject`
//! path. The opening of the database, the reader's errors and the rules for the stored
//! titles are shared with BLASTP, BLASTN and TBLASTN (`algorithm::common::database`).

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::Result;

use super::input::FastaRecord;
use super::runtime::BatchResults;
use super::subjects::Subjects;
use crate::algorithm::common::database::{self as shared, PreliminaryReadError, StoredProteins};
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::blast_xml::HitTaxonomy;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{SeqDb, SeqType};

/// The program that LOSAT's refusals name.
const PROGRAM: &str = "BLASTX";

/// The longest subject BLASTX searches in one piece (NCBI's `MAX_DBSEQ_LEN`).
const MAX_DBSEQ_LEN: usize = 5_000_000;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-119
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
// 'O', 'J'};
// ```
/// The letter of each NCBIstdaa byte.
const NCBISTDAA_TO_AMINOACID: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`BlastxArgs::database_name`: the `-db` value without white space at its ends).
pub(crate) struct BlastxDatabase {
    db: SeqDb,
    name: String,
    /// The taxonomy database of the run and the taxids its volumes have looked up (XML2 and
    /// JSON only), kept across the query batches of the run.
    taxdb: shared::TaxdbRun,
}

impl BlastxDatabase {
    /// Opens the database of `-db name` (`algorithm::common::database::open_database`)
    /// with the command line's errors.
    pub(crate) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        let db = shared::open_database(name, SeqType::Protein, settings, PROGRAM)?;
        Ok(Self {
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
        })
    }

    /// The `-db` value as NCBI keeps it.
    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    /// Whether `taxdb.bti` resolves on NCBI's search path, so that the XML2 and JSON formats
    /// write the scientific name of a taxid (`SeqDb::taxdb_found`).
    pub(crate) fn taxdb_found(&self) -> bool {
        self.db.taxdb_found()
    }

    /// Whether a later query batch of an XML2 or JSON run could still refuse a taxonomy
    /// database name or record (`algorithm::common::database::taxdb_may_refuse`: every
    /// sequence's taxid looked up once, before the first byte).
    pub(crate) fn taxdb_may_refuse(&self) -> bool {
        shared::taxdb_may_refuse(&self.db)
    }

    /// What the `Database:` lines of outfmt 0 show (title, date, real totals).
    pub(crate) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db, &self.name)
    }

    // A database without sequences (makeblastdb builds none; a volume whose index holds no
    // sequence): NCBI opens it, writes the outfmt 0 prolog, reads the first query batch and
    // fails in that batch's preliminary search, where its reader cannot map the sequence
    // file (BLASTP's `-db`, fixture `LOSAT/tests/seqdb/db/empty`).
    /// Reads every sequence of the database (`read_protein_stores`). The error (a reader
    /// error of NCBI's preliminary search, LOSAT's refusal of a sequence, or the empty
    /// database's) belongs to the first query batch with a valid context, which the caller
    /// fails with it.
    pub(crate) fn read_stores(&self) -> Result<StoredProteins> {
        if self.db.num_seqs() == 0 {
            return Err(PreliminaryReadError("Database memory map file error").into());
        }
        let stores = shared::read_protein_stores(&self.db, PROGRAM)?;
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
        // ```c
        // /** Split subject sequences if longer than this */
        // #define MAX_DBSEQ_LEN 5000000
        // ```
        // NCBI searches a longer subject in overlapping parts, which BLASTX does not port (the
        // same limit stops `-subject` in `search::search_core`); the refusal belongs to the
        // first batch with a valid context, as the reader's errors.
        if let Some(ordinal) = (0..stores.len()).find(|&i| stores.length(i) > MAX_DBSEQ_LEN) {
            return Err(too_long(ordinal, stores.length(ordinal)));
        }
        Ok(stores)
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): `read_stores` without the store
    /// (`common::database::x_check_protein_stores`): the same reads and checks with the same
    /// errors, the residues of the longest sequence (`StoredProteins::max_length`) when they pass.
    pub(crate) fn x_check_stores(&self) -> Result<usize> {
        if self.db.num_seqs() == 0 {
            return Err(PreliminaryReadError("Database memory map file error").into());
        }
        let mut longest = 0;
        let mut first_too_long = None;
        shared::x_check_protein_stores(&self.db, PROGRAM, |ordinal, length| {
            longest = longest.max(length);
            if length > MAX_DBSEQ_LEN && first_too_long.is_none() {
                first_too_long = Some((ordinal, length));
            }
        })?;
        if let Some((ordinal, length)) = first_too_long {
            return Err(too_long(ordinal as usize, length));
        }
        Ok(longest)
    }

    /// EXPERIMENT (LOSAT_X_BXDBSTREAM): the opened database.
    pub(crate) fn x_seqdb(&self) -> &SeqDb {
        &self.db
    }
}

/// LOSAT's refusal of database sequence `ordinal` of `length` residues, longer than
/// `MAX_DBSEQ_LEN`.
fn too_long(ordinal: usize, length: usize) -> anyhow::Error {
    anyhow::anyhow!(
        "database sequence {ordinal} has {length} residues; a subject longer than {MAX_DBSEQ_LEN} residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX"
    )
}

/// The reported sequences of one query batch, by rank (their order of ordinal ids): the
/// records the report writers index (`gnl|BL_ORD_ID|n`, the stored title, the letters of
/// the stored bytes), the stored titles, and the tabular ids when a tabular format is
/// written.
#[derive(Default)]
pub(crate) struct ReportedSubjects {
    records: Vec<FastaRecord>,
    titles: Vec<Arc<[u8]>>,
    tabular_ids: Vec<String>,
    /// The ordinal id of each rank, ascending (the Seq-ids of the Seq-annots,
    /// `general { db "BL_ORD_ID", tag id <oid> }`).
    ordinals: Vec<u32>,
    /// By rank: the Blast-def-line taxid and what the taxonomy database says of it, and the
    /// taxonomy database's messages of the batch (XML2 and JSON with `-db` only; empty
    /// otherwise, which reads as no taxid).
    taxids: Vec<i32>,
    taxonomy: Vec<HitTaxonomy>,
    taxdb_notes: shared::TaxdbNotes,
}

impl ReportedSubjects {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:286-291
    // ```c++
    //                 BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		                            fmt_args->GetFormattedOutputChoice());
    //             	ITERATE(CSearchResultSet, result, *results) {
    //                	    formatter.PrintOneResultSet(**result, query_batch);
    //             	}
    // ```
    /// The reported sequences of the batches `batches` (one query batch: the results of
    /// one search), with their hit lists renamed from ordinal ids to ranks. The titles are
    /// read before any report of the batch (NCBI fetches the batch's sequences first),
    /// query by query in hit-list order; NCBI's reader error stops the run at this batch
    /// (`HeaderReadError`), and a title LOSAT does not show as NCBI is its refusal
    /// (`check_reported_title`; with `tabular`, `tabular_subject_id`; with `pairwise`,
    /// `check_shown_protein_title`).
    pub(crate) fn for_batch(
        database: &BlastxDatabase,
        subjects: Subjects<'_>,
        batches: &mut [BatchResults],
        gapped: bool,
        tabular: bool,
        pairwise: bool,
        taxids: bool,
    ) -> Result<Self> {
        // Every listed subject gets a record (the tabular writer looks its list's subject
        // up before it filters the HSPs); only the reported ones have their titles read.
        let mut listed: BTreeMap<usize, Option<Arc<[u8]>>> = BTreeMap::new();
        for batch in batches.iter() {
            for (qidx, lists) in batch.queries.iter().enumerate() {
                for list in lists {
                    let reported = super::report::reported(list, gapped);
                    if matches!(listed.get(&list.oid), Some(Some(_)))
                        || (!reported && listed.contains_key(&list.oid))
                    {
                        continue;
                    }
                    if !reported {
                        listed.insert(list.oid, None);
                        continue;
                    }
                    let ordinal = u32::try_from(list.oid)?;
                    let stored = shared::read_title(
                        &database.db,
                        ordinal,
                        batch.query_ordinal + qidx,
                        PROGRAM,
                    )?;
                    shared::check_reported_title(stored.replaced_bytes, ordinal, PROGRAM)?;
                    listed.insert(list.oid, Some(Arc::from(stored.title)));
                }
            }
        }
        let mut reported = ReportedSubjects::default();
        for (&oid, title) in &listed {
            let ordinal = u32::try_from(oid)?;
            let shown = title.is_some();
            let title = title.clone().unwrap_or_else(|| Arc::from(&b""[..]));
            if pairwise && shown {
                shared::check_shown_protein_title(&title, ordinal, PROGRAM)?;
            }
            if tabular {
                // An unreported subject's id is never written.
                let id = if shown {
                    shared::tabular_subject_id(&title, ordinal, true, PROGRAM)?
                } else {
                    Vec::new()
                };
                // A stored title is a VisibleString (checked above), so its id is ASCII.
                reported
                    .tabular_ids
                    .push(String::from_utf8_lossy(&id).into_owned());
            }
            let encoded = subjects.read_encoded(oid)?;
            reported.records.push(FastaRecord {
                internal_id: format!("gnl|BL_ORD_ID|{ordinal}"),
                title: String::from_utf8_lossy(&title).into_owned(),
                sequence: encoded
                    .iter()
                    .map(|&residue| NCBISTDAA_TO_AMINOACID[usize::from(residue)])
                    .collect(),
                lowercase_masks: Vec::new(),
                warnings: Vec::new(),
                source: None,
            });
            reported.titles.push(title);
            reported.ordinals.push(ordinal);
        }
        // Ranks keep the order of the ordinal ids, so every tie broken by subject index
        // is broken as before.
        let rank: BTreeMap<usize, usize> = listed
            .keys()
            .enumerate()
            .map(|(rank, &oid)| (oid, rank))
            .collect();
        for batch in batches.iter_mut() {
            for lists in &mut batch.queries {
                for list in lists.iter_mut() {
                    list.oid = rank[&list.oid];
                }
            }
        }
        if taxids {
            // The taxid of each reported sequence, read apart from the titles for the XML2 and
            // JSON formats only (`algorithm::common::database::read_reported_taxids`); the
            // sequences with a title are the reported ones.
            let shown: std::collections::BTreeSet<usize> = listed
                .values()
                .enumerate()
                .filter(|(_, title)| title.is_some())
                .map(|(rank, _)| rank)
                .collect();
            // A scientific name NCBI cannot write, or a taxonomy database LOSAT cannot read as NCBI
            // does, is LOSAT's refusal of the XML2 document (a run that holds its output drops it).
            (reported.taxids, reported.taxonomy, reported.taxdb_notes) =
                shared::read_reported_taxids(
                    &database.db,
                    &database.taxdb,
                    &reported.ordinals,
                    &shown,
                    PROGRAM,
                )
                .map_err(super::blast_xml::XmlRefusal::of_taxonomy)?;
        }
        Ok(reported)
    }

    /// The records by rank.
    pub(crate) fn records(&self) -> &[FastaRecord] {
        &self.records
    }

    /// The ordinal id of rank `rank`.
    pub(crate) fn ordinal(&self, rank: usize) -> u32 {
        self.ordinals[rank]
    }

    /// The Blast-def-line taxid of the sequence of rank `rank` (0 for none, and for a run that
    /// does not read taxids).
    pub(crate) fn taxid(&self, rank: usize) -> i32 {
        self.taxids.get(rank).copied().unwrap_or(0)
    }

    /// What the taxonomy database says of that taxid (`HitTaxonomy::Unnamed` when nothing).
    pub(crate) fn taxonomy(&self, rank: usize) -> HitTaxonomy {
        self.taxonomy.get(rank).cloned().unwrap_or_default()
    }

    /// NCBI's messages about the taxonomy database that come before the `Hit`s' own, for the
    /// batch.
    pub(crate) fn taxdb_notes(&self) -> &shared::TaxdbNotes {
        &self.taxdb_notes
    }

    /// The tabular subject id of rank `rank` (only with `tabular` in `for_batch`).
    pub(crate) fn tabular_id(&self, rank: usize) -> &str {
        &self.tabular_ids[rank]
    }

    /// The stored title of rank `rank` (`stitle`).
    pub(crate) fn title(&self, rank: usize) -> &str {
        &self.records[rank].title
    }

    /// The stored titles by rank (the outfmt 0 alignment headings).
    pub(crate) fn stored_titles(&self) -> &[Arc<[u8]>] {
        &self.titles
    }

    /// The description-table labels by rank: the title as NCBI's defline generator makes it
    /// for the description table (as BLASTP's, BLASTN's and TBLASTN's `-db` reports).
    pub(crate) fn descriptions(&self) -> Vec<Arc<str>> {
        self.titles
            .iter()
            .map(|title| {
                let label = match crate::report::defline::generate_defline(title, true, true) {
                    Ok(defline) => defline.text,
                    Err(crate::report::defline::UnknownEncoding) => b"Unknown".to_vec(),
                };
                Arc::from(String::from_utf8_lossy(&label).as_ref())
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/db")
            .join(name)
            .to_str()
            .expect("fixture path is UTF-8")
            .to_string()
    }

    fn open(name: &str) -> BlastxDatabase {
        BlastxDatabase::open(&fixture(name), &ApplicationSettings::default())
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    #[test]
    fn the_stores_hold_the_bytes_the_reader_reads_in_ordinal_order() {
        // LOSAT/tests/seqdb/expected/prot_mv.info.tsv: 10 sequences, 801 residues, in six
        // volumes (the stores span them in ordinal order).
        let database = open("prot_mv");
        let stores = database.read_stores().unwrap();
        assert_eq!((stores.len(), stores.total_length()), (10, 801));
        let subjects = Subjects::Database(&stores);
        let mut letters = 0;
        for ordinal in 0..stores.len() {
            let block = database.db.read_sequence(ordinal as u32).unwrap();
            let stored = block.protein(0);
            assert_eq!(&*subjects.encoded(ordinal), stored, "{ordinal}");
            let with_sentinel = subjects.encoded_with_sentinel(ordinal);
            assert_eq!(&with_sentinel[..stored.len()], stored, "{ordinal}");
            assert_eq!(with_sentinel.last(), Some(&0), "{ordinal}");
            assert_eq!(subjects.length(ordinal), stored.len());
            letters += stored.len();
        }
        assert_eq!(letters, 801);
        assert_eq!(
            (subjects.total_length().unwrap(), subjects.shortest_length()),
            (801, 10)
        );
    }

    #[test]
    fn a_record_of_the_stored_letters_encodes_to_the_stored_bytes() {
        // Every residue byte 1..=27 (O = 26 included) reads back as itself through its
        // letter: the records of the reported sequences and `-subject` records alike.
        let bytes: Vec<u8> = (1..28).collect();
        let record = FastaRecord {
            internal_id: "gnl|BL_ORD_ID|0".into(),
            title: String::new(),
            sequence: bytes
                .iter()
                .map(|&residue| NCBISTDAA_TO_AMINOACID[usize::from(residue)])
                .collect(),
            lowercase_masks: Vec::new(),
            warnings: Vec::new(),
            source: None,
        };
        let records = [record];
        assert_eq!(&*Subjects::Fasta(&records).encoded(0), &bytes[..]);
    }

    #[test]
    fn an_empty_database_fails_as_ncbis_preliminary_search() {
        let error = open("empty").read_stores().unwrap_err();
        let error = error
            .downcast::<PreliminaryReadError>()
            .expect("NCBI's error of the preliminary search");
        assert_eq!(error.0, "Database memory map file error");
    }
}
