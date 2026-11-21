//! The BLAST database input of dustmasker and segmasker: `CMaskBDBReader`
//! (`c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp`), the sequences of a database in OID
//! order, each as `CSeqDB::GetBioseq` makes it, with the errors of the database reader as
//! NCBI's exceptions.
//!
//! # API
//!
//! - [`MaskBdbReader::open`]`(name, molecule, registry, title_printed, program)`:
//!   `CMaskBDBReader(name, is_nucl)` (`new CSeqDB(name, eNucleotide or eProtein)`).
//! - [`MaskBdbReader::next_record`]`()`: `GetNextSequence`, `Ok(None)` after the last OID.
//! - [`BdbRecord`]: the Seq-ids of the first definition line, the title joined from all of them
//!   (`CSeqDB::GetBioseq`), and the residues: IUPACna letters of a nucleotide
//!   (`eCoding_Iupac`), NCBIstdaa codes of a protein (`eCoding_Ncbi`; [`iupacaa_letters`] gives
//!   their `eCoding_Iupac` letters).
//! - [`BdbError`]: NCBI's exception (`Ncbi`, with the place NCBI throws it from where LOSAT
//!   knows it: [`NcbiException::site`] is [`ThrowSite::generated`] otherwise) or LOSAT's refusal
//!   (`Refused`, a noun phrase that the program completes).

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::db_subject::bioseq_facts;
use crate::report::defline::bioseq_title_from_source;
use crate::seqdb::alias::SearchPath;
use crate::seqdb::{SeqDb, SeqDbError, SeqType};
use crate::seqid::ParsedSeqId;

/// One sequence of the database (`CBioseq`).
pub(crate) struct BdbRecord {
    /// The Seq-ids of the first definition line.
    pub ids: Vec<ParsedSeqId>,
    /// The title (`CSeqDB::GetBioseq`'s title descriptor).
    pub title: Vec<u8>,
    /// The name of the organism of the Bioseq's first source, from which `GenerateDefline` makes
    /// an empty title and rewrites the end of a protein title ending in `]`
    /// (`report::defline::BioseqFacts::taxname`); looked up only for such a title that the
    /// writer prints.
    pub taxname: Option<std::sync::Arc<[u8]>>,
    /// Whether the Bioseq has a source: a definition line with a positive taxid
    /// (`report::defline::BioseqFacts::has_source`); asked only for a title as `taxname` (else
    /// `false`, which writes the same title).
    pub has_source: bool,
    /// IUPACna letters of a nucleotide, NCBIstdaa codes of a protein.
    pub residues: Vec<u8>,
    /// The `ERR_POST` texts of the taxonomy database that NCBI's `GetBioseq` writes while it
    /// reads the sequence ([`MaskBdbReader::next_record`]), to be posted before the sequence is
    /// written (`seqmasks_io::run::post_read_texts`).
    pub texts: Vec<String>,
}

/// Why the reader gave no record.
#[derive(Debug)]
pub(crate) enum BdbError {
    /// NCBI's exception.
    Ncbi(NcbiException),
    /// LOSAT's refusal: a whole message (`algorithm::common::database::refusal`), or a noun
    /// phrase (`Unsupported`) that the program completes with "is not supported by LOSAT's
    /// `<program>`".
    Refused(String),
    Unsupported(String),
}

/// `CMaskBDBReader`: the sequences of a database in OID order. Every OID of the volumes, also
/// those outside the alias files' OID ranges (`SeqDb::open_every_oid_in`).
///
/// NCBI reference (598d8ae6): c++/include/objtools/seqmasks_io/mask_bdb_reader.hpp:55-61
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp:41-55
pub(crate) struct MaskBdbReader {
    db: SeqDb,
    molecule: SeqType,
    oid: u32,
    /// Whether the writer prints the title (`interval`, `acclist`, `fasta`; the `seqloc_*` and
    /// `maskinfo_*` writers take only the Seq-id).
    title_printed: bool,
    /// The program's name in LOSAT's refusals.
    program: &'static str,
    /// The taxonomy database of the search path (`CSeqDBTaxInfo`, opened at the first lookup) and
    /// the lookups of the volumes' caches, each `GetBioseq` a load.
    taxdb: crate::algorithm::common::database::TaxdbRun,
    /// What `source_taxname` looked up, as NCBI's volume keeps it (`seqdb::taxdb::TaxCache`).
    source_names: crate::seqdb::taxdb::TaxCache,
}

/// Whether the writer of `outfmt` prints the title of a sequence (`CMaskWriter::IdToString`).
pub(crate) fn prints_title(outfmt: &str) -> bool {
    matches!(outfmt, "interval" | "acclist" | "fasta")
}

impl MaskBdbReader {
    /// `new CSeqDB(name, molecule)` on the search path of the environment and `registry` (the
    /// `[BLAST] BLASTDB` of the configuration: `ApplicationSettings::blastdb`, `Err` when LOSAT
    /// does not read the entry as NCBI does). `program` names LOSAT's program in a refusal.
    pub(crate) fn open(
        name: &str,
        molecule: SeqType,
        registry: Result<Option<&str>, &str>,
        title_printed: bool,
        program: &'static str,
    ) -> Result<Self, BdbError> {
        let registry = registry.map_err(|reason| {
            BdbError::Refused(
                crate::algorithm::common::database::refusal(reason, program).to_string(),
            )
        })?;
        let path = SearchPath::from_environment(registry)
            .map_err(|error| database_error(error, program))?;
        // `GetNumOIDs`, not the OIDs that the alias files' ranges include.
        let db = SeqDb::open_every_oid_in(name, molecule, &path)
            .map_err(|error| database_error(error, program))?;
        Ok(MaskBdbReader {
            db,
            molecule,
            oid: 0,
            title_printed,
            program,
            taxdb: crate::algorithm::common::database::TaxdbRun::sequential(),
            source_names: Default::default(),
        })
    }

    /// `GetNextSequence`: the sequence of the next OID, `None` after the last.
    pub(crate) fn next_record(&mut self) -> Result<Option<BdbRecord>, BdbError> {
        if self.oid >= self.db.num_seqs() {
            return Ok(None);
        }
        let oid = self.oid;
        self.oid += 1;
        let program = self.program;
        let error = |error| database_error(error, program);
        // `GetBioseq` builds NCBI's OID list first (`SeqDb::oid_list_ready`).
        self.db.oid_list_ready().map_err(error)?;
        let set = self.db.defline_set(oid).map_err(error)?;
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err(BdbError::Unsupported(format!(
                "sequence {oid} of the database: a definition line with bytes NCBI reads as '#' \
                 (it writes an error each time it reads them)"
            )));
        }
        let ids = set
            .deflines
            .first()
            .map(|defline| defline.ids.clone())
            .unwrap_or_default();
        let title = set.bioseq_title();
        let residues = match self.molecule {
            SeqType::Nucleotide => self.db.nucleotide_iupac(oid).map_err(error)?,
            SeqType::Protein => self
                .db
                .read_sequence(oid)
                .map_err(error)?
                .protein(0)
                .to_vec(),
        };
        // `GetBioseq` gives no Bioseq for a sequence without letters, before its lookups
        // (seqdbvol.cpp:1380-1386), and the writer prints none. A sequence without a positive
        // taxid has no lookup (`TaxdbRun::bioseq_texts`).
        let mut texts = Vec::new();
        if !residues.is_empty() {
            texts = self
                .taxdb
                .bioseq_texts(&self.db, oid, &set, None)
                .map_err(|error| lookup_refusal(oid, error))?;
        }
        // NCBI makes the title of a sequence without one from its BioSource (the taxonomy
        // names of a definition line with a taxid; seqdbvol.cpp:1100-1188), and rewrites the end
        // of a protein title ending in `]` from it (`x_AdjustProteinTitleSuffix`,
        // `report::defline::adjust_protein_title_suffix`, which runs for any source: probe of
        // line RG T5c, `T/t5c/out/seg`, `>alpha [Bacillus sp.]` of `alpha, partial [Bacillus
        // sp.]`, `>a` of `a, [ ]` without a name), for a Bioseq it makes (one with letters).
        // The title is asked first: a printed title that NCBI neither makes from the source nor
        // rewrites the end of (the reader's usual one) asks nothing of the source.
        let source_shown = self.title_printed
            && !residues.is_empty()
            && (bioseq_title_from_source(&title, || bioseq_facts(&ids))
                || (self.molecule == SeqType::Protein
                    && crate::report::defline::protein_title_suffix_adjustable(&title)));
        let has_source = source_shown && set.first_source_taxid(None).is_some();
        let taxname = if has_source {
            self.source_taxname(oid, &set)?
        } else {
            None
        };
        Ok(Some(BdbRecord {
            ids,
            title,
            taxname,
            has_source,
            residues,
            texts,
        }))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1126-1155
    // The first source names the organism (`DeflineSet::first_source_taxid`,
    // `seqdb::taxdb::SourceNames::taxname`). The load of the sequence (`TaxdbRun::bioseq_texts`,
    // in `next_record` before this) made the lookups and has their texts; the name is asked again
    // here, through the cache of the lookups (`source_names`, NCBI's `m_TaxCache`), and its texts
    // are dropped. A lookup that faults or whose names NCBI reads past its mapping of `taxdb.btd`
    // is refused at the first that meets it (the cache keeps none of those).
    /// The name of the organism of the first source of the sequence `oid`, whose definition lines
    /// are `set`, from which NCBI makes its empty title or the end of its protein title.
    fn source_taxname(
        &self,
        oid: u32,
        set: &crate::seqid::DeflineSet,
    ) -> Result<Option<std::sync::Arc<[u8]>>, BdbError> {
        let Some(taxid) = set.first_source_taxid(None) else {
            return Ok(None);
        };
        let names = self
            .source_names
            .names(taxid, &mut Vec::new(), |messages| {
                self.taxdb.lookup(&self.db, taxid, messages)?.into_names()
            })
            .map_err(|error| lookup_refusal(oid, error))?;
        Ok(names.taxname().cloned())
    }
}

/// LOSAT's refusal of the sequence `oid` whose taxonomy lookup failed with `error`: a record of an
/// empty `taxdb.btd` (NCBI maps it as a null pointer and faults, PD-LOSAT-NCBI-DEFECTS rule 3), a
/// title made from names that NCBI reads past its mapping of `taxdb.btd` (its bytes, or whether
/// NCBI faults, depend on its memory map), or a taxonomy database that does not open or read.
fn lookup_refusal(oid: u32, error: SeqDbError) -> BdbError {
    BdbError::Unsupported(format!(
        "sequence {oid} of the database, whose taxid's lookup in the taxonomy database fails \
         ({error}),"
    ))
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1785-1813
/// The `eCoding_Iupac` letters of a protein's NCBIstdaa codes (`CSeqVector`'s IUPACaa): the gap
/// (0) and the stop (25) are the byte 255, as NCBI's table has them (segmasker 2.17.0 writes
/// them so with `-infmt blastdb -outfmt fasta`, fixture `db_ambig_fasta`).
pub(crate) fn iupacaa_letters(stdaa: &[u8]) -> Vec<u8> {
    crate::sequence::stdaa_letters::protein_string_letters(stdaa)
}

/// An error of the database reader: NCBI's exception (`SeqDbError::ncbi_exception`, or the
/// `CSeqDBException` of `SeqDbError::ncbi_message`, with its place where LOSAT knows it,
/// `seqdb::throw_site`), or LOSAT's refusal.
fn database_error(error: SeqDbError, program: &str) -> BdbError {
    // An NCBI exception of the reader itself (an alias value that does not convert), reported
    // as the application reports it (blastdbcmd's `database_error`).
    if let Some(exception) = error.ncbi_exception() {
        return BdbError::Ncbi(exception.clone());
    }
    match error.ncbi_request_message() {
        Some(text) => {
            let site = crate::seqdb::throw_site(&text)
                .unwrap_or(ThrowSite::generated("CSeqDBException::eInvalid"));
            BdbError::Ncbi(NcbiException::new(site, text))
        }
        None => BdbError::Refused(
            crate::algorithm::common::database::refusal(&error.to_string(), program).to_string(),
        ),
    }
}
