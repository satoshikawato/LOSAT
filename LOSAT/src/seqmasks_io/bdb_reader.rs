//! The BLAST database input of dustmasker and segmasker: `CMaskBDBReader`
//! (`c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp`), the sequences of a database in OID
//! order, each as `CSeqDB::GetBioseq` makes it, with the errors of the database reader as
//! NCBI's exceptions.
//!
//! # API
//!
//! - [`MaskBdbReader::open`]`(name, molecule, registry, title_printed, program)`:
//!   `CMaskBDBReader` of `name` and `is_nucl` (a new `CSeqDB` of `name`, nucleotide or protein).
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
use crate::seqdb::{NucleotideAmbiguity, SeqDb, SeqDbError, SeqType};
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
    /// `sequence::GetTitle` of the Bioseq (`sequence::title::get_title`) of a sequence with
    /// letters, when the reader was asked for it (`MaskBdbReader::ask_get_title`: windowmasker's
    /// `-ids` matched as text, makembindex's `-idmap`). An empty title is made from the Bioseq's
    /// first source, which it has for a definition line with a positive taxid whatever the taxonomy
    /// database holds (seqdbvol.cpp:1126-1176).
    pub get_title: Option<Vec<u8>>,
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
    /// The search path the database was opened on (`volume_paths`).
    path: SearchPath,
    /// Whether the records carry `sequence::GetTitle` (`BdbRecord::get_title`).
    get_title: bool,
    /// Whether the last record's read of its header may have met a UTF8String for a VisibleString,
    /// for which NCBI posts its warning 803.10 (`take_utf8_header`): it met one, or it failed (the
    /// record is read again to see whether one came before the failure).
    utf8_header: bool,
}

/// Whether the writer of `outfmt` prints the title of a sequence (`CMaskWriter::IdToString`).
pub(crate) fn prints_title(outfmt: &str) -> bool {
    matches!(outfmt, "interval" | "acclist" | "fasta")
}

impl MaskBdbReader {
    /// A new `CSeqDB` of `name` (the bytes of NCBI's `std::string`) and the molecule, on the
    /// search path of the environment and `registry` (the `[BLAST] BLASTDB` of the
    /// configuration: `ApplicationSettings::blastdb`). `program` names LOSAT's program in a
    /// refusal.
    pub(crate) fn open(
        name: &[u8],
        molecule: SeqType,
        registry: Option<&[u8]>,
        title_printed: bool,
        program: &'static str,
    ) -> Result<Self, BdbError> {
        let path = SearchPath::from_environment(registry);
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
            path,
            get_title: false,
            utf8_header: false,
        })
    }

    /// From now on the records carry `sequence::GetTitle` of their Bioseq
    /// (`BdbRecord::get_title`).
    pub(crate) fn ask_get_title(&mut self) {
        self.get_title = true;
    }

    /// Whether the read of the last record's header may have met a UTF8String for a
    /// VisibleString (`SeqDb::take_utf8_post_at` finds its member and whether the run posts it
    /// now): makembindex's post level shows NCBI's warning 803.10, the maskers' does not.
    pub(crate) fn take_utf8_header(&mut self) -> bool {
        std::mem::take(&mut self.utf8_header)
    }

    /// The OID that `next_record` reads next (makembindex's `CSequenceIStreamBlastDB::oid_`).
    pub(crate) fn next_oid(&self) -> u32 {
        self.oid
    }

    /// `CSequenceIStreamBlastDB::putback()`: the last OID is read again (none before OID 0).
    // NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/sequence_istream_bdb.cpp:138-139
    pub(crate) fn putback(&mut self) {
        self.oid = self.oid.saturating_sub(1);
    }

    /// The database (its mask algorithms and mask data: makembindex's `-db_mask`).
    pub(crate) fn db(&self) -> &SeqDb {
        &self.db
    }

    /// `FindVolumePaths` (with its flag true) of the database opened as `name`: the paths of its
    /// volumes (without extension), each once, in NCBI's order.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:929-966
    pub(crate) fn volume_paths(&self, name: &[u8]) -> Result<Vec<Vec<u8>>, BdbError> {
        SeqDb::volume_and_alias_paths(name, self.molecule, &self.path)
            .map(|(volumes, _)| volumes)
            .map_err(|error| database_error(error, self.program))
    }

    /// The exception or refusal of `error`, as `next_record` reports it.
    pub(crate) fn database_error(&self, error: SeqDbError) -> BdbError {
        database_error(error, self.program)
    }

    /// `GetNextSequence`: the sequence of the next OID, `None` after the last.
    pub(crate) fn next_record(&mut self) -> Result<Option<BdbRecord>, BdbError> {
        Ok(self.next_record_of(true)?.map(|(record, _)| record))
    }

    /// [`MaskBdbReader::next_record`] of a nucleotide database without the letters (`residues`
    /// is empty): the length and the ambiguity runs of the sequence come with the record, so
    /// that the caller reads the letters by ranges ([`SeqDb::nucleotide_iupac_range_with`],
    /// makembindex: the whole sequence is never held as letters).
    pub(crate) fn next_record_by_ranges(
        &mut self,
    ) -> Result<Option<(BdbRecord, NucleotideAmbiguity)>, BdbError> {
        Ok(self
            .next_record_of(false)?
            .map(|(record, ambiguity)| (record, ambiguity.expect("a nucleotide database"))))
    }

    /// `next_record` with the letters (`letters`), or with the length and the ambiguity runs of
    /// a nucleotide instead.
    fn next_record_of(
        &mut self,
        letters: bool,
    ) -> Result<Option<(BdbRecord, Option<NucleotideAmbiguity>)>, BdbError> {
        self.utf8_header = false;
        if self.oid >= self.db.num_seqs() {
            return Ok(None);
        }
        let oid = self.oid;
        self.oid += 1;
        let program = self.program;
        let error = |error| database_error(error, program);
        // `GetBioseq` builds NCBI's OID list first (`SeqDb::build_oid_list`).
        self.db.build_oid_list().map_err(error)?;
        let set = self.db.defline_set(oid);
        self.utf8_header = set.as_ref().map_or(true, |set| set.utf8);
        let set = set.map_err(error)?;
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
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1362-1366
        // `GetBioseq` gives no Bioseq for a sequence whose every definition line the alias files'
        // `OID_MASK_TYPE` erased (line RG D1b-2), before it reads the letters: the record of a
        // sequence without letters, which segmasker and windowmasker skip (probe S2_seg_walk of
        // line RG D1b), dustmasker stops at (probes dust_omnx_int, dust_ntall_fasta of D1b-2) and
        // makembindex refuses.
        if set.deflines.is_empty() {
            let record = BdbRecord {
                ids: Vec::new(),
                title: Vec::new(),
                taxname: None,
                has_source: false,
                get_title: None,
                residues: Vec::new(),
                texts: Vec::new(),
            };
            let ambiguity =
                (self.molecule == SeqType::Nucleotide && !letters).then(NucleotideAmbiguity::none);
            return Ok(Some((record, ambiguity)));
        }
        let ids = set
            .deflines
            .first()
            .map(|defline| defline.ids.clone())
            .unwrap_or_default();
        let title = set.bioseq_title();
        let (residues, ambiguity) = match self.molecule {
            SeqType::Nucleotide if letters => (self.db.nucleotide_iupac(oid).map_err(error)?, None),
            SeqType::Nucleotide => (
                Vec::new(),
                Some(self.db.nucleotide_ambiguity(oid).map_err(error)?),
            ),
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1380-1386
            // NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp:228-232
            // `GetBioseq` returns no Bioseq when the signed 32-bit protein length is below 1
            // (`0` or `>= 2^31` residues), after opening `.psq` via `x_GetSequence`.
            SeqType::Protein => {
                if self.db.protein_length_i32(oid).map_err(error)? < 1 {
                    self.db.check_protein_range(oid).map_err(error)?;
                    (Vec::new(), None)
                } else {
                    (
                        self.db
                            .read_sequence(oid)
                            .map_err(error)?
                            .protein(0)
                            .to_vec(),
                        None,
                    )
                }
            }
        };
        let has_letters = ambiguity
            .as_ref()
            .map_or(!residues.is_empty(), |ambiguity| ambiguity.length() > 0);
        // `GetBioseq` gives no Bioseq for a sequence without letters, before its lookups
        // (seqdbvol.cpp:1380-1386), and the writer prints none. A sequence without a positive
        // taxid has no lookup (`TaxdbRun::bioseq_texts`).
        let mut texts = Vec::new();
        if has_letters {
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
            && has_letters
            && (bioseq_title_from_source(&title, || bioseq_facts(&ids))
                || (self.molecule == SeqType::Protein
                    && crate::report::defline::protein_title_suffix_adjustable(&title)));
        let any_source = set.first_source_taxid(None).is_some();
        let has_source = source_shown && any_source;
        let taxname = if has_source {
            self.source_taxname(oid, &set)?
        } else {
            None
        };
        // `GetBioseq` gives no Bioseq for a sequence without letters, before its lookups: none
        // is made, and the callers skip the record (windowmasker) or refuse it (makembindex).
        let get_title = if self.get_title && has_letters {
            Some(crate::sequence::title::get_title(
                &ids,
                &title,
                any_source,
                || self.source_taxname(oid, &set),
            )?)
        } else {
            None
        };
        Ok(Some((
            BdbRecord {
                ids,
                title,
                taxname,
                has_source,
                get_title,
                residues,
                texts,
            },
            ambiguity,
        )))
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

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1785-1818
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
        None => BdbError::Refused(crate::algorithm::common::database::losat_reader_error(
            &error, program,
        )),
    }
}
