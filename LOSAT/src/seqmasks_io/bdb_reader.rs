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
use crate::report::db_subject::bioseq_ids;
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
    /// IUPACna letters of a nucleotide, NCBIstdaa codes of a protein.
    pub residues: Vec<u8>,
}

/// Why the reader gave no record.
#[derive(Debug)]
pub(crate) enum BdbError {
    /// NCBI's exception.
    Ncbi(NcbiException),
    /// LOSAT's refusal: a whole message (`algorithm::common::database::refusal`), or a noun
    /// phrase (`Unsupported`) that the program completes with "is not supported by LOSAT's
    /// <program>".
    Refused(String),
    Unsupported(String),
}

/// `CMaskBDBReader`: the sequences of a database in OID order. Every OID of the volumes, also
/// those outside the alias files' OID ranges (`SeqDb::open_every_oid_in`).
///
/// NCBI reference (598d8ae6): c++/include/objtools/seqmasks_io/mask_bdb_reader.hpp:55-61
/// ```c++
///     CMaskBDBReader( const string & dbname, bool is_nucl = true)
///         : CMaskReader( cin ),
///           seqdb_( new CSeqDB( dbname, (is_nucl
///                                        ? CSeqDB::eNucleotide
///                                        : CSeqDB::eProtein)) ),
///           oid_( 0 )
///     {}
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_bdb_reader.cpp:41-55
/// ```c++
/// CRef< CSeq_entry > CMaskBDBReader::GetNextSequence()
/// {
///     CRef< CSeq_entry > entry( null );
///
///     if( oid_ < seqdb_->GetNumOIDs() ) {
///         CRef< CBioseq > seq = seqdb_->GetBioseq( oid_++ );
///         entry.Reset( new CSeq_entry );
///         if( seq.NotEmpty() ) entry->SetSeq( *seq );
///     }
///
///     return entry;
/// }
/// ```
pub(crate) struct MaskBdbReader {
    db: SeqDb,
    molecule: SeqType,
    oid: u32,
    /// Whether the writer prints the title (`interval`, `acclist`, `fasta`; the `seqloc_*` and
    /// `maskinfo_*` writers take only the Seq-id).
    title_printed: bool,
    /// The program's name in LOSAT's refusals.
    program: &'static str,
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
        // NCBI makes the title of a sequence without one from its BioSource (the taxonomy
        // names of a definition line with a taxid; seqdbvol.cpp:1055-1160).
        let has_source = set
            .deflines
            .iter()
            .any(|defline| defline.taxid.is_some_and(|taxid| taxid > 0));
        if self.title_printed && has_source && bioseq_title_from_source(&title, &bioseq_ids(&ids)) {
            return Err(BdbError::Unsupported(format!(
                "sequence {oid} of the database, which has no title and a taxid (NCBI makes its \
                 title from the taxonomy),"
            )));
        }
        let residues = match self.molecule {
            SeqType::Nucleotide => {
                let length = self.db.nucleotide_length(oid).map_err(error)?;
                self.db
                    .nucleotide_iupac_range(oid, 0, length)
                    .map_err(error)?
            }
            SeqType::Protein => self
                .db
                .read_sequence(oid)
                .map_err(error)?
                .protein(0)
                .to_vec(),
        };
        Ok(Some(BdbRecord {
            ids,
            title,
            residues,
        }))
    }
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1785-1813
// ```c++
// const Uint1 CStdaaToIupacaa::scm_Table[32] = {
//     255, 'A', 'B', 'C', 'D', 'E', 'F', 'G',
//     'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
//     'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z',
//     'U', 255, 'O', 'J', 255, 255, 255, 255
// };
// ```
/// The `eCoding_Iupac` letters of a protein's NCBIstdaa codes (`CSeqVector`'s IUPACaa): the gap
/// (0) and the stop (25) are the byte 255, as NCBI's table has them (segmasker 2.17.0 writes
/// them so with `-infmt blastdb -outfmt fasta`, fixture `db_ambig_fasta`).
pub(crate) fn iupacaa_letters(stdaa: &[u8]) -> Vec<u8> {
    crate::sequence::stdaa_letters::protein_string_letters(stdaa)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:444-450
// ```c++
// const char* CException::GetErrCodeString(void) const
// {
//     switch (GetErrCode()) {
//     case eUnknown: return "eUnknown";
//     default:       return "eInvalid";
//     }
// }
// ```
// (`CSeqDBException`, seqdbcommon.hpp:73-100, does not name its codes, so every one is reported
// as `eInvalid`: `(CSeqDBException::eInvalid) No alias or index file found ...` for the
// `eFileErr` of seqdbalias.cpp:311-321, fixture db_missing.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:311-322
// ```c++
//             ostringstream oss;
//             oss << "No alias or index file found for " << p_or_n
//                 << " database [" << m_DBList[i].GetBasePathS()
//                 << "] in search path [" << search_path << "]. "
//                 << "Please verify the spelling of the BLAST "
//                 << "database and its molecule type.";
//
//             string msg(oss.str());
//
//             NCBI_THROW(CSeqDBException,
//                        eFileErr,
//                        msg);
// ```
// (`NCBI_THROW` is on line 320; its function is `BLASTDB::ncbi::CSeqDBAliasNode::x_ResolveNames`,
// the module `BLASTDB` before it: segmasker 2.17.0's `e.what()`, fixture db_missing.)
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:740-752, 921-925
// ```c++
//     if (m_DBList.empty()) {
//         ...
//             situation = "passed to CSeqDB::CSeqDB().";
//         ...
//             situation = string("found in alias file [")
//                 + this_name.GetBasePathS() + "].";
//         ...
//         NCBI_THROW(CSeqDBException,
//                    eArgErr,
//                    string("No database names were ") + situation);
// ...
//         oss << "Could not find volume or alias file ("
//             << m_DBList[i].GetBasePathS() << ") referenced in alias file ("
//             << this_name.GetBasePathS() << ").";
//
//         NCBI_THROW(CSeqDBException, eFileErr, oss.str());
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:155-164 (the
// constructor of `new CSeqDB(dbname, seqtype)`)
// ```c++
// CSeqDB::CSeqDB(const string & dbname,
//                ESeqType       seqtype,
//                CSeqDBGiList * gi_list,
//                bool           use_atlas_lock)
//
// {
//     if (dbname.size() == 0) {
//         NCBI_THROW(CSeqDBException,
//                    eArgErr,
//                    "Database name is required.");
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:124-128
// ```c++
//     if (! m_File.Open(CSeqDB_Path(m_FileName))) {
//         //m_Atlas.Unlock(locked);
//
//         string msg = string("Error: File (") + m_FileName + ") not found.";
//
//         NCBI_THROW(CSeqDBException, eFileErr, msg);
// ```
// (Each of these texts has this one place in NCBI's database reader; segmasker 2.17.0 names the
// lines 750, 925, 162 and 128 and these functions, fixtures db_alias_*, db_empty_name,
// db_space_name; review of line RG item K2, L2.)
/// The places of NCBI's throws of a database error, by the start of its text.
const DATABASE_ERROR_SITES: [(&str, ThrowSite); 5] = [
    (
        "No alias or index file found for ",
        seqdb_site(
            SEQDBALIAS,
            320,
            "BLASTDB::ncbi::CSeqDBAliasNode::x_ResolveNames",
        ),
    ),
    (
        "Could not find volume or alias file (",
        seqdb_site(
            SEQDBALIAS,
            925,
            "BLASTDB::ncbi::CSeqDBAliasNode::x_ExpandAliases",
        ),
    ),
    (
        "No database names were ",
        seqdb_site(
            SEQDBALIAS,
            750,
            "BLASTDB::ncbi::CSeqDBAliasNode::x_ExpandAliases",
        ),
    ),
    (
        "Database name is required.",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdb.cpp",
            162,
            "BLASTDB::ncbi::CSeqDB::CSeqDB",
        ),
    ),
    (
        "Error: File (",
        seqdb_site(
            "src/objtools/blast/seqdb_reader/seqdbfile.cpp",
            128,
            "BLASTDB::ncbi::CSeqDBExtFile::CSeqDBExtFile",
        ),
    ),
];

/// NCBI's alias reader.
const SEQDBALIAS: &str = "src/objtools/blast/seqdb_reader/seqdbalias.cpp";

/// The place of a `CSeqDBException` (its codes are all reported as `eInvalid`).
const fn seqdb_site(file: &'static str, line: u32, function: &'static str) -> ThrowSite {
    ThrowSite {
        file,
        line,
        kind: "CSeqDBException::eInvalid",
        function: Some(function),
    }
}

/// An error of the database reader: NCBI's exception (`SeqDbError::ncbi_exception`, or the
/// `CSeqDBException` of `SeqDbError::ncbi_message`, with its place where LOSAT knows it), or
/// LOSAT's refusal.
fn database_error(error: SeqDbError, program: &str) -> BdbError {
    // An NCBI exception of the reader itself (an alias value that does not convert), reported
    // as the application reports it (blastdbcmd's `database_error`).
    if let Some(exception) = error.ncbi_exception() {
        return BdbError::Ncbi(exception.clone());
    }
    match error.ncbi_message() {
        Some(text) => {
            let site = DATABASE_ERROR_SITES
                .iter()
                .find(|(start, _)| text.starts_with(start))
                .map_or(
                    ThrowSite::generated("CSeqDBException::eInvalid"),
                    |&(_, site)| site,
                );
            BdbError::Ncbi(NcbiException::new(site, text))
        }
        None => BdbError::Refused(
            crate::algorithm::common::database::refusal(&error.to_string(), program).to_string(),
        ),
    }
}
