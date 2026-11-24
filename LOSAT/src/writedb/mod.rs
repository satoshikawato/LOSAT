//! Writer of NCBI BLAST databases (format version 5, makeblastdb's default,
//! and version 4) and the `losat makeblastdb` command.
//!
//! NCBI's writer is `c++/src/objtools/blast/seqdb_writer/` (CWriteDB,
//! CBuildDatabase) under the application `c++/src/app/blastdb/makeblastdb.cpp`.
//! The files written are those makeblastdb writes for a FASTA input without
//! parsed Seq-ids: the index file (`.pin` / `.nin`), the sequence file
//! (`.psq` / `.nsq`), the header file (`.phr` / `.nhr`), an alias file
//! (`.pal` / `.nal`) when the database has more than one volume, the JSON
//! metadata file (`.pjs` / `.njs`), and for version 5 the LMDB files
//! (`.pdb` / `.ndb` with the volume tables, `.ptf` / `.ntf` with the tax ids)
//! and the tax id lookup files (`.pot` / `.not`, `.pto` / `.nto`).
//!
//! With `-parse_seqids` and `-hash_index`, the volumes also get their ISAM
//! files (`.psi`/`.psd` of the Seq-ids in version 4, `.pni`/`.pnd` of the GIs,
//! `.phi`/`.phd` of the sequence hashes), their OID to GI file (`.pog`), in
//! version 5 the acc2oid keys of the Seq-ids (in `.pdb` / `.ndb`) and the OID to
//! Seq-ids file (`.pos` / `.nos`), and the duplicate check of their Seq-ids
//! (`volume`, `isam`, `id_set`); `-taxid` and `-taxid_map` (which requires
//! `-parse_seqids`) put tax ids on the definition lines and into the version 5
//! tax id files (`taxid_set`).
//!
//! Mask data (`-mask_data`, `-mask_id`, `-mask_desc`): the library writes the mask column of a
//! volume (the files `.?aa`, `.?ab`, `.?ac`, version 4 and 5 alike) from the mask files that
//! `mask_info` reads (ASN.1 text or binary, as dustmasker, segmasker, windowmasker and
//! convert2blastmask write them) and `mask_data` turns into ranges per sequence, with or without
//! `-parse_seqids`. With `-gi_mask` (and `-gi_mask_name`) the masks go instead to the GI-indexed
//! files of each algorithm (`gi_mask`), and the database gets an alias file that names them.
//! `losat makeblastdb` refuses no option of makeblastdb 2.17.0 by name today (no row of
//! `ncbi_app::args::refusal`); `-oid_masks` writes the OID mask files (`oid_list`);
//! `-metadata_output_prefix` names the files of the metadata file
//! (`metadata::metadata_user_path`).
//!
//! `-input_type blastdb` copies existing databases (`raw_source`, `CRawSeqDBSource`): each
//! sequence in the form of the sequence file with its definition lines (and their PIG), the
//! ranges of its masking algorithms and the blobs of its other columns, through the same
//! `BuildDatabase` and `WriteDb` as a FASTA input; with `-randomize`, in the order of NCBI's
//! shuffle (`crate::utils::std_shuffle`).

pub(crate) mod alias_file;
pub(crate) mod binary_list;
pub mod bioseq;
pub mod build_db;
pub(crate) mod column;
pub mod convert;
pub(crate) mod cout;
pub mod data_size;
pub mod deflines;
pub(crate) mod dir_entry;
pub mod fasta_source;
pub mod files;
pub mod format_guess;
pub(crate) mod gi_mask;
pub mod header;
pub(crate) mod id_set;
pub(crate) mod isam;
pub(crate) mod lmdb;
pub mod local_offset;
pub mod makeblastdb;
pub(crate) mod mask_data;
pub(crate) mod mask_info;
#[allow(dead_code)]
pub(crate) mod mdb;
pub(crate) mod metadata;
pub(crate) mod ofstream;
pub(crate) mod oid_list;
pub mod platform;
pub(crate) mod raw_source;
pub(crate) mod seqidlist_writer;
pub mod taxid_set;
pub mod time;
pub mod volume;
mod windows_zone;
pub mod writer;

#[cfg(test)]
mod tests;

use crate::wording::{DATABASE_CREATION_ERROR, DATABASE_ERROR, OPTIONS_ERROR};
use std::path::PathBuf;

use crate::ncbi_exception::NcbiException;

pub use platform::Platform;

/// Error of the database writer.
#[derive(Debug)]
pub enum WriteDbError {
    /// NCBI's `CWriteDBException`: the text makeblastdb prints after
    /// "BLAST Database creation error: ".
    Creation(String),
    /// `Creation` whose text holds bytes that are not UTF-8 (a Seq-id label of the
    /// duplicate check), which makeblastdb prints as they are.
    CreationBytes(Vec<u8>),
    /// LOSAT's refusal of a version 5 database with parsed Seq-ids in which the sequence of OID
    /// `oid` (first id `name`, its FASTA label) has GIs only and lies between sequences with
    /// other ids: NCBI's writer aborts at its close (`writer::gi_only_between_keys_error`).
    /// makeblastdb reports it as a creation error.
    GiOnlyBetweenKeys { oid: u32, name: Vec<u8> },
    /// NCBI's `std::invalid_argument` with its `what()`: `BuildLMDBFileName` of an empty
    /// database name (`writer::build_lmdb_file_name`).
    InvalidArgument(&'static str),
    /// LOSAT's refusal of given definition lines whose first line has no Seq-id in a database
    /// without parsed Seq-ids: NCBI's writer puts the `gnl|BL_ORD_ID` id past the end of the
    /// empty list (undefined behaviour). makeprofiledb gives such lines (a query without ids);
    /// makeblastdb's lines always have an id.
    OrdinalIdWithoutIds,
    /// NCBI's `CInvalidDataException`: the text makeblastdb prints after
    /// "BLAST options error: ".
    Options(String),
    /// An exception of NCBI's C++ Toolkit that makeblastdb reports through the `CException`
    /// clause of its general handler: "Error: " and the exception's report
    /// (`CException::what()`), exit status 255.
    Exception(NcbiException),
    /// NCBI's `CSeqDBException` (the LMDB environment that could not be opened, `CBlastEnv`):
    /// the text makeblastdb prints after "BLAST Database error: ", exit status 2. At the close
    /// it is a creation error, as `EndBuild` makes every exception of `Close`.
    Database(String),
    /// LOSAT's own refusals (an input on which makeblastdb crashes, does not end or cannot be
    /// followed), printed after "Error: " with exit status 255, and the errors of the LMDB
    /// layer after its environment is open (`lmdb::error`, NCBI's `CSeqDBException` texts of
    /// `CWriteDB_LMDB`), which occur only in the close and so become creation errors; and the
    /// `what()` of a `std::exception` that the general handler reports the same way (the
    /// `std::length_error` of a mask file's string longer than `std::string` holds).
    Other(String),
    /// LOSAT's refusal of an input that makeblastdb reads, raised inside the build (mask data
    /// whose ranges NCBI orders by memory addresses, `mask_data::range_set`, at the close of a
    /// build that would otherwise succeed; a mask column created after the first volume): "losat
    /// makeblastdb: " and the message, exit status 1, like `MakedbError::Unsupported`. Both end
    /// the build as an error of NCBI's build at that point does (`EndBuild` with `erase` true: its
    /// files are removed); a registration after the first volume (makeblastdb registers every
    /// algorithm before the input) is refused before anything changes and leaves the build open.
    Unsupported(String),
    /// An allocation failed (`std::bad_alloc`): makeblastdb prints "BLAST ran
    /// out of memory".
    OutOfMemory,
    /// The operating system failed to open an input file, or to write a data file of
    /// LMDB (`mdb`). The other output files of the database are not tested, as in
    /// makeblastdb (`ofstream`).
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-259
// NCBI's writer does not test whether an output file opened, a write succeeded or a
// rename took place, so makeblastdb has no message for these failures, and LOSAT's writer
// does not test them either (`ofstream`, `dir_entry`). An input file that cannot be opened
// or read, and an error of LMDB's own writes, LOSAT reports in the form of makeblastdb's
// "other" errors ("Error: " and exit status 255, `makeblastdb::MakedbError::exit_code`),
// with the path and the system's message; they are not errors of the database's contents
// ("BLAST Database creation error").
impl WriteDbError {
    /// A creation error with the message `message`, whose bytes are kept when they are not UTF-8.
    pub(crate) fn creation(message: Vec<u8>) -> Self {
        match String::from_utf8(message) {
            Ok(text) => WriteDbError::Creation(text),
            Err(error) => WriteDbError::CreationBytes(error.into_bytes()),
        }
    }

    /// The message makeblastdb prints for the error (without the line feed that ends it), as
    /// bytes: a name or a word of the input that is not UTF-8 keeps its bytes, as in NCBI's
    /// `std::string`.
    pub fn message(&self) -> Vec<u8> {
        let text = |prefix: &str, message: &[u8]| [prefix.as_bytes(), message].concat();
        // LOSAT's prefixes (`crate::wording`) and a space, as NCBI's `BLAST ... error: `.
        let losat = |prefix: &str, message: &[u8]| [prefix.as_bytes(), b" ", message].concat();
        match self {
            WriteDbError::Creation(message) => losat(DATABASE_CREATION_ERROR, message.as_bytes()),
            WriteDbError::CreationBytes(message) => losat(DATABASE_CREATION_ERROR, message),
            WriteDbError::GiOnlyBetweenKeys { oid, name } => losat(
                DATABASE_CREATION_ERROR,
                &writer::gi_only_between_keys_message(*oid, name),
            ),
            WriteDbError::InvalidArgument(what) => text("Error: ", what.as_bytes()),
            WriteDbError::OrdinalIdWithoutIds => text(
                "losat makeblastdb: ",
                b"a first definition line without a Seq-id, whose gnl|BL_ORD_ID id NCBI's writer \
                  puts past the end of the empty list",
            ),
            WriteDbError::Options(message) => losat(OPTIONS_ERROR, message.as_bytes()),
            WriteDbError::Database(message) => losat(DATABASE_ERROR, message.as_bytes()),
            // `what()` ends with a line feed and the message gets one more: the report ends
            // with an empty line.
            WriteDbError::Exception(e) => text("Error: ", &e.what()),
            WriteDbError::Other(message) => text("Error: ", message.as_bytes()),
            WriteDbError::Unsupported(message) => text("losat makeblastdb: ", message.as_bytes()),
            WriteDbError::OutOfMemory => crate::wording::OUT_OF_MEMORY.as_bytes().to_vec(),
            WriteDbError::Io { path, source } => [
                b"Error: ".as_slice(),
                &crate::utils::os_bytes::path_bytes(path),
                format!(": {source}").as_bytes(),
            ]
            .concat(),
        }
    }
}

impl std::fmt::Display for WriteDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.message()))
    }
}

impl std::error::Error for WriteDbError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            WriteDbError::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}
