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
//! `-parse_seqids`. `losat makeblastdb` refuses `-gi_mask`, `-gi_mask_name`, `-oid_masks` and
//! the other options of `makeblastdb::UNSUPPORTED_OPTIONS`.

pub mod bioseq;
pub mod build_db;
pub mod convert;
pub mod data_size;
pub mod deflines;
pub(crate) mod dir_entry;
pub mod fasta_source;
pub mod files;
pub mod format_guess;
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
pub mod platform;
pub mod taxid_set;
pub mod time;
pub mod volume;
pub mod writer;

#[cfg(test)]
mod tests;

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
    /// `CWriteDB_LMDB`), which occur only in the close and so become creation errors.
    Other(String),
    /// LOSAT's refusal of an input that makeblastdb reads, raised inside the build (mask data
    /// whose ranges NCBI orders by memory addresses, `mask_data::range_set`): "losat
    /// makeblastdb: " and the message, exit status 1, like `MakedbError::Unsupported`.
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
// ```c++
//     catch (const CException& e) {                                           \
//         LOG_POST(Error << "Error: " << e.what());                           \
//         exit_code = BLAST_UNKNOWN_ERROR;                                    \
//     }                                                                       \
//     catch (const std::ios::failure&) {                                      \
//         LOG_POST(Error << "BLAST failed to write output");                  \
//         exit_code = BLAST_OUTPUT_ERROR;                                     \
//     }                                                                       \
//     catch (const std::bad_alloc&) {                                         \
//         LOG_POST(Error << "BLAST ran out of memory");                       \
//         exit_code = BLAST_OUT_OF_MEMORY;                                    \
//     }                                                                       \
// ```
// NCBI's writer does not test whether an output file opened, a write succeeded or a
// rename took place, so makeblastdb has no message for these failures, and LOSAT's writer
// does not test them either (`ofstream`, `dir_entry`). An input file that cannot be opened
// or read, and an error of LMDB's own writes, LOSAT reports in the form of makeblastdb's
// "other" errors ("Error: " and exit status 255, `makeblastdb::MakedbError::exit_code`),
// with the path and the system's message; they are not errors of the database's contents
// ("BLAST Database creation error").
impl std::fmt::Display for WriteDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WriteDbError::Creation(message) => {
                write!(f, "BLAST Database creation error: {message}")
            }
            WriteDbError::CreationBytes(message) => write!(
                f,
                "BLAST Database creation error: {}",
                String::from_utf8_lossy(message)
            ),
            WriteDbError::Options(message) => write!(f, "BLAST options error: {message}"),
            WriteDbError::Database(message) => write!(f, "BLAST Database error: {message}"),
            // `what()` ends with a line feed and the message gets one more: the report ends
            // with an empty line.
            WriteDbError::Exception(e) => write!(f, "Error: {}", e.what()),
            WriteDbError::Other(message) => write!(f, "Error: {message}"),
            WriteDbError::Unsupported(message) => write!(f, "losat makeblastdb: {message}"),
            WriteDbError::OutOfMemory => f.write_str("BLAST ran out of memory"),
            WriteDbError::Io { path, source } => {
                write!(f, "Error: {}: {source}", path.display())
            }
        }
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
