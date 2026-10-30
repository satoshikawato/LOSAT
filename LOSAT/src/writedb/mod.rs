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
//! With `-parse_seqids` (library only until the version 5 keys are written,
//! `MakedbArgs::parse_seqids`) and `-hash_index`, the volumes also get their ISAM
//! files (`.psi`/`.psd` of the Seq-ids in version 4, `.pni`/`.pnd` of the GIs,
//! `.phi`/`.phd` of the sequence hashes), their OID to GI file (`.pog`) and the
//! duplicate check of their Seq-ids (`volume`, `isam`, `id_set`); `-taxid` and
//! `-taxid_map` (library only, it requires `-parse_seqids`) put tax ids on the
//! definition lines and into the version 5 tax id files (`taxid_set`).
//!
//! Not written: the version 5 acc2oid keys and `.pos` of `-parse_seqids`, and
//! mask data.

pub mod bioseq;
pub mod build_db;
pub mod convert;
pub mod data_size;
pub mod deflines;
pub mod fasta_source;
pub mod files;
pub mod format_guess;
pub mod header;
pub(crate) mod id_set;
pub(crate) mod isam;
// The LMDB files of version 5. Parts that only `-parse_seqids` and the tax id options
// reach are not called yet.
#[allow(dead_code)]
pub(crate) mod lmdb;
pub mod local_offset;
pub mod makeblastdb;
#[allow(dead_code)]
pub(crate) mod mdb;
pub(crate) mod metadata;
pub mod platform;
pub mod taxid_set;
pub mod time;
pub mod volume;
pub mod writer;

#[cfg(test)]
mod tests;

use std::path::PathBuf;

pub use platform::Platform;
pub use writer::WriteDb;

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
    /// Another NCBI exception that ends makeblastdb: the text it prints after
    /// "Error: ".
    Other(String),
    /// An allocation failed (`std::bad_alloc`): makeblastdb prints "BLAST ran
    /// out of memory".
    OutOfMemory,
    /// The operating system failed to create or write a file of the
    /// database, or to open the input file.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// The operating system failed to give a file of a volume its final
    /// name.
    Rename {
        from: PathBuf,
        to: PathBuf,
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
// NCBI's writer does not test whether a file opened, a write succeeded or a
// rename took place, so makeblastdb has no message for these failures. LOSAT
// reports them in the form of makeblastdb's "other" errors ("Error: " and
// exit status 255, `makeblastdb::MakedbError::exit_code`), with the path and
// the system's message; they are not errors of the database's contents
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
            WriteDbError::Other(message) => write!(f, "Error: {message}"),
            WriteDbError::OutOfMemory => f.write_str("BLAST ran out of memory"),
            WriteDbError::Io { path, source } => {
                write!(f, "Error: {}: {source}", path.display())
            }
            WriteDbError::Rename { from, to, source } => {
                write!(f, "Error: {} -> {}: {source}", from.display(), to.display())
            }
        }
    }
}

impl std::error::Error for WriteDbError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            WriteDbError::Io { source, .. } | WriteDbError::Rename { source, .. } => Some(source),
            _ => None,
        }
    }
}
