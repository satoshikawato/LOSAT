//! Writer of NCBI BLAST databases (the version 4 volume layout that
//! `makeblastdb -blastdb_version 4` writes) and the `LOSAT makedb` command.
//!
//! NCBI's writer is `c++/src/objtools/blast/seqdb_writer/` (CWriteDB,
//! CBuildDatabase) under the application `c++/src/app/blastdb/makeblastdb.cpp`.
//! The files written are the ones `crate::seqdb` reads: the index file
//! (`.pin` / `.nin`), the sequence file (`.psq` / `.nsq`), the header file
//! (`.phr` / `.nhr`), and an alias file (`.pal` / `.nal`) when the database
//! has more than one volume.
//!
//! Not written: parsed Seq-ids and their ISAM files (`-parse_seqids`),
//! taxonomy data, the LMDB files of version 5, and the JSON metadata file
//! (`.pjs` / `.njs`). Mask data (`-mask_data`: the mask column files `.?aa`,
//! `.?ab`, `.?ac`) can be written by the library (`mask_data`), but
//! `LOSAT makedb` does not offer the option yet.

pub mod bioseq;
pub mod build_db;
pub mod convert;
pub mod data_size;
pub mod fasta_source;
pub mod files;
pub mod format_guess;
pub mod header;
pub mod local_offset;
pub mod makeblastdb;
pub mod mask_data;
pub mod mask_info;
pub mod time;
pub mod volume;
pub mod writer;

#[cfg(test)]
mod tests;

use std::path::PathBuf;

pub use writer::WriteDb;

/// Error of the database writer.
#[derive(Debug)]
pub enum WriteDbError {
    /// NCBI's `CWriteDBException`: the text makeblastdb prints after
    /// "BLAST Database creation error: ".
    Creation(String),
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
