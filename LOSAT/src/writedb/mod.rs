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
//! Not written: parsed Seq-ids and their ISAM files (`-parse_seqids`), mask
//! data, taxonomy data, the LMDB files of version 5, and the JSON metadata
//! file (`.pjs` / `.njs`).

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
    /// The operating system failed to create, write, rename or remove a
    /// file.
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
}

impl std::fmt::Display for WriteDbError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WriteDbError::Creation(message) => {
                write!(f, "BLAST Database creation error: {message}")
            }
            WriteDbError::Options(message) => write!(f, "BLAST options error: {message}"),
            WriteDbError::Other(message) => write!(f, "Error: {message}"),
            WriteDbError::Io { path, source } => {
                write!(
                    f,
                    "BLAST Database creation error: {}: {source}",
                    path.display()
                )
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
