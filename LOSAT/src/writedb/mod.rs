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
pub mod fasta_source;
