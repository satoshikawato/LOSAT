//! LOSAT: a pure-Rust reimplementation of NCBI BLAST+ for local searches, with the command
//! line's output byte for byte.
//!
//! The public surface of this crate is the module [`api`] (`docs/losat_library_surface_plan.md`):
//! a search of FASTA records in memory ([`api::read_fasta`], [`api::validate`],
//! [`api::run_collect`]), with the command line's options and the command line's reports.
//!
//! Every other module is the command line or the engine: NCBI's algorithms, ported function by
//! function. The modules that the command line, the web adapters or the tests name are `pub` but
//! hidden from this documentation and outside the compatibility promise of the version number:
//! their items change as the port advances. The others are `pub(crate)`.

#[doc(hidden)]
pub mod algorithm;
#[doc(hidden)]
pub mod common;
pub(crate) mod sequence;
#[doc(hidden)]
pub mod utils;

#[doc(hidden)]
pub mod config;
#[doc(hidden)]
pub mod report;
#[doc(hidden)]
pub mod stats;

// NCBI api/, blastinput/ and core/ layers.
pub mod api;
#[doc(hidden)]
pub mod blastinput;
#[doc(hidden)]
pub mod core;
// NCBI objtools/blast/seqdb_reader (CSeqDB): reader of BLAST database volumes.
#[doc(hidden)]
pub mod seqdb;
// NCBI objects/seqloc (CSeq_id), objects/blastdb (CBlast_def_line_set) and the id helpers
// of objtools/align_format: the Seq-ids of database sequences.
pub(crate) mod seqid;
// NCBI objtools/blast/seqdb_writer (CWriteDB) and app/blastdb/makeblastdb.cpp:
// writer of BLAST database volumes.
#[doc(hidden)]
pub mod writedb;
// NCBI app/blastdb/blastdbcmd.cpp: the BLAST database client (its -info).
#[doc(hidden)]
pub mod blastdbcmd;

#[cfg(target_arch = "wasm32")]
pub(crate) mod web_api;

// NCBI blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
#[doc(hidden)]
pub mod cli;

// The type of `api::FastaRecord`: private, so that rustdoc inlines its re-exports (record.rs).
mod record;
