#![allow(warnings, clippy::all)]

pub mod algorithm;
pub mod common;
pub mod sequence;
pub mod utils;

pub mod align;
pub mod config;
pub mod post;
pub mod report;
pub mod seed;
pub mod stats;

// NCBI api/, blastinput/ and core/ layers.
pub mod api;
pub mod blastinput;
pub mod core;
// NCBI objtools/blast/seqdb_reader (CSeqDB): reader of BLAST database volumes.
pub mod seqdb;
// NCBI objects/seqloc (CSeq_id), objects/blastdb (CBlast_def_line_set) and the id helpers
// of objtools/align_format: the Seq-ids of database sequences.
pub mod seqid;
// NCBI objtools/blast/seqdb_writer (CWriteDB) and app/blastdb/makeblastdb.cpp:
// writer of BLAST database volumes.
pub mod writedb;
// NCBI app/blastdb/blastdbcmd.cpp: the BLAST database client (its -info).
pub mod blastdbcmd;

#[cfg(target_arch = "wasm32")]
pub mod web_api;

// NCBI blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
pub mod cli;
