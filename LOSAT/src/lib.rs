//! LOSAT: a pure-Rust reimplementation of NCBI BLAST+ for local searches, with the command
//! line's output byte for byte.
//!
//! The public surface of this crate is the module [`api`] (`docs/losat_library_surface_plan.md`):
//! a search of FASTA records in memory ([`api::read_fasta`], [`api::validate`],
//! [`api::run_collect`]), with the command line's options and the command line's reports.
//!
//! Every other module is the command line or the engine: NCBI's algorithms, ported function by
//! function. They are `pub(crate)`: their items change as the port advances. A few hidden items are
//! `pub` for the crate's own clients, outside the compatibility promise of the version number: the
//! `losat` binary calls the function `main_entry`; LOSAT Web's adapter names the engine through
//! `api` and the module `web_adapter`, reads and sets the field `FastaRecord::source` (what the
//! BLAST archive writes besides the letters) and calls the methods of `FastaRecord`'s hidden `impl`
//! (`believed`, `shown_id`, `subject_tabular_id` and the others there); the tests use the module
//! `wording`.

// The NCBI BLAST+ release that LOSAT's output follows. This is the one place in the code that
// names the release: `-version` says "output compatible with blast 2.17.0", and `crate::wording`
// maps NCBI's version statement (`BLASTN 2.17.0+` in the freezes) to LOSAT's (`BLASTN (LOSAT
// 0.2.0)`, the Owner 2026-11-21 「Q1=B」). NCBI's frozen `-help` texts of every program that has
// one (the search programs' `help.txt`, `blastdbcmd/help.txt`, `writedb/makeblastdb_help.txt`)
// and its `-xmlhelp` files carry it too and are printed through that mapping.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:81-83
// NCBI reference (598d8ae6): c++/include/algo/blast/api/version.hpp:53-62
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:60-64 (report headers)
// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:732-739 (`-version`)
//
// `ncbi_blast_release!()`: NCBI's package version (`GetPackageVersion().Print()`, "2.17.0"),
// which LOSAT's `-version` names after `output compatible with blast` (LOSAT's own `Package:`
// version is `CARGO_PKG_VERSION`). A macro, so that `concat!` builds constants with it.
macro_rules! ncbi_blast_release {
    () => {
        "2.17.0"
    };
}

/// NCBI's package version (`GetPackageVersion().Print()`, "2.17.0"), which LOSAT's `-version`
/// names after `output compatible with blast` (LOSAT's own `Package:` version is
/// `CARGO_PKG_VERSION`).
pub(crate) const NCBI_BLAST_RELEASE: &str = ncbi_blast_release!();
/// What NCBI's `CBlastVersion::Print()` writes ("2.17.0+"): the version in NCBI's report headers,
/// XML and JSON version, archive ident, SAM's `@PG VN` and `-version`, where LOSAT writes its own
/// (`crate::wording::LOSAT_VERSION`, `REPORT_VERSION`); the mapping of NCBI's freezes reads it.
pub(crate) const NCBI_BLAST_VERSION: &str = concat!(ncbi_blast_release!(), "+");

pub(crate) mod algorithm;
pub(crate) mod common;
pub(crate) mod sequence;
pub(crate) mod utils;

pub(crate) mod config;
pub(crate) mod report;
pub(crate) mod stats;

// NCBI api/, blastinput/ and core/ layers.
pub mod api;
pub(crate) mod blastinput;
pub(crate) mod core;
// NCBI objtools/blast/seqdb_reader (CSeqDB): reader of BLAST database volumes.
pub(crate) mod seqdb;
// NCBI objects/seqloc (CSeq_id), objects/blastdb (CBlast_def_line_set) and the id helpers
// of objtools/align_format: the Seq-ids of database sequences.
pub(crate) mod seqid;

/// The mask readers and writers of NCBI's masking tools (`dustmasker`, `segmasker`,
/// `windowmasker`, `convert2blastmask`); makembindex reads its sequences through them too.
pub(crate) mod seqmasks_io;
// NCBI objtools/blast/seqdb_writer (CWriteDB) and c++/src/app/blastdb/makeblastdb.cpp:
// writer of BLAST database volumes.
pub(crate) mod writedb;
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp: the BLAST database client (its -info).
pub(crate) mod blastdbcmd;
// NCBI reference (598d8ae6): c++/src/app/dustmasker/dust_mask_app.cpp: the masker of low-complexity nucleotide regions.
pub(crate) mod dustmasker;
// NCBI reference (598d8ae6): c++/src/app/segmasker/segmasker.cpp: the masker of low-complexity protein regions.
pub(crate) mod segmasker;
// NCBI reference (598d8ae6): c++/src/app/winmasker/win_mask_app.cpp: the masker of over-represented nucleotide
// regions by unit counts (WindowMasker).
pub(crate) mod windowmasker;
// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp: the masking information of lower case FASTA input.
pub(crate) mod convert2blastmask;
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex: the megablast database index (CDbIndex) and its sequence streams.
pub(crate) mod dbindex;
// NCBI reference (598d8ae6): c++/src/algo/blast/dbindex/makeindex/mkindex_app.cpp: the writer of megablast indices.
pub(crate) mod makembindex;
// NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp: the writer of RPS databases of profiles.
pub(crate) mod makeprofiledb;
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdb_aliastool.cpp: the writer of alias files and of binary
// GI and seqid lists.
pub(crate) mod blastdb_aliastool;
// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp and c++/src/corelib/ncbiargs.cpp: the application and argument reader of
// blastdbcmd, dustmasker, segmasker, windowmasker, convert2blastmask, makembindex, makeprofiledb and
// blastdb_aliastool.
pub(crate) mod ncbi_app;
// NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp: the reports of NCBI's exceptions (`CException::what()`).
pub(crate) mod ncbi_exception;
// LOSAT's own words in place of NCBI's names and toolkit internals in messages, and the one
// mapping of NCBI's frozen messages to LOSAT's (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 12):
// hidden, outside the public API, for the tests and LOSAT Web's adapter.
#[doc(hidden)]
pub mod wording;
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp: the number conversions of `NStr`.
pub(crate) mod ncbi_str;
// NCBI reference (598d8ae6): c++/src/corelib/version.cpp: what `-version` and the full versions print (`CVersionAPI`).
pub(crate) mod ncbi_version;
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp: what `-xmlhelp` prints (`CArgDescriptions::PrintUsageXml`).
pub(crate) mod ncbi_xmlhelp;
// NCBI's rpsbproc 0.5.1 (RpsbProc/src/rpsbproc_main.cpp of the 2025-05-13 package): the CDD
// domain annotation of rpsblast's archives.
pub(crate) mod rpsbproc;

#[cfg(target_arch = "wasm32")]
pub(crate) mod web_api;
// The words of a request of LOSAT Web's ABI v1 (frozen), read as v1 reads them.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod web_v1_argv;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
pub(crate) mod cli;
// The body of the `losat` binary (`src/main.rs`): the one item of this crate that the binary names.
#[doc(hidden)]
pub use cli::main_entry::main_entry;

// The type of `api::FastaRecord`: private, so that rustdoc inlines its re-exports (record.rs).
mod record;

// The engine items that LOSAT Web's adapter (`web/adapter`) names: one hidden façade, outside the
// public API (web_adapter.rs).
#[doc(hidden)]
pub mod web_adapter;
