//! LOSAT: a pure-Rust reimplementation of NCBI BLAST+ for local searches, with the command
//! line's output byte for byte.
//!
//! The public surface of this crate is the module [`api`] (`docs/losat_library_surface_plan.md`):
//! a search of FASTA records in memory ([`api::read_fasta`], [`api::validate`],
//! [`api::run_collect`]), with the command line's options and the command line's reports.
//!
//! Every other module is the command line or the engine: NCBI's algorithms, ported function by
//! function. They are `pub(crate)`: their items change as the port advances. Two hidden items are
//! `pub` for the crate's own two clients, outside the compatibility promise of the version number:
//! the `losat` binary calls the function `main_entry`, and LOSAT Web's adapter names the engine
//! through `api` and the module `web_adapter`.

// The NCBI BLAST+ release that LOSAT's output follows, and the spellings of it that the output
// carries. This is the one place in the code that names the release; every report, archive and
// `-version` text takes it from here. `blastdbcmd/help.txt` is NCBI's frozen `-help` text and
// carries it too; the test `the_frozen_help_names_the_version_that_losat_writes` ties the two.
//
// The macros produce literals so that `concat!` can build `&'static str` constants such as
// "BLASTX 2.17.0+" (`crate::ncbi_blast_version!()`). Other modules call them by path, through the
// `pub(crate) use` lines below, which therefore have to follow the definitions.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:81-83
// ```c
// NCBI_XBLAST_EXPORT const int   kBlastMajorVersion = 2;
// NCBI_XBLAST_EXPORT const int   kBlastMinorVersion = 17;
// NCBI_XBLAST_EXPORT const int   kBlastPatchVersion = 0;
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/version.hpp:53-62
// ```cpp
// class CBlastVersion : public CVersionInfo {
// public:
//     CBlastVersion()
//         : CVersionInfo(kBlastMajorVersion,
//                        kBlastMinorVersion,
//                        kBlastPatchVersion) {}
//     virtual string Print(void) const {
//         return CVersionInfo::Print() + "+";
//     }
// };
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:60-64 (report headers)
// ```c++
// string CBlastFormatUtil::BlastGetVersion(const string program)
// {
//     string program_uc = program;
//     return NStr::ToUpper(program_uc) + " " + blast::CBlastVersion().Print();
// }
// ```
// NCBI reference (598d8ae6): c++/src/corelib/version.cpp:732-739 (`-version`)
// ```cpp
// if (flags & fVersionInfo) {
//     os << appname << ": " << m_VersionInfo->Print() << endl;
// }
// ...
//     os << " Package: " << GetPackageName() << ' '
//        << GetPackageVersion().Print() << ", build "
// ```
//
// `ncbi_blast_release!()`: NCBI's package version (`GetPackageVersion().Print()`, "2.17.0"),
// which LOSAT's `-version` names after `output compatible with blast` (LOSAT's own `Package:`
// version is `CARGO_PKG_VERSION`).
macro_rules! ncbi_blast_release {
    () => {
        "2.17.0"
    };
}
// `ncbi_blast_version!()`: what `CBlastVersion::Print()` writes ("2.17.0+"), used in the report
// headers ("# BLASTN 2.17.0+"), the XML and JSON version, the archive ident, SAM's `@PG VN` and
// the first line of `-version`.
macro_rules! ncbi_blast_version {
    () => {
        concat!($crate::ncbi_blast_release!(), "+")
    };
}
pub(crate) use ncbi_blast_release;
pub(crate) use ncbi_blast_version;

/// NCBI's package version (`GetPackageVersion().Print()`, "2.17.0"), which LOSAT's `-version`
/// names after `output compatible with blast` (LOSAT's own `Package:` version is
/// `CARGO_PKG_VERSION`).
pub(crate) const NCBI_BLAST_RELEASE: &str = ncbi_blast_release!();
/// What NCBI's `CBlastVersion::Print()` writes: the version in report headers, the XML and JSON
/// version, the archive ident, SAM's `@PG VN` and the first line of `-version`.
pub(crate) const NCBI_BLAST_VERSION: &str = ncbi_blast_version!();

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

/// The mask readers and writers of NCBI's masking tools (`dustmasker`, `segmasker`).
pub(crate) mod seqmasks_io;
// NCBI objtools/blast/seqdb_writer (CWriteDB) and app/blastdb/makeblastdb.cpp:
// writer of BLAST database volumes.
pub(crate) mod writedb;
// NCBI app/blastdb/blastdbcmd.cpp: the BLAST database client (its -info).
pub(crate) mod blastdbcmd;
// NCBI app/dustmasker/dust_mask_app.cpp: the masker of low-complexity nucleotide regions.
pub(crate) mod dustmasker;
// NCBI app/segmasker/segmasker.cpp: the masker of low-complexity protein regions.
pub(crate) mod segmasker;
// NCBI corelib/ncbiapp.cpp and corelib/ncbiargs.cpp: the application and argument reader of
// blastdbcmd, dustmasker and segmasker.
pub(crate) mod ncbi_app;
// NCBI corelib/ncbiexpt.cpp: the reports of NCBI's exceptions (`CException::what()`).
pub(crate) mod ncbi_exception;
// NCBI corelib/ncbistr.cpp: the number conversions of `NStr`.
pub(crate) mod ncbi_str;
// NCBI corelib/version.cpp: what `-version` and the full versions print (`CVersionAPI`).
pub(crate) mod ncbi_version;
// NCBI's rpsbproc 0.5.1 (RpsbProc/src/rpsbproc_main.cpp of the 2025-05-13 package): the CDD
// domain annotation of rpsblast's archives.
pub(crate) mod rpsbproc;

#[cfg(target_arch = "wasm32")]
pub(crate) mod web_api;

// NCBI blastinput/blast_args.cpp:332-349: shared string-valued filtering arguments.
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
