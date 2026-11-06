//! The library API: a search of FASTA records in memory, with the command line's options and
//! the command line's output bytes.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/api/ (`CLocalBlast`, `CSearchResultSet`); the plan
//! is `docs/losat_library_surface_plan.md`.
//!
//! # Use
//!
//! ```no_run
//! use LOSAT::api::{read_fasta, run_collect, validate, Inputs, Program, Role};
//!
//! # fn main() -> Result<(), LOSAT::api::Error> {
//! let query = read_fasta(Program::Blastp, Role::Query, b">q\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n")?;
//! let subject = read_fasta(Program::Blastp, Role::Subject, b">s\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n")?;
//! let search = validate(Program::Blastp, &["-evalue", "1e-5", "-num_threads", "2"])?;
//! let outcome = run_collect(&search, Inputs::new(&query, &subject), &["6", "0"])?;
//! // outcome.reports[0].bytes: the outfmt 6 rows; outcome.hsps: one record per HSP;
//! // outcome.query_stats: the search space and Karlin-Altschul values of each query.
//! for query in outcome.queries() {
//!     for subject in query.subjects() {
//!         println!("{} {}: {} HSPs", query.q_idx, subject.s_idx, subject.hsps.len());
//!     }
//! }
//! # Ok(()) }
//! ```
//!
//! A BLAST database in place of the subject records (the command line's `-db`, found on the
//! `BLASTDB` search path):
//!
//! ```no_run
//! # use LOSAT::api::{read_inputs, run_collect, validate, Error, Inputs, Program};
//! # fn use_it(fasta: &[u8]) -> Result<(), Error> {
//! let search = validate(Program::Blastn, &["-taxids", "9606", "-evalue", "1e-5"])?;
//! let (query, _) = read_inputs(&search, fasta, b"")?;
//! let outcome = run_collect(&search, Inputs::database(&query, "core_nt"), &["6"])?;
//! # Ok(()) }
//! ```
//!
//! Many queries, one of NCBI's query batches in memory at a time:
//!
//! ```no_run
//! # use LOSAT::api::{run_batches, Error, Inputs, Search};
//! # fn use_it(search: &Search, inputs: Inputs<'_>) -> Result<(), Error> {
//! let batches = run_batches(search, inputs)?;
//! let total = batches.query_count();
//! for batch in batches {
//!     let batch = batch?;
//!     println!("{} of {total} queries, {} HSPs", batch.query_indices.end, batch.hsps.len());
//!     // Leaving the loop stops the search: no further batch is searched.
//! }
//! # Ok(()) }
//! ```
//!
//! # Contract
//!
//! - The options are one argv, parsed by the command line's parser ([`validate`]). The inputs
//!   (`-query`, `-subject`, `-db`), the output file and the output format are not options here:
//!   the records come from [`read_fasta`] in an [`Inputs`], a database by its name
//!   ([`Inputs::database`]), the formats are arguments of [`run_collect`] and [`run`](run()).
//!   [`validate`] holds an argv for either kind of subject; the search of the inputs fails with
//!   the parser's message on an option of the other kind (`-subject_loc` excludes a database,
//!   `-taxids` needs one).
//! - Each report of [`run_collect`] is byte for byte what the command line writes to `-out` for
//!   that `-outfmt` alone, with the inputs' names where the command line shows its file names.
//!   [`Report`] tells where each HSP lies in its bytes. The HSP records ([`HspRecord`]) are the
//!   engine's values: each number that names a `-outfmt` column is the column's value for the
//!   HSP's row, not rounded. Text to show is read from a report, never formatted here.
//! - An E-value is a value of the whole search, not of the pair: the same pair has another
//!   E-value among other subjects. [`QueryStats`] holds what fixed it, and `-dbsize` or
//!   `-searchsp` fix it across searches.
//! - A search writes nothing to the process's standard output or error and never ends the
//!   process; the command line's warnings come back as `Outcome::diagnostics`, its errors as
//!   [`Error`] with the command line's message.
//! - `-num_threads` is an option like any other, checked as on the command line. A search of
//!   two or more threads builds its own thread pool for the call and cannot start from inside
//!   a Rayon worker; one of one thread can. Searches may run on several threads of a process at
//!   once; each has its own pool and its own outputs.
//! - The environment variables that the command line reads (`LOSAT_X_*` performance switches,
//!   which change no output; the diagnostic variables of `docs/agents/diagnostics.md`) are read
//!   here too. A database search reads what the command line reads to find a database (`BLASTDB`,
//!   `.ncbirc`) and the query batch and archive settings (`BATCH_SIZE`, `ARCHIVE_FORMAT`, ...);
//!   a search of records reads no NCBI setting.
//! - [`run_batches`] searches one of NCBI's query batches at a time: the host reports progress
//!   between batches and stops the search by not asking for the next one. Inside a batch only
//!   the host can stop it, by ending the thread or the process.
//! - The records are read before the search, so a reader error is an [`Error::Input`] before
//!   any report, where the command line has reported the batches before the error.
//! - Outfmt 13 and 14 write side files next to their report, so they come through
//!   [`run_to_file`] only; a report in memory has no place for them (the in-memory search under
//!   the API, `local_blast`, takes a receiver of them, `OutputSink::Files`, which LOSAT Web gives).
//! - [`Program::Psiblast`] is the search of LOSAT's `psiblast`: its first round, and the formats
//!   it writes. The API writes no PSSM file, so it refuses `-out_ascii_pssm`,
//!   `-save_pssm_after_last_round` and `-save_each_pssm` ([`validate`]).
//! - [`run_makeblastdb`] and [`run_blastdbcmd`] are the two database programs: the command line's
//!   words, with the caller's writers for its standard streams.
//!
//! # Items
//!
//! - The program: [`Program`], [`describe`](describe()), [`Description`], [`Parameter`],
//!   [`ValueKind`]
//! - The inputs: [`read_input_file`], [`read_fasta`], [`read_inputs`], [`Role`], [`FastaRecord`],
//!   [`Inputs`]
//! - The options: [`validate`], [`Search`], [`API_OWNED`]
//! - The search: [`run_collect`], [`run`](run()), [`run_to_file`], [`run_batches`]
//! - The results: [`Outcome`], [`Batches`], [`Batch`], [`Report`], [`HspRecord`], [`QueryStats`],
//!   [`Strand`], [`QueryHits`], [`SubjectHits`]
//! - The database programs: [`run_makeblastdb`], [`run_blastdbcmd`]
//! - The failures: [`Error`]
//! - The engine's version: [`VERSION`]

mod batches;
mod describe;
mod error;
mod input;
// The shared local search entry of the command line and the web ABIs: each program's `run_local`
// and the routing of its outputs. It speaks the engine's types, so it stays outside the API.
#[doc(hidden)]
pub mod local_blast;
mod program;
mod run;
mod search;
mod tools;
mod view;

pub use crate::blastinput::fasta_reader::FastaRecord;
pub use batches::{run_batches, Batch, Batches};
pub use describe::{describe, Description, Parameter, ValueKind};
pub use error::Error;
pub use input::{read_fasta, read_input_file, read_inputs, Inputs, Role};
pub use program::Program;
pub use run::{run, run_collect, run_to_file, HspRecord, Outcome, QueryStats, Report, Strand};
pub use search::{validate, Search, API_OWNED};
pub use tools::{run_blastdbcmd, run_makeblastdb};
pub use view::{QueryHits, SubjectHits};

// NCBI reference (598d8ae6): c++/include/algo/blast/api/version.hpp:49-57
// ```c++
// /// Keeps track of the version of the BLAST engine in the NCBI C++ toolkit.
// /// Used to perform run-time version checks
// ///
// /// For reference, please refer to http://apr.apache.org/versioning.html
// class CBlastVersion : public CVersionInfo {
// public:
//     CBlastVersion()
//         : CVersionInfo(kBlastMajorVersion,
//                        kBlastMinorVersion,
// ```
/// The version of the engine: what `losat --version` prints after the name.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
