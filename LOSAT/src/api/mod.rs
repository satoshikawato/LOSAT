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
//! # Ok(()) }
//! ```
//!
//! # Contract
//!
//! - The options are one argv, parsed by the command line's parser ([`validate`]). The inputs
//!   (`-query`, `-subject`), the output file and the output format are not options here: the
//!   records come from [`read_fasta`] in an [`Inputs`], the formats are arguments of
//!   [`run_collect`] and [`run`](run()).
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
//!   here too.
//! - Not in this API: a database subject (`-db`), `psiblast`, `makeblastdb`, `blastdbcmd`,
//!   stopping a search (the host stops the thread or process).
//!
//! # Items
//!
//! - The program: [`Program`], [`describe`](describe()), [`Description`], [`Parameter`]
//! - The inputs: [`read_fasta`], [`Role`], [`FastaRecord`], [`Inputs`]
//! - The options: [`validate`], [`Search`], [`API_OWNED`]
//! - The search: [`run_collect`], [`run`](run()), [`Outcome`], [`Report`], [`HspRecord`],
//!   [`QueryStats`], [`Strand`]
//! - The failures: [`Error`]

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

pub use crate::blastinput::fasta_reader::FastaRecord;
pub use describe::{describe, Description, Parameter};
pub use error::Error;
pub use input::{read_fasta, Inputs, Role};
pub use program::Program;
pub use run::{run, run_collect, HspRecord, Outcome, QueryStats, Report, Strand};
pub use search::{validate, Search, API_OWNED};
