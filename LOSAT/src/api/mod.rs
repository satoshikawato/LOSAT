//! The library API: a search of FASTA records in memory, with the command line's options and
//! the command line's output bytes.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/ (`CLocalBlast`, `CSearchResultSet`); the plan
//! is `docs/losat_library_surface_plan.md`.
//!
//! # Use
//!
//! ```no_run
//! use losat::api::{read_fasta, run_collect, validate, Inputs, Program, Role};
//!
//! # fn main() -> Result<(), losat::api::Error> {
//! // Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
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
//! # use losat::api::{read_inputs, run_collect, validate, Error, Inputs, Program};
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
//! # use losat::api::{run_batches, Error, Inputs, Search};
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
//!   HSP's row, not rounded, and each id the column's bytes ([`Id`]). Text to show is read
//!   from a report, never formatted here.
//! - An E-value is a value of the whole search, not of the pair: the same pair has another
//!   E-value among other subjects. [`QueryStats`] holds what fixed it, and `-dbsize` or
//!   `-searchsp` fix it across searches.
//! - A search writes nothing to the process's standard output or error and never ends the
//!   process; the command line's warnings come back as `Outcome::diagnostics`, its errors as
//!   [`Error`], whose message is made from the command line's text for the same failure by the
//!   rules that [`Error`] states (a search's NCBI message keeps its newlines, and an input's, a
//!   tool's and that of a start that ends the program lose their last ones; LOSAT's own refusals
//!   have no `Error: `; text that a search writes before the failure is not in it).
//! - `-num_threads` is an option like any other, checked as on the command line. A search of
//!   two or more threads builds its own thread pool for the call and cannot start from inside
//!   a Rayon worker; one of one thread can. Searches may run on several threads of a process at
//!   once; each has its own pool and its own outputs.
//! - The environment variables that the command line reads (`LOSAT_X_*` performance switches,
//!   which change no output; the diagnostic variables of `docs/agents/diagnostics.md`) are read
//!   here too. Each search call runs NCBI's application start at its beginning, as the command
//!   line runs it at its own: it reads the registry files (`.ncbirc`, the program's `.ini` or the
//!   file of `-conffile`) and NCBI's environment settings, which the search reads where the
//!   command line reads them: to find a database (`BLASTDB`), and where NCBI's library looks one
//!   up itself (blastn's `-window_masker_taxid`, the taxonomy database of the organism report of
//!   outfmt 18 and of the plural taxonomy names (`sscinames`, `scomnames`, `sblastnames`,
//!   `sskingdoms`), and the file of a `-matrix` that NCBI's tables do not have: tblastx, blastx
//!   `-ungapped`). Every search reads the query batch and archive settings (`BATCH_SIZE`,
//!   `ARCHIVE_FORMAT`, ...) as the command line does, of a database or of records.
//! - The start's messages about the registry files go where the command line writes them (the
//!   call's diagnostics stand for its standard error), `-logfile`'s file gets what the command line
//!   writes there (with tblastx, the messages and the errors of its run too), and a start that ends
//!   the command line's program is an [`Error::Search`] with its text. The registry and environment
//!   settings that the command line refuses after its start (a usage report setting that NCBI
//!   aborts on, an entry or a variable that LOSAT does not take) are the call's [`Error::Search`]
//!   then, with the command line's text after `Error: `, whatever the call's inputs. With `-dryrun`
//!   ([`Search::dry_run`]) a call ends after the start, the settings and the parse of the argv with
//!   the call's kind of subject (`-subject_loc` with a database fails as without `-dryrun`): no
//!   record or `-outfmt` value is read, no report is written. The versions and `-xmlhelp`, whose
//!   texts only the command line writes, are refused ([`validate`]); the tools print the versions.
//!   [`validate`] parses the argv without the start, which runs with each search: an argv that the
//!   parser refuses makes no `-logfile` file (the command line's start has made it before its
//!   parser's error), and one whose start ends the command line's program before its arguments are
//!   read (a missing `-conffile` file) and that the parser refuses too (two `-conffile`) gets the
//!   parser's error there. As any error of [`run_collect`] and [`run_batches`], one of the start
//!   carries its message alone: the start's messages before it (a usage report's) were in the
//!   call's diagnostics, which an error does not return; [`run()`] has written them to its writer.
//! - [`run_batches`] searches one of NCBI's query batches at a time: the host reports progress
//!   between batches and stops the search by not asking for the next one. Inside a batch only
//!   the host can stop it, by ending the thread or the process. A run that ends with exit status
//!   0 and no message in a batch (a database header record that ends before NCBI's read enters
//!   its elements) ends with that batch: it carries no HSPs and its diagnostics.
//! - The records are read before the search. [`read_inputs`] reads every query record or fails
//!   with the reader's error ([`Error::Input`]) before any report. [`read_input_records`] keeps
//!   the query records read before a reader error and the error ([`QueryEnd`]), and a search of
//!   them ([`InputRecords::inputs`], [`Inputs::query_end`]) ends with that error where the
//!   command line does: in the query batch that reads the record, after the reports of the
//!   batches before (and, with `-out`, after making the file, as [`run_to_file`] makes it where
//!   the command line makes `-out`). A subject's reader error fails the reading, as the command
//!   line reads every subject before the first query.
//! - Outfmt 13 and 14 write side files next to their report, so they come through
//!   [`run_to_file`] only; a report in memory has no place for them (the in-memory search under
//!   the API, `local_blast`, takes a receiver of them, `OutputSink::Files`, which LOSAT Web gives).
//! - [`Program::Psiblast`] is the search of LOSAT's `psiblast`, with its rounds
//!   (`-num_iterations`) and the formats it writes: the reports hold every round as the command
//!   line writes them, each HSP record its round ([`HspRecord::round`]) and each query the
//!   statistics of each round and whether its rounds converged ([`QueryRounds`],
//!   [`QueryHits::rounds`]). Its PSSM files (`-out_pssm`, `-out_ascii_pssm`, with
//!   `-save_pssm_after_last_round` and `-save_each_pssm`) are written where the command line
//!   writes them, with the argv's names: `-save_each_pssm` numbers them over the call, and over
//!   the batches of [`run_batches`] as over the command line's run; a call of [`run`](run()) or
//!   [`run_collect`] whose formats take two searches and that writes them is refused.
//! - A PSSM query of psiblast and tblastn (the command line's `-in_pssm`) is an input:
//!   [`Inputs::pssm`] and [`Inputs::pssm_database`] name the file, which the search reads; its one
//!   query is the PSSM's sequence, query 0 of the results, in one batch of [`run_batches`].
//! - A query alignment of psiblast (the command line's `-in_msa`) is an input the same way:
//!   [`Inputs::msa`] and [`Inputs::msa_database`] name the alignment file, which the search reads
//!   as the command line reads it, errors included; its one query is the alignment's master,
//!   query 0 of the results, in one batch of [`run_batches`]. `-msa_master_idx` and
//!   `-ignore_msa_master` are options of the argv ([`validate`] takes them; a search of other
//!   inputs fails with the parser's message, as the command line does without `-in_msa`).
//! - The API writes no search strategy file either: every program refuses
//!   `-export_search_strategy` ([`validate`]); the command line writes it.
//!   A search strategy file is an input of the API: [`Search::read_strategy`] reads its bytes (the
//!   command line's `-import_search_strategy FILE`, through the same import), [`Inputs::strategy`]
//!   gives the queries and the subjects (or the database) that it carries, and
//!   `-import_search_strategy` in an argv is refused ([`validate`]).
//! - [`run_makeblastdb`] and [`run_blastdbcmd`] are the two database programs, and
//!   [`run_dustmasker`] and [`run_segmasker`] the maskers of low-complexity regions (nucleotide,
//!   protein), [`run_windowmasker`] the masker of over-represented nucleotide regions,
//!   [`run_convert2blastmask`] the converter of lower case masks, [`run_makembindex`] the
//!   writer of megablast indices and [`run_blast_formatter`] the formatter of BLAST archives: the
//!   command line's words, with the caller's reader and writers for its standard streams.
//!
//! # Items
//!
//! - The program: [`Program`], [`describe`](describe()), [`Description`], [`Parameter`],
//!   [`ValueKind`]
//! - The inputs: [`read_input_file`], [`read_fasta`], [`read_inputs`], [`read_input_records`],
//!   [`Role`], [`FastaRecord`], [`Inputs`], [`InputRecords`], [`QueryEnd`], [`QueryError`]
//! - The options: [`validate`], [`Search`], [`API_OWNED`], [`Strategy`]
//! - The search: [`run_collect`], [`run`](run()), [`run_to_file`], [`run_batches`]
//! - The results: [`Outcome`], [`Batches`], [`Batch`], [`Report`], [`HspRecord`], [`Id`],
//!   [`QueryStats`], [`QueryRounds`], [`Strand`], [`QueryHits`], [`SubjectHits`]
//! - The tools: [`run_makeblastdb`], [`run_blastdbcmd`], [`run_dustmasker`], [`run_segmasker`],
//!   [`run_windowmasker`], [`run_convert2blastmask`], [`run_makembindex`],
//!   [`run_blast_formatter`]
//! - The failures: [`Error`]
//! - The engine's version: [`VERSION`]

mod batches;
mod describe;
mod error;
mod id;
pub(crate) mod input;
mod parts;
// The shared local search entry of the command line and the web ABIs: each program's `run_local`
// and the routing of its outputs. It speaks the engine's types, so it stays outside the API.
pub(crate) mod local_blast;
mod program;
mod run;
mod search;
pub(crate) mod start;
pub(crate) mod strategy;
mod tools;
mod view;

pub use crate::record::FastaRecord;
pub use batches::{run_batches, Batch, Batches};
pub use describe::{describe, Description, Parameter, ValueKind};
pub use error::Error;
pub use id::Id;
pub use input::{
    read_fasta, read_input_file, read_input_records, read_inputs, InputRecords, Inputs, QueryEnd,
    QueryError, Role,
};
pub use program::Program;
pub use run::{
    run, run_collect, run_to_file, HspRecord, Outcome, QueryRounds, QueryStats, Report, Strand,
};
pub(crate) use search::file_options;
pub use search::{validate, Search, API_OWNED};
pub use strategy::Strategy;
pub use tools::{
    run_blast_formatter, run_blastdbcmd, run_convert2blastmask, run_dustmasker, run_makeblastdb,
    run_makembindex, run_segmasker, run_windowmasker,
};
pub use view::{QueryHits, SubjectHits};

// NCBI reference (598d8ae6): c++/include/algo/blast/api/version.hpp:49-57
/// The version of the engine: what `losat --version` prints after the name.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
