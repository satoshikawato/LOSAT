//! A search strategy file as an input of a search (`-import_search_strategy` of the command line,
//! line RG item F2-4).
//!
//! The command line reads the file named by `-import_search_strategy` in `RecoverSearchStrategy`
//! and the engine's import (`import::apply`) turns it into the arguments of the strategy's
//! program, with the queries and the subjects (or the database's name) that the strategy
//! carries. The library API and LOSAT Web ABI v2 give the strategy as bytes and come through the
//! same import: one engine path (`apply`), then the program's `run_local` over records, as for
//! any other search.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:606-628
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3837
#![warn(missing_docs, unused)]

use std::path::Path;

use super::error::Error;
use super::input::read_search_inputs;
use super::program::Program;
use super::search::Search;
use crate::cli::Commands;
use crate::record::FastaRecord;

/// What a strategy file gives a search: the arguments as the import made them, and the queries
/// and subjects that the strategy carries as FASTA bytes.
#[derive(Debug)]
pub struct Applied {
    /// The arguments of the strategy's program with the strategy's options set (the engine's
    /// import): the search's program's, or another's where the search's program runs another
    /// program's strategy as the command line does.
    pub command: Commands,
    /// The strategy's queries (NCBI's temporary FASTA file), empty when it has none.
    pub queries: Vec<u8>,
    /// The strategy's subject sequences as FASTA text, empty for a database or none.
    pub subjects: Vec<u8>,
    /// The strategy's database name, when its subject is one.
    pub database: Option<String>,
}

/// The file that the argv's `-import_search_strategy` names, as the parser kept it.
pub fn import_name(command: &Commands) -> Option<&Path> {
    match command {
        Commands::Blastp(args) => args.import_search_strategy.as_deref(),
        Commands::Tblastn(args) => args.import_search_strategy.as_deref(),
        Commands::Blastn(args) => args.import_search_strategy.as_deref(),
        Commands::Tblastx(args) => args.import_search_strategy.as_deref(),
        Commands::Blastx(args) => args.import_search_strategy.as_deref(),
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => None,
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:606-628
/// `RecoverSearchStrategy` of the arguments `command` (parsed with an `-import_search_strategy`
/// of any name) over the strategy `bytes`: the engine's import (`import::apply`), the one the
/// command line runs at the start of its run. The options the command line may give besides
/// (`-outfmt`, `-num_threads`, `-query_loc` ...) are those that the strategy's import takes from
/// them, as there. A strategy of another program is that program's arguments where the command
/// line runs it so (blastn and blastx running a TBLASTX strategy, tblastx a BLASTX or BLASTN one).
/// The strategy's queries and subjects move out of the arguments into `Applied`: the search
/// reads them as records and the arguments keep no copy.
///
/// Errors are the command line's messages (after the warnings it posts before them): a file
/// that is not a strategy is NCBI's error that the strategy cannot be read, and a strategy of
/// another program that the command line does not run or with a PSSM query is refused as there.
pub fn apply(command: Commands, bytes: &[u8]) -> Result<Applied, Error> {
    use crate::report::archive::read::import;
    #[cfg(test)]
    tests::APPLIED.with(|count| count.set(count.get() + 1));
    if !import::imports(&command) {
        return Err(Error::Input(
            b"-import_search_strategy is not supported by this program of LOSAT".to_vec(),
        ));
    }
    // The warnings that the command line posts before a failure (to its standard error) come
    // first in the error's text, as the two come on a terminal.
    let mut posted = Vec::new();
    let command = import::apply(command, Some(bytes), &mut posted).map_err(|error| {
        posted.extend_from_slice(format!("{error:#}").as_bytes());
        Error::Input(posted)
    })?;
    let text = |held: Option<&std::sync::Arc<import::StrategyText>>| {
        held.map_or(Ok(Vec::new()), |text| text.take())
            .map_err(|error| Error::Input(format!("{error:#}").into_bytes()))
    };
    let run = import_run(&command);
    let queries = text(run.and_then(|run| run.queries.as_ref()))?;
    let subjects = text(run.and_then(|run| run.subjects.as_ref()))?;
    let database = database_of(&command);
    Ok(Applied {
        command,
        queries,
        subjects,
        database,
    })
}

/// The program whose arguments `command` holds: after an import, the strategy's program
/// (`Applied::command`).
pub fn program_of(command: &Commands) -> Option<Program> {
    Some(match command {
        Commands::Blastp(_) => Program::Blastp,
        Commands::Tblastn(_) => Program::Tblastn,
        Commands::Blastn(_) => Program::Blastn,
        Commands::Tblastx(_) => Program::Tblastx,
        Commands::Blastx(_) => Program::Blastx,
        _ => return None,
    })
}

/// The import run state of the applied arguments (the files and the warnings of the strategy).
fn import_run(command: &Commands) -> Option<&crate::report::archive::read::import::ImportRun> {
    match command {
        Commands::Blastp(args) => args.import.as_ref(),
        Commands::Tblastn(args) => args.import.as_ref(),
        Commands::Blastn(args) => args.import.as_ref(),
        Commands::Tblastx(args) => args.import.as_ref(),
        Commands::Blastx(args) => args.import.as_ref(),
        _ => None,
    }
}

fn import_run_mut(
    command: &mut Commands,
) -> Option<&mut crate::report::archive::read::import::ImportRun> {
    match command {
        Commands::Blastp(args) => args.import.as_mut(),
        Commands::Tblastn(args) => args.import.as_mut(),
        Commands::Blastn(args) => args.import.as_mut(),
        Commands::Tblastx(args) => args.import.as_mut(),
        Commands::Blastx(args) => args.import.as_mut(),
        _ => None,
    }
}

/// The query range that the strategy gives (`-query_loc` of the command line's import comes from
/// the strategy's `QueryRange`, `Start::run`): the search's range where the argv has none.
pub(crate) fn query_range(
    command: &Commands,
) -> Option<crate::blastinput::seq_range::SequenceRange> {
    import_run(command).and_then(|run| run.query_range)
}

/// The database that the applied arguments search (`-db`), when the strategy's subject is one.
fn database_of(command: &Commands) -> Option<String> {
    let db = match command {
        Commands::Blastp(args) => args.db.as_deref(),
        Commands::Tblastn(args) => args.db.as_deref(),
        Commands::Blastn(args) => args.db.as_deref(),
        Commands::Tblastx(args) => args.db.as_deref(),
        Commands::Blastx(args) => args.db.as_deref(),
        _ => None,
    };
    db.map(|name| name.to_string_lossy().into_owned())
}

/// Sets the arguments' `-out` (the file of `run_to_file`).
fn set_out(command: &mut Commands, path: &Path) {
    let path = Some(path.to_path_buf());
    match command {
        Commands::Blastp(args) => args.out = path,
        Commands::Tblastn(args) => args.out = path,
        Commands::Blastn(args) => args.out = path,
        Commands::Tblastx(args) => args.out = path,
        Commands::Blastx(args) => args.out = path,
        _ => {}
    }
}

/// The warnings that the import posts before the program's own (`RecoverSearchStrategy`'s and
/// the strategy's own, before the arguments' options are checked), taken from the arguments so
/// that the search, which writes them itself in one of its entries (tblastn's), does not write
/// them twice. A search of records writes them once at its start, first. The warnings for the
/// options that the strategy replaces (`x_IssueWarningsForIgnoredOptions`) are the search's own,
/// after the formatting handler's.
pub fn take_warnings(command: &mut Commands) -> Vec<u8> {
    let Some(run) = import_run_mut(command) else {
        return Vec::new();
    };
    let lines = run.warning_lines();
    run.warnings.clear();
    lines
}

/// A search strategy read for a search (`Search::read_strategy`): the options that the
/// command line's import takes from the search's argv and from the strategy, and the queries and
/// subjects (or the database) that the strategy carries. [`Inputs::strategy`](super::Inputs::strategy)
/// searches them; the search is the one that read the strategy.
#[derive(Debug)]
pub struct Strategy {
    program: Program,
    argv: Vec<String>,
    /// The arguments of the engine's import of the strategy file, made once when the strategy is
    /// read: each of the search's calls takes a copy (`command`). The sequence texts are not here
    /// (they became the records below), and each copy has the import's warnings, which a call
    /// takes for itself (`take_warnings`). Nothing changes them after the read, so a call that
    /// panics leaves them whole: the strategy stays unwind-safe, as [`Inputs`](super::Inputs) and
    /// [`Batches`](super::Batches) that refer to it are, although the arguments hold trait objects
    /// (the window masker's counts source).
    applied: std::panic::AssertUnwindSafe<Commands>,
    queries: Vec<FastaRecord>,
    /// How the reading of the strategy's queries ended (`InputRecords::query_end`), which the
    /// search meets where the command line's query batches meet it.
    query_end: super::QueryEnd,
    subjects: Vec<FastaRecord>,
    database: Option<String>,
}

impl Strategy {
    /// The records of the strategy's queries.
    pub(crate) fn queries(&self) -> &[FastaRecord] {
        &self.queries
    }

    /// How the reading of the strategy's queries ended.
    pub(crate) fn query_end(&self) -> &super::QueryEnd {
        &self.query_end
    }

    /// The records of the strategy's subjects, none for a database.
    pub(crate) fn subjects(&self) -> &[FastaRecord] {
        &self.subjects
    }

    /// The strategy's database name.
    pub(crate) fn database(&self) -> Option<&str> {
        self.database.as_deref()
    }

    /// The arguments of `search`, which must be the search that read this strategy.
    pub(crate) fn command(&self, search: &Search) -> Result<Commands, Error> {
        if search.program() != self.program || search.argv() != self.argv {
            return Err(Error::Argv(
                "the search strategy was read for another search: read it with the search that runs it"
                    .into(),
            ));
        }
        // The import was made once, at the reading of the strategy (`Search::read_strategy`): the
        // arguments of a call are a copy of it, with the import's warnings for the call, and no
        // sequence text, as the call searches the records read once. The search's `-out` is the
        // only thing of the call that the copy takes.
        let mut command = self.applied.0.clone();
        if let Some(path) = search.out_path() {
            set_out(&mut command, path);
        }
        Ok(command)
    }
}

impl Search {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:128-133
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3836
    /// Reads a search strategy file for this search: the command line's `-import_search_strategy
    /// FILE`, with the file's bytes. The strategy's options become the search's, set on the
    /// program's defaults in the file's order; of this search's argv only the options that NCBI's
    /// import takes from the command line count (the formatting options, `-num_threads`,
    /// `-query_loc`, `-parse_deflines`, `-lcase_masking`, `-strand` ...), and the others are
    /// ignored with NCBI's warnings, which the search writes to its diagnostics first. The
    /// strategy's queries and subjects (or its database) are the inputs of
    /// [`Inputs::strategy`](super::Inputs::strategy); a strategy's subjects are searched once
    /// with the argv's `-num_threads`, and the output is that of `-num_threads 1`, as on the
    /// command line.
    ///
    /// A strategy of another program runs as that program where the command line runs it so
    /// (blastn and blastx with a TBLASTX strategy, tblastx with a BLASTX or BLASTN one): the
    /// records are read as that program reads them.
    ///
    /// The error is the command line's: a file that is not a strategy, a strategy of a program the
    /// strategy of which is not run by this one, a strategy the engine does not run yet (a PSSM
    /// query, standard XML, a query by Seq-id), or `-import_search_strategy`'s refusal for
    /// psiblast.
    pub fn read_strategy(&self, strategy: &[u8]) -> Result<Strategy, Error> {
        // No call has started (`validate` reads its argv so too): the posts of an argument error
        // are the error's.
        let command = self.strategy_command(
            None,
            &crate::blastinput::ncbi_environment::DiagTarget::Stderr,
        )?;
        let applied = apply(command, strategy)?;
        let records = read_search_inputs(
            program_of(&applied.command).unwrap_or(self.program()),
            &applied.command,
            None,
            true,
            &applied.queries,
            &applied.subjects,
        )?;
        Ok(Strategy {
            program: self.program(),
            argv: self.argv().to_vec(),
            applied: std::panic::AssertUnwindSafe(applied.command),
            queries: records.query,
            query_end: records.query_end,
            subjects: records.subject,
            database: applied.database,
        })
    }
}

#[cfg(test)]
#[path = "strategy_tests.rs"]
mod tests;
