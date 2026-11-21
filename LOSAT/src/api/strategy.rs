//! A search strategy file as an input of a search (`-import_search_strategy` of the command line,
//! line RG item F2-4).
//!
//! The command line reads the file named by `-import_search_strategy` in `RecoverSearchStrategy`
//! and each program's `apply_strategy` turns it into the program's arguments, with the queries
//! and the subjects (or the database's name) that the strategy carries. The library API and LOSAT
//! Web ABI v2 give the strategy as bytes and come through the same `apply_strategy`: one engine
//! path (`apply`), then the program's `run_local` over records, as for any other search.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:529-552
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3800-3836
#![warn(missing_docs, unused)]

use std::path::Path;

use super::error::Error;
use super::input::read_search_inputs;
use super::program::Program;
use super::search::Search;
use crate::cli::Commands;
use crate::record::FastaRecord;

/// What a strategy file gives a search: the arguments as `apply_strategy` made them, and the
/// queries and subjects that the strategy carries as FASTA bytes.
#[derive(Debug)]
pub struct Applied {
    /// The program's arguments with the strategy's options set (the engine's `apply_strategy`).
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

/// Gives the arguments the strategy's `bytes` for their `apply_strategy` to read, in place of the
/// file that `-import_search_strategy` names (the name stays as the parser kept it, for the
/// messages that quote it), and the running program's NCBI name for the warnings.
fn set_source(command: &mut Commands, bytes: &[u8]) {
    use crate::report::archive::read::import::ImportRun;
    let source = std::sync::Arc::new(bytes.to_vec());
    let run = |app: &'static str| ImportRun {
        source: Some(source.clone()),
        app,
        ..ImportRun::default()
    };
    match command {
        Commands::Blastp(args) => args.import = Some(run("blastp")),
        Commands::Tblastn(args) => args.import = Some(run("tblastn")),
        Commands::Blastn(args) => args.import = Some(run("blastn")),
        Commands::Tblastx(args) => args.import = Some(run("tblastx")),
        Commands::Blastx(args) => args.import = Some(run("blastx")),
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {}
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:128-133
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:529-552
/// `RecoverSearchStrategy` of the arguments `command` (parsed with an `-import_search_strategy`
/// of any name) over the strategy `bytes`: the engine's `apply_strategy` of the program, which
/// the command line runs at the start of its run. The options the command line may give
/// besides (`-outfmt`, `-num_threads`, `-query_loc` ...) are those that the strategy's import
/// takes from them, as there.
///
/// Errors are the command line's messages (after the warnings it posts before them): a file that is not a strategy is NCBI's `Failed to
/// read search strategy`, a strategy of another program or with a PSSM query is refused as on the
/// command line.
pub fn apply(mut command: Commands, bytes: &[u8]) -> Result<Applied, Error> {
    set_source(&mut command, bytes);
    // The warnings that the command line posts before a failure (to its standard error) come
    // first in the error's text, as the two come on a terminal.
    let mut posted = Vec::new();
    let command = apply_one(command, &mut posted).map_err(|error| {
        posted.extend_from_slice(format!("{error:#}").as_bytes());
        Error::Input(posted)
    })?;
    let text =
        |held: Option<&std::sync::Arc<Vec<u8>>>| held.map(|text| text.to_vec()).unwrap_or_default();
    let run = import_run(&command);
    let (queries, subjects) = (
        text(run.and_then(|run| run.queries.as_ref())),
        text(run.and_then(|run| run.subjects.as_ref())),
    );
    let database = database_of(&command);
    Ok(Applied {
        command,
        queries,
        subjects,
        database,
    })
}

/// The program's `apply_strategy` (one engine path for the command line, the library and Web).
fn apply_one(command: Commands, diagnostics: &mut dyn std::io::Write) -> anyhow::Result<Commands> {
    use crate::algorithm::{blastn, blastp, blastx, tblastn, tblastx};
    Ok(match command {
        Commands::Blastp(args) => {
            Commands::Blastp(blastp::archive::apply_strategy(args, diagnostics)?)
        }
        Commands::Tblastn(args) => {
            Commands::Tblastn(tblastn::archive::apply_strategy(args, diagnostics)?)
        }
        Commands::Blastn(args) => {
            Commands::Blastn(blastn::archive::apply_strategy(args, diagnostics)?)
        }
        Commands::Tblastx(args) => {
            Commands::Tblastx(tblastx::archive::apply_strategy(args, diagnostics)?)
        }
        Commands::Blastx(args) => {
            Commands::Blastx(blastx::archive::apply_strategy(args, diagnostics)?)
        }
        other => {
            anyhow::bail!("-import_search_strategy is not available for this program: {other:?}")
        }
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

/// The query range that the strategy gives (`-query_loc` of the command line's import comes from the
/// strategy's `QueryRange`, `Start::query_range`): the search's range where the argv has none.
pub(crate) fn query_range(
    command: &Commands,
) -> Option<crate::blastinput::seq_range::SequenceRange> {
    import_run(command).and_then(|run| run.query_range)
}

/// The database that the applied arguments search (`-db`), when the strategy's subject is one.
fn database_of(command: &Commands) -> Option<String> {
    match command {
        Commands::Blastp(args) => args.db.clone(),
        Commands::Tblastn(args) => args.db.clone(),
        Commands::Blastn(args) => args.db.clone(),
        Commands::Tblastx(args) => args.db.clone(),
        Commands::Blastx(args) => args.db.clone(),
        _ => None,
    }
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
    /// The strategy file's bytes, applied again for each of the search's calls (`command`).
    bytes: Vec<u8>,
    queries: Vec<FastaRecord>,
    subjects: Vec<FastaRecord>,
    database: Option<String>,
}

impl Strategy {
    /// The records of the strategy's queries.
    pub(crate) fn queries(&self) -> &[FastaRecord] {
        &self.queries
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
        // The engine's `apply_strategy` over the file again: `Commands` has no copy, and the
        // arguments of a call are the import's of that call (its warnings, its files).
        let mut command = apply(search.strategy_command(None)?, &self.bytes)?.command;
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
    /// [`Inputs::strategy`](super::Inputs::strategy); a strategy with subjects is searched with
    /// one thread whatever `-num_threads` says, as on the command line.
    ///
    /// The error is the command line's: a file that is not a strategy, a strategy of a program the
    /// strategy of which is not run by this one, a strategy the engine does not run yet (a PSSM
    /// query, standard XML, a query by Seq-id), or `-import_search_strategy`'s refusal for
    /// psiblast.
    pub fn read_strategy(&self, strategy: &[u8]) -> Result<Strategy, Error> {
        let command = self.strategy_command(None)?;
        let applied = apply(command, strategy)?;
        let (queries, subjects) = read_search_inputs(
            self.program(),
            &applied.command,
            None,
            true,
            &applied.queries,
            &applied.subjects,
        )?;
        Ok(Strategy {
            program: self.program(),
            argv: self.argv().to_vec(),
            bytes: strategy.to_vec(),
            queries,
            subjects,
            database: applied.database,
        })
    }
}
