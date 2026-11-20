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

use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

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

/// Points the arguments' `-import_search_strategy` at `path`.
fn set_import_name(command: &mut Commands, path: PathBuf) {
    match command {
        Commands::Blastp(args) => args.import_search_strategy = Some(path),
        Commands::Tblastn(args) => args.import_search_strategy = Some(path),
        Commands::Blastn(args) => args.import_search_strategy = Some(path),
        Commands::Tblastx(args) => args.import_search_strategy = Some(path),
        Commands::Blastx(args) => args.import_search_strategy = Some(path),
        Commands::Psiblast(_)
        | Commands::Rpsblast(_)
        | Commands::Rpstblastn(_)
        | Commands::Makedb(_)
        | Commands::Rpsbproc(_) => {}
    }
}

/// A strategy file made of bytes, which the engine's reader opens by its name; it is removed
/// when this value is dropped.
///
/// No NCBI counterpart. The reader of the strategy takes a path (`report::archive::read::import`),
/// so the bytes go to a file for the one call of `apply_strategy`. Where there is no file system
/// (the WASI targets, whose browser host preopens no directory) there is no such file: the reader
/// has to take the bytes before the Web can import a strategy (the note at `apply`).
struct StrategyFile {
    path: PathBuf,
}

impl StrategyFile {
    #[cfg(not(target_os = "wasi"))]
    fn create(bytes: &[u8]) -> Result<Self, Error> {
        use std::io::Write;
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        #[allow(clippy::disallowed_methods)] // not WASI: `temp_dir()` panics on the WASI targets
        let directory = std::env::temp_dir();
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|since| since.as_nanos())
            .unwrap_or(0);
        let path = directory.join(format!(
            ".losat-strategy-{}-{stamp}-{}",
            crate::utils::process_id::process_id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let made = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .and_then(|mut file| file.write_all(bytes));
        match made {
            Ok(()) => Ok(Self { path }),
            Err(error) => {
                let _ = std::fs::remove_file(&path);
                Err(Error::Search(format!(
                    "the search strategy could not be written to a temporary file: {error}"
                )))
            }
        }
    }

    #[cfg(target_os = "wasi")]
    fn create(_bytes: &[u8]) -> Result<Self, Error> {
        Err(Error::Search(
            "-import_search_strategy is not available on this target: it has no file system for the strategy"
                .to_string(),
        ))
    }
}

impl Drop for StrategyFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

/// The bytes of an imported file (NCBI's temporary FASTA file of the strategy's queries or
/// subjects), from its start.
fn bytes_of(file: &std::fs::File) -> Result<Vec<u8>, Error> {
    let mut reader = file;
    let mut bytes = Vec::new();
    reader
        .seek(SeekFrom::Start(0))
        .and_then(|_| reader.read_to_end(&mut bytes))
        .map_err(|error| {
            Error::Search(format!("the strategy's sequences cannot be read: {error}"))
        })?;
    Ok(bytes)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:128-133
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:529-552
/// `RecoverSearchStrategy` of the arguments `command` (parsed with an `-import_search_strategy`
/// of any name) over the strategy `bytes`: the engine's `apply_strategy` of the program, which
/// the command line runs at the start of its run. The options the command line may give
/// besides (`-outfmt`, `-num_threads`, `-query_loc` ...) are those that the strategy's import
/// takes from them, as there.
///
/// Errors are the command line's messages: a file that is not a strategy is NCBI's `Failed to
/// read search strategy`, a strategy of another program or with a PSSM query is refused as on the
/// command line.
pub fn apply(mut command: Commands, bytes: &[u8]) -> Result<Applied, Error> {
    let file = StrategyFile::create(bytes)?;
    set_import_name(&mut command, file.path.clone());
    let mut command = apply_one(command).map_err(|error| Error::Input(format!("{error:#}")))?;
    let (queries, subjects) = {
        let run = import_run(&command);
        let read = |file: Option<&std::sync::Arc<std::fs::File>>| match file {
            Some(file) => bytes_of(file),
            None => Ok(Vec::new()),
        };
        (
            read(run.and_then(|run| run.queries.as_ref()))?,
            read(run.and_then(|run| run.subjects.as_ref()))?,
        )
    };
    let database = database_of(&command);
    // BLASTX's search of records reads the name of the subject set from `-subject`, where the
    // command line's run reads the strategy's subjects with no name (the empty name: `User
    // specified sequence set.`, `report::pairwise::scan_mode_definition`); the other programs'
    // searches of records take their labels from the call.
    if let Commands::Blastx(args) = &mut command {
        if args.subject.is_none() && args.db.is_none() && !subjects.is_empty() {
            args.subject = Some(PathBuf::new());
        }
    }
    drop(file);
    Ok(Applied {
        command,
        queries,
        subjects,
        database,
    })
}

/// The program's `apply_strategy` (one engine path for the command line, the library and Web).
fn apply_one(command: Commands) -> anyhow::Result<Commands> {
    use crate::algorithm::{blastn, blastp, blastx, tblastn, tblastx};
    Ok(match command {
        Commands::Blastp(args) => Commands::Blastp(blastp::archive::apply_strategy(args)?),
        Commands::Tblastn(args) => Commands::Tblastn(tblastn::archive::apply_strategy(args)?),
        Commands::Blastn(args) => Commands::Blastn(blastn::archive::apply_strategy(args)?),
        Commands::Tblastx(args) => Commands::Tblastx(tblastx::archive::apply_strategy(args)?),
        Commands::Blastx(args) => Commands::Blastx(blastx::archive::apply_strategy(args)?),
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
/// `x_IssueWarningsForIgnoredOptions`'), taken from the arguments so that the search, which
/// writes the first of them itself in one of its entries (tblastn's), does not write them twice.
/// A search of records writes them once at its start, first.
pub fn take_warnings(command: &mut Commands) -> Vec<u8> {
    let Some(run) = import_run_mut(command) else {
        return Vec::new();
    };
    let mut lines = run.warning_lines();
    lines.extend(run.option_warning_lines());
    run.warnings.clear();
    run.option_warnings.clear();
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
                    .to_string(),
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
