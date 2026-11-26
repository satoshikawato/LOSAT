//! The import of a run's search strategy, made once per run (`-import_search_strategy`, line RG
//! item F2-4, review 2's R2-N3).
//!
//! The engine's import (`losat::web_adapter::import_strategy`) renders every query and subject
//! of the strategy as FASTA text. LOSAT Web reads the arguments of an argv again at each of its
//! entries (`start`, `validate`, `register_inputs`, `run`) and in several places of each, so that
//! an import for each reading was a full render of the strategy's sequences each time. The import
//! of a run is held here instead, where the run's state is: the arguments that it made, and the
//! texts, which only `register_inputs` reads and which it takes (`Supplied::take`) so that the
//! adapter holds them no longer than the reading of the records.
//!
//! The import is a function of the argv's arguments and the strategy file's bytes: an entry is
//! the import of one argv and one registered file, and a reading of the same argv and file takes
//! a copy of its arguments (they carry the import's warnings, which each use takes for itself,
//! `take_strategy_warnings`). The held import is the latest one (a run that starts replaces it),
//! and the end of the run (`end_run`) drops it. It is the thread's: the host calls the ABI from
//! one thread.

use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Weak};

use losat::web_adapter::Commands;

use crate::run::Program;

/// The queries and the subjects (FASTA bytes) that the strategy carries and the argv does not
/// override with its own `-query` or `-subject`: empty for an input that the strategy does not
/// give.
pub struct Texts {
    pub queries: Vec<u8>,
    pub subjects: Vec<u8>,
}

/// The import of one argv over one registered strategy file.
struct Imported {
    words: Vec<String>,
    /// The registered file that was imported: an entry is for that file alone (a released file
    /// has no strong holder, and a file registered again is another one).
    file: Weak<Vec<u8>>,
    /// The arguments as the parser made them, before the import: what the texts are made from
    /// again, for a second reading of the inputs (`Supplied::take`).
    unapplied: Commands,
    program: Program,
    /// The arguments of the import, without the sequence texts.
    command: Commands,
    /// The texts of the import, until `register_inputs` takes them.
    texts: RefCell<Option<Texts>>,
}

thread_local! {
    static HELD: RefCell<Option<Rc<Imported>>> = const { RefCell::new(None) };
    /// The imports made on this thread.
    #[cfg(test)]
    pub static IMPORTS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// The engine's import of `command` over the strategy file `file`: the arguments of the
/// strategy's program, and the texts of its queries and subjects.
fn import_once(command: Commands, file: &[u8]) -> Result<(Commands, Texts), String> {
    #[cfg(test)]
    IMPORTS.with(|count| count.set(count.get() + 1));
    let applied =
        losat::web_adapter::import_strategy(command, file).map_err(|error| error.to_string())?;
    Ok((
        applied.command,
        Texts {
            queries: applied.queries,
            subjects: applied.subjects,
        },
    ))
}

/// The import of the argv `words`, whose arguments `command` (parsed, with the strategy named)
/// name the registered strategy file `file`, over the program `program` of the argv: the held
/// one where it is for the same argv and file, else made now (and the held one dropped before, so
/// that the texts of two imports are never held together).
fn imported(
    words: &[&str],
    program: Program,
    command: Commands,
    file: &Arc<Vec<u8>>,
) -> Result<Rc<Imported>, String> {
    let held = HELD.with(|held| held.borrow().clone());
    if let Some(entry) = held {
        let same_file = entry
            .file
            .upgrade()
            .is_some_and(|held| Arc::ptr_eq(&held, file));
        if same_file
            && entry
                .words
                .iter()
                .map(String::as_str)
                .eq(words.iter().copied())
        {
            return Ok(entry);
        }
    }
    end_run();
    let unapplied = command.clone();
    let (command, texts) = import_once(command, file)?;
    // The strategy's program, which reads the inputs and runs (`losat::web_adapter::program_of`).
    let program = match losat::web_adapter::program_of(&command) {
        Some(applied) => Program::parse(applied.name())?,
        None => program,
    };
    let entry = Rc::new(Imported {
        words: words.iter().map(|word| word.to_string()).collect(),
        file: Arc::downgrade(file),
        unapplied,
        program,
        command,
        texts: RefCell::new(Some(texts)),
    });
    HELD.with(|held| *held.borrow_mut() = Some(entry.clone()));
    Ok(entry)
}

/// The run of the held import ended (or its argv failed): the import and the texts it still holds
/// are dropped.
pub fn end_run() {
    let dropped = HELD.with(|held| held.borrow_mut().take());
    drop(dropped);
}

/// What the argv's import supplies to its run: the arguments and the program of the import, and
/// the strategy's inputs.
pub struct Supplied {
    imported: Rc<Imported>,
}

impl Supplied {
    /// The import of `words` over `file` (`imported`), once for the run.
    pub fn of(
        words: &[&str],
        program: Program,
        command: Commands,
        file: &Arc<Vec<u8>>,
    ) -> Result<Self, String> {
        imported(words, program, command, file).map(|imported| Self { imported })
    }

    /// The program that the strategy's arguments are of.
    pub fn program(&self) -> Program {
        self.imported.program
    }

    /// A copy of the arguments of the import, for one use: each use takes the import's warnings
    /// for itself.
    pub fn command(&self) -> Commands {
        self.imported.command.clone()
    }

    /// The texts of the strategy's queries and subjects, for the reading of the run's inputs
    /// (`register_inputs`): taken from the import, which then holds none. A second reading of
    /// the same import makes them again.
    pub fn take(&self) -> Result<Texts, String> {
        if let Some(texts) = self.imported.texts.borrow_mut().take() {
            return Ok(texts);
        }
        let file = self.imported.file.upgrade().ok_or_else(|| {
            "the search strategy file was released before the run read its sequences".to_string()
        })?;
        import_once(self.imported.unapplied.clone(), &file).map(|(_, texts)| texts)
    }
}
