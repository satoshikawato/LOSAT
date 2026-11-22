//! A search from Python: the inputs are read and the search runs with the GIL released.
//!
//! Python cannot stop a thread inside its process, so a search stops where `losat::api` lets a
//! caller stop: between two of NCBI's query batches. `Batches` returns to Python after each
//! batch, and `run` without reports looks for a signal between them (Ctrl-C).

use std::borrow::Cow;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::Mutex;

use arrow_array::RecordBatch;
use losat::api::{self, FastaRecord, HspRecord, Inputs, Program, QueryEnd, Role, Search};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::{PyBytes, PyList, PySlice, PyTuple};
use self_cell::self_cell;

use crate::exception;
use crate::table::{of_records, HspTable};

/// A FASTA input: its bytes, or the name of a file.
#[derive(FromPyObject)]
pub enum Source {
    Bytes(PyBackedBytes),
    Path(PathBuf),
}

impl Source {
    fn bytes(&self, role: Role) -> Result<Cow<'_, [u8]>, api::Error> {
        match self {
            Source::Bytes(bytes) => Ok(Cow::Borrowed(bytes)),
            Source::Path(path) => api::read_input_file(role, path).map(Cow::Owned),
        }
    }
}

/// An input and the name that the reports show for it.
type NamedSource = (Source, String);

/// The query of a search: a FASTA input and its name, or a PSSM file (`-in_pssm` of psiblast and
/// tblastn, `api::Inputs::pssm`), or an alignment file (`-in_msa` of psiblast, `api::Inputs::msa`),
/// which the search reads. `_api.py` gives an alignment as the dict `{"alignment": path}`, a PSSM
/// as the path alone.
#[derive(FromPyObject)]
pub enum NamedQuery {
    Alignment {
        #[pyo3(item("alignment"))]
        path: PathBuf,
    },
    Fasta(NamedSource),
    Pssm(PathBuf),
}

/// The file that holds the query of a search, in place of query records.
enum QueryFile {
    Pssm(PathBuf),
    Alignment(PathBuf),
}

/// The subject of a search: a FASTA input and its name, or no input and the name of a database.
type NamedSubject = (Option<Source>, String);

/// The search of `program` with `argv`, checked by the command line's parser.
pub fn checked(program: &str, argv: &[String]) -> Result<Search, api::Error> {
    let words: Vec<&str> = argv.iter().map(String::as_str).collect();
    api::validate(Program::parse(program)?, &words)
}

/// A checked search and its records.
pub struct Loaded {
    search: Search,
    /// The PSSM or alignment file of a query that a file holds; the query records are then none.
    file: Option<QueryFile>,
    /// The query records read before the end of the query input (`query_end`).
    query: Vec<FastaRecord>,
    /// How the reading of the query input ended after `query` (`api::read_input_records`): a
    /// search of the records meets a reader error, or blank lines where a record should be, in
    /// the query batch that reads it, after the batches before, as the command line does.
    query_end: QueryEnd,
    subject: Vec<FastaRecord>,
    query_name: String,
    subject_name: String,
    /// Whether `subject_name` is a database.
    database: bool,
}

impl Loaded {
    fn read(
        program: &str,
        argv: &[String],
        query: NamedQuery,
        (subject, subject_name): NamedSubject,
    ) -> Result<Self, api::Error> {
        let search = checked(program, argv)?;
        let database = subject.is_none();
        let (query, query_name, file) = match query {
            NamedQuery::Fasta((query, query_name)) => (Some(query), query_name, None),
            NamedQuery::Pssm(path) => (None, String::new(), Some(QueryFile::Pssm(path))),
            NamedQuery::Alignment { path } => {
                (None, String::new(), Some(QueryFile::Alignment(path)))
            }
        };
        // A dry run reads no record, as the command line's does not (`Search::dry_run`).
        if search.dry_run() {
            return Ok(Self {
                search,
                file,
                query: Vec::new(),
                query_end: QueryEnd::Input,
                subject: Vec::new(),
                query_name,
                subject_name,
                database,
            });
        }
        // A PSSM or alignment query is read by the search, as the command line's `-in_pssm` and
        // `-in_msa`.
        let query_bytes = match &query {
            Some(query) => query.bytes(Role::Query)?,
            None => Cow::Borrowed(&[][..]),
        };
        let subject_bytes = match &subject {
            Some(subject) => subject.bytes(Role::Subject)?,
            None => Cow::Borrowed(&[][..]),
        };
        // A query record that the reader refuses ends the search in its batch, not here: the
        // command line reads its queries one batch at a time (line WQA's D10). A subject's
        // error fails here, as the command line reads every subject before the first query.
        let records = api::read_input_records(&search, &query_bytes, &subject_bytes)?;
        Ok(Self {
            search,
            file,
            query: records.query,
            query_end: records.query_end,
            subject: records.subject,
            query_name,
            subject_name,
            database,
        })
    }

    fn inputs(&self) -> Inputs<'_> {
        let inputs = match (self.database, self.file.as_ref()) {
            (true, Some(QueryFile::Pssm(pssm))) => Inputs::pssm_database(pssm, &self.subject_name),
            (false, Some(QueryFile::Pssm(pssm))) => Inputs::pssm(pssm, &self.subject),
            (true, Some(QueryFile::Alignment(msa))) => {
                Inputs::msa_database(msa, &self.subject_name)
            }
            (false, Some(QueryFile::Alignment(msa))) => Inputs::msa(msa, &self.subject),
            (true, None) => {
                Inputs::database(&self.query, &self.subject_name).query_end(&self.query_end)
            }
            (false, None) => Inputs::new(&self.query, &self.subject).query_end(&self.query_end),
        };
        inputs.named(&self.query_name, &self.subject_name)
    }
}

/// The statistics of a query, as `QueryStats` of `_api.py` takes them.
type StatsRow = (i64, i32, f64, f64, f64);

fn stats_rows(stats: &[Option<api::QueryStats>]) -> Vec<Option<StatsRow>> {
    stats
        .iter()
        .map(|stats| {
            stats.map(|stats| {
                (
                    stats.eff_space,
                    stats.hsp_len,
                    stats.lambda,
                    stats.kappa,
                    stats.entropy,
                )
            })
        })
        .collect()
}

/// The rounds of a query, as `QueryRounds` of `_api.py` takes them: each round's statistics and
/// whether the rounds converged.
type RoundsRow = (Vec<Option<StatsRow>>, bool);

fn rounds_rows(rounds: &[api::QueryRounds]) -> Vec<RoundsRow> {
    rounds
        .iter()
        .map(|rounds| (stats_rows(&rounds.stats), rounds.converged))
        .collect()
}

/// A list of byte ranges of a report, as Python slices of its bytes.
fn slices<'py>(
    py: Python<'py>,
    ranges: &[Option<std::ops::Range<usize>>],
) -> PyResult<Bound<'py, PyList>> {
    let slices = ranges.iter().map(|range| {
        range
            .as_ref()
            .map(|range| PySlice::new(py, range.start as isize, range.end as isize, 1))
    });
    PyList::new(py, slices)
}

/// `(reports, hsps, query stats, rounds, diagnostics)` of a search, which `_api.py` turns into an
/// `Outcome`. The reports need one call of the engine (`run_collect`). Without reports the
/// search runs one query batch at a time, which gives the same records (`run_batches`), and a
/// signal stops it between two batches.
#[pyfunction]
pub fn run<'py>(
    py: Python<'py>,
    program: &str,
    argv: Vec<String>,
    query: NamedQuery,
    subject: NamedSubject,
    formats: Vec<String>,
) -> PyResult<Bound<'py, PyTuple>> {
    let loaded = py
        .detach(|| Loaded::read(program, &argv, query, subject))
        .map_err(exception)?;
    if formats.is_empty() {
        return run_typed(py, &loaded);
    }
    let formats: Vec<&str> = formats.iter().map(String::as_str).collect();
    let (outcome, table) = py
        .detach(|| {
            let outcome = api::run_collect(&loaded.search, loaded.inputs(), &formats)?;
            let table = of_records(&outcome.hsps);
            Ok((outcome, table))
        })
        .map_err(exception)?;
    let reports = outcome
        .reports
        .iter()
        .map(|report| {
            Ok((
                report.outfmt.as_str(),
                PyBytes::new(py, &report.bytes),
                slices(py, &report.spans)?,
                slices(py, &report.subject_headings)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    (
        reports,
        HspTable::new(table),
        stats_rows(&outcome.query_stats),
        rounds_rows(&outcome.rounds),
        PyBytes::new(py, &outcome.diagnostics),
    )
        .into_pyobject(py)
}

/// The report of `outfmt` written to the file `path`, as the command line's `-out` writes it, with
/// the side files of outfmt 13 and 14 next to it; returns the bytes that the command line writes
/// to its standard error.
#[pyfunction]
pub fn run_to_file<'py>(
    py: Python<'py>,
    program: &str,
    argv: Vec<String>,
    query: NamedQuery,
    subject: NamedSubject,
    outfmt: &str,
    path: PathBuf,
) -> PyResult<Bound<'py, PyBytes>> {
    let diagnostics = py
        .detach(|| {
            let loaded = Loaded::read(program, &argv, query, subject)?;
            let mut diagnostics = Vec::new();
            api::run_to_file(
                &loaded.search,
                loaded.inputs(),
                outfmt,
                &path,
                &mut diagnostics,
            )?;
            Ok(diagnostics)
        })
        .map_err(exception)?;
    Ok(PyBytes::new(py, &diagnostics))
}

/// A tool of the library API (`api::run_makeblastdb` and the others).
type Tool = fn(&[&str], &mut dyn Read, &mut dyn Write, &mut dyn Write) -> Result<(), api::Error>;

/// The bytes that a tool writes to standard output and to standard error.
type ToolBytes<'py> = (Bound<'py, PyBytes>, Bound<'py, PyBytes>);

/// `tool` with the command line's words and `input` as its standard input, without the GIL:
/// what it writes to standard output and to standard error.
fn run_tool<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
    tool: Tool,
) -> PyResult<ToolBytes<'py>> {
    let (output, diagnostics) = py
        .detach(|| {
            let words: Vec<&str> = argv.iter().map(String::as_str).collect();
            let (mut output, mut diagnostics) = (Vec::new(), Vec::new());
            tool(&words, &mut &*input, &mut output, &mut diagnostics)?;
            Ok((output, diagnostics))
        })
        .map_err(exception)?;
    Ok((PyBytes::new(py, &output), PyBytes::new(py, &diagnostics)))
}

/// `makeblastdb` with the command line's words: `input` is its standard input. Returns what it
/// writes to standard output and to standard error.
#[pyfunction]
pub fn makeblastdb<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_makeblastdb)
}

/// `blastdbcmd` with the command line's words: `input` is its standard input (a batch file `-`).
/// Returns what it writes to standard output and to standard error.
#[pyfunction]
pub fn blastdbcmd<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_blastdbcmd)
}

/// `dustmasker` with the command line's words: `input` is its standard input (`-in -`, and
/// without `-in`). Returns what it writes to standard output and to standard error.
#[pyfunction]
pub fn dustmasker<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_dustmasker)
}

/// `segmasker` with the command line's words: `input` is its standard input (`-in -`, and
/// without `-in`). Returns what it writes to standard output and to standard error.
#[pyfunction]
pub fn segmasker<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_segmasker)
}

/// `windowmasker` with the command line's words: `input` is its standard input (`-in -`, and
/// without `-in`). Returns what it writes to standard output and to standard error.
#[pyfunction]
pub fn windowmasker<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_windowmasker)
}

/// `convert2blastmask` with the command line's words: `input` is its standard input (`-in -`,
/// and without `-in`). Returns what it writes to standard output and to standard error.
#[pyfunction]
pub fn convert2blastmask<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_convert2blastmask)
}

/// `makembindex` with the command line's words: `input` is its standard input (FASTA input of
/// `-old_style_index true` without `-input`). Returns what it writes to standard output and to
/// standard error; the index files are written where the command line writes them.
#[pyfunction]
pub fn makembindex<'py>(
    py: Python<'py>,
    argv: Vec<String>,
    input: PyBackedBytes,
) -> PyResult<ToolBytes<'py>> {
    run_tool(py, argv, input, api::run_makembindex)
}

/// `run` without reports.
fn run_typed<'py>(py: Python<'py>, loaded: &Loaded) -> PyResult<Bound<'py, PyTuple>> {
    let mut batches = api::run_batches(&loaded.search, loaded.inputs()).map_err(exception)?;
    let mut hsps: Vec<HspRecord> = Vec::new();
    let mut stats = Vec::new();
    let mut rounds = Vec::new();
    let mut diagnostics = Vec::new();
    loop {
        py.check_signals()?;
        let Some(batch) = py.detach(|| batches.next()) else {
            break;
        };
        let mut batch = batch.map_err(exception)?;
        hsps.append(&mut batch.hsps);
        stats.append(&mut batch.query_stats);
        rounds.append(&mut batch.rounds);
        diagnostics.append(&mut batch.diagnostics);
    }
    let table = py.detach(|| of_records(&hsps));
    let reports: Vec<Bound<'py, PyAny>> = Vec::new();
    (
        reports,
        HspTable::new(table),
        stats_rows(&stats),
        rounds_rows(&rounds),
        PyBytes::new(py, &diagnostics),
    )
        .into_pyobject(py)
}

type BatchIterator<'a> = api::Batches<'a>;

self_cell!(
    /// A `Batches` of the API with the search and the records that it borrows, so that it can
    /// live in a Python object.
    struct Running {
        owner: Loaded,

        #[covariant]
        dependent: BatchIterator,
    }
);

/// `(first query, end of the queries, hsps, query stats, rounds, diagnostics)` of a batch.
type BatchRow = (
    usize,
    usize,
    RecordBatch,
    Vec<Option<StatsRow>>,
    Vec<RoundsRow>,
    Vec<u8>,
);

/// The search of one query batch at a time. Each `next` searches a batch with the GIL released;
/// dropping the object stops the search, and a signal stops it before the next batch.
#[pyclass(frozen, module = "pylosat._pylosat")]
pub struct Batches {
    running: Mutex<Running>,
    query_count: usize,
}

#[pymethods]
impl Batches {
    #[new]
    fn new(
        py: Python<'_>,
        program: &str,
        argv: Vec<String>,
        query: NamedQuery,
        subject: NamedSubject,
    ) -> PyResult<Self> {
        let running = py
            .detach(|| {
                let loaded = Loaded::read(program, &argv, query, subject)?;
                Running::try_new(loaded, |loaded| {
                    api::run_batches(&loaded.search, loaded.inputs())
                })
            })
            .map_err(exception)?;
        let query_count = running.borrow_dependent().query_count();
        Ok(Self {
            running: Mutex::new(running),
            query_count,
        })
    }

    /// How many queries the inputs hold: the records, or the one sequence of a PSSM file.
    #[getter]
    fn query_count(&self) -> usize {
        self.query_count
    }

    fn __iter__(this: PyRef<'_, Self>) -> PyRef<'_, Self> {
        this
    }

    fn __next__<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyTuple>>> {
        py.check_signals()?;
        let batch = py
            .detach(|| {
                let mut running = self.running.lock().expect("a batch panicked");
                running
                    .with_dependent_mut(|_, batches| batches.next())
                    .map(|batch| batch.map(batch_row))
            })
            .transpose()
            .map_err(exception)?;
        batch
            .map(|(start, end, table, stats, rounds, diagnostics)| {
                let diagnostics = PyBytes::new(py, &diagnostics);
                (start, end, HspTable::new(table), stats, rounds, diagnostics).into_pyobject(py)
            })
            .transpose()
    }
}

fn batch_row(batch: api::Batch) -> BatchRow {
    (
        batch.query_indices.start,
        batch.query_indices.end,
        of_records(&batch.hsps),
        stats_rows(&batch.query_stats),
        rounds_rows(&batch.rounds),
        batch.diagnostics,
    )
}
