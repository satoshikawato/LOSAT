//! The harness shared by the `run_local` tests of every program that compare a search with the
//! `losat` binary: the temporary FASTA files, and a search of files through `LOSAT::api` (the
//! library's entry, which calls the same `run_local_*` functions as the command line), with the
//! reports, the warnings and the error that the command line gives. The tests of the entry's own
//! contract (its sinks, its observer, its refusals of hand-set arguments) are in the crate
//! (`LOSAT/src/api/local_blast_tests/`).
#![allow(
    dead_code,
    clippy::duplicate_mod,
    reason = "each test file that declares this module uses a different subset of it, and the blastp and tblastn suites load it through several #[path] files"
)]

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use LOSAT::api::{
    read_fasta, read_input_file, read_inputs, run, run_collect, run_to_file, validate, Error,
    FastaRecord, Inputs, Outcome, Program, Role,
};

/// The sequence of the first record of a FASTA file in `LOSAT/tests/fasta`.
pub fn fixture_sequence(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    read_records(&path)
        .into_iter()
        .next()
        .expect("fixture record")
        .seq()
        .to_vec()
}

/// The `-out` file of the CLI run of `outfmt` for the query `query`: the query path with the
/// extension `cli<outfmt>.out`, where every character of `outfmt` that is not an ASCII letter
/// or digit (a space, `|`, `=`) is written `_`, so the name is a valid file name on Windows too.
pub fn cli_out_path(query: &Path, outfmt: &str) -> PathBuf {
    let tag: String = outfmt
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    query.with_extension(format!("cli{tag}.out"))
}

/// A path in the temporary directory that no other call of this process gets. The tests of
/// one binary run at the same time and remove their files when they end; the process id and
/// the clock alone give two of them the same path where the clock ticks in microseconds
/// (macOS), so a counter of the calls is part of the name. A suite binary holds one copy of this
/// module per file that declares it, each with its own counter, so the module path is part of
/// the name too.
pub fn temp_path(name: &str) -> PathBuf {
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "losat_run_local_{}_{}_{nanos}_{}_{name}",
        module_path!().replace("::", "-"),
        std::process::id(),
        CALLS.fetch_add(1, Ordering::Relaxed)
    ))
}

/// A temporary FASTA file, removed when dropped.
pub struct TempFasta(pub PathBuf);

impl TempFasta {
    /// Writes `records` (id, sequence) with 60 residues per line.
    pub fn new(name: &str, records: &[(&str, &[u8])]) -> Self {
        let path = temp_path(name);
        let mut text = Vec::new();
        for (id, sequence) in records {
            text.extend_from_slice(format!(">{id}\n").as_bytes());
            for line in sequence.chunks(60) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
        }
        std::fs::write(&path, text).expect("write temporary FASTA");
        Self(path)
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// The records of a fixture of either molecule (the protein reader takes every nucleotide letter).
pub fn read_records(path: &Path) -> Vec<FastaRecord> {
    read_fasta(
        Program::Blastp,
        Role::Query,
        &std::fs::read(path).expect("FASTA"),
    )
    .expect("records")
}

/// The files of a search, shown in the reports by the names `query` and `subject` (the
/// command line's `-query` and `-subject` arguments, which a report shows as they are given).
#[derive(Clone, Copy)]
pub struct Files<'a> {
    pub program: Program,
    pub query: &'a Path,
    pub subject: &'a Path,
    /// The names the reports show, where they are not the paths.
    pub names: Option<(&'a str, &'a str)>,
}

impl<'a> Files<'a> {
    pub fn new(program: Program, query: &'a Path, subject: &'a Path) -> Self {
        Self {
            program,
            query,
            subject,
            names: None,
        }
    }

    /// The search of `options` (the command line's, without the inputs and the output) reading
    /// the files as the command line reads them, and `with` its records.
    fn with_records<T>(
        &self,
        options: &[&str],
        with: impl FnOnce(&LOSAT::api::Search, Inputs<'_>) -> Result<T, Error>,
    ) -> Result<T, Error> {
        let search = validate(self.program, options)?;
        let (query, subject) = read_inputs(
            &search,
            &read_input_file(Role::Query, self.query)?,
            &read_input_file(Role::Subject, self.subject)?,
        )?;
        let (query_name, subject_name) = match self.names {
            Some(names) => names,
            None => (
                &*self.query.display().to_string(),
                &*self.subject.display().to_string(),
            ),
        };
        with(
            &search,
            Inputs::new(&query, &subject).named(query_name, subject_name),
        )
    }

    /// `run_collect` of `formats`: the reports, the records and the warnings.
    pub fn try_collect(&self, formats: &[&str], options: &[&str]) -> Result<Outcome, Error> {
        self.with_records(options, |search, inputs| {
            run_collect(search, inputs, formats)
        })
    }

    /// `try_collect` of a search that must succeed.
    pub fn collect(&self, formats: &[&str], options: &[&str]) -> Outcome {
        self.try_collect(formats, options)
            .unwrap_or_else(|error| panic!("{formats:?} {options:?}: {error}"))
    }

    /// `run_to_file`: the report of `outfmt` in the file `path` (outfmt 13 and 14 write their side
    /// files next to it), and the warnings.
    pub fn write_file(
        &self,
        outfmt: &str,
        options: &[&str],
        path: &Path,
    ) -> (Result<(), Error>, Vec<u8>) {
        let mut diagnostics = Vec::new();
        let result = self.with_records(options, |search, inputs| {
            run_to_file(search, inputs, outfmt, path, &mut diagnostics)
        });
        (result, diagnostics)
    }

    /// `run`: one writer per format, and what the writers hold when the search fails.
    pub fn write(&self, formats: &[&str], options: &[&str]) -> Written {
        let mut outputs = vec![Vec::new(); formats.len()];
        let mut diagnostics = Vec::new();
        let result = self.with_records(options, |search, inputs| {
            let mut reports: Vec<(&str, &mut (dyn std::io::Write + Send))> = formats
                .iter()
                .copied()
                .zip(outputs.iter_mut())
                .map(|(outfmt, writer)| (outfmt, writer as &mut (dyn std::io::Write + Send)))
                .collect();
            run(search, inputs, &mut reports, &mut diagnostics)
        });
        Written {
            result,
            outputs,
            diagnostics,
        }
    }
}

/// What `Files::write` gave: the result of the search and what the writers hold.
pub struct Written {
    pub result: Result<(), Error>,
    pub outputs: Vec<Vec<u8>>,
    pub diagnostics: Vec<u8>,
}

/// The `index`-th tab-separated field of a tabular row.
pub fn field(row: &[u8], index: usize) -> String {
    let row = std::str::from_utf8(row).expect("utf-8 row").trim_end();
    row.split('\t').nth(index).expect("field").to_string()
}
