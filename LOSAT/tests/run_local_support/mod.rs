//! The harness shared by the `run_local` tests of every program that compare a search with the
//! `losat` binary: the temporary FASTA files, and a search of files through `losat::api` (the
//! library's entry, which calls the same `run_local_*` functions as the command line), with the
//! reports, the warnings and the error that the command line gives. The tests of the entry's own
//! contract (its sinks, its observer, its refusals of hand-set arguments) are in the crate
//! (`LOSAT/src/api/local_blast_tests/`).
#![allow(
    dead_code,
    clippy::duplicate_mod,
    reason = "each test file that declares this module uses a different subset of it, and the blastp and tblastn suites load it through several #[path] files"
)]

use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use losat::api::{
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

/// The outcome of `losat <program> -query <query> -subject <subject> <args>` with the variables
/// `env` set and the batch size and `OLD_FSC` removed otherwise: the exit status, stdout and
/// stderr.
pub fn subject_cli(
    program: &str,
    query: &Path,
    subject: &Path,
    args: &[&str],
    env: &[(&str, &str)],
) -> (Option<i32>, String, String) {
    let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_losat"));
    command
        .arg(program)
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .args(args)
        .env_remove("LOSAT_TIMING")
        .env_remove("BATCH_SIZE")
        .env_remove("OLD_FSC")
        // `cargo test --all-features` builds the probe of blastx's pool workers, which needs a log.
        .env(
            "LOSAT_BLASTX_WORKER_LOG",
            std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
                .unwrap_or_else(|| temp_path("subject_cli_worker.log").into_os_string()),
        )
        .envs(env.iter().copied());
    let output = command.output().expect("run LOSAT CLI");
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
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
        with: impl FnOnce(&losat::api::Search, Inputs<'_>) -> Result<T, Error>,
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

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     // Add tab in front of field, except for the first field.
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
/// The ranges of `Report::spans` of report `index`, one per HSP of `Outcome::hsps` in the order of
/// the hit list; a range is missing only where the report does not show the HSP.
fn spans_of(outcome: &Outcome, index: usize, context: &str) -> Vec<Range<usize>> {
    let report = &outcome.reports[index];
    assert_eq!(report.spans.len(), outcome.hsps.len(), "{context}");
    report
        .spans
        .iter()
        .enumerate()
        .map(|(hsp, span)| {
            span.clone()
                .unwrap_or_else(|| panic!("{context}: HSP {hsp} has no range in report {index}"))
        })
        .collect()
}

/// Checks the ranges of the tabular reports of `outcome`: the outfmt 6 rows (at position `fmt6`
/// of the request) are in HSP index order and tile the output, one line each; with `fmt7`, the
/// outfmt 7 rows are its non-comment lines and equal the outfmt 6 rows. Returns the number of
/// rows.
pub fn assert_tabular_ranges(
    outcome: &Outcome,
    fmt6: usize,
    fmt7: Option<usize>,
    context: &str,
) -> usize {
    let (out6, rows6) = (
        &outcome.reports[fmt6].bytes,
        spans_of(outcome, fmt6, context),
    );
    let mut position = 0;
    for (hsp, row) in rows6.iter().enumerate() {
        assert_eq!(row.start, position, "{context}: outfmt 6 HSP {hsp} start");
        assert_eq!(
            out6[row.clone()]
                .iter()
                .filter(|&&byte| byte == b'\n')
                .count(),
            1,
            "{context}: outfmt 6 HSP {hsp} is one line"
        );
        position = row.end;
    }
    assert_eq!(
        position,
        out6.len(),
        "{context}: outfmt 6 rows tile the output"
    );

    if let Some(fmt7) = fmt7 {
        let (out7, rows7) = (
            &outcome.reports[fmt7].bytes,
            spans_of(outcome, fmt7, context),
        );
        for (hsp, (row6, row7)) in rows6.iter().zip(&rows7).enumerate() {
            assert_eq!(
                &out6[row6.clone()],
                &out7[row7.clone()],
                "{context}: outfmt 7 HSP {hsp}"
            );
        }
        let data_lines = out7
            .split_inclusive(|&byte| byte == b'\n')
            .filter(|line| !line.starts_with(b"#"))
            .count();
        assert_eq!(data_lines, rows6.len(), "{context}: outfmt 7 data lines");
    }
    rows6.len()
}

/// Checks the ranges of outfmt 0, 6 and 7 of `outcome` (at the given positions of the request)
/// against its hit records: `assert_tabular_ranges`; the same index has the same coordinates in
/// the hit record and in its outfmt 6 row; each outfmt 0 section is one score block with the raw
/// score and the query and subject coordinates of its HSP, inside the report block of its query;
/// and the subject headings: one per subject of each query, marked with its first HSP and ending
/// where its section starts.
pub fn assert_observer_ranges(outcome: &Outcome, [fmt0, fmt6, fmt7]: [usize; 3], context: &str) {
    let hits = &outcome.hsps;
    assert!(!hits.is_empty(), "{context}: the fixture must produce hits");
    let rows = assert_tabular_ranges(outcome, fmt6, Some(fmt7), context);
    assert_eq!(rows, hits.len(), "{context}: outfmt 6 rows");
    let rows6 = spans_of(outcome, fmt6, context);
    for (hit, row) in hits.iter().zip(&rows6) {
        assert_eq!(
            [6, 7, 8, 9].map(|index| field(&outcome.reports[fmt6].bytes[row.clone()], index)),
            [hit.q_start, hit.q_end, hit.s_start, hit.s_end].map(|value| value.to_string()),
            "{context}: outfmt 6 HSP {} coordinates",
            hit.index
        );
    }

    let out0 = &outcome.reports[fmt0].bytes;
    let sections = spans_of(outcome, fmt0, context);
    // Byte offsets of the "Query=" lines that open each query's report block.
    let query_blocks: Vec<usize> = out0
        .windows(7)
        .enumerate()
        .filter(|&(offset, window)| {
            window == b"Query= " && (offset == 0 || out0[offset - 1] == b'\n')
        })
        .map(|(offset, _)| offset)
        .collect();
    let mut blocks_of_queries = std::collections::BTreeMap::new();
    for (hit, section) in hits.iter().zip(&sections) {
        let hsp = hit.index;
        let text = std::str::from_utf8(&out0[section.clone()]).expect("utf-8 section");
        assert!(
            text.starts_with(" Score ="),
            "{context}: section {hsp}: {text:.20?}"
        );
        assert_eq!(
            text.matches(" Score =").count(),
            1,
            "{context}: section {hsp}"
        );
        assert!(
            text.contains(&format!("bits ({}),", hit.raw_score)),
            "{context}: raw score of HSP {hsp}"
        );
        assert_eq!(
            alignment_span(text, "Query"),
            (hit.q_start, hit.q_end),
            "{context}: query coordinates of HSP {hsp}"
        );
        assert_eq!(
            alignment_span(text, "Sbjct"),
            (hit.s_start, hit.s_end),
            "{context}: subject coordinates of HSP {hsp}"
        );
        // The section lies in one report block, the block of its query: the block of a query
        // is the same for all its HSPs and later for a later query (a query that the search
        // skips has no block, so the index of the block is not always the index of the query),
        // and where every query has a block, it is the query's own.
        let block = query_blocks.partition_point(|&offset| offset < section.start);
        assert!(block > 0, "{context}: HSP {hsp} before its query block");
        let block = block - 1;
        if let Some(&next) = query_blocks.get(block + 1) {
            assert!(
                section.end <= next,
                "{context}: HSP {hsp} after its query block"
            );
        }
        if query_blocks.len() == outcome.query_stats.len() {
            assert_eq!(
                block, hit.q_idx,
                "{context}: HSP {hsp} in the block of its query"
            );
        }
        let (earlier, mut later) = (
            blocks_of_queries.range(..hit.q_idx),
            blocks_of_queries.range(hit.q_idx + 1..),
        );
        assert!(
            earlier.last().is_none_or(|(_, &b)| b < block)
                && later.next().is_none_or(|(_, &b)| b > block),
            "{context}: HSP {hsp}: the blocks of the queries are in the order of the queries"
        );
        assert_eq!(
            *blocks_of_queries.entry(hit.q_idx).or_insert(block),
            block,
            "{context}: HSP {hsp}: one block for all HSPs of a query"
        );
    }

    // One heading per subject of each query, just before the section of its first HSP; the
    // other HSPs of the subject have the same heading.
    let headings = &outcome.reports[fmt0].subject_headings;
    assert_eq!(headings.len(), hits.len(), "{context}: headings");
    for (index, hit) in hits.iter().enumerate() {
        let heading = headings[index]
            .clone()
            .unwrap_or_else(|| panic!("{context}: HSP {index} has no heading"));
        let first =
            index == 0 || (hit.q_idx, hit.s_idx) != (hits[index - 1].q_idx, hits[index - 1].s_idx);
        if first {
            let text = std::str::from_utf8(&out0[heading.clone()]).expect("utf-8 heading");
            assert!(
                text.starts_with("> "),
                "{context}: heading of HSP {index}: {text:.20?}"
            );
            assert_eq!(
                heading.end, sections[index].start,
                "{context}: heading of HSP {index} ends at its section"
            );
        } else {
            assert_eq!(
                headings[index],
                headings[index - 1],
                "{context}: HSP {index} has the heading of its subject"
            );
        }
    }
}

/// The first and last coordinates on the `label` ("Query" or "Sbjct") lines of one
/// outfmt 0 alignment section.
fn alignment_span(section: &str, label: &str) -> (usize, usize) {
    let numbers: Vec<Vec<usize>> = section
        .lines()
        .filter(|line| line.starts_with(label))
        .map(|line| {
            line.split_whitespace()
                .filter_map(|word| word.parse().ok())
                .collect()
        })
        .collect();
    let first = *numbers
        .first()
        .and_then(|line| line.first())
        .expect("first coordinate");
    let last = *numbers
        .last()
        .and_then(|line| line.last())
        .expect("last coordinate");
    (first, last)
}

/// The build tree that NCBI's reports name: that of the official Linux x64 release of BLAST+
/// 2.17.0 (`ncbi_exception::NCBI_BUILD_ROOT` in the crate).
pub const NCBI_BUILD_ROOT: &str = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";

/// NCBI's report of the `CStringException` that `NStr::StringToInt` (flags 0) throws for `value`
/// at `line` of ncbistr.cpp and position `pos`, as `CATCH_ALL` writes it (`Error: ` and `what()`,
/// exit status 255; blast_app_util.hpp:248-251): line 640 is the `int` range check of
/// `StringToInt`, 852 the overflow of `StringToInt8`, 861 no digit and 870 a byte after the
/// digits.
pub fn string_to_int_exception(value: &str, line: u32, pos: usize) -> String {
    let (function, kind) = match line {
        640 => ("StringToInt", "int, overflow"),
        852 => ("StringToInt8", "Int8, overflow"),
        _ => ("StringToInt8", "Int8"),
    };
    format!(
        "Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/corelib/ncbistr.cpp\", line {line}: Error: (CStringException::eConvert) ncbi::NStr::{function}() - Cannot convert string '{value}' to {kind} (m_Pos = {pos})\n\n"
    )
}
