//! The shared BLASTX entry `run_local` (which `LOSAT::api` calls) writes outfmt 0, 6 and 7 from
//! one search. Every format must be byte-identical to the CLI run that requests only that format
//! and to the frozen NCBI BLAST+ 2.17.0 bytes (`outfmt0_manifest.tsv`), the warnings must be
//! written once, as the CLI writes them to standard error, and the api's ranges must be those of
//! every HSP in every format. The records are those of NCBI's reader (`fasta_reader`), which
//! LOSAT Web registers; the CLI reads with BLASTX's own reader, so the two readers are compared
//! over the FASTA-input fixtures too. The entry's own contract (its observer, its refusals, its
//! typed errors) is checked in the crate, `LOSAT/src/api/local_blast_tests/blastx.rs`.

mod run_local_support;

use std::path::{Path, PathBuf};
use std::process::Command;

use run_local_support::{
    assert_observer_ranges, cli_out_path, fixture_sequence, temp_path, Files, Written,
};
use LOSAT::api::{Error, Outcome, Program};

const FORMATS: [&str; 3] = ["0", "6", "7"];

/// The LOSAT crate's directory, where the manifests' paths start.
fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The input paths of a run, as given to the CLI (relative to the crate's directory or
/// absolute); the api shows them as its display names.
struct Inputs {
    query: String,
    subject: String,
}

impl Inputs {
    fn new(query: impl AsRef<Path>, subject: impl AsRef<Path>) -> Self {
        Self {
            query: query.as_ref().display().to_string(),
            subject: subject.as_ref().display().to_string(),
        }
    }

    fn path(name: &str) -> PathBuf {
        crate_dir().join(name)
    }

    /// `with` the files of the search, shown by the names that the CLI is given.
    fn with_files<T>(&self, with: impl FnOnce(Files<'_>) -> T) -> T {
        let (query, subject) = (Self::path(&self.query), Self::path(&self.subject));
        let mut files = Files::new(Program::Blastx, &query, &subject);
        files.names = Some((&self.query, &self.subject));
        with(files)
    }

    /// The api's search of `formats` with the options `extra`: the reports, the records and the
    /// warnings.
    fn collect(&self, formats: &[&str], extra: &[&str]) -> Outcome {
        self.with_files(|files| files.collect(formats, extra))
    }

    /// The api's search of `formats` with the options `extra`, with what the writers hold when it
    /// fails.
    fn write(&self, formats: &[&str], extra: &[&str]) -> Written {
        self.with_files(|files| files.write(formats, extra))
    }

    /// The CLI run of one format: its exit status, output file and standard error.
    fn cli(&self, outfmt: &str, extra: &[&str]) -> (bool, Vec<u8>, Vec<u8>) {
        let out = cli_out_path(&temp_path("blastx_cli.fna"), outfmt);
        let output = Command::new(env!("CARGO_BIN_EXE_losat"))
            .current_dir(crate_dir())
            .env_remove("BATCH_SIZE")
            // `cargo test --all-features` builds the probe of the pools' workers, which needs a log.
            .env(
                "LOSAT_BLASTX_WORKER_LOG",
                std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
                    .unwrap_or_else(|| temp_path("blastx_worker.log").into_os_string()),
            )
            .arg("blastx")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(&self.subject)
            .args(extra)
            .arg("-outfmt")
            .arg(outfmt)
            .arg("-out")
            .arg(&out)
            .output()
            .expect("run LOSAT CLI");
        let bytes = std::fs::read(&out).unwrap_or_default();
        let _ = std::fs::remove_file(&out);
        (output.status.success(), bytes, output.stderr)
    }
}

/// A query of several reader batches (more than 10,002 letters in all): a record without letters
/// and a long record, then records with the reader's warnings, a record with a lower-case island
/// and a query too short to search; the subjects hold a record without letters.
struct Batches {
    query: PathBuf,
    subject: PathBuf,
}

impl Batches {
    fn new() -> Self {
        let genome = fixture_sequence("MjeNMV.fasta");
        // The first batch: a record without letters and one of 10,100 letters.
        let mut query = b">empty record\n".to_vec();
        query.extend_from_slice(b">long record of 10100 letters\n");
        for line in genome[0..10_100].chunks(70) {
            query.extend_from_slice(line);
            query.push(b'\n');
        }
        // The reader's warnings (invalid residues, hyphens) in the second batch.
        query.extend_from_slice(
            &std::fs::read(Inputs::path(
                "tests/fixtures/fasta_input/batch_late_warn.nuc.fa",
            ))
            .expect("fixture"),
        );
        query.extend_from_slice(b">lower case island\n");
        query.extend_from_slice(&genome[20_000..20_300]);
        query.extend(genome[20_300..20_600].iter().map(u8::to_ascii_lowercase));
        query.extend_from_slice(&genome[20_600..21_200]);
        query.extend_from_slice(b"\n>short\nAC\n");
        let query_path = temp_path("blastx_batches_query.fna");
        std::fs::write(&query_path, query).expect("write query");
        let mut subject =
            std::fs::read(Inputs::path("tests/fasta/MjeNMV.faa")).expect("subject fixture");
        subject.extend_from_slice(b">no letters\n");
        subject.extend_from_slice(
            &std::fs::read(Inputs::path("tests/fixtures/fasta_input/base.s.prot.faa"))
                .expect("fixture"),
        );
        let subject_path = temp_path("blastx_batches_subject.faa");
        std::fs::write(&subject_path, subject).expect("write subject");
        Self {
            query: query_path,
            subject: subject_path,
        }
    }

    fn inputs(&self) -> Inputs {
        Inputs::new(&self.query, &self.subject)
    }
}

impl Drop for Batches {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.query);
        let _ = std::fs::remove_file(&self.subject);
    }
}

/// Checks that one api call of outfmt 0, 6 and 7 equals the single-format CLI runs and,
/// for inputs with hits, that its ranges and hits agree with the reports.
fn assert_matches_cli(inputs: &Inputs, extra: &[&str], hits: bool, context: &str) {
    let together = inputs.collect(&FORMATS, extra);
    if hits {
        assert_observer_ranges(&together, [0, 1, 2], context);
    } else {
        assert!(together.hsps.is_empty(), "{context}: no hits");
        assert!(
            together
                .reports
                .iter()
                .all(|report| report.spans.iter().all(Option::is_none)),
            "{context}: no ranges"
        );
    }
    for (index, outfmt) in FORMATS.iter().enumerate() {
        let (ok, output, stderr) = inputs.cli(outfmt, extra);
        assert!(ok, "{context}: CLI outfmt {outfmt} failed");
        assert!(
            together.reports[index].bytes == output,
            "{context}: outfmt {outfmt} differs from the CLI"
        );
        assert_eq!(
            String::from_utf8_lossy(&together.diagnostics),
            String::from_utf8_lossy(&stderr),
            "{context}: the warnings of three formats equal the CLI's of outfmt {outfmt}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-297
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
// ...
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ...
//             	ITERATE(CSearchResultSet, result, *results) {
//                	    formatter.PrintOneResultSet(**result, query_batch);
//             	}
// ...
//         formatter.PrintEpilog(opt);
// ```
// The half of the test of this name in the crate (`local_blast_tests/blastx.rs`) checks the
// observer's ranges against the entry's hit records; this half compares the reports and the
// warnings with the CLI runs, and the api's ranges with its hit records.
#[test]
fn every_format_of_one_search_matches_the_cli_run_of_that_format() {
    let batches = Batches::new();
    // (inputs, whether they have hits)
    let cases = [
        (
            Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa"),
            true,
        ),
        (batches.inputs(), true),
        (
            Inputs::new("tests/fasta/small_test.fasta", "tests/fasta/PajaWSV.faa"),
            false,
        ),
    ];
    for (inputs, hits) in &cases {
        for threads in ["1", "2", "4"] {
            for extra in [&[][..], &["-lcase_masking"][..]] {
                let extra = [extra, &["-num_threads", threads]].concat();
                let context = format!("{} {extra:?}", inputs.query);
                assert_matches_cli(inputs, &extra, *hits, &context);
            }
        }
    }
    // The second batch holds the reader's warnings: they come after the reports of the first
    // batch's queries (its query warnings), as the CLI writes them.
    let run = batches.inputs().collect(&["6"], &[]);
    let diagnostics = String::from_utf8_lossy(&run.diagnostics);
    let empty = diagnostics
        .find("Query_1 empty record")
        .expect("empty query");
    let hyphens = diagnostics.find("Hyphens").expect("reader warning");
    assert!(empty < hyphens, "{diagnostics}");
    assert!(diagnostics.contains("Subject sequence contains no data"));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2973
// ```c++
//     if (args.Exist(kArgNumDescriptions) && args[kArgNumDescriptions]) {
//        m_NumDescriptions = args[kArgNumDescriptions].AsInteger();
//     }
// ```
// The pairwise report's counts differ by format: formats that need different hit lists are
// searched by separate calls (LOSAT Web's search groups), which may share one pool.
// The half of the test of this name in the crate checks that the entry refuses formats that need
// different hit list sizes; this half compares the api, which splits them into searches itself,
// with the CLI runs.
#[test]
fn display_counts_of_the_pairwise_report() {
    let inputs = Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa");
    let extra = [
        "-num_descriptions",
        "3",
        "-num_alignments",
        "2",
        "-num_threads",
        "2",
    ];
    let (pairwise, tabular) = (
        inputs.collect(&["0"], &extra),
        inputs.collect(&["6", "7"], &extra),
    );
    let (ok, output, stderr) = inputs.cli("0", &extra);
    assert!(ok);
    assert!(pairwise.reports[0].bytes == output, "outfmt 0");
    assert_eq!(pairwise.diagnostics, stderr, "outfmt 0 warnings");
    for (index, outfmt) in ["6", "7"].iter().enumerate() {
        let (ok, output, stderr) = inputs.cli(outfmt, &extra);
        assert!(ok);
        assert!(tabular.reports[index].bytes == output, "outfmt {outfmt}");
        assert_eq!(tabular.diagnostics, stderr, "outfmt {outfmt} warnings");
    }
    assert!(
        String::from_utf8_lossy(&tabular.diagnostics).contains("-num_descriptions is ignored"),
        "the tabular formats' warning"
    );
    // All three formats in one call: the api searches the pairwise format and the tabular ones
    // apart, and the warnings are those of the tabular formats' search.
    let all = inputs.collect(&FORMATS, &extra);
    assert!(
        all.reports[0].bytes == pairwise.reports[0].bytes,
        "outfmt 0"
    );
    assert!(all.reports[1].bytes == tabular.reports[0].bytes, "outfmt 6");
    assert!(all.reports[2].bytes == tabular.reports[1].bytes, "outfmt 7");
    assert_eq!(all.diagnostics, tabular.diagnostics, "warnings");
}

/// The rows of `outfmt0_manifest.tsv` of BLASTX with outfmt 0, 6 or 7 (6 and 7 with custom fields
/// too, W10) and without `-parse_deflines` (LOSAT Web refuses it): (fixture id, query, subject,
/// extra arguments, outfmt).
fn manifest_rows() -> Vec<(String, String, String, Vec<String>, String)> {
    let text =
        std::fs::read_to_string(crate_dir().join("tests/outfmt0_manifest.tsv")).expect("manifest");
    let mut lines = text.lines().filter(|line| !line.starts_with('#'));
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let column = |name: &str| header.iter().position(|&h| h == name).expect(name);
    let (id, program, task, query, subject, extra, outfmt) = (
        column("fixture_id"),
        column("program"),
        column("task"),
        column("query"),
        column("subject"),
        column("extra_args"),
        column("outfmt"),
    );
    lines
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .filter(|row| {
            row[program] == "blastx"
                && matches!(row[outfmt].split(' ').next(), Some("0" | "6" | "7"))
        })
        .filter(|row| !row[extra].contains("-parse_deflines"))
        .map(|row| {
            (
                row[id].to_string(),
                row[query].to_string(),
                row[subject].to_string(),
                (!row[task].is_empty())
                    .then(|| ["-task", row[task]])
                    .into_iter()
                    .flatten()
                    .chain(row[extra].split_whitespace())
                    .map(str::to_string)
                    .collect(),
                row[outfmt].to_string(),
            )
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
// ```c++
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
// ```
#[test]
fn formats_of_one_search_equal_the_frozen_ncbi_rows() {
    let rows = manifest_rows();
    assert!(rows.len() >= 30, "{} BLASTX rows", rows.len());
    for (id, query, subject, extra, outfmt) in rows {
        let extra: Vec<&str> = extra.iter().map(String::as_str).collect();
        // The row's format first: its warnings are the diagnostics.
        let mut formats = vec![outfmt.as_str()];
        formats.extend(FORMATS.iter().filter(|&&f| f != outfmt));
        let inputs = Inputs::new(&query, &subject);
        let run = inputs.collect(&formats, &extra);
        let expected = std::fs::read(crate_dir().join(format!("tests/fixtures/outfmt0/{id}.out")))
            .expect("frozen NCBI output");
        assert!(run.reports[0].bytes == expected, "{id}: differs from NCBI");
        if let Ok(stderr) =
            std::fs::read(crate_dir().join(format!("tests/fixtures/outfmt0/{id}.err")))
        {
            assert_eq!(
                String::from_utf8_lossy(&run.diagnostics),
                String::from_utf8_lossy(&stderr),
                "{id}: warnings"
            );
        }
    }
}

/// Whether the CLI's standard error holds the api's error message: the api joins an error and its
/// causes by `: `, the CLI writes `Error: ` and the error, then `Caused by:` and each cause on a
/// line of its own, so every part of the message comes in the standard error, in order.
fn error_is_in_stderr(message: &str, stderr: &[u8]) -> bool {
    let stderr = String::from_utf8_lossy(stderr);
    let mut from = 0;
    message
        .split(": ")
        .all(|part| match stderr[from..].find(part) {
            Some(at) => {
                from += at + part.len();
                true
            }
            None => false,
        })
}

/// The FASTA-input fixtures that BLASTX can read: each nucleotide query file against the base
/// protein subjects, and the base nucleotide queries against each protein subject file.
fn reader_cases() -> Vec<(String, String)> {
    let directory = "tests/fixtures/fasta_input";
    let mut names: Vec<String> = std::fs::read_dir(crate_dir().join(directory))
        .expect("fixtures")
        .map(|entry| {
            entry
                .expect("entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    names.sort();
    let mut cases = Vec::new();
    for name in names {
        let path = format!("{directory}/{name}");
        if name.ends_with(".q.nuc.fa") {
            cases.push((path, format!("{directory}/base.s.prot.faa")));
        } else if name.ends_with(".s.prot.faa") && name != "base.s.prot.faa" {
            cases.push((format!("{directory}/base.q.nuc.fa"), path));
        }
    }
    cases
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:328-330
// ```c++
//     flags += (m_ReadProteins
//               ? CFastaReader::fAssumeProt
//               : CFastaReader::fAssumeNuc);
// ```
// LOSAT Web registers the records of NCBI's reader (`fasta_reader`), which `run_local` turns into
// BLASTX's records; the CLI reads with BLASTX's reader. Both must give the same report and
// messages, and fail where the other fails, except where the readers differ (`ReaderDifference`).
#[test]
fn the_records_of_ncbis_reader_search_as_blastxs_reader_reads_them() {
    let cases = reader_cases();
    assert!(cases.len() > 150, "{} cases", cases.len());
    let mut differences = Vec::new();
    // Cases whose report has rows: the comparison is not one of empty reports only.
    let mut with_rows = 0;
    for (query, subject) in cases {
        for extra in [&[][..], &["-lcase_masking"][..]] {
            let inputs = Inputs::new(&query, &subject);
            let (ok, output, stderr) = inputs.cli("6", extra);
            let run = inputs.write(&["6"], extra);
            // The api reads the files before any report: a reader error comes with no output.
            if let Err(Error::Input(error)) = &run.result {
                // A first line that NCBI may fetch as a Seq-id: LOSAT Web's `register` refuses
                // it (NCBI's reader with its data loaders); BLASTX's reader reads the line.
                if query.contains("/seqid3.") || subject.contains("/seqid3.") {
                    assert!(error.contains("data loader"), "{query} {subject}: {error}");
                    continue;
                }
                if ok {
                    differences.push(format!("{query} {subject} {extra:?}: reader {error}"));
                }
                continue;
            }
            let Written {
                result,
                outputs,
                diagnostics,
            } = run;
            if let Some(difference) = ReaderDifference::of(&query) {
                difference.check(ok, &stderr, &result, &diagnostics, &query);
                continue;
            }
            with_rows += usize::from(!output.is_empty());
            let same = match &result {
                Ok(()) => ok && outputs[0] == output && diagnostics == stderr,
                Err(error) => {
                    !ok && diagnostics == stderr[..diagnostics.len().min(stderr.len())]
                        && error_is_in_stderr(error.message(), &stderr)
                }
            };
            if !same {
                differences.push(format!(
                    "{query} {subject} {extra:?}: CLI ok={ok} stderr={:?}; the api {:?} diagnostics={:?}",
                    String::from_utf8_lossy(&stderr),
                    result.as_ref().map_err(|e| e.message()),
                    String::from_utf8_lossy(&diagnostics),
                ));
            }
        }
    }
    assert!(differences.is_empty(), "{}", differences.join("\n"));
    assert!(with_rows > 200, "{with_rows} cases with rows");
}

/// The environment variable that makes `waves_child` run (`queued_batches_as_waves`).
const WAVES_CHILD: &str = "SX_RUN_LOCAL_WAVES_CHILD";

// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:345-357
// ```c
// 		while (master_node.Processing()) {
// 			if (!input.AtEOF()) {
// 			 	if (!master_node.IsFull()) {
// ...
// 					int num_q = input.GetQueryBatch(qb, q_index);
// ```
// With 2+ threads the batch loop queues the reader's batches; `LOSAT_X_BXFLOW=0` (read once per
// process) searches them in waves (`x_flush`) instead of the flow (`x_flow`), so a child process
// of this test binary runs the comparison of several batches with it.
#[test]
fn queued_batches_as_waves() {
    let output = Command::new(std::env::current_exe().expect("test binary"))
        .args([
            "run_local_blastx::waves_child",
            "--exact",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(WAVES_CHILD, "1")
        .env("LOSAT_X_BXFLOW", "0")
        .env_remove("BATCH_SIZE")
        .output()
        .expect("run the child");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success(),
        "{stdout}\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(stdout.contains("waves compared"), "{stdout}");
}

/// The child of `queued_batches_as_waves`: without `WAVES_CHILD` it does nothing.
#[test]
fn waves_child() {
    if std::env::var_os(WAVES_CHILD).is_none() {
        return;
    }
    let batches = Batches::new();
    for threads in ["2", "4"] {
        assert_matches_cli(
            &batches.inputs(),
            &["-num_threads", threads],
            true,
            &format!("waves n{threads}"),
        );
    }
    println!("waves compared");
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1361-1365
// ```c++
//
//     // apply source mods *after* figuring out mol type
//     ITERATE(vector<SLineTextAndLoc>, title_ci, m_CurrentSeqTitles) {
//         ParseTitle(*title_ci, pMessageListener);
//     }
// ```
// NCBI posts a record's title warning when it assembles the record, after the warnings of its
// data lines. NCBI BLAST+ 2.17.0 (`blastx -query q.fna -subject s.faa -outfmt 6`, the inputs
// below) wrote this standard error and the rows that LOSAT writes, from LOSAT Web's records and on
// the command line.
#[test]
fn a_records_title_warning_follows_its_data_warnings() {
    let query = temp_path("blastx_title_order.fna");
    let subject = temp_path("blastx_title_order.faa");
    std::fs::write(
        &query,
        b">q1 x ACGTACGTACGTACGTACGTACGT\nATGCTGTGTCGCATATTTAAG-GACCACTTGCATCAGGAGTGGCGCGCACCC7AATATAATT\nCTGCGGTGCCACGACGTCCTCCCCTCACGGCGTTGTGCCAACGACTCGGACATGTGGACA\n",
    )
    .expect("write query");
    std::fs::write(
        &subject,
        b">s1 t ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN\nMLCRIFKDHLHQEW-RAPNIILRCHDVLPSRRCAN1DSDMWTYDMFNCDELHKGDKAHEATQ\n",
    )
    .expect("write subject");
    let inputs = Inputs::new(&query, &subject);
    let run = inputs.collect(&["6"], &[]);
    let (ok, output, cli_stderr) = inputs.cli("6", &[]);
    let _ = std::fs::remove_file(&query);
    let _ = std::fs::remove_file(&subject);
    assert!(ok && !output.is_empty() && run.reports[0].bytes == output);
    let ncbi = "CFastaReader: Hyphens are invalid and will be ignored around line 2
FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 36
FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?
CFastaReader: Hyphens are invalid and will be ignored around line 2
FASTA-Reader: Ignoring invalid residues at position(s): On line 2: 53
FASTA-Reader: Title ends with at least 20 valid nucleotide characters.  Was the sequence accidentally put in the title line?
";
    assert_eq!(String::from_utf8_lossy(&run.diagnostics), ncbi);
    assert_eq!(String::from_utf8_lossy(&cli_stderr), ncbi);
}

/// Where BLASTX's reader (the CLI) and NCBI's reader (`fasta_reader`, LOSAT Web's records) read a
/// FASTA-input fixture differently. Both end when BLASTX's command line reads with NCBI's reader.
enum ReaderDifference {
    /// A query file of comment lines only: NCBI's `Empty CBlastQueryVector` (exit 3), which
    /// BLASTX's CLI gives; the records hold no query, which every program's `run_local` reports
    /// as NCBI's `Query is Empty!`.
    NoRecords,
}

impl ReaderDifference {
    fn of(query: &str) -> Option<Self> {
        let stem = Path::new(query)
            .file_name()
            .map(|name| {
                name.to_string_lossy()
                    .split('.')
                    .next()
                    .unwrap_or("")
                    .to_string()
            })
            .unwrap_or_default();
        (stem == "file_comment_only").then_some(Self::NoRecords)
    }

    fn check(
        &self,
        cli_ok: bool,
        cli_stderr: &[u8],
        result: &Result<(), Error>,
        diagnostics: &[u8],
        context: &str,
    ) {
        let cli_stderr = String::from_utf8_lossy(cli_stderr);
        assert!(!cli_ok, "{context}: the CLI fails");
        match self {
            Self::NoRecords => {
                assert_eq!(cli_stderr, "BLAST engine error: Empty CBlastQueryVector\n");
                assert!(result.is_ok(), "{context}");
                assert_eq!(diagnostics, b"Warning: [blastx] Query is Empty!\n");
            }
        }
    }
}

/// A query whose first batch ends on its first record only when the batch counts the whole record:
/// 10,050 letters (9,950 in the range `101-20000`), then 200 (100 in the range) and 50 (skipped).
/// Counting the cut letters would put the first two in one batch and leave the third alone, an
/// empty last batch.
fn whole_record_batch_query() -> PathBuf {
    let genome = fixture_sequence("MjeNMV.fasta");
    let mut query = Vec::new();
    for (name, letters) in [
        ("long", &genome[0..10_050]),
        ("middle", &genome[20_000..20_200]),
        ("short", &genome[30_000..30_050]),
    ] {
        query.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in letters.chunks(70) {
            query.extend_from_slice(line);
            query.push(b'\n');
        }
    }
    let path = temp_path("blastx_whole_record_batch.fna");
    std::fs::write(&path, query).expect("write query");
    path
}

/// A subject file of records with the reader's warnings (`batch_late_warn.prot.faa`) and a record
/// of 10 letters between them, after the first warning.
fn short_subject_between_warnings() -> PathBuf {
    let text = std::fs::read_to_string(Inputs::path(
        "tests/fixtures/fasta_input/batch_late_warn.prot.faa",
    ))
    .expect("fixture");
    let at = text.find(">q5").expect("q5");
    let subject = format!("{}>short subject\nMLCRIFKDHL\n{}", &text[..at], &text[at..]);
    let path = temp_path("blastx_short_subject.faa");
    std::fs::write(&path, subject).expect("write subject");
    path
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:138-170
// ```c++
//     TSeqPos size_read = 0;
//
//     while (size_read < GetBatchSize()) {
// ...
//         catch (const exception&) {
//             continue; //SB-2307. ignore well formed, not found accession
//         }
// ...
//             size_read += sequence::GetLength(loc->GetInt().GetId(),
//                                              q->GetScope());
// ```
/// `-query_loc` and `-subject_loc` cut the records as the command line cuts them: the reports of
/// the three formats from one search equal the CLI's of each format, with queries skipped or
/// without letters in the range, subjects without letters, lower-case masks, the minus strand and
/// queries of several batches, at 1, 2 and 4 threads. The half of the test of this name in the crate
/// checks the observer's ranges of the same cases; this half compares the reports and the
/// warnings with the CLI runs, and the api's ranges with its hit records.
#[test]
fn ranges_cut_the_records_as_the_command_line_does() {
    let batches = Batches::new();
    let whole_record_batch = whole_record_batch_query();
    let bslice = || Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa");
    let lowercase_subjects = || {
        Inputs::new(
            "tests/fasta/outfmt0/sml_bx_subj_q.fna",
            "tests/fasta/outfmt0/sml_bx_subj_s.faa",
        )
    };
    // (inputs, ranges and options, whether they have hits)
    let cases = [
        (bslice(), &["-query_loc", "1001-3000"][..], true),
        // The first record ends before the range: it is not searched.
        (bslice(), &["-query_loc", "5200-5600"][..], true),
        // The first record has no letters in the range.
        (bslice(), &["-query_loc", "5001-5600"][..], true),
        (
            bslice(),
            &["-query_loc", "1002-3000", "-strand", "minus"][..],
            true,
        ),
        (bslice(), &["-subject_loc", "10-70"][..], true),
        // Subjects without letters in the range.
        (bslice(), &["-subject_loc", "77-300"][..], true),
        (
            bslice(),
            &["-query_loc", "1001-3000", "-subject_loc", "10-70"][..],
            true,
        ),
        (
            Inputs::new(
                "tests/fasta/outfmt0/x5_bx_q_lc_mix.fna",
                "tests/db_search/bslice.faa",
            ),
            &["-lcase_masking", "-query_loc", "501-2800"][..],
            true,
        ),
        // Several reader batches: the records without letters and the short query are skipped,
        // the batch counts the long record whole.
        (batches.inputs(), &["-query_loc", "50-10050"][..], true),
        // The batch counts the whole record: the first record fills the first batch.
        (
            Inputs::new(&whole_record_batch, "tests/fasta/MjeNMV.faa"),
            &["-query_loc", "101-20000"][..],
            true,
        ),
        // Lower-case runs in the subjects, cut by the subject range.
        (
            lowercase_subjects(),
            &["-lcase_masking", "-subject_loc", "50-150"][..],
            true,
        ),
        (
            lowercase_subjects(),
            &[
                "-lcase_masking",
                "-subject_loc",
                "100-300",
                "-query_loc",
                "1-900",
            ][..],
            true,
        ),
        (
            batches.inputs(),
            &["-query_loc", "50-10050", "-lcase_masking"][..],
            true,
        ),
    ];
    for (inputs, ranges, hits) in &cases {
        for threads in ["1", "2", "4"] {
            let extra = [ranges, &["-num_threads", threads][..]].concat();
            let context = format!("{} {extra:?}", inputs.query);
            assert_matches_cli(inputs, &extra, *hits, &context);
        }
    }
    let _ = std::fs::remove_file(&whole_record_batch);
}

/// The errors of the ranges and the records come where the command line's come, after the same
/// output and warnings: a subject whose interval starts past its end stops the reading of the
/// subjects after its warnings, and a batch whose queries are all skipped is NCBI's `Empty
/// CBlastQueryVector`, the first batch before any report, a later one after the reports of the
/// earlier batches. The half of the test of this name in the crate checks that the entry's error
/// is the typed `NativeError`; this half compares the partial output, the warnings and the message
/// with the CLI runs (`Files::write`: the writers keep what was written before the error).
#[test]
fn errors_of_the_ranges_follow_the_command_line() {
    let batches = Batches::new();
    let short_subject = short_subject_between_warnings();
    let cases = [
        (
            Inputs::new("tests/db_search/bslice_q.fna", &short_subject),
            &["-subject_loc", "30-50"][..],
            "Invalid from coordinate",
        ),
        (
            Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa"),
            &["-query_loc", "6000-7000"][..],
            "Empty CBlastQueryVector",
        ),
        // The long record fills the first batch; the later records are all skipped.
        (
            batches.inputs(),
            &["-query_loc", "10001-10050"][..],
            "Empty CBlastQueryVector",
        ),
    ];
    for (inputs, ranges, message) in &cases {
        for threads in ["1", "2", "4"] {
            let extra = [ranges, &["-num_threads", threads][..]].concat();
            let context = format!("{} {extra:?}", inputs.query);
            let Written {
                result,
                outputs,
                diagnostics,
            } = inputs.write(&FORMATS, &extra);
            let error = result.expect_err(&context);
            assert!(
                matches!(error, Error::Search(_)) && error.message().contains(message),
                "{context}: {error}"
            );
            for (index, outfmt) in FORMATS.iter().enumerate() {
                let (ok, output, stderr) = inputs.cli(outfmt, &extra);
                assert!(!ok, "{context}: the CLI of outfmt {outfmt} fails too");
                assert!(
                    outputs[index] == output,
                    "{context}: outfmt {outfmt} differs from the CLI"
                );
                assert_eq!(
                    String::from_utf8_lossy(
                        &[&diagnostics[..], error.message().as_bytes()].concat()
                    ),
                    String::from_utf8_lossy(&stderr),
                    "{context}: outfmt {outfmt}'s warnings and error"
                );
            }
        }
    }
    let _ = std::fs::remove_file(&short_subject);
}

// NCBI reference: c++/src/app/blast/blast_app_util.cpp:206-210,732-737 and
// c++/src/app/blast/blastx_app.cpp:283-288: NCBI searches each `-subject` sequence on its own
// with BL2SEQ_LEGACY, and converts PRE_FETCH_SEQS_LIMIT with NStr::StringToInt before each
// batch's report (its CStringException names its build's source files), except for the archive
// of outfmt 11, which it writes without the prefetch. The -db cases are the db.bxenv.* rows.
#[test]
fn environment_that_ncbi_reads_around_the_report_is_rejected_where_it_reads_it() {
    let inputs = Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa");
    let run = |outfmt: &str, variable: Option<(&str, &str)>| {
        let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
        command
            .current_dir(crate_dir())
            .env_remove("BATCH_SIZE")
            .env_remove("PRE_FETCH_SEQS_LIMIT")
            .env_remove("BL2SEQ_LEGACY")
            .env(
                "LOSAT_BLASTX_WORKER_LOG",
                std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
                    .unwrap_or_else(|| temp_path("blastx_worker.log").into_os_string()),
            )
            .args([
                "blastx",
                "-query",
                &inputs.query,
                "-subject",
                &inputs.subject,
            ])
            .args(["-outfmt", outfmt]);
        if let Some((name, value)) = variable {
            command.env(name, value);
        }
        command.output().expect("run LOSAT CLI")
    };
    for (variable, value, outfmt, text) in [
        ("PRE_FETCH_SEQS_LIMIT", "abc", "6", "has the value \"abc\""),
        ("PRE_FETCH_SEQS_LIMIT", "", "0", "has the value \"\""),
        (
            "BL2SEQ_LEGACY",
            "1",
            "6",
            "the environment variable BL2SEQ_LEGACY",
        ),
        (
            "BL2SEQ_LEGACY",
            "",
            "11",
            "the environment variable BL2SEQ_LEGACY",
        ),
    ] {
        let output = run(outfmt, Some((variable, value)));
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{variable}: {output:?}");
        assert!(
            output.stdout.is_empty()
                && stderr.contains(text)
                && stderr.contains("not supported by LOSAT's BLASTX"),
            "{variable} {outfmt}: {stderr}"
        );
    }
    // The archive, and an integer in any format, are written as without the variable.
    for (value, outfmt) in [("abc", "11"), ("", "11"), ("5", "6")] {
        let plain = run(outfmt, None);
        let output = run(outfmt, Some(("PRE_FETCH_SEQS_LIMIT", value)));
        assert!(
            plain.status.success() && !plain.stdout.is_empty(),
            "{plain:?}"
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{value} {outfmt}: {output:?}"
        );
        assert!(
            output.stdout == plain.stdout && output.stderr == plain.stderr,
            "{value} {outfmt}"
        );
    }
}

/// The environment variable that makes `environment_child` run: the file it writes.
const ENVIRONMENT_CHILD: &str = "DBID_RUN_LOCAL_ENVIRONMENT_CHILD";

// `run_local` (LOSAT Web, the library) refuses the environment of the CLI's `-subject` search
// (`environment_that_ncbi_reads_around_the_report_is_rejected_where_it_reads_it`), and an
// integer PRE_FETCH_SEQS_LIMIT leaves its outputs as they are without it. The variables are
// read from the process environment, so a child process of this test binary runs `run_local`
// with each of them.
#[test]
fn run_local_refuses_the_environment_as_the_cli() {
    let child = |variable: Option<(&str, &str)>| {
        let out = temp_path("blastx_environment_child.out");
        let mut command = Command::new(std::env::current_exe().expect("test binary"));
        command
            .args([
                "run_local_blastx::environment_child",
                "--exact",
                "--nocapture",
                "--test-threads",
                "1",
            ])
            .env(ENVIRONMENT_CHILD, &out)
            .env_remove("BATCH_SIZE")
            .env_remove("PRE_FETCH_SEQS_LIMIT")
            .env_remove("BL2SEQ_LEGACY");
        if let Some((name, value)) = variable {
            command.env(name, value);
        }
        let output = command.output().expect("run the child");
        assert!(
            output.status.success(),
            "{variable:?}: {}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let written = std::fs::read(&out).expect("the child's file");
        let _ = std::fs::remove_file(&out);
        written
    };
    let plain = child(None);
    assert!(
        plain.starts_with(b"ok\n") && plain.len() > 1000,
        "{plain:?}"
    );
    for (variable, value, text) in [
        ("PRE_FETCH_SEQS_LIMIT", "abc", "has the value \"abc\""),
        ("PRE_FETCH_SEQS_LIMIT", "", "has the value \"\""),
        (
            "BL2SEQ_LEGACY",
            "1",
            "the environment variable BL2SEQ_LEGACY",
        ),
    ] {
        let written = String::from_utf8(child(Some((variable, value)))).expect("text");
        assert!(
            written.starts_with("refused: ")
                && written.contains(text)
                && written.contains("not supported by LOSAT's BLASTX"),
            "{variable}={value:?}: {written}"
        );
    }
    assert!(child(Some(("PRE_FETCH_SEQS_LIMIT", "5"))) == plain);
}

/// The child of `run_local_refuses_the_environment_as_the_cli`: without `ENVIRONMENT_CHILD` it
/// does nothing. It writes `ok` and every format's output and the diagnostics, or `refused:` and
/// the error, after checking that nothing was written.
#[test]
fn environment_child() {
    let Some(out) = std::env::var_os(ENVIRONMENT_CHILD) else {
        return;
    };
    let inputs = Inputs::new("tests/db_search/bslice_q.fna", "tests/db_search/bslice.faa");
    let Written {
        result,
        outputs,
        diagnostics,
    } = inputs.write(&FORMATS, &[]);
    let written = match result {
        Ok(()) => [b"ok\n".to_vec(), outputs.concat(), diagnostics].concat(),
        Err(error) => {
            assert!(
                outputs.iter().all(Vec::is_empty) && diagnostics.is_empty(),
                "{error}"
            );
            format!("refused: {error}\n").into_bytes()
        }
    };
    std::fs::write(out, written).expect("write the child's file");
}
