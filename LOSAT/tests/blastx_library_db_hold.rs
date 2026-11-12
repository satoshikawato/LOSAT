//! BLASTX through the library API with `-db`: a run in which a later query batch can still refuse
//! the taxonomy database (a scientific name NCBI's XML stream cannot write) or a subject of a
//! database with parsed Seq-ids holds every format of the call until its epilog, as the command
//! line holds its BLAST XML, XML2 or JSON document, so that the refusal follows no byte (line NP,
//! series SM; the command line's hold: `tests/outfmt16_blastx.rs`, DECISIONS 21 and 32 of line
//! OF). The side files of outfmt 13 and 14 are part of what is held.
//!
//! The library runs are `examples/library_command_line` (the command line's arguments through
//! `LOSAT::api::run` and `run_to_file`; a failure ends with status 1), each next to the command
//! line's run of the same arguments, and `LOSAT::api::run` with several formats in one call.
//! `BATCH_SIZE=1` makes a query batch of each record; a first record without a hit makes the
//! refusal one of a later batch.
//!
//! No NCBI counterpart for the refusals (LOSAT's own); what a held run writes when nothing is
//! refused is NCBI BLAST+ 2.17.0's output (`tests/fixtures/outfmt16/mixtax_blastx_badname_*`).

use std::path::{Path, PathBuf};
use std::process::Command;

fn engine() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

/// Without the caller's experiment switches (`LOSAT_X_*`) and NCBI's variables that change a
/// search (`tests/search_environment.rs`); `LOSAT_BLASTX_WORKER_LOG` stays, which a build with
/// every feature needs once it starts a pool. `HOME` is an empty directory: no `.ncbirc` names a
/// taxonomy database.
fn clear_environment(command: &mut Command, home: &Path) {
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || name.starts_with("NCBI")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT"
                    | "BATCH_SIZE"
                    | "BL2SEQ_LEGACY"
                    | "BLASTDB"
                    | "CHUNK_SIZE"
                    | "OLD_FSC"
                    | "OVERLAP_CHUNK_SIZE"
                    | "PRE_FETCH_SEQS_LIMIT"
            )
        {
            command.env_remove(name);
        }
    }
    command.env("HOME", home).env("USERPROFILE", home);
}

/// A temporary directory of this test, removed at its end.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat-blastx-library-db-hold-{tag}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("home")).expect("create the temporary directory");
        TempDir(dir)
    }

    fn write(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, text).expect("write an input");
        path
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// An example of the library API, built (with the profile of the binary) when a single test
/// target was asked for.
fn example_binary(name: &str) -> PathBuf {
    let bin_dir = Path::new(env!("CARGO_BIN_EXE_losat"))
        .parent()
        .unwrap()
        .to_path_buf();
    let example = bin_dir
        .join("examples")
        .join(format!("{name}{}", std::env::consts::EXE_SUFFIX));
    if !example.exists() {
        let profile: &[&str] = if bin_dir.ends_with("release") {
            &["--release"]
        } else {
            &["--profile", "test"]
        };
        let status = Command::new(env!("CARGO"))
            .args(["build", "--example", name])
            .args(profile)
            .arg("--target-dir")
            .arg(bin_dir.parent().unwrap())
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .status()
            .expect("build the example");
        assert!(status.success(), "cargo build --example {name}");
    }
    example
}

struct Run {
    status: Option<i32>,
    stdout: Vec<u8>,
    stderr: String,
    /// The files that the run left in its directory, by name, with their bytes: the master file
    /// of outfmt 13 and 14 (`o.json`, `o.xml`) and their side files (`o_<n>.json`, `o_<n>.xml`).
    files: Vec<(String, Vec<u8>)>,
}

impl Run {
    fn side_files(&self) -> Vec<&str> {
        self.files
            .iter()
            .map(|(name, _)| name.as_str())
            .filter(|name| name.starts_with("o_"))
            .collect()
    }

    /// Nothing but the empty `-out` file that both runs create before the search.
    fn wrote_nothing(&self) -> bool {
        self.stdout.is_empty()
            && self
                .files
                .iter()
                .all(|(name, bytes)| name.starts_with("o.") && bytes.is_empty())
    }
}

/// One BLASTX search of `query` against the database `db` (from `LOSAT/`), by the command line
/// (`library` false) or the library example, with `-outfmt outfmt` (13 and 14 write to `o.json`
/// and `o.xml` in a directory of the run's own) at one thread.
fn search(
    dir: &TempDir,
    library: bool,
    db: &str,
    query: &Path,
    outfmt: &str,
    env: &[(&str, &str)],
    more: &[&str],
) -> Run {
    let mut command = if library {
        Command::new(example_binary("library_command_line"))
    } else {
        Command::new(env!("CARGO_BIN_EXE_losat"))
    };
    clear_environment(&mut command, &dir.0.join("home"));
    let out = dir.0.join(if library { "library" } else { "cli" });
    let _ = std::fs::remove_dir_all(&out);
    std::fs::create_dir_all(&out).expect("the run's directory");
    command
        .current_dir(engine())
        .arg("blastx")
        .arg("-query")
        .arg(query)
        .args(["-db", db, "-num_threads", "1", "-outfmt", outfmt])
        .args(more);
    match outfmt {
        "13" => command.arg("-out").arg(out.join("o.json")),
        "14" => command.arg("-out").arg(out.join("o.xml")),
        _ => &mut command,
    };
    for (name, value) in env {
        command.env(name, value);
    }
    let output = command.output().expect("the search runs");
    let mut files: Vec<(String, Vec<u8>)> = std::fs::read_dir(&out)
        .expect("list the directory")
        .map(|entry| {
            let entry = entry.expect("an entry");
            (
                entry.file_name().to_string_lossy().into_owned(),
                std::fs::read(entry.path()).expect("read a file of the run"),
            )
        })
        .collect();
    files.sort();
    Run {
        status: output.status.code(),
        stdout: output.stdout,
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        files,
    }
}

/// The library's run ends as the command line's run of the same arguments: the same bytes on
/// standard output and in the files, and the command line's message (the example writes LOSAT's
/// own error without the `Error: ` before it, and ends every failure with status 1).
fn same_as_command_line(what: &str, library: &Run, cli: &Run) {
    assert_eq!(
        library.status == Some(0),
        cli.status == Some(0),
        "{what}: status {:?}, the command line's {:?}; stderr {}",
        library.status,
        cli.status,
        library.stderr
    );
    assert_eq!(
        String::from_utf8_lossy(&library.stdout),
        String::from_utf8_lossy(&cli.stdout),
        "{what}: standard output"
    );
    assert_eq!(
        library.files.len(),
        cli.files.len(),
        "{what}: files {:?}, the command line's {:?}",
        library
            .files
            .iter()
            .map(|(name, _)| name)
            .collect::<Vec<_>>(),
        cli.files.iter().map(|(name, _)| name).collect::<Vec<_>>()
    );
    for ((name, bytes), (cli_name, cli_bytes)) in library.files.iter().zip(&cli.files) {
        assert_eq!(name, cli_name, "{what}: file names");
        assert_eq!(
            String::from_utf8_lossy(bytes),
            String::from_utf8_lossy(cli_bytes),
            "{what}: {name}"
        );
    }
    let cli_message: String = cli
        .stderr
        .lines()
        .map(|line| match line.strip_prefix("Error: ") {
            Some(own) if !own.starts_with("[blastx]") => own,
            _ => line,
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        library.stderr.trim_end(),
        cli_message.trim_end(),
        "{what}: standard error"
    );
}

/// LOSAT's refusal, after no byte and no side file.
fn refused(what: &str, run: &Run, text: &str) {
    assert_eq!(run.status, Some(1), "{what}: status, stderr {}", run.stderr);
    assert!(run.stderr.contains(text), "{what}: stderr {}", run.stderr);
    assert!(
        run.wrote_nothing(),
        "{what}: wrote {} bytes and the files {:?}",
        run.stdout.len(),
        run.files
            .iter()
            .map(|(name, bytes)| (name, bytes.len()))
            .collect::<Vec<_>>()
    );
}

/// 300 bases that hit nothing in the fixtures' databases: the first query batch of a run whose
/// refusal comes in a later batch.
const NO_HIT: &str =
    "GATTACAGGCTTCAGTCCATAGCGTAACTGGACCTTGAGCATCGGATACTTGCAGCTAAGTCCGATAGGCTTACCGATAGC\
TTAGGCATCAGTCGATACGGATTCAGCTAGGACTTCAGATCGGCATTACGGATACCTGATCGGATTACAGCTAGGTCCATA\
GCTAGGATTCAGCATACGGTTACCAGATCGGTACATTGCAGGCTAACTGGATCAGTCCGATTAGGCATACCGTTAGGCATC\
AGGCTTACAGTCGGATACCTAGGATTCAGGCATACTGGATCAGGCTAACCGTTAGCATCGG";

/// `>blhfill0 no hit` (`NO_HIT`) before the records of the query file `name` (from `LOSAT/`),
/// and `last` after them.
fn after_no_hit(dir: &TempDir, name: &str, last: &str) -> PathBuf {
    let records = std::fs::read_to_string(engine().join(name)).expect("the query fixture");
    let stem = Path::new(name).file_name().unwrap().to_string_lossy();
    dir.write(
        &format!("f.{stem}"),
        &format!(">blhfill0 no hit\n{NO_HIT}\n{}\n{last}", records.trim_end()),
    )
}

const TAXDB_REFUSAL: &str =
    "a taxonomy database scientific name with a character outside the US-ASCII range";
const PARSED_REFUSAL: &str = "by which NCBI BLAST+'s formatters do not find the sequence";
const EACH: (&str, &str) = ("BATCH_SIZE", "1");

/// A scientific name of the taxonomy database that the XML2 and JSON documents cannot write
/// (`tests/seqdb/taxdb/high_bytes_utf8`: taxid 9606), of a sequence that the second query batch
/// is the first to report (`db/sprot_tax`, `qnx_late.fna`: two batches of 10,003 bases): the
/// library, like the command line, writes nothing, and leaves no side file of the first batch.
/// Outfmt 5 asks no scientific name and writes both batches.
#[test]
fn a_taxonomy_refusal_of_a_later_batch_follows_no_byte_of_the_library() {
    let dir = TempDir::new("taxdb");
    let taxdb = engine().join("tests/seqdb/taxdb/high_bytes_utf8");
    let env = [("BLASTDB", taxdb.to_str().unwrap())];
    let db = "tests/fixtures/outfmt16/db/sprot_tax";
    let query = engine().join("tests/fixtures/outfmt16/qnx_late.fna");
    for outfmt in ["13", "14", "15", "16"] {
        let what = format!("outfmt {outfmt}");
        let library = search(&dir, true, db, &query, outfmt, &env, &[]);
        refused(&what, &library, TAXDB_REFUSAL);
        let cli = search(&dir, false, db, &query, outfmt, &env, &[]);
        refused(&format!("{what}, the command line"), &cli, TAXDB_REFUSAL);
        same_as_command_line(&what, &library, &cli);
    }
    let library = search(&dir, true, db, &query, "5", &env, &[]);
    assert_eq!(library.status, Some(0), "outfmt 5: {}", library.stderr);
    assert_eq!(
        String::from_utf8_lossy(&library.stdout)
            .matches("<Iteration>")
            .count(),
        2
    );
    let cli = search(&dir, false, db, &query, "5", &env, &[]);
    same_as_command_line("outfmt 5", &library, &cli);
}

/// The hold changes no byte of a run that nothing refuses: `db/sprot_mixtax` has a sequence, never
/// hit, with the taxid whose name would be refused, so the run holds its output (DECISIONS 32 of
/// line OF) and writes what NCBI BLAST+ 2.17.0 writes (`mixtax_blastx_badname_{15,16}`), through
/// the library too.
#[test]
fn a_held_library_run_that_nothing_refuses_writes_ncbis_bytes() {
    let dir = TempDir::new("mixtax");
    let taxdb = engine().join("tests/seqdb/taxdb/high_bytes_utf8");
    let env = [("BLASTDB", taxdb.to_str().unwrap())];
    let fixtures = engine().join("tests/fixtures/outfmt16");
    for outfmt in ["15", "16"] {
        let frozen = std::fs::read(fixtures.join(format!("mixtax_blastx_badname_{outfmt}.out")))
            .expect("NCBI's bytes");
        // The frozen run names the database as it was given from the fixtures' directory.
        let mut command = Command::new(example_binary("library_command_line"));
        clear_environment(&mut command, &dir.0.join("home"));
        let output = command
            .current_dir(&fixtures)
            .args(["blastx", "-query", "qnx_late.fna", "-db", "db/sprot_mixtax"])
            .args(["-num_threads", "1", "-outfmt", outfmt])
            .env(env[0].0, env[0].1)
            .output()
            .expect("the search runs");
        assert_eq!(
            output.status.code(),
            Some(0),
            "outfmt {outfmt}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&frozen),
            "outfmt {outfmt}"
        );
    }
    // Outfmt 13 and 14: the held side files and the master file are the command line's.
    let db = "tests/fixtures/outfmt16/db/sprot_mixtax";
    let query = fixtures.join("qnx_late.fna");
    for outfmt in ["13", "14"] {
        let what = format!("outfmt {outfmt}");
        let library = search(&dir, true, db, &query, outfmt, &env, &[]);
        assert_eq!(library.status, Some(0), "{what}: {}", library.stderr);
        assert_eq!(library.side_files().len(), 2, "{what}: a file per query");
        let cli = search(&dir, false, db, &query, outfmt, &env, &[]);
        same_as_command_line(&what, &library, &cli);
    }
}

/// A database with parsed Seq-ids (`ids/refuse_v5`): the second query batch reports a sequence
/// whose Seq-align id NCBI's formatters do not find, which LOSAT refuses. The library holds the
/// first batch's reports and its side file, and writes nothing, as the command line.
#[test]
fn a_parsed_seqid_refusal_of_a_later_batch_follows_no_byte_of_the_library() {
    let dir = TempDir::new("parsed");
    let db = "tests/db_search/ids/refuse_v5";
    let query = after_no_hit(&dir, "tests/db_search/ids/refuse_disp.fna", "");
    for outfmt in ["5", "13", "14", "15", "16"] {
        let what = format!("outfmt {outfmt}");
        let more = ["-evalue", "1e-20"];
        let library = search(&dir, true, db, &query, outfmt, &[EACH], &more);
        refused(&what, &library, PARSED_REFUSAL);
        let cli = search(&dir, false, db, &query, outfmt, &[EACH], &more);
        refused(&format!("{what}, the command line"), &cli, PARSED_REFUSAL);
        same_as_command_line(&what, &library, &cli);
    }
}

/// A held run that one of NCBI's own errors stops writes what NCBI wrote before it: the reports
/// of the batches before (a last batch without letters fails NCBI's set-up), and with
/// `BATCH_SIZE=0` the prolog's bytes before `Empty CBlastQueryVector`. The database has
/// parsed Seq-ids (`ids/prot_v5`), so every such run is held.
#[test]
fn a_held_library_run_that_ncbis_error_stops_writes_the_batches_before() {
    let dir = TempDir::new("ncbi-error");
    let db = "tests/db_search/ids/prot_v5";
    let query = after_no_hit(
        &dir,
        "tests/db_search/ids/ids_q.fna",
        ">blhlast9 no letters\n",
    );
    for outfmt in ["5", "13", "14", "15", "16"] {
        let what = format!("no letters, outfmt {outfmt}");
        let library = search(&dir, true, db, &query, outfmt, &[EACH], &[]);
        assert!(
            library.stderr.contains("Sequence contains no data"),
            "{what}: {}",
            library.stderr
        );
        assert!(!library.wrote_nothing(), "{what}: the batches before");
        if matches!(outfmt, "13" | "14") {
            assert!(!library.side_files().is_empty(), "{what}: side files");
        }
        let cli = search(&dir, false, db, &query, outfmt, &[EACH], &[]);
        assert_eq!(cli.status, Some(3), "{what}: {}", cli.stderr);
        same_as_command_line(&what, &library, &cli);

        let what = format!("BATCH_SIZE=0, outfmt {outfmt}");
        let zero = [("BATCH_SIZE", "0")];
        let library = search(&dir, true, db, &query, outfmt, &zero, &[]);
        assert!(
            library.stderr.contains("Empty CBlastQueryVector"),
            "{what}: {}",
            library.stderr
        );
        let cli = search(&dir, false, db, &query, outfmt, &zero, &[]);
        assert_eq!(cli.status, Some(3), "{what}: {}", cli.stderr);
        same_as_command_line(&what, &library, &cli);
        // NCBI writes the start of an XML2 or JSON document in its prolog, and that of BLAST XML
        // with the first query's report, which this run never reaches.
        assert_eq!(
            !library.stdout.is_empty(),
            matches!(outfmt, "15" | "16"),
            "{what}: the start of the document"
        );
    }
}

/// One call of `LOSAT::api::run` with several formats: when a later batch refuses a subject with
/// parsed Seq-ids, no format has written a byte (the tabular rows of the first batch are held
/// with the document), whichever format comes first. The first query (more than 10,002 bases, a
/// batch of its own under the default batch size) has no hit; the second is refused.
#[test]
fn a_library_call_with_several_formats_holds_them_all() {
    use LOSAT::api::{read_inputs, run, validate, Inputs, Program};
    // The call runs in this process: with one of the variables that the other tests clear for
    // their runs, its batches or its database's taxonomy are not those of this test.
    if let Some((name, _)) = std::env::vars_os().find(|(name, _)| {
        let name = name.to_string_lossy();
        name.starts_with("NCBI")
            || matches!(&*name, "BATCH_SIZE" | "BLASTDB" | "PRE_FETCH_SEQS_LIMIT")
    }) {
        eprintln!("skipped: {} is set", name.to_string_lossy());
        return;
    }
    let mut filler = String::from(NO_HIT);
    while filler.len() < 10_003 {
        filler.push_str(NO_HIT);
    }
    let refused_record =
        std::fs::read_to_string(engine().join("tests/db_search/ids/refuse_disp.fna"))
            .expect("the query fixture");
    let query = format!(">blhfill0 no hit\n{filler}\n{refused_record}");
    let db = engine().join("tests/db_search/ids/refuse_v5");
    let search = validate(Program::Blastx, &["-evalue", "1e-20"]).expect("blastx");
    let (queries, _) = read_inputs(&search, query.as_bytes(), b"").expect("the records");
    for formats in [["6", "16"], ["16", "6"], ["5", "7"]] {
        let mut sinks: Vec<Vec<u8>> = vec![Vec::new(); formats.len()];
        let mut diagnostics = Vec::new();
        let result = {
            let mut reports: Vec<(&str, &mut (dyn std::io::Write + Send))> = formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(&format, sink)| (format, sink as &mut (dyn std::io::Write + Send)))
                .collect();
            run(
                &search,
                Inputs::database(&queries, db.to_str().unwrap()),
                &mut reports,
                &mut diagnostics,
            )
        };
        let message = result.expect_err("the refusal").message().to_string();
        assert!(message.contains(PARSED_REFUSAL), "{formats:?}: {message}");
        assert!(
            sinks.iter().all(Vec::is_empty),
            "{formats:?}: wrote {:?} bytes",
            sinks.iter().map(Vec::len).collect::<Vec<_>>()
        );
    }
}
