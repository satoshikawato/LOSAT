#![allow(warnings, clippy::all)]

//! Release-facing CLI smoke tests.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_path(stem: &str, extension: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "losat_cli_{stem}_{}_{}.{}",
        std::process::id(),
        nanos,
        extension
    ))
}

fn losat_binary() -> PathBuf {
    if let Some(path) = option_env!("CARGO_BIN_EXE_losat") {
        return PathBuf::from(path);
    }

    let mut path = std::env::current_exe().expect("current test executable");
    path.pop();
    if path.file_name().is_some_and(|name| name == "deps") {
        path.pop();
    }
    path.push(format!("losat{}", std::env::consts::EXE_SUFFIX));
    path
}

fn clean_losat_command() -> Command {
    let mut command = Command::new(losat_binary());
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
            command.env_remove(name);
        }
    }
    command
}

fn run_losat(args: &[&str]) -> Output {
    clean_losat_command()
        .args(args)
        .output()
        .expect("run LOSAT")
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blastn_args.cpp:48-53
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blastp_args.cpp:47-52
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/tblastx_args.cpp:47-52
// ```c
// m_ClientId = string(kProgram) + " " + CBlastVersion().Print();
// m_ClientId = kProgram + " " + CBlastVersion().Print();
// ```
#[test]
fn cli_help_and_version_are_available_without_inputs() {
    for args in [
        vec!["--help"],
        vec!["blastn", "--help"],
        vec!["tblastx", "--help"],
        vec!["blastp", "--help"],
    ] {
        let output = run_losat(&args);
        assert!(
            output.status.success(),
            "{args:?} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "{args:?} should not write help to stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            String::from_utf8_lossy(&output.stdout).contains("USAGE\n"),
            "{args:?} help should open with NCBI's USAGE heading"
        );
    }

    let output = run_losat(&["--version"]);
    assert!(
        output.status.success(),
        "--version failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        format!("losat {}", env!("CARGO_PKG_VERSION"))
    );
    assert!(
        output.stderr.is_empty(),
        "--version should not write stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

// No NCBI counterpart: LOSAT_STARTUP_TRACE is a LOSAT diagnostic variable; it does not change any
// value NCBI computes.
// CI's `--all-features` run checks the feature side of this test (the trace line, no notice); a
// default `cargo test` checks the other (the notice line, no trace). `ci_fast_regressions.py`
// checks the notice on the default release binary.
#[test]
fn startup_trace_needs_the_diagnostics_feature() {
    let plain = clean_losat_command()
        .arg("--version")
        .output()
        .expect("run LOSAT");
    let traced = clean_losat_command()
        .arg("--version")
        .env("LOSAT_STARTUP_TRACE", "1")
        .output()
        .expect("run LOSAT");
    assert!(plain.status.success() && traced.status.success());
    assert_eq!(plain.stdout, traced.stdout);
    assert!(
        plain.stderr.is_empty(),
        "--version should not write stderr: {}",
        String::from_utf8_lossy(&plain.stderr)
    );
    let expected = if cfg!(feature = "diagnostics") {
        "[startup] enter main\n"
    } else {
        "LOSAT_STARTUP_TRACE is ignored: this losat was built without the diagnostics feature \
         (docs/agents/diagnostics.md, \"Diagnostic build\")\n"
    };
    assert_eq!(String::from_utf8_lossy(&traced.stderr), expected);
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
// ```c
//     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
//         (*arg)->ExtractAlgorithmOptions(args, opts);
//     }
//
//     m_IsUngapped = !opts.GetGappedMode();
//     try { retval->Validate(); }
// ```
// NCBI's handlers read the subjects and open the query before the options are checked, so a
// missing file is reported first; LOSAT's rejection of an option that NCBI runs comes where
// NCBI checks the options.
#[test]
fn unsupported_blastp_option_fails_after_ncbis_file_checks() {
    let missing_query = temp_path("missing_query", "faa");
    let missing_subject = temp_path("missing_subject", "faa");
    let output = clean_losat_command()
        .args(["blastp", "-query"])
        .arg(&missing_query)
        .arg("-subject")
        .arg(&missing_subject)
        .args(["-word_size", "7"])
        .output()
        .expect("run unsupported blastp");
    assert!(!output.status.success(), "a missing subject should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("not supported by LOSAT's BLASTP"),
        "the subject is read before the options are checked: {stderr}"
    );

    let query = temp_path("word7_query", "faa");
    let subject = temp_path("word7_subject", "faa");
    std::fs::write(&query, ">q\nMKTAYIAKQRQISFVKSHFSRQ\n").unwrap();
    std::fs::write(&subject, ">s\nMKTAYIAKQRQISFVKSHFSRQ\n").unwrap();
    let output = clean_losat_command()
        .args(["blastp", "-query"])
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .args(["-word_size", "7"])
        .output()
        .expect("run unsupported blastp");
    let _ = std::fs::remove_file(&query);
    let _ = std::fs::remove_file(&subject);
    assert!(!output.status.success(), "unsupported blastp should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(
            "-word_size 7 (the compressed lookup table) is not supported by LOSAT's BLASTP"
        ),
        "unsupported error should name program, option, and reason: {stderr}"
    );
    assert!(output.stdout.is_empty(), "nothing is written: {stderr}");
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
// ```c
// CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
// CRef<CBlastInput> input(new CBlastInput(fasta));
// CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
// sequences = input->GetAllSeqs(*scope);
// ```
#[test]
fn missing_input_file_reports_ncbi_error() {
    let missing_query = temp_path("missing_query", "fa");
    let missing_subject = temp_path("missing_subject", "fa");
    let output = clean_losat_command()
        .args(["blastn", "-outfmt", "6", "-query"])
        .arg(&missing_query)
        .arg("-subject")
        .arg(&missing_subject)
        .output()
        .expect("run blastn with missing input");

    assert!(!output.status.success(), "missing input should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    // NCBI opens the subject first (blastn_args.cpp:63-70), with NCBI's message
    // (ncbiargs.cpp:615-619).
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr,
        format!(
            "Command line argument error: Argument \"subject\". File is not accessible:  `{}'\n",
            missing_subject.display()
        ),
        "missing input error should be NCBI's: {stderr}"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:751-758
// ```c++
//     if (bad >= good / 3  &&
//         (len_to_check > 3  ||  good == 0  ||  bad > good))
//     {
//         FASTA_ERROR( LineNumber(),
//             "CFastaReader: Near line " << LineNumber()
//             << ", there's a line that doesn't look like plausible data, "
//             "but it's not marked as defline or comment.",
//             CObjReaderParseException::eFormat);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
// ```c++
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
// A comment line before the first defline is skipped; a first data line that does not
// look like sequence data stops BLASTN with NCBI's reader error.
#[test]
fn malformed_fasta_reports_parse_error() {
    let query = temp_path("malformed_query", "fa");
    let subject = temp_path("valid_subject", "fa");
    fs::write(&query, "!a comment line\n@@@@ACGT\n").expect("write malformed FASTA");
    fs::write(&subject, ">s\nACGTACGTACGT\n").expect("write subject FASTA");

    let output = clean_losat_command()
        .args(["blastn", "-outfmt", "6", "-query"])
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .output()
        .expect("run blastn with malformed FASTA");

    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);

    assert_eq!(output.status.code(), Some(1), "malformed FASTA should fail");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n"
    );
}

/// Gzip members of `members`, each with bgzip's extra field, concatenated as `bgzip` and
/// `cat a.gz b.gz` write them.
fn gzip_members(members: &[&[u8]]) -> Vec<u8> {
    use std::io::Write;
    let mut bytes = Vec::new();
    for member in members {
        let mut encoder = flate2::GzBuilder::new()
            .extra(b"BC\x02\x00\x00\x00".to_vec())
            .write(&mut bytes, flate2::Compression::default());
        encoder.write_all(member).unwrap();
        encoder.finish().unwrap();
    }
    bytes
}

const GZ_QUERY: &str = ">q1 first\nACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGTACCTTAGCAGTCGATCGACGTTGCAAGGCTTAACCGTATGCG\n>q2 second\nGGCATGCATGCATTTAGGCCAAGTACGATCGATCGACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGG\n";
const GZ_SUBJECT: &str = ">s1 the subject\nTTTTACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGTACCTTAGCAGTCGATCGACGTTGCAAGGCTTAACCGTATGCGAAAAGGCATGCATGCATTTAGGCCAAGTACGATCGATCGACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGCCCC\n";

fn blastn_output(query: &std::path::Path, subject: &std::path::Path, outfmt: &str) -> Output {
    clean_losat_command()
        .args(["blastn", "-outfmt", outfmt, "-query"])
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .output()
        .expect("run blastn")
}

// No NCBI counterpart (line GZ, DECISIONS 1 and 3; `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
// exception 8): a `-query` or `-subject` named `.gz` reads as the plain file, in any case of
// the suffix; the name the user gave is the one the reports show.
#[test]
fn gzip_query_and_subject_read_as_the_plain_files() {
    let query = temp_path("gz_query", "fa");
    let subject = temp_path("gz_subject", "fa");
    let query_gz = temp_path("gz_query", "fa.gz");
    let subject_gz = temp_path("gz_subject", "fna.GZ");
    fs::write(&query, GZ_QUERY).unwrap();
    fs::write(&subject, GZ_SUBJECT).unwrap();
    fs::write(&query_gz, gzip_members(&[GZ_QUERY.as_bytes()])).unwrap();
    fs::write(&subject_gz, gzip_members(&[GZ_SUBJECT.as_bytes()])).unwrap();

    let plain = blastn_output(&query, &subject, "6");
    let inflated = blastn_output(&query_gz, &subject_gz, "6");
    let plain_7 = blastn_output(&query, &subject, "7");
    let inflated_7 = blastn_output(&query_gz, &subject_gz, "7");

    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);
    let _ = fs::remove_file(&query_gz);
    let _ = fs::remove_file(&subject_gz);

    assert_eq!(
        plain.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&plain.stderr)
    );
    assert_eq!(
        inflated.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&inflated.stderr)
    );
    assert!(plain.stdout.len() > 0, "the plain run has rows");
    assert_eq!(inflated.stdout, plain.stdout);
    assert_eq!(inflated.stderr, plain.stderr);
    let plain_7 = String::from_utf8_lossy(&plain_7.stdout).into_owned();
    let inflated_7 = String::from_utf8_lossy(&inflated_7.stdout).into_owned();
    let input_line = |text: &str, path: &std::path::Path| {
        let line = format!(
            "# Database: User specified sequence set (Input: {})",
            path.display()
        );
        assert!(text.contains(&line), "{text}");
        text.replace(
            &line,
            "# Database: User specified sequence set (Input: SUBJECT)",
        )
    };
    assert_eq!(
        input_line(&inflated_7, &subject_gz),
        input_line(&plain_7, &subject)
    );
}

// No NCBI counterpart (line GZ, DECISIONS 4): the members of a concatenated gzip file (bgzip's
// blocks, each with its extra field, and its empty end block) are all read; a reader error
// inside the inflated file keeps the line number of the plain file.
#[test]
fn concatenated_gzip_members_are_read_whole_with_the_plain_line_numbers() {
    let subject = temp_path("members_subject", "fa");
    let query_gz = temp_path("members_query", "fa.gz");
    fs::write(&subject, GZ_SUBJECT).unwrap();
    let (first, second) = GZ_QUERY.split_at(GZ_QUERY.find(">q2").unwrap());
    fs::write(
        &query_gz,
        gzip_members(&[first.as_bytes(), second.as_bytes(), b""]),
    )
    .unwrap();
    let both = blastn_output(&query_gz, &subject, "6");
    let bad_gz = temp_path("members_bad", "fa.gz");
    fs::write(
        &bad_gz,
        gzip_members(&[b"!a comment line\n", b"@@@@ACGT\n"]),
    )
    .unwrap();
    let bad = blastn_output(&bad_gz, &subject, "6");
    let _ = fs::remove_file(&subject);
    let _ = fs::remove_file(&query_gz);
    let _ = fs::remove_file(&bad_gz);

    assert_eq!(
        both.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&both.stderr)
    );
    let rows = String::from_utf8_lossy(&both.stdout).into_owned();
    assert!(rows.lines().any(|row| row.starts_with("q1\t")), "{rows}");
    assert!(rows.lines().any(|row| row.starts_with("q2\t")), "{rows}");
    assert_eq!(bad.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&bad.stderr),
        "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n"
    );
}

// No NCBI counterpart (line GZ, DECISIONS 5): a `.gz` that is not gzip data is refused before
// the search, with exit status 1, and `-` is read as it comes (never inflated).
#[test]
fn a_gz_name_without_gzip_data_is_refused_before_the_search() {
    let query_gz = temp_path("text_query", "fa.gz");
    let subject = temp_path("text_subject", "fa");
    fs::write(&query_gz, GZ_QUERY).unwrap();
    fs::write(&subject, GZ_SUBJECT).unwrap();
    let refused = blastn_output(&query_gz, &subject, "6");
    let _ = fs::remove_file(&query_gz);
    let _ = fs::remove_file(&subject);
    assert_eq!(refused.status.code(), Some(1));
    assert_eq!(refused.stdout, b"");
    let stderr = String::from_utf8_lossy(&refused.stderr).into_owned();
    assert!(stderr.starts_with("Error: -query: '"), "{stderr}");
    assert!(
        stderr.contains("' is named .gz but is not gzip data ("),
        "{stderr}"
    );
    assert!(
        stderr
            .trim_end()
            .ends_with("); LOSAT reads a file named .gz through gzip"),
        "{stderr}"
    );
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050 (an `-in` name that
// ends in `.gz`, in any case, is read through a gzip stream; the title defaults to the name as
// given). The database built from the inflated file is the one built from the plain file;
// `.bz2` and `.zst` stay refused (line GZ, DECISIONS 2).
#[test]
fn makeblastdb_reads_a_gz_input_as_makeblastdb_does() {
    let dir = temp_path("makedb_gz", "dir");
    fs::create_dir_all(&dir).unwrap();
    let plain = dir.join("s.fna");
    let gz = dir.join("s.fna.GZ");
    let query = dir.join("q.fna");
    fs::write(&plain, GZ_SUBJECT).unwrap();
    fs::write(&gz, gzip_members(&[GZ_SUBJECT.as_bytes()])).unwrap();
    fs::write(&query, GZ_QUERY).unwrap();
    let build = |input: &std::path::Path, out: &str| {
        clean_losat_command()
            .current_dir(&dir)
            .args(["makeblastdb", "-dbtype", "nucl", "-in"])
            .arg(input.file_name().unwrap())
            .args(["-out", out])
            .output()
            .expect("run makeblastdb")
    };
    let from_plain = build(&plain, "db_plain");
    let from_gz = build(&gz, "db_gz");
    assert_eq!(
        from_plain.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&from_plain.stderr)
    );
    assert_eq!(
        from_gz.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&from_gz.stderr)
    );
    assert_eq!(from_gz.stderr, from_plain.stderr);
    let log = String::from_utf8_lossy(&from_gz.stdout).into_owned();
    assert!(log.contains("New DB title:  s.fna.GZ\n"), "{log}");
    assert!(log.contains("added 1 sequences in"), "{log}");
    for extension in ["nsq", "nhr"] {
        assert_eq!(
            fs::read(dir.join(format!("db_gz.{extension}"))).unwrap(),
            fs::read(dir.join(format!("db_plain.{extension}"))).unwrap(),
            "{extension}"
        );
    }
    let search = |db: &str| {
        clean_losat_command()
            .current_dir(&dir)
            .args(["blastn", "-outfmt", "6", "-query", "q.fna", "-db", db])
            .output()
            .expect("run blastn -db")
    };
    let plain_rows = search("db_plain");
    let gz_rows = search("db_gz");
    assert_eq!(
        plain_rows.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&plain_rows.stderr)
    );
    assert!(plain_rows.stdout.len() > 0);
    assert_eq!(gz_rows.stdout, plain_rows.stdout);

    let text_gz = dir.join("text.fna.gz");
    fs::write(&text_gz, GZ_SUBJECT).unwrap();
    let refused = build(&text_gz, "db_text");
    assert_eq!(refused.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&refused.stderr).into_owned();
    assert!(
        stderr
            .contains("losat makeblastdb: -in: 'text.fna.gz' is named .gz but is not gzip data ("),
        "{stderr}"
    );
    assert!(
        !dir.join("db_text.nsq").exists(),
        "the failed build leaves no files"
    );

    let bz2 = dir.join("s.fna.bz2");
    fs::write(&bz2, GZ_SUBJECT).unwrap();
    let refused = build(&bz2, "db_bz2");
    assert_eq!(refused.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&refused.stderr).into_owned();
    assert!(stderr.contains("s.fna.bz2: makeblastdb reads a file of this name as compressed; losat makeblastdb reads uncompressed FASTA and gzip (.gz) only"), "{stderr}");
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:140-152
// ```c++
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
// ```
// NCBI reads the queries one batch at a time: a reader error in a later batch comes after
// the reports of the batches before, without the epilog (oracle BI b_q_b2_bad_7), and in
// the first batch before any report.
#[test]
fn a_reader_error_in_a_later_query_batch_follows_the_earlier_reports() {
    let query = temp_path("later_batch_query", "fa");
    let subject = temp_path("later_batch_subject", "fa");
    let sequence = "ACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGTACCTTAGCAGTCGATCG".repeat(5);
    let mut lines = String::from(">q1\n");
    for chunk in sequence.as_bytes().chunks(60) {
        lines.push_str(std::str::from_utf8(chunk).unwrap());
        lines.push('\n');
    }
    lines.push_str(">q2\n@@@@\n");
    fs::write(&query, &lines).expect("write query FASTA");
    fs::write(&subject, format!(">s1\n{sequence}\n")).expect("write subject FASTA");
    let run = |batch_size: Option<&str>| {
        let mut command = clean_losat_command();
        command.args(["blastn", "-outfmt", "7", "-query"]);
        command.arg(&query).arg("-subject").arg(&subject);
        match batch_size {
            Some(size) => command.env("BATCH_SIZE", size),
            None => command.env_remove("BATCH_SIZE"),
        };
        command.output().expect("run blastn")
    };
    let later = run(Some("300"));
    let first = run(None);
    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);
    let error = "BLAST query error: CFastaReader: Near line 8, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n";
    for output in [&later, &first] {
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(String::from_utf8_lossy(&output.stderr), error);
    }
    let report = String::from_utf8_lossy(&later.stdout);
    assert!(
        report.starts_with("# BLASTN 2.17.0+\n# Query: q1\n")
            && !report.contains("# BLAST processed"),
        "{report}"
    );
    assert!(first.stdout.is_empty(), "{:?}", first.stdout);
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
// ```c
// #if _BLAST_DEBUG
// arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
//                  true);
// #endif /* _BLAST_DEBUG */
// ```
#[test]
fn normal_tblastx_run_keeps_stderr_clean_without_debug_env() {
    let query = temp_path("tblastx_query", "fa");
    let subject = temp_path("tblastx_subject", "fa");
    let out = temp_path("tblastx_out", "txt");
    fs::write(&query, ">q\nATGATGATGATGATGATGATGATGATGATG\n").expect("write query FASTA");
    fs::write(&subject, ">s\nATGATGATGATGATGATGATGATGATGATG\n").expect("write subject FASTA");

    let output = clean_losat_command()
        .args(["tblastx", "-query"])
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .args(["-seg", "no", "-outfmt", "6", "-num_threads", "1", "-out"])
        .arg(&out)
        .output()
        .expect("run tblastx");

    let output_text = fs::read_to_string(&out).unwrap_or_default();
    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);
    let _ = fs::remove_file(&out);

    assert!(
        output.status.success(),
        "tblastx smoke failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stderr.is_empty(),
        "normal tblastx run should keep stderr clean: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !output_text.is_empty(),
        "tblastx smoke input should produce at least one hit"
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-977
// ```c++
//             } else if ( NStr::strcmp(argv[i], s_ArgVersion) == 0 ) {
//                 delete[] v;
//                 // Print VERSION
//                 cout << GetFullVersion().Print( appname, CVersionAPI::fVersionInfo | CVersionAPI::fPackageShort );
//                 diag_context.DiscardMessages();
//                 return 0;
// ```
// NCBI reads `-version` before any other argument: anywhere before `--`, in the value place
// of another option, and next to arguments that are invalid.
#[test]
fn ncbi_version_is_two_lines_for_every_program() {
    let package = format!(
        " Package: LOSAT {}, output compatible with blast 2.17.0\n",
        env!("CARGO_PKG_VERSION")
    );
    for (name, printed) in [
        ("blastn", "blastn"),
        ("blastp", "blastp"),
        ("blastx", "blastx"),
        ("tblastn", "tblastn"),
        ("tblastx", "tblastx"),
        ("psiblast", "psiblast"),
        ("makeblastdb", "makeblastdb"),
    ] {
        let expected = format!("{printed}: 2.17.0+\n{package}");
        for words in [
            vec![name, "-version"],
            vec![name, "-out", "-version"],
            vec![name, "-no_such_option", "-version"],
            vec![name, "-evalue", "x", "-version", "-bogus"],
        ] {
            let output = run_losat(&words);
            assert!(output.status.success(), "{words:?}: {:?}", output.status);
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                expected,
                "{words:?}"
            );
            assert!(output.stderr.is_empty(), "{words:?}");
        }
        // After `--` the word is a positional argument, which no program has.
        let output = run_losat(&[name, "--", "-version"]);
        assert!(!output.status.success(), "{name} -- -version");
    }
    // The top level `--version` is clap's.
    let output = run_losat(&["--version"]);
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        format!("losat {}", env!("CARGO_PKG_VERSION"))
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:76-77
// ```c++
// static const char* s_AutoHelp     = "h";
// static const char* s_AutoHelpFull = "help";
// ```
#[test]
fn short_help_h_succeeds_for_every_program() {
    for name in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "makeblastdb",
    ] {
        for words in [vec![name, "-h"]] {
            let output = run_losat(&words);
            assert!(output.status.success(), "{words:?}: {:?}", output.status);
            assert!(output.stderr.is_empty(), "{words:?}");
            let text = String::from_utf8_lossy(&output.stdout);
            assert!(text.starts_with("USAGE\n"), "{words:?}: {text}");
            assert!(text.contains("\nDESCRIPTION\n"), "{words:?}: {text}");
            assert!(!text.contains("error:"), "{words:?}: {text}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:217-219
// ```c++
//     arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
//                   "Application to create BLAST databases, version "
//                   + CBlastVersion().Print());
// ```
// `losat makeblastdb` writes the database; its old name `makedb` is an unknown command.
#[test]
fn makeblastdb_writes_a_database_and_makedb_is_unknown() {
    let input = temp_path("mkdb_in", "faa");
    fs::write(
        &input,
        ">a first\nMKTAYIAKQRQISFVKSHFSRQ\n>b second\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let base = temp_path("makeblastdb", "db");
    let (input_text, base_text) = (input.to_str().unwrap(), base.to_str().unwrap());
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
    ]);
    assert!(output.status.success(), "{:?}", output);
    for extension in ["phr", "psq", "pin"] {
        let path = format!("{base_text}.{extension}");
        assert!(!fs::read(&path).unwrap().is_empty(), ".{extension}");
        fs::remove_file(path).unwrap();
    }
    // `makedb` is not a command, with or without a valid command line.
    for words in [
        vec![
            "makedb", "-in", input_text, "-dbtype", "prot", "-out", base_text,
        ],
        vec!["makedb", "-help"],
        vec!["makedb", "-version"],
    ] {
        let output = run_losat(&words);
        assert!(!output.status.success(), "{words:?}");
        assert!(output.stdout.is_empty(), "{words:?}");
        assert!(
            !std::path::Path::new(&format!("{base_text}.psq")).exists(),
            "{words:?}"
        );
    }
    fs::remove_file(&input).unwrap();
    // `-help` works, and the top level lists makeblastdb and not makedb.
    let output = run_losat(&["makeblastdb", "-help"]);
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stdout).contains("-dbtype <molecule_type>"));
    let top = String::from_utf8_lossy(&run_losat(&["-help"]).stdout).into_owned();
    assert!(
        top.contains("makeblastdb") && !top.contains("makedb"),
        "{top}"
    );
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1203-1209
// ```c++
//     bool parse_seqids = x_ShouldParseSeqIds();
//     bool hash_index = args["hash_index"];
//     bool use_gi_mask = args["gi_mask"];
//
//     CWriteDB::TIndexType indexing = CWriteDB::eNoIndex;
//     indexing |= (hash_index ? CWriteDB::eAddHash : 0);
//     indexing |= (parse_seqids ? CWriteDB::eFullIndex : 0);
// ```
// `-parse_seqids`, `-hash_index` and `-taxid_map` on the command line: the id files of version 5
// (`.pos`, `.pog`, the hash index) and of version 4 (`.psi`/`.psd`), the argument rules (exit
// status 2 from LOSAT's parser), a version 5 input makeblastdb aborts on (exit status 1), and
// `[BLAST] MAX_SEQID_LENGTH` refused with `-parse_seqids` (makedb-v5 DECISIONS 12). The bytes are
// checked against makeblastdb in `writedb::tests::parse_seqids_databases_are_makeblastdb_s`.
#[test]
fn makeblastdb_parse_seqids_hash_index_and_taxid_map() {
    let input = temp_path("mkdb_parse_in", "faa");
    fs::write(
        &input,
        ">ref|NP_000001.1| first\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|second b\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let map = temp_path("mkdb_parse_map", "txt");
    fs::write(&map, "NP_000001.1 9606\nsecond 10090\n").unwrap();
    let base = temp_path("mkdb_parse", "db");
    let (input_text, map_text, base_text) = (
        input.to_str().unwrap(),
        map.to_str().unwrap(),
        base.to_str().unwrap(),
    );
    let remove = |extensions: &[&str]| {
        for extension in extensions {
            let _ = fs::remove_file(format!("{base_text}.{extension}"));
        }
    };
    const V5: [&str; 12] = [
        "pdb", "phd", "phi", "phr", "pin", "pjs", "pog", "pos", "pot", "psq", "ptf", "pto",
    ];
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
        "-hash_index",
        "-taxid_map",
        map_text,
    ]);
    assert!(output.status.success(), "{output:?}");
    for extension in V5 {
        let path = format!("{base_text}.{extension}");
        assert!(!fs::read(&path).unwrap().is_empty(), ".{extension}");
    }
    let metadata = fs::read_to_string(format!("{base_text}.pjs")).unwrap();
    assert!(metadata.contains("\"number-of-taxids\": 2"), "{metadata}");
    remove(&V5);
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
        "-blastdb_version",
        "4",
    ]);
    assert!(output.status.success(), "{output:?}");
    for extension in ["psi", "psd", "pog"] {
        assert!(
            std::path::Path::new(&format!("{base_text}.{extension}")).exists(),
            ".{extension}"
        );
    }
    remove(&["phr", "pin", "pjs", "pog", "psd", "psi", "psq"]);
    // `-taxid_map` requires `-parse_seqids` and excludes `-taxid`: LOSAT's parser, exit 2.
    for words in [
        vec!["-taxid_map", map_text],
        vec!["-parse_seqids", "-taxid_map", map_text, "-taxid", "1"],
    ] {
        let mut argv = vec![
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
        ];
        argv.extend(words.iter().copied());
        let output = run_losat(&argv);
        assert_eq!(output.status.code(), Some(2), "{words:?}: {output:?}");
        assert!(!std::path::Path::new(&format!("{base_text}.psq")).exists());
    }
    // A GI-only sequence between sequences with other ids: makeblastdb 2.17.0 aborts on version 5.
    fs::write(
        &input,
        ">ref|NP_1.1| a\nMKTAYIAKQR\n>gi|42 b\nMKTAYIAKQR\n>ref|NP_2.1| c\nMKTAYIAKQR\n",
    )
    .unwrap();
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
    ]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(String::from_utf8_lossy(&output.stderr).starts_with(
        "BLAST Database creation error: Input id list not in ascending oid order: the sequence gi|42 (OID 1)"
    ));
    assert!(!std::path::Path::new(&format!("{base_text}.psq")).exists());
    remove(&V5);
    // `[BLAST] MAX_SEQID_LENGTH` is refused with `-parse_seqids` only.
    let refused = clean_losat_command()
        .env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "100")
        .args([
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
            "-parse_seqids",
            "-blastdb_version",
            "4",
        ])
        .output()
        .unwrap();
    assert!(!refused.status.success(), "{refused:?}");
    assert!(
        String::from_utf8_lossy(&refused.stderr).contains("[BLAST] MAX_SEQID_LENGTH (100)"),
        "{refused:?}"
    );
    assert!(!std::path::Path::new(&format!("{base_text}.psq")).exists());
    let accepted = clean_losat_command()
        .env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "100")
        .args([
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
        ])
        .output()
        .unwrap();
    assert!(accepted.status.success(), "{accepted:?}");
    remove(&V5);
    // A value makeblastdb ignores (StringToUInt cannot convert it, or 0) is accepted with
    // -parse_seqids (review B N4; makeblastdb 2.17.0: pr2/fix/n4probe).
    for ignored in ["abc", "0", ""] {
        let accepted = clean_losat_command()
            .env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", ignored)
            .args([
                "makeblastdb",
                "-in",
                input_text,
                "-dbtype",
                "prot",
                "-out",
                base_text,
                "-parse_seqids",
                "-blastdb_version",
                "4",
            ])
            .output()
            .unwrap();
        assert!(accepted.status.success(), "{ignored:?}: {accepted:?}");
        remove(&V5);
    }
    // `-taxid_map -` reads the map from the standard input (CArg_InputFile maps `-` to cin;
    // review B S1).
    fs::write(
        &input,
        ">ref|NP_000001.1| first\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|second b\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let mut child = clean_losat_command()
        .args([
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
            "-parse_seqids",
            "-taxid_map",
            "-",
        ])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    {
        use std::io::Write as _;
        let mut stdin = child.stdin.take().unwrap();
        stdin
            .write_all(b"NP_000001.1 9606\nsecond 10090\n")
            .unwrap();
    }
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "{output:?}");
    let metadata = fs::read_to_string(format!("{base_text}.pjs")).unwrap();
    assert!(metadata.contains("\"number-of-taxids\": 2"), "{metadata}");
    remove(&V5);
    // `-taxid_map ""`: "File is not accessible" without a name part (s_ArgExptMsg; review B
    // N1), exit status 1. A directory is read as an empty map (a read error ends the reading;
    // review B S2).
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
        "-taxid_map",
        "",
    ]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Command line argument error: Argument \"taxid_map\". File is not accessible\n"
    );
    remove(&V5);
    #[cfg(target_os = "linux")]
    {
        let directory = std::env::temp_dir();
        let output = run_losat(&[
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
            "-parse_seqids",
            "-taxid_map",
            directory.to_str().unwrap(),
        ]);
        assert!(output.status.success(), "{output:?}");
        assert!(output.stderr.is_empty(), "{output:?}");
        let metadata = fs::read_to_string(format!("{base_text}.pjs")).unwrap();
        assert!(!metadata.contains("number-of-taxids"), "{metadata}");
        remove(&V5);
    }
    fs::remove_file(&input).unwrap();
    fs::remove_file(&map).unwrap();
}

const MT_PROTEIN: &str = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV";
const MT_NUCLEOTIDE: &str = "ATGAAAACCGCCTATATTGCCAAACAGCGTCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGTCAGCTGGAAGAACGTCTGGGCCTGATTGAAGTGCAGGCCCCGATTCTGAGCCGTGTGGGCGATGGCACCCAGGATAACCTGAGCGGCGCCGAAAAAGCCGTGCAGGTGAAAGTGAAAGCCCTGCCGGATGCCCAGTTTGAAGTGGTG";

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3176-3181
// ```c++
//     	arg_desc.SetConstraint(kArgMTMode,
//                                new CArgAllowValuesBetween(0, 2, true));
//     	arg_desc.SetDependency(kArgMTMode,
//                                CArgDescriptions::eRequires,
//                                kArgNumThreads);
// ```
// `-mt_mode` has no effect on LOSAT's output (CLI DECISIONS 5).
#[test]
fn mt_mode_does_not_change_the_output() {
    let protein = temp_path("mt_protein", "faa");
    let nucleotide = temp_path("mt_nucleotide", "fna");
    fs::write(&protein, format!(">p\n{MT_PROTEIN}\n")).unwrap();
    fs::write(&nucleotide, format!(">n\n{MT_NUCLEOTIDE}\n")).unwrap();
    let (protein, nucleotide) = (protein.to_str().unwrap(), nucleotide.to_str().unwrap());
    for (program, query, subject) in [
        ("blastn", nucleotide, nucleotide),
        ("blastp", protein, protein),
        ("blastx", nucleotide, protein),
        ("tblastn", protein, nucleotide),
    ] {
        let base = [
            program, "-query", query, "-subject", subject, "-outfmt", "6",
        ];
        let run = |words: &[&str]| {
            let mut command = clean_losat_command();
            // A build with the diagnostics feature needs this log for two or more threads.
            if cfg!(feature = "diagnostics") {
                command.env("LOSAT_BLASTX_WORKER_LOG", temp_path("mt_worker", "log"));
            }
            command.args(words).output().expect("run LOSAT")
        };
        let plain = run(&base);
        assert!(plain.status.success(), "{program}: {:?}", plain);
        assert!(!plain.stdout.is_empty(), "{program}");
        for mode in ["0", "1", "2"] {
            let mut words = base.to_vec();
            words.extend(["-mt_mode", mode, "-num_threads", "2"]);
            let output = run(&words);
            assert!(output.status.success(), "{program} {mode}: {:?}", output);
            assert_eq!(output.stdout, plain.stdout, "{program} -mt_mode {mode}");
        }
        // Out of range, and without `-num_threads`: argument errors.
        for extra in [
            &["-mt_mode", "3", "-num_threads", "2"][..],
            &["-mt_mode", "-1", "-num_threads", "2"],
            &["-mt_mode", "1"],
        ] {
            let mut words = base.to_vec();
            words.extend(extra);
            let output = run(&words);
            assert!(!output.status.success(), "{program} {extra:?}");
            assert!(output.stdout.is_empty(), "{program} {extra:?}");
        }
    }
    // NCBI's tblastx and psiblast have no -mt_mode.
    for program in ["tblastx", "psiblast"] {
        let output = run_losat(&[
            program,
            "-query",
            protein,
            "-subject",
            protein,
            "-mt_mode",
            "1",
            "-num_threads",
            "2",
        ]);
        assert!(!output.status.success(), "{program}");
    }
    fs::remove_file(protein).unwrap();
    fs::remove_file(nucleotide).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3354-3358
// ```c++
//     arg_desc.AddDefaultKey(
//             kArgUseIndex, "boolean",
//             "Use MegaBLAST database index",
//             CArgDescriptions::eBoolean, NStr::BoolToString(kDfltArgUseIndex));
// ```
#[test]
fn use_index_false_is_accepted_and_true_is_refused() {
    let nucleotide = temp_path("use_index", "fna");
    fs::write(&nucleotide, format!(">n\n{MT_NUCLEOTIDE}\n")).unwrap();
    let nucleotide = nucleotide.to_str().unwrap();
    let base = [
        "blastn", "-query", nucleotide, "-subject", nucleotide, "-outfmt", "6",
    ];
    let plain = run_losat(&base);
    assert!(plain.status.success() && !plain.stdout.is_empty());
    for value in ["false", "FALSE", "f", "no", "0"] {
        let mut words = base.to_vec();
        words.extend(["-use_index", value]);
        let output = run_losat(&words);
        assert!(output.status.success(), "-use_index {value}: {:?}", output);
        assert_eq!(output.stdout, plain.stdout, "-use_index {value}");
    }
    for value in ["true", "t", "1", "maybe"] {
        let mut words = base.to_vec();
        words.extend(["-use_index", value]);
        let output = run_losat(&words);
        assert!(!output.status.success(), "-use_index {value}");
        assert!(output.stdout.is_empty(), "-use_index {value}");
    }
    let mut words = base.to_vec();
    words.extend(["-use_index", "true"]);
    let error = String::from_utf8_lossy(&run_losat(&words).stderr).into_owned();
    assert!(
        error.contains("-use_index true (the megablast index) is not supported by LOSAT's BLASTN"),
        "{error}"
    );
    fs::remove_file(nucleotide).unwrap();
}
