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
            String::from_utf8_lossy(&output.stdout).contains("Usage:"),
            "{args:?} help should include clap usage text"
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
        ("makedb", "makeblastdb"),
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
        "makedb",
        "makeblastdb",
    ] {
        for words in [vec![name, "-h"]] {
            let output = run_losat(&words);
            assert!(output.status.success(), "{words:?}: {:?}", output.status);
            assert!(output.stderr.is_empty(), "{words:?}");
            let text = String::from_utf8_lossy(&output.stdout);
            assert!(text.contains("Usage:"), "{words:?}: {text}");
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
// `losat makeblastdb` is `losat makedb`: the same files from the same input.
#[test]
fn makeblastdb_writes_the_files_makedb_writes() {
    let input = temp_path("mkdb_in", "faa");
    fs::write(
        &input,
        ">a first\nMKTAYIAKQRQISFVKSHFSRQ\n>b second\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let input_text = input.to_str().unwrap();
    let mut outputs = Vec::new();
    for name in ["makedb", "makeblastdb"] {
        let base = temp_path(name, "db");
        let base_text = base.to_str().unwrap();
        let output = run_losat(&[
            name, "-in", input_text, "-dbtype", "prot", "-out", base_text,
        ]);
        assert!(output.status.success(), "{name}: {:?}", output);
        // The log names the database's own path; compare it with the path removed.
        // The log also holds the time and the seconds of the run.
        let log = String::from_utf8_lossy(&output.stdout)
            .replace(base_text, "<db>")
            .lines()
            .filter(|line| !line.contains("current time") && !line.contains(" seconds"))
            .collect::<Vec<_>>()
            .join("\n");
        let mut files = Vec::new();
        for extension in ["phr", "psq", "pin"] {
            let path = format!("{base_text}.{extension}");
            files.push((extension, fs::read(&path).unwrap()));
            fs::remove_file(path).unwrap();
        }
        outputs.push((log.replace(input_text, "<in>"), files));
    }
    fs::remove_file(&input).unwrap();
    assert_eq!(outputs[0].0, outputs[1].0);
    for ((extension, first), (_, second)) in outputs[0].1.iter().zip(&outputs[1].1) {
        // The index file holds the creation date (README, section 5).
        if *extension == "pin" {
            assert_eq!(first.len(), second.len());
        } else {
            assert_eq!(first, second, ".{extension}");
        }
    }
    // `-help` works under both names, and the top level lists makeblastdb.
    for name in ["makedb", "makeblastdb"] {
        let output = run_losat(&[name, "-help"]);
        assert!(output.status.success());
        assert!(String::from_utf8_lossy(&output.stdout).contains("-dbtype <molecule_type>"));
    }
    let top = String::from_utf8_lossy(&run_losat(&["-help"]).stdout).into_owned();
    assert!(top.contains("makeblastdb"), "{top}");
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
