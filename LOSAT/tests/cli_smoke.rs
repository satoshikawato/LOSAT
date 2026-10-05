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
    if let Some(path) = option_env!("CARGO_BIN_EXE_LOSAT") {
        return PathBuf::from(path);
    }

    let mut path = std::env::current_exe().expect("current test executable");
    path.pop();
    if path.file_name().is_some_and(|name| name == "deps") {
        path.pop();
    }
    path.push(format!("LOSAT{}", std::env::consts::EXE_SUFFIX));
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

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/readers/fasta.cpp:391-396
// ```c
// NCBI_THROW2(CObjReaderParseException, eNoDefline,
//             "CFastaReader: Input doesn't start with"
//             " a defline or comment around line " + NStr::NumericToString(lineNum),
//              lineNum);
// ```
#[test]
fn malformed_fasta_reports_parse_error() {
    let query = temp_path("malformed_query", "fa");
    let subject = temp_path("valid_subject", "fa");
    fs::write(&query, "!not a FASTA defline\nACGTACGT\n").expect("write malformed FASTA");
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

    assert!(!output.status.success(), "malformed FASTA should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("failed to read query FASTA"),
        "malformed FASTA error should preserve query role: {stderr}"
    );
    assert!(
        stderr.contains("Expected > at record start."),
        "malformed FASTA error should preserve parser reason: {stderr}"
    );
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
