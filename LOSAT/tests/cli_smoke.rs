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

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:275-285
// ```c
// if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
//     opt.SetLookupTableType(eCompressedAaLookupTable);
//     opt.SetWordThreshold(19.3);
//     if (args[kArgWordSize].AsInteger() > 5) {
//         opt.SetWordThreshold(21.0);
//     }
// }
// ```
#[test]
fn unsupported_blastp_option_fails_before_file_io() {
    let missing_query = temp_path("missing_query", "faa");
    let missing_subject = temp_path("missing_subject", "faa");
    let output = clean_losat_command()
        .args(["blastp", "-q"])
        .arg(&missing_query)
        .arg("-s")
        .arg(&missing_subject)
        .args(["--word-size", "7"])
        .output()
        .expect("run unsupported blastp");

    assert!(!output.status.success(), "unsupported blastp should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unsupported pure-Rust blastp word_size=7"),
        "unsupported error should name program, option, and reason: {stderr}"
    );
    assert!(
        !stderr.contains("failed to open query FASTA"),
        "unsupported option should fail before file I/O: {stderr}"
    );
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
// ```c
// CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
// CRef<CBlastInput> input(new CBlastInput(fasta));
// CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
// sequences = input->GetAllSeqs(*scope);
// ```
#[test]
fn missing_input_file_reports_explicit_query_error() {
    let missing_query = temp_path("missing_query", "fa");
    let missing_subject = temp_path("missing_subject", "fa");
    let output = clean_losat_command()
        .args(["blastn", "-q"])
        .arg(&missing_query)
        .arg("-s")
        .arg(&missing_subject)
        .output()
        .expect("run blastn with missing input");

    assert!(!output.status.success(), "missing input should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("failed to open query FASTA"),
        "missing input error should name query FASTA: {stderr}"
    );
    assert!(
        stderr.contains(&missing_query.display().to_string()),
        "missing input error should include path: {stderr}"
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
        .args(["blastn", "-q"])
        .arg(&query)
        .arg("-s")
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
        .args(["tblastx", "-q"])
        .arg(&query)
        .arg("-s")
        .arg(&subject)
        .args(["--seg=false", "--outfmt", "6", "-n", "1", "-o"])
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
