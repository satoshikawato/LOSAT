//! Engine-level tests for tblastx runs.

use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use LOSAT::algorithm::tblastx::{run, TblastxArgs};

fn temp_path(stem: &str, extension: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "losat_tblastx_{stem}_{}_{}.{}",
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

#[test]
fn db_gencode_controls_local_subject_search_translation() {
    let query = temp_path("query", "fa");
    let subject = temp_path("subject", "fa");
    let out = temp_path("out", "txt");

    fs::write(&query, ">q\nTGATGATGATGATGATGATGATGATGA\n").expect("write query FASTA");
    fs::write(&subject, ">s\nTGATGATGATGATGATGATGATGATGA\n").expect("write subject FASTA");

    let args = TblastxArgs {
        query: query.clone(),
        subject: subject.clone(),
        evalue: 10.0,
        threshold: 13,
        word_size: 3,
        num_threads: 1,
        out: Some(out.clone()),
        query_gencode: 4,
        db_gencode: 4,
        max_target_seqs: 500,
        seg: false,
        seg_window: 12,
        seg_locut: 2.2,
        seg_hicut: 2.5,
        window_size: 40,
        outfmt: "6".to_string(),
        culling_limit: 0,
    };

    run(args).expect("tblastx run");
    let output = fs::read_to_string(&out).expect("read tblastx output");

    let has_full_length_w_hit = output.lines().any(|line| {
        let fields: Vec<_> = line.split('\t').collect();
        fields.len() >= 12
            && fields[0] == "q"
            && fields[1] == "s"
            && fields[3] == "9"
            && fields[6] == "1"
            && fields[7] == "27"
            && fields[8] == "1"
            && fields[9] == "27"
    });

    let _ = fs::remove_file(query);
    let _ = fs::remove_file(subject);
    let _ = fs::remove_file(out);

    assert!(
        has_full_length_w_hit,
        "--db-gencode 4 should translate local subject TGA codons as W during search"
    );
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
// ```c
// #if _BLAST_DEBUG
// arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
//                  true);
// arg_desc.AddFlag("remote_verbose",
//                  "Produce verbose output for remote searches", true);
// #endif /* _BLAST_DEBUG */
// ```
#[test]
fn tblastx_output_filter_debug_stderr_is_env_gated() {
    let query = temp_path("debug_query", "fa");
    let subject = temp_path("debug_subject", "fa");
    let out_normal = temp_path("debug_normal", "txt");
    let out_debug = temp_path("debug_enabled", "txt");

    fs::write(&query, ">q\nATGATGATGATGATGATGATGATGATGATG\n").expect("write query FASTA");
    fs::write(&subject, ">s\nATGATGATGATGATGATGATGATGATGATG\n").expect("write subject FASTA");

    let base_args = [
        "tblastx",
        "-q",
        query.to_str().expect("query path UTF-8"),
        "-s",
        subject.to_str().expect("subject path UTF-8"),
        "--seg=false",
        "--outfmt",
        "6",
        "-n",
        "1",
    ];

    let normal = Command::new(losat_binary())
        .args(base_args)
        .args(["-o", out_normal.to_str().expect("normal output path UTF-8")])
        .env_remove("LOSAT_DEBUG_OUTPUT_FILTER")
        .env_remove("LOSAT_DEBUG_HSP_SAVING")
        .env_remove("LOSAT_DIAGNOSTICS")
        .env_remove("LOSAT_TIMING")
        .output()
        .expect("run normal tblastx");
    assert!(
        normal.status.success(),
        "normal tblastx run failed: {}",
        String::from_utf8_lossy(&normal.stderr)
    );
    let normal_stderr = String::from_utf8_lossy(&normal.stderr);
    assert!(
        !normal_stderr.contains("[DEBUG OUTPUT_FILTER]"),
        "normal run should not emit output-filter debug lines: {normal_stderr}"
    );
    assert!(
        !fs::read_to_string(&out_normal)
            .expect("read normal output")
            .is_empty(),
        "test input should produce at least one hit"
    );

    let debug = Command::new(losat_binary())
        .args(base_args)
        .args(["-o", out_debug.to_str().expect("debug output path UTF-8")])
        .env("LOSAT_DEBUG_OUTPUT_FILTER", "1")
        .env_remove("LOSAT_DEBUG_HSP_SAVING")
        .env_remove("LOSAT_DIAGNOSTICS")
        .env_remove("LOSAT_TIMING")
        .output()
        .expect("run debug tblastx");
    assert!(
        debug.status.success(),
        "debug tblastx run failed: {}",
        String::from_utf8_lossy(&debug.stderr)
    );
    let debug_stderr = String::from_utf8_lossy(&debug.stderr);
    assert!(
        debug_stderr.contains("[DEBUG OUTPUT_FILTER]"),
        "LOSAT_DEBUG_OUTPUT_FILTER=1 should emit output-filter debug lines: {debug_stderr}"
    );

    let _ = fs::remove_file(query);
    let _ = fs::remove_file(subject);
    let _ = fs::remove_file(out_normal);
    let _ = fs::remove_file(out_debug);
}
