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

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
// ```c
// status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
//                                dbseq_chunk_overlap);
// BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
//         hit_params->options, &hsp_list);
// Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
// status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list, ...);
// ```
#[test]
fn tblastx_parallel_chunks_matches_serial_on_real_chunk_boundaries() {
    if !cfg!(debug_assertions) {
        return;
    }

    let query =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/truncated/AP027133_50k.fasta");
    let subject = query.clone();
    let out_serial = temp_path("parallel_chunks_serial", "txt");
    let out_parallel = temp_path("parallel_chunks_enabled", "txt");
    let chunk_size_aa = "600";

    let base_args = [
        "tblastx",
        "-q",
        query.to_str().expect("query path UTF-8"),
        "-s",
        subject.to_str().expect("subject path UTF-8"),
        "--seg=false",
        "--outfmt",
        "6",
    ];

    let serial = Command::new(losat_binary())
        .args(base_args)
        .args([
            "-n",
            "1",
            "-o",
            out_serial.to_str().expect("serial output path UTF-8"),
        ])
        .env("LOSAT_TBLASTX_TEST_CHUNK_SIZE", chunk_size_aa)
        .env_remove("LOSAT_TBLASTX_PARALLEL_CHUNKS")
        .env_remove("LOSAT_DIAGNOSTICS")
        .env_remove("LOSAT_TIMING")
        .output()
        .expect("run serial chunked tblastx");
    assert!(
        serial.status.success(),
        "serial chunked tblastx failed: {}",
        String::from_utf8_lossy(&serial.stderr)
    );

    let parallel = Command::new(losat_binary())
        .args(base_args)
        .args([
            "-n",
            "2",
            "-o",
            out_parallel.to_str().expect("parallel output path UTF-8"),
        ])
        .env("LOSAT_TBLASTX_TEST_CHUNK_SIZE", chunk_size_aa)
        .env("LOSAT_TBLASTX_PARALLEL_CHUNKS", "1")
        .env_remove("LOSAT_DIAGNOSTICS")
        .env_remove("LOSAT_TIMING")
        .output()
        .expect("run parallel chunked tblastx");
    assert!(
        parallel.status.success(),
        "parallel chunked tblastx failed: {}",
        String::from_utf8_lossy(&parallel.stderr)
    );

    let serial_output = fs::read_to_string(&out_serial).expect("read serial output");
    let parallel_output = fs::read_to_string(&out_parallel).expect("read parallel output");
    assert!(
        !serial_output.is_empty(),
        "real AP027133 50k self-search should produce tblastx hits"
    );
    assert!(
        output_spans_forced_chunk_boundary(&serial_output, 600),
        "test output should include a subject hit spanning a forced chunk boundary"
    );
    assert_eq!(
        serial_output, parallel_output,
        "LOSAT_TBLASTX_PARALLEL_CHUNKS must preserve serial chunked output"
    );

    let _ = fs::remove_file(out_serial);
    let _ = fs::remove_file(out_parallel);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:83-123
// ```c
// for (s = s_first; s <= s_last; s++) {
//     ...
//     offset_pairs[i + totalhits].qs_offsets.s_off = s_off;
// }
// ```
//
// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:614
// ```c
// /* increment the offset in the diagonal array */
// Blast_ExtendWordExit(ewp, subject->length);
// ```
#[test]
fn tblastx_parallel_scan_chunks_matches_normal_on_ap027133_50k() {
    if !cfg!(debug_assertions) {
        return;
    }

    let query =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/truncated/AP027133_50k.fasta");
    let subject = query.clone();
    assert_parallel_scan_chunks_match_normal(
        &query,
        &subject,
        &["--query-gencode", "1", "--db-gencode", "1", "--seg=false"],
        true,
    );
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:452-584
// ```c
// aux_struct->WordFinder(..., init_hitlist, ...);
// BLAST_GetUngappedHSPList(init_hitlist, query_info, subject,
//         hit_params->options, &hsp_list);
// Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
// ```
#[test]
fn tblastx_parallel_scan_chunks_matches_normal_on_mjenmv() {
    if !cfg!(debug_assertions) {
        return;
    }

    let query = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/MjeNMV.fasta");
    let subject = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/MelaMJNV.fasta");
    assert_parallel_scan_chunks_match_normal(
        &query,
        &subject,
        &["--query-gencode", "1", "--db-gencode", "1", "--seg=true"],
        false,
    );
}

fn assert_parallel_scan_chunks_match_normal(
    query: &PathBuf,
    subject: &PathBuf,
    extra_args: &[&str],
    require_boundary_crossing_hsp: bool,
) {
    let out_normal = temp_path("parallel_scan_chunks_normal", "txt");
    let out_scan = temp_path("parallel_scan_chunks_enabled", "txt");
    let chunk_size_aa = "600";

    let mut base_args = vec![
        "tblastx",
        "-q",
        query.to_str().expect("query path UTF-8"),
        "-s",
        subject.to_str().expect("subject path UTF-8"),
        "--outfmt",
        "6",
    ];
    base_args.extend_from_slice(extra_args);

    let normal = Command::new(losat_binary())
        .args(&base_args)
        .args([
            "-n",
            "1",
            "-o",
            out_normal.to_str().expect("normal output path UTF-8"),
        ])
        .env_remove("LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS")
        .env_remove("LOSAT_TBLASTX_PARALLEL_CHUNKS")
        .env_remove("LOSAT_TBLASTX_TEST_SCAN_CHUNK_SIZE")
        .env_remove("LOSAT_TBLASTX_TEST_CHUNK_SIZE")
        .env_remove("LOSAT_DIAGNOSTICS")
        .env_remove("LOSAT_TIMING")
        .output()
        .expect("run normal tblastx");
    assert!(
        normal.status.success(),
        "normal tblastx failed: {}",
        String::from_utf8_lossy(&normal.stderr)
    );

    let scan = Command::new(losat_binary())
        .args(&base_args)
        .args([
            "-n",
            "2",
            "-o",
            out_scan.to_str().expect("scan output path UTF-8"),
        ])
        .env("LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS", "1")
        .env("LOSAT_TBLASTX_TEST_SCAN_CHUNK_SIZE", chunk_size_aa)
        .env_remove("LOSAT_TBLASTX_PARALLEL_CHUNKS")
        .env_remove("LOSAT_TBLASTX_TEST_CHUNK_SIZE")
        .env_remove("LOSAT_DIAGNOSTICS")
        .env_remove("LOSAT_TIMING")
        .output()
        .expect("run scan-chunk tblastx");
    assert!(
        scan.status.success(),
        "scan-chunk tblastx failed: {}",
        String::from_utf8_lossy(&scan.stderr)
    );

    let normal_output = fs::read_to_string(&out_normal).expect("read normal output");
    let scan_output = fs::read_to_string(&out_scan).expect("read scan output");
    assert!(
        !normal_output.is_empty(),
        "test input should produce tblastx hits"
    );
    if require_boundary_crossing_hsp {
        assert!(
            output_spans_forced_chunk_boundary(&normal_output, 600),
            "test output should include a subject hit spanning a forced scan boundary"
        );
    }
    assert_eq!(
        normal_output, scan_output,
        "LOSAT_TBLASTX_PARALLEL_SCAN_CHUNKS must preserve normal LOSAT output"
    );

    let _ = fs::remove_file(out_normal);
    let _ = fs::remove_file(out_scan);
}

fn output_spans_forced_chunk_boundary(output: &str, chunk_size_aa: usize) -> bool {
    let chunk_size_nt = chunk_size_aa * 3;
    output.lines().any(|line| {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 10 {
            return false;
        }
        let Ok(s_start) = fields[8].parse::<usize>() else {
            return false;
        };
        let Ok(s_end) = fields[9].parse::<usize>() else {
            return false;
        };
        let lo = s_start.min(s_end);
        let hi = s_start.max(s_end);
        let first_boundary = lo.div_ceil(chunk_size_nt) * chunk_size_nt;
        first_boundary > lo && first_boundary < hi
    })
}
