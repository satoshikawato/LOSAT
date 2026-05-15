//! Engine-level tests for tblastx runs.

use std::fs;
use std::path::PathBuf;
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
