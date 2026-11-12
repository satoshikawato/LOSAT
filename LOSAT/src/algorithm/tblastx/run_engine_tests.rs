//! Engine-level tests for tblastx runs.

use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::algorithm::tblastx::{run, TblastxArgs};

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
        subject: Some(subject.clone()),
        db: None,
        db_filter: Default::default(),
        dbsize: None,
        show_gis: false,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        evalue: 10.0,
        threshold: None,
        matrix: None,
        word_size: None,
        num_threads: 1,
        out: Some(out.clone()),
        query_gencode: 4,
        db_gencode: 4,
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        report_counts: Default::default(),
        max_hsps: None,
        qcov_hsp_perc: None,
        seg: "no".to_string(),
        soft_masking: None,
        lcase_masking: false,
        window_size: None,
        outfmt: "6".to_string(),
        culling_limit: None,
        query_loc: None,
        subject_loc: None,
        archive: None,
        tabular_taxdb: None,
        first_query_ordinal: 0,
        matrix_file: None,
    };

    run(
        args,
        crate::blastinput::ncbi_environment::ApplicationSettings::default(),
    )
    .expect("tblastx run");
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
        "-db_gencode 4 should translate local subject TGA codons as W during search"
    );
}
