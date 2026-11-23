//! Engine-level tests for tblastx runs.

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::algorithm::tblastx::{run, TblastxArgs};

fn temp_path(stem: &str, extension: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    crate::utils::test_dir::temp_dir().join(format!(
        "losat_tblastx_{stem}_{}_{}.{}",
        crate::utils::test_dir::process_id(),
        nanos,
        extension
    ))
}

/// The arguments of a local search of `query` against `subject` with the translation table 4 of
/// both (`-query_gencode 4 -db_gencode 4`), written to `out`.
fn local_args(query: &Path, subject: &Path, out: &Path) -> TblastxArgs {
    TblastxArgs {
        query: query.to_path_buf(),
        subject: Some(subject.to_path_buf()),
        db: None,
        db_filter: Default::default(),
        db_mask: Default::default(),
        formatting: Default::default(),
        dbsize: None,
        show_gis: false,
        html: false,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        evalue: 10.0,
        threshold: None,
        matrix: None,
        word_size: None,
        num_threads: 1,
        out: Some(out.to_path_buf()),
        query_gencode: 4,
        db_gencode: 4,
        max_intron_length: None,
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        report_counts: Default::default(),
        html_links: None,
        max_hsps: None,
        qcov_hsp_perc: None,
        seg: "no".to_string(),
        soft_masking: None,
        sum_stats: None,
        lcase_masking: false,
        strand: "both".to_string(),
        xdrop_ungap: None,
        window_size: None,
        outfmt: "6".to_string(),
        culling_limit: None,
        best_hit_overhang: None,
        best_hit_score_edge: None,
        subject_besthit: false,
        query_loc: None,
        subject_loc: None,
        archive: None,
        export_search_strategy: None,
        import_search_strategy: None,
        import: None,
        given: Default::default(),
        strategy: None,
        tabular_taxdb: None,
        subject_taxdb: None,
        taxdb_setting: None,
        first_query_ordinal: 0,
        one_query_batch: false,
        scope_refusal: None,
        matrix_file: None,
        command_line: None,
    }
}

#[test]
fn db_gencode_controls_local_subject_search_translation() {
    let query = temp_path("query", "fa");
    let subject = temp_path("subject", "fa");
    let out = temp_path("out", "txt");

    fs::write(&query, ">q\nTGATGATGATGATGATGATGATGATGA\n").expect("write query FASTA");
    fs::write(&subject, ">s\nTGATGATGATGATGATGATGATGATGA\n").expect("write subject FASTA");

    let args = local_args(&query, &subject, &out);

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:137 and
// c++/src/algo/blast/core/aa_ungapped.c:214-235 (a window of 0 or less takes the one-hit word
// finder): an imported search strategy can carry a negative `WindowSize`, which NCBI runs.
#[test]
fn a_negative_window_runs_the_one_hit_word_finder() {
    let query = temp_path("negwin_query", "fa");
    let subject = temp_path("negwin_subject", "fa");
    fs::write(&query, ">q\nTGATGATGATGATGATGATGATGATGA\n").expect("write query FASTA");
    fs::write(&subject, ">s\nTGATGATGATGATGATGATGATGATGA\n").expect("write subject FASTA");

    let search = |window_size: i32| {
        let out = temp_path("negwin_out", "txt");
        let mut args = local_args(&query, &subject, &out);
        args.window_size = Some(window_size);
        run(
            args,
            crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .expect("tblastx run");
        let output = fs::read_to_string(&out).expect("read tblastx output");
        let _ = fs::remove_file(out);
        output
    };
    let has_full_length_w_hit = |output: &str| {
        output.lines().any(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            fields.len() >= 12
                && fields[0] == "q"
                && fields[1] == "s"
                && fields[3] == "9"
                && fields[6] == "1"
                && fields[7] == "27"
                && fields[8] == "1"
                && fields[9] == "27"
        })
    };
    // -1: the word at 0 of a fresh table is skipped, the one at 1 on the same diagonal is not.
    assert!(has_full_length_w_hit(&search(-1)));
    // Any other value runs (a window of -2^31 is NCBI's undefined `Int4` sums, here exact): the
    // one-hit scan of the diagonals that a table of one cell holds all together.
    for window_size in [-40, -1_000_000, i32::MIN] {
        let _ = search(window_size);
    }

    let _ = fs::remove_file(query);
    let _ = fs::remove_file(subject);
}
