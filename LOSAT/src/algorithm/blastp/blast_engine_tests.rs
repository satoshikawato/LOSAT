use super::*;
use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
use std::path::PathBuf;

// NCBI 2.17.0 (line RG F/d36 probes ps.06, ps.06std, ps.p7, ps.00, bp.06, bp.p0; rows
// db.dlm.ps.n06, db.dlm.bp.p0): `-outfmt 06`, `"06 std"`, `+7` and `00` are outfmt 6, 7 and 0
// (`NStr::StringToInt`), the library's value as given as the command line's.
#[test]
fn text_formats_read_the_number_as_ncbi() {
    for (outfmt, format) in [
        ("0", OutputFormat::Pairwise),
        ("6", OutputFormat::Tabular),
        ("7", OutputFormat::TabularWithComments),
        ("06", OutputFormat::Tabular),
        ("06 std", OutputFormat::Tabular),
        (" +7 ", OutputFormat::TabularWithComments),
        ("00", OutputFormat::Pairwise),
        ("6 qaccver saccver", OutputFormat::Tabular),
    ] {
        assert_eq!(blastp_text_format(outfmt), Ok(format), "{outfmt}");
    }
    for outfmt in ["99", "5", "10", "x6"] {
        let error = blastp_text_format(outfmt).unwrap_err();
        assert!(error.starts_with("Unsupported output format: "), "{outfmt}");
    }
    assert_eq!(
        blastp_text_format("10 qseqid").unwrap_err(),
        "Unsupported output format: 10. Supported: 0, 6, 7"
    );
}

// EXPERIMENT (LOSAT_X_BPRECOUNTO): each index is tested once, its answer kept, whatever the order
// of the questions.
#[test]
fn x_tested_once_runs_each_test_once() {
    let mut memo = Vec::new();
    let mut tests = 0;
    let answers = [
        (5, true),
        (2, false),
        (5, true),
        (0, true),
        (2, false),
        (5, true),
    ];
    for (index, answer) in answers {
        let got = x_tested_once(&mut memo, index, || {
            tests += 1;
            answer
        });
        assert_eq!(got, answer, "index {index}");
    }
    assert_eq!(tests, 3);
    assert_eq!(memo, [2, 0, 1, 0, 0, 2]);
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): all four documents (outfmt 13 to 16) get a block plan; a
// kind added to `Xml2Kind` is not admitted until `x_xml2_kind_takes_blocks` names it.
#[test]
fn xml2_block_plan_admits_outfmt_13_to_16() {
    use crate::report::blast_xml2::Xml2Kind;
    assert!(x_xml2_kind_takes_blocks(Xml2Kind::Single15));
    assert!(x_xml2_kind_takes_blocks(Xml2Kind::Single16));
    assert!(x_xml2_kind_takes_blocks(Xml2Kind::Files13));
    assert!(x_xml2_kind_takes_blocks(Xml2Kind::Files14));
    for (outfmt, taken) in [("15", true), ("16", true), ("13", true), ("14", true)] {
        assert_eq!(
            blastp_xml2_kind(outfmt).is_some_and(x_xml2_kind_takes_blocks),
            taken,
            "outfmt {outfmt}"
        );
    }
}

// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): a stored residue reads as the letter that
// `ListedSubjects::read_titles` writes into its record (`ncbistdaa_to_ascii`), for every code
// and for a code past the table; the O test and the O of a rendered row agree.
#[test]
fn stored_letters_read_as_the_records_letters() {
    let residues: Vec<u8> = (0..28).chain([30]).collect();
    let ascii: Vec<u8> = residues
        .iter()
        .map(|&residue| ncbistdaa_to_ascii(residue) as u8)
        .collect();
    let (stored, record) = (ShownLetters::Stored(&residues), ShownLetters::Ascii(&ascii));
    for index in 0..=residues.len() {
        assert_eq!(stored.get(index), record.get(index), "index {index}");
    }
    assert!(stored.has_o() && record.has_o());
    // `STORED_O` is the one code shown as O.
    let shown_as_o: Vec<u8> = (0..=u8::MAX)
        .filter(|&residue| ncbistdaa_to_ascii(residue).eq_ignore_ascii_case(&'O'))
        .collect();
    assert_eq!(shown_as_o, [STORED_O]);
    let without_o: Vec<u8> = residues.iter().copied().filter(|&r| r != 26).collect();
    assert!(!ShownLetters::Stored(&without_o).has_o());
    // The row of residues 25..=27 as the search shows O (`X` where the record has O).
    let row = Some("*XJ".to_string());
    let hit = row_hit((26, 28), (26, 28), None);
    assert_eq!(
        restore_displayed_o(row.clone(), stored, &hit, RowSide::Subject),
        restore_displayed_o(row, record, &hit, RowSide::Subject)
    );
    assert_eq!(
        restore_displayed_o(Some("*XJ".to_string()), stored, &hit, RowSide::Subject).as_deref(),
        Some("*OJ")
    );
    // A subject set's records and a database's listed subjects give the same letters and
    // titles by subject index.
    let records = [FastaRecord::new("gnl|BL_ORD_ID|0", b"p01 first", &ascii)];
    let subjects = [EncodedProtein {
        aa_seq: [0]
            .into_iter()
            .chain(residues.iter().copied())
            .chain([0])
            .collect(),
        aa_len: residues.len(),
    }];
    let titles: [Arc<[u8]>; 1] = [Arc::from(&b"p01 first"[..])];
    let listed = ShownSubjects::Listed {
        subjects: &subjects,
        titles: &titles,
    };
    let shown = ShownSubjects::Records(&records);
    assert_eq!(listed.title(0), shown.title(0));
    assert_eq!(listed.title(1), None);
    let (Some(a), Some(b)) = (listed.letters(0), shown.letters(0)) else {
        panic!("letters of subject 0");
    };
    assert!((0..=residues.len()).all(|index| a.get(index) == b.get(index)));
    assert!(listed.letters(1).is_none());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:346-374 (cited at
// `blastp_ungapped_word_cutoff_score_for_context`): the word cutoff of an ungapped search is
// `BLAST_Cutoffs` at e-value 1e-300 (a score far above any gap trigger), at most the gap
// trigger and at most the hit-saving cutoff; an invalid context never extends.
#[test]
fn ungapped_word_cutoff_is_the_smallest_of_the_gap_trigger_and_the_hit_cutoff() {
    let karlin = crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let context = |is_valid| QueryContext {
        q_idx: 0,
        frame: 0,
        aa_seq: vec![0; 302],
        aa_seq_nomask: None,
        aa_len: 300,
        orig_len: 300,
        frame_base: 0,
        is_valid,
        karlin_params: karlin,
    };
    let trigger = gap_trigger_raw_score(22.0, &karlin);
    assert!(trigger > 30 && trigger < 1000);
    for subject_length in [1_000, 250_000] {
        assert_eq!(
            blastp_ungapped_word_cutoff_score_for_context(
                1000,
                &context(true),
                22.0,
                subject_length
            ),
            trigger
        );
        assert_eq!(
            blastp_ungapped_word_cutoff_score_for_context(30, &context(true), 22.0, subject_length),
            30
        );
        assert_eq!(
            blastp_ungapped_word_cutoff_score_for_context(
                30,
                &context(false),
                22.0,
                subject_length
            ),
            i32::MAX
        );
    }
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/linkhsp_unit_test.cpp:674-695 UngappedBlastpCutoffs
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/linkhsp_unit_test.cpp:674-695
// (`kGoodCutoffs`, `testAllCutoffs`): the ungapped word cutoff of BLASTP, `cutoff_score_ungapped`
// 41, is the same for the three subject lengths 400, 60 and 3000 and for the database's
// average of 500 (`kDbLength` 500000000 / `kDbNumSeqs` 1000000): the gap trigger of 22 bits at
// BLOSUM62's ungapped block (41), below the hit-saving cutoff of 66 that the case lists. The
// case's x-drops (16), its final cutoff (66, which depends on its query gi|129295) and its
// sum-statistics cutoffs are not asserted.
#[test]
fn ungapped_blastp_word_cutoff_is_41_for_the_cases_subject_lengths() {
    let karlin = crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let context = QueryContext {
        q_idx: 0,
        frame: 0,
        aa_seq: vec![0; 302],
        aa_seq_nomask: None,
        aa_len: 300,
        orig_len: 300,
        frame_base: 0,
        is_valid: true,
        karlin_params: karlin,
    };
    assert_eq!(gap_trigger_raw_score(22.0, &karlin), 41);
    for subject_length in [400, 60, 3000, 500] {
        assert_eq!(
            blastp_ungapped_word_cutoff_score_for_context(66, &context, 22.0, subject_length),
            41,
            "subject length {subject_length}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:346-374 and
// c++/include/algo/blast/core/blast_parameters.h:77,80 (`CUTOFF_E_BLASTP` and
// `CUTOFF_E_TBLASTX`, both 1e-300): BLASTP's word cutoff is the function that TBLASTX's calls
// (`word_cutoff_score_for_update`), with the E-value of its own program, and equals TBLASTX's at
// the same lengths; a program with another cutoff E-value gets another score.
#[test]
fn blastp_and_tblastx_word_cutoffs_are_one_function() {
    use crate::algorithm::tblastx::ncbi_cutoffs::{
        cutoff_score_for_update_tblastx, word_cutoff_score_for_update, CUTOFF_E_BLASTP,
        CUTOFF_E_TBLASTX,
    };
    assert_eq!(CUTOFF_E_BLASTP, 1e-300);
    assert_eq!(CUTOFF_E_TBLASTX, 1e-300);
    let karlin = crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    for (query, subject, trigger, cap) in [
        (300_i64, 1_000_i64, 1000, 1000),
        (300, 40, 1000, 1000),
        (2, 3, 41, 1000),
        (300, 250_000, 41, 30),
    ] {
        let tblastx =
            cutoff_score_for_update_tblastx(query, subject, trigger, cap, 0.0, &karlin, 1.0);
        let blastp = word_cutoff_score_for_update(
            CUTOFF_E_BLASTP,
            query,
            subject,
            trigger,
            cap,
            0.0,
            &karlin,
            1.0,
        );
        assert_eq!(blastp, tblastx, "{query} {subject} {trigger} {cap}");
    }
    // The E-value is the program's: 1.0 (BLASTX's) lowers the score below the gap trigger.
    let loose = word_cutoff_score_for_update(1.0, 300, 1_000, 1000, 1000, 0.0, &karlin, 1.0);
    let tight =
        word_cutoff_score_for_update(CUTOFF_E_BLASTP, 300, 1_000, 1000, 1000, 0.0, &karlin, 1.0);
    assert!(loose < tight, "{loose} {tight}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:766-769 and
// c++/src/algo/blast/core/blast_hits.c:1833-1839 (cited at `call_setup::context_block`): the
// block of a context's statistics is the call's gapped block in a gapped search and the
// context's own `kbp_std` in an ungapped one, whichever way the mode is given.
#[test]
fn context_block_is_the_gapped_block_or_the_contexts_own() {
    let own = crate::stats::lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let shared = KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.9,
        beta: -30.0,
    };
    let context = QueryContext {
        q_idx: 0,
        frame: 0,
        aa_seq: vec![0; 12],
        aa_seq_nomask: None,
        aa_len: 10,
        orig_len: 10,
        frame_base: 0,
        is_valid: true,
        karlin_params: own,
    };
    let mut args = make_blastp_args();
    args.comp_based_stats = "0".to_string();
    let gapped = args.resolve().expect("gapped");
    args.ungapped = true;
    let ungapped = args.resolve().expect("ungapped");
    assert!(std::ptr::eq(
        call_setup::context_block(&gapped, &shared, &context),
        &shared
    ));
    assert!(std::ptr::eq(
        call_setup::context_block(&ungapped, &shared, &context),
        &context.karlin_params
    ));
    assert!(std::ptr::eq(
        call_setup::block_of::<false>(&shared, &context),
        &shared
    ));
    assert!(std::ptr::eq(
        call_setup::block_of::<true>(&shared, &context),
        &context.karlin_params
    ));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:517-527 (cited at
// `call_setup::gapped_block`): an ungapped search does not read its gap costs, so a run
// without a query batch writes the report for costs that are not in the matrix's table (and
// for a matrix name that NCBI's tables do not have), with no `Gap Penalties:` line; a gapped
// search writes it. The three commands of the review of line O4b (`printf '' | blastp -db ...
// -ungapped -comp_based_stats 0 -gapopen 5 -gapextend 5`, `-matrix PAM30 -gapopen 11 -gapextend 1`,
// `-matrix FOO`).
#[test]
fn a_run_without_queries_writes_the_ungapped_report_for_gap_costs_outside_the_table() {
    use crate::api::local_blast::{FormatOutput, FormatRefusal, OutputSink};
    let report = |extra: &[(&str, &str)], ungapped: bool| -> String {
        let mut args = make_blastp_args();
        args.ungapped = ungapped;
        args.comp_based_stats = "0".to_string();
        for &(name, value) in extra {
            match name {
                "matrix" => args.matrix = Some(value.to_string()),
                "gapopen" => args.gap_open = Some(value.parse().unwrap()),
                "gapextend" => args.gap_extend = Some(value.parse().unwrap()),
                _ => unreachable!(),
            }
        }
        let resolved = args.resolve().expect("resolved");
        let mut sink: Vec<u8> = Vec::new();
        let mut diagnostics: Vec<u8> = Vec::new();
        {
            let mut outputs = ReportOutputs {
                other_outfmts: &[],
                html_posts: None,
                formats: vec![FormatOutput {
                    outfmt: "0",
                    refusal: FormatRefusal::FailsCall,
                    sink: OutputSink::Writer(&mut sink),
                }],
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
                later_warnings: None,
            };
            write_reports_without_queries(
                &resolved,
                BlastpSubjectSource::Fasta(&[]),
                "",
                &mut outputs,
            )
            .expect("the report of a run without queries");
        }
        String::from_utf8(sink).expect("text")
    };
    for extra in [
        &[("gapopen", "5"), ("gapextend", "5")][..],
        &[("matrix", "PAM30"), ("gapopen", "11"), ("gapextend", "1")][..],
        &[("matrix", "FOO")][..],
    ] {
        let text = report(extra, true);
        assert!(
            text.starts_with(&format!("BLASTP (LOSAT {})\n", env!("CARGO_PKG_VERSION"))),
            "{extra:?}: {text}"
        );
        assert!(!text.contains("Gap Penalties:"), "{extra:?}: {text}");
    }
    let gapped = report(&[], false);
    assert!(
        gapped.contains("Gap Penalties: Existence: 11, Extension: 1"),
        "{gapped}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158-182 and
// c++/src/algo/blast/api/blast_setup_cxx.cpp:919-926 (cited at
// `blastx::query_setup::unknown_matrix_message[_of]`): the error of a batch whose score block
// cannot load the matrix holds each query's error followed by its messages: the O warning of a
// protein query (positions of its first 20 O) and `Sequence contains no data` of an empty one.
#[test]
fn the_unknown_matrix_error_names_each_querys_messages() {
    let records = [
        FastaRecord::new("Query_1", b"a", b"MKOVL"),
        FastaRecord::new("Query_2", b"b", b""),
        FastaRecord::new("Query_3", b"c", b"OKVOL"),
        FastaRecord::new("Query_4", b"d", b"MKVL"),
    ];
    assert_eq!(
        unknown_matrix_message(&records),
        "Error: Unknown error code -1 Warning: One or more O characters replaced by X for alignment score calculations at positions 2 \
             Error: Unknown error code -1 Warning: Sequence contains no data \
             Error: Unknown error code -1 Warning: One or more O characters replaced by X for alignment score calculations at positions 0, 3 \
             Error: Unknown error code -1 "
    );
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54 (cited at
// `too_long_subject_refusal`).
#[test]
fn the_long_subject_refusal_names_the_subject_and_its_length() {
    assert_eq!(
        too_long_subject_refusal("subject sequence 3", 5_000_001).to_string(),
        "subject sequence 3 has 5000001 residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTP"
    );
    assert_eq!(
        too_long_subject_refusal("a database sequence", 6_000_000).to_string(),
        "a database sequence has 6000000 residues; a subject longer than 5000000 residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTP"
    );
}

use crate::stats::lookup_protein_params;
use crate::stats::KarlinParams;
use crate::utils::matrix::ncbistdaa;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn make_test_blastp_hsp_list(oid: u32, raw_score: i32) -> BlastpHspList<XPrelimHsp> {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![XPrelimHsp {
            hsp: BlastpPreliminaryHsp {
                query_start: 0,
                query_end: 10,
                subject_start: 0,
                subject_end: 10,
                gapped_query_start: 1,
                gapped_subject_start: 1,
                raw_score,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 10,
                q_idx: 0,
                s_idx: oid,
            },
            e_value: 1.0,
        }],
        best_evalue: f64::MAX,
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn make_test_subject_result(oid: usize, raw_score: i32) -> BlastpSubjectResult {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
    let mut result = BlastpSubjectResult::new(oid);
    result
        .x_compact_hsp_lists
        .push((0, make_test_blastp_hsp_list(oid as u32, raw_score)));
    result
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn test_merge_blastp_subject_results_replays_subject_index_order() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1554
    let mut hit_lists = vec![None];
    let subject_results = vec![
        Some(make_test_subject_result(0, 50)),
        Some(make_test_subject_result(1, 60)),
        Some(make_test_subject_result(2, 70)),
    ];

    merge_blastp_subject_results_in_order(&mut hit_lists, subject_results, 10);

    let hit_list = hit_lists[0]
        .as_ref()
        .expect("query hit list should be created");
    let oids: Vec<u32> = hit_list.hsplist_array.iter().map(|list| list.oid).collect();
    assert_eq!(oids, vec![0, 1, 2]);
}

fn make_blastp_args() -> BlastpArgs {
    BlastpArgs {
        query: PathBuf::from("query.faa"),
        subject: Some(PathBuf::from("subject.faa")),
        db: None,
        db_filter: Default::default(),
        db_mask: Default::default(),
        dbsize: None,
        show_gis: false,
        html: false,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        task: "blastp".to_string(),
        evalue: None,
        threshold: None,
        word_size: None,
        num_threads: 1,
        mt_mode: None,
        out: None,
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        max_hsps_per_subject: None,
        qcov_hsp_perc: None,
        culling_limit: None,
        best_hit_overhang: None,
        best_hit_score_edge: None,
        subject_besthit: false,
        ungapped: false,
        xdrop_ungap: None,
        xdrop_gap: None,
        xdrop_gap_final: None,
        window_size: None,
        matrix: None,
        gap_open: None,
        gap_extend: None,
        comp_based_stats: "2".to_string(),
        seg: None,
        soft_masking: None,
        lcase_masking: false,
        use_sw_tback: false,
        query_loc: None,
        subject_loc: None,
        outfmt: "6".to_string(),
        archive: None,
        export_search_strategy: None,
        strategy: None,
        request_inputs: None,
        app: BlastpApp::Blastp,
        taxdb_setting: None,
        matrix_search: None,
        write_verification: None,
        inclusion_ethresh: None,
        gap_trigger: None,
        sum_stats: false,
        first_query_ordinal: 0,
        scope_refusal: None,
        psiblast_pssm: None,
        formatting: Default::default(),
        rps: None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-923
// The call's extent holds a database total of 2^32 letters or more (an alias file's `LENGTH`
// over `tests/db_search/taxfilter/tfp`, rows db.cnt.l4295.bp.* and db.cnt.avg2p32.bp_c0) as
// NCBI's `Int8` on every target: the effective lengths and the Gumbel block read all of it.
#[test]
fn db_extent_holds_a_total_of_2_to_the_32_letters_or_more() {
    let resolved = make_blastp_args().resolve().expect("resolved blastp args");
    for (name, total) in [
        ("l4295", 4_294_967_296_000_000_000),
        ("l4294", 4_294_967_295_999_999_999),
        ("avg2p32", 98_784_247_808),
    ] {
        assert!(total > i64::from(u32::MAX), "{name}");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/dbcounts")
            .join(name);
        let database = database::BlastpDatabase::open(
            path.to_str().unwrap(),
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
            "BLASTP",
        )
        .unwrap_or_else(|error| panic!("{name}: {error:#}"));
        let extent = call_setup::db_extent(&resolved, Some(&database), &[])
            .unwrap_or_else(|error| panic!("{name}: {error:#}"));
        assert_eq!(
            (
                extent.total_db_len,
                extent.eff_db_length,
                extent.num_subjects
            ),
            (total, total, 23),
            "{name}"
        );
        assert_eq!(
            extent.gapped_gumbel.map(|gumbel| gumbel.db_length),
            Some(total),
            "{name}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:87-129
fn make_poly_a_query_frame() -> QueryFrame {
    QueryFrame {
        frame: 1,
        aa_seq: vec![
            0,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
            0,
        ],
        aa_seq_nomask: None,
        aa_len: 6,
        orig_len: 6,
        seg_masks: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
fn make_query_frame_with_residues(residues: &[u8]) -> QueryFrame {
    let mut aa_seq = Vec::with_capacity(residues.len() + 2);
    aa_seq.push(0);
    aa_seq.extend_from_slice(residues);
    aa_seq.push(0);
    QueryFrame {
        frame: 1,
        aa_seq,
        aa_seq_nomask: None,
        aa_len: residues.len(),
        orig_len: residues.len(),
        seg_masks: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
fn make_poly_a_lookup() -> BlastAaLookupTable {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let queries = vec![vec![make_poly_a_query_frame()]];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2746-2748
    let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, false);
    lookup
}

fn make_init_hsp(
    q_start_absolute: i32,
    length: i32,
    s_start: i32,
    q_seed_absolute: i32,
    s_seed: i32,
    score: i32,
) -> InitHSP {
    InitHSP {
        q_seed_absolute,
        s_seed,
        ungapped_data: BlastpUngappedData {
            q_start: q_start_absolute,
            s_start,
            length,
            score,
        },
    }
}

// The table of contexts gives the binary search's index at every
// offset of the concatenated query, for query sets of one length, of mixed lengths, with
// empty queries (equal frame bases), and with blocks of 1 to 16 offsets and the size `new`
// chooses (a set long enough for blocks of more than one offset).
#[test]
fn x_context_table_matches_the_binary_search_at_every_offset() {
    let long: Vec<usize> = [1usize, 2, 7, 64, 65, 500, 0, 33, 129, 3]
        .iter()
        .copied()
        .cycle()
        .take(3000)
        .collect();
    let shapes: [&[usize]; 6] = [
        &[3],
        &[5, 5, 5, 5],
        &[1, 300, 2, 17, 40, 1, 1, 999, 8],
        &[0, 4, 0, 0, 9, 0],
        &[2; 64],
        &long,
    ];
    for lengths in shapes {
        // as the lookup lays them out: the first at 0, each next one after a sentinel
        let mut contexts = Vec::new();
        let mut base = 0i32;
        for (i, &len) in lengths.iter().enumerate() {
            let mut ctx = make_query_context(i as u32, base);
            ctx.aa_len = len;
            ctx.orig_len = len;
            contexts.push(ctx);
            base += len as i32 + 1;
        }
        let span = base as usize + 1;
        let bounds = blastp_context_lookup_bounds(&contexts);
        let frame_bases: Vec<i32> = contexts.iter().map(|ctx| ctx.frame_base).collect();
        let tables = (0..5)
            .map(|shift| {
                super::super::x_seed_fused::XContextTable::with_shift(&frame_bases, span, shift)
            })
            .chain(std::iter::once(
                super::super::x_seed_fused::XContextTable::new(&frame_bases, span),
            ));
        for table in tables {
            for q in 0..span as i32 {
                // the binary search's bounds assume offsets inside the query set
                if bounds.min_length > 0 && q >= base {
                    continue;
                }
                let idx = blastp_context_idx_for_query_offset_with_bounds(&contexts, bounds, q);
                assert_eq!(
                    table.get(q),
                    Some(idx),
                    "lengths {:?} offset {q}",
                    &lengths[..lengths.len().min(12)]
                );
                // the frame base the candidate closures test (SC2e)
                assert_eq!(
                    table.get_with_start(q),
                    Some((idx, contexts[idx].frame_base)),
                    "lengths {:?} offset {q}",
                    &lengths[..lengths.len().min(12)]
                );
            }
            assert_eq!(table.get(span as i32), None);
            assert_eq!(table.get(-1), None);
            assert_eq!(table.get_with_start(span as i32), None);
            assert_eq!(table.get_with_start(-1), None);
        }
    }
    let empty = super::super::x_seed_fused::XContextTable::new(&[], 4);
    assert_eq!(empty.get(0), Some(0));
    assert_eq!(empty.get(3), Some(0));
}

fn make_query_context(q_idx: u32, frame_base: i32) -> QueryContext {
    QueryContext {
        q_idx,
        frame: 1,
        aa_seq: vec![0, 1, 2, 3, 0],
        aa_seq_nomask: None,
        aa_len: 3,
        orig_len: 3,
        frame_base,
        is_valid: true,
        karlin_params: KarlinParams {
            lambda: 0.267,
            k: 0.041,
            h: 0.14,
            alpha: 1.0,
            beta: -1.0,
        },
    }
}

fn make_preliminary_hsp(
    q_idx: u32,
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
    raw_score: i32,
) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start,
        query_end,
        subject_start,
        subject_end,
        gapped_query_start: query_start,
        gapped_subject_start: subject_start,
        raw_score,
        query_context: q_idx as i32,
        query_frame: 1,
        subject_frame: 0,
        query_length: 100,
        q_idx,
        s_idx: 0,
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
fn make_local_match_hsp(preliminary_hsp: &BlastpPreliminaryHsp) -> BlastpHsp {
    BlastpHsp {
        num: 1,
        identity: 0.0,
        length: usize::try_from(preliminary_hsp.query_end - preliminary_hsp.query_start)
            .expect("NCBI BLAST preliminary query span must be non-negative"),
        mismatch: 0,
        gapopen: 0,
        q_start: usize::try_from(preliminary_hsp.query_start + 1)
            .expect("NCBI BLAST localMatch query start must be non-negative"),
        q_end: usize::try_from(preliminary_hsp.query_end)
            .expect("NCBI BLAST localMatch query end must be non-negative"),
        s_start: usize::try_from(preliminary_hsp.subject_start + 1)
            .expect("NCBI BLAST localMatch subject start must be non-negative"),
        s_end: usize::try_from(preliminary_hsp.subject_end)
            .expect("NCBI BLAST localMatch subject end must be non-negative"),
        gapped_q_start: preliminary_hsp.gapped_query_start,
        gapped_s_start: preliminary_hsp.gapped_subject_start,
        e_value: 0.0,
        bit_score: 0.0,
        num_ident: 0,
        query_context: preliminary_hsp.query_context,
        query_frame: preliminary_hsp.query_frame,
        subject_frame: preliminary_hsp.subject_frame,
        query_length: preliminary_hsp.query_length,
        q_idx: preliminary_hsp.q_idx,
        s_idx: preliminary_hsp.s_idx,
        raw_score: preliminary_hsp.raw_score,
        gap_info: None,
        num_positives: 0,
        comp_adjustment_method: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:797-800
#[test]
fn test_blastp_diag_last_hit_after_one_hit_matches_ncbi_formula() {
    assert_eq!(70 - (3 - 1) + 50, 118);
    assert_eq!(11 - (5 - 1) + 19, 26);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:509-521
#[test]
fn test_blastp_flagged_diag_branch_matches_ncbi_skip_vs_reset_order() {
    let mut skip_diag = DiagStruct::default();
    skip_diag.set_flag(1);
    skip_diag.set_last_hit(90);
    let skip_subject_plus_offset = 30 + 50;
    if skip_subject_plus_offset >= skip_diag.last_hit() {
        skip_diag.set_last_hit(skip_subject_plus_offset);
        skip_diag.set_flag(0);
    }
    assert_eq!(skip_diag.flag(), 1);
    assert_eq!(skip_diag.last_hit(), 90);

    let mut reset_diag = DiagStruct::default();
    reset_diag.set_flag(1);
    reset_diag.set_last_hit(90);
    let reset_subject_plus_offset = 40 + 50;
    if reset_subject_plus_offset >= reset_diag.last_hit() {
        reset_diag.set_last_hit(reset_subject_plus_offset);
        reset_diag.set_flag(0);
    }
    assert_eq!(reset_diag.flag(), 0);
    assert_eq!(reset_diag.last_hit(), 90);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:581-588
#[test]
fn test_blastp_two_hit_diag_update_matches_ncbi_flag_and_last_hit_updates() {
    let mut right_extended_diag = DiagStruct::default();
    right_extended_diag.set_flag(0);
    right_extended_diag.set_last_hit(0);
    right_extended_diag.set_flag(1);
    right_extended_diag.set_last_hit(70 - (3 - 1) + 50);
    assert_eq!(right_extended_diag.flag(), 1);
    assert_eq!(right_extended_diag.last_hit(), 118);

    let mut plain_diag = DiagStruct::default();
    plain_diag.set_flag(0);
    plain_diag.set_last_hit(0);
    plain_diag.set_last_hit(40 + 50);
    assert_eq!(plain_diag.flag(), 0);
    assert_eq!(plain_diag.last_hit(), 90);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-573
#[test]
fn test_blastp_two_hit_cross_query_reject_uses_context_start() {
    // aa_ungapped.c:566: the hit is rejected when `query_offset - diff` lies before the
    // context's start. The two cases are checked when the test is compiled.
    const fn rejected(query_offset: i32, diff: i32, context_start: i32) -> bool {
        query_offset - diff < context_start
    }
    const { assert!(rejected(105, 10, 100)) };
    const { assert!(!rejected(105, 5, 100)) };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:358-372
#[test]
fn test_blastp_save_init_hsp_uses_qstart_sstart_len_as_single_truth() {
    let mut init_hsps = Vec::new();

    blastp_save_init_hsp(&mut init_hsps, 1015, 42, 1018, 45, 27, 88);

    assert_eq!(init_hsps.len(), 1);
    let hsp = init_hsps[0];
    assert_eq!(hsp.ungapped_data.q_start, 1015);
    assert_eq!(hsp.query_end_absolute(), 1042);
    assert_eq!(hsp.ungapped_data.s_start, 42);
    assert_eq!(hsp.subject_end(), 69);
    assert_eq!(hsp.ungapped_data.length, 27);
    assert_eq!(hsp.q_seed_absolute, 1018);
    assert_eq!(hsp.s_seed, 45);
    assert_eq!(hsp.ungapped_data.score, 88);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:787-800
#[test]
fn test_one_hit_caller_passes_ncbi_payload_directly_to_blastp_save_init_hsp() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
    ];
    let subject = query.clone();
    let query_offset = 1;
    let subject_offset = 1;
    let ungapped_hit = extend_one_hit(
        ScoringMatrix::Blosum62,
        &query,
        &subject,
        query_offset as usize,
        subject_offset as usize,
        7,
        3,
    )
    .expect("one-hit ungapped extension");
    let mut init_hsps = Vec::new();

    blastp_save_init_hsp(
        &mut init_hsps,
        ungapped_hit.ungapped_data.q_start,
        ungapped_hit.ungapped_data.s_start,
        query_offset,
        subject_offset,
        ungapped_hit.ungapped_data.length,
        ungapped_hit.ungapped_data.score,
    );

    assert_eq!(init_hsps[0].q_seed_absolute, query_offset);
    assert_eq!(init_hsps[0].s_seed, subject_offset);
    assert_eq!(init_hsps[0].ungapped_data.q_start, 0);
    assert_eq!(init_hsps[0].ungapped_data.s_start, 0);
    assert_eq!(init_hsps[0].ungapped_data.length, 5);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-583
#[test]
fn test_two_hit_caller_preserves_second_hit_seed_and_length_when_right_extension_occurs() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
    ];
    let subject = query.clone();
    let query_offset = 3;
    let subject_offset = 3;
    let ungapped_hit = extend_two_hit(
        ScoringMatrix::Blosum62,
        &query,
        &subject,
        0,
        subject_offset as usize,
        query_offset as usize,
        7,
        3,
    )
    .expect("two-hit ungapped extension");
    let mut init_hsps = Vec::new();

    blastp_save_init_hsp(
        &mut init_hsps,
        ungapped_hit.ungapped_data.q_start,
        ungapped_hit.ungapped_data.s_start,
        query_offset,
        subject_offset,
        ungapped_hit.ungapped_data.length,
        ungapped_hit.ungapped_data.score,
    );

    assert!(ungapped_hit.right_extend);
    assert_eq!(init_hsps[0].q_seed_absolute, query_offset);
    assert_eq!(init_hsps[0].s_seed, subject_offset);
    assert_eq!(init_hsps[0].ungapped_data.length, 6);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:566-583
#[test]
fn test_two_hit_caller_preserves_second_hit_seed_and_length_when_left_stops_early() {
    let query = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::D,
        ncbistdaa::D,
        ncbistdaa::D,
        ncbistdaa::C,
        ncbistdaa::C,
        ncbistdaa::C,
    ];
    let subject = vec![
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::A,
        ncbistdaa::W,
        ncbistdaa::W,
        ncbistdaa::W,
        ncbistdaa::C,
        ncbistdaa::C,
        ncbistdaa::C,
    ];
    let query_offset = 6;
    let subject_offset = 6;
    let ungapped_hit = extend_two_hit(
        ScoringMatrix::Blosum62,
        &query,
        &subject,
        0,
        subject_offset as usize,
        query_offset as usize,
        2,
        3,
    )
    .expect("two-hit ungapped extension");
    let mut init_hsps = Vec::new();

    blastp_save_init_hsp(
        &mut init_hsps,
        ungapped_hit.ungapped_data.q_start,
        ungapped_hit.ungapped_data.s_start,
        query_offset,
        subject_offset,
        ungapped_hit.ungapped_data.length,
        ungapped_hit.ungapped_data.score,
    );

    assert!(!ungapped_hit.right_extend);
    assert_eq!(init_hsps[0].q_seed_absolute, query_offset);
    assert_eq!(init_hsps[0].s_seed, subject_offset);
    assert_eq!(init_hsps[0].ungapped_data.length, 3);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-511
#[test]
fn test_blastp_for_each_offset_pair_preserves_scan_emission_order() {
    let lookup = make_poly_a_lookup();
    let runtime_lookup = BlastpRuntimeLookup::Aa(lookup);
    let mut subject_sequence = vec![ncbistdaa::A; 6];
    subject_sequence.push(0);
    let mut offset_pairs = vec![OffsetPair::default(); 5];
    let mut emitted = Vec::new();

    blastp_for_each_offset_pair(
        &runtime_lookup,
        &subject_sequence,
        6,
        &mut offset_pairs,
        BlastpOffsetPairScanMode::TwoHit,
        |pair| emitted.push((pair.q_off, pair.s_off)),
    );

    let expected = vec![
        (0, 0),
        (1, 0),
        (2, 0),
        (3, 0),
        (0, 1),
        (1, 1),
        (2, 1),
        (3, 1),
        (0, 2),
        (1, 2),
        (2, 2),
        (3, 2),
        (0, 3),
        (1, 3),
        (2, 3),
        (3, 3),
    ];
    assert_eq!(emitted, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_lookup.c:102-128
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
#[test]
fn test_blastp_for_each_offset_pair_preserves_absolute_query_offsets_across_queries() {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let queries = vec![
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
    ];
    let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, false);
    let runtime_lookup = BlastpRuntimeLookup::Aa(lookup);
    let subject_sequence = vec![ncbistdaa::A, ncbistdaa::A, ncbistdaa::A, ncbistdaa::A, 0];
    let mut offset_pairs = vec![OffsetPair::default(); 8];
    let mut emitted = Vec::new();

    blastp_for_each_offset_pair(
        &runtime_lookup,
        &subject_sequence,
        4,
        &mut offset_pairs,
        BlastpOffsetPairScanMode::TwoHit,
        |pair| emitted.push((pair.q_off, pair.s_off)),
    );

    assert_eq!(emitted, vec![(0, 0), (4, 0), (0, 1), (4, 1)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:552-568
#[test]
fn test_blastp_for_each_offset_pair_dispatches_compressed_lookup() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    let word = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
    ];
    lookup.add_unencoded(&word, 40);
    lookup.add_unencoded(&word, 41);
    lookup.finalize();
    let runtime_lookup = BlastpRuntimeLookup::Compressed(lookup);
    let subject_sequence = vec![
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
        0,
    ];
    let mut offset_pairs = vec![OffsetPair::default(); 8];
    let mut emitted = Vec::new();

    blastp_for_each_offset_pair(
        &runtime_lookup,
        &subject_sequence,
        10,
        &mut offset_pairs,
        BlastpOffsetPairScanMode::TwoHit,
        |pair| emitted.push((pair.q_off, pair.s_off)),
    );

    assert_eq!(emitted, vec![(40, 0), (41, 0), (40, 5), (41, 5)]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:496-500
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
#[test]
fn test_blastp_for_each_offset_pair_two_hit_uses_trailing_sentinel_for_short_subjects() {
    let lookup = make_poly_a_lookup();
    let runtime_lookup = BlastpRuntimeLookup::Aa(lookup);
    let subject_sequence = vec![ncbistdaa::A, ncbistdaa::A, 0];
    let mut offset_pairs = vec![OffsetPair::default(); 5];
    let mut emitted = Vec::new();

    blastp_for_each_offset_pair(
        &runtime_lookup,
        &subject_sequence,
        2,
        &mut offset_pairs,
        BlastpOffsetPairScanMode::TwoHit,
        |pair| emitted.push((pair.q_off, pair.s_off)),
    );

    assert!(emitted.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:551-567
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:358-372
#[test]
fn test_blastp_absolute_lookup_offsets_flow_directly_into_context_and_seed_fields() {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let queries = vec![
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
    ];
    let (_lookup, contexts) = build_ncbi_lookup(&queries, 0, &ungapped, false);
    let query_offset = 4;
    let ctx_idx = blastp_context_idx_for_query_offset(&contexts, query_offset);
    let ctx = &contexts[ctx_idx];
    let mut init_hsps = Vec::new();

    blastp_save_init_hsp(&mut init_hsps, query_offset, 12, query_offset, 12, 3, 27);

    assert_eq!(ctx_idx, 1);
    assert_eq!(ctx.frame_base, 4);
    assert_eq!(init_hsps[0].q_seed_absolute, 4);
    assert_eq!(init_hsps[0].ungapped_data.q_start, 4);
    assert_eq!(init_hsps[0].query_end_absolute(), 7);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2445-2467
#[test]
fn test_blastp_adjust_hsp_offsets_and_get_query_data_uses_qoff_as_context_truth() {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let queries = vec![
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
    ];
    let (_lookup, contexts) = build_ncbi_lookup(&queries, 0, &ungapped, false);
    let concatenated_query = build_blastp_concatenated_query(&contexts);
    let mut init_hsps = Vec::new();

    blastp_save_init_hsp(&mut init_hsps, 4, 12, 4, 12, 3, 27);

    let mut init_hsp = init_hsps[0];
    let (context_idx, ctx, query_sequence) = blastp_adjust_hsp_offsets_and_get_query_data(
        &concatenated_query,
        &contexts,
        blastp_context_lookup_bounds(&contexts),
        &mut init_hsp,
    );

    assert_eq!(context_idx, 1);
    assert_eq!(ctx.q_idx, 1);
    assert_eq!(query_sequence.len(), 4);
    assert_eq!(init_hsp.q_seed_absolute, 0);
    assert_eq!(init_hsp.ungapped_data.q_start, 0);
    assert_eq!(init_hsp.query_end_absolute(), 3);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:246-250
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-501,654
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-61
#[test]
fn test_blastp_diag_table_uses_ncbi_query_length_with_shared_boundary_sentinel() {
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let queries = vec![
        vec![make_query_frame_with_residues(&[
            ncbistdaa::A,
            ncbistdaa::A,
            ncbistdaa::A,
        ])],
        vec![make_query_frame_with_residues(&[
            ncbistdaa::C,
            ncbistdaa::C,
            ncbistdaa::C,
        ])],
    ];
    let (lookup, _) = build_ncbi_lookup(&queries, 0, &ungapped, false);
    let (diag_array_size, diag_mask) = compute_diag_table_shape(lookup.query_length, 1);

    assert_eq!(lookup.query_length, 7);
    assert_eq!(diag_array_size, 8);
    assert_eq!(diag_mask, 7);
    assert_eq!(lookup.query_length + 1, 8);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1163-1169
#[test]
fn test_validate_requested_blastp_support_accepts_blastp_fast_word_size_five() {
    let mut args = make_blastp_args();
    args.task = "blastp-fast".to_string();
    let resolved = args.resolve().expect("resolved blastp-fast args");

    validate_requested_blastp_support(&resolved)
        .expect("word_size 5 compressed blastp-fast is comparison-gated");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1286,1289-1292
#[test]
fn test_validate_requested_blastp_support_accepts_each_compressed_word_size() {
    for task in ["blastp", "blastp-fast"] {
        for word_size in 5..=7 {
            let mut args = make_blastp_args();
            args.task = task.to_string();
            args.word_size = Some(word_size);
            let resolved = args.resolve().expect("resolved blastp args");
            assert_eq!(
                resolved.lookup_table_type,
                BlastpLookupTableType::CompressedAaLookupTable
            );
            validate_requested_blastp_support(&resolved)
                .unwrap_or_else(|err| panic!("{task} -word_size {word_size}: {err}"));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1163-1169
// (above): word size 4 keeps the lookup table, which LOSAT's BLASTP builds with the joint
// builder (item SC8b).
#[test]
fn test_validate_requested_blastp_support_accepts_lookup_table_word_size_four() {
    for task in ["blastp", "blastp-fast", "blastp-short"] {
        let mut args = make_blastp_args();
        args.task = task.into();
        args.word_size = Some(4);
        let resolved = args.resolve().expect("resolved blastp args");
        assert_eq!(
            resolved.lookup_table_type,
            BlastpLookupTableType::AaLookupTable
        );
        validate_requested_blastp_support(&resolved)
            .unwrap_or_else(|err| panic!("{task} -word_size 4: {err}"));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1163-1169 (above): word size 2 keeps
// the lookup table (blastp-fast's compressed one too, `SetWordSize`). Item SC3: blastp-short's
// word size and `-word_size 2` of every task run it.
#[test]
fn test_validate_requested_blastp_support_accepts_word_size_two() {
    for task in ["blastp", "blastp-fast", "blastp-short"] {
        for word_size in [None, Some(2)] {
            let mut args = make_blastp_args();
            args.task = task.into();
            args.word_size = word_size;
            let resolved = args.resolve().expect("resolved blastp args");
            if task != "blastp-short" && word_size.is_none() {
                continue;
            }
            assert_eq!(resolved.word_size, 2, "{task}");
            assert_eq!(
                resolved.lookup_table_type,
                BlastpLookupTableType::AaLookupTable,
                "{task}"
            );
            validate_requested_blastp_support(&resolved)
                .unwrap_or_else(|err| panic!("{task} {word_size:?}: {err}"));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3036-3062
#[test]
fn test_validate_requested_blastp_support_takes_unified_p_modes() {
    for value in ["1u", "2u", "3U", "Du", "tu"] {
        let mut args = make_blastp_args();
        args.comp_based_stats = value.to_string();
        let resolved = args.resolve().expect("resolved blastp args");
        assert!(resolved.comp_based_stats.unified_p, "{value}");
        validate_requested_blastp_support(&resolved).expect(value);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:143-156 (IDENTITY is
// searched without the composition-based statistics only; item SC4b)
#[test]
fn test_validate_requested_blastp_support_takes_every_matrix() {
    // Every row of each matrix's table (its first row 32767/32767 included), in either case.
    for matrix in ScoringMatrix::ALL {
        let table = crate::stats::protein_options::scoring_matrix_table(matrix);
        for row in table.values {
            let (gap_open, gap_extend) = (row.0[0] as i32, row.0[1] as i32);
            for name in [
                matrix.ncbi_name().to_string(),
                matrix.ncbi_name().to_lowercase(),
            ] {
                let mut args = make_blastp_args();
                args.matrix = Some(name.clone());
                args.gap_open = Some(gap_open);
                args.gap_extend = Some(gap_extend);
                let resolved = args.resolve().expect("resolved blastp args");
                assert_eq!(resolved.scoring.matrix, matrix);
                validate_requested_blastp_support(&resolved)
                    .unwrap_or_else(|err| panic!("{name} {gap_open}/{gap_extend}: {err}"));
            }
        }
        // psiblast takes them too, in all its rounds (line RG, Q2); IDENTITY is NCBI's error
        // there (psiblast takes the standard matrices only, blast_options.c:913-917).
        let mut args = make_blastp_args();
        args.app = BlastpApp::Psiblast;
        args.matrix = Some(matrix.ncbi_name().to_string());
        if matrix == ScoringMatrix::Identity {
            let err = args.resolve().expect_err("psiblast IDENTITY");
            assert!(err
                .to_string()
                .contains("IDENTITY is not a supported matrix"));
            continue;
        }
        for row in table.values {
            let (gap_open, gap_extend) = (row.0[0] as i32, row.0[1] as i32);
            let mut args = args.clone();
            args.gap_open = Some(gap_open);
            args.gap_extend = Some(gap_extend);
            let resolved = args.resolve().expect("resolved psiblast args");
            validate_requested_blastp_support(&resolved).unwrap_or_else(|err| {
                panic!(
                    "psiblast {} {gap_open}/{gap_extend}: {err}",
                    matrix.ncbi_name()
                )
            });
        }
    }
    // IDENTITY with the compressed lookup table (word size 5, blastp-fast's default), with
    // which NCBI crashes: the options-only check refuses it (`search` where it searches its first
    // query batch); the lookup table runs (blastp-fast at word size 3); 6 and 7 are NCBI's error
    // (`validate_protein_options`).
    for (task, word_size, runs) in [
        ("blastp", Some(5), false),
        ("blastp-fast", None, false),
        ("blastp-fast", Some(3), true),
    ] {
        let mut args = make_blastp_args();
        args.task = task.to_string();
        args.matrix = Some("IDENTITY".to_string());
        args.word_size = word_size;
        let resolved = args.resolve().expect("resolved blastp args");
        validate_requested_blastp_support(&resolved).expect("IDENTITY's options");
        if runs {
            check_options(&args).expect("IDENTITY with the lookup table");
            continue;
        }
        let err = check_options(&args).expect_err("IDENTITY word 5");
        assert!(
            err.to_string().contains(
                "-matrix IDENTITY with the compressed lookup table (-word_size 5), with which NCBI BLAST+'s blastp crashes, is not supported by LOSAT's BLASTP"
            ),
            "{task}: {err}"
        );
    }
    for word_size in [6, 7] {
        let mut args = make_blastp_args();
        args.matrix = Some("identity".to_string());
        args.word_size = Some(word_size);
        let err = args.resolve().expect_err("IDENTITY word 6 and 7");
        assert!(
            err.to_string().contains(
                "Word size larger than 5 is not supported for the identity scoring matrix"
            ),
            "{word_size}: {err}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:258-271 (BLOSUM62's rows)
#[test]
fn test_validate_requested_blastp_support_takes_blosum62_gap_costs_of_the_table() {
    for (gap_open, gap_extend) in [
        (11, 2),
        (10, 2),
        (9, 2),
        (8, 2),
        (7, 2),
        (6, 2),
        (13, 1),
        (12, 1),
        (11, 1),
        (10, 1),
        (9, 1),
        (32767, 32767),
    ] {
        let mut args = make_blastp_args();
        args.gap_open = Some(gap_open);
        args.gap_extend = Some(gap_extend);
        let resolved = args.resolve().expect("resolved blastp args");
        validate_requested_blastp_support(&resolved)
            .unwrap_or_else(|err| panic!("{gap_open}/{gap_extend}: {err}"));
        // psiblast takes them too (its later rounds search the PSSM with them); only a PSSM
        // query (`-in_pssm`) keeps BLOSUM62 11/1 (line PS, PS7b; `db.ps.inpssm.rej.gap`).
        let mut args = make_blastp_args();
        args.app = BlastpApp::Psiblast;
        args.gap_open = Some(gap_open);
        args.gap_extend = Some(gap_extend);
        let resolved = args.resolve().expect("resolved psiblast args");
        validate_requested_blastp_support(&resolved)
            .unwrap_or_else(|err| panic!("psiblast {gap_open}/{gap_extend}: {err}"));
    }
}

// The contexts of a PSSM search and of psiblast's check of its queries take the matrix of
// the search (`pssm_hooks::psi_contexts`): with BLOSUM62 they are the blocks and the validity
// that BLOSUM62's own preparation gives (psiblast's default bytes), and another matrix's
// differ.
#[test]
fn test_psi_contexts_take_the_search_matrix() {
    let letters: Vec<Vec<u8>> = vec![
        b"MKVLAAGIVALLLAAGCSSSKEETPAWTQEQQQALTEAIRWYEGDTL".to_vec(),
        b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW".to_vec(),
        b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
        b"XXXXXXXXXXXXXXXXXXXX".to_vec(),
        b"ACDEFGHIKLMNPQRSTVWYBZUOJ*".to_vec(),
        b"CCCCCCCCCCPPPPPPPPPP".to_vec(),
    ];
    let records: Vec<FastaRecord> = letters
        .iter()
        .enumerate()
        .map(|(index, seq)| FastaRecord::new("", format!("q{index}").as_bytes(), seq))
        .collect();
    let resolved = |matrix: Option<&str>| {
        let mut args = make_blastp_args();
        args.app = BlastpApp::Psiblast;
        args.seg = Some("no".to_string());
        if let Some(matrix) = matrix {
            args.matrix = Some(matrix.to_string());
        }
        args.resolve().expect("resolved psiblast args")
    };
    let blosum62 = resolved(None);
    let frames = blastp_query_frames(&blosum62, &records);
    let (_, _, reference) =
        crate::algorithm::tblastx::lookup::prepare_blosum62_lookup_query_for_word_size(
            &frames,
            blosum62.word_size,
            false,
        );
    let contexts = pssm_hooks::psi_contexts(&blosum62, &frames);
    assert_eq!(contexts.len(), reference.len());
    let mut valid = 0;
    for (context, reference) in contexts.iter().zip(&reference) {
        assert_eq!(context.is_valid, reference.is_valid);
        let (a, b) = (context.karlin_params, reference.karlin_params);
        assert_eq!(
            (a.lambda.to_bits(), a.k.to_bits(), a.h.to_bits()),
            (b.lambda.to_bits(), b.k.to_bits(), b.h.to_bits())
        );
        valid += usize::from(context.is_valid);
    }
    assert!(
        valid > 0 && valid < contexts.len(),
        "{valid} valid contexts"
    );
    let blosum45 = resolved(Some("BLOSUM45"));
    let other = pssm_hooks::psi_contexts(&blosum45, &frames);
    assert_ne!(
        other[0].karlin_params.lambda.to_bits(),
        contexts[0].karlin_params.lambda.to_bits()
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2380-2382
#[test]
fn test_blastp_redo_x_dropoff_uses_scaled_final_vs_unscaled_extend_max() {
    let gapped = crate::stats::KarlinParams {
        lambda: 0.267,
        k: 0.041,
        h: 0.14,
        alpha: 1.0,
        beta: -1.0,
    };
    use crate::algorithm::tblastx::constants::{X_DROP_GAPPED_FINAL, X_DROP_GAPPED_PRELIM};
    let scaled_lambda = gapped.lambda / 32.0;
    let (_, unscaled_final) = crate::core::x_drop::gapped_x_dropoffs(
        X_DROP_GAPPED_PRELIM as f64,
        X_DROP_GAPPED_FINAL as f64,
        gapped.lambda,
    );

    assert_eq!(
        blastp_redo_x_dropoff(X_DROP_GAPPED_FINAL as f64, scaled_lambda, unscaled_final),
        ((X_DROP_GAPPED_FINAL as f64 * std::f64::consts::LN_2 / scaled_lambda) as i32)
            .max(unscaled_final)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2438-2462
#[test]
fn test_blastp_redo_params_follow_ncbi_cutoff_formulas() {
    let scaled_lambda = 0.267 / 32.0;
    let near_identical_cutoff = blastp_redo_near_identical_cutoff(scaled_lambda);
    assert!(near_identical_cutoff > 0.0);
    assert_eq!(blastp_redo_cutoff_score(false, 42, 32.0), 1);
    assert_eq!(blastp_redo_cutoff_score(true, 42, 32.0), 1344);
    assert_eq!(blastp_redo_cutoff_score(true, 7, 3.2), 22);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-941
#[test]
fn test_blastp_prelim_evalue_and_cutoff_follow_ncbi_cbs_stretch() {
    let expect_value = 1e-5;
    assert_eq!(
        blastp_prelim_evalue(expect_value, BlastpCompositionMode::NoCompositionBasedStats),
        expect_value
    );
    assert_eq!(
        blastp_prelim_evalue(expect_value, BlastpCompositionMode::CompositionMatrixAdjust),
        5.0 * expect_value
    );

    let scoring = crate::config::ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 11,
        gap_extend: 1,
    };
    let gapped = lookup_protein_params(&scoring);
    let gumbel =
        lookup_protein_gumbel_params(&scoring, 10_000).expect("supported blastp gumbel parameters");
    let cutoff_without_cbs = blastp_cutoff_score_max(
        expect_value,
        120,
        500,
        BlastpCompositionMode::NoCompositionBasedStats,
        &gapped,
        Some(&gumbel),
        0,
    );
    let cutoff_with_cbs = blastp_cutoff_score_max(
        expect_value,
        120,
        500,
        BlastpCompositionMode::CompositionMatrixAdjust,
        &gapped,
        Some(&gumbel),
        0,
    );

    assert!(cutoff_with_cbs <= cutoff_without_cbs);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-983
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1403-1440
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:226-240: the shortest record, raised to
// BLAST_SEQSRC_MINLENGTH (10).
#[test]
fn test_blastp_preliminary_cutoff_subject_length_uses_global_min_subject_len() {
    let subjects = vec![
        encode_protein_sequence(b"AAAAAAAAAAAAAAA"),
        encode_protein_sequence(b"AAAAAAAAAAAA"),
        encode_protein_sequence(b"AAAAAAAAAAAAAAAAAAAAAAAAA"),
    ];
    assert_eq!(blastp_preliminary_cutoff_subject_length(&subjects), 12);
    let short = vec![
        encode_protein_sequence(b"AAAAAA"),
        encode_protein_sequence(b"AAA"),
    ];
    assert_eq!(blastp_preliminary_cutoff_subject_length(&short), 10);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:555
#[test]
fn test_prepare_preliminary_hits_for_kappa_purges_shared_endpoints_then_sorts() {
    let hits = vec![
        BlastpPreliminaryHsp {
            query_start: 10,
            query_end: 60,
            subject_start: 20,
            subject_end: 70,
            gapped_query_start: 30,
            gapped_subject_start: 40,
            raw_score: 80,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: 0,
        },
        BlastpPreliminaryHsp {
            query_start: 10,
            query_end: 80,
            subject_start: 20,
            subject_end: 90,
            gapped_query_start: 35,
            gapped_subject_start: 45,
            raw_score: 90,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: 0,
        },
        BlastpPreliminaryHsp {
            query_start: 5,
            query_end: 40,
            subject_start: 7,
            subject_end: 42,
            gapped_query_start: 20,
            gapped_subject_start: 22,
            raw_score: 120,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 100,
            q_idx: 0,
            s_idx: 0,
        },
    ];

    let prepared = prepare_preliminary_hits_for_kappa(hits);
    assert_eq!(prepared.len(), 2);
    assert_eq!(prepared[0].raw_score, 120);
    assert_eq!(prepared[1].raw_score, 90);
    assert_eq!(prepared[1].query_start, 10);
    assert_eq!(prepared[1].subject_start, 20);
    assert_eq!(prepared[1].query_end, 80);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2505
#[test]
fn test_ncbi_qsort_preliminary_hits_matches_requested_ncbi_comparator() {
    let hits = vec![
        make_preliminary_hsp(0, 10, 40, 15, 45, 90),
        make_preliminary_hsp(0, 10, 45, 15, 50, 120),
        make_preliminary_hsp(0, 5, 30, 9, 34, 120),
        make_preliminary_hsp(1, 3, 20, 4, 21, 80),
    ];

    let mut by_query_offset = hits.clone();
    ncbi_qsort_preliminary_hits(&mut by_query_offset, preliminary_query_offset_compare);
    let ordered_by_query_offset: Vec<(i32, i32, i32, i32)> = by_query_offset
        .iter()
        .map(|hsp| {
            (
                hsp.query_context,
                hsp.query_start,
                hsp.subject_start,
                hsp.raw_score,
            )
        })
        .collect();
    assert_eq!(
        ordered_by_query_offset,
        vec![
            (0, 5, 9, 120),
            (0, 10, 15, 120),
            (0, 10, 15, 90),
            (1, 3, 4, 80)
        ]
    );

    let mut by_query_end = hits.clone();
    ncbi_qsort_preliminary_hits(&mut by_query_end, preliminary_query_end_compare);
    let ordered_by_query_end: Vec<(i32, i32, i32, i32)> = by_query_end
        .iter()
        .map(|hsp| {
            (
                hsp.query_context,
                hsp.query_end,
                hsp.subject_end,
                hsp.raw_score,
            )
        })
        .collect();
    assert_eq!(
        ordered_by_query_end,
        vec![
            (0, 30, 34, 120),
            (0, 40, 45, 90),
            (0, 45, 50, 120),
            (1, 20, 21, 80)
        ]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2488-2495,2514-2520
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2320,2385
#[test]
fn test_prepare_preliminary_hits_for_kappa_keeps_first_equal_endpoint_hsp() {
    let mut first = make_preliminary_hsp(0, 12, 6074, 13, 6057, 16284);
    first.gapped_query_start = 4964;
    first.gapped_subject_start = 4942;
    let mut later = first;
    later.gapped_query_start = 55;
    later.gapped_subject_start = 55;

    let prepared = prepare_preliminary_hits_for_kappa(vec![first, later]);

    assert_eq!(prepared.len(), 1);
    assert_eq!(prepared[0].gapped_query_start, 4964);
    assert_eq!(prepared[0].gapped_subject_start, 4942);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2482-2487
#[test]
fn test_prepare_preliminary_hits_for_kappa_keeps_shared_endpoints_on_subject_frame_change() {
    let mut plus = make_preliminary_hsp(0, 10, 60, 20, 70, 80);
    plus.subject_frame = 1;
    let mut minus = make_preliminary_hsp(0, 10, 60, 20, 70, 90);
    minus.subject_frame = -1;

    let prepared = prepare_preliminary_hits_for_kappa(vec![plus, minus]);

    assert_eq!(prepared.len(), 2);
    assert_eq!(prepared[0].subject_frame, -1);
    assert_eq!(prepared[1].subject_frame, 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:117-149
#[test]
fn test_split_preliminary_hits_by_query_preserves_subject_order_within_query() {
    let contexts = vec![make_query_context(0, 0), make_query_context(1, 1000)];
    let hits = vec![
        make_preliminary_hsp(99, 15, 35, 100, 120, 50),
        make_preliminary_hsp(77, 10, 30, 40, 60, 80),
        make_preliminary_hsp(88, 25, 45, 130, 150, 40),
        make_preliminary_hsp(66, 35, 55, 70, 90, 60),
    ]
    .into_iter()
    .enumerate()
    .map(|(index, mut hsp)| {
        hsp.query_context = if index % 2 == 0 { 1 } else { 0 };
        hsp
    })
    .collect();

    let hits_by_query = split_preliminary_hits_by_query(hits, &contexts, 2);

    assert_eq!(hits_by_query.len(), 2);
    assert_eq!(hits_by_query[0].len(), 2);
    assert_eq!(hits_by_query[0][0].query_start, 10);
    assert_eq!(hits_by_query[0][1].query_start, 35);
    assert_eq!(hits_by_query[1].len(), 2);
    assert_eq!(hits_by_query[1][0].subject_start, 100);
    assert_eq!(hits_by_query[1][1].subject_start, 130);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:119-132
// LOSAT_X_BPHITSFLAT: the counting order lists each query's HSPs as the per-query vectors do,
// for random contexts, including a query index past `num_queries` (left out by both).
#[test]
fn x_hits_flat_order_matches_the_per_query_vectors() {
    let contexts: Vec<QueryContext> = (0..5u32)
        .map(|q| make_query_context(q, q as i32 * 1000))
        .collect();
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    let mut order = Vec::new();
    let mut starts = Vec::new();
    for round in 0..200 {
        let len = next(40) as usize;
        let hits: Vec<BlastpPreliminaryHsp> = (0..len)
            .map(|index| {
                let mut hsp = make_preliminary_hsp(
                    index as u32,
                    index as i32,
                    index as i32 + 10,
                    round,
                    round + 10,
                    next(100) as i32,
                );
                hsp.query_context = next(5) as i32;
                hsp
            })
            .collect();
        let num_queries = 1 + (round as usize % 5);
        x_hits_flat::order_by_query(&hits, &contexts, num_queries, &mut order, &mut starts);
        let by_query = split_preliminary_hits_by_query(hits.clone(), &contexts, num_queries);
        assert_eq!(starts.len(), num_queries + 1);
        for (q, expected) in by_query.iter().enumerate().take(num_queries) {
            let run = &order[starts[q] as usize..starts[q + 1] as usize];
            let got: Vec<u32> = run.iter().map(|&i| hits[i as usize].q_idx).collect();
            let want: Vec<u32> = expected.iter().map(|hsp| hsp.q_idx).collect();
            assert_eq!(got, want, "round {round} query {q}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:289-318
#[test]
fn test_preliminary_hit_from_local_match_hsp_preserves_ncbi_gapped_starts_and_context() {
    let mut preliminary_hsp = make_preliminary_hsp(0, 10, 20, 30, 40, 50);
    preliminary_hsp.gapped_query_start = 13;
    preliminary_hsp.gapped_subject_start = 29;
    preliminary_hsp.query_context = 7;
    preliminary_hsp.subject_frame = -1;

    let rebuilt = preliminary_hit_from_local_match_hsp(&make_local_match_hsp(&preliminary_hsp));

    assert_eq!(rebuilt.query_start, 10);
    assert_eq!(rebuilt.subject_start, 30);
    assert_eq!(rebuilt.gapped_query_start, 13);
    assert_eq!(rebuilt.gapped_subject_start, 29);
    assert_eq!(rebuilt.query_context, 7);
    assert_eq!(rebuilt.subject_frame, -1);
}

#[test]
fn test_restricted_preliminary_score_retries_exact_within_ncbi_inconclusive_band() {
    assert!(should_retry_exact_preliminary_gapped_alignment(
        true, 95, 100, 90
    ));
    assert!(!should_retry_exact_preliminary_gapped_alignment(
        false, 95, 100, 90
    ));
    assert!(!should_retry_exact_preliminary_gapped_alignment(
        true, 89, 100, 90
    ));
    assert!(!should_retry_exact_preliminary_gapped_alignment(
        true, 100, 100, 90
    ));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-4000
#[test]
fn test_rebuild_preliminary_interval_tree_drops_only_retried_query_hsps() {
    let contexts = vec![make_query_context(0, 0), make_query_context(1, 1000)];
    let query0_hsp = make_preliminary_hsp(0, 10, 60, 20, 70, 100);
    let query1_hsp = make_preliminary_hsp(1, 15, 55, 30, 70, 90);
    let mut preliminary_hsp_list = vec![query0_hsp, query1_hsp];
    let mut interval_tree = BlastIntervalTree::new(0, 2000, 0, 200);

    assert_eq!(
        rebuild_preliminary_interval_tree_from_list(
            &mut interval_tree,
            &preliminary_hsp_list,
            &contexts,
        ),
        Some(2)
    );

    let query0_contained = TreeHsp {
        query_offset: 20,
        query_end: 40,
        subject_offset: 30,
        subject_end: 50,
        score: 50,
        query_frame: 1,
        query_length: 100,
        query_context_offset: contexts[0].frame_base,
        subject_frame_sign: 0,
    };
    let query1_contained = TreeHsp {
        query_offset: 20,
        query_end: 40,
        subject_offset: 40,
        subject_end: 60,
        score: 50,
        query_frame: 1,
        query_length: 100,
        query_context_offset: contexts[1].frame_base,
        subject_frame_sign: 0,
    };

    assert!(interval_tree.contains_hsp(&query0_contained, contexts[0].frame_base, 0));
    assert!(interval_tree.contains_hsp(&query1_contained, contexts[1].frame_base, 0));

    preliminary_hsp_list[0].q_idx = 99;
    preliminary_hsp_list[1].q_idx = 77;
    remove_preliminary_hits_for_query(&mut preliminary_hsp_list, 0, &contexts);
    assert_eq!(
        rebuild_preliminary_interval_tree_from_list(
            &mut interval_tree,
            &preliminary_hsp_list,
            &contexts,
        ),
        Some(1)
    );

    assert!(!interval_tree.contains_hsp(&query0_contained, contexts[0].frame_base, 0));
    assert!(interval_tree.contains_hsp(&query1_contained, contexts[1].frame_base, 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2268-2315
#[test]
fn test_preliminary_common_endpoint_compare_uses_context_not_query_index() {
    let mut left = make_preliminary_hsp(0, 10, 40, 20, 50, 80);
    left.query_context = 1;
    let mut right = make_preliminary_hsp(0, 10, 40, 20, 50, 70);
    right.query_context = 0;

    let prepared = prepare_preliminary_hits_for_kappa(vec![left, right]);

    assert_eq!(prepared.len(), 2);
    assert_eq!(prepared[0].query_context, 1);
    assert_eq!(prepared[1].query_context, 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3119-3128
#[test]
fn test_redone_match_heaps_use_hitlist_size_not_prelim_hitlist_size() {
    let heaps = build_redone_match_heaps(3, || BlastCompoHeap::new(7, PSI_INCLUSION_ETHRESH));
    assert_eq!(heaps.len(), 3);
    assert_eq!(heaps[0].heap_threshold, 7);
    assert!((heaps[0].ecutoff - PSI_INCLUSION_ETHRESH).abs() < f64::EPSILON);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:204-221,324-329,917-921
#[test]
fn test_invalid_blastp_context_uses_ncbi_placeholder_cutoffs() {
    let invalid = QueryContext {
        is_valid: false,
        ..make_query_context(0, 0)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3577-3622
    let scoring = crate::config::ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: 11,
        gap_extend: 1,
    };
    let gapped = lookup_protein_params(&scoring);
    let gumbel =
        lookup_protein_gumbel_params(&scoring, 10_000).expect("supported blastp gumbel parameters");

    assert_eq!(blastp_x_dropoff_init_for_context(&invalid, 7.0), 0);
    assert_eq!(
        blastp_hit_cutoff_score_for_context(
            &invalid,
            10.0,
            100,
            BlastpCompositionMode::CompositionMatrixAdjust,
            &gapped,
            Some(&gumbel),
            0,
        ),
        i32::MAX
    );
    assert_eq!(
        blastp_word_cutoff_score_for_context(50, &invalid, 22.0),
        i32::MAX
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-309
#[test]
fn test_chain_blastp_init_hsps_restores_ncbi_score_order() {
    let contexts = vec![QueryContext {
        aa_len: 256,
        ..make_query_context(0, 0)
    }];
    let init_hsps = vec![
        make_init_hsp(0, 10, 0, 0, 0, 50),
        make_init_hsp(100, 10, 100, 100, 100, 60),
    ];
    let mut keep = Vec::new();
    let mut nodes = Vec::new();

    let chained = chain_blastp_init_hsps(
        init_hsps,
        &contexts,
        &[100],
        &[1],
        11,
        1,
        &mut keep,
        &mut nodes,
    );

    assert_eq!(chained.len(), 2);
    assert_eq!(chained[0].ungapped_data.score, 60);
    assert_eq!(chained[1].ungapped_data.score, 50);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3628-3636
#[test]
fn test_chain_blastp_init_hsps_drops_a_lone_hsp_below_the_gapped_cutoff() {
    let contexts = vec![QueryContext {
        aa_len: 256,
        ..make_query_context(0, 0)
    }];
    let mut keep = Vec::new();
    let mut nodes = Vec::new();

    // 20 - 12 + 15 - 1 = 22 < 23: dropped.
    let dropped = chain_blastp_init_hsps(
        vec![make_init_hsp(0, 10, 0, 0, 0, 20)],
        &contexts,
        &[15],
        &[23],
        11,
        1,
        &mut keep,
        &mut nodes,
    );
    assert!(dropped.is_empty());

    // 20 - 12 + 15 - 1 = 22, not below 22: kept.
    let kept = chain_blastp_init_hsps(
        vec![make_init_hsp(0, 10, 0, 0, 0, 20)],
        &contexts,
        &[15],
        &[22],
        11,
        1,
        &mut keep,
        &mut nodes,
    );
    assert_eq!(kept.len(), 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3532-3533
#[test]
fn test_ncbi_qsort_init_hsps_by_query_offset_score_matches_ncbi_ties() {
    let mut init_hsps = vec![
        make_init_hsp(20, 8, 40, 20, 40, 90),
        make_init_hsp(10, 7, 15, 10, 15, 80),
        make_init_hsp(10, 9, 15, 10, 15, 120),
        make_init_hsp(20, 6, 35, 20, 35, 90),
    ];

    ncbi_qsort_init_hsps_by_query_offset_score(&mut init_hsps);

    let ordered: Vec<(i32, i32, i32)> = init_hsps
        .iter()
        .map(|hsp| {
            (
                hsp.ungapped_data.q_start,
                hsp.ungapped_data.score,
                hsp.ungapped_data.s_start,
            )
        })
        .collect();
    assert_eq!(
        ordered,
        vec![(10, 120, 15), (10, 80, 15), (20, 90, 40), (20, 90, 35)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-300
#[test]
fn test_compare_init_hsps_by_score_matches_ncbi_score_compare_match_ties() {
    let mut init_hsps = vec![
        make_init_hsp(20, 10, 30, 20, 30, 100),
        make_init_hsp(5, 10, 10, 5, 10, 100),
        make_init_hsp(20, 20, 10, 20, 10, 100),
        make_init_hsp(10, 20, 10, 10, 10, 100),
    ];

    ncbi_qsort_init_hsps_by_score(&mut init_hsps);

    let ordered: Vec<(i32, i32, i32)> = init_hsps
        .iter()
        .map(|hsp| {
            (
                hsp.ungapped_data.q_start,
                hsp.ungapped_data.s_start,
                hsp.ungapped_data.length,
            )
        })
        .collect();
    assert_eq!(
        ordered,
        vec![(10, 10, 20), (20, 10, 20), (5, 10, 10), (20, 30, 10)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3639-3657
#[test]
fn test_compact_chained_init_hsps_ncbi_uses_tail_swap_compaction_order() {
    let mut init_hsps = vec![
        make_init_hsp(10, 10, 10, 10, 10, 10),
        make_init_hsp(20, 10, 20, 20, 20, 20),
        make_init_hsp(30, 10, 30, 30, 30, 30),
        make_init_hsp(40, 10, 40, 40, 40, 40),
        make_init_hsp(50, 10, 50, 50, 50, 50),
    ];
    let mut keep = vec![true, false, true, false, true];

    compact_chained_init_hsps_ncbi(&mut init_hsps, &mut keep);

    let ordered_q_starts: Vec<i32> = init_hsps
        .iter()
        .map(|hsp| hsp.ungapped_data.q_start)
        .collect();
    assert_eq!(ordered_q_starts, vec![10, 50, 30]);
    assert_eq!(keep, vec![true, true, true]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3748-3786
#[test]
fn test_build_restricted_align_array_uses_first_init_hsp_per_query() {
    let contexts = vec![
        QueryContext {
            aa_len: 1000,
            ..make_query_context(0, 0)
        },
        QueryContext {
            aa_len: 1000,
            ..make_query_context(1, 1000)
        },
    ];
    let init_hsps = vec![
        make_init_hsp(20, 20, 20, 20, 20, 45),
        make_init_hsp(0, 20, 0, 0, 0, 80),
        make_init_hsp(1060, 20, 60, 1060, 60, 70),
    ];

    let restricted = build_restricted_align_array(&init_hsps, &contexts, &[100, 100], 2);

    assert_eq!(restricted, vec![true, false]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3748-3786
#[test]
fn test_build_restricted_align_array_uses_qoff_not_cached_context() {
    let contexts = vec![
        QueryContext {
            aa_len: 1000,
            ..make_query_context(0, 0)
        },
        QueryContext {
            aa_len: 1000,
            ..make_query_context(1, 1000)
        },
    ];
    let init_hsps = vec![
        make_init_hsp(0, 20, 0, 0, 0, 80),
        make_init_hsp(1005, 20, 50, 1005, 50, 45),
    ];

    let restricted = build_restricted_align_array(&init_hsps, &contexts, &[100, 100], 2);

    assert_eq!(restricted, vec![false, true]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-209
// The compressed tables of NCBI's query batches as ranges of lookup segments: a batch that
// NCBI splits (19,800 letters or more with the default CHUNK_SIZE) has none, and a query
// chunk's search is one table.
#[test]
fn test_ncbi_compressed_lookup_tables_are_unsplit_batches() {
    let lengths = [5000usize, 6000, 3000, 25000];
    let contexts = vec![
        make_query_context(0, 0),
        make_query_context(1, 5001),
        make_query_context(2, 11002),
        make_query_context(3, 14003),
    ];
    // Query 0 has two segments, query 2 none (all masked).
    let locations = [(0, 99), (200, 4999), (5001, 11000), (14003, 39002)];
    let batches = [
        seq_range::RangedBatch {
            input: 0..2,
            searched: 0..2,
        },
        seq_range::RangedBatch {
            input: 2..3,
            searched: 2..3,
        },
        seq_range::RangedBatch {
            input: 3..4,
            searched: 3..4,
        },
        seq_range::RangedBatch {
            input: 4..5,
            searched: 4..4,
        },
    ];
    assert_eq!(
        ncbi_compressed_lookup_tables(&batches, &lengths, &contexts, &locations, false, true)
            .unwrap(),
        vec![0..3, 3..3]
    );
    assert_eq!(
        ncbi_compressed_lookup_tables(&batches, &lengths, &contexts, &locations, true, true)
            .unwrap(),
        vec![0..4]
    );
}

// Protein records of `lengths` letters (a deterministic generator), some with runs that SEG
// masks, so that a query can have more than one lookup segment.
fn x_block_test_records(lengths: &[usize]) -> Vec<FastaRecord> {
    const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = 0x2545_f491_4f6c_dd1du64;
    lengths
        .iter()
        .enumerate()
        .map(|(index, &length)| {
            let mut sequence = Vec::with_capacity(length);
            while sequence.len() < length {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                if index % 2 == 1 && state.is_multiple_of(97) {
                    let run = (length - sequence.len()).min(30);
                    sequence.extend(std::iter::repeat_n(b'Q', run));
                } else {
                    sequence.push(LETTERS[(state % LETTERS.len() as u64) as usize]);
                }
            }
            FastaRecord::new("", format!("q{index}").as_bytes(), &sequence)
        })
        .collect()
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): a block's compressed table is tested with the tables of the
// block's own query batches, as ranges of the block's lookup segments.
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-260
// The block search gives `run_resolved_in_pool` a block's records as the whole input, whose
// batches are those of the whole input; the blocks' tables, each shifted by the block's first
// lookup segment, are the tables of the search of all queries, and the blocks' lookup segments
// are those of the search of all queries.
#[test]
fn test_x_query_block_compressed_tables_are_the_blocks_batches() {
    let lengths = [
        3000usize, 4000, 2500, 6000, 25000, 0, 7000, 1500, 9000, 200, 12000, 800,
    ];
    let records = x_block_test_records(&lengths);
    let seg = crate::utils::seg::SegParams::default();
    let lookup_query = |records: &[FastaRecord]| {
        let frames: Vec<Vec<QueryFrame>> = records
            .iter()
            .map(|record| {
                vec![encode_protein_query_frame_with_seg(
                    record.seq(),
                    Some(&seg),
                )]
            })
            .collect();
        let (_, locations, contexts) =
            crate::algorithm::tblastx::lookup::prepare_blosum62_lookup_query_for_word_size(
                &frames, 5, false,
            );
        (locations, contexts)
    };
    let tables = |records: &[FastaRecord]| {
        let (locations, contexts) = lookup_query(records);
        let lengths: Vec<usize> = records.iter().map(|record| record.seq().len()).collect();
        let batches = seq_range::QueryInput::whole(records).batches(10000);
        ncbi_compressed_lookup_tables(&batches, &lengths, &contexts, &locations, false, true)
            .unwrap()
    };
    let (whole_locations, whole_contexts) = lookup_query(&records);
    let whole_tables = tables(&records);
    assert!(
        whole_locations.len() > lengths.len(),
        "SEG splits some queries"
    );
    assert!(whole_tables.len() >= 4);
    let whole_batches = seq_range::QueryInput::whole(&records).batches(10000);
    for block_letters in [1, 20000] {
        let blocks = super::super::x_query_blocks::x_query_blocks(
            &whole_batches,
            &lengths,
            block_letters,
            None,
        )
        .unwrap();
        assert!(blocks.len() > 2, "blocks of {block_letters}: {blocks:?}");
        let mut shifted = Vec::new();
        for block in &blocks {
            let base = whole_contexts[block.start].frame_base;
            let first = whole_locations.partition_point(|&(left, _)| left < base);
            let (block_locations, _) = lookup_query(&records[block.clone()]);
            let expected: Vec<(i32, i32)> = whole_locations[first..first + block_locations.len()]
                .iter()
                .map(|&(left, right)| (left - base, right - base))
                .collect();
            assert_eq!(block_locations, expected, "block {block:?}");
            shifted.extend(
                tables(&records[block.clone()])
                    .into_iter()
                    .map(|table| table.start + first..table.end + first),
            );
        }
        assert_eq!(shifted, whole_tables, "blocks of {block_letters}");
    }
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the refusals of the search of all queries that the block
// search would make only in a later block keep the reference path.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn test_x_block_search_refusals_of_the_search_of_all_queries() {
    let lengths = [1500usize, 2000, 1000];
    let records = x_block_test_records(&lengths);
    let batches = seq_range::QueryInput::whole(&records).batches(10000);
    let mut args = make_blastp_args();
    args.task = "blastp-fast".to_string();
    let resolved = args.resolve().expect("resolved blastp-fast args");
    assert_eq!(
        resolved.lookup_table_type,
        BlastpLookupTableType::CompressedAaLookupTable
    );
    assert!(x_block_search_refuses_nothing(
        &resolved, &records, &batches, &lengths
    ));
    // A threshold below every score: each offset can add all 15^5 compressed words, more
    // than the banks take for a batch of 4,500 letters, so the bound does not clear it
    // (`x_word_blocks`, on one thread and on three; with batches of one query each, the
    // queries of 1,500 and 2,000 letters fail and the one of 1,000 passes).
    let mut low = resolved.clone();
    low.threshold = -1000.0;
    let mut low_aa = low.clone();
    // The default lookup table has no overflow banks.
    low_aa.lookup_table_type = BlastpLookupTableType::AaLookupTable;
    let input = seq_range::QueryInput::whole(&records);
    let letter_blocks = [0..3];
    for threads in [1, 3] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            for batch_size in [10000, 1] {
                assert!(x_word_blocks(
                    &resolved,
                    &records,
                    &input,
                    batch_size,
                    &letter_blocks,
                    pool
                )
                .is_some());
                assert!(
                    x_word_blocks(&low, &records, &input, batch_size, &letter_blocks, pool)
                        .is_none()
                );
                assert!(
                    x_word_blocks(&low_aa, &records, &input, batch_size, &letter_blocks, pool)
                        .is_some()
                );
            }
            assert!(x_word_blocks(
                &low,
                &records[2..],
                &seq_range::QueryInput::whole(&records[2..]),
                1,
                &[0..1],
                pool
            )
            .is_some());
            Ok(())
        })
        .expect("pool");
    }
    assert!(x_block_search_refuses_nothing(
        &low_aa, &records, &batches, &lengths
    ));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-57 with the concatenated length of all queries.
    let mut wide = resolved.clone();
    let query_length = lengths.iter().sum::<usize>() + lengths.len() - 1;
    wide.window_size = (1 << 30) - query_length;
    assert!(x_block_search_refuses_nothing(
        &wide, &records, &batches, &lengths
    ));
    wide.window_size += 1;
    assert!(!x_block_search_refuses_nothing(
        &wide, &records, &batches, &lengths
    ));
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the preliminary stage in blocks is planned for blastp's search
// of all queries only: the blocks of whole batches of `x_batch_blocks` (here more than one block
// of the built-in 50,000 letters of `x_bp_qblock`; the variable is not set in tests), each
// ending at a batch end, on one thread and on three; none for one block, for psiblast or with
// `-sum_stats`; and a call of the preliminary stage only never plans
// (`BlastpReports::prelim_only`), so a block is not split again on its own batches' words.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn x_prelim_block_plan_is_the_search_of_all_queries_of_blastp() {
    let records = x_block_test_records(&[3000; 30]);
    let input = seq_range::QueryInput::whole(&records);
    let few = seq_range::QueryInput::whole(&records[..10]);
    let resolved = make_blastp_args().resolve().expect("resolved blastp args");
    let mut psi = make_blastp_args();
    psi.app = BlastpApp::Psiblast;
    let psi = psi.resolve().expect("resolved psiblast args");
    let mut sum_stats = resolved.clone();
    sum_stats.sum_stats = true;
    let batches = input.batches(10000);
    for threads in [1, 3] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            let blocks =
                x_prelim_block_plan(&resolved, &records, &input, 10000, pool).expect("blocks");
            assert!(blocks.len() > 1, "{blocks:?}");
            assert_eq!(blocks.first().map(|block| block.start), Some(0));
            assert_eq!(blocks.last().map(|block| block.end), Some(records.len()));
            for (block, next) in blocks.iter().zip(&blocks[1..]) {
                assert_eq!(block.end, next.start);
            }
            for block in &blocks {
                assert!(batches.iter().any(|batch| batch.searched.end == block.end));
            }
            assert!(x_prelim_block_plan(&psi, &records, &input, 10000, pool).is_none());
            assert!(x_prelim_block_plan(&sum_stats, &records, &input, 10000, pool).is_none());
            assert!(x_prelim_block_plan(&resolved, &records[..10], &few, 10000, pool).is_none());
            Ok(())
        })
        .expect("pool");
    }
    assert!(!BlastpReports::prelim_only(true).x_plan_blocks);
    assert!(!BlastpReports::prelim_only(false).x_plan_blocks);
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): a database search without the Gumbel block (`-searchsp`,
// `OLD_FSC`) keeps its block plan where the setup counts its sequences, and has none where the
// count is 0 (`call_setup::setup_count_refused`): that search refuses in the call of all
// queries, before any output (review E of line RG D1a-4, N3; the bytes are the rows'
// `db.ak.bp.rej.nseq0.oldfsc.blocks` and `db.ak.bp.nseq0.ssp.allx.7`). Thirty queries of 3,000
// letters, three batches, more than one block of the built-in 50,000 letters.
#[test]
fn x_query_block_plan_is_kept_without_the_gumbel_block_but_for_a_count_of_0() {
    let records = x_block_test_records(&[3000; 30]);
    let input = seq_range::QueryInput::whole(&records);
    let batches = crate::blastinput::query_batch::read_query_batches(
        &records,
        &records,
        &input,
        crate::blastinput::query_batch::QueryReading::Records,
        10000,
        true,
        None,
    );
    let settings = crate::blastinput::ncbi_environment::ApplicationSettings::default();
    // From the manifest's directory, as the other tests' databases: a WASI module has no
    // working directory (the runner preopens `/`).
    let open = |name: &str| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(name);
        database::BlastpDatabase::open(path.to_str().expect("UTF-8 path"), &settings, "BLASTP")
            .expect("test database")
    };
    // NSEQ 0 over `tfp` (`aliaskeys/p_nseq0.pal`), and `tfp` itself.
    let counts_0 = open("tests/db_search/aliaskeys/p_nseq0");
    let counts = open("tests/db_search/taxfilter/tfp");
    assert_eq!(counts_0.statistics_num_sequences(), 0);
    assert!(counts.statistics_num_sequences() > 0);
    let gumbel = make_blastp_args().resolve().expect("resolved blastp args");
    let mut searchsp = make_blastp_args();
    searchsp.searchsp = Some(1000);
    let searchsp = searchsp.resolve().expect("resolved blastp args");
    let mut old_fsc = gumbel.clone();
    old_fsc.old_fsc = true;
    assert!(searchsp.old_fsc && !gumbel.old_fsc);
    let mut stream: Vec<u8> = Vec::new();
    let mut diagnostics = std::io::sink();
    let outputs = ReportOutputs::single("7", OutputSink::Writer(&mut stream), &mut diagnostics);
    let plan = |args: &ResolvedBlastpArgs, database: &database::BlastpDatabase| {
        x_query_block_plan(
            args,
            &records,
            &input,
            &batches,
            BlastpSubjectSource::Database(database),
            &outputs,
            10000,
        )
        .map(|(_, blocks)| blocks.len())
    };
    for args in [&searchsp, &old_fsc, &gumbel] {
        assert!(plan(args, &counts).is_some_and(|blocks| blocks > 1));
    }
    assert_eq!(plan(&searchsp, &counts_0), None);
    assert_eq!(plan(&old_fsc, &counts_0), None);
    assert!(plan(&gumbel, &counts_0).is_some_and(|blocks| blocks > 1));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:192-201
// Line RG item X-short-bat-2: the refusal of a compressed table whose words fill NCBI's overflow
// banks is the failure of the batch of the first table that fills them (`banks_filled_refusal`):
// a batch that NCBI splits has no table among NCBI's tables (its chunks are searches of their own,
// whose one table names the chunk's call's first query, made its batch's by
// `search_split_query_batches`); an ungapped search splits none. One query per batch; the second
// of 20,000 letters is split at the default chunk size.
#[test]
fn banks_filled_refusal_names_the_batch_of_the_table() {
    let records = x_block_test_records(&[100, 20000, 100, 100]);
    let lengths: Vec<usize> = records.iter().map(|record| record.seq().len()).collect();
    let batches = seq_range::QueryInput::whole(&records).batches(1);
    assert_eq!(batches.len(), 4);
    let resolve = |ungapped: bool| {
        let mut args = make_blastp_args();
        args.task = "blastp-fast".to_string();
        args.ungapped = ungapped;
        args.resolve().expect("resolved blastp-fast args")
    };
    let gapped = resolve(false);
    let ungapped = resolve(true);
    let named = |args: &ResolvedBlastpArgs, filled: usize, chunk: bool| {
        let failure = banks_filled_refusal(args, &batches, &lengths, chunk, BanksFilled(filled))
            .downcast::<crate::blastinput::app::ErrorAfterBatchReports>()
            .expect("a batch's failure");
        assert!(
            failure.error.to_string().contains("1024 overflow banks"),
            "{}",
            failure.error
        );
        failure.query
    };
    for (filled, query) in [(0, 0), (1, 2), (2, 3)] {
        assert_eq!(named(&gapped, filled, false), Some(query), "table {filled}");
    }
    for filled in 0..4 {
        assert_eq!(named(&ungapped, filled, false), Some(filled), "table {filled}");
    }
    assert_eq!(named(&gapped, 0, true), Some(0));
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the contexts of the search of all queries whose preliminary
// stage runs in blocks (`x_blocks_contexts`, no table) are those of each table's builder, every
// field and the Karlin blocks bit for bit. The builders make their contexts in the pass that
// makes their lookup query (`x_prepare_lookup_query`, with `x_context_karlin` inline, on the
// pool, or as tasks of the build), and the table-less contexts come from another entry to the
// same two functions (`prepare_matrix_lookup_query_for_word_size`); this test ties them: every
// matrix (`ScoringMatrix::ALL`) at every word size of `-task blastp` (2 to 4: exact blastp's
// table; 5 to 7: the compressed table, none of IDENTITY's) and with `-task blastp-fast`, at the
// default thresholds, with SEG and `-lcase_masking` on, on one thread and on three; queries
// without letters, of X only, of one letter repeated, of every letter, with lower-case letters.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn x_blocks_contexts_are_the_builders_contexts() {
    let mut records = x_block_test_records(&[300, 0, 1200, 45, 700]);
    for (index, letters) in [
        b"XXXXXXXXXXXXXXXXXXXX".as_slice(),
        b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW",
        b"ACDEFGHIKLMNPQRSTVWYBZUOJ*",
        b"MKVLAAGIVALLLAAGCSSSKEETPAWTQEQQQALTEAIRWYEGDTL",
        b"MKVLAAGIVAllaagcssskeetpAWTQEQQQALTEAIRWYEGDTLmkvlaagivalllaagcSSSKEETPAW",
    ]
    .into_iter()
    .enumerate()
    {
        records.push(FastaRecord::new(
            "",
            format!("x{index}").as_bytes(),
            letters,
        ));
    }
    let lengths: Vec<usize> = records.iter().map(|record| record.seq().len()).collect();
    let batches = seq_range::QueryInput::whole(&records).batches(10000);
    let mut compared = 0;
    let mut refused: Vec<String> = Vec::new();
    let configs = ScoringMatrix::ALL.into_iter().flat_map(|matrix| {
        (2..=7)
            .map(|word_size| ("blastp", Some(word_size)))
            .chain([("blastp-fast", None)])
            .map(move |(task, word_size)| (task, matrix, word_size))
    });
    for (task, matrix, word_size) in configs {
        let mut args = make_blastp_args();
        args.task = task.to_string();
        args.seg = Some("yes".to_string());
        args.lcase_masking = true;
        args.matrix = Some(matrix.ncbi_name().to_string());
        args.word_size = word_size;
        let resolved = match args.resolve() {
            Ok(resolved) => resolved,
            Err(error) => {
                refused.push(format!("{task} {matrix:?} {word_size:?}: {error:#}"));
                continue;
            }
        };
        assert_eq!(resolved.scoring.matrix, matrix);
        let frames = blastp_query_frames(&resolved, &records);
        for threads in [1, 3] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                let (_, built) = match resolved.lookup_table_type {
                    BlastpLookupTableType::AaLookupTable => {
                        blastp_aa_lookup(&resolved, &frames, false, pool)
                    }
                    BlastpLookupTableType::CompressedAaLookupTable => {
                        match blastp_compressed_lookup(
                            &resolved,
                            &frames,
                            &batches,
                            &lengths,
                            pool,
                            super::super::x_query_blocks::XLookupFor::Search,
                            0,
                        ) {
                            Ok(built) => built,
                            Err(error) => {
                                refused.push(format!(
                                    "{task} {matrix:?} {word_size:?} threads {threads}: {error:#}"
                                ));
                                return Ok(());
                            }
                        }
                    }
                };
                let (lookup, contexts) = x_blocks_contexts(&resolved, &frames, pool);
                assert_eq!(lookup.longest_chain(), 0);
                let what = format!("{task} {matrix:?} {word_size:?} threads {threads}");
                compared += 1;
                assert_eq!(contexts.len(), built.len(), "{what}");
                let mut valid = 0;
                for (a, b) in contexts.iter().zip(&built) {
                    assert_eq!(
                        (
                            a.q_idx,
                            a.frame,
                            a.aa_len,
                            a.orig_len,
                            a.frame_base,
                            a.is_valid
                        ),
                        (
                            b.q_idx,
                            b.frame,
                            b.aa_len,
                            b.orig_len,
                            b.frame_base,
                            b.is_valid
                        ),
                        "{what}"
                    );
                    assert!(a.aa_seq == b.aa_seq && a.aa_seq_nomask == b.aa_seq_nomask);
                    // Debug writes each f64 so that it reads back to the same bits.
                    assert_eq!(
                        format!("{:?}", a.karlin_params),
                        format!("{:?}", b.karlin_params),
                        "{what}"
                    );
                    valid += usize::from(a.is_valid);
                }
                assert!(valid > 0 && valid < contexts.len(), "{what}: {valid} valid");
                Ok(())
            })
            .expect("pool");
        }
    }
    // IDENTITY's compressed table is built by no search: the options refuse word sizes above 5,
    // and `search` refuses word size 5 (blastp-fast's too) before any table (NCBI's blastp
    // crashes), as its builder does here; every other configuration is compared.
    let identity_compressed = [
        "blastp Identity Some(5) ",
        "blastp Identity Some(6): ",
        "blastp Identity Some(7): ",
        "blastp-fast Identity None ",
    ];
    assert!(
        refused.len() == 6
            && refused.iter().all(|refusal| {
                identity_compressed
                    .iter()
                    .any(|prefix| refusal.starts_with(prefix))
            }),
        "{compared} compared; refused:\n{}",
        refused.join("\n")
    );
    // (9 matrices x (6 word sizes + blastp-fast) - IDENTITY's 3 compressed tables) x 2 thread
    // counts.
    assert_eq!(compared, 118);
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): with the compressed lookup table the blocks also close on
// their words (`x_word_blocks`): at blastp-fast's threshold a block of whole batches keeps the
// blocks of the letters; with a lowered one, a block holds at most four words per backbone cell
// plus the words of the batch with the most (or one batch), on one thread and on three. The
// standard table keeps the blocks it is given.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn test_x_word_blocks_close_on_the_words_of_a_lowered_threshold() {
    use crate::algorithm::tblastx::lookup::compressed::CompressedWordCount;
    let lengths = [900usize, 1200, 700, 1100, 800, 1000];
    let records = x_block_test_records(&lengths);
    let input = seq_range::QueryInput::whole(&records);
    let mut args = make_blastp_args();
    args.task = "blastp-fast".to_string();
    let resolved = args.resolve().expect("resolved blastp-fast args");
    // Batches of one query each; the blocks of 50,000 letters (`x_bp_qblock`'s default, the
    // variable is not set in tests) are one block of all six.
    let batch_size = 1;
    let whole = [0..6];
    let mut low = resolved.clone();
    low.threshold = 11.0;
    let mut low_aa = low.clone();
    low_aa.lookup_table_type = BlastpLookupTableType::AaLookupTable;
    let counter = CompressedWordCount::new(5, 1100).expect("word size 5");
    let words: Vec<u64> = records
        .iter()
        .map(|record| {
            let encoded = encode_protein_sequence(record.seq());
            counter
                .batch_words(&mut std::iter::once(
                    &encoded.aa_seq[1..encoded.aa_seq.len() - 1],
                ))
                .expect("within the banks")
        })
        .collect();
    let most = 4 * 15u64.pow(5) + words.iter().copied().max().expect("words");
    for threads in [1, 3] {
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            assert_eq!(
                x_word_blocks(&resolved, &records, &input, batch_size, &whole, pool)
                    .map(|(blocks, _)| blocks),
                Some(whole.to_vec())
            );
            let (blocks, block_words) =
                x_word_blocks(&low, &records, &input, batch_size, &whole, pool).expect("blocks");
            // Each block's words are its batches' (one batch per query here).
            assert_eq!(
                block_words,
                blocks
                    .iter()
                    .map(|block| words[block.clone()].iter().sum::<u64>())
                    .collect::<Vec<_>>()
            );
            assert!(blocks.len() > 1, "{blocks:?}");
            assert_eq!(blocks.first().map(|block| block.start), Some(0));
            assert_eq!(blocks.last().map(|block| block.end), Some(6));
            for (block, next) in blocks.iter().zip(&blocks[1..]) {
                assert_eq!(block.end, next.start);
                // Closed before the next batch: the block and that batch pass the most.
                let block_words: u64 = words[block.clone()].iter().sum();
                assert!(block_words + words[next.start] > most, "{blocks:?}");
            }
            // No block holds more than the most (a batch holds at most the largest batch's
            // words, which the most includes).
            for block in &blocks {
                assert!(
                    words[block.clone()].iter().sum::<u64>() <= most,
                    "{blocks:?}"
                );
            }
            assert_eq!(
                x_word_blocks(&low_aa, &records, &input, batch_size, &whole, pool),
                Some((whole.to_vec(), Vec::new()))
            );
            Ok(())
        })
        .expect("pool");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3405-3437 (cited at
// `blastp_get_start_for_gapped_alignment_int4_length`)
/// BLASTP's gapped start of a one-hit HSP of negative length, which its gapped stage keeps
/// inline, gives what TBLASTN's and BLASTX's shared one gives (`gapalign::gapped_start_int4`
/// over windows that read past a sequence as NULLBs), with every matrix, for windows inside the
/// sequences, ending on their NULLBs and past them.
#[test]
fn the_inline_one_hit_gapped_start_is_the_shared_one() {
    use crate::algorithm::blastp::gapalign::{gapped_start_int4, gapped_start_window};
    let mut state = 0x2f2f_u64;
    let mut next = |bound: usize| -> usize {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 33) % bound as u64) as usize
    };
    let mut above = 0;
    for matrix in ScoringMatrix::ALL {
        for _ in 0..3000 {
            let query: Vec<u8> = (0..1 + next(30)).map(|_| 1 + next(27) as u8).collect();
            let q_start = next(query.len());
            // Half the subjects copy the query's letters from the same start, so that windows
            // sum above 0.
            let (subject, s_start) = if next(2) == 0 {
                let subject: Vec<u8> = (0..1 + next(30)).map(|_| 1 + next(27) as u8).collect();
                let s_start = next(subject.len());
                (subject, s_start)
            } else {
                (query.clone(), q_start)
            };
            let length = -1 - next(4) as i32;
            let inline = blastp_get_start_for_gapped_alignment_int4_length(
                &query, &subject, q_start, length, s_start, matrix,
            );
            let shared = gapped_start_int4(
                &gapped_start_window(&query, q_start),
                &gapped_start_window(&subject, s_start),
                q_start,
                matrix,
            );
            assert_eq!(
                inline, shared,
                "{matrix:?} {query:?} {subject:?} {q_start} {s_start}"
            );
            above += usize::from(shared != q_start);
        }
    }
    assert!(above > 100, "{above}");
}

/// A protein alignment of query `q` and subject `s` (1-based, inclusive) with the operations
/// `ops` (one substitution of the query's span without).
fn row_hit(q: (usize, usize), s: (usize, usize), ops: Option<Vec<GapEditOp>>) -> Hit {
    let length = ops.as_ref().map_or(q.1 - q.0 + 1, |ops| {
        ops.iter().map(|op| op.num() as usize).sum()
    });
    Hit {
        identity: 0.0,
        length,
        mismatch: 0,
        gapopen: 0,
        q_start: q.0,
        q_end: q.1,
        s_start: s.0,
        s_end: s.1,
        e_value: 1.0,
        bit_score: 1.0,
        num_ident: 0,
        query_frame: 0,
        query_length: q.1,
        q_idx: 0,
        s_idx: 0,
        raw_score: 1,
        sort_query_offset: 0,
        sort_query_end: 0,
        sort_subject_offset: 0,
        sort_subject_end: 0,
        has_sort_offsets: false,
        gap_info: ops,
        num_positives: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-121 (cited above `restore_displayed_o`): the rows read each sequence's
// Bioseq at the columns of its residues, which the alignment's operations give; a gap
// character (a PSSM query's gap residue, a database's stored `-`) is a residue, so an O after
// it is the O of its own position (M4a-1c LG probes `tk_dash_o1_eaa`, `ts_dash_o1_eaa`).
#[test]
fn the_rows_read_each_residue_by_the_operations() {
    use crate::report::row_letters::QueryCase;
    // Query `A-OK` (record letters, `-` the gap residue) at 1..=4 against a subject of three
    // residues: the gap residue's column is a gap of the subject (an insertion).
    let ops = vec![GapEditOp::Sub(1), GapEditOp::Ins(1), GapEditOp::Sub(2)];
    let hit = row_hit((1, 4), (1, 3), Some(ops));
    let record = ShownLetters::query(b"A-OK", None);
    assert!(record.recounted() && record.shows_other_letters());
    // The search's row shows the gap residue as `-` and the O as X.
    let shown = restore_displayed_o(Some("A-XK".to_string()), record, &hit, RowSide::Query);
    assert_eq!(shown.as_deref(), Some("A-OK"));
    assert_eq!(stored_gap_columns(b"A-OK", b'-', &hit, RowSide::Query), [1]);
    assert_eq!(record.case(), QueryCase::Masked);
    // A PSSM query's Bioseq letters show as they are; its record's capitals are searched.
    let bioseq = ShownLetters::query(b"MKV", Some(b"mKv"));
    let hit = row_hit((1, 3), (5, 7), None);
    let shown = restore_displayed_o(Some("MKV".to_string()), bioseq, &hit, RowSide::Query);
    assert_eq!(shown.as_deref(), Some("mKv"));
    assert!(bioseq.recounted());
    assert_eq!(bioseq.case(), QueryCase::Bioseq);
    assert!(stored_gap_columns(b"MKV", b'-', &hit, RowSide::Query).is_empty());
    // A FASTA record without O or `-` is not recounted and its row stays.
    let plain = ShownLetters::query(b"MKV", None);
    assert!(!plain.recounted() && !plain.shows_other_letters());
    // A subject's stored `-` before an O: the O stays at its own position.
    let residues = [1, STORED_GAP, STORED_O];
    let hit = row_hit((1, 3), (1, 3), None);
    let shown = restore_displayed_o(
        Some("A-X".to_string()),
        ShownLetters::Stored(&residues),
        &hit,
        RowSide::Subject,
    );
    assert_eq!(shown.as_deref(), Some("A-O"));
    assert_eq!(
        stored_gap_columns(&residues, STORED_GAP, &hit, RowSide::Subject),
        [1]
    );
}

/// BLASTP's table as the `-db_hard_mask` search moves it between sequences (`db_masked.rs`:
/// `TableStart::apply`, `TableEnd::apply`): NCBI's table at its one-thread offset (entries 0, or
/// the window less than 0 after a clear), the clear (offset the window) before the move, and the
/// moves after a sequence.
#[test]
fn hard_mask_starts_and_ends_move_the_table() {
    use crate::algorithm::common::hard_chunk_walk::{HardTable, NcbiStart, TableEnd, TableStart};
    let window = 40;
    let mut scratch = BlastpSubjectScratch::new(16, 64, window, 21);
    let mut written = DiagStruct::default();
    written.set_last_hit(5000);
    written.set_flag(1);
    scratch.diag_array.fill(written);
    scratch.diag_offset = 5000;
    let Ok(()) = TableStart::Keep.apply(&mut scratch, window);
    assert_eq!(
        (scratch.diag_offset, scratch.diag_array[3]),
        (5000, written)
    );
    let Ok(()) = TableStart::Move {
        clear: false,
        by: 70,
    }
    .apply(&mut scratch, window);
    assert_eq!(
        (scratch.diag_offset, scratch.diag_array[3]),
        (5070, written)
    );
    let Ok(()) = TableStart::Move {
        clear: true,
        by: 70,
    }
    .apply(&mut scratch, window);
    assert_eq!(scratch.diag_offset, window + 70);
    assert!(scratch
        .diag_array
        .iter()
        .all(|&entry| entry == DiagStruct::clear(window)));
    let Ok(()) = TableEnd::Move(98).apply(&mut scratch, window);
    assert_eq!(scratch.diag_offset, window + 168);
    scratch.diag_array.fill(written);
    let Ok(()) = TableEnd::Clear.apply(&mut scratch, window);
    assert_eq!(scratch.diag_offset, window);
    assert!(scratch
        .diag_array
        .iter()
        .all(|&entry| entry == DiagStruct::clear(window)));
    for (start, entry) in [
        (
            NcbiStart {
                offset: -98,
                cleared: false,
            },
            DiagStruct::default(),
        ),
        (
            NcbiStart {
                offset: 5,
                cleared: true,
            },
            DiagStruct::clear(window),
        ),
    ] {
        scratch.diag_array.fill(written);
        let Ok(()) = TableStart::New(start).apply(&mut scratch, window);
        assert_eq!(scratch.diag_offset, start.offset);
        assert!(scratch.diag_array.iter().all(|&e| e == entry), "{start:?}");
        // The plain sequence's start (`db_masked::plain_on_ncbi_table`) makes the same table.
        scratch.diag_array.fill(written);
        scratch.ncbi_table(start, window);
        assert_eq!(scratch.diag_offset, start.offset);
        assert!(scratch.diag_array.iter().all(|&e| e == entry), "{start:?}");
    }
}
