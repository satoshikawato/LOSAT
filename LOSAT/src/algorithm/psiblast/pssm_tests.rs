//! Tests of psiblast's PSSM options, input, ASCII PSSM and files (`pssm.rs`).

use super::*;
use crate::core::psi::PsiPssm;

fn native(error: &anyhow::Error) -> (i32, String) {
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .expect("a NativeError");
    (native.exit, native.message.clone())
}

fn options(save_last: bool, save_each: bool, ascii: Option<&str>) -> PsiblastPssmOptions {
    PsiblastPssmOptions {
        save_pssm_after_last_round: save_last,
        save_each_pssm: save_each,
        out_ascii_pssm: ascii.map(str::to_string),
        pseudocount: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1283-1303 (quoted at
// `PsiblastPssmOptions::check`); NCBI 2.17.0 `tk.lastonly` and `tk.eachonly`.
#[test]
fn saving_a_pssm_needs_a_pssm_file() {
    let error = options(true, false, None).check().unwrap_err();
    assert_eq!(
        native(&error),
        (
            1,
            "BLAST query/options error: save_pssm_after_last_round option requires out_pssm or out_ascii_pssm\nPlease refer to the BLAST+ user manual.\n"
                .to_string()
        )
    );
    // The last round's check comes first.
    let error = options(true, true, None).check().unwrap_err();
    assert!(native(&error)
        .1
        .contains("save_pssm_after_last_round option"));
    let error = options(false, true, None).check().unwrap_err();
    assert_eq!(
        native(&error),
        (
            1,
            "BLAST query/options error: save_each_pssm option requires out_pssm or out_ascii_pssm\nPlease refer to the BLAST+ user manual.\n"
                .to_string()
        )
    );
    assert!(options(true, true, Some("a.txt")).check().is_ok());
    assert!(options(false, false, Some("a.txt")).check().is_ok());
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-321 (quoted at
// `PsiblastPssmOptions::request`): without `-save_pssm_after_last_round` one round makes no PSSM.
#[test]
fn only_the_last_rounds_pssm_with_an_ascii_file_is_a_request() {
    assert_eq!(options(false, false, Some("a.txt")).request(), None);
    assert_eq!(options(false, true, Some("a.txt")).request(), None);
    assert_eq!(options(true, false, None).request(), None);
    assert_eq!(
        options(true, true, Some("-")).request(),
        Some(PssmRequest {
            ascii_path: "-".to_string(),
            save_each: true,
            pseudocount: 0,
        })
    );
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1223-1227 testRejectUnsupportedMatrix
// The matrix without frequency ratios that LOSAT knows is IDENTITY (NCBI's case names a matrix
// "TEST", which LOSAT's parser does not accept).
#[test]
fn a_matrix_without_frequency_ratios_has_no_score_block() {
    assert!(score_block(ScoringMatrix::Identity).is_err());
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1239-1245 testPsiLowerCaseMatrix
#[test]
fn a_lower_case_matrix_name_has_frequency_ratios() {
    let matrix: ScoringMatrix = "blosum62".parse().unwrap();
    let sbp = score_block(matrix).unwrap();
    assert_eq!(sbp.bit_scale_factor, 2);
    assert_eq!(
        score_block(ScoringMatrix::Blosum45)
            .unwrap()
            .bit_scale_factor,
        3
    );
    // BLOSUM62's ideal Lambda (`kbp_ideal`), as blastp's engine computes it.
    let ideal = crate::core::blast_stat::compute_blosum62_ideal_karlin_params().unwrap();
    assert_eq!(sbp.ideal_lambda, ideal.lambda);
}

fn alignment(
    subject: u32,
    e_value: f64,
    q_start: usize,
    s_start: usize,
    q_length: usize,
    ops: Option<&[GapEditOp]>,
) -> PssmAlignment<'_> {
    PssmAlignment {
        subject,
        e_value,
        q_start,
        s_start,
        q_length,
        ops,
    }
}

fn row_letters(row: &[PsiMsaCell]) -> Vec<Option<u8>> {
    row.iter()
        .map(|cell| cell.is_aligned.then_some(cell.letter))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:283-323 and 363-408
// (quoted at `input_msa` and `process_denseg`).
#[test]
fn the_msa_has_a_row_per_qualifying_subject_in_seq_align_order() {
    let query = [1, 2, 3, 4, 5, 6];
    let subjects: [&[u8]; 3] = [&[10, 11, 12, 13, 14, 15, 16], &[20, 21, 22], &[7, 8, 9, 10]];
    let subject = |s_idx: u32| subjects[s_idx as usize].to_vec();
    let gapped = [
        GapEditOp::Sub(2),
        GapEditOp::Ins(1),
        GapEditOp::Del(2),
        GapEditOp::Sub(1),
    ];
    let alignments = [
        // Subject 0: two HSPs, the second does not overwrite aligned cells.
        alignment(0, 1e-10, 0, 1, 0, Some(&gapped)),
        alignment(0, 1e-5, 1, 0, 3, None),
        // Subject 1 is not below the threshold, but moves `last_sid` on.
        alignment(1, 1.0, 0, 0, 3, None),
        alignment(2, 0.0, 3, 0, 3, None),
        // Subject 0 again, after another subject: the index moves past the rows NCBI allocates
        // (one per qualifying subject), where NCBI writes out of bounds and LOSAT drops it.
        alignment(0, 1e-3, 0, 0, 1, None),
    ];
    let msa = input_msa(&query, &alignments, &subject, 0.002);
    assert_eq!(msa.rows.len(), 3);
    assert_eq!(row_letters(&msa.rows[0]), query.map(Some).to_vec());
    // HSP 1 at query 0: 11 12, a gap in the subject at query 2, the subject's 13 14 skipped,
    // 15 at query 3; HSP 2 fills query 1..4 only where nothing is aligned (none left).
    assert_eq!(
        row_letters(&msa.rows[1]),
        vec![Some(11), Some(12), Some(0), Some(15), None, None]
    );
    // Subject 2's HSP comes after subject 1's: the next row, query 3..6 = 7 8 9.
    assert_eq!(
        row_letters(&msa.rows[2]),
        vec![None, None, None, Some(7), Some(8), Some(9)]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:529-547
// (`GetLowestEvalue`, quoted at `input_msa`): the e-value counts as at most
// BLAST_EXPECT_VALUE (10), and a NaN as 10.
#[test]
fn an_e_value_from_10_on_counts_as_10() {
    let query = [1, 2];
    let subject = |_: u32| vec![3, 4];
    for (e_value, threshold, rows) in [
        (9.0, 10.0, 2),
        (9.0, 9.0, 1),
        (10.0, 10.0, 1),
        (50.0, 20.0, 2),
        (f64::NAN, 20.0, 2),
        (f64::NAN, 10.0, 1),
        (0.0, 1e-300, 2),
    ] {
        let msa = input_msa(
            &query,
            &[alignment(0, e_value, 0, 0, 2, None)],
            &subject,
            threshold,
        );
        assert_eq!(msa.rows.len(), rows, "{e_value} {threshold}");
    }
}

#[test]
fn fixed_formats_as_the_c_library() {
    for (value, precision, text) in [
        (0.125, 2, "0.12"),
        (0.375, 2, "0.38"),
        (-0.001, 2, "-0.00"),
        // The double nearest 0.13695 is below it (0.136949999...).
        (0.13695, 4, "0.1369"),
        (0.03125, 4, "0.0312"),
        (0.09375, 4, "0.0938"),
        (0.04100, 4, "0.0410"),
        (f64::NAN, 2, "nan"),
        (-f64::NAN, 2, "-nan"),
        (f64::INFINITY, 4, "inf"),
        (f64::NEG_INFINITY, 4, "-inf"),
    ] {
        assert_eq!(fixed(value, precision), text, "{value}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:191-327 (quoted at
// `print_ascii_pssm`).
#[test]
fn the_ascii_pssm_has_ncbis_columns() {
    let mut scores = vec![[BLAST_SCORE_MIN; BLASTAA_SIZE]; 2];
    let mut freqs = vec![[0.0; BLASTAA_SIZE]; 2];
    for (c, &residue) in RESIDUE_ORDER.iter().enumerate() {
        scores[0][residue] = c as i32 - 3;
        scores[1][residue] = -(c as i32);
        freqs[0][residue] = if c == 0 { 0.995 } else { 0.0 };
        freqs[1][residue] = if c == 1 { 0.004 } else { 0.0 };
    }
    // A residue at BLAST_SCORE_MIN (`-I`, no percentage), and -12 makes the width 4.
    scores[1][RESIDUE_ORDER[19]] = BLAST_SCORE_MIN;
    scores[1][RESIDUE_ORDER[18]] = -12;
    let pssm = PsiPssm {
        scores,
        information_content: vec![0.5, -0.004],
        weighted_residue_freqs: freqs,
        gapless_column_weights: vec![1.234, 9.0],
        num_matching_seqs: vec![2, 2],
    };
    let query = [ncbistdaa::A, ncbistdaa::X];
    let text = String::from_utf8(print_ascii_pssm(
        &pssm,
        &query,
        PssmAncillary {
            ungapped: (0.13697, 0.31998),
            gapped: (0.041, 0.267),
        },
    ))
    .unwrap();
    let mut expected = String::from("\nLast position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n         ");
    for letter in "ARNDCQEGHILKMFPSTWYV".chars() {
        expected.push_str(&format!("{letter:>4}"));
    }
    for letter in "ARNDCQEGHILKMFPSTWYV".chars() {
        expected.push_str(&format!("   {letter}"));
    }
    expected.push_str("\n    1 A  ");
    for c in 0..20 {
        expected.push_str(&format!("{:>4}", c - 3));
    }
    expected.push(' ');
    expected.push_str(" 100");
    for _ in 1..20 {
        expected.push_str("   0");
    }
    expected.push_str("  0.50 1.23");
    expected.push_str("\n    2 X  ");
    for c in 0..18 {
        expected.push_str(&format!("{:>4}", -c));
    }
    expected.push_str(" -12-I  ");
    for _ in 0..19 {
        expected.push_str("   0");
    }
    expected.push_str("  -0.00     0.00");
    expected.push_str("\n\n                      K         Lambda\n");
    expected.push_str("Standard Ungapped    0.1370     0.3200\n");
    expected.push_str("Standard Gapped      0.0410     0.2670\n");
    expected.push_str("PSI Ungapped         0.1370     0.3200\n");
    expected.push_str("PSI Gapped           0.0410     0.2670\n");
    assert_eq!(text, expected);
}

/// A directory of its own for a test's files.
fn scratch(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("losat-pssm-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn request(path: &std::path::Path, save_each: bool) -> PssmRequest {
    PssmRequest {
        ascii_path: path.to_str().unwrap().to_string(),
        save_each,
        pseudocount: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:50-66 (quoted at
// `write_ascii_pssm_files`); NCBI 2.17.0 `two.last`, `pre.nohit.keep`, `two.each.last`,
// `err.baddir`.
#[test]
fn each_pssm_replaces_the_file_or_gets_a_numbered_one() {
    let dir = scratch("files");
    let file = dir.join("a.txt");
    let pssms: [&[u8]; 2] = [b"first\n", b"second\n"];
    write_ascii_pssm_files(&request(&file, false), pssms.into_iter());
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    // No PSSM: the file stays.
    write_ascii_pssm_files(&request(&file, false), std::iter::empty());
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    // Versions: F.1, F.2, and F is not touched.
    write_ascii_pssm_files(&request(&file, true), pssms.into_iter());
    assert_eq!(std::fs::read(dir.join("a.txt.1")).unwrap(), b"first\n");
    assert_eq!(std::fs::read(dir.join("a.txt.2")).unwrap(), b"second\n");
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    // A file that cannot be opened is not reported.
    let missing = dir.join("nodir").join("a.txt");
    write_ascii_pssm_files(&request(&missing, false), pssms.into_iter());
    assert!(!missing.exists());
    std::fs::remove_dir_all(&dir).unwrap();
}
