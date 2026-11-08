//! Tests of psiblast's PSSM options, input, ASCII PSSM and files (`pssm.rs`).

use super::*;
use crate::core::psi::{PsiDiagnostics, PsiPssm};

/// BLOSUM62 11/1's gapped block (`Blast_KarlinBlkGappedCalc`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:268
// ```c
//     {11, 1, (double) INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.602800, 43.636200},
// ```
const BLOSUM62_11_1: KarlinBlk = KarlinBlk {
    lambda: 0.267,
    k: 0.041,
    h: 0.14,
};

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
        out_pssm: None,
        out_ascii_pssm: ascii.map(str::to_string),
        pseudocount: 0,
        in_pssm: None,
        rounds: None,
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
    // `-out_pssm` satisfies both checks too.
    let asn_only = PsiblastPssmOptions {
        out_pssm: Some("p.asn".to_string()),
        ..options(true, true, None)
    };
    assert!(asn_only.check().is_ok());
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-321 (quoted at
// `PsiblastPssmOptions::request`): without `-save_pssm_after_last_round` one round makes no PSSM.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:591-608 (quoted at
// `PsiDiagnosticsRequest`): the frequency ratios always, the rest with `-out_ascii_pssm`.
#[test]
fn only_the_last_rounds_pssm_with_a_pssm_file_is_a_request() {
    assert_eq!(options(false, false, Some("a.txt")).request(), None);
    assert_eq!(options(false, true, Some("a.txt")).request(), None);
    assert_eq!(options(true, false, None).request(), None);
    let request = options(true, true, Some("-")).request().unwrap();
    assert_eq!(
        request,
        PssmRequest {
            asn_path: None,
            ascii_path: Some("-".to_string()),
            save_each: true,
            pseudocount: 0,
        }
    );
    assert_eq!(request.diagnostics(), PsiDiagnosticsRequest::new_ex(true));
    let with_asn = PsiblastPssmOptions {
        out_pssm: Some("p.asn".to_string()),
        ..options(true, false, None)
    };
    let request = with_asn.request().unwrap();
    assert_eq!(request.asn_path.as_deref(), Some("p.asn"));
    assert_eq!(request.ascii_path, None);
    let diagnostics = request.diagnostics();
    assert!(diagnostics.frequency_ratios);
    assert!(
        !diagnostics.information_content && !diagnostics.sigma && !diagnostics.num_matching_seqs
    );
    let no_asn_flag = PsiblastPssmOptions {
        out_pssm: Some("p.asn".to_string()),
        ..options(false, false, None)
    };
    assert_eq!(no_asn_flag.request(), None);
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1223-1227 testRejectUnsupportedMatrix
// The matrix without frequency ratios that LOSAT knows is IDENTITY (NCBI's case names a matrix
// "TEST", which LOSAT's parser does not accept).
#[test]
fn a_matrix_without_frequency_ratios_has_no_score_block() {
    assert!(score_block(ScoringMatrix::Identity, BLOSUM62_11_1).is_err());
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1239-1245 testPsiLowerCaseMatrix
#[test]
fn a_lower_case_matrix_name_has_frequency_ratios() {
    let matrix: ScoringMatrix = "blosum62".parse().unwrap();
    let sbp = score_block(matrix, BLOSUM62_11_1).unwrap();
    assert_eq!(sbp.bit_scale_factor, 2);
    assert_eq!(
        score_block(ScoringMatrix::Blosum45, BLOSUM62_11_1)
            .unwrap()
            .bit_scale_factor,
        3
    );
    // BLOSUM62's ideal Lambda (`kbp_ideal`), as blastp's engine computes it.
    let ideal = crate::stats::compute_blosum62_ideal_karlin_params().unwrap();
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
    ];
    let msa = input_msa(&query, &alignments, &subject, 0.002).unwrap();
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

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:310-323 (quoted at
// `input_msa`): `msa_index` moves on after a subject that is not counted, and an HSP whose index
// passes the rows NCBI allocates is written out of bounds (`_ASSERT` only). LOSAT refuses it.
#[test]
fn an_alignment_past_the_last_row_is_refused() {
    let query = [1, 2, 3];
    let subject = |_: u32| vec![4, 5, 6];
    // Subject 1 is not counted but comes first: subject 0's HSP takes index 2 of 2 rows.
    let alignments = [
        alignment(1, 1.0, 0, 0, 3, None),
        alignment(0, 1e-10, 0, 0, 3, None),
    ];
    let error = input_msa(&query, &alignments, &subject, 0.002).unwrap_err();
    assert!(error.downcast_ref::<crate::cli::NativeError>().is_none());
    assert_eq!(
        error.to_string(),
        "an alignment past the last row of the PSSM's multiple alignment (a subject included by -inclusion_ethresh after one that is not, where NCBI BLAST+'s psiblast writes out of bounds) is not supported by LOSAT's PSIBLAST"
    );
    // A counted subject first keeps every later counted subject in its row.
    let alignments = [
        alignment(0, 1e-10, 0, 0, 3, None),
        alignment(1, 1.0, 0, 0, 3, None),
        alignment(2, 1e-10, 0, 0, 3, None),
    ];
    let msa = input_msa(&query, &alignments, &subject, 0.002).unwrap();
    assert_eq!(msa.rows.len(), 3);
    assert_eq!(row_letters(&msa.rows[2]), vec![Some(4), Some(5), Some(6)]);
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
        )
        .unwrap();
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
        gapped: BLOSUM62_11_1,
        ungapped: BLOSUM62_11_1,
        diagnostics: PsiDiagnostics {
            information_content: Some(vec![0.5, -0.004]),
            weighted_residue_freqs: Some(freqs),
            gapless_column_weights: Some(vec![1.234, 9.0]),
            num_matching_seqs: Some(vec![2, 2]),
            ..PsiDiagnostics::default()
        },
    };
    let query = [ncbistdaa::A, ncbistdaa::X];
    let text = String::from_utf8(print_ascii_pssm(&pssm, &query, ANCILLARY)).unwrap();
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

const TK_1_40: &[u8] = b"MQLILSHHLTMAGRVELVTGPMFAGKSTYLKNIYQQENGG";

fn pssm_request(pseudocount: i32) -> PssmRequest {
    PssmRequest {
        asn_path: None,
        ascii_path: Some("-".to_string()),
        save_each: false,
        pseudocount,
    }
}

const ANCILLARY: PssmAncillary = PssmAncillary {
    ungapped: KarlinBlk {
        lambda: 0.31998,
        k: 0.13697,
        h: 0.4,
    },
    gapped: BLOSUM62_11_1,
    psi_ungapped: KarlinBlk {
        lambda: 0.31998,
        k: 0.13697,
        h: 0.4,
    },
    psi_gapped: BLOSUM62_11_1,
};

const SCORING: PssmScoring<'static> = PssmScoring {
    matrix: ScoringMatrix::Blosum62,
    matrix_name: "BLOSUM62",
    gap_open: 11,
    gap_extend: 1,
};

/// `query_pssm` of a query of `letters` named `Query_1`.
fn run_query_pssm(
    request: &PssmRequest,
    inclusion_ethresh: f64,
    letters: &[u8],
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Vec<u8>,
) -> Result<QueryPssm> {
    query_pssm(
        request,
        SCORING,
        BLOSUM62_11_1,
        inclusion_ethresh,
        PssmQueryRecord {
            letters,
            id: "Query_1",
            title: b"a title",
        },
        alignments,
        subject,
        ANCILLARY,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:894-903 and
// psi_pssm_input.cpp:166-171 (quoted at `query_pssm`); NCBI 2.17.0 `db.ps.pssm.o.0`: an O of
// the query record is X in the search's buffer, which is the PSSM's query, so the ASCII PSSM
// prints X and the X branch of the gapless weight (`    0.00`) at that position.
#[test]
fn an_o_in_the_query_prints_as_x() {
    // TK 1-40 with O at 30 (as `tk_o.faa`); the subject is TK 1-40 with three substitutions
    // (92.5% identity: not purged as near-identical), so column 30 has two sequences.
    let mut record = TK_1_40.to_vec();
    record[29] = b'O';
    let mut subject_letters = TK_1_40.to_vec();
    for (position, letter) in [(4, b'A'), (14, b'I'), (24, b'S')] {
        subject_letters[position] = letter;
    }
    let subject = |_: u32| ncbistdaa_letters(&subject_letters);
    let alignments = [alignment(0, 1e-20, 0, 0, 40, None)];
    let query = crate::algorithm::blastp::encoding::encode_protein_sequence(&record);
    let query = &query.aa_seq[1..=query.aa_len];
    assert_eq!(query[29], ncbistdaa::X);
    let msa = input_msa(query, &alignments, &subject, 0.002).unwrap();
    let sbp = score_block(ScoringMatrix::Blosum62, BLOSUM62_11_1).unwrap();
    let pssm =
        psi::create_pssm_with_diagnostics(&msa, 0, &sbp, PsiDiagnosticsRequest::new_ex(true))
            .unwrap();
    assert!(pssm.diagnostics.num_matching_seqs.as_ref().unwrap()[29] > 1);
    let text = String::from_utf8(
        run_query_pssm(&pssm_request(0), 0.002, &record, &alignments, &subject)
            .unwrap()
            .ascii
            .unwrap(),
    )
    .unwrap();
    assert_eq!(text.as_bytes(), print_ascii_pssm(&pssm, query, ANCILLARY));
    let line = text
        .lines()
        .find(|line| line.starts_with("   30 "))
        .unwrap();
    assert!(line.starts_with("   30 X  "), "{line}");
    // `"  " info " "` then the X branch's `"    0.00"` (the other branch has no padding).
    assert!(line.ends_with("     0.00"), "{line}");
    // A column with an O nowhere keeps its letter and its weight.
    let line = text
        .lines()
        .find(|line| line.starts_with("   31 "))
        .unwrap();
    assert!(line.starts_with("   31 K  "), "{line}");
    assert!(!line.ends_with("     0.00"), "{line}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1577-1600 and
// c++/src/algo/blast/api/pssm_engine.cpp:131-148 (quoted at `query_pssm`); NCBI 2.17.0
// `db.ps.pssm.pcneg.7` and `db.ps.pssm.inc0.6`: the PSSM engine's options checks, the pseudocount
// first, then the inclusion threshold.
#[test]
fn the_pssm_engine_checks_the_pseudocount_then_the_inclusion_threshold() {
    let subject = |_: u32| ncbistdaa_letters(TK_1_40);
    let alignments = [alignment(0, 1e-20, 0, 0, 40, None)];
    let run = |pseudocount: i32, inclusion_ethresh: f64| {
        run_query_pssm(
            &pssm_request(pseudocount),
            inclusion_ethresh,
            TK_1_40,
            &alignments,
            &subject,
        )
    };
    let pseudocount = (
        1,
        "BLAST options error: IPssmInputData returns invalid PSIBlastOptions: Pseudo count must be greater than or equal to 0\n"
            .to_string(),
    );
    let inclusion = (
        1,
        "BLAST options error: IPssmInputData returns invalid PSIBlastOptions: Inclusion threshold must be greater than 0\n"
            .to_string(),
    );
    assert_eq!(native(&run(-1, 0.002).unwrap_err()), pseudocount);
    assert_eq!(native(&run(-1, 0.0).unwrap_err()), pseudocount);
    assert_eq!(native(&run(0, 0.0).unwrap_err()), inclusion);
    assert_eq!(native(&run(5, -1e-5).unwrap_err()), inclusion);
    assert!(run(0, 1e-300).is_ok());
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
        asn_path: None,
        ascii_path: Some(path.to_str().unwrap().to_string()),
        save_each,
        pseudocount: 0,
    }
}

fn asn_request(path: &std::path::Path, save_each: bool) -> PssmRequest {
    PssmRequest {
        asn_path: Some(path.to_str().unwrap().to_string()),
        ascii_path: None,
        save_each,
        pseudocount: 0,
    }
}

/// The PSSM files of `request` for the queries' PSSMs `asn` and `ascii` (one entry per query),
/// written after the "reports".
fn write_files(request: &PssmRequest, asn: &[&[u8]], ascii: &[&[u8]]) -> Result<()> {
    let mut files = PssmFiles::new(request);
    for query in 0..asn.len().max(ascii.len()) {
        files.add(
            asn.get(query).copied(),
            ascii.get(query).map(|pssm| pssm.to_vec()),
        )?;
    }
    files.commit()
}

/// The names in `dir`, sorted.
fn names(dir: &std::path::Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    names
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:50-66 (quoted at
// `OutputFileReset`); NCBI 2.17.0 `two.last`, `pre.nohit.keep`, `two.each.last`,
// `err.baddir`.
#[test]
fn each_pssm_replaces_the_file_or_gets_a_numbered_one() {
    let dir = scratch("files");
    let file = dir.join("a.txt");
    let pssms: [&[u8]; 2] = [b"first\n", b"second\n"];
    write_files(&request(&file, false), &[], &pssms).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    // No PSSM: the file stays.
    write_files(&request(&file, false), &[], &[]).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    // Versions: F.1, F.2, and F is not touched.
    write_files(&request(&file, true), &[], &pssms).unwrap();
    assert_eq!(std::fs::read(dir.join("a.txt.1")).unwrap(), b"first\n");
    assert_eq!(std::fs::read(dir.join("a.txt.2")).unwrap(), b"second\n");
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    // An ASCII file that cannot be opened is not reported.
    let missing = dir.join("nodir").join("a.txt");
    write_files(&request(&missing, false), &[], &pssms).unwrap();
    assert!(!missing.exists());
    // The same name for both files: each query's ASN.1 file, then its ASCII one.
    let both = PssmRequest {
        asn_path: Some(file.to_str().unwrap().to_string()),
        ..request(&file, false)
    };
    write_files(&both, &[b"asn 1\n", b"asn 2\n"], &pssms).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"second\n");
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:319-321 (quoted at
// `run_resolved_in_pool`'s end): a run that fails before `SavePssmToFile` leaves every PSSM file
// as it was (NCBI 2.17.0, `-out /dev/full`: the earlier file stays, none is made).
#[test]
fn staged_asn_files_change_no_target_before_the_commit() {
    let dir = scratch("staged");
    let file = dir.join("p.asn");
    // Dropped (a run that fails before its reports): no file is made, none is left.
    let mut files = PssmFiles::new(&asn_request(&file, false));
    files.add(Some(b"pssm 1\n"), None).unwrap();
    files.add(Some(b"pssm 2\n"), None).unwrap();
    // The second PSSM's staging file replaces the first one's.
    assert_eq!(names(&dir).len(), 1);
    assert!(!file.exists());
    drop(files);
    assert!(names(&dir).is_empty());
    // An existing file stays until the commit, which replaces it.
    std::fs::write(&file, b"previous file\n").unwrap();
    let mut files = PssmFiles::new(&asn_request(&file, false));
    files.add(Some(b"pssm 1\n"), None).unwrap();
    files.add(Some(b"pssm 2\n"), None).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"previous file\n");
    drop(files);
    assert_eq!(std::fs::read(&file).unwrap(), b"previous file\n");
    assert_eq!(names(&dir), ["p.asn"]);
    let mut files = PssmFiles::new(&asn_request(&file, false));
    files.add(Some(b"pssm 1\n"), None).unwrap();
    files.add(Some(b"pssm 2\n"), None).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"previous file\n");
    files.commit().unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"pssm 2\n");
    assert_eq!(names(&dir), ["p.asn"]);
    // Numbered files: each staged, F.1 truncated in place by the commit, F.2 made.
    std::fs::write(dir.join("p.asn.1"), b"previous file\n").unwrap();
    let mut files = PssmFiles::new(&asn_request(&file, true));
    files.add(Some(b"pssm 1\n"), None).unwrap();
    files.add(Some(b"pssm 2\n"), None).unwrap();
    assert_eq!(names(&dir).len(), 4);
    assert_eq!(
        std::fs::read(dir.join("p.asn.1")).unwrap(),
        b"previous file\n"
    );
    assert!(!dir.join("p.asn.2").exists());
    files.commit().unwrap();
    assert_eq!(names(&dir), ["p.asn", "p.asn.1", "p.asn.2"]);
    assert_eq!(std::fs::read(dir.join("p.asn.1")).unwrap(), b"pssm 1\n");
    assert_eq!(std::fs::read(dir.join("p.asn.2")).unwrap(), b"pssm 2\n");
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:50-66 (quoted at
// `OutputFileReset`): version 0 removes the file and makes it again, a numbered file is
// truncated (`ofstream`), which keeps its inode and mode.
#[cfg(unix)]
#[test]
fn the_commit_makes_or_truncates_the_file_as_ncbi() {
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    let dir = scratch("modes");
    let fresh = dir.join("fresh");
    std::fs::write(&fresh, b"").unwrap();
    let fresh_mode = std::fs::metadata(&fresh).unwrap().mode();
    std::fs::remove_file(&fresh).unwrap();
    let file = dir.join("p.asn");
    for name in ["p.asn", "p.asn.1"] {
        std::fs::write(dir.join(name), b"previous file\n").unwrap();
        std::fs::set_permissions(dir.join(name), std::fs::Permissions::from_mode(0o600)).unwrap();
    }
    let numbered = std::fs::metadata(dir.join("p.asn.1")).unwrap();
    write_files(&asn_request(&file, false), &[b"pssm\n"], &[]).unwrap();
    assert_eq!(std::fs::metadata(&file).unwrap().mode(), fresh_mode);
    write_files(&asn_request(&file, true), &[b"pssm\n"], &[]).unwrap();
    let truncated = std::fs::metadata(dir.join("p.asn.1")).unwrap();
    assert_eq!(truncated.ino(), numbered.ino());
    assert_eq!(truncated.mode(), numbered.mode());
    assert_eq!(std::fs::read(dir.join("p.asn.1")).unwrap(), b"pssm\n");
    std::fs::remove_dir_all(&dir).unwrap();
}

// NCBI 2.17.0 `err.baddir.p`: an -out_pssm file that cannot be opened stops the run after that
// query's report (`PssmFiles::add`); `-out_pssm /dev/full`: one that cannot be written (the review
// probe of PS3).
#[test]
fn an_asn_file_that_cannot_be_written_is_refused_when_staged() {
    let dir = scratch("refused");
    let refused = |path: &std::path::Path, save_each: bool| -> String {
        let mut files = PssmFiles::new(&asn_request(path, save_each));
        files.add(Some(b"pssm\n"), None).unwrap_err().to_string()
    };
    let cannot = "a -out_pssm file that cannot be written (";
    assert!(refused(&dir.join("nodir").join("p.asn"), false).starts_with(cannot));
    std::fs::create_dir(dir.join("d")).unwrap();
    assert!(refused(&dir.join("d"), false).starts_with(cannot));
    assert_eq!(names(&dir), ["d"]);
    // Earlier queries' files are written by the commit after a refused one.
    let file = dir.join("p.asn");
    let mut files = PssmFiles::new(&asn_request(&file, true));
    files.add(Some(b"pssm 1\n"), None).unwrap();
    std::fs::create_dir(dir.join("p.asn.2")).unwrap();
    assert!(files
        .add(Some(b"pssm 2\n"), None)
        .unwrap_err()
        .to_string()
        .starts_with(cannot));
    files.commit().unwrap();
    assert_eq!(std::fs::read(dir.join("p.asn.1")).unwrap(), b"pssm 1\n");
    assert_eq!(names(&dir), ["d", "p.asn.1", "p.asn.2"]);
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(dir.join("missing"), dir.join("link")).unwrap();
        assert_eq!(
            refused(&dir.join("link"), false),
            format!(
                "an -out_pssm file that is a symbolic link to a missing file ({}) is not supported by LOSAT's PSIBLAST",
                dir.join("link").display()
            )
        );
    }
    #[cfg(target_os = "linux")]
    {
        assert!(refused(std::path::Path::new("/dev/full"), false)
            .starts_with("a -out_pssm file that cannot be written (/dev/full: "));
        write_files(
            &asn_request(std::path::Path::new("/dev/null"), false),
            &[b"pssm\n"],
            &[],
        )
        .unwrap();
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// A file that cannot be opened for writing: version 0 removes it when its directory allows
// (`CDirEntry::RemoveEntry`), a numbered file is opened in place and fails.
#[cfg(unix)]
#[test]
fn a_read_only_file_is_replaced_only_with_version_0() {
    use std::os::unix::fs::PermissionsExt;
    let dir = scratch("readonly");
    let file = dir.join("p.asn");
    for name in ["p.asn", "p.asn.1"] {
        std::fs::write(dir.join(name), b"previous file\n").unwrap();
        std::fs::set_permissions(dir.join(name), std::fs::Permissions::from_mode(0o444)).unwrap();
    }
    // Root opens any file for writing.
    if std::fs::OpenOptions::new().write(true).open(&file).is_err() {
        let mut files = PssmFiles::new(&asn_request(&file, true));
        assert!(files
            .add(Some(b"pssm\n"), None)
            .unwrap_err()
            .to_string()
            .starts_with("a -out_pssm file that cannot be written ("));
        drop(files);
        assert_eq!(
            std::fs::read(dir.join("p.asn.1")).unwrap(),
            b"previous file\n"
        );
    }
    write_files(&asn_request(&file, false), &[b"pssm\n"], &[]).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"pssm\n");
    assert_eq!(names(&dir), ["p.asn", "p.asn.1"]);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The re-review of PS3 (RS1): the staging file's name is short, so a target name that NCBI can
// make (up to NAME_MAX, 255 bytes here) is staged; a longer one fails as NCBI's open does
// (`ENAMETOOLONG` from the target's `stat`).
#[cfg(unix)]
#[test]
fn a_long_target_name_is_staged_and_a_too_long_one_refused() {
    let dir = scratch("longname");
    let long = dir.join("L".repeat(245));
    write_files(&asn_request(&long, false), &[b"pssm\n"], &[]).unwrap();
    assert_eq!(std::fs::read(&long).unwrap(), b"pssm\n");
    let mut files = PssmFiles::new(&asn_request(&dir.join("L".repeat(300)), false));
    assert!(files
        .add(Some(b"pssm\n"), None)
        .unwrap_err()
        .to_string()
        .starts_with("a -out_pssm file that cannot be written ("));
    drop(files);
    assert_eq!(names(&dir), ["L".repeat(245)]);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The re-review of PS3 (RS2): an existing file that opens for writing in a directory that takes no
// new file is written in place after the reports, as NCBI's version 0 (which cannot remove its
// name) and numbered files are (NCBI 2.17.0 `rodir`); nothing changes before the commit.
#[cfg(unix)]
#[test]
fn an_existing_file_in_a_directory_that_takes_no_file_is_written_in_place() {
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    let dir = scratch("rodir");
    let ro = dir.join("ro");
    std::fs::create_dir(&ro).unwrap();
    for name in ["p.asn", "p.asn.1"] {
        std::fs::write(ro.join(name), b"previous file\n").unwrap();
    }
    let inode = |name: &str| std::fs::metadata(ro.join(name)).unwrap().ino();
    let (before, before_1) = (inode("p.asn"), inode("p.asn.1"));
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o555)).unwrap();
    // Root makes files in any directory: then the files are staged as elsewhere.
    let in_place = std::fs::File::create(ro.join("probe")).is_err();
    let file = ro.join("p.asn");
    let mut files = PssmFiles::new(&asn_request(&file, false));
    files.add(Some(b"pssm 1\n"), None).unwrap();
    files.add(Some(b"pssm 2\n"), None).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"previous file\n");
    drop(files);
    assert_eq!(std::fs::read(&file).unwrap(), b"previous file\n");
    write_files(&asn_request(&file, false), &[b"pssm 1\n", b"pssm 2\n"], &[]).unwrap();
    write_files(&asn_request(&file, true), &[b"pssm 1\n"], &[]).unwrap();
    assert_eq!(std::fs::read(&file).unwrap(), b"pssm 2\n");
    assert_eq!(std::fs::read(ro.join("p.asn.1")).unwrap(), b"pssm 1\n");
    assert_eq!(inode("p.asn.1"), before_1);
    if in_place {
        assert_eq!(inode("p.asn"), before);
    }
    std::fs::set_permissions(&ro, std::fs::Permissions::from_mode(0o755)).unwrap();
    let _ = std::fs::remove_file(ro.join("probe"));
    assert_eq!(names(&ro), ["p.asn", "p.asn.1"]);
    std::fs::remove_dir_all(&dir).unwrap();
}

// The re-review of PS3 (RN1): `/dev/null` keeps no PSSM bytes; another device keeps them for the
// commit. RN2: `EPERM` from the test open (an immutable or append-only file) is not removable.
#[cfg(target_os = "linux")]
#[test]
fn dev_null_keeps_no_bytes_and_eperm_is_immutable() {
    let mut files = PssmFiles::new(&asn_request(std::path::Path::new("/dev/null"), false));
    files.add(Some(b"pssm\n"), None).unwrap();
    assert!(matches!(files.pending[0].asn, Some((_, Staged::Null))));
    files.commit().unwrap();
    let mut files = PssmFiles::new(&asn_request(std::path::Path::new("/dev/zero"), false));
    files.add(Some(b"pssm\n"), None).unwrap();
    assert!(matches!(files.pending[0].asn, Some((_, Staged::Device(_)))));
    // Linux's EPERM and EACCES.
    assert!(immutable(&std::io::Error::from_raw_os_error(1)));
    assert!(!immutable(&std::io::Error::from_raw_os_error(13)));
}

// The review of PS3 (B2): NCBI writes a PSSM file that is also a report stream among the reports'
// bytes, in an order of their buffers; refused before the search.
#[test]
fn a_pssm_file_that_is_a_report_stream_is_refused() {
    let dir = scratch("alias");
    let out = dir.join("out.txt");
    std::fs::write(&out, b"").unwrap();
    let check = |request: &PssmRequest, queries: usize| {
        check_pssm_files(request, queries, Some(&out)).map_err(|error| error.to_string())
    };
    assert_eq!(
        check(&asn_request(&out, false), 1),
        Err(format!(
            "an -out_pssm file that is also the -out file ({}) is not supported by LOSAT's PSIBLAST",
            out.display()
        ))
    );
    assert!(check(&request(&out, false), 1)
        .unwrap_err()
        .starts_with("an -out_ascii_pssm file that is also the -out file"));
    // With -save_each_pssm, every numbered name up to the number of queries.
    let base = dir.join("p.asn");
    std::fs::write(dir.join("p.asn.2"), b"").unwrap();
    assert!(check(&asn_request(&base, true), 2).is_ok());
    assert!(check(&asn_request(&base, false), 2).is_ok());
    let numbered = check_pssm_files(&asn_request(&base, true), 2, Some(&dir.join("p.asn.2")));
    assert!(numbered
        .unwrap_err()
        .to_string()
        .starts_with("an -out_pssm file that is also the -out file"));
    assert!(check_pssm_files(&asn_request(&base, true), 1, Some(&dir.join("p.asn.2"))).is_ok());
    // `/dev/null` shows no order: it is no report stream (`report_stream_of`).
    assert!(check(&asn_request(std::path::Path::new("/dev/null"), false), 1).is_ok());
    #[cfg(target_os = "linux")]
    {
        use std::os::unix::fs::MetadataExt;
        let null = std::fs::metadata("/dev/null").unwrap().rdev();
        let is_null = |name: &str| std::fs::metadata(name).is_ok_and(|file| file.rdev() == null);
        if !is_null("/dev/stdout") {
            let stdout = check(&asn_request(std::path::Path::new("/dev/stdout"), false), 1);
            assert!(stdout
                .unwrap_err()
                .starts_with("an -out_pssm file that is also the standard output (/dev/stdout)"));
        }
        if !is_null("/dev/stderr") {
            let stderr = check(&request(std::path::Path::new("/dev/stderr"), false), 1);
            assert!(stderr
                .unwrap_err()
                .starts_with("an -out_ascii_pssm file that is also the standard "));
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
