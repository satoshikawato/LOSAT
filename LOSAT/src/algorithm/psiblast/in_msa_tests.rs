//! Tests of `-in_msa`'s port of `CPsiBlastInputClustalW` and of the use of the PSSM engine on it.
//! The expected answers are NCBI BLAST+ 2.17.0's (the probes of line PS item PS6, named after
//! each case).
use super::*;
use crate::utils::matrix::aa_char_to_ncbistdaa;

fn rows(rows: &[&str]) -> Vec<Vec<u8>> {
    rows.iter().map(|row| row.as_bytes().to_vec()).collect()
}

fn codes(letters: &str) -> Vec<u8> {
    letters
        .bytes()
        .map(|letter| aa_char_to_ncbistdaa(letter.to_ascii_uppercase()))
        .collect()
}

fn message(error: &anyhow::Error) -> (i32, String) {
    let native = error
        .downcast_ref::<crate::cli::NativeError>()
        .expect("a NativeError");
    (
        native.exit,
        String::from_utf8(native.message.clone()).expect("a UTF-8 message"),
    )
}

// `x_ExtractQueryFromMsa`: the master is the first row, or the `-msa_master_idx`-th (1-based); one
// past the rows is an options error (`o.midx.v[4]`: "between 1 and 3").
#[test]
fn the_master_is_a_row_without_its_gaps() {
    let aln = rows(&["MK-VL", "mkAVL", "-KAVL"]);
    assert_eq!(extract_query(&aln, None).unwrap(), codes("MKVL"));
    assert_eq!(extract_query(&aln, Some(1)).unwrap(), codes("MKVL"));
    // Lower case is read as upper case.
    assert_eq!(extract_query(&aln, Some(2)).unwrap(), codes("MKAVL"));
    assert_eq!(extract_query(&aln, Some(3)).unwrap(), codes("KAVL"));
    let error = extract_query(&aln, Some(4)).unwrap_err();
    assert_eq!(
        message(&error),
        (
            1,
            "LOSAT options error: Invalid master sequence index, please use a value between 1 and 3\n"
                .to_string()
        )
    );
    // `*` is a letter of the alphabet and has a code.
    assert_eq!(
        extract_query(&rows(&["M*K", "MAK"]), None).unwrap(),
        codes("M*K")
    );
}

// `x_ValidateQueryInMsa`: the first row whose letters start with the master's is swapped with row
// 0, so a master of row 2 whose letters row 0 does not have goes first.
#[test]
fn the_master_is_swapped_to_the_first_row() {
    let mut aln = rows(&["MKVLA", "TKVLA", "MKV-A"]);
    let query = extract_query(&aln, Some(2)).unwrap();
    validate_query(&mut aln, &query).unwrap();
    assert_eq!(aln, rows(&["TKVLA", "MKVLA", "MKV-A"]));
    // A row that goes on after the query matches too (the check is of a prefix), and the first
    // row that matches is the one swapped, not the master's own.
    let mut aln = rows(&["MKVLA", "MKV--", "TKVLA"]);
    let query = extract_query(&aln, Some(2)).unwrap();
    validate_query(&mut aln, &query).unwrap();
    assert_eq!(aln[0], b"MKVLA");
    // Selenocysteine: `U` of a row matches an `X` of the query.
    let mut aln = rows(&["MKVLA", "MKULA"]);
    let query = codes("MKXLA");
    validate_query(&mut aln, &query).unwrap();
    // No row starts with the master's letters: the options error.
    let mut aln = rows(&["MKVLA", "TKVLA"]);
    let error = validate_query(&mut aln, &codes("MKVLAA")).unwrap_err();
    assert_eq!(
        message(&error),
        (
            1,
            "LOSAT options error: No sequence in the multiple sequence alignment provided matches the query sequence\n"
                .to_string()
        )
    );
}

// `x_CopyQueryToMsa` and `x_ExtractAlignmentData`: a lower-case residue of the master is not
// aligned; the flanking gaps and the internal gaps of 10 or more of another row are not.
#[test]
fn the_engine_msa_marks_unaligned_cells() {
    let master = "MKvLAGTEDWQAS-HGFRKLMNPQ";
    let long = "MK----------LAGTEDWQA---";
    let short = "-KVLAGTEDW---AHGFRKLMNPQ";
    assert_eq!((master.len(), long.len(), short.len()), (24, 24, 24));
    let aln = rows(&[master, long, short]);
    let query = extract_query(&aln, None).unwrap();
    let msa = engine_msa(&QueryRows::new(&aln, &query).unwrap());
    assert_eq!(msa.query_length, 23);
    assert_eq!(msa.rows.len(), 3);
    // The master: only the lower-case `v` is not aligned.
    let unaligned: Vec<usize> = (0..23).filter(|&i| !msa.rows[0][i].is_aligned).collect();
    assert_eq!(unaligned, vec![2]);
    // Row 1: a gap run of 10 inside, and a gap run of 3 at the end (flanking: unaligned).
    let unaligned: Vec<usize> = (0..23).filter(|&i| !msa.rows[1][i].is_aligned).collect();
    assert_eq!(unaligned, (2..12).chain([20, 21, 22]).collect::<Vec<_>>());
    // Row 2: a leading gap (unaligned), and an internal gap of 3 (aligned, letter 0).
    assert!(!msa.rows[2][0].is_aligned);
    assert_eq!(msa.rows[2][0].letter, 0);
    assert!(msa.rows[2][10].is_aligned && msa.rows[2][10].letter == 0);
    assert!(msa.rows[2][1].is_aligned);
}

// `x_ExtractQueryForPssm`: the first Bioseq of the master's length and letters in the file's
// order (`e.same13.m3`); a master with `*` has none (IUPACaa has no `*`: `r.ch.q.star`), and a
// nucleotide Bioseq of that length that comes first ends the run (`e.len1`).
#[test]
fn the_query_bioseq_is_the_first_equal_one() {
    let read = |rows: &[(&str, &str)]| {
        let mut text = String::from("CLUSTAL W (1.83) multiple sequence alignment\n\n\n");
        for (id, letters) in rows {
            text.push_str(&format!("{id:<8}{letters}\n"));
        }
        text.push_str(&format!("        {}\n\n", " ".repeat(rows[0].1.len())));
        read_alignment(text.as_bytes()).unwrap()
    };
    let aln = read(&[("a1", "MKVLAGT"), ("b2", "MKVIAGT"), ("a3", "MKVLAGT")]);
    let query = codes("MKVLAGT");
    let found = query_bioseq(&aln.bioseqs, &query).unwrap().unwrap();
    assert_eq!(found.letters, b"MKVLAGT");
    assert_eq!(found.ids, aln.bioseqs[0].ids);
    let star = read(&[("a1", "MK*LAGT"), ("b2", "MKVIAGT")]);
    assert!(query_bioseq(&star.bioseqs, &codes("MK*LAGT"))
        .unwrap()
        .is_none());
    let dna = read(&[("q1", "T"), ("s1", "T")]);
    let error = query_bioseq(&dna.bioseqs, &codes("T")).unwrap_err();
    // NCBI's report (PS6 corpus `e.len1`, standard error after `Error: `) in LOSAT's words
    // (`crate::wording::losat_text`, the one mapping of the freeze): the message.
    let ncbi = format!(
        "NCBI C++ Exception:\n    T0 \"{}/src/util/sequtil/sequtil_convert_imp.cpp\", line 225: Error: (CSeqUtilException::Invalid coding) ncbi::CSeqConvert_imp::Convert() - Unknown conversion.\n",
        crate::ncbi_exception::NCBI_BUILD_ROOT
    );
    assert_eq!(
        error.what(),
        crate::wording::losat_bytes(ncbi.as_bytes()).as_ref()
    );
    assert_eq!(error.what(), b"Unknown conversion.\n");
}

fn scoring(matrix_name: &str, gap_open: i32, gap_extend: i32) -> MsaScoring<'_> {
    MsaScoring {
        matrix_name,
        gap_open,
        gap_extend,
        save_ascii_pssm: false,
    }
}

/// The error of `score_block` for the options and the query's letters.
fn score_block_error(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    query: &str,
) -> (i32, String) {
    let error = score_block(&scoring(matrix_name, gap_open, gap_extend), &codes(query))
        .err()
        .unwrap();
    message(&error)
}

const TK_MASTER: &str = "TMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEV";

// What NCBI's score block says of the options' matrix and gap costs (`o.ok.matrix`, `b.gap.bad`).
#[test]
fn the_score_block_reports_the_matrix_and_the_gap_costs() {
    let query = codes(TK_MASTER);
    assert!(score_block(&scoring("BLOSUM62", 11, 1), &query).is_ok());
    assert!(score_block(&scoring("blosum80", 10, 1), &query).is_ok());
    assert_eq!(
        score_block_error("NOPE", 11, 1, TK_MASTER),
        (3, "LOSAT engine error: Unknown error code -1\n".to_string())
    );
    assert_eq!(
        score_block_error("blosum62", 5, 5, TK_MASTER),
        (
            3,
            "LOSAT engine error: Gap existence and extension values of 5 and 5 not supported for BLOSUM62\n"
                .to_string()
        )
    );
    // IDENTITY has no frequency ratios: no score block for the engine.
    assert!(score_block(&scoring("IDENTITY", 15, 2), &query)
        .unwrap()
        .1
        .is_none());
}

// S1 of the PS6 review: a gap cost of 0 is "not given" to the PSSM engine, which then takes its
// default (11 for the existence, 1 for the extension), and so is an empty matrix name (BLOSUM62).
// NCBI BLAST+ 2.17.0, review probes `z.go0`, `z.ge0`, `z.go0.ge0` (the engine goes on),
// `z.go0.ge5`, `z.go5.ge0`, `z.b45.go0.ge0`, `z.pam30.go0.ge0` (the engine's own message, with the
// values it took), and `PS/logs/fix6/probe_mx` for the matrix name.
#[test]
fn a_gap_cost_of_zero_is_the_engines_default() {
    let query = codes(TK_MASTER);
    for (gap_open, gap_extend) in [(0, 1), (11, 0), (0, 0)] {
        let (matrix, block) = score_block(&scoring("BLOSUM62", gap_open, gap_extend), &query)
            .unwrap_or_else(|error| panic!("{gap_open}/{gap_extend}: {error}"));
        assert_eq!(matrix, ScoringMatrix::Blosum62);
        // The gapped block is the one of 11/1.
        let wanted = score_block(&scoring("BLOSUM62", 11, 1), &query)
            .unwrap()
            .1
            .unwrap();
        let block = block.unwrap();
        assert_eq!(block.gap_psi, wanted.gap_psi, "{gap_open}/{gap_extend}");
        assert_eq!(
            block.gap_std_k.to_bits(),
            wanted.gap_std_k.to_bits(),
            "{gap_open}/{gap_extend}"
        );
    }
    // BLOSUM80 has 11/1 too (`z.b80.go0.ge0`: the engine goes on).
    assert!(score_block(&scoring("BLOSUM80", 0, 0), &query).is_ok());
    for (matrix_name, gap_open, gap_extend, said) in [
        ("BLOSUM62", 0, 5, "11 and 5 not supported for BLOSUM62"),
        ("BLOSUM62", 5, 0, "5 and 1 not supported for BLOSUM62"),
        ("BLOSUM45", 0, 0, "11 and 1 not supported for BLOSUM45"),
        ("BLOSUM45", 0, 2, "11 and 2 not supported for BLOSUM45"),
        ("PAM30", 0, 0, "11 and 1 not supported for PAM30"),
        ("", 5, 5, "5 and 5 not supported for BLOSUM62"),
    ] {
        assert_eq!(
            score_block_error(matrix_name, gap_open, gap_extend, TK_MASTER),
            (
                3,
                format!("LOSAT engine error: Gap existence and extension values of {said}\n")
            ),
            "{matrix_name} {gap_open}/{gap_extend}"
        );
    }
    // An empty matrix name is BLOSUM62 to the engine.
    let (matrix, block) = score_block(&scoring("", 11, 1), &query).unwrap();
    assert_eq!(matrix, ScoringMatrix::Blosum62);
    assert!(block.is_some());
}

// S2 of the PS6 review: with gap costs that the matrix does not have, the engine's error is its
// first message, which for a master whose ungapped parameters cannot be computed is the ungapped
// calculation's (NCBI BLAST+ 2.17.0, review probes `k.X20.badgap`, `k.star20.badgap`,
// `o.allx.badgap`, `z.x20.go0.ge5`, `z.x20.b45.go5.ge5`: exit status 3 and this line).
#[test]
fn the_ungapped_message_comes_before_the_gap_costs_one() {
    let ungapped = (
        3,
        "LOSAT engine error: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options\n"
            .to_string(),
    );
    let all_x = "X".repeat(20);
    assert_eq!(score_block_error("BLOSUM62", 5, 5, &all_x), ungapped);
    assert_eq!(score_block_error("BLOSUM62", 0, 5, &all_x), ungapped);
    assert_eq!(score_block_error("BLOSUM45", 5, 5, &all_x), ungapped);
    assert_eq!(
        score_block_error("BLOSUM62", 5, 5, &"*".repeat(20)),
        ungapped
    );
    // The text has one owner (the search's own report of such a query).
    assert_eq!(
        ungapped.1,
        format!(
            "LOSAT engine error: {}\n",
            crate::report::query_warnings::INVALID_QUERY_MESSAGE
        )
    );
    // With gap costs that the matrix has there is no error here: the engine goes on without the
    // query's block (`ungapped_karlin_valid`), where NCBI crashes and LOSAT refuses the run.
    assert!(score_block(&scoring("BLOSUM62", 11, 1), &codes(&all_x)).is_ok());
    assert!(score_block(&scoring("BLOSUM62", 0, 0), &codes(&all_x)).is_ok());
    // A master with letters whose parameters can be computed keeps the gap costs' message.
    assert_eq!(
        score_block_error("BLOSUM62", 5, 5, "XXXXXXXXXXXXXXXXMK"),
        (
            3,
            "LOSAT engine error: Gap existence and extension values of 5 and 5 not supported for BLOSUM62\n"
                .to_string()
        )
    );
}

// A query of X alone leaves NCBI's ungapped Karlin-Altschul block invalid (`r.ch.q.allX`,
// `x.in.allO`).
#[test]
fn a_query_of_x_alone_has_no_ungapped_karlin_block() {
    let matrix = ScoringMatrix::Blosum62;
    assert!(ungapped_karlin_valid(matrix, &codes("MKVLAGTEDWQASHGFRK")));
    assert!(!ungapped_karlin_valid(matrix, &codes(&"X".repeat(60))));
    assert!(!ungapped_karlin_valid(matrix, &codes(&"O".repeat(60))));
    assert!(ungapped_karlin_valid(matrix, &codes("XXXXXXXXXXXXXXXXMK")));
}

// The PSSM engine on the alignment of `ncbi/x.ina.tk.plain` and `x.ina.tk.ign` (NCBI 2.17.0's
// `-out_ascii_pssm` of the same three-row Clustal file, BLOSUM62, 11/1): the scores of the first
// six and the last four positions, `-ignore_msa_master` being the engine's
// `nsg_compatibility_mode`. The letters are `ARNDCQEGHILKMFPSTWYV`.
const TK_ALIGNMENT: &str = "CLUSTAL W (1.83) multiple sequence alignment\n\n\nq1      TMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEV\ns1      TNIGKIHLVIGPMFAGKSSTIIKIYQKKIKMQKRCVFIKHVKDNRYGCELTETVTHSGQT\ns2      -MTATIHMILGPMFSGKSTTLMKLIKHYKISNAKCMLITHTLDKRYTNGN-QVSSHDGEK\n         .  :  :  ****:***  :  :  :  :  : *:: *     ** .     .* *  \n\n";

const TK_PLAIN: [(usize, [i32; 20]); 10] = [
    (
        0,
        [
            0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -3, -2, 0,
        ],
    ),
    (
        1,
        [
            -1, -1, 3, -1, -2, 0, -2, -2, -1, 0, 1, -1, 5, -1, -3, -1, -1, -2, -2, 0,
        ],
    ),
    (
        2,
        [
            2, -2, -1, -2, -1, -1, -2, -1, -2, 2, 0, -1, 0, -2, -2, 1, 3, -3, -2, 1,
        ],
    ),
    (
        3,
        [
            2, -2, -1, -2, -2, -2, -2, 5, -2, -3, -3, -1, -2, -3, -2, 0, -1, -3, -3, -2,
        ],
    ),
    (
        4,
        [
            -1, 4, 0, -1, -2, 1, 0, -2, -1, -2, -2, 3, -1, -3, -2, 0, 3, -3, -2, -1,
        ],
    ),
    (
        5,
        [
            -1, -3, -3, -3, -1, -3, -3, -4, -4, 4, 1, -3, 1, 0, -3, -2, -1, -3, -1, 3,
        ],
    ),
    (
        56,
        [
            2, -1, 0, 4, -2, 0, 0, -1, -1, -2, -3, -1, -2, -3, -1, 2, 0, -3, -2, -2,
        ],
    ),
    (
        57,
        [
            0, -3, -1, -2, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -3, -3, -4,
        ],
    ),
    (
        58,
        [
            -1, 0, 0, 1, -4, 4, 5, -2, 0, -3, -3, 1, -2, -4, -1, 0, -1, -3, -2, -3,
        ],
    ),
    (
        59,
        [
            0, 0, -1, -2, -2, 0, -1, -2, -2, 0, -1, 2, 0, -2, -2, 0, 3, -3, -2, 2,
        ],
    ),
];

const TK_IGNORE: [(usize, [i32; 20]); 10] = [
    (
        0,
        [
            0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0,
        ],
    ),
    (
        1,
        [
            -1, -1, 4, 0, -2, 0, -1, -2, 0, 0, 0, -1, 4, -1, -3, 0, 0, -3, -2, -1,
        ],
    ),
    (
        2,
        [
            -1, -2, -1, -2, -1, -2, -2, -3, -3, 3, 1, -2, 0, -1, -2, 0, 3, -3, -2, 2,
        ],
    ),
    (
        3,
        [
            3, -2, -1, -2, -1, -1, -2, 4, -2, -3, -3, -1, -2, -3, -2, 1, -1, -3, -3, -1,
        ],
    ),
    (
        4,
        [
            0, 1, 0, -1, -2, 1, 0, -2, -1, -2, -2, 3, -1, -3, -1, 1, 3, -3, -2, -1,
        ],
    ),
    (
        5,
        [
            -1, -3, -4, -4, -1, -3, -4, -4, -4, 5, 2, -3, 1, 0, -3, -3, -1, -3, -2, 3,
        ],
    ),
    (
        56,
        [
            0, -1, 1, 5, -2, 0, 1, -1, -1, -3, -3, -1, -2, -3, -1, 3, 0, -4, -3, -3,
        ],
    ),
    (
        57,
        [
            0, -3, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -3, -3, -4,
        ],
    ),
    (
        58,
        [
            -1, 1, 0, 1, -4, 5, 4, -2, 0, -3, -3, 1, -1, -4, -1, 0, -1, -3, -2, -3,
        ],
    ),
    (
        59,
        [
            0, 1, 0, -1, -2, 1, 0, -2, -1, -2, -2, 3, -1, -3, -1, 1, 3, -3, -2, -1,
        ],
    ),
];

/// `pssm_from_msa` of the Clustal text `alignment` written to a file of a directory named after
/// `name`.
fn pssm_of(
    name: &str,
    alignment: &str,
    master_index: Option<i32>,
    ignore_master: bool,
    scoring: &MsaScoring<'_>,
) -> Result<crate::report::serial::pssm::PssmWithParameters> {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-in-msa-{}-{name}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("m.aln");
    std::fs::write(&path, alignment).unwrap();
    let input = MsaInput {
        path: path.clone(),
        master_index,
        ignore_master,
    };
    let pssm = pssm_from_msa(&input, scoring, &mut std::io::sink());
    let _ = std::fs::remove_dir_all(&dir);
    pssm
}

fn tk_scores(ignore_master: bool) -> (Vec<i32>, usize) {
    let pssm = pssm_of(
        &format!("tk-{ignore_master}"),
        TK_ALIGNMENT,
        None,
        ignore_master,
        &scoring("BLOSUM62", 11, 1),
    )
    .unwrap();
    (
        pssm.pssm.final_data.unwrap().scores,
        pssm.pssm.num_rows as usize,
    )
}

fn check_tk(expected: &[(usize, [i32; 20]); 10], ignore_master: bool) {
    let (scores, rows) = tk_scores(ignore_master);
    assert_eq!(rows, 28);
    assert_eq!(scores.len(), 28 * 60);
    for (column, wanted) in expected {
        let got: Vec<i32> = "ARNDCQEGHILKMFPSTWYV"
            .bytes()
            .map(|letter| scores[column * 28 + usize::from(aa_char_to_ncbistdaa(letter))])
            .collect();
        assert_eq!(
            &got[..],
            &wanted[..],
            "column {column}, ignore {ignore_master}"
        );
    }
}

#[test]
fn the_engine_scores_the_alignment_as_ncbi_does() {
    check_tk(&TK_PLAIN, false);
}

#[test]
fn the_ignored_master_changes_the_scores_as_in_ncbis_nsg_mode() {
    check_tk(&TK_IGNORE, true);
}

// B1 of the PS6 review (DECISIONS 75.1 of line PS): the review's `p.pfx.10of12`. The master of
// `-msa_master_idx 3` (10 letters) is the start of row 1 (12 letters), which is the first row that
// starts with the master's letters and so the one that NCBI's alignment is filled from.
const PREFIX_ALIGNMENT: &str = "CLUSTAL W (1.83) multiple sequence alignment\n\n\np1      TMAGRVELVTGP\np2      TNIGKIHLVIGP\np3      TMAGRVELVT--\n        *::*:::**:  \n\n";

// The first row after the swap can have more letters than the query; the engine's alignment is
// made only of rows whose first row has the query's letters, one a cell.
#[test]
fn the_engine_alignment_takes_no_first_row_longer_than_the_query() {
    let mut aln = rows(&["TMAGRVELVTGP", "TNIGKIHLVIGP", "TMAGRVELVT--"]);
    let query = extract_query(&aln, Some(3)).unwrap();
    assert_eq!(query.len(), 10);
    // Row 1 is found (it starts with the master's letters) and stays in front.
    assert_eq!(validate_query(&mut aln, &query).unwrap(), 0);
    assert_eq!(aln[0], b"TMAGRVELVTGP");
    assert!(QueryRows::new(&aln, &query).is_none());
    // The same letters in a row with gaps, and a row longer by one letter.
    assert!(QueryRows::new(&rows(&["TMAGRV-ELVTG-", "TNIGKI-HLVIGP"]), &query).is_none());
    // No row at all, and a first row with fewer letters (which `validate_query` never leaves).
    assert!(QueryRows::new(&[], &query).is_none());
    assert!(QueryRows::new(&rows(&["TMAGRVELV---"]), &query).is_none());
    // The master's own row, and a row above it with the same letters (lower case, other gaps).
    let own = rows(&["TMAGRVELVT--", "TNIGKIHLVIGP"]);
    let aligned = QueryRows::new(&own, &query).unwrap();
    let msa = engine_msa(&aligned);
    assert_eq!((msa.query_length, msa.rows.len()), (10, 2));
    assert!(msa.rows.iter().all(|cells| cells.len() == 10));
    assert!(QueryRows::new(&rows(&["tmag-rvelvt-", "TNIGKIHLVIGP"]), &query).is_some());
}

// The run of such a master is refused where NCBI fills its alignment, after the score block's
// errors (`w.pfx.badgap`: NCBI's gap costs message), and the refusal says what to change; with
// the master above the longer row the same alignment gives a PSSM.
#[test]
fn a_master_that_starts_a_longer_row_above_it_is_refused() {
    let error = pssm_of(
        "prefix",
        PREFIX_ALIGNMENT,
        Some(3),
        false,
        &scoring("BLOSUM62", 11, 1),
    )
    .err()
    .unwrap();
    assert!(error.downcast_ref::<crate::cli::NativeError>().is_none());
    assert_eq!(
        error.to_string(),
        "-in_msa with -msa_master_idx 3, a master sequence that is the start of the longer sequence of row 1 above it, is not supported by LOSAT's PSIBLAST (NCBI BLAST+ 2.17.0 reads and writes past its alignment buffers there); put the master sequence above row 1 in the alignment file"
    );
    // `-ignore_msa_master` changes nothing of it.
    let ignored = pssm_of(
        "prefix-ign",
        PREFIX_ALIGNMENT,
        Some(3),
        true,
        &scoring("BLOSUM62", 11, 1),
    )
    .err()
    .unwrap();
    assert_eq!(ignored.to_string(), error.to_string());
    // NCBI's errors of the score block come first.
    let gap = pssm_of(
        "prefix-gap",
        PREFIX_ALIGNMENT,
        Some(3),
        false,
        &scoring("BLOSUM62", 5, 5),
    )
    .err()
    .unwrap();
    assert_eq!(
        message(&gap),
        (
            3,
            "LOSAT engine error: Gap existence and extension values of 5 and 5 not supported for BLOSUM62\n"
                .to_string()
        )
    );
    // The other masters of the file are read, and so is this one from above the longer row.
    for master in [None, Some(1), Some(2)] {
        let pssm = pssm_of(
            &format!("prefix-{master:?}"),
            PREFIX_ALIGNMENT,
            master,
            false,
            &scoring("BLOSUM62", 11, 1),
        )
        .unwrap();
        assert_eq!(pssm.pssm.num_columns, 12);
    }
    let moved = PREFIX_ALIGNMENT.replace(
        "p1      TMAGRVELVTGP\np2      TNIGKIHLVIGP\np3      TMAGRVELVT--\n",
        "p3      TMAGRVELVT--\np1      TMAGRVELVTGP\np2      TNIGKIHLVIGP\n",
    );
    assert_ne!(moved, PREFIX_ALIGNMENT);
    let pssm = pssm_of(
        "prefix-moved",
        &moved,
        None,
        false,
        &scoring("BLOSUM62", 11, 1),
    )
    .unwrap();
    assert_eq!(pssm.pssm.num_columns, 10);
}

// N1 of the PS6 review: the two refusals of a query Seq-id name what the user gave, the
// alignment's master sequence, not a PSSM (review probes `i.r1.gim`, `i.r1.gnl.ord`).
#[test]
fn the_id_refusals_name_the_master_sequence() {
    for (id, said) in [
        (
            "gim|12345",
            "the query Seq-id gim|12345 of the alignment's master sequence is a Geninfo import id, which NCBI BLAST+ cannot resolve; this is not supported by LOSAT's PSIBLAST",
        ),
        (
            "gnl|BL_ORD_ID|3",
            "the query Seq-id gnl|BL_ORD_ID|3 of the alignment's master sequence names a database ordinal; this is not supported by LOSAT's PSIBLAST",
        ),
    ] {
        let alignment = TK_ALIGNMENT.replace("q1      TMAG", &format!("{id} TMAG"));
        assert_ne!(alignment, TK_ALIGNMENT);
        let dir = crate::utils::test_dir::temp_dir().join(format!(
            "losat-in-msa-{}-ids-{}",
            crate::utils::test_dir::process_id(),
            id.len()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("m.aln");
        std::fs::write(&path, alignment).unwrap();
        let input = MsaInput {
            path,
            master_index: None,
            ignore_master: false,
        };
        let error = read_in_msa(
            &input,
            &scoring("BLOSUM62", 11, 1),
            crate::algorithm::pssm_search::in_pssm::PssmReader::PSIBLAST,
            &mut std::io::sink(),
        )
        .err()
        .unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(error.to_string(), said);
    }
}
