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
    (native.exit, native.message.clone())
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
            "BLAST options error: Invalid master sequence index, please use a value between 1 and 3\n"
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
            "BLAST options error: No sequence in the multiple sequence alignment provided matches the query sequence\n"
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
    let msa = engine_msa(&aln, &query);
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
    assert!(error
        .what()
        .contains("(CSeqUtilException::Invalid coding) ncbi::CSeqConvert_imp::Convert() - Unknown conversion."));
}

// What NCBI's score block says of the options' matrix and gap costs (`o.ok.matrix`, `b.gap.bad`).
#[test]
fn the_score_block_reports_the_matrix_and_the_gap_costs() {
    let scoring = |matrix_name, gap_open, gap_extend| MsaScoring {
        matrix_name,
        gap_open,
        gap_extend,
        save_ascii_pssm: false,
    };
    assert!(score_block(&scoring("BLOSUM62", 11, 1)).is_ok());
    assert!(score_block(&scoring("blosum80", 10, 1)).is_ok());
    let error = score_block(&scoring("NOPE", 11, 1)).err().unwrap();
    assert_eq!(
        message(&error),
        (3, "BLAST engine error: Unknown error code -1\n".to_string())
    );
    let error = score_block(&scoring("blosum62", 5, 5)).err().unwrap();
    assert_eq!(
        message(&error),
        (
            3,
            "BLAST engine error: Gap existence and extension values of 5 and 5 not supported for BLOSUM62\n"
                .to_string()
        )
    );
    // IDENTITY has no frequency ratios: no score block for the engine.
    assert!(score_block(&scoring("IDENTITY", 15, 2))
        .unwrap()
        .1
        .is_none());
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

fn tk_scores(ignore_master: bool) -> (Vec<i32>, usize) {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-in-msa-{}-{ignore_master}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("tk.aln");
    std::fs::write(&path, TK_ALIGNMENT).unwrap();
    let input = MsaInput {
        path: path.clone(),
        master_index: None,
        ignore_master,
    };
    let scoring = MsaScoring {
        matrix_name: "BLOSUM62",
        gap_open: 11,
        gap_extend: 1,
        save_ascii_pssm: false,
    };
    let pssm = pssm_from_msa(&input, &scoring).unwrap();
    let _ = std::fs::remove_dir_all(&dir);
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
