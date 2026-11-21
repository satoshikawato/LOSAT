use super::*;

#[test]
fn residue_probabilities_leave_x_out() {
    // A A C X: X is not counted, A is 2/3 and C 1/3.
    let probabilities = residue_probabilities(&[1, 1, 3, 21]);
    assert_eq!(probabilities[1], 2.0 / 3.0);
    assert_eq!(probabilities[3], 1.0 / 3.0);
    assert_eq!(probabilities[21], 0.0);
    assert_eq!(residue_probabilities(&[21, 21]), [0.0; BLASTAA_SIZE]);
}

// A profile scaled by 2 whose scores read on the query's composition have BLOSUM62's
// Lambda at scale 2 rescales to itself; the X column and the columns past the alphabet
// follow NCBI's rules.
#[test]
fn rescaling_keeps_x_and_marks_the_columns_past_the_alphabet() {
    let mut row = [0i32; BLASTAA_SIZE];
    for (j, score) in row.iter_mut().enumerate() {
        *score = if j % 3 == 0 { 8 } else { -6 };
    }
    row[0] = BLAST_SCORE_MIN;
    row[21] = 99;
    let rows = vec![row; 4];
    let mut probabilities = [0.0; BLASTAA_SIZE];
    for (j, p) in probabilities.iter_mut().enumerate().skip(1).take(20) {
        *p = if j % 2 == 0 { 0.06 } else { 0.04 };
    }
    let mut out = Vec::new();
    assert!(rps_rescale_pssm(100.0, &probabilities, &rows, 26, &mut out).unwrap());
    assert_eq!(out.len(), 4);
    assert_eq!(out[0][0], BLAST_SCORE_MIN);
    assert_eq!(out[0][21], 99);
    assert_eq!((out[0][26], out[0][27]), (BLAST_SCORE_MIN, BLAST_SCORE_MIN));
    // Every scaled score has the sign of its source.
    assert!(out[0][3] > 0 && out[0][1] < 0);
    // No negative expected score: no Lambda, NULL.
    let positive = vec![[5i32; BLASTAA_SIZE]; 2];
    assert!(!rps_rescale_pssm(1.0, &probabilities, &positive, 28, &mut out).unwrap());
}
