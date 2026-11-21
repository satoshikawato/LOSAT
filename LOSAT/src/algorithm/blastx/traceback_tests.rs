use super::*;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3248-3330
#[test]
fn offsets_match_pinned_ncbi_boundary_oracle() {
    let expected = include_str!("../../../../docs/evidence/losatx_stage_d/offsets_expected.tsv");
    let mut count = 0;
    for line in expected.lines() {
        let f: Vec<i32> = line.split('\t').map(|v| v.parse().unwrap()).collect();
        let mode = f[0];
        let n = f[1] as usize;
        let q = vec![1u8; 100];
        let mut s = vec![if mode == 0 { 1 } else { 21 }; 100];
        let end = 7 + n + if mode == 3 { 12 } else { 0 };
        if mode == 2 {
            s[7 + n / 2..7 + n].fill(1);
        }
        if mode == 3 && n > 11 {
            s[end - 11..end].fill(1);
        }
        let h = PreliminaryHsp {
            context: 0,
            frame: 1,
            score: 0,
            q_start: 3,
            q_end: 3 + n as i32,
            q_gapped_start: 0,
            s_start: 7,
            s_end: end as i32,
            s_gapped_start: 0,
        };
        let result = get_offsets_for_gapped_alignment(&q, &s, &h, ScoringMatrix::Blosum62);
        assert_eq!(
            result,
            if f[2] != 0 { Some((f[3], f[4])) } else { None },
            "mode {mode} length {n}"
        );
        count += 1;
    }
    assert_eq!(count, 204);
}
