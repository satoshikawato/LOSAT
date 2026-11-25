use super::*;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:53-62
fn entry(row: &[i32]) -> Entry {
    Entry {
        id: 0,
        oid: row[3] as usize,
        merit: 0,
        hsp: Hsp {
            hsp: super::super::preliminary::PreliminaryHsp {
                context: 0,
                frame: 1,
                score: row[2],
                q_start: row[0],
                q_end: row[1],
                q_gapped_start: row[0],
                s_start: row[4],
                s_end: row[4] + 1,
                s_gapped_start: row[4],
            },
            num: 1,
            evalue: 0.0,
            bit_score: 0.0,
            identity: 0,
            positive: 0,
            edit_script: Vec::new(),
            composition_method: 0,
            query_left_out: 0,
        },
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c:79-120
#[test]
fn pinned_c_dominance_boundaries() {
    let golden = include_str!("../../../tests/unit/blastx_d_boundary_expected.tsv");
    let mut count = 0;
    for line in golden.lines().filter(|l| l.starts_with("C\t")) {
        let row: Vec<i32> = line
            .split('\t')
            .skip(1)
            .map(|s| s.parse().unwrap())
            .collect();
        assert_eq!(
            dominates_source(&entry(&row[..5]), &entry(&row[5..10])),
            row[10] != 0,
            "{line}"
        );
        count += 1;
    }
    assert_eq!(count, 12);
}
