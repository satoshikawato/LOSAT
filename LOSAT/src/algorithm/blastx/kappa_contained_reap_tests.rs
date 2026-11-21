use super::*;

#[test]
fn final_reap_matches_direct_pinned_c_boundaries() {
    let text = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/unit/blastx_contained_reap_expected_v020.tsv"
    ));
    let mut cases: std::collections::BTreeMap<&str, Vec<(RedoneHsp, i32)>> =
        std::collections::BTreeMap::new();
    for line in text
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
    {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 10);
        let n: Vec<i32> = fields[1..]
            .iter()
            .map(|x| x.parse().expect("C integer field"))
            .collect();
        let entries = cases.entry(fields[0]).or_default();
        assert_eq!(n[0] as usize, entries.len(), "C PRE slot order");
        entries.push((
            RedoneHsp {
                hsp: PreliminaryHsp {
                    context: n[1] as usize,
                    score: n[2],
                    frame: n[3] as i8,
                    q_start: n[4],
                    q_end: n[5],
                    q_gapped_start: 0,
                    s_start: n[6],
                    s_end: n[7],
                    s_gapped_start: 0,
                },
                edit_script: Vec::new(),
                composition_method: 0,
                identity_matrix: None,
            },
            n[8],
        ));
    }
    assert_eq!(cases.len(), 21);
    for (name, entries) in cases {
        let mut expected: Vec<(i32, PreliminaryHsp)> = entries
            .iter()
            .filter(|(_, slot)| *slot >= 0)
            .map(|(h, slot)| (*slot, h.hsp.clone()))
            .collect();
        expected.sort_by_key(|(slot, _)| *slot);
        let input = entries.into_iter().map(|(h, _)| h).collect();
        let actual = reap_contained_hsps(input);
        let actual_hsps: Vec<PreliminaryHsp> = actual.into_iter().map(|h| h.hsp).collect();
        let expected_hsps: Vec<PreliminaryHsp> = expected.into_iter().map(|(_, h)| h).collect();
        assert_eq!(actual_hsps, expected_hsps, "direct C boundary case {name}");
    }
}
