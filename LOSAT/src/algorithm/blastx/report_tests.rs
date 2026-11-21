use super::query_endpoints;
// NCBI reference (598d8ae6): c++/src/objects/seqfeat/Genetic_code_table.cpp:195-209
// Expected values are the verbatim C++ function over all 26 BLASTX
// genetic-code tables and 4,096 ncbi4na triplets, not Rust output.
#[test]
fn display_codons_match_pinned_cpp_all_26_codes() {
    let expected = include_str!("../../../tests/unit/blastx_stage_e_display_codon_expected.tsv");
    let mut last_code = 0;
    let mut code = crate::utils::genetic_code::GeneticCode::from_id(1);
    let mut count = 0;
    for row in expected.lines() {
        let values: Vec<_> = row.split('\t').collect();
        let id: u8 = values[0].parse().unwrap();
        if id != last_code {
            code = crate::utils::genetic_code::GeneticCode::from_id(id);
            last_code = id;
        }
        let masks = [
            values[1].parse().unwrap(),
            values[2].parse().unwrap(),
            values[3].parse().unwrap(),
        ];
        assert_eq!(
            super::display_codon(masks, &code),
            values[4].as_bytes()[0],
            "{row}"
        );
        count += 1;
    }
    assert_eq!(count, 106496);
}
#[test]
fn translated_endpoints_match_independent_ncbi_cpp_boundaries() {
    for row in include_str!("../../../tests/unit/blastx_stage_e_coordinates_expected.tsv").lines() {
        let f: Vec<i32> = row.split('\t').map(|x| x.parse().unwrap()).collect();
        assert_eq!(
            query_endpoints(f[0] as usize, f[1] as i8, f[2], f[3]),
            (f[4], f[5]),
            "{row}"
        );
    }
}
