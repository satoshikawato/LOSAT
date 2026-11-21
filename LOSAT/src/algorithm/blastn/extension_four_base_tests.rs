use super::build_query_four_base_bytes;

#[test]
fn every_blastna_word_and_starting_phase_matches_c_uint1_expression() {
    for word in 0..65_536u32 {
        let bases = [word >> 12, (word >> 8) & 15, (word >> 4) & 15, word & 15];
        let expected = ((bases[0] << 6) | (bases[1] << 4) | (bases[2] << 2) | bases[3]) as u8;
        for phase in 0..4 {
            let mut query = vec![15; phase];
            query.extend(bases.map(|base| base as u8));
            assert_eq!(build_query_four_base_bytes(&query)[phase], expected);
        }
    }
    for len in 0..8 {
        let query = vec![15; len];
        assert_eq!(
            build_query_four_base_bytes(&query),
            vec![255; len.saturating_sub(3)]
        );
    }
}
