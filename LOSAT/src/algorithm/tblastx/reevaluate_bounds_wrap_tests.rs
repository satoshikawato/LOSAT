use super::{get_num_identities_and_positives_ungapped, reevaluate_ungapped_hit_ncbi_translated};

#[test]
fn an_offset_near_usize_max_is_refused() {
    let seq = [1u8; 40];
    // `q_off + hsp_len` wraps to 1; 16 letters take the AVX2 loop where it is detected.
    for len in [2, 16] {
        let q_off = usize::MAX - len + 2;
        assert_eq!(
            reevaluate_ungapped_hit_ncbi_translated(&seq, &seq, q_off, 0, len, 0),
            None
        );
        assert_eq!(
            reevaluate_ungapped_hit_ncbi_translated(&seq, &seq, 0, usize::MAX, len, 0),
            None
        );
        assert_eq!(
            get_num_identities_and_positives_ungapped(&seq, &seq, usize::MAX, 0, len),
            None
        );
    }
    assert_eq!(
        reevaluate_ungapped_hit_ncbi_translated(&seq, &seq, 0, 0, 41, 0),
        None
    );
    assert!(get_num_identities_and_positives_ungapped(&seq, &seq, 24, 24, 16).is_some());
}
