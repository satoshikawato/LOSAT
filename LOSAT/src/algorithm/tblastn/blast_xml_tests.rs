use super::*;

fn masks_of(intervals: &[(usize, usize)], length: usize, offset: usize) -> Vec<(i32, i32)> {
    let mut out = Vec::new();
    query_masks(intervals, length, offset, &mut out);
    out.iter()
        .map(|mask| {
            assert_eq!(mask.frame, 0);
            (mask.from, mask.to)
        })
        .collect()
}

#[test]
fn adjacent_intervals_stay_two_masks_and_overlapping_ones_join() {
    // A lower-case run [96, 121) and a SEG interval [121, 153) of a `-query_loc 21-195`
    // query (fuzz seed 2, case_200159): NCBI keeps both (blast_filter.c:993), and its XML
    // mask loop then leaves the letter at 141 (0-based) unmasked when 140 is the
    // second-last column of a segment.
    assert_eq!(
        masks_of(&[(121, 153), (96, 121)], 175, 20),
        vec![(116, 140), (141, 172)]
    );
    // One shared position joins them (inclusive 120 > 120 is false: still apart).
    assert_eq!(
        masks_of(&[(96, 121), (120, 153)], 175, 20),
        vec![(116, 140), (140, 172)]
    );
    // Two shared positions join them.
    assert_eq!(
        masks_of(&[(96, 122), (120, 153)], 175, 20),
        vec![(116, 172)]
    );
    // The whole searched interval is not a mask.
    assert_eq!(
        masks_of(&[(0, 100), (100, 175)], 175, 20),
        vec![(20, 119), (120, 194)]
    );
    assert_eq!(masks_of(&[(0, 175)], 175, 20), Vec::<(i32, i32)>::new());
}

#[test]
fn report_intervals_keep_the_lowercase_run_apart_from_seg() {
    // Without SEG, the lower-case runs alone, unjoined.
    let query = b"ACDEFghiklMNPQrsT";
    assert_eq!(
        super::super::search_seed::tblastn_report_mask_intervals(query, None, true),
        vec![(5, 10), (14, 16)]
    );
    assert!(
        super::super::search_seed::tblastn_report_mask_intervals(query, None, false).is_empty()
    );
}
