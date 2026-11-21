use super::build_two_stage_lookup;

fn hits(ascending_cells: bool) -> Vec<u32> {
    // "ACGTACGT" repeated: the 4-letter word ACGT occurs at offsets 0, 4, 8, ...
    // of the first context and of the second (offset 21 in the block).
    let unit = [0u8, 1, 2, 3];
    let context: Vec<u8> = unit.iter().cycle().take(20).copied().collect();
    let lookup = build_two_stage_lookup(
        &[&context, &context],
        &[0, 21],
        8,
        4,
        &[Vec::new(), Vec::new()],
        None,
        0,
        40,
        ascending_cells,
    );
    let acgt = 0b00_01_10_11u64;
    let mut found = Vec::new();
    lookup.for_each_hit(acgt, |q_off_1| found.push(q_off_1 - 1));
    found
}

#[test]
fn small_and_standard_cells_list_offsets_in_indexing_order() {
    let ascending = hits(true);
    assert_eq!(ascending, vec![0, 4, 8, 12, 16, 21, 25, 29, 33, 37]);
    let mut newest_first = hits(false);
    assert_eq!(newest_first.first(), Some(&37));
    newest_first.reverse();
    assert_eq!(newest_first, ascending);
}
