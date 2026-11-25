use super::*;

#[test]
fn merging_joins_overlapping_and_touching_intervals() {
    let mut dest = vec![(0, 9), (30, 39)];
    merge_mask_info(&mut dest, &[(10, 12), (20, 25), (35, 50)]);
    assert_eq!(dest, vec![(0, 12), (20, 25), (30, 50)]);
    let mut dest = vec![];
    merge_mask_info(&mut dest, &[(1, 2)]);
    assert_eq!(dest, vec![(1, 2)]);
    let mut dest = vec![(1, 2)];
    merge_mask_info(&mut dest, &[]);
    assert_eq!(dest, vec![(1, 2)]);
}

#[test]
fn a_sequence_without_masks_is_dusted_whole() {
    let data = b"ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC";
    let mut masker = SDustMasker::new(20);
    let whole = masker.mask(data, &[]);
    let direct = DustMasker::new(20, 64, 1)
        .mask_range(data, 0, data.len(), &mut NcbiRandom::default())
        .into_iter()
        .map(|interval| (interval.start as u32, (interval.end - 1) as u32))
        .collect::<Vec<_>>();
    assert_eq!(whole, direct);
    assert!(!whole.is_empty());
}
