use super::*;

#[test]
fn volume_numbers_have_at_least_two_digits() {
    assert_eq!(num_digits(1), 2);
    assert_eq!(num_digits(99), 2);
    assert_eq!(num_digits(100), 3);
    assert_eq!(num_digits(1234), 4);
}

#[test]
fn dir_entry_name_is_the_part_after_the_last_separator() {
    assert_eq!(dir_entry_name(b"in/pmv.00"), b"pmv.00");
    assert_eq!(dir_entry_name(b"pmv.00"), b"pmv.00");
    assert_eq!(dir_entry_name(b"a/b/"), b"");
    assert_eq!(dir_entry_name(b""), b"");
}

#[test]
fn filter_keys_are_ncbis() {
    assert_eq!(FilterType::None.key(), "");
    assert_eq!(FilterType::GiList.key(), "GILIST");
    assert_eq!(FilterType::SeqIdList.key(), "SEQIDLIST");
    assert_eq!(FilterType::TaxIdList.key(), "TAXIDLIST");
}
