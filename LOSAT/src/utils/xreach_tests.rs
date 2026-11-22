use super::*;

#[test]
fn index_finds_every_name_of_the_set_then_the_graduated_ones() {
    for (i, name) in STRICT_SET.iter().chain(GRADUATED).enumerate() {
        assert_eq!(index(name), i);
    }
}
