use super::*;

#[test]
fn old_style_names_split_and_parse_as_ncbi() {
    assert_eq!(split_old_style_names("a"), ["a"]);
    assert_eq!(split_old_style_names(" a  b "), ["", "a", "b"]);
    assert_eq!(split_old_style_names(""), [""]);
    let max = u64::MAX;
    for (name, parsed) in [
        ("oidx", ("oidx", 0, 99)),
        ("oidx,,1,3", ("oidx", 1, 3)),
        ("oidx,x,1x,+2", ("oidx", 1, 2)),
        ("oidx,,-1", ("oidx", max, 99)),
        ("oidx,,0,-2", ("oidx", 0, max - 1)),
        ("oidx,,", ("oidx", 0, 99)),
        ("oidx,,,5", ("oidx", 0, 5)),
        ("oidx,", ("oidx", 0, 99)),
        (",,7", ("", 7, 99)),
        ("o,,4,5,6", ("o", 4, 5)),
    ] {
        let (base, first, last) = parse_old_style_name(name);
        assert_eq!((base.as_str(), first, last), parsed, "{name}");
    }
}
