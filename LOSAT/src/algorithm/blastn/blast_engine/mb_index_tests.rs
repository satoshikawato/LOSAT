use super::*;

#[test]
fn old_style_names_split_and_parse_as_ncbi() {
    assert_eq!(split_old_style_names(b"a"), [b"a".to_vec()]);
    assert_eq!(
        split_old_style_names(b" a  b "),
        [b"".to_vec(), b"a".to_vec(), b"b".to_vec()]
    );
    assert_eq!(split_old_style_names(b""), [b"".to_vec()]);
    let max = u64::MAX;
    for (name, parsed) in [
        (&b"oidx"[..], (&b"oidx"[..], 0, 99)),
        (b"oidx,,1,3", (b"oidx", 1, 3)),
        (b"oidx,x,1x,+2", (b"oidx", 1, 2)),
        (b"oidx,,-1", (b"oidx", max, 99)),
        (b"oidx,,0,-2", (b"oidx", 0, max - 1)),
        (b"oidx,,", (b"oidx", 0, 99)),
        (b"oidx,,,5", (b"oidx", 0, 5)),
        (b"oidx,", (b"oidx", 0, 99)),
        (b",,7", (b"", 7, 99)),
        (b"o,,4,5,6", (b"o", 4, 5)),
        (b"d\xe9/idx,,1,2", (b"d\xe9/idx", 1, 2)),
    ] {
        let (base, first, last) = parse_old_style_name(name);
        assert_eq!((base.as_slice(), first, last), parsed);
    }
}
