use super::*;

#[test]
fn a_name_that_is_not_utf8_keeps_its_bytes() {
    let path = path_from_bytes(b"l/d\xe9.nin".to_vec());
    assert_eq!(path_bytes(&path).as_ref(), b"l/d\xe9.nin");
    assert_eq!(
        path_bytes(&path.with_file_name(path_from_bytes(b"d\xe9.pdb".to_vec()))).as_ref(),
        b"l/d\xe9.pdb"
    );
    assert_eq!(os_str_bytes(OsStr::new("a")).as_ref(), b"a");
}
