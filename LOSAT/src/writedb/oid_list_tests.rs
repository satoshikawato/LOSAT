use super::*;

/// The bits of makeblastdb 2.17.0's `db.pxm` of 10 OIDs whose models are 0, 2, 4 and 9
/// (probe `p_ps1`: `0000000957bf`), written to a file of its own.
#[test]
fn the_mask_file_is_makeblastdb_s() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-oidlist-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let base = dir.join("db").to_string_lossy().into_owned();
    let mut list = OidList::exclude_model(base.as_bytes(), true, 0, 0);
    for oid in [0, 2, 4, 9] {
        list.add_oid(oid);
    }
    list.close(10).unwrap();
    let written = std::fs::read(crate::utils::os_bytes::path_of_bytes(list.filename())).unwrap();
    assert_eq!(written, [0, 0, 0, 9, 0x57, 0xBF]);
    let mut empty = OidList::exclude_model(dir.join("e").to_string_lossy().as_bytes(), false, 0, 0);
    empty.close(0).unwrap();
    assert!(!crate::utils::os_bytes::path_of_bytes(empty.filename()).exists());
    std::fs::remove_dir_all(&dir).unwrap();
}
