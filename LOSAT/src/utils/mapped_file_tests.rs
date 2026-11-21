use super::*;

#[test]
fn a_mapped_file_has_the_files_bytes_and_aligned_views() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-mapped-file-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("words.bin");
    let words: Vec<i32> = (0..70).map(|i| i * 3 - 7).collect();
    let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_ne_bytes()).collect();
    std::fs::write(&path, &bytes).unwrap();
    let mapped = MappedFile::open(&path).unwrap();
    assert_eq!(mapped.bytes(), &bytes[..]);
    assert_eq!(mapped.len(), 280);
    assert_eq!(mapped.i32s(8, 3).unwrap(), &words[2..5]);
    assert!(mapped.i32s(2, 1).is_none());
    assert!(mapped.i32s(276, 2).is_none());
    let rows = mapped.rows::<28>(0, 2).unwrap();
    assert_eq!(rows[1][0], words[28]);
    assert!(mapped.rows::<28>(0, 3).is_none());
    let empty = dir.join("empty.bin");
    std::fs::write(&empty, b"").unwrap();
    assert_eq!(MappedFile::open(&empty).unwrap().len(), 0);
    assert!(MappedFile::open(&dir.join("missing")).is_err());
    let _ = std::fs::remove_dir_all(&dir);
}
