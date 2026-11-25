use super::*;

// Probes L01-L05, J160-J178 of line BDC.
#[test]
fn the_format_is_read_as_ncbi_s_loop_reads_it() {
    let format = ListFormat::parse(b"%f %p").unwrap();
    assert_eq!(format.offsets, [0, 3]);
    assert_eq!(format.types, b"fp");
    let format = ListFormat::parse(b"%%f %f %%").unwrap();
    assert_eq!(format.spec, b"%f %f %");
    assert_eq!(format.offsets, [3]);
    assert_eq!(format.replace(&[b"db/x".to_vec()]), b"%f db/x %");
    let format = ListFormat::parse(b"%f %").unwrap();
    assert_eq!(format.types, b"f\0");
    assert_eq!(format.replace(&[b"a".to_vec(), b"b".to_vec()]), b"a b");
    // After `%%` the loop steps over the byte that follows the second `%`.
    assert!(ListFormat::parse(b"%%%f").is_err());
    assert!(ListFormat::parse(b"abc").is_err());
    assert!(ListFormat::parse(b"").is_err());
    assert!(ListFormat::parse(b"%%").is_err());
}

#[test]
fn the_extensions_follow_the_version() {
    assert_eq!(file_extensions(false).len(), 20);
    let v5 = file_extensions(true);
    assert_eq!(v5.len(), 25);
    assert!(!v5.contains(&"si"));
    assert!(v5.contains(&"tf-lock"));
}
