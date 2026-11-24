//! Tests of the OID list file reader (`oid_list.rs`): NCBI's errors and the bits it copies for a
//! volume (probes G and H of line RG D1b).

use std::path::PathBuf;

use super::*;

/// A directory that is removed when the test ends.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = crate::utils::test_dir::temp_dir().join(format!(
            "losat-seqdb-oidlist-{tag}-{}-{}",
            crate::utils::test_dir::process_id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    /// Writes `bytes` to `file` and returns its path.
    fn write(&self, file: &str, bytes: &[u8]) -> String {
        let path = self.0.join(file);
        std::fs::write(&path, bytes).unwrap();
        path.to_str().unwrap().to_string()
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A file of the header `last_oid` and the bytes of bits `bits`.
fn list_file(last_oid: u32, bits: &[u8]) -> Vec<u8> {
    [last_oid.to_be_bytes().as_slice(), bits].concat()
}

/// The OIDs set in `words`.
fn oids(words: &[u64]) -> Vec<u32> {
    (0..words.len() as u32 * 64)
        .filter(|&oid| words[(oid / 64) as usize] & (1 << (oid % 64)) != 0)
        .collect()
}

#[test]
fn a_path_that_is_not_a_file_is_ncbis_memory_map_error() {
    let dir = TempDir::new("missing");
    let missing = dir.0.join("absent.msk").to_str().unwrap().to_string();
    let error = read(missing.as_bytes()).unwrap_err();
    assert!(matches!(&error, OidListError::NotAFile(path) if *path == missing.as_bytes()));
    let error = error.to_error();
    let text = format!("To be memory mapped the file must exist: '{missing}'");
    assert!(
        matches!(&error, SeqDbError::MemoryMap(message) if message == text.as_bytes()),
        "{error}"
    );
    assert_eq!(
        error.ncbi_message(),
        Some(format!("Database memory map file error - {text}").into_bytes())
    );
    // blastdbcmd's request errors print the exception's own message.
    assert_eq!(error.ncbi_request_message(), Some(text.into_bytes()));
    // A directory is not a file.
    let directory = dir.0.to_str().unwrap().to_string();
    assert!(matches!(
        read(directory.as_bytes()).unwrap_err(),
        OidListError::NotAFile(path) if path == directory.as_bytes()
    ));
}

#[test]
fn files_without_a_byte_of_bits_are_refused() {
    let dir = TempDir::new("short");
    let empty = dir.write("empty.msk", b"");
    assert!(matches!(
        read(empty.as_bytes()).unwrap_err(),
        OidListError::Empty(_)
    ));
    assert!(matches!(
        read(empty.as_bytes()).unwrap_err().to_error(),
        SeqDbError::NcbiCrashes { .. }
    ));
    for len in 1..=4u64 {
        let path = dir.write("short.msk", &vec![0u8; len as usize]);
        let error = read(path.as_bytes()).unwrap_err();
        assert!(
            matches!(&error, OidListError::Short(at, n) if at.as_path() == std::path::Path::new(&path) && *n == len)
        );
        assert!(matches!(
            error.to_error(),
            SeqDbError::Integrity { message, .. } if message.contains(&format!("holds {len} bytes"))
        ));
    }
    let path = dir.write("five.msk", &list_file(0, &[0x80]));
    assert_eq!(
        oids(&read(path.as_bytes()).unwrap().volume_words(1).unwrap()),
        [0]
    );
}

#[test]
fn bits_are_msb_first_and_volume_relative() {
    let dir = TempDir::new("bits");
    let path = dir.write("m.msk", &list_file(9, &[0b1010_0000, 0b0100_0000]));
    let file = read(path.as_bytes()).unwrap();
    assert_eq!(oids(&file.volume_words(10).unwrap()), [0, 2, 9]);
    // The bits of the last byte past the volume's OIDs are not the volume's.
    assert_eq!(oids(&file.volume_words(2).unwrap()), [0]);
    assert_eq!(oids(&file.volume_words(3).unwrap()), [0, 2]);
    // A bigger volume reads the bytes up to the header's 32-bit word: zeros past the file's
    // end (probe H2).
    assert_eq!(oids(&file.volume_words(100).unwrap()), [0, 2, 9]);
    assert_eq!(file.volume_words(100).unwrap().len(), 2);
}

#[test]
fn the_header_bounds_the_bytes_read() {
    let dir = TempDir::new("header");
    // Last OID 9: (10 + 31) / 32 = 1 word, 4 bytes; the fifth byte is not read (probe H1).
    let path = dir.write("w.msk", &list_file(9, &[0, 0, 0, 0x01, 0xFF]));
    let file = read(path.as_bytes()).unwrap();
    assert_eq!(oids(&file.volume_words(40).unwrap()), [31]);
    // Last OID 32: 2 words, 8 bytes.
    let path = dir.write("w2.msk", &list_file(32, &[0, 0, 0, 0x01, 0xFF]));
    assert_eq!(
        oids(&read(path.as_bytes()).unwrap().volume_words(40).unwrap()),
        [31, 32, 33, 34, 35, 36, 37, 38, 39]
    );
    // A header of 0xFFFFFFFF counts 0 OIDs in `Uint4`: no byte is read (probe H3).
    let path = dir.write("ff.msk", &list_file(u32::MAX, &[0xFF; 8]));
    assert!(oids(&read(path.as_bytes()).unwrap().volume_words(64).unwrap()).is_empty());
}

#[test]
fn a_read_past_the_last_mapped_page_is_refused() {
    let dir = TempDir::new("page");
    // Five bytes map one page of 4096: 4092 bytes of bits after the header are inside it
    // (zeros past the file, probe H7), one more is not.
    let path = dir.write("p.msk", &list_file(0x7FFF_FFFF, &[0x80]));
    let file = read(path.as_bytes()).unwrap();
    assert_eq!(oids(&file.volume_words(4092 * 8).unwrap()), [0]);
    let error = file.volume_words(4092 * 8 + 1).unwrap_err();
    assert!(matches!(
        &error,
        OidListError::PastPage { path: at, read: 4093, len: 5 } if at.as_path() == std::path::Path::new(&path)
    ));
    assert!(matches!(error.to_error(), SeqDbError::NcbiCrashes { .. }));
}
