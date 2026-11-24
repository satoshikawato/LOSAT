use super::*;
use crate::dbindex::{test_letters, test_volume, Options};

fn fasta() -> String {
    format!(
        ">s1\n{}\n>s2\n{}\n",
        test_letters(300, 7),
        test_letters(120, 11)
    )
}

#[test]
fn a_legacy_volume_reads_with_its_fixed_stride_and_hint() {
    let mut options = Options::default_options();
    options.hkey_width = 8;
    let (_dir, path) = test_volume("legacy", &fasta(), &options);
    let volume = IndexVolume::open(&path).unwrap();
    assert_eq!(
        (volume.hkey_width, volume.stride, volume.ws_hint),
        (8, STRIDE, LEGACY_WS_HINT)
    );
    assert_eq!((volume.start, volume.stop), (0, 2));
    assert_eq!(volume.hash_table().len(), (1 << 16) + 1);
    assert_eq!(volume.offsets()[0], 0);
    // One local id of the two sequences, one chunk each, and the sentinel chunk.
    assert_eq!(volume.lid_map().len(), 4);
    assert_eq!(volume.chunks().len(), 3);
    assert_eq!(volume.subjects()[..2], [1, 2]);
    assert_eq!(volume.min_offset(), 64);
}

#[test]
fn a_volume_of_the_other_format_reads_its_stride_and_hint() {
    let mut options = Options::default_options();
    options.legacy = false;
    options.hkey_width = 9;
    options.stride = 3;
    options.ws_hint = 16;
    let (_dir, path) = test_volume("new", &fasta(), &options);
    let volume = IndexVolume::open(&path).unwrap();
    assert_eq!(
        (volume.hkey_width, volume.stride, volume.ws_hint),
        (9, 3, 16)
    );
    assert_eq!(volume.min_offset(), 16);
}

#[test]
fn a_volume_outside_its_sections_is_refused() {
    let mut options = Options::default_options();
    options.hkey_width = 8;
    let (_dir, path) = test_volume("bad", &fasta(), &options);
    let bytes = std::fs::read(&path).unwrap();
    let write = |name: &str, bytes: &[u8]| {
        let bad = path.with_file_name(name);
        std::fs::write(&bad, bytes).unwrap();
        IndexVolume::open(&bad).err()
    };
    // Cut short: the packed letters (or a section) end past the file.
    assert!(write("short.idx", &bytes[..bytes.len() - 8]).is_some());
    assert_eq!(
        write("head.idx", &bytes[..20]).as_deref(),
        Some("shorter than its header")
    );
    // A version NCBI rejects.
    let mut version = bytes.clone();
    version[0] = 7;
    assert_eq!(
        write("version.idx", &version).as_deref(),
        Some("wrong index version")
    );
    // A local id whose chunks run past the chunk starts.
    let volume = IndexVolume::open(&path).unwrap();
    let lid_at = volume.lid_map.at * 4;
    drop(volume);
    let mut lid = bytes.clone();
    lid[lid_at + 4..lid_at + 8].copy_from_slice(&1000u32.to_ne_bytes());
    assert_eq!(
        write("lid.idx", &lid).as_deref(),
        Some("a local id outside its sections")
    );
}
