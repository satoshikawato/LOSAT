//! Tests of the GI-mask reader on the files that NCBI BLAST+ 2.17.0's makeblastdb wrote
//! (`tests/seqdb/gi_mask/db/big*`: 1, 512, 513, 1024, 1025 and 1,100 GIs, one data volume, and
//! 1,100 GIs in four data volumes).

use super::*;

fn temp_dir(tag: &str) -> PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-gi-mask-{tag}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/gi_mask/db")
        .join(name)
}

fn be(bytes: &[u8], at: usize) -> u32 {
    u32::from_be_bytes(bytes[at..at + 4].try_into().unwrap())
}

fn le(bytes: &[u8], at: usize) -> u32 {
    u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap())
}

/// Every `(gi, vol, off)` of an offset file, read page by page (512 entries a page) in the byte
/// order `read`.
fn entries(offsets: &[u8], num_gis: usize, read: fn(&[u8], usize) -> u32) -> Vec<(u32, u32, u32)> {
    let mut out = Vec::new();
    for page in 0..num_gis.div_ceil(512) {
        let count = (num_gis - page * 512).min(512);
        let begin = page * 512 * 12;
        for i in 0..count {
            out.push((
                read(offsets, begin + 4 * i),
                read(offsets, begin + 4 * count + 8 * i),
                read(offsets, begin + 4 * count + 8 * i + 4),
            ));
        }
    }
    out
}

/// The ranges of the record at `off` of a data file, in the byte order `read`.
fn record(data: &[u8], off: usize, read: fn(&[u8], usize) -> u32) -> Vec<(u32, u32)> {
    let n = read(data, off) as usize;
    (0..n)
        .map(|i| (read(data, off + 4 + 8 * i), read(data, off + 8 + 8 * i)))
        .collect()
}

/// Each GI of the offset file gets the record that the big-endian copy of the files (`.gmo`,
/// `.gmd`, which LOSAT does not read on a little-endian host) holds for it, through the index and
/// the pages; the numbers just outside every GI that is not in the list get none.
#[cfg(target_endian = "little")]
#[test]
fn every_gi_of_ncbi_files_finds_its_record() {
    for (name, volumes) in [
        ("big1", 1),
        ("big512", 1),
        ("big513", 1),
        ("big1024", 1),
        ("big1025", 1),
        ("big1100", 1),
        ("bigmax", 4),
    ] {
        let dir = fixture(name);
        let search = SearchPath::new(&dir, None, None).unwrap();
        let files = GiMaskFiles::open("big", &search).unwrap();
        assert_eq!(files.data.len(), volumes, "{name}");
        let index = std::fs::read(dir.join("big.gmi")).unwrap();
        let num_gis = be(&index, 24) as usize;
        let big_endian = entries(&std::fs::read(dir.join("big.gmo")).unwrap(), num_gis, be);
        let little_endian = entries(&std::fs::read(dir.join("big.gno")).unwrap(), num_gis, le);
        assert_eq!(big_endian, little_endian, "{name}");
        let data_name = |vol: u32| {
            if volumes == 1 {
                "big.gmd".to_string()
            } else {
                format!("big.{:02}.gmd", vol)
            }
        };
        let gis: std::collections::BTreeSet<u32> = big_endian.iter().map(|e| e.0).collect();
        for &(gi, vol, off) in &big_endian {
            let data = std::fs::read(dir.join(data_name(vol))).unwrap();
            let mut got = Vec::new();
            files.ranges_into(gi, &mut got).unwrap();
            assert_eq!(got, record(&data, off as usize, be), "{name} gi {gi}");
            assert!(!got.is_empty(), "{name} gi {gi}");
            for other in [gi.wrapping_sub(1), gi.wrapping_add(1)] {
                if !gis.contains(&other) {
                    let mut none = Vec::new();
                    files.ranges_into(other, &mut none).unwrap();
                    assert!(none.is_empty(), "{name} gi {other}");
                }
            }
        }
    }
}

/// `s_BinarySearch` as NCBI wrote it: the first and last keys are tested first, then the middle
/// entries; with equal keys the entry that the halving reaches is the answer (survey D7); a key
/// between two keys gives the entry before it, one outside them -1.
#[test]
fn binary_search_is_ncbi_s() {
    let search = |keys: &[u32], key: u32| {
        binary_search(keys.len() as i64, key, |i| Ok(keys[i as usize])).unwrap()
    };
    assert_eq!(search(&[1, 2, 2, 2, 3], 2), (true, 2));
    assert_eq!(search(&[2, 2, 2, 2, 3], 2), (true, 0));
    assert_eq!(search(&[1, 2, 2, 2, 2], 2), (true, 4));
    assert_eq!(search(&[10, 20, 30, 40], 25), (false, 1));
    assert_eq!(search(&[10, 20, 30, 40], 9), (false, -1));
    assert_eq!(search(&[10, 20, 30, 40], 41), (false, -1));
    assert_eq!(search(&[10], 10), (true, 0));
    assert_eq!(search(&[10, 20], 15), (false, 0));
}

/// The fields of a damaged index file: a format version other than 1 is NCBI's text; a file too
/// short for its fields is LOSAT's integrity error, and one whose description and date do not
/// read is an error (NCBI's blob text, or LOSAT's integrity error for an empty one).
#[test]
fn damaged_index_files() {
    let dir = temp_dir("damaged");
    let good = std::fs::read(fixture("big1100").join("big.gni")).unwrap();
    let search = SearchPath::new(&dir, None, None).unwrap();
    let open = |index: &[u8]| {
        for (from, to) in [("big.gno", "x.gno"), ("big.gnd", "x.gnd")] {
            std::fs::copy(fixture("big1100").join(from), dir.join(to)).unwrap();
        }
        std::fs::write(dir.join("x.gni"), index).unwrap();
        GiMaskFiles::open("x", &search).map(|_| ()).unwrap_err()
    };
    let mut version = good.clone();
    version[3] = 2;
    assert_eq!(
        open(&version).ncbi_message().as_deref(),
        Some(b"Gi-mask file uses unknown format_version.".as_slice())
    );
    assert!(matches!(open(&good[..20]), SeqDbError::Integrity { .. }));
    let mut empty_desc = good.clone();
    empty_desc[32] = 0;
    let error = open(&empty_desc);
    assert!(
        matches!(
            &error,
            SeqDbError::Integrity { .. } | SeqDbError::Format { .. }
        ),
        "{error}"
    );
    std::fs::remove_file(dir.join("x.gno")).unwrap();
    std::fs::write(dir.join("x.gni"), &good).unwrap();
    assert_eq!(
        GiMaskFiles::open("x", &search)
            .map(|_| ())
            .unwrap_err()
            .ncbi_message()
            .as_deref(),
        Some(b"Could not open gi-mask index files.".as_slice())
    );
}

/// A data volume that is not found is the memory-map error NCBI's lease gives the empty name
/// (survey D8): "To be memory mapped the file must exist: ''".
#[test]
fn missing_data_volume_is_the_memory_map_error() {
    let dir = temp_dir("novol");
    for file in [
        "big.gni",
        "big.gno",
        "big.00.gnd",
        "big.02.gnd",
        "big.03.gnd",
    ] {
        std::fs::copy(fixture("bigmax").join(file), dir.join(file)).unwrap();
    }
    let search = SearchPath::new(&dir, None, None).unwrap();
    let error = GiMaskFiles::open("big", &search).map(|_| ()).unwrap_err();
    assert!(
        matches!(&error, SeqDbError::MemoryMap(text) if text == b"To be memory mapped the file must exist: ''"),
        "{error}"
    );
    assert_eq!(
        error.ncbi_message().as_deref(),
        Some(
            b"Database memory map file error - To be memory mapped the file must exist: ''"
                .as_slice()
        )
    );
}
