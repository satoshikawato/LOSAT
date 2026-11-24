//! Tests of the GI-mask reader on the files that NCBI BLAST+ 2.17.0's makeblastdb wrote
//! (`tests/seqdb/gi_mask/db/big*`: 1, 512, 513, 1024, 1025 and 1,100 GIs, one data volume, and
//! 1,100 GIs in four data volumes; `tests/seqdb/gi_mask_read/db/sv4` and `big` with their `.nog`
//! cut or edited), and of the reads across chunks and windows.

use super::*;

fn temp_dir(tag: &str) -> crate::utils::test_dir::TempDir {
    crate::utils::test_dir::TempDir::new(&format!(
        "losat-gi-mask-{tag}-{}",
        crate::utils::test_dir::process_id()
    ))
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
        let search = SearchPath::new(&dir, None, None);
        let files = GiMaskFiles::open(b"big", &search).unwrap();
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
        // One lookup and one window for every lookup (in the order of the offset file, not of
        // the data file), and a new lookup and window for each.
        let mut lookup = files.lookup();
        let mut window = DataWindow::default();
        for &(gi, vol, off) in &big_endian {
            let data = std::fs::read(dir.join(data_name(vol))).unwrap();
            let mut got = Vec::new();
            lookup.ranges_into(gi, &mut got, &mut window).unwrap();
            assert_eq!(got, record(&data, off as usize, be), "{name} gi {gi}");
            assert!(!got.is_empty(), "{name} gi {gi}");
            let mut alone = Vec::new();
            files
                .lookup()
                .ranges_into(gi, &mut alone, &mut DataWindow::default())
                .unwrap();
            assert_eq!(alone, got, "{name} gi {gi}");
            for other in [gi.wrapping_sub(1), gi.wrapping_add(1)] {
                if !gis.contains(&other) {
                    let mut none = Vec::new();
                    lookup.ranges_into(other, &mut none, &mut window).unwrap();
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
    let search = SearchPath::new(&dir, None, None);
    let open = |index: &[u8]| {
        for (from, to) in [("big.gno", "x.gno"), ("big.gnd", "x.gnd")] {
            std::fs::copy(fixture("big1100").join(from), dir.join(to)).unwrap();
        }
        std::fs::write(dir.join("x.gni"), index).unwrap();
        GiMaskFiles::open(b"x", &search).map(|_| ()).unwrap_err()
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
        GiMaskFiles::open(b"x", &search)
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
    let search = SearchPath::new(&dir, None, None);
    let error = GiMaskFiles::open(b"big", &search).map(|_| ()).unwrap_err();
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

/// A read that crosses the 64 KiB chunks of a `MappedFile` (the index, offset and `.nog` files),
/// or the windows of a reader of a data file (`DataWindow`), gives the file's bytes (review S3 of
/// M4c-2: every fixture is smaller than a chunk).
#[test]
fn reads_across_chunks_and_windows_are_the_file_s_bytes() {
    let dir = temp_dir("chunks");
    let c = CHUNK as usize;
    let bytes: Vec<u8> = (0..3 * c + 100).map(|i| (i * 7 + i / 251) as u8).collect();
    let path = dir.join("f");
    std::fs::write(&path, &bytes).unwrap();
    let file = MappedFile::open(&path, INDEX_FILE).unwrap();
    for start in (c - 12..c + 12).chain(2 * c - 12..2 * c + 12) {
        for len in 1..=24 {
            let mut out = vec![0u8; len];
            file.read_at(start as i64, &mut out).unwrap();
            assert_eq!(out, bytes[start..start + len], "{start} {len}");
        }
    }
    // The numbers of the arrays through one view (the chunk it holds, another chunk, back) and
    // through a view each: in one chunk, across a chunk's end, the last ones, outside.
    let mut view = file.view();
    for start in (c - 6..c + 2)
        .chain(2 * c - 6..2 * c + 2)
        .chain([10, 2 * c + 8, c + 8, 11])
        .chain(bytes.len() - 6..bytes.len() - 3)
    {
        let want: [u8; 4] = bytes[start..start + 4].try_into().unwrap();
        assert_eq!(
            view.u32(start as i64).unwrap(),
            u32::from_ne_bytes(want),
            "{start}"
        );
        assert_eq!(
            file.view().u32(start as i64).unwrap(),
            u32::from_ne_bytes(want),
            "{start}"
        );
        assert_eq!(
            view.mapped::<4>(start as i64).unwrap(),
            Some(want),
            "{start}"
        );
        assert_eq!(
            file.view().mapped::<4>(start as i64).unwrap(),
            Some(want),
            "{start}"
        );
    }
    for start in [bytes.len() as i64 - 3, -1, -(CHUNK as i64), i64::MIN] {
        for error in [view.u32(start), file.view().u32(start)] {
            let error = error.unwrap_err();
            assert!(matches!(&error, SeqDbError::Format { .. }), "{error}");
        }
    }
    // One read over three chunks, the last bytes, and reads outside the file.
    let mut out = vec![0u8; c + 40];
    file.read_at(c as i64 - 20, &mut out).unwrap();
    assert_eq!(out, bytes[c - 20..2 * c + 20]);
    let last = bytes.len() - 7;
    let mut out = [0u8; 7];
    file.read_at(last as i64, &mut out).unwrap();
    assert_eq!(out, bytes[last..]);
    for (start, len) in [(last as i64, 8), (-1, 2)] {
        let error = file.read_at(start, &mut vec![0u8; len]).unwrap_err();
        assert!(
            matches!(&error, SeqDbError::Format { message, .. } if message.contains(INDEX_FILE)),
            "{error}"
        );
    }
    // A data file through one window: held bytes, a read past the window, a read longer than a
    // window, and the bytes of another file at the same offsets.
    let data = DbFile::open(&path).unwrap();
    let mut window = DataWindow::default();
    for (start, len) in [
        (c - 4, 8),
        (c + 4, 8),
        (2 * c - 2, 4),
        (10, c + 50),
        (c - 30, 30),
        (3 * c + 90, 10),
    ] {
        let got = window.get((7, 0), &data, start as u64, len as u64).unwrap();
        assert_eq!(got, &bytes[start..start + len], "{start} {len}");
    }
    let other: Vec<u8> = bytes.iter().map(|b| !b).collect();
    std::fs::write(dir.join("g"), &other).unwrap();
    let other_file = DbFile::open(&dir.join("g")).unwrap();
    for of in [(8, 0), (7, 1)] {
        let got = window.get(of, &other_file, 3 * c as u64 + 90, 10).unwrap();
        assert_eq!(got, &other[3 * c + 90..3 * c + 100], "{of:?}");
        window.get((7, 0), &data, 3 * c as u64 + 90, 10).unwrap();
    }
}

/// The masks of `oids` (all with `None`) of the GI-mask database `db` copied from
/// `tests/seqdb/gi_mask_read/db/<db>` with its `.nog` changed by `edit`, read through the search
/// path of the copy's directory (where the mask files are), for the algorithm 0.
fn nog_masks(
    tag: &str,
    db: &str,
    edit: impl FnOnce(&Path),
    oids: Option<&[u32]>,
) -> Vec<Result<Vec<(u32, u32)>, SeqDbError>> {
    let dir = temp_dir(tag);
    let from = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/gi_mask_read/db")
        .join(db);
    for entry in std::fs::read_dir(&from).unwrap() {
        let entry = entry.unwrap();
        std::fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
    }
    edit(&dir.join("db.nog"));
    let search = SearchPath::new(&dir, None, None);
    let db = super::super::SeqDb::open_in(b"db", SeqType::Nucleotide, &search).unwrap();
    let all: Vec<u32> = (0..db.num_oids()).collect();
    oids.unwrap_or(&all)
        .iter()
        .map(|&oid| {
            let at = db.locate_oid(oid).unwrap();
            db.mask_ranges_at(at, 0).map(|found| match found {
                super::super::mask::MaskRanges::Ranges(ranges) => ranges,
                other => panic!("{other:?}"),
            })
        })
        .collect()
}

fn unwrapped(masks: Vec<Result<Vec<(u32, u32)>, SeqDbError>>) -> Vec<Vec<(u32, u32)>> {
    masks.into_iter().map(Result::unwrap).collect()
}

fn cut(len: u64) -> impl FnOnce(&Path) {
    move |path: &Path| {
        std::fs::OpenOptions::new()
            .write(true)
            .open(path)
            .unwrap()
            .set_len(len)
            .unwrap()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:52-71
/// The `.nog` as NCBI's mapping reads it (review S4 and N12 of M4c-2; the NCBI runs of the same
/// files are the rows `bdc_sv4n*` of `tests/seqdb/gi_mask_read`): cut inside its header, every
/// sequence takes the GI of its definition lines, as without the file; cut inside its entries, a
/// sequence whose entry is cut off reads GI 0 (or the bytes left, then zeros); a count below the
/// volume's sends the sequences after it to their definition lines; an empty file and an entry
/// past the last 4 KiB page of the mapping are LOSAT's refusals.
#[test]
fn the_nog_reads_as_ncbi_s_mapping_of_it() {
    let normal = unwrapped(nog_masks("nog-normal", "sv4", |_| {}, None));
    assert_eq!(normal.len(), 12);
    assert!(normal[0] != normal[1] && !normal[0].is_empty() && !normal[2].is_empty());
    for (tag, edit) in [
        ("nog-12", Box::new(cut(12)) as Box<dyn FnOnce(&Path)>),
        ("nog-15", Box::new(cut(15))),
        (
            "nog-none",
            Box::new(|path: &Path| std::fs::remove_file(path).unwrap()),
        ),
    ] {
        let got = unwrapped(nog_masks(tag, "sv4", edit, None));
        assert_eq!(got, normal, "{tag}");
    }
    // Cut at 42: the entries of OIDs 0 and 1, then two bytes of OID 2's.
    let got = unwrapped(nog_masks("nog-42", "sv4", cut(42), None));
    assert_eq!(got[..2], normal[..2]);
    assert!(got[2..].iter().all(Vec::is_empty), "{got:?}");
    // A count of 3, with the entries of OIDs 0 and 1 swapped.
    let swap = |path: &Path| {
        let mut bytes = std::fs::read(path).unwrap();
        bytes[12..16].copy_from_slice(&3u32.to_be_bytes());
        let first: [u8; 4] = bytes[32..36].try_into().unwrap();
        bytes.copy_within(36..40, 32);
        bytes[36..40].copy_from_slice(&first);
        std::fs::write(path, bytes).unwrap();
    };
    let got = unwrapped(nog_masks("nog-few", "sv4", swap, None));
    assert_eq!((&got[0], &got[1]), (&normal[1], &normal[0]));
    assert_eq!(got[2..], normal[2..]);
    // Empty: NCBI maps a null pointer and crashes.
    for found in nog_masks("nog-0", "sv4", cut(0), Some(&[0, 5])) {
        let error = found.unwrap_err();
        assert!(
            error.ncbi_crashes() && error.to_string().contains(NOG_FILE),
            "{error}"
        );
    }
    // `big` (1,100 OIDs) cut to 100 bytes: OIDs 17 to 1015 read zeros, OID 1016 is past the page.
    let got = nog_masks(
        "nog-page",
        "big",
        cut(100),
        Some(&[16, 17, 1015, 1016, 1099]),
    );
    assert!(!got[0].as_ref().unwrap().is_empty());
    assert!(got[1].as_ref().unwrap().is_empty() && got[2].as_ref().unwrap().is_empty());
    for error in got[3..].iter().map(|found| found.as_ref().unwrap_err()) {
        assert!(
            matches!(error, SeqDbError::Unsupported { message, .. }
                if message.contains("which has 100 bytes, outside the 4096 bytes of NCBI")),
            "{error}"
        );
    }
}
