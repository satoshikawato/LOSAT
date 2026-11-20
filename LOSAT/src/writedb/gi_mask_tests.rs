//! Tests of the GI-mask files (`super`): the files of makeblastdb 2.17.0 rebuilt from what they
//! hold (`tests/seqdb/gi_mask`, `make_gi_mask_fixtures.py`), and the refusal of two algorithms
//! whose files are the same files. The tests run no NCBI program.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use super::{files_place, GiMask, GiMasks};
use crate::writedb::tests::{seqdb_fixtures, TempDir};
use crate::writedb::WriteDbError;

const EXTENSIONS: [&str; 6] = ["gmd", "gnd", "gmo", "gno", "gmi", "gni"];

fn dates() -> [Vec<u8>; 2] {
    [
        b"11/20/2026 15:00:00".to_vec(),
        b"11/20/2026 15:00:01".to_vec(),
    ]
}

/// `CBlastDbBlob::x_ReadVarInt`: the most significant bits first, 7 bits a byte with 0x80 for
/// "more", the last byte 6 bits and the sign 0x40.
fn read_var_int(bytes: &[u8], at: &mut usize) -> u64 {
    let mut value = 0u64;
    loop {
        let byte = bytes[*at];
        *at += 1;
        if byte & 0x80 != 0 {
            value = (value << 7) | u64::from(byte & 0x7F);
        } else {
            assert_eq!(byte & 0x40, 0, "a negative length");
            return (value << 6) | u64::from(byte & 0x3F);
        }
    }
}

/// The fields of an index file's header that the rebuild needs.
struct Header {
    num_vols: u32,
    num_gis: usize,
    desc: Vec<u8>,
    date: Vec<u8>,
}

fn header(bytes: &[u8]) -> Header {
    let int = |i: usize| u32::from_be_bytes(bytes[4 * i..4 * i + 4].try_into().unwrap());
    let mut at = 32;
    let mut string = || {
        let length = read_var_int(bytes, &mut at) as usize;
        let text = bytes[at..at + length].to_vec();
        at += length;
        text
    };
    let desc = string();
    let date = string();
    Header {
        num_vols: int(1),
        num_gis: int(6) as usize,
        desc,
        date,
    }
}

/// The entries of a big-endian offset file: GI, data volume and offset, page by page.
fn offset_entries(bytes: &[u8], num_gis: usize) -> Vec<(i64, i32, i32)> {
    let int = |at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().unwrap());
    let mut entries = Vec::with_capacity(num_gis);
    let mut page = 0;
    let mut left = num_gis;
    while left > 0 {
        let n = left.min(512);
        for k in 0..n {
            let offsets = page + 4 * n + 8 * k;
            entries.push((
                i64::from(int(page + 4 * k)),
                int(offsets) as i32,
                int(offsets + 4) as i32,
            ));
        }
        page += 12 * n;
        left -= n;
    }
    assert_eq!(page, bytes.len(), "the offset file holds more than its GIs");
    entries
}

/// The record at `offset` of a big-endian data file.
fn record(data: &[u8], offset: usize) -> Vec<(u32, u32)> {
    let int = |at: usize| u32::from_be_bytes(data[at..at + 4].try_into().unwrap());
    (0..int(offset) as usize)
        .map(|k| (int(offset + 4 + 8 * k), int(offset + 8 + 8 * k)))
        .collect()
}

/// Every file under `dir`, as paths relative to it with `/`.
fn files_under(dir: &Path) -> Vec<String> {
    let mut files = Vec::new();
    let mut stack = vec![PathBuf::new()];
    while let Some(rel) = stack.pop() {
        for entry in std::fs::read_dir(dir.join(&rel)).unwrap() {
            let entry = entry.unwrap();
            let path = rel.join(entry.file_name());
            if entry.file_type().unwrap().is_dir() {
                stack.push(path);
            } else {
                files.push(path.to_string_lossy().replace('\\', "/"));
            }
        }
    }
    files.sort();
    files
}

/// The files of the GI mask `name` among `files`: `name.<ext>` and `name.NN.<ext>`.
fn mask_files<'a>(files: &'a [String], name: &str) -> Vec<&'a String> {
    files
        .iter()
        .filter(|file| {
            file.strip_prefix(name)
                .and_then(|rest| rest.strip_prefix('.'))
                .is_some_and(|rest| {
                    let (volume, extension) = rest.rsplit_once('.').unwrap_or(("", rest));
                    EXTENSIONS.contains(&extension)
                        && (volume.is_empty() || volume.bytes().all(|b| b.is_ascii_digit()))
                })
        })
        .collect()
}

/// The cases of `cases.tsv` whose files LOSAT writes (not those it refuses): the name and the
/// value of `-max_file_sz` (0 without it: every file of these cases is far below the default).
fn fixture_cases(dir: &Path) -> Vec<(String, u64)> {
    let table = std::fs::read_to_string(dir.join("cases.tsv")).unwrap();
    table
        .lines()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let max = fields
                .iter()
                .position(|&word| word == "-max_file_sz")
                .map_or(0, |at| fields[at + 1].parse().unwrap());
            (fields[1] == "same").then(|| (fields[0].to_string(), max))
        })
        .collect()
}

/// The GI mask `name` of the case in `fixture`, written again from the records, GIs, description
/// and dates its files hold; returns the number of GIs. The records go in in the order of their
/// volume and offset, which is the order makeblastdb wrote them.
fn rebuild(fixture: &Path, name: &str, max_file_size: u64, out: &TempDir) -> usize {
    let read = |file: String| std::fs::read(fixture.join(file)).unwrap();
    let head = header(&read(format!("{name}.gmi")));
    let le_date = header(&read(format!("{name}.gni"))).date;
    let entries = offset_entries(&read(format!("{name}.gmo")), head.num_gis);
    let mut records: BTreeMap<(i32, i32), Vec<i64>> = BTreeMap::new();
    for (gi, volume, offset) in entries {
        records.entry((volume, offset)).or_default().push(gi);
    }
    let base = out.base(name);
    if let Some(parent) = Path::new(&base).parent() {
        std::fs::create_dir_all(parent).unwrap();
    }
    let mut mask = GiMask::new(&base, &head.desc, max_file_size, [head.date, le_date]);
    let mut data = BTreeMap::new();
    for ((volume, offset), gis) in &records {
        let bytes = data.entry(*volume).or_insert_with(|| match head.num_vols {
            1 => read(format!("{name}.gmd")),
            _ => read(format!("{name}.{volume:02}.gmd")),
        });
        mask.add_gi_mask(gis, &record(bytes, *offset as usize));
    }
    mask.close();
    assert_eq!(mask.name(), base);
    head.num_gis
}

#[test]
fn gi_masks_rebuild_makeblastdb_s_files_byte_for_byte() {
    let dir = seqdb_fixtures().join("gi_mask");
    let mut rebuilt = Vec::new();
    for (case, max_file_size) in fixture_cases(&dir) {
        let fixture = dir.join("db").join(&case);
        if !fixture.is_dir() {
            continue;
        }
        let expected = files_under(&fixture);
        let names: Vec<&str> = expected
            .iter()
            .filter_map(|file| file.strip_suffix(".gmi"))
            .collect();
        for name in names {
            let out = TempDir::new("gi-mask-rebuild");
            let num_gis = rebuild(&fixture, name, max_file_size, &out);
            let written = files_under(&out.0);
            let want = mask_files(&expected, name);
            assert_eq!(
                written.iter().collect::<Vec<_>>(),
                want,
                "{case}: files of {name}"
            );
            for file in want {
                assert!(
                    std::fs::read(out.0.join(file)).unwrap()
                        == std::fs::read(fixture.join(file)).unwrap(),
                    "{case}: {file} differs from makeblastdb's"
                );
            }
            rebuilt.push(format!("{case}/{name}:{num_gis}"));
        }
    }
    // The page layouts of the survey: 1, 512, 513, 1024, 1025 and 1100 GIs; several data volumes
    // (a first record over the limit, records over the limit); header pads of 19 and 71 bytes of
    // description; duplicate GIs of version 4; a GI over 2^31.
    for want in [
        "big1/big:1",
        "big512/big:512",
        "big513/big:513",
        "big1024/big:1024",
        "big1025/big:1025",
        "big1100/big:1100",
        "bigmax/big:1100",
        "bigmax4/big:1100",
        "d02_10/m:",
        "d02_20/m:",
        "d03/abcdefgh:",
        "d04/nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn:",
        "d01/dg:",
        "g01/ad:",
        "m01/m1:",
        "m01/m2:",
        "x02/gx:",
    ] {
        assert!(
            rebuilt.iter().any(|done| done.starts_with(want)),
            "{want} was not rebuilt: {rebuilt:?}"
        );
    }
}

#[test]
fn an_algorithm_without_gis_writes_nothing_and_loses_its_name() {
    let out = TempDir::new("gi-mask-empty");
    let mut masks = GiMasks::default();
    masks
        .add(2, &out.base("used"), b"dust:", 0, dates())
        .unwrap();
    masks
        .add(3, &out.base("unused"), b"dust:x", 0, dates())
        .unwrap();
    // `SetMaskData` passes no sequence without a GI (`super::super::mask_data`).
    masks.add_gi_mask(2, &[5, 3], &[(1, 4)]);
    // An algorithm that was never registered is ignored, as `m_MaskAlgoMap` has no entry for it.
    masks.add_gi_mask(9, &[7], &[(1, 4)]);
    masks.close();
    let names: Vec<&str> = masks.names().collect();
    assert_eq!(names, [out.base("used").as_str(), ""]);
    let files = files_under(&out.0);
    assert_eq!(
        files,
        ["used.gmd", "used.gmi", "used.gmo", "used.gnd", "used.gni", "used.gno"]
    );
    let gmo = std::fs::read(out.0.join("used.gmo")).unwrap();
    // The GIs sorted, then their volume and offset: one record at offset 0.
    assert_eq!(
        gmo,
        [0, 0, 0, 3, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );
    let gno = std::fs::read(out.0.join("used.gno")).unwrap();
    assert_eq!(
        gno,
        [3, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );
}

fn assert_refused(result: Result<(), WriteDbError>, names: &str) {
    match result {
        Err(WriteDbError::Unsupported(text)) => assert!(
            text.starts_with(&format!(
                "-gi_mask: two masking algorithms would write their GI-mask files under one name \
                 ({names})"
            )),
            "{text}"
        ),
        other => panic!("{names}: {other:?}"),
    }
}

#[test]
fn two_algorithms_with_the_same_files_are_refused() {
    let out = TempDir::new("gi-mask-place");
    std::fs::create_dir(out.0.join("sub")).unwrap();
    let x = out.base("x");
    let add = |masks: &mut GiMasks, id: i32, name: &str| masks.add(id, name, b"dust:", 0, dates());

    let mut masks = GiMasks::default();
    add(&mut masks, 2, &x).unwrap();
    assert_refused(add(&mut masks, 3, &x), &format!("\"{x}\""));
    let dotted = format!("{}/./x", out.0.display());
    assert_refused(
        add(&mut masks, 3, &dotted),
        &format!("\"{x}\" and \"{dotted}\""),
    );
    let up = format!("{}/sub/../x", out.0.display());
    assert_refused(add(&mut masks, 3, &up), &format!("\"{x}\" and \"{up}\""));
    // Other files: another name, a volume name, another directory; and a directory that does not
    // exist, where nothing is written (twice).
    for (id, name) in [
        (4, out.base("y")),
        (5, out.base("x.00")),
        (6, out.base("sub/x")),
        (7, out.base("missing/x")),
        (8, out.base("missing/x")),
    ] {
        add(&mut masks, id, &name).unwrap();
    }
    // Case is part of the name on Linux; Windows and macOS fold it.
    let upper = out.base("X");
    if cfg!(any(windows, target_os = "macos")) {
        assert_refused(
            add(&mut masks, 9, &upper),
            &format!("\"{x}\" and \"{upper}\""),
        );
    } else {
        add(&mut masks, 9, &upper).unwrap();
    }
    // Nothing is written before the first record.
    assert_eq!(files_under(&out.0), Vec::<String>::new());

    // Relative names are read from the current directory.
    let mut relative = GiMasks::default();
    add(&mut relative, 2, "x").unwrap();
    assert_refused(add(&mut relative, 3, "./x"), "\"x\" and \"./x\"");
    assert_eq!(files_place("x"), files_place("./x"));
    assert_eq!(files_place(&out.base("missing/x")), None);
}

#[cfg(unix)]
#[test]
fn a_directory_reached_through_a_link_is_the_same_directory() {
    let out = TempDir::new("gi-mask-link");
    std::os::unix::fs::symlink(&out.0, out.0.join("link")).unwrap();
    let mut masks = GiMasks::default();
    masks.add(2, &out.base("x"), b"dust:", 0, dates()).unwrap();
    let linked = out.base("link/x");
    assert_refused(
        masks.add(3, &linked, b"dust:", 0, dates()),
        &format!("\"{}\" and \"{linked}\"", out.base("x")),
    );
}
