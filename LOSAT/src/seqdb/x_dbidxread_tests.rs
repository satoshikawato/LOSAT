//! EXPERIMENT (LOSAT_X_DBIDXREAD): the index read from the file by position gives what the held
//! arrays give: the same values, the same errors, the same blocks and the same definition lines.

use std::path::{Path, PathBuf};

use super::file::DbFile;
use super::index::IndexFile;
use super::volume::{held_block_end, Volume};
use super::*;

fn fixture_db_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db")
}

/// Every index file of the fixtures, with its sequence type.
fn index_files() -> Vec<(PathBuf, SeqType)> {
    let mut files: Vec<(PathBuf, SeqType)> = std::fs::read_dir(fixture_db_dir())
        .unwrap()
        .filter_map(|entry| {
            let path = entry.unwrap().path();
            match path.extension().and_then(|ext| ext.to_str()) {
                Some("pin") => Some((path, SeqType::Protein)),
                Some("nin") => Some((path, SeqType::Nucleotide)),
                _ => None,
            }
        })
        .collect();
    files.sort_by(|a, b| a.0.cmp(&b.0));
    assert!(files.len() >= 10, "fixture index files: {files:?}");
    files
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-dbidxread-{tag}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The values an index gives, through its accessors, or the error text.
fn summary(index: &Result<IndexFile, SeqDbError>) -> String {
    let index = match index {
        Ok(index) => index,
        Err(error) => return format!("error: {error}"),
    };
    let n = index.num_oids;
    let mut text = format!(
        "{} {:?} {:?} {:?} {:?} {:?} {} {} {}\n",
        index.version,
        index.seqtype,
        index.volume_number,
        index.title,
        index.lmdb_file,
        index.date,
        n,
        index.volume_length,
        index.max_length
    );
    text += &format!(
        "hdr {:?}\nseq {:?}\namb {:?}\n",
        index.header_offsets(0, n).unwrap(),
        index.sequence_offsets(0, n).unwrap(),
        index.ambiguity_offsets(0, n).unwrap()
    );
    for oid in 0..n {
        text += &format!(
            "{oid}: {:?} {:?} {:?} {}\n",
            index.header_range(oid).unwrap(),
            index.sequence_range(oid).unwrap(),
            index.ambiguity_range(oid).unwrap(),
            index.sequence_start(oid).unwrap()
        );
    }
    text
}

/// `IndexFile::open` of `bytes` (written to `path`) with checks of `chunk` entries, with and
/// without the shadow, against `IndexFile::parse`.
fn check_open(bytes: &[u8], path: &Path, seqtype: SeqType, chunks: &[usize]) {
    std::fs::write(path, bytes).unwrap();
    let expected = summary(&IndexFile::parse(bytes, seqtype, path));
    for &chunk in chunks {
        for shadow in [false, true] {
            let file = DbFile::open(path).unwrap();
            let got = summary(&IndexFile::open_checked(file, seqtype, shadow, chunk));
            assert_eq!(
                got,
                expected,
                "{} ({} bytes), chunk {chunk}, shadow {shadow}",
                path.display(),
                bytes.len()
            );
        }
    }
}

#[test]
fn open_reads_what_parse_reads() {
    let dir = TempDir::new("same");
    for (source, seqtype) in index_files() {
        let bytes = std::fs::read(&source).unwrap();
        let path = dir.0.join(source.file_name().unwrap());
        check_open(&bytes, &path, seqtype, &[1, 2, 3, 7, 1 << 18]);
        // Opened as the other sequence type: the same error.
        let other = match seqtype {
            SeqType::Protein => SeqType::Nucleotide,
            SeqType::Nucleotide => SeqType::Protein,
        };
        check_open(&bytes, &path, other, &[1 << 18]);
    }
}

#[test]
fn open_gives_the_errors_of_parse() {
    let dir = TempDir::new("errors");
    for (source, seqtype) in index_files() {
        let bytes = std::fs::read(&source).unwrap();
        let path = dir.0.join(source.file_name().unwrap());
        // Cut short at every length.
        for cut in 0..bytes.len() {
            check_open(&bytes[..cut], &path, seqtype, &[2, 1 << 18]);
        }
        // Every byte set to 0, to 0xFF and raised by one: versions, types, lengths, counts and
        // every entry of the arrays (backwards, equal, no sentinel, ambiguity outside).
        for at in 0..bytes.len() {
            for value in [0u8, 0xFF, bytes[at].wrapping_add(1)] {
                let mut damaged = bytes.clone();
                damaged[at] = value;
                check_open(&damaged, &path, seqtype, &[1, 3, 1 << 18]);
            }
        }
        // Two faults at once: the errors keep `parse`'s order.
        let index = IndexFile::parse(&bytes, seqtype, &source).unwrap();
        let n = index.num_oids as usize;
        if n < 2 {
            continue;
        }
        let entries = n + 1;
        let arrays = match seqtype {
            SeqType::Protein => 2,
            SeqType::Nucleotide => 3,
        };
        let start = bytes.len() - 4 * entries * arrays;
        let entry = |array: usize, i: usize| start + 4 * (array * entries + i);
        let word =
            |bytes: &[u8], at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().unwrap());
        for (a, b) in [(0, 1), (1, 0), (1, 2), (2, 1), (0, 2)] {
            if a >= arrays || b >= arrays {
                continue;
            }
            for (i, j) in [(n, 1), (1, n), (n - 1, n)] {
                let mut damaged = bytes.clone();
                // Entry i of array a equal to its predecessor, entry j of array b set to 0.
                let previous = word(&damaged, entry(a, i - 1));
                damaged[entry(a, i)..entry(a, i) + 4].copy_from_slice(&previous.to_be_bytes());
                damaged[entry(b, j)..entry(b, j) + 4].copy_from_slice(&0u32.to_be_bytes());
                check_open(&damaged, &path, seqtype, &[1, 2, 1 << 18]);
            }
        }
    }
}

/// The volumes of the fixture databases, each opened with the arrays held and read.
fn volume_pairs() -> Vec<(String, Volume, Volume)> {
    index_files()
        .into_iter()
        .filter_map(|(path, seqtype)| {
            let base = path.with_extension("");
            let held = Volume::open_as(&base, seqtype, XDbIdxRead::Off);
            let read = Volume::open_as(&base, seqtype, XDbIdxRead::On);
            assert_eq!(
                held.as_ref().err().map(ToString::to_string),
                read.as_ref().err().map(ToString::to_string),
                "{}",
                path.display()
            );
            Some((path.display().to_string(), held.ok()?, read.ok()?))
        })
        .collect()
}

#[test]
fn blocks_end_where_the_held_arrays_end_them() {
    for (name, held, read) in volume_pairs() {
        assert!(held.index().holds_offsets() && !read.index().holds_offsets());
        let n = held.num_oids();
        let offsets = held.index().sequence_offsets(0, n).unwrap();
        let lead = match held.index().seqtype {
            SeqType::Protein => 1,
            SeqType::Nucleotide => 0,
        };
        let total = u64::from(offsets[n as usize]);
        let mut budgets: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 8, 16, 64, 100, 200, 1 << 20, u64::MAX];
        budgets.extend(offsets.iter().flat_map(|&o| {
            let o = u64::from(o);
            [o.saturating_sub(1), o, o + 1]
        }));
        budgets.push(total);
        for first in 0..n {
            for &max_bytes in &budgets {
                let expected = held_block_end(&offsets, first, max_bytes, lead);
                assert_eq!(held.block_end(first, max_bytes).unwrap(), expected);
                assert_eq!(
                    read.block_end(first, max_bytes).unwrap(),
                    expected,
                    "{name}"
                );
                for window in [1, 2, 3, 1024] {
                    assert_eq!(
                        read.x_read_block_end(first, max_bytes, window).unwrap(),
                        expected,
                        "{name} from {first} within {max_bytes}, window {window}"
                    );
                }
                assert_eq!(
                    held.block_bytes(first, expected).unwrap(),
                    read.block_bytes(first, expected).unwrap()
                );
            }
        }
    }
}

#[test]
fn blocks_read_alike() {
    for (name, held, read) in volume_pairs() {
        let n = held.num_oids();
        for first in 0..n {
            for end in first + 1..=n {
                let a = held.read_block(first, end, 7).unwrap();
                let b = read.read_block(first, end, 7).unwrap();
                assert_eq!(format!("{a:?}"), format!("{b:?}"), "{name} {first}..{end}");
            }
            assert_eq!(
                format!("{:?}", held.check_ambiguity(first, first)),
                format!("{:?}", read.check_ambiguity(first, first)),
                "{name} {first}"
            );
            assert_eq!(
                held.header_bytes(first).unwrap(),
                read.header_bytes(first).unwrap()
            );
        }
    }
}

/// The definition lines that `visit_deflines` hands out for `first..end`, with pieces of
/// `piece_bytes`, and the number of calls before it stopped (`stop_at`: the call that returns
/// false).
fn visited(
    volume: &Volume,
    first: u32,
    end: u32,
    piece_bytes: u64,
    stop_at: Option<u32>,
) -> (bool, Vec<String>) {
    let mut seen = Vec::new();
    let all = volume.visit_deflines(first, end, piece_bytes, &mut |oid, deflines| {
        seen.push(format!(
            "{oid}: {:?}",
            deflines.map_err(|error| error.to_string())
        ));
        Some(oid) != stop_at
    });
    (all, seen)
}

fn one_by_one(volume: &Volume, first: u32, end: u32, stop_at: Option<u32>) -> (bool, Vec<String>) {
    let mut seen = Vec::new();
    let mut all = true;
    for oid in first..end {
        seen.push(format!(
            "{oid}: {:?}",
            volume.deflines(oid).map_err(|error| error.to_string())
        ));
        if Some(oid) == stop_at {
            all = false;
            break;
        }
    }
    (all, seen)
}

#[test]
fn definition_lines_read_in_ranges_are_those_read_one_by_one() {
    for (name, _, read) in volume_pairs() {
        let n = read.num_oids();
        for first in 0..n {
            for end in first + 1..=n {
                for piece_bytes in [0, 1, 40, 100, 1 << 20] {
                    for stop_at in [None, Some(first), Some(end - 1)] {
                        assert_eq!(
                            visited(&read, first, end, piece_bytes, stop_at),
                            one_by_one(&read, first, end, stop_at),
                            "{name} {first}..{end}, pieces of {piece_bytes}, stop at {stop_at:?}"
                        );
                    }
                }
            }
        }
    }
}

/// A copy of the fixture database `name` (its volume files `exts`) in `dir`.
fn copy_database(name: &str, exts: &[&str], dir: &Path) -> PathBuf {
    for ext in exts {
        let file = format!("{name}.{ext}");
        std::fs::copy(fixture_db_dir().join(&file), dir.join(&file)).unwrap();
    }
    dir.join(name)
}

#[test]
fn damaged_headers_read_in_ranges_as_one_by_one() {
    let dir = TempDir::new("headers");
    let base = copy_database("prot_v4", &["pin", "psq", "phr"], &dir.0);
    let phr = base.with_extension("phr");
    let bytes = std::fs::read(&phr).unwrap();
    // The header file cut short (records pass its end), a damaged byte in every record, and no
    // header file.
    let mut damages: Vec<Option<Vec<u8>>> = vec![None];
    for cut in [0, 1, bytes.len() / 3, bytes.len() / 2, bytes.len() - 1] {
        damages.push(Some(bytes[..cut].to_vec()));
    }
    for at in (0..bytes.len()).step_by(7) {
        let mut damaged = bytes.clone();
        damaged[at] ^= 0xA5;
        damages.push(Some(damaged));
    }
    for damage in damages {
        match &damage {
            Some(bytes) => std::fs::write(&phr, bytes).unwrap(),
            None => std::fs::remove_file(&phr).unwrap(),
        }
        let volume = Volume::open_as(&base, SeqType::Protein, XDbIdxRead::On).unwrap();
        let n = volume.num_oids();
        for (first, end) in [(0, n), (0, 1), (n / 2, n), (n - 1, n)] {
            for piece_bytes in [0, 50, 1 << 20] {
                assert_eq!(
                    visited(&volume, first, end, piece_bytes, None),
                    one_by_one(&volume, first, end, None),
                    "{first}..{end}, pieces of {piece_bytes}, {} header bytes",
                    damage.as_ref().map_or(0, Vec::len)
                );
            }
        }
        if damage.is_none() {
            std::fs::write(&phr, &bytes).unwrap();
        }
    }
}

#[test]
fn database_ranges_cover_every_sequence_once() {
    for name in ["prot_v4", "prot_mv", "nucl_v5", "nucl_mv", "empty"] {
        let seqtype = if name.starts_with("nucl") {
            SeqType::Nucleotide
        } else {
            SeqType::Protein
        };
        let db = SeqDb::open(fixture_db_dir().join(name).to_str().unwrap(), seqtype).unwrap();
        for max_seqs in [0, 1, 2, 3, 1024] {
            let ranges = db.header_ranges(max_seqs);
            let mut next = 0;
            for range in &ranges {
                assert_eq!(range.first_ordinal, next, "{name}");
                assert!(!range.is_empty() && range.len() <= max_seqs.max(1));
                assert!(range.end <= db.volumes()[range.volume].num_oids());
                next = range.end_ordinal();
                let mut titles = Vec::new();
                assert!(db.visit_deflines(*range, |ordinal, deflines| {
                    titles.push((ordinal, deflines.unwrap()));
                    true
                }));
                let expected: Vec<_> = (range.first_ordinal..range.end_ordinal())
                    .map(|ordinal| (ordinal, db.deflines(ordinal).unwrap()))
                    .collect();
                assert_eq!(titles, expected, "{name}");
                let mut checked = Vec::new();
                assert!(db.visit_ambiguity(*range, |ordinal, result| {
                    checked.push((ordinal, result.is_ok()));
                    true
                }));
                let expected: Vec<_> = (range.first_ordinal..range.end_ordinal())
                    .map(|ordinal| (ordinal, db.check_ambiguity(ordinal).is_ok()))
                    .collect();
                assert_eq!(checked, expected, "{name}");
            }
            assert_eq!(next, db.num_seqs(), "{name}");
        }
    }
}

#[test]
fn damaged_ambiguity_data_check_in_ranges_as_one_by_one() {
    let dir = TempDir::new("ambiguity");
    let base = copy_database("nucl_v4", &["nin", "nsq", "nhr"], &dir.0);
    let nsq = base.with_extension("nsq");
    let bytes = std::fs::read(&nsq).unwrap();
    let clean = Volume::open_as(&base, SeqType::Nucleotide, XDbIdxRead::Off).unwrap();
    let n = clean.num_oids();
    let mut damages = Vec::new();
    for oid in 0..n {
        let (at, next) = clean.index().ambiguity_range(oid).unwrap().unwrap();
        if at < next {
            // The entry count raised by one, and the data cut off inside the sequence's data.
            let mut damaged = bytes.clone();
            damaged[at as usize + 3] = damaged[at as usize + 3].wrapping_add(1);
            damages.push(damaged);
            damages.push(bytes[..next as usize - 1].to_vec());
        }
    }
    drop(clean);
    assert!(damages.len() >= 2, "the fixture has ambiguity data");
    for damaged in damages {
        std::fs::write(&nsq, &damaged).unwrap();
        let read = Volume::open_as(&base, SeqType::Nucleotide, XDbIdxRead::On).unwrap();
        for (first, end) in [(0, n), (1, n), (0, 1), (n - 1, n)] {
            let mut seen = Vec::new();
            read.visit_ambiguity(first, end, 100 + first, &mut |oid, result| {
                seen.push(format!(
                    "{oid}: {:?}",
                    result.map_err(|error| error.to_string())
                ));
                true
            });
            let expected: Vec<String> = (first..end)
                .map(|oid| {
                    format!(
                        "{oid}: {:?}",
                        read.check_ambiguity(oid, 100 + oid)
                            .map_err(|error| error.to_string())
                    )
                })
                .collect();
            assert_eq!(seen, expected, "{first}..{end}");
        }
    }
}
