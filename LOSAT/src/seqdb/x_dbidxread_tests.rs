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

/// `IndexFile::open` of `bytes` (written to `path`) with checks of `chunk` entries, against
/// `IndexFile::parse`.
fn check_open(bytes: &[u8], path: &Path, seqtype: SeqType, chunks: &[usize]) {
    std::fs::write(path, bytes).unwrap();
    let expected = summary(&IndexFile::parse(bytes, seqtype, path));
    for &chunk in chunks {
        let file = DbFile::open(path).unwrap();
        let got = summary(&IndexFile::open_checked(file, seqtype, chunk));
        assert_eq!(
            got,
            expected,
            "{} ({} bytes), chunk {chunk}",
            path.display(),
            bytes.len()
        );
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
fn a_header_is_read_alike_alone_and_in_a_range() {
    let ids = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids");
    for (path, seqtype) in [
        (fixture_db_dir().join("prot_v4"), SeqType::Protein),
        (fixture_db_dir().join("prot_mv"), SeqType::Protein),
        (fixture_db_dir().join("nucl_v5"), SeqType::Nucleotide),
        (fixture_db_dir().join("nucl_mv"), SeqType::Nucleotide),
        (ids.join("nucl_v4"), SeqType::Nucleotide),
        (ids.join("nucl_v5"), SeqType::Nucleotide),
        (ids.join("njoined"), SeqType::Nucleotide),
        (ids.join("prot_v4"), SeqType::Protein),
        (ids.join("refuse_v5"), SeqType::Protein),
    ] {
        let db = SeqDb::open(path.to_str().unwrap(), seqtype).unwrap();
        let parsed = check_headers(&db, &path);
        assert_eq!(parsed > 0, path.starts_with(&ids), "{path:?}");
    }
    // A taxid list leaves out sequences and definition lines (`SeqDb::restrict_taxids`).
    for list in [
        TaxidList::Positive([9606].into_iter().collect()),
        TaxidList::Negative([9606].into_iter().collect()),
    ] {
        let path = ids.join("nucl_v5");
        let mut db = SeqDb::open(path.to_str().unwrap(), SeqType::Nucleotide).unwrap();
        let whole = db.num_seqs();
        db.restrict_taxids(&list).unwrap();
        assert!(db.num_seqs() > 0 && db.num_seqs() < whole, "{list:?}");
        assert!(check_headers(&db, &path) > 0, "{list:?}");
    }
    // A seqid list too. The negative list keeps a sequence without any of its lines here, which
    // both reads refuse alike (`no_line_kept`).
    for negative in [false, true] {
        let path = ids.join("nucl_v5");
        let mut db = SeqDb::open(path.to_str().unwrap(), SeqType::Nucleotide).unwrap();
        let whole = db.num_seqs();
        let file =
            super::seqidlist::read_si_list(b"AF000002.3\nAJ000003.1\n", &mut Vec::new()).unwrap();
        let list = DbList::Seqids(SeqidList { negative, file });
        db.restrict(&list, &mut Vec::new()).unwrap();
        assert!(db.num_seqs() > 0 && db.num_seqs() < whole, "{negative}");
        if !negative {
            assert!(check_headers(&db, &path) > 0);
            continue;
        }
        let mut refused = 0;
        for range in db.header_ranges(1024) {
            assert!(db.visit_sequence_headers(range, |ordinal, header| {
                let walked = header.err().map(|error| error.to_string());
                refused += usize::from(walked.is_some());
                let alone = db.read_header(ordinal).err();
                assert_eq!(walked, alone.map(|error| error.to_string()), "{ordinal}");
                true
            }));
        }
        assert!(refused > 0);
    }
}

/// Checks `visit_sequence_headers` and `read_header` of every sequence of `db` against
/// `deflines` and `defline_set`; returns how many records with parsed Seq-ids it saw.
fn check_headers(db: &SeqDb, path: &Path) -> usize {
    let mut parsed = 0;
    for max_seqs in [1, 3, 1024] {
        let mut headers = Vec::new();
        for range in db.header_ranges(max_seqs) {
            assert!(db.visit_sequence_headers(range, |ordinal, header| {
                headers.push((ordinal, header.unwrap()));
                true
            }));
            // `visit` ends the walk.
            let mut calls = 0;
            assert!(!db.visit_sequence_headers(range, |_, _| {
                calls += 1;
                false
            }));
            assert_eq!(calls, 1, "{path:?}");
        }
        let expected: Vec<_> = (0..db.num_seqs())
            .map(|ordinal| {
                // A record with parsed Seq-ids is what `deflines` does not read.
                let header = match db.deflines(ordinal) {
                    Err(error) if error.is_parsed_seqids() => {
                        parsed += 1;
                        SequenceHeader::Parsed(db.defline_set(ordinal).unwrap())
                    }
                    deflines => SequenceHeader::Ordinal(deflines.unwrap()),
                };
                assert_eq!(db.read_header(ordinal).unwrap(), header, "{path:?}");
                (ordinal, header)
            })
            .collect();
        assert_eq!(headers, expected, "{path:?}");
    }
    parsed
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
            // EXPERIMENT (LOSAT_X_BPDBKEEP): the ranges made one at a time from their index.
            assert_eq!(db.x_header_range_count(max_seqs), ranges.len(), "{name}");
            for (index, range) in ranges.iter().enumerate() {
                assert_eq!(db.x_header_range(max_seqs, index), *range, "{name} {index}");
            }
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

/// Writes `words` at entry `entry` of array `array` (0 header, 1 sequence, 2 ambiguity) of the
/// index file `path` in place (the same file, as a writer that edits it during a run would).
fn overwrite_entries(path: &Path, arrays: usize, num_oids: usize, array: usize, words: &[u32]) {
    use std::io::{Seek, SeekFrom, Write};
    let len = std::fs::metadata(path).unwrap().len() as usize;
    let at = len - 4 * (num_oids + 1) * arrays + 4 * (num_oids + 1) * array;
    let mut file = std::fs::OpenOptions::new().write(true).open(path).unwrap();
    file.seek(SeekFrom::Start(at as u64)).unwrap();
    let bytes: Vec<u8> = words.iter().flat_map(|word| word.to_be_bytes()).collect();
    file.write_all(&bytes).unwrap();
}

#[test]
fn an_index_file_changed_after_the_open_gives_errors_not_panics() {
    let dir = TempDir::new("changed");
    for (name, exts, seqtype, arrays) in [
        ("prot_v4", ["pin", "psq", "phr"], SeqType::Protein, 2),
        ("nucl_v4", ["nin", "nsq", "nhr"], SeqType::Nucleotide, 3),
    ] {
        let base = copy_database(name, &exts, &dir.0);
        let index = base.with_extension(exts[0]);
        let pristine = std::fs::read(&index).unwrap();
        let num_oids = Volume::open_as(&base, seqtype, XDbIdxRead::Off)
            .unwrap()
            .num_oids() as usize;
        assert!(num_oids >= 3, "{name}");
        let descending: Vec<u32> = (0..=num_oids as u32).rev().collect();
        let zeros = vec![0u32; num_oids + 1];
        let mut changes: Vec<(usize, &[u32])> =
            vec![(0, &descending), (1, &descending), (1, &zeros)];
        if arrays == 3 {
            changes.push((2, &zeros[..num_oids]));
        }
        for (array, words) in changes {
            std::fs::write(&index, &pristine).unwrap();
            let volume = Volume::open_as(&base, seqtype, XDbIdxRead::On).unwrap();
            overwrite_entries(&index, arrays, num_oids, array, words);
            let n = num_oids as u32;
            let label = format!("{name}, array {array} set to {:?}..", &words[..2]);
            let block = volume.read_block(0, n, 0);
            if array != 0 {
                let error = block.expect_err(&label).to_string();
                assert!(error.contains("the file changed"), "{label}: {error}");
            }
            for first in 0..n {
                let _ = volume.block_end(first, 64);
                let _ = volume.block_bytes(first, n);
                let _ = volume.read_block(first, first + 1, first);
            }
            let mut visited = 0;
            volume.visit_deflines(0, n, 1 << 20, &mut |_, _| {
                visited += 1;
                true
            });
            assert_eq!(visited, n, "{label}");
            if array == 0 {
                // Header offsets read again in a range break the order: the range is read one
                // sequence at a time, and every sequence gets the error of its own read.
                let error = volume
                    .index()
                    .header_offsets(0, n)
                    .expect_err(&label)
                    .to_string();
                assert!(error.contains("the file changed"), "{label}: {error}");
            }
        }
    }
}
