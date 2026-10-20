//! Tests of the mask-column reader (`column.rs`, `mask.rs`).
//!
//! The fixtures `db/mask_*` were written by makeblastdb 2.17.0 with `-mask_data`
//! (`make_mask_fixtures.py`); `expected/mask_*.algorithms.txt` and `expected/mask_*.masks.tsv`
//! are what NCBI's blastdbcmd reads back from them. Damaged and hand-made column files are
//! built in a temporary directory from the same databases.

use std::path::{Path, PathBuf};

use super::mask::{ncbi_split, string_to_int};
use super::*;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb")
}

fn db_name(name: &str) -> String {
    fixture_dir()
        .join("db")
        .join(name)
        .to_str()
        .expect("fixture path is UTF-8")
        .to_string()
}

fn read_expected(file: &str) -> String {
    std::fs::read_to_string(fixture_dir().join("expected").join(file))
        .unwrap_or_else(|error| panic!("read {file}: {error}"))
}

/// A directory that is removed when the test ends.
struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-seqdb-mask-{tag}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }

    fn path(&self, file: &str) -> PathBuf {
        self.0.join(file)
    }

    fn name(&self, db: &str) -> String {
        self.0.join(db).to_str().unwrap().to_string()
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn copy_files(db: &str, exts: &[&str], dir: &TempDir) {
    for ext in exts {
        let file = format!("{db}.{ext}");
        std::fs::copy(fixture_dir().join("db").join(&file), dir.path(&file)).unwrap();
    }
}

/// Copy of the files of `mask_nucl_v4`: the sequence files and the column `a`.
fn copy_v4(exts: &[&str]) -> TempDir {
    let dir = TempDir::new("v4");
    copy_files("mask_nucl_v4", exts, &dir);
    dir
}

const V4_ALL: [&str; 6] = ["nin", "nsq", "nhr", "naa", "nab", "nac"];

fn patch(path: &Path, at: usize, bytes: &[u8]) {
    let mut data = std::fs::read(path).unwrap();
    data[at..at + bytes.len()].copy_from_slice(bytes);
    std::fs::write(path, data).unwrap();
}

fn truncate(path: &Path, len: usize) {
    let mut data = std::fs::read(path).unwrap();
    data.truncate(len);
    std::fs::write(path, data).unwrap();
}

fn open_masked(name: &str) -> Result<SeqDb, SeqDbError> {
    SeqDb::open_with_masks(name, SeqType::Nucleotide)
}

fn format_message(error: &SeqDbError) -> &str {
    match error {
        SeqDbError::Format { message, .. } => message,
        other => panic!("expected a format error, got {other}"),
    }
}

/// The text of a `Ranges` result as `blastdbcmd -outfmt %m<algo>` prints it.
fn ranges_text(result: &MaskRanges) -> Option<String> {
    match result {
        MaskRanges::Ranges(ranges) if ranges.is_empty() => Some("none".to_string()),
        MaskRanges::Ranges(ranges) => Some(
            ranges
                .iter()
                .map(|(start, end)| format!("{start}-{end};"))
                .collect(),
        ),
        MaskRanges::AlgorithmNotInVolume { .. } => None,
    }
}

/// The databases of `expected/`: name and sequence type.
const DATABASES: [(&str, SeqType); 9] = [
    ("mask_nucl_v4", SeqType::Nucleotide),
    ("mask_nucl_v5", SeqType::Nucleotide),
    ("mask_nucl_mv", SeqType::Nucleotide),
    ("mask_nucl_ids", SeqType::Nucleotide),
    ("mask_nucl_z", SeqType::Nucleotide),
    ("mask_nucl_y", SeqType::Nucleotide),
    ("mask_nucl_opts", SeqType::Nucleotide),
    ("mask_prot_v4", SeqType::Protein),
    ("mask_remap", SeqType::Nucleotide),
];

/// The ids that the table of `expected/<db>.algorithms.txt` lists.
fn expected_ids(text: &str) -> Vec<i32> {
    text.lines()
        .skip(4)
        .map(|line| line.split_whitespace().next().unwrap().parse().unwrap())
        .collect()
}

#[test]
fn descriptions_and_ids_equal_ncbi_for_every_fixture() {
    for (name, seqtype) in DATABASES {
        let db = SeqDb::open_with_masks(&db_name(name), seqtype)
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        let algorithms = db.mask_algorithms().unwrap();
        let expected = read_expected(&format!("{name}.algorithms.txt"));
        assert_eq!(algorithms.descriptions().unwrap(), expected, "{name}");
        assert_eq!(algorithms.ids(), expected_ids(&expected), "{name}");
        assert!(!algorithms.ids().is_empty(), "{name}");
    }
}

/// The rows of `expected/<db>.masks.tsv`: algorithm, ordinal id, ranges.
fn expected_rows(name: &str) -> Vec<(i32, u32, String)> {
    read_expected(&format!("{name}.masks.tsv"))
        .lines()
        .map(|line| {
            let mut fields = line.split('\t');
            let algo = fields.next().unwrap().parse().unwrap();
            let oid = fields.next().unwrap().parse().unwrap();
            let ranges = fields.next().unwrap().to_string();
            assert!(fields.next().is_none(), "{line}");
            (algo, oid, ranges)
        })
        .collect()
}

/// Name of the volume of the sequence `ordinal`, without directory.
fn volume_name(db: &SeqDb, ordinal: u32) -> String {
    let (volume, _) = db.find_volume(ordinal).unwrap();
    db.volumes()[volume]
        .base()
        .file_name()
        .unwrap()
        .to_string_lossy()
        .into_owned()
}

#[test]
fn mask_ranges_equal_ncbi_for_every_fixture_row() {
    for (name, seqtype) in DATABASES {
        let db = SeqDb::open_with_masks(&db_name(name), seqtype).unwrap();
        let rows = expected_rows(name);
        assert!(!rows.is_empty(), "{name}");
        let mut not_in_volume: Vec<(i32, String)> = Vec::new();
        let mut compared = 0;
        for (algo, oid, expected) in &rows {
            let result = db.mask_ranges(*oid, *algo).unwrap();
            match ranges_text(&result) {
                Some(text) => {
                    assert_eq!(&text, expected, "{name}: algorithm {algo}, sequence {oid}");
                    compared += 1;
                }
                None => {
                    let pair = (*algo, volume_name(&db, *oid));
                    if !not_in_volume.contains(&pair) {
                        not_in_volume.push(pair);
                    }
                }
            }
        }
        // Only mask_remap has volumes that lack an algorithm that other volumes have.
        if name == "mask_remap" {
            assert!(compared > 0);
        } else {
            assert!(not_in_volume.is_empty(), "{name}: {not_in_volume:?}");
            assert_eq!(compared, rows.len());
        }
    }
}

/// `mask_remap` is an alias of five volumes (NCBI's order: ids, plain, v4, y, z) with
/// different ids for the same algorithms. In a volume that has a mask column but not the
/// algorithm, NCBI's rows come from its stale cache (for example `101 37 0-100;` in volume
/// z); LOSAT answers `AlgorithmNotInVolume` there and equals the file everywhere else.
#[test]
fn mask_remap_algorithm_not_in_volume_pairs() {
    let db = SeqDb::open_with_masks(&db_name("mask_remap"), SeqType::Nucleotide).unwrap();
    let volumes: Vec<String> = (0..db.num_seqs())
        .step_by(8)
        .map(|ordinal| volume_name(&db, ordinal))
        .collect();
    assert_eq!(
        volumes,
        [
            "mask_nucl_ids",
            "mask_nucl_plain",
            "mask_nucl_v4",
            "mask_nucl_y",
            "mask_nucl_z"
        ]
    );
    let mut pairs: Vec<(i32, String)> = Vec::new();
    for (algo, oid, expected) in expected_rows("mask_remap") {
        let result = db.mask_ranges(oid, algo).unwrap();
        match result {
            MaskRanges::AlgorithmNotInVolume { volume } => {
                let name = volume.file_name().unwrap().to_string_lossy().into_owned();
                assert_eq!(name, volume_name(&db, oid));
                // A volume without a mask column never gets here.
                assert_ne!(name, "mask_nucl_plain");
                let pair = (algo, name);
                if !pairs.contains(&pair) {
                    pairs.push(pair);
                }
            }
            MaskRanges::Ranges(_) => {
                assert_eq!(ranges_text(&result).unwrap(), expected, "{algo} {oid}");
            }
        }
    }
    pairs.sort();
    // Each volume knows the algorithms its own metadata lists:
    // ids: 100 (dustX) and 101 (wmY); v4: 11 and 30; y: 101 (wmY); z: 102 (wmZ).
    let owners: [(&str, &[i32]); 4] = [
        ("mask_nucl_ids", &[100, 101]),
        ("mask_nucl_v4", &[11, 30]),
        ("mask_nucl_y", &[101]),
        ("mask_nucl_z", &[102]),
    ];
    let mut expected_pairs: Vec<(i32, String)> = Vec::new();
    for (volume, known) in owners {
        for algo in [11, 30, 100, 101, 102] {
            if !known.contains(&algo) {
                expected_pairs.push((algo, volume.to_string()));
            }
        }
    }
    // Rows of an algorithm that the volume lacks are `AlgorithmNotInVolume` only for the
    // sequences that have mask data, so the pairs seen are a subset of the lacking ones, and
    // every volume that has any mask data shows at least one pair.
    for pair in &pairs {
        assert!(expected_pairs.contains(pair), "{pair:?}");
    }
    assert_eq!(
        pairs,
        [
            (11, "mask_nucl_ids".to_string()),
            (11, "mask_nucl_y".to_string()),
            (11, "mask_nucl_z".to_string()),
            (30, "mask_nucl_ids".to_string()),
            (30, "mask_nucl_y".to_string()),
            (30, "mask_nucl_z".to_string()),
            (100, "mask_nucl_v4".to_string()),
            (100, "mask_nucl_y".to_string()),
            (100, "mask_nucl_z".to_string()),
            (101, "mask_nucl_v4".to_string()),
            (101, "mask_nucl_z".to_string()),
            (102, "mask_nucl_ids".to_string()),
            (102, "mask_nucl_v4".to_string()),
            (102, "mask_nucl_y".to_string()),
        ]
    );
}

#[test]
fn unknown_algorithm_and_names_of_algorithms() {
    let db = SeqDb::open_with_masks(&db_name("mask_nucl_ids"), SeqType::Nucleotide).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    let table = read_expected("mask_nucl_ids.algorithms.txt");
    let error = algorithms.details(77).unwrap_err();
    assert!(matches!(&error, SeqDbError::Range(_)), "{error}");
    assert_eq!(
        error.to_string(),
        format!("Filtering algorithm ID 77 is not supported.\n{table}")
    );
    assert_eq!(
        algorithms.details(100).unwrap(),
        MaskAlgorithm {
            id: 100,
            program: "100".to_string(),
            program_name: "dustX".to_string(),
            options: "window=64; level=20; linker=1".to_string(),
        }
    );
    assert_eq!(algorithms.id_of("dustX"), Some(100));
    assert_eq!(algorithms.id_of("wmY"), Some(101));
    assert_eq!(algorithms.id_of("nothing"), None);

    // A description of two items is its own name: the key is the whole text.
    let db = SeqDb::open_with_masks(&db_name("mask_nucl_v4"), SeqType::Nucleotide).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    assert_eq!(algorithms.id_of("30:"), Some(30));
    assert_eq!(
        algorithms.details(30).unwrap(),
        MaskAlgorithm {
            id: 30,
            program: "30".to_string(),
            program_name: "windowmasker".to_string(),
            options: String::new(),
        }
    );
    assert_eq!(algorithms.details(11).unwrap().program_name, "dust");

    let db = SeqDb::open_with_masks(&db_name("mask_remap"), SeqType::Nucleotide).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    assert_eq!(algorithms.id_of("dustX"), Some(100));
    assert_eq!(algorithms.id_of("wmY"), Some(101));
    assert_eq!(algorithms.id_of("wmZ"), Some(102));

    // A sequence outside the database.
    assert!(matches!(
        db.mask_ranges(db.num_seqs(), 100),
        Err(SeqDbError::Range(message)) if message == "OID not found"
    ));
}

#[test]
fn search_open_still_rejects_masked_databases() {
    for (name, seqtype) in DATABASES {
        let error = SeqDb::open(&db_name(name), seqtype).unwrap_err();
        assert!(
            matches!(&error, SeqDbError::Unsupported { message, .. } if message.contains("column data")),
            "{name}: {error}"
        );
        assert!(
            SeqDb::open_with_masks(&db_name(name), seqtype).is_ok(),
            "{name}"
        );
    }
}

#[test]
fn a_database_without_columns_has_no_masks() {
    let db = SeqDb::open(&db_name("mask_nucl_plain"), SeqType::Nucleotide).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    assert!(algorithms.ids().is_empty());
    assert_eq!(algorithms.descriptions().unwrap(), "");
    assert_eq!(algorithms.id_of("dust"), None);
    assert_eq!(
        db.mask_ranges(0, 11).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
    // Also when opened to keep columns.
    let db = open_masked(&db_name("mask_nucl_plain")).unwrap();
    assert!(db.mask_algorithms().unwrap().ids().is_empty());
    assert_eq!(
        db.mask_ranges(7, 30).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
    assert!(db.volumes()[0].columns().unwrap().is_empty());
}

#[test]
fn the_column_of_mask_nucl_v4_is_read_as_ncbi_wrote_it() {
    let db = open_masked(&db_name("mask_nucl_v4")).unwrap();
    let columns = db.volumes()[0].columns().unwrap();
    assert_eq!(columns.len(), 1);
    let column = &columns[0];
    assert_eq!(column.title(), b"BlastDb/MaskData");
    assert_eq!(column.num_oids(), 8);
    assert_eq!(column.data_length(), 236);
    assert_eq!(column.date(), b"10/20/2026 12:54:21");
    let metadata: Vec<(&[u8], &[u8])> = column
        .metadata()
        .iter()
        .map(|(key, value)| (key.as_slice(), value.as_slice()))
        .collect();
    assert_eq!(
        metadata,
        [
            (
                b"11".as_slice(),
                b"10:window=64; level=20; linker=1".as_slice()
            ),
            (b"30".as_slice(), b"30:".as_slice())
        ]
    );
}

// Damaged index files. Offsets are those of mask_nucl_v4.naa (hexdump in the comments of
// `make_mask_fixtures.py`'s output): 32 fixed bytes, the title at 32, the date, metadata at
// 0x45, offsets from 0x78.

#[test]
fn damaged_column_index_gives_ncbis_messages() {
    let cases: [(&str, usize, &[u8], &str); 5] = [
        (
            "format_version",
            3,
            &[2],
            "Column file uses unknown format_version.",
        ),
        (
            "column_type",
            7,
            &[2],
            "Column file uses unknown data type.",
        ),
        (
            "offset_size",
            11,
            &[8],
            "Column file uses unsupported offset size.",
        ),
        // The metadata start moves one byte: the title and date no longer end there.
        (
            "metadata_start",
            27,
            &[0x46],
            "CSeqDBColumn: File format error.",
        ),
        // A negative metadata count (the single byte 0x42 is -2).
        (
            "metadata_count",
            0x45,
            &[0x42],
            "CSeqDBColumn: File format error.",
        ),
    ];
    for (tag, at, bytes, message) in cases {
        let dir = copy_v4(&V4_ALL);
        patch(&dir.path("mask_nucl_v4.naa"), at, bytes);
        let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
        let error = db.mask_algorithms().unwrap_err();
        assert_eq!(format_message(&error), message, "{tag}");
        // The next call tries again and fails alike.
        assert_eq!(format_message(&db.mask_algorithms().unwrap_err()), message);
    }
}

#[test]
fn damaged_metadata_gives_ncbis_messages() {
    // The second key "30" (at 0x6b) becomes the first key "11".
    let dir = copy_v4(&V4_ALL);
    let index = dir.path("mask_nucl_v4.naa");
    let bytes = std::fs::read(&index).unwrap();
    assert_eq!(&bytes[0x6b..0x6d], b"30");
    assert_eq!(bytes[0x77], 0);
    patch(&index, 0x6b, b"11");
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "CSeqDBColumn: Error; duplicate metadata key."
    );

    // The NUL that ends the padding is gone.
    let dir = copy_v4(&V4_ALL);
    patch(&dir.path("mask_nucl_v4.naa"), 0x77, b"#");
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "CBlastDbBlob::ReadString: Unterminated string."
    );

    // A string cut by the end of the metadata (the value of "30" grows past it).
    let dir = copy_v4(&V4_ALL);
    patch(&dir.path("mask_nucl_v4.naa"), 0x6d, &[0x20]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "CBlastDbBlob::x_ReadRaw: hit end of data"
    );

    // The metadata count states more pairs than the metadata holds.
    let dir = copy_v4(&V4_ALL);
    patch(&dir.path("mask_nucl_v4.naa"), 0x45, &[0x05]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "CBlastDbBlob::x_ReadRaw: hit end of data"
    );

    // A varint that never ends.
    let dir = copy_v4(&V4_ALL);
    let index = dir.path("mask_nucl_v4.naa");
    let mut data = std::fs::read(&index).unwrap();
    for byte in &mut data[0x45..0x78] {
        *byte = 0x82;
    }
    std::fs::write(&index, data).unwrap();
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "CBlastDbBlob::ReadVarInt: eof while reading integer."
    );
}

#[test]
fn a_column_with_the_wrong_number_of_sequences_is_refused() {
    let dir = copy_v4(&V4_ALL);
    patch(&dir.path("mask_nucl_v4.naa"), 15, &[7]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "Error: column has wrong #oids [7 vs 8]."
    );
}

#[test]
fn two_columns_with_one_title_are_refused() {
    let dir = copy_v4(&V4_ALL);
    std::fs::copy(dir.path("mask_nucl_v4.naa"), dir.path("mask_nucl_v4.nba")).unwrap();
    std::fs::copy(dir.path("mask_nucl_v4.nac"), dir.path("mask_nucl_v4.nbc")).unwrap();
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "Error: duplicate column title [BlastDb/MaskData]."
    );

    // Both faults in one column: the text about the number of sequences wins.
    patch(&dir.path("mask_nucl_v4.nba"), 15, &[7]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(
        format_message(&db.mask_algorithms().unwrap_err()),
        "Error: column has wrong #oids [7 vs 8]."
    );
}

#[test]
fn column_data_files_are_little_endian_or_refused() {
    // An index file without any data file is skipped.
    let dir = copy_v4(&["nin", "nsq", "nhr", "naa"]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert!(db.mask_algorithms().unwrap().ids().is_empty());
    assert_eq!(
        db.mask_ranges(0, 11).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );

    // Only the big-endian data file: NCBI would read its ranges as little-endian.
    let dir = copy_v4(&["nin", "nsq", "nhr", "naa", "nab"]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    let error = db.mask_algorithms().unwrap_err();
    assert!(
        matches!(&error, SeqDbError::Unsupported { message, .. }
            if message.contains("no little-endian data file mask_nucl_v4.nac")),
        "{error}"
    );

    // Only the little-endian data file is used; the big-endian one is not read.
    let dir = copy_v4(&["nin", "nsq", "nhr", "naa", "nac"]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert_eq!(db.mask_algorithms().unwrap().ids(), [11, 30]);
}

#[test]
fn damaged_index_and_data_files_are_format_errors() {
    // An index file shorter than its fixed fields, than the start of its offset array.
    for (len, text) in [
        (20, "fewer than its 32 fixed bytes"),
        (100, "fewer than the start of its offset array"),
    ] {
        let dir = copy_v4(&V4_ALL);
        truncate(&dir.path("mask_nucl_v4.naa"), len);
        let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
        let error = db.mask_algorithms().unwrap_err();
        assert!(format_message(&error).contains(text), "{len}: {error}");
    }

    // The offset array cut after the metadata is read: the last sequence has no entry.
    let dir = copy_v4(&V4_ALL);
    truncate(&dir.path("mask_nucl_v4.naa"), 0x78 + 4 * 8);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert!(db.mask_algorithms().is_ok());
    assert!(db.mask_ranges(6, 11).is_ok());
    let error = db.mask_ranges(7, 11).unwrap_err();
    assert!(
        format_message(&error)
            .contains("offset array of the column ends before the entry of sequence 7"),
        "{error}"
    );

    // A blob range past the end of the data file.
    let dir = copy_v4(&V4_ALL);
    let data = dir.path("mask_nucl_v4.nac");
    let len = std::fs::metadata(&data).unwrap().len() as usize;
    truncate(&data, len - 4);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    assert!(db.mask_ranges(0, 11).is_ok());
    let error = db.mask_ranges(7, 11).unwrap_err();
    assert!(
        format_message(&error).contains("passes the end of the file"),
        "{error}"
    );

    // The data range of a sequence ends before it starts; or starts below 0.
    let dir = copy_v4(&V4_ALL);
    patch(
        &dir.path("mask_nucl_v4.naa"),
        0x78 + 4 * 3,
        &0x1000i32.to_be_bytes(),
    );
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    let error = db.mask_ranges(3, 11).unwrap_err();
    assert!(
        format_message(&error).contains("ends (140) before it starts (4096)"),
        "{error}"
    );
    patch(
        &dir.path("mask_nucl_v4.naa"),
        0x78 + 4 * 3,
        &(-4i32).to_be_bytes(),
    );
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    let error = db.mask_ranges(3, 11).unwrap_err();
    assert!(format_message(&error).contains("starts at -4"), "{error}");

    // Data for no sequences; a negative, zero or too small start of the metadata or of the
    // offset array.
    for (at, bytes, text) in [
        (
            12usize,
            [0u8, 0, 0, 0],
            "states 236 data bytes for no sequences",
        ),
        (24, [0xff, 0xff, 0xff, 0xff], "negative metadata start"),
        (24, [0, 0, 0, 0], "metadata start of 0"),
        (
            28,
            [0, 0, 0, 0x40],
            "starts its offset array (64) before its metadata (69)",
        ),
        (28, [0, 0, 0, 0x45], "has no metadata"),
    ] {
        let dir = copy_v4(&V4_ALL);
        patch(&dir.path("mask_nucl_v4.naa"), at, &bytes);
        let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
        let error = db.mask_algorithms().unwrap_err();
        assert!(format_message(&error).contains(text), "{text}: {error}");
    }

    // A date of length 0 (a title of length 0 is made below).
    let dir = copy_v4(&V4_ALL);
    patch(&dir.path("mask_nucl_v4.naa"), 0x31, &[0]);
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    let error = db.mask_algorithms().unwrap_err();
    assert!(format_message(&error).contains("empty date"), "{error}");
}

// Hand-made columns.

fn varint(value: i64) -> Vec<u8> {
    assert!(value >= 0);
    let mut groups = Vec::new();
    let mut rest = value >> 6;
    while rest > 0 {
        groups.push((rest & 0x7f) as u8 | 0x80);
        rest >>= 7;
    }
    groups.reverse();
    groups.push((value & 0x3f) as u8);
    groups
}

fn var_string(text: &[u8]) -> Vec<u8> {
    let mut out = varint(text.len() as i64);
    out.extend_from_slice(text);
    out
}

/// A mask blob: `(algorithm id, ranges)` per algorithm, ranges as little-endian pairs.
fn mask_blob(entries: &[(i32, &[(u32, u32)])]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&(entries.len() as i32).to_be_bytes());
    for (algo, ranges) in entries {
        out.extend_from_slice(&algo.to_be_bytes());
        out.extend_from_slice(&(ranges.len() as i32).to_be_bytes());
        for (start, end) in *ranges {
            out.extend_from_slice(&start.to_le_bytes());
            out.extend_from_slice(&end.to_le_bytes());
        }
    }
    out
}

/// Writes `<base>.n<letter>a` and `<base>.n<letter>c` the way makeblastdb does.
fn write_column(
    dir: &TempDir,
    base: &str,
    letter: char,
    title: &str,
    metadata: &[(&str, &str)],
    blobs: &[Vec<u8>],
) {
    let mut head = Vec::new();
    head.extend_from_slice(&var_string(title.as_bytes()));
    head.extend_from_slice(&var_string(b"10/20/2026 12:00:00"));
    let meta_data_start = 32 + head.len();
    let mut meta = varint(metadata.len() as i64);
    for (key, value) in metadata {
        meta.extend_from_slice(&var_string(key.as_bytes()));
        meta.extend_from_slice(&var_string(value.as_bytes()));
    }
    // Pad with '#' and one NUL to a multiple of 8 in the file.
    let end = meta_data_start + meta.len();
    let pad = 8 - end % 8;
    meta.extend(std::iter::repeat_n(b'#', pad - 1));
    meta.push(0);
    let offset_array_start = meta_data_start + meta.len();

    let mut data = Vec::new();
    let mut offsets = vec![0i32];
    for blob in blobs {
        data.extend_from_slice(blob);
        offsets.push(data.len() as i32);
    }

    let mut index = Vec::new();
    for value in [1i32, 1, 4, blobs.len() as i32] {
        index.extend_from_slice(&value.to_be_bytes());
    }
    index.extend_from_slice(&(data.len() as i64).to_be_bytes());
    index.extend_from_slice(&(meta_data_start as i32).to_be_bytes());
    index.extend_from_slice(&(offset_array_start as i32).to_be_bytes());
    index.extend_from_slice(&head);
    index.extend_from_slice(&meta);
    for offset in offsets {
        index.extend_from_slice(&offset.to_be_bytes());
    }
    std::fs::write(dir.path(&format!("{base}.n{letter}a")), index).unwrap();
    std::fs::write(dir.path(&format!("{base}.n{letter}c")), data).unwrap();
}

/// A copy of mask_nucl_v4's sequence files (8 sequences) with a column made here.
fn made_column(metadata: &[(&str, &str)], blobs: &[Vec<u8>]) -> (TempDir, String) {
    let dir = TempDir::new("made");
    copy_files("mask_nucl_v4", &["nin", "nsq", "nhr"], &dir);
    write_column(
        &dir,
        "mask_nucl_v4",
        'a',
        "BlastDb/MaskData",
        metadata,
        blobs,
    );
    let name = dir.name("mask_nucl_v4");
    (dir, name)
}

fn eight_empty_blobs() -> Vec<Vec<u8>> {
    vec![Vec::new(); 8]
}

#[test]
fn a_made_column_reads_like_a_fixture() {
    let mut blobs = eight_empty_blobs();
    blobs[2] = mask_blob(&[(7, &[(3, 9), (20, 25)]), (8, &[(1, 2)])]);
    blobs[5] = mask_blob(&[(8, &[(10, 11)])]);
    let (_dir, name) = made_column(&[("7", "20:window=12"), ("8", "10:")], &blobs);
    let db = open_masked(&name).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    assert_eq!(algorithms.ids(), [7, 8]);
    assert_eq!(
        db.mask_ranges(2, 7).unwrap(),
        MaskRanges::Ranges(vec![(3, 9), (20, 25)])
    );
    assert_eq!(
        db.mask_ranges(2, 8).unwrap(),
        MaskRanges::Ranges(vec![(1, 2)])
    );
    assert_eq!(
        db.mask_ranges(5, 7).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
    assert_eq!(
        db.mask_ranges(5, 8).unwrap(),
        MaskRanges::Ranges(vec![(10, 11)])
    );
    // An empty blob: no ranges, whatever the algorithm.
    assert_eq!(
        db.mask_ranges(0, 7).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
    assert_eq!(
        db.mask_ranges(0, 99).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
    // A blob with data in a volume that does not know the algorithm.
    assert!(matches!(
        db.mask_ranges(2, 99).unwrap(),
        MaskRanges::AlgorithmNotInVolume { .. }
    ));
    let text = algorithms.descriptions().unwrap();
    assert!(text.contains(&format!("\n{:<13}{:<40}window=12", 7, "seg")));
    assert!(text.contains(&format!("\n{:<13}{:<40}default options used", 8, "dust")));
}

#[test]
fn an_empty_title_is_refused() {
    let dir = TempDir::new("title");
    copy_files("mask_nucl_v4", &["nin", "nsq", "nhr"], &dir);
    write_column(
        &dir,
        "mask_nucl_v4",
        'a',
        "",
        &[("1", "10:")],
        &eight_empty_blobs(),
    );
    let db = open_masked(&dir.name("mask_nucl_v4")).unwrap();
    let error = db.mask_algorithms().unwrap_err();
    assert!(format_message(&error).contains("empty title"), "{error}");
}

#[test]
fn a_non_integer_algorithm_key_is_refused() {
    let (_dir, name) = made_column(&[("x", "10:")], &eight_empty_blobs());
    let db = open_masked(&name).unwrap();
    let error = db.mask_algorithms().unwrap_err();
    assert!(
        format_message(&error).contains("\"x\" is not an integer"),
        "{error}"
    );
    // Keys that NStr::StringToInt takes.
    let (_dir, name) = made_column(&[("+7", "10:"), ("-3", "20:")], &eight_empty_blobs());
    let db = open_masked(&name).unwrap();
    assert_eq!(db.mask_algorithms().unwrap().ids(), [-3, 7]);
}

#[test]
fn equal_descriptions_in_a_volume_are_refused() {
    // Item 3 of a four-item description is compared, not the whole text.
    let (_dir, name) = made_column(
        &[("1", "5:a:same:x"), ("2", "6:b:same:y")],
        &eight_empty_blobs(),
    );
    let db = open_masked(&name).unwrap();
    let error = db.mask_algorithms().unwrap_err();
    let message = format_message(&error);
    assert!(
        message.starts_with("Error: volume (")
            && message.ends_with(") mask data has duplicates value (6:b:same:y)"),
        "{error}"
    );
    assert!(message.contains("mask_nucl_v4"), "{error}");

    let (_dir, name) = made_column(&[("1", "10:"), ("2", "10:")], &eight_empty_blobs());
    let db = open_masked(&name).unwrap();
    assert!(format_message(&db.mask_algorithms().unwrap_err()).ends_with("duplicates value (10:)"));
}

#[test]
fn stored_descriptions_that_cannot_be_read_are_refused() {
    // Unknown program number; text that is no number; three and one items; none.
    for (description, text) in [
        ("99:opts", "program 99 of the mask algorithm description"),
        ("seg:opts", "is not an integer"),
        ("10:a:b", "Error in stored mask algorithm description data."),
        ("10", "Error in stored mask algorithm description data."),
        ("", "Error in stored mask algorithm description data."),
        (
            "1:2:3:4:5",
            "Error in stored mask algorithm description data.",
        ),
    ] {
        let (_dir, name) = made_column(&[("4", description)], &eight_empty_blobs());
        let db = open_masked(&name).unwrap();
        let algorithms = db.mask_algorithms().unwrap();
        assert_eq!(algorithms.ids(), [4], "{description}");
        let error = algorithms.details(4).unwrap_err();
        assert!(
            format_message(&error).contains(text),
            "{description}: {error}"
        );
        let error = algorithms.descriptions().unwrap_err();
        assert!(
            format_message(&error).contains(text),
            "{description}: {error}"
        );
        // An unknown id names the table, which cannot be built.
        let error = algorithms.details(5).unwrap_err();
        assert!(
            format_message(&error).contains(text),
            "{description}: {error}"
        );
    }
}

#[test]
fn descriptions_of_four_items_and_wide_columns() {
    let long_name = "n".repeat(45);
    let description = format!("tool\u{1}x:o1\u{1}y:{long_name}:text");
    let (_dir, name) = made_column(&[("3", &description)], &eight_empty_blobs());
    let db = open_masked(&name).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    // A program that is not a number is the name; ':' comes back from 0x01.
    assert_eq!(
        algorithms.details(3).unwrap(),
        MaskAlgorithm {
            id: 3,
            program: "tool:x".to_string(),
            program_name: long_name.clone(),
            options: "o1:y".to_string(),
        }
    );
    assert_eq!(algorithms.id_of(&long_name), Some(3));
    // setw(40) does not cut a longer text.
    assert_eq!(
        algorithms.descriptions().unwrap(),
        format!(
            "\nAvailable filtering algorithms applied to database sequences:\n\n\
             {:<13}{:<40}{:<40}\n{:<13}{:<40}{:<40}\n",
            "Algorithm ID", "Algorithm name", "Algorithm options", 3, "tool:x", "o1:y"
        )
    );
}

#[test]
fn mask_blobs_that_end_early_are_format_errors() {
    let mut blobs = eight_empty_blobs();
    // Entry count 2 but only one entry (the algorithm sought is the missing second one).
    blobs[0] = mask_blob(&[(8, &[(1, 2)])]);
    blobs[0][3] = 2;
    // A count of ranges that passes the end of the blob.
    blobs[1] = mask_blob(&[(7, &[(1, 2)])]);
    blobs[1][11] = 9;
    // A negative count of ranges.
    blobs[2] = mask_blob(&[(7, &[(1, 2)])]);
    blobs[2][8..12].copy_from_slice(&(-1i32).to_be_bytes());
    // The entry of another algorithm is skipped past the end; the count then ends the loop.
    blobs[3] = mask_blob(&[(8, &[(1, 2)])]);
    blobs[3][11] = 100;
    // Two entries, the first skipped past the end: the read of the second fails.
    blobs[4] = mask_blob(&[(8, &[(1, 2)]), (7, &[(3, 4)])]);
    blobs[4][11] = 100;
    // A negative number of entries: nothing is read.
    blobs[5] = mask_blob(&[(7, &[(1, 2)])]);
    blobs[5][0..4].copy_from_slice(&(-1i32).to_be_bytes());
    // A blob of fewer than four bytes.
    blobs[6] = vec![0, 0];
    let (_dir, name) = made_column(&[("7", "10:")], &blobs);
    let db = open_masked(&name).unwrap();
    let end = "CBlastDbBlob::x_ReadRaw: hit end of data";
    for oid in [0, 1, 2, 4, 6] {
        let error = db.mask_ranges(oid, 7).unwrap_err();
        assert_eq!(format_message(&error), end, "sequence {oid}");
    }
    assert_eq!(
        db.mask_ranges(3, 7).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
    assert_eq!(
        db.mask_ranges(5, 7).unwrap(),
        MaskRanges::Ranges(Vec::new())
    );
}

#[test]
fn volumes_with_different_ids_for_one_algorithm_are_remapped() {
    // Two volumes: id 30 means "30:" in the first and "10:" in the second, which also has
    // "30:" as id 31. The table keeps 30 for "30:" and gives "10:" the synthetic id 100.
    let dir = TempDir::new("remap");
    for volume in ["mask_nucl_mv.00", "mask_nucl_mv.01"] {
        copy_files(volume, &["nin", "nsq", "nhr"], &dir);
    }
    copy_files("mask_nucl_mv", &["nal"], &dir);
    let mut first = vec![Vec::new(); 4];
    first[1] = mask_blob(&[(30, &[(1, 2)])]);
    let mut second = vec![Vec::new(); 4];
    second[0] = mask_blob(&[(30, &[(5, 6)]), (31, &[(7, 8)])]);
    write_column(
        &dir,
        "mask_nucl_mv.00",
        'a',
        "BlastDb/MaskData",
        &[("30", "30:")],
        &first,
    );
    write_column(
        &dir,
        "mask_nucl_mv.01",
        'a',
        "BlastDb/MaskData",
        &[("30", "10:"), ("31", "30:")],
        &second,
    );
    let db = open_masked(&dir.name("mask_nucl_mv")).unwrap();
    let algorithms = db.mask_algorithms().unwrap();
    assert_eq!(algorithms.ids(), [30, 100]);
    assert_eq!(algorithms.details(100).unwrap().program_name, "dust");
    assert_eq!(
        db.mask_ranges(1, 30).unwrap(),
        MaskRanges::Ranges(vec![(1, 2)])
    );
    // The first volume has 4 sequences: sequence 4 is the first of the second volume.
    assert_eq!(
        db.mask_ranges(4, 30).unwrap(),
        MaskRanges::Ranges(vec![(7, 8)])
    );
    assert_eq!(
        db.mask_ranges(4, 100).unwrap(),
        MaskRanges::Ranges(vec![(5, 6)])
    );
    // The first volume does not know 100.
    assert!(matches!(
        db.mask_ranges(1, 100).unwrap(),
        MaskRanges::AlgorithmNotInVolume { .. }
    ));
}

#[test]
fn the_other_volumes_columns_are_opened_first() {
    // A fault in the second volume's column is reported for a request in the first volume
    // (x_GetColumnId opens the columns of every volume).
    let dir = TempDir::new("order");
    for volume in ["mask_nucl_mv.00", "mask_nucl_mv.01"] {
        copy_files(volume, &["nin", "nsq", "nhr", "naa", "nab", "nac"], &dir);
    }
    copy_files("mask_nucl_mv", &["nal"], &dir);
    patch(&dir.path("mask_nucl_mv.01.naa"), 3, &[9]);
    let db = open_masked(&dir.name("mask_nucl_mv")).unwrap();
    let error = db.mask_ranges(0, 30).unwrap_err();
    assert_eq!(
        format_message(&error),
        "Column file uses unknown format_version."
    );
}

#[test]
fn integers_and_splits_follow_ncbi() {
    for (text, value) in [
        ("0", Some(0)),
        ("30", Some(30)),
        ("+30", Some(30)),
        ("-5", Some(-5)),
        ("007", Some(7)),
        ("2147483647", Some(i32::MAX)),
        ("-2147483648", Some(i32::MIN)),
        ("2147483648", None),
        ("-2147483649", None),
        ("99999999999999999999999999", None),
        ("", None),
        ("+", None),
        ("-", None),
        (" 5", None),
        ("5 ", None),
        ("5x", None),
        ("0x10", None),
        ("1,000", None),
        ("--5", None),
        // The text ends at a NUL, as the C string does.
        ("5\u{0}x", Some(5)),
        ("\u{0}5", None),
    ] {
        assert_eq!(string_to_int(text), value, "{text:?}");
    }
    let split: [(&str, &[&str]); 9] = [
        ("", &[]),
        ("30", &["30"]),
        ("30:", &["30", ""]),
        (":", &["", ""]),
        (":a", &["", "a"]),
        ("a::b", &["a", "", "b"]),
        ("a:b:c:d", &["a", "b", "c", "d"]),
        ("a:b:c:d:", &["a", "b", "c", "d", ""]),
        ("::", &["", "", ""]),
    ];
    for (text, items) in split {
        assert_eq!(ncbi_split(text), items, "{text:?}");
    }
}
