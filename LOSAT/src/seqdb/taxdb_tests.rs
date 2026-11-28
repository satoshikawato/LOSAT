//! Tests against taxdb files written by `make_taxdb.py` (the fixtures under
//! `tests/seqdb/taxdb/<variant>/`, with what NCBI BLAST+ 2.17.0 returned
//! for them in `tests/seqdb/taxdb/expected.tsv`) and against files built
//! here for the cases NCBI cannot be asked about.

use super::super::file::MAP_PAGE;
use super::*;
use crate::utils::test_dir::TempDir;
use std::path::Path;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb")
}

/// A search path with the empty current directory it was made with, which is removed with it.
struct Searching {
    search: SearchPath,
    _cwd: TempDir,
}

impl std::ops::Deref for Searching {
    type Target = SearchPath;

    fn deref(&self) -> &SearchPath {
        &self.search
    }
}

/// A search path whose only useful directory is `dir` (`BLASTDB`), run from an empty current directory.
fn path_of(dir: &Path) -> Searching {
    let cwd = TempDir::unique("taxdb-cwd");
    let search = SearchPath::new(&cwd, Some(dir.as_os_str()), None);
    Searching { search, _cwd: cwd }
}

fn open_fixture(variant: &str) -> Result<TaxDb, SeqDbError> {
    TaxDb::open(&path_of(&fixture_dir().join(variant)))
}

fn be(v: u32) -> [u8; 4] {
    v.to_be_bytes()
}

/// Builds a taxdb: `rows` are (taxid, record), laid out back to back in the order given.
fn build(rows: &[(i32, &[u8])]) -> (Vec<u8>, Vec<u8>) {
    let mut index = Vec::new();
    index.extend(be(MAGIC));
    index.extend(be(rows.len() as u32));
    index.extend([0u8; 16]);
    let mut data = Vec::new();
    for (taxid, record) in rows {
        index.extend(be(*taxid as u32));
        index.extend(be(data.len() as u32));
        data.extend(*record);
    }
    (index, data)
}

fn write_db(index: &[u8], data: &[u8]) -> TempDir {
    let dir = TempDir::unique("taxdb-db");
    std::fs::write(dir.join("taxdb.bti"), index).unwrap();
    std::fs::write(dir.join("taxdb.btd"), data).unwrap();
    dir
}

/// Writes a taxdb and opens it: the directory is returned first, so that it is dropped after
/// the open database.
fn open_built(index: &[u8], data: &[u8]) -> (TempDir, TaxDb) {
    let dir = write_db(index, data);
    let db = TaxDb::open(&path_of(&dir)).unwrap();
    (dir, db)
}

fn names_of(db: &TaxDb, taxid: i32) -> (Option<TaxNames>, Vec<String>) {
    let mut messages = Vec::new();
    let names = db.names(taxid, &mut messages).unwrap();
    (names, messages)
}

// seqdbtax.cpp:258-259 (`info.taxid` is set to `tax_id` once the search finds the taxid, before the
// lookup can still fail): blastdbcmd's `-tax_info` writes empty names then, `N/A` only for a
// taxid the search does not find (probe of line RG T5b, `T/t5b/out/extra`).
#[test]
fn a_lookup_tells_a_taxid_not_found_from_one_found_without_names() {
    let lookup = |variant: &str, taxid: i32| {
        let mut messages = Vec::new();
        let found = open_fixture(variant)
            .unwrap()
            .lookup(taxid, &mut messages)
            .unwrap();
        (found, messages)
    };
    assert!(matches!(lookup("good", 9606), (TaxLookup::Named(_), m) if m.is_empty()));
    for taxid in [150, 0, 77777, 100000] {
        assert_eq!(lookup("good", taxid), (TaxLookup::NotInIndex, Vec::new()));
    }
    // A record without three tabs, an empty record of the last entry, the offset error.
    assert_eq!(lookup("tabs0", 9606), (TaxLookup::Unnamed, Vec::new()));
    assert_eq!(lookup("tabs1", 9606), (TaxLookup::Unnamed, Vec::new()));
    assert_eq!(
        lookup("offset_last_at_end", 99990),
        (TaxLookup::Unnamed, Vec::new())
    );
    assert_eq!(
        lookup("offset_last_past_end", 99990),
        (TaxLookup::Unnamed, vec![MSG_OFFSET_ERROR.to_string()])
    );
    // A database NCBI treats as missing finds no taxid.
    assert_eq!(
        lookup("bad_magic", 9606),
        (TaxLookup::NotInIndex, Vec::new())
    );
    assert_eq!(lookup("no_bti", 9606), (TaxLookup::NotInIndex, Vec::new()));
}

fn names_of_lazy(lazy: &LazyTaxDb, taxid: i32) -> (Option<TaxNames>, Vec<String>) {
    let mut messages = Vec::new();
    let names = lazy.names(taxid, &mut messages).unwrap();
    (names, messages)
}

fn four(names: &TaxNames) -> [&[u8]; 4] {
    [
        &names.scientific_name,
        &names.common_name,
        &names.blast_name,
        names
            .super_kingdom
            .as_deref()
            .expect("a kingdom inside the mapping"),
    ]
}

#[test]
fn good_database_finds_every_entry_and_nothing_else() {
    let db = open_fixture("good").unwrap();
    assert!(!db.is_missing());
    assert!(db.open_messages().is_empty());
    let expected: [(i32, [&[u8]; 4]); 7] = [
        (
            562,
            [b"Escherichia coli", b"", b"enterobacteria", b"Bacteria"],
        ),
        (7777, [b"Quote \"x\", y; z", b"a,b", b"c\"d", b"K;L"]),
        (9606, [b"Homo sapiens", b"human", b"primates", b"Eukaryota"]),
        (
            10090,
            [b"Mus musculus", b"house mouse", b"rodents", b"Eukaryota"],
        ),
        (
            12908,
            [
                b"unclassified sequences",
                b"",
                b"unclassified",
                b"unclassified",
            ],
        ),
        (32630, [b"synthetic construct", b"", b"-", b"-"]),
        // The last entry: the record ends at the end of taxdb.btd.
        (
            99990,
            [
                b"Last entry sp.",
                b"last common",
                b"lastblast",
                b"Eukaryota",
            ],
        ),
    ];
    for (taxid, fields) in expected {
        let (names, messages) = names_of(&db, taxid);
        let names = names.unwrap_or_else(|| panic!("taxid {taxid}"));
        assert_eq!(names.taxid, taxid);
        assert_eq!(four(&names), fields, "taxid {taxid}");
        assert!(messages.is_empty());
    }
    // Between two entries, below the first, above the last, zero, negative.
    for taxid in [
        777,
        8000,
        50000,
        3,
        561,
        99991,
        100000,
        0,
        -1,
        i32::MIN,
        i32::MAX,
    ] {
        assert_eq!(names_of(&db, taxid).0, None, "taxid {taxid}");
    }
}

#[test]
fn every_position_of_every_size_is_found() {
    // The search loop for index sizes 1..=40, every entry found, every gap not found.
    for n in 1..=40i32 {
        let rows: Vec<(i32, Vec<u8>)> = (0..n)
            .map(|i| (10 * i + 5, format!("s{i}\tc\tb\tk{i}").into_bytes()))
            .collect();
        let borrowed: Vec<(i32, &[u8])> = rows.iter().map(|(t, r)| (*t, r.as_slice())).collect();
        let (index, data) = build(&borrowed);
        let (_dir, db) = open_built(&index, &data);
        for i in 0..n {
            let names = names_of(&db, 10 * i + 5)
                .0
                .unwrap_or_else(|| panic!("n={n} i={i}"));
            assert_eq!(names.scientific_name, format!("s{i}").into_bytes());
            assert_eq!(names.super_kingdom, Ok(format!("k{i}").into_bytes()));
            assert_eq!(names_of(&db, 10 * i + 6).0, None);
            assert_eq!(names_of(&db, 10 * i + 4).0, None);
        }
    }
}

#[test]
fn missing_and_unusable_databases_find_nothing() {
    // (variant, messages written when it is opened)
    let cases: [(&str, &[&str]); 11] = [
        ("no_bti", &[]),
        ("no_btd", &[]),
        ("index_short", &[]),
        ("index_header_only", &[]),
        ("bad_magic", &[MSG_WRONG_MAGIC]),
        ("bad_magic_le", &[MSG_WRONG_MAGIC]),
        (
            "count_larger",
            &["SeqDB: Taxid metadata indicates (8) entries but file has room for (7)."],
        ),
        (
            "count_smaller",
            &["SeqDB: Taxid metadata indicates (6) entries but file has room for (7)."],
        ),
        (
            "count_zero",
            &["SeqDB: Taxid metadata indicates (0) entries but file has room for (7)."],
        ),
        (
            "index_trailing_8",
            &["SeqDB: Taxid metadata indicates (7) entries but file has room for (8)."],
        ),
        ("nonexistent", &[]),
    ];
    for (variant, messages) in cases {
        let db = open_fixture(variant).unwrap();
        assert!(db.is_missing(), "{variant}");
        assert_eq!(db.open_messages(), messages, "{variant}");
        assert_eq!(names_of(&db, 9606), (None, vec![]), "{variant}");
    }
}

#[test]
fn index_sizes_at_the_edges() {
    // Exactly header + one entry is valid; seven stray bytes are ignored (integer division).
    let db = open_fixture("index_one_entry").unwrap();
    assert!(!db.is_missing());
    assert!(names_of(&db, 9606).0.is_some());
    assert_eq!(names_of(&db, 9607).0, None);
    let db = open_fixture("index_trailing").unwrap();
    assert!(!db.is_missing());
    assert!(names_of(&db, 99990).0.is_some());
}

#[test]
fn empty_data_file_is_refused_only_where_a_record_is_read() {
    let db = open_fixture("btd_empty").unwrap();
    assert!(!db.is_missing());
    let mut messages = Vec::new();
    // NCBI: segmentation fault.
    assert!(matches!(
        db.names(9606, &mut messages),
        Err(SeqDbError::Unsupported { .. })
    ));
    // NCBI: not found; and the offset error of the last entry.
    assert_eq!(names_of(&db, 777), (None, vec![]));
    assert_eq!(names_of(&db, 100000), (None, vec![]));
    assert_eq!(
        names_of(&db, 99990),
        (None, vec![MSG_OFFSET_ERROR.to_string()])
    );
}

/// `std::env::join_paths`, which WASI lacks (it returns an error there): the directories join
/// with ":" as `path_splitter` splits them off Windows.
fn joined_paths(dirs: [&std::path::Path; 2]) -> std::ffi::OsString {
    #[cfg(not(target_os = "wasi"))]
    {
        std::env::join_paths(dirs).unwrap()
    }
    #[cfg(target_os = "wasi")]
    {
        let mut joined = dirs[0].as_os_str().to_os_string();
        joined.push(":");
        joined.push(dirs[1]);
        joined
    }
}

#[test]
fn the_index_path_decides_and_the_data_file_is_its_sibling() {
    // The first directory with a taxdb.bti wins, even when its taxdb.btd is missing.
    let good = fixture_dir().join("good");
    let nobtd = fixture_dir().join("no_btd");
    let joined = joined_paths([&nobtd, &good]);
    let cwd = TempDir::unique("taxdb-cwd");
    let search = SearchPath::new(&cwd, Some(joined.as_os_str()), None);
    assert!(TaxDb::open(&search).unwrap().is_missing());
    let joined = joined_paths([&fixture_dir().join("no_bti"), &good]);
    let search = SearchPath::new(&cwd, Some(joined.as_os_str()), None);
    assert!(!TaxDb::open(&search).unwrap().is_missing());
    // The current directory comes before BLASTDB.
    let search = SearchPath::new(&good, Some(nobtd.as_os_str()), None);
    assert!(!TaxDb::open(&search).unwrap().is_missing());
    // The registry value comes after BLASTDB.
    let search = SearchPath::new(
        &cwd,
        Some(nobtd.as_os_str()),
        good.to_str().map(str::as_bytes),
    );
    assert!(TaxDb::open(&search).unwrap().is_missing());
    let search = SearchPath::new(&cwd, None, good.to_str().map(str::as_bytes));
    assert!(!TaxDb::open(&search).unwrap().is_missing());
    // No directory at all.
    let search = SearchPath::new(&cwd, None, None);
    assert!(TaxDb::open(&search).unwrap().is_missing());
}

#[test]
fn records_with_other_than_four_fields() {
    let rows: [(i32, &[u8]); 9] = [
        (1, b"a\tb\tc\td"),
        (2, b"a\tb\tc"),       // two tabs: no kingdom
        (3, b"a\tb\tc\t"),     // empty kingdom
        (4, b"a"),             // no tab
        (5, b"a\tb\tc\td\te"), // kingdom takes the rest
        (6, b"\t\t\tx"),       // all other fields empty
        (7, b"a\tb"),
        (8, b"a\tb\tc\td\n"), // a terminator belongs to the kingdom
        (9, b"a\tb\tc\td\0"),
    ];
    let (index, data) = build(&rows);
    let (_dir, db) = open_built(&index, &data);
    assert_eq!(four(&names_of(&db, 1).0.unwrap()), [b"a", b"b", b"c", b"d"]);
    assert_eq!(names_of(&db, 2).0, None);
    assert_eq!(names_of(&db, 3).0, None);
    assert_eq!(names_of(&db, 4).0, None);
    assert_eq!(
        four(&names_of(&db, 5).0.unwrap()),
        [b"a".as_slice(), b"b", b"c", b"d\te"]
    );
    assert_eq!(
        four(&names_of(&db, 6).0.unwrap()),
        [b"".as_slice(), b"", b"", b"x"]
    );
    assert_eq!(names_of(&db, 7).0, None);
    assert_eq!(four(&names_of(&db, 8).0.unwrap())[3], b"d\n");
    // The last entry: the record ends at the end of the file.
    assert_eq!(four(&names_of(&db, 9).0.unwrap())[3], b"d\0");
}

#[test]
fn last_entry_record_shapes() {
    for (record, found) in [
        (&b"a\tb\tc\td"[..], true),
        (b"a\tb\tc\t", false),
        (b"a\tb\tc", false),
        (b"a", false),
        (b"", false),
    ] {
        let (index, data) = build(&[(1, b"x\ty\tz\tw"), (2, record)]);
        let (_dir, db) = open_built(&index, &data);
        assert_eq!(names_of(&db, 2).0.is_some(), found, "{record:?}");
        assert!(names_of(&db, 1).0.is_some());
    }
}

#[test]
fn last_offset_past_the_end_is_an_error_message_and_not_found() {
    let (mut index, data) = build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th")]);
    let at = 24 + 8 + 4; // offset of the second entry
    index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 10));
    let (_dir, db) = open_built(&index, &data);
    let (names, messages) = names_of(&db, 2);
    assert_eq!(names, None);
    assert_eq!(messages, [MSG_OFFSET_ERROR]);
    // The message comes at every lookup.
    assert_eq!(names_of(&db, 2).1, [MSG_OFFSET_ERROR]);
    // An offset equal to the size is an empty record: not found, no message.
    index[at..at + 4].copy_from_slice(&be(data.len() as u32));
    let (_dir, db) = open_built(&index, &data);
    assert_eq!(names_of(&db, 2), (None, vec![]));
}

#[test]
fn the_taxid_of_the_offset_error_is_the_last_entry_s_where_its_lookup_writes_it() {
    // Every variant: the taxid for which `blastdbcmd -tax_info` of NCBI BLAST+ 2.17.0 wrote
    // the message (`expected.tsv`), none where it wrote it for no taxid.
    let text = std::fs::read_to_string(fixture_dir().join("expected.tsv")).unwrap();
    let mut variants: Vec<String> = std::fs::read_dir(fixture_dir())
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_dir())
        .map(|path| path.file_name().unwrap().to_string_lossy().to_string())
        .collect();
    variants.sort();
    let mut with_message = 0;
    for variant in &variants {
        let Ok(db) = open_fixture(variant) else {
            continue;
        };
        let ncbi: Vec<i32> = text
            .lines()
            .map(|line| line.splitn(4, '\t').collect::<Vec<&str>>())
            .filter(|cols| {
                cols[0] == variant.as_str()
                    && cols.get(3).is_some_and(|m| m.contains(MSG_OFFSET_ERROR))
            })
            .map(|cols| cols[1].parse().unwrap())
            .collect();
        assert!(ncbi.len() <= 1, "{variant}");
        assert_eq!(
            db.offset_error_taxid().unwrap(),
            ncbi.first().copied(),
            "{variant}"
        );
        with_message += ncbi.len();
    }
    // offset_last_past_end, btd_truncated and btd_empty.
    assert_eq!(with_message, 3);
    // The last entry (50) begins past the end, but the search for 50 ends on entry 2 (90 >
    // 50, then 20 < 50): no lookup writes the message.
    let rows: [(i32, &[u8]); 5] = [
        (10, b"a\tb\tc\td"),
        (20, b"a\tb\tc\td"),
        (90, b"a\tb\tc\td"),
        (30, b"a\tb\tc\td"),
        (50, b"a\tb\tc\td"),
    ];
    let (mut index, data) = build(&rows);
    let at = 24 + 8 * 4 + 4;
    index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 10));
    let (_dir, db) = open_built(&index, &data);
    assert_eq!(names_of(&db, 50), (None, vec![]));
    assert_eq!(db.offset_error_taxid().unwrap(), None);
    // Sorted, the same entry writes it.
    let rows: [(i32, &[u8]); 3] = [(10, b"a\tb\tc\td"), (20, b"a\tb\tc\td"), (50, b"x")];
    let (mut index, data) = build(&rows);
    let at = 24 + 8 * 2 + 4;
    index[at..at + 4].copy_from_slice(&be(data.len() as u32 + 1));
    let (_dir, db) = open_built(&index, &data);
    assert_eq!(db.offset_error_taxid().unwrap(), Some(50));
}

#[test]
fn record_ending_before_it_begins_is_not_found() {
    // Entry 1 starts after entry 2: the size of 1 is negative.
    let (mut index, data) = build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
    let off1 = 24 + 4;
    let off2 = 24 + 8 + 4;
    index[off1..off1 + 4].copy_from_slice(&be(8));
    index[off2..off2 + 4].copy_from_slice(&be(0));
    let (_dir, db) = open_built(&index, &data);
    assert_eq!(names_of(&db, 1), (None, vec![]));
    // Entry 2 is now [0, 16): both records, the kingdom takes the rest.
    let names = names_of(&db, 2).0.unwrap();
    assert_eq!(names.scientific_name, b"a");
    assert_eq!(
        names.super_kingdom,
        Ok(b"d"
            .to_vec()
            .into_iter()
            .chain(*b"e\tf\tg\th")
            .collect::<Vec<u8>>())
    );
}

#[test]
fn record_reaching_past_the_end_reads_page_padding_and_then_memory_past_the_mapping() {
    let (mut index, data) = build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
    let off2 = 24 + 8 + 4;
    let off3 = 24 + 16 + 4;
    // Entry 2 ends 10 bytes after the file: zero padding belongs to its kingdom.
    index[off3..off3 + 4].copy_from_slice(&be(data.len() as u32 + 10));
    let (_dir, db) = open_built(&index, &data);
    let names = names_of(&db, 2).0.unwrap();
    let mut kingdom = b"h".to_vec();
    kingdom.extend(*b"i\tj\tk\tl");
    kingdom.extend([0u8; 10]);
    assert_eq!(names.super_kingdom, Ok(kingdom));
    // Entry 3 (the last) begins after the end: NCBI's message.
    assert_eq!(names_of(&db, 3), (None, vec![MSG_OFFSET_ERROR.to_string()]));
    // Past the page: the three names before the third tab inside the mapping, a kingdom
    // NCBI reads past it. Entry 1 is fine.
    index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 1));
    let (_dir, db) = open_built(&index, &data);
    let past = PastMapping {
        path: db.files.as_ref().unwrap().data.path().to_path_buf(),
        taxid: 2,
        begin: 7,
        end: MAP_PAGE + 1,
        len: data.len() as u64,
    };
    let mut messages = Vec::new();
    assert_eq!(
        db.lookup(2, &mut messages).unwrap(),
        TaxLookup::Named(TaxNames {
            taxid: 2,
            scientific_name: b"e".to_vec(),
            common_name: b"f".to_vec(),
            blast_name: b"g".to_vec(),
            super_kingdom: Err(past.clone()),
        })
    );
    assert!(messages.is_empty());
    assert!(names_of(&db, 1).0.is_some());
    // A record past the mapping without its third tab inside it: no defined result, which
    // only a caller that writes the names refuses.
    index[off2..off2 + 4].copy_from_slice(&be(data.len() as u32 - 2));
    let (_dir, db) = open_built(&index, &data);
    let past = PastMapping {
        path: db.files.as_ref().unwrap().data.path().to_path_buf(),
        begin: data.len() as u64 - 2,
        ..past
    };
    assert_eq!(
        db.lookup(2, &mut messages).unwrap(),
        TaxLookup::PastMapping(past.clone())
    );
    assert!(matches!(
        db.names(2, &mut messages),
        Err(SeqDbError::Unsupported { .. })
    ));
    assert!(messages.is_empty());
    // A record that begins past the mapping.
    index[off2..off2 + 4].copy_from_slice(&be(MAP_PAGE as u32));
    index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 9));
    let (_dir, db) = open_built(&index, &data);
    assert!(matches!(
        db.lookup(2, &mut messages).unwrap(),
        TaxLookup::PastMapping(PastMapping {
            begin: 4096,
            end: 4105,
            ..
        })
    ));
    // Neither faults for certain, nor does the pass over the index say so; an empty
    // taxdb.btd does.
    assert!(!db.any_lookup_faults().unwrap());
    assert!(!db.maps_null());
    let empty = open_fixture("btd_empty").unwrap();
    assert!(empty.any_lookup_faults().unwrap());
    assert!(empty.maps_null());
}

#[test]
fn unsorted_index_follows_the_search_loop() {
    // A port of the loop gives these results for an unsorted index; the same loop in NCBI
    // finds an entry only where its probes happen to land.
    let rows: [(i32, &[u8]); 5] = [
        (10, b"a\tb\tc\td10"),
        (50, b"a\tb\tc\td50"),
        (20, b"a\tb\tc\td20"),
        (40, b"a\tb\tc\td40"),
        (60, b"a\tb\tc\td60"),
    ];
    let (index, data) = build(&rows);
    let (_dir, db) = open_built(&index, &data);
    let found: Vec<i32> = [10, 20, 40, 50, 60]
        .into_iter()
        .filter(|t| names_of(&db, *t).0.is_some())
        .collect();
    // 50: new=2 (20) -> low=2, new=3 (40) -> low=3, new=3 again: 50 > 40 so new=4 (60), which is not 50.
    assert_eq!(found, [10, 20, 40, 60]);
    // High-bit taxids are negative Int4 values, sorted after 99990 in the file: out of range.
    let db = open_fixture("taxid_high_bit").unwrap();
    assert_eq!(names_of(&db, 9606).0, None);
    assert_eq!(names_of(&db, i32::MIN + 5).0, None);
}

#[test]
fn names_are_bytes() {
    let db = open_fixture("high_bytes_utf8").unwrap();
    assert_eq!(
        names_of(&db, 9606).0.unwrap().scientific_name,
        "Caf\u{e9} sapiens".as_bytes()
    );
    let db = open_fixture("high_bytes_latin1").unwrap();
    let names = names_of(&db, 10090).0.unwrap();
    assert_eq!(names.common_name, b"mus\xe9");
    assert_eq!(names.super_kingdom.as_deref(), Ok(&b"Euk\xffaryota"[..]));
}

#[test]
fn lazy_handle_opens_once_and_hands_over_the_open_messages_once() {
    let dir = fixture_dir().join("bad_magic");
    let searching = path_of(&dir);
    let lazy = LazyTaxDb::new(searching.search.clone());
    assert!(!lazy.is_open());
    let mut messages = Vec::new();
    assert_eq!(lazy.names(9606, &mut messages).unwrap(), None);
    assert!(lazy.is_open());
    assert_eq!(messages, [MSG_WRONG_MAGIC]);
    assert_eq!(lazy.names(9606, &mut messages).unwrap(), None);
    assert_eq!(messages, [MSG_WRONG_MAGIC]);
    // A missing database: opened at the first name, silent, finds nothing.
    let searching = path_of(&fixture_dir().join("nonexistent"));
    let lazy = LazyTaxDb::new(searching.search.clone());
    assert!(!lazy.is_open());
    for _ in 0..2 {
        assert_eq!(names_of_lazy(&lazy, 9606), (None, vec![]));
    }
    assert!(lazy.is_open());
}

#[test]
fn err_post_line_is_what_blast_prints() {
    assert_eq!(
        err_post_line("blastp", MSG_WRONG_MAGIC),
        "Error: [blastp] Error: Tax database file has wrong magic number."
    );
}

/// `make_expected.py`'s escaping: bytes outside 0x20..0x7e and the backslash as `\xNN`.
fn escape(bytes: &[u8]) -> String {
    bytes
        .iter()
        .map(|&b| {
            if (0x20..0x7f).contains(&b) && b != b'\\' {
                (b as char).to_string()
            } else {
                format!("\\x{b:02x}")
            }
        })
        .collect()
}

/// Every (variant, taxid) of `expected.tsv`: what `blastdbcmd -tax_info` of NCBI BLAST+ 2.17.0
/// printed (`%T|%S|%L|%B|%K`, escaped), `N/A` for a taxid it did not find, `CRASH` for a
/// signal;
/// then the lines it wrote to standard error.
#[test]
fn matches_blastdbcmd_tax_info_on_every_variant() {
    let text = std::fs::read_to_string(fixture_dir().join("expected.tsv")).unwrap();
    let mut checked = 0;
    for line in text.lines().filter(|l| !l.is_empty()) {
        let cols: Vec<&str> = line.splitn(4, '\t').collect();
        let (variant, taxid, expect) = (cols[0], cols[1].parse::<i32>().unwrap(), cols[2]);
        let ncbi_messages: Vec<&str> = cols[3].split('\u{1f}').filter(|m| !m.is_empty()).collect();
        let what = format!("{variant} taxid {taxid}");
        let opened = open_fixture(variant);
        let mut messages: Vec<String> = Vec::new();
        let result = match &opened {
            Ok(db) => {
                messages.extend(db.open_messages().iter().cloned());
                db.names(taxid, &mut messages)
            }
            Err(error) => Err(SeqDbError::Range(error.to_string())),
        };
        match (expect, &result) {
            ("CRASH", Err(_)) => continue,
            ("N/A", Ok(None)) => {}
            (expect, Ok(Some(names))) => {
                let mut raw = format!("{}|", names.taxid).into_bytes();
                for (i, field) in four(names).iter().enumerate() {
                    if i > 0 {
                        raw.push(b'|');
                    }
                    raw.extend(*field);
                }
                assert_eq!(escape(&raw), expect, "{what}");
            }
            (expect, other) => panic!("{what}: NCBI {expect:?}, got {other:?}"),
        }
        assert_eq!(messages, ncbi_messages, "{what}: messages");
        checked += 1;
    }
    assert!(checked > 300, "{checked}");
}

/// What [`TaxDb::any_lookup_may_fail`] must agree with: the lookup of every taxid of the
/// index, failing (`Err`) or giving a part of its record from past NCBI's mapping of
/// `taxdb.btd` (`TaxLookup::PastMapping`, a kingdom `Err`); `names` then refuses the first.
fn brute_force_may_fail(db: &TaxDb) -> bool {
    let Some(files) = &db.files else {
        return false;
    };
    (0..files.count).any(|i| {
        let taxid = TaxDb::entry(files, i).unwrap().0;
        !matches!(
            db.lookup(taxid, &mut Vec::new()),
            Ok(TaxLookup::NotInIndex
                | TaxLookup::Unnamed
                | TaxLookup::Named(TaxNames {
                    super_kingdom: Ok(_),
                    ..
                }))
        )
    })
}

/// Whether the index of `db` is in the layout the pass relies on (strictly increasing taxids,
/// non-decreasing offsets).
fn regular_layout(db: &TaxDb) -> bool {
    let Some(files) = &db.files else {
        return true;
    };
    let entries: Vec<(i32, u32)> = (0..files.count)
        .map(|i| TaxDb::entry(files, i).unwrap())
        .collect();
    entries
        .windows(2)
        .all(|pair| pair[0].0 < pair[1].0 && pair[0].1 <= pair[1].1)
}

/// The pass over the index answers as the lookups do on every fixture variant: never `false`
/// where the lookup of a taxid of the index, of a taxid next to one, or of the extremes fails;
/// and exactly the lookups' answer on the taxids of the index when the index is sorted. Of the
/// variants, `btd_empty` fails (every record past the mapped end); a record cut by the end of
/// the file inside its last page does not (`btd_truncated`, `offset_next_past_end`).
#[test]
fn the_pass_over_the_index_agrees_with_the_failing_lookups_on_every_variant() {
    let (mut failing, mut clean) = (Vec::new(), 0);
    let mut dirs: Vec<_> = std::fs::read_dir(fixture_dir())
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_dir())
        .collect();
    dirs.sort();
    for dir in dirs {
        let variant = dir.file_name().unwrap().to_string_lossy().to_string();
        let Ok(db) = open_fixture(&variant) else {
            continue;
        };
        let pass = db.any_lookup_faults().unwrap();
        let fails = |taxid: i32| db.lookup(taxid, &mut Vec::new()).is_err();
        let indexed: Vec<i32> = db.files.as_ref().map_or(Vec::new(), |files| {
            (0..files.count)
                .map(|i| TaxDb::entry(files, i).unwrap().0)
                .collect()
        });
        let near = indexed
            .iter()
            .flat_map(|&taxid| [taxid.wrapping_sub(1), taxid.wrapping_add(1)])
            .chain([i32::MIN, -1, 0, 1, i32::MAX]);
        assert!(
            pass || !indexed.iter().copied().chain(near).any(fails),
            "{variant}: the pass says no, a lookup fails"
        );
        if regular_layout(&db) {
            assert_eq!(pass, indexed.iter().copied().any(fails), "{variant}");
        }
        if pass {
            failing.push(variant);
        } else {
            clean += 1;
        }
    }
    assert!(
        failing.iter().any(|variant| variant == "btd_empty") && clean >= 20,
        "{failing:?} {clean}"
    );
}

/// The one pass over the files answers as the lookup of every taxid does on every fixture
/// variant: never `false` where a lookup refuses, and exactly the same when the index has the
/// layout the pass relies on.
#[test]
fn one_pass_over_the_files_agrees_with_the_lookups_on_every_variant() {
    let mut refusing = 0;
    let mut clean = 0;
    let mut dirs: Vec<_> = std::fs::read_dir(fixture_dir())
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_dir())
        .collect();
    dirs.sort();
    for dir in dirs {
        let variant = dir.file_name().unwrap().to_string_lossy().to_string();
        let Ok(db) = open_fixture(&variant) else {
            continue;
        };
        let brute = brute_force_may_fail(&db);
        let walk = db.any_lookup_may_fail().unwrap();
        assert!(
            walk || !brute,
            "{variant}: the pass says no, a lookup refuses"
        );
        if regular_layout(&db) {
            assert_eq!(walk, brute, "{variant}");
        }
        if brute {
            refusing += 1;
        } else {
            clean += 1;
        }
    }
    // `btd_empty` refuses (an empty taxdb.btd maps no page); a scientific name of any bytes is
    // no refusal (the two `high_bytes_*` variants are clean); the test of the index's blocks
    // (`the_pass_counts_entries_and_reads_past_one_index_block`) refuses past the first.
    assert!(refusing >= 1 && clean >= 20, "{refusing} {clean}");
}

#[test]
fn the_pass_counts_entries_and_reads_past_one_index_block() {
    assert_eq!(open_fixture("good").unwrap().entry_count(), 7);
    assert_eq!(open_fixture("no_btd").unwrap().entry_count(), 0);
    assert_eq!(open_fixture("bad_magic").unwrap().entry_count(), 0);
    // 70,000 entries (more than one block of the index); in the damaged copy the record of an
    // entry near the end reaches past NCBI's mapping of taxdb.btd (the two entries after it
    // begin past the file's last page), whose lookup gives its kingdom from past the mapping.
    let n = 70_000u32;
    let mut data = Vec::new();
    let mut offsets = Vec::new();
    for i in 0..n {
        offsets.push(data.len() as u32);
        data.extend(format!("name {i}\tc\tb\tK\n").as_bytes());
    }
    let past = (data.len().div_ceil(MAP_PAGE as usize) * MAP_PAGE as usize + 4096) as u32;
    let mut damaged = offsets.clone();
    damaged[n as usize - 2] = past;
    damaged[n as usize - 1] = past;
    for (tag, offsets, refuses) in [("walk", &offsets, false), ("walkpast", &damaged, true)] {
        let dir = TempDir::unique(&format!("taxdb-{tag}"));
        let mut index = Vec::new();
        index.extend(MAGIC.to_be_bytes());
        index.extend(n.to_be_bytes());
        index.extend([0u8; 16]);
        for (i, offset) in offsets.iter().enumerate() {
            index.extend((i as i32 + 1).to_be_bytes());
            index.extend(offset.to_be_bytes());
        }
        std::fs::write(dir.join("taxdb.bti"), &index).unwrap();
        std::fs::write(dir.join("taxdb.btd"), &data).unwrap();
        let db = TaxDb::open(&path_of(&dir)).unwrap();
        assert_eq!(db.entry_count(), u64::from(n));
        assert_eq!(db.any_lookup_may_fail().unwrap(), refuses, "{tag}");
        assert_eq!(brute_force_may_fail(&db), refuses, "{tag}");
    }
}

/// `TaxCache::names` of `taxid` through `db`, with the number of lookups it made.
fn cached(
    cache: &TaxCache,
    db: &TaxDb,
    taxid: i32,
) -> (
    Result<std::sync::Arc<SourceNames>, SeqDbError>,
    Vec<String>,
    u32,
) {
    let mut messages = Vec::new();
    let mut asked = 0;
    let names = cache.names(taxid, &mut messages, |messages| {
        asked += 1;
        db.names(taxid, messages)
    });
    (names, messages, asked)
}

// NCBI's `CSeqDBIntCache` (c++/include/objtools/blast/seqdb_reader/impl/seqdbgeneral.hpp:1184-1196): a slot per `taxid & 255` that holds the
// last taxid looked up there. LOSAT keeps only the answers that a lookup gives again without
// a text: those that wrote one, and the refusals, are asked again.
#[test]
fn the_cache_keeps_answers_without_texts_in_a_slot_per_taxid() {
    // 257 shares the slot of 1; 300 is in no entry (found nothing); 600 is the last entry,
    // whose offset is past the end (the offset error at each lookup).
    let (mut index, data) = build(&[
        (1, b"Bacillus\tb\tc\td"),
        (257, b"\tb\tc\td"),
        (600, b"e\tf\tg\th"),
    ]);
    let last = 24 + 16 + 4;
    index[last..last + 4].copy_from_slice(&be(data.len() as u32 + 10));
    let (_dir, db) = open_built(&index, &data);
    let cache = TaxCache::default();
    let (first, messages, asked) = cached(&cache, &db, 1);
    let first = first.unwrap();
    assert_eq!((messages.len(), asked), (0, 1));
    assert_eq!(
        first.taxname().map(|name| &name[..]),
        Some(&b"Bacillus"[..])
    );
    let (again, _, asked) = cached(&cache, &db, 1);
    assert!(std::sync::Arc::ptr_eq(&first, &again.unwrap()));
    assert_eq!(asked, 0);
    // An empty scientific name is no `m_Taxname`; the names are still there.
    let (other, _, asked) = cached(&cache, &db, 257);
    let other = other.unwrap();
    assert_eq!(
        (other.taxname(), other.names().is_some(), asked),
        (None, true, 1)
    );
    // 257 took the slot: 1 is asked again.
    assert_eq!(cached(&cache, &db, 1).2, 1);
    // Nothing found is kept.
    let (none, _, asked) = cached(&cache, &db, 300);
    assert!(none.unwrap().names().is_none());
    assert_eq!(asked, 1);
    assert_eq!(cached(&cache, &db, 300).2, 0);
    // The offset error comes at each lookup: never kept.
    for _ in 0..2 {
        let (names, messages, asked) = cached(&cache, &db, 600);
        assert!(names.unwrap().names().is_none());
        assert_eq!((messages, asked), (vec![MSG_OFFSET_ERROR.to_string()], 1));
    }
    // A refused lookup (a record that lies past NCBI's mapping of taxdb.btd, whose names are
    // not defined) is asked again.
    let (mut index, data) = build(&[(1, b"a\tb\tc\td"), (2, b"e\tf\tg\th"), (3, b"i\tj\tk\tl")]);
    let off2 = 24 + 8 + 4;
    let off3 = 24 + 16 + 4;
    index[off2..off2 + 4].copy_from_slice(&be(MAP_PAGE as u32));
    index[off3..off3 + 4].copy_from_slice(&be(MAP_PAGE as u32 + 9));
    let (_dir, db) = open_built(&index, &data);
    let cache = TaxCache::default();
    for _ in 0..2 {
        let (names, _, asked) = cached(&cache, &db, 2);
        assert!(matches!(names, Err(SeqDbError::Unsupported { .. })));
        assert_eq!(asked, 1);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:154-161
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:199
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbtax.cpp:217-220
/// A `taxdb.bti` or `taxdb.btd` that does not open (mode 000) is `SeqDbError::TaxdbThrows`, with
/// `CMemoryFile`'s exception, and the lazy handle gives it again at every lookup, as NCBI's
/// taxonomy file object throws at every lookup (line WQA D18b, 2026-11-21). An index too short
/// for one entry is missing before it is opened (NCBI takes its length first), and a data file
/// is not opened where the index has a wrong magic number. Skipped for a user that reads every
/// file.
#[cfg(unix)]
#[test]
fn a_taxonomy_file_that_does_not_open_throws_at_every_lookup() {
    use std::os::unix::fs::PermissionsExt;
    let mode = |path: &Path, mode: u32| {
        std::fs::set_permissions(path, std::fs::Permissions::from_mode(mode)).unwrap()
    };
    let good = fixture_dir().join("good");
    for file in ["taxdb.bti", "taxdb.btd"] {
        let dir = TempDir::unique("taxdb-unopened");
        for name in ["taxdb.bti", "taxdb.btd"] {
            std::fs::copy(good.join(name), dir.join(name)).unwrap();
        }
        let path = dir.join(file);
        mode(&path, 0o000);
        if std::fs::File::open(&path).is_ok() {
            // A user that reads every file (the superuser).
            mode(&path, 0o644);
            return;
        }
        let search = path_of(&dir);
        // The file in the directory of the search path, `BLASTDB` with its links followed
        // (seqdbcommon.cpp:378-392, `crate::utils::test_dir::ncbi_dir`).
        let named = crate::utils::test_dir::ncbi_dir(&dir).join(file);
        let expect = |error: SeqDbError| match error {
            SeqDbError::TaxdbThrows {
                path: failed,
                reason,
                exception,
            } => {
                assert_eq!(failed, named, "{file}");
                assert_eq!(reason, "Permission denied", "{file}");
                let message = format!(
                    "CMemoryFile: Cannot memory map file '{}':Permission denied",
                    named.display()
                );
                assert_eq!(exception.get_msg(), message.as_bytes(), "{file}");
                let error = SeqDbError::TaxdbThrows {
                    path: failed,
                    reason,
                    exception,
                };
                // Not NCBI's database error: the exception, exit status 255 where nothing
                // catches it; blastdbcmd's request writes its message.
                assert_eq!(error.ncbi_message(), None, "{file}");
                assert!(error.ncbi_exception().is_some(), "{file}");
                assert_eq!(
                    error.ncbi_request_message(),
                    Some(message.into_bytes()),
                    "{file}"
                );
                assert_eq!(
                    error.to_string(),
                    format!("{} cannot be opened: Permission denied", named.display())
                );
            }
            other => panic!("{file}: {other:?}"),
        };
        match TaxDb::open(&search) {
            Err(error) => expect(error),
            Ok(_) => panic!("{file}: opened"),
        }
        let lazy = LazyTaxDb::new(search.search.clone());
        for _ in 0..2 {
            let mut messages = Vec::new();
            expect(lazy.names(9606, &mut messages).expect_err(file));
            assert!(messages.is_empty(), "{file}");
        }
        mode(&path, 0o644);
    }
    // An index too short for one entry: missing, though it does not open.
    let dir = TempDir::unique("taxdb-short-unopened");
    std::fs::write(dir.join("taxdb.bti"), [0u8; 8]).unwrap();
    std::fs::copy(good.join("taxdb.btd"), dir.join("taxdb.btd")).unwrap();
    mode(&dir.join("taxdb.bti"), 0o000);
    let db = TaxDb::open(&path_of(&dir)).unwrap();
    assert!(db.is_missing());
    mode(&dir.join("taxdb.bti"), 0o644);
    // An index of a wrong magic number: missing with NCBI's text, the data file never opened.
    let dir = TempDir::unique("taxdb-bad-magic-unopened");
    let bad = fixture_dir().join("bad_magic");
    std::fs::copy(bad.join("taxdb.bti"), dir.join("taxdb.bti")).unwrap();
    std::fs::copy(bad.join("taxdb.btd"), dir.join("taxdb.btd")).unwrap();
    mode(&dir.join("taxdb.btd"), 0o000);
    let db = TaxDb::open(&path_of(&dir)).unwrap();
    assert!(db.is_missing());
    assert_eq!(db.open_messages(), [MSG_WRONG_MAGIC]);
    mode(&dir.join("taxdb.btd"), 0o644);
}

#[test]
fn read_paths_are_the_places_that_the_run_looks_at_up_to_the_one_found() {
    // Line WQA D35 (WQA DECISIONS 62): the taxonomy files that a run checks or opens, as the
    // native open trace shows them: `taxdb.bti` in each directory of the search path up to the
    // first that holds one, the `taxdb.btd` beside that one, and `taxonomy4blast.sqlite3` in the
    // same way for a taxonomy filter that is widened.
    use crate::seqdb::DbUse;
    let cwd = TempDir::unique("taxdb-read-cwd");
    let dir = TempDir::unique("taxdb-read-dir");
    let search = SearchPath::new(&cwd, Some(dir.as_os_str()), None);
    let bytes = |path: PathBuf| path.into_os_string().into_encoded_bytes();
    let at = |root: &Path, name: &str| bytes(root.join(name));
    let search_use = DbUse::of_program("blastp");
    let mut every_taxid = DbUse::of_program("blastp");
    every_taxid.note_every_taxid();
    let mut filtered = DbUse::of_program("blastp");
    filtered.note_option("taxids");
    let mut kept = filtered;
    kept.note_flag("no_taxid_expansion");

    // No taxonomy file anywhere: a search of a database looks in every directory, whatever its
    // report; a `-subject` search only for a report that looks every taxid up.
    let nowhere = vec![at(&cwd, "taxdb.bti"), at(&dir, "taxdb.bti")];
    assert_eq!(read_paths(&search, search_use, true), nowhere);
    assert_eq!(read_paths(&search, DbUse::TOOL, true), {
        let mut with_tree = nowhere.clone();
        with_tree.extend([at(&cwd, "taxonomy4blast.sqlite3"), at(&dir, "taxonomy4blast.sqlite3")]);
        with_tree
    });
    assert!(read_paths(&search, search_use, false).is_empty());
    assert_eq!(read_paths(&search, every_taxid, false), nowhere);

    // The database in `dir` (`BLASTDB`): the working directory's name is looked at first (a
    // `-logfile taxdb.bti` would be found there), then `dir`'s, and its data file is opened.
    std::fs::write(dir.join("taxdb.bti"), [0u8; 32]).unwrap();
    std::fs::write(dir.join("taxdb.btd"), [0u8; 32]).unwrap();
    let in_dir = vec![
        at(&cwd, "taxdb.bti"),
        at(&dir, "taxdb.bti"),
        at(&dir, "taxdb.btd"),
    ];
    assert_eq!(read_paths(&search, search_use, true), in_dir);
    assert_eq!(read_paths(&search, every_taxid, false), in_dir);

    // A widened taxonomy filter looks for `taxonomy4blast.sqlite3` the same way; with
    // `-no_taxid_expansion` it does not.
    std::fs::write(dir.join("taxonomy4blast.sqlite3"), [0u8; 32]).unwrap();
    let mut with_tree = in_dir.clone();
    with_tree.extend([
        at(&cwd, "taxonomy4blast.sqlite3"),
        at(&dir, "taxonomy4blast.sqlite3"),
    ]);
    assert_eq!(read_paths(&search, filtered, true), with_tree);
    assert_eq!(read_paths(&search, kept, true), in_dir);

    // A database in the working directory hides `dir`'s: the run looks no further, and a name
    // in `dir` is not one that it reads.
    std::fs::write(cwd.join("taxdb.bti"), [0u8; 32]).unwrap();
    assert_eq!(
        read_paths(&search, search_use, true),
        vec![at(&cwd, "taxdb.bti"), at(&cwd, "taxdb.btd")]
    );
}
