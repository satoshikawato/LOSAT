use super::*;
use std::collections::HashMap;
use std::path::PathBuf;

/// Directory of the fixtures (small databases made with Python's `sqlite3` from the
/// schema text of NCBI's real file) and of their `.ref` files (NCBI's query run by SQLite).
fn fx_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxsqlite")
}

fn fnv(ids: &[i32]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for i in ids {
        for b in format!("{},", i).bytes() {
            h = (h ^ u64::from(b)).wrapping_mul(0x100000001b3);
        }
    }
    h
}

fn load(name: &str) -> Vec<u8> {
    std::fs::read(fx_dir().join(name)).unwrap_or_else(|e| panic!("fixture {}: {}", name, e))
}

/// Compare every node of a fixture with the reference written by gen_fx.py.
fn check_fixture(name: &str) {
    let data = load(&format!("{}.sqlite3", name));
    let db = TaxonomySqlite::open(&data).unwrap_or_else(|e| panic!("{}: open: {}", name, e));
    let refs = std::fs::read_to_string(fx_dir().join(format!("{}.ref", name))).unwrap();
    let mut n = 0;
    let mut full: HashMap<i64, Vec<i32>> = HashMap::new();
    for line in refs.lines() {
        let f: Vec<&str> = line.split(' ').collect();
        if f[0] == "wide" {
            continue;
        }
        if f[0] == "full" {
            full.insert(
                f[1].parse().unwrap(),
                f[2..]
                    .iter()
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse().unwrap())
                    .collect(),
            );
            continue;
        }
        let t: i64 = f[0].parse().unwrap();
        // taxids beyond i32 cannot be asked through the i32 API
        if t > i64::from(i32::MAX) {
            continue;
        }
        if f[1] == "LOOP" {
            // taxid 1 is its own parent: SQLite's query never ends; see `root_in_a_parent_cycle`
            continue;
        }
        let mut got = db
            .descendants(t as i32)
            .unwrap_or_else(|e| panic!("{} taxid {}: {}", name, t, e));
        got.sort_unstable();
        let (cnt, h): (usize, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        assert_eq!((got.len(), fnv(&got)), (cnt, h), "{} taxid {}", name, t);
        if let Some(want) = full.get(&t) {
            assert_eq!(&got, want, "{} taxid {}", name, t);
        }
        n += 1;
    }
    assert!(n > 50, "{}: only {} nodes compared", name, n);
    // taxid <= 0 and absent taxids
    assert_eq!(db.descendants(0).unwrap(), Vec::<i32>::new());
    assert_eq!(db.descendants(-5).unwrap(), Vec::<i32>::new());
    assert_eq!(db.descendants(i32::MAX).unwrap(), Vec::<i32>::new());
    println!("{}: {} nodes equal to sqlite3", name, n);
}

#[test]
fn fixture_page_512() {
    check_fixture("p512");
}

#[test]
fn fixture_utf16() {
    check_fixture("utf16");
}

/// Taxid 1 is its own parent (as in NCBI's file). NCBI's query for it never ends; LOSAT
/// skips visited taxids (DI-8 DECISIONS 10), so the result is every other node reachable from 1.
#[test]
fn root_in_a_parent_cycle() {
    for name in ["p512", "utf16"] {
        let data = load(&format!("{name}.sqlite3"));
        let db = TaxonomySqlite::open(&data).unwrap();
        let refs = std::fs::read_to_string(fx_dir().join(format!("{name}.ref"))).unwrap();
        // The reference lists every node of the fixture (the "1 LOOP" line is the root); the
        // nodes that 1 does not reach are 2 (parent NULL), 3 (parent 0) and the two absent ids.
        let mut want: Vec<i32> = Vec::new();
        let mut listed = 0;
        for l in refs.lines() {
            let f: Vec<&str> = l.split(' ').collect();
            if f[0] == "wide" || f[0] == "full" {
                continue;
            }
            listed += 1;
            let t: i64 = f[0].parse().unwrap();
            if [1, 2, 3, 12345, 999_999_999].contains(&t) {
                continue;
            }
            // sqlite3_column_int keeps the low 32 bits
            want.push(t as i32);
        }
        want.sort_unstable();
        let mut got = db.descendants(1).unwrap();
        got.sort_unstable();
        assert_eq!(got.len(), want.len(), "{name}");
        assert_eq!(got, want, "{name}");
        assert_eq!(listed, want.len() + 5, "{name}");
        assert!(!got.contains(&1));
    }
}

/// Two nodes that are each other's parent (20 and 21), and the root: the walk ends and
/// returns what NCBI's query would return before it started to repeat itself.
#[test]
fn parent_cycles_end() {
    let data = load("cycle.sqlite3");
    let db = TaxonomySqlite::open(&data).unwrap();
    let d = |t: i32| {
        let mut v = db.descendants(t).unwrap();
        v.sort_unstable();
        v
    };
    assert_eq!(d(1), vec![10, 11]);
    assert_eq!(d(10), vec![11]);
    assert_eq!(d(20), vec![21, 22]);
    assert_eq!(d(21), vec![20, 22]);
    assert_eq!(d(22), Vec::<i32>::new());
    assert_eq!(d(30), Vec::<i32>::new());
    assert_eq!(d(31), Vec::<i32>::new());
}

#[test]
fn wide_taxid_truncates_like_getint() {
    // gen_fx.py stores taxid 5_000_000_000 as a child of the node named in the
    // `wide` line of the reference, and 5_000_000_001 as its child. GetInt keeps
    // the low 32 bits: 5_000_000_000 - 2^32 = 705_032_704. (The reference equality
    // of check_fixture already covers it; this names the values.)
    let data = load("p512.sqlite3");
    let db = TaxonomySqlite::open(&data).unwrap();
    let refs = std::fs::read_to_string(fx_dir().join("p512.ref")).unwrap();
    let wide: i32 = refs
        .lines()
        .find_map(|l| l.strip_prefix("wide "))
        .unwrap()
        .parse()
        .unwrap();
    let got = db.descendants(wide).unwrap();
    assert!(
        got.contains(&705_032_704) && got.contains(&705_032_705),
        "{:?}",
        got
    );
}

#[test]
fn sanity_check_failures() {
    let data = load("noindex.sqlite3");
    let e = TaxonomySqlite::open_named(&data, "/db/taxonomy4blast.sqlite3")
        .err()
        .unwrap()
        .to_string();
    assert_eq!(
        e,
        "Database '/db/taxonomy4blast.sqlite3' does not have index TaxidInfoCompositeIdx_parent. Please run the following command or contact your system administrator to install it:\nupdate_blastdb.pl --decompress taxdb"
    );
    let e = TaxonomySqlite::open(&data).err().unwrap().to_string();
    assert!(
        e.starts_with("Database 'taxonomy4blast.sqlite3' does not have index "),
        "{}",
        e
    );
    // an empty file is an empty database: both missing, the index is reported first
    let empty = load("empty.bin");
    let e = TaxonomySqlite::open(&empty).err().unwrap().to_string();
    assert!(
        e.contains("does not have index TaxidInfoCompositeIdx_parent"),
        "{}",
        e
    );
    // the table row removed from sqlite_master, the index there
    let d = load("notable.sqlite3");
    let e = TaxonomySqlite::open(&d).err().unwrap().to_string();
    assert!(
        e.contains("does not have table TaxidInfo. Please run"),
        "{}",
        e
    );
}

#[test]
fn not_a_database() {
    let d = load("notdb.bin");
    let e = TaxonomySqlite::open(&d).err().unwrap().to_string();
    assert!(e.ends_with("[26] file is not a database"), "{}", e);
    let e = TaxonomySqlite::open(&b"SQLite format 3\0"[..])
        .err()
        .unwrap()
        .to_string();
    assert!(e.contains("not a database"), "{}", e);
}

#[test]
fn wal_refused() {
    let mut d = load("p512.sqlite3");
    d[18] = 2;
    d[19] = 2;
    let e = TaxonomySqlite::open(&d).err().unwrap().to_string();
    assert!(e.contains("WAL"), "{}", e);
}

#[test]
fn page_in_memory_and_vec_sources() {
    // the two impls of ReadAt answer the same
    let v = load("p512.sqlite3");
    let a = TaxonomySqlite::open(&v).unwrap();
    let b = TaxonomySqlite::open(v.as_slice()).unwrap();
    let first: i32 = {
        let refs = std::fs::read_to_string(fx_dir().join("p512.ref")).unwrap();
        refs.lines()
            .filter(|l| {
                !l.starts_with("full")
                    && !l.starts_with("wide")
                    && !l.contains("LOOP")
                    && !l.contains(" 0 ")
            })
            .map(|l| l.split(' ').next().unwrap().parse::<i64>().unwrap())
            .find(|&t| t < i64::from(i32::MAX))
            .unwrap() as i32
    };
    assert_eq!(a.descendants(first).unwrap(), b.descendants(first).unwrap());
    assert!(!a.descendants(first).unwrap().is_empty());
}

/// Random byte-damaged copies per fixture (the source test used 1500; kept short for a debug build).
const DAMAGE_ROUNDS: usize = 200;

/// Damaged files must give `Err` (or a result), never a panic or a hang.
#[test]
fn damaged_files_never_panic() {
    let mut seed: u64 = 0x1234_5678_9abc_def0;
    let mut rnd = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for name in ["p512.sqlite3", "utf16.sqlite3"] {
        let orig = load(name);
        // the root, a node without children and the largest subtrees of the fixture
        let refs =
            std::fs::read_to_string(fx_dir().join(name.replace(".sqlite3", ".ref"))).unwrap();
        let mut big: Vec<(usize, i32)> = refs
            .lines()
            .filter(|l| !l.starts_with("wide") && !l.starts_with("full"))
            .filter_map(|l| {
                let f: Vec<&str> = l.split(' ').collect();
                Some((
                    f.get(1)?.parse().ok()?,
                    i32::try_from(f[0].parse::<i64>().ok()?).ok()?,
                ))
            })
            .collect();
        big.sort_unstable_by(|a, b| b.cmp(a));
        let mut ids: Vec<i32> = vec![1, 2, 705_032_704];
        ids.extend(big.iter().take(3).map(|b| b.1));
        // truncations
        for cut in [
            0usize,
            50,
            99,
            100,
            101,
            511,
            1000,
            orig.len() / 3,
            orig.len() / 2,
            orig.len() - 1,
        ] {
            let d = &orig[..cut.min(orig.len())];
            if let Ok(db) = TaxonomySqlite::open(d) {
                for &t in &ids {
                    let _ = db.descendants(t);
                }
            }
        }
        // random byte damage, mostly in the first pages and the cell pointers
        let mut opens = 0;
        for round in 0..DAMAGE_ROUNDS {
            let mut d = orig.clone();
            for _ in 0..(1 + round % 6) {
                let at = if rnd() % 3 == 0 {
                    (rnd() % 160) as usize
                } else {
                    (rnd() % d.len() as u64) as usize
                };
                d[at] = rnd() as u8;
            }
            if let Ok(db) = TaxonomySqlite::open(&d) {
                opens += 1;
                for &t in &ids {
                    let _ = db.descendants(t);
                }
            }
        }
        assert!(opens > 100, "{}: only {} damaged files opened", name, opens);
    }
}

// --- the real file ------------------------------------------------------------

/// Positioned reads of a real file (pread), the way the integration will use it.
struct FileSrc(std::fs::File, u64);

impl ReadAt for FileSrc {
    fn len(&self) -> u64 {
        self.1
    }
    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
        ReadAt::read_exact_at(&self.0, buf, offset)
    }
}

/// Needs the real NCBI file: `LOSAT_TEST_TAX_SQLITE=/path/to/taxonomy4blast.sqlite3` (3.4 million
/// rows) and `LOSAT_TEST_TAX_SQLITE_EXPECTED`, the file that `tests/seqdb/tax_sqlite_ref.py`
/// writes from it with Python's `sqlite3`: a line `T: sorted descendants` for 212 taxids (NCBI's
/// query run by SQLite) and for taxid 1, whose NCBI query never ends (a walk of the whole tree).
#[test]
#[ignore = "needs LOSAT_TEST_TAX_SQLITE (the real taxonomy4blast.sqlite3) and LOSAT_TEST_TAX_SQLITE_EXPECTED"]
fn real_file_equals_sqlite3() {
    let path = std::env::var("LOSAT_TEST_TAX_SQLITE").expect("LOSAT_TEST_TAX_SQLITE");
    let expected =
        std::env::var("LOSAT_TEST_TAX_SQLITE_EXPECTED").expect("LOSAT_TEST_TAX_SQLITE_EXPECTED");
    let t0 = std::time::Instant::now();
    let f = std::fs::File::open(&path).unwrap();
    let file = FileSrc(f.try_clone().unwrap(), f.metadata().unwrap().len());
    let db = TaxonomySqlite::open_named(&file, &path).unwrap();
    println!("real file: {} bytes, open {:?}", file.1, t0.elapsed());
    let refs = std::fs::read_to_string(&expected).unwrap();
    let mut total = 0usize;
    let mut count = 0usize;
    for line in refs.lines() {
        let (t, rest) = line.split_once(':').unwrap();
        let t: i32 = t.parse().unwrap();
        let want: Vec<i32> = rest
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        let t1 = std::time::Instant::now();
        // root: SQLite never finishes; LOSAT returns the whole tree
        let mut got = db
            .descendants(t)
            .unwrap_or_else(|e| panic!("taxid {}: {}", t, e));
        let el = t1.elapsed();
        got.sort_unstable();
        assert_eq!(got, want, "taxid {}", t);
        total += got.len();
        count += 1;
        if got.len() > 20000 || [1, 9606, 561, 10090].contains(&t) {
            println!("taxid {:>8}: {:>8} descendants in {:?}", t, got.len(), el);
        }
    }
    assert_eq!(count, 213, "{expected}: lines");
    println!(
        "{} taxids equal to sqlite3 ({} descendants in all)",
        count, total
    );
}

/// The bytes of `data`, whose reads from `from` fail as the system's do
/// ([`ReadError::Unreadable`]).
struct FailsFrom<'a> {
    data: &'a [u8],
    from: u64,
}

impl ReadAt for FailsFrom<'_> {
    fn len(&self) -> u64 {
        self.data.len() as u64
    }

    fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), ReadError> {
        if offset + buf.len() as u64 > self.from {
            return Err(ReadError::Unreadable(crate::seqdb::SeqDbError::Read {
                path: PathBuf::from("taxonomy4blast.sqlite3"),
                source: std::io::Error::from_raw_os_error(5),
            }));
        }
        ReadAt::read_exact_at(self.data, buf, offset)
    }
}

// NCBI reference (598d8ae6): c++/src/db/sqlite/sqlitewrapp.cpp:649-655 (cited at
// `QueryReadFailure`)
/// A read that the system fails after the sanity check, while the query steps, is
/// `ReadError::Unreadable` of `descendants`, for which the callers write NCBI's texts
/// (`QueryReadFailure`: NCBI BLAST+ 2.17.0 with the reads of `tests/db_search/taxfilter`'s file
/// failing with EIO from offset 4096, line WQA D18b round 2, probe L3); one that fails at the
/// sanity check is `Unreadable` of `open`, which the callers take as a file that is not there.
#[test]
fn a_read_that_fails_in_the_query_gives_ncbis_sqlite_texts() {
    let data = std::fs::read(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/taxfilter/taxonomy4blast.sqlite3"),
    )
    .unwrap();
    let late = FailsFrom {
        data: &data,
        from: 4096,
    };
    let db = TaxonomySqlite::open(&late).expect("the sanity check reads page 1 alone");
    assert!(matches!(
        db.descendants(9606),
        Err(ReadError::Unreadable(crate::seqdb::SeqDbError::Read { .. }))
    ));
    let early = FailsFrom {
        data: &data,
        from: 100,
    };
    assert!(matches!(
        TaxonomySqlite::open(&early),
        Err(ReadError::Unreadable(_))
    ));
    let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
    let place = |line: u32, kind: &str, function: &str| {
        format!(
            "    T0 \"{root}/src/db/sqlite/sqlitewrapp.cpp\", line {line}: Error: \
             (CSQLITE_Exception::{kind}) ncbi::{function}() - "
        )
    };
    let ncbi_finalize = format!(
        "Error: (1010.3) [blastp] [Error in x_Finalize()] Exception: NCBI C++ Exception:\n\
         {}Error from sqlite3: [266] disk I/O error\n\
         {}Cannot finalize statement\n\n",
        place(193, "eStmtFinalize", "s_ProcessErrorCode"),
        place(684, "eStmtFinalize", "CSQLITE_Statement::x_Finalize"),
    );
    assert_eq!(
        String::from_utf8(QueryReadFailure::finalize_post("blastp")).unwrap(),
        crate::wording::losat_text(&ncbi_finalize)
    );
    let ncbi_step = format!(
        "NCBI C++ Exception:\n\
         {}Error from sqlite3: [266] disk I/O error\n\
         {}Error stepping through statement results\n",
        place(193, "eStmtStep", "s_ProcessErrorCode"),
        place(762, "eStmtStep", "CSQLITE_Statement::Step"),
    );
    assert_eq!(
        String::from_utf8(QueryReadFailure::step().what()).unwrap(),
        crate::wording::losat_text(&ncbi_step)
    );
}
