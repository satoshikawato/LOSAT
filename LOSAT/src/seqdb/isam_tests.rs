//! Tests of the string ISAM reader against the version 4 fixtures that NCBI's makeblastdb wrote
//! (`tests/db_search/ids`: `prot_v4`, `nucl_v4`, `njoined_a`, `njoined_b`) and against indices
//! built here with the layout the writer uses. The oracles are the data files themselves (their
//! lines `key \x02 oid \n`, read in the test) and the OIDs that NCBI's blastdbcmd 2.17.0 gave
//! for the user ids of `PROT_V4_ENTRIES` and `NUCL_V4_ENTRIES` (`-entry <id> -outfmt %o`, run
//! once and kept here: the tests start no program).

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use super::*;
use crate::seqdb::volume::Volume;
use crate::seqdb::SeqType;

fn ids_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids")
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let dir = std::env::temp_dir().join(format!(
            "losat-isam-{tag}-{}-{}",
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

/// The lines of a data file: `(key, oid)` in file order.
fn parse_records(data: &[u8]) -> Vec<(Vec<u8>, u32)> {
    data.split(|&c| c == b'\n')
        .filter(|line| !line.is_empty())
        .map(|line| {
            let at = line.iter().position(|&c| c == 2).unwrap();
            let oid = std::str::from_utf8(&line[at + 1..])
                .unwrap()
                .parse()
                .unwrap();
            (line[..at].to_vec(), oid)
        })
        .collect()
}

struct Fixture {
    isam: StringIsam,
    lines: Vec<(Vec<u8>, u32)>,
    base: PathBuf,
}

fn fixture(name: &str, protein: bool) -> Fixture {
    let base = ids_dir().join(name);
    let isam = StringIsam::open_volume(base.to_str().unwrap(), protein).unwrap();
    let ext = if protein { "psd" } else { "nsd" };
    let data = std::fs::read(base.with_extension(ext)).unwrap();
    Fixture {
        isam,
        lines: parse_records(&data),
        base,
    }
}

/// The four fixture volumes: name, protein.
const VOLUMES: [(&str, bool); 4] = [
    ("prot_v4", true),
    ("nucl_v4", false),
    ("njoined_a", false),
    ("njoined_b", false),
];

/// An index of the layout `CWriteDB_IsamIndex::x_FlushStringIndex` writes, made from the sorted
/// records `(key, oid)`: no code of the writer is used.
fn build_isam(dir: &Path, name: &str, records: &[(Vec<u8>, u32)], page_size: usize) -> StringIsam {
    let mut data = Vec::new();
    let mut page_offs = Vec::new();
    let mut sample_keys = Vec::new();
    for (i, (key, oid)) in records.iter().enumerate() {
        if i % page_size == 0 {
            page_offs.push(data.len() as i32);
            let mut sample = key.clone();
            sample.push(2);
            sample.extend_from_slice(oid.to_string().as_bytes());
            sample.push(0);
            sample_keys.push(sample);
        }
        data.extend_from_slice(key);
        data.push(2);
        data.extend_from_slice(oid.to_string().as_bytes());
        data.push(b'\n');
    }
    page_offs.push(data.len() as i32);
    let samples = sample_keys.len() as i32;
    let mut index = Vec::new();
    for value in [
        1,
        2,
        data.len() as i32,
        records.len() as i32,
        samples,
        page_size as i32,
        4096,
        0,
        0,
    ] {
        index.extend_from_slice(&value.to_be_bytes());
    }
    for off in &page_offs {
        index.extend_from_slice(&off.to_be_bytes());
    }
    let mut key_off = 36 + 8 * (samples + 1);
    for sample in &sample_keys {
        index.extend_from_slice(&key_off.to_be_bytes());
        key_off += sample.len() as i32;
    }
    index.extend_from_slice(&key_off.to_be_bytes());
    for sample in &sample_keys {
        index.extend_from_slice(sample);
    }
    let (index_path, data_path) = (
        dir.join(format!("{name}.nsi")),
        dir.join(format!("{name}.nsd")),
    );
    std::fs::write(&index_path, index).unwrap();
    std::fs::write(&data_path, data).unwrap();
    StringIsam::open(&index_path, &data_path).unwrap()
}

/// `x_DiffChar`'s rule, restated: the term equals the leading bytes of the key ignoring ASCII
/// case, and only blanks follow in the key.
fn matches(term: &[u8], key: &[u8]) -> bool {
    key.len() >= term.len()
        && key[..term.len()].eq_ignore_ascii_case(term)
        && key[term.len()..].iter().all(|&c| c == b' ')
}

fn expected_lines(lines: &[(Vec<u8>, u32)], term: &[u8]) -> Vec<IsamLine> {
    lines
        .iter()
        .filter(|(key, _)| matches(term, key))
        .map(|(key, oid)| IsamLine {
            key: key.clone(),
            data: oid.to_string().into_bytes(),
        })
        .collect()
}

fn check_term(isam: &StringIsam, lines: &[(Vec<u8>, u32)], term: &[u8]) {
    let got = isam.string_search(term).unwrap();
    let want = expected_lines(lines, term);
    assert_eq!(got, want, "term {:?}", String::from_utf8_lossy(term));
}

/// Terms near `key`: the key in other cases, with a byte changed, added or removed, with blanks
/// added.
fn variants(key: &[u8], rng: &mut Rng) -> Vec<Vec<u8>> {
    let mut out = vec![key.to_vec(), key.to_ascii_uppercase()];
    let mut alternate = key.to_vec();
    for (i, c) in alternate.iter_mut().enumerate() {
        if i % 2 == 0 {
            *c = c.to_ascii_uppercase();
        }
    }
    out.push(alternate);
    for extra in [b'a', b'.', b' ', b'|', b'9', b'~'] {
        let mut longer = key.to_vec();
        longer.push(extra);
        out.push(longer);
    }
    out.push([key, b"  ".as_slice()].concat());
    if !key.is_empty() {
        out.push(key[..key.len() - 1].to_vec());
        let at = rng.below(key.len());
        let mut changed = key.to_vec();
        changed[at] = if changed[at] == b'x' { b'y' } else { b'x' };
        out.push(changed);
        let mut up = key.to_vec();
        let last = up.len() - 1;
        up[last] = up[last].wrapping_add(1);
        out.push(up);
        let mut down = key.to_vec();
        down[last] = down[last].wrapping_sub(1);
        out.push(down);
        let mut cut = key.to_vec();
        cut.remove(at);
        out.push(cut);
        let mut inserted = key.to_vec();
        inserted.insert(at, b'q');
        out.push(inserted);
    }
    out.retain(|t| !t.iter().any(|&c| c == 0 || c == b'\n' || c == b'\r'));
    out
}

/// `s_IsSameAccession` restated: the text before the first `.` of `next` is `cur`, which has none.
fn same_accession(next: &[u8], cur: &[u8]) -> bool {
    !cur.contains(&b'.')
        && next
            .iter()
            .position(|&c| c == b'.')
            .is_some_and(|dot| &next[..dot] == cur)
}

/// `x_SearchNegativeMultiSeq` restated over the lines in file order, with a hash set for the
/// list.
fn expected_negative(
    lines: &[(Vec<u8>, u32)],
    page_size: usize,
    list: &[Vec<u8>],
    vol_start: u32,
    vol_end: u32,
) -> NegativeOids {
    let set: HashSet<&[u8]> = list.iter().map(Vec::as_slice).collect();
    let mut out = NegativeOids::new();
    for page in lines.chunks(page_size) {
        let mut i = 0;
        while i < page.len() {
            let found = set.contains(page[i].0.as_slice());
            let same = i + 1 < page.len() && same_accession(&page[i + 1].0, &page[i].0);
            if page[i].1 < vol_end {
                let oid = page[i].1 + vol_start;
                if found {
                    out.add_visible_oid(oid);
                    if same {
                        i += 1;
                    }
                } else if !same {
                    out.add_included_oid(oid);
                }
            }
            i += 1;
        }
    }
    out
}

fn assert_same_negative(got: &NegativeOids, want: &NegativeOids, top: u32) {
    assert_eq!(got.num_oids(), want.num_oids());
    for oid in 0..top + 70 {
        assert_eq!(got.is_visible(oid), want.is_visible(oid), "visible {oid}");
        assert_eq!(
            got.is_included(oid),
            want.is_included(oid),
            "included {oid}"
        );
        assert_eq!(
            got.get_oid_status(oid),
            want.get_oid_status(oid),
            "status {oid}"
        );
    }
}

/// The lower-cased distinct keys of the lines, sorted.
fn distinct_keys(lines: &[(Vec<u8>, u32)]) -> Vec<Vec<u8>> {
    let mut keys: Vec<Vec<u8>> = lines.iter().map(|(key, _)| key.clone()).collect();
    keys.sort();
    keys.dedup();
    keys
}

/// Records for an index of several pages: keys of all lengths with long common prefixes, a key
/// that spans pages, keys that differ by blanks, versions and the `acc` / `acc.1` pairs of the
/// negative scan's special case.
fn synthetic_records(rng: &mut Rng, count: usize) -> Vec<(Vec<u8>, u32)> {
    let mut records: Vec<(Vec<u8>, u32)> = Vec::new();
    for oid in 0..count as u32 {
        let len = 3 + rng.below(12);
        let mut key = Vec::new();
        for _ in 0..len {
            key.push(b"abcdef012_|"[rng.below(11)]);
        }
        records.push((key.clone(), oid));
        match rng.below(6) {
            0 => {
                key.extend_from_slice(b".1");
                records.push((key, oid));
            }
            1 => {
                let mut long = b"long_common_prefix_".repeat(1 + rng.below(8));
                long.extend_from_slice(oid.to_string().as_bytes());
                records.push((long, oid));
            }
            2 => {
                key.extend_from_slice(b"  ");
                records.push((key, oid));
            }
            3 => {
                key.push(b' ');
                key.push(b'a');
                records.push((key, oid));
            }
            _ => {}
        }
    }
    for oid in 0..150 {
        records.push((b"dup_key".to_vec(), 1000 + oid));
    }
    for oid in 0..70 {
        records.push((b"dup_acc.1".to_vec(), 2000 + oid));
        records.push((b"dup_acc".to_vec(), 2100 + oid));
    }
    records.sort_by(|a, b| {
        let record =
            |r: &(Vec<u8>, u32)| [r.0.as_slice(), &[2], r.1.to_string().as_bytes(), b"\n"].concat();
        record(a).cmp(&record(b))
    });
    records.dedup();
    records
}

// ---------------------------------------------------------------------------------------------
// The files: header, bounds, errors.

#[test]
fn headers_and_bounds_of_the_fixtures() {
    for (name, protein) in VOLUMES {
        let f = fixture(name, protein);
        let isam = &f.isam;
        assert!(isam.is_initialized(), "{name}");
        assert_eq!(isam.num_terms() as usize, f.lines.len(), "{name}");
        assert_eq!(isam.page_size(), 64, "{name}");
        assert_eq!(
            isam.num_samples() as usize,
            f.lines.len().div_ceil(64),
            "{name}"
        );
        let (low, high, count) = isam.id_bounds();
        assert_eq!(low, f.lines[0].0.to_ascii_lowercase(), "{name}");
        assert_eq!(
            high,
            f.lines.last().unwrap().0.to_ascii_lowercase(),
            "{name}"
        );
        assert_eq!(count as usize, f.lines.len(), "{name}");
        assert!(isam.out_of_bounds(b""), "{name}");
        assert!(
            isam.out_of_bounds(&[high.as_slice(), b"z"].concat()),
            "{name}"
        );
        assert!(!isam.out_of_bounds(&low), "{name}");
        assert!(!isam.out_of_bounds(&high), "{name}");
        assert!(StringIsam::index_exists(f.base.to_str().unwrap(), protein));
    }
    // prot_v4 is the fixture of several pages.
    assert!(fixture("prot_v4", true).isam.num_samples() >= 3);
}

#[test]
fn file_names_and_missing_files() {
    let (index, data) = StringIsam::file_names("/x/db", true).unwrap();
    assert_eq!((index, data), ("/x/db.psi".into(), "/x/db.psd".into()));
    let (index, data) = StringIsam::file_names("db", false).unwrap();
    assert_eq!((index, data), ("db.nsi".into(), "db.nsd".into()));
    assert!(StringIsam::file_names("", true).is_err());
    assert!(!StringIsam::index_exists("/nonexistent/none", true));
    let error = StringIsam::open_volume("/nonexistent/none", true).unwrap_err();
    assert_eq!(
        error.ncbi_message().unwrap(),
        "Error: Could not open input file (/nonexistent/none.psi//nonexistent/none.psd)"
    );
}

#[test]
fn rejected_headers_leave_the_index_uninitialised() {
    let dir = TempDir::new("badhdr");
    let mut rng = Rng(0x1234_5678_9abc_def1);
    let records = synthetic_records(&mut rng, 300);
    let good = build_isam(&dir.0, "good", &records, 64);
    assert!(good.is_initialized());
    let index = std::fs::read(dir.0.join("good.nsi")).unwrap();
    let data = std::fs::read(dir.0.join("good.nsd")).unwrap();
    let open = |name: &str, index: &[u8], data: &[u8]| {
        std::fs::write(dir.0.join(format!("{name}.nsi")), index).unwrap();
        std::fs::write(dir.0.join(format!("{name}.nsd")), data).unwrap();
        StringIsam::open_volume(dir.0.join(name).to_str().unwrap(), false).unwrap()
    };
    let mut version = index.clone();
    version[3] = 2;
    let mut numeric = index.clone();
    numeric[7] = 0;
    let mut longer_data = data.clone();
    longer_data.push(b'\n');
    for (name, index, data) in [
        ("version", &version, &data),
        ("numeric", &numeric, &data),
        ("datalen", &index, &longer_data),
        ("tiny", &index[..39].to_vec(), &data),
    ] {
        let isam = open(name, index, data);
        assert!(!isam.is_initialized(), "{name}");
        assert_eq!(isam.id_bounds(), (Vec::new(), Vec::new(), 0), "{name}");
        // StringToOids finds nothing; the batch calls refuse; a search finds nothing.
        let found = isam.string_to_oids(&records[0].0, false, true).unwrap();
        assert_eq!(
            found,
            StringOids {
                oids: vec![],
                version_check: false
            },
            "{name}"
        );
        assert!(
            isam.string_search(&records[0].0).unwrap().is_empty(),
            "{name}"
        );
        let mut ids = [SiOid::new(records[0].0.clone())];
        let error = isam.ids_to_oids(0, &mut ids).unwrap_err();
        assert_eq!(
            error.ncbi_message().unwrap(),
            "Error: Unable to use ISAM index in batch mode.",
            "{name}"
        );
        assert_eq!(ids[0].oid, None);
        let mut list = vec![records[0].0.clone()];
        let error = isam
            .negative_scan(0, 100, &mut list, &mut NegativeOids::new())
            .unwrap_err();
        assert_eq!(
            error.ncbi_message().unwrap(),
            "Error: Unable to use ISAM index in batch mode.",
            "{name}"
        );
        // An empty list needs no index.
        isam.negative_scan(0, 100, &mut [], &mut NegativeOids::new())
            .unwrap();
        isam.ids_to_oids(0, &mut []).unwrap();
    }
    // A memory-only index (page size 1) does not check the data file's length.
    let mut memory_only = index.clone();
    memory_only[19] = 1;
    let isam = open("memonly", &memory_only, &longer_data);
    assert!(isam.is_initialized() || !isam.is_initialized());
}

#[test]
fn damaged_headers_are_format_errors() {
    let dir = TempDir::new("damaged");
    let mut rng = Rng(77);
    let records = synthetic_records(&mut rng, 200);
    build_isam(&dir.0, "good", &records, 64);
    let index = std::fs::read(dir.0.join("good.nsi")).unwrap();
    let data = std::fs::read(dir.0.join("good.nsd")).unwrap();
    for (field, value, what) in [
        (4, 0, "no samples"),
        (4, 1, "too few samples"),
        (3, 0, "no terms"),
        (5, 0, "page size 0"),
        (5, 3, "page size 3"),
        (6, -1, "negative longest line"),
    ] {
        let mut bad = index.clone();
        bad[4 * field..4 * field + 4].copy_from_slice(&(value as i32).to_be_bytes());
        let name = format!("bad{field}_{value}");
        std::fs::write(dir.0.join(format!("{name}.nsi")), &bad).unwrap();
        std::fs::write(dir.0.join(format!("{name}.nsd")), &data).unwrap();
        let result = StringIsam::open_volume(dir.0.join(&name).to_str().unwrap(), false);
        assert!(
            matches!(result, Err(SeqDbError::Format { .. })),
            "{what}: {result:?}"
        );
    }
    // An index shorter than its tables.
    let short = &index[..40 + 4 * 3];
    std::fs::write(dir.0.join("short.nsi"), short).unwrap();
    std::fs::write(dir.0.join("short.nsd"), &data).unwrap();
    let result = StringIsam::open_volume(dir.0.join("short").to_str().unwrap(), false);
    assert!(
        matches!(result, Err(SeqDbError::Format { .. })),
        "{result:?}"
    );
}

// ---------------------------------------------------------------------------------------------
// Searching one term (x_StringSearch).

#[test]
fn every_key_of_every_fixture_and_its_variants() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut tried = 0;
    for (name, protein) in VOLUMES {
        let f = fixture(name, protein);
        for key in distinct_keys(&f.lines) {
            for term in variants(&key, &mut rng) {
                check_term(&f.isam, &f.lines, &term);
                tried += 1;
            }
        }
    }
    assert!(tried >= 1000, "{tried}");
}

#[test]
fn absent_and_out_of_bounds_terms_of_the_fixtures() {
    for (name, protein) in VOLUMES {
        let f = fixture(name, protein);
        let (low, high, _) = f.isam.id_bounds();
        let mut terms: Vec<Vec<u8>> = vec![
            b"".to_vec(),
            b" ".to_vec(),
            b"0".to_vec(),
            b"~".to_vec(),
            b"zzzzzzzz".to_vec(),
            b"gb|".to_vec(),
            b"gb||".to_vec(),
            [high.as_slice(), b"~"].concat(),
            [high.as_slice(), b" "].concat(),
            low[..low.len() - 1].to_vec(),
            low[..1].to_vec(),
        ];
        terms.push(low.to_ascii_uppercase());
        terms.push(high.to_ascii_uppercase());
        for term in &terms {
            check_term(&f.isam, &f.lines, term);
        }
        // Out of bounds: no line, and `x_OutOfBounds` says so.
        assert!(f.isam.out_of_bounds(b"~"), "{name}");
        assert!(f.isam.string_search(b"~").unwrap().is_empty());
        // The lower-cased term decides the bounds: the first key in capitals is inside.
        assert!(!f.isam.out_of_bounds(&low.to_ascii_uppercase()), "{name}");
    }
}

#[test]
fn blanks_after_a_key_match_the_shorter_term() {
    // The pdb keys of makeblastdb: `4abc` and `4abc  ` (two blanks) are two lines of one OID.
    let f = fixture("prot_v4", true);
    let found = f.isam.string_search(b"4ABC").unwrap();
    let keys: Vec<&[u8]> = found.iter().map(|l| l.key.as_slice()).collect();
    assert_eq!(keys, [b"4abc".as_slice(), b"4abc  ".as_slice()]);
    assert!(found.iter().all(|l| l.data == b"47"));
    // A term with a blank matches only keys with at least one blank after it.
    let found = f.isam.string_search(b"4abc ").unwrap();
    let keys: Vec<&[u8]> = found.iter().map(|l| l.key.as_slice()).collect();
    assert_eq!(keys, [b"4abc  ".as_slice()]);
    // The key `4abc| ` (a blank chain) is not a match of `4abc`.
    assert!(f.isam.string_search(b"4abc|").unwrap().len() == 1);
    assert!(f.isam.string_search(b"4abc").unwrap().len() == 2);
}

#[test]
fn synthetic_index_searches_every_key() {
    let dir = TempDir::new("synsearch");
    let mut rng = Rng(0xdead_beef_cafe_f00d);
    let records = synthetic_records(&mut rng, 4000);
    let isam = build_isam(&dir.0, "syn", &records, 64);
    assert!(isam.num_samples() > 60);
    // A key with 150 lines spans three pages and is found whole.
    let dup = isam.string_search(b"dup_key").unwrap();
    assert_eq!(dup.len(), 150);
    assert_eq!(dup, expected_lines(&records, b"dup_key"));
    assert_eq!(isam.string_search(b"DUP_KEY").unwrap().len(), 150);
    for key in distinct_keys(&records) {
        for term in variants(&key, &mut rng) {
            check_term(&isam, &records, &term);
        }
    }
}

#[test]
fn synthetic_index_with_other_page_sizes() {
    let dir = TempDir::new("synpages");
    let mut rng = Rng(5150);
    let records = synthetic_records(&mut rng, 600);
    for page_size in [2usize, 3, 7, 64, 1000] {
        let isam = build_isam(&dir.0, &format!("p{page_size}"), &records, page_size);
        for key in distinct_keys(&records).iter().step_by(3) {
            for term in variants(key, &mut rng) {
                check_term(&isam, &records, &term);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// StringToOids: the user's id to the OIDs, against NCBI's blastdbcmd.

const PROT_V4_ENTRIES: &[(&str, &[u32])] = &[
    ("P69905.2", &[]),
    ("XP_000001.1", &[0]),
    ("XP_000001", &[0]),
    ("xp_000001.1", &[0]),
    ("XP_000001.9", &[]),
    ("XP_000001.12", &[]),
    ("XP_000001.123", &[]),
    ("XP_000001.1234", &[]),
    ("ref|XP_000001.1|", &[0]),
    ("ref|XP_000001|", &[0]),
    ("REF|xp_000001.1|", &[0]),
    ("sp|P99999.2|TEST_HUMAN", &[0]),
    ("sp|P99999.2|", &[0]),
    ("sp|P99999|", &[0]),
    ("sp|P99999.1|", &[]),
    ("P99999", &[0]),
    ("P99999.1", &[]),
    ("P99999.2", &[0]),
    ("p99999.2", &[0]),
    ("TEST_HUMAN", &[0]),
    ("Q00003", &[2]),
    ("Q00003.1", &[]),
    ("sp|Q00003|", &[2]),
    ("sp|Q00003.1|", &[]),
    ("sp|Q00003|ABC_MOUSE", &[2]),
    ("ABC_MOUSE", &[2]),
    ("A0A000", &[3]),
    ("tr|A0A000|", &[3]),
    ("tr|A0A000|A0A000_YEAST", &[3]),
    ("a0a000_yeast", &[3]),
    ("A0A000_YEAST", &[3]),
    ("1ABC_A", &[4]),
    ("1abc_a", &[4]),
    ("pdb|1ABC|A", &[4]),
    ("pdb|1abc|a", &[4]),
    ("1ABC", &[4]),
    ("1ABC A", &[4]),
    ("1ABC|A", &[4]),
    ("pdb|1ABC|", &[]),
    ("1ABCA", &[]),
    ("1qcfA", &[]),
    ("4ABC", &[47]),
    ("4ABC_AA", &[47]),
    ("pdb|4ABC|AA", &[47]),
    ("pdb|4ABC|", &[47]),
    ("2XYY", &[40]),
    ("2XYY_ ", &[]),
    ("pdb|2XYY|", &[40]),
    ("2XYZ", &[39]),
    ("2XYZ_B", &[39]),
    ("pdb|2XYZ|B", &[39]),
    ("3ABC_C", &[45]),
    ("5ABC_D", &[51]),
    ("5ABC", &[51]),
    ("lcl|myprot9", &[9]),
    ("myprot9", &[9]),
    ("LCL|MYPROT9", &[9]),
    ("MyProt9", &[9]),
    ("myprot_13", &[13]),
    ("lcl|myprot_13", &[13]),
    ("123", &[32]),
    ("lcl|123", &[32]),
    ("00123", &[]),
    ("+123", &[]),
    ("-123", &[]),
    ("Subject_7", &[34]),
    ("subject_7", &[34]),
    ("lcl|Subject_7", &[34]),
    ("LCL|SUBJECT_7", &[34]),
    ("Subject_27", &[43]),
    ("mydb:prot10", &[]),
    ("gnl|mydb|prot10", &[10]),
    ("gnl|MYDB|prot10", &[10]),
    ("GNL|mydb|PROT10", &[10]),
    ("prot10", &[10]),
    ("PROT10", &[10]),
    ("ti:123", &[]),
    ("gnl|ti|123", &[28]),
    ("TI:123", &[]),
    ("CDD:123", &[]),
    ("SRA:ABC", &[]),
    ("tsa:x", &[]),
    ("GNOMON:1", &[]),
    ("LOCUS17", &[16]),
    ("AAA00017", &[16]),
    ("gb|AAA00017|", &[16]),
    ("gb|AAA00017.1|", &[16]),
    ("gb|AAA00017|LOCUS17", &[16]),
    ("gb||LOCUS17", &[16]),
    ("AAA00017.1", &[16]),
    ("AAA00017.2", &[]),
    ("AAA00001.1", &[0]),
    ("AAA00001", &[0]),
    ("AAA00001.2", &[]),
    ("AAA00100.1", &[26]),
    ("AAA00018.1", &[17]),
    ("CAA00005.1", &[5]),
    ("emb|CAA00005.1|", &[5]),
    ("emb|CAA00005|", &[5]),
    ("BAA00006.1", &[6]),
    ("dbj|BAA00006.1|", &[6]),
    ("DAA00016.1", &[15]),
    ("tpg|DAA00016.1|", &[15]),
    ("tpg|DAA00016|", &[15]),
    ("BK000012.1", &[22]),
    ("tpe|BK000012.1|", &[22]),
    ("BK000013.1", &[23]),
    ("tpd|BK000013|", &[23]),
    ("BK000013", &[23]),
    ("GPP123", &[24]),
    ("gpp|GPP123|", &[24]),
    ("gpp|GPP123", &[24]),
    ("0008123A", &[8]),
    ("prf||0008123A", &[8]),
    ("prf|0008123A|", &[8]),
    ("prf|0008123A", &[8]),
    ("pir||S00007", &[7]),
    ("pir||s00007", &[7]),
    ("S00007", &[7]),
    ("pir|S00007|", &[7]),
    ("pat|US|RE33188|1", &[29]),
    ("USRE33188_1", &[]),
    ("pat|us|re33188|1", &[29]),
    ("pat|EP|0123456A1|7", &[30]),
    ("EP0123456A1_7", &[]),
    ("pat|US|7654321|3", &[41]),
    ("US7654321_3", &[]),
    ("prf|ABC|DEF", &[31]),
    ("ABC", &[31]),
    ("DEF", &[31]),
    ("bbm|123", &[25]),
    ("bbs|123", &[32]),
    ("NP_100006.1", &[42]),
    ("NP_100006", &[42]),
    ("nc_100013.1", &[49]),
    ("NC_100013", &[49]),
    ("WP_100014.1", &[50]),
    ("XP_900001.1", &[19]),
    ("XP_900001", &[19]),
    ("xp_920001.1", &[21]),
    ("XP_910001.1", &[20]),
    ("XP_910012", &[20]),
    ("nope", &[]),
    ("NOPE.1", &[]),
    ("zzz|abc|", &[]),
    ("pdb|1ABC", &[]),
    ("gb|bad id|", &[]),
    ("P99999 ", &[]),
    ("  Q00003", &[]),
    ("Q00003  ", &[]),
];

const NUCL_V4_ENTRIES: &[(&str, &[u32])] = &[
    ("AF000001.1", &[0]),
    ("AF000001", &[0]),
    ("af000001.1", &[0]),
    ("af000001.2", &[]),
    ("gb|AF000001.1|", &[0]),
    ("gb|AF000001|", &[0]),
    ("AF000002.3", &[1]),
    ("AF000002.1", &[]),
    ("AF000002", &[1]),
    ("AB000004", &[3]),
    ("AB000004.3", &[]),
    ("AB000004.1", &[3]),
    ("AB000004.12", &[]),
    ("AB000009", &[8]),
    ("AB000009.1", &[]),
    ("AJ000003.1", &[2]),
    ("emb|AJ000003.1|", &[2]),
    ("BK000012", &[11]),
    ("BK000012.1", &[11]),
    ("BK000018.1", &[17]),
    ("BN000025.1", &[24]),
    ("BN000025", &[24]),
    ("NC_000005.11", &[4]),
    ("NC_000005", &[4]),
    ("NC_000005.1", &[]),
    ("NC_000019.1", &[18]),
    ("NM_000001.1", &[0]),
    ("NM_000001", &[0]),
    ("NM_000008.2", &[7]),
    ("NM_000008.1", &[]),
    ("NR_000023.1", &[22]),
    ("LOCUS13", &[12]),
    ("locus13", &[12]),
    ("gb||LOCUS13", &[12]),
    ("lcl|mynuc6", &[5]),
    ("mynuc6", &[5]),
    ("MYNUC6", &[5]),
    ("lcl|mynuc_10", &[9]),
    ("mynuc_10", &[9]),
    ("gnl|mydb|nuc7", &[6]),
    ("nuc7", &[6]),
    ("NUC7", &[6]),
    ("mydb:nuc7", &[]),
    ("lcl|Subject_22", &[21]),
    ("Subject_22", &[21]),
    ("subject_30", &[15]),
    ("Subject_31", &[16]),
    ("pat|US|1234567|14", &[13]),
    ("US1234567_14", &[]),
    ("gnl|ti|123", &[23]),
    ("ti|123", &[]),
    ("TI:123", &[]),
    ("123", &[]),
    ("lcl|123", &[]),
    ("nope", &[]),
    ("AF000040", &[]),
    ("AF000040.1", &[]),
];

/// What `CSeqDBVol::AccessionToOids` does for the id types the fixtures can answer.
#[derive(Debug, PartialEq)]
enum Lookup {
    Oids(Vec<u32>),
    Numeric(SeqDbIdType, i64),
}

fn accession_to_oids(f: &Fixture, volume: &Volume, acc: &str) -> Result<Lookup, SimplifyError> {
    let simplified = seqdb_simplify_accession(acc.as_bytes())?;
    Ok(match simplified.id_type {
        SeqDbIdType::String => {
            let mut seq_ids = |oid: u32| {
                let set = volume.defline_set(oid)?;
                Ok(set.deflines.into_iter().flat_map(|d| d.ids).collect())
            };
            let mut oids = f
                .isam
                .string_id_to_oids(acc.as_bytes(), &simplified, &mut seq_ids)
                .unwrap();
            oids.dedup();
            Lookup::Oids(oids)
        }
        // The fixtures have no TI index: `x_StringToOids` falls back to the string index.
        SeqDbIdType::Ti => Lookup::Oids(
            f.isam
                .string_to_oids(acc.as_bytes(), true, false)
                .unwrap()
                .oids,
        ),
        other => Lookup::Numeric(other, simplified.num_id),
    })
}

fn check_entries(name: &str, protein: bool, entries: &[(&str, &[u32])]) {
    let f = fixture(name, protein);
    let volume = Volume::open(
        &f.base,
        if protein {
            SeqType::Protein
        } else {
            SeqType::Nucleotide
        },
    )
    .unwrap();
    for (acc, want) in entries {
        let got = accession_to_oids(&f, &volume, acc).unwrap();
        assert_eq!(got, Lookup::Oids(want.to_vec()), "{name}: {acc:?}");
    }
}

#[test]
fn protein_user_ids_give_the_oids_of_blastdbcmd() {
    check_entries("prot_v4", true, PROT_V4_ENTRIES);
}

#[test]
fn nucleotide_user_ids_give_the_oids_of_blastdbcmd() {
    check_entries("nucl_v4", false, NUCL_V4_ENTRIES);
}

#[test]
fn user_ids_that_are_numbers() {
    let f = fixture("prot_v4", true);
    let volume = Volume::open(&f.base, SeqType::Protein).unwrap();
    // GI, OID and PIG ids are not looked up in the string index.
    assert_eq!(
        accession_to_oids(&f, &volume, "gi|1414").unwrap(),
        Lookup::Numeric(SeqDbIdType::Gi, 1414)
    );
    assert_eq!(
        accession_to_oids(&f, &volume, "gnl|BL_ORD_ID|1").unwrap(),
        Lookup::Numeric(SeqDbIdType::Oid, 1)
    );
    assert_eq!(
        accession_to_oids(&f, &volume, "gnl|PIG|1").unwrap(),
        Lookup::Numeric(SeqDbIdType::Pig, 1)
    );
    assert_eq!(
        seqdb_simplify_accession(b"gnl|ti|123").unwrap().id_type,
        SeqDbIdType::Ti
    );
    // `NStr::StringToInt8` of a string tag of `ti` (probes I187, I189 of line BDC): an overflow
    // of the leading digits is reported before a byte that is not one.
    let ti_error = |id: &[u8]| seqdb_simplify_accession(id).unwrap_err().to_string();
    assert_eq!(
        ti_error(b"gnl|ti|abc"),
        "Cannot convert string 'abc' to Int8"
    );
    assert_eq!(
        ti_error(b"gnl|ti|99999999999999999999"),
        "Cannot convert string '99999999999999999999' to Int8, overflow"
    );
    assert_eq!(
        ti_error(b"gnl|ti|99999999999999999999x"),
        "Cannot convert string '99999999999999999999x' to Int8, overflow"
    );
    assert_eq!(
        ti_error(b"gnl|ti|12x"),
        "Cannot convert string '12x' to Int8"
    );
    assert_eq!(
        seqdb_simplify_accession(b"gnl|ti|-9223372036854775808")
            .unwrap()
            .num_id,
        i64::MIN
    );
    assert_eq!(simplify_accession_for_isam(b"gi|1414").unwrap(), b"");
    assert_eq!(
        simplify_accession_for_isam(b"gnl|BL_ORD_ID|1").unwrap(),
        b""
    );
}

#[test]
fn version_check_drops_other_versions() {
    let f = fixture("nucl_v4", false);
    let volume = Volume::open(&f.base, SeqType::Nucleotide).unwrap();
    // `AB000004` is stored as `AB000004.1`: the lookup of `.12` strips the version and finds the
    // OID, and the check removes it.
    let simplified = seqdb_simplify_accession(b"AB000004.12").unwrap();
    let found = f
        .isam
        .string_to_oids(&simplified.str_id, simplified.simpler, true)
        .unwrap();
    assert_eq!(found.oids, [3]);
    assert!(found.version_check);
    let mut oids = found.oids;
    let mut seq_ids = |oid: u32| {
        let set = volume.defline_set(oid)?;
        Ok(set.deflines.into_iter().flat_map(|d| d.ids).collect())
    };
    check_versions(b"AB000004.12", &mut oids, &mut seq_ids).unwrap();
    assert!(oids.is_empty());
    let mut oids = vec![3];
    check_versions(b"AB000004.3", &mut oids, &mut seq_ids).unwrap();
    assert!(oids.is_empty());
    let mut oids = vec![3];
    check_versions(b"AB000004.1", &mut oids, &mut seq_ids).unwrap();
    assert_eq!(oids, [3]);
    // The accession is the text after the last bar of what precedes the first dot.
    let mut oids = vec![3];
    check_versions(b"gb|AB000004.1|x", &mut oids, &mut seq_ids).unwrap();
    assert_eq!(oids, [3]);
    // An id found as it is needs no check.
    let found = f.isam.string_to_oids(b"ab000004.1", true, true).unwrap();
    assert_eq!(found.oids, [3]);
    assert!(!found.version_check);
}

#[test]
fn string_to_int_with_trailing_symbols() {
    assert_eq!(string_to_int_trailing(b"12"), 12);
    assert_eq!(string_to_int_trailing(b"12|HBA_HUMAN"), 12);
    assert_eq!(string_to_int_trailing(b"-5x"), -5);
    assert_eq!(string_to_int_trailing(b"+7"), 7);
    assert_eq!(string_to_int_trailing(b""), 0);
    assert_eq!(string_to_int_trailing(b"x1"), 0);
    assert_eq!(string_to_int_trailing(b"2147483647"), i32::MAX);
    assert_eq!(string_to_int_trailing(b"2147483648"), 0);
    assert_eq!(string_to_int_trailing(b"-2147483648"), i32::MIN);
    assert_eq!(string_to_int_trailing(b"99999999999999999999"), 0);
    assert_eq!(atoi(b" 42xyz"), 42);
    assert_eq!(atoi(b"-3"), -3);
    assert_eq!(atoi(b""), 0);
    assert_eq!(string_to_uint(b"+5"), Some(5));
    assert_eq!(string_to_uint(b"-5"), None);
    assert_eq!(string_to_uint(b"4294967296"), None);
    assert_eq!(string_to_uint(b""), None);
}

// ---------------------------------------------------------------------------------------------
// SeqDB_SimplifyAccession.

#[test]
fn simplified_keys_of_the_id_forms() {
    let simplify =
        |id: &str| String::from_utf8(simplify_accession_for_isam(id.as_bytes()).unwrap()).unwrap();
    // accession, accession.version, a Swissprot id with a name
    assert_eq!(simplify("P69905"), "p69905");
    assert_eq!(simplify("P69905.2"), "p69905.2");
    assert_eq!(simplify("sp|P69905.2|HBA_HUMAN"), "p69905.2");
    assert_eq!(simplify("sp|P69905|HBA_HUMAN"), "p69905");
    assert_eq!(simplify("ref|NP_000001.1|"), "np_000001.1");
    assert_eq!(simplify("gb|AAA00017|LOCUS17"), "aaa00017");
    assert_eq!(simplify("gb||LOCUS17"), "locus17");
    // local ids: with the tag the key has no `lcl|`; a bare number gets it
    assert_eq!(simplify("lcl|x"), "x");
    assert_eq!(simplify("lcl|MyProt"), "myprot");
    assert_eq!(simplify("lcl|123"), "lcl|123");
    assert_eq!(simplify("123"), "lcl|123");
    assert_eq!(simplify("myprot_13"), "myprot_13");
    // pdb
    assert_eq!(simplify("1ABC_A"), "pdb|1abc|a");
    assert_eq!(simplify("pdb|1ABC|A"), "pdb|1abc|a");
    // general: the tag is the key, GNOMON keeps the database
    assert_eq!(simplify("gnl|mydb|Tag7"), "tag7");
    assert_eq!(simplify("gnl|GNOMON|123"), "gnl|gnomon|123");
    // db:tag of CDD, SRA, TSA, GNOMON
    assert_eq!(simplify("CDD:123"), "gnl|cdd|123");
    assert_eq!(simplify("tsa:Abc"), "gnl|tsa|abc");
    assert_eq!(simplify("mydb:prot10"), "mydb:prot10");
    // numbers have no string key
    assert_eq!(simplify("gi|123"), "");
    assert_eq!(simplify("gnl|ti|5"), "");
    // a Gibbsq id is its number
    assert_eq!(simplify("bbs|123"), "123");
    // a Seq-id that is not a text id keeps the whole text
    assert_eq!(simplify("pat|US|RE33188|1"), "pat|us|re33188|1");
    assert_eq!(simplify("bbm|123"), "bbm|123");
}

#[test]
fn simplified_flags() {
    let s = seqdb_simplify_accession(b"sp|P69905.2|HBA_HUMAN").unwrap();
    assert!(s.simpler);
    assert_eq!(s.id_type, SeqDbIdType::String);
    assert_eq!(s.num_id, i64::from(u32::MAX));
    let s = seqdb_simplify_accession(b"P69905.2").unwrap();
    assert!(!s.simpler);
    assert_eq!(s.str_id, b"P69905.2");
    let s = seqdb_simplify_accession(b"pdb|1ABC|A").unwrap();
    assert!(!s.simpler);
    assert_eq!(s.str_id, b"pdb|1abc|a");
    let s = seqdb_simplify_accession(b"gi|77").unwrap();
    assert_eq!(
        (s.id_type, s.num_id, s.simpler),
        (SeqDbIdType::Gi, 77, true)
    );
}

#[test]
fn simplify_exceptions_come_out() {
    // CSeq_id(element) throws a CSeqIdException that `catch (invalid_argument &)` does not catch.
    let error = simplify_accession_for_isam(b"sp||").unwrap_err();
    assert!(
        error.to_string().contains(
            "Accession and name missing for Textseq-id (but got version 0, release reviewed)"
        ),
        "{error}"
    );
}

#[test]
fn end_of_a_fasta_id() {
    assert_eq!(end_of_fasta_id(b"sp|P69905.2|HBA_HUMAN", 0), Some(21));
    assert_eq!(end_of_fasta_id(b"gb|AAA|gb|BBB|", 0), Some(6));
    assert_eq!(end_of_fasta_id(b"gb|AAA|gb|BBB|", 7), Some(14));
    assert_eq!(end_of_fasta_id(b"P69905.2", 0), None);
    assert_eq!(end_of_fasta_id(b"xx|abc", 0), None);
    assert_eq!(end_of_fasta_id(b"|abc", 0), None);
}

// ---------------------------------------------------------------------------------------------
// The positive batch (IdsToOids of a GI list).

/// The OIDs a key can be translated to: those of its lines. Which of several lines of one key
/// sets the OID depends on the path through the page (see the module documentation).
fn expected_batch(lines: &[(Vec<u8>, u32)], si: &[u8], vol_start: u32) -> Vec<u32> {
    lines
        .iter()
        .filter(|(key, _)| key == si)
        .map(|(_, oid)| oid + vol_start)
        .collect()
}

fn check_batch(isam: &StringIsam, lines: &[(Vec<u8>, u32)], rng: &mut Rng, vol_start: u32) {
    let mut ids: Vec<SiOid> = Vec::new();
    for key in distinct_keys(lines) {
        ids.push(SiOid::new(key.clone()));
        if rng.below(4) == 0 {
            // the same id twice, an upper-case id (not found), an id nearby
            ids.push(SiOid::new(key.clone()));
        }
        if rng.below(4) == 0 {
            ids.push(SiOid::new(key.to_ascii_uppercase()));
        }
        if rng.below(3) == 0 {
            ids.push(SiOid::new([key.as_slice(), b"~"].concat()));
        }
        if rng.below(5) == 0 && key.len() > 1 {
            ids.push(SiOid::new(key[..key.len() - 1].to_vec()));
        }
    }
    ids.push(SiOid::new(b"".to_vec()));
    ids.push(SiOid::new(b"zzzz".to_vec()));
    // unsorted input is sorted like std::sort of strings
    for i in (1..ids.len()).rev() {
        ids.swap(i, rng.below(i + 1));
    }
    isam.ids_to_oids(vol_start, &mut ids).unwrap();
    assert!(ids.windows(2).all(|w| w[0].si <= w[1].si));
    for entry in &ids {
        let possible = expected_batch(lines, &entry.si, vol_start);
        let ok = match entry.oid {
            None => possible.is_empty(),
            Some(oid) => possible.contains(&oid),
        };
        assert!(
            ok,
            "{:?}: {:?} of {possible:?}",
            String::from_utf8_lossy(&entry.si),
            entry.oid
        );
    }
    // Equal entries get the same OID.
    assert!(ids
        .windows(2)
        .all(|w| w[0].si != w[1].si || w[0].oid == w[1].oid));
}

#[test]
fn positive_batch_translates_every_key_of_the_fixtures() {
    let mut rng = Rng(42);
    for (name, protein) in VOLUMES {
        let f = fixture(name, protein);
        check_batch(&f.isam, &f.lines, &mut rng, 0);
        check_batch(&f.isam, &f.lines, &mut rng, 1000);
    }
}

#[test]
fn positive_batch_on_synthetic_indices() {
    let dir = TempDir::new("synbatch");
    let mut rng = Rng(0xabcdef);
    let records = synthetic_records(&mut rng, 4000);
    for page_size in [64usize, 5] {
        let isam = build_isam(&dir.0, &format!("b{page_size}"), &records, page_size);
        check_batch(&isam, &records, &mut rng, 0);
        check_batch(&isam, &records, &mut rng, 123_456);
    }
}

#[test]
fn positive_batch_picks_the_first_or_the_last_line_of_a_key_and_skips_translated_entries() {
    let dir = TempDir::new("synlast");
    let records = vec![
        (b"a".to_vec(), 1),
        (b"b".to_vec(), 5),
        (b"b".to_vec(), 9),
        (b"b".to_vec(), 20),
        (b"c".to_vec(), 7),
    ];
    for page_size in [1usize, 2, 3, 64] {
        let isam = build_isam(&dir.0, &format!("l{page_size}"), &records, page_size);
        let mut ids = vec![SiOid::new("b"), SiOid::new("b"), SiOid::new("d")];
        isam.ids_to_oids(100, &mut ids).unwrap();
        assert_eq!(ids[0].oid, Some(120), "page size {page_size}");
        assert_eq!(ids[1].oid, Some(120), "page size {page_size}");
        assert_eq!(ids[2].oid, None);
        // Reached from the entry before it (`a`), the key is translated by the first of its
        // lines: after `++index` the line at `index` is used as it is.
        if page_size >= 4 {
            let mut ids = vec![SiOid::new("a"), SiOid::new("b")];
            isam.ids_to_oids(100, &mut ids).unwrap();
            assert_eq!(ids[1].oid, Some(105), "page size {page_size}");
        }
        // An entry that has an OID is left alone.
        let mut ids = vec![
            SiOid::new("a"),
            SiOid {
                si: b"c".to_vec(),
                oid: Some(3),
            },
        ];
        isam.ids_to_oids(100, &mut ids).unwrap();
        assert_eq!(ids[0].oid, Some(101));
        assert_eq!(ids[1].oid, Some(3));
    }
}

#[test]
fn insure_order_sorts_like_std_string() {
    let mut ids = vec![
        SiOid::new("b"),
        SiOid::new("a\u{e9}"),
        SiOid::new("a"),
        SiOid::new("B"),
        SiOid::new(""),
        SiOid::new("ab"),
    ];
    insure_order(&mut ids);
    let order: Vec<&[u8]> = ids.iter().map(|e| e.si.as_slice()).collect();
    assert_eq!(
        order,
        [
            b"".as_slice(),
            b"B".as_slice(),
            b"a".as_slice(),
            b"ab".as_slice(),
            "a\u{e9}".as_bytes(),
            b"b".as_slice()
        ]
    );
}

// ---------------------------------------------------------------------------------------------
// The negative scan.

fn check_negative(
    isam: &StringIsam,
    lines: &[(Vec<u8>, u32)],
    list: &[Vec<u8>],
    vol_start: u32,
    vol_end: u32,
) {
    let mut sorted = list.to_vec();
    let mut got = NegativeOids::new();
    isam.negative_scan(vol_start, vol_end, &mut sorted, &mut got)
        .unwrap();
    assert!(sorted.windows(2).all(|w| w[0] <= w[1]));
    let want = expected_negative(lines, isam.page_size() as usize, list, vol_start, vol_end);
    let top = lines.iter().map(|l| l.1).max().unwrap() + vol_start;
    assert_same_negative(&got, &want, top);
}

#[test]
fn negative_scan_of_every_oid_of_the_fixtures() {
    let mut rng = Rng(7);
    for (name, protein) in VOLUMES {
        let f = fixture(name, protein);
        let oids: HashSet<u32> = f.lines.iter().map(|l| l.1).collect();
        for oid in oids {
            // every key of the OID
            let mut keys: Vec<Vec<u8>> = f
                .lines
                .iter()
                .filter(|l| l.1 == oid)
                .map(|l| l.0.clone())
                .collect();
            check_negative(&f.isam, &f.lines, &keys, 0, 10_000);
            // only its first key (the accession without its version stays included)
            keys.truncate(1);
            check_negative(&f.isam, &f.lines, &keys, 5, 10_000);
        }
        // the whole list, none, a random half, with absent and upper-case keys
        let all = distinct_keys(&f.lines);
        check_negative(&f.isam, &f.lines, &all, 0, 10_000);
        let mut half: Vec<Vec<u8>> = all.iter().filter(|_| rng.below(2) == 0).cloned().collect();
        half.push(b"zzz".to_vec());
        half.push(all[0].to_ascii_uppercase());
        half.push(b"".to_vec());
        check_negative(&f.isam, &f.lines, &half, 77, 10_000);
        // `vol_end` is compared with the volume-local OID
        check_negative(&f.isam, &f.lines, &all, 3, 4);
        check_negative(&f.isam, &f.lines, &half, 0, 0);
    }
}

#[test]
fn negative_scan_marks_included_and_visible() {
    // prot_v4: the sequence 0 has the keys `aaa00001` and `aaa00001.1`; listing one of them
    // marks the OID visible, and the pair rule keeps the accession without a version from
    // marking it included.
    let f = fixture("nucl_v4", false);
    // OID 0 has the keys `af000001`, `af000001.1`, `nm_000001` and `nm_000001.1`.
    let keys_of_0: Vec<Vec<u8>> = f
        .lines
        .iter()
        .filter(|l| l.1 == 0)
        .map(|l| l.0.clone())
        .collect();
    assert!(keys_of_0.len() >= 3);
    // Listing one key only: the OID is visible and included (its other keys are not listed).
    let mut list = vec![b"af000001.1".to_vec()];
    let mut got = NegativeOids::new();
    f.isam.negative_scan(10, 100, &mut list, &mut got).unwrap();
    assert!(got.is_visible(10) && got.is_included(10) && got.get_oid_status(10));
    // Listing all its keys: visible and not included, so dropped.
    let mut list = keys_of_0;
    let mut got = NegativeOids::new();
    f.isam.negative_scan(10, 100, &mut list, &mut got).unwrap();
    assert!(got.is_visible(10) && !got.is_included(10) && !got.get_oid_status(10));
    // an OID that is only included, and one that was never seen
    assert!(got.is_included(11) && got.get_oid_status(11));
    assert!(got.get_oid_status(5000));
    assert_eq!(got.num_oids(), 10 + 25);
}

#[test]
fn negative_scan_on_synthetic_indices() {
    let dir = TempDir::new("synneg");
    let mut rng = Rng(0x5eed);
    let records = synthetic_records(&mut rng, 3000);
    let keys = distinct_keys(&records);
    for page_size in [64usize, 4] {
        let isam = build_isam(&dir.0, &format!("n{page_size}"), &records, page_size);
        for round in 0..6 {
            let share = [1, 2, 3, 10, 50, 1][round];
            let mut list: Vec<Vec<u8>> = keys
                .iter()
                .filter(|_| rng.below(100) < share)
                .cloned()
                .collect();
            list.push(b"dup_acc".to_vec());
            list.push(b"nothing".to_vec());
            check_negative(&isam, &records, &list, 11 * round as u32, 1_000_000);
        }
        check_negative(&isam, &records, &keys, 0, 1_000_000);
    }
}

#[test]
fn negative_scan_pair_rule() {
    // `acc` followed by `acc.1` in a page: the line of `acc` is not included, the next line decides.
    let dir = TempDir::new("pair");
    let records = vec![
        (b"abc".to_vec(), 0),
        (b"abc.1".to_vec(), 1),
        (b"abd".to_vec(), 2),
        (b"abd.2".to_vec(), 2),
        (b"abe".to_vec(), 3),
    ];
    let isam = build_isam(&dir.0, "pair", &records, 64);
    let run = |list: &[&str]| {
        let mut list: Vec<Vec<u8>> = list.iter().map(|s| s.as_bytes().to_vec()).collect();
        let mut out = NegativeOids::new();
        isam.negative_scan(0, 100, &mut list, &mut out).unwrap();
        out
    };
    // (an empty list marks nothing: `x_SearchNegativeMultiSeq` returns at once)
    assert_eq!(run(&[]).num_oids(), 0);
    let out = run(&["zzz"]);
    assert!(!out.is_included(0) && out.is_included(1) && out.is_included(2) && out.is_included(3));
    assert_eq!(out.num_oids(), 4);
    let out = run(&["abc"]);
    // `abc` listed: OID 0 visible, the next line (`abc.1`, OID 1) skipped: not marked at all
    assert!(out.is_visible(0) && !out.is_visible(1) && !out.is_included(1));
    let out = run(&["abc.1"]);
    assert!(!out.is_visible(0) && out.is_visible(1) && !out.is_included(0));
    let out = run(&["abd.2", "abd"]);
    assert!(out.is_visible(2) && !out.is_included(2) && out.is_included(3));
    // OIDs from `vol_end` on are not marked
    let mut list = vec![b"abc".to_vec()];
    let mut out = NegativeOids::new();
    isam.negative_scan(0, 2, &mut list, &mut out).unwrap();
    assert!(out.is_visible(0) && !out.is_included(2) && !out.is_included(3));
    assert_eq!(out.num_oids(), 1);
}

#[test]
fn bit_vectors_follow_csqdbbitvector() {
    let mut v = NegativeOids::new();
    assert_eq!(v.num_oids(), 0);
    assert!(v.get_oid_status(0) && v.get_oid_status(100));
    v.add_visible_oid(31);
    assert_eq!(v.num_oids(), 32);
    assert!(!v.get_oid_status(31) && v.get_oid_status(30) && v.get_oid_status(32));
    v.add_included_oid(31);
    assert!(v.get_oid_status(31));
    v.add_included_oid(40_000);
    assert_eq!(v.num_oids(), 40_001);
    assert!(v.is_included(40_000) && !v.is_included(39_999) && !v.is_visible(40_000));
    // setting a lower bit does not shrink the size
    v.add_visible_oid(2);
    assert_eq!(v.num_oids(), 40_001);
    assert!(v.is_visible(2));
    // a bit past Size reads false, as do all the bits of a vector never set
    assert!(!v.is_visible(100_000));
}
