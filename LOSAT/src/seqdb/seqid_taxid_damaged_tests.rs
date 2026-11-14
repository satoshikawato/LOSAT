//! Damaged copies of the files that map ids and taxids to OIDs, and of the taxonomy name files.
//!
//! The first three tests build the damage that once drove a reader to allocate or loop far past
//! what the file holds: the volume tables of the LMDB file of version 5 (`tax_files::volumes`),
//! the string ISAM index of version 4 (`StringIsam::load_data`) and the taxonomy SQLite file
//! (`TaxonomySqlite::open_named`). Each case was a process abort (a failed allocation) or a loop
//! that did not end before the fix; now each is an error.
//!
//! The sweeps after them damage, in the form of `damaged_files_tests.rs` (whose generator, cases
//! and temporary directory they share), the string (`.?si`, `.?sd`) and numeric (`.?ni`, `.?nd`)
//! ISAM files of version 4, the lookup files of two LMDB databases of version 5 (`.?db`, `.?tf`,
//! `.?to`, `.?ot`, `.?os`) and the taxonomy name files (`taxdb.bti`, `taxdb.btd`). After each damage
//! the files are read through the lookups of `-seqidlist`, `-negative_seqidlist`, `-taxids`,
//! `-negative_taxids`, `blastdbcmd -entry` and `-entry_batch`, the search's lookup of the reported
//! Seq-ids and the taxonomy columns: directly and, where a database is needed, through [`SeqDb`].
//! Any result but a panic or a hang is accepted. Every value read (header fields included) and
//! every error text goes into a digest; each target must have seeded rounds in which every reader
//! opened and the digest differs from the intact copy's (the damage reached a reader that did not
//! refuse it at opening), and the restored files must read as the intact ones did. The loops over
//! OIDs and ordinal ids stop at `MAX_OIDS`. The taxonomy SQLite file is swept by
//! `damaged_files_never_panic` (`tax_sqlite.rs`), and the LMDB pages themselves by
//! `corrupt_files_give_errors_not_panics` (`lmdb_read.rs`).

use std::collections::BTreeSet;
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};

use super::alias::SearchPath;
use super::damaged_files_tests::{cases, Rng, TempDir};
use super::filter::{DbList, SeqidList, TaxidList};
use super::isam::{seqdb_simplify_accession, NegativeOids, SiOid, StringIsam};
use super::isam_num::NumericIsam;
use super::lmdb_read::Env;
use super::seqidlist::{parse_list_id, read_si_list};
use super::tax_files::{self, OidSeqIds, OidTaxIds, TaxError};
use super::tax_sqlite::TaxonomySqlite;
use super::taxdb::TaxDb;
use super::{SeqDb, SeqDbError, SeqType};

fn lmdb_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/lmdb")
}

fn ids_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids")
}

fn taxsqlite_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxsqlite")
}

/// Sets `md_entries` (offset 32 of the 48-byte `MDB_db` record) of the named database `name`
/// at every match of its node (`flags, key size, key, record`) in `file`, live or stale.
fn set_entries(file: &mut [u8], name: &[u8], entries: u64) -> usize {
    let mut node = vec![2u8, 0];
    node.extend_from_slice(&(name.len() as u16).to_le_bytes());
    node.extend_from_slice(name);
    let mut found = 0;
    let mut at = 0;
    while let Some(i) = file[at..]
        .windows(node.len())
        .position(|w| w == node.as_slice())
    {
        let record = at + i + node.len();
        file[record + 32..record + 40].copy_from_slice(&entries.to_le_bytes());
        found += 1;
        at = record;
    }
    found
}

/// `md_entries` of both volume tables set to a count that the file cannot hold. Before the fix the
/// count of 2^28 ran 84 s up to 4.2 GB, then aborted on an allocation of 8 GiB (under an 8 GB
/// address-space limit); 2^32 - 1 entries would ask for about 128 GiB.
#[test]
fn volume_tables_with_more_entries_than_the_file_holds_are_errors() {
    for name in ["p5.pdb", "mv.pdb"] {
        let mut file = std::fs::read(lmdb_dir().join(name)).unwrap();
        let intact = Env::open(&file).unwrap();
        let volumes = tax_files::volumes(&intact).unwrap();
        assert!(!volumes.is_empty());
        for entries in [1u64 << 28, u64::from(u32::MAX)] {
            assert!(set_entries(&mut file, b"volname", entries) > 0, "{name}");
            assert!(set_entries(&mut file, b"volinfo", entries) > 0, "{name}");
            let env = Env::open(&file).unwrap();
            let error = tax_files::volumes(&env).unwrap_err();
            assert_eq!(
                error,
                format!("Volinfo error : {entries} volumes"),
                "{name}"
            );
            assert!(!tax_files::is_ncbi_volume_error(&error), "{name}: {error}");
        }
    }
}

/// `num_terms` and `page_size` of a string ISAM index rewritten to numbers that pass the header
/// checks (`(samples - 1) * page < terms <= samples * page`) but that the data file cannot hold.
/// Before the fix the first page read asked for 11.4 GB at once (an abort under an 8 GB
/// address-space limit).
#[test]
fn string_isam_pages_larger_than_the_data_file_are_errors() {
    let dir = TempDir::new("ids-isam");
    let index = std::fs::read(ids_dir().join("prot_v4.psi")).unwrap();
    let samples = u32::from_be_bytes(index[16..20].try_into().unwrap());
    assert_eq!(samples, 3, "the fixture's sample count");
    let page: u32 = 715_000_000;
    let terms: u32 = 2 * page + 1;
    let mut damaged = index.clone();
    damaged[12..16].copy_from_slice(&terms.to_be_bytes());
    damaged[20..24].copy_from_slice(&page.to_be_bytes());
    let (idx, dat) = (dir.0.join("d.psi"), dir.0.join("d.psd"));
    std::fs::write(&idx, &damaged).unwrap();
    std::fs::copy(ids_dir().join("prot_v4.psd"), &dat).unwrap();
    let isam = StringIsam::open(&idx, &dat).unwrap();
    // The batch functions read a whole page of `page_size` records (`load_data`).
    let mut ids = [
        SiOid::new(b"xp_000011.2".to_vec()),
        SiOid::new(b"0008123a".to_vec()),
    ];
    let error = isam.ids_to_oids(0, &mut ids).unwrap_err().to_string();
    assert!(error.contains("string ISAM data page"), "{error}");
    let mut negative = NegativeOids::new();
    let mut list = vec![b"xp_000011.2".to_vec()];
    let error = isam
        .negative_scan(0, 10, &mut list, &mut negative)
        .unwrap_err()
        .to_string();
    assert!(error.contains("string ISAM data page"), "{error}");
    // The intact index gives the OID of the same key.
    let intact = StringIsam::open(
        &ids_dir().join("prot_v4.psi"),
        &ids_dir().join("prot_v4.psd"),
    )
    .unwrap();
    let mut ids = [SiOid::new(b"xp_000011.2".to_vec())];
    intact.ids_to_oids(0, &mut ids).unwrap();
    assert_eq!(ids[0].oid, Some(11));
}

/// A taxonomy SQLite file whose `sqlite_master` tree is a chain of interior pages, each pointing
/// at the next page from every cell and from its right-most pointer (a tree of 58 * 72^5 leaves
/// on 7 pages), and whose header claims 2^32 - 1 pages.
fn branching_sqlite() -> Vec<u8> {
    let mut b = std::fs::read(taxsqlite_dir().join("p512.sqlite3")).unwrap();
    let ps = 512usize;
    assert_eq!(u16::from_be_bytes([b[16], b[17]]) as usize, ps);
    b[28..32].copy_from_slice(&u32::MAX.to_be_bytes());
    let counter: [u8; 4] = b[24..28].try_into().unwrap();
    b[92..96].copy_from_slice(&counter);
    for pg in 1..=7usize {
        let base = (pg - 1) * ps;
        let h = if pg == 1 { 100 } else { 0 };
        b[base + h..base + ps].fill(0);
        let page = &mut b[base..base + ps];
        if pg == 7 {
            page[h] = 0x0d;
            page[h + 5..h + 7].copy_from_slice(&(ps as u16).to_be_bytes());
            continue;
        }
        page[h] = 0x05;
        let k = (ps - h - 12) / 7;
        let mut content = ps;
        for i in 0..k {
            content -= 5;
            page[content..content + 4].copy_from_slice(&((pg + 1) as u32).to_be_bytes());
            page[content + 4] = 1;
            let ptr = h + 12 + 2 * i;
            page[ptr..ptr + 2].copy_from_slice(&(content as u16).to_be_bytes());
        }
        page[h + 3..h + 5].copy_from_slice(&(k as u16).to_be_bytes());
        page[h + 5..h + 7].copy_from_slice(&(content as u16).to_be_bytes());
        page[h + 8..h + 12].copy_from_slice(&((pg + 1) as u32).to_be_bytes());
    }
    b
}

/// Before the fix the walk of `sqlite_master` did not end: the tree is not cyclic, so the depth
/// limit is never reached, and the page budget (the header's page count plus 64) was 2^32 + 63.
/// SQLite refuses the file (`database disk image is malformed`).
#[test]
fn sqlite_header_page_count_past_the_file_is_an_error() {
    let file = branching_sqlite();
    let error = TaxonomySqlite::open(&file).err().expect("must not open");
    assert!(error.contains("malformed"), "{error}");
    // The header's count alone, as SQLite 3.51.1 reads it: 248 pages on a file of 247 pages is
    // malformed; on 247 pages and 100 bytes (248 pages rounded up) the file opens.
    let mut file = std::fs::read(taxsqlite_dir().join("p512.sqlite3")).unwrap();
    assert_eq!(file.len(), 247 * 512);
    let counter: [u8; 4] = file[24..28].try_into().unwrap();
    file[92..96].copy_from_slice(&counter);
    file[28..32].copy_from_slice(&247u32.to_be_bytes());
    assert!(TaxonomySqlite::open(&file).is_ok());
    file[28..32].copy_from_slice(&248u32.to_be_bytes());
    let error = TaxonomySqlite::open(&file).err().expect("248 pages on 247");
    assert!(error.contains("malformed"), "{error}");
    file.extend_from_slice(&[0; 100]);
    assert!(TaxonomySqlite::open(&file).is_ok());
}

/// The most OIDs read from a damaged OID-to-taxids or OID-to-seqids file (besides its last OID,
/// the one past it and `u32::MAX`), and the most ordinal ids of a restricted database listed. The
/// fixtures hold at most 14 OIDs (`tx`).
const MAX_OIDS: u64 = 16;

/// The text of an error, with the forms in which LOSAT reports it.
trait ErrorText {
    fn text(&self) -> String;
}

impl ErrorText for SeqDbError {
    fn text(&self) -> String {
        format!(
            "{self} {:?} {:?}",
            self.ncbi_message(),
            self.ncbi_search_message()
        )
    }
}

impl ErrorText for String {
    fn text(&self) -> String {
        self.clone()
    }
}

impl ErrorText for TaxError {
    fn text(&self) -> String {
        format!("{self} {} {}", self.message("db"), self.detail())
    }
}

/// What the readers of one case gave: how many opened, and a digest (`DefaultHasher::new`, whose
/// keys are fixed) of every value read and every error text.
#[derive(Default)]
struct Trace {
    opened: usize,
    digest: DefaultHasher,
}

impl Trace {
    /// An opened reader, or `None` after adding the error to the digest.
    fn open<T, E: ErrorText>(&mut self, result: Result<T, E>) -> Option<T> {
        match result {
            Ok(value) => {
                self.opened += 1;
                Some(value)
            }
            Err(error) => {
                error.text().hash(&mut self.digest);
                None
            }
        }
    }

    /// A reader that may be absent (`Ok(None)`, which is not counted as opened).
    fn open_some<T, E: ErrorText>(&mut self, result: Result<Option<T>, E>) -> Option<T> {
        match result {
            Ok(None) => {
                "absent".hash(&mut self.digest);
                None
            }
            result => self.open(result).flatten(),
        }
    }

    /// A value read, or the error, added to the digest.
    fn read<T: Debug, E: ErrorText>(&mut self, result: Result<T, E>) {
        match result {
            Ok(value) => self.note(&value),
            Err(error) => error.text().hash(&mut self.digest),
        }
    }

    fn note<T: Debug + ?Sized>(&mut self, value: &T) {
        format!("{value:?}").hash(&mut self.digest);
    }

    fn finish(&self) -> (usize, u64) {
        (self.opened, self.digest.finish())
    }
}

/// Damages each file of `targets` in `dir` in turn (the others intact) and runs `read` on the
/// directory after each damage. Fails listing every case that panicked, every target with no
/// seeded round in which all `opens` readers opened and gave something other than the intact
/// files, and every target whose restored file reads otherwise than the intact one.
fn sweep(
    tag: &str,
    dir: &TempDir,
    targets: &[&str],
    opens: usize,
    read: impl Fn(&Path, &mut Trace),
) {
    let read_all = || {
        let mut trace = Trace::default();
        read(&dir.0, &mut trace);
        trace.finish()
    };
    let intact = read_all();
    assert_eq!(
        intact.0, opens,
        "{tag}: the readers of the intact files must open"
    );
    let mut problems = Vec::new();
    for target in targets {
        let path = dir.0.join(target);
        let orig = std::fs::read(&path).unwrap();
        let mut rng = Rng::seeded(&format!("{tag}/{target}"));
        let (mut seeded_opened, mut seeded_reached) = (0, 0);
        for case in cases(&orig, &mut rng) {
            std::fs::write(&path, &case.bytes).unwrap();
            match catch_unwind(AssertUnwindSafe(&read_all)) {
                Err(_) => problems.push(format!("{target}: {} panicked", case.describe())),
                Ok((opened, digest)) if case.round.is_some() && opened == opens => {
                    seeded_opened += 1;
                    seeded_reached += usize::from(digest != intact.1);
                }
                Ok(_) => {}
            }
        }
        std::fs::write(&path, &orig).unwrap();
        if seeded_reached == 0 {
            problems.push(format!(
                "{target}: the seeded damage did not reach the readers past their opening \
                 (opened {seeded_opened}, read otherwise {seeded_reached})"
            ));
        }
        if read_all() != intact {
            problems.push(format!("{target}: the restored file reads otherwise"));
        }
    }
    assert!(
        problems.is_empty(),
        "{tag}: {} problems:\n{}",
        problems.len(),
        problems.join("\n")
    );
}

/// The ids of the lines of a string ISAM data file (`key 0x02 oid 0x0a`): about eight spread over
/// the file and the last, then ids before the first and after the last line, an accession with a
/// version the file does not hold, and one (of `prot_v4`) with and without its version.
fn string_keys(data: &[u8]) -> Vec<Vec<u8>> {
    let lines: Vec<&[u8]> = data
        .split(|&byte| byte == b'\n')
        .filter_map(|line| line.split(|&byte| byte == 2).next())
        .filter(|key| !key.is_empty())
        .collect();
    let step = lines.len().div_ceil(8).max(1);
    let mut keys: Vec<Vec<u8>> = lines.iter().step_by(step).map(|key| key.to_vec()).collect();
    keys.push(lines.last().unwrap().to_vec());
    for key in [
        &b"0"[..],
        b"zzzz",
        b"aaa00001.9",
        b"xp_000011.2",
        b"xp_000011",
    ] {
        keys.push(key.to_vec());
    }
    keys
}

/// The GIs of a numeric ISAM data file of 4-byte keys (big-endian pairs of an unsigned GI and an
/// OID).
fn file_gis(data: &[u8]) -> Vec<i64> {
    data.chunks_exact(8)
        .map(|pair| i64::from(u32::from_be_bytes(pair[..4].try_into().unwrap())))
        .collect()
}

/// What the ISAM sweep of a database looks up.
struct IsamQueries {
    /// Ids of the string ISAM ([`string_keys`]).
    keys: Vec<Vec<u8>>,
    /// The GIs of the numeric ISAM, each with its neighbours, and numbers outside every file.
    gis: Vec<i64>,
    /// Up to 5 of the `keys`, two ids that neither file holds and 4 of the file's GIs (as `gi|N`),
    /// of those that parse as the ids of a seqid list.
    list: Vec<Vec<u8>>,
    /// The OID count of the volume, for the negative scan.
    vol_end: u32,
}

impl IsamQueries {
    fn new(keys: Vec<Vec<u8>>, file_gis: &[i64], vol_end: u32) -> Self {
        let mut gis = vec![0, -5, i64::from(u32::MAX) + 1, i64::MAX];
        gis.extend(file_gis.iter().flat_map(|&gi| [gi - 1, gi, gi + 1]));
        let parses = |id: &Vec<u8>| parse_list_id(id).is_ok();
        let mut list: Vec<Vec<u8>> = keys
            .iter()
            .filter(|id| parses(id))
            .take(5)
            .cloned()
            .collect();
        let absent = [&b"zzzz"[..], b"aaa00001.9"].map(<[u8]>::to_vec);
        list.extend(absent.into_iter().filter(parses));
        let gi_ids = file_gis.iter().map(|gi| format!("gi|{gi}").into_bytes());
        list.extend(gi_ids.filter(parses).take(4));
        Self {
            keys,
            gis,
            list,
            vol_end,
        }
    }
}

/// The string and numeric ISAMs of the version 4 database `stem` in `dir`, read directly and
/// through the lookups of [`SeqDb`]: `blastdbcmd -entry` with a GI or an accession, the search's
/// lookup of the reported Seq-ids, and `-seqidlist` and `-negative_seqidlist` (`restrict`).
fn read_isams(dir: &Path, stem: &str, seqtype: SeqType, queries: &IsamQueries, trace: &mut Trace) {
    let IsamQueries {
        keys,
        gis,
        list,
        vol_end,
    } = queries;
    let letter = match seqtype {
        SeqType::Protein => 'p',
        SeqType::Nucleotide => 'n',
    };
    let file = |ext: &str| dir.join(format!("{stem}.{letter}{ext}"));
    if let Some(isam) = trace.open(StringIsam::open(&file("si"), &file("sd"))) {
        trace.note(&isam.id_bounds());
        for key in keys {
            trace.read(isam.string_search(key));
            trace.read(isam.string_to_oids(key, false, true));
            trace.read(isam.string_to_oids(key, true, false));
        }
        let mut ids: Vec<SiOid> = keys.iter().map(|key| SiOid::new(key.clone())).collect();
        trace.read(isam.ids_to_oids(0, &mut ids));
        trace.note(&ids);
        let mut negative_ids = keys.to_vec();
        let mut negative = NegativeOids::new();
        trace.read(isam.negative_scan(0, *vol_end, &mut negative_ids, &mut negative));
        trace.note(&negative);
        // The search programs look ids up in the index kept in memory (`lookup.rs`).
        if let Some(kept) = trace.open(isam.keep_index()) {
            for key in keys {
                trace.read(kept.string_search(key));
            }
        }
    }
    if let Some(isam) = trace.open_some(NumericIsam::open_volume(&dir.join(stem), seqtype)) {
        for &gi in gis {
            trace.read(isam.id_to_oid(gi));
        }
    }
    let name = dir.join(stem);
    let Some(db) = trace.open(SeqDb::open(name.to_str().unwrap(), seqtype)) else {
        return;
    };
    let Some(lookup) = trace.open(db.id_lookup()) else {
        return;
    };
    for &gi in gis {
        trace.read(lookup.gi_to_oid(gi));
    }
    for key in keys {
        if let Ok(simplified) = seqdb_simplify_accession(key) {
            trace.read(lookup.accession_to_oids_v4(key, &simplified));
        }
    }
    read_seqid_lookups(&db, list, trace);
    read_seqid_lists(name.to_str().unwrap(), seqtype, list, trace);
}

/// The search's lookup of the Seq-ids of the reported sequences (`seqid_to_oids`).
fn read_seqid_lookups(db: &SeqDb, ids: &[Vec<u8>], trace: &mut Trace) {
    let Some(lookup) = trace.open(db.seqid_lookup()) else {
        return;
    };
    for id in ids {
        if let Ok(id) = parse_list_id(id) {
            trace.read(lookup.seqid_to_oids(&id, true));
        }
    }
}

/// `-seqidlist` and `-negative_seqidlist`: the database `name` restricted to a text list of `ids`,
/// and the OIDs of its first ordinal ids.
fn read_seqid_lists(name: &str, seqtype: SeqType, ids: &[Vec<u8>], trace: &mut Trace) {
    let mut text = Vec::new();
    for id in ids {
        text.extend_from_slice(id);
        text.push(b'\n');
    }
    for negative in [false, true] {
        let mut warnings = Vec::new();
        let file = read_si_list(&text, &mut warnings).unwrap();
        let list = DbList::Seqids(SeqidList { negative, file });
        let Ok(mut db) = SeqDb::open(name, seqtype) else {
            continue;
        };
        let mut warnings = Vec::new();
        trace.read(db.restrict(&list, &mut warnings));
        trace.note(&warnings);
        let ordinals = db.num_seqs().min(MAX_OIDS as u32);
        let oids: Vec<u32> = (0..ordinals).map(|ordinal| db.oid(ordinal).get()).collect();
        trace.note(&oids);
    }
}

fn sweep_isams(stem: &str, seqtype: SeqType) {
    let dir = TempDir::new(&format!("ids-{stem}"));
    dir.copy_fixtures(&ids_dir(), &[stem]);
    let letter = match seqtype {
        SeqType::Protein => 'p',
        SeqType::Nucleotide => 'n',
    };
    let ext = |ext: &str| format!("{stem}.{letter}{ext}");
    let queries = IsamQueries::new(
        string_keys(&std::fs::read(dir.0.join(ext("sd"))).unwrap()),
        &file_gis(&std::fs::read(dir.0.join(ext("nd"))).unwrap()),
        SeqDb::open(&dir.name(stem), seqtype).unwrap().num_seqs(),
    );
    assert!(
        queries.list.len() >= 6,
        "{stem}: list ids {:?}",
        queries.list
    );
    let targets = [ext("si"), ext("sd"), ext("ni"), ext("nd")];
    let targets: Vec<&str> = targets.iter().map(String::as_str).collect();
    // The string ISAM, its kept index, the numeric ISAM, the database and its two lookups.
    sweep(stem, &dir, &targets, 6, |dir, trace| {
        read_isams(dir, stem, seqtype, &queries, trace)
    });
}

#[test]
fn damaged_isam_files_never_panic() {
    sweep_isams("prot_v4", SeqType::Protein);
    sweep_isams("nucl_v4", SeqType::Nucleotide);
}

/// The accessions of `p5.pdb` (`tests/seqdb/lmdb/expected/p5.pdb.txt`) with and without their
/// versions, a key of other case and keys it does not hold.
const P5_ACCESSIONS: [&[u8]; 12] = [
    b"AAA00001",
    b"AAA00001.1",
    b"AAA00001.2",
    b"BBB00001.1",
    b"CCC00001",
    b"myprot9",
    b"MYPROT9",
    b"NP_000002",
    b"P12345.2",
    b"XP_000009.1",
    b"aaa00001",
    b"ZZZ00001",
];

/// Sets of the taxids of `p5` (`expected/p5.tax.txt`: 0, 7227, 9606, 10090), with a taxid it
/// does not hold.
fn p5_tax_sets() -> Vec<BTreeSet<i32>> {
    [
        &[9606][..],
        &[10090],
        &[7227, 9606],
        &[0],
        &[1],
        &[0, 7227, 9606, 10090],
    ]
    .iter()
    .map(|set| set.iter().copied().collect())
    .collect()
}

/// The first `MAX_OIDS` OIDs of a file that claims `num_oids`, its last OID, the one past it and
/// `u32::MAX`.
fn oids_to_read(num_oids: u64) -> impl Iterator<Item = u32> {
    let past = u32::try_from(num_oids).unwrap_or(u32::MAX);
    (0..num_oids.min(MAX_OIDS) as u32).chain([past.saturating_sub(1), past, u32::MAX])
}

/// The lookup files of the LMDB database `p5`, read by the `tax_files` functions on bytes:
/// accessions to OIDs (`.pdb`), taxids to OIDs (`.ptf`, `.pto`, and `.pot` for a negative list),
/// OIDs to taxids (`.pot`) and OIDs to seqids (`.pos`). `num_oids` is the intact count (a damaged
/// volume table can claim 2^32 OIDs, whose bitmap `get_oids_for_tax_ids` would allocate).
fn read_lmdb_files(dir: &Path, num_oids: u32, trace: &mut Trace) {
    let file = |ext: &str| std::fs::read(dir.join(format!("p5.{ext}"))).unwrap();
    let (pdb, ptf, pto, pot, pos) = (
        file("pdb"),
        file("ptf"),
        file("pto"),
        file("pot"),
        file("pos"),
    );
    if let Some(env) = trace.open(Env::open(&pdb)) {
        trace.read(tax_files::volumes(&env));
        for acc in P5_ACCESSIONS {
            trace.read(tax_files::acc_first_oid(&env, acc));
            trace.read(tax_files::acc_oids(&env, acc));
        }
    }
    if let Some(ntf) = trace.open(Env::open(&ptf)) {
        trace.read(tax_files::get_db_tax_ids(&ntf));
        for set in p5_tax_sets() {
            trace.read(tax_files::get_oids_for_tax_ids(&ntf, &pto, num_oids, &set));
            trace.read(tax_files::negative_tax_ids_to_oids(
                &ntf, &pto, &pot, num_oids, &set,
            ));
        }
    }
    if let Some(lookup) = trace.open(OidTaxIds::open(&pot)) {
        trace.note(&lookup.num_oids());
        for oid in oids_to_read(lookup.num_oids()) {
            trace.read(lookup.taxids(oid));
        }
    }
    if let Some(lookup) = trace.open(OidSeqIds::open(&pos)) {
        trace.note(&lookup.num_oids());
        for oid in oids_to_read(lookup.num_oids()) {
            trace.read(lookup.seqids(oid));
        }
    }
}

#[test]
fn damaged_lmdb_lookup_files_never_panic() {
    let dir = TempDir::new("lmdb-p5");
    dir.copy_fixtures(&lmdb_dir(), &["p5"]);
    let pdb = std::fs::read(dir.0.join("p5.pdb")).unwrap();
    let volumes = tax_files::volumes(&Env::open(&pdb).unwrap()).unwrap();
    let num_oids = u32::try_from(tax_files::total_oids(&volumes)).unwrap();
    assert_eq!(num_oids, 7, "the OIDs of p5");
    let targets = ["p5.pdb", "p5.ptf", "p5.pto", "p5.pot", "p5.pos"];
    // The two LMDB environments and the two OID files.
    sweep("p5", &dir, &targets, 4, |dir, trace| {
        read_lmdb_files(dir, num_oids, trace)
    });
}

/// The accessions of `tx` (`tests/seqdb/lmdb/tx.map`), with and without versions, and one it
/// does not hold.
const TX_ACCESSIONS: [&[u8]; 5] = [
    b"AAA00001.1",
    b"AAA00002",
    b"BBB00001.1",
    b"aaa00003.1",
    b"ZZZ00001",
];

/// The taxids of `tx` (100 to 500), with one it does not hold.
fn tx_tax_sets() -> Vec<BTreeSet<i32>> {
    [&[100][..], &[200], &[100, 200], &[500], &[999]]
        .iter()
        .map(|set| set.iter().copied().collect())
        .collect()
}

/// The database `tx` (version 5) opened by [`SeqDb`] and read through the lookups of
/// `blastdbcmd -entry` and `-entry_batch` (`accession_to_oids`, `accessions_to_oids`),
/// `blastdbcmd -tax_info` (`db_tax_ids`), the search's lookup of the reported Seq-ids,
/// `-seqidlist` and `-negative_seqidlist` (`restrict`, which reads `tx.pos` for a negative list),
/// and `-taxids` and `-negative_taxids` (`restrict_taxids`).
fn read_tx(dir: &Path, trace: &mut Trace) {
    let name = dir.join("tx");
    let name = name.to_str().unwrap();
    let ids: Vec<Vec<u8>> = TX_ACCESSIONS.iter().map(|acc| acc.to_vec()).collect();
    if let Some(db) = trace.open(SeqDb::open(name, SeqType::Protein)) {
        if let Some(lookup) = trace.open(db.id_lookup()) {
            for acc in TX_ACCESSIONS {
                trace.read(lookup.accession_to_oids(acc));
            }
            let mut oids = Vec::new();
            trace.read(lookup.accessions_to_oids(&ids, &mut oids));
            trace.note(&oids);
        }
        if let Some(lookup) = trace.open(db.tax_id_lookup()) {
            trace.read(lookup.db_tax_ids());
            for set in tx_tax_sets() {
                trace.read(lookup.tax_ids_to_oids(&set));
            }
        }
        read_seqid_lookups(&db, &ids, trace);
    }
    read_seqid_lists(name, SeqType::Protein, &ids, trace);
    for set in [&[100][..], &[100, 200], &[999]] {
        let set: BTreeSet<i32> = set.iter().copied().collect();
        for list in [TaxidList::Positive(set.clone()), TaxidList::Negative(set)] {
            let Ok(mut db) = SeqDb::open(name, SeqType::Protein) else {
                continue;
            };
            trace.read(db.restrict_taxids(&list));
            let ordinals = db.num_seqs().min(MAX_OIDS as u32);
            let oids: Vec<u32> = (0..ordinals).map(|ordinal| db.oid(ordinal).get()).collect();
            trace.note(&oids);
        }
    }
}

#[test]
fn damaged_lookup_files_of_a_database_never_panic() {
    let dir = TempDir::new("lmdb-tx");
    dir.copy_fixtures(&lmdb_dir(), &["tx"]);
    let targets = ["tx.pdb", "tx.ptf", "tx.pto", "tx.pot", "tx.pos"];
    // The database and its three lookups.
    sweep("tx", &dir, &targets, 4, read_tx);
}

/// The taxids of `taxdb/good` (`taxdb.bti`) and taxids it does not hold.
const TAXDB_TAXIDS: [i32; 11] = [
    562,
    7777,
    9606,
    10090,
    12908,
    32630,
    99990,
    0,
    1,
    -1,
    i32::MAX,
];

/// The taxonomy name files on the search path, read by the lookups of the taxonomy columns
/// (`names`) and of the checks before a search (`offset_error_taxid`, `any_lookup_may_*`).
fn read_taxdb(dir: &Path, trace: &mut Trace) {
    let search = SearchPath::new(dir, Some(dir.as_os_str()), None).unwrap();
    match TaxDb::open(&search) {
        // A database that NCBI treats as missing.
        Ok(db) if db.is_missing() => trace.note(db.open_messages()),
        result => {
            let Some(db) = trace.open(result) else {
                return;
            };
            trace.note(db.open_messages());
            trace.note(&db.entry_count());
            trace.read(db.offset_error_taxid());
            trace.read(db.any_lookup_may_fail());
            trace.read(db.any_lookup_may_refuse(&|name: &[u8]| name.contains(&b'\t')));
            for tax_id in TAXDB_TAXIDS {
                let mut messages = Vec::new();
                trace.read(db.names(tax_id, &mut messages));
                trace.note(&messages);
            }
        }
    }
}

#[test]
fn damaged_taxdb_files_never_panic() {
    let dir = TempDir::new("taxdb");
    dir.copy_fixtures(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/good"),
        &["taxdb"],
    );
    sweep("taxdb", &dir, &["taxdb.bti", "taxdb.btd"], 1, read_taxdb);
}
