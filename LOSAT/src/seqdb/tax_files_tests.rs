use super::*;
use std::fs::File;

/// Fixture databases (made by NCBI makeblastdb 2.17.0, see `tests/seqdb/README.md`).
const FX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/seqdb/lmdb");
/// Values written by `dump.py` (a Python reference reader) and by NCBI programs.
const EXP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/seqdb/lmdb/expected");

/// Path prefix of the NCBI-distributed `18S_fungal_sequences` database
/// (`LOSAT_TEST_18S_DB`, without extension) for the `#[ignore]` tests.
fn s18() -> String {
    std::env::var("LOSAT_TEST_18S_DB").expect("LOSAT_TEST_18S_DB")
}

/// Path of an `expected/<name>.txt` file: the repository fixtures, else the directory named by
/// `LOSAT_TEST_18S_EXPECTED` (the dumps of the 18S database are too big to commit).
fn exp_path(name: &str) -> String {
    let repo = format!("{EXP}/{name}.txt");
    if std::path::Path::new(&repo).exists() {
        return repo;
    }
    let dir = std::env::var("LOSAT_TEST_18S_EXPECTED").unwrap_or_default();
    format!("{dir}/{name}.txt")
}

/// Values dumped by `gen/expect_tax.py` with `dump.py`.
#[derive(Default)]
struct Expected {
    vols: Vec<(String, u32)>,
    /// taxid -> OIDs (ascending)
    t: Vec<(i32, Vec<u32>)>,
    /// oid -> taxids
    o: Vec<Vec<i32>>,
    /// oid -> seqids
    s: Vec<Vec<String>>,
}

fn csv<T: std::str::FromStr>(s: &str) -> Vec<T>
where
    T::Err: std::fmt::Debug,
{
    s.split(',')
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().unwrap())
        .collect()
}

fn expected(name: &str) -> Expected {
    let text = std::fs::read_to_string(exp_path(&format!("{name}.tax"))).unwrap();
    let mut e = Expected::default();
    for line in text.lines() {
        let f: Vec<&str> = line.splitn(5, ' ').collect();
        match f[0] {
            "V" => e.vols.push((f[1].to_string(), f[2].parse().unwrap())),
            "T" => e.t.push((f[1].parse().unwrap(), csv(f[3]))),
            "O" => {
                assert_eq!(f[1].parse::<usize>().unwrap(), e.o.len());
                e.o.push(csv(f.get(2).copied().unwrap_or("")));
            }
            "S" => {
                assert_eq!(f[1].parse::<usize>().unwrap(), e.s.len());
                let ids = f.get(2).copied().unwrap_or("");
                e.s.push(
                    ids.split('|')
                        .filter(|x| !x.is_empty())
                        .map(String::from)
                        .collect(),
                );
            }
            "X" => {}
            _ => panic!("bad line {line}"),
        }
    }
    e
}

/// The five files of one LMDB entry.
struct Files {
    ndb: File,
    ntf: File,
    nto: File,
    not: File,
    nos: Option<File>,
}

fn files(prefix: &str, n: char) -> Files {
    let o = |ext: &str| File::open(format!("{prefix}.{n}{ext}")).unwrap();
    Files {
        ndb: o("db"),
        ntf: o("tf"),
        nto: o("to"),
        not: o("ot"),
        nos: File::open(format!("{prefix}.{n}os")).ok(),
    }
}

fn set(v: &[i32]) -> BTreeSet<i32> {
    v.iter().copied().collect()
}

/// (path prefix, p/n, expected-file name)
fn all_dbs() -> Vec<(String, char, &'static str)> {
    vec![
        (format!("{FX}/n5"), 'n', "n5"),
        (format!("{FX}/p5"), 'p', "p5"),
        (format!("{FX}/mv"), 'p', "mv"),
        (format!("{FX}/dup"), 'p', "dup"),
        (format!("{FX}/tx"), 'p', "tx"),
    ]
}

fn check_volumes(prefix: &str, n: char, name: &str) {
    {
        let f = files(prefix, n);
        let env = Env::open(&f.ndb).unwrap();
        let vols = volumes(&env).unwrap();
        let want = expected(name).vols;
        let got: Vec<(String, u32)> = vols
            .iter()
            .map(|(a, b)| (String::from_utf8(a.clone()).unwrap(), *b))
            .collect();
        assert_eq!(got, want, "{name}");
        // the sum is the OID count of the .not file
        let oidtax = OidTaxIds::open(&f.not).unwrap();
        assert_eq!(total_oids(&vols), oidtax.num_oids(), "{name}");
    }
}

#[test]
fn volumes_match_dump_py() {
    for (prefix, n, name) in all_dbs() {
        check_volumes(&prefix, n, name);
    }
    // literals from the notes: mv has 4 volumes 8,8,8,6
    let f = files(&format!("{FX}/mv"), 'p');
    let vols = volumes(&Env::open(&f.ndb).unwrap()).unwrap();
    assert_eq!(
        vols,
        vec![
            (b"mv.00".to_vec(), 8),
            (b"mv.01".to_vec(), 8),
            (b"mv.02".to_vec(), 8),
            (b"mv.03".to_vec(), 6)
        ]
    );
}

#[test]
fn volumes_of_400_volume_db() {
    let f = File::open(format!("{FX}/mvbig.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    let vols = volumes(&env).unwrap();
    assert_eq!(vols.len(), 400);
    for (i, (name, _)) in vols.iter().enumerate() {
        assert_eq!(name, format!("mvbig.{i:03}").as_bytes());
    }
    // counts from dump.py's volinfo dump
    let text = std::fs::read_to_string(exp_path("mvbig.pdb")).unwrap();
    let mut in_info = false;
    let mut want = Vec::new();
    for l in text.lines() {
        if l.starts_with("DB ") {
            in_info = l.starts_with("DB volinfo ");
        } else if in_info && l.starts_with("K ") {
            let v = l.split(' ').nth(2).unwrap();
            want.push(u32::from_le_bytes([
                u8::from_str_radix(&v[0..2], 16).unwrap(),
                u8::from_str_radix(&v[2..4], 16).unwrap(),
                u8::from_str_radix(&v[4..6], 16).unwrap(),
                u8::from_str_radix(&v[6..8], 16).unwrap(),
            ]));
        }
    }
    assert_eq!(vols.iter().map(|v| v.1).collect::<Vec<_>>(), want);
}

#[test]
fn missing_tables_give_ncbi_texts() {
    // a .ptf has no volume or accession table; a .pdb has no taxid table
    let ptf = File::open(format!("{FX}/p5.ptf")).unwrap();
    let env = Env::open(&ptf).unwrap();
    assert_eq!(
        volumes(&env).unwrap_err().to_string(),
        "DB contains no vol info."
    );
    assert_eq!(
        acc_first_oid(&env, b"x").unwrap_err().to_string(),
        "DB contains no accession info."
    );
    assert_eq!(
        acc_oids(&env, b"x").unwrap_err().to_string(),
        NO_ACCESSION_INFO
    );
    let pdb = File::open(format!("{FX}/p5.pdb")).unwrap();
    let env = Env::open(&pdb).unwrap();
    let nto: Vec<u8> = Vec::new();
    let e = get_oids_for_tax_ids(&env, &nto, 7, &set(&[1])).unwrap_err();
    assert!(matches!(e, TaxError::NoTaxInfo), "{e:?}");
    assert_eq!(e.message(b"p5.ptf"), b"No taxonomy info found in p5.ptf");
    assert!(matches!(
        get_db_tax_ids(&env).unwrap_err(),
        TaxError::NoTaxInfo
    ));
    assert_eq!(
        accessions_lookup_error(b"p5.pdb"),
        b"Accessions to Oids lookup error in p5.pdb"
    );
    let e = TaxError::Lookup("detail".into());
    assert_eq!(
        e.message(b"x.ntf"),
        b"Taxonomy Id to Oids lookup error in x.ntf"
    );
    assert_eq!(e.detail(), "detail");
    // volinfo/volname with different counts and a missing volinfo: synthetic LMDB is covered in lmdb_read
}

#[test]
fn accession_lookups() {
    let f = File::open(format!("{FX}/p5.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    assert_eq!(acc_first_oid(&env, b"AAA00001").unwrap(), Some(2));
    assert_eq!(acc_oids(&env, b"AAA00001").unwrap(), vec![2, 3]);
    assert_eq!(acc_first_oid(&env, b"AAA00001.2").unwrap(), Some(3));
    assert_eq!(acc_first_oid(&env, b"aaa00001").unwrap(), None);
    assert_eq!(acc_oids(&env, b"nokey").unwrap(), Vec::<u32>::new());
    // 700 OIDs under one key (sub-database): memcmp order of the LE bytes, NCBI fails here
    let f = File::open(format!("{FX}/dup.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    let oids = acc_oids(&env, b"DUP00001").unwrap();
    assert_eq!(oids.len(), 700);
    assert_eq!(&oids[..7], &[0, 256, 512, 1, 257, 513, 2]);
    assert_eq!(acc_first_oid(&env, b"DUP00001").unwrap(), Some(0));
    // global OIDs over volumes (notes: MVA00017.1 -> 17, volume 2 holds 16..23)
    let f = File::open(format!("{FX}/mv.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    assert_eq!(acc_first_oid(&env, b"MVA00017.1").unwrap(), Some(17));
}

/// `acc_case_variant_oids`: the other keys that differ only in the case of letters (p5's
/// `lcl|myprot9` and the upper-case copy makeblastdb writes, OID 4), none for a key alone in
/// its letters, and none before byte `from`; the answers kept in the map give the same.
#[test]
fn accession_case_variants() {
    let f = File::open(format!("{FX}/p5.pdb")).unwrap();
    let env = Env::open(&f).unwrap();
    // The second pass reads the prefixes that the first kept.
    let mut known = rustc_hash::FxHashMap::default();
    let db = acc2oid(&env).unwrap();
    for _ in 0..2 {
        for (acc, from, want) in [
            (&b"myprot9"[..], 0, vec![4]),
            (b"MYPROT9", 0, vec![4]),
            (b"MyProt9", 0, vec![4, 4]),
            (b"MyProt9", 2, vec![]),
            (b"MyProt9", 1, vec![4]),
            (b"aaa00001", 0, vec![2, 3]),
            (b"AAA00001", 0, vec![]),
            (b"np_000002.1", 0, vec![1]),
            (b"zzz", 0, vec![]),
            (b"", 0, vec![]),
        ] {
            let got = acc_case_variant_oids(&env, &db, acc, from, &mut known).unwrap();
            assert_eq!(got, want, "{}", String::from_utf8_lossy(acc));
        }
    }
}

/// `case_variant_keys` asks once for each letter where nothing differs in case, follows a
/// letter of the other case only where a prefix is found, and gives the variants in byte order.
#[test]
fn case_variant_keys_follow_found_prefixes_only() {
    let keys: Vec<&[u8]> = vec![b"pir||ABc1", b"pir||AbC1", b"pir||Abc1", b"pir||abc1"];
    let mut asked = Vec::new();
    let found = case_variant_keys(b"pir||Abc1", 4, |prefix: &[u8]| {
        asked.push(prefix.to_vec());
        Ok::<_, ()>(keys.iter().any(|key| key.starts_with(prefix)))
    })
    .unwrap();
    assert_eq!(
        found,
        vec![
            b"pir||ABc1".to_vec(),
            b"pir||AbC1".to_vec(),
            b"pir||abc1".to_vec()
        ]
    );
    // Never the own letters' prefixes, nor beyond a prefix not found.
    assert!(!asked.iter().any(|prefix| prefix.as_slice() == b"pir||Ab"));
    assert!(!asked
        .iter()
        .any(|prefix| prefix.starts_with(b"pir||aB") && prefix.len() > b"pir||aB".len()));
    let mut asked = 0;
    let none = case_variant_keys(b"gb_LOCUS9", 0, |_: &[u8]| {
        asked += 1;
        Ok::<_, ()>(false)
    })
    .unwrap();
    assert!(none.is_empty());
    assert_eq!(asked, 7, "one question for each letter");
}

/// Every acc2oid key of the dump of `path` (`name` names the `expected` file).
fn check_accession_keys(path: &str, name: &str) {
    {
        let f = File::open(path).unwrap();
        let env = Env::open(&f).unwrap();
        let text = std::fs::read_to_string(exp_path(name)).unwrap();
        let mut in_acc = false;
        let mut n = 0;
        for l in text.lines() {
            if l.starts_with("DB ") {
                in_acc = l.starts_with("DB acc2oid ");
            } else if in_acc && l.starts_with("K ") {
                let p: Vec<&str> = l.split(' ').collect();
                let key: Vec<u8> = (0..p[1].len() / 2)
                    .map(|i| u8::from_str_radix(&p[1][2 * i..2 * i + 2], 16).unwrap())
                    .collect();
                let want: Vec<u32> = p[2]
                    .split(',')
                    .map(|v| {
                        u32::from_str_radix(
                            &format!("{}{}{}{}", &v[6..8], &v[4..6], &v[2..4], &v[0..2]),
                            16,
                        )
                        .unwrap()
                    })
                    .collect();
                assert_eq!(
                    acc_first_oid(&env, &key).unwrap(),
                    Some(want[0]),
                    "{name} {}",
                    String::from_utf8_lossy(&key)
                );
                assert_eq!(
                    acc_oids(&env, &key).unwrap(),
                    want,
                    "{name} {}",
                    String::from_utf8_lossy(&key)
                );
                n += 1;
            }
        }
        assert!(n > 0, "{name}");
    }
}

/// Every acc2oid key of the dump (all fixtures).
#[test]
fn every_accession_key_matches_dump_py() {
    for name in ["p5.pdb", "n5.ndb", "mv.pdb", "dup.pdb", "thr501.pdb"] {
        check_accession_keys(&format!("{FX}/{name}"), name);
    }
}

#[test]
fn n5_literals_from_the_notes() {
    let f = files(&format!("{FX}/n5"), 'n');
    let ntf = Env::open(&f.ntf).unwrap();
    let ndb = Env::open(&f.ndb).unwrap();
    let num = total_oids(&volumes(&ndb).unwrap()) as u32;
    assert_eq!(num, 5);
    // n5.nto: 562 -> {2,4}; 7227 -> {0,3,4}; 10090 -> {1}
    let r = |t: &[i32]| get_oids_for_tax_ids(&ntf, &f.nto, num, &set(t)).unwrap();
    assert_eq!(r(&[562]).oids, vec![2, 4]);
    assert_eq!(r(&[7227]).oids, vec![0, 3, 4]);
    assert_eq!(r(&[10090]).oids, vec![1]);
    let both = r(&[7227, 562, 9606]);
    assert_eq!(both.oids, vec![0, 2, 3, 4]);
    assert_eq!(both.tax_ids_found, vec![562, 7227]);
    assert_eq!(r(&[9606]), TaxOids::default());
    // n5.not: oid 4 has two taxids (two deflines)
    assert_eq!(taxids_of_oid(&f.not, 0).unwrap(), vec![7227]);
    assert_eq!(taxids_of_oid(&f.not, 4).unwrap(), vec![562, 7227]);
    assert!(taxids_of_oid(&f.not, 5).is_err());
    // n5.nos: OID 4 = DD000001.1, Y00002.1
    let nos = f.nos.as_ref().unwrap();
    assert_eq!(
        seqids_of_oid(nos, 4).unwrap(),
        vec![b"DD000001.1".to_vec(), b"Y00002.1".to_vec()]
    );
    assert_eq!(seqids_of_oid(nos, 0).unwrap(), vec![b"AB000001.1".to_vec()]);
    assert!(seqids_of_oid(nos, 5).is_err());
    // negative: {562}: OID 2 only (OID 4 also has 7227); {562, 7227}: all four
    let neg = |t: &[i32]| negative_tax_ids_to_oids(&ntf, &f.nto, &f.not, num, &set(t)).unwrap();
    assert_eq!(neg(&[562]).oids, vec![2]);
    assert_eq!(neg(&[562, 7227]).oids, vec![0, 2, 3, 4]);
    // size rule: 2 requested ids, found {562}; OID 4 has {562, 7227} (2 <= 2) but 7227 not found
    let n = neg(&[562, 1]);
    assert_eq!((n.oids, n.tax_ids_found), (vec![2], vec![562]));
    assert_eq!(neg(&[1]), TaxOids::default());
}

/// Every OID list, taxid list and seqid list of every fixture against the dump.py values.
fn check_lookup_files(prefix: &str, n: char, name: &str) {
    {
        let f = files(prefix, n);
        let ex = expected(name);
        let ndb = Env::open(&f.ndb).unwrap();
        let ntf = Env::open(&f.ntf).unwrap();
        let num = total_oids(&volumes(&ndb).unwrap()) as u32;
        // GetDBTaxIds
        let mut keys = get_db_tax_ids(&ntf).unwrap();
        let want_keys: Vec<i32> = ex.t.iter().map(|t| t.0).collect();
        assert_eq!(keys, want_keys, "{name} key order");
        // one taxid at a time
        for (tax, oids) in &ex.t {
            let r = get_oids_for_tax_ids(&ntf, &f.nto, num, &set(&[*tax])).unwrap();
            assert_eq!(&r.oids, oids, "{name} taxid {tax}");
            assert_eq!(r.tax_ids_found, vec![*tax]);
        }
        // all taxids at once: the union (a big list takes the bit-vector branch)
        keys.sort_unstable();
        let all = get_oids_for_tax_ids(&ntf, &f.nto, num, &keys.iter().copied().collect()).unwrap();
        let mut union: Vec<u32> = ex.t.iter().flat_map(|t| t.1.iter().copied()).collect();
        union.sort_unstable();
        union.dedup();
        assert_eq!(all.oids, union, "{name} union");
        assert_eq!(all.tax_ids_found, keys);
        // per OID
        let oidtax = OidTaxIds::open(&f.not).unwrap();
        assert_eq!(oidtax.num_oids(), ex.o.len() as u64);
        for (oid, want) in ex.o.iter().enumerate() {
            assert_eq!(
                &oidtax.taxids(oid as u32).unwrap(),
                want,
                "{name} taxids of {oid}"
            );
        }
        if let Some(nos) = &f.nos {
            let seq = OidSeqIds::open(nos).unwrap();
            assert_eq!(seq.num_oids(), ex.s.len() as u64);
            for (oid, want) in ex.s.iter().enumerate() {
                let got: Vec<String> = seq
                    .seqids(oid as u32)
                    .unwrap()
                    .into_iter()
                    .map(|b| String::from_utf8(b).unwrap())
                    .collect();
                assert_eq!(&got, want, "{name} seqids of {oid}");
            }
        }
    }
}

/// Every OID list, taxid list and seqid list of every fixture against the dump.py values.
#[test]
fn lookup_files_match_dump_py() {
    for (prefix, n, name) in all_dbs() {
        check_lookup_files(&prefix, n, name);
    }
}

fn parse_oracle(name: &str) -> Vec<(char, Vec<i32>, Option<Vec<u32>>)> {
    let text = std::fs::read_to_string(exp_path(name)).unwrap();
    text.lines()
        .map(|l| {
            let p: Vec<&str> = l.splitn(3, ' ').collect();
            let rest = p.get(2).copied().unwrap_or("");
            let oids = if rest.starts_with("Error") {
                None
            } else {
                Some(
                    rest.split_whitespace()
                        .map(|x| x.parse().unwrap())
                        .collect(),
                )
            };
            (p[0].chars().next().unwrap(), csv(p[1]), oids)
        })
        .collect()
}

/// NCBI itself (blastdbcmd -taxids for the positive rule, blastp -negative_taxids with one
/// query per OID for the negative rule), on a DB with OIDs that have 1, 2, 3 taxids, taxid 0
/// and a taxid listed twice. `gen/oracle.sh` wrote the file.
#[test]
fn positive_and_negative_rules_match_ncbi_tx_fixture() {
    let f = files(&format!("{FX}/tx"), 'p');
    let ntf = Env::open(&f.ntf).unwrap();
    let num = total_oids(&volumes(&Env::open(&f.ndb).unwrap()).unwrap()) as u32;
    assert_eq!(num, 14);
    let mut checked = (0, 0);
    for (kind, taxids, oids) in parse_oracle("tx_oracle") {
        let want_set = set(&taxids);
        match kind {
            'P' => {
                let r = get_oids_for_tax_ids(&ntf, &f.nto, num, &want_set).unwrap();
                // NCBI's blastdbcmd errors "Taxonomy ID(s) not found" when nothing matched
                assert_eq!(Some(r.oids), oids.or(Some(vec![])), "P {taxids:?}");
                checked.0 += 1;
            }
            'N' => {
                // the oracle lists the OIDs still searched; the port returns the excluded ones
                let excluded = negative_tax_ids_to_oids(&ntf, &f.nto, &f.not, num, &want_set)
                    .unwrap()
                    .oids;
                let included: Vec<u32> = (0..num).filter(|o| !excluded.contains(o)).collect();
                assert_eq!(included, oids.unwrap(), "N {taxids:?}");
                checked.1 += 1;
            }
            _ => panic!(),
        }
    }
    assert_eq!(checked, (16, 16));
}

/// Needs the NCBI-distributed `18S_fungal_sequences` database (`LOSAT_TEST_18S_DB`, path without
/// extension). The expected OIDs are what `blastdbcmd -taxids` printed.
#[test]
#[ignore = "needs LOSAT_TEST_18S_DB"]
fn positive_rule_matches_ncbi_on_distributed_18s() {
    let f = files(&s18(), 'n');
    let ntf = Env::open(&f.ntf).unwrap();
    let num = total_oids(&volumes(&Env::open(&f.ndb).unwrap()).unwrap()) as u32;
    let rows = parse_oracle("18S_oracle");
    assert!(rows.len() >= 10);
    for (_, taxids, oids) in rows {
        let r = get_oids_for_tax_ids(&ntf, &f.nto, num, &set(&taxids)).unwrap();
        assert_eq!(r.oids, oids.unwrap_or_default(), "{taxids:?}");
    }
}

/// The negative rule re-derived from the dump.py tables on random sets, for every fixture
/// (independent code: candidates from the T table, per-OID lists from the O table).
fn check_negative_rule(prefix: &str, n: char, name: &str, next: &mut impl FnMut() -> usize) {
    {
        let f = files(prefix, n);
        let ex = expected(name);
        let ntf = Env::open(&f.ntf).unwrap();
        let num = ex.o.len() as u32;
        for round in 0..60 {
            let size = 1 + next() % 5;
            let mut req = BTreeSet::new();
            for _ in 0..size {
                // mostly present taxids, sometimes an absent one
                if next().is_multiple_of(5) {
                    req.insert(1_000_000 + (next() % 50) as i32);
                } else {
                    req.insert(ex.t[next() % ex.t.len()].0);
                }
            }
            let found: BTreeSet<i32> =
                ex.t.iter()
                    .map(|t| t.0)
                    .filter(|t| req.contains(t))
                    .collect();
            let mut cand: Vec<u32> =
                ex.t.iter()
                    .filter(|t| found.contains(&t.0))
                    .flat_map(|t| t.1.iter().copied())
                    .collect();
            cand.sort_unstable();
            cand.dedup();
            let want: Vec<u32> = cand
                .into_iter()
                .filter(|&o| {
                    let l = &ex.o[o as usize];
                    l.len() <= req.len() && l.iter().all(|t| found.contains(t))
                })
                .collect();
            let got = negative_tax_ids_to_oids(&ntf, &f.nto, &f.not, num, &req).unwrap();
            assert_eq!(got.oids, want, "{name} round {round} {req:?}");
            assert_eq!(got.tax_ids_found, found.iter().copied().collect::<Vec<_>>());
        }
    }
}

#[test]
fn negative_rule_on_random_sets() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    for (prefix, n, name) in all_dbs() {
        check_negative_rule(&prefix, n, name, &mut next);
    }
}

#[test]
fn corrupt_lookup_files_are_errors() {
    let f = files(&format!("{FX}/n5"), 'n');
    let ntf = Env::open(&f.ntf).unwrap();
    // offsets past the end of the .nto
    let short: Vec<u8> = std::fs::read(format!("{FX}/n5.nto")).unwrap()[..20].to_vec();
    let e = get_oids_for_tax_ids(&ntf, &short, 5, &set(&[10090])).unwrap_err();
    assert!(
        matches!(&e, TaxError::Data(m) if m.contains("beyond the end")),
        "{e:?}"
    );
    // a list that starts inside the file but runs past its end
    let e = get_oids_for_tax_ids(&ntf, &short, 5, &set(&[7227])).unwrap_err();
    assert!(
        matches!(&e, TaxError::Data(m) if m.contains("runs past the end")),
        "{e:?}"
    );
    assert!(get_oids_for_tax_ids(&ntf, &Vec::<u8>::new(), 5, &set(&[562])).is_err());
    // an OID at or above the OID count (undefined behaviour in NCBI)
    let e = get_oids_for_tax_ids(&ntf, &f.nto, 4, &set(&[562])).unwrap_err();
    assert!(
        matches!(&e, TaxError::Data(m) if m.contains("OID 4 is not below the 4")),
        "{e:?}"
    );
    // a huge count
    let mut bad = std::fs::read(format!("{FX}/n5.nto")).unwrap();
    bad[0..4].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(get_oids_for_tax_ids(&ntf, &bad, 5, &set(&[562])).is_err());
    // .not: empty file, truncated, bad header, ends going backwards
    assert_eq!(
        OidTaxIds::open(&Vec::<u8>::new())
            .err()
            .unwrap()
            .to_string(),
        "Failed to open oid-to-taxids lookup file"
    );
    assert_eq!(
        OidSeqIds::open(&Vec::<u8>::new())
            .err()
            .unwrap()
            .to_string(),
        "Failed to open oid-to-seqid lookup file"
    );
    let not = std::fs::read(format!("{FX}/n5.not")).unwrap();
    assert!(OidTaxIds::open(&not[..12]).is_err());
    let t = OidTaxIds::open(&not[..not.len() - 4]).unwrap();
    assert!(t.taxids(4).is_err());
    assert!(t.taxids(3).is_ok());
    let mut b = not.clone();
    b[8..16].copy_from_slice(&5u64.to_le_bytes());
    b[16..24].copy_from_slice(&1u64.to_le_bytes());
    assert!(OidTaxIds::open(&b)
        .unwrap()
        .taxids(1)
        .unwrap_err()
        .to_string()
        .contains("go backwards"));
    let mut b = not.clone();
    b[0..8].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(OidTaxIds::open(&b).is_err());
    // the negative rule reports a bad .not as Data
    let e =
        negative_tax_ids_to_oids(&ntf, &f.nto, &not[..40].to_vec(), 5, &set(&[562])).unwrap_err();
    assert!(matches!(e, TaxError::Data(_)), "{e:?}");
}

#[test]
fn seqid_list_forms() {
    // 3 OIDs: [two ids], [empty], [one id of 300 bytes written as 0xFF + u32 length]
    let long = vec![b'x'; 300];
    let mut data: Vec<u8> = Vec::new();
    data.extend([3, b'a', b'b', b'c', 2, b'd', b'e']);
    let end0 = data.len() as u64;
    let end1 = end0;
    data.push(0xFF);
    data.extend(300u32.to_le_bytes());
    data.extend(&long);
    let end2 = data.len() as u64;
    let mut file: Vec<u8> = Vec::new();
    file.extend(3u64.to_le_bytes());
    for e in [end0, end1, end2] {
        file.extend(e.to_le_bytes());
    }
    file.extend(&data);
    let s = OidSeqIds::open(&file).unwrap();
    assert_eq!(s.seqids(0).unwrap(), vec![b"abc".to_vec(), b"de".to_vec()]);
    assert_eq!(s.seqids(1).unwrap(), Vec::<Vec<u8>>::new());
    assert_eq!(s.seqids(2).unwrap(), vec![long]);
    assert!(s.seqids(3).is_err());
    // an id that runs past its OID's end
    let mut bad = file.clone();
    bad[8..16].copy_from_slice(&5u64.to_le_bytes());
    assert!(OidSeqIds::open(&bad)
        .unwrap()
        .seqids(0)
        .unwrap_err()
        .to_string()
        .contains("runs past"));
    // truncated file
    assert!(OidSeqIds::open(&file[..file.len() - 10])
        .unwrap()
        .seqids(2)
        .is_err());
    // a 0xFF length cut off
    let mut cut = file.clone();
    cut[24..32].copy_from_slice(&(end0 + 3).to_le_bytes());
    cut[16..24].copy_from_slice(&end0.to_le_bytes());
    assert!(OidSeqIds::open(&cut)
        .unwrap()
        .seqids(2)
        .unwrap_err()
        .to_string()
        .contains("truncated"));
}

/// Needs the NCBI-distributed `18S_fungal_sequences` database (`LOSAT_TEST_18S_DB`, path without
/// extension).
#[test]
#[ignore = "needs LOSAT_TEST_18S_DB"]
fn nos_of_distributed_18s_has_uncovered_tail() {
    // The NCBI-distributed 18S .nos has 12 bytes after the last end: the last OID reads empty.
    let f = File::open(format!("{}.nos", s18())).unwrap();
    let s = OidSeqIds::open(&f).unwrap();
    let n = s.num_oids() as u32;
    assert_eq!(s.seqids(n - 1).unwrap(), Vec::<Vec<u8>>::new());
    assert_eq!(s.seqids(0).unwrap(), vec![b"NG_065155.1".to_vec()]);
}

/// Volumes, every lookup list, every accession key and the negative rule of the NCBI-distributed
/// `18S_fungal_sequences` database against the `dump.py` outputs. Needs `LOSAT_TEST_18S_DB` (path
/// without extension) and `LOSAT_TEST_18S_EXPECTED` (directory with `18S.tax.txt` and
/// `18S_fungal_sequences.ndb.txt`).
#[test]
#[ignore = "needs LOSAT_TEST_18S_DB and LOSAT_TEST_18S_EXPECTED"]
fn ncbi_distributed_18s_lookup_files_match_dump_py() {
    let prefix = s18();
    check_volumes(&prefix, 'n', "18S");
    check_lookup_files(&prefix, 'n', "18S");
    check_accession_keys(&format!("{prefix}.ndb"), "18S_fungal_sequences.ndb");
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    check_negative_rule(&prefix, 'n', "18S", &mut next);
}

#[test]
fn many_oids_in_one_list_cross_the_read_blocks() {
    // 20000 OIDs in one list (more than one 16384 block), unsorted input OIDs with a duplicate
    let n = 20000u32;
    let mut nto: Vec<u8> = Vec::new();
    nto.extend((n + 1).to_le_bytes());
    for i in (0..n).rev() {
        nto.extend(i.to_le_bytes());
    }
    nto.extend(7u32.to_le_bytes());
    // a minimal .ntf is needed: borrow p5.ptf's layout is not possible, so build the lookup by
    // hand through the dumped fixture dup.ptf (taxid 0 -> offset 0)
    let f = File::open(format!("{FX}/dup.ptf")).unwrap();
    let ntf = Env::open(&f).unwrap();
    let r = get_oids_for_tax_ids(&ntf, &nto, n, &set(&[0])).unwrap();
    assert_eq!(r.oids, (0..n).collect::<Vec<_>>());
    assert_eq!(r.tax_ids_found, vec![0]);
}

/// Read back a version 5 database written by LOSAT's own writer (`writedb::lmdb`): five proteins
/// in two volumes, two taxids on some, none on one (the writer files it under taxid 0, as
/// makeblastdb does for a defline without a map entry).
#[test]
fn reads_back_a_database_written_by_the_losat_writer() {
    use crate::writedb::lmdb::{lmdb_file_name, LmdbEnvSettings, WriteDbLmdb, WriteDbTaxId};
    use crate::writedb::platform::{LmdbPlatform, Platform};

    let linux = Platform::LINUX_X86_64;
    let platform = LmdbPlatform::LINUX_X86_64;
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-tax-files-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // (accession, tax ids as written, tax ids as read back)
    let oids: [(&str, Vec<i32>); 5] = [
        ("NP_000001", vec![9606, 10090]),
        ("NP_000002", vec![9606]),
        ("NP_000003", vec![]),
        ("NP_000004", vec![562, 10090, 65536]),
        ("NP_000005", vec![562]),
    ];
    let pdb = dir.join("db.pdb");
    let ptf = lmdb_file_name(&pdb, "tf");
    let settings = LmdbEnvSettings::default();
    let mut lmdb =
        WriteDbLmdb::new(&pdb, linux.default_lmdb_map_size, &platform, &settings).unwrap();
    let mut tax =
        WriteDbTaxId::new(&ptf, linux.default_taxid_map_size, &platform, &settings).unwrap();
    for (oid, (acc, taxids)) in oids.iter().enumerate() {
        lmdb.insert_entry(acc.as_bytes(), oid as u32, false)
            .unwrap();
        lmdb.insert_entry(format!("{acc}.1").as_bytes(), oid as u32, true)
            .unwrap();
        tax.insert_entries(taxids, oid as u32);
    }
    lmdb.insert_volumes_info(&[b"db.00".as_slice(), b"db.01"], &[3, 2])
        .unwrap();
    lmdb.close().unwrap();
    tax.close().unwrap();

    let open = |ext: &str| File::open(dir.join(format!("db.{ext}"))).unwrap();
    let (ndb, ntf, nto, not, nos) = (
        open("pdb"),
        open("ptf"),
        open("pto"),
        open("pot"),
        open("pos"),
    );
    let db = Env::open(&ndb).unwrap();
    let tf = Env::open(&ntf).unwrap();

    let vols = volumes(&db).unwrap();
    assert_eq!(vols, vec![(b"db.00".to_vec(), 3), (b"db.01".to_vec(), 2)]);
    let num = total_oids(&vols) as u32;
    assert_eq!(num, 5);

    // what was written, as a read-back expects it: no tax ids = taxid 0
    let written: Vec<Vec<i32>> = oids
        .iter()
        .map(|o| if o.1.is_empty() { vec![0] } else { o.1.clone() })
        .collect();
    let mut all: Vec<i32> = written.iter().flatten().copied().collect();
    all.sort_unstable();
    all.dedup();
    let mut keys = get_db_tax_ids(&tf).unwrap();
    keys.sort_unstable();
    assert_eq!(keys, all);
    for &t in &all {
        let want: Vec<u32> = (0..num)
            .filter(|&o| written[o as usize].contains(&t))
            .collect();
        let got = get_oids_for_tax_ids(&tf, &nto, num, &set(&[t])).unwrap();
        assert_eq!((got.oids, got.tax_ids_found), (want, vec![t]), "taxid {t}");
    }
    let both = get_oids_for_tax_ids(&tf, &nto, num, &set(&[10090, 562, 7])).unwrap();
    assert_eq!(
        (both.oids, both.tax_ids_found),
        (vec![0, 3, 4], vec![562, 10090])
    );
    for (oid, want) in written.iter().enumerate() {
        assert_eq!(&taxids_of_oid(&not, oid as u32).unwrap(), want, "OID {oid}");
    }
    assert!(taxids_of_oid(&not, num).is_err());
    // negative rule: excluding 9606 drops OID 1 only (OID 0 also has 10090)
    assert_eq!(
        negative_tax_ids_to_oids(&tf, &nto, &not, num, &set(&[9606]))
            .unwrap()
            .oids,
        vec![1]
    );
    // accessions and seqids
    assert_eq!(acc_first_oid(&db, b"NP_000004.1").unwrap(), Some(3));
    assert_eq!(acc_oids(&db, b"NP_000002").unwrap(), vec![1]);
    for (oid, (acc, _)) in oids.iter().enumerate() {
        assert_eq!(
            seqids_of_oid(&nos, oid as u32).unwrap(),
            vec![format!("{acc}.1").into_bytes()],
            "OID {oid}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}
