//! `taxdb_may_refuse`: the answer of the lookup of every sequence's taxid, with the
//! databases of `tests/fixtures/outfmt16/db` (`sprot_tax`: two sequences, taxid 9606;
//! `sprot_mixtax`: those two with taxid 10090 and a third with taxid 9606) and the taxonomy
//! databases of `tests/seqdb/taxdb`.
use super::*;
use crate::seqdb::alias::SearchPath;
use crate::seqdb::SeqType;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
}

/// A database and the empty current directory it was opened from (its search path keeps it for
/// the files looked for later), which is removed with it.
struct Opened {
    db: SeqDb,
    cwd: PathBuf,
}

impl std::ops::Deref for Opened {
    type Target = SeqDb;

    fn deref(&self) -> &SeqDb {
        &self.db
    }
}

impl Drop for Opened {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.cwd);
    }
}

/// The database `name` of `tests/fixtures/outfmt16/db` opened with `BLASTDB=taxdb`
/// (`None`: not set), from an empty current directory.
fn database(name: &str, taxdb: Option<&Path>) -> Opened {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat_prepass_cwd_{}_{}",
        crate::utils::test_dir::process_id(),
        NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
    let path = root().join("fixtures/outfmt16/db").join(name);
    let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Protein, &search).unwrap();
    Opened { db, cwd }
}

fn variant(name: &str) -> PathBuf {
    root().join("seqdb/taxdb").join(name)
}

/// An empty current directory of the test `tag`, which the test removes at its end.
fn empty_cwd(tag: &str) -> PathBuf {
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat_{tag}_cwd_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    cwd
}

/// The subject of each sequence of a database with parsed Seq-ids, by ordinal id.
fn parsed_subjects(db: &SeqDb) -> Vec<crate::report::db_subject::ParsedSubject> {
    (0..db.num_seqs())
        .map(|ordinal| {
            let set = db.defline_set(ordinal).unwrap();
            crate::report::db_subject::ParsedSubject::new(set).unwrap()
        })
        .collect()
}

/// The database `name` of `tests/db_search` opened with `BLASTDB=taxdb` (`None`: not set),
/// from the empty current directory `cwd`.
fn db_search(name: &str, taxdb: Option<&Path>, cwd: &Path) -> SeqDb {
    let search = SearchPath::new(cwd, taxdb.map(Path::as_os_str), None).unwrap();
    let path = root().join("db_search").join(name);
    SeqDb::open_in(path.to_str().unwrap(), SeqType::Protein, &search).unwrap()
}

#[test]
fn a_reported_sequence_whose_seq_align_id_finds_another_record_first() {
    // tests/db_search/ids/make_ids.py (line DBID N1): the second database of an alias holds
    // the record of the first, which NCBI shows for it (`pdup_ab`), or the record with another
    // title, which it shows in its place (`pdup_ac`).
    let cwd = empty_cwd("found_instead");
    for (name, want) in [
        ("ids/pdup_ab", None),
        (
            "ids/pdup_ac",
            Some("lookup gives database sequence 0 first"),
        ),
        (
            "ids/pfound_v5",
            Some("lookup gives database sequence 0 first"),
        ),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db);
        let ReportedHeader::Parsed(subject) = read_header(&headers, 1, 0, "BLASTP").unwrap() else {
            panic!("{name}: a sequence with parsed Seq-ids");
        };
        match (
            subject.fetched(false, true, crate::report::db_subject::HitLoad::Loaded),
            want,
        ) {
            (Ok(()), None) => {}
            (Err(error), Some(text)) => assert!(error.contains(text), "{name}: {error}"),
            (got, _) => panic!("{name}: {got:?}"),
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1447-1455 (quoted
// at `SeqidLookup::seqid_to_oids`): without `multi`, the lookup returns at the first OID.
/// Line DBID v4late (review 3 of line DBID loadord, H3): in a list of version 4 volumes whose
/// later volume's string ISAM index does not open (mode 000), NCBI's lookup of an id that an
/// earlier volume holds returns there (NCBI BLAST+ 2.17.0: exit status 0, probe
/// `blastdb-ids-20261021/R/v4late`), where LOSAT's list of every holder fails: the other
/// holders are unknown, with what the lookup gives (the subject, the same record, another one).
/// An id that only the later volume can hold, or any id where that volume comes first, opens
/// it: NCBI's stop. What the prefetch reads (review 1 of the item, H1 and H2): what the lookup
/// gives, whose definition line `GetBioseq` may reject (NCBI then stops), or nothing where a
/// holder that LOSAT cannot list was read before: with a taxonomy database, the same texts of
/// both, or LOSAT's refusal.
#[cfg(unix)]
#[test]
fn a_later_volume_whose_isam_does_not_open_leaves_the_other_holders_unknown() {
    use crate::report::db_subject::FoundInstead;
    use std::os::unix::fs::PermissionsExt;
    // Made empty: a file of mode 000 left by a run before cannot be copied over.
    let cwd = empty_cwd("v4late");
    std::fs::remove_dir_all(&cwd).unwrap();
    std::fs::create_dir_all(&cwd).unwrap();
    let source = root().join("db_search/multidb");
    // `l5a`: a copy of l5e, whose name sorts before l5c's (the volumes of a list are in the order
    // of their names, as NCBI's).
    for (stem, from) in [
        ("l5c", "l5c"),
        ("l5e", "l5e"),
        ("l5a", "l5e"),
        ("l5b", "l5b"),
        ("l6o", "l6o"),
        ("l6u", "l6u"),
    ] {
        for ext in ["phr", "pin", "pog", "psd", "psi", "psq"] {
            let to = cwd.join(format!("{stem}.{ext}"));
            std::fs::copy(source.join(format!("{from}.{ext}")), &to).unwrap();
            if ext == "psi" && !matches!(stem, "l5b" | "l5c") {
                std::fs::set_permissions(&to, std::fs::Permissions::from_mode(0o000)).unwrap();
            }
        }
    }
    if std::fs::File::open(cwd.join("l5e.psi")).is_ok() {
        // A user that reads every file (root): nothing fails to open.
        std::fs::remove_dir_all(&cwd).unwrap();
        return;
    }
    let unknown = "for which LOSAT cannot list the other sequences that hold it, as a file of a volume after the first one that holds it, which NCBI BLAST+'s lookup does not read, fails (Error: Could not open input file (";
    let open = |list: &str, taxdb: Option<&str>| {
        let names = list
            .split(' ')
            .map(|stem| cwd.join(stem).to_str().unwrap().to_owned())
            .collect::<Vec<_>>()
            .join(" ");
        let taxdb = taxdb.map(|variant| root().join("seqdb/taxdb").join(variant));
        let search = SearchPath::new(&cwd, taxdb.as_deref().map(Path::as_os_str), None).unwrap();
        SeqDb::open_in(&names, SeqType::Protein, &search).unwrap()
    };
    let subject_of = |db: &SeqDb, ordinal: u32| {
        let headers = HeaderReader::new(db);
        match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
            ReportedHeader::Parsed(subject) => subject,
            _ => panic!("{ordinal}: a sequence with parsed Seq-ids"),
        }
    };
    for (list, ordinal, instead, stops) in [
        // The subject's own record first, then the same record, then another one.
        ("l5c l5e", 0, None, false),
        ("l5c l5e", 1, None, false),
        ("l5c l6o", 1, Some(0), false),
        // `ref|NP_000004.1|`, which l5c does not hold, and l5c's id behind l5a.
        ("l5c l6u", 1, None, true),
        ("l5c l5a", 1, None, true),
    ] {
        let db = open(list, None);
        let subject = subject_of(&db, ordinal);
        if stops {
            assert!(
                subject.prefetch_stops().is_some_and(
                    |message| message.starts_with(b"Error: Could not open input file (")
                ),
                "{list} {ordinal}: {:?}",
                subject.found()
            );
            continue;
        }
        assert!(subject.prefetch_stops().is_none(), "{list} {ordinal}");
        match subject.found() {
            Some(FoundInstead::HoldersUnknown { instead: got, .. }) => assert_eq!(
                got.map(|other| other.ordinal()),
                instead,
                "{list} {ordinal}"
            ),
            got => panic!("{list} {ordinal}: {got:?}"),
        }
        let refusal = subject
            .fetched(false, false, crate::report::db_subject::HitLoad::Loaded)
            .unwrap_err();
        assert!(refusal.contains(unknown), "{list} {ordinal}: {refusal}");
        assert_eq!(
            refusal.contains("lookup gives database sequence 0 first"),
            instead.is_some(),
            "{list} {ordinal}: {refusal}"
        );
        assert!(subject.pssm_read(false).is_err(), "{list} {ordinal}");
    }
    // l5b, which the lookup of l5e's `ref|NP_000002.1|` gives first, holds it after
    // `ref|NP_000999.1|` in one definition line: NCBI's prefetch reads it and stops.
    let db = open("l5b l5e", None);
    let subject = subject_of(&db, 1);
    assert!(matches!(
        subject.found(),
        Some(FoundInstead::HoldersUnknown { instead: Some(other), .. }) if other.ordinal() == 0
    ));
    let thrown = prefetch_reads(&db, 1, &subject, false, "BLASTP")
        .unwrap()
        .unwrap_err();
    assert!(
        thrown.contains("oid headers do not contain target gi/seq_id"),
        "{thrown}"
    );
    // The subject's own blob (taxid 9606) or nothing: the same texts with a whole taxonomy
    // database, other texts with one of a wrong magic number (its texts of opening).
    for (taxdb, refused) in [("good", false), ("bad_magic", true)] {
        let db = open("l5c l5e", Some(taxdb));
        for ordinal in [0, 1] {
            match prefetch_reads(&db, ordinal, &subject_of(&db, ordinal), true, "BLASTP").unwrap() {
                Ok(blob) => {
                    assert!(!refused, "{taxdb} {ordinal}");
                    assert_eq!(blob, Some(ordinal), "{taxdb} {ordinal}");
                }
                Err(reason) => {
                    assert!(refused, "{taxdb} {ordinal}: {reason}");
                    assert!(reason.contains(unknown), "{taxdb} {ordinal}: {reason}");
                }
            }
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn the_same_record_with_other_taxids_is_shown_apart_by_the_reports_of_taxids_only() {
    // tests/db_search/ids/make_ids.py (line DBID N1b): `pdup_at` holds `gb|EEE00005.1|`
    // without a taxid then with taxid 9606, `pdup_tu` the other way, `pdup_no` twice with
    // taxid 99990 (the last entry of every variant of tests/seqdb/taxdb, whose lookup writes
    // the offset error of `offset_last_past_end`), `pdupe_at` an empty title without a taxid
    // then with 9606; `taxdb`, the variant on BLASTDB. `without`, `with`: the refusal without
    // the reports of taxids (`fetched(false)`) and with them.
    let cwd = empty_cwd("other_taxids");
    let first = |n: u32| format!("lookup gives database sequence {n} first");
    let taxids = |n: u32| {
        format!("database sequence {n} also holds with the same definition lines and letters but other taxids (")
    };
    let title = |n: u32| {
        format!("database sequence {n} also holds with the same definition lines and letters but other taxids, from which")
    };
    let texts = |n: u32| {
        format!(
            "database sequence {n} also holds with the same definition lines and letters: where"
        )
    };
    let both = |text: String| (Some(text.clone()), Some(text));
    for (name, taxdb, ordinal, (without, with)) in [
        ("ids/pdup_at", None, 0, (None, Some(taxids(1)))),
        ("ids/pdup_at", None, 1, (None, Some(taxids(0)))),
        ("ids/pdup_tu", None, 0, (None, Some(taxids(1)))),
        ("ids/pdup_tu", None, 1, (None, Some(taxids(0)))),
        ("ids/pdup_at", Some("good"), 1, (None, Some(taxids(0)))),
        // Texts of opening at the run's first lookup, which only the copy with a taxid makes.
        ("ids/pdup_at", Some("bad_magic"), 0, both(texts(1))),
        ("ids/pdup_at", Some("bad_magic"), 1, both(texts(0))),
        ("ids/pdup_tu", Some("count_larger"), 1, both(texts(0))),
        // The offset error of a taxid that neither copy holds.
        (
            "ids/pdup_at",
            Some("offset_last_past_end"),
            1,
            (None, Some(taxids(0))),
        ),
        // The same taxids: the texts of opening come at the same lookup, the offset error
        // once per volume.
        ("ids/pdup_no", None, 1, (None, None)),
        ("ids/pdup_no", Some("good"), 1, (None, None)),
        ("ids/pdup_no", Some("bad_magic"), 1, (None, None)),
        ("ids/pdup_no", Some("count_larger"), 0, (None, None)),
        (
            "ids/pdup_no",
            Some("offset_last_past_end"),
            1,
            both(texts(0)),
        ),
        (
            "ids/pdup_no",
            Some("offset_last_past_end"),
            0,
            both(texts(1)),
        ),
        // The same record without taxids looks nothing up.
        ("ids/pdup_ab", Some("offset_last_past_end"), 1, (None, None)),
        // An empty title made from the taxids where the taxonomy database is found, which the
        // reports of generated titles refuse (below).
        ("ids/pdupe_at", None, 0, (None, Some(taxids(1)))),
        ("ids/pdupe_at", Some("good"), 0, (None, Some(taxids(1)))),
        ("ids/pdupe_at", Some("good"), 1, (None, Some(taxids(0)))),
        // A taxonomy database that names no taxid: NCBI makes no title from them.
        ("ids/pdupe_at", Some("no_btd"), 0, (None, Some(taxids(1)))),
        (
            "ids/pdupe_at",
            Some("empty_fields"),
            0,
            (None, Some(taxids(1))),
        ),
        // One that writes a text at the run's first lookup, which only the copy with a taxid
        // makes.
        ("ids/pdupe_at", Some("bad_magic"), 0, both(texts(1))),
        // Another title is another record whatever the taxids.
        ("ids/pdup_ac", None, 1, both(first(0))),
    ] {
        let taxdb = taxdb.map(|variant| root().join("seqdb/taxdb").join(variant));
        let db = db_search(name, taxdb.as_deref(), &cwd);
        let headers = HeaderReader::new(&db);
        let ReportedHeader::Parsed(subject) = read_header(&headers, ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
        };
        for (taxonomy, want) in [(false, &without), (true, &with)] {
            match (
                subject.fetched(taxonomy, true, crate::report::db_subject::HitLoad::Loaded),
                want,
            ) {
                (Ok(()), None) => {}
                (Err(error), Some(text)) => {
                    assert!(error.contains(text.as_str()), "{name} {ordinal}: {error}")
                }
                (got, _) => panic!("{name} {ordinal} {taxdb:?} {taxonomy}: {got:?}"),
            }
        }
        // Line DBID N1b L1: outfmt 0 and the tabular word of a replaced id
        // (`ParsedSubject::generated_title_fetched`).
        let shown =
            name == "ids/pdupe_at" && taxdb.as_deref().is_some_and(|dir| dir.ends_with("good"));
        match subject.generated_title_fetched() {
            Ok(()) => assert!(!shown, "{name} {ordinal} {taxdb:?}"),
            Err(error) => assert!(
                shown && error.contains(&title(1 - ordinal)),
                "{name} {ordinal} {taxdb:?}: {error}"
            ),
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn a_seq_align_id_of_any_kind_is_left_to_the_lookup() {
    // tests/db_search/ids/make_look_fasta.py (line DBID N2b): `None`, NCBI's formatters find the
    // sequence itself; else the text of the refusal.
    let cwd = empty_cwd("any_kind");
    let other = |n: u32| format!("lookup gives database sequence {n} first");
    let none = "formatters do not find the sequence".to_string();
    let held = |n: u32| format!("which database sequence {n} also holds");
    for (name, ordinal, want) in [
        ("ids/plook_v5", 1, None),
        ("ids/plook_v5", 2, None),
        ("ids/plook_v5", 3, None),
        ("ids/plook_v5", 4, Some(none.clone())),
        ("ids/plook_v5", 5, None),
        ("ids/plook_v5", 6, Some(other(5))),
        ("ids/plook_v5", 8, Some(other(7))),
        ("ids/plook_v5", 9, None),
        ("ids/plook_v5", 10, None),
        ("ids/plook_v5", 11, Some(other(4))),
        ("ids/plook_v5", 12, Some(none.clone())),
        ("ids/plook_v5", 13, None),
        ("ids/plook_v4", 1, Some(none.clone())),
        ("ids/plook_v4", 2, Some(none.clone())),
        ("ids/plook_v4", 3, None),
        ("ids/plook_v4", 4, Some(none.clone())),
        ("ids/plook_v4", 5, Some(other(6))),
        ("ids/plook_v4", 9, None),
        ("ids/plook_v4", 11, Some(other(4))),
        ("ids/plook_v4", 12, Some(none.clone())),
        ("ids/plook_vz", 16, None),
        // Another sequence holds the Seq-align id (review 1 of N2b): NCBI's scope shows it
        // where both are reported.
        ("ids/pshare_v5", 0, Some(held(1))),
        ("ids/pshare_v5", 1, None),
        ("ids/pshare_v4", 0, Some(held(1))),
        ("ids/pshare_ab", 0, Some(held(1))),
        // Ids that differ in case alone (review 2 of N2b), and a later copy of a record whose
        // id a third sequence also holds.
        ("ids/pshare_v5", 4, Some(held(5))),
        ("ids/pshare_v4", 4, Some(held(5))),
        ("ids/pshare_cd", 0, Some(held(2))),
        ("ids/pshare_cd", 1, Some(held(3))),
        ("ids/pdup_abc", 0, Some(held(2))),
        ("ids/pdup_abc", 1, Some(held(2))),
        (
            "ids/pshare_v5",
            2,
            Some("whose lookup in the database fails".to_string()),
        ),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db);
        let ReportedHeader::Parsed(subject) = read_header(&headers, ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
        };
        assert!(!subject.unfetchable(), "{name} {ordinal}");
        match (
            subject.fetched(false, true, crate::report::db_subject::HitLoad::Loaded),
            want,
        ) {
            (Ok(()), None) => {}
            (Err(error), Some(text)) => {
                assert!(error.contains(&text), "{name} {ordinal}: {error}")
            }
            (got, _) => panic!("{name} {ordinal}: {got:?}"),
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn the_prefetch_reads_a_subject_by_its_seq_align_id_as_ncbis_loader_does() {
    // Line DBID P3 (probes of `blastdb-ids-20261021/R/intmaxpf`): where no report reads the
    // Bioseq, `prefetch_reads` without texts (no `taxdb.bti` on this search path).
    let cwd = empty_cwd("prefetch_reads");
    let stops = "and then stops (oid headers do not contain target gi/seq_id)";
    for (name, ordinal, want) in [
        // `pdb|1ABC|A` of pshare_c and `pdb|1abc|A` of pshare_d: the scope can read the
        // other, whose definition lines do not hold the id (NCBI exit 2).
        ("ids/pshare_cd", 0, Err(stops)),
        // `lcl|555`, which sequence 1 also holds in a second definition line: read as it or
        // as itself, nothing to refuse without texts.
        ("ids/pshare_v5", 0, Ok(None)),
        ("ids/pshare_v5", 1, Ok(Some(1))),
        // `gnl|BL_ORD_ID|0` stored with a text tag.
        (
            "ids/pshare_v5",
            2,
            Err("whose lookup in the database fails"),
        ),
        // An id that the lookup finds nowhere, and one that it finds first on sequence 5,
        // none of whose ids compares equal: nothing is read.
        ("ids/plook_v5", 4, Ok(None)),
        ("ids/plook_v5", 6, Ok(None)),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db);
        let ReportedHeader::Parsed(subject) = read_header(&headers, ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
        };
        match (
            prefetch_reads(&db, ordinal, &subject, false, "BLASTP").unwrap(),
            want,
        ) {
            (Ok(got), Ok(want)) => assert_eq!(got, want, "{name} {ordinal}"),
            (Err(reason), Err(text)) => {
                assert!(reason.contains(text), "{name} {ordinal}: {reason}")
            }
            (got, _) => panic!("{name} {ordinal}: {got:?}"),
        }
        // A lookup that fails is the batch's warning (`batch_prefetch`), not a refusal.
        let check = check_prefetched_parsed(&db, &subject, ordinal, "BLASTP");
        assert_eq!(
            check.is_err(),
            want.is_err() && !subject.prefetch_fails(),
            "{name} {ordinal}"
        );
        assert_eq!(
            subject.prefetch_fails(),
            want == Err("whose lookup in the database fails"),
            "{name} {ordinal}"
        );
    }
    // `lcl|777` before `gi|777`: the formatters do not find the Bioseq, which the prefetch
    // reads as itself; it fails where it asks for the length (a minus-strand alignment).
    let search = SearchPath::new(&cwd, None, None).unwrap();
    let path = root().join("db_search/ids/nrefuse_v5");
    let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
    let headers = HeaderReader::new(&db);
    let ReportedHeader::Parsed(subject) = read_header(&headers, 0, 0, "TBLASTX").unwrap() else {
        panic!("nrefuse_v5 0: a sequence with parsed Seq-ids");
    };
    assert!(subject.unfetchable());
    assert_eq!(
        prefetch_reads(&db, 0, &subject, true, "TBLASTX").unwrap(),
        Ok(Some(0))
    );
    let run = TaxdbRun::default();
    let warning = b"Warning: [tblastx] Error pre-fetching sequence data\n".to_vec();
    let prefetch = |flipped: bool| {
        let subject = PrefetchSubject {
            ordinal: 0,
            parsed: Some(&*subject),
            flipped,
        };
        batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [subject],
            "TBLASTX",
            &stop_queries(0),
        )
        .unwrap()
        .text()
    };
    assert_eq!(prefetch(false), b"");
    assert_eq!(prefetch(true), warning);
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn a_batch_whose_prefetch_fails_writes_ncbis_warning_after_its_texts() {
    // tests/db_search/ids/make_ids.py (line DBID pfwarn): `ntag_t5` 0 `gnl|BL_ORD_ID|0`
    // stored with a text tag (its lookup fails), 1 `gb|FFF00092.1|`; `nlook_t4` 1 `bbs|123`,
    // which a version 4 database does not look up (no Bioseq), 0 `gb|AF000070.1|`; every
    // sequence with taxid 99990, whose lookup in `offset_last_past_end` writes the offset
    // error.
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-batch-prefetch-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let taxdb = root().join("seqdb/taxdb/offset_last_past_end");
    let open = |name: &str, taxdb: Option<&Path>| {
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
        let path = root().join("db_search/ids").join(name);
        SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap()
    };
    let offset = "Error: [blastn] Error: Offset error at end of taxdb file.\n";
    let warning = "Warning: [blastn] Error pre-fetching sequence data\n";
    let text = |posts: crate::report::diag_log::Posts| String::from_utf8(posts.text()).unwrap();
    for taxdb in [None, Some(taxdb.as_path())] {
        let db = open("ntag_t5", taxdb);
        let headers = HeaderReader::new(&db);
        let parsed = |ordinal: u32| match read_header(&headers, ordinal, 0, "BLASTN").unwrap() {
            ReportedHeader::Parsed(subject) => subject,
            _ => panic!("ntag_t5 {ordinal}: a sequence with parsed Seq-ids"),
        };
        let (failed, other) = (parsed(0), parsed(1));
        assert!(failed.prefetch_fails() && !other.prefetch_fails());
        let run = TaxdbRun::default();
        let alone = batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [PrefetchSubject {
                ordinal: 0,
                parsed: Some(&*failed),
                flipped: false,
            }],
            "BLASTN",
            &stop_queries(0),
        );
        assert_eq!(text(alone.unwrap()), warning);
        // With the other subject, whose load writes the offset error before the failed
        // lookup or not at all: refused where the taxonomy database writes it.
        let both = batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [
                PrefetchSubject {
                    ordinal: 0,
                    parsed: Some(&*failed),
                    flipped: false,
                },
                PrefetchSubject {
                    ordinal: 1,
                    parsed: Some(&*other),
                    flipped: false,
                },
            ],
            "BLASTN",
            &stop_queries(0),
        );
        match taxdb {
            None => assert_eq!(text(both.unwrap()), warning),
            Some(_) => {
                let refused = both.unwrap_err().to_string();
                assert!(
                    refused.contains("in an order that LOSAT does not reproduce"),
                    "{refused}"
                );
            }
        }
        // A later batch of the run, whose lookup of the other subject wrote its text before:
        // nothing to place.
        if taxdb.is_some() {
            let run = TaxdbRun::default();
            let first = batch_prefetch(
                &db,
                &run,
                std::iter::empty(),
                [PrefetchSubject {
                    ordinal: 1,
                    parsed: Some(&*other),
                    flipped: false,
                }],
                "BLASTN",
                &stop_queries(0),
            );
            assert_eq!(text(first.unwrap()), offset);
            let later = batch_prefetch(
                &db,
                &run,
                std::iter::empty(),
                [
                    PrefetchSubject {
                        ordinal: 1,
                        parsed: Some(&*other),
                        flipped: false,
                    },
                    PrefetchSubject {
                        ordinal: 0,
                        parsed: Some(&*failed),
                        flipped: false,
                    },
                ],
                "BLASTN",
                &stop_queries(0),
            );
            assert_eq!(text(later.unwrap()), warning);
        }
    }
    // A minus-strand Seq-align of a subject without a Bioseq fails after every load: the
    // texts of the batch, then the warning.
    let db = open("nlook_t4", Some(&taxdb));
    let headers = HeaderReader::new(&db);
    let parsed = |ordinal: u32| match read_header(&headers, ordinal, 0, "BLASTN").unwrap() {
        ReportedHeader::Parsed(subject) => subject,
        _ => panic!("nlook_t4 {ordinal}: a sequence with parsed Seq-ids"),
    };
    let (control, bbs) = (parsed(0), parsed(1));
    let run = TaxdbRun::default();
    let prefetch = batch_prefetch(
        &db,
        &run,
        std::iter::empty(),
        [
            PrefetchSubject {
                ordinal: 1,
                parsed: Some(&*bbs),
                flipped: true,
            },
            PrefetchSubject {
                ordinal: 0,
                parsed: Some(&*control),
                flipped: false,
            },
        ],
        "BLASTN",
        &stop_queries(0),
    );
    assert_eq!(text(prefetch.unwrap()), format!("{offset}{warning}"));
    let run = TaxdbRun::default();
    let plus = batch_prefetch(
        &db,
        &run,
        std::iter::empty(),
        [PrefetchSubject {
            ordinal: 1,
            parsed: Some(&*bbs),
            flipped: false,
        }],
        "BLASTN",
        &stop_queries(0),
    );
    assert_eq!(text(plus.unwrap()), "");
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn the_loads_before_a_stopping_lookup_follow_ncbis_order_of_seq_id_handles() {
    // Line DBID loadord (tests/db_search/make_multidb.py --loadord-only): `ldp` holds
    // `ref|NP_000002.1|` (packed 2) and `ref|NP_000900.1|`; the unparsed `qp1` is
    // `gnl|BL_ORD_ID|0` beside it in `lda ldp`, `|1` in `ldb ldp` and `|2` in `ldp ldu`.
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-load-order-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let open = |names: &[&str]| {
        let search = SearchPath::new(&cwd, None, None).unwrap();
        let list: Vec<String> = names
            .iter()
            .map(|name| {
                root()
                    .join("db_search/multidb")
                    .join(name)
                    .to_str()
                    .unwrap()
                    .to_owned()
            })
            .collect();
        SeqDb::open_in(&list.join(" "), SeqType::Protein, &search).unwrap()
    };
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db);
            (0..db.num_seqs())
                .map(
                    |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                        ReportedHeader::Parsed(subject) => Some(subject),
                        _ => None,
                    },
                )
                .collect()
        };
    let accession =
        |parsed: &[Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
         number: &str| {
            parsed
                .iter()
                .position(|header| {
                    header.as_ref().is_some_and(|subject| {
                        subject.prefetch_stops().is_some()
                            && String::from_utf8_lossy(&subject.alignment_fasta()).contains(number)
                    })
                })
                .unwrap() as u32
        };
    for (names, want) in [
        (["lda", "ldp"], LoadBeforeStop::No),
        (["ldb", "ldp"], LoadBeforeStop::Yes),
    ] {
        let db = open(&names);
        let parsed = headers_of(&db);
        let stop = accession(&parsed, "NP_000002.1");
        let unparsed = parsed.iter().rposition(Option::is_none).unwrap() as u32;
        let subjects = [
            PrefetchSubject::unflipped(unparsed, None),
            PrefetchSubject::unflipped(stop, parsed[stop as usize].as_deref()),
        ];
        assert_eq!(
            load_before_stop(&db, &subjects).at_first_stop(),
            [want, LoadBeforeStop::No],
            "{names:?}"
        );
        let before = load_before_stop(&db, &subjects).at_first_stop();
        // No load before the batch: the run has reported nothing (its first batch, or batches
        // whose subjects have no parsed Seq-ids), and the scope holds none of its ids.
        assert_eq!(
            stop_lookup_unsure(
                &db,
                &TaxdbRun::default(),
                &subjects,
                &before,
                &stop_queries(0),
                stop,
                parsed[stop as usize].as_ref().unwrap()
            )
            .unwrap(),
            None
        );
    }
    // The same packed number 2 (the handles' addresses decide), one Seq-align id of two
    // sequences, more than 200 ids, and two stops (the smaller number decides).
    let db = open(&["ldp", "ldu"]);
    let parsed = headers_of(&db);
    let stop = accession(&parsed, "NP_000002.1");
    let other = accession(&parsed, "NP_000900.1");
    let tie = PrefetchSubject::unflipped(2, None);
    let stopping = PrefetchSubject::unflipped(stop, parsed[stop as usize].as_deref());
    assert!(matches!(
        load_before_stop(&db, &[tie, stopping]).at_first_stop()[0],
        LoadBeforeStop::Unknown(why) if why.contains("same packed number")
    ));
    let unstopped = (**parsed[stop as usize].as_ref().unwrap())
        .clone()
        .found_instead(None);
    let copy = PrefetchSubject::unflipped(other, Some(&unstopped));
    assert!(matches!(
        load_before_stop(&db, &[stopping, copy]).at_first_stop()[1],
        LoadBeforeStop::Unknown(why) if why.contains("the same Seq-align id")
    ));
    // `gnl|BL_ORD_ID|3` and on (packed 3 and more, after the stop): 200 distinct ids in one
    // map, 201 in pieces.
    let unparsed: Vec<PrefetchSubject<'_>> = (3..204)
        .map(|ordinal| PrefetchSubject::unflipped(ordinal, None))
        .collect();
    let mut many = unparsed.clone();
    many.push(stopping);
    assert!(load_before_stop(&db, &many).at_first_stop()[..201]
        .iter()
        .all(|before| matches!(
            before,
            LoadBeforeStop::Unknown(why) if why.contains("more than 200")
        )));
    let mut few = unparsed[..199].to_vec();
    few.push(stopping);
    assert!(load_before_stop(&db, &few)
        .at_first_stop()
        .iter()
        .all(|before| *before == LoadBeforeStop::No));
    let later = PrefetchSubject::unflipped(other, parsed[other as usize].as_deref());
    assert_eq!(
        load_before_stop(&db, &[later, tie, stopping]).at_first_stop(),
        [
            LoadBeforeStop::No,
            load_before_stop(&db, &[tie, stopping]).at_first_stop()[0],
            LoadBeforeStop::No
        ]
    );
    // NCBI may never look the stopping id up (reviews 1 and 2 of line DBID loadord, as
    // DECISIONS 60 (6) to (8) decided them): `ldg`'s `gi|1|ref|NP_000002.1|` holds it, and its
    // load registers it in NCBI's data loader: refused where the batch loads that sequence
    // before the stop (its gi's packed number 1 is smaller), or where the run loaded it
    // before (`TaxdbRun::note_loaded`: a batch, query or round before; review 2's
    // `B.loaded`), or after calls that nothing recorded (the library's batches); not where no
    // batch loaded it (review 2, M1: `B.never`, `P.first.it2`).
    let unsure = |db: &SeqDb,
                  run: &TaxdbRun,
                  subjects: &[PrefetchSubject<'_>],
                  queries: PrefetchQueries<'_>,
                  stop: u32,
                  parsed: &crate::report::db_subject::ParsedSubject| {
        stop_lookup_unsure(
            db,
            run,
            subjects,
            &load_before_stop(db, subjects).at_first_stop(),
            &queries,
            stop,
            parsed,
        )
        .unwrap()
        .map(|unsure| unsure.why)
    };
    let db = open(&["lda", "ldg", "ldp"]);
    let parsed = headers_of(&db);
    let stop = accession(&parsed, "NP_000002.1");
    let holder = parsed
        .iter()
        .position(|header| {
            header
                .as_ref()
                .is_some_and(|subject| subject.prefetch_stops().is_none())
        })
        .unwrap() as u32;
    let stopping = PrefetchSubject::unflipped(stop, parsed[stop as usize].as_deref());
    let held = PrefetchSubject::unflipped(holder, parsed[holder as usize].as_deref());
    let unparsed = PrefetchSubject::unflipped(0, None);
    assert!(parsed[0].is_none());
    let stop_parsed = parsed[stop as usize].as_ref().unwrap();
    let why = unsure(
        &db,
        &TaxdbRun::default(),
        &[held, stopping],
        stop_queries(0),
        stop,
        stop_parsed,
    )
    .unwrap();
    assert!(why.contains("a load of the batch's prefetch"), "{why}");
    assert_eq!(
        unsure(
            &db,
            &TaxdbRun::default(),
            &[stopping],
            stop_queries(0),
            stop,
            stop_parsed
        ),
        None
    );
    // A run that reported the holder before (a batch, query or round before; the record of
    // a run that can stop: `lda` has no accession table).
    let loaded = TaxdbRun::default();
    loaded.note_loaded(&db, || [unparsed, held]);
    assert_eq!(loaded.loaded(), [holder]);
    let why = unsure(
        &db,
        &loaded,
        &[stopping],
        stop_queries(1),
        stop,
        stop_parsed,
    )
    .unwrap();
    assert!(why.contains("the run loaded that sequence"), "{why}");
    // A run that reported only sequences without parsed Seq-ids: nothing recorded, NCBI looks
    // the stopping id up (M1).
    let never = TaxdbRun::default();
    never.note_loaded(&db, || [unparsed]);
    assert!(never.loaded().is_empty());
    assert_eq!(
        unsure(&db, &never, &[stopping], stop_queries(1), stop, stop_parsed),
        None
    );
    // A call of the library's batches after others, whose reports nothing recorded.
    let why = unsure(
        &db,
        &never,
        &[stopping],
        stop_queries(1).after_unrecorded_calls(true),
        stop,
        stop_parsed,
    )
    .unwrap();
    assert!(
        why.contains("whose reports LOSAT has no record of"),
        "{why}"
    );
    // A copy of the stopping sequence with the same Seq-align id is the same lookup: no
    // batch before can have loaded it without stopping there (`ldc`, a copy of `ldp`'s
    // record; kit case `copy.21b` of line DBID mixdb2's review 1, the rows `.copy.6`).
    let db = open(&["ldp", "ldc", "lda"]);
    let parsed = headers_of(&db);
    let stop = accession(&parsed, "NP_000002.1");
    let copy = parsed
        .iter()
        .enumerate()
        .rposition(|(ordinal, header)| {
            ordinal != stop as usize
                && header.as_ref().is_some_and(|subject| {
                    String::from_utf8_lossy(&subject.alignment_fasta()).contains("NP_000002.1")
                })
        })
        .unwrap() as u32;
    let stopping = PrefetchSubject::unflipped(stop, parsed[stop as usize].as_deref());
    let copied = TaxdbRun::default();
    copied.note_loaded(&db, || {
        [PrefetchSubject::unflipped(
            copy,
            parsed[copy as usize].as_deref(),
        )]
    });
    assert_eq!(copied.loaded(), [copy]);
    assert_eq!(
        unsure(
            &db,
            &copied,
            &[stopping],
            stop_queries(1),
            stop,
            parsed[stop as usize].as_ref().unwrap()
        ),
        None
    );
    // Review 4 of line DBID loadord, M1: the run's record is made only where a prefetch of
    // the run can stop, and the gate is asked first: a run that cannot stop (`ldp` alone, a
    // version 5 database with the accession table) never makes the subjects it is given.
    let closed = open(&["ldp"]);
    let run = TaxdbRun::default();
    run.note_loaded(&closed, || -> [PrefetchSubject<'_>; 0] {
        panic!("the subjects of a run that cannot stop were made")
    });
    assert!(!run.prefetch_can_stop(&closed));
    assert!(run.loaded().is_empty());
    // (The gate is decided once for a run: another run for another database.)
    let open_gate = open(&["lda", "ldp"]);
    let run = TaxdbRun::default();
    let mut made = 0;
    run.note_loaded(&open_gate, || {
        made += 1;
        std::iter::empty()
    });
    assert!(run.prefetch_can_stop(&open_gate) && made == 1);
    // Review 4, L3: `ldd` and `ldh` hold one Seq-align id, gi|9000, and the stopping
    // accession: a load of the later one (recorded) by that id loads the first (NCBI's
    // `x_GetOid` takes the lookup's first OID), so the sequence that NCBI may give for the
    // stopping id instead is `ldd`, never `ldh` itself.
    let db = open(&["lda", "ldd", "ldh", "ldp"]);
    let parsed = headers_of(&db);
    let stop = accession(&parsed, "NP_000002.1");
    assert_eq!(db.oid(stop).get(), 3);
    let recorded = TaxdbRun::default();
    recorded.note_loaded(&db, || {
        [PrefetchSubject::unflipped(2, parsed[2].as_deref())]
    });
    assert_eq!(recorded.loaded(), [2]);
    let stopping = [PrefetchSubject::unflipped(
        stop,
        parsed[stop as usize].as_deref(),
    )];
    let given = stop_lookup_unsure(
        &db,
        &recorded,
        &stopping,
        &[LoadBeforeStop::No],
        &stop_queries(1),
        stop,
        parsed[stop as usize].as_ref().unwrap(),
    )
    .unwrap()
    .unwrap();
    assert_eq!(given.instead, Some(vec![1]), "{}", given.why);
    assert!(given.why.contains("(OID 1)"), "{}", given.why);
    // `ldq`'s `lcl|Query_1`: NCBI's scope finds it where it holds a query with that id
    // (`QueryScope`, the owner of the queries' ids), never in the database.
    let db = open(&["lda", "ldq"]);
    let parsed = headers_of(&db);
    let stop = parsed
        .iter()
        .position(|header| {
            header
                .as_ref()
                .is_some_and(|s| s.prefetch_stops().is_some())
        })
        .unwrap() as u32;
    let stopping = PrefetchSubject::unflipped(stop, parsed[stop as usize].as_deref());
    let stop_parsed = parsed[stop as usize].as_ref().unwrap();
    let scope_refuses = |scope: &dyn Fn() -> crate::blastinput::query_batch::QueryScope| {
        unsure(
            &db,
            &TaxdbRun::default(),
            &[stopping],
            PrefetchQueries::new(0, scope),
            stop,
            stop_parsed,
        )
        .inspect(|why| assert!(why.contains("NCBI's scope holds a query"), "{why}"))
        .is_some()
    };
    use crate::blastinput::query_batch::QueryScope;
    let records = |fasta: &[u8], parse_deflines: bool| {
        let config = crate::blastinput::fasta_reader::ReaderConfig::query("BLASTP", true, false)
            .with_parse_deflines(parse_deflines);
        let mut source =
            crate::blastinput::fasta_reader::FastaInputSource::from_bytes(fasta, config);
        crate::blastinput::fasta_reader::read_all(&mut source, &mut |_| Ok(())).unwrap()
    };
    let letters = "A".repeat(150);
    let two = records(
        format!(">qa\n{letters}\n>qb\n{letters}\n").as_bytes(),
        false,
    );
    // The batch of query 0 (its generated id `Query_1`, `FastaRecord::local_id`), of both, or
    // of query 1 alone (the scope of a later batch, outside outfmt 5, which keeps the batches
    // before: `ResetScopeHistory`).
    assert!(scope_refuses(&|| QueryScope::of_shown_records(&two[..1])));
    assert!(scope_refuses(&|| QueryScope::of_shown_records(&two)));
    assert!(!scope_refuses(&|| QueryScope::of_shown_records(&two[1..])));
    // A record that `-query_loc` drops (40 letters, `-query_loc 50-140`) keeps its number and
    // stays in its batch's scope (review 2, H2): the batch's input records, not its searched
    // queries, make the scope.
    let short = records(
        format!(">qa\n{}\n>qb\n{letters}\n", "A".repeat(40)).as_bytes(),
        false,
    );
    let cut = crate::blastinput::seq_range::cut_queries(
        &short,
        &crate::blastinput::seq_range::SequenceRange { from: 50, to: 140 },
    );
    let batch = cut.input.batches(u32::MAX).remove(0);
    assert_eq!((batch.input.clone(), batch.searched.clone()), (0..2, 0..1));
    assert_eq!(cut.records[0].local_id, "Query_2");
    assert!(scope_refuses(&|| QueryScope::of_shown_records(
        &short[batch.input.clone()]
    )));
    // With `-parse_deflines` the scope holds the ids of the deflines, and no generated id
    // (review 2, M2: NCBI looks `lcl|Query_1` up and stops).
    let believed = records(
        format!(">lcl|qa\n{letters}\n>lcl|qb\n{letters}\n").as_bytes(),
        true,
    );
    assert!(believed[0].believed().is_some());
    assert!(!scope_refuses(&|| QueryScope::of_shown_records(&believed)));
    // The PSSM's query of `-in_pssm` (review 2, H4): its Bioseq's ids, `lcl|Query_1` or
    // `lcl|myprot`.
    let local = |id: &str| {
        crate::seqid::ParsedSeqId::Local(crate::seqid::model::ObjectId::Str(id.as_bytes().to_vec()))
    };
    assert!(scope_refuses(&|| QueryScope::of_bioseq_ids(&[local(
        "Query_1"
    )])));
    assert!(!scope_refuses(&|| QueryScope::of_bioseq_ids(&[local(
        "myprot"
    )])));
    // Local ids of any case are one handle (`lcl|QUERY_1`); a leading zero is another id.
    assert!(scope_refuses(&|| QueryScope::of_bioseq_ids(&[local(
        "QUERY_1"
    )])));
    assert!(!scope_refuses(&|| QueryScope::of_bioseq_ids(&[local(
        "Query_01"
    )])));
    std::fs::remove_dir_all(&cwd).unwrap();
}

// Review 3 of line DBID loadord, H1 (its kit's `R.*` cases): NCBI's loads end at the first
// lookup that throws, in the order of the Seq-id handles (`LoadSequencesToScope` rethrows the
// `CSeqDBException`); a later stopping id that NCBI's scope or data loader would find without
// a lookup is never reached. A stop before it that NCBI may skip leaves the loads up to the
// stop that NCBI surely makes, and the sequence that the skipped lookup is given instead, of
// an unknown place: refused only where one of them writes a text.
#[test]
fn a_batch_with_two_stopping_lookups_ends_at_the_first_that_ncbi_makes() {
    // `multidb` (make_multidb.py --loadord-only): `lda` (no accession table) beside `ldp`'s
    // `ref|NP_000002.1|` (packed 2) and `ref|NP_000900.1|` (packed 900), `ldq`'s
    // `lcl|Query_1` (not packed), `ldg`'s `gi|1|ref|NP_000002.1|` and the unparsed `ldu`.
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-two-stops-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let open = |names: &[&str], taxdb: Option<&Path>| {
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
        let list: Vec<String> = names
            .iter()
            .map(|name| {
                root()
                    .join("db_search/multidb")
                    .join(name)
                    .to_str()
                    .unwrap()
                    .to_owned()
            })
            .collect();
        SeqDb::open_in(&list.join(" "), SeqType::Protein, &search).unwrap()
    };
    type Headers = [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>];
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db);
            (0..db.num_seqs())
                .map(
                    |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                        ReportedHeader::Parsed(subject) => Some(subject),
                        _ => None,
                    },
                )
                .collect()
        };
    let find = |parsed: &Headers, text: &str, stops: bool| {
        parsed
            .iter()
            .position(|header| {
                header.as_ref().is_some_and(|subject| {
                    subject.prefetch_stops().is_some() == stops
                        && String::from_utf8_lossy(&subject.alignment_fasta()).contains(text)
                })
            })
            .unwrap() as u32
    };
    let stopped = |result: Result<crate::report::diag_log::Posts>| -> (usize, String) {
        let failure = result.unwrap_err().downcast::<HeaderReadError>().unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (failure.query, String::from_utf8(native.message).unwrap())
    };
    let refused = |result: Result<crate::report::diag_log::Posts>| -> String {
        let refusal = result.unwrap_err();
        assert!(!refusal.is::<HeaderReadError>(), "{refusal}");
        refusal.to_string()
    };
    let error = "BLAST Database error: DB contains no accession info.\n".to_owned();
    // 1. The first stop in NCBI's order is one that it makes: `ref|NP_000002.1|` before
    // `lcl|Query_1`, which the scope holds (kit `R.q1`), in either order of the batch.
    let db = open(&["lda", "ldp", "ldq"], None);
    let parsed = headers_of(&db);
    let s2 = find(&parsed, "NP_000002.1", true);
    let q1 = find(&parsed, "Query_1", true);
    fn subject(parsed: &Headers, ordinal: u32) -> PrefetchSubject<'_> {
        PrefetchSubject::unflipped(ordinal, parsed[ordinal as usize].as_deref())
    }
    let local = |id: &str| {
        crate::seqid::ParsedSeqId::Local(crate::seqid::model::ObjectId::Str(id.as_bytes().to_vec()))
    };
    let scope = || crate::blastinput::query_batch::QueryScope::of_bioseq_ids(&[local("Query_1")]);
    let queries = PrefetchQueries::new(0, &scope);
    for subjects in [
        [subject(&parsed, q1), subject(&parsed, s2)],
        [subject(&parsed, s2), subject(&parsed, q1)],
    ] {
        assert_eq!(
            stopped(batch_prefetch(
                &db,
                &TaxdbRun::default(),
                std::iter::empty(),
                subjects,
                "BLASTP",
                &queries
            )),
            (0, error.clone())
        );
    }
    // `lcl|Query_1` alone: no stop that NCBI surely makes, refused (NCBI goes on).
    let why = refused(batch_prefetch(
        &db,
        &TaxdbRun::default(),
        std::iter::empty(),
        [subject(&parsed, q1)],
        "BLASTP",
        &queries,
    ));
    assert!(
        why.contains("may not look up: NCBI's scope holds a query"),
        "{why}"
    );
    // 2. The first stop may be skipped: a batch before loaded `gi|1|ref|NP_000002.1|` (the
    // run's record), whose load registered the accession, so NCBI may look
    // `ref|NP_000900.1|` up next and stop there (kit `R.rev`). With `bad_magic`, the loads
    // write the texts of opening the taxonomy database in a run that has not opened it.
    let bad_magic = root().join("seqdb/taxdb/bad_magic");
    for taxdb in [None, Some(bad_magic.as_path())] {
        let db = open(&["lda", "ldg", "ldp", "ldu"], taxdb);
        let parsed = headers_of(&db);
        let s2 = find(&parsed, "NP_000002.1", true);
        let s900 = find(&parsed, "NP_000900.1", true);
        let holder = find(&parsed, "gi|1", false);
        // `ldu`'s `gnl|BL_ORD_ID|4` (packed 4): between the two stops.
        let between = parsed.iter().rposition(Option::is_none).unwrap() as u32;
        assert_eq!(db.oid(between).get(), 4);
        let both = [subject(&parsed, s900), subject(&parsed, s2)];
        let with_between = [
            subject(&parsed, s2),
            PrefetchSubject::unflipped(between, None),
            subject(&parsed, s900),
        ];
        // The run that loaded the holder and looked its taxonomy up (its batch's prefetch).
        let looked_up = || {
            let run = TaxdbRun::default();
            run.note_loaded(&db, || [subject(&parsed, holder)]);
            let opened = prefetch_messages(
                &db,
                &run,
                [(holder, parsed[holder as usize].as_deref())],
                "BLASTP",
            )
            .unwrap();
            assert_eq!(!opened.is_empty(), taxdb.is_some());
            run
        };
        // The holder's load by the skipped id, again: nothing written, NCBI's error.
        assert_eq!(
            stopped(batch_prefetch(
                &db,
                &looked_up(),
                std::iter::empty(),
                both,
                "BLASTP",
                &stop_queries(1)
            )),
            (1, error.clone())
        );
        // A run whose record holds the holder but which has opened nothing (a holder loaded
        // by a call before, the library's): that load, or the load between the stops, writes
        // the texts of opening, which NCBI writes only where it skips the first lookup.
        let fresh = || {
            let run = TaxdbRun::default();
            run.note_loaded(&db, || [subject(&parsed, holder)]);
            run
        };
        for (subjects, what) in [
            (
                &both[..],
                "the sequence that NCBI's BLAST database data loader gives for a stopping id",
            ),
            (
                &with_between[..],
                "its Seq-id handle comes after the first stopping one's",
            ),
        ] {
            let result = batch_prefetch(
                &db,
                &fresh(),
                std::iter::empty(),
                subjects.iter().copied(),
                "BLASTP",
                &stop_queries(1),
            );
            match taxdb {
                None => assert_eq!(stopped(result), (1, error.clone())),
                Some(_) => {
                    let why = refused(result);
                    assert!(
                        why.contains("makes the taxonomy database write a text"),
                        "{why}"
                    );
                    assert!(why.contains(what), "{why}");
                }
            }
        }
        // 3. No stop that NCBI surely makes: `ref|NP_000002.1|` alone, refused.
        let why = refused(batch_prefetch(
            &db,
            &looked_up(),
            std::iter::empty(),
            [subject(&parsed, s2)],
            "BLASTP",
            &stop_queries(1),
        ));
        assert!(why.contains("the run loaded that sequence"), "{why}");
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:756-758 (cited at
// `report::query_warnings::write_posts_to`: a cut `Iteration` makes none of its loads)
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3299 and
// c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-396 (cited at `load_before_stop` and
// `prefetch_stop`: the loads in the order of the handles, up to the first stopping lookup)
/// Line RG K4h's E2 in the programs that model a stopping batch before they write the BLAST XML
/// reports before it (review 5 of line DBID loadord, M1): the report of query 0, which NCBI's
/// stream cut before its hits, left the texts of opening the taxonomy database, which the first
/// load before the stop that asks the database writes before NCBI's error (NCBI 2.17.0,
/// review 5 `E2.sure`, `E2.skip`); a load after the stop is not made (`E2.after`: the error
/// alone); one of unknown place refuses the batch; a load of the reports or prefetches between
/// takes them first.
#[test]
fn the_texts_that_a_cut_report_left_go_to_the_first_asking_load_before_the_stop() {
    // `multidb` (make_multidb.py --loadord-only, --loadord2-only, --loadord3-only,
    // --loadord4-only; taxid 9606): `ldb`'s unparsed OID 1 (`gnl|BL_ORD_ID|1`, before
    // `ldp`'s stopping `ref|NP_000002.1|`), `lda`'s OID 0 (`|0`, after it), `ldu`'s OID 2
    // (`|2`, the same packed number), `ldd`'s `gi|9000|ref|NP_000002.1|` (a holder of the
    // stopping accession) and `lds`'s stopping `ref|NP_000003.1|`.
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-e2-placed-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let bad_magic = root().join("seqdb/taxdb/bad_magic");
    let open = |names: &[&str]| {
        let search = SearchPath::new(&cwd, Some(bad_magic.as_os_str()), None).unwrap();
        let list: Vec<String> = names
            .iter()
            .map(|name| {
                root()
                    .join("db_search/multidb")
                    .join(name)
                    .to_str()
                    .unwrap()
                    .to_owned()
            })
            .collect();
        SeqDb::open_in(&list.join(" "), SeqType::Protein, &search).unwrap()
    };
    type Headers = [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>];
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db);
            (0..db.num_seqs())
                .map(
                    |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                        ReportedHeader::Parsed(subject) => Some(subject),
                        _ => None,
                    },
                )
                .collect()
        };
    let find = |parsed: &Headers, text: &str| {
        parsed
            .iter()
            .position(|header| {
                header.as_ref().is_some_and(|subject| {
                    String::from_utf8_lossy(&subject.alignment_fasta()).contains(text)
                })
            })
            .unwrap() as u32
    };
    fn subject(parsed: &Headers, ordinal: u32) -> PrefetchSubject<'_> {
        PrefetchSubject::unflipped(ordinal, parsed[ordinal as usize].as_deref())
    }
    let stopped = |result: Result<crate::report::diag_log::Posts>| -> (usize, String) {
        let failure = result.unwrap_err().downcast::<HeaderReadError>().unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (failure.query, String::from_utf8(native.message).unwrap())
    };
    let error = "BLAST Database error: DB contains no accession info.\n".to_owned();
    // The model of the cut report of query 0 in `run`: its reads of `read` (each sequence's
    // first load in the run; the first writes the texts of opening).
    let cut_report = |db: &SeqDb, run: &TaxdbRun, read: &[u32]| {
        let reads = subject_read_texts(
            db,
            run,
            read.iter().map(|&ordinal| (ordinal, ordinal, None)),
            false,
            "BLASTP",
        )
        .unwrap();
        assert!(reads.iter().all(|read| read.load.is_some()));
        reads
    };
    // The batch of query `query` after the outfmt 5 reports `before`.
    fn after_reports<'a>(query: usize, before: &'a XmlReportsBefore<'a>) -> PrefetchQueries<'a> {
        PrefetchQueries::new(query, &crate::blastinput::query_batch::QueryScope::default)
            .after_xml_reports(Some(before))
    }
    let none = |_: usize| None;
    let cut = |query: usize| query == 0;
    let uncut = |_: usize| false;
    // 1. Before the stop (`E2.sure`): `|1` loads before `ref|NP_000002.1|` and writes the
    // texts of opening that query 0's report left.
    let db = open(&["ldb", "ldp"]);
    let parsed = headers_of(&db);
    let s2 = find(&parsed, "NP_000002.1");
    assert_eq!(db.oid(1).get(), 1);
    let batch = [PrefetchSubject::unflipped(1, None), subject(&parsed, s2)];
    let run = TaxdbRun::default();
    let reads = vec![cut_report(&db, &run, &[1])];
    let open_text = String::from_utf8(reads[0][0].text.clone()).unwrap();
    assert!(open_text.contains("wrong magic number"), "{open_text}");
    let before = XmlCuts::of(true, &cut)
        .reports_before(&reads, &none)
        .unwrap();
    assert_eq!(
        stopped(batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            batch,
            "BLASTP",
            &after_reports(1, &before)
        )),
        (1, format!("{open_text}{error}"))
    );
    // The report not cut wrote them itself; without outfmt 5 first, the same.
    let made = XmlCuts::of(true, &uncut)
        .reports_before(&reads, &none)
        .unwrap();
    for queries in [
        after_reports(1, &made),
        PrefetchQueries::new(1, &crate::blastinput::query_batch::QueryScope::default),
    ] {
        assert_eq!(
            stopped(batch_prefetch(
                &db,
                &run,
                std::iter::empty(),
                batch,
                "BLASTP",
                &queries
            )),
            (1, error.clone())
        );
    }
    // A batch between whose prefetch loads `|1` (its posts before query 1's report) takes
    // them first: the stop of query 2 writes NCBI's error alone.
    let mut between = crate::report::diag_log::Posts::new();
    between.loads(subject_read_texts(&db, &run, [(1, 1, None)], false, "BLASTP").unwrap());
    let prefetched = |query: usize| (query == 1).then_some(&between);
    let taken = XmlCuts::of(true, &cut)
        .reports_before(&reads, &prefetched)
        .unwrap();
    assert_eq!(
        stopped(batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            batch,
            "BLASTP",
            &after_reports(2, &taken)
        )),
        (2, error.clone())
    );
    // 2. After the stop (`E2.after`): `|0` loads after it, which NCBI never makes; the texts
    // are not written.
    let db = open(&["lda", "ldp"]);
    let parsed = headers_of(&db);
    let s2 = find(&parsed, "NP_000002.1");
    assert_eq!(db.oid(0).get(), 0);
    let run = TaxdbRun::default();
    let reads = vec![cut_report(&db, &run, &[0])];
    let before = XmlCuts::of(true, &cut)
        .reports_before(&reads, &none)
        .unwrap();
    assert_eq!(
        stopped(batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [PrefetchSubject::unflipped(0, None), subject(&parsed, s2)],
            "BLASTP",
            &after_reports(1, &before)
        )),
        (1, error.clone())
    );
    // 3. Of unknown place: `|2` has the stop's packed number 2, so NCBI's order is that of
    // the handles' addresses; where it loads first it writes the texts: refused.
    let db = open(&["ldp", "ldu"]);
    let parsed = headers_of(&db);
    let s2 = find(&parsed, "NP_000002.1");
    let run = TaxdbRun::default();
    let reads = vec![cut_report(&db, &run, &[2])];
    let before = XmlCuts::of(true, &cut)
        .reports_before(&reads, &none)
        .unwrap();
    let refusal = batch_prefetch(
        &db,
        &run,
        std::iter::empty(),
        [PrefetchSubject::unflipped(2, None), subject(&parsed, s2)],
        "BLASTP",
        &after_reports(1, &before),
    )
    .unwrap_err();
    assert!(!refusal.is::<HeaderReadError>(), "{refusal}");
    let why = refusal.to_string();
    assert!(
        why.contains("same packed number")
            && why.contains("for what a BLAST XML report before it did not"),
        "{why}"
    );
    // 4. One skipped stop: the run's record holds the holder `gi|9000|...` (a load before
    // the cut report's, here set by hand: the cut report's own reads load nothing,
    // `XmlCuts`, review 6 H2), so `ref|NP_000002.1|` may be skipped and `ref|NP_000003.1|`
    // is the sure stop; `|1` before both writes the texts that the cut report left, the
    // holder's load by the skipped id writes nothing more.
    let db = open(&["ldb", "ldd", "ldp", "lds"]);
    let parsed = headers_of(&db);
    let s2 = find(&parsed, "ref|NP_000002.1|");
    let s3 = find(&parsed, "NP_000003.1");
    let holder = find(&parsed, "gi|9000");
    assert!(parsed[s2 as usize]
        .as_ref()
        .unwrap()
        .prefetch_stops()
        .is_some());
    let run = TaxdbRun::default();
    let reads = vec![subject_read_texts(
        &db,
        &run,
        [
            (1, 1, None),
            (holder, holder, parsed[holder as usize].as_deref()),
        ],
        false,
        "BLASTP",
    )
    .unwrap()];
    run.note_loaded(&db, || [subject(&parsed, holder)]);
    let open_text = String::from_utf8(reads[0][0].text.clone()).unwrap();
    let before = XmlCuts::of(true, &cut)
        .reports_before(&reads, &none)
        .unwrap();
    assert_eq!(
        stopped(batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [
                subject(&parsed, s3),
                PrefetchSubject::unflipped(1, None),
                subject(&parsed, s2),
            ],
            "BLASTP",
            &after_reports(1, &before)
        )),
        (1, format!("{open_text}{error}"))
    );
    std::fs::remove_dir_all(&cwd).unwrap();
}

// Review 6 of line DBID loadord (its kit's `n6.ps.H.*` and `n6.*.K.cut.hold.5`): the run's
// record of what NCBI's data loader loaded for the batches before a stopping batch. H1: a
// batch whose prefetch `PRE_FETCH_SEQS_LIMIT` turns off loads what its reports read after
// the cut to `-max_target_seqs` (`BatchLoads::Read`, also where its rows read every
// subject), not psiblast's traceback beyond the cut, which only its prefetch loads
// (`BatchLoads::Prefetched`). H2: an outfmt 5 report that NCBI's XML stream cut before its
// hits loads nothing (`XmlCuts`, the cut that the texts it left ask too), so the holder that
// its model reads leaves the stopping accession sure.
#[test]
fn the_record_before_a_stop_holds_only_what_ncbis_reports_and_prefetches_loaded() {
    use crate::report::tabular::SubjectFetch::Every;
    // `multidb` (make_multidb.py --loadord-only): `ldb`'s unparsed OID 1 (`gnl|BL_ORD_ID|1`),
    // `ldd`'s `gi|9000|ref|NP_000002.1|` (a holder of the stopping accession), `ldp`'s
    // stopping `ref|NP_000002.1|` (no accession table); taxid 9606, BLASTDB=bad_magic.
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-record-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let bad_magic = root().join("seqdb/taxdb/bad_magic");
    let search = SearchPath::new(&cwd, Some(bad_magic.as_os_str()), None).unwrap();
    let list: Vec<String> = ["ldb", "ldd", "ldp"]
        .iter()
        .map(|name| {
            root()
                .join("db_search/multidb")
                .join(name)
                .to_str()
                .unwrap()
                .to_owned()
        })
        .collect();
    let db = SeqDb::open_in(&list.join(" "), SeqType::Protein, &search).unwrap();
    let headers = HeaderReader::new(&db);
    let parsed: Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> = (0..db
        .num_seqs())
        .map(
            |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                ReportedHeader::Parsed(subject) => Some(subject),
                _ => None,
            },
        )
        .collect();
    let find = |text: &str| {
        parsed
            .iter()
            .position(|header| {
                header.as_ref().is_some_and(|subject| {
                    String::from_utf8_lossy(&subject.alignment_fasta()).contains(text)
                })
            })
            .unwrap() as u32
    };
    let (s2, holder) = (find("ref|NP_000002.1|"), find("gi|9000"));
    let subject =
        |ordinal: u32| PrefetchSubject::unflipped(ordinal, parsed[ordinal as usize].as_deref());
    // Query 0 (batch 0..1): its prefetched set (psiblast's traceback) holds `|1` and the
    // holder; its report shows `shown`.
    struct Query {
        traceback: Vec<u32>,
        shown: Vec<u32>,
    }
    let record = |query: &Query, limit: Option<i32>, cuts: XmlCuts<'_>| -> Vec<u32> {
        LoadedBefore {
            hits: std::slice::from_ref(query),
            spans: &[0..1],
            rows: Every,
            cuts,
        }
        .subjects_at(
            limit,
            |hits, span| hits.get(span),
            |_, hits| hits.iter().map(|query| query.traceback.len()).sum(),
            |_, hits, loads| {
                hits.iter()
                    .flat_map(|query| match loads {
                        BatchLoads::Prefetched => query.traceback.clone(),
                        BatchLoads::Read(_) => query.shown.clone(),
                    })
                    .collect::<Vec<u32>>()
            },
        )
        .collect()
    };
    let error = "BLAST Database error: DB contains no accession info.\n".to_owned();
    // The stop of query 1's batch (`ref|NP_000002.1|`) after the record `loaded`.
    let stop = |loaded: &[u32], queries: &PrefetchQueries<'_>| {
        batch_prefetch(
            &db,
            &TaxdbRun::default(),
            loaded.iter().map(|&ordinal| subject(ordinal)),
            [subject(s2)],
            "BLASTP",
            queries,
        )
        .unwrap_err()
    };
    let native = |error: anyhow::Error| -> (usize, String) {
        let failure = error.downcast::<HeaderReadError>().unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (failure.query, String::from_utf8(native.message).unwrap())
    };
    let refused = |error: anyhow::Error| {
        assert!(!error.is::<HeaderReadError>(), "{error}");
        let why = error.to_string();
        assert!(why.contains("the run loaded that sequence"), "{why}");
    };
    // H1: `-max_target_seqs 1` cuts the report to `|1`; the prefetch off (limit 1 < 2
    // Seq-aligns), its rows read every subject they show: NCBI never loads the holder, looks
    // `ref|NP_000002.1|` up and stops.
    let mts1 = Query {
        traceback: vec![1, holder],
        shown: vec![1],
    };
    assert_eq!(record(&mts1, Some(1), XmlCuts::NONE), [1]);
    assert_eq!(
        native(stop(
            &record(&mts1, Some(1), XmlCuts::NONE),
            &stop_queries(1)
        )),
        (1, error.clone())
    );
    // The prefetch on loads the whole traceback, the holder too: NCBI gives the holder for
    // the accession and goes on, which LOSAT refuses.
    assert_eq!(record(&mts1, None, XmlCuts::NONE), [1, holder]);
    refused(stop(&record(&mts1, None, XmlCuts::NONE), &stop_queries(1)));
    // H2: query 0's outfmt 5 report shows both, but NCBI's XML stream cuts its Iteration
    // before its hits: it loads neither, so the stop is sure; the texts of opening that the
    // model of its reads left go to no load before the stop (NCBI writes its error alone).
    let both = Query {
        traceback: vec![1, holder],
        shown: vec![1, holder],
    };
    let cut = |query: usize| query == 0;
    let cuts = XmlCuts::of(true, &cut);
    assert!(record(&both, Some(1), cuts).is_empty());
    let run = TaxdbRun::default();
    let reads = vec![subject_read_texts(
        &db,
        &run,
        [
            (1, 1, None),
            (holder, holder, parsed[holder as usize].as_deref()),
        ],
        false,
        "BLASTP",
    )
    .unwrap()];
    assert!(String::from_utf8_lossy(&reads[0][0].text).contains("wrong magic number"));
    let none = |_: usize| None;
    let before = cuts.reports_before(&reads, &none).unwrap();
    let queries = PrefetchQueries::new(1, &crate::blastinput::query_batch::QueryScope::default)
        .after_xml_reports(Some(&before));
    assert_eq!(
        native(stop(&record(&both, Some(1), cuts), &queries)),
        (1, error.clone())
    );
    // The report written (no outfmt 5 first: no cut, no reports before to ask) loads the
    // holder: refused as NCBI goes on.
    let uncut = XmlCuts::of(false, &cut);
    assert!(uncut.reports_before(&reads, &none).is_none());
    assert_eq!(record(&both, Some(1), uncut), [1, holder]);
    refused(stop(&record(&both, Some(1), uncut), &stop_queries(1)));
    // A cut does not touch a batch whose prefetch ran: it loaded both before the report.
    assert_eq!(record(&both, None, cuts), [1, holder]);
    std::fs::remove_dir_all(&cwd).unwrap();
}

/// The batch whose stop names the query `query` and whose scope holds no query (a run whose
/// queries have none of the ids of the database).
fn stop_queries(query: usize) -> PrefetchQueries<'static> {
    PrefetchQueries::new(query, &crate::blastinput::query_batch::QueryScope::default)
}

#[test]
fn a_lookup_that_fails_on_a_file_stops_the_batchs_prefetch_as_ncbis() {
    // Line DBID mixdb2: `multidb/na` (-parse_seqids) beside `multidb/nnp` (no accession
    // table), whose parsed Seq-ids' lookups fail with NCBI's `DB contains no accession
    // info.`, and `taxfilter/tfu_n` (-input_type blastdb: the table of its ordinal ids; taxids
    // 7955 and 9606), whose loads look the taxonomy up.
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-prefetch-stop-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let open = |names: &[&str], taxdb: Option<&Path>| {
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None).unwrap();
        let list: Vec<String> = names
            .iter()
            .map(|name| root().join(name).to_str().unwrap().to_owned())
            .collect();
        SeqDb::open_in(&list.join(" "), SeqType::Nucleotide, &search).unwrap()
    };
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db);
            (0..db.num_seqs())
                .map(
                    |ordinal| match read_header(&headers, ordinal, 0, "BLASTN").unwrap() {
                        ReportedHeader::Parsed(subject) => Some(subject),
                        _ => None,
                    },
                )
                .collect()
        };
    fn subject(
        ordinal: u32,
        parsed: &Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>,
    ) -> PrefetchSubject<'_> {
        PrefetchSubject {
            ordinal,
            parsed: parsed.as_deref(),
            flipped: false,
        }
    }
    let stopped = |result: Result<()>| -> (usize, String) {
        let failure = result.unwrap_err().downcast::<HeaderReadError>().unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (failure.query, String::from_utf8(native.message).unwrap())
    };
    let error = "BLAST Database error: DB contains no accession info.\n".to_owned();
    let names = [
        "db_search/multidb/na",
        "db_search/multidb/nnp",
        "db_search/taxfilter/tfu_n",
    ];
    let bad_magic = root().join("seqdb/taxdb/bad_magic");
    for taxdb in [None, Some(bad_magic.as_path())] {
        let db = open(&names, taxdb);
        let parsed = headers_of(&db);
        let stop = parsed
            .iter()
            .position(|header| {
                header
                    .as_ref()
                    .is_some_and(|s| s.prefetch_stops().is_some())
            })
            .unwrap() as u32;
        // A sequence of tfu_n (the last volume by name): its load looks its taxid up.
        let other = db.num_seqs() - 1;
        assert!(parsed[other as usize].is_none());
        let open = String::from_utf8(
            prefetch_messages(&db, &TaxdbRun::default(), [(other, None)], "BLASTN").unwrap(),
        )
        .unwrap();
        assert_eq!(!open.is_empty(), taxdb.is_some());
        let both = || [subject(stop, &parsed[stop as usize]), subject(other, &None)];
        // The stopping subject alone: NCBI's error, of the batch's first query.
        let run = TaxdbRun::default();
        let alone = [subject(stop, &parsed[stop as usize])];
        assert_eq!(
            stopped(
                batch_prefetch(
                    &db,
                    &run,
                    std::iter::empty(),
                    alone,
                    "BLASTN",
                    &stop_queries(3)
                )
                .map(drop)
            ),
            (3, error.clone())
        );
        assert_eq!(
            stopped(batch_prefetch_stop(&db, &run, alone, "BLASTN", &stop_queries(3)).map(drop)),
            (3, error.clone())
        );
        // With another subject whose load writes the texts of opening the taxonomy database:
        // its `gnl|BL_ORD_ID|n` (packed) loads before the stopping `lcl|` id (not packed) in
        // NCBI's order, so its texts come before the error.
        let stop_id = parsed[stop as usize].as_ref().unwrap().alignment_id();
        assert_eq!(
            crate::seqid::mask_handle::mask_handle_packed(stop_id),
            0,
            "{stop_id:?}"
        );
        let run = TaxdbRun::default();
        match taxdb {
            None => {
                assert_eq!(
                    stopped(
                        batch_prefetch(
                            &db,
                            &run,
                            std::iter::empty(),
                            both(),
                            "BLASTN",
                            &stop_queries(0)
                        )
                        .map(drop)
                    )
                    .1,
                    error
                );
                assert_eq!(
                    stopped(
                        batch_prefetch_stop(&db, &run, both(), "BLASTN", &stop_queries(0))
                            .map(drop)
                    )
                    .1,
                    error
                );
                assert_eq!(
                    stopped(
                        held_batch_prefetch_stop(
                            &db,
                            &TaxdbRun::default(),
                            std::iter::empty(),
                            std::iter::empty(),
                            both(),
                            "BLASTN",
                            &stop_queries(0)
                        )
                        .map(drop)
                    )
                    .1,
                    error
                );
            }
            Some(_) => {
                // Each in a run that has looked nothing up: a stopped prefetch has recorded
                // its loads in its run.
                for stop in [
                    batch_prefetch(
                        &db,
                        &TaxdbRun::default(),
                        std::iter::empty(),
                        both(),
                        "BLASTN",
                        &stop_queries(0),
                    )
                    .map(drop),
                    batch_prefetch_stop(
                        &db,
                        &TaxdbRun::default(),
                        both(),
                        "BLASTN",
                        &stop_queries(0),
                    )
                    .map(drop),
                    held_batch_prefetch_stop(
                        &db,
                        &TaxdbRun::default(),
                        std::iter::empty(),
                        std::iter::empty(),
                        both(),
                        "BLASTN",
                        &stop_queries(0),
                    )
                    .map(drop),
                ] {
                    assert_eq!(stopped(stop), (0, format!("{open}{error}")));
                }
                // A batch after one whose prefetch opened it and looked the other subject up
                // (the formats whose prefetch writes the texts, `batch_prefetch`): its load
                // writes nothing now, so NCBI's error.
                let run = TaxdbRun::default();
                assert!(!batch_prefetch(
                    &db,
                    &run,
                    std::iter::empty(),
                    [subject(other, &None)],
                    "BLASTN",
                    &stop_queries(0)
                )
                .unwrap()
                .is_empty());
                assert_eq!(
                    stopped(
                        batch_prefetch(
                            &db,
                            &run,
                            std::iter::empty(),
                            both(),
                            "BLASTN",
                            &stop_queries(1)
                        )
                        .map(drop)
                    ),
                    (1, error.clone())
                );
                // The same after a batch whose XML2 or JSON reports looked it up: in the run
                // that holds their lookups (BLASTX, `batch_prefetch_stop`), and in the record
                // made where the batch stops of the subjects that the batches before report
                // (`held_batch_prefetch_stop`; with no batch before, refused above).
                let run = TaxdbRun::default();
                read_reported_taxonomy(
                    &db,
                    &run,
                    &[other],
                    &std::collections::BTreeSet::from([0]),
                    &[None::<&crate::report::db_subject::ParsedSubject>],
                    "BLASTN",
                )
                .unwrap();
                assert_eq!(
                    stopped(
                        batch_prefetch_stop(&db, &run, both(), "BLASTN", &stop_queries(1))
                            .map(drop)
                    ),
                    (1, error.clone())
                );
                assert_eq!(
                    stopped(
                        held_batch_prefetch_stop(
                            &db,
                            &TaxdbRun::default(),
                            [subject(other, &None)],
                            std::iter::empty(),
                            both(),
                            "BLASTN",
                            &stop_queries(1)
                        )
                        .map(drop)
                    ),
                    (1, error.clone())
                );
            }
        }
        // No stopping subject: nothing for the XML2 and JSON formats.
        batch_prefetch_stop(
            &db,
            &run,
            [subject(other, &None)],
            "BLASTN",
            &stop_queries(0),
        )
        .map(drop)
        .unwrap();
    }
    // Another subject whose lookup throws on its id, which the prefetch catches and which ends
    // its loads: whether the stop comes first is the order's, refused. (A text tag of
    // `ids/ntag_t5` is the second volume beside `multidb/nnp`, which the ordinal-tag rule of
    // line DBID ordsl refuses first: the throw is set on a copy of the stopping subject.)
    let db = open(&names[..2], None);
    let parsed = headers_of(&db);
    let stop = parsed
        .iter()
        .position(|header| {
            header
                .as_ref()
                .is_some_and(|s| s.prefetch_stops().is_some())
        })
        .unwrap();
    let throws = Some(std::sync::Arc::new(
        (**parsed[stop].as_ref().unwrap())
            .clone()
            .found_instead(Some(crate::report::db_subject::FoundInstead::ThrowsOnId(
                "Invalid choice selection".into(),
            ))),
    ));
    assert!(throws.as_ref().unwrap().prefetch_fails());
    let run = TaxdbRun::default();
    let other = parsed.len() as u32 - 1;
    for refused in [
        batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [subject(stop as u32, &parsed[stop]), subject(other, &throws)],
            "BLASTN",
            &stop_queries(0),
        )
        .map(drop),
        batch_prefetch_stop(
            &db,
            &run,
            [subject(other, &throws), subject(stop as u32, &parsed[stop])],
            "BLASTN",
            &stop_queries(0),
        )
        .map(drop),
    ] {
        let refused = refused.unwrap_err();
        assert!(!refused.is::<HeaderReadError>());
        assert!(
            refused.to_string().contains("whose lookup throws"),
            "{refused}"
        );
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn a_weak_match_is_refused_where_ncbi_may_read_the_reported_sequences_one_by_one() {
    // tests/db_search/ids/make_ids.py (line DBID N2c): `pweak_v5` and `pweak_v4` 0
    // `sp|P00041.1|NAMEW`, 1 `gb|FFF00040|`, 2 `pir|PIRW1|`, 3 `sp||NAMEX`, 4 `gb|FFF00041.1|`,
    // 5 `sp||NAMEY`, 6 `sp|P00042.1|NAMEY`, 7 `gb|FFF00043|`, 8 `gb|FFF00043.1|`, 9
    // `sp|P00044.1|NAMEZ`, 10 `sp||NAMEZ`. With the prefetch on, NCBI writes each sequence's
    // own report (in version 4 the lookup of 10 gives 9 first: refused either way, not listed).
    let cwd = empty_cwd("weak_match");
    // `None`: fetched; `Some(Some(n))`: another sequence n matches weakly; `Some(None)`: weak
    // keys that version 5 does not look up.
    let unlisted = Some(None);
    for (name, ordinal, want) in [
        ("ids/pweak_v5", 0, unlisted),
        ("ids/pweak_v5", 1, None),
        ("ids/pweak_v5", 2, unlisted),
        ("ids/pweak_v5", 3, unlisted),
        ("ids/pweak_v5", 4, None),
        ("ids/pweak_v5", 6, unlisted),
        ("ids/pweak_v5", 7, Some(Some(8))),
        ("ids/pweak_v5", 8, None),
        ("ids/pweak_v4", 0, None),
        ("ids/pweak_v4", 1, None),
        ("ids/pweak_v4", 2, None),
        ("ids/pweak_v4", 3, None),
        ("ids/pweak_v4", 5, Some(Some(6))),
        ("ids/pweak_v4", 6, Some(Some(5))),
        ("ids/pweak_v4", 7, Some(Some(8))),
        ("ids/pweak_v4", 9, Some(Some(10))),
    ] {
        let db = db_search(name, None, &cwd);
        for prefetch_off in [false, true] {
            let headers = HeaderReader {
                prefetch_off,
                ..HeaderReader::new(&db)
            };
            let ReportedHeader::Parsed(subject) =
                read_header(&headers, ordinal, 0, "BLASTP").unwrap()
            else {
                panic!("{name} {ordinal}: no parsed Seq-ids");
            };
            // A batch whose prefetch is on loads the subject by its own id.
            assert!(
                subject
                    .fetched(false, false, crate::report::db_subject::HitLoad::Loaded)
                    .is_ok(),
                "{name} {ordinal}"
            );
            match (
                subject.fetched(false, true, crate::report::db_subject::HitLoad::Loaded),
                want.filter(|_| prefetch_off),
            ) {
                (Ok(()), None) => {}
                (Err(error), Some(Some(other))) => assert!(
                    error.contains(&format!("match an id of database sequence {other} weakly")),
                    "{name} {ordinal}: {error}"
                ),
                (Err(error), Some(None)) => assert!(
                    error.contains("a version 5 database does not look up"),
                    "{name} {ordinal}: {error}"
                ),
                (got, want) => {
                    panic!("{name} {ordinal} {prefetch_off}: {got:?}, want {want:?}")
                }
            }
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn a_stored_ordinal_id_of_another_number_is_not_fetched() {
    // tests/db_search/ids/make_ids.py (line DBID N2a): sequences whose only Seq-id is
    // `gnl|BL_ORD_ID|n`. With n the sequence's number in its volume it reads as a sequence
    // without parsed Seq-ids (in `pord_vx`, sequence 6 is the second of `pord_x`); another n is
    // a parsed Seq-id, n plus the first ordinal id of the volume, by which NCBI's formatters
    // find no sequence or another one.
    let cwd = empty_cwd("stored_ordinal_id");
    // `Ok(true)`: read as a sequence without parsed Seq-ids; `Ok(false)`: a parsed subject
    // that the reports fetch; `Err(id)`: refused, with its Seq-align id.
    for (name, ordinal, want) in [
        ("ids/pord_v5", 0, Ok(false)),
        ("ids/pord_v5", 2, Ok(true)),
        ("ids/pord_v5", 3, Err("gnl|BL_ORD_ID|6")),
        ("ids/pord_v5", 4, Err("gnl|BL_ORD_ID|1")),
        ("ids/pord_v4", 2, Ok(true)),
        ("ids/pord_v4", 3, Err("gnl|BL_ORD_ID|6")),
        ("ids/pord_vx", 3, Err("gnl|BL_ORD_ID|6")),
        ("ids/pord_vx", 5, Ok(false)),
        ("ids/pord_vx", 6, Ok(true)),
        ("ids/pord_vx", 7, Err("gnl|BL_ORD_ID|12")),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db);
        match (read_header(&headers, ordinal, 0, "BLASTP").unwrap(), want) {
            (ReportedHeader::Stored(_), Ok(true)) => {}
            (ReportedHeader::Parsed(subject), Ok(false)) => {
                assert!(
                    subject
                        .fetched(false, true, crate::report::db_subject::HitLoad::Loaded)
                        .is_ok(),
                    "{name} {ordinal}"
                )
            }
            (ReportedHeader::Parsed(subject), Err(id)) => {
                let error = subject
                    .fetched(false, true, crate::report::db_subject::HitLoad::Loaded)
                    .unwrap_err();
                assert!(
                    error.contains(&format!("the Seq-align id {id},")),
                    "{name} {ordinal}: {error}"
                );
            }
            (got, want) => panic!("{name} {ordinal}: {got:?}, want {want:?}"),
        }
    }
    // The first sequence of `nord_v5` holds `gnl|BL_ORD_ID|3`: the database has parsed Seq-ids.
    let search = SearchPath::new(&cwd, None, None).unwrap();
    for name in ["nord_v5", "nord_v4"] {
        let path = root().join("db_search/ids").join(name);
        let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
        assert!(db.has_parsed_seqids().unwrap(), "{name}");
        assert!(matches!(
            db.read_header(1).unwrap(),
            crate::seqdb::SequenceHeader::Ordinal(_)
        ));
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn an_empty_stored_title_is_made_from_the_name_of_its_taxid() {
    use crate::report::defline::SourceTaxname;
    // tests/db_search/ids/make_ids.py (line DBID H3): `pbare_t` and `nbare_t` hold one
    // sequence without parsed Seq-ids, without a title, with taxid 9606, which the good
    // taxonomy database names, `empty_fields` with an empty scientific name, `no_btd` and
    // `bad_magic` not at all, and whose lookup in `btd_empty` reads past the mapped end of
    // the file (NCBI BLAST+ crashes); tests/fixtures/outfmt16/db/other_tax9606 has titles
    // with that taxid.
    let cwd = empty_cwd("stored_taxonomy");
    let homo = || SourceTaxname::Named(std::sync::Arc::from(&b"Homo sapiens"[..]));
    for (name, seq_type, taxdb, want) in [
        (
            "db_search/ids/pbare_t",
            SeqType::Protein,
            Some("good"),
            homo(),
        ),
        (
            "db_search/ids/pbare_t",
            SeqType::Protein,
            Some("no_btd"),
            SourceTaxname::Unnamed,
        ),
        (
            "db_search/ids/pbare_t",
            SeqType::Protein,
            None,
            SourceTaxname::Unnamed,
        ),
        (
            "db_search/ids/nbare_t",
            SeqType::Nucleotide,
            Some("good"),
            homo(),
        ),
        (
            "db_search/ids/nbare_t",
            SeqType::Nucleotide,
            Some("bad_magic"),
            SourceTaxname::Unnamed,
        ),
        (
            "db_search/ids/nbare_t",
            SeqType::Nucleotide,
            Some("empty_fields"),
            SourceTaxname::Unnamed,
        ),
        (
            "db_search/ids/nbare_t",
            SeqType::Nucleotide,
            Some("btd_empty"),
            SourceTaxname::Refused(""),
        ),
        (
            "fixtures/outfmt16/db/other_tax9606",
            SeqType::Protein,
            Some("good"),
            SourceTaxname::Unnamed,
        ),
    ] {
        let taxdb = taxdb.map(variant);
        let search = SearchPath::new(&cwd, taxdb.as_deref().map(Path::as_os_str), None).unwrap();
        let path = root().join(name);
        let db = SeqDb::open_in(path.to_str().unwrap(), seq_type, &search).unwrap();
        let headers = HeaderReader::new(&db);
        let ReportedHeader::Stored(stored) = read_header(&headers, 0, 0, "BLASTP").unwrap() else {
            panic!("{name}: a sequence without parsed Seq-ids");
        };
        let refused = matches!(want, SourceTaxname::Refused(_));
        if refused {
            assert!(
                matches!(stored.source, SourceTaxname::Refused(_)),
                "{name} {taxdb:?}"
            );
        } else {
            assert_eq!(stored.source, want, "{name} {taxdb:?}");
        }
        assert_eq!(
            stored
                .check_shown_source(|| crate::seqdb::Oid::of_unrestricted(0), "BLASTP")
                .is_err(),
            refused,
            "{name} {taxdb:?}"
        );
        // TBLASTX's check before the search (`report_stream`) for outfmt 0 and for a report
        // that shows no title (outfmt 5 of a nucleotide database): only a source that LOSAT
        // refuses is refused.
        if seq_type == SeqType::Nucleotide {
            let pool = crate::utils::threading::SearchPool::x_serial();
            let shown = |pairwise| ParsedShown {
                pairwise,
                ..Default::default()
            };
            let found = db.taxdb_found();
            let refusal = |pairwise| {
                reports_may_refuse(
                    &db,
                    shown(pairwise),
                    false,
                    found,
                    false,
                    None,
                    &pool,
                    "TBLASTX",
                )
            };
            assert_eq!(
                refusal(true),
                refused.then_some(ReportRefusal::Headers),
                "{taxdb:?}"
            );
            assert_eq!(refusal(false), None, "{taxdb:?}");
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn scientific_names_never_hold() {
    // Names with bytes >= 0x80 (`high_bytes_utf8`: 9606 has one; `high_bytes_latin1`: 10090's
    // common name and kingdom): NCBI's streams write them (outfmt 16 throws at its own
    // level, `report::serial::xml`), so nothing is refused.
    for name in ["sprot_tax", "sprot_mixtax"] {
        for taxonomy in ["good", "high_bytes_utf8", "high_bytes_latin1"] {
            let db = database(name, Some(&variant(taxonomy)));
            assert!(db.taxdb_found());
            assert!(!taxdb_may_refuse(&db), "{name} {taxonomy}");
        }
    }
}

#[test]
fn a_record_past_the_last_page_holds_where_a_written_name_comes_from_past_it() {
    // The record of 10090 (entry 3) ends at the offset of entry 4: 5000 bytes, past the
    // 4096-byte page that ends NCBI's mapping of the small `taxdb.btd`. From 105 its three
    // tabs are inside the mapping (its scientific name too): nothing holds. From 4096 the
    // record lies wholly past the mapping: its names, which XML2 and JSON write, hold (the
    // check does not tell the formats apart: a run whose loads write no name holds too).
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_prepass_btd_{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::copy(variant("good").join("taxdb.btd"), dir.join("taxdb.btd")).unwrap();
    let mut index = std::fs::read(variant("good").join("taxdb.bti")).unwrap();
    let next = 24 + 8 * 4 + 4;
    index[next..next + 4].copy_from_slice(&5000u32.to_be_bytes());
    for (own, holds) in [(105u32, false), (4096, true)] {
        let at = 24 + 8 * 3 + 4;
        index[at..at + 4].copy_from_slice(&own.to_be_bytes());
        std::fs::write(dir.join("taxdb.bti"), &index).unwrap();
        assert_eq!(
            taxdb_may_refuse(&database("sprot_mixtax", Some(&dir))),
            holds,
            "{own}"
        );
        // 9606 (entry 2) is fine: the database with that taxid alone does not hold.
        assert!(!taxdb_may_refuse(&database("sprot_tax", Some(&dir))));
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_missing_or_damaged_taxdb_that_ncbi_survives_does_not_hold() {
    // No `taxdb.bti` on the search path: the formats write no name.
    let db = database("sprot_mixtax", None);
    assert!(!db.taxdb_found());
    assert!(!taxdb_may_refuse(&db));
    // NCBI turns a wrong magic number into a message and N/A (`expected.tsv`).
    let db = database("sprot_mixtax", Some(&variant("bad_magic")));
    assert!(!taxdb_may_refuse(&db));
}

#[test]
fn the_two_sides_agree_and_the_taxonomy_side_alone_clears_a_good_taxdb() {
    use crate::seqdb::taxdb::LazyTaxDb;
    for (taxonomy, expected) in [
        ("good", Some(false)),
        ("high_bytes_latin1", Some(false)),
        ("high_bytes_utf8", Some(false)),
        // The record of 9606 runs past the end of taxdb.btd but not past its mapped page, and
        // 10090's ends before it begins: no lookup fails (`TaxDb::any_lookup_may_fail`; the
        // index's offsets do not increase, which the pass's former stricter copy answered
        // with a hold, review 2 of K4h L1).
        ("offset_next_past_end", Some(false)),
        // An empty taxdb.btd maps no page: a lookup reads past it.
        ("btd_empty", Some(true)),
        ("bad_magic", Some(false)),
        // Not checked: whatever the pass says, the answer is the database side's.
        ("unsorted_adjacent", None),
    ] {
        for name in ["sprot", "sprot_tax", "sprot_mixtax"] {
            let db = database(name, Some(&variant(taxonomy)));
            let taxdb = LazyTaxDb::new(db.search_path().unwrap().clone());
            let side = taxdb_side(&taxdb);
            if let Some(side) = side {
                // The taxonomy side is a superset of the database side.
                assert!(side || !database_side(&db, &taxdb), "{taxonomy} {name}");
                if !side {
                    assert!(!taxdb_may_refuse(&db), "{taxonomy} {name}");
                }
            }
            // Only a record past NCBI's mapping of taxdb.btd may refuse.
            if name == "sprot_mixtax" && expected.is_some() {
                assert_eq!(side, expected, "{taxonomy}");
            }
            // The database side alone and the dispatch give one answer.
            let alone = LazyTaxDb::new(db.search_path().unwrap().clone());
            assert_eq!(
                taxdb_may_refuse(&db),
                database_side(&db, &alone),
                "{taxonomy} {name}"
            );
        }
    }
}

/// Times the taxonomy side on a real NCBI `taxdb` (`LOSAT_TEST_NCBI_TAXDB=<dir with taxdb.bti
/// and taxdb.btd>`): `cargo test --release`-like timing of one pass, and whether any record
/// of it could be refused. Ignored: it needs the 190 MB download of OF5b's fix.
#[test]
#[ignore]
fn real_ncbi_taxdb_walk() {
    use crate::seqdb::taxdb::LazyTaxDb;
    let Some(dir) = std::env::var_os("LOSAT_TEST_NCBI_TAXDB") else {
        return;
    };
    let cwd = crate::utils::test_dir::temp_dir();
    let search = SearchPath::new(&cwd, Some(&dir), None).unwrap();
    let taxdb = LazyTaxDb::new(search);
    let entries = taxdb.entry_count().unwrap();
    for _ in 0..3 {
        let start = std::time::Instant::now();
        let answer = taxdb_side(&taxdb);
        eprintln!(
            "real taxdb: {entries} entries, taxonomy side {answer:?} in {:?}",
            start.elapsed()
        );
        assert_eq!(answer, Some(false));
    }
}

/// `TaxdbRun::bioseq_texts` of a sequential run (dustmasker's and segmasker's reader read each
/// sequence once, blastdbcmd's FASTA formatter twice) on `tests/db_search/t5b/pmv`
/// (make_t5b.py): the OIDs 4, 6, 7 (99990 then 150) and 9 (the second volume) look 99990 up,
/// 5 and 7 take its slot with 150. The counts and places are NCBI 2.17.0's (probe of line RG
/// T5b): the maskers 3 offset errors (OID 7's 99990 hits the slot that 6 filled), blastdbcmd 4
/// (the second read of 7 misses after its own 150).
#[test]
fn a_sequential_run_follows_the_volume_cache_at_every_read() {
    use crate::seqdb::taxdb::{MSG_OFFSET_ERROR, MSG_WRONG_MAGIC};
    let cwd = empty_cwd("sequential");
    let reads = |taxdb: Option<&Path>, times: usize| {
        let db = db_search("t5b/pmv", taxdb, &cwd);
        let run = TaxdbRun::sequential();
        let mut texts = Vec::new();
        for oid in 0..db.num_seqs() {
            let set = db.defline_set(oid).unwrap();
            for _ in 0..times {
                for text in run.bioseq_texts(&db, oid, &set, None)? {
                    texts.push((oid, text));
                }
            }
        }
        Ok::<_, SeqDbError>((texts, run.texts.get().map(|texts| texts.offset_error_taxid)))
    };
    let at = |oids: &[u32], text: &str| -> Vec<(u32, String)> {
        oids.iter().map(|&oid| (oid, text.to_string())).collect()
    };
    let offset = variant("offset_last_past_end");
    assert_eq!(
        reads(Some(&offset), 1).unwrap(),
        (at(&[4, 6, 9], MSG_OFFSET_ERROR), Some(Some(99990)))
    );
    assert_eq!(
        reads(Some(&offset), 2).unwrap(),
        (at(&[4, 6, 7, 9], MSG_OFFSET_ERROR), Some(Some(99990)))
    );
    // The texts of opening at the first lookup (OID 2), then nothing to look up.
    for (name, text) in [
        ("bad_magic", MSG_WRONG_MAGIC.to_string()),
        (
            "count_zero",
            "SeqDB: Taxid metadata indicates (0) entries but file has room for (7).".to_string(),
        ),
    ] {
        assert_eq!(
            reads(Some(&variant(name)), 2).unwrap(),
            (at(&[2], &text), Some(None))
        );
    }
    // A good taxonomy database: no text; none on the path: nothing opened.
    assert_eq!(
        reads(Some(&variant("good")), 2).unwrap(),
        (Vec::new(), Some(None))
    );
    assert_eq!(reads(None, 2).unwrap(), (Vec::new(), None));
    // An empty taxdb.btd (NCBI maps it as a null pointer): the first lookup that reads a
    // record, of 9606 (OID 2), faults.
    assert!(matches!(
        reads(Some(&variant("btd_empty")), 1),
        Err(SeqDbError::Unsupported { .. })
    ));
    // A database without the offset error whose record of 99990 (OIDs 4, 6, 7, 9) reaches
    // past NCBI's mapping of taxdb.btd (to the next entry's offset, 5000): NCBI reads memory
    // past the mapping, and its loads write nothing of it, so they go on as NCBI's where it
    // completes; the names of 99990 have no defined result (no tab inside the mapping).
    let built = cwd.join("taxdb_99990_past_page");
    std::fs::create_dir_all(&built).unwrap();
    let mut index = Vec::new();
    for word in [
        0x8739u32, 4, 0, 0, 0, 0, 9606, 0, 99990, 20, 99991, 5000, 99992, 30,
    ] {
        index.extend(word.to_be_bytes());
    }
    let mut data = b"Homo sapiens\th\tp\tE".to_vec();
    data.resize(100, b'x');
    std::fs::write(built.join("taxdb.bti"), index).unwrap();
    std::fs::write(built.join("taxdb.btd"), data).unwrap();
    let db = db_search("t5b/pmv", Some(&built), &cwd);
    let run = TaxdbRun::sequential();
    for oid in 0..4 {
        let set = db.defline_set(oid).unwrap();
        assert_eq!(
            run.bioseq_texts(&db, oid, &set, None).unwrap(),
            Vec::<String>::new()
        );
    }
    assert_eq!(
        run.texts.get().map(|texts| texts.offset_error_taxid),
        Some(None)
    );
    for oid in 4..db.num_seqs() {
        let set = db.defline_set(oid).unwrap();
        assert_eq!(
            run.bioseq_texts(&db, oid, &set, None).unwrap(),
            Vec::<String>::new()
        );
    }
    assert!(matches!(
        run.lookup(&db, 99990, &mut Vec::new()),
        Ok(crate::seqdb::taxdb::TaxLookup::PastMapping(_))
    ));
    assert!(matches!(
        run.lookup(&db, 9606, &mut Vec::new()),
        Ok(crate::seqdb::taxdb::TaxLookup::Named(_))
    ));
    std::fs::remove_dir_all(&cwd).unwrap();
}

/// `PrefetchTexts::of`, the one producer of what a run's lookups can write, for the search's
/// prefetch and the readers' first load alike: `tests/db_search/t5b/taxdb/dup_last_btd_empty`
/// (make_t5b.py), an empty taxdb.btd whose index holds its last taxid 99990 twice, so that the
/// lookup of 99990 ends on the entry before the last one and reads a record (NCBI's fault,
/// no offset error). Both runs refuse at their first lookup that reads a record, where NCBI
/// 2.17.0 faults (SIGSEGV in all six search programs, probe of line RG T5b review 2).
#[test]
fn the_search_and_the_readers_refuse_a_null_mapping_whose_last_lookup_faults() {
    let cwd = empty_cwd("null_map_last_faults");
    let taxdb = root().join("db_search/t5b/taxdb/dup_last_btd_empty");
    let db = db_search("t5b/pmv", Some(&taxdb), &cwd);
    let want = PrefetchTexts {
        offset_error_taxid: None,
        maps_null: true,
        opening: false,
    };
    let search = TaxdbRun::default();
    // The database stores parsed Seq-ids: each subject with its parsed definition lines.
    let parsed = parsed_subjects(&db);
    let all = (0..db.num_seqs()).map(|ordinal| (ordinal, Some(&parsed[ordinal as usize])));
    let refusal = prefetch_messages(&db, &search, all, "BLASTP").unwrap_err();
    assert!(
        refusal
            .to_string()
            .contains("reaches past the mapped end of the file (0 bytes)"),
        "{refusal}"
    );
    assert_eq!(search.texts.get(), Some(&want));
    let readers = TaxdbRun::sequential();
    let mut first_refused = None;
    for oid in 0..db.num_seqs() {
        let set = db.defline_set(oid).unwrap();
        if let Err(error) = readers.bioseq_texts(&db, oid, &set, None) {
            assert!(matches!(error, SeqDbError::Unsupported { .. }), "{error}");
            first_refused = Some(oid);
            break;
        }
    }
    // The OIDs 0 and 1 have no taxid; OID 2 is 9606, whose record (68..105) faults.
    assert_eq!(first_refused, Some(2));
    assert_eq!(readers.texts.get(), Some(&want));
    std::fs::remove_dir_all(&cwd).unwrap();
}

/// `prefetch_messages`: the open texts at the run's first lookup, a volume's lookup once per
/// taxid over the batches of a run where no volume holds another taxid of the slot of NCBI's
/// cache that the offset error's taxid takes, nothing without a taxonomy database or a taxid.
#[test]
fn the_prefetch_writes_each_text_once_without_a_cache_collision() {
    let text = |db: &SeqDb, run: &TaxdbRun, ordinals: &[u32]| -> String {
        let subjects = ordinals.iter().map(|&ordinal| (ordinal, None));
        String::from_utf8(prefetch_messages(db, run, subjects, "BLASTP").unwrap()).unwrap()
    };
    let db = database("sprot_mixtax", Some(&variant("bad_magic")));
    let run = TaxdbRun::default();
    assert_eq!(
        text(&db, &run, &[0, 1, 2, 0]),
        "Error: [blastp] Error: Tax database file has wrong magic number.\n"
    );
    assert_eq!(text(&db, &run, &[2, 1]), "");
    // No taxonomy database on the search path, or one that writes nothing.
    for taxdb in [None, Some(variant("good"))] {
        let db = database("sprot_mixtax", taxdb.as_deref());
        assert_eq!(text(&db, &TaxdbRun::default(), &[0, 1, 2]), "");
    }
    // A database without taxids does not open the taxonomy database.
    let db = database("sprot", Some(&variant("bad_magic")));
    assert_eq!(text(&db, &TaxdbRun::default(), &[0, 1]), "");
    let cwd = empty_cwd("prefetch");
    // The last entry's offset error (taxid 99990 of x16_taxid_99990), once per volume; the
    // two volumes of tm_2vol (copies of x16_taxid_99990) write it once each.
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    for (name, volumes) in [("x16_taxid_99990", 1), ("tm_2vol", 2)] {
        let db = db_search(name, Some(&variant("offset_last_past_end")), &cwd);
        let run = TaxdbRun::default();
        let all: Vec<u32> = (0..db.num_seqs()).collect();
        assert_eq!(text(&db, &run, &all), offset.repeat(volumes), "{name}");
        assert_eq!(text(&db, &run, &all), "", "{name}");
    }
    // tm_cache151 (parsed Seq-ids): 99990 (st1, st3) and 151 (st2), which takes another slot
    // than 99990's (99990 & 255 = 150): once, whatever the order.
    let db = db_search("tm_cache151", Some(&variant("offset_last_past_end")), &cwd);
    let parsed = parsed_subjects(&db);
    for order in [[0, 1, 2], [2, 1, 0], [1, 0, 2]] {
        let run = TaxdbRun::default();
        for expected in [offset, ""] {
            let subjects = order
                .iter()
                .map(|&ordinal| (ordinal, Some(&parsed[ordinal as usize])));
            let text = prefetch_messages(&db, &run, subjects, "BLASTP").unwrap();
            assert_eq!(String::from_utf8(text).unwrap(), expected, "{order:?}");
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

/// `prefetch_messages`: a run whose volume holds the offset error's taxid (99990, slot 150 of
/// NCBI's 256) and another taxid of that slot (150) is refused at its first prefetch, before
/// its subjects are looked at, and at every later one; the taxonomy databases that write no
/// offset error, and the same database with 151, are not.
#[test]
fn a_cache_collision_with_the_offset_error_is_refused_at_the_first_prefetch() {
    // The databases have parsed Seq-ids: the subjects with their definition lines.
    let prefetch = |db: &SeqDb, run: &TaxdbRun, ordinals: &[u32]| {
        let parsed = parsed_subjects(db);
        let subjects = ordinals
            .iter()
            .map(|&ordinal| (ordinal, Some(&parsed[ordinal as usize])))
            .collect::<Vec<_>>();
        prefetch_messages(db, run, subjects, "BLASTP")
    };
    let cwd = empty_cwd("cache_collision");
    let db = db_search("tm_cache150", Some(&variant("offset_last_past_end")), &cwd);
    let run = TaxdbRun::default();
    assert!(run.prefetch_refused(&db));
    for ordinals in [&[][..], &[0], &[1, 2]] {
        let error = prefetch(&db, &run, ordinals).unwrap_err();
        assert!(error.is::<TaxonomyRefusal>(), "{error:#}");
        let text = format!("{error:#}");
        assert!(text.contains("taxid 99990"), "{text}");
        assert!(text.contains("taxid 150,"), "{text}");
        assert!(text.contains("256-slot taxonomy cache"), "{text}");
        assert!(text.contains("LOSAT's BLASTP"), "{text}");
    }
    // Nothing was looked up: the run's own handle is not open.
    assert!(run.lazy.get().is_none());
    for taxdb in ["good", "bad_magic", "high_bytes_utf8"] {
        let db = db_search("tm_cache150", Some(&variant(taxdb)), &cwd);
        let run = TaxdbRun::default();
        assert!(!run.prefetch_refused(&db), "{taxdb}");
        assert!(prefetch(&db, &run, &[0, 1, 2]).is_ok(), "{taxdb}");
    }
    let db = db_search("tm_cache150", None, &cwd);
    assert!(!TaxdbRun::default().prefetch_refused(&db));
    let db = db_search("tm_cache151", Some(&variant("offset_last_past_end")), &cwd);
    assert!(!TaxdbRun::default().prefetch_refused(&db));
    // The two volumes of tm_2vol hold 99990 alone.
    let db = db_search("tm_2vol", Some(&variant("offset_last_past_end")), &cwd);
    assert!(!TaxdbRun::default().prefetch_refused(&db));
    std::fs::remove_dir_all(&cwd).unwrap();
}

/// `read_reported_taxonomy` (XML2 and JSON, line DBID L1d): NCBI fetches every reported
/// subject, so a run whose volume has the cache collision is refused where a reported subject
/// has a positive taxid (`nz_rep`: a1 99990, b1 150, c1 99990, d1 without a taxid, e1 -5;
/// `nz_repo`: a1 to d1 without parsed Seq-ids), and not for d1 or e1 alone or for nothing
/// reported; the taxonomy databases without the offset error read the same subjects.
#[test]
fn the_xml2_lookups_ask_for_the_cache_refusal_where_a_subject_has_a_taxid() {
    let cwd = empty_cwd("xml2_cache");
    let open = |taxdb: &str| {
        let search = SearchPath::new(&cwd, Some(variant(taxdb).as_os_str()), None).unwrap();
        let path = root().join("db_search/nz_rep");
        SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap()
    };
    let read = |db: &SeqDb, run: &TaxdbRun, ranks: &[usize]| {
        let parsed: Vec<_> = parsed_subjects(db)
            .into_iter()
            .map(|subject| Some(std::sync::Arc::new(subject)))
            .collect();
        let ordinals: Vec<u32> = (0..db.num_seqs()).collect();
        let reported: std::collections::BTreeSet<usize> = ranks.iter().copied().collect();
        read_reported_taxonomy(db, run, &ordinals, &reported, &parsed, "BLASTN").map(|_| ())
    };
    let db = open("offset_last_past_end");
    for ranks in [&[0][..], &[1], &[3, 2], &[0, 1, 2, 3]] {
        let run = TaxdbRun::default();
        let error = read(&db, &run, ranks).unwrap_err();
        assert!(error.is::<TaxonomyRefusal>(), "{ranks:?} {error:#}");
        let text = format!("{error:#}");
        assert!(text.contains("256-slot taxonomy cache"), "{text}");
        assert!(text.contains("LOSAT's BLASTN"), "{text}");
        // Refused before any lookup: the run's own handle is not open.
        assert!(run.lazy.get().is_none());
    }
    // Nothing reported, d1 (no taxid), e1 (taxid -5: a volume looks a positive taxid up only).
    for ranks in [&[][..], &[3], &[4], &[3, 4]] {
        assert!(read(&db, &TaxdbRun::default(), ranks).is_ok(), "{ranks:?}");
    }
    for taxdb in ["good", "bad_magic"] {
        assert!(
            read(&open(taxdb), &TaxdbRun::default(), &[0, 1, 2, 3, 4]).is_ok(),
            "{taxdb}"
        );
    }
    // nz_repo: a1 to d1 without parsed Seq-ids (the taxid of the definition line).
    let taxdb = variant("offset_last_past_end");
    let search = SearchPath::new(&cwd, Some(taxdb.as_os_str()), None).unwrap();
    let path = root().join("db_search/nz_repo");
    let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
    let ordinals: Vec<u32> = (0..db.num_seqs()).collect();
    for (ranks, refused) in [(&[0][..], true), (&[1, 3], true), (&[3], false)] {
        let reported: std::collections::BTreeSet<usize> = ranks.iter().copied().collect();
        let run = TaxdbRun::default();
        let no_parsed: &[Option<&crate::report::db_subject::ParsedSubject>] = &[];
        let result = read_reported_taxonomy(&db, &run, &ordinals, &reported, no_parsed, "BLASTN");
        assert_eq!(result.is_err(), refused, "{ranks:?}");
        if let Err(error) = result {
            assert!(error.is::<TaxonomyRefusal>(), "{ranks:?} {error:#}");
        }
    }
    std::fs::remove_dir_all(&cwd).unwrap();
}

/// `VolumeLookups::offset_slot`: the slot of the offset error's taxid L = 99990 (slot 150) of
/// NCBI's cache of each volume, in the order of the loads (review 4 of line DBID pssmread):
/// the first load of a sequence writes L's text at each lookup of L that misses the slot (twice
/// for L, 150, L in one load), a read again of a sequence whose load looks L up writes nothing
/// where the slot holds L whether NCBI's scope keeps it or loads it again, and is refused where
/// loading it again would write the text; a load after a read again that may have put 150
/// there is refused; other slots and other volumes are apart.
#[test]
fn the_slot_of_the_offset_errors_taxid_follows_the_loads() {
    const L: i32 = 99990;
    let fetch = |ordinal, taxid, again| SlotFetch {
        ordinal,
        taxid,
        again,
    };
    // H1 (cc_m): one load of L, 150, L; again: its own 150 is between.
    let mut record = VolumeLookups::default();
    assert_eq!(
        record.offset_slot(0, 0, &[L, 150, L], L),
        Ok(vec![true, false, true])
    );
    assert_eq!(
        record.offset_slot(0, 0, &[L, 150, L], L),
        Err(SlotUncertain {
            ordinal: 0,
            again: true,
            after: Some(fetch(0, 150, true)),
        })
    );
    // M1 (b) (cc_m2): A (L, 150) and B (L), first loads: two texts. Read again, A may put 150
    // there (Either), and B's read again is not certain.
    let mut record = VolumeLookups::default();
    assert_eq!(
        record.offset_slot(0, 0, &[L, 150], L),
        Ok(vec![true, false])
    );
    assert_eq!(record.offset_slot(0, 1, &[L], L), Ok(vec![true]));
    assert_eq!(
        record.offset_slot(0, 0, &[L, 150], L),
        Ok(vec![false, false])
    );
    assert_eq!(
        record.offset_slot(0, 1, &[L], L),
        Err(SlotUncertain {
            ordinal: 1,
            again: true,
            after: Some(fetch(0, 150, true)),
        })
    );
    // A first load after such a read again (C of 99990 + 256 * 2 is not L).
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![true]));
    assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
    assert_eq!(record.offset_slot(0, 2, &[L], L), Ok(vec![true]));
    assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 3, &[9606, L + 512, L], L),
        Ok(vec![false, false, true])
    );
    assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 4, &[L], L),
        Err(SlotUncertain {
            ordinal: 4,
            again: false,
            after: Some(fetch(1, 150, true)),
        })
    );
    // M1 (a) (cc_ax): 150 never looked up: L's reads again write nothing.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![true]));
    for _ in 0..3 {
        assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![false]));
    }
    // cc_ac: A (L), C (150), then A again: refused, after C's fetch.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &[L], L), Ok(vec![true]));
    assert_eq!(record.offset_slot(0, 1, &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 0, &[L], L),
        Err(SlotUncertain {
            ordinal: 0,
            again: true,
            after: Some(fetch(1, 150, false)),
        })
    );
    // C, then A in another volume (cc_c1 cc_a1), then both again; taxids of other slots.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(1, 0, &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 1, &[L, 9606], L),
        Ok(vec![true, false])
    );
    assert_eq!(record.offset_slot(1, 0, &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 1, &[L, 9606], L),
        Ok(vec![false, false])
    );
    assert_eq!(
        record.offset_slot(0, 5, &[9606, 151], L),
        Ok(vec![false, false])
    );
    assert!(!record.slot_loads.contains(&5));
}

/// `pssm_load_texts` on `tm_cache150` (st1 99990, st2 150, st3 99990, parsed Seq-ids, one
/// volume) with the offset error of `offset_last_past_end`: the PSSM's loads in their order
/// write NCBI's texts, and a read again that NCBI's scope decides is refused, naming the fetch
/// of 150; with 151 (`tm_cache151`) the read again writes nothing.
#[test]
fn the_pssms_loads_write_the_offset_error_as_the_slot_holds_it() {
    let cwd = empty_cwd("pssm_slot");
    let offset = "Error: [psiblast] Error: Offset error at end of taxdb file.\n";
    let loads = |name: &str, ordinals: &[u32]| -> Result<String> {
        let db = db_search(name, Some(&variant("offset_last_past_end")), &cwd);
        let parsed = parsed_subjects(&db);
        let run = TaxdbRun::default();
        let mut text = Vec::new();
        for &ordinal in ordinals {
            let subject = &parsed[ordinal as usize];
            let lookups = blob_taxids(&db, ordinal, Some(subject), ordinal, "PSIBLAST")?;
            text.extend(pssm_load_texts(
                &db, &run, None, ordinal, &lookups, "PSIBLAST",
            )?);
        }
        Ok(String::from_utf8(text).unwrap())
    };
    assert_eq!(loads("tm_cache150", &[0, 1, 2]).unwrap(), offset.repeat(2));
    assert_eq!(loads("tm_cache150", &[0, 2, 1]).unwrap(), offset);
    assert_eq!(loads("tm_cache150", &[0, 0, 2]).unwrap(), offset);
    let error = loads("tm_cache150", &[0, 1, 0]).unwrap_err();
    assert!(error.is::<TaxonomyRefusal>(), "{error:#}");
    let text = format!("{error:#}");
    for part in [
        "taxid 99990, the last entry of the taxonomy database",
        "reads database sequence 0 of the volume",
        "again after the run's fetch of database sequence 1 with taxid 150",
        "256-slot taxonomy cache",
        "LOSAT's PSIBLAST",
    ] {
        assert!(text.contains(part), "{part}: {text}");
    }
    assert_eq!(loads("tm_cache151", &[0, 1, 0, 2]).unwrap(), offset);
    std::fs::remove_dir_all(&cwd).unwrap();
}

#[test]
fn a_database_without_taxids_does_not_hold() {
    let db = database("sprot", Some(&variant("high_bytes_utf8")));
    assert!(db.taxdb_found());
    assert!(!taxdb_may_refuse(&db));
}

/// NCBI's `CATCH_ALL` takes the `CEofException` of a record that ends before its read enters
/// the set's elements as a `CIOException` other than `eFlush`: no post, exit status 0 (RG-T5d
/// fix 1, probe le8); the `CSerialException` of a record cut inside is `Error: what()`, 255.
#[test]
fn a_header_record_cut_before_its_elements_ends_the_search_with_exit_status_0() {
    let path = std::path::Path::new("t.phr");
    let native = |record: &[u8]| {
        let error = crate::seqdb::header::decode_full(record, path).unwrap_err();
        let error = database_error(error, "BLASTP");
        let native = error.downcast_ref::<NativeError>().expect("NCBI's error");
        (native.exit, native.message.is_empty())
    };
    assert_eq!(native(&[0x30]), (0, true));
    assert_eq!(native(&[0x3f]), (0, true));
    assert_eq!(native(&[0x30, 0x80]), (255, false));
    assert_eq!(native(&[0x30, 0x80, 0x30, 0x80, 0xbf]), (255, false));
}
