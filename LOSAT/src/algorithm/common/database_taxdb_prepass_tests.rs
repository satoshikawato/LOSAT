//! `taxdb_may_refuse`: the answer of the lookup of every sequence's taxid, with the
//! databases of `tests/fixtures/outfmt16/db` (`sprot_tax`: two sequences, taxid 9606;
//! `sprot_mixtax`: those two with taxid 10090 and a third with taxid 9606) and the taxonomy
//! databases of `tests/seqdb/taxdb`.
use super::*;
use crate::seqdb::alias::SearchPath;
use crate::seqdb::SeqType;
use crate::utils::test_dir::TempDir;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
}

/// A database and the empty current directory it was opened from (its search path keeps it for
/// the files looked for later), which is removed with it.
struct Opened {
    db: SeqDb,
    _cwd: TempDir,
}

impl std::ops::Deref for Opened {
    type Target = SeqDb;

    fn deref(&self) -> &SeqDb {
        &self.db
    }
}

/// The database `name` of `tests/fixtures/outfmt16/db` opened with `BLASTDB=taxdb`
/// (`None`: not set), from an empty current directory.
fn database(name: &str, taxdb: Option<&Path>) -> Opened {
    let cwd = TempDir::unique("prepass-cwd");
    let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None);
    let path = root().join("fixtures/outfmt16/db").join(name);
    let db = SeqDb::open_in(path.to_str().unwrap().as_bytes(), SeqType::Protein, &search).unwrap();
    Opened { db, _cwd: cwd }
}

fn variant(name: &str) -> PathBuf {
    root().join("seqdb/taxdb").join(name)
}

/// An empty current directory of the test `tag`, of its own for each call, which is removed
/// with what it holds when it is dropped, at the end of the test, passed or failed.
fn empty_cwd(tag: &str) -> TempDir {
    TempDir::unique(&format!("{tag}-cwd"))
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
    let search = SearchPath::new(cwd, taxdb.map(Path::as_os_str), None);
    let path = root().join("db_search").join(name);
    SeqDb::open_in(path.to_str().unwrap().as_bytes(), SeqType::Protein, &search).unwrap()
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
        let headers = HeaderReader::new(&db, false);
        let ReportedHeader::Parsed(subject) = read_header(&headers, 1, 0, "BLASTP").unwrap() else {
            panic!("{name}: a sequence with parsed Seq-ids");
        };
        match (
            subject.fetched(false, crate::report::db_subject::HitLoad::Loaded),
            want,
        ) {
            (Ok(()), None) => {}
            (Err(error), Some(text)) => assert!(error.contains(text), "{name}: {error}"),
            (got, _) => panic!("{name}: {got:?}"),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1447-1455 (quoted
// at `SeqidLookup::seqid_to_oids`): without `multi`, the lookup returns at the first OID.
/// Line DBID v4late (review 3 of line DBID loadord, H3): in a list of version 4 volumes whose
/// later volume's string ISAM index does not open (mode 000), NCBI's lookup of an id that an
/// earlier volume holds returns there (NCBI BLAST+ 2.17.0: exit status 0, probe
/// `blastdb-ids-20261021/R/v4late`), and so does LOSAT's lookup (`HeaderReader::resolution`,
/// which lists no other holder: NCBI's places hold only what the run loaded, line DBID hr-1):
/// what the lookup gives (the subject, the same record, another one) is the answer. An id that
/// only the later volume can hold, or any id where that volume comes first, opens it: NCBI's
/// stop. What the prefetch reads (review 1 of the item, H1 and H2): what the lookup gives, whose
/// definition line `GetBioseq` may reject (NCBI then stops).
#[cfg(unix)]
#[test]
fn a_later_volume_whose_isam_does_not_open_is_not_read_where_the_lookup_returns_before_it() {
    use crate::report::db_subject::FoundInstead;
    use std::os::unix::fs::PermissionsExt;
    // A directory of its own, made empty: a file of mode 000 left by a run before cannot be
    // copied over.
    let cwd = empty_cwd("v4late");
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
        return;
    }
    let open = |list: &str, taxdb: Option<&str>| {
        let names = list
            .split(' ')
            .map(|stem| cwd.join(stem).to_str().unwrap().to_owned())
            .collect::<Vec<_>>()
            .join(" ");
        let taxdb = taxdb.map(|variant| root().join("seqdb/taxdb").join(variant));
        let search = SearchPath::new(&cwd, taxdb.as_deref().map(Path::as_os_str), None);
        SeqDb::open_in(names.as_bytes(), SeqType::Protein, &search).unwrap()
    };
    let subject_of = |db: &SeqDb, ordinal: u32| {
        let headers = HeaderReader::new(db, false);
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
        match (subject.found(), instead) {
            (None, None) => {}
            (Some(FoundInstead::Sequence(got)), Some(want)) => {
                assert_eq!(got.ordinal(), want, "{list} {ordinal}")
            }
            (got, _) => panic!("{list} {ordinal}: {got:?}"),
        }
        let fetched = subject.fetched(false, crate::report::db_subject::HitLoad::Loaded);
        match instead {
            None => assert!(fetched.is_ok(), "{list} {ordinal}: {fetched:?}"),
            Some(_) => assert!(
                fetched
                    .unwrap_err()
                    .contains("lookup gives database sequence 0 first"),
                "{list} {ordinal}"
            ),
        }
        assert!(subject.pssm_read().is_ok(), "{list} {ordinal}");
    }
    // l5b, which the lookup of l5e's `ref|NP_000002.1|` gives first, holds it after
    // `ref|NP_000999.1|` in one definition line: NCBI's prefetch reads it and stops.
    let db = open("l5b l5e", None);
    let subject = subject_of(&db, 1);
    assert!(matches!(
        subject.found(),
        Some(FoundInstead::Sequence(other)) if other.ordinal() == 0
    ));
    let thrown = prefetch_reads(&db, 1, &subject, false, false, "BLASTP")
        .unwrap()
        .unwrap_err();
    assert!(
        thrown.contains("oid headers do not contain target gi/seq_id"),
        "{thrown}"
    );
    // The subject's own blob (taxid 9606), whatever the taxonomy database: no other holder is
    // listed, so the prefetch reads what the lookup gives.
    for taxdb in ["good", "bad_magic"] {
        let db = open("l5c l5e", Some(taxdb));
        for ordinal in [0, 1] {
            let blob = prefetch_reads(
                &db,
                ordinal,
                &subject_of(&db, ordinal),
                true,
                false,
                "BLASTP",
            )
            .unwrap()
            .unwrap_or_else(|reason| panic!("{taxdb} {ordinal}: {reason}"));
            assert_eq!(blob, Some(ordinal), "{taxdb} {ordinal}");
        }
    }
}

#[test]
fn the_same_record_with_other_taxids_is_shown_apart_by_the_reports_of_taxids_only() {
    // tests/db_search/ids/make_ids.py (line DBID N1b): `pdup_at` holds `gb|EEE00005.1|`
    // without a taxid then with taxid 9606, `pdup_tu` the other way, `pdup_no` twice with
    // taxid 99990 (the last entry of every variant of tests/seqdb/taxdb, whose lookup writes
    // the offset error of `offset_last_past_end`), `pdupe_at` an empty title without a taxid
    // then with 9606; `taxdb`, the variant on BLASTDB. `without`, `with`: the refusal without
    // the reports of taxids (`fetched(false)`) and with them. Sequence 0, which the lookup of
    // the shared id gives first, is refused for nothing: every load of that id gives its blob,
    // also after sequence 1 (line DBID hr-1, `ScopeWalk`; the static scan of the other holders
    // that refused it is gone), so the refusals are sequence 1's.
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
        ("ids/pdup_at", None, 0, (None, None)),
        ("ids/pdup_at", None, 1, (None, Some(taxids(0)))),
        ("ids/pdup_tu", None, 0, (None, None)),
        ("ids/pdup_tu", None, 1, (None, Some(taxids(0)))),
        ("ids/pdup_at", Some("good"), 1, (None, Some(taxids(0)))),
        // Texts of opening at the run's first lookup, which only the copy with a taxid makes.
        ("ids/pdup_at", Some("bad_magic"), 0, (None, None)),
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
        ("ids/pdup_no", Some("offset_last_past_end"), 0, (None, None)),
        // The same record without taxids looks nothing up.
        ("ids/pdup_ab", Some("offset_last_past_end"), 1, (None, None)),
        // An empty title made from the taxids where the taxonomy database is found, which the
        // reports of generated titles refuse (below).
        ("ids/pdupe_at", None, 0, (None, None)),
        ("ids/pdupe_at", Some("good"), 0, (None, None)),
        ("ids/pdupe_at", Some("good"), 1, (None, Some(taxids(0)))),
        // A taxonomy database that names no taxid: NCBI makes no title from them.
        ("ids/pdupe_at", Some("no_btd"), 1, (None, Some(taxids(0)))),
        (
            "ids/pdupe_at",
            Some("empty_fields"),
            1,
            (None, Some(taxids(0))),
        ),
        // One that writes a text at the run's first lookup, which only the copy with a taxid
        // makes.
        ("ids/pdupe_at", Some("bad_magic"), 0, (None, None)),
        // Another title is another record whatever the taxids.
        ("ids/pdup_ac", None, 1, both(first(0))),
    ] {
        let taxdb = taxdb.map(|variant| root().join("seqdb/taxdb").join(variant));
        let db = db_search(name, taxdb.as_deref(), &cwd);
        let headers = HeaderReader::new(&db, false);
        let ReportedHeader::Parsed(subject) = read_header(&headers, ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
        };
        for (taxonomy, want) in [(false, &without), (true, &with)] {
            match (
                subject.fetched(taxonomy, crate::report::db_subject::HitLoad::Loaded),
                want,
            ) {
                (Ok(()), None) => {}
                (Err(error), Some(text)) => {
                    assert!(error.contains(text.as_str()), "{name} {ordinal}: {error}")
                }
                (got, _) => panic!("{name} {ordinal} {taxdb:?} {taxonomy}: {got:?}"),
            }
        }
        // Read after the other copy, the subject's load gives what its lookup gives.
        let ReportedHeader::Parsed(other) =
            read_header(&headers, 1 - ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {}: a sequence with parsed Seq-ids", 1 - ordinal);
        };
        let mut walk = crate::report::db_subject::ScopeWalk::new(
            crate::report::db_subject::LoadsBefore::none(),
            Some(false),
        );
        walk.load(&other);
        assert_eq!(
            walk.load(&subject),
            crate::report::db_subject::HitLoad::Loaded,
            "{name} {ordinal} {taxdb:?}"
        );
        // Line DBID N1b L1: outfmt 0 and the tabular word of a replaced id
        // (`ParsedSubject::generated_title_fetched`).
        let shown = name == "ids/pdupe_at"
            && ordinal == 1
            && taxdb.as_deref().is_some_and(|dir| dir.ends_with("good"));
        match subject.generated_title_fetched() {
            Ok(()) => assert!(!shown, "{name} {ordinal} {taxdb:?}"),
            Err(error) => assert!(
                shown && error.contains(&title(1 - ordinal)),
                "{name} {ordinal} {taxdb:?}: {error}"
            ),
        }
    }
}

#[test]
fn a_seq_align_id_of_any_kind_is_left_to_the_lookup() {
    // tests/db_search/ids/make_look_fasta.py (line DBID N2b): `None`, NCBI's formatters find the
    // sequence itself; else the text of the refusal.
    let cwd = empty_cwd("any_kind");
    let other = |n: u32| format!("lookup gives database sequence {n} first");
    let none = "formatters do not find the sequence".to_string();
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
        ("ids/pshare_v5", 1, None),
        (
            "ids/pshare_v5",
            2,
            Some("whose lookup in the database throws on the id itself".to_string()),
        ),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db, false);
        let ReportedHeader::Parsed(subject) = read_header(&headers, ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
        };
        assert!(!subject.unfetchable(), "{name} {ordinal}");
        match (
            subject.fetched(false, crate::report::db_subject::HitLoad::Loaded),
            want,
        ) {
            (Ok(()), None) => {}
            (Err(error), Some(text)) => {
                assert!(error.contains(&text), "{name} {ordinal}: {error}")
            }
            (got, _) => panic!("{name} {ordinal}: {got:?}"),
        }
    }
    // Another sequence holds the Seq-align id (review 1 of N2b; ids that differ in case alone,
    // review 2): NCBI's places give it where it was loaded before the subject's load (line DBID
    // hr-1, `ScopeWalk`), and the reports refuse the subject with its words; loaded first, the
    // subject is its own. A third sequence of the id of two copies of a record (`pdup_abc`,
    // another title) loads, by its own id, the first copy that the lookup gives: no load of the
    // copies gives it (`None`).
    use crate::report::db_subject::{LoadsBefore, ScopeWalk};
    let mut differ = Vec::new();
    for (name, ordinal, other, holder) in [
        ("ids/pshare_v5", 0, 1, Some(1)),
        ("ids/pshare_v4", 0, 1, Some(1)),
        ("ids/pshare_ab", 0, 1, Some(1)),
        ("ids/pshare_v5", 4, 5, Some(5)),
        ("ids/pshare_v4", 4, 5, Some(5)),
        ("ids/pshare_cd", 0, 2, Some(2)),
        ("ids/pshare_cd", 1, 3, Some(3)),
        ("ids/pdup_abc", 0, 2, None),
        ("ids/pdup_abc", 1, 2, None),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db, false);
        let parsed = |ordinal: u32| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
            ReportedHeader::Parsed(subject) => subject,
            _ => panic!("{name} {ordinal}: a sequence with parsed Seq-ids"),
        };
        let (subject, before) = (parsed(ordinal), parsed(other));
        let first = subject.fetched(false, crate::report::db_subject::HitLoad::Loaded);
        if first.is_err() {
            differ.push(format!("{name} {ordinal} read first: {first:?}"));
        }
        let mut walk = ScopeWalk::new(LoadsBefore::none(), Some(false));
        walk.load(&before);
        let load = walk.load(&subject);
        match (load, holder) {
            (crate::report::db_subject::HitLoad::Shown(shown), Some(holder))
                if shown.holder.map(|holder| holder.ordinal()) == Some(holder) =>
            {
                let error = subject.fetched(false, load).unwrap_err();
                let words = format!("which database sequence {holder} also holds");
                if !error.contains(&words) {
                    differ.push(format!("{name} {ordinal} after {other}: {error}"));
                }
            }
            (crate::report::db_subject::HitLoad::Loaded, None) => {}
            _ => differ.push(format!("{name} {ordinal} after {other}: {load:?}")),
        }
    }
    assert!(differ.is_empty(), "{}", differ.join("\n"));
}

#[test]
fn the_prefetch_reads_a_subject_by_its_seq_align_id_as_ncbis_loader_does() {
    // Line DBID P3 (probes of `blastdb-ids-20261021/R/intmaxpf`): where no report reads the
    // Bioseq, `prefetch_reads` without texts (no `taxdb.bti` on this search path).
    let cwd = empty_cwd("prefetch_reads");
    let stops = "and then stops (oid headers do not contain target gi/seq_id)";
    for (name, ordinal, want) in [
        // `pdb|1ABC|A` of pshare_c, whose lookup gives it first, and `lcl|555`, which sequence
        // 1 also holds in a second definition line: read as itself where nothing loaded before
        // gives another (below, the walk where one did).
        ("ids/pshare_cd", 0, Ok(Some(0))),
        ("ids/pshare_v5", 0, Ok(Some(0))),
        ("ids/pshare_v5", 1, Ok(Some(1))),
        // `gnl|BL_ORD_ID|0` stored with a text tag.
        (
            "ids/pshare_v5",
            2,
            Err("whose lookup in the database throws on the id itself"),
        ),
        // An id that the lookup finds nowhere, and one that it finds first on sequence 5,
        // none of whose ids compares equal: nothing is read.
        ("ids/plook_v5", 4, Ok(None)),
        ("ids/plook_v5", 6, Ok(None)),
    ] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db, false);
        let ReportedHeader::Parsed(subject) = read_header(&headers, ordinal, 0, "BLASTP").unwrap()
        else {
            panic!("{name} {ordinal}: a sequence with parsed Seq-ids");
        };
        match (
            prefetch_reads(&db, ordinal, &subject, false, false, "BLASTP").unwrap(),
            want,
        ) {
            (Ok(got), Ok(want)) => assert_eq!(got, want, "{name} {ordinal}"),
            (Err(reason), Err(text)) => {
                assert!(reason.contains(text), "{name} {ordinal}: {reason}")
            }
            (got, _) => panic!("{name} {ordinal}: {got:?}"),
        }
        // A lookup that fails is the batch's warning (`batch_prefetch`), not a refusal.
        let check = check_prefetched_parsed(
            &db,
            &subject,
            crate::report::db_subject::HitLoad::Loaded,
            ordinal,
            false,
            "BLASTP",
        );
        assert_eq!(
            check.is_err(),
            want.is_err() && !subject.prefetch_fails(),
            "{name} {ordinal}"
        );
        assert_eq!(
            subject.prefetch_fails(),
            want == Err("whose lookup in the database throws on the id itself"),
            "{name} {ordinal}"
        );
    }
    // `pdb|1ABC|A` of pshare_c and `pdb|1abc|A` of pshare_d, prefetched together: the bulk load
    // gives one the other's blob, whose definition lines do not hold its id, and NCBI stops
    // there (exit 2; `ScopeWalk::bulk`).
    {
        use crate::report::db_subject::{HitLoad, LoadsBefore, ScopeWalk};
        let db = db_search("ids/pshare_cd", None, &cwd);
        let headers = HeaderReader::new(&db, false);
        let parsed: Vec<_> = [0, 2]
            .map(
                |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                    ReportedHeader::Parsed(subject) => (ordinal, subject),
                    _ => panic!("pshare_cd {ordinal}: a sequence with parsed Seq-ids"),
                },
            )
            .into();
        let loads = ScopeWalk::new(LoadsBefore::none(), Some(false)).bulk(
            parsed
                .iter()
                .map(|(ordinal, subject)| (*ordinal as usize, &**subject)),
        );
        let shown: Vec<_> = parsed
            .iter()
            .filter(|(ordinal, _)| matches!(loads[&(*ordinal as usize)], HitLoad::Shown(_)))
            .collect();
        assert_eq!(shown.len(), 1, "{loads:?}");
        let (ordinal, subject) = shown[0];
        let load = loads[&(*ordinal as usize)];
        let error = check_prefetched_parsed(&db, subject, load, *ordinal, false, "BLASTP")
            .expect_err(&format!("pshare_cd {ordinal}: {load:?}"));
        assert!(format!("{error:#}").contains(stops), "{error:#}");
    }
    // `lcl|777` before `gi|777`: the formatters do not find the Bioseq, which the prefetch
    // reads as itself; it fails where it asks for the length (a minus-strand alignment).
    let search = SearchPath::new(&cwd, None, None);
    let path = root().join("db_search/ids/nrefuse_v5");
    let db = SeqDb::open_in(
        path.to_str().unwrap().as_bytes(),
        SeqType::Nucleotide,
        &search,
    )
    .unwrap();
    let headers = HeaderReader::new(&db, false);
    let ReportedHeader::Parsed(subject) = read_header(&headers, 0, 0, "TBLASTX").unwrap() else {
        panic!("nrefuse_v5 0: a sequence with parsed Seq-ids");
    };
    assert!(subject.unfetchable());
    assert_eq!(
        prefetch_reads(&db, 0, &subject, true, false, "TBLASTX").unwrap(),
        Ok(Some(0))
    );
    let run = TaxdbRun::default();
    let warning = b"Warning: [tblastx] Error pre-fetching sequence data\n".to_vec();
    let prefetch = |flipped: bool| {
        let subject = PrefetchSubject {
            ordinal: 0,
            parsed: Some(&*subject),
            flipped,
            load: None,
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
}

#[test]
fn a_batch_whose_prefetch_fails_writes_ncbis_warning_after_its_texts() {
    // tests/db_search/ids/make_ids.py (line DBID pfwarn): `ntag_t5` 0 `gnl|BL_ORD_ID|0`
    // stored with a text tag (its lookup fails), 1 `gb|FFF00092.1|`; `nlook_t4` 1 `bbs|123`,
    // which a version 4 database does not look up (no Bioseq), 0 `gb|AF000070.1|`; every
    // sequence with taxid 99990, whose lookup in `offset_last_past_end` writes the offset
    // error.
    let cwd = TempDir::unique("batch-prefetch");
    let taxdb = root().join("seqdb/taxdb/offset_last_past_end");
    let open = |name: &str, taxdb: Option<&Path>| {
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None);
        let path = root().join("db_search/ids").join(name);
        SeqDb::open_in(
            path.to_str().unwrap().as_bytes(),
            SeqType::Nucleotide,
            &search,
        )
        .unwrap()
    };
    let offset = "Error: [blastn] Error: Offset error at end of taxdb file.\n";
    let warning = "Warning: [blastn] Error pre-fetching sequence data\n";
    let text = |posts: crate::report::diag_log::Posts| String::from_utf8(posts.text()).unwrap();
    for taxdb in [None, Some(taxdb.as_path())] {
        let db = open("ntag_t5", taxdb);
        let headers = HeaderReader::new(&db, false);
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
                load: None,
            }],
            "BLASTN",
            &stop_queries(0),
        );
        assert_eq!(text(alone.unwrap()), warning);
        // With the other subject, whose Seq-id handle (packed 92) NCBI's map puts before the
        // failed lookup's (a text tag, not packed): its load writes the offset error before the
        // warning (line NP item O4b2-10; LOSAT's refusal before it).
        let both = batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [
                PrefetchSubject {
                    ordinal: 0,
                    parsed: Some(&*failed),
                    flipped: false,
                    load: None,
                },
                PrefetchSubject {
                    ordinal: 1,
                    parsed: Some(&*other),
                    flipped: false,
                    load: None,
                },
            ],
            "BLASTN",
            &stop_queries(0),
        );
        match taxdb {
            None => assert_eq!(text(both.unwrap()), warning),
            Some(_) => assert_eq!(text(both.unwrap()), format!("{offset}{warning}")),
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
                    load: None,
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
                        load: None,
                    },
                    PrefetchSubject {
                        ordinal: 0,
                        parsed: Some(&*failed),
                        flipped: false,
                        load: None,
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
    let headers = HeaderReader::new(&db, false);
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
                load: None,
            },
            PrefetchSubject {
                ordinal: 0,
                parsed: Some(&*control),
                flipped: false,
                load: None,
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
            load: None,
        }],
        "BLASTN",
        &stop_queries(0),
    );
    assert_eq!(text(plus.unwrap()), "");
}

#[test]
fn the_loads_before_a_stopping_lookup_follow_ncbis_order_of_seq_id_handles() {
    // Line DBID loadord (tests/db_search/make_multidb.py --loadord-only): `ldp` holds
    // `ref|NP_000002.1|` (packed 2) and `ref|NP_000900.1|`; the unparsed `qp1` is
    // `gnl|BL_ORD_ID|0` beside it in `lda ldp`, `|1` in `ldb ldp` and `|2` in `ldp ldu`.
    let cwd = TempDir::unique("load-order");
    let open = |names: &[&str]| {
        let search = SearchPath::new(&cwd, None, None);
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
        SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Protein, &search).unwrap()
    };
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db, false);
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
    loaded.note_loaded(0, || [unparsed, held]);
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
    never.note_loaded(0, || [unparsed]);
    assert!(never.loaded().is_empty());
    assert_eq!(
        unsure(&db, &never, &[stopping], stop_queries(1), stop, stop_parsed),
        None
    );
    // A call of the library's batches after others, whose reports nothing recorded: its
    // reader lists the other holders from the database (`HeaderReader::unrecorded_holders`).
    let mut after_calls = TaxdbRun::default();
    after_calls.note_calls_before(true);
    let ReportedHeader::Parsed(unrecorded_stop) =
        read_header(&HeaderReader::new(&db, true), stop, 0, "BLASTP").unwrap()
    else {
        panic!("the stopping subject has parsed Seq-ids");
    };
    let stopping_after = PrefetchSubject::unflipped(stop, Some(&unrecorded_stop));
    let why = unsure(
        &db,
        &after_calls,
        &[stopping_after],
        stop_queries(1),
        stop,
        &unrecorded_stop,
    )
    .unwrap();
    assert!(
        why.contains("whose reports LOSAT has no record of"),
        "{why}"
    );
    // A copy of the stopping sequence with the same Seq-align id is the same lookup: no
    // batch before can have loaded it without stopping there (`ldc`, a copy of `ldp`'s
    // record; kit case `copy.21b` of line DBID mixdb2's review 1, the rows `.copy.6`), so its
    // record registers nothing.
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
    copied.note_loaded(0, || {
        [PrefetchSubject::unflipped(
            copy,
            parsed[copy as usize].as_deref(),
        )]
    });
    assert!(copied.loaded().is_empty());
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
    // Line DBID hr-1: the run's record is made in every run (the walks of the reports' loads
    // read it, not only a prefetch's stop), each record's subjects made once; whether a
    // prefetch of the run can stop is still asked of the database (BLASTX's scopes).
    let closed = open(&["ldp"]);
    let run = TaxdbRun::default();
    let mut made = 0;
    run.note_loaded(0, || {
        made += 1;
        std::iter::empty::<PrefetchSubject<'_>>()
    });
    assert!(!run.prefetch_can_stop(&closed) && made == 1);
    assert!(run.loaded().is_empty());
    let open_gate = open(&["lda", "ldp"]);
    assert!(TaxdbRun::default().prefetch_can_stop(&open_gate));
    // Review 4, L3: `ldd` and `ldh` hold one Seq-align id, gi|9000, and the stopping
    // accession: a load of the later one (recorded) by that id loads the first (NCBI's
    // `x_GetOid` takes the lookup's first OID), so the record registers `ldd`, and the
    // sequence that NCBI may give for the stopping id instead is `ldd`, never `ldh` itself.
    let db = open(&["lda", "ldd", "ldh", "ldp"]);
    let parsed = headers_of(&db);
    let stop = accession(&parsed, "NP_000002.1");
    assert_eq!(db.oid(stop).get(), 3);
    let recorded = TaxdbRun::default();
    recorded.note_loaded(0, || [PrefetchSubject::unflipped(2, parsed[2].as_deref())]);
    assert_eq!(recorded.loaded(), [1]);
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
    let cwd = TempDir::unique("two-stops");
    let open = |names: &[&str], taxdb: Option<&Path>| {
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None);
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
        SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Protein, &search).unwrap()
    };
    type Headers = [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>];
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db, false);
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
        let failure = result
            .unwrap_err()
            .downcast::<ErrorAfterBatchReports>()
            .unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (
            failure.query.expect("the stop names its query"),
            String::from_utf8(native.message).unwrap(),
        )
    };
    let refused = |result: Result<crate::report::diag_log::Posts>| -> String {
        let refusal = result.unwrap_err();
        assert!(!refusal.is::<ErrorAfterBatchReports>(), "{refusal}");
        refusal.to_string()
    };
    let error = "LOSAT database error: DB contains no accession info.\n".to_owned();
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
            run.note_loaded(0, || [subject(&parsed, holder)]);
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
            run.note_loaded(0, || [subject(&parsed, holder)]);
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
    let cwd = TempDir::unique("e2-placed");
    let bad_magic = root().join("seqdb/taxdb/bad_magic");
    let open = |names: &[&str]| {
        let search = SearchPath::new(&cwd, Some(bad_magic.as_os_str()), None);
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
        SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Protein, &search).unwrap()
    };
    type Headers = [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>];
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db, false);
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
        let failure = result
            .unwrap_err()
            .downcast::<ErrorAfterBatchReports>()
            .unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (
            failure.query.expect("the stop names its query"),
            String::from_utf8(native.message).unwrap(),
        )
    };
    let error = "LOSAT database error: DB contains no accession info.\n".to_owned();
    // The model of the cut report of query 0 in `run`: its reads of `read` (each sequence's
    // first load in the run; the first writes the texts of opening).
    let cut_report = |db: &SeqDb, run: &TaxdbRun, read: &[u32]| {
        let reads = subject_read_texts(
            db,
            run,
            read.iter().map(|&ordinal| (ordinal, ordinal, None)),
            RowLookups::None,
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
    between
        .loads(subject_read_texts(&db, &run, [(1, 1, None)], RowLookups::None, "BLASTP").unwrap());
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
    assert!(!refusal.is::<ErrorAfterBatchReports>(), "{refusal}");
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
        RowLookups::None,
        "BLASTP",
    )
    .unwrap()];
    run.note_loaded(0, || [subject(&parsed, holder)]);
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
    let cwd = TempDir::unique("record");
    let bad_magic = root().join("seqdb/taxdb/bad_magic");
    let search = SearchPath::new(&cwd, Some(bad_magic.as_os_str()), None);
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
    let db = SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Protein, &search).unwrap();
    let headers = HeaderReader::new(&db, false);
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
            format: PrefetchFormat::of(6, crate::report::archive::options::RequestKind::Blastp),
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
    let error = "LOSAT database error: DB contains no accession info.\n".to_owned();
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
        let failure = error.downcast::<ErrorAfterBatchReports>().unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (
            failure.query.expect("the stop names its query"),
            String::from_utf8(native.message).unwrap(),
        )
    };
    let refused = |error: anyhow::Error| {
        assert!(!error.is::<ErrorAfterBatchReports>(), "{error}");
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
        RowLookups::None,
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
    let cwd = TempDir::unique("prefetch-stop");
    let open = |names: &[&str], taxdb: Option<&Path>| {
        let search = SearchPath::new(&cwd, taxdb.map(Path::as_os_str), None);
        let list: Vec<String> = names
            .iter()
            .map(|name| root().join(name).to_str().unwrap().to_owned())
            .collect();
        SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Nucleotide, &search).unwrap()
    };
    let headers_of =
        |db: &SeqDb| -> Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>> {
            let headers = HeaderReader::new(db, false);
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
            load: None,
        }
    }
    let stopped = |result: Result<()>| -> (usize, String) {
        let failure = result
            .unwrap_err()
            .downcast::<ErrorAfterBatchReports>()
            .unwrap();
        let native = failure.error.downcast::<NativeError>().unwrap();
        assert_eq!(native.exit, 2);
        (
            failure.query.expect("the stop names its query"),
            String::from_utf8(native.message).unwrap(),
        )
    };
    let error = "LOSAT database error: DB contains no accession info.\n".to_owned();
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
            }
            Some(_) => {
                // In a run that has looked nothing up: the texts of the load before the stop,
                // then the error.
                let stop = batch_prefetch(
                    &db,
                    &TaxdbRun::default(),
                    std::iter::empty(),
                    both(),
                    "BLASTN",
                    &stop_queries(0),
                )
                .map(drop);
                assert_eq!(stopped(stop), (0, format!("{open}{error}")));
                // A batch after one whose prefetch opened it and looked the other subject up
                // (every format that prefetches, the XML2 and JSON ones too, line DBID taxpf-1):
                // its load writes nothing now, so NCBI's error.
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
            }
        }
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
                crate::ncbi_exception::NcbiException::new(
                    crate::seqdb::isam::OBJECT_ID_GET_ID,
                    "Invalid choice selection",
                ),
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
        batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            [subject(other, &throws), subject(stop as u32, &parsed[stop])],
            "BLASTN",
            &stop_queries(0),
        )
        .map(drop),
    ] {
        let refused = refused.unwrap_err();
        assert!(!refused.is::<ErrorAfterBatchReports>());
        assert!(
            refused.to_string().contains("whose lookup throws"),
            "{refused}"
        );
    }
}

#[test]
fn a_weak_match_is_shown_where_ncbi_reads_the_holder_one_by_one_before() {
    use crate::report::db_subject::{HitLoad, LoadsBefore, ParsedSubject, ScopeWalk, ShownBy};
    // tests/db_search/ids/make_ids.py (line DBID N2c): `pweak_v5` and `pweak_v4` 0
    // `sp|P00041.1|NAMEW`, 1 `gb|FFF00040|`, 2 `pir|PIRW1|`, 3 `sp||NAMEX`, 4 `gb|FFF00041.1|`,
    // 5 `sp||NAMEY`, 6 `sp|P00042.1|NAMEY`, 7 `gb|FFF00043|`, 8 `gb|FFF00043.1|`, 9
    // `sp|P00044.1|NAMEZ`, 10 `sp||NAMEZ`. With the prefetch on (the bulk load, no scope's
    // index), NCBI writes each sequence's own report; where a report reads one by one a
    // sequence that the subject's Seq-align id matches weakly before it, in the batch's scope,
    // the scope gives that sequence (line DBID hr-1: the walk answers it among the loaded
    // sequences, in both versions; `ScopeWalk`).
    let cwd = empty_cwd("weak_match");
    // (subject, the sequence that its Seq-align id matches weakly)
    let weak = [(5, 6), (6, 5), (7, 8), (9, 10), (10, 9)];
    for name in ["ids/pweak_v5", "ids/pweak_v4"] {
        let db = db_search(name, None, &cwd);
        let headers = HeaderReader::new(&db, false);
        let subjects: Vec<std::sync::Arc<ParsedSubject>> = (0..db.num_seqs())
            .map(
                |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                    ReportedHeader::Parsed(subject) => subject,
                    _ => panic!("{name} {ordinal}: no parsed Seq-ids"),
                },
            )
            .collect();
        let after = |before: &[u32], subject: u32| {
            let mut walk = ScopeWalk::new(LoadsBefore::none(), Some(false));
            walk.enter(Some(0));
            for &other in before {
                walk.load(&subjects[other as usize]);
            }
            walk.load(&subjects[subject as usize])
        };
        for subject in 0..db.num_seqs() {
            let holder = weak
                .iter()
                .find(|&&(of, _)| of == subject)
                .map(|&(_, holder)| holder);
            // The subject read first, and the prefetch's bulk load of every subject.
            let first = after(&[], subject);
            assert!(
                !matches!(first, HitLoad::Shown(_)),
                "{name} {subject}: {first:?}"
            );
            let bulk = ScopeWalk::new(LoadsBefore::none(), Some(false)).bulk(
                subjects
                    .iter()
                    .enumerate()
                    .map(|(key, subject)| (key, &**subject)),
            );
            assert!(
                !matches!(bulk.get(&(subject as usize)), Some(HitLoad::Shown(_))),
                "{name} {subject}: {:?}",
                bulk.get(&(subject as usize))
            );
            // Every other sequence read one by one before it.
            let others: Vec<u32> = (0..db.num_seqs())
                .filter(|&other| other != subject)
                .collect();
            match (after(&others, subject), holder) {
                (HitLoad::Shown(shown), Some(holder)) => {
                    assert_eq!(
                        shown.holder.map(|holder| holder.ordinal()),
                        Some(holder),
                        "{name} {subject}"
                    );
                    assert_eq!(shown.by, ShownBy::Weak, "{name} {subject}");
                    let error = subjects[subject as usize]
                        .fetched(false, HitLoad::Shown(shown))
                        .unwrap_err();
                    assert!(
                        error
                            .contains(&format!("match an id of database sequence {holder} weakly")),
                        "{name} {subject}: {error}"
                    );
                }
                (load, None) => assert!(
                    !matches!(load, HitLoad::Shown(_)),
                    "{name} {subject}: {load:?}"
                ),
                (load, Some(holder)) => panic!("{name} {subject}: {load:?}, want {holder}"),
            }
        }
    }
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
        let headers = HeaderReader::new(&db, false);
        match (read_header(&headers, ordinal, 0, "BLASTP").unwrap(), want) {
            (ReportedHeader::Stored(_), Ok(true)) => {}
            (ReportedHeader::Parsed(subject), Ok(false)) => {
                assert!(
                    subject
                        .fetched(false, crate::report::db_subject::HitLoad::Loaded)
                        .is_ok(),
                    "{name} {ordinal}"
                )
            }
            (ReportedHeader::Parsed(subject), Err(id)) => {
                let error = subject
                    .fetched(false, crate::report::db_subject::HitLoad::Loaded)
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
    let search = SearchPath::new(&cwd, None, None);
    for name in ["nord_v5", "nord_v4"] {
        let path = root().join("db_search/ids").join(name);
        let db = SeqDb::open_in(
            path.to_str().unwrap().as_bytes(),
            SeqType::Nucleotide,
            &search,
        )
        .unwrap();
        assert!(db.has_parsed_seqids().unwrap(), "{name}");
        assert!(matches!(
            db.read_header(1).unwrap(),
            crate::seqdb::SequenceHeader::Ordinal(_)
        ));
    }
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
        let search = SearchPath::new(&cwd, taxdb.as_deref().map(Path::as_os_str), None);
        let path = root().join(name);
        let db = SeqDb::open_in(path.to_str().unwrap().as_bytes(), seq_type, &search).unwrap();
        let headers = HeaderReader::new(&db, false);
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
                    (false, false),
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
    let dir = TempDir::unique("prepass-btd");
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
    let search = SearchPath::new(&cwd, Some(&dir), None);
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
}

/// `prefetch_messages` (the reads of a report, in their order): a run whose volume holds the offset
/// error's taxid (99990, slot 150 of NCBI's 256) and another taxid of that slot (150) follows the
/// slot through the reads (line NP item O4b2-5; refused at its first read until then): st1 writes
/// the offset error, st2 (150) takes the slot, st3 asks again; the taxonomy databases that write no
/// offset error, and the same database with 151, have no collision.
#[test]
fn a_cache_collision_with_the_offset_error_follows_the_slot_through_the_reads() {
    // The databases have parsed Seq-ids: the subjects with their definition lines.
    let prefetch = |db: &SeqDb, run: &TaxdbRun, ordinals: &[u32]| {
        let parsed = parsed_subjects(db);
        let subjects = ordinals
            .iter()
            .map(|&ordinal| (ordinal, Some(&parsed[ordinal as usize])))
            .collect::<Vec<_>>();
        prefetch_messages(db, run, subjects, "BLASTP")
    };
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    let cwd = empty_cwd("cache_collision");
    let db = db_search("tm_cache150", Some(&variant("offset_last_past_end")), &cwd);
    let run = TaxdbRun::default();
    assert!(run.cache_may_refuse(&db));
    for (ordinals, expected) in [(&[][..], ""), (&[0], offset), (&[1, 2], offset)] {
        let text = String::from_utf8(prefetch(&db, &run, ordinals).unwrap()).unwrap();
        assert_eq!(text, expected, "{ordinals:?}");
    }
    // Read again: NCBI keeps the three sequences, so nothing loads.
    assert_eq!(prefetch(&db, &run, &[0, 1, 2]).unwrap(), b"");
    for taxdb in ["good", "bad_magic", "high_bytes_utf8"] {
        let db = db_search("tm_cache150", Some(&variant(taxdb)), &cwd);
        let run = TaxdbRun::default();
        assert!(!run.cache_may_refuse(&db), "{taxdb}");
        assert!(prefetch(&db, &run, &[0, 1, 2]).is_ok(), "{taxdb}");
    }
    let db = db_search("tm_cache150", None, &cwd);
    assert!(!TaxdbRun::default().cache_may_refuse(&db));
    let db = db_search("tm_cache151", Some(&variant("offset_last_past_end")), &cwd);
    assert!(!TaxdbRun::default().cache_may_refuse(&db));
    // The two volumes of tm_2vol hold 99990 alone.
    let db = db_search("tm_2vol", Some(&variant("offset_last_past_end")), &cwd);
    assert!(!TaxdbRun::default().cache_may_refuse(&db));
}

/// `unordered_boundary_misses` against every order of up to seven loads of the four kinds, from a
/// slot that holds the offset error's taxid or another one (line NP item O4b2-3).
#[test]
fn unordered_misses_are_every_order() {
    // The kinds as (first lookup is `L`, last lookup is `L`).
    const KINDS: [(bool, bool); 4] = [(true, true), (true, false), (false, true), (false, false)];
    fn orders(left: &mut [usize; 4], prev_l: bool, misses: usize, out: &mut Vec<usize>) {
        if left.iter().all(|&n| n == 0) {
            out.push(misses);
            return;
        }
        for kind in 0..4 {
            if left[kind] == 0 {
                continue;
            }
            left[kind] -= 1;
            let (first_l, last_l) = KINDS[kind];
            orders(left, last_l, misses + usize::from(first_l && !prev_l), out);
            left[kind] += 1;
        }
    }
    for total in 2..=7usize {
        for ll in 0..=total {
            for lo in 0..=total - ll {
                for ol in 0..=total - ll - lo {
                    let oo = total - ll - lo - ol;
                    for held_last in [false, true] {
                        let mut all = Vec::new();
                        orders(&mut [ll, lo, ol, oo], held_last, 0, &mut all);
                        let want = (*all.iter().min().unwrap(), *all.iter().max().unwrap());
                        let loads = UnorderedLoads { ll, lo, ol, oo };
                        assert_eq!(
                            unordered_boundary_misses(held_last, loads),
                            want,
                            "{loads:?} held_last={held_last}"
                        );
                    }
                }
            }
        }
    }
}

/// A batch's prefetch on a database whose volume holds taxid 99990 (the offset error of
/// `offset_last_past_end`) and 150 (its slot of NCBI's cache) follows the slot in NCBI's order of
/// the loads, and refuses only where the addresses of the Seq-id handles decide the texts (line
/// NP item O4b2-3; the local ids of `nz_pc` and `tm_cache150` are one key of NCBI's map): `nz_pc`
/// sZ 150, sI 562, sM 10090, sN 99990 writes the offset error once in any order; `tm_cache150`
/// st1 99990, st2 150, st3 99990 once or twice; st1 and st2 once, and then the slot holds either.
#[test]
fn a_prefetch_follows_ncbis_order_of_its_loads_in_the_cache_slot() {
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    let cwd = empty_cwd("prefetch_slot");
    let taxdb = variant("offset_last_past_end");
    let prefetch = |db: &SeqDb,
                    run: &TaxdbRun,
                    parsed: &[crate::report::db_subject::ParsedSubject],
                    ordinals: &[u32]| {
        batch_prefetch(
            db,
            run,
            std::iter::empty(),
            ordinals.iter().map(|&ordinal| {
                PrefetchSubject::unflipped(ordinal, Some(&parsed[ordinal as usize]))
            }),
            "BLASTP",
            &stop_queries(0),
        )
        .map(|posts| String::from_utf8(posts.text()).unwrap())
    };
    let db = db_search("nz_pc", Some(&taxdb), &cwd);
    let parsed = parsed_subjects(&db);
    for ordinals in [&[0, 1, 2, 3][..], &[3, 2, 1, 0], &[3, 0], &[0, 3]] {
        let run = TaxdbRun::default();
        assert!(run.cache_may_refuse(&db));
        assert_eq!(
            prefetch(&db, &run, &parsed, ordinals).unwrap(),
            offset,
            "{ordinals:?}"
        );
    }
    // A later batch that loads sN again after sZ: NCBI keeps every sequence of a database of
    // four (`keeps_every_load`), so it loads nothing (the probe of item O4b2-3,
    // `A.nzpc.nz_q.*.limu.b1`: one offset error).
    let run = TaxdbRun::default();
    assert_eq!(prefetch(&db, &run, &parsed, &[0, 3]).unwrap(), offset);
    assert_eq!(prefetch(&db, &run, &parsed, &[0, 3]).unwrap(), "");
    assert_eq!(prefetch(&db, &run, &parsed, &[3, 0, 1]).unwrap(), "");

    let db = db_search("tm_cache150", Some(&taxdb), &cwd);
    let parsed = parsed_subjects(&db);
    let refused = prefetch(&db, &TaxdbRun::default(), &parsed, &[0, 1, 2]).unwrap_err();
    assert!(refused.is::<TaxonomyRefusal>(), "{refused:#}");
    let text = format!("{refused:#}");
    for part in [
        "the prefetch of a query batch reads database sequences",
        "taxids 99990 and 150",
        "same slot of NCBI BLAST+'s 256-slot taxonomy cache",
        "addresses of their Seq-id handles",
        "LOSAT's BLASTP",
    ] {
        assert!(text.contains(part), "{part}: {text}");
    }
    // One batch each: the order of the batches.
    let run = TaxdbRun::default();
    assert_eq!(prefetch(&db, &run, &parsed, &[0]).unwrap(), offset);
    assert_eq!(prefetch(&db, &run, &parsed, &[1]).unwrap(), "");
    assert_eq!(prefetch(&db, &run, &parsed, &[2]).unwrap(), offset);
    // st1 and st2 in one batch: once, and the slot then holds 99990 or 150 by their order. The
    // same batch again loads nothing (NCBI keeps the three sequences).
    let run = TaxdbRun::default();
    assert_eq!(prefetch(&db, &run, &parsed, &[1, 0]).unwrap(), offset);
    assert_eq!(prefetch(&db, &run, &parsed, &[0, 1]).unwrap(), "");
    let later = format!("{:#}", prefetch(&db, &run, &parsed, &[2]).unwrap_err());
    assert!(
        later.contains("after the prefetch of a query batch read database sequences"),
        "{later}"
    );
    // Without the offset error, nothing to follow.
    for name in ["good", "bad_magic"] {
        let db = db_search("tm_cache150", Some(&variant(name)), &cwd);
        assert!(
            prefetch(&db, &TaxdbRun::default(), &parsed, &[0, 1, 2]).is_ok(),
            "{name}"
        );
    }
}

/// `read_reported_taxonomy` (XML2 and JSON) makes no load of a reported subject: those formats
/// prefetch as the others (`prefetches_subjects`), so the batches' prefetch or the reports' reads
/// make the loads before this call, and refuse what NCBI's order of them does not decide
/// (`prefetch_slot_plan`, `load_read_texts`; the merge of line NP item O4b2-6 into line DBID item
/// pf7). On a volume with the cache collision (`nz_rep`: a1 99990, b1 150, c1 99990, d1 without a
/// taxid, e1 -5; `nz_repo`: a1 to d1 without parsed Seq-ids) it reads the `Hit`s' lookups alone,
/// as on the taxonomy databases without the offset error.
#[test]
fn the_xml2_lookups_leave_the_cache_slot_to_the_loads() {
    let cwd = empty_cwd("xml2_cache");
    let open = |taxdb: &str| {
        let search = SearchPath::new(&cwd, Some(variant(taxdb).as_os_str()), None);
        let path = root().join("db_search/nz_rep");
        SeqDb::open_in(
            path.to_str().unwrap().as_bytes(),
            SeqType::Nucleotide,
            &search,
        )
        .unwrap()
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
    for taxdb in ["offset_last_past_end", "good", "bad_magic"] {
        let db = open(taxdb);
        for ranks in [&[][..], &[0], &[1], &[3, 2], &[3], &[4], &[0, 1, 2, 3, 4]] {
            let result = read(&db, &TaxdbRun::default(), ranks);
            assert!(result.is_ok(), "{taxdb} {ranks:?} {result:?}");
        }
    }
    // nz_repo: a1 to d1 without parsed Seq-ids (the taxid of the definition line).
    let taxdb = variant("offset_last_past_end");
    let search = SearchPath::new(&cwd, Some(taxdb.as_os_str()), None);
    let path = root().join("db_search/nz_repo");
    let db = SeqDb::open_in(
        path.to_str().unwrap().as_bytes(),
        SeqType::Nucleotide,
        &search,
    )
    .unwrap();
    let ordinals: Vec<u32> = (0..db.num_seqs()).collect();
    for ranks in [&[0][..], &[1, 3], &[3]] {
        let reported: std::collections::BTreeSet<usize> = ranks.iter().copied().collect();
        let no_parsed: &[Option<&crate::report::db_subject::ParsedSubject>] = &[];
        let result = read_reported_taxonomy(
            &db,
            &TaxdbRun::default(),
            &ordinals,
            &reported,
            no_parsed,
            "BLASTN",
        );
        assert!(result.is_ok(), "{ranks:?}");
    }
}

/// `VolumeLookups::offset_slot`: the slot of the offset error's taxid L = 99990 (slot 150) of
/// NCBI's cache of each volume, in the order of the loads (review 4 of line DBID pssmread):
/// the first load of a sequence writes L's text at each lookup of L that misses the slot (twice
/// for L, 150, L in one load), a read again of a sequence whose load looks L up writes nothing
/// where the slot holds L whether NCBI's scope keeps it or loads it again, and is refused where
/// loading it again would write the text; a load after a read again that may have put 150
/// there is refused; other slots and other volumes are apart.
/// Where NCBI keeps every sequence that the run loads (a database of no more sequences than its
/// data source's cache holds, `keeps_every_load`), a read again loads nothing: no lookup, the slot
/// as it was (the probe of item O4b2-3, `A.nzpc.nz_q.*.limu.b1`: a later batch's prefetch of the
/// same subjects writes no second offset error with four heap layouts).
#[test]
fn a_read_again_loads_nothing_where_ncbi_keeps_every_sequence() {
    const L: i32 = 99990;
    let mut record = VolumeLookups {
        kept: true,
        ..VolumeLookups::default()
    };
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![true, false])
    );
    assert_eq!(record.offset_slot(0, 1, &by(1), &[L], L), Ok(vec![true]));
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![false, false])
    );
    assert_eq!(record.offset_slot(0, 1, &by(1), &[L], L), Ok(vec![false]));
    assert_eq!(record.offset_slot(0, 2, &by(2), &[150], L), Ok(vec![false]));
    assert_eq!(record.offset_slot(0, 3, &by(3), &[L], L), Ok(vec![true]));
    let after_three = (record.offset_slots.clone(), record.slot_blobs.clone());
    // A read by another Seq-align id is another blob of NCBI's data loader (the sequence with
    // that id's handle, review 1 of item O4b2-5, H1), a load of its own where no Bioseq of the
    // sequence's blobs holds that id (item O4b2-7): sequence 2 by id 4 puts 150 there again.
    assert_eq!(record.offset_slot(0, 2, &by(4), &[150], L), Ok(vec![false]));
    // Sequence 3 again by its own id loads nothing; a first load of `L` after it misses.
    assert_eq!(record.offset_slot(0, 3, &by(3), &[L], L), Ok(vec![false]));
    assert_eq!(record.offset_slot(0, 5, &by(5), &[L], L), Ok(vec![true]));
    // Where the Bioseq of sequence 2's blob holds id 4, NCBI's scope may give it without a load:
    // the slot holds either, and the first load of `L` after it is not certain.
    (record.offset_slots, record.slot_blobs) = after_three;
    assert_eq!(
        record.offset_slot(0, 2, &found_in(4, 2), &[150], L),
        Ok(vec![false])
    );
    assert_eq!(record.offset_slot(0, 3, &by(3), &[L], L), Ok(vec![false]));
    assert!(record.offset_slot(0, 5, &by(5), &[L], L).is_err());
    // A read by an id whose handle LOSAT does not know (NCBI makes one at each read) is refused.
    assert_eq!(
        record.offset_slot(0, 6, &BlobName::default(), &[L], L),
        Err(SlotUncertain {
            ordinal: 6,
            again: false,
            after: Some(SlotFetch {
                ordinal: 2,
                taxid: 150,
                again: true,
                maybe: false,
            }),
            order: None,
            unnamed: true,
        })
    );
}

/// The blob name of a read by the Seq-align id `id` ([`BlobName`]): its handle, and a scope key
/// that its Bioseq holds and by which it finds such a Bioseq.
fn by(id: u32) -> BlobName {
    BlobName {
        handle: Some(id.to_le_bytes().to_vec()),
        finds: (vec![u64::from(id)], Vec::new()),
        holds: BioseqKeys {
            scope: vec![u64::from(id)],
            ..BioseqKeys::default()
        },
    }
}

/// [`by`] of `id`, which also finds the Bioseq of a read by `other`.
fn found_in(id: u32, other: u32) -> BlobName {
    let mut name = by(id);
    name.finds.0.push(u64::from(other));
    name
}

/// `VolumeLookups::unordered_slot`: loads of one key whose order LOSAT does not know, a load that
/// NCBI may make or not counted both ways (review 1 of item O4b2-3, F2).
#[test]
fn a_group_of_loads_takes_both_cases_of_a_load_made_or_not() {
    const L: i32 = 99990;
    let load = |ordinal: u32, blob: u32, lookups: &[i32], name: BlobName| SlotLoad {
        ordinal,
        blob,
        volume: 0,
        lookups: lookups.to_vec(),
        order: 7,
        name,
        pins: false,
        maybe: false,
    };
    // Batch 1: sZ (150) alone; batch 2: sZ again and sQ (150): without `kept`, sZ's blob is
    // loaded or not, and no case writes a text; the slot holds 150 in every case.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &by(0), &[150], L), Ok(vec![false]));
    let (z, q) = (load(0, 0, &[150], by(0)), load(1, 1, &[150], by(1)));
    let misses = record.unordered_slot(0, &[&z, &q], L).unwrap();
    assert_eq!(misses, vec![(0, vec![false]), (1, vec![false])]);
    assert_eq!(record.offset_slots[&0].holds, SlotHolds::Other);
    // A (L) and then A again with B (L): A's second load, made or not, comes after L in any
    // order (the slot holds L): no text.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &by(0), &[L], L), Ok(vec![true]));
    let (a, b) = (load(0, 0, &[L], by(0)), load(1, 1, &[L], by(1)));
    let misses = record.unordered_slot(0, &[&a, &b], L).unwrap();
    assert_eq!(misses, vec![(0, vec![false]), (1, vec![false])]);
    assert_eq!(record.offset_slots[&0].holds, SlotHolds::Last);
    // A (L, 150) and then A again with B (L): A's second load leaves 150 before B or not, which
    // decides B's miss: refused.
    let mut record = VolumeLookups::default();
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![true, false])
    );
    let (a, b) = (load(0, 0, &[L, 150], by(0)), load(1, 1, &[L], by(1)));
    assert!(record.unordered_slot(0, &[&a, &b], L).is_err());
    // Where NCBI keeps every sequence, A's second load is none for certain: B misses after A's 150.
    let mut record = VolumeLookups {
        kept: true,
        ..VolumeLookups::default()
    };
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![true, false])
    );
    let misses = record.unordered_slot(0, &[&a, &b], L).unwrap();
    assert_eq!(misses, vec![(0, vec![false, false]), (1, vec![true])]);
    // Two loads of one sequence by two ids in one group, neither found in the other's Bioseq: two
    // loads, whose order decides where each looks its lines up in another order (L, 150 and 150,
    // L: once or twice), and not where both look L up (once).
    let mut record = VolumeLookups::default();
    let (one, two) = (load(0, 0, &[L, 150], by(0)), load(1, 0, &[150, L], by(1)));
    assert!(record.unordered_slot(0, &[&one, &two], L).is_err());
    let (one, two) = (load(0, 0, &[L], by(0)), load(1, 0, &[L], by(1)));
    let misses = record.unordered_slot(0, &[&one, &two], L).unwrap();
    assert_eq!(misses, vec![(0, vec![true]), (1, vec![false])]);
    // The second id found in the first load's Bioseq: in one prefetch a load of its own all the
    // same (NCBI's prefetch asks its scope only for ids resolved before it; review 1 of item
    // O4b2-7, M2).
    let mut record = VolumeLookups::default();
    let two = load(1, 0, &[L], found_in(1, 0));
    let misses = record.unordered_slot(0, &[&one, &two], L).unwrap();
    assert_eq!(misses, vec![(0, vec![true]), (1, vec![false])]);
    assert!(record.slot_blobs[&0][1].sure);
}

/// [`by`] of `id`, whose Bioseq also holds the Seq-ids of `others` (by their handles).
fn holding(id: u32, others: &[u32]) -> BlobName {
    let mut name = by(id);
    name.holds.handles = std::iter::once(id)
        .chain(others.iter().copied())
        .map(|id| id.to_le_bytes().to_vec())
        .collect();
    name
}

/// `VolumeLookups::prefetch_slots`: one prefetch loads every handle of a sequence, also one that
/// the Bioseq of its first load holds (NCBI's prefetch asks its scope only for the ids that it
/// resolved before, then its loader; review 1 of item O4b2-7, M2), each for certain: within the
/// same live scope a succeeding prefetch pins its blobs so a second prefetch and a report read by
/// an exact handle of the first Bioseq load nothing; after `reset_scope`, a later batch's prefetch
/// of both loads nothing where NCBI keeps every sequence and is refused otherwise, and a report
/// read by id 3 loads fresh because the cleared scope no longer indexes the first Bioseq.
#[test]
fn a_prefetch_loads_each_handle_of_a_sequence() {
    const L: i32 = 99990;
    let loads = || {
        [(0, 1, holding(0, &[1, 3])), (1, 2, by(1))].map(|(ordinal, order, name)| SlotLoad {
            ordinal,
            blob: 0,
            volume: 0,
            lookups: vec![L, 150],
            order,
            name,
            pins: true,
            maybe: false,
        })
    };
    for kept in [false, true] {
        let mut record = VolumeLookups {
            kept,
            ..VolumeLookups::default()
        };
        let misses = record.prefetch_slots(&loads(), |_| false, L, None).unwrap();
        assert_eq!(misses[&0], vec![true, false]);
        assert_eq!(misses[&1], vec![true, false]);
        assert!(record.prefetch_from.is_none());
        assert!(record.slot_blobs[&0]
            .iter()
            .all(|blob| blob.sure && blob.in_scope && blob.pinned));
        let after = (record.offset_slots.clone(), record.slot_blobs.clone());
        // Within the same live scope (e.g. outfmt 5 keeping scope), pinned blobs load nothing.
        let same_scope = record.prefetch_slots(&loads(), |_| false, L, None).unwrap();
        assert_eq!(same_scope[&0], vec![false, false]);
        assert_eq!(same_scope[&1], vec![false, false]);
        // A report read in the same live scope by exact handle 3 hits the first Bioseq's handles.
        (record.offset_slots, record.slot_blobs) = after.clone();
        assert_eq!(record.blob_read(0, &by(3)), BlobRead::Kept);
        assert_eq!(
            record.offset_slot(0, 0, &by(3), &[L, 150], L),
            Ok(vec![false, false])
        );
        assert_eq!(record.slot_blobs[&0].len(), 2);
        // After a batch scope reset, `in_scope` and `pinned` clear: only `kept` keeps the blobs.
        (record.offset_slots, record.slot_blobs) = after;
        record.reset_scope();
        let after_reset = (record.offset_slots.clone(), record.slot_blobs.clone());
        let again = record.prefetch_slots(&loads(), |_| false, L, None);
        if kept {
            let again = again.unwrap();
            assert_eq!(again[&0], vec![false, false]);
            assert_eq!(again[&1], vec![false, false]);
        } else {
            assert!(again.is_err());
        }
        // After reset, id 3 is no longer in live scope: fresh load (misses after 150).
        (record.offset_slots, record.slot_blobs) = after_reset;
        assert_eq!(record.blob_read(0, &by(3)), BlobRead::Fresh);
        assert_eq!(
            record.offset_slot(0, 0, &by(3), &[L, 150], L),
            Ok(vec![true, false])
        );
    }
}

/// The bounds of NCBI's object manager: a database of at most 10 sequences keeps every load, of
/// at most 20 loses none at the reports' reads; an environment variable of its parameters leaves
/// both to LOSAT's refusals (review 1 of item O4b2-5, M2).
#[test]
fn the_object_managers_bounds_follow_its_sizes() {
    let cwd = empty_cwd("om_bounds");
    let three = db_search("tm_cache150", Some(&variant("offset_last_past_end")), &cwd);
    let twelve = db_search("tm_cache_gi", Some(&variant("offset_last_past_end")), &cwd);
    if default_object_manager() {
        assert!(keeps_every_load(&three));
        assert!(!reads_may_load_again(&twelve));
    }
    assert!(!keeps_every_load(&twelve));
    let var = |name: &str| {
        (
            std::ffi::OsString::from(name),
            std::ffi::OsString::from("0"),
        )
    };
    assert!(object_manager_default_in([var("PATH")]));
    assert!(!object_manager_default_in([var("OBJMGR_BLOB_CACHE")]));
    assert!(!object_manager_default_in([var(
        "objmgr_scope_autorelease_size"
    )]));
}

#[test]
fn the_slot_of_the_offset_errors_taxid_follows_the_loads() {
    const L: i32 = 99990;
    let fetch = |ordinal, taxid, again| SlotFetch {
        ordinal,
        taxid,
        again,
        maybe: false,
    };
    // H1 (cc_m): one load of L, 150, L; again: its own 150 is between.
    let mut record = VolumeLookups::default();
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150, L], L),
        Ok(vec![true, false, true])
    );
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150, L], L),
        Err(SlotUncertain {
            ordinal: 0,
            again: true,
            after: Some(fetch(0, 150, true)),
            order: None,
            unnamed: false,
        })
    );
    // M1 (b) (cc_m2): A (L, 150) and B (L), first loads: two texts. Read again, A may put 150
    // there (Either), and B's read again is not certain.
    let mut record = VolumeLookups::default();
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![true, false])
    );
    assert_eq!(record.offset_slot(0, 1, &by(1), &[L], L), Ok(vec![true]));
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![false, false])
    );
    assert_eq!(
        record.offset_slot(0, 1, &by(1), &[L], L),
        Err(SlotUncertain {
            ordinal: 1,
            again: true,
            after: Some(fetch(0, 150, true)),
            order: None,
            unnamed: false,
        })
    );
    // A first load after such a read again (C of 99990 + 256 * 2 is not L).
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &by(0), &[L], L), Ok(vec![true]));
    assert_eq!(record.offset_slot(0, 1, &by(1), &[150], L), Ok(vec![false]));
    assert_eq!(record.offset_slot(0, 2, &by(2), &[L], L), Ok(vec![true]));
    assert_eq!(record.offset_slot(0, 1, &by(1), &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 3, &by(3), &[9606, L + 512, L], L),
        Ok(vec![false, false, true])
    );
    assert_eq!(record.offset_slot(0, 1, &by(1), &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 4, &by(4), &[L], L),
        Err(SlotUncertain {
            ordinal: 4,
            again: false,
            after: Some(fetch(1, 150, true)),
            order: None,
            unnamed: false,
        })
    );
    // M1 (a) (cc_ax): 150 never looked up: L's reads again write nothing.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &by(0), &[L], L), Ok(vec![true]));
    for _ in 0..3 {
        assert_eq!(record.offset_slot(0, 0, &by(0), &[L], L), Ok(vec![false]));
    }
    // cc_ac: A (L), C (150), then A again: refused, after C's fetch.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(0, 0, &by(0), &[L], L), Ok(vec![true]));
    assert_eq!(record.offset_slot(0, 1, &by(1), &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 0, &by(0), &[L], L),
        Err(SlotUncertain {
            ordinal: 0,
            again: true,
            after: Some(fetch(1, 150, false)),
            order: None,
            unnamed: false,
        })
    );
    // C, then A in another volume (cc_c1 cc_a1), then both again; taxids of other slots.
    let mut record = VolumeLookups::default();
    assert_eq!(record.offset_slot(1, 0, &by(0), &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 1, &by(1), &[L, 9606], L),
        Ok(vec![true, false])
    );
    assert_eq!(record.offset_slot(1, 0, &by(0), &[150], L), Ok(vec![false]));
    assert_eq!(
        record.offset_slot(0, 1, &by(1), &[L, 9606], L),
        Ok(vec![false, false])
    );
    assert_eq!(
        record.offset_slot(0, 5, &by(5), &[9606, 151], L),
        Ok(vec![false, false])
    );
    assert!(!record.slot_blobs.contains_key(&5));
}

/// `pssm_load_texts` on `tm_cache150` (st1 99990, st2 150, st3 99990, parsed Seq-ids, one
/// volume) with the offset error of `offset_last_past_end`: the PSSM's loads in their order
/// write NCBI's texts, a read again loads nothing where NCBI keeps every sequence (three), and one
/// that NCBI's scope decides (`tm_cache_gi`, twelve) is refused, naming the fetch of 150; with 151
/// (`tm_cache151`) the read again writes nothing.
#[test]
fn the_pssms_loads_write_the_offset_error_as_the_slot_holds_it() {
    let cwd = empty_cwd("pssm_slot");
    let offset = "Error: [psiblast] Error: Offset error at end of taxdb file.\n";
    let loads = |name: &str, ordinals: &[u32], reset_at: Option<usize>| -> Result<String> {
        let db = db_search(name, Some(&variant("offset_last_past_end")), &cwd);
        let parsed = parsed_subjects(&db);
        let run = TaxdbRun::default();
        let mut text = Vec::new();
        for (index, &ordinal) in ordinals.iter().enumerate() {
            if reset_at == Some(index) {
                run.reset_scope();
            }
            let subject = &parsed[ordinal as usize];
            let lookups = blob_taxids(&db, ordinal, Some(subject), ordinal, "PSIBLAST")?;
            text.extend(pssm_load_texts(
                &db,
                &run,
                None,
                (ordinal, ordinal),
                &lookups,
                "PSIBLAST",
            )?);
        }
        Ok(String::from_utf8(text).unwrap())
    };
    assert_eq!(
        loads("tm_cache150", &[0, 1, 2], None).unwrap(),
        offset.repeat(2)
    );
    assert_eq!(loads("tm_cache150", &[0, 2, 1], None).unwrap(), offset);
    assert_eq!(loads("tm_cache150", &[0, 0, 2], None).unwrap(), offset);
    // NCBI keeps every sequence of a database of three across resets: the read again loads nothing.
    assert_eq!(loads("tm_cache150", &[0, 1, 0], Some(2)).unwrap(), offset);
    // Of twelve (`tm_cache_gi`: 0 to 5 hold 99990, 6 to 11 hold 150), within one query's live
    // scope (no `ResetScopeHistory`) 12 <= 20 fits in the scope queue + DS cache; after a query
    // reset, NCBI may have released it from the 10-slot DS cache.
    if default_object_manager() {
        assert_eq!(loads("tm_cache_gi", &[0, 6, 0], None).unwrap(), offset);
    }
    let error = loads("tm_cache_gi", &[0, 6, 0], Some(2)).unwrap_err();
    assert!(error.is::<TaxonomyRefusal>(), "{error:#}");
    let text = format!("{error:#}");
    for part in [
        "taxid 99990, the last entry of the taxonomy database",
        "reads database sequence 0 of the volume",
        "again after the run's fetch of database sequence 6 with taxid 150",
        "256-slot taxonomy cache",
        "LOSAT's PSIBLAST",
    ] {
        assert!(text.contains(part), "{part}: {text}");
    }
    assert_eq!(
        loads("tm_cache151", &[0, 1, 0, 2], Some(2)).unwrap(),
        offset
    );
}

#[test]
fn a_database_without_taxids_does_not_hold() {
    let db = database("sprot", Some(&variant("high_bytes_utf8")));
    assert!(db.taxdb_found());
    assert!(!taxdb_may_refuse(&db));
}

/// Review 3 of line DBID pfoff-4, M1 and L2: the sequences that NCBI may give for a Seq-align id
/// whose lookup fails on a file, without that lookup, after the calls of the run before this one
/// (`HeaderReader::unrecorded_holders`, read with `unrecorded`), in lists
/// of `tests/db_search/multidb` beside `pU`, whose unparsed volume makes the lookup of a text id
/// fail (NCBI BLAST+ 2.17.0 holds or throws as rows `db.pfo.bp.*.held.*` of
/// `db_search_manifest/blastp.pfo.tsv` show). The data loader gives a sequence it loaded that
/// holds the id's own handle in any of its definition lines (`Lines::Every`); the scope one that
/// holds an id it takes for it, or matches it weakly, in the line its load made the Bioseq of (the
/// line of the sequence's own Seq-align id, `Lines::Loaded`). Holders by their Seq-align ids.
#[test]
fn the_scope_holds_the_line_of_a_load_and_the_loader_every_line() {
    use crate::report::db_subject::{Holders, ParsedSubject};
    let cwd = empty_cwd("unrecorded_holders");
    let open = |names: &str| {
        let search = SearchPath::new(&cwd, None, None);
        let list: Vec<String> = names
            .split(' ')
            .map(|name| {
                let path = root().join("db_search/multidb").join(name);
                path.to_str().unwrap().to_owned()
            })
            .collect();
        SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Protein, &search).unwrap()
    };
    let subject = |db: &SeqDb, ordinal: u32| -> Option<std::sync::Arc<ParsedSubject>> {
        match read_header(&HeaderReader::new(db, true), ordinal, 0, "BLASTP").unwrap() {
            ReportedHeader::Parsed(subject) => Some(subject),
            _ => None,
        }
    };
    let id = |subject: &ParsedSubject| String::from_utf8(subject.alignment_fasta()).unwrap();
    let names = |db: &SeqDb, holders: &Holders| match holders {
        Holders::Listed(listed) => listed
            .iter()
            .map(|holder| id(&subject(db, holder.ordinal()).unwrap()))
            .collect::<Vec<_>>(),
        Holders::Unknown(why) => panic!("unknown holders: {why}"),
    };
    // (list, the failing Seq-align id, the loader's holders, the scope's holders)
    let cases: [(&str, &str, &[&str], &[&str]); 7] = [
        // A weak holder (gi|5001|ref|NP_100001.1| for the unversioned id) and one of another type
        // (gi|5004|emb|NP_100001.1| for ref): the scope's only.
        ("vpP zzgpP pU", "ref|NP_100001|", &[], &["gi|5001"]),
        ("epP pP pU", "ref|NP_100001.1|", &[], &["gi|5004"]),
        // A holder of the id's own handle (gi|7022|gnl|TDB|TAG1): both places.
        ("gnpF zgpH pU", "gnl|TDB|TAG1", &["gi|7022"], &["gi|7022"]),
        // zmdpH: gi|7051|ref|XP_900001.1| then gi|7052|ref|NP_100001.1|, loaded by its first
        // line (its Seq-align id gi|7051): the loader's holder of ref|NP_100001.1| by the second
        // line, the scope's of nothing (review 3's families mdP, mdV and mdE).
        ("pP zmdpH pU", "ref|NP_100001.1|", &["gi|7051"], &[]),
        ("vpP zmdpH pU", "ref|NP_100001|", &[], &[]),
        ("enpF zmdpH pU", "emb|NP_100001.1|", &[], &[]),
        // zmgpH: ref|XP_900002.1| then gi|7053|ref|NP_100001.1|, loaded by its second line (the
        // line of its Seq-align id gi|7053): the scope's weak holder by that line.
        ("vpP zmgpH pU", "ref|NP_100001|", &[], &["gi|7053"]),
    ];
    for (list, failing, loader, scope) in cases {
        let db = open(list);
        let failed = (0..db.num_seqs())
            .filter_map(|ordinal| subject(&db, ordinal))
            .find(|subject| id(subject) == failing)
            .unwrap_or_else(|| panic!("{list}: no {failing}"));
        assert!(
            matches!(
                failed.found(),
                Some(crate::report::db_subject::FoundInstead::Failed(_))
            ),
            "{list}: {:?}",
            failed.found()
        );
        let holders = failed.unrecorded();
        assert_eq!(names(&db, &holders.loader), loader, "{list}: loader");
        assert_eq!(names(&db, &holders.scope), scope, "{list}: scope");
    }
}

#[test]
fn a_prefetch_with_a_failed_lookup_takes_a_load_of_the_cache_slot_as_made_or_not() {
    // Line NP item O4b2-3, review 1 F1: NCBI's prefetch ends its loads at a Seq-align id whose
    // lookup fails, so whether it loads sZ (150) before it, which puts 150 in the slot of 99990
    // without a text, is not known; a later load of 99990 would write the offset error or not.
    // Line NP item O4b2-10, review F2: so the slot holds 99990 or 150 after it, which only a later
    // lookup of 99990 shows (refused there); where it held no 99990 before, either leaves no 99990.
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    let warning = "Warning: [blastp] Error pre-fetching sequence data\n";
    let cwd = empty_cwd("prefetch_slot_stop");
    let taxdb = variant("offset_last_past_end");
    let db = db_search("nz_pc", Some(&taxdb), &cwd);
    let parsed = parsed_subjects(&db);
    let throws =
        parsed[1]
            .clone()
            .found_instead(Some(crate::report::db_subject::FoundInstead::ThrowsOnId(
                crate::ncbi_exception::NcbiException::new(
                    crate::seqdb::isam::OBJECT_ID_GET_ID,
                    "Invalid choice selection",
                ),
            )));
    assert!(throws.prefetch_fails());
    let run = TaxdbRun::default();
    let batch = |subjects: Vec<PrefetchSubject<'_>>| {
        batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            subjects,
            "BLASTP",
            &stop_queries(0),
        )
        .map(|posts| String::from_utf8(posts.text()).unwrap())
    };
    assert_eq!(
        batch(vec![PrefetchSubject::unflipped(3, Some(&parsed[3]))]).unwrap(),
        offset
    );
    // sM (10090) takes another slot: the batch writes NCBI's warning alone.
    assert_eq!(
        batch(vec![
            PrefetchSubject::unflipped(1, Some(&throws)),
            PrefetchSubject::unflipped(2, Some(&parsed[2])),
        ])
        .unwrap(),
        warning
    );
    let tie = || {
        vec![
            PrefetchSubject::unflipped(1, Some(&throws)),
            PrefetchSubject::unflipped(0, Some(&parsed[0])),
        ]
    };
    let fresh = TaxdbRun::default();
    let prefetch = |run: &TaxdbRun, subjects: Vec<PrefetchSubject<'_>>| {
        batch_prefetch(
            &db,
            run,
            std::iter::empty(),
            subjects,
            "BLASTP",
            &stop_queries(0),
        )
        .map(|posts| String::from_utf8(posts.text()).unwrap())
    };
    assert_eq!(prefetch(&fresh, tie()).unwrap(), warning);
    assert_eq!(
        fresh.volume_lookups(&db).offset_slots[&0].holds,
        SlotHolds::Other
    );
    // After a load of sN (99990): 99990 or 150, and a later first load of 99990 is refused.
    let after = TaxdbRun::default();
    assert_eq!(
        prefetch(
            &after,
            vec![PrefetchSubject::unflipped(3, Some(&parsed[3]))]
        )
        .unwrap(),
        offset
    );
    assert_eq!(prefetch(&after, tie()).unwrap(), warning);
    let mut lookups = after.volume_lookups(&db);
    assert_eq!(lookups.offset_slots[&0].holds, SlotHolds::Either);
    let other = BlobName {
        handle: Some(b"another handle".to_vec()),
        ..BlobName::default()
    };
    let refused = slot_uncertain_error(
        &db,
        0,
        99990,
        lookups
            .offset_slot(0, 3, &other, &[99990], 99990)
            .unwrap_err(),
    )
    .to_string();
    assert!(
        refused.contains(
            "loaded database sequence 0 with taxid 150 before that failure or not at all"
        ),
        "{refused}"
    );
}

/// Line NP item O4b2-10 (`tests/db_search/ids/make_ids.py --pfworder-only`, probe `O4b2/p10` of
/// NCBI 2.17.0): a lookup that throws on its id ends the prefetch's loads in NCBI's order of the
/// batch's Seq-id handles. Every sequence has taxid 99990 (the offset error of
/// `offset_last_past_end`) but the last of `ntagc_t5` and `ptagc_t5` (150, the same slot of a
/// volume's taxonomy cache: a colliding volume). `ntag_t5`: `gb|FFF00092.1|` (packed 92) loads
/// before `gnl|BL_ORD_ID|0` (a text tag, not packed); `ntagx_t5`: `gb|FFF00094.1|` (94) after
/// `gnl|BL_ORD_ID|x002` (packed 2); `ntagc_t5`: `gb|FFF00010.1|` (10) before `gnl|BL_ORD_ID|x050`
/// (50), `gb|FFF00090.1|` (90, taxid 150) after it; `ptagx_t5` and `ptagc_t5` the same with
/// x002 and x051. `pshare_t5`: `lcl|555` and the text tag are both unpacked (one key, ordered by
/// address), `gb|FFF00032.1|` (32) before.
#[test]
fn a_failed_lookup_ends_the_prefetchs_loads_in_ncbis_order_of_seq_id_handles() {
    use LoadBeforeStop::{No, Unknown, Yes};
    let cwd = empty_cwd("pfworder");
    let taxdb = variant("offset_last_past_end");
    let open = |name: &str, kind: SeqType| {
        let search = SearchPath::new(&cwd, Some(taxdb.as_os_str()), None);
        let path = root().join("db_search/ids").join(name);
        SeqDb::open_in(path.to_str().unwrap().as_bytes(), kind, &search).unwrap()
    };
    let parsed = |db: &SeqDb| -> Vec<std::sync::Arc<crate::report::db_subject::ParsedSubject>> {
        let headers = HeaderReader::new(db, false);
        (0..db.num_seqs())
            .map(
                |ordinal| match read_header(&headers, ordinal, 0, "BLASTN").unwrap() {
                    ReportedHeader::Parsed(subject) => subject,
                    _ => panic!("{ordinal}: a sequence with parsed Seq-ids"),
                },
            )
            .collect()
    };
    let offset =
        |program: &str| format!("Error: [{program}] Error: Offset error at end of taxdb file.\n");
    let warning =
        |program: &str| format!("Warning: [{program}] Error pre-fetching sequence data\n");
    for (name, kind, program, failed, want, text) in [
        ("ntag_t5", SeqType::Nucleotide, "blastn", 0, vec![Yes], true),
        (
            "ntagx_t5",
            SeqType::Nucleotide,
            "blastn",
            0,
            vec![No],
            false,
        ),
        (
            "ntagc_t5",
            SeqType::Nucleotide,
            "blastn",
            1,
            vec![Yes, No],
            true,
        ),
        ("ptagx_t5", SeqType::Protein, "blastp", 0, vec![No], false),
        (
            "ptagc_t5",
            SeqType::Protein,
            "blastp",
            1,
            vec![Yes, No],
            true,
        ),
    ] {
        let db = open(name, kind);
        let parsed = parsed(&db);
        assert!(parsed[failed as usize].prefetch_fails(), "{name}");
        let loaded: Vec<(u32, Option<&crate::report::db_subject::ParsedSubject>)> = (0..db
            .num_seqs())
            .filter(|&ordinal| ordinal != failed)
            .map(|ordinal| (ordinal, Some(&*parsed[ordinal as usize])))
            .collect();
        let failing = [PrefetchSubject::unflipped(
            failed,
            Some(&*parsed[failed as usize]),
        )];
        assert_eq!(
            loads_before_failure(&db, &loaded, &failing, None),
            want,
            "{name}"
        );
        // The batch in any order: the loads before the failure write the offset error, then the
        // warning; those after it write nothing (on `*tagc_t5` also no step through the slot).
        let subjects = failing.iter().copied().chain(
            loaded
                .iter()
                .rev()
                .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed)),
        );
        let posts = batch_prefetch(
            &db,
            &TaxdbRun::default(),
            std::iter::empty(),
            subjects,
            &program.to_ascii_uppercase(),
            &stop_queries(0),
        )
        .unwrap();
        let expected = if text {
            offset(program) + &warning(program)
        } else {
            warning(program)
        };
        assert_eq!(String::from_utf8(posts.text()).unwrap(), expected, "{name}");
    }
    // `pshare_t5`: the load of `lcl|555` may come before the failed lookup or not at all: refused
    // where it writes the offset error, not where a load before it (`gb|FFF00032.1|`) wrote it.
    let db = open("pshare_t5", SeqType::Protein);
    let parsed = parsed(&db);
    let (local, failed, control) = (0, 2, 3);
    assert!(parsed[failed].prefetch_fails());
    let subject =
        |ordinal: usize| PrefetchSubject::unflipped(ordinal as u32, Some(&*parsed[ordinal]));
    let pair = [
        (local as u32, Some(&*parsed[local])),
        (control as u32, Some(&*parsed[control])),
    ];
    assert!(matches!(
        loads_before_failure(&db, &pair, &[subject(failed)], None)[..],
        [Unknown(why), Yes] if why == ADDRESS_ORDER
    ));
    let batch = |run: &TaxdbRun, subjects: &[usize]| {
        batch_prefetch(
            &db,
            run,
            std::iter::empty(),
            subjects.iter().map(|&ordinal| subject(ordinal)),
            "BLASTP",
            &stop_queries(0),
        )
        .map(|posts| String::from_utf8(posts.text()).unwrap())
    };
    let refused = batch(&TaxdbRun::default(), &[failed, local])
        .unwrap_err()
        .to_string();
    assert!(
        refused.contains("in an order that LOSAT does not reproduce")
            && refused
                .contains("database sequence 0, which the prefetch may load before or not at all")
            && refused.contains("makes the taxonomy database write a text"),
        "{refused}"
    );
    assert_eq!(
        batch(&TaxdbRun::default(), &[failed, local, control]).unwrap(),
        offset("blastp") + &warning("blastp")
    );
}

/// The databases of `tests/db_search/ids` (`make_ids.py --pfworder-only`) of line NP item O4b2-10,
/// opened with `BLASTDB=offset_last_past_end` from the empty current directory `cwd`, and the subject
/// of each sequence as the reports' reader makes it (a text tag's lookup fails:
/// `ParsedSubject::prefetch_fails`).
fn pfworder_db(
    cwd: &Path,
    name: &str,
    kind: SeqType,
) -> (
    SeqDb,
    Vec<std::sync::Arc<crate::report::db_subject::ParsedSubject>>,
) {
    let taxdb = variant("offset_last_past_end");
    let search = SearchPath::new(cwd, Some(taxdb.as_os_str()), None);
    let path = root().join("db_search/ids").join(name);
    let db = SeqDb::open_in(path.to_str().unwrap().as_bytes(), kind, &search).unwrap();
    let headers = HeaderReader::new(&db, false);
    let parsed = (0..db.num_seqs())
        .map(
            |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                ReportedHeader::Parsed(subject) => subject,
                _ => panic!("{ordinal}: a sequence with parsed Seq-ids"),
            },
        )
        .collect();
    (db, parsed)
}

/// A copy of `subject` whose lookup throws on its id where the prefetch loads it
/// (`ParsedSubject::prefetch_fails`).
fn throwing(
    subject: &crate::report::db_subject::ParsedSubject,
) -> crate::report::db_subject::ParsedSubject {
    subject
        .clone()
        .found_instead(Some(crate::report::db_subject::FoundInstead::ThrowsOnId(
            crate::ncbi_exception::NcbiException::new(
                crate::seqdb::isam::OBJECT_ID_GET_ID,
                "Invalid choice selection",
            ),
        )))
}

/// Line NP item O4b2-10, review F1 (`make_ids.py --pfworder-only`, probe `O4b2/impl10/r2` of NCBI
/// 2.17.0 in four heap layouts): the blobs that a prefetch's loads lock before a lookup that throws
/// go to NCBI's data source's cache of 10, not to the scope, released from the largest key down.
/// `pslot11_t5`: `gb|FFF00090.1|` (ordinal 0, taxid 150), `gb|FFF00011.1|` to `gb|FFF00001.1|`
/// (ordinals 1 to 11, keys 11 to 1, taxid 99990), the failing `gnl|BL_ORD_ID|x050` (key 50);
/// `pslot12_t5` and `pslot10_t5` with 12 and 10 loads before the failure; `pslotc_t5` the control
/// without a failing lookup. `pslotr_t5` holds key 11 at ordinal 0 and `pslotm_t5` key 1 there, key 11
/// at ordinal 11. A read again of a dropped blob loads it (its lookups run again: NCBI's second offset
/// error after the 150 load), one that the cache kept loads nothing.
#[test]
fn a_failed_lookups_loads_go_to_the_data_sources_cache_of_ten() {
    let cwd = empty_cwd("pfwcache");
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    let warning = "Warning: [blastp] Error pre-fetching sequence data\n";
    let prefetch =
        |db: &SeqDb,
         run: &TaxdbRun,
         parsed: &[std::sync::Arc<crate::report::db_subject::ParsedSubject>]| {
            batch_prefetch(
                db,
                run,
                std::iter::empty(),
                (0..db.num_seqs()).rev().map(|ordinal| {
                    PrefetchSubject::unflipped(ordinal, Some(&*parsed[ordinal as usize]))
                }),
                "BLASTP",
                &stop_queries(0),
            )
            .map(|posts| String::from_utf8(posts.text()).unwrap())
        };
    let read = |db: &SeqDb, run: &TaxdbRun, ordinal: u32| {
        run.volume_lookups(db)
            .blob_read(ordinal, &blob_name(db, ordinal, ordinal, "BLASTP").unwrap())
    };
    for (name, dropped, kept, fails) in [
        ("pslot11_t5", 1..2, 2..12, true),
        ("pslot12_t5", 1..3, 3..13, true),
        ("pslot10_t5", 1..1, 1..11, true),
        ("pslotc_t5", 1..1, 1..12, false),
        ("pslotr_t5", 0..1, 2..12, true),
        ("pslotm_t5", 11..12, 2..11, true),
    ] {
        let (db, parsed) = pfworder_db(&cwd, name, SeqType::Protein);
        let failing = db.num_seqs() - 1;
        assert_eq!(parsed[failing as usize].prefetch_fails(), fails, "{name}");
        let run = TaxdbRun::default();
        let expected = if fails {
            offset.to_owned() + warning
        } else {
            offset.to_owned()
        };
        assert_eq!(prefetch(&db, &run, &parsed).unwrap(), expected, "{name}");
        for ordinal in dropped {
            assert_eq!(
                read(&db, &run, ordinal),
                BlobRead::Fresh,
                "{name}: {ordinal}"
            );
        }
        for ordinal in kept {
            assert_eq!(
                read(&db, &run, ordinal),
                BlobRead::Kept,
                "{name}: {ordinal}"
            );
        }
        // The cache keeps them to the end of the batch only: NCBI may drop them after
        // (`ResetScopeHistory` releases the scope's blobs to the same cache).
        let mut lookups = run.volume_lookups(&db);
        lookups.reset_scope();
        let name_of = blob_name(&db, 5, 5, "BLASTP").unwrap();
        assert_eq!(lookups.blob_read(5, &name_of), BlobRead::Again, "{name}");
    }
}

/// Line NP item O4b2-10, review F2 (`ptie_t5`, `ptie150_t5` of `make_ids.py --pfworder-only`,
/// NCBI 2.17.0 SAME in four heap layouts): `gb|FFF00050.1|` (ordinal 1) has the packed number of the
/// failing `gnl|BL_ORD_ID|x050` (ordinal 2), so NCBI loads it before the failure or not at all, by
/// the addresses of their handles; `gb|FFF00010.1|` (0, taxid 99990) before them, `gb|FFF00090.1|`
/// (3, taxid 150) after. With taxid 99990 the load finds 99990 in the slot either way: no refusal,
/// the slot holds 99990. With taxid 150 it leaves 99990 or 150, which only a later lookup of
/// 99990 would show: no refusal there, and that lookup is refused.
#[test]
fn a_load_of_the_failures_key_is_made_or_not() {
    use LoadBeforeStop::{No, Unknown, Yes};
    const L: i32 = 99990;
    let cwd = empty_cwd("pfwtie");
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    let warning = "Warning: [blastp] Error pre-fetching sequence data\n";
    for (name, holds) in [
        ("ptie_t5", SlotHolds::Last),
        ("ptie150_t5", SlotHolds::Either),
    ] {
        let (db, parsed) = pfworder_db(&cwd, name, SeqType::Protein);
        assert!(parsed[2].prefetch_fails(), "{name}");
        let loaded: Vec<(u32, Option<&crate::report::db_subject::ParsedSubject>)> = [0, 1, 3]
            .into_iter()
            .map(|ordinal| (ordinal, Some(&*parsed[ordinal as usize])))
            .collect();
        let failing = [PrefetchSubject::unflipped(2, Some(&*parsed[2]))];
        assert!(matches!(
            loads_before_failure(&db, &loaded, &failing, None)[..],
            [Yes, Unknown(why), No] if why == ADDRESS_ORDER
        ));
        let run = TaxdbRun::default();
        let posts = batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            (0..4).map(|ordinal| {
                PrefetchSubject::unflipped(ordinal, Some(&*parsed[ordinal as usize]))
            }),
            "BLASTP",
            &stop_queries(0),
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(posts.text()).unwrap(),
            offset.to_owned() + warning
        );
        let mut lookups = run.volume_lookups(&db);
        assert_eq!(lookups.offset_slots[&0].holds, holds, "{name}");
        // NCBI may have made the tie's load: a read again of it is a load or none.
        assert_eq!(
            lookups.blob_read(1, &blob_name(&db, 1, 1, "BLASTP").unwrap()),
            BlobRead::Again,
            "{name}"
        );
        if holds == SlotHolds::Either {
            let other = BlobName {
                handle: Some(b"another handle".to_vec()),
                ..BlobName::default()
            };
            let cause = lookups.offset_slot(0, 0, &other, &[L], L).unwrap_err();
            assert_eq!(
                cause.after,
                Some(SlotFetch {
                    ordinal: 1,
                    taxid: 150,
                    again: false,
                    maybe: true,
                })
            );
            let refusal = slot_uncertain_error(&db, 0, L, cause).to_string();
            assert!(
                refusal.contains(
                    "loaded database sequence 1 with taxid 150 before that failure or not at all"
                ) && refusal.contains("as that prefetch loaded database sequence 1 or not"),
                "{refusal}"
            );
        }
    }
}

/// Line NP item O4b2-10, re-review R2-1 and R2-2 (`make_ids.py --pfworder-only`, probe
/// `O4b2/impl10/r3` of NCBI 2.17.0 in four heap layouts): every blob that a load locks before the
/// failed lookup counts among those that NCBI's data source releases after the loads of larger keys,
/// whatever its taxids; a load of the failure's key that NCBI makes or not, after which 10 or more
/// are released, is held nowhere either way. `pslotmix_t5`: `gb|FFF00090.1|` (ordinal 0, taxid 150),
/// `gb|FFF00011.1|` to `gb|FFF00002.1|` (ordinals 1 to 10, keys 11 to 2, taxid 99990),
/// `gb|FFF00001.1|` (11, key 1, taxid 9606), the failing `gnl|BL_ORD_ID|x050`: 10 releases follow
/// key 11, which is dropped. `pslotmix12_t5`: `gb|FFF00011.1|` (1) of taxid 9606, released first,
/// then keys 10 to 1 (ordinals 2 to 11, 99990): all kept. `pslottie_t5`: keys 10 to 1 (ordinals 1 to
/// 10, 99990), `gb|FFF00050.1|` (11, 99990) of the failure's key, then the failure: the tie is
/// loaded again at its read, made or not. Where the scope may hold a sequence that the records of
/// the slot do not name (a read since its reset), the load of another slot may ask no data source:
/// key 11 is a load or none at a read again.
#[test]
fn every_load_before_a_failed_lookup_counts_among_the_releases_of_the_data_sources_cache() {
    let cwd = empty_cwd("pfwmix");
    let offset = "Error: [blastp] Error: Offset error at end of taxdb file.\n";
    let warning = "Warning: [blastp] Error pre-fetching sequence data\n";
    let prefetch =
        |db: &SeqDb,
         run: &TaxdbRun,
         parsed: &[std::sync::Arc<crate::report::db_subject::ParsedSubject>]| {
            batch_prefetch(
                db,
                run,
                std::iter::empty(),
                (0..db.num_seqs()).rev().map(|ordinal| {
                    PrefetchSubject::unflipped(ordinal, Some(&*parsed[ordinal as usize]))
                }),
                "BLASTP",
                &stop_queries(0),
            )
            .map(|posts| String::from_utf8(posts.text()).unwrap())
        };
    let read = |db: &SeqDb, run: &TaxdbRun, ordinal: u32| {
        run.volume_lookups(db)
            .blob_read(ordinal, &blob_name(db, ordinal, ordinal, "BLASTP").unwrap())
    };
    for (name, fresh, kept) in [
        ("pslotmix_t5", 1..2, 2..11),
        ("pslotmix12_t5", 1..2, 2..12),
        ("pslottie_t5", 11..12, 1..11),
    ] {
        let (db, parsed) = pfworder_db(&cwd, name, SeqType::Protein);
        assert!(parsed[12].prefetch_fails(), "{name}");
        let run = TaxdbRun::default();
        assert_eq!(
            prefetch(&db, &run, &parsed).unwrap(),
            offset.to_owned() + warning,
            "{name}"
        );
        for ordinal in fresh {
            assert_eq!(
                read(&db, &run, ordinal),
                BlobRead::Fresh,
                "{name}: {ordinal}"
            );
        }
        for ordinal in kept {
            assert_eq!(
                read(&db, &run, ordinal),
                BlobRead::Kept,
                "{name}: {ordinal}"
            );
        }
    }
    // A read of the live scope before the prefetch: the load of taxid 9606 may be resolved there.
    let (db, parsed) = pfworder_db(&cwd, "pslotmix_t5", SeqType::Protein);
    let run = TaxdbRun::default();
    run.note_scope_read();
    prefetch(&db, &run, &parsed).unwrap();
    assert_eq!(read(&db, &run, 1), BlobRead::Again);
    assert_eq!(read(&db, &run, 2), BlobRead::Kept);
    // Until the scope's reset, which a run that keeps its scope does not make.
    for keeps_scope in [false, true] {
        let mut lookups = VolumeLookups {
            keeps_scope,
            unrecorded_reads: true,
            ..VolumeLookups::default()
        };
        lookups.reset_scope();
        assert_eq!(lookups.unrecorded_reads, keeps_scope);
        lookups.offset_slot(0, 0, &by(0), &[150], 99990).unwrap();
        assert!(lookups.unrecorded_reads);
    }
}

/// `released_after` (re-review of item O4b2-10, R2-1): of the loads of a failed prefetch, those that
/// NCBI's data source releases after one for certain (released, of a smaller key) and at most (not
/// resolved by the scope, of a key not larger or not known).
#[test]
fn the_releases_after_a_load_count_every_load_by_its_key() {
    use FailedRelease::{Maybe, Released, Scope, Unheld};
    let load = |ordinal: u32, order: Option<u64>| FailedLoad {
        ordinal,
        order,
        surely: true,
    };
    let made = [
        load(0, Some(5)),
        load(1, Some(3)),
        load(2, Some(3)),
        load(3, Some(7)),
        load(4, None),
        load(5, Some(1)),
        load(6, Some(5)),
        load(7, Some(2)),
    ];
    let releases = [
        Released, Released, Maybe, Released, Released, Scope, Released, Unheld,
    ];
    // Of the key 5: keys 3, 3 and 2 below it (one released for certain), 5 (by address), the
    // unknown one.
    assert_eq!(released_after(&made, &releases, 0), (1, 5));
    // Of an unknown key: every load that the scope does not resolve, none for certain.
    assert_eq!(released_after(&made, &releases, 4), (0, 6));
    // Of the key 7: keys 5, 5 and 3 released for certain.
    assert_eq!(released_after(&made, &releases, 3), (3, 6));
}

/// Line NP item O4b2-10, review F6: of several lookups of a batch that throw on their ids, the first
/// in NCBI's map ends the loads (`pslot11_t5` with `gb|FFF00005.1|`, ordinal 7, key 5, made to throw
/// beside the text tag of key 50): keys 4 to 1 (ordinals 8 to 11) are loaded before it, keys 11 to 6
/// and `gb|FFF00090.1|` not.
#[test]
fn the_first_of_several_failing_lookups_ends_the_loads() {
    use LoadBeforeStop::{No, Yes};
    let cwd = empty_cwd("pfwmany");
    let (db, parsed) = pfworder_db(&cwd, "pslot11_t5", SeqType::Protein);
    let five = throwing(&parsed[7]);
    let failing = [
        PrefetchSubject::unflipped(12, Some(&*parsed[12])),
        PrefetchSubject::unflipped(7, Some(&five)),
    ];
    let ordinals: Vec<u32> = (0..12).filter(|&ordinal| ordinal != 7).collect();
    let loaded: Vec<(u32, Option<&crate::report::db_subject::ParsedSubject>)> = ordinals
        .iter()
        .map(|&ordinal| (ordinal, Some(&*parsed[ordinal as usize])))
        .collect();
    let want: Vec<LoadBeforeStop> = ordinals
        .iter()
        .map(|&ordinal| if ordinal >= 8 { Yes } else { No })
        .collect();
    assert_eq!(loads_before_failure(&db, &loaded, &failing, None), want);
    // In either order of the batch: the first lookup of 99990 writes the offset error.
    for reversed in [false, true] {
        let mut subjects: Vec<PrefetchSubject<'_>> = loaded
            .iter()
            .map(|&(ordinal, parsed)| PrefetchSubject::unflipped(ordinal, parsed))
            .chain(failing)
            .collect();
        if reversed {
            subjects.reverse();
        }
        let run = TaxdbRun::default();
        let posts = batch_prefetch(
            &db,
            &run,
            std::iter::empty(),
            subjects,
            "BLASTP",
            &stop_queries(0),
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(posts.text()).unwrap(),
            "Error: [blastp] Error: Offset error at end of taxdb file.\nWarning: [blastp] Error pre-fetching sequence data\n"
        );
        // Four loads, all kept by the data source's cache; those after the failure never loaded.
        let lookups = run.volume_lookups(&db);
        for ordinal in 8..12 {
            assert_eq!(
                lookups.blob_read(
                    ordinal,
                    &blob_name(&db, ordinal, ordinal, "BLASTP").unwrap()
                ),
                BlobRead::Kept
            );
        }
        assert!(!lookups.slot_blobs.contains_key(&1));
    }
}

/// Line NP item O4b2-10, review F3 and F6: NCBI's `GetBioseqHandles` counts the distinct Seq-id
/// handles of every Seq-align id of a batch (`CSortedSeq_ids`) before it asks the scope for any, so a
/// subject that a query resolves (which the prefetch's loads skip) counts: 200 ids in one map, 201 in
/// pieces of 100 sorted otherwise (`MANY_IDS`), whose loads' places to a failed lookup LOSAT does not
/// know; nor where a Seq-align id has no handle that LOSAT knows (`UNNAMED_MAP`).
#[test]
fn the_map_of_a_batch_counts_every_seq_align_id_of_it() {
    use LoadBeforeStop::{No, Unknown};
    let cwd = empty_cwd("pfwmap");
    let search = SearchPath::new(&cwd, None, None);
    let list: Vec<String> = ["ldp", "ldu"]
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
    let db = SeqDb::open_in(list.join(" ").as_bytes(), SeqType::Protein, &search).unwrap();
    let headers = HeaderReader::new(&db, false);
    let stop = (0..db.num_seqs())
        .find_map(
            |ordinal| match read_header(&headers, ordinal, 0, "BLASTP").unwrap() {
                ReportedHeader::Parsed(subject)
                    if String::from_utf8_lossy(&subject.alignment_fasta())
                        .contains("NP_000002.1") =>
                {
                    Some((ordinal, subject))
                }
                _ => None,
            },
        )
        .expect("ref|NP_000002.1|");
    let throws = throwing(&stop.1);
    assert!(throws.prefetch_fails());
    let failing = PrefetchSubject::unflipped(stop.0, Some(&throws));
    // `gnl|BL_ORD_ID|3` and on: packed 3 and more, after the failing packed 2.
    let unparsed: Vec<PrefetchSubject<'_>> = (3..204)
        .map(|ordinal| PrefetchSubject::unflipped(ordinal, None))
        .collect();
    let batch = |subjects: &[PrefetchSubject<'_>]| {
        let mut all = subjects.to_vec();
        all.push(failing);
        batch_map(&db, &all)
    };
    assert_eq!(batch(&unparsed[..199]), None);
    assert_eq!(batch(&unparsed[..200]), Some(MANY_IDS));
    // One handle twice counts once.
    let mut twice = unparsed[..199].to_vec();
    twice.push(unparsed[0]);
    assert_eq!(batch(&twice), None);
    // The prefetch's loads of 199 of them: after the failure where the map is known; not known
    // where the batch has one more id, which its query resolves.
    let loaded: Vec<(u32, Option<&crate::report::db_subject::ParsedSubject>)> = unparsed[..199]
        .iter()
        .map(|subject| (subject.ordinal, None))
        .collect();
    assert!(
        loads_before_failure(&db, &loaded, &[failing], batch(&unparsed[..199]))
            .iter()
            .all(|before| *before == No)
    );
    for map in [batch(&unparsed[..200]), Some(UNNAMED_MAP)] {
        assert!(loads_before_failure(&db, &loaded, &[failing], map)
            .iter()
            .all(|before| *before == Unknown(map.unwrap())));
    }
    assert_eq!(
        map_unknown(&[(2, Some(vec![1])), (3, None)], 1),
        Some(UNNAMED_MAP)
    );
    assert_eq!(map_unknown(&[(2, Some(vec![1]))], 201), Some(MANY_IDS));
    assert_eq!(map_unknown(&[(2, Some(vec![1]))], 200), None);
}

/// Line NP item O4b2-10, review F4 and F6: NCBI posts its warning of a UTF8String header (803.10) once
/// per process (`ERR_POST_X_ONCE`), at the first load that reads one. A load that the prefetch may
/// make before a failed lookup or not at all is refused where it would be that load, not where a
/// load that NCBI surely makes reads one with a smaller key of a map that LOSAT knows (`ntagc_t5`:
/// `gb|FFF00010.1|` key 10 before the failing text tag, `gb|FFF00090.1|` key 90).
#[test]
fn a_load_that_may_not_be_made_is_refused_where_it_would_post_the_utf8_warning() {
    let cwd = empty_cwd("pfwutf8");
    let (db, parsed) = pfworder_db(&cwd, "ntagc_t5", SeqType::Nucleotide);
    db.note_scope_utf8();
    assert!(db.scope_utf8_loads());
    let utf8 = |ordinal: usize| {
        (*parsed[ordinal])
            .clone()
            .scope_utf8(Some(crate::seqdb::DbSequence::of_unrestricted(
                ordinal as u32,
            )))
    };
    let (ten, ninety) = (utf8(0), utf8(2));
    let stop = (1, "gnl|BL_ORD_ID|x050".to_owned());
    let loaded = [(0, Some(&ten)), (2, Some(&ninety))];
    let check = |unknown: [Option<&'static str>; 2], map: Option<&'static str>| {
        unknown_utf8_load(&db, &stop, &loaded, &unknown, map, "BLASTN")
            .map_err(|refusal| refusal.to_string())
    };
    // The load of key 10 surely posts it first.
    assert!(check([None, Some(ADDRESS_ORDER)], None).is_ok());
    for (unknown, map, sequence) in [
        ([Some(ADDRESS_ORDER), None], None, 0),
        ([Some(ADDRESS_ORDER), Some(ADDRESS_ORDER)], None, 0),
        ([None, Some(MANY_IDS)], Some(MANY_IDS), 2),
    ] {
        let refusal = check(unknown, map).unwrap_err();
        assert!(
            refusal.contains(&format!(
                "database sequence {sequence}, which the prefetch may load before or not at all"
            )) && refusal.contains("reads by its lookup a header with a UTF8String"),
            "{refusal}"
        );
    }
    // Where the load that NCBI surely makes reads no such header, the other one posts the warning
    // where NCBI makes it; a load that reads none writes none.
    let plain = [(0, Some(&ten)), (2, Some(&*parsed[2]))];
    assert!(unknown_utf8_load(
        &db,
        &stop,
        &plain,
        &[Some(ADDRESS_ORDER), None],
        None,
        "BLASTN"
    )
    .is_err());
    assert!(unknown_utf8_load(
        &db,
        &stop,
        &plain,
        &[None, Some(ADDRESS_ORDER)],
        None,
        "BLASTN"
    )
    .is_ok());
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

/// `map_unknown`: NCBI's map of a batch's Seq-id handles is unknown with more than 200 distinct
/// ids, and with any id whose handle LOSAT does not know, whatever the count (NCBI makes one
/// handle for each of its Seq-aligns; re-review of item O4b2-7, R2-H1).
#[test]
fn the_map_of_a_batch_with_an_unknown_handle_is_unknown() {
    let named = |n: u64| (n, Some(n.to_le_bytes().to_vec()));
    let keys: Vec<_> = (0..200).map(named).collect();
    assert_eq!(map_unknown(&keys, 200), None);
    let keys: Vec<_> = (0..201).map(named).collect();
    assert_eq!(map_unknown(&keys, 201), Some(MANY_IDS));
    let keys = [named(1), (u64::MAX, None)];
    assert_eq!(map_unknown(&keys, 1), Some(UNNAMED_MAP));
}

/// Scope lifetime across batches (`VolumeLookups::reset_scope`, `blob_read`, `keeps_scope`):
/// - after a batch reset (non-outfmt 5), another handle previously held by a Bioseq (`holds` /
///   `found_in`) is no longer in live scope and loads `Fresh` (fixes R2-L1);
/// - in the same live scope, a report read by an exact held handle hits `Kept` without recording a
///   synthetic blob, a report read by a scope/weak key hits `Again`, and a prefetch by a held
///   handle hits `Again`;
/// - when `kept` is true (<= 10 seqs), reading the same handle after a reset re-attaches its blob
///   into the new live scope (`in_scope = true`).
fn scope_slot_load(ordinal: u32, blob: u32, order: u64, name: BlobName, pins: bool) -> SlotLoad {
    SlotLoad {
        ordinal,
        blob,
        volume: 0,
        lookups: vec![99990, 150],
        order,
        name,
        pins,
        maybe: false,
    }
}

fn prefetch_one(
    lookups: &mut VolumeLookups,
    ordinal: u32,
    blob: u32,
    order: u64,
    name: BlobName,
    pins: bool,
) -> std::result::Result<rustc_hash::FxHashMap<u32, Vec<bool>>, (usize, SlotUncertain)> {
    lookups.prefetch_slots(
        &[scope_slot_load(ordinal, blob, order, name, pins)],
        |_| false,
        99990,
        None,
    )
}

#[test]
fn a_batch_reset_clears_live_scope_matches_unless_keeps_scope() {
    const L: i32 = 99990;
    let mut record = VolumeLookups {
        kept: true,
        in_scope_kept: true,
        ..VolumeLookups::default()
    };
    let b1 = prefetch_one(&mut record, 0, 0, 1, holding(0, &[3]), true).unwrap();
    assert_eq!(b1[&0], vec![true, false]);
    let after_b1 = (record.offset_slots.clone(), record.slot_blobs.clone());
    // Same live scope: exact handle -> Kept on report read; scope key -> Again; prefetch -> Again.
    assert_eq!(record.blob_read(0, &by(3)), BlobRead::Kept);
    assert_eq!(record.blob_read(0, &found_in(4, 0)), BlobRead::Again);
    assert!(prefetch_one(&mut record, 1, 0, 2, by(3), true).is_err());
    // Reset scope (non-outfmt 5): batch 2 prefetch by id 3 loads Fresh -> second offset error.
    (record.offset_slots, record.slot_blobs) = after_b1;
    record.reset_scope();
    assert_eq!(record.blob_read(0, &by(3)), BlobRead::Fresh);
    assert_eq!(record.blob_read(0, &found_in(4, 0)), BlobRead::Fresh);
    let b2 = prefetch_one(&mut record, 1, 0, 2, by(3), true).unwrap();
    assert_eq!(b2[&1], vec![true, false]);
    // Re-attachment: after reset, reading id 0 (Kept via DS cache) re-attaches it to live scope,
    // so a same-batch report read of id 3 hits Kept again.
    let mut reattach = VolumeLookups {
        kept: true,
        in_scope_kept: true,
        ..VolumeLookups::default()
    };
    assert_eq!(
        reattach.offset_slot(0, 0, &holding(0, &[3]), &[L, 150], L),
        Ok(vec![true, false])
    );
    reattach.reset_scope();
    assert_eq!(reattach.blob_read(0, &by(3)), BlobRead::Fresh);
    assert_eq!(reattach.blob_read(0, &by(0)), BlobRead::Kept);
    assert_eq!(
        reattach.offset_slot(0, 0, &by(0), &[L, 150], L),
        Ok(vec![false, false])
    );
    assert_eq!(reattach.blob_read(0, &by(3)), BlobRead::Kept);
    // With `keeps_scope = true` (outfmt 5), `reset_scope` leaves `in_scope` and `pinned` intact.
    let mut xml = VolumeLookups {
        keeps_scope: true,
        ..VolumeLookups::default()
    };
    let m1 = prefetch_one(&mut xml, 0, 0, 1, holding(0, &[3]), true).unwrap();
    assert_eq!(m1[&0], vec![true, false]);
    xml.reset_scope();
    assert_eq!(xml.blob_read(0, &by(0)), BlobRead::Kept);
    assert_eq!(xml.blob_read(0, &by(3)), BlobRead::Kept);
    let m2 = prefetch_one(&mut xml, 0, 0, 1, by(0), true).unwrap();
    assert_eq!(m2[&0], vec![false, false]);
}

/// `outfmt 5` pinning vs non-pinning prefetches, 11-20 `in_scope_kept` lifetime,
/// `TaxdbRun::restart` preserving `keeps_scope`, and `FormatThrow::keeps_scope`.
#[test]
fn outfmt_five_pins_prefetched_blobs_across_batches_and_restart_preserves_keeps_scope() {
    use crate::report::db_subject::FormatThrow;
    const L: i32 = 99990;
    // >20 DB / OBJMGR_* (`kept = false, in_scope_kept = false`), outfmt 5 (`keeps_scope = true`):
    // a failing/stopping prefetch (`pins = false`) does NOT pin; a succeeding one (`pins = true`)
    // pins across `reset_scope`.
    let mut unpinned = VolumeLookups {
        keeps_scope: true,
        ..VolumeLookups::default()
    };
    let m = prefetch_one(&mut unpinned, 0, 0, 1, by(0), false).unwrap();
    assert_eq!(m[&0], vec![true, false]);
    unpinned.reset_scope();
    assert_eq!(unpinned.blob_read(0, &by(0)), BlobRead::Again);
    assert!(prefetch_one(&mut unpinned, 0, 0, 1, by(0), true).is_err());
    // 11-20 DB (`kept = false, in_scope_kept = true`): report reads survive `reset_scope` only
    // when `keeps_scope = true`.
    for keeps_scope in [false, true] {
        let mut record = VolumeLookups {
            in_scope_kept: true,
            keeps_scope,
            ..VolumeLookups::default()
        };
        assert_eq!(
            record.offset_slot(0, 0, &by(0), &[L, 150], L),
            Ok(vec![true, false])
        );
        record.reset_scope();
        assert_eq!(
            record.blob_read(0, &by(0)),
            if keeps_scope {
                BlobRead::Kept
            } else {
                BlobRead::Again
            }
        );
    }
    // `TaxdbRun::restart` preserves `keeps_scope`.
    let mut run = TaxdbRun::default();
    run.set_keeps_scope(true);
    run.restart();
    assert!(run.volume_lookups.lock().unwrap().keeps_scope);
    // `FormatThrow::outfmt_keeps_scope`: true for outfmt 5 alone (derived directly from `ReportFormat::Xml`).
    assert!(FormatThrow::outfmt_keeps_scope("5"));
    assert!(FormatThrow::outfmt_keeps_scope("5 "));
    for outfmt in [
        "0",
        "6",
        "6 qseqid evalue",
        "7",
        "8",
        "11",
        "15",
        "16",
        "17",
        "18",
        "invalid",
    ] {
        assert!(!FormatThrow::outfmt_keeps_scope(outfmt), "{outfmt}");
    }
}

/// Audit invariants for scope lifetime:
/// - `blob_read` (`!in_prefetch`) requires every `in_scope` exact-handle holder to be `blob_kept`
///   (`x_FindBestTSEInIndex` picks the earliest matching TSE);
/// - `unordered_slot` re-attaches a `Kept` blob to live scope via `touch_same_handle`;
/// - `TaxdbRun::reset_scope` fast path vs initialized `offset_error_taxid`, and `restart`
///   preserving `keeps_scope`;
/// - `FormatThrow::scope_flags` and `FormatThrow::is_xml2`;
/// - `prefetch_slot_plan` sets `pinned` on `tm_gi21` (`keeps_scope = true`) when every segment
///   resolves to a non-empty `CBioseq_Handle`.
#[test]
fn scope_lifetime_audit_invariants_hold() {
    use crate::report::db_subject::FormatThrow as Throw;
    const L: i32 = 99990;

    // 1. Mixed `in_scope` blobs holding handle 3 on a >20 DB (`kept = false, in_scope_kept = false`,
    //    `keeps_scope = true`): an earlier unpinned blob (`pins = false`) followed by a later
    //    pinned blob (`pins = true`) must give `Again`, not `Kept`.
    let mut mixed = VolumeLookups {
        keeps_scope: true,
        ..VolumeLookups::default()
    };
    prefetch_one(&mut mixed, 0, 0, 1, holding(0, &[3]), false).unwrap();
    prefetch_one(&mut mixed, 0, 0, 1, holding(1, &[3]), true).unwrap();
    assert_eq!(mixed.blob_read(0, &by(3)), BlobRead::Again);

    // When both `in_scope` holders are pinned, `blob_read` gives `Kept`.
    let mut both_pinned = VolumeLookups {
        keeps_scope: true,
        ..VolumeLookups::default()
    };
    prefetch_one(&mut both_pinned, 0, 0, 1, holding(0, &[3]), true).unwrap();
    prefetch_one(&mut both_pinned, 0, 0, 1, holding(1, &[3]), true).unwrap();
    assert_eq!(both_pinned.blob_read(0, &by(3)), BlobRead::Kept);

    // When the earlier unpinned blob is cleared by `reset_scope` (`keeps_scope = false`), only the
    // new pinned `in_scope` blob matches, giving `Kept`.
    let mut after_reset = VolumeLookups::default();
    prefetch_one(&mut after_reset, 0, 0, 1, holding(0, &[3]), false).unwrap();
    after_reset.reset_scope();
    prefetch_one(&mut after_reset, 0, 0, 1, holding(1, &[3]), true).unwrap();
    assert_eq!(after_reset.blob_read(0, &by(3)), BlobRead::Kept);

    // 2. `unordered_slot` re-attaches a `Kept` blob to live scope after `reset_scope`.
    let mut unordered = VolumeLookups {
        kept: true,
        in_scope_kept: true,
        ..VolumeLookups::default()
    };
    prefetch_one(&mut unordered, 0, 0, 1, holding(0, &[3]), true).unwrap();
    unordered.reset_scope();
    assert_eq!(unordered.blob_read(0, &by(3)), BlobRead::Fresh);
    // Two loads with the same `order` key enter `unordered_slot`: blob 0 (`by(0)`) is `Kept` and
    // must re-attach (`in_scope = true`), while blob 1 (`by(1)`) is `Fresh`.
    let reattached = unordered
        .prefetch_slots(
            &[
                scope_slot_load(0, 0, 5, by(0), true),
                SlotLoad {
                    ordinal: 1,
                    blob: 1,
                    volume: 0,
                    lookups: vec![150],
                    order: 5,
                    name: by(1),
                    pins: true,
                    maybe: false,
                },
            ],
            |_| false,
            L,
            None,
        )
        .unwrap();
    assert_eq!(reattached[&0], vec![false, false]);
    assert_eq!(unordered.blob_read(0, &by(3)), BlobRead::Kept);

    // 3. `TaxdbRun::reset_scope` fast path and `restart` preserving `keeps_scope`.
    let mut kept_run = TaxdbRun::default();
    kept_run.set_keeps_scope(true);
    kept_run.set_separate_xml2(true);
    kept_run.note_calls_before(true);
    kept_run
        .volume_lookups
        .lock()
        .unwrap()
        .offset_slot(0, 0, &by(0), &[L, 150], L)
        .unwrap();
    assert!(kept_run.unrecorded_calls());
    kept_run.restart();
    assert!(kept_run.volume_lookups.lock().unwrap().keeps_scope);
    assert!(kept_run
        .volume_lookups
        .lock()
        .unwrap()
        .slot_blobs
        .is_empty());

    // Fast path: without `texts` initialized, `reset_scope` is a no-op; once `offset_error_taxid`
    // is `Some(L)`, `reset_scope` clears `in_scope` when `keeps_scope = false`.
    let fast = TaxdbRun::default();
    fast.volume_lookups
        .lock()
        .unwrap()
        .offset_slot(0, 0, &by(0), &[L], L)
        .unwrap();
    fast.reset_scope();
    assert!(fast.volume_lookups.lock().unwrap().slot_blobs[&0][0].in_scope);
    let _ = fast.texts.set(PrefetchTexts {
        offset_error_taxid: Some(L),
        ..PrefetchTexts::default()
    });
    fast.reset_scope();
    assert!(!fast.volume_lookups.lock().unwrap().slot_blobs[&0][0].in_scope);

    // 4. `FormatThrow::scope_flags` and `FormatThrow::is_xml2`.
    assert_eq!(Throw::scope_flags(Some("5")), (true, true));
    assert_eq!(Throw::scope_flags(Some("5 ")), (true, true));
    assert_eq!(Throw::scope_flags(Some("6")), (false, true));
    assert_eq!(Throw::scope_flags(Some("15")), (false, false));
    assert_eq!(Throw::scope_flags(None::<&str>), (false, false));
    for outfmt in ["13", "14", "15", "16", "15 "] {
        assert!(Throw::is_xml2(outfmt), "{outfmt}");
    }
    for outfmt in ["0", "5", "6", "7", "8", "11", "17", "18"] {
        assert!(!Throw::is_xml2(outfmt), "{outfmt}");
    }

    // 5. `prefetch_slot_plan` pins `tm_gi21` blobs under `keeps_scope = true` across `reset_scope`.
    let cwd = empty_cwd("pin_gi21");
    let taxdb = variant("offset_last_past_end");
    let db = db_search("tm_gi21", Some(&taxdb), &cwd);
    let parsed = parsed_subjects(&db);
    let mut xml_run = TaxdbRun::default();
    xml_run.set_keeps_scope(true);
    let run_batch = |query: usize| {
        batch_prefetch(
            &db,
            &xml_run,
            std::iter::empty(),
            [0u32, 1]
                .into_iter()
                .map(|o| PrefetchSubject::unflipped(o, Some(&parsed[o as usize]))),
            "BLASTP",
            &stop_queries(query),
        )
        .unwrap()
    };
    assert!(!run_batch(0).text().is_empty());
    xml_run.reset_scope();
    {
        let lookups = xml_run.volume_lookups.lock().unwrap();
        assert!(lookups.slot_blobs[&0][0].pinned);
        assert!(lookups.slot_blobs[&0][0].in_scope);
    }
    assert!(run_batch(1).text().is_empty());
}

/// Secondary XML2/JSON (`separate_xml2 = true`) shares the run's taxonomy handle, open state,
/// non-collision volume lookups, and `calls_before`, while refusing cache-colliding volumes
/// (`refuse_cache_collision`).
#[test]
fn separate_xml2_shares_open_and_non_collision_state_and_refuses_collision() {
    let cwd = empty_cwd("sep_xml2_shared");
    // 1. Opening-texts taxdb (`bad_magic`): primary prefetch opens `lazy` once, so secondary
    //    `read_reported_taxonomy` gets empty `notes.open` (opening texts are not repeated).
    let db_open = db_search("tm_cache151", Some(&variant("bad_magic")), &cwd);
    let parsed_open = parsed_subjects(&db_open);
    let mut run_open = TaxdbRun::default();
    run_open.set_separate_xml2(true);
    run_open.note_calls_before(true);
    let first_open = batch_prefetch(
        &db_open,
        &run_open,
        std::iter::empty(),
        [PrefetchSubject::unflipped(0, Some(&parsed_open[0]))],
        "BLASTP",
        &stop_queries(0),
    )
    .unwrap();
    assert!(!first_open.text().is_empty());
    let ordinals: Vec<u32> = (0..db_open.num_seqs()).collect();
    let reported: std::collections::BTreeSet<usize> = [0].into_iter().collect();
    let wrapped_open: Vec<_> = parsed_open.iter().map(Some).collect();
    let (_, _, _, notes_open) = read_reported_taxonomy(
        &db_open,
        &run_open,
        &ordinals,
        &reported,
        &wrapped_open,
        "BLASTP",
    )
    .unwrap();
    assert!(notes_open.open.is_empty());
    assert!(run_open.unrecorded_calls());

    // 2. Non-collision `offset_last_past_end` (`tm_cache151`): primary prefetch records `Last`, so
    //    secondary `read_reported_taxonomy` succeeds and repeats no opening texts.
    let db_nc = db_search("tm_cache151", Some(&variant("offset_last_past_end")), &cwd);
    let parsed_nc = parsed_subjects(&db_nc);
    let mut run_nc = TaxdbRun::default();
    run_nc.set_separate_xml2(true);
    let first_nc = batch_prefetch(
        &db_nc,
        &run_nc,
        std::iter::empty(),
        [PrefetchSubject::unflipped(0, Some(&parsed_nc[0]))],
        "BLASTP",
        &stop_queries(0),
    )
    .unwrap();
    assert!(!first_nc.text().is_empty());
    let wrapped_nc: Vec<_> = parsed_nc.iter().map(Some).collect();
    let (_, _, _, notes_nc) = read_reported_taxonomy(
        &db_nc,
        &run_nc,
        &[0, 1, 2],
        &[0usize].into_iter().collect(),
        &wrapped_nc,
        "BLASTP",
    )
    .unwrap();
    assert!(notes_nc.open.is_empty());

    // 3. Collision `offset_last_past_end` (`tm_cache_gi`): `separate_xml2 = true` refuses
    //    `read_reported_taxonomy` with `refuse_cache_collision` without affecting primary scope.
    let db_col = db_search("tm_cache_gi", Some(&variant("offset_last_past_end")), &cwd);
    let parsed_col = parsed_subjects(&db_col);
    let mut run_col = TaxdbRun::default();
    run_col.set_keeps_scope(true);
    run_col.set_separate_xml2(true);
    let all_ordinals: Vec<u32> = (0..db_col.num_seqs()).collect();
    let all_reported: std::collections::BTreeSet<usize> = (0..db_col.num_seqs() as usize).collect();
    let wrapped_col: Vec<_> = parsed_col.iter().map(Some).collect();
    let primary_batch = |query: usize| {
        batch_prefetch(
            &db_col,
            &run_col,
            std::iter::empty(),
            all_ordinals
                .iter()
                .map(|&o| PrefetchSubject::unflipped(o, Some(&parsed_col[o as usize]))),
            "BLASTP",
            &stop_queries(query),
        )
        .unwrap()
    };
    assert!(!primary_batch(0).text().is_empty());
    run_col.reset_scope();
    let sec_err = read_reported_taxonomy(
        &db_col,
        &run_col,
        &all_ordinals,
        &all_reported,
        &wrapped_col,
        "BLASTP",
    )
    .err()
    .expect("secondary XML2/JSON on colliding volume must refuse");
    assert!(
        sec_err.to_string().contains("256-slot taxonomy cache"),
        "unexpected error: {sec_err}"
    );
    // Primary `volume_lookups` still keeps its scope from batch 1 -> primary batch 2 writes nothing.
    assert!(primary_batch(1).text().is_empty());
}

/// `prefetch_slot_plan` pins only the prefix of segments preceding the first empty
/// `CBioseq_Handle` in `CSeq_id::CompareOrdered` order (`CanResolveRange`, seq_map.cpp:1085-1097),
/// not packed handle order (`load_places`).
#[test]
fn prefetch_slot_plan_pins_only_prefix_before_first_empty_handle_in_compare_ordered_order() {
    use crate::report::db_subject::{FoundInstead, ParsedSubject};
    use crate::seqid::defline_set::{BlastDefline, DeflineSet};
    use crate::seqid::model::{ParsedSeqId, PdbSeqId};
    const L: i32 = 99990;

    // 1. Direct `VolumeLookups::prefetch_slots` with mixed per-load `pins` (including inside
    //    `unordered_slot` when two loads share an order key).
    let mut mixed = VolumeLookups {
        keeps_scope: true,
        ..VolumeLookups::default()
    };
    let loads = [
        scope_slot_load(0, 0, 5, by(0), true),
        SlotLoad {
            ordinal: 1,
            blob: 1,
            volume: 0,
            lookups: vec![150],
            order: 5,
            name: by(1),
            pins: false,
            maybe: false,
        },
    ];
    mixed.prefetch_slots(&loads, |_| false, L, None).unwrap();
    assert!(mixed.slot_blobs[&0][0].pinned);
    assert!(!mixed.slot_blobs[&1][0].pinned);

    // 2. Real `prefetch_slot_plan` on `tm_gi21` (`gi|1` at ordinal 0, `gi|5` at ordinal 4, `gi|9`
    //    at ordinal 8; all with taxids in slot 150).
    let cwd = empty_cwd("partial_pin");
    let db = db_search("tm_gi21", Some(&variant("offset_last_past_end")), &cwd);
    let parsed = parsed_subjects(&db);
    let pinned_blob = |run: &TaxdbRun, ordinal: u32| {
        run.volume_lookups.lock().unwrap().slot_blobs[&ordinal][0].pinned
    };
    let make_empty = |id: ParsedSeqId| {
        ParsedSubject::new(DeflineSet::from(vec![BlastDefline {
            ids: vec![id],
            ..BlastDefline::default()
        }]))
        .unwrap()
        .found_instead(Some(FoundInstead::Nothing))
    };

    // Case A: empty `gi|5` splits `gi|1` (< `gi|5`, pinned) and `gi|9` (> `gi|5`, unpinned).
    let empty_gi5 = parsed[4].clone().found_instead(Some(FoundInstead::Nothing));
    let mut run_a = TaxdbRun::default();
    run_a.set_keeps_scope(true);
    let loaded_a = [
        (0u32, Some(&parsed[0])),
        (4u32, Some(&empty_gi5)),
        (8u32, Some(&parsed[8])),
    ];
    prefetch_slot_plan(&db, &run_a, &loaded_a, None, true, None, "BLASTP").unwrap();
    assert!(pinned_blob(&run_a, 0));
    assert!(!pinned_blob(&run_a, 8));

    // Case B: empty `bbs|123` (`e_Gibbsq`, unpacked -> `load_places.order` is `u64::MAX`) vs
    // `gi|1` (`e_Gi`, packed): under `CSeq_id::CompareOrdered`, `bbs|123` precedes `gi|1`, so
    // `gi|1` must not be pinned.
    let empty_bbs = make_empty(ParsedSeqId::Gibbsq(123));
    let mut run_b = TaxdbRun::default();
    run_b.set_keeps_scope(true);
    let loaded_b = [(0u32, Some(&parsed[0])), (1u32, Some(&empty_bbs))];
    prefetch_slot_plan(&db, &run_b, &loaded_b, None, true, None, "BLASTP").unwrap();
    assert!(!pinned_blob(&run_b, 0));

    // Case C: empty `pdb|1ABC|A` (`e_Pdb`, unpacked) vs `gi|1` (`e_Gi`): `gi|1` precedes
    // `pdb|1ABC|A` under `CSeq_id::CompareOrdered`, so `gi|1` is pinned.
    let empty_pdb = make_empty(ParsedSeqId::Pdb(PdbSeqId {
        mol: b"1ABC".to_vec(),
        chain: Some(b'A' as i32),
        rel: None,
        chain_id: None,
    }));
    let mut run_c = TaxdbRun::default();
    run_c.set_keeps_scope(true);
    let loaded_c = [(0u32, Some(&parsed[0])), (1u32, Some(&empty_pdb))];
    prefetch_slot_plan(&db, &run_c, &loaded_c, None, true, None, "BLASTP").unwrap();
    assert!(pinned_blob(&run_c, 0));

    // Case D: `Holders::Unknown` (even on a taxid-free subject that `prefetch_reads` does not
    // refuse) causes `prefetch_pin_stops` to return `None`, leaving the whole batch unpinned.
    let unknown = ParsedSubject::new(DeflineSet::from(vec![BlastDefline {
        ids: vec![ParsedSeqId::Gi(99)],
        ..BlastDefline::default()
    }]))
    .unwrap()
    .unrecorded_holders(crate::report::db_subject::ScopeHolders {
        loader: crate::report::db_subject::Holders::Unknown("unreadable volume".to_string()),
        scope: crate::report::db_subject::Holders::default(),
    });
    let mut run_d = TaxdbRun::default();
    run_d.set_keeps_scope(true);
    let loaded_d = [(0u32, Some(&parsed[0])), (4u32, Some(&unknown))];
    assert_eq!(
        prefetch_pin_stops(&db, &loaded_d, true, "BLASTP").unwrap(),
        None
    );
    prefetch_slot_plan(&db, &run_d, &loaded_d, None, true, None, "BLASTP").unwrap();
    assert!(!pinned_blob(&run_d, 0));
}

#[test]
fn one_by_one_query_reads_stops_at_throw_and_collapses_joined_copies() {
    use crate::report::db_subject::{
        AfterThrow, FormatThrow, HitLoads, LoadsBefore, ParsedSubject, RunThrow,
    };
    use crate::seqid::defline_set::{BlastDefline, DeflineSet};
    use crate::seqid::model::ParsedSeqId;

    let cwd = empty_cwd("one_by_one_query_reads");
    let throwing_db = db_search("t5/ppig1", None, &cwd);
    let headers = HeaderReader::new(&throwing_db, false);
    let ReportedHeader::Parsed(throwing) = read_header(&headers, 0, 0, "BLASTP").unwrap() else {
        panic!("t5/ppig1 sequence 0 must be a parsed subject");
    };
    let make_gi = |gi: i64| {
        ParsedSubject::new(DeflineSet::from(vec![BlastDefline {
            ids: vec![ParsedSeqId::Gi(gi)],
            ..BlastDefline::default()
        }]))
        .unwrap()
    };
    // Subjects: 0 = gi|1, 1 = gi|1 (adjacent copy of 0), 2 = throwing (t5/ppig1), 3 = gi|3.
    let subjects = [make_gi(1), make_gi(1), (*throwing).clone(), make_gi(3)];
    let parsed = |rank: usize| Some(&subjects[rank]);

    // Query 0 hits [0, 1, 2(throws), 3]; Query 1 hits [3], so run-wide `loads.reached(3)` is
    // true even though Query 0 stopped at rank 2 under `AfterThrow::EndsReport`.
    let loads = HitLoads::of(
        LoadsBefore::default(),
        FormatThrow::Writes(RunThrow::every_row(AfterThrow::EndsReport)),
        [
            (true, vec![Some(0), Some(1), Some(2), Some(3)]),
            (true, vec![Some(3)]),
        ],
        parsed,
        |_, _| 1,
    );
    assert!(loads.reached(3));
    assert!(loads.thrown(2).is_some());

    // XML2/JSON (`joins_copies = true`, `ends_on_throw = true`): collapses adjacent copy 1 and
    // stops at throwing rank 2 without a phantom read of rank 3.
    assert_eq!(
        one_by_one_query_reads([0, 1, 2, 3], parsed, &loads, true, true),
        vec![0, 2]
    );
    // BLAST XML (`joins_copies = true`, `ends_on_throw = false`): collapses adjacent copy 1 and
    // continues to rank 3.
    assert_eq!(
        one_by_one_query_reads([0, 1, 2, 3], parsed, &loads, true, false),
        vec![0, 2, 3]
    );
    // Tabular (`joins_copies = false`, `ends_on_throw = false`): reads every reached entry.
    assert_eq!(
        one_by_one_query_reads([0, 1, 2, 3], parsed, &loads, false, false),
        vec![0, 1, 2, 3]
    );
}
