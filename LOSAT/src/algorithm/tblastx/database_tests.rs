use super::*;

// DECISIONS 3 Q8 of line DBID: NCBI's choice of the Seq-ids of a database with parsed
// Seq-ids depends on SEQ_ID_PREFER_ACCESSION_OVER_GI and SEQ_ID_AVOID_GI.
#[test]
fn a_seq_id_variable_refuses_a_database_with_parsed_seqids_only() {
    let path = |name: &str| {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(name)
            .to_str()
            .expect("fixture path is UTF-8")
            .to_string()
    };
    let settings = ApplicationSettings {
        seq_id_rank_variable: Some("SEQ_ID_PREFER_ACCESSION_OVER_GI"),
        ..ApplicationSettings::default()
    };
    let parsed = path("tests/db_search/ids/nucl_v5");
    TblastxDatabase::open(&parsed, &ApplicationSettings::default())
        .expect("a database with parsed Seq-ids opens");
    let Err(error) = TblastxDatabase::open(&parsed, &settings) else {
        panic!("a database with parsed Seq-ids opened with SEQ_ID_PREFER_ACCESSION_OVER_GI")
    };
    assert_eq!(
        format!("{error:#}"),
        "the environment variable SEQ_ID_PREFER_ACCESSION_OVER_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's TBLASTX"
    );
    TblastxDatabase::open(&path("tests/seqdb/db/nucl_v4"), &settings)
        .expect("a database without parsed Seq-ids opens with SEQ_ID_PREFER_ACCESSION_OVER_GI");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-911
// An alias file's `LENGTH` of 2^32 letters or more over the 16 sequences of
// `tests/db_search/taxfilter/tfn` (`tests/db_search/dbcounts`, row db.cnt.avg2p32.tx) is
// NCBI's `Int8` on every target, wasm32's 32-bit `usize` included: the database opens with
// its whole total.
#[test]
fn totals_of_2_to_the_32_letters_or_more_are_held_whole() {
    for (name, total) in [
        ("l4295", 4_294_967_296_000_000_000),
        ("l4294", 4_294_967_295_999_999_999),
        ("avg2p32", 68_719_476_736),
    ] {
        assert!(total > i64::from(u32::MAX), "{name}");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/dbcounts")
            .join(name);
        let database =
            TblastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
                .unwrap_or_else(|error| panic!("{name}: {error:#}"));
        let numbers = database.numbers(None);
        assert_eq!(
            (numbers.total, numbers.count, numbers.link_total),
            (total, 16, total),
            "{name}"
        );
    }
}

/// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
/// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> TblastxDatabase {
    use crate::seqdb::alias::SearchPath;
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
    let search =
        SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None).unwrap();
    let path = tests.join("db_search").join(name);
    let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
    TblastxDatabase {
        parsed: db.has_parsed_seqids().unwrap(),
        db,
        name: name.to_string(),
        taxdb: shared::TaxdbRun::default(),
    }
}

/// A reported hit of the batch's first query on the database sequence `ordinal`.
fn hit_on(ordinal: u32) -> TblastxHsp {
    TblastxHsp {
        hit: crate::common::Hit {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 30,
            s_start: 1,
            s_end: 30,
            e_value: 1e-9,
            bit_score: 50.0,
            num_ident: 10,
            query_frame: 1,
            query_length: 0,
            q_idx: 0,
            s_idx: ordinal,
            raw_score: 100,
            sort_query_offset: 0,
            sort_query_end: 10,
            sort_subject_offset: 0,
            sort_subject_end: 10,
            has_sort_offsets: true,
            gap_info: None,
            num_positives: 10,
        },
        subject_frame: 1,
        num: 1,
    }
}

// DBID item T: a batch's prefetch looks up the sequences that its reported hits name, once
// per volume and taxid over the batches of a run, and opens the taxonomy database at the
// run's first lookup; nothing without `taxdb.bti`. `x16_ntaxid_99990` (two sequences of
// taxid 99990) with `offset_last_past_end` (the last entry, 99990, is past the end of
// `taxdb.btd`), `x16_ntaxid` (taxid 9606) with `bad_magic`.
#[test]
fn the_prefetch_looks_up_the_reported_sequences_once_per_run() {
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-tblastx-db-prefetch-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let text = |database: &TblastxDatabase, ordinals: &[u32]| -> String {
        let ordered: Vec<TblastxHsp> = ordinals.iter().map(|&ordinal| hit_on(ordinal)).collect();
        let mut titles = DatabaseTitles::new();
        read_titles(database, &ordered, 0, &mut titles).unwrap();
        String::from_utf8(prefetch_messages(database, &ordered, &titles).unwrap()).unwrap()
    };
    let database = open_with_taxdb("x16_ntaxid_99990", Some("offset_last_past_end"), &cwd);
    assert!(database.taxdb_found());
    // A batch without reported hits fetches nothing: the taxonomy database is not opened.
    assert_eq!(text(&database, &[]), "");
    // Two hits on sequence 1 and one on sequence 0, all of taxid 99990: one lookup.
    let offset = "Error: [tblastx] Error: Offset error at end of taxdb file.\n";
    assert_eq!(text(&database, &[1, 1, 0]), offset);
    // The same volume and taxid in a later batch.
    assert_eq!(text(&database, &[0, 1]), "");
    let database = open_with_taxdb("x16_ntaxid", Some("bad_magic"), &cwd);
    assert_eq!(
        text(&database, &[1, 0]),
        "Error: [tblastx] Error: Tax database file has wrong magic number.\n"
    );
    let database = open_with_taxdb("x16_ntaxid", None, &cwd);
    assert!(!database.taxdb_found());
    assert_eq!(text(&database, &[0, 1]), "");
    std::fs::remove_dir_all(&cwd).unwrap();
}

// Line DBID, ID4t fix 2b: the loop writes each query batch's reports when it takes the batch
// only where no report of any sequence of the database can be refused (`report_stream`): a
// sequence that NCBI's formatters do not find by its Seq-align id (nrefuse_v5) refuses the
// tabular formats and outfmt 0, and two sequences whose Seq-align ids match (njoined, in two
// volumes) refuse outfmt 0 only, read on one thread and on four.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn a_report_of_any_sequence_that_could_be_refused_is_found() {
    let open = |name: &str| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids")
            .join(name);
        TblastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    };
    let tabular = ParsedShown {
        tabular: true,
        ..ParsedShown::default()
    };
    let pairwise = ParsedShown {
        pairwise: true,
        joined: true,
        ..ParsedShown::default()
    };
    let answers = |pool: &crate::utils::threading::SearchPool<'_>| {
        ["nucl_v5", "nucl_v4", "nrefuse_v5", "njoined"].map(|name| {
            let database = open(name);
            (
                database.reports_may_refuse(tabular, false, None, pool),
                database.reports_may_refuse(pairwise, false, None, pool),
            )
        })
    };
    use shared::ReportRefusal::{Headers, Joined};
    let serial = answers(&crate::utils::threading::SearchPool::x_serial());
    assert_eq!(
        serial,
        [
            (None, None),
            (None, None),
            (Some(Headers), Some(Headers)),
            (None, Some(Joined)),
        ]
    );
    let parallel =
        crate::utils::threading::with_search_pool(4, "TBLASTX", |pool| Ok(answers(pool))).unwrap();
    assert_eq!(serial, parallel);
}

// Line DBID, DECISIONS 19: a database with parsed Seq-ids is beyond a bound of fewer
// sequences than its volumes hold (`report_stream::Off::LargeParsed`), also where a taxid list
// keeps fewer (the range walks decode the excluded ones of every range they read); a database
// without parsed Seq-ids is never beyond it (read whatever its size, as blastn's).
#[test]
fn a_database_with_parsed_seqids_is_beyond_a_bound_below_its_volumes_sequences() {
    let open = |name: &str| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join(name);
        TblastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    };
    for name in ["db_search/ids/nucl_v5", "db_search/ids/nucl_v4"] {
        let database = open(name);
        let count = database.db.num_seqs();
        assert!(count > 1, "{name}");
        assert!(!database.parsed_beyond(count), "{name}");
        assert!(database.parsed_beyond(count - 1), "{name}");
    }
    let mut restricted = open("db_search/ids/nucl_v5");
    let count = restricted.db.num_seqs();
    restricted
        .db
        .restrict_taxids(&crate::seqdb::TaxidList::Positive(
            [9606].into_iter().collect(),
        ))
        .unwrap();
    let included = restricted.db.num_seqs();
    assert!(included > 0 && included < count);
    assert!(!restricted.parsed_beyond(count));
    assert!(restricted.parsed_beyond(included));
    let stored = open("seqdb/db/nucl_v4");
    assert!(stored.db.num_seqs() > 0);
    assert!(!stored.parsed_beyond(0));
}
