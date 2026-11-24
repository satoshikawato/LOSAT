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
    TblastnDatabase::open(&parsed, &ApplicationSettings::default())
        .expect("a database with parsed Seq-ids opens");
    let Err(error) = TblastnDatabase::open(&parsed, &settings) else {
        panic!("a database with parsed Seq-ids opened with SEQ_ID_PREFER_ACCESSION_OVER_GI")
    };
    assert_eq!(
        format!("{error:#}"),
        "the environment variable SEQ_ID_PREFER_ACCESSION_OVER_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's TBLASTN"
    );
    TblastnDatabase::open(&path("tests/seqdb/db/nucl_v4"), &settings)
        .expect("a database without parsed Seq-ids opens with SEQ_ID_PREFER_ACCESSION_OVER_GI");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-911
// An alias file's `LENGTH` of 2^32 letters or more over the 16 sequences of
// `tests/db_search/taxfilter/tfn` (`tests/db_search/dbcounts`) is NCBI's `Int8` on every
// target, wasm32's 32-bit `usize` included: the database opens with its whole total.
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
            TblastnDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
                .unwrap_or_else(|error| panic!("{name}: {error:#}"));
        let numbers = database.numbers();
        assert_eq!(
            (numbers.total_length, numbers.num_seqs),
            (total, 16),
            "{name}"
        );
    }
}

/// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
/// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> TblastnDatabase {
    use crate::seqdb::alias::SearchPath;
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
    let search = SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None);
    let path = tests.join("db_search").join(name);
    TblastnDatabase {
        db: {
            let mut db = SeqDb::open_in(
                path.to_str().unwrap().as_bytes(),
                SeqType::Nucleotide,
                &search,
            )
            .unwrap();
            db.read_name_records(path.to_str().unwrap(), false).unwrap();
            db
        },
        name: name.to_string(),
        taxdb: shared::TaxdbRun::default(),
        parsed_seqids: true,
    }
}

/// The final hit list of one query whose lists name the sequences `lists` by ordinal id, each
/// with one HSP of the score given (0: no Seq-align, `KappaResultHitList::drop_zero_score_hsps`).
fn hitlist(lists: &[(i32, i32)]) -> KappaResultHitList {
    use super::super::search_gapped::GappedHsp;
    use super::super::stage_d_linking::{LinkedHsp, LinkedHspList};
    use super::super::stage_d_results::KappaResultList;
    let mut hitlist = KappaResultHitList::new(500).unwrap();
    for &(oid, score) in lists {
        let hsp = LinkedHsp {
            context: 0,
            hsp: GappedHsp {
                frame: 1,
                score,
                q_start: 0,
                q_end: 1,
                q_gapped_start: 0,
                s_start: 0,
                s_end: 1,
                s_gapped_start: 0,
            },
            num: 1,
            evalue: 1e-9,
            source_index: 0,
        };
        let list = KappaResultList {
            oid,
            hsps: LinkedHspList {
                hsps: vec![hsp],
                best_evalue: 1e-9,
            },
            payloads: Vec::new(),
        };
        assert!(hitlist.update(list).unwrap().is_none());
    }
    hitlist
}

// Reports that read no letters get the lengths of the reported sequences and no store; the
// others get both, by rank.
#[test]
fn the_reported_sequences_are_read_only_for_reports_of_letters() {
    use super::super::x_subject_store::SubjectLetters;
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-tblastn-db-reported-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let database = open_with_taxdb("x16_ntaxid", None, &cwd);
    let stream = database.x_stream(1, None).unwrap();
    let of = |letters: bool| {
        let mut hitlists = vec![hitlist(&[(1, 50), (0, 40)])];
        let mut titles = DatabaseTitles::new();
        read_titles(
            &database,
            &hitlists,
            0,
            (Default::default(), true),
            &mut titles,
        )
        .unwrap();
        reported_subjects(
            &database,
            &titles,
            &mut hitlists,
            SubjectLetters::Stream(&stream),
            letters,
            &[],
            (&[], false),
            None,
            &[],
            true,
        )
        .unwrap()
    };
    let whole = [
        stream.read_length(0).unwrap(),
        stream.read_length(1).unwrap(),
    ];
    let reported = of(false);
    assert_eq!(reported.lengths, whole);
    assert!(reported.stores.is_empty());
    let letters = reported.letters();
    assert!(matches!(letters, SubjectLetters::Lengths(_)));
    assert_eq!(letters.count(), 2);
    assert_eq!(letters.length(1), whole[1]);
    let reported = of(true);
    assert_eq!(reported.lengths, whole);
    assert_eq!(reported.stores.len(), 2);
    assert_eq!(reported.stores[0].len(), whole[0]);
    assert_eq!(reported.stores[1].len(), whole[1]);
    let letters = reported.letters();
    assert!(matches!(letters, SubjectLetters::Packed(_)));
    assert_eq!(letters.length(1), whole[1]);
    std::fs::remove_dir_all(&cwd).unwrap();
}

// DBID item T: a batch's prefetch looks up the sequences that its reports show (an HSP with a
// Seq-align), once per volume and taxid over the batches of a run, and opens the taxonomy
// database at the run's first lookup; nothing without `taxdb.bti`. `x16_ntaxid_99990` (two
// sequences of taxid 99990) with `offset_last_past_end` (the last entry, 99990, is past the
// end of `taxdb.btd`), `x16_ntaxid` (taxid 9606) with `bad_magic`.
#[test]
fn the_prefetch_looks_up_the_reported_sequences_once_per_run() {
    use super::super::x_subject_store::SubjectLetters;
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-tblastn-db-prefetch-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    // The reports of one batch of `queries` (the lists of each query), as the run makes them:
    // the titles, the reported sequences by rank, then the batch's prefetch.
    let text = |database: &TblastnDatabase, queries: &[&[(i32, i32)]]| -> String {
        let stream = database.x_stream(1, None).unwrap();
        let mut hitlists: Vec<KappaResultHitList> =
            queries.iter().map(|&lists| hitlist(lists)).collect();
        let mut titles = DatabaseTitles::new();
        read_titles(
            database,
            &hitlists,
            0,
            (Default::default(), true),
            &mut titles,
        )
        .unwrap();
        let reported = reported_subjects(
            database,
            &titles,
            &mut hitlists,
            SubjectLetters::Stream(&stream),
            true,
            &[],
            (&[], false),
            None,
            &[],
            true,
        )
        .unwrap();
        let messages = database.prefetch_messages(&reported, &hitlists, SubjectFetch::Every);
        String::from_utf8(messages.unwrap()).unwrap()
    };
    let database = open_with_taxdb("x16_ntaxid_99990", Some("offset_last_past_end"), &cwd);
    assert!(database.taxdb_found());
    // A sequence whose HSPs have no Seq-align is not fetched: the taxonomy database is not
    // opened.
    assert_eq!(text(&database, &[&[(0, 0)]]), "");
    let offset = "Error: [tblastn] Error: Offset error at end of taxdb file.\n";
    assert_eq!(text(&database, &[&[(0, 0)], &[(1, 50)]]), offset);
    // The same volume and taxid in a later batch.
    assert_eq!(text(&database, &[&[(0, 50), (1, 50)]]), "");
    let database = open_with_taxdb("x16_ntaxid", Some("bad_magic"), &cwd);
    assert_eq!(
        text(&database, &[&[(1, 50), (0, 40)]]),
        "Error: [tblastn] Error: Tax database file has wrong magic number.\n"
    );
    let database = open_with_taxdb("x16_ntaxid", None, &cwd);
    assert!(!database.taxdb_found());
    assert_eq!(text(&database, &[&[(0, 50), (1, 40)]]), "");
    std::fs::remove_dir_all(&cwd).unwrap();
}
