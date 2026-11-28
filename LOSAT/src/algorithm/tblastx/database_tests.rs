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
    TblastxDatabase::open(path("tests/seqdb/db/nucl_v4"), &settings)
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
    let search = SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None);
    let path = tests.join("db_search").join(name);
    let db = SeqDb::open_in(
        path.to_str().unwrap().as_bytes(),
        SeqType::Nucleotide,
        &search,
    )
    .unwrap();
    TblastxDatabase {
        parsed: db.has_parsed_seqids().unwrap(),
        db,
        name: name.as_bytes().to_vec(),
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
        read_titles(database, &ordered, 0, Default::default(), &mut titles).unwrap();
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
// volumes) refuse outfmt 0 only; where a batch may read its subjects one by one, the scope's
// index can give either for the other's id (a weak match; line DBID hr-1, the keys of NCBI's
// places, `ParsedSubject::place_keys`), which refuses both. Read on one thread and on four.
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
        [
            ("nucl_v5", false),
            ("nucl_v4", false),
            ("nrefuse_v5", false),
            ("njoined", false),
            ("njoined", true),
        ]
        .map(|(name, one_by_one)| {
            let database = open(name);
            (
                database.reports_may_refuse(tabular, (false, one_by_one), None, pool),
                database.reports_may_refuse(pairwise, (false, one_by_one), None, pool),
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
            (Some(Headers), Some(Headers)),
        ]
    );
    let parallel =
        crate::utils::threading::with_search_pool(4, "TBLASTX", |pool| Ok(answers(pool))).unwrap();
    assert_eq!(serial, parallel);
}

// Review of D1c-1c-b PR 2b: a database sequence with a definition line without a Seq-id
// (`tests/db_search/idlist/pn`) only refuses tabular formats that print `sallacc`
// (`stops_at_lineless`), not default or `sacc`/`staxid` field lists (`stops_at_failed_lookup`).
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn a_lineless_defline_only_refuses_report_streaming_when_sallacc_is_shown() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/idlist/pn");
    let database = TblastxDatabase::open(path.to_str().unwrap(), &ApplicationSettings::default())
        .unwrap_or_else(|error| panic!("pn: {error:#}"));
    let without_sallacc = ParsedShown {
        tabular: true,
        stops_at_failed_lookup: true,
        ..ParsedShown::default()
    };
    let with_sallacc = ParsedShown {
        tabular: true,
        stops_at_lineless: true,
        ..ParsedShown::default()
    };
    let answers = |pool: &crate::utils::threading::SearchPool<'_>| {
        (
            database.reports_may_refuse(without_sallacc, (false, false), None, pool),
            database.reports_may_refuse(with_sallacc, (false, false), None, pool),
        )
    };
    let serial = answers(&crate::utils::threading::SearchPool::x_serial());
    assert_eq!(serial, (None, Some(shared::ReportRefusal::Headers)));
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

// Review 2 of line DBID pfoff-3, L5: TBLASTX's model of what one query's pairwise report reads
// (`pairwise_read_order` of the query's hits in report order) is the writer's group of the
// rows that `report::pairwise_hits` makes of those hits (`x_pairwise_par::writer_reads`), for
// every sort and count.
#[test]
fn the_pairwise_reads_are_those_of_the_writers_rows() {
    use crate::algorithm::tblastx::report::{pairwise_hits, seqalign_identities, SubjectLetters};
    use crate::blastinput::fasta_reader::FastaRecord;
    use crate::blastinput::seq_range::Placements;
    use crate::utils::genetic_code::GeneticCode;
    let hsp =
        |s_idx: u32, raw_score: i32, e_value: f64, offset: usize, num_ident: usize| TblastxHsp {
            hit: crate::common::Hit {
                identity: 0.0,
                length: 10,
                mismatch: 0,
                gapopen: 0,
                q_start: 1,
                q_end: 30,
                s_start: 1,
                s_end: 30,
                e_value,
                bit_score: raw_score as f64 / 2.0,
                num_ident,
                query_frame: 1,
                query_length: 0,
                q_idx: 0,
                s_idx,
                raw_score,
                sort_query_offset: offset,
                sort_query_end: offset + 10,
                sort_subject_offset: offset,
                sort_subject_end: offset + 10,
                has_sort_offsets: true,
                gap_info: None,
                num_positives: 0,
            },
            subject_frame: 1,
            num: 1,
        };
    let ordered = [
        hsp(0, 50, 1e-9, 4, 3),
        hsp(0, 30, 1e-6, 0, 9),
        hsp(1, 60, 1e-12, 2, 1),
        hsp(2, 20, 1e-3, 1, 8),
        hsp(2, 40, 1e-7, 6, 5),
    ];
    // Long enough for the residues of every HSP (offsets up to 16 residues).
    let letters = b"ATGGCTAAACTGGAAGGTCCATTTAGCCGAATGGCTAAACTGGAAGGTCCATTTAGCCGAATGGCT";
    let queries = [FastaRecord::new("q", b"q", letters)];
    let subjects: Vec<FastaRecord> = (0..3)
        .map(|_| FastaRecord::new("s", b"s", letters))
        .collect();
    let code = GeneticCode::from_id(1);
    let placements = Placements::default();
    let identities = seqalign_identities(&ordered);
    let rows = pairwise_hits(
        &ordered,
        Some(&identities),
        &queries,
        &subjects,
        SubjectLetters::Records(&subjects),
        &[super::super::blast_engine::TblastxQueryStats::default()],
        &code,
        &code,
        &placements,
        &placements,
        Default::default,
        None,
    );
    let titles = DatabaseTitles::new();
    for hits in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
        for hsps in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
            for (descriptions, alignments) in [(0, 0), (1, 2), (2, 1), (500, 500)] {
                let reading = crate::report::pairwise::PairwiseReading {
                    sort: crate::report::hit_sort::ReportSort { hits, hsps },
                    shown: crate::report::db_subject::Shown {
                        descriptions,
                        alignments,
                    },
                };
                assert_eq!(
                    pairwise_read_order(&ordered, &titles, reading),
                    crate::report::x_pairwise_par::writer_reads(
                        &rows,
                        |_, _| false,
                        reading,
                        true,
                        None
                    ),
                    "{reading:?}"
                );
            }
        }
    }
}
