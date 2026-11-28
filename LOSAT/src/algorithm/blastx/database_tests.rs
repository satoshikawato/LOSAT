use super::*;

fn fixture(name: &str) -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/db")
        .join(name)
        .to_str()
        .expect("fixture path is UTF-8")
        .to_string()
}

fn open(name: &str) -> BlastxDatabase {
    BlastxDatabase::open(fixture(name), &ApplicationSettings::default())
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
}

// DECISIONS 3 Q8 of line DBID: NCBI's choice of the Seq-ids of a database with parsed
// Seq-ids depends on SEQ_ID_PREFER_ACCESSION_OVER_GI and SEQ_ID_AVOID_GI.
#[test]
fn a_seq_id_variable_refuses_a_database_with_parsed_seqids_only() {
    assert!(open("prot_ids_v4").parsed_seqids());
    assert!(!open("prot_v4").parsed_seqids());
    let settings = ApplicationSettings {
        seq_id_rank_variable: Some("SEQ_ID_AVOID_GI"),
        ..ApplicationSettings::default()
    };
    let Err(error) = BlastxDatabase::open(fixture("prot_ids_v4"), &settings) else {
        panic!("a database with parsed Seq-ids opened with SEQ_ID_AVOID_GI")
    };
    assert_eq!(
        format!("{error:#}"),
        "the environment variable SEQ_ID_AVOID_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's BLASTX"
    );
    BlastxDatabase::open(fixture("prot_v4"), &settings)
        .expect("a database without parsed Seq-ids opens with SEQ_ID_AVOID_GI");
}

#[test]
fn a_record_of_the_stored_letters_encodes_to_the_searched_bytes() {
    // A record of the stored letters encodes every residue byte 1..=27 back to itself but
    // O (26): the search reads a `-subject` O as X (`subjects::search_residue`), and a
    // database's stored O as it is (`read_encoded`).
    let bytes: Vec<u8> = (1..28).collect();
    let searched: Vec<u8> = bytes
        .iter()
        .map(|&residue| if residue == 26 { 21 } else { residue })
        .collect();
    let record = FastaRecord {
        internal_id: "gnl|BL_ORD_ID|0".into(),
        title: Vec::new(),
        sequence: bytes
            .iter()
            .map(|&residue| NCBISTDAA_TO_AMINOACID[usize::from(residue)])
            .collect(),
        lowercase_masks: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let records = [record];
    assert_eq!(&*Subjects::Fasta(&records).encoded(0), &searched[..]);
}

/// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
/// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> BlastxDatabase {
    use crate::seqdb::alias::SearchPath;
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
    let search = SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None);
    let path = tests.join("db_search").join(name);
    BlastxDatabase {
        db: SeqDb::open_in(path.to_str().unwrap().as_bytes(), SeqType::Protein, &search).unwrap(),
        name: name.as_bytes().to_vec(),
        taxdb: shared::TaxdbRun::default(),
        parsed_seqids: false,
        program: PROGRAM,
    }
}

/// The results of a query batch of one query whose hit lists name the sequences `lists` by
/// ordinal id, each with one HSP of the score given (0: no Seq-align in a gapped search).
fn one_query(lists: &[(usize, i32)]) -> Vec<BatchResults> {
    use super::super::preliminary::PreliminaryHsp;
    use super::super::results::{Hsp, HspList};
    let query = FastaRecord {
        internal_id: "Query_1".into(),
        title: Vec::new(),
        sequence: b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec(),
        lowercase_masks: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let hsp = |score: i32| Hsp {
        hsp: PreliminaryHsp {
            context: 0,
            frame: 1,
            score,
            q_start: 0,
            q_end: 3,
            q_gapped_start: 0,
            s_start: 0,
            s_end: 1,
            s_gapped_start: 0,
        },
        num: 1,
        evalue: 1e-9,
        bit_score: 50.0,
        identity: 1,
        positive: 1,
        edit_script: Vec::new(),
        composition_method: 0,
        query_left_out: 0,
    };
    vec![BatchResults {
        query_ordinal: 0,
        in_report_order: false,
        queries: vec![lists
            .iter()
            .map(|&(oid, score)| HspList {
                query_index: 0,
                oid,
                hsps: vec![hsp(score)],
                best_evalue: 1e-9,
            })
            .collect()],
        prepared: super::super::query_setup::prepare_queries(
            &[query],
            &super::super::stage_c_tests::options(),
        )
        .unwrap(),
        parameters: Vec::new(),
        search_skipped: false,
        request_space: None,
        unlisted: Vec::new(),
        unlisted_subject_row_throws: false,
        shown: Vec::new(),
        gapped_k: None,
        fetched_utf8: None,
    }]
}

// Line RG D1a-4 (fix round 7, review F N2): the one owner of the formatter's cut of a query's
// lists, which every report that NCBI cuts reads (`BatchResults::shown_lists`) and the loads
// of `ReportedSubjects::for_batch`: `shown` (RPSTBLASTN's count over several volumes, which
// joins following lists of one Seq-align id) wins over the cut, the cut keeps at most the
// lists there are, and `usize::MAX` (every other search) keeps every list.
#[test]
fn the_cut_of_a_querys_lists_reads_shown_before_the_cut() {
    let mut batches = one_query(&[(0, 90), (1, 80), (2, 70), (3, 60)]);
    let batch = &mut batches[0];
    let oids = |lists: &[super::super::results::HspList]| {
        lists.iter().map(|list| list.oid).collect::<Vec<_>>()
    };
    assert_eq!(oids(batch.cut_lists(0, usize::MAX)), [0, 1, 2, 3]);
    assert_eq!(oids(batch.cut_lists(0, 2)), [0, 1]);
    assert_eq!(oids(batch.cut_lists(0, 9)), [0, 1, 2, 3]);
    assert!(batch.cut_lists(0, 0).is_empty());
    batch.shown = vec![3];
    assert_eq!(oids(batch.cut_lists(0, 2)), [0, 1, 2]);
    assert_eq!(oids(batch.cut_lists(0, usize::MAX)), [0, 1, 2]);
    batch.shown = vec![9];
    assert_eq!(oids(batch.cut_lists(0, 2)), [0, 1, 2, 3]);
}

// DBID item T: a batch's prefetch looks up the sequences that its reports show (an HSP with a
// Seq-align), once per volume and taxid over the batches of a run, and opens the taxonomy
// database at the run's first lookup; nothing without `taxdb.bti`. `x16_taxid_99990` (two
// sequences of taxid 99990) with `offset_last_past_end` (the last entry, 99990, is past the
// end of `taxdb.btd`), `x16_taxid` (taxid 9606) with `bad_magic`.
#[test]
fn the_prefetch_looks_up_the_reported_sequences_once_per_run() {
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-db-prefetch-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let text = |database: &BlastxDatabase, lists: &[(usize, i32)]| -> String {
        let longest = database.x_check_stores().unwrap();
        let stream = super::super::x_db_stream::XStreamDb::new(&database.db, longest);
        let mut batches = one_query(lists);
        let titles = ReportedSubjects::read_titles(
            database,
            &batches,
            true,
            (false, Default::default()),
            Default::default(),
        )
        .unwrap();
        let reported = ReportedSubjects::for_batch(
            database,
            Subjects::Stream(super::super::x_db_stream::XStream {
                db: &stream,
                block: None,
            }),
            &mut batches,
            titles,
            true,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            crate::report::tabular::PrintedIds::default(),
            false,
            crate::report::db_subject::ReportRead {
                one_by_one: false,
                batch: Some(0),
                prefetch_fails: false,
            },
            usize::MAX,
            (crate::report::db_subject::FormatThrow::Unread, None),
            false,
        )
        .unwrap();
        let scope = crate::blastinput::query_batch::QueryScope::default;
        let messages = reported.batch_prefetch(
            database,
            &batches,
            true,
            &shared::PrefetchQueries::new(0, &scope),
        );
        String::from_utf8(messages.unwrap().text()).unwrap()
    };
    let database = open_with_taxdb("x16_taxid_99990", Some("offset_last_past_end"), &cwd);
    assert!(database.taxdb_found());
    // A sequence whose HSPs have no Seq-align is not fetched: the taxonomy database is not
    // opened.
    assert_eq!(text(&database, &[(0, 0)]), "");
    let offset = "Error: [blastx] Error: Offset error at end of taxdb file.\n";
    assert_eq!(text(&database, &[(0, 0), (1, 50)]), offset);
    // The same volume and taxid in a later batch.
    assert_eq!(text(&database, &[(0, 50), (1, 50)]), "");
    let database = open_with_taxdb("x16_taxid", Some("bad_magic"), &cwd);
    assert_eq!(
        text(&database, &[(1, 50), (0, 40)]),
        "Error: [blastx] Error: Tax database file has wrong magic number.\n"
    );
    let database = open_with_taxdb("x16_taxid", None, &cwd);
    assert!(!database.taxdb_found());
    assert_eq!(text(&database, &[(0, 50), (1, 40)]), "");
    std::fs::remove_dir_all(&cwd).unwrap();
}

// DBID P2: where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off, each query's report
// fetches the sequences it shows: the messages are kept by query ordinal (the writers write
// them after the query's warnings), once per volume and taxid over the run; a report whose
// fields read no subject fetches nothing. `x16_taxid_99990` with `offset_last_past_end`, as
// above; query 0 shows nothing (its HSP has no Seq-align), query 1 shows both sequences.
#[test]
fn with_the_prefetch_off_each_query_fetches_what_it_shows() {
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-db-prefetch-off-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    // `identities`: the engine's identities of every HSP.
    let fetched = |rows: SubjectFetch, identities: usize| -> [Vec<u8>; 3] {
        let database = open_with_taxdb("x16_taxid_99990", Some("offset_last_past_end"), &cwd);
        let longest = database.x_check_stores().unwrap();
        let stream = super::super::x_db_stream::XStreamDb::new(&database.db, longest);
        let mut batches = one_query(&[(0, 0)]);
        let mut second = one_query(&[(1, 50), (0, 40)]);
        second[0].query_ordinal = 1;
        batches.append(&mut second);
        let titles = ReportedSubjects::read_titles(
            &database,
            &batches,
            true,
            (false, Default::default()),
            Default::default(),
        )
        .unwrap();
        let mut reported = ReportedSubjects::for_batch(
            &database,
            Subjects::Stream(super::super::x_db_stream::XStream {
                db: &stream,
                block: None,
            }),
            &mut batches,
            titles,
            true,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            crate::report::tabular::PrintedIds::default(),
            false,
            crate::report::db_subject::ReportRead {
                one_by_one: false,
                batch: Some(0),
                prefetch_fails: false,
            },
            usize::MAX,
            (crate::report::db_subject::FormatThrow::Unread, None),
            false,
        )
        .unwrap();
        for hsp in batches
            .iter_mut()
            .flat_map(|batch| batch.queries.iter_mut().flatten())
            .flat_map(|list| list.hsps.iter_mut())
        {
            hsp.identity = identities;
        }
        reported
            .fetch_by_query(
                &database,
                &batches,
                true,
                (
                    rows,
                    crate::algorithm::common::database::PrefetchFormat::of(
                        6,
                        crate::report::archive::options::RequestKind::Blastx,
                    ),
                ),
                false,
                crate::algorithm::common::database::RowLookups::None,
                None,
                (usize::MAX, false, false),
            )
            .unwrap();
        [
            reported.query_fetches(0).to_vec(),
            reported.query_fetches(1).to_vec(),
            reported.take_query_fetches(),
        ]
    };
    let offset = b"Error: [blastx] Error: Offset error at end of taxdb file.\n".to_vec();
    let none = [Vec::new(), Vec::new(), Vec::new()];
    let second = [Vec::new(), offset.clone(), offset];
    assert_eq!(fetched(SubjectFetch::Every, 1), second);
    assert_eq!(fetched(SubjectFetch::Every, 0), second);
    assert_eq!(fetched(SubjectFetch::Never, 1), none);
    // L1: rows that read the subject only for its identities (`nident` alone) read none of
    // HSPs without identities.
    assert_eq!(fetched(SubjectFetch::WithIdentities, 1), second);
    assert_eq!(fetched(SubjectFetch::WithIdentities, 0), none);
    std::fs::remove_dir_all(&cwd).unwrap();
}

// RG-KB 9(b) review (PR #1637 small finding 3): in pairwise mode, a later subject's
// `StoredTitle::check_shown_source` refusal must stay in the deferred `refused` slot so it
// neither overtakes an earlier `check_reported_title` refusal nor overtakes a later damaged
// header record's `ErrorAfterBatchReports`.
#[test]
fn read_titles_does_not_let_shown_source_refusal_overtake_earlier_refusal_or_header_error() {
    use crate::seqdb::alias::SearchPath;
    let cwd = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-db-title-order-{}-{:?}",
        crate::utils::test_dir::process_id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&cwd).unwrap();
    let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let bare_v4 = cwd.join("bare_v4");
    let bare_args = crate::writedb::makeblastdb::MakedbArgs {
        input: std::ffi::OsString::from("-"),
        input_type: "fasta".to_string(),
        dbtype: "prot".to_string(),
        title: Some(std::ffi::OsString::from("bare_v4")),
        out: Some(bare_v4.as_os_str().to_os_string()),
        blastdb_version: 4,
        max_file_sz: std::ffi::OsString::from("1GB"),
        metadata_output_prefix: None,
        oid_masks: None,
        logfile: None,
        parse_seqids: false,
        hash_index: false,
        taxid: Some(9606),
        taxid_map: None,
        mask_data: None,
        mask_id: None,
        mask_desc: None,
        gi_mask: false,
        gi_mask_name: None,
        randomize: false,
    };
    let now = crate::writedb::time::LocalTime {
        year: 2026,
        month: 11,
        day: 25,
        hour: 12,
        minute: 0,
        second: 0,
        unix_seconds: 1_795_575_600,
    };
    let (mut log, mut messages) = (Vec::new(), Vec::new());
    let mut fasta: &[u8] = b">\nMDIVWSTFSARIQFAQYVPP\n";
    let status = crate::writedb::makeblastdb::run_with(
        &bare_args,
        &ApplicationSettings::default(),
        Some(&now),
        crate::writedb::Platform::host(),
        &mut fasta,
        &mut log,
        &mut messages,
    );
    assert_eq!(status, 0, "{}", String::from_utf8_lossy(&messages));
    let taxdb = tests.join("seqdb/taxdb/btd_empty");
    let search = SearchPath::new(&cwd, Some(taxdb.as_os_str()), None);
    let spec = format!(
        "{} {}",
        tests.join("db_search/tchar/tc_dmg_p").to_str().unwrap(),
        bare_v4.to_str().unwrap(),
    );
    let database = BlastxDatabase {
        db: SeqDb::open_in(spec.as_bytes(), SeqType::Protein, &search).unwrap(),
        name: spec.into_bytes(),
        taxdb: shared::TaxdbRun::default(),
        parsed_seqids: false,
        program: PROGRAM,
    };
    let read = |lists: &[(usize, i32)]| {
        let batches = one_query(lists);
        ReportedSubjects::read_titles(
            &database,
            &batches,
            true,
            (true, shared::TitleReads::Reports),
            Default::default(),
        )
    };
    // Multi-volume specs sort volumes by path, so locate bare_v4 (length 20) vs tc_dmg_p
    // (lengths 60, 60, 80: OID +0 has a 0xE9 title byte, OID +1 has a damaged header).
    assert_eq!(database.db.num_seqs(), 4);
    let bare_oid = (0..4usize)
        .find(|&oid| database.db.protein_length(oid as u32).unwrap() == 20)
        .expect("bare_v4 sequence");
    let dmg_base = if bare_oid == 0 { 1 } else { 0 };
    let bad_title_oid = dmg_base;
    let damaged_oid = dmg_base + 1;
    let err_source_only = read(&[(bare_oid, 90)]).err().expect("shown-source refusal");
    assert!(
        err_source_only
            .downcast_ref::<crate::blastinput::app::ErrorAfterBatchReports>()
            .is_none(),
        "expected shown-source refusal, got ErrorAfterBatchReports: {err_source_only:#}"
    );
    let err_title_first = read(&[(bad_title_oid, 90), (bare_oid, 80)])
        .err()
        .expect("bad-title refusal");
    assert!(
        format!("{err_title_first:#}")
            .contains("has a stored title with a byte outside ' ' to '~'"),
        "expected bad-title refusal first, got: {err_title_first:#}"
    );
    let err_damaged_after = read(&[(bare_oid, 90), (damaged_oid, 80)])
        .err()
        .expect("damaged-header ErrorAfterBatchReports");
    assert!(
        err_damaged_after
            .downcast_ref::<crate::blastinput::app::ErrorAfterBatchReports>()
            .is_some(),
        "expected damaged-header ErrorAfterBatchReports to win over shown-source refusal, got: {err_damaged_after:#}"
    );
    std::fs::remove_dir_all(&cwd).unwrap();
}
