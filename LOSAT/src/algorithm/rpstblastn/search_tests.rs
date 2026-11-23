use super::*;
use crate::algorithm::blastp::hsp::BlastpHspList;
use crate::common::GapEditOp;

fn options(extra: &[&str]) -> ResolvedOptions {
    let mut argv = vec!["losat", "rpstblastn", "-query", "q.fna", "-db", "cdd"];
    argv.extend_from_slice(&["-comp_based_stats", "0"]);
    argv.extend_from_slice(extra);
    let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
    let crate::cli::Commands::Rpstblastn(args) = cli.command else {
        unreachable!()
    };
    args.into_blastx_args().unwrap().resolve().unwrap()
}

fn record(id: &str, letters: &[u8]) -> FastaRecord {
    FastaRecord {
        internal_id: id.to_string(),
        title: Vec::new(),
        sequence: letters.to_vec(),
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    }
}

/// 60 letters of a coding sequence (20 codons) whose six frames all have letters.
const DNA: &[u8] = b"ATGGCTAAAGGTGAACGTCTGATCGCTGTTGCTAAAGAACTGGGTCTGACCGCTAAAGAA";

// blast_engine.c:805-813 (cited at `RpsFrame`): the frames of each query in NCBI's context
// order, their letters at BLASTX's offsets between two sentinels, and the frames of the strand
// that `-strand` leaves out without letters.
#[test]
fn the_frames_of_a_query_are_blastxs_contexts() {
    let options = options(&[]);
    let records = [record("a", DNA), record("b", &DNA[..31])];
    let prepared = prepare_unit_queries(&records, None, &options).unwrap();
    let parameters = score_block(&prepared);
    let masked = vec![None; prepared.contexts.len()];
    let queries = frame_queries(&prepared, &parameters, &masked, true);
    assert_eq!(queries.len(), 2);
    let frames: Vec<i8> = queries[0].frames.iter().map(|f| f.frame).collect();
    assert_eq!(frames, [1, 2, 3, -1, -2, -3]);
    let lengths: Vec<usize> = queries[0].frames.iter().map(RpsFrame::length).collect();
    assert_eq!(lengths, [20, 19, 19, 20, 19, 19]);
    let lengths: Vec<usize> = queries[1].frames.iter().map(RpsFrame::length).collect();
    assert_eq!(lengths, [10, 10, 9, 10, 10, 9]);
    for (query, frame_query) in queries.iter().enumerate() {
        for (i, frame) in frame_query.frames.iter().enumerate() {
            let context = &prepared.contexts[query * NUM_FRAMES + i];
            assert_eq!(frame.aa_seq.first(), Some(&0));
            assert_eq!(frame.aa_seq.last(), Some(&0));
            assert_eq!(
                &frame.aa_seq[1..frame.aa_seq.len() - 1],
                &prepared.sequence_start[context.offset + 1..context.offset + 1 + context.length]
            );
            assert!(frame.valid);
            assert_eq!(frame.aa_seq_nomask, None);
            assert_eq!(frame_to_context(i32::from(frame.frame)), i);
        }
    }
    // The first frame translates the coding sequence (BLASTX's genetic code 1).
    let translated: Vec<u8> = queries[0].frames[0].aa_seq[1..21]
        .iter()
        .map(|&letter| b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ"[usize::from(letter)])
        .collect();
    assert_eq!(translated, b"MAKGERLIAVAKELGLTAKE");
    let plus = options_and_frames(&["-strand", "plus"]);
    assert_eq!(plus, [20, 19, 19, 0, 0, 0]);
    let minus = options_and_frames(&["-strand", "minus"]);
    assert_eq!(minus, [0, 0, 0, 20, 19, 19]);
}

fn options_and_frames(extra: &[&str]) -> Vec<usize> {
    let options = options(extra);
    let records = [record("a", DNA)];
    let prepared = prepare_unit_queries(&records, None, &options).unwrap();
    let parameters = score_block(&prepared);
    let masked = vec![None; prepared.contexts.len()];
    let queries = frame_queries(&prepared, &parameters, &masked, true);
    for frame in &queries[0].frames {
        assert_eq!(frame.valid, frame.length() > 0);
    }
    queries[0].frames.iter().map(RpsFrame::length).collect()
}

// blast_engine.c:813 (cited at `RpsFrame`): the last frame with letters.
#[test]
fn the_preliminary_stat_length_is_the_last_frame_with_letters() {
    let options_both = options(&[]);
    let records = [record("a", &DNA[..59])];
    let prepared = prepare_unit_queries(&records, None, &options_both).unwrap();
    let parameters = score_block(&prepared);
    let masked = vec![None; prepared.contexts.len()];
    let queries = frame_queries(&prepared, &parameters, &masked, true);
    // Frames of 59 letters: +1 19, +2 19, +3 19, -1 19, -2 19, -3 19 ((59 - i) / 3).
    assert_eq!(queries[0].stat_length(), queries[0].frames[5].length());
    let options_plus = options(&["-strand", "plus"]);
    let prepared = prepare_unit_queries(&[record("b", DNA)], None, &options_plus).unwrap();
    let parameters = score_block(&prepared);
    let masked = vec![None; prepared.contexts.len()];
    let queries = frame_queries(&prepared, &parameters, &masked, true);
    assert_eq!(queries[0].stat_length(), 19);
    // A query of 4 letters: one letter in frames +1, -1 only (frames +2, +3: (4 - 1) / 3 =
    // 1, (4 - 2) / 3 = 0 letters).
    let prepared = prepare_unit_queries(&[record("c", b"ACGT")], None, &options_both).unwrap();
    let parameters = score_block(&prepared);
    let masked = vec![None; prepared.contexts.len()];
    let queries = frame_queries(&prepared, &parameters, &masked, true);
    let lengths: Vec<usize> = queries[0].frames.iter().map(RpsFrame::length).collect();
    assert_eq!(lengths, [1, 1, 0, 1, 1, 0]);
    assert_eq!(queries[0].stat_length(), 1);
}

// blast_setup_cxx.cpp:72-77 and blast_query_info.c:379-381 (cited at
// `RpsQuery::translation_length`): each context starts one sentinel after the previous one,
// an empty context takes no byte, and the buffer's length is the last context's offset plus
// its length; BLASTX's layout (`prepare_unit_queries`) is NCBI's.
#[test]
fn the_translation_buffer_has_ncbis_layout() {
    let layout = |extra: &[&str], letters: &[u8]| -> (Vec<(usize, usize)>, usize) {
        let options = options(extra);
        let prepared = prepare_unit_queries(&[record("a", letters)], None, &options).unwrap();
        let parameters = score_block(&prepared);
        let masked = vec![None; prepared.contexts.len()];
        let queries = frame_queries(&prepared, &parameters, &masked, true);
        let first = prepared.contexts[0].offset;
        let contexts: Vec<(usize, usize)> = prepared.contexts[..NUM_FRAMES]
            .iter()
            .map(|c| (c.offset - first, c.length))
            .collect();
        let (offset, length) = contexts[NUM_FRAMES - 1];
        assert_eq!(
            queries[0].translation_length(),
            offset + length,
            "{extra:?}"
        );
        (contexts, queries[0].translation_length())
    };
    // Both strands: six frames of 20, 19, 19 letters.
    let (both, length) = layout(&[], DNA);
    assert_eq!(
        both,
        [(0, 20), (21, 19), (41, 19), (61, 20), (82, 19), (102, 19)]
    );
    assert_eq!(length, 121);
    // -strand plus: the minus frames are empty and start where the last plus frame ends.
    let (plus, length) = layout(&["-strand", "plus"], DNA);
    assert_eq!(
        plus,
        [(0, 20), (21, 19), (41, 19), (61, 0), (61, 0), (61, 0)]
    );
    assert_eq!(length, 61);
    // -strand minus: the plus frames take no byte.
    let (minus, length) = layout(&["-strand", "minus"], DNA);
    assert_eq!(minus, [(0, 0), (0, 0), (0, 0), (0, 20), (21, 19), (41, 19)]);
    assert_eq!(length, 60);
    // A query of 4 letters: frames +3 and -3 are empty.
    let (short, length) = layout(&[], b"ACGT");
    assert_eq!(short, [(0, 1), (2, 1), (4, 0), (4, 1), (6, 1), (8, 0)]);
    assert_eq!(length, 8);
}

fn hsp(frame: i32, score: i32, q: (usize, usize), s: (usize, usize)) -> BlastpHsp {
    BlastpHsp {
        num: 1,
        identity: 50.0,
        length: q.1 - q.0 + 2,
        mismatch: 1,
        gapopen: 1,
        q_start: q.0,
        q_end: q.1,
        s_start: s.0,
        s_end: s.1,
        gapped_q_start: (q.0 + 2) as i32,
        gapped_s_start: (s.0 + 2) as i32,
        e_value: 1e-5,
        bit_score: 30.5,
        num_ident: 7,
        query_context: 0,
        query_frame: frame,
        subject_frame: 0,
        query_length: 100,
        q_idx: 0,
        s_idx: 0,
        raw_score: score,
        // The query's columns (1-based closed), with one deletion after the fifth.
        gap_info: Some(match (q.1 - q.0 + 1) as u32 {
            columns if columns > 5 => vec![
                GapEditOp::Sub(5),
                GapEditOp::Del(1),
                GapEditOp::Sub(columns - 5),
            ],
            columns => vec![GapEditOp::Sub(columns)],
        }),
        num_positives: 9,
        comp_adjustment_method: 0,
    }
}

// The conversion of the RPS traceback's final HSPs (1-based closed coordinates after the role
// swap, `traceback::swap_roles`) into BLASTX's (0-based half-open, context = query * 6 +
// the frame's context).
#[test]
fn rps_hsps_become_blastx_hsps() {
    let converted = blastx_hsp(&hsp(-2, 41, (11, 30), (3, 23)), 2);
    assert_eq!(
        converted.hsp,
        PreliminaryHsp {
            context: 2 * 6 + 4,
            frame: -2,
            score: 41,
            q_start: 10,
            q_end: 30,
            q_gapped_start: 13,
            s_start: 2,
            s_end: 23,
            s_gapped_start: 5,
        }
    );
    assert_eq!(
        (
            converted.evalue,
            converted.bit_score,
            converted.identity,
            converted.positive,
            converted.num,
            converted.composition_method
        ),
        (1e-5, 30.5, 7, 9, 0, 0)
    );
    assert_eq!(
        converted.edit_script,
        [GapEditOp::Sub(5), GapEditOp::Del(1), GapEditOp::Sub(15)]
    );
    // Every frame's context.
    for (frame, context) in [(1, 0), (2, 1), (3, 2), (-1, 3), (-2, 4), (-3, 5)] {
        assert_eq!(
            blastx_hsp(&hsp(frame, 1, (1, 2), (1, 2)), 0).hsp.context,
            context
        );
    }
    let mut method = hsp(1, 1, (1, 2), (1, 2));
    method.comp_adjustment_method = 1;
    assert_eq!(blastx_hsp(&method, 0).composition_method, 1);
    // `found_end` FALSE (blast_gapalign.c:4650-4653, cited at `Hsp::seqalign_ends`): the
    // HSP's ends stop before its starts while its edit script runs on. The HSP keeps them; its
    // Seq-align ends where the script does.
    let mut short = hsp(-2, 41, (11, 30), (3, 23));
    short.q_end = 10;
    short.s_end = 2;
    let converted = blastx_hsp(&short, 0);
    assert_eq!((converted.hsp.q_end, converted.hsp.s_end), (10, 2));
    assert_eq!(converted.seqalign_ends(), (30, 23));
}

// NCBI sorts a list's HSPs before it makes their Seq-aligns (`Blast_HSPListSortByEvalue`,
// blast_seqalign.cpp:1577), on the HSPs' own ends. Of two HSPs with the same e-value, score
// and starts, the one whose profile end is higher comes first (`ScoreCompareHSPs`): here the
// HSP that ends on profile letter 22, before the `found_end` FALSE HSP whose Seq-align ends
// on letter 23 but whose own end is letter 2.
#[test]
fn the_sort_before_the_seqaligns_reads_the_hsps_own_ends() {
    let mut short = hsp(-2, 41, (11, 30), (3, 23));
    short.q_end = 10;
    short.s_end = 2;
    let mut whole = hsp(-2, 41, (11, 30), (3, 22));
    whole.gap_info = Some(vec![GapEditOp::Sub(20)]);
    for order in [[&short, &whole], [&whole, &short]] {
        let mut list = HspList {
            query_index: 0,
            oid: 0,
            hsps: order.iter().map(|h| blastx_hsp(h, 0)).collect(),
            best_evalue: 0.0,
        };
        list.sort_evalue();
        let ends: Vec<_> = list.hsps.iter().map(Hsp::seqalign_ends).collect();
        assert_eq!(ends, [(30, 22), (30, 23)]);
    }
}

// A query's hit list: one BLASTX list per profile with HSPs, in the hit list's order, named
// by the profile's ordinal id; the best e-value refreshed (`HspList::refresh`).
#[test]
fn a_hit_list_becomes_blastx_lists_in_its_order() {
    let list = |oid: u32, evalues: &[f64]| BlastpHspList {
        oid,
        query_index: 0,
        hsps: evalues
            .iter()
            .map(|&e| BlastpHsp {
                num: 1,
                e_value: e,
                ..hsp(1, 10, (1, 9), (1, 9))
            })
            .collect(),
        best_evalue: 0.0,
    };
    let hit_list = BlastpHitList {
        hsplist_count: 3,
        hsplist_max: 3,
        worst_evalue: 0.0,
        low_score: 0,
        heapified: false,
        hsplist_array: vec![
            list(5, &[1e-3, 1e-9]),
            list(1, &[]),
            list(2, &[0.5]),
            list(9, &[1e-30]),
        ],
        hsplist_current: 0,
        num_hits: 0,
    };
    let lists = blastx_lists(Some(hit_list), 3);
    let oids: Vec<usize> = lists.iter().map(|l| l.oid).collect();
    assert_eq!(oids, [5, 2]);
    assert!(lists.iter().all(|l| l.query_index == 3));
    assert_eq!(lists[0].best_evalue, 1e-9);
    assert_eq!(lists[0].hsps.len(), 2);
    assert_eq!(lists[0].hsps[0].hsp.context, 3 * 6);
    assert!(blastx_lists(None, 0).is_empty());
}

// local_blast.cpp:176-181 (cited at `search_window`): a batch without a valid frame is not
// searched; its results are BLASTX's skipped batch.
#[test]
fn a_batch_without_a_valid_frame_is_not_searched() {
    let options = options(&[]);
    let records = [record("n", b"NNNNNN"), record("a", b"AC")];
    let prepared = prepare_unit_queries(&records, None, &options).unwrap();
    let parameters = score_block(&prepared);
    assert!(!parameters.iter().any(|p| p.valid));
    let results = unsearched_batch(&prepared, parameters);
    assert!(results.search_skipped);
    assert_eq!(results.queries.len(), 2);
    assert!(results.queries.iter().all(Vec::is_empty));
    assert_eq!(results.parameters.len(), 12);
}

// The refusals before the reports (`for_database`).
#[test]
fn sum_statistics_and_an_ungapped_search_are_refused_before_the_reports() {
    let mut search = TranslatedRpsSearch {
        search: RpsSearch::default(),
        ungapped: false,
        several_volumes: false,
        costs: None,
        strategy_db_length: None,
        units: None,
        strategy: false,
        setup: None,
    };
    assert!(search.check_options().is_ok());
    search.search.sum_stats = true;
    assert!(search
        .check_options()
        .unwrap_err()
        .to_string()
        .contains("-sum_stats true"));
    search.search.sum_stats = false;
    search.ungapped = true;
    assert!(search
        .check_options()
        .unwrap_err()
        .to_string()
        .contains("crashes after an ungapped search"));
    assert!(!search.merges_volumes());
}

// D1 of PS12s: a window's events are replayed in the reader's order, each warning written and
// each batch reported, up to the first batch that failed, whose error ends the replay after the
// events before it (the warnings and batches after it are not written, as NCBI's loop never
// reads them); a window without a failure writes everything.
#[test]
fn a_window_replays_the_readers_events_up_to_the_first_failed_batch() {
    let events = vec![
        RpsEvent::Warn(b"w0\n".to_vec()),
        RpsEvent::Batch(vec![record("a", DNA)]),
        RpsEvent::Batch(vec![record("b", DNA), record("c", DNA)]),
        RpsEvent::Warn(b"w1\n".to_vec()),
        RpsEvent::Batch(vec![record("d", DNA)]),
        RpsEvent::Warn(b"w2\n".to_vec()),
        RpsEvent::Batch(vec![record("e", DNA)]),
    ];
    // The log of a replay: warnings, each batch's request (`take:` and its first query) and
    // report; `take_fails` names the first query of a batch whose request fails.
    let run = |searched: Vec<Option<Result<u32>>>,
               take_fails: Option<&str>|
     -> (Vec<String>, Option<String>) {
        let log = std::cell::RefCell::new(Vec::new());
        let result = replay(
            &events,
            searched,
            &mut |bytes| {
                log.borrow_mut()
                    .push(String::from_utf8_lossy(bytes).trim().to_string());
                Ok(())
            },
            &mut |queries| {
                let first = queries[0].internal_id.as_str();
                if take_fails == Some(first) {
                    anyhow::bail!("request {first} failed");
                }
                log.borrow_mut().push(format!("take:{first}"));
                Ok(())
            },
            &mut |queries, results| {
                let ids: Vec<&str> = queries.iter().map(|q| q.internal_id.as_str()).collect();
                log.borrow_mut()
                    .push(format!("{}={results}", ids.join("+")));
                Ok(())
            },
        );
        (
            log.into_inner(),
            result.err().map(|error| error.to_string()),
        )
    };
    let lines =
        |lines: &[&str]| -> Vec<String> { lines.iter().map(|line| line.to_string()).collect() };
    let all = run(
        vec![Some(Ok(1)), Some(Ok(2)), Some(Ok(3)), Some(Ok(4))],
        None,
    );
    let expected = lines(&[
        "w0", "take:a", "a=1", "take:b", "b+c=2", "w1", "take:d", "d=3", "w2", "take:e", "e=4",
    ]);
    assert_eq!(all, (expected, None));
    // A failing batch has its request; the batches after it have none.
    let failed = run(
        vec![
            Some(Ok(1)),
            Some(Ok(2)),
            Some(Err(anyhow::anyhow!("batch 3 failed"))),
            None,
        ],
        None,
    );
    assert_eq!(
        failed,
        (
            lines(&["w0", "take:a", "a=1", "take:b", "b+c=2", "w1", "take:d"]),
            Some("batch 3 failed".to_string())
        )
    );
    let first = run(
        vec![
            Some(Err(anyhow::anyhow!("batch 1 failed"))),
            None,
            None,
            None,
        ],
        None,
    );
    assert_eq!(
        first,
        (lines(&["w0", "take:a"]), Some("batch 1 failed".to_string()))
    );
    // A batch's failing request comes after the batches before it and before the batch's
    // own error.
    let request = run(
        vec![
            Some(Ok(1)),
            Some(Ok(2)),
            Some(Err(anyhow::anyhow!("batch 3 failed"))),
            None,
        ],
        Some("d"),
    );
    assert_eq!(
        request,
        (
            lines(&["w0", "take:a", "a=1", "take:b", "b+c=2", "w1"]),
            Some("request d failed".to_string())
        )
    );
}

/// The files of the fixture database `name` (tests/db_search/rps/NAME.*.gz) inflated into a
/// directory of this test but those whose names end with one of `left_out`, and the database
/// opened.
fn rps_database(name: &str, left_out: &[&str], test: &str) -> (std::path::PathBuf, BlastxDatabase) {
    use std::io::Read;
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/rps");
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-rpstbn-{test}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    for entry in std::fs::read_dir(&fixtures).unwrap() {
        let file = entry.unwrap().file_name().into_string().unwrap();
        let Some(inflated) = file.strip_suffix(".gz") else {
            continue;
        };
        let Some(extension) = inflated.strip_prefix(name) else {
            continue;
        };
        if !extension.starts_with('.') || left_out.contains(&extension) {
            continue;
        }
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(std::fs::read(fixtures.join(&file)).unwrap().as_slice())
            .read_to_end(&mut bytes)
            .unwrap();
        std::fs::write(dir.join(inflated), bytes).unwrap();
    }
    let database = BlastxDatabase::open(
        dir.join(name).to_str().unwrap(),
        &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
    )
    .unwrap();
    (dir, database)
}

/// `search_window` of `batches` with the RPS options `search` on `database`, serially.
fn window_errors(
    search: &TranslatedRpsSearch,
    database: &BlastxDatabase,
    batches: &[&[FastaRecord]],
) -> Vec<Option<std::result::Result<(), String>>> {
    let options = options(&[]);
    let volumes = RpsVolumeSet::of(database.x_seqdb());
    let setup = RpsVolumeSetup::new(PROGRAM);
    let units = RpsUnits::of(database.x_seqdb(), 0, 0, PROGRAM).unwrap();
    search
        .search_window(
            batches,
            &options,
            database,
            &volumes,
            &units,
            &setup,
            &crate::utils::threading::SearchPool::x_serial(),
            &RpsWorkers::new(1),
        )
        .into_iter()
        .map(|outcome| outcome.map(|result| result.map(|_| ()).map_err(|e| e.to_string())))
        .collect()
}

/// The record `w45` of nwin.fna (build_nq.py: a back-translated window of a profile's
/// fragment, which hits `rps_mini`).
fn hit_record() -> FastaRecord {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta/rps/nwin.fna"),
    )
    .unwrap();
    let block = text
        .split('>')
        .find(|block| block.starts_with("w45 "))
        .unwrap();
    let letters: String = block.lines().skip(1).collect();
    record("w45", letters.as_bytes())
}

// B1 of the review of PS12s: a volume's refusal met by a later batch of the window does not
// end the volume loop for the batches before it (without a valid frame), which set up every
// volume as NCBI's loop of one batch at a time does and meet NCBI's error of a later volume's
// missing file first (rpsblast_local.cpp:384-399: one `CLocalBlast` per volume, whose set-up
// maps the RPS files; probe e1 of the review: NCBI rc 3).
#[test]
fn a_volumes_refusal_leaves_the_earlier_batches_to_meet_later_volumes() {
    let (dir, database) = rps_database("rps_mini3v", &[".01.loo"], "b1");
    let mut search = TranslatedRpsSearch {
        search: RpsSearch::default(),
        ungapped: false,
        several_volumes: true,
        costs: None,
        strategy_db_length: None,
        units: None,
        strategy: false,
        setup: None,
    };
    search.search.x_drop_gapped_bits = Some(1e9);
    let n_only = [record("n", &[b'N'; 210])];
    let valid = [hit_record()];
    let outcome = window_errors(&search, &database, &[&n_only, &valid]);
    let Some(Err(error)) = &outcome[0] else {
        panic!("{outcome:?}");
    };
    assert!(
        error.contains("Cannot memory map RPS-BLAST database file")
            && error.contains("rps_mini3v.01.loo"),
        "{error}"
    );
    assert!(outcome[1].is_none(), "{outcome:?}");
    // With every file the X-drop's refusal is the valid batch's, after the first batch's
    // results.
    let (whole_dir, whole) = rps_database("rps_mini3v", &[], "b1w");
    let outcome = window_errors(&search, &whole, &[&n_only, &valid]);
    assert!(matches!(outcome[0], Some(Ok(()))), "{outcome:?}");
    assert!(
        matches!(&outcome[1], Some(Err(error)) if error.contains("an -xdrop_gap of 1000000000.0 bits")),
        "{outcome:?}"
    );
    let _ = std::fs::remove_dir_all(dir);
    let _ = std::fs::remove_dir_all(whole_dir);
}

// Line RG D1a-4 fix round 2: NCBI's cut of a merged list to the hit list size counts following
// profiles of one Seq-align id as one subject (`PruneSeqalign`; `report_ordinal_ids`'s
// `formatter_kept`), and its reports find a reported profile that the range leaves out by its
// Seq-align id among the included ones (`found_as_copy`). `-db "ex12 dup/rps_mini3v.01"`:
// rps_mini3v.00 (gnl|CDD|101), rps_mini3v.01 (gnl|CDD|102) and a copy of rps_mini3v.01; the
// alias's `LAST_OID 2` is the range of the whole -db (the copy out), `LAST_OID 1` leaves both
// gnl|CDD|102 out. With the cut at 1, a profile of the reported Seq-id past the cut is
// reported: as its included copy (NCBI's rows of both, exit 0: probe dup6), or LOSAT's refusal
// where no included profile holds the id (NCBI's `Unknown` or null pointer exception: probe
// dup4_y1); after another profile it leaves the list past the cut, its HSPs counted.
#[test]
fn a_profile_of_the_reported_seq_id_past_the_cut_is_reported() {
    let (dir, _) = rps_database("rps_mini3v", &[], "dupcut");
    std::fs::create_dir_all(dir.join("dup")).unwrap();
    for entry in std::fs::read_dir(&dir).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        if name.starts_with("rps_mini3v.01.") {
            std::fs::copy(dir.join(&name), dir.join("dup").join(&name)).unwrap();
        }
    }
    let open = |last_oid: u32, copy_dir: &str| {
        std::fs::write(
            dir.join(format!("ex{last_oid}.pal")),
            format!("TITLE t\nDBLIST rps_mini3v.00 rps_mini3v.01\nLAST_OID {last_oid}\n"),
        )
        .unwrap();
        let name = format!(
            "{} {}",
            dir.join(format!("ex{last_oid}")).display(),
            dir.join(copy_dir).join("rps_mini3v.01").display()
        );
        BlastxDatabase::open(
            &name,
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .unwrap()
    };
    let list = |oid: u32| BlastpHspList {
        oid,
        query_index: 0,
        hsps: vec![hsp(1, 10, (1, 9), (1, 9))],
        best_evalue: 0.0,
    };
    let hit_list = |oids: &[u32]| BlastpHitList {
        hsplist_count: oids.len(),
        hsplist_max: oids.len(),
        worst_evalue: 0.0,
        low_score: 0,
        heapified: false,
        hsplist_array: oids.iter().map(|&oid| list(oid)).collect(),
        hsplist_current: 0,
        num_hits: 0,
    };
    let oids = |lists: &[Option<BlastpHitList>]| -> Vec<u32> {
        lists[0]
            .as_ref()
            .unwrap()
            .hsplist_array
            .iter()
            .map(|l| l.oid)
            .collect()
    };
    // The OIDs: the -db names in order, rps_mini3v.00, .01 and the copy.
    let (other, twin, copy) = (0, 1, 2);
    let database = open(2, "dup");
    let db = database.x_seqdb();
    assert_eq!(
        (0..3)
            .map(|oid| db.included_ordinal(oid))
            .collect::<Vec<_>>(),
        [Some(0), Some(1), None]
    );
    let mut joined = vec![Some(hit_list(&[twin, copy, other]))];
    assert_eq!(
        report_ordinal_ids(db, &mut joined, 1, PROGRAM, true).unwrap(),
        None
    );
    assert_eq!(oids(&joined), [1, 1, 0]);
    let mut apart = vec![Some(hit_list(&[twin, other, copy]))];
    assert_eq!(
        report_ordinal_ids(db, &mut apart, 1, PROGRAM, true).unwrap(),
        Some(vec![1])
    );
    assert_eq!(oids(&apart), [1, 0]);
    let database = open(1, "dup");
    let db = database.x_seqdb();
    let mut unfound = vec![Some(hit_list(&[twin, copy, other]))];
    let error = report_ordinal_ids(db, &mut unfound, 1, PROGRAM, true).unwrap_err();
    assert!(
        error.to_string().contains(&format!(
            "a hit of the profile of OID {twin} outside the RPS-BLAST database's alias files' FIRST_OID, LAST_OID and OIDLIST (NCBI BLAST+ 2.17.0's report does not find it"
        )),
        "{error}"
    );
    drop(database);
    // Review C of D1a-4 (SF2 d, N4): where the included profile of the Seq-id is another
    // record (another title: `duph`) or holds other letters (`dupl`), NCBI shows that
    // profile's names or letters with the HSPs of the one out of the range: LOSAT's refusal
    // (`Found::Other`), not the copy's rows.
    for (copy_dir, file, at, byte) in [
        ("dupl", "rps_mini3v.01.psq", 1, 0x0a),
        ("duph", "rps_mini3v.01.phr", 0x0f, b'R'),
    ] {
        std::fs::create_dir_all(dir.join(copy_dir)).unwrap();
        for entry in std::fs::read_dir(dir.join("dup")).unwrap() {
            let name = entry.unwrap().file_name();
            std::fs::copy(dir.join("dup").join(&name), dir.join(copy_dir).join(&name)).unwrap();
        }
        let path = dir.join(copy_dir).join(file);
        let mut bytes = std::fs::read(&path).unwrap();
        assert_ne!(bytes[at], byte, "{copy_dir}");
        bytes[at] = byte;
        std::fs::write(&path, bytes).unwrap();
        let database = open(2, copy_dir);
        let db = database.x_seqdb();
        let mut found_other = vec![Some(hit_list(&[twin, copy, other]))];
        let error = report_ordinal_ids(db, &mut found_other, 1, PROGRAM, true).unwrap_err();
        assert!(
            error.to_string().contains(&format!(
                "a hit of the profile of OID {copy} outside the RPS-BLAST database's alias files' FIRST_OID, LAST_OID and OIDLIST, whose Seq-id NCBI BLAST+ 2.17.0's report finds in another record or other letters (the profile of OID {twin})"
            )),
            "{copy_dir}: {error}"
        );
    }
    let _ = std::fs::remove_dir_all(dir);
}

// S1 of the review of PS12s: in a batch, a query's refusal of its cutoffs comes before an
// earlier query's search error (base's and NCBI's order, blast_engine.c:1372-1374). `w45`'s
// traceback reads the consensus letters, whose file is missing (NCBI's error); `p_tgg`'s frames
// have a Lambda whose ungapped X-drop of 9e6 bits leaves 32 bits (LOSAT's refusal).
#[test]
fn a_batchs_cutoff_refusal_comes_before_an_earlier_querys_search_error() {
    let (dir, database) = rps_database("rps_mini", &[".psq"], "s1");
    let mut search = TranslatedRpsSearch {
        search: RpsSearch::default(),
        ungapped: false,
        several_volumes: false,
        costs: None,
        strategy_db_length: None,
        units: None,
        strategy: false,
        setup: None,
    };
    search.search.x_drop_ungapped_bits = Some(9e6);
    let repeat = b"TGGTGTTGC".repeat(20);
    let both = [hit_record(), record("p_tgg", &repeat)];
    let outcome = window_errors(&search, &database, &[&both]);
    assert!(
        matches!(&outcome[0], Some(Err(error)) if error.contains("an -xdrop_ungap of 9000000.0 bits")),
        "{outcome:?}"
    );
    // In batches of their own the hit record's batch fails first, with NCBI's error.
    let hit = [hit_record()];
    let tgg = [record("p_tgg", &repeat)];
    let outcome = window_errors(&search, &database, &[&hit, &tgg]);
    assert!(
        matches!(&outcome[0], Some(Err(error)) if error.contains("rps_mini.psq")),
        "{outcome:?}"
    );
    assert!(outcome[1].is_none(), "{outcome:?}");
    let _ = std::fs::remove_dir_all(dir);
}

// blast_results.cpp:82-104 (cited at `footer_gapped_k`): only the first query's footer reads
// `kbp_gap[0]`, which the batch's last redo wrote (the highest profile: the stream reads the
// lists by ordinal id), and only where frame +1 is valid.
#[test]
fn only_the_first_footer_shows_the_redos_k() {
    let marks = [Some(2), None, Some(0)];
    let karlin_k = [0.1, 0.2, 0.3];
    let k = footer_gapped_k(&marks, &karlin_k, 0.041, true);
    assert_eq!(k, [1.2 * 0.3, 0.041, 0.041]);
    assert_eq!(footer_gapped_k(&marks, &karlin_k, 0.041, false), [0.041; 3]);
    assert_eq!(
        footer_gapped_k(&[None, None], &karlin_k, 0.041, true),
        [0.041; 2]
    );
}
