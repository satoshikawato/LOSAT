use super::*;

// NCBI BLAST+ 2.17.0, `tblastn -ungapped -comp_based_stats 0 -matrix FOO` with the queries of
// tests/fasta/outfmt0/sc1_q_mix.faa (an empty record, one with O at 77 and 97, one of 480
// letters, one of 32 X) against tests/db_search/tn_amb_s.fna, line O4c: each query's error,
// then the warnings of its set-up; no Karlin-Altschul message, which comes after the score
// block.
#[test]
fn an_unknown_matrix_fails_each_query_with_the_warnings_of_its_setup() {
    let mut with_o = vec![b'A'; 100];
    with_o[77] = b'O';
    with_o[97] = b'o';
    let queries = [Vec::new(), with_o, vec![b'M'; 480], vec![b'X'; 32]];
    assert_eq!(
        SetupEngineError::unknown_matrix(&queries).0,
        "Error: Unknown error code -1 Warning: Sequence contains no data \
         Error: Unknown error code -1 Warning: One or more O characters replaced by X for \
         alignment score calculations at positions 77, 97 \
         Error: Unknown error code -1 Error: Unknown error code -1 "
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-922 and :716-735
// (cited at `SubjectSetNumbers`)
// A total of 2^32 letters or more is NCBI's `Int8` on every target, wasm32's 32-bit `usize`
// included: a database's (an alias file's `LENGTH`, rows db.cnt.l4295.*, db.cnt.l4294.* and
// db.cnt.avg2p32.*), a FASTA set's sum of long records, and `-dbsize`. The effective lengths
// and the Gumbel block read it whole.
#[test]
fn subject_set_numbers_hold_totals_of_2_to_the_32_letters_or_more() {
    let lengths = [3_000_000_000usize, 3_000_000_000];
    let subjects = SubjectLetters::Lengths(&lengths);
    for total in [
        4_294_967_296_000_000_000i64,
        4_294_967_295_999_999_999,
        68_719_476_736,
    ] {
        assert!(total > i64::from(u32::MAX));
        let numbers = SubjectSetNumbers {
            database: Some(DatabaseNumbers {
                total_length: total,
                num_seqs: 16,
                gumbel_db_length: total,
                engine_db_length: total,
                min_length: 10,
                request_count: 16,
                average: (total / 16) as i32,
            }),
            ..SubjectSetNumbers::default()
        };
        assert_eq!(numbers.total_nt(&subjects), total);
        assert_eq!(numbers.dbsize_nt(), None);
        assert_eq!(numbers.eff_length_nt(total), total);
        assert_eq!(numbers.gumbel_db_length(total), total / 3);
    }
    // A FASTA set: the sum of its records.
    let fasta = SubjectSetNumbers::default();
    assert_eq!(fasta.total_nt(&subjects), 6_000_000_000);
    assert_eq!(fasta.eff_length_nt(6_000_000_000), 6_000_000_000);
    assert_eq!(fasta.gumbel_db_length(6_000_000_000), 2_000_000_000);
    // A `-dbsize` of 2^32 or more replaces the total.
    let dbsize = SubjectSetNumbers {
        dbsize: Some(4_294_967_296_000_000_000),
        ..fasta
    };
    assert_eq!(dbsize.dbsize_nt(), Some(4_294_967_296_000_000_000));
    assert_eq!(
        dbsize.eff_length_nt(6_000_000_000),
        4_294_967_296_000_000_000
    );
    assert_eq!(
        dbsize.gumbel_db_length(6_000_000_000),
        1_431_655_765_333_333_333
    );
}

// `BlastSeqSrcGetAvgSeqLen` (cited at `average_subject_length`), the length of the ungapped
// search's word cutoff: a FASTA set's total over its number of records in NCBI's `Uint4`
// (truncated; 0 for no record), a database's average as the database gives it.
#[test]
fn average_subject_length_is_the_sets_or_the_databases() {
    let fasta = SubjectSetNumbers::default();
    assert_eq!(
        fasta.average_subject_length(&SubjectLetters::Lengths(&[300, 400, 0])),
        233
    );
    assert_eq!(
        fasta.average_subject_length(&SubjectLetters::Lengths(&[])),
        0
    );
    // An average of 4e9 letters: NCBI's `Uint4`, given back as the source's `Int4`.
    assert_eq!(
        fasta.average_subject_length(&SubjectLetters::Lengths(&[4_000_000_000, 4_000_000_000])),
        -294_967_296
    );
    let database = SubjectSetNumbers {
        database: Some(DatabaseNumbers {
            total_length: 700,
            num_seqs: 2,
            gumbel_db_length: 700,
            engine_db_length: 700,
            min_length: 300,
            request_count: 2,
            average: 351,
        }),
        ..SubjectSetNumbers::default()
    };
    assert_eq!(
        database.average_subject_length(&SubjectLetters::Lengths(&[300, 400])),
        351
    );
}

// The search reads the subjects as slices of the records' residues (`tblastn_app.rs`).
fn subject_slices(subjects: &[Vec<u8>]) -> Vec<&[u8]> {
    subjects.iter().map(Vec::as_slice).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2777-2827;
// c++/src/algo/blast/core/blast_parameters.c:208-223,335-347:
// Points the active Karlin block to `kbp_std`; initial ungapped x-dropoff uses that context's Lambda.
// c++/src/algo/blast/core/aa_ungapped.c:570-590:
// For BDT63307, NCBI's query-specific x-drop is 15; using the matrix
// average's 16 changes which gapped seed survives into Kappa traceback.
#[test]
fn query_specific_word_xdrop_preserves_ncbi_bdt63307_seed() {
    let base = std::path::Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fasta"));
    // The residues of NCBI's reader, as TBLASTN reads the query and the subjects.
    let residues = |name: &str, config| {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource};
        let file = std::fs::File::open(base.join(name)).unwrap();
        read_all(&mut FastaInputSource::from_file(file, config), &mut |_| {
            Ok(())
        })
        .unwrap()
        .into_iter()
        .map(|record| record.sequence)
        .collect::<Vec<Vec<u8>>>()
    };
    use crate::blastinput::fasta_reader::ReaderConfig;
    let queries = residues("AvCLPV.faa", ReaderConfig::query("TBLASTN", true, false));
    let subjects = residues(
        "AvCLPV.fasta",
        ReaderConfig::subject("TBLASTN", false, false),
    );
    let seg = SegParams::default();
    let mut trace = StageDBoundaryTrace::default();
    let results = run_local_search(
        &queries[33..34],
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        LocalStageDScoring::default(),
        Some(&mut trace),
    )
    .unwrap();
    assert_eq!(trace.parameters.as_ref().unwrap().cutoffs[0].word_xdrop, 15);
    assert_eq!(
        trace.parameters.as_ref().unwrap().cutoffs[0].word_cutoff,
        36
    );
    let hsp = trace
        .preliminary
        .iter()
        .flat_map(|(_, hsps)| hsps)
        .map(|(_, hsp)| hsp)
        .find(|hsp| hsp.frame == -1 && hsp.s_start == 76752 && hsp.s_end == 76934)
        .unwrap();
    assert_eq!((hsp.q_gapped_start, hsp.s_gapped_start), (173, 76909));
    assert!(!results[0].lists().is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
// The pinned NCBI no-hit probe records context valid=0 for this query.
#[test]
fn all_invalid_query_skips_subject_search_and_returns_empty_hitlist() {
    let queries = vec![vec![b'W'; 60]];
    let subjects = vec![vec![b'A'; 600]];
    let seg = SegParams::default();
    let mut trace = StageDBoundaryTrace::default();
    let results = run_local_search(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        LocalStageDScoring::default(),
        Some(&mut trace),
    )
    .unwrap();
    assert_eq!(trace.query_validity, [false]);
    assert!(trace.preliminary.is_empty());
    assert_eq!(results.len(), 1);
    assert!(results[0].lists().is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:728-729,804-813;
// c++/src/algo/blast/core/blast_util.c:508-527,1070-1101:
// stat_length retains the raw nucleotide length until a nonempty frame
// replaces it, and the final translated context is frame -3.
#[test]
fn initial_stat_length_uses_last_nonempty_translated_frame() {
    for (nt_length, expected) in [
        (0, 0),
        (1, 1),
        (2, 2),
        (3, 1),
        (4, 1),
        (5, 1),
        (6, 1),
        (7, 1),
        (8, 2),
        (360, 119),
        (362, 120),
    ] {
        assert_eq!(initial_translated_stat_length(nt_length), expected);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905;
// core/blast_kappa.c:3525-3565; core/blast_traceback.c:717-719:
// compare the actual Rust C-to-D list and Kappa incoming order with the
// saved function-entry snapshots, before comparing retained results.
fn assert_local_stage_d_boundaries(
    trace: &StageDBoundaryTrace,
    calls: &str,
    kappa: Option<&str>,
    do_sum_stats: bool,
) {
    let initial_phase = if do_sum_stats {
        "link_before"
    } else {
        "evalue_before"
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905; c++/src/algo/blast/core/blast_traceback.c:1644-1707:
    // composition redo begins only after all preliminary lists were
    // written. In the ordinary path, the first posttraceback E-value
    // input has recomputed positive num_ident; initial HSPs have zero.
    let boundary = if kappa.is_some() {
        calls
            .lines()
            .find(|line| line.starts_with("D_CALL\t") && line.contains("\tredo\t"))
            .and_then(|line| calls.find(line))
            .expect("NCBI redo call")
    } else {
        calls
            .lines()
            .find(|line| {
                let f: Vec<_> = line.split('\t').collect();
                f.len() > 16
                    && f[0] == "D_HSP"
                    && f[2] == initial_phase
                    && f[11].parse::<i32>().unwrap() > 0
            })
            .and_then(|line| {
                let event = line.split('\t').nth(1).unwrap();
                calls.find(&format!("D_CALL\t{event}\t"))
            })
            .unwrap_or(calls.len())
    };
    let preliminary_calls = &calls[..boundary];
    let actual: Vec<_> = trace
        .preliminary
        .iter()
        .flat_map(|(_, hsps)| hsps.iter())
        .map(|(context, hsp)| {
            (
                *context as i32,
                hsp.frame as i32,
                hsp.q_start,
                hsp.q_end,
                hsp.s_start,
                hsp.s_end,
                hsp.score,
                hsp.s_gapped_start,
            )
        })
        .collect();
    let expected: Vec<_> = preliminary_calls
        .lines()
        .filter_map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            if f.len() > 16 && f[0] == "D_HSP" && f[2] == initial_phase {
                Some((
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap(),
                    f[16].parse().unwrap(),
                ))
            } else {
                None
            }
        })
        .collect();
    assert_eq!(actual, expected, "Stage C to D HSP input/order/count");
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905; c++/src/algo/blast/core/link_hsps.c:1765-1810:
    // each allocated subject list calls link with raw nucleotide length,
    // or direct E-values with translated length, before preliminary reap.
    let initial_call = if do_sum_stats { "link" } else { "evalue" };
    let call_rows: Vec<Vec<_>> = preliminary_calls
        .lines()
        .filter_map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (f.len() > 4 && f[0] == "D_CALL" && f[2] == initial_call).then_some(f)
        })
        .collect();
    assert_eq!(
        call_rows.len(),
        trace.preliminary.len(),
        "Stage C to D initial function call count"
    );
    for (row, (oid, _)) in call_rows.iter().zip(&trace.preliminary) {
        let nt_length = trace.subject_nt_lengths[*oid];
        assert_eq!(
            row[4].parse::<usize>().unwrap(),
            if do_sum_stats {
                nt_length
            } else {
                initial_translated_stat_length(nt_length)
            },
            "Stage C to D subject length at initial call"
        );
        if do_sum_stats {
            let link = trace.parameters.as_ref().unwrap().link.unwrap();
            assert_eq!(row[6].parse::<i32>().unwrap(), link.longest_intron);
            assert_eq!(
                row[7].parse::<f64>().unwrap().to_bits(),
                link.gap_decay_rate.to_bits()
            );
            assert_eq!(row[8].parse::<i32>().unwrap(), link.gap_size);
            assert_eq!(row[9].parse::<i32>().unwrap(), link.cutoff_small_gap);
        }
    }
    if let Some(params) = &trace.parameters {
        for line in calls.lines().filter(|line| line.starts_with("D_CONTEXT\t")) {
            let f: Vec<_> = line.split('\t').collect();
            let context: usize = f[2].parse().unwrap();
            assert_eq!(
                params.lengths[context].eff_searchsp,
                f[4].parse().unwrap(),
                "D effective search space"
            );
            assert_eq!(
                params.lengths[context].length_adjustment,
                f[5].parse().unwrap(),
                "D length adjustment"
            );
        }
    }
    if let Some(kappa) = kappa {
        let actual: Vec<_> = trace
            .redo
            .iter()
            .enumerate()
            .flat_map(|(event, (context, _, hsps))| {
                hsps.iter().enumerate().map(move |(index, hsp)| {
                    (
                        event,
                        index,
                        hsp.score,
                        *context,
                        hsp.frame as i32,
                        hsp.q_start,
                        hsp.q_end,
                        hsp.q_gapped_start,
                        hsp.s_start,
                        hsp.s_end,
                        hsp.s_gapped_start,
                    )
                })
            })
            .collect();
        let expected: Vec<_> = kappa
            .lines()
            .filter(|line| line.starts_with("K_TRACE_PRELIM\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].parse().unwrap(),
                    f[2].parse().unwrap(),
                    f[3].parse().unwrap(),
                    f[4].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[6].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap(),
                    f[11].parse().unwrap(),
                )
            })
            .collect();
        assert_eq!(actual, expected, "Kappa redo incoming HSPs/order");
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:717-719: stat_length is passed
        // to the posttraceback E-value call. D_HSP field 16 is instead
        // subject.gapped_start (ncbi_d_call_trace.c:113-118).
        let mut seen = std::collections::HashSet::new();
        let mut current_length = None;
        let mut expected = Vec::new();
        for line in calls.lines() {
            let f: Vec<_> = line.split('\t').collect();
            if f.len() > 4
                && f[0] == "D_CALL"
                && f[2] == if do_sum_stats { "link" } else { "evalue" }
            {
                current_length = Some(f[4].parse::<i32>().unwrap());
            } else if f.len() > 16 && f[0] == "D_HSP" && f[2] == "bits_after" && seen.insert(f[1]) {
                expected.push((f[4].parse::<usize>().unwrap(), current_length.unwrap()));
            }
        }
        let actual: Vec<_> = trace
            .posttraceback_lengths
            .iter()
            .map(|(context, _, length)| (*context, *length))
            .collect();
        assert_eq!(
            actual, expected,
            "posttraceback translated stat_length/call order"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3736,2494-2515;
// c++/src/algo/blast/core/blast_traceback.c:1763-1776:
// The 112-subject comparison-only local run supplies natural Stage C
// candidates and two retained HSPs after Kappa and the result post-pipes.
#[test]
fn natural_112_subject_local_path_matches_retained_ncbi_fields() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/"
    );
    let read_fasta = |name: &str| -> Vec<Vec<u8>> {
        let mut records: Vec<Vec<u8>> = Vec::new();
        for line in std::fs::read_to_string(format!("{root}{name}"))
            .unwrap()
            .lines()
        {
            if line.starts_with('>') {
                records.push(Vec::new());
            } else {
                records.last_mut().unwrap().extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("subjects.fna");
    let seg = SegParams::default();
    let mut boundary = StageDBoundaryTrace::default();
    let results = run_local_search(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10.0,
            max_target_seqs: 2,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        LocalStageDScoring::default(),
        Some(&mut boundary),
    )
    .unwrap();
    let calls =
        std::fs::read_to_string(format!("{root}natural_c_d_early_20260925/calls.tsv")).unwrap();
    let kappa = std::fs::read_to_string(format!("{root}result_order_20260925/ncbi.trace")).unwrap();
    assert_local_stage_d_boundaries(&boundary, &calls, Some(&kappa), true);
    assert_eq!(results.len(), 1);
    let hitlists = results[0].lists();
    assert_eq!(hitlists.iter().map(|h| h.oid).collect::<Vec<_>>(), [10, 11]);
    let report =
        std::fs::read_to_string(format!("{root}report_payload_20260925/report_fields.tsv"))
            .unwrap();
    for (list, row) in hitlists.iter().zip(report.lines()) {
        let f: Vec<_> = row.split('\t').collect();
        assert_eq!(list.payloads.len(), 1);
        let hsp = &list.hsps.hsps[0];
        let payload = &list.payloads[0];
        assert_eq!(f[1], format!("weak_{}", list.oid + 7));
        assert_eq!(hsp.hsp.score, f[2].parse().unwrap());
        assert_eq!(payload.num_ident, f[3].parse().unwrap());
        assert_eq!(payload.num_positives, f[4].parse().unwrap());
        assert_eq!(payload.align_length, f[5].parse().unwrap());
        assert_eq!(payload.mismatches, f[6].parse().unwrap());
        assert_eq!(payload.gap_letters, f[7].parse().unwrap());
        assert_eq!(payload.gap_opens, f[8].parse().unwrap());
        assert_eq!(hsp.hsp.q_start + 1, f[9].parse().unwrap());
        assert_eq!(hsp.hsp.q_end, f[10].parse().unwrap());
        assert_eq!(
            3 * hsp.hsp.s_start + i32::from(hsp.hsp.frame),
            f[11].parse().unwrap()
        );
        assert_eq!(
            3 * hsp.hsp.s_end + i32::from(hsp.hsp.frame) - 1,
            f[12].parse().unwrap()
        );
        assert_eq!(hsp.hsp.frame, f[13].parse().unwrap());
        assert!(!payload.edit_script.is_empty());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3577-3736;
// c++/src/algo/blast/core/blast_traceback.c:1763-1776:
// Three query contexts share one local nucleotide subject; the third
// contains only X and NCBI leaves it invalid with no retained HSPs.
#[test]
fn natural_multi_query_local_path_matches_retained_ncbi_scores_and_coordinates() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let read_fasta = |name: &str| -> Vec<Vec<u8>> {
        let mut records: Vec<Vec<u8>> = Vec::new();
        for line in std::fs::read_to_string(format!("{root}{name}"))
            .unwrap()
            .lines()
        {
            if line.starts_with('>') {
                records.push(Vec::new());
            } else {
                records.last_mut().unwrap().extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("subjects.fna");
    let seg = SegParams::default();
    let mut boundary = StageDBoundaryTrace::default();
    let results = run_local_search(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        LocalStageDScoring::default(),
        Some(&mut boundary),
    )
    .unwrap();
    let calls = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.trace"
    ))
    .unwrap();
    let kappa = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/multi_query_20260924_default.tsv"
    )).unwrap();
    assert_local_stage_d_boundaries(&boundary, &calls, Some(&kappa), true);
    let expected = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/run_20260924/multi_query_20260924_default.out"
    ))
    .unwrap();
    let mut rows = Vec::new();
    for (query_index, hitlist) in results.iter().enumerate() {
        for list in hitlist.lists() {
            for hsp in &list.hsps.hsps {
                let frame = i32::from(hsp.hsp.frame);
                let subject_len = subjects[list.oid as usize].len() as i32;
                let (sstart, send) = if frame > 0 {
                    (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                } else {
                    (
                        subject_len - 3 * hsp.hsp.s_start + frame + 1,
                        subject_len - 3 * hsp.hsp.s_end + frame + 2,
                    )
                };
                rows.push((
                    query_index,
                    hsp.hsp.score,
                    hsp.hsp.q_start + 1,
                    hsp.hsp.q_end,
                    sstart,
                    send,
                    frame,
                ));
            }
        }
    }
    let expected_rows: Vec<_> = expected
        .lines()
        .map(|row| {
            let f: Vec<_> = row.split('\t').collect();
            let query_index = match f[0] {
                "full" => 0,
                "internal70" => 1,
                _ => 2,
            };
            (
                query_index,
                f[2].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
            )
        })
        .collect();
    assert_eq!(rows, expected_rows);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2510;
    // c++/src/algo/blast/core/blast_hits.c:3243-3297:
    // Compare the unformatted per-HSP values after the full local path.
    let trace = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/multi_query_20260924_default.tsv")).unwrap();
    let heap_rows: std::collections::HashMap<_, _> = trace
        .lines()
        .filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t"))
        .map(|row| {
            let f: Vec<_> = row.split('\t').collect();
            (
                (
                    f[7].parse::<usize>().unwrap(),
                    f[8].parse::<i8>().unwrap(),
                    f[9].parse::<i32>().unwrap(),
                    f[10].parse::<i32>().unwrap(),
                    f[11].parse::<i32>().unwrap(),
                    f[12].parse::<i32>().unwrap(),
                ),
                f,
            )
        })
        .collect();
    for (context, hitlist) in results.iter().enumerate() {
        for list in hitlist.lists() {
            for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                let key = (
                    context,
                    hsp.hsp.frame,
                    hsp.hsp.q_start,
                    hsp.hsp.q_end,
                    hsp.hsp.s_start,
                    hsp.hsp.s_end,
                );
                let row = &heap_rows[&key];
                assert_eq!(hsp.hsp.score, row[3].parse().unwrap());
                assert_eq!(
                    payload.bit_score.to_bits(),
                    row[4].parse::<f64>().unwrap().to_bits()
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    row[5].parse::<f64>().unwrap().to_bits()
                );
                assert_eq!(payload.num_ident, row[6].parse().unwrap());
            }
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:800-812;
// core/blast_engine.c:1460-1466; core/blast_kappa.c:3687-3736:
// The comparison-only API looks up genetic code 32 via `FindGeneticCode` with the CLI-calibrated
// local seqSrc state and gives one retained HSP with full-precision values.
#[test]
fn code32_local_path_matches_cli_calibrated_api_oracle() {
    let input = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_a/fixtures/"
    );
    let oracle =
        concat!(env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/code32_local_api_20260925_cli_calibrated/code32.trace");
    let read = |name: &str| -> Vec<u8> {
        std::fs::read_to_string(format!("{input}{name}"))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(str::bytes)
            .collect()
    };
    let queries = vec![read("query.faa")];
    let subjects = vec![read("subject_code32.fna")];
    let seg = SegParams::default();
    let results = run_local_mode2_sum_stats(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 32,
            expect_value: 10.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
    )
    .unwrap();
    let trace = std::fs::read_to_string(oracle).unwrap();
    let row: Vec<_> = trace
        .lines()
        .find(|row| row.starts_with("K_TRACE_HEAP_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].lists().len(), 1);
    let list = &results[0].lists()[0];
    assert_eq!(list.hsps.hsps.len(), 1);
    let hsp = &list.hsps.hsps[0];
    let payload = &list.payloads[0];
    assert_eq!(hsp.hsp.score, row[3].parse().unwrap());
    assert_eq!(
        payload.bit_score.to_bits(),
        row[4].parse::<f64>().unwrap().to_bits()
    );
    assert_eq!(
        hsp.evalue.to_bits(),
        row[5].parse::<f64>().unwrap().to_bits()
    );
    assert_eq!(payload.num_ident, row[6].parse().unwrap());
    assert_eq!(hsp.hsp.frame, row[8].parse().unwrap());
    assert_eq!(
        (
            hsp.hsp.q_start,
            hsp.hsp.q_end,
            hsp.hsp.s_start,
            hsp.hsp.s_end
        ),
        (
            row[9].parse().unwrap(),
            row[10].parse().unwrap(),
            row[11].parse().unwrap(),
            row[12].parse().unwrap()
        )
    );
    assert!(!payload.edit_script.is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:200-251,
// 1481-1501,1644-1707,1763-1776:
// without composition statistics, ordinary traceback recomputes E-values,
// reaps, assigns bit scores, and inserts HSP lists before result sorting.
#[test]
fn ordinary_traceback_controls_match_local_ncbi_full_precision_fields() {
    for case in [
        "seg_hard_query_20260924",
        "multi_query_20260924",
        "run_20260923",
    ] {
        let input_root = format!(
            "{}/../docs/evidence/tlosan_stage_c/{case}",
            env!("CARGO_MANIFEST_DIR")
        );
        let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
            let mut records: Vec<(String, Vec<u8>)> = Vec::new();
            for line in std::fs::read_to_string(format!("{input_root}/{name}"))
                .unwrap()
                .lines()
            {
                if let Some(id) = line.strip_prefix('>') {
                    records.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                } else {
                    records.last_mut().unwrap().1.extend(line.bytes());
                }
            }
            records
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
        let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
        let seg = SegParams::default();
        let results = run_local_mode0_no_sum_stats(
            &query_sequences,
            &subject_slices(&subject_sequences),
            LocalStageDProfile {
                max_intron_length: 0,
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 500,
                query_records: &[],
                subject_set: Default::default(),
                engine_positives: false,
                culling_limit: 0,
                best_hit: None,
                subject_besthit: false,
                pssm: None,
                fetch_excluded: None,
                query_cov_hsp_perc: 0.0,
            },
        )
        .unwrap();
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_d/run_20260924/"
        );
        let trace = std::fs::read_to_string(format!("{root}{case}_control.trace")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("D_HSP\t") && line.contains("\tbits_after\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    (
                        f[4].parse::<usize>().unwrap(),
                        f[5].parse::<i8>().unwrap(),
                        f[6].parse::<i32>().unwrap(),
                        f[7].parse::<i32>().unwrap(),
                        f[8].parse::<i32>().unwrap(),
                        f[9].parse::<i32>().unwrap(),
                    ),
                    (
                        f[10].parse::<i32>().unwrap(),
                        f[11].parse::<i32>().unwrap(),
                        f[12].parse::<i32>().unwrap(),
                        f[13].parse::<f64>().unwrap().to_bits(),
                        f[14].parse::<f64>().unwrap().to_bits(),
                    ),
                )
            })
            .collect();
        let mut actual: Vec<_> = results
            .iter()
            .enumerate()
            .flat_map(|(query_index, hitlist)| {
                hitlist.lists().iter().flat_map(move |list| {
                    list.hsps
                        .hsps
                        .iter()
                        .zip(&list.payloads)
                        .map(move |(hsp, payload)| {
                            (
                                (
                                    query_index,
                                    hsp.hsp.frame,
                                    hsp.hsp.q_start,
                                    hsp.hsp.q_end,
                                    hsp.hsp.s_start,
                                    hsp.hsp.s_end,
                                ),
                                (
                                    hsp.hsp.score,
                                    payload.num_ident,
                                    hsp.num,
                                    hsp.evalue.to_bits(),
                                    payload.bit_score.to_bits(),
                                ),
                            )
                        })
                })
            })
            .collect();
        let mut expected = expected;
        actual.sort_unstable();
        expected.sort_unstable();
        assert_eq!(actual, expected, "{case} full-precision ordinary traceback");
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1776;
        // c++/src/algo/blast/core/blast_hits.c:3383-3400:
        // positive local TotLen sorts result lists by best E-value.
        // Here the saved diagnostic outfmt-6 columns witness Stage D
        // order, coordinates and alignment length; Stage E byte output
        // remains gated.
        let expected_output = std::fs::read_to_string(format!("{root}{case}_control.out")).unwrap();
        let expected_order: Vec<_> = expected_output
            .lines()
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[0].to_owned(),
                    f[1].to_owned(),
                    f[2].parse::<i32>().unwrap(),
                    f[5].parse::<i32>().unwrap(),
                    f[6].parse::<i32>().unwrap(),
                    f[7].parse::<i32>().unwrap(),
                    f[8].parse::<i32>().unwrap(),
                    f[9].parse::<i8>().unwrap(),
                    f[10].parse::<usize>().unwrap(),
                )
            })
            .collect();
        let mut actual_order = Vec::new();
        for (query_index, hitlist) in results.iter().enumerate() {
            for list in hitlist.lists() {
                for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                    let frame = i32::from(hsp.hsp.frame);
                    let subject_len = subjects[list.oid as usize].1.len() as i32;
                    let (sstart, send) = if frame > 0 {
                        (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                    } else {
                        (
                            subject_len - 3 * hsp.hsp.s_start + frame + 1,
                            subject_len - 3 * hsp.hsp.s_end + frame + 2,
                        )
                    };
                    actual_order.push((
                        queries[query_index].0.clone(),
                        subjects[list.oid as usize].0.clone(),
                        hsp.hsp.score,
                        hsp.hsp.q_start + 1,
                        hsp.hsp.q_end,
                        sstart,
                        send,
                        hsp.hsp.frame,
                        payload.align_length,
                    ));
                }
            }
        }
        assert_eq!(
            actual_order, expected_order,
            "{case} result order/report-ready fields"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:198-250,717-719;
// core/blast_kappa.c:3670-3733; core/blast_hits.c:3243-3297:
// every retained local Stage C FASTA fixture must preserve its final HSP
// values and order through the ordinary or composition result branch.
#[test]
fn remaining_five_local_fixture_families_match_ncbi_stage_d() {
    let croot = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/remaining_local_20260925/run_20260925/"
    );
    for case in [
        "ambiguity_20260923",
        "lowercase_20260923",
        "multi_hsp_20260924",
        "multi_hsp_no_fence_20260924",
        "six_frame_merged_20260924",
    ] {
        let read_fasta = |file: &str| -> Vec<(String, Vec<u8>)> {
            let mut rows = Vec::new();
            for line in std::fs::read_to_string(format!("{croot}{case}/{file}"))
                .unwrap()
                .lines()
            {
                if let Some(id) = line.strip_prefix('>') {
                    rows.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                } else {
                    rows.last_mut().unwrap().1.extend(line.bytes());
                }
            }
            rows
        };
        let queries = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
        let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
        let mut scoring = LocalStageDScoring::default();
        if case == "multi_hsp_no_fence_20260924" {
            scoring.gap_xdrop_bits = 5.0;
            scoring.final_xdrop_bits = 5.0;
        }
        for mode in [0, 2] {
            let name = format!("{case}.mode{mode}");
            let profile = LocalStageDProfile {
                max_intron_length: 0,
                seg: None,
                soft_masking: false,
                mask_lowercase: case == "lowercase_20260923",
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 500,
                query_records: &[],
                subject_set: Default::default(),
                engine_positives: false,
                culling_limit: 0,
                best_hit: None,
                subject_besthit: false,
                pssm: None,
                fetch_excluded: None,
                query_cov_hsp_perc: 0.0,
            };
            let mut boundary = StageDBoundaryTrace::default();
            let results = run_local_search(
                &query_sequences,
                &subject_slices(&subject_sequences),
                profile,
                if mode == 2 {
                    BlastCompoAdjustMode::CompositionMatrixAdjust
                } else {
                    BlastCompoAdjustMode::NoCompositionBasedStats
                },
                mode == 2,
                scoring,
                Some(&mut boundary),
            )
            .unwrap();
            assert_eq!(results.len(), 1, "{name} query count");
            let mut actual = Vec::new();
            let mut ordered = Vec::new();
            for list in results[0].lists() {
                for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                    actual.push((
                        hsp.hsp.frame,
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                        hsp.hsp.score,
                        payload.num_ident,
                        hsp.num,
                        hsp.evalue.to_bits(),
                        payload.bit_score.to_bits(),
                    ));
                    let frame = i32::from(hsp.hsp.frame);
                    let slen = subjects[list.oid as usize].1.len() as i32;
                    let (sstart, send) = if frame > 0 {
                        (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                    } else {
                        (
                            slen - 3 * hsp.hsp.s_start + frame + 1,
                            slen - 3 * hsp.hsp.s_end + frame + 2,
                        )
                    };
                    ordered.push((
                        (queries[0].0.clone(), subjects[list.oid as usize].0.clone()),
                        (
                            hsp.hsp.score,
                            payload.report_num_ident,
                            payload.align_length,
                            payload.report_mismatches,
                            payload.gap_letters,
                            payload.gap_opens,
                        ),
                        (
                            hsp.hsp.q_start + 1,
                            hsp.hsp.q_end,
                            sstart,
                            send,
                            hsp.hsp.frame,
                        ),
                    ));
                }
            }
            let trace = std::fs::read_to_string(format!(
                "{root}{name}.{}",
                if mode == 0 {
                    "calls.trace"
                } else {
                    "kappa.trace"
                }
            ))
            .unwrap();
            let calls = std::fs::read_to_string(format!("{root}{name}.calls.trace")).unwrap();
            assert_local_stage_d_boundaries(
                &boundary,
                &calls,
                if mode == 2 { Some(&trace) } else { None },
                mode == 2,
            );
            let mut expected: Vec<_> = trace
                .lines()
                .filter(|line| {
                    if mode == 0 {
                        line.starts_with("D_HSP\t") && line.contains("\tbits_after\t")
                    } else {
                        line.starts_with("K_TRACE_HEAP_HSP\t")
                    }
                })
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    if mode == 0 {
                        (
                            f[5].parse().unwrap(),
                            f[6].parse().unwrap(),
                            f[7].parse().unwrap(),
                            f[8].parse().unwrap(),
                            f[9].parse().unwrap(),
                            f[10].parse().unwrap(),
                            f[11].parse().unwrap(),
                            f[12].parse().unwrap(),
                            f[13].parse::<f64>().unwrap().to_bits(),
                            f[14].parse::<f64>().unwrap().to_bits(),
                        )
                    } else {
                        (
                            f[8].parse().unwrap(),
                            f[9].parse().unwrap(),
                            f[10].parse().unwrap(),
                            f[11].parse().unwrap(),
                            f[12].parse().unwrap(),
                            f[3].parse().unwrap(),
                            f[6].parse().unwrap(),
                            f[13].parse().unwrap(),
                            f[5].parse::<f64>().unwrap().to_bits(),
                            f[4].parse::<f64>().unwrap().to_bits(),
                        )
                    }
                })
                .collect();
            actual.sort_unstable();
            expected.sort_unstable();
            assert_eq!(actual, expected, "{name} full-precision internal HSPs");
            let expected_order: Vec<_> =
                std::fs::read_to_string(format!("{root}{name}.report.tsv"))
                    .unwrap()
                    .lines()
                    .map(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        (
                            (f[0].to_owned(), f[1].to_owned()),
                            (
                                f[2].parse().unwrap(),
                                f[5].parse().unwrap(),
                                f[7].parse().unwrap(),
                                f[8].parse().unwrap(),
                                f[9].parse().unwrap(),
                                f[10].parse().unwrap(),
                            ),
                            (
                                f[11].parse().unwrap(),
                                f[12].parse().unwrap(),
                                f[13].parse().unwrap(),
                                f[14].parse().unwrap(),
                                f[15].parse().unwrap(),
                            ),
                        )
                    })
                    .collect();
            assert_eq!(ordered, expected_order, "{name} order/report fields");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:283-307,870-905;
// core/blast_traceback.c:198-250,717-719,1481-1501;
// core/blast_kappa.c:3525-3736:
// generated masked and unmasked translated chunk boundaries use the same
// local Stage C search then ordinary or composition Stage D result path.
#[test]
fn four_generated_stage_c_chunk_families_match_ncbi_stage_d() {
    let croot = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/extended_chunks_20260925/run_20260925/"
    );
    let source = std::fs::read_to_string(format!("{croot}run_20260923/subjects.fna")).unwrap();
    let insert = source
        .split('>')
        .find(|record| record.starts_with("plus1\n"))
        .unwrap()
        .lines()
        .skip(1)
        .flat_map(str::bytes)
        .collect::<Vec<_>>();
    assert_eq!(insert.len(), 362);
    let query = std::fs::read_to_string(format!("{croot}run_20260923/query.faa"))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(str::bytes)
        .collect::<Vec<_>>();
    for case in [
        "long_chunk_20260924",
        "masked_chunk_boundary_20260924",
        "no_range_middle_20260924",
        "no_range_two_hits_20260924",
    ] {
        let mut subject = if case == "no_range_two_hits_20260924" {
            let first = 4_998_450;
            let second = 10_000_050;
            let mut seq = b"ATG".repeat(first);
            seq.extend(&insert);
            seq.push(b'A');
            seq.extend(b"ATG".repeat(second - first - 121));
            seq.extend(&insert);
            seq.extend(b"ATG".repeat(250));
            seq
        } else if case == "no_range_middle_20260924" {
            let mut seq = b"ATG".repeat(10_000_050);
            seq.extend(&insert);
            seq.extend(b"ATG".repeat(250));
            seq
        } else {
            let mut seq = b"ATG".repeat(4_999_950);
            seq.extend(&insert);
            seq.extend(b"ATG".repeat(250));
            seq
        };
        if case.starts_with("no_range") {
            subject[4_999_000 * 3..10_000_000 * 3].make_ascii_lowercase();
        } else if case == "masked_chunk_boundary_20260924" {
            subject[4_999_900 * 3..5_000_000 * 3].make_ascii_lowercase();
        }
        for mode in [0, 2] {
            let name = format!("{case}.mode{mode}");
            let profile = LocalStageDProfile {
                max_intron_length: 0,
                seg: None,
                soft_masking: false,
                mask_lowercase: case != "long_chunk_20260924",
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 500,
                query_records: &[],
                subject_set: Default::default(),
                engine_positives: false,
                culling_limit: 0,
                best_hit: None,
                subject_besthit: false,
                pssm: None,
                fetch_excluded: None,
                query_cov_hsp_perc: 0.0,
            };
            let mut boundary = StageDBoundaryTrace::default();
            let results = run_local_search(
                std::slice::from_ref(&query),
                &[subject.as_slice()],
                profile,
                if mode == 2 {
                    BlastCompoAdjustMode::CompositionMatrixAdjust
                } else {
                    BlastCompoAdjustMode::NoCompositionBasedStats
                },
                mode == 2,
                LocalStageDScoring::default(),
                Some(&mut boundary),
            )
            .unwrap();
            let calls = std::fs::read_to_string(format!("{root}{name}.calls.trace")).unwrap();
            let kappa = if mode == 2 {
                Some(std::fs::read_to_string(format!("{root}{name}.kappa.trace")).unwrap())
            } else {
                None
            };
            assert_local_stage_d_boundaries(&boundary, &calls, kappa.as_deref(), mode == 2);
            let actual: Vec<_> = results[0]
                .lists()
                .iter()
                .flat_map(|list| {
                    list.hsps
                        .hsps
                        .iter()
                        .zip(&list.payloads)
                        .map(|(hsp, payload)| {
                            let frame = i32::from(hsp.hsp.frame);
                            let slen = subject.len() as i32;
                            let (sstart, send) = if frame > 0 {
                                (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                            } else {
                                (
                                    slen - 3 * hsp.hsp.s_start + frame + 1,
                                    slen - 3 * hsp.hsp.s_end + frame + 2,
                                )
                            };
                            (
                                hsp.hsp.score,
                                payload.report_num_ident,
                                payload.align_length,
                                payload.report_mismatches,
                                payload.gap_letters,
                                payload.gap_opens,
                                hsp.hsp.q_start + 1,
                                hsp.hsp.q_end,
                                sstart,
                                send,
                                hsp.hsp.frame,
                                hsp.evalue.to_bits(),
                                payload.bit_score.to_bits(),
                            )
                        })
                })
                .collect();
            let report = std::fs::read_to_string(format!("{root}{name}.report.tsv")).unwrap();
            let expected: Vec<_> = report
                .lines()
                .map(|line| {
                    let f: Vec<_> = line.split('\t').collect();
                    (
                        f[2].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                        f[11].parse().unwrap(),
                        f[12].parse().unwrap(),
                        f[13].parse().unwrap(),
                        f[14].parse().unwrap(),
                        f[15].parse().unwrap(),
                    )
                })
                .collect();
            assert_eq!(actual.len(), expected.len(), "{name} retained count");
            for (
                (score, ident, length, mismatch, gaps, gapopen, q0, q1, s0, s1, frame, e, bits),
                expected,
            ) in actual.into_iter().zip(expected)
            {
                assert_eq!(
                    (score, ident, length, mismatch, gaps, gapopen, q0, q1, s0, s1, frame),
                    expected,
                    "{name} report values/order"
                );
                let trace = std::fs::read_to_string(format!(
                    "{root}{name}.{}",
                    if mode == 0 {
                        "calls.trace"
                    } else {
                        "kappa.trace"
                    }
                ))
                .unwrap();
                let row = trace
                    .lines()
                    .find(|line| {
                        let f: Vec<_> = line.split('\t').collect();
                        if mode == 0 {
                            f.len() > 14
                                && f[0] == "D_HSP"
                                && f[2] == "bits_after"
                                && f[10].parse::<i32>().unwrap() == score
                                && f[5].parse::<i8>().unwrap() == frame
                        } else {
                            f.len() > 13
                                && f[0] == "K_TRACE_HEAP_HSP"
                                && f[3].parse::<i32>().unwrap() == score
                                && f[8].parse::<i8>().unwrap() == frame
                        }
                    })
                    .unwrap();
                let f: Vec<_> = row.split('\t').collect();
                let (ecol, bcol) = if mode == 0 { (13, 14) } else { (5, 4) };
                assert_eq!(
                    e,
                    f[ecol].parse::<f64>().unwrap().to_bits(),
                    "{name} E-value"
                );
                assert_eq!(
                    bits,
                    f[bcol].parse::<f64>().unwrap().to_bits(),
                    "{name} bits"
                );
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:747-850,1446-1467;
// core/blast_traceback.c:198-250,717-719:
// the long six-frame subject enters the same local posttraceback call with
// translated statistical length and query-indexed effective search spaces.
#[test]
fn long_multi_query_local_stage_d_matches_ncbi() {
    let croot = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/long_subject_20260925/run_20260925/"
    );
    let query_text =
        std::fs::read_to_string(format!("{croot}long_multi_query_20260924/query.faa")).unwrap();
    let mut queries = Vec::new();
    let mut query_names = Vec::new();
    for line in query_text.lines() {
        if let Some(id) = line.strip_prefix('>') {
            query_names.push(id.split_whitespace().next().unwrap().to_owned());
            queries.push(Vec::new());
        } else {
            queries.last_mut().unwrap().extend(line.bytes());
        }
    }
    let source = std::fs::read_to_string(format!("{croot}run_20260923/subjects.fna")).unwrap();
    let plus1: Vec<u8> = source
        .split('>')
        .find_map(|record| {
            let mut lines = record.lines();
            (lines.next()?.split_whitespace().next()? == "plus1")
                .then(|| lines.flat_map(str::bytes).collect())
        })
        .unwrap();
    let mut subject = Vec::with_capacity(15_000_962);
    for _ in 0..4_999_950 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(&plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 15_000_962);
    let mut boundary = StageDBoundaryTrace::default();
    let results = run_local_search(
        &queries,
        &[subject.as_slice()],
        LocalStageDProfile {
            max_intron_length: 0,
            seg: None,
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10000.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
        BlastCompoAdjustMode::NoCompositionBasedStats,
        false,
        LocalStageDScoring::default(),
        Some(&mut boundary),
    )
    .unwrap();
    let mut actual = Vec::new();
    let mut ordered = Vec::new();
    for (qidx, hitlist) in results.iter().enumerate() {
        for list in hitlist.lists() {
            for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                actual.push((
                    qidx,
                    hsp.hsp.frame,
                    hsp.hsp.q_start,
                    hsp.hsp.q_end,
                    hsp.hsp.s_start,
                    hsp.hsp.s_end,
                    hsp.hsp.score,
                    payload.num_ident,
                    hsp.evalue.to_bits(),
                    payload.bit_score.to_bits(),
                ));
                ordered.push((
                    (query_names[qidx].clone(), "chunk_edge".to_owned()),
                    (
                        hsp.hsp.score,
                        payload.report_num_ident,
                        payload.align_length,
                        payload.report_mismatches,
                        payload.gap_letters,
                        payload.gap_opens,
                    ),
                    (
                        hsp.hsp.q_start + 1,
                        hsp.hsp.q_end,
                        3 * hsp.hsp.s_start + 1,
                        3 * hsp.hsp.s_end,
                        hsp.hsp.frame,
                    ),
                ));
            }
        }
    }
    let trace = std::fs::read_to_string(format!("{root}calls.trace")).unwrap();
    assert_local_stage_d_boundaries(&boundary, &trace, None, false);
    let mut expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("D_HSP\t") && line.contains("\tbits_after\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
                f[10].parse().unwrap(),
                f[11].parse().unwrap(),
                f[13].parse::<f64>().unwrap().to_bits(),
                f[14].parse::<f64>().unwrap().to_bits(),
            )
        })
        .collect();
    actual.sort_unstable();
    expected.sort_unstable();
    assert_eq!(actual, expected, "long subject full-precision final HSPs");
    let expected_order: Vec<_> = std::fs::read_to_string(format!("{root}report.tsv"))
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                (f[0].to_owned(), f[1].to_owned()),
                (
                    f[2].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap(),
                ),
                (
                    f[11].parse().unwrap(),
                    f[12].parse().unwrap(),
                    f[13].parse().unwrap(),
                    f[14].parse().unwrap(),
                    f[15].parse().unwrap(),
                ),
            )
        })
        .collect();
    assert_eq!(
        ordered, expected_order,
        "long subject final order/report fields"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2547-2605,4091-4171;
// core/blast_traceback.c:503-536,717-719: BLOSUM45 gap-14/2 uses
// its own score block and Karlin E-values with no Spouge gbp.
#[test]
fn alternate_blosum45_word2_local_stage_d_matches_ncbi() {
    let input = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
    );
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/alternate_matrix_20260925/run_20260925/"
    );
    let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        for line in std::fs::read_to_string(format!("{input}{name}"))
            .unwrap()
            .lines()
        {
            if let Some(id) = line.strip_prefix('>') {
                records.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
            } else {
                records.last_mut().unwrap().1.extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("subjects.fna");
    let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
    let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
    let mut boundary = StageDBoundaryTrace::default();
    let results = run_local_search(
        &query_sequences,
        &subject_slices(&subject_sequences),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: None,
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10000.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
        BlastCompoAdjustMode::NoCompositionBasedStats,
        false,
        LocalStageDScoring {
            matrix: ScoringMatrix::Blosum45,
            gap_open: 14,
            gap_extend: 2,
            word_size: 2,
            threshold: 16.0,
            window: 60,
            gap_xdrop_bits: 15.0,
            final_xdrop_bits: 25.0,
            word_xdrop_bits: 7.0,
            gapped: true,
            unknown_matrix: false,
            smith_waterman: false,
        },
        Some(&mut boundary),
    )
    .unwrap();
    let mut actual = Vec::new();
    let mut ordered = Vec::new();
    for (qidx, hitlist) in results.iter().enumerate() {
        for list in hitlist.lists() {
            for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                actual.push((
                    qidx,
                    hsp.hsp.frame,
                    hsp.hsp.q_start,
                    hsp.hsp.q_end,
                    hsp.hsp.s_start,
                    hsp.hsp.s_end,
                    hsp.hsp.score,
                    payload.num_ident,
                    hsp.evalue.to_bits(),
                    payload.bit_score.to_bits(),
                ));
                let frame = i32::from(hsp.hsp.frame);
                let slen = subjects[list.oid as usize].1.len() as i32;
                let (sstart, send) = if frame > 0 {
                    (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                } else {
                    (
                        slen - 3 * hsp.hsp.s_start + frame + 1,
                        slen - 3 * hsp.hsp.s_end + frame + 2,
                    )
                };
                ordered.push((
                    (
                        queries[qidx].0.clone(),
                        subjects[list.oid as usize].0.clone(),
                    ),
                    (
                        hsp.hsp.score,
                        payload.report_num_ident,
                        payload.align_length,
                        payload.report_mismatches,
                        payload.gap_letters,
                        payload.gap_opens,
                    ),
                    (
                        hsp.hsp.q_start + 1,
                        hsp.hsp.q_end,
                        sstart,
                        send,
                        hsp.hsp.frame,
                    ),
                ));
            }
        }
    }
    let trace = std::fs::read_to_string(format!("{root}calls.trace")).unwrap();
    assert_local_stage_d_boundaries(&boundary, &trace, None, false);
    let mut expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("D_HSP\t") && line.contains("\tbits_after\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[9].parse().unwrap(),
                f[10].parse().unwrap(),
                f[11].parse().unwrap(),
                f[13].parse::<f64>().unwrap().to_bits(),
                f[14].parse::<f64>().unwrap().to_bits(),
            )
        })
        .collect();
    actual.sort_unstable();
    expected.sort_unstable();
    assert_eq!(actual, expected, "BLOSUM45 full-precision retained HSPs");
    let expected_order: Vec<_> = std::fs::read_to_string(format!("{root}report.tsv"))
        .unwrap()
        .lines()
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                (f[0].to_owned(), f[1].to_owned()),
                (
                    f[2].parse().unwrap(),
                    f[5].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap(),
                ),
                (
                    f[11].parse().unwrap(),
                    f[12].parse().unwrap(),
                    f[13].parse().unwrap(),
                    f[14].parse().unwrap(),
                    f[15].parse().unwrap(),
                ),
            )
        })
        .collect();
    assert_eq!(
        ordered, expected_order,
        "BLOSUM45 final order and report fields"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:198-250,1486-1501;
// c++/src/algo/blast/core/blast_kappa.c:411-427,3670-3701:
// compositionBasedStats chooses redo while do_sum_stats independently
// chooses BLAST_LinkHsps or direct Blast_HSPListGetEvalues.
#[test]
fn independent_composition_and_sum_stats_options_match_ncbi() {
    let stage_c = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/option_cross_20260925/run_20260925/"
    );
    let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        for line in std::fs::read_to_string(format!("{stage_c}{name}"))
            .unwrap()
            .lines()
        {
            if let Some(id) = line.strip_prefix('>') {
                records.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
            } else {
                records.last_mut().unwrap().1.extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("subjects.fna");
    let query_sequences: Vec<_> = queries.iter().map(|(_, seq)| seq.clone()).collect();
    let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
    let seg = SegParams::default();
    for (name, mode) in [("mode0_sumtrue", 0), ("mode2_sumfalse", 2)] {
        let profile = LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10.0,
            max_target_seqs: 500,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        };
        let mut boundary = StageDBoundaryTrace::default();
        let result = run_local_search(
            &query_sequences,
            &subject_slices(&subject_sequences),
            profile,
            if mode == 2 {
                BlastCompoAdjustMode::CompositionMatrixAdjust
            } else {
                BlastCompoAdjustMode::NoCompositionBasedStats
            },
            mode == 0,
            LocalStageDScoring::default(),
            Some(&mut boundary),
        )
        .unwrap();
        let mut actual = Vec::new();
        let mut ordered = Vec::new();
        for (query_index, hitlist) in result.iter().enumerate() {
            for list in hitlist.lists() {
                for (hsp, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                    let frame = i32::from(hsp.hsp.frame);
                    let subject_len = subjects[list.oid as usize].1.len() as i32;
                    let (sstart, send) = if frame > 0 {
                        (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
                    } else {
                        (
                            subject_len - 3 * hsp.hsp.s_start + frame + 1,
                            subject_len - 3 * hsp.hsp.s_end + frame + 2,
                        )
                    };
                    actual.push((
                        query_index,
                        hsp.hsp.frame,
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                        hsp.hsp.score,
                        hsp.evalue.to_bits(),
                    ));
                    ordered.push((
                        (
                            queries[query_index].0.clone(),
                            subjects[list.oid as usize].0.clone(),
                        ),
                        (
                            hsp.hsp.score,
                            payload.report_num_ident,
                            payload.align_length,
                            payload.report_mismatches,
                            payload.gap_letters,
                            payload.gap_opens,
                        ),
                        (
                            hsp.hsp.q_start + 1,
                            hsp.hsp.q_end,
                            sstart,
                            send,
                            hsp.hsp.frame,
                        ),
                    ));
                }
            }
        }
        let trace = std::fs::read_to_string(format!(
            "{root}{name}.{}",
            if mode == 0 {
                "calls.trace"
            } else {
                "kappa.trace"
            }
        ))
        .unwrap();
        let calls = std::fs::read_to_string(format!("{root}{name}.calls.trace")).unwrap();
        assert_local_stage_d_boundaries(
            &boundary,
            &calls,
            if mode == 2 { Some(&trace) } else { None },
            mode == 0,
        );
        let mut expected: Vec<_> = trace
            .lines()
            .filter(|line| {
                if mode == 0 {
                    line.starts_with("D_HSP\t") && line.contains("\tbits_after\t")
                } else {
                    line.starts_with("K_TRACE_RESULT_HSP\t")
                }
            })
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                if mode == 0 {
                    (
                        f[4].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                        f[13].parse::<f64>().unwrap().to_bits(),
                    )
                } else {
                    (
                        usize::from(f[7] == "70"),
                        f[10].parse().unwrap(),
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[4].parse().unwrap(),
                        f[5].parse::<f64>().unwrap().to_bits(),
                    )
                }
            })
            .collect();
        actual.sort_unstable();
        expected.sort_unstable();
        assert_eq!(
            actual, expected,
            "{name} exact retained score/coords/E-value"
        );
        let expected_order: Vec<_> = std::fs::read_to_string(format!("{root}{name}.report.tsv"))
            .unwrap()
            .lines()
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    (f[0].to_owned(), f[1].to_owned()),
                    (
                        f[2].parse().unwrap(),
                        f[5].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap(),
                        f[10].parse().unwrap(),
                    ),
                    (
                        f[11].parse().unwrap(),
                        f[12].parse().unwrap(),
                        f[13].parse().unwrap(),
                        f[14].parse().unwrap(),
                        f[15].parse().unwrap(),
                    ),
                )
            })
            .collect();
        assert_eq!(
            ordered, expected_order,
            "{name} final order/report numeric fields"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-654;
// c++/src/algo/blast/core/blast_traceback.c:380-391,583-605,
// 1481-1501; c++/src/algo/blast/core/blast_kappa.c:3690-3736:
// hard/soft SEG and query lowercase flags choose the working-query bytes
// before Stage C, then each composition option follows its source branch.
#[test]
fn six_masking_fixtures_match_both_local_stage_d_option_paths() {
    let cases = [
        ("seg_cross_traceback_20260924", true, false, false),
        ("seg_soft_traceback_20260924", true, true, false),
        ("lcase_query_20260924", false, false, true),
        ("lcase_soft_query_20260924", false, true, true),
        ("seg_lcase_overlap_20260924", true, false, true),
        ("seg_lcase_overlap_soft_20260924", true, true, true),
    ];
    let read = |path: &str| -> Vec<u8> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(str::bytes)
            .collect()
    };
    let seg = SegParams::default();
    for (case, use_seg, soft, lowercase) in cases {
        let input = format!(
            "{}/../docs/evidence/tlosan_stage_c/{case}",
            env!("CARGO_MANIFEST_DIR")
        );
        let queries = vec![read(&format!("{input}/query.faa"))];
        let subjects = vec![read(&format!("{input}/subjects.fna"))];
        for mode in [0, 2] {
            let profile = LocalStageDProfile {
                max_intron_length: 0,
                seg: if use_seg { Some(&seg) } else { None },
                soft_masking: soft,
                mask_lowercase: lowercase,
                genetic_code: 1,
                expect_value: 10000.0,
                max_target_seqs: 500,
                query_records: &[],
                subject_set: Default::default(),
                engine_positives: false,
                culling_limit: 0,
                best_hit: None,
                subject_besthit: false,
                pssm: None,
                fetch_excluded: None,
                query_cov_hsp_perc: 0.0,
            };
            let mut boundary = StageDBoundaryTrace::default();
            let results = run_local_search(
                &queries,
                &subject_slices(&subjects),
                profile,
                if mode == 2 {
                    BlastCompoAdjustMode::CompositionMatrixAdjust
                } else {
                    BlastCompoAdjustMode::NoCompositionBasedStats
                },
                mode == 2,
                LocalStageDScoring::default(),
                Some(&mut boundary),
            )
            .unwrap();
            assert_eq!(results.len(), 1);
            assert_eq!(results[0].lists().len(), 1, "{case} mode {mode}");
            let list = &results[0].lists()[0];
            assert_eq!(list.hsps.hsps.len(), 1, "{case} mode {mode}");
            let hsp = &list.hsps.hsps[0];
            let payload = &list.payloads[0];
            let root = concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../docs/evidence/tlosan_stage_d/masking_options_20260925/run_20260925/"
            );
            let stem = format!("{case}.mode{mode}");
            let report = std::fs::read_to_string(format!("{root}{stem}.report.tsv")).unwrap();
            let f: Vec<_> = report.trim_end().split('\t').collect();
            assert_eq!(hsp.hsp.score, f[2].parse().unwrap(), "{stem} score");
            assert_eq!(
                (
                    payload.report_num_ident,
                    payload.align_length,
                    payload.report_mismatches,
                    payload.gap_letters,
                    payload.gap_opens
                ),
                (
                    f[5].parse().unwrap(),
                    f[7].parse().unwrap(),
                    f[8].parse().unwrap(),
                    f[9].parse().unwrap(),
                    f[10].parse().unwrap()
                ),
                "{stem} report fields"
            );
            let frame = i32::from(hsp.hsp.frame);
            let subject_len = subjects[0].len() as i32;
            let (sstart, send) = if frame > 0 {
                (3 * hsp.hsp.s_start + frame, 3 * hsp.hsp.s_end + frame - 1)
            } else {
                (
                    subject_len - 3 * hsp.hsp.s_start + frame + 1,
                    subject_len - 3 * hsp.hsp.s_end + frame + 2,
                )
            };
            assert_eq!(
                (
                    hsp.hsp.q_start + 1,
                    hsp.hsp.q_end,
                    sstart,
                    send,
                    hsp.hsp.frame
                ),
                (
                    f[11].parse().unwrap(),
                    f[12].parse().unwrap(),
                    f[13].parse().unwrap(),
                    f[14].parse().unwrap(),
                    f[15].parse().unwrap()
                ),
                "{stem} coordinates"
            );
            let trace = std::fs::read_to_string(format!(
                "{root}{stem}.{}",
                if mode == 0 {
                    "calls.trace"
                } else {
                    "kappa.trace"
                }
            ))
            .unwrap();
            let calls = std::fs::read_to_string(format!("{root}{stem}.calls.trace")).unwrap();
            assert_local_stage_d_boundaries(
                &boundary,
                &calls,
                if mode == 2 { Some(&trace) } else { None },
                mode == 2,
            );
            let row: Vec<_> = trace
                .lines()
                .find(|line| {
                    if mode == 0 {
                        line.starts_with("D_HSP\t") && line.contains("\tbits_after\t")
                    } else {
                        line.starts_with("K_TRACE_HEAP_HSP\t")
                    }
                })
                .unwrap()
                .split('\t')
                .collect();
            let (bits_col, e_col) = if mode == 0 { (14, 13) } else { (4, 5) };
            assert_eq!(
                payload.bit_score.to_bits(),
                row[bits_col].parse::<f64>().unwrap().to_bits(),
                "{stem} full-precision bits"
            );
            assert_eq!(
                hsp.evalue.to_bits(),
                row[e_col].parse::<f64>().unwrap().to_bits(),
                "{stem} full-precision E-value"
            );
            assert!(!payload.edit_script.is_empty(), "{stem} owned edit script");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/gc.prt:105-357;
// c++/src/algo/blast/core/blast_engine.c:762-774,1460-1466;
// c++/src/algo/blast/core/blast_kappa.c:3690-3736,2494-2515:
// the selected FindGeneticCode table enters subject translation before
// preliminary HSPs, redo, heap insertion and final hitlist updates.
// Each saved comparison uses the CLI-calibrated local API oracle, including
// code 32, and independently runs the identical input under code 1.
#[test]
fn all_27_codes_and_code1_controls_match_local_ncbi_result_traces() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/all_codes_20260925/"
    );
    let fixtures = std::fs::read_to_string(format!("{root}fixtures.tsv")).unwrap();
    let read = |name: &str| -> Vec<u8> {
        std::fs::read_to_string(format!("{root}fixtures/{name}"))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(str::bytes)
            .collect()
    };
    let seg = SegParams::default();
    let mut compared = 0;
    for row in fixtures.lines().skip(1) {
        let code: u8 = row.split('\t').next().unwrap().parse().unwrap();
        let queries = vec![read(&format!("code{code}.faa"))];
        let subjects = vec![read(&format!("code{code}.fna"))];
        for (selected, label) in [(1, "control"), (code, "selected")] {
            let mut boundary = StageDBoundaryTrace::default();
            let results = run_local_search(
                &queries,
                &subject_slices(&subjects),
                LocalStageDProfile {
                    max_intron_length: 0,
                    seg: Some(&seg),
                    soft_masking: false,
                    mask_lowercase: false,
                    genetic_code: selected,
                    expect_value: 10.0,
                    max_target_seqs: 500,
                    query_records: &[],
                    subject_set: Default::default(),
                    engine_positives: false,
                    culling_limit: 0,
                    best_hit: None,
                    subject_besthit: false,
                    pssm: None,
                    fetch_excluded: None,
                    query_cov_hsp_perc: 0.0,
                },
                BlastCompoAdjustMode::CompositionMatrixAdjust,
                true,
                LocalStageDScoring::default(),
                Some(&mut boundary),
            )
            .unwrap();
            let trace = std::fs::read_to_string(format!(
                "{root}run_20260925/code{code}.{label}.kappa.trace"
            ))
            .unwrap();
            let calls =
                std::fs::read_to_string(format!("{root}run_20260925/code{code}.{label}.d.trace"))
                    .unwrap();
            assert_local_stage_d_boundaries(&boundary, &calls, Some(&trace), true);
            let expected: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_RESULT_HSP\t"))
                .collect();
            let actual: Vec<_> = results[0]
                .lists()
                .iter()
                .flat_map(|list| list.hsps.hsps.iter().zip(&list.payloads))
                .collect();
            let report = std::fs::read_to_string(format!(
                "{root}run_20260925/code{code}.{label}.report.tsv"
            ))
            .unwrap();
            let report_rows: Vec<_> = report.lines().collect();
            let fields = std::fs::read_to_string(format!(
                "{root}fields_rules_20260925/code{code}.{label}.fields.tsv"
            ))
            .unwrap();
            let field_rows: Vec<_> = fields.lines().collect();
            let rule_trace = std::fs::read_to_string(format!(
                "{root}fields_rules_20260925/code{code}.{label}.rules.trace"
            ))
            .unwrap();
            let rule_rows: Vec<Vec<_>> = rule_trace
                .lines()
                .map(|line| line.split('\t').collect())
                .collect();
            assert_eq!(
                actual.len(),
                expected.len(),
                "code {code} {label} result count"
            );
            assert_eq!(
                actual.len(),
                report_rows.len(),
                "code {code} {label} report count"
            );
            assert_eq!(
                actual.len(),
                field_rows.len(),
                "code {code} {label} custom report count"
            );
            let heap: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
                .collect();
            assert_eq!(heap.len(), expected.len(), "code {code} {label} heap count");
            let returns: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_RETURN\t"))
                .collect();
            let edits: Vec<_> = trace
                .lines()
                .filter(|line| line.starts_with("K_TRACE_EDIT\t"))
                .collect();
            assert_eq!(
                returns.len(),
                edits.len(),
                "code {code} {label} redo scripts"
            );
            for ((((hsp, payload), expected), reported), custom_reported) in actual
                .into_iter()
                .zip(expected)
                .zip(report_rows)
                .zip(field_rows)
            {
                let f: Vec<_> = expected.split('\t').collect();
                let report_fields: Vec<_> = reported.split('\t').collect();
                let custom_fields: Vec<_> = custom_reported.split('\t').collect();
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:971-1021,
                // 1090-1095: the diagnostic field list exposes every
                // report count on the same retained HSP ordering.
                assert_eq!(custom_fields.len(), 16, "code {code} {label} custom width");
                assert_eq!(
                    hsp.hsp.score,
                    custom_fields[2].parse().unwrap(),
                    "code {code} {label} custom score"
                );
                assert_eq!(
                    payload.report_num_ident,
                    custom_fields[5].parse().unwrap(),
                    "code {code} {label} report identities"
                );
                assert_eq!(
                    payload.align_length,
                    custom_fields[7].parse().unwrap(),
                    "code {code} {label} custom length"
                );
                assert_eq!(
                    payload.report_mismatches,
                    custom_fields[8].parse().unwrap(),
                    "code {code} {label} custom mismatches"
                );
                assert_eq!(
                    payload.gap_letters,
                    custom_fields[9].parse().unwrap(),
                    "code {code} {label} gap letters"
                );
                assert_eq!(
                    payload.gap_opens,
                    custom_fields[10].parse().unwrap(),
                    "code {code} {label} custom gap opens"
                );
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:971-1021,
                // 1090-1095: report counts are recomputed from the aligned
                // sequence strings after the result HSPs are retained.
                assert_eq!(report_fields.len(), 12, "code {code} {label} report width");
                assert_eq!(
                    payload.align_length,
                    report_fields[3].parse().unwrap(),
                    "code {code} {label} alignment length"
                );
                assert_eq!(
                    payload.report_mismatches,
                    report_fields[4].parse().unwrap(),
                    "code {code} {label} mismatches"
                );
                assert_eq!(
                    payload.gap_opens,
                    report_fields[5].parse().unwrap(),
                    "code {code} {label} gap opens"
                );
                let reported_identity: f64 = report_fields[2].parse().unwrap();
                let actual_identity =
                    100.0 * payload.report_num_ident as f64 / payload.align_length as f64;
                assert!(
                    (actual_identity - reported_identity).abs() <= 0.0005,
                    "code {code} {label} percent identity"
                );
                let coords = (
                    hsp.hsp.q_start,
                    hsp.hsp.q_end,
                    hsp.hsp.s_start,
                    hsp.hsp.s_end,
                );
                assert_eq!(
                    hsp.hsp.score,
                    f[4].parse().unwrap(),
                    "code {code} {label} score"
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    f[5].parse::<f64>().unwrap().to_bits(),
                    "code {code} {label} E-value"
                );
                assert_eq!(
                    coords,
                    (
                        f[6].parse().unwrap(),
                        f[7].parse().unwrap(),
                        f[8].parse().unwrap(),
                        f[9].parse().unwrap()
                    ),
                    "code {code} {label} coordinates/order"
                );
                assert_eq!(
                    (hsp.hsp.frame, hsp.num),
                    (f[10].parse().unwrap(), f[11].parse().unwrap()),
                    "code {code} {label} frame/link count"
                );
                let heap_row: Vec<_> = heap
                    .iter()
                    .map(|line| line.split('\t').collect::<Vec<_>>())
                    .find(|row| {
                        row[3].parse::<i32>().unwrap() == hsp.hsp.score
                            && (
                                row[9].parse::<i32>().unwrap(),
                                row[10].parse::<i32>().unwrap(),
                                row[11].parse::<i32>().unwrap(),
                                row[12].parse::<i32>().unwrap(),
                            ) == coords
                    })
                    .unwrap_or_else(|| panic!("code {code} {label} missing heap HSP"));
                assert_eq!(
                    payload.bit_score.to_bits(),
                    heap_row[4].parse::<f64>().unwrap().to_bits(),
                    "code {code} {label} bits"
                );
                assert_eq!(
                    payload.num_ident,
                    heap_row[6].parse().unwrap(),
                    "code {code} {label} identities"
                );
                let matching_scripts: Vec<_> = returns
                    .iter()
                    .zip(&edits)
                    .filter(|(returned, edit)| {
                        let r: Vec<_> = returned.split('\t').collect();
                        let e: Vec<_> = edit.split('\t').collect();
                        assert_eq!(r[2], e[1], "code {code} {label} redo event identity");
                        if e[2] == "-1"
                            || (
                                r[5].parse::<i32>().unwrap(),
                                r[6].parse::<i32>().unwrap(),
                                r[7].parse::<i32>().unwrap(),
                                r[8].parse::<i32>().unwrap(),
                            ) != coords
                        {
                            return false;
                        }
                        let script: Vec<_> = e[3..]
                            .iter()
                            .map(|op| {
                                let (kind, len) = op.split_once(':').unwrap();
                                let len = len.parse().unwrap();
                                match kind {
                                    "3" => crate::common::GapEditOp::Sub(len),
                                    "0" => crate::common::GapEditOp::Del(len),
                                    "6" => crate::common::GapEditOp::Ins(len),
                                    _ => panic!("unknown NCBI edit operation"),
                                }
                            })
                            .collect();
                        payload.edit_script == script
                    })
                    .collect();
                assert_eq!(
                    matching_scripts.len(),
                    1,
                    "code {code} {label} unique redo script owner"
                );
                let returned: Vec<_> = matching_scripts[0].0.split('\t').collect();
                // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:103-124;
                // core/blast_kappa.c:325-342: BlastCompo_AlignmentNew stores
                // the raw matrix rule, then conversion moves it with the HSP.
                let rules: Vec<_> = rule_rows
                    .iter()
                    .filter(|row| {
                        row[1] == returned[4]
                            && row[4] == returned[5]
                            && row[5] == returned[6]
                            && row[6] == returned[7]
                            && row[7] == returned[8]
                            && row[8].parse::<i32>().unwrap() == i32::from(hsp.hsp.frame)
                    })
                    .collect();
                assert_eq!(rules.len(), 1, "code {code} {label} unique raw matrix rule");
                assert_eq!(
                    payload.matrix_adjust_rule as i32,
                    rules[0][2].parse().unwrap(),
                    "code {code} {label} matrix rule owner"
                );
            }
            compared += 1;
        }
    }
    assert_eq!(compared, 54);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-358,3687-3736;
// c++/src/algo/blast/core/blast_hits.c:3243-3297:
// A broad 112-target local run retains every significant redone HSP,
// including both distinct scripts in the two-HSP OID 5 list.
#[test]
fn broad_112_subject_results_keep_every_ncbi_report_payload() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/kappa_heap_rejection_20260925/"
    );
    let read_fasta = |name: &str| -> Vec<Vec<u8>> {
        let mut records: Vec<Vec<u8>> = Vec::new();
        for line in std::fs::read_to_string(format!("{root}{name}"))
            .unwrap()
            .lines()
        {
            if line.starts_with('>') {
                records.push(Vec::new());
            } else {
                records.last_mut().unwrap().extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("subjects.fna");
    let seg = SegParams::default();
    let results = run_local_mode2_sum_stats(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: Some(&seg),
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10.0,
            max_target_seqs: 112,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
    )
    .unwrap();
    let expected = std::fs::read_to_string(format!(
        "{root}report_payload_20260925/report_fields_all.tsv"
    ))
    .unwrap();
    let actual: Vec<_> = results[0]
        .lists()
        .iter()
        .flat_map(|list| {
            list.hsps
                .hsps
                .iter()
                .zip(&list.payloads)
                .map(move |(hsp, payload)| {
                    (
                        list.oid,
                        hsp.hsp.score,
                        payload.num_ident,
                        payload.num_positives,
                        payload.align_length,
                        payload.mismatches,
                        payload.gap_letters,
                        payload.gap_opens,
                        hsp.hsp.frame,
                    )
                })
        })
        .collect();
    let expected: Vec<_> = expected
        .lines()
        .map(|row| {
            let f: Vec<_> = row.split('\t').collect();
            (
                if let Some(id) = f[1].strip_prefix("weak_") {
                    id.parse::<i32>().unwrap() - 7
                } else {
                    f[1].strip_prefix("background_")
                        .unwrap()
                        .parse::<i32>()
                        .unwrap()
                        + 12
                },
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[6].parse().unwrap(),
                f[7].parse().unwrap(),
                f[8].parse().unwrap(),
                f[13].parse().unwrap(),
            )
        })
        .collect();
    assert_eq!(actual, expected);
    assert_eq!(
        results[0]
            .lists()
            .iter()
            .map(|list| list.payloads.len())
            .sum::<usize>(),
        19
    );
    assert_eq!(
        results[0]
            .lists()
            .iter()
            .find(|list| list.oid == 5)
            .unwrap()
            .payloads
            .len(),
        2
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2510;
    // c++/src/algo/blast/core/blast_hits.c:1437-1456,3243-3297:
    // The max-112 run retains each of the 12 alignments recorded by the
    // independent max-2 redo trace, including the one max-2 heap rejected.
    let redo_trace =
        std::fs::read_to_string(format!("{root}result_order_20260925/ncbi.trace")).unwrap();
    let returns: Vec<_> = redo_trace
        .lines()
        .filter(|line| line.starts_with("K_TRACE_RETURN\t"))
        .collect();
    let edits: Vec<_> = redo_trace
        .lines()
        .filter(|line| line.starts_with("K_TRACE_EDIT\t"))
        .collect();
    assert_eq!((returns.len(), edits.len()), (12, 12));
    let final_scripts: Vec<_> = results[0]
        .lists()
        .iter()
        .flat_map(|list| {
            list.hsps
                .hsps
                .iter()
                .zip(&list.payloads)
                .map(|(hsp, payload)| {
                    (
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end,
                        payload.edit_script.as_slice(),
                    )
                })
        })
        .collect();
    for (returned, edit) in returns.iter().zip(edits) {
        let r: Vec<_> = returned.split('\t').collect();
        let e: Vec<_> = edit.split('\t').collect();
        assert_eq!(r[2], e[1], "redo event identity");
        let script: Vec<_> = e[3..]
            .iter()
            .map(|op| {
                let (kind, len) = op.split_once(':').unwrap();
                let len = len.parse().unwrap();
                match kind {
                    "3" => crate::common::GapEditOp::Sub(len),
                    "0" => crate::common::GapEditOp::Del(len),
                    "6" => crate::common::GapEditOp::Ins(len),
                    _ => panic!("unexpected NCBI edit operation"),
                }
            })
            .collect();
        let coords = (
            r[5].parse::<i32>().unwrap(),
            r[6].parse::<i32>().unwrap(),
            r[7].parse::<i32>().unwrap(),
            r[8].parse::<i32>().unwrap(),
        );
        assert!(
            final_scripts
                .iter()
                .any(|&(q0, q1, s0, s1, ops)| (q0, q1, s0, s1) == coords && ops == script),
            "NCBI redo event {} lost its owned script",
            r[2]
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-905;
// c++/src/algo/blast/core/blast_kappa.c:3525-3736,2494-2515;
// c++/src/algo/blast/core/blast_traceback.c:1763-1776:
// These local -subject traces cover thirteen subjects and hard SEG with
// the same setup, Kappa and result-post-pipe order.
#[test]
fn thirteen_subject_and_hard_seg_local_results_match_exact_ncbi_heap_fields() {
    for (case, expected_count) in [("run_20260923", 11), ("seg_hard_query_20260924", 1)] {
        let input_root = format!(
            "{}/../docs/evidence/tlosan_stage_c/{case}",
            env!("CARGO_MANIFEST_DIR")
        );
        let read_fasta = |name: &str| -> Vec<(String, Vec<u8>)> {
            let mut records: Vec<(String, Vec<u8>)> = Vec::new();
            for line in std::fs::read_to_string(format!("{input_root}/{name}"))
                .unwrap()
                .lines()
            {
                if let Some(id) = line.strip_prefix('>') {
                    records.push((id.split_whitespace().next().unwrap().to_owned(), Vec::new()));
                } else {
                    records.last_mut().unwrap().1.extend(line.bytes());
                }
            }
            records
        };
        let query = read_fasta("query.faa");
        let subjects = read_fasta("subjects.fna");
        let query_sequences: Vec<_> = query.iter().map(|(_, seq)| seq.clone()).collect();
        let subject_sequences: Vec<_> = subjects.iter().map(|(_, seq)| seq.clone()).collect();
        let seg = SegParams::default();
        let results = run_local_mode2_sum_stats(
            &query_sequences,
            &subject_slices(&subject_sequences),
            LocalStageDProfile {
                max_intron_length: 0,
                seg: Some(&seg),
                soft_masking: false,
                mask_lowercase: false,
                genetic_code: 1,
                expect_value: 10.0,
                max_target_seqs: 500,
                query_records: &[],
                subject_set: Default::default(),
                engine_positives: false,
                culling_limit: 0,
                best_hit: None,
                subject_besthit: false,
                pssm: None,
                fetch_excluded: None,
                query_cov_hsp_perc: 0.0,
            },
        )
        .unwrap();
        let trace_path = format!(
            "{}/../docs/evidence/tlosan_stage_d/kappa_result_order_20260925/{case}_default.tsv",
            env!("CARGO_MANIFEST_DIR")
        );
        let trace = std::fs::read_to_string(trace_path).unwrap();
        let expected: std::collections::HashMap<(i32, usize), Vec<&str>> = trace
            .lines()
            .filter(|line| line.starts_with("K_TRACE_HEAP_HSP\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                ((f[1].parse().unwrap(), f[2].parse().unwrap()), f)
            })
            .collect();
        assert_eq!(expected.len(), expected_count);
        let lists = results[0].lists();
        let actual_count: usize = lists.iter().map(|list| list.hsps.hsps.len()).sum();
        assert_eq!(actual_count, expected_count, "{case} HSP count");
        for list in lists {
            for (index, (hsp, payload)) in list.hsps.hsps.iter().zip(&list.payloads).enumerate() {
                let row = &expected[&(list.oid, index)];
                assert_eq!(hsp.hsp.score, row[3].parse().unwrap(), "{case} score");
                assert_eq!(
                    payload.bit_score.to_bits(),
                    row[4].parse::<f64>().unwrap().to_bits(),
                    "{case} bits"
                );
                assert_eq!(
                    hsp.evalue.to_bits(),
                    row[5].parse::<f64>().unwrap().to_bits(),
                    "{case} E-value"
                );
                assert_eq!(
                    payload.num_ident,
                    row[6].parse().unwrap(),
                    "{case} identity"
                );
                assert_eq!(hsp.hsp.frame, row[8].parse().unwrap(), "{case} frame");
                assert_eq!(
                    (
                        hsp.hsp.q_start,
                        hsp.hsp.q_end,
                        hsp.hsp.s_start,
                        hsp.hsp.s_end
                    ),
                    (
                        row[9].parse().unwrap(),
                        row[10].parse().unwrap(),
                        row[11].parse().unwrap(),
                        row[12].parse().unwrap()
                    ),
                    "{case} coordinates"
                );
            }
        }
        let output_path = format!(
            "{}/../docs/evidence/tlosan_stage_d/run_20260924/{case}_default.out",
            env!("CARGO_MANIFEST_DIR")
        );
        let output = std::fs::read_to_string(output_path).unwrap();
        let actual_order: Vec<_> = lists
            .iter()
            .map(|list| subjects[list.oid as usize].0.as_str())
            .collect();
        let expected_order: Vec<_> = output
            .lines()
            .map(|line| line.split('\t').nth(1).unwrap())
            .collect();
        assert_eq!(actual_order, expected_order, "{case} result order");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1437-1456,3243-3297;
// c++/src/algo/blast/core/blast_kappa.c:305-358,3687-3736:
// Natural cross-frame sum-statistics HSPs leave Kappa in score order
// 238,189,135 but hitlist replacement E-value-sorts them 238,135,189.
#[test]
fn natural_evalue_sort_moves_the_owned_kappa_payload() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/natural_positive_20260925/"
    );
    let read_fasta = |name: &str| -> Vec<Vec<u8>> {
        let mut records: Vec<Vec<u8>> = Vec::new();
        for line in std::fs::read_to_string(format!("{root}{name}"))
            .unwrap()
            .lines()
        {
            if line.starts_with('>') {
                records.push(Vec::new());
            } else {
                records.last_mut().unwrap().extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("evalue_subjects.fna");
    let results = run_local_mode2_sum_stats(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: None,
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10000.0,
            max_target_seqs: 1,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
    )
    .unwrap();
    let lists = results[0].lists();
    assert_eq!(lists.len(), 1);
    assert_eq!(lists[0].oid, 0);
    assert_eq!(
        lists[0]
            .hsps
            .hsps
            .iter()
            .map(|h| h.hsp.score)
            .collect::<Vec<_>>(),
        [238, 135, 189]
    );
    let trace =
        std::fs::read_to_string(format!("{root}run_20260925/evalue_sort.kappa.trace")).unwrap();
    let heap: std::collections::HashMap<i32, Vec<&str>> = trace
        .lines()
        .filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t0\t"))
        .map(|row| {
            let f: Vec<_> = row.split('\t').collect();
            (f[3].parse().unwrap(), f)
        })
        .collect();
    for (hsp, payload) in lists[0].hsps.hsps.iter().zip(&lists[0].payloads) {
        let row = &heap[&hsp.hsp.score];
        assert_eq!(
            hsp.evalue.to_bits(),
            row[5].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            payload.bit_score.to_bits(),
            row[4].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(payload.num_ident, row[6].parse().unwrap());
        assert_eq!(
            (
                hsp.hsp.q_start,
                hsp.hsp.q_end,
                hsp.hsp.s_start,
                hsp.hsp.s_end
            ),
            (
                row[9].parse().unwrap(),
                row[10].parse().unwrap(),
                row[11].parse().unwrap(),
                row[12].parse().unwrap()
            )
        );
        let script_len = match hsp.hsp.score {
            238 => 45,
            135 | 189 => 35,
            _ => unreachable!(),
        };
        assert_eq!(
            payload.edit_script,
            vec![crate::common::GapEditOp::Sub(script_len)]
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:224-273,3658-3689:
// The pinned local fixture naturally redos two same-frame alignments;
// the smaller one is contained, leaving one postredo link input.
#[test]
fn natural_postredo_containment_keeps_surviving_script_and_report_fields() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/natural_positive_20260925/"
    );
    let read = |name: &str| -> Vec<u8> {
        std::fs::read_to_string(format!("{root}{name}"))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(str::bytes)
            .collect()
    };
    let queries = vec![read("containment_query.faa")];
    let subjects = vec![read("containment_subject.fna")];
    let results = run_local_mode2_sum_stats(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: None,
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10000.0,
            max_target_seqs: 1,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
    )
    .unwrap();
    let lists = results[0].lists();
    assert_eq!(lists.len(), 1);
    assert_eq!(lists[0].hsps.hsps.len(), 1);
    let hsp = &lists[0].hsps.hsps[0];
    let payload = &lists[0].payloads[0];
    let trace =
        std::fs::read_to_string(format!("{root}containment_run_20260925/kappa.trace")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .find(|row| row.starts_with("K_TRACE_HEAP_HSP\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(hsp.hsp.score, expected[3].parse().unwrap());
    assert_eq!(
        payload.bit_score.to_bits(),
        expected[4].parse::<f64>().unwrap().to_bits()
    );
    assert_eq!(
        hsp.evalue.to_bits(),
        expected[5].parse::<f64>().unwrap().to_bits()
    );
    assert_eq!(payload.num_ident, expected[6].parse().unwrap());
    assert_eq!(
        (
            hsp.hsp.frame,
            hsp.hsp.q_start,
            hsp.hsp.q_end,
            hsp.hsp.s_start,
            hsp.hsp.s_end
        ),
        (
            expected[8].parse().unwrap(),
            expected[9].parse().unwrap(),
            expected[10].parse().unwrap(),
            expected[11].parse().unwrap(),
            expected[12].parse().unwrap()
        )
    );
    assert_eq!(
        payload.edit_script,
        vec![
            crate::common::GapEditOp::Sub(18),
            crate::common::GapEditOp::Ins(5),
            crate::common::GapEditOp::Sub(97),
        ]
    );
    assert_eq!(payload.matrix_adjust_rule as i32, 4);
    let report =
        std::fs::read_to_string(format!("{root}containment_run_20260925/report.tsv")).unwrap();
    let f: Vec<_> = report.trim_end().split('\t').collect();
    assert_eq!(
        (
            payload.num_ident,
            payload.num_positives,
            payload.align_length,
            payload.mismatches,
            payload.gap_letters,
            payload.gap_opens
        ),
        (
            f[5].parse().unwrap(),
            f[6].parse().unwrap(),
            f[7].parse().unwrap(),
            f[8].parse().unwrap(),
            f[9].parse().unwrap(),
            f[10].parse().unwrap()
        )
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/compo_heap.c:330-391;
// c++/src/algo/blast/core/blast_kappa.c:3719-3755:
// Two naturally generated weak subjects make the second Kappa record
// replace the first at heap threshold one, both above inclusion 0.002.
#[test]
fn natural_composition_heap_replacement_keeps_new_hsp_payloads() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_d/natural_positive_20260925/"
    );
    let read_fasta = |name: &str| -> Vec<Vec<u8>> {
        let mut records: Vec<Vec<u8>> = Vec::new();
        for line in std::fs::read_to_string(format!("{root}{name}"))
            .unwrap()
            .lines()
        {
            if line.starts_with('>') {
                records.push(Vec::new());
            } else {
                records.last_mut().unwrap().extend(line.bytes());
            }
        }
        records
    };
    let queries = read_fasta("query.faa");
    let subjects = read_fasta("heap_subjects.fna");
    let results = run_local_mode2_sum_stats(
        &queries,
        &subject_slices(&subjects),
        LocalStageDProfile {
            max_intron_length: 0,
            seg: None,
            soft_masking: false,
            mask_lowercase: false,
            genetic_code: 1,
            expect_value: 10000.0,
            max_target_seqs: 1,
            query_records: &[],
            subject_set: Default::default(),
            engine_positives: false,
            culling_limit: 0,
            best_hit: None,
            subject_besthit: false,
            pssm: None,
            fetch_excluded: None,
            query_cov_hsp_perc: 0.0,
        },
    )
    .unwrap();
    let lists = results[0].lists();
    assert_eq!(lists.len(), 1);
    assert_eq!(lists[0].oid, 0);
    let trace = std::fs::read_to_string(format!("{root}run_20260925/heap_replacement.kappa.trace"))
        .unwrap();
    assert!(trace
        .lines()
        .any(|row| row.starts_with("K_TRACE_HEAP_INSERT_RETURN\t0\t") && row.ends_with("\t1")));
    let rows: Vec<_> = trace
        .lines()
        .filter(|row| row.starts_with("K_TRACE_HEAP_HSP\t0\t"))
        .map(|row| row.split('\t').collect::<Vec<_>>())
        .collect();
    assert_eq!(lists[0].hsps.hsps.len(), rows.len());
    for ((hsp, payload), row) in lists[0].hsps.hsps.iter().zip(&lists[0].payloads).zip(rows) {
        assert_eq!(hsp.hsp.score, row[3].parse().unwrap());
        assert_eq!(
            hsp.evalue.to_bits(),
            row[5].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(
            payload.bit_score.to_bits(),
            row[4].parse::<f64>().unwrap().to_bits()
        );
        assert_eq!(payload.num_ident, row[6].parse().unwrap());
        assert!(!payload.edit_script.is_empty());
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:697
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:499
// The budget of a chunk: the largest file, at most 2^30 bytes, over four times the threads,
// plus one. For `-subject` the file is ceil(bases / 4), the `.nsq` of `makeblastdb`.
#[cfg(feature = "parallel")]
#[test]
fn x_chunk_budget_is_the_slice_over_four_times_the_threads() {
    for (threads, largest_file, bases) in [
        // the fly genome, 143.7 Mb: a .nsq of 35.9 MB
        (8, 35_925_000, 1_122_657),
        (2, 35_925_000, 4_490_626),
        // the human genome, 3,100 Mb: 775 MB, below the slice limit
        (8, 775_000_000, 24_218_751),
        // above 2^30 bytes the slice is 2^30
        (8, 2_000_000_000, 33_554_433),
        (32, 1 << 30, 8_388_609),
        (8, 0, 1),
        (0, 100, 26),
        (usize::MAX, 1 << 30, 1),
    ] {
        assert_eq!(
            x_oid_chunk_bases(threads, largest_file),
            bases,
            "{threads} threads, largest file {largest_file}"
        );
    }
    for (lengths, size) in [
        (&[][..], 0),
        (&[0][..], 0),
        (&[1][..], 1),
        (&[4][..], 1),
        (&[5][..], 2),
        (&[3, 3][..], 2),
        (&[1_000; 240][..], 60_000),
    ] {
        let subjects: Vec<Vec<u8>> = lengths.iter().map(|&n| vec![b'A'; n]).collect();
        assert_eq!(
            x_subject_file_size(SubjectLetters::Records(&subject_slices(&subjects))),
            size,
            "{lengths:?}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:699-707
// A chunk is at least one subject, then the next ones while the bases stay within the budget
// (a sum equal to the budget is in); a subject longer than the budget is a chunk of its own,
// and empty subjects join a chunk that is not over its budget.
#[cfg(feature = "parallel")]
#[test]
fn x_chunks_hold_consecutive_subjects_within_the_budget() {
    let set =
        |lengths: &[usize]| -> Vec<Vec<u8>> { lengths.iter().map(|&n| vec![b'A'; n]).collect() };
    for (lengths, begin, bases, end) in [
        // the budget's edge
        (&[10, 10, 10][..], 0, 20, 2),
        (&[10, 10, 10][..], 0, 19, 1),
        (&[10, 10, 10][..], 0, 30, 3),
        (&[10, 10, 10][..], 1, 20, 3),
        // one subject over the budget, first or after shorter ones
        (&[50, 1, 1][..], 0, 20, 1),
        (&[5, 50, 5][..], 0, 20, 1),
        (&[5, 50, 5][..], 1, 20, 2),
        (&[5, 50, 5][..], 2, 20, 3),
        // empty subjects
        (&[0, 0, 0][..], 0, 1, 3),
        (&[30, 0, 0][..], 0, 20, 1),
        (&[10, 0, 10, 0, 1][..], 0, 20, 4),
        (&[0, 30][..], 0, 20, 1),
        // no subject at `begin`
        (&[10][..], 1, 20, 1),
    ] {
        assert_eq!(
            x_oid_chunk_end(
                SubjectLetters::Records(&subject_slices(&set(lengths))),
                begin,
                bases
            ),
            end,
            "{lengths:?} from {begin}, {bases} bases"
        );
    }
    // `x_all_chunks` takes every chunk from the first subject on, in OID order.
    let subjects = set(&[50, 5, 5, 5, 5, 50, 5]);
    assert_eq!(
        x_all_chunks(SubjectLetters::Records(&subject_slices(&subjects)), 12),
        [0..1, 1..3, 3..5, 5..6, 6..7]
    );
    assert!(x_all_chunks(SubjectLetters::Records(&[]), 12).is_empty());
    // 240 subjects of 1,000 bases at 4 threads (tests/tblastn_oid_chunk.rs): 3,751 bases, three
    // subjects per chunk, 80 chunks in all.
    let subjects = set(&[1_000; 240]);
    let bases = x_oid_chunk_bases(
        4,
        x_subject_file_size(SubjectLetters::Records(&subject_slices(&subjects))),
    );
    assert_eq!(bases, 3_751);
    let all = x_all_chunks(SubjectLetters::Records(&subject_slices(&subjects)), bases);
    assert_eq!(all.len(), 80);
    assert_eq!(all[..2], [0..3, 3..6]);
    assert_eq!(all[79], 237..240);
}

// `XChunkCut` cuts, one length after the other, the chunks that `x_oid_chunk_end` cuts from
// each chunk's first subject (NCBI's loop), on the lengths of every edge above and on random
// ones; `x_all_chunks` of a `-db` stream cuts the index's lengths and keeps them.
#[cfg(feature = "parallel")]
#[test]
fn x_chunk_cut_gives_the_chunks_of_each_first_subject() {
    let reference = |lengths: &[usize], bases: u64| -> Vec<u32> {
        let records: Vec<Vec<u8>> = lengths.iter().map(|&n| vec![b'A'; n]).collect();
        let slices = subject_slices(&records);
        let subjects = SubjectLetters::Records(&slices);
        let mut ends = Vec::new();
        let mut start = 0;
        while start < lengths.len() {
            start = x_oid_chunk_end(subjects, start, bases);
            ends.push(start as u32);
        }
        ends
    };
    let cut = |lengths: &[usize], bases: u64| -> Vec<u32> {
        let mut cut = XChunkCut::new(bases);
        for &length in lengths {
            cut.push(length as u64);
        }
        cut.ends()
    };
    let mut state = 0x2026_1103_u64;
    let mut next = |below: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % below
    };
    let mut cases: Vec<(Vec<usize>, u64)> = vec![
        (vec![], 12),
        (vec![10, 10, 10], 20),
        (vec![10, 10, 10], 19),
        (vec![50, 1, 1], 20),
        (vec![5, 50, 5], 20),
        (vec![0, 0, 0], 1),
        (vec![30, 0, 0], 20),
        (vec![10, 0, 10, 0, 1], 20),
        (vec![0, 30], 20),
    ];
    for _ in 0..500 {
        let count = next(40) as usize;
        let lengths = (0..count)
            .map(|_| match next(4) {
                0 => 0,
                1 => next(5) as usize,
                _ => next(300) as usize,
            })
            .collect();
        cases.push((lengths, 1 + next(400)));
    }
    for (lengths, bases) in cases {
        assert_eq!(
            cut(&lengths, bases),
            reference(&lengths, bases),
            "{lengths:?} {bases}"
        );
    }
    // A sum past `u64` starts the next chunk.
    let mut big = XChunkCut::new(u64::MAX);
    for length in [u64::MAX, 1, 0] {
        big.push(length);
    }
    assert_eq!(big.ends(), [1, 3]);
    // A `-db` stream: the chunks of the index's lengths, cut once for each budget.
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/db/nucl_mv");
    let db = crate::seqdb::SeqDb::open(path.to_str().unwrap(), crate::seqdb::SeqType::Nucleotide)
        .unwrap();
    let stream = crate::algorithm::common::x_db_stream::XDbStream::new(&db, 4, "TBLASTN").unwrap();
    let mut approx = Vec::new();
    stream.x_approx_lengths(|length| approx.push(length as usize));
    for bases in [1, 50, 500, 5_000, u64::MAX] {
        let chunks = x_all_chunks(SubjectLetters::Stream(&stream), bases);
        let ends: Vec<u32> = chunks.iter().map(|chunk| chunk.end as u32).collect();
        assert_eq!(ends, reference(&approx, bases), "{bases}");
        assert_eq!(chunks.first().map(|chunk| chunk.start), Some(0));
        assert_eq!(chunks.last().map(|chunk| chunk.end), Some(stream.count()));
        assert_eq!(x_all_chunks(SubjectLetters::Stream(&stream), bases), chunks);
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-289
// NCBI searches every batch with the same code, whatever their number. A search of fewer than
// four batches keeps no subject plans (XSubjectSide::KEEP_FROM_BATCHES), and that does not
// move it to the other loop: x_search_batches takes the batches whatever their number, at one
// thread too, and gives the report of run_local_for_report_threads, with the plans kept (four
// batches) and without (one, two, three).
#[test]
fn x_search_batches_takes_the_batches_whatever_their_number() {
    // Test protein: made from the first residues of E. coli AsnA (RefSeq NP_418200.1) and other letters.
    let queries = vec![b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ".to_vec(); 4];
    let subjects: Vec<Vec<u8>> = (0..2usize)
        .map(|n| {
            (0..3_000usize)
                .map(|i| b"ACGT"[(i * 7 + i / 3 + i / 11 + n) % 4])
                .collect()
        })
        .collect();
    let seg = SegParams::default();
    let profile = |_: std::ops::Range<usize>| LocalStageDProfile {
        max_intron_length: 0,
        seg: Some(&seg),
        soft_masking: false,
        mask_lowercase: false,
        genetic_code: 1,
        expect_value: 10.0,
        max_target_seqs: 500,
        query_records: &[],
        subject_set: Default::default(),
        engine_positives: false,
        culling_limit: 0,
        best_hit: None,
        subject_besthit: false,
        pssm: None,
        fetch_excluded: None,
        query_cov_hsp_perc: 0.0,
    };
    for ranges in [
        vec![0..1],
        vec![0..1, 1..2],
        vec![0..1, 1..2, 2..3],
        vec![0..1, 1..2, 2..3, 3..4],
    ] {
        let searched = x_search_batches(
            &queries,
            SubjectLetters::Records(&subject_slices(&subjects)),
            ranges.clone(),
            profile,
            BlastCompoAdjustMode::CompositionMatrixAdjust,
            true,
            LocalStageDScoring::default(),
            1,
        )
        .unwrap();
        assert_eq!(searched.as_ref().map(Vec::len), Some(ranges.len()));
        for (range, (results, parameters, _, validity)) in searched.into_iter().flatten() {
            let (reference, reference_parameters, _, reference_validity) =
                run_local_for_report_threads(
                    &queries[range.clone()],
                    SubjectLetters::Records(&subject_slices(&subjects)),
                    profile(range.clone()),
                    BlastCompoAdjustMode::CompositionMatrixAdjust,
                    true,
                    LocalStageDScoring::default(),
                    1,
                )
                .unwrap();
            assert!(ranges.contains(&range));
            assert_eq!(parameters, reference_parameters);
            assert_eq!(validity, reference_validity);
            assert_eq!(results.len(), reference.len());
            for (hitlist, expected) in results.iter().zip(&reference) {
                assert_eq!(hitlist.lists().len(), expected.lists().len());
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:627-631
// NCBI translates a subject for each batch and frees it after the subject. Batches searched
// side by side on five subjects at two threads (five subjects hold more than two of the
// longest, so only the side-by-side layout keeps the plans, x_searched_side_by_side): each
// batch's report is the one it has alone.
#[cfg(feature = "parallel")]
#[test]
fn x_search_batches_side_by_side_gives_each_batch_its_report() {
    let queries = vec![b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ".to_vec(); 4];
    let subjects: Vec<Vec<u8>> = (0..5usize)
        .map(|n| {
            (0..3_000usize)
                .map(|i| b"ACGT"[(i * 7 + i / 3 + i / 11 + n) % 4])
                .collect()
        })
        .collect();
    let lengths = || subjects.iter().map(Vec::len);
    assert!(
        !super::super::search_seed::XSubjectSide::kept_within_threads(
            lengths().max().unwrap(),
            lengths().sum(),
            2
        )
    );
    assert!(
        x_searched_side_by_side(SubjectLetters::Records(&subject_slices(&subjects)), 4, 2)
            || !x_tbn_batch_rule()
    );
    assert_eq!(
        x_batches_side_by_side(lengths().sum(), subjects.len(), 4, 2),
        4
    );
    let seg = SegParams::default();
    let profile = |_: std::ops::Range<usize>| LocalStageDProfile {
        max_intron_length: 0,
        seg: Some(&seg),
        soft_masking: false,
        mask_lowercase: false,
        genetic_code: 1,
        expect_value: 10.0,
        max_target_seqs: 500,
        query_records: &[],
        subject_set: Default::default(),
        engine_positives: false,
        culling_limit: 0,
        best_hit: None,
        subject_besthit: false,
        pssm: None,
        fetch_excluded: None,
        query_cov_hsp_perc: 0.0,
    };
    let ranges = vec![0..1, 1..2, 2..3, 3..4];
    let searched = x_search_batches(
        &queries,
        SubjectLetters::Records(&subject_slices(&subjects)),
        ranges.clone(),
        profile,
        BlastCompoAdjustMode::CompositionMatrixAdjust,
        true,
        LocalStageDScoring::default(),
        2,
    )
    .unwrap();
    let searched = searched.expect("two threads search the batches on the pool");
    assert_eq!(
        searched
            .iter()
            .map(|(range, _)| range.clone())
            .collect::<Vec<_>>(),
        ranges
    );
    for (range, (results, parameters, _, validity)) in searched {
        let (reference, reference_parameters, _, reference_validity) =
            run_local_for_report_threads(
                &queries[range.clone()],
                SubjectLetters::Records(&subject_slices(&subjects)),
                profile(range.clone()),
                BlastCompoAdjustMode::CompositionMatrixAdjust,
                true,
                LocalStageDScoring::default(),
                1,
            )
            .unwrap();
        assert_eq!(parameters, reference_parameters);
        assert_eq!(validity, reference_validity);
        assert_eq!(results.len(), reference.len());
        for (hitlist, expected) in results.iter().zip(&reference) {
            let lists = |hitlist: &KappaResultHitList| {
                hitlist
                    .lists()
                    .iter()
                    .map(|list| (list.oid, list.payloads.len()))
                    .collect::<Vec<_>>()
            };
            assert_eq!(lists(hitlist), lists(expected));
            assert_eq!(
                hitlist.worst_evalue().to_bits(),
                expected.worst_evalue().to_bits()
            );
        }
    }
}

// EXPERIMENT (LOSAT_X_TBNKEEPFEW): subjects searched by batches side by side or on own pools keep
// their plans; batches all one after another, or one thread, do not.
#[cfg(feature = "parallel")]
#[test]
fn x_plans_are_kept_for_batches_side_by_side() {
    let subjects: Vec<Vec<u8>> = (0..100usize)
        .map(|n| (0..3_000usize).map(|i| b"ACGT"[(i * 7 + n) % 4]).collect())
        .collect();
    let slices = subject_slices(&subjects);
    let side_by_side = |batches, threads| {
        x_searched_side_by_side(SubjectLetters::Records(&slices), batches, threads)
    };
    if !x_tbn_batch_rule() {
        return;
    }
    assert!(side_by_side(8, 4));
    assert!(side_by_side(7, 4));
    assert!(side_by_side(8, 8));
    assert!(side_by_side(16, 8));
    assert!(!side_by_side(7, 8));
    assert!(!side_by_side(8, 1));
    // A set of more than 24 Mb of many subjects is searched one batch after another.
    let large: Vec<&[u8]> = std::iter::repeat_n(&subjects[0][..], 8_001).collect();
    assert!(!x_searched_side_by_side(
        SubjectLetters::Records(&large),
        16,
        4
    ));
    // Ten long subjects (above 24 Mb) at 8 threads: four batches, none side by side, on four
    // own pools (natively).
    let long: Vec<u8> = (0..2_500_001usize).map(|i| b"ACGT"[i % 4]).collect();
    let ten: Vec<&[u8]> = std::iter::repeat_n(&long[..], 10).collect();
    assert_eq!(x_side_by_side_batches(25_000_010, 10, 4, 8, true), 0);
    assert_eq!(
        x_searched_side_by_side(SubjectLetters::Records(&ten), 4, 8),
        x_tbn_own_pools()
    );
}

// EXPERIMENT (LOSAT_X_TBNBATCHRULE): the full waves of a small subject set or of fewer than
// two subjects per thread; none of a large set of many subjects or of a search on one thread.
#[cfg(feature = "parallel")]
#[test]
fn x_side_by_side_takes_the_full_waves_of_a_small_or_few_subject_set() {
    let small = X_SIDE_BY_SIDE_MAX_BASES;
    assert_eq!(x_batches_side_by_side(small, 1_000, 10, 8), 8);
    assert_eq!(x_batches_side_by_side(small, 1_000, 16, 8), 16);
    assert_eq!(x_batches_side_by_side(small, 1_000, 4, 3), 3);
    assert_eq!(x_batches_side_by_side(small, 1_000, 5, 8), 0);
    assert_eq!(x_batches_side_by_side(small, 1_000, 64, 16), 64);
    assert_eq!(x_batches_side_by_side(small, 1, 10, 1), 0);
    assert_eq!(x_batches_side_by_side(small + 1, 32, 64, 16), 0);
    assert_eq!(x_batches_side_by_side(usize::MAX, 1_000, 64, 16), 0);
    assert_eq!(x_batches_side_by_side(small + 1, 31, 40, 16), 32);
    assert_eq!(x_batches_side_by_side(small + 1, 32, 40, 16), 0);
    assert_eq!(x_batches_side_by_side(usize::MAX, 3, 20, 8), 16);
}

// EXPERIMENT (LOSAT_X_TBNOWNPOOLS): two to threads - 1 batches left (threads / 2 above
// X_SIDE_BY_SIDE_MAX_BASES) go to pools of their own, all or none, and only with fewer than
// two subjects per thread; the threads are shared out with the remainder first.
#[cfg(feature = "parallel")]
#[test]
fn x_own_pools_take_the_batches_left_of_few_subjects() {
    let small = X_SIDE_BY_SIDE_MAX_BASES;
    assert_eq!(x_batches_on_own_pools(16, small, 1, 32), 16);
    assert_eq!(x_batches_on_own_pools(5, small, 1, 8), 5);
    assert_eq!(x_batches_on_own_pools(7, small, 1, 8), 7);
    assert_eq!(x_batches_on_own_pools(8, small, 1, 8), 0);
    assert_eq!(x_batches_on_own_pools(2, small, 1, 3), 2);
    assert_eq!(x_batches_on_own_pools(1, small, 1, 8), 0);
    assert_eq!(x_batches_on_own_pools(0, small, 1, 8), 0);
    assert_eq!(x_batches_on_own_pools(5, small + 1, 3, 16), 5);
    assert_eq!(x_batches_on_own_pools(8, small + 1, 3, 16), 8);
    assert_eq!(x_batches_on_own_pools(9, small + 1, 3, 16), 0);
    assert_eq!(x_batches_on_own_pools(5, small + 1, 3, 8), 0);
    assert_eq!(x_batches_on_own_pools(16, small, 63, 32), 16);
    assert_eq!(x_batches_on_own_pools(16, small, 64, 32), 0);
    assert_eq!(x_batches_on_own_pools(5, small, 5_000, 16), 0);
    let shares: Vec<_> = (0..5).map(|i| x_own_pool_threads(16, 5, i)).collect();
    assert_eq!(shares, vec![4, 3, 3, 3, 3]);
    assert_eq!(
        (0..16)
            .map(|i| x_own_pool_threads(32, 16, i))
            .sum::<usize>(),
        32
    );
    assert_eq!(
        (0..7).map(|i| x_own_pool_threads(16, 7, i)).sum::<usize>(),
        16
    );
}

/// A fixed xorshift sequence.
#[cfg(feature = "parallel")]
struct XorShift(u64);

#[cfg(feature = "parallel")]
impl XorShift {
    fn below(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }
}

// EXPERIMENT (LOSAT_X_TBNREDOLPT): the order names each wanted match once, by HSPs largest
// first, ties in index order, and no match early termination ended.
#[cfg(feature = "parallel")]
#[test]
fn the_redo_lpt_order_takes_the_wanted_matches_largest_first() {
    assert!(x_redo_lpt_order(&[], |_| 0).is_empty());
    assert!(x_redo_lpt_order(&[false, false], |_| 1).is_empty());
    let hsps = [3, 9, 3, 0, 9, 4];
    assert_eq!(
        x_redo_lpt_order(&[true; 6], |index| hsps[index]),
        [1, 4, 5, 0, 2, 3]
    );
    assert_eq!(
        x_redo_lpt_order(&[true, false, true, true, true, false], |index| hsps[index]),
        [4, 0, 2, 3]
    );
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    for len in [1, 2, 7, 64, 256] {
        // Few distinct HSP counts (many ties), and many.
        for distinct in [1, 3, 100_000] {
            let hsps: Vec<usize> = (0..len).map(|_| rng.below(distinct)).collect();
            let wanted: Vec<bool> = (0..len).map(|_| rng.below(4) != 0).collect();
            let order = x_redo_lpt_order(&wanted, |index| hsps[index]);
            let mut named = order.clone();
            named.sort_unstable();
            let expected: Vec<usize> = (0..len).filter(|&index| wanted[index]).collect();
            assert_eq!(named, expected, "the wanted matches, once each");
            for pair in order.windows(2) {
                let (a, b) = (pair[0], pair[1]);
                assert!(
                    hsps[a] > hsps[b] || (hsps[a] == hsps[b] && a < b),
                    "match {a} ({} HSPs) before match {b} ({} HSPs)",
                    hsps[a],
                    hsps[b]
                );
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_TBNREDOLPT): the scope returns what the `par_iter` path returns for the
// same chunk (`Some(run(match))` where wanted, `None` elsewhere, by index), runs each wanted
// match once, and gives each task one worker of its own; also when a match runs nested
// parallel work, where a thread that waits in it may run another task's loop.
#[cfg(feature = "parallel")]
#[test]
fn the_redo_lpt_scope_returns_the_par_iter_results_by_index() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    for threads in [2, 3, 7, 17] {
        crate::utils::threading::with_search_pool(threads, "tblastn", |pool| {
            for matches in [0_usize, 1, 2, 5, 40, 256] {
                let hsps: Vec<usize> = (0..matches).map(|_| rng.below(50)).collect();
                let wanted: Vec<bool> = (0..matches).map(|_| rng.below(5) != 0).collect();
                let order = x_redo_lpt_order(&wanted, |index| hsps[index]);
                // Each worker records the matches it ran.
                let workers: Vec<std::sync::Mutex<Vec<usize>>> =
                    (0..threads).map(|_| Default::default()).collect();
                let run = |index: usize| -> (usize, u64) {
                    let nested: u64 = (0..64_u64).into_par_iter().map(|x| x * index as u64).sum();
                    (index * 3 + hsps[index], nested)
                };
                let results = x_redo_lpt(pool, &workers, matches, &order, |ran, index| {
                    assert!(rayon::current_thread_index().is_some(), "on the pool");
                    ran.push(index);
                    run(index)
                });
                let par_iter: Vec<Option<(usize, u64)>> = pool.install(|| {
                    (0..matches)
                        .into_par_iter()
                        .zip(wanted.par_iter())
                        .map(|(index, &wanted)| wanted.then(|| run(index)))
                        .collect()
                });
                assert_eq!(results, par_iter, "threads {threads} matches {matches}");
                let ran: Vec<Vec<usize>> = workers
                    .into_iter()
                    .map(|worker| worker.into_inner().unwrap())
                    .collect();
                let tasks = threads.min(order.len());
                assert!(
                    ran[tasks..].iter().all(Vec::is_empty),
                    "only the first {tasks} workers ran"
                );
                let mut all: Vec<usize> = ran.concat();
                all.sort_unstable();
                let mut expected = order.clone();
                expected.sort_unstable();
                assert_eq!(all, expected, "each wanted match once");
                // One worker runs its matches in the order's order.
                for worker in &ran {
                    let position = |index: &usize| order.iter().position(|o| o == index);
                    assert!(worker.windows(2).all(|w| position(&w[0]) < position(&w[1])));
                }
            }
            Ok(())
        })
        .expect("a pool");
    }
}

// EXPERIMENT (LOSAT_X_TBNREDOLPT): a match that panics: the panic comes out of the scope (the
// other tasks finish first, so nothing hangs) and the same pool computes the next chunk, with
// workers made again as the search would make them.
#[cfg(feature = "parallel")]
#[test]
fn a_panicking_redo_match_propagates_and_the_pool_stays_usable() {
    for threads in [2, 3, 7, 17] {
        crate::utils::threading::with_search_pool(threads, "tblastn", |pool| {
            let matches = 40_usize;
            let order: Vec<usize> = (0..matches).collect();
            let run = |panicking: bool| {
                let workers: Vec<std::sync::Mutex<()>> =
                    (0..threads).map(|_| Default::default()).collect();
                x_redo_lpt(pool, &workers, matches, &order, |_, index| {
                    assert!(
                        !(panicking && (index == 5 || index == 30)),
                        "match {index} panics"
                    );
                    index * 3
                })
            };
            let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(true)));
            assert!(caught.is_err(), "threads {threads}: the panic was lost");
            let expected: Vec<Option<usize>> = (0..matches).map(|index| Some(index * 3)).collect();
            assert_eq!(run(false), expected, "threads {threads}");
            Ok(())
        })
        .expect("a pool");
    }
}
