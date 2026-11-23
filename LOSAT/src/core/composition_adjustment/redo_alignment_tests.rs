use super::BlastCompoAlignmentContext::PreliminaryHspIndex;
use super::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:80,644-751;
// c++/src/algo/blast/core/blast_kappa.c:1498-1525:
// Compare ordered joined Rust windows to the pinned partial-translation calls.
#[test]
fn translated_kappa_windows_match_ncbi_partial_translation_ranges() {
    for (case, query_lengths) in [
        ("seg_hard_query_20260924_default", vec![160]),
        ("multi_query_20260924_default", vec![120, 70, 120]),
    ] {
        let trace = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_d/kappa_mode2_20260924/{case}.tsv",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let enter = trace
            .lines()
            .find(|line| line.starts_with("K_CALL\t0\tredo_enter\t"))
            .unwrap();
        let enter_fields: Vec<_> = enter.split('\t').collect();
        let subject_nt_length: i32 = enter_fields[7].parse().unwrap();
        let input_rows: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("K_ALIGN\t0\tincoming\t"))
            .collect();
        assert_eq!(input_rows.len(), enter_fields[3].parse::<usize>().unwrap());
        let mut input: Option<Box<BlastCompoAlignment>> = None;
        for row in input_rows.into_iter().rev() {
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(f[6], "-1");
            input = Some(Box::new(BlastCompoAlignment {
                score: f[5].parse().unwrap(),
                matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
                query_index: f[7].parse().unwrap(),
                query_start: f[8].parse().unwrap(),
                query_end: f[9].parse().unwrap(),
                match_start: f[10].parse().unwrap(),
                match_end: f[11].parse().unwrap(),
                frame: f[12].parse().unwrap(),
                context: None,
                next: input,
            }));
        }
        let actual =
            windows_from_translated_aligns(&input, &query_lengths, subject_nt_length, 200, false)
                .unwrap();
        let mut expected = Vec::new();
        let mut expected_translation = Vec::new();
        let mut pending_translation = false;
        for row in trace.lines() {
            if row.starts_with("K_CALL\t0\tredo_return\t") {
                break;
            }
            if row.starts_with("K_TRANSLATED\t") {
                if pending_translation {
                    expected_translation.push(row.split('\t').nth(2).unwrap().to_string());
                    pending_translation = false;
                }
                continue;
            }
            if !row.contains("\tpartial_translation\t") {
                continue;
            }
            let f: Vec<_> = row.split('\t').collect();
            assert_eq!(f[5], "0");
            let value = (
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[6].parse::<i32>().unwrap(),
            );
            if expected.last() != Some(&value) {
                expected.push(value);
                pending_translation = true;
            }
        }
        let observed: Vec<_> = actual
            .iter()
            .map(|window| {
                let aa_length = window.subject_range.end - window.subject_range.begin;
                (
                    3 * aa_length + window.subject_range.context.abs() - 1,
                    window.subject_range.context,
                    aa_length,
                )
            })
            .collect();
        assert_eq!(observed, expected, "{case}");
        assert_eq!(expected_translation.len(), actual.len(), "{case}");
        let subject_fasta = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_c/{}/subjects.fna",
            env!("CARGO_MANIFEST_DIR"),
            case.strip_suffix("_default").unwrap()
        ))
        .unwrap();
        let subject_nt: Vec<u8> = subject_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        assert_eq!(subject_nt.len(), subject_nt_length as usize, "{case}");
        let genetic_code = GeneticCode::from_id(1);
        // EXPERIMENT (LOSAT_X_TBNSTORE): the windows of the packed subject are the same
        let store =
            crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(&subject_nt)
                .expect("subject letters");
        for (window, expected_hex) in actual.iter().zip(expected_translation) {
            let translated = translated_subject_window_sequence(
                SubjectNt::Letters(&subject_nt),
                &window.subject_range,
                &genetic_code,
            )
            .unwrap();
            assert_eq!(
                translated_subject_window_sequence(
                    SubjectNt::Packed(&store),
                    &window.subject_range,
                    &genetic_code,
                )
                .unwrap(),
                translated,
                "{case}: packed {:?}",
                window.subject_range
            );
            let observed_hex = translated[1..translated.len() - 1]
                .iter()
                .map(|value| format!("{value:02x}"))
                .collect::<String>();
            assert_eq!(
                observed_hex, expected_hex,
                "{case}: {:?}",
                window.subject_range
            );
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:906-938:
        let first_window = &actual[0];
        let first_align = first_window.align.as_deref().unwrap();
        let translated = translated_subject_window_sequence(
            SubjectNt::Letters(&subject_nt),
            &first_window.subject_range,
            &genetic_code,
        )
        .unwrap();
        let observed =
            translated_subject_composition(&translated, &first_window.subject_range, first_align)
                .unwrap();
        let composition_trace = std::fs::read_to_string(format!(
            "{}/../docs/evidence/tlosan_stage_d/kappa_composition_matrix_scores_20260924/{case}.tsv",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let expected_row = composition_trace
            .lines()
            .find(|line| line.starts_with("K_COMP\t0\tsubject\t"))
            .unwrap();
        let fields: Vec<_> = expected_row.split('\t').collect();
        assert_eq!(
            observed.num_true_amino_acids,
            fields[3].parse().unwrap(),
            "{case}"
        );
        for (index, (value, bits)) in observed.prob.iter().zip(&fields[4..]).enumerate() {
            assert_eq!(
                value.to_bits(),
                u64::from_str_radix(bits, 16).unwrap(),
                "{case}: composition index {index}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1417-1455,1530-1548:
#[test]
fn translated_kappa_seg_bytes_match_ncbi_mask_calls() {
    let trace = std::fs::read_to_string(format!(
        "{}/../docs/evidence/tlosan_stage_d/kappa_seg_20260924/multi_query_20260924_default.tsv",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap();
    let mut raw = std::collections::HashMap::new();
    let decode = |hex: &str| -> Vec<u8> {
        hex.as_bytes()
            .chunks_exact(2)
            .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
            .collect()
    };
    for row in trace.lines() {
        let fields: Vec<_> = row.split('\t').collect();
        if fields[0] == "K_SEG_RAW" {
            let residues = decode(fields[5]);
            assert_eq!(residues.len(), fields[4].parse::<usize>().unwrap());
            raw.insert(fields[1].parse::<usize>().unwrap(), residues);
        }
    }
    let mut compared = 0;
    for row in trace.lines() {
        let fields: Vec<_> = row.split('\t').collect();
        if fields[0] != "K_SEG_MASKED" {
            continue;
        }
        let event = fields[1].parse::<usize>().unwrap();
        let mut translated = Vec::with_capacity(raw[&event].len() + 2);
        translated.push(0);
        translated.extend_from_slice(&raw[&event]);
        translated.push(0);
        assert_eq!(translated.len() - 2, fields[2].parse::<usize>().unwrap());
        assert!(mask_translated_subject_seg(&mut translated).unwrap());
        assert_eq!(
            &translated[1..translated.len() - 1],
            decode(fields[3]),
            "SEG event {event}"
        );
        compared += 1;
    }
    assert_eq!(compared, 7);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1172-1208;
// core/blast_kappa.c:1488-1552,1684-1704
#[test]
fn translated_get_range_order_and_seg_match_ncbi_saved_calls() {
    fn decode(hex: &str) -> Vec<u8> {
        hex.as_bytes()
            .chunks_exact(2)
            .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
            .collect()
    }
    for case in ["seg_hard_query_20260924", "multi_query_20260924"] {
        let root = format!("{}/../docs/evidence", env!("CARGO_MANIFEST_DIR"));
        let mode_trace = std::fs::read_to_string(format!(
            "{root}/tlosan_stage_d/kappa_mode2_20260924/{case}_default.tsv"
        ))
        .unwrap();
        let seg_trace = std::fs::read_to_string(format!(
            "{root}/tlosan_stage_d/kappa_seg_20260924/{case}_default.tsv"
        ))
        .unwrap();
        let query_fasta =
            std::fs::read_to_string(format!("{root}/tlosan_stage_c/{case}/query.faa")).unwrap();
        let subject_fasta =
            std::fs::read_to_string(format!("{root}/tlosan_stage_c/{case}/subjects.fna")).unwrap();
        let mut queries: Vec<Vec<u8>> = Vec::new();
        for line in query_fasta.lines() {
            if line.starts_with('>') {
                queries.push(Vec::new());
            } else {
                queries
                    .last_mut()
                    .unwrap()
                    .extend(line.bytes().map(aa_char_to_ncbistdaa));
            }
        }
        let query_lengths: Vec<_> = queries.iter().map(|q| q.len() as i32).collect();
        let subject_nt: Vec<_> = subject_fasta
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.bytes())
            .collect();
        let matching_seq =
            BlastCompoMatchingSequence::new_translated(0, SubjectNt::Letters(&subject_nt), 1)
                .unwrap();
        let mut raw_events = Vec::new();
        let mut masked_events = std::collections::HashMap::new();
        for line in seg_trace.lines() {
            let fields: Vec<_> = line.split('\t').collect();
            match fields[0] {
                "K_SEG_RAW" => raw_events.push((
                    fields[1].parse::<usize>().unwrap(),
                    fields[2].parse::<i32>().unwrap(),
                    fields[3].parse::<i32>().unwrap(),
                    decode(fields[5]),
                )),
                "K_SEG_MASKED" => {
                    masked_events.insert(fields[1].parse::<usize>().unwrap(), decode(fields[3]));
                }
                _ => {}
            }
        }
        let mut observed_calls = 0;
        let mut observed_masks = 0;
        for enter in mode_trace
            .lines()
            .filter(|line| line.contains("\tredo_enter\t"))
        {
            let fields: Vec<_> = enter.split('\t').collect();
            let redo_id = fields[1];
            let lambda = fields[4].parse::<f64>().unwrap();
            let nt_length = fields[7].parse::<i32>().unwrap();
            assert_eq!(nt_length, matching_seq.length);
            let input_rows: Vec<_> = mode_trace
                .lines()
                .filter(|line| line.starts_with(&format!("K_ALIGN\t{redo_id}\tincoming\t")))
                .collect();
            assert_eq!(input_rows.len(), fields[3].parse::<usize>().unwrap());
            let mut input: Option<Box<BlastCompoAlignment>> = None;
            for row in input_rows.into_iter().rev() {
                let f: Vec<_> = row.split('\t').collect();
                input = Some(Box::new(BlastCompoAlignment {
                    score: f[5].parse().unwrap(),
                    matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
                    query_index: f[7].parse().unwrap(),
                    query_start: f[8].parse().unwrap(),
                    query_end: f[9].parse().unwrap(),
                    match_start: f[10].parse().unwrap(),
                    match_end: f[11].parse().unwrap(),
                    frame: f[12].parse().unwrap(),
                    context: None,
                    next: input,
                }));
            }
            let windows = windows_from_translated_aligns(
                &input,
                &query_lengths,
                nt_length,
                WINDOW_BORDER,
                false,
            )
            .unwrap();
            let params = BlastRedoAlignParams {
                matrix_info: crate::core::composition_adjustment::adjust_scores::build_matrix_info(
                    crate::config::ScoringMatrix::Blosum62,
                    0.3176,
                )
                .unwrap(),
                gapping_params: BlastCompoGappingParams {
                    gap_open: 11,
                    gap_extend: 1,
                    x_dropoff: 10,
                    context: Cell::new(None),
                },
                compo_adjust_mode: BlastCompoAdjustMode::CompositionMatrixAdjust,
                alphsize: BLASTAA_SIZE as i32,
                composition_test_index: 0,
                log_k: 0.0,
                score_divisor: 32.0,
                smith_waterman: false,
                near_identical_cutoff: 1.74 * LOCAL_LN2 / lambda,
                position_based: false,
                ccat_query_length: *query_lengths.iter().max().unwrap(),
                query_is_translated: false,
                subject_is_translated: true,
                cutoff_score: 1,
                do_link_hsps: false,
            };
            for window in &windows {
                let query_index = usize::try_from(window.query_range.context).unwrap();
                let query_data = BlastCompoSequenceData::from_ncbistdaa(&queries[query_index]);
                let query_info = BlastCompoQueryInfo {
                    origin: 0,
                    seq: query_data.clone(),
                    composition: BlastAminoAcidComposition::empty(),
                    eff_search_space: 1.0,
                    words: Some(build_query_word_hashes(query_data.data())),
                };
                let mut old_near_identical = false;
                let mut maybe_biased = true;
                let mut hsp_index = 0;
                let mut align = window.align.as_deref();
                while let Some(current) = align {
                    let near_identical = if hsp_index == 0 || maybe_biased {
                        preliminary_test_near_identical(
                            &query_info,
                            window,
                            current,
                            params.near_identical_cutoff,
                        )
                    } else {
                        old_near_identical
                    };
                    if hsp_index == 0 || (maybe_biased && near_identical != old_near_identical) {
                        let event = &raw_events[observed_calls];
                        assert_eq!(event.0, observed_calls, "{case}: event index");
                        let raw = translated_subject_window_sequence(
                            SubjectNt::Letters(&subject_nt),
                            &window.subject_range,
                            &GeneticCode::from_id(1),
                        )
                        .unwrap();
                        assert_eq!(
                            event.1,
                            3 * (window.subject_range.end - window.subject_range.begin)
                                + window.subject_range.context.abs()
                                - 1
                        );
                        assert_eq!(event.2, window.subject_range.context);
                        assert_eq!(
                            &raw[1..raw.len() - 1],
                            event.3,
                            "{case}: raw event {}",
                            event.0
                        );
                        let result = translated_subject_get_range(
                            &matching_seq,
                            &window.subject_range,
                            &query_data,
                            &window.query_range,
                            query_info.words.as_deref(),
                            window.align.as_deref().unwrap(),
                            near_identical,
                            maybe_biased,
                            &params,
                        )
                        .unwrap();
                        assert_eq!(
                            result.query.data(),
                            query_data.data(),
                            "{case}: query {query_index}"
                        );
                        let expected = masked_events.get(&event.0).unwrap_or(&event.3);
                        assert_eq!(
                            result.subject.data(),
                            expected,
                            "{case}: subject event {}",
                            event.0
                        );
                        if masked_events.contains_key(&event.0) {
                            assert!(result.subject_maybe_biased);
                            observed_masks += 1;
                        }
                        maybe_biased = result.subject_maybe_biased;
                        observed_calls += 1;
                    }
                    old_near_identical = near_identical;
                    hsp_index += 1;
                    align = current.next.as_deref();
                }
            }
        }
        assert_eq!(
            observed_calls,
            raw_events.len(),
            "{case}: get-range call order"
        );
        assert_eq!(
            observed_masks,
            masked_events.len(),
            "{case}: SEG call count"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1101-1146,1172-1269;
// core/blast_kappa.c:1475-1552,1896-1957
#[test]
fn translated_redo_core_uses_ncbi_hard_seg_mode2_inputs() {
    fn trace_alignment(row: &str) -> Box<BlastCompoAlignment> {
        let f: Vec<_> = row.split('\t').collect();
        blast_compo_alignment_new(
            f[5].parse().unwrap(),
            EMatrixAdjustRule::DontAdjustMatrix,
            f[8].parse().unwrap(),
            f[9].parse().unwrap(),
            f[7].parse().unwrap(),
            f[10].parse().unwrap(),
            f[11].parse().unwrap(),
            f[12].parse().unwrap(),
            None,
        )
    }
    fn assert_mode2_redo(
        incoming_align: &BlastCompoAlignment,
        rule: EMatrixAdjustRule,
        adjusted: Option<&AdjustedProteinMatrix>,
        query: &BlastCompoSequenceData,
        query_range: &BlastCompoSequenceRange,
        _ccat_query_length: i32,
        subject: &BlastCompoSequenceData,
        subject_range: &BlastCompoSequenceRange,
        _full_subject_length: i32,
        _params: &BlastRedoAlignParams,
    ) -> Result<Option<Box<BlastCompoAlignment>>> {
        assert_eq!(rule, EMatrixAdjustRule::UserSpecifiedRelEntropy);
        let matrix = adjusted.unwrap();
        assert_eq!(
            (
                matrix.scores[1][1],
                matrix.scores[1][16],
                matrix.scores[22][22]
            ),
            (149, -28, 202)
        );
        assert_eq!(query.length, 160);
        assert_eq!(
            query_range,
            &BlastCompoSequenceRange {
                begin: 0,
                end: 160,
                context: 0
            }
        );
        assert_eq!(&query.data()[..40], &[21u8; 40]);
        assert_eq!(subject.length, 120);
        assert_eq!(subject_range.context, 1);
        assert_eq!(subject.buffer[0], 0);
        assert_eq!(subject.buffer[subject.buffer.len() - 1], 0);
        Ok(Some(alignment_copy(incoming_align)))
    }
    let root = format!("{}/../docs/evidence", env!("CARGO_MANIFEST_DIR"));
    let mode_trace = std::fs::read_to_string(format!(
        "{root}/tlosan_stage_d/kappa_mode2_20260924/seg_hard_query_20260924_default.tsv"
    ))
    .unwrap();
    let incoming = Some(trace_alignment(
        mode_trace
            .lines()
            .find(|line| line.starts_with("K_ALIGN\t0\tincoming\t"))
            .unwrap(),
    ));
    let subject_fasta = std::fs::read_to_string(format!(
        "{root}/tlosan_stage_c/seg_hard_query_20260924/subjects.fna"
    ))
    .unwrap();
    let subject_nt: Vec<_> = subject_fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.bytes())
        .collect();
    let matching_seq =
        BlastCompoMatchingSequence::new_translated(0, SubjectNt::Letters(&subject_nt), 1).unwrap();
    let query_fasta = std::fs::read_to_string(format!(
        "{root}/tlosan_stage_c/seg_hard_query_20260924/query.faa"
    ))
    .unwrap();
    let mut query: Vec<_> = query_fasta
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.bytes().map(aa_char_to_ncbistdaa))
        .collect();
    assert_eq!(query.len(), 160);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:614-638 and saved hard-SEG query composition:
    // hard query filtering replaces the first 40 K residues with X.
    query[..40].fill(21);
    let composition = read_aa_composition(&query);
    let comp_trace = std::fs::read_to_string(format!(
        "{root}/tlosan_stage_d/kappa_composition_matrix_scores_20260924/seg_hard_query_20260924_default.tsv"
    )).unwrap();
    let expected_query = comp_trace
        .lines()
        .find(|line| line.starts_with("K_COMP\t0\tquery\t"))
        .unwrap();
    let expected: Vec<_> = expected_query.split('\t').collect();
    assert_eq!(
        composition.num_true_amino_acids,
        expected[3].parse().unwrap()
    );
    for (value, bits) in composition.prob.iter().zip(&expected[4..]) {
        assert_eq!(value.to_bits(), u64::from_str_radix(bits, 16).unwrap());
    }
    let query_info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(&query),
        composition,
        eff_search_space: 1.0,
        words: Some(build_query_word_hashes(&query)),
    };
    let params = BlastRedoAlignParams {
        matrix_info: crate::core::composition_adjustment::adjust_scores::build_matrix_info(
            crate::config::ScoringMatrix::Blosum62,
            0.0099251861761165822,
        )
        .unwrap(),
        gapping_params: BlastCompoGappingParams {
            gap_open: 11,
            gap_extend: 1,
            x_dropoff: 10,
            context: Cell::new(None),
        },
        compo_adjust_mode: BlastCompoAdjustMode::CompositionMatrixAdjust,
        alphsize: BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: 0.0,
        score_divisor: 32.0,
        smith_waterman: false,
        near_identical_cutoff: 1.74 * LOCAL_LN2 / 0.0083437500000000005,
        position_based: false,
        ccat_query_length: 160,
        query_is_translated: false,
        subject_is_translated: true,
        cutoff_score: 320,
        do_link_hsps: false,
    };
    let callbacks = BlastRedoAlignCallbacks {
        calc_lambda: None,
        get_range: translated_subject_get_range,
        redo_one_alignment: assert_mode2_redo,
    };
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    let result = blast_redo_one_match_with_workspace_queries_and_matrix(
        &incoming,
        &params,
        &matching_seq,
        &[query_info],
        0.0083437500000000005,
        &callbacks,
        &mut workspace,
        &mut None,
    )
    .unwrap();
    assert_eq!(result.lambda_ratio.unwrap().to_bits(), 1.0f64.to_bits());
    assert_eq!(result.alignments_by_query.len(), 1);
    assert_eq!(
        result.alignments_by_query[0]
            .as_ref()
            .unwrap()
            .matrix_adjust_rule,
        EMatrixAdjustRule::DontAdjustMatrix
    );
}

fn make_align(
    score: i32,
    q_start: i32,
    q_end: i32,
    s_start: i32,
    s_end: i32,
    context: usize,
) -> Box<BlastCompoAlignment> {
    blast_compo_alignment_new(
        score,
        EMatrixAdjustRule::DontAdjustMatrix,
        q_start,
        q_end,
        0,
        s_start,
        s_end,
        0,
        Some(PreliminaryHspIndex(context)),
    )
}

#[test]
fn test_distinct_alignments_sort_orders_ncbi_style() {
    let mut head = Some(make_align(50, 20, 40, 10, 30, 0));
    if let Some(node) = head.as_mut() {
        node.next = Some(make_align(60, 30, 50, 10, 31, 1));
        if let Some(node2) = node.next.as_mut() {
            node2.next = Some(make_align(60, 25, 45, 9, 31, 2));
        }
    }

    distinct_alignments_sort(&mut head);
    let values = list_to_vec(head);
    assert_eq!(values[0].context, Some(PreliminaryHspIndex(2)));
    assert_eq!(values[1].context, Some(PreliminaryHspIndex(1)));
    assert_eq!(values[2].context, Some(PreliminaryHspIndex(0)));
}

#[test]
fn test_with_distinct_ends_keeps_best_shared_endpoint() {
    let mut old = Some(make_align(80, 10, 30, 100, 120, 0));
    let mut new_align = Some(make_align(70, 10, 32, 100, 122, 1));
    with_distinct_ends(&mut new_align, &mut old, true);
    let values = list_to_vec(old);
    assert_eq!(values.len(), 1);
    assert_eq!(values[0].context, Some(PreliminaryHspIndex(0)));
}

#[test]
fn test_with_distinct_ends_similar_endpoint_requires_same_frame() {
    let mut old = Some(make_align(80, 10, 30, 100, 120, 0));
    old.as_mut().unwrap().frame = 1;

    let mut new_align = Some(make_align(70, 12, 32, 102, 122, 1));
    new_align.as_mut().unwrap().frame = -1;

    with_distinct_ends(&mut new_align, &mut old, false);

    let values = list_to_vec(old);
    assert_eq!(values.len(), 2);
    assert_eq!(values[0].context, Some(PreliminaryHspIndex(1)));
    assert_eq!(values[1].context, Some(PreliminaryHspIndex(0)));
}

#[test]
fn test_is_contained_requires_score_tolerance() {
    let mut existing = Some(make_align(100, 10, 50, 100, 140, 0));
    let mut inner = make_align(90, 20, 40, 110, 130, 1);
    inner.next = None;
    assert!(is_contained(&inner, &existing, 0.267));
    existing.as_mut().unwrap().score = 80;
    assert!(!is_contained(&inner, &existing, 0.267));
}

#[test]
fn test_windows_from_protein_aligns_preserves_incoming_order_after_reverse_only() {
    let mut aligns = Some(make_align(10, 0, 10, 0, 10, 0));
    if let Some(node) = aligns.as_mut() {
        node.next = Some(make_align(20, 5, 15, 20, 30, 1));
    }
    let query_info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(&[1, 2, 3, 4]),
        composition: BlastAminoAcidComposition::empty(),
        eff_search_space: 1.0,
        words: None,
    };

    let windows = windows_from_protein_aligns(&aligns, &query_info, 100).unwrap();
    assert_eq!(windows.len(), 1);
    assert_eq!(windows[0].subject_range.begin, 0);
    assert_eq!(windows[0].subject_range.end, 100);
    let mut contexts = Vec::new();
    let mut current = windows[0].align.as_deref();
    while let Some(node) = current {
        contexts.push(node.context.clone());
        current = node.next.as_deref();
    }
    assert_eq!(
        contexts,
        vec![Some(PreliminaryHspIndex(0)), Some(PreliminaryHspIndex(1))]
    );
}

#[test]
fn test_windows_from_protein_aligns_separates_queries_by_index() {
    let mut aligns = Some(make_align(10, 0, 10, 0, 10, 0));
    if let Some(node) = aligns.as_mut() {
        node.next = Some(make_align(20, 5, 15, 20, 30, 1));
        node.next.as_mut().unwrap().query_index = 1;
        node.next.as_mut().unwrap().next = Some(make_align(24, 4, 16, 40, 50, 2));
        node.next
            .as_mut()
            .unwrap()
            .next
            .as_mut()
            .unwrap()
            .query_index = 1;
    }

    let query_info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(&[1, 2, 3, 4]),
        composition: BlastAminoAcidComposition::empty(),
        eff_search_space: 1.0,
        words: None,
    };

    let windows = windows_from_protein_aligns(&aligns, &query_info, 100).unwrap();
    assert_eq!(windows.len(), 2);
    assert_eq!(windows[0].query_range.context, 0);
    assert_eq!(windows[1].query_range.context, 1);
    assert_eq!(windows[0].hspcnt, 1);
    assert_eq!(windows[1].hspcnt, 2);

    let mut query1_contexts = Vec::new();
    let mut current = windows[1].align.as_deref();
    while let Some(node) = current {
        query1_contexts.push(node.context.clone());
        current = node.next.as_deref();
    }
    assert_eq!(
        query1_contexts,
        vec![Some(PreliminaryHspIndex(1)), Some(PreliminaryHspIndex(2))]
    );
}

#[test]
fn test_blast_redo_one_match_keeps_distinct_alignments_per_query_index() {
    let mut aligns = Some(make_align(80, 0, 40, 0, 40, 0));
    if let Some(node) = aligns.as_mut() {
        node.query_index = 2;
        node.next = Some(make_align(60, 10, 30, 50, 70, 1));
        node.next.as_mut().unwrap().query_index = 2;
    }

    let query_info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(&[1u8; 100]),
        composition: BlastAminoAcidComposition::empty(),
        eff_search_space: 1.0,
        words: None,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2939-2955
    let matching_residues = vec![1u8; 100];
    let matching_seq = BlastCompoMatchingSequence::new(0, &matching_residues);
    let params = BlastRedoAlignParams {
        matrix_info: crate::core::composition_adjustment::adjust_scores::build_matrix_info(
            crate::config::ScoringMatrix::Blosum62,
            0.3176,
        )
        .unwrap(),
        gapping_params: BlastCompoGappingParams {
            gap_open: 11,
            gap_extend: 1,
            x_dropoff: 10,
            context: Cell::new(None),
        },
        compo_adjust_mode: BlastCompoAdjustMode::NoCompositionBasedStats,
        alphsize: BLASTAA_SIZE as i32,
        composition_test_index: 0,
        log_k: 0.0,
        score_divisor: 1.0,
        smith_waterman: false,
        near_identical_cutoff: 0.0,
        position_based: false,
        ccat_query_length: 100,
        query_is_translated: false,
        subject_is_translated: false,
        cutoff_score: 1,
        do_link_hsps: false,
    };
    fn test_get_range<'seq>(
        _matching_seq: &BlastCompoMatchingSequence<'seq>,
        subject_range: &BlastCompoSequenceRange,
        orig_query: &BlastCompoSequenceData,
        query_range: &BlastCompoSequenceRange,
        _query_words: Option<&[u64]>,
        _align: &BlastCompoAlignment,
        _should_test_identical: bool,
        _subject_maybe_biased: bool,
        _params: &BlastRedoAlignParams,
    ) -> Result<BlastRedoRangeResult> {
        Ok(BlastRedoRangeResult {
            query: orig_query.slice_range(query_range),
            subject: BlastCompoSequenceData::from_ncbistdaa(&vec![
                1u8;
                (subject_range.end - subject_range.begin)
                    as usize
            ]),
            subject_maybe_biased: false,
            x_subject_composition: None,
        })
    }

    fn test_redo_one_alignment(
        incoming_align: &BlastCompoAlignment,
        matrix_adjust_rule: EMatrixAdjustRule,
        _adjusted_matrix: Option<&AdjustedProteinMatrix>,
        _query_data: &BlastCompoSequenceData,
        _query_range: &BlastCompoSequenceRange,
        _ccat_query_length: i32,
        _subject_data: &BlastCompoSequenceData,
        _subject_range: &BlastCompoSequenceRange,
        _full_subject_length: i32,
        _params: &BlastRedoAlignParams,
    ) -> Result<Option<Box<BlastCompoAlignment>>> {
        Ok(Some(blast_compo_alignment_new(
            incoming_align.score,
            matrix_adjust_rule,
            incoming_align.query_start,
            incoming_align.query_end,
            incoming_align.query_index,
            incoming_align.match_start,
            incoming_align.match_end,
            incoming_align.frame,
            incoming_align.context.clone(),
        )))
    }

    let callbacks = BlastRedoAlignCallbacks {
        calc_lambda: None,
        get_range: test_get_range,
        redo_one_alignment: test_redo_one_alignment,
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
    let mut workspace = BlastCompositionWorkspace::new_blosum62();
    let result = blast_redo_one_match_with_workspace_queries_and_matrix(
        &aligns,
        &params,
        &matching_seq,
        std::slice::from_ref(&query_info),
        0.267,
        &callbacks,
        &mut workspace,
        &mut None,
    )
    .unwrap();

    assert_eq!(result.alignments_by_query.len(), 3);
    assert!(result.alignments_by_query[0].is_none());
    assert!(result.alignments_by_query[1].is_none());
    let mut contexts = Vec::new();
    let mut current = result.alignments_by_query[2].as_deref();
    while let Some(node) = current {
        contexts.push(node.context.clone());
        current = node.next.as_deref();
    }
    assert_eq!(
        contexts,
        vec![Some(PreliminaryHspIndex(1)), Some(PreliminaryHspIndex(0))]
    );
}

#[test]
fn test_preliminary_test_near_identical_uses_cutoff_branch() {
    let query_info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(&[1u8; 60]),
        composition: BlastAminoAcidComposition::empty(),
        eff_search_space: 1.0,
        words: None,
    };
    let window = window_info_new(0, 60, 0, 0, 60, 0, Some(make_align(240, 0, 59, 0, 59, 0)));
    let align = window.align.as_deref().unwrap();
    assert!(preliminary_test_near_identical(
        &query_info,
        &window,
        align,
        3.5
    ));
    assert!(!preliminary_test_near_identical(
        &query_info,
        &window,
        align,
        5.0
    ));
}

#[test]
fn test_location_compare_windows_sorts_by_query_then_subject() {
    let left = window_info_new(0, 100, 0, 0, 50, 0, None);
    let right = window_info_new(10, 100, 0, 0, 50, 0, None);
    assert_eq!(location_compare_windows(&left, &right), Ordering::Less);
}

#[test]
fn test_build_query_word_hashes_matches_ncbi_shift_hash() {
    let query = [1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let hashes = build_query_word_hashes(&query);
    assert_eq!(hashes.len(), query.len() - NEAR_IDENTICAL_WORD_SIZE + 1);
    assert_eq!(hashes[0], get_hash(&query, NEAR_IDENTICAL_WORD_SIZE));
    assert_eq!(hashes[1], get_hash(&query[1..], NEAR_IDENTICAL_WORD_SIZE));
    assert_eq!(hashes[2], 0);
}

#[test]
fn test_find_num_identical_accepts_short_hash_slice_for_9aa_window() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1160-1188
    //
    // A 9-aa slice has only two stored 8-mer hashes; NCBI still scans it.
    let query = [1u8, 3, 4, 5, 6, 7, 8, 9, 10];
    let hashes = build_query_word_hashes(&query);
    assert_eq!(hashes.len(), 2);
    assert_eq!(
        find_num_identical(&query[..], &hashes[..], &query, NEAR_IDENTICAL_MAX_SHIFT),
        query.len()
    );
}

#[test]
fn test_test_near_identical_accepts_high_identity_alignment() {
    let residues = vec![1u8; 80];
    let query = BlastCompoSequenceData::from_ncbistdaa(&residues);
    let subject = BlastCompoSequenceData::from_ncbistdaa(&residues);
    let align = make_align(80, 0, 80, 0, 80, 0);
    let words = build_query_word_hashes(query.data());
    assert!(test_near_identical(
        &subject,
        0,
        &query,
        0,
        Some(&words),
        &align,
    ));
}

#[test]
fn test_test_near_identical_rejects_low_identity_alignment() {
    let query = BlastCompoSequenceData::from_ncbistdaa(&[1u8; 80]);
    let subject = BlastCompoSequenceData::from_ncbistdaa(&[2u8; 80]);
    let align = make_align(20, 0, 80, 0, 80, 0);
    let words = build_query_word_hashes(query.data());
    assert!(!test_near_identical(
        &subject,
        0,
        &query,
        0,
        Some(&words),
        &align,
    ));
}

#[test]
fn test_test_near_identical_keeps_fraction_check_when_middle_slice_has_no_words() {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1314-1333
    //
    // A short unresolved middle slice still flows to the final fraction
    // check; NCBI does not return FALSE just because there are no 8-mer
    // hashes left to inspect.
    let mut query_residues = vec![1u8; 48];
    query_residues.extend_from_slice(&[3, 4, 5, 6]);
    query_residues.extend_from_slice(&[1u8; 48]);

    let mut subject_residues = vec![1u8; 48];
    subject_residues.extend_from_slice(&[7, 8, 9, 10]);
    subject_residues.extend_from_slice(&[1u8; 48]);

    let query = BlastCompoSequenceData::from_ncbistdaa(&query_residues);
    let subject = BlastCompoSequenceData::from_ncbistdaa(&subject_residues);
    let align = make_align(96, 0, 100, 0, 100, 0);

    assert!(test_near_identical(&subject, 0, &query, 0, None, &align));
}
// LOSAT_X_BXENDMEM: the test on the residues gives the table's answer: letters from small
// alphabets (so that words repeat), ranges not starting at 0 (the table is read at the range's
// offsets, as NCBI reads it), and subjects with runs of letter 0 (the table's last entry).
#[test]
fn near_identical_on_residues_equals_the_table() {
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move |n: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % n
    };
    let mut true_answers = 0;
    for round in 0..20_000 {
        let letters = [2, 3, 28][round % 3];
        let context: Vec<u8> = (0..4 + next(120)).map(|_| next(letters) as u8).collect();
        let begin = next(1 + context.len() as u64 / 2) as usize;
        let end = begin + next(1 + (context.len() - begin) as u64) as usize;
        if end == begin {
            continue;
        }
        let orig = BlastCompoSequenceData::from_ncbistdaa(&context);
        let query = BlastCompoSequenceData::from_ncbistdaa(&context[begin..end]);
        let mut subject: Vec<u8> = (0..4 + next(120))
            .map(|_| {
                if next(2) == 0 {
                    context[next(context.len() as u64) as usize]
                } else {
                    next(letters) as u8
                }
            })
            .collect();
        if next(4) == 0 && subject.len() > 20 {
            subject[2..14].fill(0);
        }
        let subject = BlastCompoSequenceData::from_ncbistdaa(&subject);
        let query_start = begin as i32 + next(1 + (end - begin) as u64 / 4) as i32;
        let query_end = end as i32 - next(1 + (end - begin) as u64 / 4) as i32;
        let match_start = next(1 + subject.data().len() as u64 / 4) as i32;
        let match_end =
            subject.data().len() as i32 - next(1 + subject.data().len() as u64 / 4) as i32;
        if query_start >= query_end || match_start >= match_end {
            continue;
        }
        let align = blast_compo_alignment_new(
            1,
            EMatrixAdjustRule::DontAdjustMatrix,
            query_start,
            query_end,
            0,
            match_start,
            match_end,
            1,
            None,
        );
        let table = build_query_word_hashes(orig.data());
        let want = test_near_identical(&subject, 0, &query, begin as i32, Some(&table), &align);
        let got = x_test_near_identical_on_residues(
            &subject,
            0,
            &query,
            begin as i32,
            orig.data(),
            &align,
        );
        assert_eq!(got, want, "round {round}");
        true_answers += usize::from(want);
    }
    assert!(true_answers > 100, "true answers {true_answers}");
}

/// The ncbistdaa letters of `text`.
fn ncbistdaa(text: &str) -> Vec<u8> {
    text.bytes()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect()
}

/// `test_near_identical_on_copy` of the profile letters `letters` against `subject` over the
/// alignment `align`: its answer and whether it read past the copy or the words.
fn walk_on_copy(letters: &[u8], subject: &[u8], align: &BlastCompoAlignment) -> (bool, bool) {
    let query = BlastCompoSequenceData::from_ncbistdaa(letters);
    let words = build_query_word_hashes(letters);
    let read_past = Cell::new(false);
    let identical = test_near_identical_on_copy(
        &BlastCompoSequenceData::from_ncbistdaa(subject),
        0,
        NearIdenticalCopy::new(&query, &words, &read_past),
        0,
        align,
    );
    (identical, read_past.get())
}

// blast_kappa.c:1259-1340 and 944-1100 (cited at `test_near_identical_on_copy`), PS10e's probes
// T_AA1_q1, Q_P10_q1 and Q_P6_q1: c1's 22 columns aligned with the user query of c1's letters,
// and the profile's copy of its `.psq` letters: 6 letters 00 00 00 00 00 02 (c1 as 22 A),
// c1's first 10 (IVLARIDDRL) and first 6 (IVLARI, no word array). NCBI's walk reads past its copy
// of L_pin + 2 bytes in each (valgrind's invalid reads; SIGSEGV for the null words of P6), as the
// copy view reports; c1's own 22 letters read nothing past it and answer as the slices' walk.
#[test]
fn the_walk_on_a_profile_copy_reports_reads_past_it() {
    let c1 = ncbistdaa("IVLARIDDRLIHGQVATTWTKE");
    let align = make_align(100, 0, 22, 0, 22, 0);
    let (identical, read_past) = walk_on_copy(&c1, &c1, &align);
    assert!(!read_past);
    let query = BlastCompoSequenceData::from_ncbistdaa(&c1);
    let words = build_query_word_hashes(&c1);
    let on_slices = test_near_identical(
        &BlastCompoSequenceData::from_ncbistdaa(&c1),
        0,
        &query,
        0,
        Some(&words),
        &align,
    );
    assert!(identical);
    assert_eq!(identical, on_slices);
    for letters in [vec![0, 0, 0, 0, 0, 2], c1[..10].to_vec(), c1[..6].to_vec()] {
        let (_, read_past) = walk_on_copy(&letters, &c1, &align);
        assert!(read_past, "{letters:?}");
    }
}

// blast_kappa.c:1161-1166 (cited at `CopyWords::missing`): a null word table (fewer than 8
// letters) stops `s_FindNumIdentical` only at offset 0; from another start NCBI reads the null
// table plus the start, which the copy view reports. A table of L - 7 words reads past it from
// word L - 7 on.
#[test]
fn the_copy_words_are_missing_only_at_their_start() {
    let read_past = Cell::new(false);
    let none = CopyWords {
        words: None,
        start: 0,
        read_past: &read_past,
    };
    assert!(none.missing(20));
    assert!(!read_past.get());
    let later = CopyWords { start: 3, ..none };
    assert!(!later.missing(20));
    assert_eq!(later.word(0), 0);
    assert!(read_past.get());
    let letters = ncbistdaa("IVLARIDDRL");
    let words = build_query_word_hashes(&letters);
    assert_eq!(words.len(), 3);
    let read_past = Cell::new(false);
    let table = CopyWords {
        words: Some(&words),
        start: 1,
        read_past: &read_past,
    };
    assert_eq!(table.word(1), words[2]);
    assert!(!read_past.get());
    table.word(2);
    assert!(read_past.get());
}

// blast_kappa.c:1671-1708 (cited at `NearIdenticalCopy`): the copy holds the letters and one zero
// byte after them (index L_pin); index L_pin + 1 on is past NCBI's block.
#[test]
fn the_copy_letters_end_after_their_zero_byte() {
    let letters = ncbistdaa("IVLARI");
    let query = BlastCompoSequenceData::from_ncbistdaa(&letters);
    let read_past = Cell::new(false);
    let copy = NearIdenticalCopy::new(&query, &[], &read_past).letters();
    let window = copy.window(2, 9).unwrap();
    assert_eq!(window.len(), 8);
    assert_eq!(window.at(0), letters[2]);
    assert_eq!(window.at(4), 0);
    assert!(!read_past.get());
    assert_eq!(window.at(5), 0);
    assert!(read_past.get());
}
