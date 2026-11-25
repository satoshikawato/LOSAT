use super::*;
use crate::algorithm::tblastn::search_gapped::{
    preliminary_stages, preliminary_stages_of, word_finder_profile,
};
use crate::sequence::packed_subject::SubjectNt;
use std::fs;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:259-310,478-487:
// subject->seq_ranges and subject->length are filled before WordFinder.
// The comparison-only probe records those exact input rows at the call.
fn ncbi_chunk_range_rows(path: &str) -> Vec<(i8, usize, Vec<(i32, i32)>)> {
    let trace = fs::read_to_string(path).unwrap();
    let mut calls: Vec<(i8, usize, Vec<(i32, i32)>)> = Vec::new();
    for line in trace.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        match fields[0] {
            "RANGE_CALL" => {
                assert_eq!(fields[1].parse::<usize>().unwrap(), calls.len());
                calls.push((
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                    Vec::new(),
                ));
            }
            "RANGE" => {
                let call = fields[1].parse::<usize>().unwrap();
                assert_eq!(call, calls.len() - 1);
                assert_eq!(fields[2].parse::<usize>().unwrap(), calls[call].2.len());
                calls[call]
                    .2
                    .push((fields[3].parse().unwrap(), fields[4].parse().unwrap()));
            }
            "SET_RANGES" | "SET_RANGE" => {}
            _ => panic!("unexpected NCBI range trace row"),
        }
    }
    calls
}

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let text = fs::read_to_string(path).unwrap();
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(id) = line.strip_prefix('>') {
            records.push((id.to_string(), Vec::new()));
        } else {
            records.last_mut().unwrap().1.extend(line.as_bytes());
        }
    }
    records
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-614:
// the selected matrix, wordsize, x-dropoff, and per-context cutoff
// determine two-hit extension before frame-local score sorting.
#[test]
fn alternate_matrix_word2_wordfinder_hsps_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let observed_cutoffs: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("PARAM\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[2].parse::<i32>().unwrap(),
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i32>().unwrap(),
            )
        })
        .collect();
    assert_eq!(observed_cutoffs, [(21, 21, 0); 6]);
    let trace = fs::read_to_string(format!("{root}wordfinder.stderr")).unwrap();
    let frame_order = [1, 2, 3, -1, -2, -3];
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            InitHsp {
                frame: frame_order[f[1].parse::<usize>().unwrap()],
                chunk_offset: 0,
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    let actual = preliminary_stages(
        &refs,
        subject,
        1,
        word_finder_profile(
            None,
            16,
            60,
            &[21, 21, 0],
            &[0, 0, i32::MAX],
            false,
            ScoringMatrix::Blosum45,
            2,
        ),
    )
    .unwrap()
    .initial;
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!("first BLOSUM45/word-2 WordFinder difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
        actual.get(first), expected.get(first), actual.len(), expected.len());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:562-583:
// BSearchContextInfo finds the query context before two-hit extension;
// the left word must remain within that context's query_offset.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234:
// `Blast_InitHitListSortByScore` fixes frame-local HSP order.
#[test]
fn multi_query_wordfinder_hsps_match_ncbi_in_order() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let trace = fs::read_to_string(format!("{root}wordfinder.stderr")).unwrap();
    let frame_order = [1, 2, 3, -1, -2, -3];
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            InitHsp {
                frame: frame_order[f[1].parse::<usize>().unwrap()],
                chunk_offset: 0,
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    let actual = preliminary_stages(
        &refs,
        subject,
        1,
        word_finder_profile(
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!(
            "first multi-query WordFinder difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
            actual.get(first),
            expected.get(first),
            actual.len(),
            expected.len()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
// The saved natural multi-query fixture has score 20 in context 0 and
// score 393 in context 1, so equality and just-above cutoffs are distinct.
#[test]
fn wordfinder_cutoff_is_selected_per_query_context_at_score_boundary() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let equal = preliminary_stages(
        &refs,
        subject,
        1,
        word_finder_profile(
            None,
            13,
            40,
            &[16, 16, 0],
            &[20, 393, i32::MAX],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    assert!(equal.iter().any(|h| h.score == 20 && h.q_seed < 121));
    assert!(equal
        .iter()
        .any(|h| h.score == 393 && (121..192).contains(&h.q_seed)));
    let above = preliminary_stages(
        &refs,
        subject,
        1,
        word_finder_profile(
            None,
            13,
            40,
            &[16, 16, 0],
            &[21, 394, i32::MAX],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    assert!(!above.iter().any(|h| h.score == 20 && h.q_seed < 121));
    assert!(!above
        .iter()
        .any(|h| h.score == 393 && (121..192).contains(&h.q_seed)));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:562-590:
// The comparison-only real call changes context 1 x-drop from 16 to 1.
#[test]
fn real_wordfinder_distinct_context_xdrop_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
    let trace =
        fs::read_to_string(format!("{root}word_xdrop_real_path_20260924/trace.tsv")).unwrap();
    let injections: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("WORD_XDROP_INJECT\t"))
        .collect();
    assert_eq!(injections.len(), 6);
    for (call, line) in injections.iter().enumerate() {
        assert_eq!(*line, format!("WORD_XDROP_INJECT\t{call}\t16\t1"));
    }
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            InitHsp {
                frame: [1, 2, 3, -1, -2, -3][call],
                chunk_offset: 0,
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    assert_eq!(expected.len(), 32);
    let actual = preliminary_stages(
        &refs,
        subject,
        1,
        word_finder_profile(
            None,
            13,
            40,
            &[16, 1, 0],
            &[0, 0, i32::MAX],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!("first distinct-context x-drop HSP difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
               actual.get(first), expected.get(first), actual.len(), expected.len());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:816-831
// Converts subject nucleotide scan ranges to protein coordinates via `CONV_NUCL2PROT_COORDINATES`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-614
#[test]
fn lowercase_mask_init_hsps_match_ncbi_scores_and_internal_offsets() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/lowercase_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subjects = read_fasta(&format!("{root}subjects.fna"));
    let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].to_string(),
                InitHsp {
                    frame: f[2].parse().unwrap(),
                    chunk_offset: 0,
                    q_seed: f[4].parse().unwrap(),
                    s_seed: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    s_start: f[7].parse().unwrap(),
                    length: f[8].parse().unwrap(),
                    score: f[9].parse().unwrap(),
                },
            )
        })
        .collect();
    let mut actual = Vec::new();
    for (id, subject) in &subjects {
        for hit in preliminary_stages(
            &[query],
            subject,
            1,
            word_finder_profile(None, 13, 40, &[16], &[0], true, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .initial
        {
            actual.push((id.clone(), hit));
        }
    }
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:575-614
// Extends two-hit seeds via `s_BlastAaExtendTwoHit`.
#[test]
fn ambiguity_spectrum_init_hsp_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/ambiguity_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            InitHsp {
                frame: f[2].parse().unwrap(),
                chunk_offset: 0,
                q_seed: f[4].parse().unwrap(),
                s_seed: f[5].parse().unwrap(),
                q_start: f[6].parse().unwrap(),
                s_start: f[7].parse().unwrap(),
                length: f[8].parse().unwrap(),
                score: f[9].parse().unwrap(),
            }
        })
        .collect();
    let actual = preliminary_stages(
        &[query],
        subject,
        1,
        word_finder_profile(None, 13, 40, &[16], &[0], false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .initial;
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234,575-614
// Populates the initial hit list via `s_BlastAaWordFinder_TwoHit`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:747-775
#[test]
fn code32_api_oracle_init_hsps_match_scores_internal_offsets_and_order() {
    let fixture = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_a/fixtures/"
    );
    let trace_root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/code32_20260923/"
    );
    let query = &read_fasta(&format!("{fixture}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{fixture}subject_code32.fna"))[0].1;
    let trace = fs::read_to_string(format!("{trace_root}ncbi_init_trace.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            InitHsp {
                frame: f[0].parse().unwrap(),
                chunk_offset: 0,
                q_seed: f[2].parse().unwrap(),
                s_seed: f[3].parse().unwrap(),
                q_start: f[4].parse().unwrap(),
                s_start: f[5].parse().unwrap(),
                length: f[6].parse().unwrap(),
                score: f[7].parse().unwrap(),
            }
        })
        .collect();
    let actual = preliminary_stages(
        &[query],
        subject,
        32,
        word_finder_profile(
            None,
            13,
            40,
            &[16],
            &[13],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234,575-614
// Extends two-hit seeds via `s_BlastAaExtendTwoHit` using the context ungapped x-dropoff.
#[test]
fn saved_init_hsps_match_ncbi_fixture_scores_and_internal_offsets() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subjects = read_fasta(&format!("{root}subjects.fna"));
    let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].to_string(),
                InitHsp {
                    frame: f[2].parse().unwrap(),
                    chunk_offset: 0,
                    q_seed: f[4].parse().unwrap(),
                    s_seed: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    s_start: f[7].parse().unwrap(),
                    length: f[8].parse().unwrap(),
                    score: f[9].parse().unwrap(),
                },
            )
        })
        .collect();
    let params = fs::read_to_string(format!("{root}ncbi_word_params.tsv")).unwrap();
    let observed: Vec<_> = params
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (f[2].parse::<i32>().unwrap(), f[3].parse::<i32>().unwrap())
        })
        .collect();
    assert!(observed.iter().all(|pair| *pair == (16, 0)));
    let mut actual = Vec::new();
    for (id, subject) in &subjects {
        for hit in preliminary_stages(
            &[query],
            subject,
            1,
            word_finder_profile(None, 13, 40, &[16], &[0], false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .initial
        {
            actual.push((id.clone(), hit));
        }
    }
    assert_eq!(actual, expected);
}
// No NCBI counterpart: the inputs and NCBI's initial HSPs of the long-subject tests below
// (docs/evidence/tlosan_stage_c/long_chunk_20260924): a 15,000,962-base subject whose frames
// have two chunks each.
fn long_chunk_case() -> (Vec<u8>, Vec<u8>, Vec<InitHsp>) {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let query = read_fasta(&format!("{root}long_chunk_20260924/query.faa"))[0]
        .1
        .clone();
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    assert_eq!(plus1.len(), 362);
    let mut subject = Vec::with_capacity(15_000_962);
    for _ in 0..4_999_950 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 15_000_962);
    let trace = fs::read_to_string(format!("{root}long_chunk_20260924/frame_chunks.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            InitHsp {
                frame: [1, 2, 3, -1, -2, -3][call / 2],
                chunk_offset: if call.is_multiple_of(2) { 0 } else { 4_999_900 },
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    (query, subject, expected)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-552:
// s_GetNextSubjectChunk precedes each WordFinder call and the hit list
// is reset before the next chunk. Subject offsets stay chunk-local.
#[test]
fn long_subject_chunk_wordfinder_hsps_match_ncbi() {
    let (query, subject, expected) = long_chunk_case();
    let query = &query;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // Retained long-fixture PARAM rows give cutoff 28.
    let actual = preliminary_stages(
        &[query],
        &subject,
        1,
        word_finder_profile(
            None,
            13,
            40,
            &[16],
            &[28],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-552,805-812 (above)
// The same subject packed, as every `-db` sequence is: with LOSAT_X_TBNUNITXLATE its plan has
// no residues (its first frame has more than 1 MiB), and each unit translates its chunk and
// 4,096 residues past it (the first chunk of a frame) or to the end sentinel (the second) when
// it runs; the initial HSPs are NCBI's either way.
#[test]
fn long_packed_subject_chunk_wordfinder_hsps_match_ncbi() {
    let (query, subject, expected) = long_chunk_case();
    let store = crate::sequence::packed_subject::PackedNucleotide::from_subject_letters(&subject)
        .expect("letters");
    let plan = super::super::search_seed::x_scan_plan(
        &[&query],
        SubjectNt::Packed(&store).into(),
        1,
        None,
        13,
        false,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &Default::default(),
        None,
        None,
    )
    .unwrap()
    .unwrap();
    assert_eq!(plan.subject.units.len(), 12);
    let lean = super::super::x_subject_translate::x_units_translate(subject.len());
    assert!(plan
        .subject
        .frames
        .iter()
        .all(|frame| frame.aa_seq.is_empty() == lean));
    drop(plan);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590 (above): cutoff 28.
    let actual = preliminary_stages_of(
        &[&query],
        SubjectNt::Packed(&store),
        1,
        word_finder_profile(
            None,
            13,
            40,
            &[16],
            &[28],
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .initial;
    assert_eq!(actual, expected);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:283-310:
// a soft range ending at the second chunk's start is included, then
// clipped to zero width before aa_ungapped.c:495-500 scans range zero.
#[test]
fn masked_chunk_boundary_wordfinder_hsp_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let query = &read_fasta(&format!("{root}masked_chunk_boundary_20260924/query.faa"))[0].1;
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    let mut subject = Vec::with_capacity(15_000_962);
    for _ in 0..4_999_950 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    subject[14_999_700..15_000_000].make_ascii_lowercase();
    // The (frame, chunk, scan ranges) units the search runs: the unit path's plan.
    let plan = super::super::search_seed::x_scan_plan(
        &[query],
        SubjectNt::Letters(&subject).into(),
        1,
        None,
        13,
        true,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &Default::default(),
        None,
        None,
    )
    .unwrap()
    .unwrap();
    let actual_ranges: Vec<_> = plan
        .subject
        .units
        .iter()
        .map(|(frame_index, chunk, ranges)| {
            (
                plan.subject.frames[*frame_index].frame,
                chunk.length,
                ranges.clone(),
            )
        })
        .collect();
    assert_eq!(
        actual_ranges,
        ncbi_chunk_range_rows(&format!(
            "{root}masked_chunk_boundary_20260924/chunk_ranges.tsv"
        ))
    );
    drop(plan);
    let trace = fs::read_to_string(format!(
        "{root}masked_chunk_boundary_20260924/frame_chunks.tsv"
    ))
    .unwrap();
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            InitHsp {
                frame: [1, 2, 3, -1, -2, -3][call / 2],
                chunk_offset: if call.is_multiple_of(2) { 0 } else { 4_999_900 },
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // Retained long-fixture PARAM rows give cutoff 28.
    let actual = preliminary_stages(
        &[query],
        &subject,
        1,
        word_finder_profile(None, 13, 40, &[16], &[28], true, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .initial;
    assert_eq!(actual, expected);
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:283-310,478-500:
// SUBJECT_SPLIT_NO_RANGE skips the middle positive-frame chunk; the
// later chunk still runs WordFinder with its chunk-local offsets.
#[test]
fn masked_no_range_middle_chunk_skips_and_later_hsp_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let query = &read_fasta(&format!("{root}no_range_middle_20260924/query.faa"))[0].1;
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    let mut subject = Vec::with_capacity(30_001_262);
    for _ in 0..10_000_050 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 30_001_262);
    subject[14_997_000..30_000_000].make_ascii_lowercase();
    let trace =
        fs::read_to_string(format!("{root}no_range_middle_20260924/frame_chunks.tsv")).unwrap();
    let observed_calls: Vec<(i8, usize)> = trace
        .lines()
        .filter(|line| line.starts_with("FRAME_CHUNK\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (f[2].parse().unwrap(), f[3].parse().unwrap())
        })
        .collect();
    // The (frame, chunk, scan ranges) units the search runs: the unit path's plan.
    let plan = super::super::search_seed::x_scan_plan(
        &[query],
        SubjectNt::Letters(&subject).into(),
        1,
        None,
        13,
        true,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &Default::default(),
        None,
        None,
    )
    .unwrap()
    .unwrap();
    let actual_ranges: Vec<_> = plan
        .subject
        .units
        .iter()
        .map(|(frame_index, chunk, ranges)| {
            (
                plan.subject.frames[*frame_index].frame,
                chunk.length,
                ranges.clone(),
            )
        })
        .collect();
    assert_eq!(
        actual_ranges
            .iter()
            .map(|(frame, length, _)| (*frame, *length))
            .collect::<Vec<_>>(),
        observed_calls
    );
    assert_eq!(
        actual_ranges,
        ncbi_chunk_range_rows(&format!("{root}no_range_middle_20260924/chunk_ranges.tsv"))
    );
    drop(plan);
    let expected: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("INIT\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            assert_eq!(f[1], "1");
            InitHsp {
                frame: 1,
                chunk_offset: 9_999_800,
                q_seed: f[3].parse().unwrap(),
                s_seed: f[4].parse().unwrap(),
                q_start: f[5].parse().unwrap(),
                s_start: f[6].parse().unwrap(),
                length: f[7].parse().unwrap(),
                score: f[8].parse().unwrap(),
            }
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // Retained long-fixture PARAM rows give cutoff 30.
    let actual = preliminary_stages(
        &[query],
        &subject,
        1,
        word_finder_profile(None, 13, 40, &[16], &[30], true, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .initial;
    assert_eq!(actual, expected);
}

/// The concatenated queries and their context offsets, as the unit path builds them.
fn x_concatenated(queries: &[&[u8]]) -> (Vec<u8>, Vec<i32>) {
    let mut sequence = Vec::new();
    let mut offsets = Vec::new();
    for query in queries {
        offsets.push(sequence.len() as i32);
        sequence.extend_from_slice(&encode_tblastn_lookup_query(query, None, false).aa_seq[1..]);
    }
    (sequence, offsets)
}

/// Every exact word of three residues the subject shares with the queries, in subject order
/// (the order of the subject scan).
fn x_word_pairs(
    query: &[u8],
    subject: &[u8],
) -> Vec<crate::algorithm::tblastx::blast_aascan::BlastOffsetPair> {
    let mut pairs = Vec::new();
    for s in 0..subject.len().saturating_sub(2) {
        for q in 0..query.len().saturating_sub(2) {
            if query[q..q + 3] == subject[s..s + 3] && query[q..q + 3].iter().all(|&r| r != 0) {
                pairs.push(crate::algorithm::tblastx::blast_aascan::BlastOffsetPair {
                    q_off: q as u32,
                    s_off: s as u32,
                });
            }
        }
    }
    pairs
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-172
// The unit path's kept state: one table carried through many units with `exit` after each, as
// NCBI's thread carries it through its subjects, gives every unit the hits of a cleared table
// (a table made afresh): from a new table, across the clear near INT4_MAX / 4, and when
// `with_table` takes it over for another query set of the same table length. The units are
// fragments of the queries with substitutions, so that two-hit extensions (and the `flag` path
// after a right extension) happen on many diagonals.
#[test]
fn x_two_hit_table_kept_across_units_gives_the_hits_of_a_cleared_table() {
    let mut seed = 0x2545_f491_u64;
    let mut next = move |bound: usize| {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((seed >> 33) as usize) % bound
    };
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    fn random(length: usize, next: &mut dyn FnMut(usize) -> usize) -> Vec<u8> {
        (0..length).map(|_| RESIDUES[next(20)]).collect()
    }
    let queries: Vec<Vec<u8>> = [150, 90, 210]
        .into_iter()
        .map(|length| random(length, &mut next))
        .collect();
    let mut units = Vec::new();
    for _ in 0..60 {
        let mut unit = random(next(40), &mut next);
        for _ in 0..1 + next(3) {
            let query = &queries[next(queries.len())];
            let start = next(query.len() - 30);
            let end = (start + 30 + next(80)).min(query.len());
            for &residue in &query[start..end] {
                unit.push(if next(100) < 15 {
                    RESIDUES[next(20)]
                } else {
                    residue
                });
            }
            unit.extend(random(next(30), &mut next));
        }
        units.push(encode_tblastn_lookup_query(&unit, None, false).aa_seq[1..=unit.len()].to_vec());
    }
    let x_dropoffs = [16, 16, 16];
    let cutoffs = [24, 24, 24];
    let window = 40;
    /// The hits of every unit in turn, each on a table made afresh by `fresh` (a cleared table).
    fn cleared_hits<'q>(
        fresh: impl Fn() -> XTwoHit<'q>,
        query: &[u8],
        units: &[Vec<u8>],
    ) -> Vec<Vec<InitHsp>> {
        units
            .iter()
            .map(|unit| {
                let mut hits = Vec::new();
                fresh().pairs(1, 0, unit, &x_word_pairs(query, unit), &mut hits);
                hits
            })
            .collect()
    }
    /// The hits of every unit in turn on one kept table (`exit` after each unit).
    fn kept_hits(two_hit: &mut XTwoHit<'_>, query: &[u8], units: &[Vec<u8>]) -> Vec<Vec<InitHsp>> {
        units
            .iter()
            .map(|unit| {
                let mut hits = Vec::new();
                two_hit.pairs(1, 0, unit, &x_word_pairs(query, unit), &mut hits);
                two_hit.exit(unit.len());
                hits
            })
            .collect()
    }
    let set_a: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let (query_a, offsets_a) = x_concatenated(&set_a);
    let new_a = |table| {
        XTwoHit::with_table(
            table,
            &query_a,
            &offsets_a,
            &x_dropoffs,
            &cutoffs,
            ScoringMatrix::Blosum62,
            window,
            3,
        )
    };
    let cleared = cleared_hits(|| new_a(None), &query_a, &units);
    let found: usize = cleared.iter().map(Vec::len).sum();
    assert!(found >= 60, "the units must give hits ({found})");

    let mut kept = new_a(None);
    assert_eq!(
        kept_hits(&mut kept, &query_a, &units),
        cleared,
        "kept table"
    );

    let mut near_limit = new_a(Some(kept.into_table()));
    let start = i32::MAX / 4 - 5_000;
    near_limit.diagonals.offset = start;
    assert_eq!(
        kept_hits(&mut near_limit, &query_a, &units),
        cleared,
        "across the clear"
    );
    assert!(
        near_limit.diagonals.offset < start,
        "the offset must have passed INT4_MAX / 4 and been cleared"
    );

    // Another query set of the same table length (the same queries in another order) takes
    // the table over.
    let set_b: Vec<&[u8]> = queries.iter().rev().map(Vec::as_slice).collect();
    let (query_b, offsets_b) = x_concatenated(&set_b);
    let new_b = |table| {
        XTwoHit::with_table(
            table,
            &query_b,
            &offsets_b,
            &x_dropoffs,
            &cutoffs,
            ScoringMatrix::Blosum62,
            window,
            3,
        )
    };
    let cleared_b = cleared_hits(|| new_b(None), &query_b, &units);
    let entries = near_limit.x_table_at();
    let table = near_limit.into_table();
    let mut taken_over = new_b(Some(table));
    assert_eq!(
        taken_over.x_table_at(),
        entries,
        "the table must be taken over"
    );
    assert_eq!(
        kept_hits(&mut taken_over, &query_b, &units),
        cleared_b,
        "taken over"
    );

    // A table of another length is not taken over.
    let (short_query, short_offsets) = x_concatenated(&set_a[..1]);
    let other = XTwoHit::with_table(
        Some(taken_over.into_table()),
        &short_query,
        &short_offsets,
        &x_dropoffs[..1],
        &cutoffs[..1],
        ScoringMatrix::Blosum62,
        window,
        3,
    );
    assert_eq!(other.diagonals.offset, window);
    assert_eq!(
        other.diagonals.entries.len(),
        (short_query.len() - 1 + window as usize).next_power_of_two()
    );

    // LOSAT_X_TBNDIAGPACK: a table of the other layout is not taken over either.
    for packed in [false, true] {
        let mut kept = XTwoHit::x_with_layout(
            packed,
            None,
            &query_b,
            &offsets_b,
            &x_dropoffs,
            &cutoffs,
            ScoringMatrix::Blosum62,
            window,
            3,
        );
        let _ = kept_hits(&mut kept, &query_b, &units);
        let at = kept.x_table_at();
        let other = XTwoHit::x_with_layout(
            !packed,
            Some(kept.into_table()),
            &query_b,
            &offsets_b,
            &x_dropoffs,
            &cutoffs,
            ScoringMatrix::Blosum62,
            window,
            3,
        );
        assert_ne!(other.x_table_at(), at, "packed {packed}: another layout");
        assert_eq!(
            other.diagonals.offset, window,
            "packed {packed}: a new table"
        );
        assert_eq!(other.diagonals.entries.packed(), !packed);
    }
}

// LOSAT_X_TBNDIAGPACK: the tables of TBLASTN's windows (40 with word size 3, 60 with 2) are
// packed when the switch is on, up to a window of 2^28 above the word size; others keep pairs.
#[test]
fn x_packs_takes_the_windows_whose_values_fit() {
    let on = x_tbn_diag_pack();
    for (window, word_size, fit) in [
        (40, 3, true),
        (60, 2, true),
        (4, 3, true),
        (1 << 28, 3, true),
        ((1 << 28) + 1, 3, false),
        (3, 3, false),
        (0, 0, false),
        (i32::MAX, 3, false),
    ] {
        assert_eq!(
            x_packs(window, word_size),
            on && fit,
            "window {window} word size {word_size}"
        );
    }
}

// `x_pairs_bf` against `pairs` on units made as in the test above (fragments
// of the queries with substitutions, so that extensions and flag-1 cells occur), on cleared
// tables, on a kept table and across the clear near INT4_MAX / 4, for windows 40 and
// wordsize + 1, with the branch-free loop on either layout of the table (LOSAT_X_TBNDIAGPACK)
// and `pairs` on `(last_hit, flag)` cells; each unit's offset pairs come in arrays of random
// lengths, as the scan fills them. The init HSPs and every cell must be equal after every unit.
// With another matrix, `x_pairs_bf_of` with the matrix's byte-table walks against `pairs` with
// the reference extension.
#[test]
fn x_pairs_bf_gives_the_hits_and_table_of_pairs() {
    let mut seed = 0x9e37_79b9_u64;
    let mut next = move |bound: usize| {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((seed >> 33) as usize) % bound
    };
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    fn random(length: usize, next: &mut dyn FnMut(usize) -> usize) -> Vec<u8> {
        (0..length).map(|_| RESIDUES[next(20)]).collect()
    }
    let queries: Vec<Vec<u8>> = [150, 90, 210]
        .into_iter()
        .map(|length| random(length, &mut next))
        .collect();
    let mut units = Vec::new();
    for _ in 0..60 {
        let mut unit = random(next(40), &mut next);
        for _ in 0..1 + next(3) {
            let query = &queries[next(queries.len())];
            let start = next(query.len() - 30);
            let end = (start + 30 + next(80)).min(query.len());
            for &residue in &query[start..end] {
                unit.push(if next(100) < 15 {
                    RESIDUES[next(20)]
                } else {
                    residue
                });
            }
            unit.extend(random(next(30), &mut next));
        }
        units.push(encode_tblastn_lookup_query(&unit, None, false).aa_seq[1..=unit.len()].to_vec());
    }
    let x_dropoffs = [16, 16, 16];
    let cutoffs = [24, 24, 24];
    let set: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let (query, offsets) = x_concatenated(&set);
    let mut found = 0;
    // too far, overlap, flag 1, candidate (the last two go to `pairs` in `x_pairs_bf`).
    let mut outcomes = [0u64; 4];
    for (matrix, window, packed) in [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Blosum45,
        ScoringMatrix::Pam250,
        ScoringMatrix::Identity,
    ]
    .into_iter()
    .flat_map(|m| [(m, 40, false), (m, 4, false), (m, 40, true), (m, 4, true)])
    {
        for (mode, start) in [
            ("cleared", None),
            ("kept", None),
            ("near the clear", Some(i32::MAX / 4 - 5_000)),
        ] {
            // The reference on `(last_hit, flag)` cells; the branch-free loop on either layout
            // (LOSAT_X_TBNDIAGPACK).
            let new_in = |packed| {
                XTwoHit::x_with_layout(
                    packed,
                    None,
                    &query,
                    &offsets,
                    &x_dropoffs,
                    &cutoffs,
                    matrix,
                    window,
                    3,
                )
            };
            let new = || {
                let mut reference = new_in(false);
                if matrix != ScoringMatrix::Blosum62 {
                    (reference.x_ungap_fast, reference.x_ungap_matrix) = (false, None);
                }
                reference
            };
            let new_bf = || {
                let branch_free = new_in(packed);
                // with the switch on (`x_tbn_ungap`), another matrix takes `x_pairs_bf_of`
                if crate::utils::xswitch::on("LOSAT_X_TBNUNGAP") {
                    assert_eq!(
                        branch_free.x_ungap_matrix().is_some(),
                        matrix != ScoringMatrix::Blosum62
                    );
                }
                branch_free
            };
            let mut reference = new();
            let mut branch_free = new_bf();
            if let Some(start) = start {
                reference.diagonals.offset = start;
                branch_free.diagonals.offset = start;
            }
            for (u, unit) in units.iter().enumerate() {
                if mode == "cleared" {
                    reference = new();
                    branch_free = new_bf();
                }
                let pairs = x_word_pairs(&query, unit);
                // The reference one hit at a time, its outcomes counted from the cell before.
                let mut expected = Vec::new();
                for pair in &pairs {
                    let mask = reference.x_table_len() - 1;
                    let (last_hit, flag) =
                        reference.x_cell((pair.q_off.wrapping_sub(pair.s_off) as usize) & mask);
                    let diff = pair.s_off as i32 - (last_hit - reference.diagonals.offset);
                    let outcome = if flag {
                        2
                    } else if diff >= window {
                        0
                    } else if diff < 3 {
                        1
                    } else {
                        3
                    };
                    outcomes[outcome] += 1;
                    reference.pairs(1, 0, unit, std::slice::from_ref(pair), &mut expected);
                }
                let mut got = Vec::new();
                let mut rest = &pairs[..];
                while !rest.is_empty() {
                    let (array, tail) = rest.split_at(1 + next(rest.len().min(50)));
                    match branch_free.x_ungap_matrix() {
                        None => branch_free.x_pairs_bf(1, 0, unit, array, &mut got),
                        Some(table) => {
                            branch_free.x_pairs_bf_of(table, 1, 0, unit, array, &mut got)
                        }
                    }
                    rest = tail;
                }
                assert_eq!(
                    got, expected,
                    "{matrix:?} window {window} packed {packed} {mode} unit {u}: init HSPs"
                );
                // After a clear (`Diagonals::finish_frame`), `Pair` cells hold NCBI's (-window,
                // false) while `Packed` cells hold (0, false); packed tables never run lead
                // subjects and both values lie at least a window behind every hit at offset 0 or
                // more.
                let reference_table: Vec<(i32, bool)> = if packed {
                    reference
                        .x_table()
                        .into_iter()
                        .map(|cell| {
                            if cell == (-window, false) {
                                (0, false)
                            } else {
                                cell
                            }
                        })
                        .collect()
                } else {
                    reference.x_table()
                };
                assert_eq!(
                    branch_free.x_table(),
                    reference_table,
                    "{matrix:?} window {window} packed {packed} {mode} unit {u}: table"
                );
                found += expected.len();
                if mode != "cleared" {
                    reference.exit(unit.len());
                    branch_free.exit(unit.len());
                }
            }
        }
    }
    assert!(found >= 800, "the units must give hits ({found})");
    assert!(
        outcomes.iter().all(|&count| count > 0),
        "outcomes {outcomes:?}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1088-1091
// EXPERIMENT (LOSAT_X_TBNUNGAP): the word finder with the byte-table walks gives, unit after
// unit on one kept table (`exit` after each), the init HSPs and the diagonal table of the word
// finder with the reference extension: on whole units and on units cut short (a chunk that ends
// inside a frame), with stop codons and X among the subject letters, and with a context whose
// X-drop the byte table does not cover (the reference runs for its candidates); with every matrix
// (BLOSUM62's walks or the matrix's table).
#[test]
fn x_byte_table_walks_give_the_reference_hits_and_table() {
    let mut seed = 0x9e37_79b9_u64;
    let mut next = move |bound: usize| {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((seed >> 33) as usize) % bound
    };
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    const SUBJECT_LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYX*";
    fn random(letters: &[u8], length: usize, next: &mut dyn FnMut(usize) -> usize) -> Vec<u8> {
        (0..length).map(|_| letters[next(letters.len())]).collect()
    }
    let queries: Vec<Vec<u8>> = [170, 80, 230]
        .into_iter()
        .map(|length| random(RESIDUES, length, &mut next))
        .collect();
    let mut units = Vec::new();
    for _ in 0..60 {
        let mut unit = random(SUBJECT_LETTERS, next(40), &mut next);
        for _ in 0..1 + next(3) {
            let query = &queries[next(queries.len())];
            let start = next(query.len() - 30);
            let end = (start + 30 + next(80)).min(query.len());
            for &residue in &query[start..end] {
                unit.push(if next(100) < 15 {
                    SUBJECT_LETTERS[next(SUBJECT_LETTERS.len())]
                } else {
                    residue
                });
            }
            unit.extend(random(SUBJECT_LETTERS, next(30), &mut next));
        }
        units.push(encode_tblastn_lookup_query(&unit, None, false).aa_seq[1..=unit.len()].to_vec());
    }
    let set: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    let (query, offsets) = x_concatenated(&set);
    let cutoffs = [24, 24, 24];
    for (matrix, x_dropoffs) in ScoringMatrix::ALL
        .into_iter()
        .flat_map(|matrix| [(matrix, [16, 16, 16]), (matrix, [16, 300, 7])])
    {
        let new = |fast: bool| {
            let mut two_hit =
                XTwoHit::with_table(None, &query, &offsets, &x_dropoffs, &cutoffs, matrix, 40, 3);
            two_hit.x_ungap_fast = fast;
            two_hit.x_ungap_matrix = (fast && matrix != ScoringMatrix::Blosum62)
                .then(|| x_ungapped_fast::matrix_table(matrix).unwrap());
            two_hit
        };
        let (mut fast, mut reference) = (new(true), new(false));
        let mut found = 0;
        for unit in &units {
            for end in [unit.len(), unit.len() * 2 / 3] {
                let unit = &unit[..end];
                let pairs = x_word_pairs(&query, unit);
                let (mut got, mut want) = (Vec::new(), Vec::new());
                fast.pairs(1, 0, unit, &pairs, &mut got);
                reference.pairs(1, 0, unit, &pairs, &mut want);
                assert_eq!(got, want, "{matrix:?} X-drops {x_dropoffs:?}: init HSPs");
                assert_eq!(
                    fast.diagonals.entries, reference.diagonals.entries,
                    "{matrix:?} X-drops {x_dropoffs:?}: diagonal table"
                );
                found += got.len();
                fast.exit(unit.len());
                reference.exit(unit.len());
            }
        }
        assert!(
            found >= 60,
            "{matrix:?} X-drops {x_dropoffs:?}: the units must give hits ({found})"
        );
    }
    // The rule: BLOSUM62 takes its own walks, another matrix the walks with its table; with
    // the switch off (`LOSAT_X_TBNUNGAP=0`) every matrix takes the reference.
    let on = crate::utils::xswitch::on("LOSAT_X_TBNUNGAP");
    let (b62, b45) = (
        x_tbn_ungap(ScoringMatrix::Blosum62),
        x_tbn_ungap(ScoringMatrix::Blosum45),
    );
    assert_eq!((b62.0, b62.1.is_some()), (on, false));
    assert_eq!((b45.0, b45.1.is_some()), (on, on));
}

/// `length` random residues from the generator state `seed`.
fn x_negative_letters(length: usize, seed: &mut u64) -> Vec<u8> {
    const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    (0..length)
        .map(|_| {
            *seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            RESIDUES[((*seed >> 33) as usize) % RESIDUES.len()]
        })
        .collect()
}

/// A subject's residues as the word finder reads them (the encoded letters, no sentinel).
fn x_negative_encoded(letters: &[u8]) -> Vec<u8> {
    encode_tblastn_lookup_query(letters, None, false).aa_seq[1..=letters.len()].to_vec()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65,109-157
// The table of a window below 0: `diag_array_length` is the smallest power of two, from 1, that is
// not below `query length + window` (so the table can be shorter than the query, and is one cell
// when the sum is 1 or less), `diag_mask` is that length less 1, `offset` starts at the window and
// the cells (`calloc`) at `last_hit` 0.
#[test]
fn a_negative_windows_table_has_ncbis_length_and_start() {
    for (query_length, window, cells) in [
        (5000, -1, 8192),
        (100, -1, 128),
        (100, -36, 64),
        (100, -37, 64),
        (100, -98, 2),
        (100, -99, 1),
        (100, -100, 1),
        (100, -5000, 1),
        (100, i32::MIN, 1),
    ] {
        let table = NegDiagonals::new(query_length, window);
        assert_eq!(
            table.cells.len(),
            cells,
            "query {query_length} window {window}"
        );
        assert_eq!(table.mask as usize, cells - 1);
        assert_eq!(table.offset, i64::from(window));
        assert_eq!(table.window, window);
        assert!(table.cells.iter().all(|&cell| cell == 0));
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:90-106,162-185
// `Blast_ExtendWordExit` adds the subject's length and the window to the offset (it falls with a
// window below the length), and from a quarter of `INT4_MAX` on sets it back to the window and every cell
// to `-window`; the offset of a window that overflows `Int4` is taken in `i64`.
#[test]
fn a_negative_windows_table_moves_its_offset_and_clears_at_ncbis_limit() {
    let mut table = NegDiagonals::new(100, -40);
    table.exit(1000);
    assert_eq!(table.offset, -40 + 1000 - 40);
    table.exit(0);
    assert_eq!(table.offset, -40 + 1000 - 40 - 40);
    table.cells[3] = 77;
    table.offset = i64::from(i32::MAX / 4) - 1;
    table.exit(10);
    assert_eq!(table.offset, i64::from(i32::MAX / 4) - 1 + 10 - 40);
    assert_eq!(table.cells[3], 77, "below the limit nothing is cleared");
    table.offset = i64::from(i32::MAX / 4);
    table.exit(10);
    assert_eq!(table.offset, -40);
    assert!(table.cells.iter().all(|&cell| cell == 40));

    let mut wide = NegDiagonals::new(100, i32::MIN);
    wide.exit(7);
    assert_eq!(wide.offset, 2 * i64::from(i32::MIN) + 7);
    wide.cells[0] = 5;
    wide.offset = i64::from(i32::MAX / 4);
    wide.exit(1);
    assert_eq!(wide.offset, i64::from(i32::MIN));
    assert_eq!(
        wide.cells[0], 0,
        "-window of INT_MIN is 2^31, whose low 31 bits are 0"
    );
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:51-55
// `last_hit` is a signed bit field of 31 bits.
#[test]
fn the_last_hit_of_a_cell_keeps_its_low_31_bits_as_a_signed_number() {
    for (value, stored) in [
        (5, 5),
        (-5, -5),
        ((1 << 30) - 1, (1 << 30) - 1),
        (1 << 30, -(1 << 30)),
        (-(1 << 30), -(1 << 30)),
        (-(1 << 30) - 1, (1 << 30) - 1),
        (1 << 31, 0),
        ((1 << 40) + 7, 7),
    ] {
        assert_eq!(last_hit_31(value), stored, "{value}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:773-803
// A table that starts with `offset` = -40 and `last_hit` 0 holds `last_hit - offset` = 40 on every
// diagonal: the one-hit word finder extends no word whose subject offset is below 40 (`diff =
// s_off - (last_hit - offset)` below 0). A call that ended before (`exit` after 100 letters, with
// no hit) leaves `last_hit - offset` = -20 on the diagonals it did not touch, so the next call
// extends from its first word on: the table carries from one call to the next.
#[test]
fn a_negative_windows_one_hit_finder_skips_the_words_before_its_window_and_carries_its_table() {
    let mut seed = 0x51ed_270b_u64;
    let query_letters = x_negative_letters(120, &mut seed);
    let (query, offsets) = x_concatenated(&[&query_letters]);
    let subject = x_negative_encoded(&query_letters[10..90]);
    let pairs = x_word_pairs(&query, &subject);
    let x_dropoffs = [16];
    let cutoffs = [24];
    let window = -40;
    let finder = XTwoHit::with_table(
        None,
        &query,
        &offsets,
        &x_dropoffs,
        &cutoffs,
        ScoringMatrix::Blosum62,
        window,
        3,
    );
    let seeds = |table: &mut NegDiagonals| -> Vec<u32> {
        let mut hits = Vec::new();
        finder
            .x_pairs_one_hit_negative(
                table,
                crate::algorithm::tblastn::score::MatrixScore,
                1,
                0,
                &subject,
                &pairs,
                &mut hits,
            )
            .unwrap();
        hits.iter().map(|hit| hit.s_seed).collect()
    };
    let mut fresh = NegDiagonals::new(query.len() - 1, window);
    let first = seeds(&mut fresh);
    assert_eq!(first.first(), Some(&40), "{first:?}");
    assert!(first.iter().all(|&s| s >= 40), "{first:?}");
    assert!(fresh.cells.iter().any(|&cell| cell != 0));

    let mut carried = NegDiagonals::new(query.len() - 1, window);
    carried.exit(100);
    assert_eq!(carried.offset, -40 + 100 - 40);
    let second = seeds(&mut carried);
    assert_eq!(second.first(), Some(&0), "{second:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:765-769
// A table shorter than the query (62 of 70 letters rounds up to 64 cells) gives the diagonals
// whose numbers differ by a multiple of its length one cell: an extension that reached offset 100
// on diagonal 100 blocks a word at offset 46 on diagonal 36.
#[test]
fn a_negative_windows_table_shorter_than_the_query_shares_cells_between_diagonals() {
    use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair;
    let mut seed = 0x2f6b_8e03_u64;
    let query_letters = x_negative_letters(70, &mut seed);
    let (query, offsets) = x_concatenated(&[&query_letters]);
    let subject = x_negative_encoded(&x_negative_letters(200, &mut seed));
    let x_dropoffs = [16];
    let cutoffs = [i32::MAX];
    let window = -8;
    let finder = XTwoHit::with_table(
        None,
        &query,
        &offsets,
        &x_dropoffs,
        &cutoffs,
        ScoringMatrix::Blosum62,
        window,
        3,
    );
    let mut table = NegDiagonals::new(query.len() - 1, window);
    assert_eq!(table.cells.len(), 64);
    let mut run = |q_off, s_off| {
        let mut hits = Vec::new();
        finder
            .x_pairs_one_hit_negative(
                &mut table,
                crate::algorithm::tblastn::score::MatrixScore,
                1,
                0,
                &subject,
                &[BlastOffsetPair { q_off, s_off }],
                &mut hits,
            )
            .unwrap();
    };
    run(0, 100);
    let written: Vec<usize> = table
        .cells
        .iter()
        .enumerate()
        .filter(|(_, &cell)| cell != 0)
        .map(|(index, _)| index)
        .collect();
    assert_eq!(written, [36]);
    let cell = table.cells[36];
    let mut run = |q_off, s_off| {
        let mut hits = Vec::new();
        finder
            .x_pairs_one_hit_negative(
                &mut table,
                crate::algorithm::tblastn::score::MatrixScore,
                1,
                0,
                &subject,
                &[BlastOffsetPair { q_off, s_off }],
                &mut hits,
            )
            .unwrap();
    };
    run(10, 46);
    assert_eq!(
        table.cells[36], cell,
        "the word at offset 46 is not extended"
    );
    assert_eq!(table.cells.iter().filter(|&&cell| cell != 0).count(), 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:42-65
// A window below 0 of any size (down to INT_MIN) runs the preliminary stage of a subject (its
// frames are units of one table) without a panic: no word reaches the cutoff, so no HSP is saved.
#[test]
fn a_negative_window_runs_the_preliminary_stage_of_a_subject() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/ambiguity_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    for window in [-1, -40, -5000, i32::MIN] {
        let stages = preliminary_stages(
            &[query],
            subject,
            1,
            word_finder_profile(
                None,
                13,
                window,
                &[16],
                &[i32::MAX],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap();
        assert!(stages.initial.is_empty(), "window {window}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:40-66,87-107,131-144,162-173
/// `HardTable` on TBLASTN's `Diagonals` and `Pair` clear to `(-window, false)`: `ncbi_table`
/// sets NCBI's offset and unwritten-cell value, `clear_table` and `finish_frame` at INT4_MAX / 4
/// clear `Pair` cells to `(-window, false)`, and `move_offset` checks Int4 overflow.
#[test]
fn diagonals_hard_table_and_pair_clear_match_ncbi() {
    use crate::algorithm::common::hard_chunk_walk::{HardTable, NcbiStart, TableStart};

    let mut table = Diagonals::new(32, 40, false);
    for (start, want_cell) in [
        (
            NcbiStart {
                offset: 536_875_512,
                cleared: false,
            },
            (0, false),
        ),
        (
            NcbiStart {
                offset: 40,
                cleared: true,
            },
            (-40, false),
        ),
    ] {
        TableStart::New(start).apply(&mut table, 40).unwrap();
        assert_eq!(table.offset, start.offset);
        assert_eq!(
            table.entries,
            XDiagCells::Pair(vec![want_cell; table.entries.len()])
        );
    }

    // Reaching INT4_MAX / 4 clears Pair cells to (-window, false) and resets offset to window.
    table.offset = i32::MAX / 4;
    table.finish_frame(5_000_000);
    assert_eq!(table.offset, 40);
    assert_eq!(
        table.entries,
        XDiagCells::Pair(vec![(-40, false); table.entries.len()])
    );

    // clear_table and move_offset (including Int4 overflow check).
    table.move_offset(100).unwrap();
    assert_eq!(table.offset, 140);
    assert!(table.move_offset(i32::MAX).is_err());
    table.clear_table(40);
    assert_eq!(
        (table.offset, &table.entries),
        (
            40,
            &XDiagCells::Pair(vec![(-40, false); table.entries.len()])
        )
    );

    // Packed clear keeps 0 (a packed table never runs a lead subject).
    let mut packed = Diagonals::new(32, 40, true);
    packed.offset = i32::MAX / 4;
    packed.finish_frame(5_000_000);
    assert_eq!(packed.offset, 40);
    assert_eq!(
        packed.entries,
        XDiagCells::Packed(vec![0; packed.entries.len()])
    );
}
