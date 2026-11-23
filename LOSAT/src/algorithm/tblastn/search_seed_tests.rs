use super::*;
use crate::algorithm::blastp::encoding::ncbistdaa_to_ascii;
use crate::algorithm::common::translated_chunks::unmasked_translated_chunks;
use crate::algorithm::tblastn::search_gapped::{preliminary_stages, scan_profile};
use crate::sequence::packed_subject::SubjectNt;
use anyhow::bail;
use std::collections::HashMap;
use std::fs;

fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let text = fs::read_to_string(path).unwrap();
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(id) = line.strip_prefix('>') {
            records.push((
                id.split_whitespace().next().unwrap().to_string(),
                Vec::new(),
            ));
        } else {
            records.last_mut().unwrap().1.extend(line.as_bytes());
        }
    }
    records
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:68-96:
// one context per protein query; query offsets advance by length + 1.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:478-505:
// scansub emits the ordered global offset pairs for each subject frame.
// This comparison uses the unmodified candidate stream from the pinned
// local -subject oracle; no final-output filter is applied to candidates.
#[test]
fn multi_query_candidate_contexts_and_order_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let refs: Vec<&[u8]> = queries
        .iter()
        .map(|(_, sequence)| sequence.as_slice())
        .collect();
    let frames: Vec<_> = refs
        .iter()
        .map(|query| vec![encode_protein_query_frame_with_seg(query, None)])
        .collect();
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (_, contexts) = build_ncbi_lookup(&frames, 13, &karlin, false);
    let actual_contexts: Vec<_> = contexts
        .iter()
        .map(|context| {
            (
                context.q_idx as usize,
                context.frame_base,
                context.aa_len,
                context.is_valid,
            )
        })
        .collect();
    let trace = fs::read_to_string(format!("{root}query_context.tsv")).unwrap();
    let expected_contexts: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("QUERY_CONTEXT\t0\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[7] == "1",
            )
        })
        .collect();
    assert_eq!(actual_contexts, expected_contexts);

    let trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let frame_order = [1, 2, 3, -1, -2, -3];
    let expected: Vec<Seed> = trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            Seed {
                frame: frame_order[f[1].parse::<usize>().unwrap()],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            }
        })
        .collect();
    let actual = preliminary_stages(
        &refs,
        subject,
        1,
        scan_profile(refs.len(), None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds;
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!(
            "first multi-query candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
            actual.get(first),
            expected.get(first),
            actual.len(),
            expected.len()
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:227-245,438-465:
// word_length and neighbor scores come from the selected options/matrix.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:478-505:
// the subject scanner emits global query offsets in its observed order.
#[test]
fn alternate_matrix_word2_candidates_and_contexts_match_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
    );
    let queries = read_fasta(&format!("{root}query.faa"));
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let frames: Vec<_> = refs
        .iter()
        .map(|query| vec![encode_protein_query_frame_with_seg(query, None)])
        .collect();
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum45);
    let (lookup, contexts) =
        build_ncbi_lookup_for_profile(&frames, 16, &karlin, false, ScoringMatrix::Blosum45, 2);
    assert_eq!(lookup.word_length, 2);
    let trace = fs::read_to_string(format!("{root}query_context.tsv")).unwrap();
    let expected_contexts: Vec<_> = trace
        .lines()
        .filter(|line| line.starts_with("QUERY_CONTEXT\t0\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse::<usize>().unwrap(),
                f[4].parse::<i32>().unwrap(),
                f[5].parse::<usize>().unwrap(),
                f[7] == "1",
            )
        })
        .collect();
    let actual_contexts: Vec<_> = contexts
        .iter()
        .map(|ctx| (ctx.q_idx as usize, ctx.frame_base, ctx.aa_len, ctx.is_valid))
        .collect();
    assert_eq!(actual_contexts, expected_contexts);
    let trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
    let frame_order = [1, 2, 3, -1, -2, -3];
    let expected: Vec<Seed> = trace
        .lines()
        .filter(|line| line.starts_with("CAND\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            Seed {
                frame: frame_order[f[1].parse::<usize>().unwrap()],
                chunk_offset: 0,
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            }
        })
        .collect();
    let actual = preliminary_stages(
        &refs,
        subject,
        1,
        scan_profile(refs.len(), None, 16, false, ScoringMatrix::Blosum45, 2),
    )
    .unwrap()
    .seeds;
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!("first BLOSUM45/word-2 candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
            actual.get(first), expected.get(first), actual.len(), expected.len());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
// Seeds `CRandom` with the nucleotide base length; ambiguous ncbi4na bases are picked once
// for the compressed ncbi2na preliminary subject.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:767-775,804-812
#[test]
fn ncbi_preliminary_subject_frames_match_all_fixture_bytes() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let records = read_fasta(&format!("{root}subjects.fna"));
    let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let mut plans = HashMap::new();
    let trace = fs::read_to_string(format!("{root}ncbi_frame_trace.tsv")).unwrap();
    let mut compared = 0;
    for line in trace.lines().skip(1) {
        let columns: Vec<_> = line.split('\t').collect();
        let subject_id = columns[1];
        let frame: i8 = columns[2].parse().unwrap();
        // The frames the search translates: the unit path's plan.
        let plan = plans.entry(subject_id).or_insert_with(|| {
            x_scan_plan(
                &[query],
                SubjectNt::Letters(by_id[subject_id]).into(),
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
            .unwrap()
        });
        let actual = plan
            .subject
            .frames
            .iter()
            .find(|f| f.frame == frame)
            .unwrap();
        let actual_hex: String = actual.aa_seq[1..actual.aa_seq.len() - 1]
            .iter()
            .map(|base| format!("{base:02x}"))
            .collect();
        assert_eq!(
            actual.aa_len,
            columns[3].parse::<usize>().unwrap(),
            "{subject_id}/{frame}"
        );
        assert_eq!(actual_hex, columns[4], "{subject_id}/{frame}");
        compared += 1;
    }
    assert_eq!(compared, 78);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:478-505
#[test]
fn ncbi_all_candidates_match_in_set_and_order() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let records = read_fasta(&format!("{root}subjects.fna"));
    let trace = fs::read_to_string(format!("{root}ncbi_candidates.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (
                fields[1].to_string(),
                fields[2].parse::<i8>().unwrap(),
                fields[4].parse::<u32>().unwrap(),
                fields[5].parse::<u32>().unwrap(),
            )
        })
        .collect();
    let mut actual = Vec::new();
    for (id, subject) in &records {
        let seeds = preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        for seed in seeds {
            actual.push((
                id.clone(),
                seed.frame,
                seed.query_offset,
                seed.subject_offset,
            ));
        }
    }
    assert_eq!(actual, expected);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:816-831
// On the first plus frame (context 0), converts subject scan ranges to amino-acid
// coordinates; on the first minus frame (context 3), reverses them for minus frames.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:51-59
#[test]
fn ncbi_lowercase_mask_candidates_match_in_set_and_order() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/lowercase_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let records = read_fasta(&format!("{root}subjects.fna"));
    let trace = fs::read_to_string(format!("{root}ncbi_candidates.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].to_string(),
                f[2].parse::<i8>().unwrap(),
                f[4].parse::<u32>().unwrap(),
                f[5].parse::<u32>().unwrap(),
            )
        })
        .collect();
    let mut actual = Vec::new();
    for (id, subject) in &records {
        for seed in preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(1, None, 13, true, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds
        {
            actual.push((
                id.clone(),
                seed.frame,
                seed.query_offset,
                seed.subject_offset,
            ));
        }
    }
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!("first lowercase candidate mismatch at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
               actual.get(first), expected.get(first), actual.len(), expected.len());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
// Seeds `CRandom` with the nucleotide base length; each ambiguous 4na value consumes one `GetRand()` draw.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
#[test]
fn ncbi_ambiguity_spectrum_frames_and_candidates_match() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/ambiguity_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
    // The frames the search translates: the unit path's plan.
    let plan = x_scan_plan(
        &[query],
        SubjectNt::Letters(subject).into(),
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
    let actual_frames = &plan.subject.frames;
    let frame_trace = fs::read_to_string(format!("{root}ncbi_frame_trace.tsv")).unwrap();
    for line in frame_trace.lines().skip(1) {
        let f: Vec<_> = line.split('\t').collect();
        let frame: i8 = f[2].parse().unwrap();
        let actual = actual_frames.iter().find(|x| x.frame == frame).unwrap();
        let hex: String = actual.aa_seq[1..actual.aa_seq.len() - 1]
            .iter()
            .map(|base| format!("{base:02x}"))
            .collect();
        assert_eq!(hex, f[4], "frame {frame}");
    }
    let trace = fs::read_to_string(format!("{root}ncbi_candidates.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[2].parse::<i8>().unwrap(),
                f[4].parse::<u32>().unwrap(),
                f[5].parse::<u32>().unwrap(),
            )
        })
        .collect();
    let actual: Vec<_> = preliminary_stages(
        &[query],
        subject,
        1,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds
    .into_iter()
    .map(|x| (x.frame, x.query_offset, x.subject_offset))
    .collect();
    assert_eq!(actual, expected);
    assert_eq!(actual.len(), 181);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1101
// New NCBI -outfmt 6 output supplies the frame and aligned subject AA.
#[test]
fn six_frame_translation_and_seed_coverage() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let records = read_fasta(&format!("{root}subjects.fna"));
    let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
    let output = fs::read_to_string(format!("{root}raw_isolation_fields.out")).unwrap();
    let mut frames_checked = 0;
    for line in output.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        let subject_id = fields[1];
        if !(subject_id.starts_with("plus") || subject_id.starts_with("minus")) {
            continue;
        }
        let expected_frame: i8 = fields[7].parse().unwrap();
        let subject = by_id[subject_id];
        // The frames the search translates: the unit path's plan.
        let plan = x_scan_plan(
            &[query],
            SubjectNt::Letters(subject).into(),
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
        let frame = plan
            .subject
            .frames
            .iter()
            .find(|f| f.frame == expected_frame)
            .unwrap();
        let translated: String = frame.aa_seq[1..frame.aa_seq.len() - 1]
            .iter()
            .map(|&aa| ncbistdaa_to_ascii(aa))
            .collect();
        assert!(translated.contains(fields[9]), "{subject_id}: {translated}");
        let seeds = preliminary_stages(
            &[query],
            subject,
            1,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        assert!(
            seeds.iter().any(|seed| seed.frame == expected_frame
                && seed.query_offset == 0
                && seed.subject_offset == 0),
            "{subject_id}"
        );
        frames_checked += 1;
    }
    assert_eq!(frames_checked, 6);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:200-234
// Collects initial hits via `s_BlastAaWordFinder_TwoHit`.
// The comparison-only C probe captures NCBI's saved seed pair. Each
// unambiguous saved pair must exist in the Rust scan for the same frame.
#[test]
fn ncbi_saved_init_seeds_exist_in_rust_word_scan() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let records = read_fasta(&format!("{root}subjects.fna"));
    let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
    let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
    let mut compared = 0;
    for line in trace.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        let subject_id = fields[1];
        if subject_id == "ambiguous" {
            continue;
        }
        let frame: i8 = fields[2].parse().unwrap();
        let q_seed: u32 = fields[4].parse().unwrap();
        let s_seed: u32 = fields[5].parse().unwrap();
        let seeds = preliminary_stages(
            &[query],
            by_id[subject_id],
            1,
            scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        assert!(
            seeds.iter().any(|seed| seed.frame == frame
                && seed.query_offset == q_seed
                && seed.subject_offset == s_seed),
            "{subject_id}: ({q_seed},{s_seed})"
        );
        compared += 1;
    }
    assert_eq!(compared, 10);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381
#[test]
fn seg_suppresses_low_complexity_query_seeds() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/run_20260923/"
    );
    let query = &read_fasta(&format!("{root}query_low.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}low_subject.fna"))[0].1;
    let unmasked = preliminary_stages(
        &[query],
        subject,
        1,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds;
    let masked = preliminary_stages(
        &[query],
        subject,
        1,
        scan_profile(
            1,
            Some(&SegParams::default()),
            13,
            false,
            ScoringMatrix::Blosum62,
            3,
        ),
    )
    .unwrap()
    .seeds;
    println!(
        "fixture=low_complexity seg_off_seeds={} seg_on_seeds={}",
        unmasked.len(),
        masked.len()
    );
    assert!(!unmasked.is_empty());
    assert!(masked.is_empty());
    assert!(fs::read(format!("{root}low_query_seg_on.out"))
        .unwrap()
        .is_empty());
    assert!(!fs::read(format!("{root}low_query_seg_off.out"))
        .unwrap()
        .is_empty());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:747-775,804-841
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:478-505
#[test]
fn code32_api_oracle_candidates_match_in_set_and_order() {
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
    let trace = fs::read_to_string(format!("{trace_root}ncbi_candidates.tsv")).unwrap();
    let expected: Vec<_> = trace
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[0].parse::<i8>().unwrap(),
                f[2].parse::<u32>().unwrap(),
                f[3].parse::<u32>().unwrap(),
            )
        })
        .collect();
    let actual: Vec<_> = preliminary_stages(
        &[query],
        subject,
        32,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds
    .into_iter()
    .map(|x| (x.frame, x.query_offset, x.subject_offset))
    .collect();
    assert_eq!(actual, expected);
    assert_eq!(actual.len(), 199);
}

// NCBI reference (598d8ae6): c++/src/objects/seqfeat/gc.prt:340-347
#[test]
fn code_32_changes_search_words_and_ambiguity_resolves() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_a/fixtures/"
    );
    let query = &read_fasta(&format!("{root}query.faa"))[0].1;
    let subject = &read_fasta(&format!("{root}subject_code32.fna"))[0].1;
    let code1 = preliminary_stages(
        &[query],
        subject,
        1,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds;
    let code32 = preliminary_stages(
        &[query],
        subject,
        32,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds;
    println!(
        "fixture=code32 code1_seeds={} code32_seeds={}",
        code1.len(),
        code32.len()
    );
    assert_ne!(code1, code32);
    assert!(code32
        .iter()
        .any(|seed| seed.frame == 1 && seed.query_offset == 0 && seed.subject_offset == 0));
    assert!(preliminary_stages(
        &[query],
        b"ATN",
        32,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3)
    )
    .is_ok());
    assert!(preliminary_stages(
        &[query],
        b"atg",
        32,
        scan_profile(1, None, 13, false, ScoringMatrix::Blosum62, 3)
    )
    .is_ok());
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:68-96:
// each protein query receives one context; offsets advance length + 1,
// and the X-only third context is invalid for the lookup.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:478-505 and
// c++/src/algo/blast/core/blast_engine.c:804-841:
// scansub emits every pair in frame/chunk order over both long chunks.
#[test]
fn long_multi_query_every_candidate_and_context_matches_ncbi() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let fixture = format!("{root}long_multi_query_20260924/");
    let queries = read_fasta(&format!("{fixture}query.faa"));
    let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
    let plus1 = read_fasta(&format!("{root}run_20260923/subjects.fna"))
        .into_iter()
        .find(|(id, _)| id == "plus1")
        .unwrap()
        .1;
    let mut subject = Vec::with_capacity(15_000_962);
    for _ in 0..4_999_950 {
        subject.extend_from_slice(b"ATG");
    }
    subject.extend_from_slice(&plus1);
    for _ in 0..250 {
        subject.extend_from_slice(b"ATG");
    }
    assert_eq!(subject.len(), 15_000_962);
    let frames: Vec<_> = refs
        .iter()
        .map(|query| vec![encode_protein_query_frame_with_seg(query, None)])
        .collect();
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (_, contexts) = build_ncbi_lookup(&frames, 13, &karlin, false);
    let actual_contexts: Vec<_> = contexts
        .iter()
        .map(|c| (c.q_idx as usize, c.frame_base, c.aa_len, c.is_valid))
        .collect();
    let gapped_trace = fs::read_to_string(format!("{fixture}gapped_events.tsv")).unwrap();
    let expected_contexts: Vec<_> = gapped_trace
        .lines()
        .filter(|line| line.starts_with("QUERY_CONTEXT\t0\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
                f[7] == "1",
            )
        })
        .collect();
    assert_eq!(actual_contexts, expected_contexts);
    assert_eq!(
        actual_contexts,
        [(0, 0, 120, true), (1, 121, 70, true), (2, 192, 120, false)]
    );
    let trace = fs::read_to_string(format!("{fixture}candidate.stderr")).unwrap();
    assert_eq!(
        trace.lines().filter(|l| l.starts_with("PARAM\t")).count(),
        12
    );
    assert_eq!(trace.lines().filter(|l| l.starts_with("END\t")).count(), 12);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:478-505:
    // scansub reads the same per-call scan
    // parameters and chunk state before returning each ordered pair.
    let params: Vec<_> = trace
        .lines()
        .filter(|l| l.starts_with("PARAM\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].parse::<usize>().unwrap(),
                f[2].parse::<i32>().unwrap(),
                f[3].parse::<i32>().unwrap(),
                f[4].parse::<i32>().unwrap(),
            )
        })
        .collect();
    assert_eq!(params.len(), 12);
    for (call, row) in params.iter().enumerate() {
        assert_eq!(*row, (call, 16, 16, 28));
    }
    let frame_order = [1, 2, 3, -1, -2, -3];
    let chunk_rows = fs::read_to_string(format!("{fixture}frame_chunks.tsv")).unwrap();
    let observed_chunks: Vec<_> = chunk_rows
        .lines()
        .filter(|line| line.starts_with("FRAME_CHUNK\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].parse::<usize>().unwrap(),
                f[2].parse::<i8>().unwrap(),
                f[3].parse::<usize>().unwrap(),
            )
        })
        .collect();
    assert_eq!(observed_chunks.len(), 12);
    for (call, row) in observed_chunks.iter().enumerate() {
        assert_eq!(
            *row,
            (
                call,
                frame_order[call / 2],
                if call % 2 == 0 { 5_000_000 } else { 420 }
            )
        );
    }
    let mut candidate_counts = [0usize; 12];
    let expected: Vec<_> = trace
        .lines()
        .filter(|l| l.starts_with("CAND\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            let call: usize = f[1].parse().unwrap();
            candidate_counts[call] += 1;
            Seed {
                frame: frame_order[call / 2],
                chunk_offset: if call.is_multiple_of(2) { 0 } else { 4_999_900 },
                query_offset: f[3].parse().unwrap(),
                subject_offset: f[4].parse().unwrap(),
            }
        })
        .collect();
    assert_eq!(expected.len(), 394);
    let ends: Vec<_> = trace
        .lines()
        .filter(|l| l.starts_with("END\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (
                f[1].parse::<usize>().unwrap(),
                f[2].parse::<usize>().unwrap(),
            )
        })
        .collect();
    assert_eq!(ends.len(), 12);
    for (call, row) in ends.iter().enumerate() {
        assert_eq!(*row, (call, candidate_counts[call]));
    }
    assert!(expected.iter().all(|seed| seed.query_offset < 192));
    let actual = preliminary_stages(
        &refs,
        &subject,
        1,
        scan_profile(refs.len(), None, 13, false, ScoringMatrix::Blosum62, 3),
    )
    .unwrap()
    .seeds;
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(left, right)| left != right)
            .unwrap_or(actual.len().min(expected.len()));
        panic!("first long multi-query candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
               actual.get(first), expected.get(first), actual.len(), expected.len());
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:804-841:
// Searches each translated frame context via `s_BlastSearchEngineOneContext`.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:495-505:
#[test]
fn long_subject_every_candidate_matches_ncbi_in_chunk_order() {
    let root = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/"
    );
    let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
    let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
    let frame_order = [1, 2, 3, -1, -2, -3];
    for (fixture, mask_lowercase) in [
        ("long_chunk_20260924", false),
        ("masked_chunk_boundary_20260924", true),
    ] {
        let query = &read_fasta(&format!("{root}{fixture}/query.faa"))[0].1;
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 15_000_962);
        if mask_lowercase {
            subject[14_999_700..15_000_000].make_ascii_lowercase();
        }
        let trace = fs::read_to_string(format!(
            "{root}ordered_candidates_20260924/{fixture}.stderr"
        ))
        .unwrap();
        assert_eq!(
            trace
                .lines()
                .filter(|line| line.starts_with("PARAM\t"))
                .count(),
            12
        );
        assert_eq!(
            trace
                .lines()
                .filter(|line| line.starts_with("END\t"))
                .count(),
            12
        );
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                Seed {
                    frame: frame_order[call / 2],
                    chunk_offset: if call.is_multiple_of(2) { 0 } else { 4_999_900 },
                    query_offset: f[3].parse().unwrap(),
                    subject_offset: f[4].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &[query.as_slice()],
            &subject,
            1,
            scan_profile(1, None, 13, mask_lowercase, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!(
                "{fixture}: first candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                actual.get(first), expected.get(first), actual.len(), expected.len()
            );
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-264,478-500,804-841:
    // s_GetNextSubjectChunk skips a middle chunk with no soft range;
    // each surviving chunk calls WordFinder in frame order.
    for (fixture, early_hit) in [
        ("no_range_middle_20260924", false),
        ("no_range_two_hits_20260924", true),
    ] {
        let query = &read_fasta(&format!("{root}{fixture}/query.faa"))[0].1;
        let mut subject = if early_hit {
            let mut nt = b"ATG".repeat(4_998_450);
            nt.extend_from_slice(plus1);
            nt.push(b'A');
            nt.extend_from_slice(&b"ATG".repeat(10_000_050 - 4_998_450 - 121));
            nt.extend_from_slice(plus1);
            nt
        } else {
            let mut nt = b"ATG".repeat(10_000_050);
            nt.extend_from_slice(plus1);
            nt
        };
        subject.extend_from_slice(&b"ATG".repeat(250));
        assert_eq!(subject.len(), 30_001_262);
        subject[14_997_000..30_000_000].make_ascii_lowercase();
        let trace = fs::read_to_string(format!(
            "{root}ordered_candidates_20260924/{fixture}.stderr"
        ))
        .unwrap();
        let frames = [1, 1, 2, 2, 3, 3, -1, -1, -1, -2, -2, -2, -3, -3, -3];
        let offsets = [
            0, 9_999_800, 0, 9_999_800, 0, 9_999_800, 0, 4_999_900, 9_999_800, 0, 4_999_900,
            9_999_800, 0, 4_999_900, 9_999_800,
        ];
        let range_trace = fs::read_to_string(format!("{root}{fixture}/chunk_ranges.tsv")).unwrap();
        let calls: Vec<_> = range_trace
            .lines()
            .filter(|line| line.starts_with("RANGE_CALL\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (f[2].parse::<i8>().unwrap(), f[3].parse::<usize>().unwrap())
            })
            .collect();
        assert_eq!(calls.len(), 15);
        for (call, (frame, length)) in calls.iter().enumerate() {
            assert_eq!(*frame, frames[call]);
            assert_eq!(
                *length,
                if offsets[call] == 9_999_800 {
                    620
                } else {
                    5_000_000
                }
            );
        }
        assert_eq!(
            trace
                .lines()
                .filter(|line| line.starts_with("PARAM\t"))
                .count(),
            15
        );
        assert_eq!(
            trace
                .lines()
                .filter(|line| line.starts_with("END\t"))
                .count(),
            15
        );
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("CAND\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                Seed {
                    frame: frames[call],
                    chunk_offset: offsets[call],
                    query_offset: f[3].parse().unwrap(),
                    subject_offset: f[4].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &[query.as_slice()],
            &subject,
            1,
            scan_profile(1, None, 13, true, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .seeds;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!(
                "{fixture}: first candidate difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                actual.get(first), expected.get(first), actual.len(), expected.len()
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-264:
// split only when offset + MAX_DBSEQ_LEN < hard-range right;
// the next offset overlaps the first chunk by DBSEQ_CHUNK_OVERLAP.
#[test]
fn unmasked_long_translated_chunks_match_ncbi_wordfinder_calls() {
    let rows = unmasked_translated_chunks(5_000_320);
    assert_eq!(
        rows,
        vec![
            TranslatedChunk::new(0, 5_000_000),
            TranslatedChunk::new(4_999_900, 420),
        ]
    );
    assert_eq!(
        unmasked_translated_chunks(5_000_000),
        vec![TranslatedChunk::new(0, 5_000_000),]
    );
    assert_eq!(
        unmasked_translated_chunks(5_000_001),
        vec![
            TranslatedChunk::new(0, 5_000_000),
            TranslatedChunk::new(4_999_900, 101),
        ]
    );
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../docs/evidence/tlosan_stage_c/long_chunk_20260924/frame_chunks.tsv"
    );
    let trace = fs::read_to_string(path).unwrap();
    let observed: Vec<(i8, usize)> = trace
        .lines()
        .filter(|line| line.starts_with("FRAME_CHUNK\t"))
        .map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (f[2].parse().unwrap(), f[3].parse().unwrap())
        })
        .collect();
    let expected: Vec<_> = [1, 2, 3, -1, -2, -3]
        .into_iter()
        .flat_map(|frame| rows.iter().map(move |chunk| (frame, chunk.length)))
        .collect();
    assert_eq!(observed, expected);
}

// No NCBI counterpart: test helper (pseudo-random nucleotide subject); it does not change any
// value NCBI computes.
fn x_test_subject(length: usize) -> Vec<u8> {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect()
}

// No NCBI counterpart: test of the plan cache (one plan, charged once, when batches race); it
// does not change any value NCBI computes.
// Subject plans: batches that reach a subject at the same time get one
// plan. Its bytes are charged to the shared limit once (`keeps`, by the
// first of them) and it is made once.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn x_subject_plan_is_made_and_charged_once_when_batches_race() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let subject = x_test_subject(9_000);
    let code = GeneticCode::try_from_id(1).unwrap();
    let side = XSubjectSide::new(3, 1, false);
    let slot = side.slot(1, 1, false).unwrap();
    let builds = AtomicUsize::new(0);
    let build = || -> Result<std::sync::Arc<XSubjectPlan>> {
        builds.fetch_add(1, Ordering::Relaxed);
        // Long enough for every thread to arrive while the first builds.
        std::thread::sleep(std::time::Duration::from_millis(50));
        Ok(std::sync::Arc::new(XSubjectPlan {
            frames: generate_frames(&subject, &code),
            units: Vec::new(),
        }))
    };
    let threads = 8;
    let barrier = std::sync::Barrier::new(threads);
    let plans: Vec<_> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    barrier.wait();
                    assert!(side.keeps(slot, subject.len()));
                    side.plan(slot, subject.len(), build).unwrap()
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    assert_eq!(builds.load(Ordering::Relaxed), 1);
    for plan in &plans[1..] {
        assert!(std::sync::Arc::ptr_eq(&plans[0], plan));
    }
    let bytes = XSubjectSide::plan_bytes(&plans[0]);
    // Three frames per strand of 3000, 2999 and 2999 residues, two sentinels each.
    assert_eq!(bytes, 2 * (3_000 + 2_999 + 2_999) + 12);
    assert_eq!(XSubjectSide::plan_bytes_of(subject.len()), bytes);
    assert_eq!(side.kept_bytes(), bytes);
    // A later batch reads the same plan and adds nothing.
    assert!(side.keeps(slot, subject.len()));
    let again = side.plan(slot, subject.len(), build).unwrap();
    assert!(std::sync::Arc::ptr_eq(&plans[0], &again));
    assert_eq!(
        (builds.load(Ordering::Relaxed), side.kept_bytes()),
        (1, bytes)
    );
    // Another gencode or masking option does not use these plans.
    assert!(side.slot(1, 4, false).is_none() && side.slot(1, 1, true).is_none());
}

// No NCBI counterpart: test of the plan cache through x_scan_plan; it does not change any value
// NCBI computes.
// The whole stage: `x_scan_plan` from eight threads on one subject.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn x_scan_plans_of_racing_batches_share_one_subject_plan() {
    let subject = x_test_subject(30_000);
    // Test protein: made from the first residues of E. coli AsnA (RefSeq NP_418200.1) and other letters.
    let query: &[u8] = b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ";
    let side = XSubjectSide::new(1, 1, false);
    let threads = 8;
    let barrier = std::sync::Barrier::new(threads);
    let plans: Vec<_> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    barrier.wait();
                    x_scan_plan(
                        &[query],
                        crate::sequence::packed_subject::SubjectNt::Letters(&subject).into(),
                        1,
                        None,
                        13,
                        false,
                        ScoringMatrix::Blosum62,
                        3,
                        &std::sync::OnceLock::new(),
                        &Default::default(),
                        Some((&side, 0)),
                        None,
                    )
                    .unwrap()
                    .unwrap()
                    .subject
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for plan in &plans[1..] {
        assert!(std::sync::Arc::ptr_eq(&plans[0], plan));
    }
    assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes(&plans[0]));
    assert!(!plans[0].units.is_empty());
}

// No NCBI counterpart: test of when the plan cache exists (XSubjectSide::for_batches); it does
// not change any value NCBI computes.
// Subject plans: a search of fewer than four query batches keeps no
// plan and every use makes its own (one batch has no later batch to read a
// kept plan; two and three give the time it would save for its memory);
// four or more batches keep the plans. The plan is the same either way.
#[test]
fn x_subject_plans_are_kept_only_from_four_batches() {
    assert_eq!(XSubjectSide::KEEP_FROM_BATCHES, 4);
    for batches in [0, 1, 2, 3] {
        assert!(XSubjectSide::for_batches(true, batches, 1, 1, false).is_none());
    }
    for batches in [4, 5, 9] {
        assert!(XSubjectSide::for_batches(true, batches, 1, 1, false).is_some());
    }
    for batches in [0, 1, 2, 3, 4, 9] {
        assert!(XSubjectSide::for_batches(false, batches, 1, 1, false).is_none());
    }
    let subject = x_test_subject(9_000);
    let query: &[u8] = b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ";
    let scan = |side: Option<&XSubjectSide>| {
        x_scan_plan(
            &[query],
            crate::sequence::packed_subject::SubjectNt::Letters(&subject).into(),
            1,
            None,
            13,
            false,
            ScoringMatrix::Blosum62,
            3,
            &std::sync::OnceLock::new(),
            &Default::default(),
            side.map(|side| (side, 0)),
            None,
        )
        .unwrap()
        .unwrap()
        .subject
    };
    // Three batches: two uses make two plans, and nothing else holds them.
    let three = XSubjectSide::for_batches(true, 3, 1, 1, false);
    let (first, second) = (scan(three.as_ref()), scan(three.as_ref()));
    assert!(!std::sync::Arc::ptr_eq(&first, &second));
    assert_eq!(std::sync::Arc::strong_count(&first), 1);
    // Four batches: the second use reads the plan the first one made.
    let four = XSubjectSide::for_batches(true, 4, 1, 1, false).unwrap();
    let (kept, read) = (scan(Some(&four)), scan(Some(&four)));
    assert!(std::sync::Arc::ptr_eq(&kept, &read));
    assert_eq!(four.kept_bytes(), XSubjectSide::plan_bytes(&kept));
    // The kept plan is the plan a use makes for itself, in every field. For subject plans
    // `aa_seq_nomask` is None and `seg_masks` is empty (and `mask_lowercase` is off here), so
    // the comparison is trivial for those two; it covers the residues, lengths and units.
    assert_eq!(
        crate::algorithm::tblastn::x_subject_translate::x_plan_fields_difference(&first, &kept),
        None
    );
    assert!(!first.units.is_empty());
}

// No NCBI counterpart: test of when LOSAT_X_TBNKEEPFEW keeps the plans; it does not change any
// value NCBI computes.
// The plans are kept when they hold no more than `threads` plans of the longest subject: one
// subject at every thread count, a few long ones at enough threads, many short ones never.
#[test]
fn x_subject_plans_are_kept_within_the_threads() {
    fn within(lengths: impl IntoIterator<Item = usize>, threads: usize) -> bool {
        let (longest, total) = lengths
            .into_iter()
            .fold((0usize, 0usize), |(longest, total), length| {
                (usize::max(longest, length), total.saturating_add(length))
            });
        XSubjectSide::kept_within_threads(longest, total, threads)
    }
    assert!(within([4_700_000usize], 1));
    assert!(within([4_700_000usize], 32));
    assert!(!within([3_000usize, 3_000], 1));
    assert!(within([3_000usize, 3_000], 2));
    assert!(within([3_000usize, 2_999, 1], 2));
    assert!(!within([3_000usize, 3_000, 1], 2));
    // Three fly chromosomes (76.9 Mb): kept at 3 or more threads (three times the longest,
    // 28.1 Mb, is 84.3 Mb), not at 1 or 2.
    let fly = [23_513_712usize, 25_286_936, 28_110_227];
    assert!(!within(fly, 1));
    assert!(!within(fly, 2));
    assert!(within(fly, 3));
    assert!(within(fly, 16));
    // 28,036 16S records of about 1,450 bases: never at the thread counts of a search.
    assert!(!within(std::iter::repeat_n(1_450usize, 28_036), 32));
    assert!(!within(std::iter::repeat_n(1_450usize, 28_036), 1_024));
    // No subjects keep nothing; no threads count as one.
    assert!(within(std::iter::empty::<usize>(), 1));
    assert!(!within([3_000usize, 3_000], 0));
    assert!(within([usize::MAX, 1], 2));
}

// No NCBI counterpart: test of the size limit of the plan cache; it does not change any value
// NCBI computes.
// A plan that does not fit leaves the total as it was; a failure is
// returned and nothing is kept for the subject.
#[test]
fn x_subject_plan_limit_counts_kept_plans_only() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let kept = AtomicUsize::new(0);
    assert!(XSubjectSide::reserve(&kept, 60, 100));
    assert!(!XSubjectSide::reserve(&kept, 60, 100));
    assert_eq!(kept.load(Ordering::Relaxed), 60);
    assert!(XSubjectSide::reserve(&kept, 40, 100));
    assert!(!XSubjectSide::reserve(&kept, 1, 100));
    assert!(XSubjectSide::reserve(&kept, 0, 100));
    assert_eq!(kept.load(Ordering::Relaxed), 100);
    assert!(!XSubjectSide::reserve(
        &AtomicUsize::new(usize::MAX),
        1,
        usize::MAX
    ));

    assert!(XSubjectSide::all_fit(100_000_000, 100));
    assert!(XSubjectSide::all_fit(
        (XSubjectSide::KEEP_LIMIT - 12) / 2,
        1
    ));
    assert!(!XSubjectSide::all_fit(
        (XSubjectSide::KEEP_LIMIT - 12) / 2 + 1,
        1
    ));
    assert!(!XSubjectSide::all_fit(usize::MAX, 2));

    let side = XSubjectSide::new(1, 1, false);
    let slot = side.slot(0, 1, false).unwrap();
    assert!(side.keeps(slot, 30_000));
    assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes_of(30_000));
    let failing = || -> Result<std::sync::Arc<XSubjectPlan>> { bail!("no plan") };
    assert_eq!(
        side.plan(slot, 30_000, failing).err().unwrap().to_string(),
        "no plan"
    );
    // The bytes of the plan that failed are given back; the next use makes (and here
    // fails) its own.
    assert_eq!(side.kept_bytes(), 0);
    assert_eq!(slot.plan.get().map(Option::is_none), Some(true));
    assert!(side.plan(slot, 30_000, failing).is_err());
}

// No NCBI counterpart: test of the length a side keeps for a `-db` sequence; it does not
// change any value NCBI computes.
// The first `keeps` of a subject keeps its length with the decision, kept or refused, for the
// search's genetic code and masking option only; a longer record's is saturated.
#[test]
fn x_side_keeps_the_length_of_its_decision() {
    let side = XSubjectSide::new(3, 1, false);
    assert_eq!(side.length(0, 1, false), None);
    assert!(side.keeps(side.slot(0, 1, false).unwrap(), 30_000));
    assert_eq!(side.length(0, 1, false), Some(30_000));
    assert_eq!(side.length(0, 2, false), None);
    assert_eq!(side.length(0, 1, true), None);
    // Refused (over the limit): the length is kept all the same; a later `keeps` keeps the
    // first decision and length.
    let length = XSubjectSide::KEEP_LIMIT;
    assert!(!side.keeps(side.slot(1, 1, false).unwrap(), length));
    assert_eq!(side.length(1, 1, false), Some(length));
    assert!(!side.keeps(side.slot(1, 1, false).unwrap(), 10));
    assert_eq!(side.length(1, 1, false), Some(length));
    assert!(side.keeps(side.slot(0, 1, false).unwrap(), 10));
    assert_eq!(side.length(0, 1, false), Some(30_000));
    let longest = i32::MAX as usize;
    side.keeps(side.slot(2, 1, false).unwrap(), longest);
    assert_eq!(side.length(2, 1, false), Some(longest));
    let side = XSubjectSide::new(1, 1, false);
    side.keeps(side.slot(0, 1, false).unwrap(), usize::MAX);
    assert_eq!(side.length(0, 1, false), Some(longest));
}

// No NCBI counterpart: test of the size a subject's plan is charged before it is made; it does
// not change any value NCBI computes.
// The bytes of a plan from the subject's length alone are the bytes of the plan made from its
// letters or its packed bases, at every length (no frame below three bases; two sentinels per
// frame), and never above the bound `all_fit` uses.
#[test]
fn x_plan_bytes_of_a_length_are_the_bytes_of_its_plan() {
    use crate::algorithm::tblastn::x_subject_translate::{x_packed_subject_plan, x_subject_plan};
    use crate::sequence::packed_subject::PackedNucleotide;
    for length in (1..=40).chain([1_000, 4_999, 100_003]) {
        let subject = x_test_subject(length);
        let bytes = XSubjectSide::plan_bytes_of(length);
        let letters = x_subject_plan(&subject, 1, false).unwrap();
        assert_eq!(
            XSubjectSide::plan_bytes(&letters),
            bytes,
            "{length} bases, letters"
        );
        let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
        let packed = x_packed_subject_plan(&store, 1).unwrap();
        assert_eq!(
            XSubjectSide::plan_bytes(&packed),
            bytes,
            "{length} bases, packed"
        );
        assert!(bytes <= 2 * length + 12, "{length} bases");
    }
    assert_eq!(XSubjectSide::plan_bytes_of(0), 0);
    assert_eq!(XSubjectSide::plan_bytes_of(2), 0);
    assert_eq!(XSubjectSide::plan_bytes_of(3), 6);
    assert_eq!(XSubjectSide::plan_bytes_of(5), 6 + 6 + 6);
    assert_eq!(XSubjectSide::plan_bytes_of(usize::MAX), usize::MAX);
    // A subject too short to be searched (x_scan_plan returns None before the plan) is not
    // charged: the pre-read decides only for a subject long enough for a lean plan, and
    // x_scan_plan decides after its length check.
    let side = XSubjectSide::new(1, 1, false);
    assert!(!x_lean_subject(Some((&side, 0)), 10, 1, false));
    assert_eq!(side.kept_bytes(), 0);
    assert!(x_scan_plan(
        &[b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ".as_slice()],
        crate::sequence::packed_subject::SubjectNt::Letters(b"ACGTACGTAC").into(),
        1,
        None,
        13,
        false,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &Default::default(),
        Some((&side, 0)),
        None,
    )
    .unwrap()
    .is_none());
    assert_eq!(side.kept_bytes(), 0);
}

// No NCBI counterpart: test of the plans a side does not keep; it does not change any value
// NCBI computes.
// A plan the side's limit refuses is made per use, as without a side: without residues when
// the subject is long (LOSAT_X_TBNUNITXLATE), with its frames when it is short; the decision
// is made once per subject, before the plan, and `x_lean_subject` reads it. A plan that fits
// is kept, with its residues, and read by the next use.
#[test]
fn x_plans_the_side_does_not_keep_are_made_per_use_as_without_a_side() {
    use crate::algorithm::tblastn::x_subject_translate::{
        x_plan_fields_difference, x_units_translate,
    };
    use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};
    let query: &[u8] = b"MKVLAAGIVGLLLAQWERTYHHNDCFPSMKTAYIAKQRQISFVKSHFSRQ";
    let scan = |side: Option<&XSubjectSide>, store: &PackedNucleotide| {
        x_scan_plan(
            &[query],
            SubjectNt::Packed(store).into(),
            1,
            None,
            13,
            false,
            ScoringMatrix::Blosum62,
            3,
            &std::sync::OnceLock::new(),
            &Default::default(),
            side.map(|side| (side, 0)),
            None,
        )
        .unwrap()
        .unwrap()
        .subject
    };
    let long = PackedNucleotide::from_subject_letters(&x_test_subject(3_200_000)).expect("letters");
    let short = PackedNucleotide::from_subject_letters(&x_test_subject(9_000)).expect("letters");
    for (store, lean) in [(&long, x_units_translate(3_200_000)), (&short, false)] {
        let length = store.len();
        // A limit of no bytes refuses every plan: each use makes its own, the plan a use
        // makes without a side in every field.
        let side = XSubjectSide::with_limit(1, 1, false, 0);
        assert_eq!(x_lean_subject(Some((&side, 0)), length, 1, false), lean);
        let (first, second) = (scan(Some(&side), store), scan(Some(&side), store));
        assert!(!std::sync::Arc::ptr_eq(&first, &second));
        assert_eq!(side.kept_bytes(), 0);
        assert!(!side.keeps(side.slot(0, 1, false).unwrap(), length));
        assert_eq!(first.frames[0].aa_seq.is_empty(), lean, "{length} bases");
        assert_eq!(x_plan_fields_difference(&first, &scan(None, store)), None);
        // The default limit keeps both: the plan has its residues and the next use reads it.
        let side = XSubjectSide::new(1, 1, false);
        assert!(!x_lean_subject(Some((&side, 0)), length, 1, false));
        let (kept, read) = (scan(Some(&side), store), scan(Some(&side), store));
        assert!(std::sync::Arc::ptr_eq(&kept, &read));
        assert!(!kept.frames[0].aa_seq.is_empty());
        assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes(&kept));
        assert_eq!(side.kept_bytes(), XSubjectSide::plan_bytes_of(length));
    }
}

// A small deterministic generator (xorshift64) for the compressed scan tests.
struct XRng(u64);

impl XRng {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % n as u64) as usize
    }
}

/// The hits of cell `index` of a compressed table, in NCBI's scan order.
fn x_cell_hits(table: &XCompressedLookup, index: usize) -> Vec<i32> {
    match table {
        XCompressedLookup::Reference(table) => {
            if table.backbone[index].num_used == 0 {
                Vec::new()
            } else {
                table.query_offsets_in_scan_order(index)
            }
        }
        XCompressedLookup::Flat(table) => table.cell(index).to_vec(),
    }
}

/// The pairs of every unit of `plan`: what `x_scan_unit` hands to its caller, and what NCBI's
/// compressed scan finds (blast_aascan.c:264-409): at each offset of each range whose word has
/// a compressed letter at every position, the hits of the word's cell (index sum c_k A^k, the
/// index the scan's preshift forms), in order.
fn x_compressed_unit_pairs(plan: &XScanPlan) -> (Vec<(u32, u32)>, Vec<(u32, u32)>) {
    let XPlanLookup::Compressed(table) = &plan.lookup else {
        panic!("a compressed plan")
    };
    let (word, alphabet, scaled) = match &**table {
        XCompressedLookup::Reference(t) => (
            t.word_length,
            t.compressed_alphabet_size,
            t.scaled_compress_table_u8,
        ),
        XCompressedLookup::Flat(t) => (
            t.word_length,
            t.compressed_alphabet_size,
            t.scaled_compress_table_u8,
        ),
    };
    let top = alphabet.pow(word as u32 - 1);
    let mut scanned = Vec::new();
    let mut expected = Vec::new();
    let mut pairs = Vec::new();
    for (unit, (frame_index, chunk, ranges)) in plan.subject.units.iter().enumerate() {
        let residues = &plan.subject.frames[*frame_index].aa_seq[1 + chunk.offset..];
        x_scan_unit(plan, unit, residues, 0, &mut pairs, |filled| {
            scanned.extend(filled.iter().map(|pair| (pair.q_off, pair.s_off)));
        });
        for (range_index, &(left, right)) in ranges.iter().enumerate() {
            let mut end = right - word;
            if range_index == 0 && end < left {
                end = left;
            }
            for s in left..=end.max(left - 1) {
                let letters: Vec<i32> = (0..word)
                    .map(|k| {
                        let letter = residues.get((s + k) as usize).copied().unwrap_or(0);
                        scaled[usize::from(letter)]
                    })
                    .collect();
                if letters.iter().any(|&c| c < 0) {
                    continue;
                }
                let index: i32 = letters
                    .iter()
                    .enumerate()
                    .map(|(k, &c)| c / top * alphabet.pow(k as u32))
                    .sum();
                for hit in x_cell_hits(table, index as usize) {
                    expected.push((hit as u32, s as u32));
                }
            }
        }
    }
    (scanned, expected)
}

fn x_random_protein(rng: &mut XRng, length: usize) -> Vec<u8> {
    const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYXBZU*";
    (0..length)
        .map(|_| LETTERS[rng.below(LETTERS.len())])
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:264-409 and
// c++/src/algo/blast/core/aa_ungapped.c:494-516
// The unit scan of the compressed table (word sizes 5 to 7) finds the hits of NCBI's scan, in
// its order, with any offset array (one that holds only the longest cell resumes at every
// cell), on the flat table and, for word size 5, on the reference table, which is the flat
// one; a threshold of 100 puts only the exact words in the table.
#[test]
fn x_compressed_unit_scan_finds_ncbis_hits_in_order() {
    use crate::algorithm::tblastx::lookup::compressed::try_build_blosum62_compressed_lookup;
    use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_assert_same_compressed_table, x_try_build_flat_compressed_lookup, XFlatBuild,
    };
    let mut rng = XRng(0x9e37_79b9_7f4a_7c15);
    let queries: Vec<Vec<u8>> = (0..3)
        .map(|i| x_random_protein(&mut rng, 40 + 23 * i))
        .collect();
    let query_refs: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
    const BASES: &[u8] = b"ACGTACGTACGTACGTACGTN";
    let mut subject: Vec<u8> = (0..2400).map(|_| BASES[rng.below(BASES.len())]).collect();
    // Fragments of the queries, back-translated (one codon for each letter, `NNN` for the
    // others) into the subject, forward and reverse-complemented, so that the exact words of
    // a threshold of 100 have hits too.
    let codon = |letter: u8| -> &'static [u8] {
        match letter {
            b'A' => b"GCT",
            b'C' => b"TGT",
            b'D' => b"GAT",
            b'E' => b"GAA",
            b'F' => b"TTT",
            b'G' => b"GGT",
            b'H' => b"CAT",
            b'I' => b"ATT",
            b'K' => b"AAA",
            b'L' => b"CTG",
            b'M' => b"ATG",
            b'N' => b"AAT",
            b'P' => b"CCG",
            b'Q' => b"CAG",
            b'R' => b"CGT",
            b'S' => b"TCT",
            b'T' => b"ACT",
            b'V' => b"GTT",
            b'W' => b"TGG",
            b'Y' => b"TAT",
            _ => b"NNN",
        }
    };
    for (index, query) in queries.iter().enumerate() {
        let mut fragment: Vec<u8> = query[5..35]
            .iter()
            .flat_map(|&l| codon(l).to_vec())
            .collect();
        if index == 1 {
            fragment.reverse();
            for base in &mut fragment {
                *base = match *base {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    other => other,
                };
            }
        }
        let at = 300 + 700 * index;
        subject.splice(at..at, fragment);
    }
    for (word, threshold) in [
        (5usize, 1100),
        (5, 2000),
        (5, 10000),
        (6, 2100),
        (6, 10000),
        (7, 2025),
        (7, 10000),
    ] {
        let frames: Vec<_> = query_refs
            .iter()
            .map(|query| vec![encode_tblastn_lookup_query(query, None, false)])
            .collect();
        let (query, locations) =
            crate::algorithm::tblastx::lookup::x_lookup_query_for_profile(&frames, word);
        let tables = [0..locations.len()];
        let XFlatBuild::Table(flat) =
            x_try_build_flat_compressed_lookup(word, threshold, &query, &locations, &tables)
                .unwrap()
        else {
            panic!("word {word} threshold {threshold}: no flat table")
        };
        if threshold == 10000 {
            // Only the exact words: one hit per query offset whose word has compressed letters.
            let hits: usize = (0..flat.starts.len() - 1).map(|i| flat.cell(i).len()).sum();
            let words: usize = locations
                .iter()
                .map(|&(from, to)| {
                    (from..=to - word as i32 + 1)
                        .filter(|&o| {
                            (0..word).all(|k| {
                                flat.scaled_compress_table_u8[usize::from(query[(o as usize) + k])]
                                    >= 0
                            })
                        })
                        .count()
                })
                .sum();
            assert_eq!(hits, words, "word {word}: exact words only");
        }
        // The side holds the flat table, so that the plan takes it whatever LOSAT_X_BPCLUT says.
        let flat = std::sync::Arc::new(XCompressedLookup::Flat(flat));
        let compressed = XCompressedSide::default();
        assert!(compressed.table.set(std::sync::Arc::clone(&flat)).is_ok());
        let mut plan = x_scan_plan(
            &query_refs,
            SubjectNt::Letters(&subject).into(),
            1,
            None,
            threshold,
            false,
            ScoringMatrix::Blosum62,
            word,
            &std::sync::OnceLock::new(),
            &compressed,
            None,
            None,
        )
        .unwrap()
        .unwrap();
        let XPlanLookup::Compressed(planned) = &plan.lookup else {
            panic!("word {word}: the plan of the compressed table")
        };
        assert!(std::sync::Arc::ptr_eq(planned, &flat));
        let mut tables = vec![flat];
        if word == 5 {
            let reference = try_build_blosum62_compressed_lookup(
                word,
                threshold,
                &query,
                &locations,
                &[0..locations.len()],
            )
            .unwrap()
            .unwrap();
            let XCompressedLookup::Flat(flat) = &*tables[0] else {
                unreachable!()
            };
            x_assert_same_compressed_table(flat, &reference);
            tables.push(std::sync::Arc::new(XCompressedLookup::Reference(reference)));
        }
        let planned_capacity = plan.pair_capacity;
        let mut first = None;
        for table in tables {
            let longest = match &*table {
                XCompressedLookup::Reference(t) => t.longest_chain,
                XCompressedLookup::Flat(t) => t.longest_chain,
            };
            plan.lookup = XPlanLookup::Compressed(table);
            for capacity in [planned_capacity, longest.max(1) as usize] {
                plan.pair_capacity = capacity;
                let (scanned, expected) = x_compressed_unit_pairs(&plan);
                assert_eq!(scanned, expected, "word {word} threshold {threshold}");
                assert!(
                    !scanned.is_empty(),
                    "word {word} threshold {threshold}: no hits"
                );
                match &first {
                    None => first = Some(scanned),
                    Some(first) => assert_eq!(&scanned, first, "word {word}"),
                }
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:496-501 and
// c++/src/algo/blast/core/blast_engine.c:305-308
// The last chunk of a soft-masked frame of more than MAX_DBSEQ_LEN letters can have one range,
// its last letter (`-lcase_masking`, a subject lowercase but its last codon). The raised word
// there runs past the frame's closing NULLB, where `subject` ends: the unit does not scan it
// (before, a debug build stopped on the read past the slice, and release read the byte after
// it). A first range of the last two letters still scans its word, which ends at the NULLB,
// and the whole frame finds every exact word of the query.
#[test]
fn x_scan_unit_does_not_scan_a_raised_word_past_the_frame() {
    let query = b"MKVLAAGIVGLLLAAGCSSHKEEAKSGGLW";
    let codon = |letter: u8| -> &'static [u8] {
        match letter {
            b'A' => b"GCT",
            b'C' => b"TGT",
            b'E' => b"GAA",
            b'G' => b"GGT",
            b'H' => b"CAT",
            b'I' => b"ATT",
            b'K' => b"AAA",
            b'L' => b"CTG",
            b'M' => b"ATG",
            b'S' => b"TCT",
            b'V' => b"GTT",
            b'W' => b"TGG",
            other => panic!("no codon for {}", other as char),
        }
    };
    // Frame +1 of the subject is the query.
    let subject: Vec<u8> = query.iter().flat_map(|&l| codon(l).to_vec()).collect();
    let mut plan = x_scan_plan(
        &[query.as_slice()],
        SubjectNt::Letters(&subject).into(),
        1,
        None,
        11,
        false,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &XCompressedSide::default(),
        None,
        None,
    )
    .unwrap()
    .unwrap();
    assert!(
        matches!(plan.lookup, XPlanLookup::Aa(_)),
        "the standard table"
    );
    let frame = plan.subject.frames[0].clone();
    assert_eq!((frame.frame, frame.aa_len), (1, query.len()));
    let length = query.len() as i32;
    let chunk = TranslatedChunk::new(0, query.len());
    plan.subject = std::sync::Arc::new(XSubjectPlan {
        frames: vec![frame],
        units: vec![
            (0, chunk, vec![(length - 1, length)]),
            (0, chunk, vec![(length - 2, length)]),
            (0, chunk, vec![(0, length)]),
        ],
    });
    // The frame from the chunk's first residue to its closing NULLB, as a unit reads it.
    let residues = &plan.subject.frames[0].aa_seq[1..];
    assert_eq!(residues.len(), query.len() + 1);
    let mut pairs = Vec::new();
    let calls: Vec<Vec<Vec<(u32, u32)>>> = (0..3)
        .map(|unit| {
            let mut calls = Vec::new();
            x_scan_unit(&plan, unit, residues, 0, &mut pairs, |filled| {
                calls.push(filled.iter().map(|p| (p.q_off, p.s_off)).collect());
            });
            calls
        })
        .collect();
    assert!(
        calls[0].is_empty(),
        "the word past the NULLB: {:?}",
        calls[0]
    );
    // The word of the last two letters and the NULLB is in no cell.
    assert_eq!(calls[1], vec![Vec::new()]);
    let whole = calls[2].concat();
    for s_off in 0..=(length - 3) as u32 {
        assert!(
            whole.contains(&(s_off, s_off)),
            "the exact word at {s_off}: {whole:?}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:494-505 and 760-765
#[test]
fn x_scan_unit_of_one_hit_does_not_scan_a_unit_whose_first_range_is_short() {
    let query = b"MKVLAAGIVGLLLAAGCSSHKEEAKSGGLW";
    let codon = |letter: u8| -> &'static [u8] {
        match letter {
            b'A' => b"GCT",
            b'C' => b"TGT",
            b'E' => b"GAA",
            b'G' => b"GGT",
            b'H' => b"CAT",
            b'I' => b"ATT",
            b'K' => b"AAA",
            b'L' => b"CTG",
            b'M' => b"ATG",
            b'S' => b"TCT",
            b'V' => b"GTT",
            b'W' => b"TGG",
            other => panic!("no codon for {}", other as char),
        }
    };
    let subject: Vec<u8> = query.iter().flat_map(|&l| codon(l).to_vec()).collect();
    let mut plan = x_scan_plan(
        &[query.as_slice()],
        SubjectNt::Letters(&subject).into(),
        1,
        None,
        11,
        false,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &XCompressedSide::default(),
        None,
        None,
    )
    .unwrap()
    .unwrap();
    let frame = plan.subject.frames[0].clone();
    let length = query.len() as i32;
    let chunk = TranslatedChunk::new(0, query.len());
    plan.subject = std::sync::Arc::new(XSubjectPlan {
        frames: vec![frame],
        units: vec![
            // A first range of two letters, then the rest.
            (0, chunk, vec![(0, 2), (5, length)]),
            // A first range of a word, then the rest.
            (0, chunk, vec![(0, 3), (5, length)]),
        ],
    });
    let residues = &plan.subject.frames[0].aa_seq[1..];
    let mut pairs = Vec::new();
    let mut scan = |unit: usize, one_hit: bool| {
        let mut words = Vec::new();
        let mut on_pairs = |filled: &[BlastOffsetPair]| {
            words.extend(filled.iter().map(|p| (p.q_off, p.s_off)));
        };
        if one_hit {
            x_scan_unit_of::<false>(&plan, unit, residues, 0, &mut pairs, &mut on_pairs);
        } else {
            x_scan_unit(&plan, unit, residues, 0, &mut pairs, &mut on_pairs);
        }
        words
    };
    // The two-hit scan raises the short first range to its word at 0 and goes on.
    let two_hit = scan(0, false);
    assert!(two_hit.contains(&(0, 0)), "{two_hit:?}");
    assert!(two_hit.contains(&(5, 5)), "{two_hit:?}");
    // The one-hit scan never starts: no word of the unit.
    assert_eq!(scan(0, true), Vec::new());
    // A first range of a word: both scans are the same.
    assert_eq!(scan(1, true), scan(1, false));
    assert!(scan(1, true).contains(&(0, 0)));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:300-306,
// c++/src/algo/blast/core/blast_aascan.c:76-77, :110 and c++/src/algo/blast/core/masksubj.inl:50-58
// Masks stored out of order (probe M4b-T2, 2026-11-23): a later range of a chunk starts before
// the chunk. The unit scan reads the frame's letters from `scan_lead` letters before the chunk
// and gives each hit at the chunk's offsets (below 0 before it, as NCBI's Uint4); a first range
// raised to one word on the frame's last letter finds no word there (the frame's NULLB), and the
// scan goes on with the later ranges.
#[test]
fn x_scan_unit_reads_a_range_before_its_chunk_at_the_chunk_s_offsets() {
    let query = b"MKVLAAGIVGLLLAAGCSSHKEEAKSGGLW";
    let codon = |letter: u8| -> &'static [u8] {
        match letter {
            b'A' => b"GCT",
            b'C' => b"TGT",
            b'E' => b"GAA",
            b'G' => b"GGT",
            b'H' => b"CAT",
            b'I' => b"ATT",
            b'K' => b"AAA",
            b'L' => b"CTG",
            b'M' => b"ATG",
            b'S' => b"TCT",
            b'V' => b"GTT",
            b'W' => b"TGG",
            other => panic!("no codon for {}", other as char),
        }
    };
    // Frame +1 of the subject is the query.
    let subject: Vec<u8> = query.iter().flat_map(|&l| codon(l).to_vec()).collect();
    let mut plan = x_scan_plan(
        &[query.as_slice()],
        SubjectNt::Letters(&subject).into(),
        1,
        None,
        11,
        false,
        ScoringMatrix::Blosum62,
        3,
        &std::sync::OnceLock::new(),
        &XCompressedSide::default(),
        None,
        None,
    )
    .unwrap()
    .unwrap();
    assert!(
        matches!(plan.lookup, XPlanLookup::Aa(_)),
        "the standard table"
    );
    let frame = plan.subject.frames[0].clone();
    let length = query.len() as i32;
    // The frame's 30 letters scanned whole, and a chunk of its last 20 (from letter 10) with the
    // ranges that masks stored out of order give it.
    let whole = TranslatedChunk::new(0, query.len());
    let chunk = TranslatedChunk::new(10, 20);
    plan.subject = std::sync::Arc::new(XSubjectPlan {
        frames: vec![frame],
        units: vec![
            (0, whole, vec![(0, length)]),
            (0, chunk, vec![(2, 6), (12, -4), (-10, 20)]),
            (0, chunk, vec![(19, 19), (19, -6), (-8, 20)]),
        ],
    });
    assert_eq!(
        plan.subject
            .units
            .iter()
            .map(|(_, _, ranges)| crate::algorithm::common::translated_chunks::scan_lead(ranges))
            .collect::<Vec<_>>(),
        vec![0, 10, 8]
    );
    let aa_seq = plan.subject.frames[0].aa_seq.clone();
    let mut pairs = Vec::new();
    let mut scan = |unit: usize, lead: usize| {
        let offset = plan.subject.units[unit].1.offset;
        let mut words = Vec::new();
        x_scan_unit(
            &plan,
            unit,
            &aa_seq[1 + offset - lead..],
            lead,
            &mut pairs,
            |filled| {
                words.extend(filled.iter().map(|p| (p.q_off, p.s_off)));
            },
        );
        words
    };
    // The whole frame's words, at the frame's offsets.
    let frame_words = scan(0, 0);
    // Those at frame offsets `offsets`, at the chunk's offsets (10 less, wrapping as Uint4).
    let at_chunk = |offsets: std::ops::RangeInclusive<u32>| -> Vec<(u32, u32)> {
        frame_words
            .iter()
            .filter(|(_, s)| offsets.contains(s))
            .map(|&(q, s)| (q, s.wrapping_sub(10)))
            .collect()
    };
    // (2, 6): the chunk's words at 2 and 3 (frame 12 and 13); (12, -4): none; (-10, 20): from 10
    // letters before the chunk (frame 0) to the chunk's last word (frame 27).
    let mut expected = at_chunk(12..=13);
    expected.extend(at_chunk(0..=27));
    assert!(expected.iter().any(|&(_, s)| (s as i32) < 0));
    assert_eq!(scan(1, 10), expected);
    // (19, 19), raised to the word at the frame's last letter: none (it holds the NULLB after the
    // frame); (19, -6): none; (-8, 20): frame 2 to 27.
    assert_eq!(scan(2, 8), at_chunk(2..=27));
}
