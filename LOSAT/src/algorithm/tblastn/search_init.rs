//! Internal NCBI TBLASTN two-hit WordFinder path for the BLOSUM62/word-3 profile.
//! The public search remains unsupported; this is a bounded Stage C diagnostic.

use super::search_seed::resolve_local_subject_ncbi2na;
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::extension::extend_two_hit;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::config::ScoringMatrix;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;
use anyhow::{ensure, Result};

// NCBI reference: c++/include/algo/blast/core/blast_extend.h:142-163
// typedef struct BlastUngappedData { Int4 q_start, s_start, length, score; } ...;
// typedef struct BlastInitHSP { BlastOffsetPair offsets;
//                              BlastUngappedData* ungapped_data; } ...;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct InitHsp {
    pub frame: i8,
    // NCBI c++/src/algo/blast/core/blast_engine.c:241-250,572-586:
    // backup.offset is retained until HSP offsets are adjusted before merge.
    pub chunk_offset: u32,
    pub q_seed: u32,
    pub s_seed: u32,
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

// NCBI reference: c++/src/algo/blast/core/blast_extend.c:46-65,162-175
// while (diag_array_length < (qlen+window_size)) diag_array_length <<= 1;
// diag_table->offset = window_size;
// if (ewp->diag_table->offset >= INT4_MAX / 4) { ... clear ... }
// else ewp->diag_table->offset += subject_length + ewp->diag_table->window;
struct Diagonals {
    entries: Vec<(i32, bool)>,
    mask: i32,
    offset: i32,
    window: i32,
}

impl Diagonals {
    fn new(query_length: usize, window: i32) -> Self {
        let size = (query_length + window as usize).next_power_of_two();
        Self {
            entries: vec![(0, false); size],
            mask: size as i32 - 1,
            offset: window,
            window,
        }
    }

    fn finish_frame(&mut self, subject_length: usize) {
        if self.offset >= i32::MAX / 4 {
            self.offset = self.window;
            self.entries.fill((0, false));
        } else {
            self.offset += subject_length as i32 + self.window;
        }
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:804-841
// for (context=first_context; context<=last_context; context++) {
//     subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//     status = s_BlastSearchEngineOneContext(...);
// }
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:478-614
// scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
// hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// for (i = 0; i < hits; ++i) { /* two-hit diagonal state and extension */ }
// Blast_ExtendWordExit(ewp, subject->length);
#[allow(dead_code)]
pub(super) fn find_blosum62_word3_init_hsps(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    window: i32,
    x_dropoff: i32,
    cutoff_score: i32,
    mask_lowercase: bool,
) -> Result<Vec<InitHsp>> {
    find_blosum62_word3_init_hsps_multi(
        &[query],
        subject,
        db_gencode,
        seg,
        threshold,
        window,
        &[x_dropoff],
        &[cutoff_score],
        mask_lowercase,
    )
}

// NCBI c++/src/algo/blast/core/aa_ungapped.c:547-583:
// curr_context = BSearchContextInfo(query_offset, query_info);
// if (query_offset - diff <
//     query_info->contexts[curr_context].query_offset) { ... continue; }
// cutoffs = word_params->cutoffs + curr_context;
// s_BlastAaExtendTwoHit(..., query, ..., query_offset, ...);
// NCBI c++/src/algo/blast/core/blast_extend.c:46-65:
// while (diag_array_length < (qlen+window_size)) diag_array_length <<= 1;
// Keep one diagonal table for the concatenated protein queries.
#[allow(dead_code)]
pub(super) fn find_blosum62_word3_init_hsps_multi(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    window: i32,
    x_dropoffs: &[i32],
    cutoff_scores: &[i32],
    mask_lowercase: bool,
) -> Result<Vec<InitHsp>> {
    find_protein_init_hsps_multi(
        queries,
        subject,
        db_gencode,
        seg,
        threshold,
        window,
        x_dropoffs,
        cutoff_scores,
        mask_lowercase,
        ScoringMatrix::Blosum62,
        3,
    )
}

// NCBI c++/src/algo/blast/core/aa_ungapped.c:516-614:
// diff < wordsize controls the second hit; s_BlastAaExtendTwoHit uses
// wordsize and score_params->matrix before saving/sorting each HSP list.
#[allow(dead_code)]
pub(super) fn find_protein_init_hsps_multi(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    window: i32,
    x_dropoffs: &[i32],
    cutoff_scores: &[i32],
    mask_lowercase: bool,
    matrix: ScoringMatrix,
    word_size: usize,
) -> Result<Vec<InitHsp>> {
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:547-583:
    // cutoffs = word_params->cutoffs + curr_context;
    // Require one cutoff per protein query for the indexed read.
    ensure!(
        x_dropoffs.len() == queries.len() && cutoff_scores.len() == queries.len(),
        "one WordFinder x-drop and cutoff per query context"
    );
    let resolved = resolve_local_subject_ncbi2na(subject)?;
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let frames = generate_frames(&resolved, &code);
    // NCBI c++/src/algo/blast/core/blast_query_info.c:246-250:
    // last context's query_offset + query_length is the sequence length;
    // each preceding protein context ends with one NULLB separator.
    let mut query_sequence = Vec::new();
    let mut context_offsets = Vec::with_capacity(queries.len());
    for query in queries {
        context_offsets.push(i32::try_from(query_sequence.len())?);
        let frame = encode_protein_query_frame_with_seg(query, seg);
        query_sequence.extend_from_slice(&frame.aa_seq[1..]);
    }
    let seeds = super::search_seed::scan_unambiguous_protein_words_multi(
        queries,
        subject,
        db_gencode,
        seg,
        threshold,
        mask_lowercase,
        matrix,
        word_size,
    )?;
    let mut diagonals = Diagonals::new(query_sequence.len().saturating_sub(1), window);
    let mut hits = Vec::new();
    let mut seed_index = 0;
    let word_size_i32 = i32::try_from(word_size)?;
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);

    for frame in frames {
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-552:
        // WordFinder, GetGappedScore, and endpoint purge run per subject chunk.
        for (chunk, _) in super::search_seed::translated_chunk_scan_ranges(
            subject,
            frame.frame,
            frame.aa_len,
            negative_first_length,
            mask_lowercase,
        )? {
            let chunk_start = hits.len();
            let subject_sequence = &frame.aa_seq[1 + chunk.offset..1 + chunk.offset + chunk.length];
            while seed_index < seeds.len()
                && seeds[seed_index].frame == frame.frame
                && seeds[seed_index].chunk_offset == chunk.offset as u32
            {
                let seed = seeds[seed_index];
                seed_index += 1;
                let q = seed.query_offset as i32;
                let s = seed.subject_offset as i32;
                // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:516-547
                // diag_coord = (query_offset - subject_offset) & diag_mask;
                // if (diag_array[diag_coord].flag) {
                //     if (subject_offset + diag_offset < last_hit) continue;
                //     last_hit = subject_offset + diag_offset; flag = 0;
                // } else {
                //     last_hit = diag_array[diag_coord].last_hit - diag_offset;
                //     diff = subject_offset - last_hit;
                //     if (diff >= window) { last_hit = subject_offset + diag_offset; continue; }
                //     if (diff < wordsize) continue;
                // }
                let index = ((q - s) & diagonals.mask) as usize;
                let (last_hit, flag) = &mut diagonals.entries[index];
                if *flag {
                    if s + diagonals.offset < *last_hit {
                        continue;
                    }
                    *last_hit = s + diagonals.offset;
                    *flag = false;
                    continue;
                }
                let previous = *last_hit - diagonals.offset;
                let diff = s - previous;
                if diff >= window {
                    *last_hit = s + diagonals.offset;
                    continue;
                }
                if diff < word_size_i32 {
                    continue;
                }
                // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:550-588
                // if (query_offset - diff < query_info->contexts[curr_context].query_offset)
                //     { last_hit = subject_offset + diag_offset; continue; }
                // score = s_BlastAaExtendTwoHit(matrix, subject, query,
                //         last_hit + wordsize, subject_offset, query_offset,
                //         cutoffs->x_dropoff, ..., wordsize, &right_extend, &s_last_off);
                // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
                // NCBI c++/src/algo/blast/core/aa_ungapped.c:562-579:
                // curr_context = BSearchContextInfo(query_offset, query_info);
                // if (query_offset - diff <
                //     query_info->contexts[curr_context].query_offset) continue;
                let context = context_offsets.partition_point(|&offset| offset <= q) - 1;
                if q - diff < context_offsets[context] {
                    *last_hit = s + diagonals.offset;
                    continue;
                }
                // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-583,1089-1155:
                // s_BlastAaExtendTwoHit(matrix, ..., wordsize, ...);
                let result = extend_two_hit(
                    matrix,
                    &query_sequence,
                    subject_sequence,
                    (previous + word_size_i32) as usize,
                    s as usize,
                    q as usize,
                    x_dropoffs[context],
                    word_size,
                );
                let Some(result) = result else {
                    continue;
                };
                // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
                // cutoffs = word_params->cutoffs + curr_context;
                // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
                if result.ungapped_data.score >= cutoff_scores[context] {
                    hits.push(InitHsp {
                        frame: frame.frame,
                        chunk_offset: seed.chunk_offset,
                        q_seed: seed.query_offset,
                        s_seed: seed.subject_offset,
                        q_start: result.ungapped_data.q_start,
                        s_start: result.ungapped_data.s_start,
                        length: result.ungapped_data.length,
                        score: result.ungapped_data.score,
                    });
                }
                // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:588-607
                // if (right_extend) {
                //     diag_array[diag_coord].flag = 1;
                //     diag_array[diag_coord].last_hit =
                //         s_last_off - (wordsize - 1) + diag_offset;
                // } else { last_hit = subject_offset + diag_offset; }
                if result.right_extend {
                    *flag = true;
                    *last_hit = result.s_last_off - (word_size_i32 - 1) + diagonals.offset;
                } else {
                    *last_hit = s + diagonals.offset;
                }
            }
            // NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234:
            // status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
            // Blast_InitHitListSortByScore(init_hitlist);
            // NCBI c++/src/algo/blast/core/blast_extend.c:273-313:
            // compare score DESC, subject start ASC, length DESC, query start ASC.
            hits[chunk_start..].sort_unstable_by(|a, b| {
                b.score
                    .cmp(&a.score)
                    .then(a.s_start.cmp(&b.s_start))
                    .then(b.length.cmp(&a.length))
                    .then(a.q_start.cmp(&b.q_start))
            });
            // NCBI c++/src/algo/blast/core/aa_ungapped.c:609-614:
            // Blast_ExtendWordExit(ewp, subject->length);
            diagonals.finish_frame(chunk.length);
        }
    }
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    // NCBI c++/src/algo/blast/core/blast_engine.c:259-310,478-487:
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

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:516-614:
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
        let actual = find_protein_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            16,
            60,
            &[21, 21, 0],
            &[0, 0, i32::MAX],
            false,
            ScoringMatrix::Blosum45,
            2,
        )
        .unwrap();
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

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:562-583:
    // BSearchContextInfo finds the query context before two-hit extension;
    // the left word must remain within that context's query_offset.
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234:
    // Blast_InitHitListSortByScore(init_hitlist) fixes frame-local HSP order.
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
        let actual = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[0, 0, i32::MAX],
            false,
        )
        .unwrap();
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

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // cutoffs = word_params->cutoffs + curr_context;
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
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
        let equal = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[20, 393, i32::MAX],
            false,
        )
        .unwrap();
        assert!(equal.iter().any(|h| h.score == 20 && h.q_seed < 121));
        assert!(equal
            .iter()
            .any(|h| h.score == 393 && (121..192).contains(&h.q_seed)));
        let above = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 16, 0],
            &[21, 394, i32::MAX],
            false,
        )
        .unwrap();
        assert!(!above.iter().any(|h| h.score == 20 && h.q_seed < 121));
        assert!(!above
            .iter()
            .any(|h| h.score == 393 && (121..192).contains(&h.q_seed)));
    }

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:562-590:
    // cutoffs = word_params->cutoffs + curr_context;
    // s_BlastAaExtendTwoHit(..., cutoffs->x_dropoff, ...);
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
        let actual = find_blosum62_word3_init_hsps_multi(
            &refs,
            subject,
            1,
            None,
            13,
            40,
            &[16, 1, 0],
            &[0, 0, i32::MAX],
            false,
        )
        .unwrap();
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

    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:816-831
    // subject->seq_ranges[i].left = CONV_NUCL2PROT_COORDINATES(...);
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:575-614
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
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
            for hit in
                find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, true).unwrap()
            {
                actual.push((id.clone(), hit));
            }
        }
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:575-614
    // score = s_BlastAaExtendTwoHit(...);
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
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
        let actual =
            find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, false).unwrap();
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:200-234,575-614
    // status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    // Blast_InitHitListSortByScore(init_hitlist);
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:747-775
    // BLAST_GetAllTranslations(..., subject->gen_code_string, ...);
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
        let actual =
            find_blosum62_word3_init_hsps(query, subject, 32, None, 13, 40, 16, 13, false).unwrap();
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:200-234,575-614
    // score = s_BlastAaExtendTwoHit(..., cutoffs->x_dropoff, ...);
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    // Blast_InitHitListSortByScore(init_hitlist);
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
            for hit in find_blosum62_word3_init_hsps(query, subject, 1, None, 13, 40, 16, 0, false)
                .unwrap()
            {
                actual.push((id.clone(), hit));
            }
        }
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:478-552:
    // s_GetNextSubjectChunk precedes each WordFinder call and the hit list
    // is reset before the next chunk. Subject offsets stay chunk-local.
    #[test]
    fn long_subject_chunk_wordfinder_hsps_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let query = &read_fasta(&format!("{root}long_chunk_20260924/query.faa"))[0].1;
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
        let trace =
            fs::read_to_string(format!("{root}long_chunk_20260924/frame_chunks.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                InitHsp {
                    frame: [1, 2, 3, -1, -2, -3][call / 2],
                    chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 28.
        let actual =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 28, false).unwrap();
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:283-310:
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
        let resolved = resolve_local_subject_ncbi2na(&subject).unwrap();
        let code = GeneticCode::try_from_id(1).unwrap();
        let frames = generate_frames(&resolved, &code);
        let negative_first_length = frames.iter().find(|f| f.frame == -1).unwrap().aa_len;
        let actual_ranges: Vec<_> = frames
            .iter()
            .flat_map(|frame| {
                super::super::search_seed::translated_chunk_scan_ranges(
                    &subject,
                    frame.frame,
                    frame.aa_len,
                    negative_first_length,
                    true,
                )
                .unwrap()
                .into_iter()
                .map(|(chunk, ranges)| (frame.frame, chunk.length, ranges))
            })
            .collect();
        assert_eq!(
            actual_ranges,
            ncbi_chunk_range_rows(&format!(
                "{root}masked_chunk_boundary_20260924/chunk_ranges.tsv"
            ))
        );
        drop(frames);
        drop(resolved);
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
                    chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 28.
        let actual =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 28, true).unwrap();
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:283-310,478-500:
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
        let resolved = resolve_local_subject_ncbi2na(&subject).unwrap();
        let code = GeneticCode::try_from_id(1).unwrap();
        let frames = generate_frames(&resolved, &code);
        let negative_first_length = frames.iter().find(|f| f.frame == -1).unwrap().aa_len;
        let actual_ranges: Vec<_> = frames
            .iter()
            .flat_map(|frame| {
                super::super::search_seed::translated_chunk_scan_ranges(
                    &subject,
                    frame.frame,
                    frame.aa_len,
                    negative_first_length,
                    true,
                )
                .unwrap()
                .into_iter()
                .map(|(chunk, ranges)| (frame.frame, chunk.length, ranges))
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
        drop(frames);
        drop(resolved);
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
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 30.
        let actual =
            find_blosum62_word3_init_hsps(query, &subject, 1, None, 13, 40, 16, 30, true).unwrap();
        assert_eq!(actual, expected);
    }
}
