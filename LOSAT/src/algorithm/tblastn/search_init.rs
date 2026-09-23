//! Internal NCBI TBLASTN two-hit WordFinder path for the BLOSUM62/word-3 profile.
//! Later gapped construction, chunk merging and traceback remain unimplemented.

use super::search_seed::{resolve_local_subject_ncbi2na, scan_unambiguous_blosum62_words};
use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::blastp::extension::extend_two_hit_blosum62;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;
use anyhow::Result;

// NCBI reference: c++/include/algo/blast/core/blast_extend.h:142-163
// typedef struct BlastUngappedData { Int4 q_start, s_start, length, score; } ...;
// typedef struct BlastInitHSP { BlastOffsetPair offsets;
//                              BlastUngappedData* ungapped_data; } ...;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct InitHsp {
    pub frame: i8,
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
    let resolved = resolve_local_subject_ncbi2na(subject)?;
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let frames = generate_frames(&resolved, &code);
    let query_frame = encode_protein_query_frame_with_seg(query, seg);
    let query_sequence = &query_frame.aa_seq[1..query_frame.aa_seq.len() - 1];
    let seeds = scan_unambiguous_blosum62_words(
        query,
        subject,
        db_gencode,
        seg,
        threshold,
        mask_lowercase,
    )?;
    let mut diagonals = Diagonals::new(query.len(), window);
    let mut hits = Vec::new();
    let mut seed_index = 0;
    const WORD_SIZE: i32 = 3;

    for frame in frames {
        let subject_sequence = &frame.aa_seq[1..frame.aa_seq.len() - 1];
        while seed_index < seeds.len() && seeds[seed_index].frame == frame.frame {
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
            if diff < WORD_SIZE {
                continue;
            }
            // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:550-588
            // if (query_offset - diff < query_info->contexts[curr_context].query_offset)
            //     { last_hit = subject_offset + diag_offset; continue; }
            // score = s_BlastAaExtendTwoHit(matrix, subject, query,
            //         last_hit + wordsize, subject_offset, query_offset,
            //         cutoffs->x_dropoff, ..., wordsize, &right_extend, &s_last_off);
            // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
            if q - diff < 0 {
                *last_hit = s + diagonals.offset;
                continue;
            }
            let result = extend_two_hit_blosum62(
                query_sequence,
                subject_sequence,
                (previous + WORD_SIZE) as usize,
                s as usize,
                q as usize,
                x_dropoff,
                WORD_SIZE as usize,
            );
            let Some(result) = result else {
                continue;
            };
            if result.ungapped_data.score >= cutoff_score {
                hits.push(InitHsp {
                    frame: frame.frame,
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
                *last_hit = result.s_last_off - (WORD_SIZE - 1) + diagonals.offset;
            } else {
                *last_hit = s + diagonals.offset;
            }
        }
        diagonals.finish_frame(frame.aa_len);
    }
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

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
}
