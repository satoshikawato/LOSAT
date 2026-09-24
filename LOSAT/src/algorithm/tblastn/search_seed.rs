//! Internal TBLASTN protein lookup and translated-subject seed stage.
//! The public search remains unsupported; this is a bounded Stage C diagnostic.

use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
use crate::algorithm::tblastx::blast_aascan::{s_blast_aa_scan_subject_one_range, BlastOffsetPair};
use crate::algorithm::tblastx::lookup::build_ncbi_lookup;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::config::ScoringMatrix;
use crate::stats::lookup_protein_params_ungapped;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;
use anyhow::{bail, Result};

// NCBI reference: c++/include/algo/blast/core/blast_def.h:141-150
// typedef union BlastOffsetPair {
//     struct { Uint4 q_off; Uint4 s_off; } qs_offsets;
// } BlastOffsetPair;
// NCBI c++/src/algo/blast/core/blast_engine.c:808-812:
// subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Seed {
    pub frame: i8,
    // NCBI c++/src/algo/blast/core/blast_engine.c:241-250:
    // backup.offset identifies this chunk before WordFinder returns local offsets.
    pub chunk_offset: u32,
    pub query_offset: u32,
    pub subject_offset: u32,
}

// NCBI c++/src/algo/blast/core/blast_engine.c:246-264:
// if (backup->offset + MAX_DBSEQ_LEN < hard_ranges[hm_index].right) {
//     subject->length = MAX_DBSEQ_LEN;
//     backup->next = backup->offset + MAX_DBSEQ_LEN - dbseq_chunk_overlap;
// } else { subject->length = hard_ranges[hm_index].right - backup->offset; }
// NCBI c++/include/algo/blast/core/blast_gapalign.h:54:
// #define MAX_DBSEQ_LEN 5000000
// NCBI c++/include/algo/blast/core/blast_hits.h:192:
// #define DBSEQ_CHUNK_OVERLAP 100
// This schedule covers the unmasked translated frame's single hard range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct TranslatedChunk {
    pub offset: usize,
    pub length: usize,
}

// NCBI c++/src/algo/blast/core/blast_engine.c:246-264:
// backup->offset = backup->next;
// while (backup->next < backup->full_range.right) { s_GetNextSubjectChunk(...); }
#[allow(dead_code)]
pub(super) fn unmasked_translated_chunks(frame_length: usize) -> Vec<TranslatedChunk> {
    const MAX_DBSEQ_LEN: usize = 5_000_000;
    const DBSEQ_CHUNK_OVERLAP: usize = 100;
    let mut chunks = Vec::new();
    let mut offset = 0;
    while offset < frame_length {
        if offset + MAX_DBSEQ_LEN < frame_length {
            chunks.push(TranslatedChunk {
                offset,
                length: MAX_DBSEQ_LEN,
            });
            offset += MAX_DBSEQ_LEN - DBSEQ_CHUNK_OVERLAP;
        } else {
            chunks.push(TranslatedChunk {
                offset,
                length: frame_length - offset,
            });
            break;
        }
    }
    chunks
}

// NCBI reference: c++/src/util/random_gen.cpp:98,227-230,287-308 and
// c++/include/util/random_gen.hpp:224-241
// static const size_t kStateOffset = 12;
// m_State[0] = m_Seed = seed;
// for (int i = 1; i < kStateSize; ++i)
//     m_State[i] = 1103515245 * m_State[i-1] + 12345;
// m_RJ = kStateOffset; m_RK = kStateSize - 1;
// for (int i = 0; i < 10 * kStateSize; ++i) GetRand();
// r = m_State[m_RK] + m_State[m_RJ--];
// m_State[m_RK--] = r;
// return r >> 1;
struct NcbiRandom {
    state: [u32; 33],
    j: usize,
    k: usize,
}

impl NcbiRandom {
    fn new(seed: u32) -> Self {
        let mut state = [0; 33];
        state[0] = seed;
        for i in 1..state.len() {
            state[i] = state[i - 1]
                .wrapping_mul(1_103_515_245)
                .wrapping_add(12_345);
        }
        let mut random = Self {
            state,
            j: 12,
            k: 32,
        };
        for _ in 0..330 {
            random.get_rand();
        }
        random
    }

    fn get_rand(&mut self) -> u32 {
        let value = self.state[self.k].wrapping_add(self.state[self.j]);
        self.state[self.k] = value;
        self.j = if self.j == 0 { 32 } else { self.j - 1 };
        self.k = if self.k == 0 { 32 } else { self.k - 1 };
        value >> 1
    }
}

// NCBI reference: c++/src/algo/blast/core/blast_encoding.c:95-103 and
// c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474,515-520
// static unsigned char ctable[16] = {0xFF,0,1,0xFF,2,0xFF,0xFF,0xFF,
//                                    3,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF};
// CRandom random(base_length);
// if (b == 0 || b == 0x0F) ncbi2na[i] = random.GetRand() & 0x3;
// else { int pick = random.GetRand() % bitcount; /* pick a set bit */ }
// The uncompressed ncbi4na sequence is retained by SetupSubjects_OMF for
// traceback reevaluation; these resolved bases are preliminary-search input.
pub(super) fn resolve_local_subject_ncbi2na(subject: &[u8]) -> Result<Vec<u8>> {
    let seed = u32::try_from(subject.len())?;
    let mut random = NcbiRandom::new(seed);
    let mut resolved = Vec::with_capacity(subject.len());
    for &base in subject {
        let mask: u8 = match base.to_ascii_uppercase() {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            b'M' => 3,
            b'R' => 5,
            b'S' => 6,
            b'V' => 7,
            b'W' => 9,
            b'Y' => 10,
            b'H' => 11,
            b'K' => 12,
            b'D' => 13,
            b'B' => 14,
            b'N' => 15,
            b'-' => 0,
            _ => bail!("TBLASTN subject base '{}' is not implemented", base as char),
        };
        let code = match mask {
            1 => 0,
            2 => 1,
            4 => 2,
            8 => 3,
            0 | 15 => random.get_rand() % 4,
            _ => {
                let mut pick = random.get_rand() % mask.count_ones();
                let mut code = 0;
                for i in 0..4 {
                    if mask & (1 << i) != 0 {
                        if pick == 0 {
                            code = i;
                            break;
                        }
                        pick -= 1;
                    }
                }
                code
            }
        };
        resolved.push(b"ACGT"[code as usize]);
    }
    Ok(resolved)
}

// NCBI reference: c++/src/algo/blast/api/blast_setup_cxx.cpp:813-826
// BlastSeqBlkSetSeqRanges(subj, (SSeqRange*) masked_ranges.get_data(),
//                          masked_ranges.size() + 1, true, eSoftSubjMasking);
// NCBI c++/src/algo/blast/blastinput/blast_fasta_input.cpp:489-502:
// apply_mask_to_both_strands = true; PackedSeqLocToMaskedQueryRegions(...);
// NCBI c++/src/algo/blast/api/blast_aux_priv.cpp:436-454:
// if (assume_both_strands) do_pos = do_neg = true;
// if (do_pos) mqr.push_back(...); if (do_neg) mqr.push_back(...);
// NCBI c++/src/algo/blast/api/blast_setup_cxx.cpp:689-707,825-826:
// both mask copies enter masked_ranges before BlastSeqBlkSetSeqRanges.
// NCBI c++/src/algo/blast/core/blast_util.c:211-215:
// tmp[0].left = 0;
// tmp[num_seq_ranges - 1].right = seq_blk->length;
// NCBI reference: c++/src/algo/blast/core/blast_engine.c:816-831
// if (context == 0) { seq_ranges[i].left = backup.seq_ranges[i].left / 3;
//                     seq_ranges[i].right = backup.seq_ranges[i].right / 3; }
// else if (context == 3) { /* reverse ranges using subject->length */ }
// NCBI reference: c++/src/algo/blast/core/masksubj.inl:51-71
// range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
// range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
fn translated_scan_ranges(
    subject: &[u8],
    frame: i8,
    frame_length: usize,
    negative_first_length: usize,
    mask_lowercase: bool,
) -> Vec<(i32, i32)> {
    if !mask_lowercase || !subject.iter().any(u8::is_ascii_lowercase) {
        return vec![(0, frame_length as i32)];
    }
    let mut nucleotide_ranges = Vec::new();
    let mut left = 0;
    let mut index = 0;
    while index < subject.len() {
        if subject[index].is_ascii_lowercase() {
            let start = index;
            while index < subject.len() && subject[index].is_ascii_lowercase() {
                index += 1;
            }
            // NCBI c++/src/algo/blast/api/blast_aux_priv.cpp:436-454:
            // if (assume_both_strands) do_pos = do_neg = true;
            // if (do_pos) mqr.push_back(...); if (do_neg) mqr.push_back(...);
            for _ in 0..2 {
                nucleotide_ranges.push((left, start));
                left = index - 1;
            }
        } else {
            index += 1;
        }
    }
    nucleotide_ranges.push((left, subject.len()));
    if frame > 0 {
        nucleotide_ranges
            .into_iter()
            .map(|(left, right)| ((left / 3) as i32, (right / 3) as i32))
            .collect()
    } else {
        let length = negative_first_length as i32;
        nucleotide_ranges
            .into_iter()
            .rev()
            .map(|(left, right)| (length - (right / 3) as i32, length - (left / 3) as i32))
            .collect()
    }
}

// NCBI c++/src/algo/blast/core/blast_engine.c:259-310:
// unsplit subjects keep backup.soft_ranges; split soft-masked subjects
// advance sm_index while right < chunk start, include ranges with left < end,
// clamp only the first left and last right, and skip SUBJECT_SPLIT_NO_RANGE.
// NCBI c++/src/algo/blast/core/blast_engine.c:261-274:
// without soft masking, each chunk has one range [0, subject->length].
pub(super) fn translated_chunk_scan_ranges(
    subject: &[u8],
    frame: i8,
    frame_length: usize,
    negative_first_length: usize,
    mask_lowercase: bool,
) -> Result<Vec<(TranslatedChunk, Vec<(i32, i32)>)>> {
    let frame_ranges = translated_scan_ranges(
        subject,
        frame,
        frame_length,
        negative_first_length,
        mask_lowercase,
    );
    let soft_masked = mask_lowercase && subject.iter().any(u8::is_ascii_lowercase);
    let mut sm_index = 0usize;
    let mut calls = Vec::new();
    for chunk in unmasked_translated_chunks(frame_length) {
        if !soft_masked {
            calls.push((chunk, vec![(0, i32::try_from(chunk.length)?)]));
            continue;
        }
        if chunk.offset == 0 && chunk.length == frame_length {
            calls.push((chunk, frame_ranges.clone()));
            continue;
        }
        let start = i32::try_from(chunk.offset)?;
        let end = i32::try_from(chunk.offset + chunk.length)?;
        let mut i = sm_index;
        while i < frame_ranges.len() && frame_ranges[i].1 < start {
            i += 1;
        }
        if i == frame_ranges.len() {
            bail!("NCBI translated soft-range sentinel is missing");
        }
        let first = i;
        while i < frame_ranges.len() && frame_ranges[i].0 < end {
            i += 1;
        }
        if i == 0 {
            bail!("NCBI translated soft-range index underflow");
        }
        sm_index = i - 1;
        if i == first {
            continue;
        }
        let mut ranges: Vec<_> = frame_ranges[first..i]
            .iter()
            .map(|&(left, right)| (left - start, right - start))
            .collect();
        ranges[0].0 = ranges[0].0.max(0);
        let last = ranges.len() - 1;
        ranges[last].1 = ranges[last].1.min(i32::try_from(chunk.length)?);
        calls.push((chunk, ranges));
    }
    Ok(calls)
}

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:747-775,804-812,835-844
// BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na,
//                          backup.full_range.right, subject->gen_code_string,
//                          &translation_buffer, &frame_offsets, NULL);
// for (context=first_context; context<=last_context; context++) {
//     subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
//     subject->sequence = translation_buffer + frame_offsets[context] + 1;
//     subject->length = frame_offsets[context+1] - frame_offsets[context] - 1;
//     status = s_BlastSearchEngineOneContext(...);
//     Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
// }
// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1023-1050
// BlastChooseProteinScanSubject(lookup_wrap);
// aux_struct->WordFinder = BlastAaWordFinder;
// NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
// scan_range[1] = subject->seq_ranges[0].left;
// scan_range[2] = subject->seq_ranges[0].right - wordsize;
// while (scan_range[1] <= scan_range[2]) {
//     hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// }
// This entry admits IUPAC DNA. NCBI resolves ambiguities during
// ncbi2na encoding (blast_objmgr_tools.cpp:428-473) and preserves ncbi4na for
// later reevaluation (blast_setup_cxx.cpp:742-745,830-853). The optional
// lowercase scan ranges follow the local -lcase_masking path above.
#[allow(dead_code)]
pub(super) fn scan_unambiguous_blosum62_words(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    mask_lowercase: bool,
) -> Result<Vec<Seed>> {
    scan_unambiguous_blosum62_words_multi(
        &[query],
        subject,
        db_gencode,
        seg,
        threshold,
        mask_lowercase,
    )
}

// NCBI c++/src/algo/blast/core/blast_query_info.c:68-96:
// retval->last_context = retval->num_queries * kNumContexts - 1;
// contexts[i].query_index = Blast_GetQueryIndexFromContext(i, program);
// NCBI c++/src/algo/blast/core/aa_ungapped.c:478-505:
// scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
// The lookup receives all protein contexts before each subject frame scan;
// scan order must retain the global query offsets from that shared lookup.
#[allow(dead_code)]
pub(super) fn scan_unambiguous_blosum62_words_multi(
    queries: &[&[u8]],
    subject: &[u8],
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
    mask_lowercase: bool,
) -> Result<Vec<Seed>> {
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:1318-1334,1429-1432
    // return word_length * 3 + 2;
    // if (seq_arg.seq->length < min_subj_seq_length) { ... continue; }
    if subject.len() < 3 * 3 + 2 {
        return Ok(Vec::new());
    }
    let resolved_subject = resolve_local_subject_ncbi2na(subject)?;
    let query_frames: Vec<_> = queries
        .iter()
        .map(|query| vec![encode_protein_query_frame_with_seg(query, seg)])
        .collect();
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (lookup, _contexts) = build_ncbi_lookup(&query_frames, threshold, &karlin, false);
    let mut seeds = Vec::new();
    let mut pairs = vec![BlastOffsetPair::default(); (lookup.longest_chain.max(1) as usize) * 1024];
    let pair_capacity = i32::try_from(pairs.len()).expect("NCBI offset array fits Int4");

    let frames = generate_frames(&resolved_subject, &code);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);
    for frame in frames {
        // NCBI c++/src/algo/blast/core/blast_engine.c:478-500:
        // SUBJECT_SPLIT_NO_RANGE skips WordFinder for this chunk.
        for (chunk, ranges) in translated_chunk_scan_ranges(
            subject,
            frame.frame,
            frame.aa_len,
            negative_first_length,
            mask_lowercase,
        )? {
            for (range_index, (left, right)) in ranges.into_iter().enumerate() {
                // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:496-501
                // scan_range[1] = subject->seq_ranges[0].left;
                // scan_range[2] = subject->seq_ranges[0].right - wordsize;
                // if (scan_range[2] < scan_range[1])
                //     scan_range[2] = scan_range[1];
                let mut end = right - lookup.word_length as i32;
                if range_index == 0 && end < left {
                    end = left;
                }
                let mut scan_range = [0, left, end];
                while scan_range[1] <= scan_range[2] {
                    let hits = s_blast_aa_scan_subject_one_range(
                        &lookup,
                        &frame.aa_seq[1 + chunk.offset..],
                        &mut pairs,
                        pair_capacity,
                        &mut scan_range,
                    );
                    for pair in pairs.iter().take(hits as usize) {
                        seeds.push(Seed {
                            frame: frame.frame,
                            chunk_offset: u32::try_from(chunk.offset)?,
                            query_offset: pair.q_off,
                            subject_offset: pair.s_off,
                        });
                    }
                }
            }
        }
    }
    Ok(seeds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::encoding::ncbistdaa_to_ascii;
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

    // NCBI c++/src/algo/blast/core/blast_query_info.c:68-96:
    // one context per protein query; query offsets advance by length + 1.
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:478-505:
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
        let actual =
            scan_unambiguous_blosum62_words_multi(&refs, subject, 1, None, 13, false).unwrap();
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

    // NCBI reference: c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
    // CRandom random(base_length); ambiguous ncbi4na bases are picked once
    // for the compressed ncbi2na preliminary subject.
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:767-775,804-812
    // BLAST_GetAllTranslations(backup.sequence, eBlastEncodingNcbi2na, ...);
    // subject->sequence = translation_buffer + frame_offsets[context] + 1;
    #[test]
    fn ncbi_preliminary_subject_frames_match_all_fixture_bytes() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let records = read_fasta(&format!("{root}subjects.fna"));
        let by_id: HashMap<_, _> = records.iter().map(|(id, dna)| (id.as_str(), dna)).collect();
        let trace = fs::read_to_string(format!("{root}ncbi_frame_trace.tsv")).unwrap();
        let mut compared = 0;
        for line in trace.lines().skip(1) {
            let columns: Vec<_> = line.split('\t').collect();
            let subject_id = columns[1];
            let frame: i8 = columns[2].parse().unwrap();
            let resolved = resolve_local_subject_ncbi2na(by_id[subject_id]).unwrap();
            let translated = generate_frames(&resolved, &GeneticCode::from_id(1));
            let actual = translated.iter().find(|f| f.frame == frame).unwrap();
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

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:478-505
    // scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    // for (i = 0; i < hits; ++i) { query_offset = offset_pairs[i].qs_offsets.q_off;
    //                            subject_offset = offset_pairs[i].qs_offsets.s_off; }
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
            let seeds =
                scan_unambiguous_blosum62_words(query, subject, 1, None, 13, false).unwrap();
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

    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:816-831
    // if (context == 0) { /* convert subject scan ranges to AA */ }
    // else if (context == 3) { /* reverse ranges for minus frames */ }
    // NCBI reference: c++/src/algo/blast/core/masksubj.inl:51-71
    // range[1] = subject->seq_ranges[range[0]].left + word_length - lut_word_length;
    // range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
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
            for seed in scan_unambiguous_blosum62_words(query, subject, 1, None, 13, true).unwrap()
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

    // NCBI reference: c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474
    // CRandom random(base_length); each ambiguous 4na value consumes one GetRand().
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:492-505
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
    #[test]
    fn ncbi_ambiguity_spectrum_frames_and_candidates_match() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/ambiguity_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let resolved = resolve_local_subject_ncbi2na(subject).unwrap();
        let actual_frames = generate_frames(&resolved, &GeneticCode::from_id(1));
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
        let actual: Vec<_> = scan_unambiguous_blosum62_words(query, subject, 1, None, 13, false)
            .unwrap()
            .into_iter()
            .map(|x| (x.frame, x.query_offset, x.subject_offset))
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 181);
    }

    // NCBI reference: c++/src/algo/blast/core/blast_util.c:1080-1101
    // frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    // offset += length + 1;
    // frame_offsets[context+1] = offset;
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
            let frame = generate_frames(subject, &GeneticCode::from_id(1))
                .into_iter()
                .find(|f| f.frame == expected_frame)
                .unwrap();
            let translated: String = frame.aa_seq[1..frame.aa_seq.len() - 1]
                .iter()
                .map(|&aa| ncbistdaa_to_ascii(aa))
                .collect();
            assert!(translated.contains(fields[9]), "{subject_id}: {translated}");
            let seeds =
                scan_unambiguous_blosum62_words(query, subject, 1, None, 13, false).unwrap();
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

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:200-234
    // status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
    // Blast_InitHitListSortByScore(init_hitlist);
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
            let seeds =
                scan_unambiguous_blosum62_words(query, by_id[subject_id], 1, None, 13, false)
                    .unwrap();
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

    // NCBI reference: c++/src/algo/blast/core/blast_filter.c:337-370
    // BlastSetUp_Filter(..., query_blk, ...);
    // query_blk->sequence_start_nomask = BlastMemDup(query_blk->sequence_start, total_length);
    #[test]
    fn seg_suppresses_low_complexity_query_seeds() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let query = &read_fasta(&format!("{root}query_low.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}low_subject.fna"))[0].1;
        let unmasked = scan_unambiguous_blosum62_words(query, subject, 1, None, 13, false).unwrap();
        let masked = scan_unambiguous_blosum62_words(
            query,
            subject,
            1,
            Some(&SegParams::default()),
            13,
            false,
        )
        .unwrap();
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

    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:747-775,804-841
    // BLAST_GetAllTranslations(..., subject->gen_code_string, ...);
    // subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:478-505
    // hits = scansub(lookup_wrap, subject, offset_pairs, array_size, scan_range);
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
        let actual: Vec<_> = scan_unambiguous_blosum62_words(query, subject, 32, None, 13, false)
            .unwrap()
            .into_iter()
            .map(|x| (x.frame, x.query_offset, x.subject_offset))
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(actual.len(), 199);
    }

    // NCBI reference: c++/src/objects/seqfeat/gc.prt:340-347
    // id 32 , ncbieaa "FFLLSSSSYYWWCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    #[test]
    fn code_32_changes_search_words_and_ambiguity_resolves() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subject_code32.fna"))[0].1;
        let code1 = scan_unambiguous_blosum62_words(query, subject, 1, None, 13, false).unwrap();
        let code32 = scan_unambiguous_blosum62_words(query, subject, 32, None, 13, false).unwrap();
        println!(
            "fixture=code32 code1_seeds={} code32_seeds={}",
            code1.len(),
            code32.len()
        );
        assert_ne!(code1, code32);
        assert!(code32
            .iter()
            .any(|seed| seed.frame == 1 && seed.query_offset == 0 && seed.subject_offset == 0));
        assert!(scan_unambiguous_blosum62_words(query, b"ATN", 32, None, 13, false).is_ok());
        assert!(scan_unambiguous_blosum62_words(query, b"atg", 32, None, 13, false).is_ok());
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:246-264:
    // split only when offset + MAX_DBSEQ_LEN < hard-range right;
    // the next offset overlaps the first chunk by DBSEQ_CHUNK_OVERLAP.
    #[test]
    fn unmasked_long_translated_chunks_match_ncbi_wordfinder_calls() {
        let rows = unmasked_translated_chunks(5_000_320);
        assert_eq!(
            rows,
            vec![
                TranslatedChunk {
                    offset: 0,
                    length: 5_000_000
                },
                TranslatedChunk {
                    offset: 4_999_900,
                    length: 420
                },
            ]
        );
        assert_eq!(
            unmasked_translated_chunks(5_000_000),
            vec![TranslatedChunk {
                offset: 0,
                length: 5_000_000
            },]
        );
        assert_eq!(
            unmasked_translated_chunks(5_000_001),
            vec![
                TranslatedChunk {
                    offset: 0,
                    length: 5_000_000
                },
                TranslatedChunk {
                    offset: 4_999_900,
                    length: 101
                },
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
}
