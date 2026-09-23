//! Internal TBLASTN protein lookup and translated-subject seed stage.
//! HSP construction and later stages are still unsupported.

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
    pub query_offset: u32,
    pub subject_offset: u32,
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
// This entry admits uppercase unambiguous DNA only. NCBI resolves ambiguities during
// ncbi2na encoding (blast_objmgr_tools.cpp:428-473) and preserves ncbi4na for
// later reevaluation (blast_setup_cxx.cpp:742-745,830-853).
#[allow(dead_code)]
pub(super) fn scan_unambiguous_blosum62_words(
    query: &[u8],
    subject: &[u8],
    db_gencode: u8,
    seg: Option<&SegParams>,
    threshold: i32,
) -> Result<Vec<Seed>> {
    if !subject
        .iter()
        .all(|base| matches!(*base, b'A' | b'C' | b'G' | b'T'))
    {
        bail!("TBLASTN ambiguous or lowercase subject encoding/masking is not implemented");
    }
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let query_frame = encode_protein_query_frame_with_seg(query, seg);
    let karlin = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let (lookup, _contexts) = build_ncbi_lookup(&[vec![query_frame]], threshold, &karlin, false);
    let mut seeds = Vec::new();
    let mut pairs = vec![BlastOffsetPair::default(); (lookup.longest_chain.max(1) as usize) * 1024];
    let pair_capacity = i32::try_from(pairs.len()).expect("NCBI offset array fits Int4");

    for frame in generate_frames(subject, &code) {
        let mut scan_range = [0, 0, frame.aa_len as i32 - lookup.word_length as i32];
        while scan_range[1] <= scan_range[2] {
            let hits = s_blast_aa_scan_subject_one_range(
                &lookup,
                &frame.aa_seq[1..],
                &mut pairs,
                pair_capacity,
                &mut scan_range,
            );
            for pair in pairs.iter().take(hits as usize) {
                seeds.push(Seed {
                    frame: frame.frame,
                    query_offset: pair.q_off,
                    subject_offset: pair.s_off,
                });
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
            let seeds = scan_unambiguous_blosum62_words(query, subject, 1, None, 13).unwrap();
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
                scan_unambiguous_blosum62_words(query, by_id[subject_id], 1, None, 13).unwrap();
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
        let unmasked = scan_unambiguous_blosum62_words(query, subject, 1, None, 13).unwrap();
        let masked =
            scan_unambiguous_blosum62_words(query, subject, 1, Some(&SegParams::default()), 13)
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

    // NCBI reference: c++/src/objects/seqfeat/gc.prt:340-347
    // id 32 , ncbieaa "FFLLSSSSYYWWCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    #[test]
    fn code_32_changes_search_words_and_ambiguity_fails_explicitly() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subject_code32.fna"))[0].1;
        let code1 = scan_unambiguous_blosum62_words(query, subject, 1, None, 13).unwrap();
        let code32 = scan_unambiguous_blosum62_words(query, subject, 32, None, 13).unwrap();
        println!(
            "fixture=code32 code1_seeds={} code32_seeds={}",
            code1.len(),
            code32.len()
        );
        assert_ne!(code1, code32);
        assert!(code32
            .iter()
            .any(|seed| seed.frame == 1 && seed.query_offset == 0 && seed.subject_offset == 0));
        assert!(scan_unambiguous_blosum62_words(query, b"ATN", 32, None, 13).is_err());
        assert!(scan_unambiguous_blosum62_words(query, b"atg", 32, None, 13).is_err());
    }
}
