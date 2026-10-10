//! One concatenated six-context BLASTX lookup and diagonal state.
use super::{
    args::ResolvedOptions, parameters::ContextParameters, query_setup::PreparedQueryBatch,
};
use crate::algorithm::blastp::extension::{extend_one_hit_blosum62, extend_two_hit_blosum62};
use crate::algorithm::tblastx::blast_aascan::s_blast_aa_scan_subject_one_range;
use crate::algorithm::tblastx::blast_aascan::BlastOffsetPair as OffsetPair;
use crate::algorithm::tblastx::lookup::compressed::{
    build_blosum62_compressed_lookup, BlastCompressedAaLookupTable,
};
use crate::algorithm::tblastx::lookup::{
    build_ncbi_lookup_from_prepared, BlastAaLookupTable, QueryContext,
};
use anyhow::{ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
// ```c
//     aux_struct->offset_pairs =
//       (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
// ```
// EXPERIMENT (LOSAT_X_BXSCAN): the seed loop on a per-thread offset-pair buffer, read in place
// (a child module, so that it uses the diagonal table's fields).
#[path = "x_seed_scan.rs"]
pub mod x_seed_scan;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
// ```c++
//         Int4 diag_array_length;
//
//         diag_table= (BLAST_DiagTable*) calloc(1, sizeof(BLAST_DiagTable));
//
//         if (diag_table)
//         {
//                 diag_array_length = 1;
//                 /* What power of 2 is just longer than the query? */
//                 while (diag_array_length < (qlen+window_size))
//                 {
//                         diag_array_length = diag_array_length << 1;
//                 }
//                 /* These are used in the word finders to shift and mask
//                 rather than dividing and taking the remainder. */
//                 diag_table->diag_array_length = diag_array_length;
//                 diag_table->diag_mask = diag_array_length-1;
//                 diag_table->multiple_hits = multiple_hits;
//                 diag_table->offset = window_size;
//                 diag_table->window = window_size;
//         }
// ```
pub struct Diagonals {
    entries: Vec<(i32, bool)>,
    offset: i32,
    window: i32,
}
impl Diagonals {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    // ```c++
    //         Int4 diag_array_length;
    //
    //         diag_table= (BLAST_DiagTable*) calloc(1, sizeof(BLAST_DiagTable));
    //
    //         if (diag_table)
    //         {
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    //                 diag_table->diag_array_length = diag_array_length;
    //                 diag_table->diag_mask = diag_array_length-1;
    //                 diag_table->multiple_hits = multiple_hits;
    //                 diag_table->offset = window_size;
    //                 diag_table->window = window_size;
    //         }
    // ```
    pub fn new(query_length: usize, window: i32) -> Result<Self> {
        let n = query_length
            .checked_add(usize::try_from(window)?)
            .and_then(|n| n.max(1).checked_next_power_of_two())
            .ok_or_else(|| anyhow::anyhow!("BLASTX diagonal length exceeds Int4"))?;
        ensure!(
            n <= i32::MAX as usize,
            "BLASTX diagonal length exceeds Int4"
        );
        Ok(Self {
            entries: vec![(0, false); n],
            offset: window,
            window,
        })
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-175
    // ```c++
    // Blast_ExtendWordExit(Blast_ExtendWord * ewp, Int4 subject_length)
    // {
    //     if (!ewp)
    //         return -1;
    //
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    //         if (ewp->hash_table->offset >= INT4_MAX / 4) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
    // ```c++
    // static Int4 s_BlastDiagClear(BLAST_DiagTable * diag)
    // {
    //     Int4 i, n;
    //     DiagStruct *diag_struct_array;
    //
    //     if (diag == NULL)
    //         return 0;
    //
    //     n = diag->diag_array_length;
    //
    //     diag->offset = diag->window;
    //
    //     diag_struct_array = diag->hit_level_array;
    //
    //     for (i = 0; i < n; i++) {
    //         diag_struct_array[i].flag = 0;
    //         diag_struct_array[i].last_hit = -diag->window;
    //         if (diag->hit_len_array) diag->hit_len_array[i] = 0;
    //     }
    //     return 0;
    // }
    // ```
    pub fn finish_subject(&mut self, length: usize) -> Result<()> {
        if self.offset >= i32::MAX / 4 {
            self.offset = self.window;
            self.entries.fill((-self.window, false));
        } else {
            self.offset = self
                .offset
                .checked_add(i32::try_from(length)?)
                .and_then(|v| v.checked_add(self.window))
                .ok_or_else(|| anyhow::anyhow!("BLASTX diagonal offset exceeds Int4"))?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:273-313
// ```c++
//  */
// static int score_compare_match(const void *v1, const void *v2)
// {
//     BlastInitHSP *h1, *h2;
//     int result = 0;
//
//     h1 = (BlastInitHSP *) v1;
//     h2 = (BlastInitHSP *) v2;
//
//     /* Check if ungapped_data substructures are initialized. If not, move
//        those array elements to the end. In reality this should never happen. */
//     if (h1->ungapped_data == NULL && h2->ungapped_data == NULL)
//         return 0;
//     else if (h1->ungapped_data == NULL)
//         return 1;
//     else if (h2->ungapped_data == NULL)
//         return -1;
//
//     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
//                                  h1->ungapped_data->score)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
//                                  h2->ungapped_data->s_start)) &&
//         0 == (result = BLAST_CMP(h2->ungapped_data->length,
//                                  h1->ungapped_data->length)) &&
//         0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
//                                  h2->ungapped_data->q_start))) {
//         result = BLAST_CMP(h2->ungapped_data->length,
//                            h1->ungapped_data->length);
//     }
//
//     return result;
// }
//
// void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
// {
//     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
//           sizeof(BlastInitHSP), score_compare_match);
// }
//
// Boolean Blast_InitHitListIsSortedByScore(BlastInitHitList * init_hitlist)
// {
// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InitHsp {
    pub q_seed: i32,
    pub s_seed: i32,
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:484-496
// ```c++
//         scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
//         wordsize = lookup->word_length;
//         use_pssm = lookup->use_pssm;
//     }
//     else {
//         BlastCompressedAaLookupTable *lookup =
//                         (BlastCompressedAaLookupTable *)(lookup_wrap->lut);
//         scansub = (TAaScanSubjectFunction)(lookup->scansub_callback);
//         wordsize = lookup->word_length;
//     }
//
//     scan_range[0] = 0;
//     scan_range[1] = subject->seq_ranges[0].left;
// ```
pub enum Lookup {
    Standard(BlastAaLookupTable),
    Compressed(BlastCompressedAaLookupTable),
}
impl Lookup {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469
    // ```c++
    //     /* create an empty backbone */
    //
    //     exact_backbone = (Int4 **) calloc(lookup->backbone_size, sizeof(Int4 *));
    //
    //     /* find all the exact matches, grouping together all offsets of identical
    //        query words. The query bias is not used here, since the next stage
    //        will need real offsets into the query sequence */
    //
    //     BlastLookupIndexQueryExactMatches(exact_backbone, lookup->word_length,
    //                                       lookup->charsize, lookup->word_length,
    //                                       query, location);
    //
    //     /* walk though the list of exact matches previously computed. Find
    //        neighboring words for entire lists at a time */
    //
    //     for (i = 0; i < lookup->backbone_size; i++) {
    //         if (exact_backbone[i] != NULL) {
    //             s_AddWordHits(lookup, matrix, query->sequence,
    //                           exact_backbone[i], query_bias, row_max);
    //             sfree(exact_backbone[i]);
    //         }
    //     }
    //
    //     sfree(exact_backbone);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1349-1351
    // ```c++
    //
    //     /* index the query and finish up */
    //
    // ```
    pub fn new(
        batch: &PreparedQueryBatch,
        parameters: &[ContextParameters],
        options: &ResolvedOptions,
    ) -> Result<Self> {
        let query = batch.sequence_start[1..].to_vec();
        if options.word_size == 5 {
            return Ok(Self::Compressed(
                build_blosum62_compressed_lookup(
                    5,
                    options.threshold,
                    &query,
                    &batch.lookup_segments,
                )
                .ok_or_else(|| anyhow::anyhow!("unsupported BLASTX compressed lookup"))?,
            ));
        }
        ensure!(
            options.word_size == 3,
            "unsupported BLASTX lookup word size"
        );
        let contexts = batch
            .contexts
            .iter()
            .zip(parameters)
            .enumerate()
            .map(|(i, (c, p))| {
                let mut aa = vec![0];
                aa.extend_from_slice(&query[c.offset..c.offset + c.length]);
                aa.push(0);
                QueryContext {
                    q_idx: c.query_index as u32,
                    f_idx: (i % 6) as u8,
                    frame: c.frame,
                    aa_seq: aa,
                    aa_seq_nomask: None,
                    aa_len: c.length,
                    orig_len: batch.original_lengths[c.query_index],
                    frame_base: c.offset as i32,
                    is_valid: p.valid,
                    karlin_params: p.ungapped,
                }
            })
            .collect();
        Ok(Self::Standard(build_ncbi_lookup_from_prepared(
            query,
            batch.lookup_segments.clone(),
            contexts,
            options.threshold as i32,
        )))
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
    // ```c++
    // Int4 GetOffsetArraySize(LookupTableWrap* lookup)
    // {
    //    Int4 offset_array_size;
    //
    //    switch (lookup->lut_type) {
    //    case eMBLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastMBLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastAaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eCompressedAaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastCompressedAaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eSmallNaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastSmallNaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eNaLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastNaLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    case eNaHashLookupTable:
    //       offset_array_size = OFFSET_ARRAY_SIZE +
    //          ((BlastNaHashLookupTable*)lookup->lut)->longest_chain;
    //       break;
    //    default:
    //       offset_array_size = OFFSET_ARRAY_SIZE;
    //       break;
    //    }
    //    return offset_array_size;
    // ```
    pub fn capacity(&self) -> usize {
        4096 + match self {
            Self::Standard(l) => l.longest_chain as usize,
            Self::Compressed(l) => l.longest_chain as usize,
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
    // ```c++
    //         wordsize = lookup->word_length;
    //     }
    //
    //     scan_range[0] = 0;
    //     scan_range[1] = subject->seq_ranges[0].left;
    //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
    //
    //     if (scan_range[2] < scan_range[1])
    //         scan_range[2] = scan_range[1];
    //
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    //         hits = scansub(lookup_wrap, subject,
    //                                   offset_pairs, array_size, scan_range);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-769
    // ```c++
    //     scan_range[0] = 0;
    //     scan_range[1] = subject->seq_ranges[0].left;
    //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
    //
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    //         hits = scansub(lookup_wrap, subject,
    //                        offset_pairs, array_size, scan_range);
    // ```
    pub fn scan(
        &self,
        subject: &[u8],
        length: usize,
        window: i32,
        mut visit: impl FnMut(&[(i32, i32)]),
    ) {
        let word = match self {
            Self::Standard(l) => l.word_length,
            Self::Compressed(l) => l.word_length,
        };
        let ranges = [(0, length as i32)];
        let mut range = [0, 0, length as i32 - word];
        if window > 0 && range[2] < range[1] {
            range[2] = range[1];
        }
        let mut pairs = vec![OffsetPair::default(); self.capacity()];
        let capacity = pairs.len() as i32;
        let mut seeds = Vec::new();
        while range[1] <= range[2] {
            seeds.clear();
            let n = match self {
                Self::Standard(l) => {
                    s_blast_aa_scan_subject_one_range(l, subject, &mut pairs, capacity, &mut range)
                }
                Self::Compressed(l) => {
                    l.scan_subject(subject, &ranges, &mut pairs, capacity, &mut range)
                }
            };
            for pair in &pairs[..n as usize] {
                seeds.push((pair.q_off as i32, pair.s_off as i32));
            }
            visit(&seeds);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-607
// ```c++
//             diag_coord = (query_offset - subject_offset) & diag_mask;
//
//             /* If the reset bit is set, an extension just happened. */
//             if (diag_array[diag_coord].flag) {
//                 /* If we've already extended past this hit, skip it. */
//                 if ((Int4) (subject_offset + diag_offset) <
//                     diag_array[diag_coord].last_hit) {
//                     continue;
//                 }
//                 /* Otherwise, start a new hit. */
//                 else {
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     diag_array[diag_coord].flag = 0;
//                 }
//             }
//             /* If the reset bit is cleared, try to start an extension. */
//             else {
//                 /* find the distance to the last hit on this diagonal */
//                 last_hit = diag_array[diag_coord].last_hit - diag_offset;
//                 diff = subject_offset - last_hit;
//
//                 if (diff >= window) {
//                     /* We are beyond the window for this diagonal; start a
//                        new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
//
//                 /* If the difference is less than the wordsize (i.e. last
//                    hit and this hit overlap), give up */
//
//                 if (diff < wordsize) {
//                     continue;
//                 }
//
//                 /* Extend this pair of hits. The extension to the left must
//                    reach the end of the first word in order for extension to
//                    the right to proceed.
//
//                    To use the cutoff and X-drop values appropriate for this
//                    extension, the query context must first be found */
//
//                 curr_context = BSearchContextInfo(query_offset, query_info);
//
//                 /* Check if the last hit hits current query. Because last_hit
//                    is never reset, it may contain a hit to the previous
//                    concatenated query */
//
//                 if (query_offset - diff <
//                     query_info->contexts[curr_context].query_offset) {
//
//                     /* there was no last hit for this diagnol; start a new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
//
//                 cutoffs = word_params->cutoffs + curr_context;
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
//                                               last_hit + wordsize,
//                                               subject_offset, query_offset,
//                                               cutoffs->x_dropoff,
//                                               &hsp_q, &hsp_s,
//                                               &hsp_len, use_pssm,
//                                               wordsize, &right_extend,
//                                               &s_last_off);
//
//                 ++hits_extended;
//
//                 /* if the hsp meets the score threshold, report it */
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                                      query_offset, subject_offset, hsp_len,
//                                      score);
//
//                 /* If an extension to the right happened, reset the last hit
//                    so that future hits to this diagonal must start over. */
//
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
//                 }
//                 /* Otherwise, make the present hit into the previous hit for
//                    this diagonal */
//                 else {
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                 }
//             }                   /* end else */
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:774-800
// ```c++
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//             diag_coord = (subject_offset - query_offset) & diag_mask;
//             diff = subject_offset -
//                 (diag_array[diag_coord].last_hit - diag_offset);
//
//             /* do an extension, but only if we have not already extended this
//                far */
//             if (diff >= 0) {
//                 Int4 curr_context = BSearchContextInfo(query_offset,
//                                                        query_info);
//                 BlastUngappedCutoffs *cutoffs = word_params->cutoffs +
//                                                         curr_context;
//                 score = s_BlastAaExtendOneHit(matrix, subject, query,
//                                               subject_offset, query_offset,
//                                               cutoffs->x_dropoff,
//                                               &hsp_q, &hsp_s, &hsp_len,
//                                               wordsize, use_pssm, &s_last_off);
//
//                 /* if the hsp meets the score threshold, report it */
//                 if (score >= cutoffs->cutoff_score) {
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
//                                      query_offset, subject_offset, hsp_len,
//                                      score);
//                 }
//                 diag_array[diag_coord].last_hit =
//                     s_last_off - (wordsize - 1) + diag_offset;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-515
// ```c++
//         scan_range[2] = scan_range[1];
//
//     while (scan_range[1] <= scan_range[2]) {
//         /* scan the subject sequence for hits */
//         hits = scansub(lookup_wrap, subject,
//                                   offset_pairs, array_size, scan_range);
//
//         totalhits += hits;
//         /* for each hit, */
//         for (i = 0; i < hits; ++i) {
//             Uint4 query_offset = offset_pairs[i].qs_offsets.q_off;
//             Uint4 subject_offset = offset_pairs[i].qs_offsets.s_off;
//
//             /* calculate the diagonal associated with this query-subject pair
//              */
//
// ```
pub fn word_finder(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    mut trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    // EXPERIMENT (LOSAT_X_SEEDBUCKET + LOSAT_X_BXSEEDBUCKET): see x_seed_bucket::blastx_mode.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    // ...
    //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
    // ...
    //                 if (score >= cutoffs->cutoff_score)
    //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    // ...
    //                 if (right_extend) {
    //                     diag_array[diag_coord].flag = 1;
    //                     diag_array[diag_coord].last_hit =
    //                         s_last_off - (wordsize - 1) + diag_offset;
    // ```
    // Dispatch point: mode 0 runs word_finder_plain (the port of this loop), mode 1 runs the
    // diagonal-bucketed variant, and the shadow mode runs both on copies of the diagonal table and
    // asserts equal hit lists and tables. Same hits, same per-diagonal order; see x_seed_bucket.
    match crate::algorithm::tblastx::x_seed_bucket::blastx_mode() {
        0 => word_finder_plain(
            batch, parameters, options, subject, lookup, diagonals, trace,
        ),
        1 => word_finder_bucketed(
            batch, parameters, options, subject, lookup, diagonals, trace,
        ),
        _ => {
            // LOSAT_X_SEEDBUCKETSHADOW: both orders, compared.
            let mut copy = Diagonals {
                entries: diagonals.entries.clone(),
                offset: diagonals.offset,
                window: diagonals.window,
            };
            let shadow = word_finder_bucketed(
                batch,
                parameters,
                options,
                subject,
                lookup,
                &mut copy,
                |_| {},
            )?;
            let reference = word_finder_plain(
                batch, parameters, options, subject, lookup, diagonals, &mut trace,
            )?;
            assert!(
                copy.offset == diagonals.offset && copy.entries == diagonals.entries,
                "LOSAT_X_SEEDBUCKETSHADOW (blastx): diagonal table differs"
            );
            assert!(
                shadow == reference,
                "LOSAT_X_SEEDBUCKETSHADOW (blastx): hit list differs ({} vs {} hits)",
                shadow.len(),
                reference.len()
            );
            crate::algorithm::tblastx::x_seed_bucket::SHADOW_CHUNKS
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            crate::algorithm::tblastx::x_seed_bucket::SHADOW_HITS
                .fetch_add(reference.len() as u64, std::sync::atomic::Ordering::Relaxed);
            Ok(reference)
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
// ...
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
// ...
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
// ```
// word_finder_plain is the original port of this loop for one subject (it is unchanged except that
// the context lookup can be moved after the diagonal tests, see LOSAT_X_BXLAZYCTX).
/// The reference word finder (the loop of NCBI's BlastAaWordFinder_TwoHit /
/// _OneHit for one subject), kept verbatim; `word_finder_bucketed` is the
/// LOSAT_X_SEEDBUCKET variant of it.
fn word_finder_plain(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    mut trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let mask = diagonals.entries.len() as i32 - 1;
    let mut hits = Vec::new();
    // No NCBI counterpart: reads the LOSAT_X_BXLAZYCTX switch once; it does not change any value
    // NCBI computes.
    let lazy_context = {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXLAZYCTX").is_some())
    };
    lookup.scan(
        subject,
        subject.len().saturating_sub(1),
        options.window_size,
        |seeds| {
            trace(seeds);
            for &(q, s) in seeds {
                // EXPERIMENT (LOSAT_X_BXLAZYCTX): the context of a seed is only
                // needed once the diagonal tests let it through; looking it up
                // has no side effect, so doing it later changes nothing.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,783-784
                // ```c
                //                 curr_context = BSearchContextInfo(query_offset, query_info);
                // ...
                //                 Int4 curr_context = BSearchContextInfo(query_offset,
                //                                                        query_info);
                // ```
                // The two-hit path looks up the context (BSearchContextInfo) only after the
                // diagonal tests have let the pair through, and the one-hit path only when diff >=
                // 0. LOSAT_X_BXLAZYCTX does the same here; the lookup has no side effect.
                let context_of = |q: i32| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                };
                let eager_context = if lazy_context {
                    usize::MAX
                } else {
                    context_of(q)
                };
                let context;
                let index = if options.window_size == 0 {
                    (s - q) & mask
                } else {
                    (q - s) & mask
                } as usize;
                let (last, flag) = &mut diagonals.entries[index];
                let (u, end, extended) = if options.window_size == 0 {
                    if s - (*last - diagonals.offset) < 0 {
                        continue;
                    }
                    context = if lazy_context {
                        context_of(q)
                    } else {
                        eager_context
                    };
                    let p = &parameters[context];
                    let Some(result) = extend_one_hit_blosum62(
                        query,
                        subject,
                        q as usize,
                        s as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, true)
                } else {
                    if *flag {
                        if s + diagonals.offset < *last {
                            continue;
                        }
                        *last = s + diagonals.offset;
                        *flag = false;
                        continue;
                    }
                    let previous = *last - diagonals.offset;
                    let diff = s - previous;
                    if diff >= options.window_size {
                        *last = s + diagonals.offset;
                        continue;
                    }
                    if diff < word {
                        continue;
                    }
                    context = if lazy_context {
                        context_of(q)
                    } else {
                        eager_context
                    };
                    let p = &parameters[context];
                    if q - diff < batch.contexts[context].offset as i32 {
                        *last = s + diagonals.offset;
                        continue;
                    }
                    let Some(result) = extend_two_hit_blosum62(
                        query,
                        subject,
                        (previous + word) as usize,
                        s as usize,
                        q as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        continue;
                    };
                    (result.ungapped_data, result.s_last_off, result.right_extend)
                };
                let p = &parameters[context];
                if u.score >= p.word_cutoff {
                    hits.push(InitHsp {
                        q_seed: q,
                        s_seed: s,
                        q_start: u.q_start,
                        s_start: u.s_start,
                        length: u.length,
                        score: u.score,
                    });
                }
                if options.window_size == 0 {
                    *last = end - (word - 1) + diagonals.offset;
                } else if extended {
                    *flag = true;
                    *last = end - (word - 1) + diagonals.offset;
                } else {
                    *last = s + diagonals.offset;
                }
            }
        },
    );
    diagonals.finish_subject(subject.len().saturating_sub(1))?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
// ...
//                 if (score >= cutoffs->cutoff_score)
//                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
// ...
//                 if (right_extend) {
//                     diag_array[diag_coord].flag = 1;
//                     diag_array[diag_coord].last_hit =
//                         s_last_off - (wordsize - 1) + diag_offset;
// ```
// word_finder_bucketed processes the same hits with the same per-hit steps, grouped by diagonal
// range (LOSAT_X_SEEDBUCKET). Each hit reads and writes only diag_array[diag_coord], so the state
// of a diagonal depends only on that diagonal's hits in scan order, which is kept.
#[allow(clippy::too_many_arguments)]
fn word_finder_bucketed(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    subject: &[u8],
    lookup: &Lookup,
    diagonals: &mut Diagonals,
    mut trace: impl FnMut(&[(i32, i32)]),
) -> Result<Vec<InitHsp>> {
    let query = &batch.sequence_start[1..];
    let word = options.word_size;
    let mask = diagonals.entries.len() as i32 - 1;
    let mut hits = Vec::new();
    // No NCBI counterpart: reads the LOSAT_X_BXLAZYCTX switch once; it does not change any value
    // NCBI computes.
    let lazy_context = {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var_os("LOSAT_X_BXLAZYCTX").is_some())
    };
    // EXPERIMENT (LOSAT_X_SEEDBUCKET): the hits of a scan window grouped by
    // diagonal range before the per-diagonal tests (see tblastx::x_seed_bucket;
    // the tests below read and write one diagonal cell per hit, so the order
    // across diagonals is free; the order within a diagonal is kept, and the
    // saved HSPs are put back into scan order before the stable sort).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:500-509
    // ```c
    //         scan_range[2] = scan_range[1];
    // ...
    //     while (scan_range[1] <= scan_range[2]) {
    //         /* scan the subject sequence for hits */
    // ...
    //
    //         totalhits += hits;
    //         /* for each hit, */
    //         for (i = 0; i < hits; ++i) {
    // ```
    // NCBI handles the hits of one scansub() call in scan order. Here the hits of a scan window are
    // collected into diagonal-range buckets first and flushed when the budget is reached; the saved
    // HSPs are put back into scan order before the sort.
    let mut buckets = crate::algorithm::tblastx::x_seed_bucket::SeedBuckets::new(
        diagonals.entries.len() as u32,
        mask as u32,
    );
    let x_budget = crate::algorithm::tblastx::x_seed_bucket::budget();
    let mut x_seq: u32 = 0;
    let mut seq_keys: Vec<u32> = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    // ...
    //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
    // ...
    //                 if (score >= cutoffs->cutoff_score)
    //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    // ...
    //                 if (right_extend) {
    //                     diag_array[diag_coord].flag = 1;
    //                     diag_array[diag_coord].last_hit =
    //                         s_last_off - (wordsize - 1) + diag_offset;
    // ```
    // x_one_hit is a copy of the loop body of word_finder_plain for one hit (the C body above),
    // expanded in the flush callbacks.
    // One hit of the NCBI loop (aa_ungapped.c:531-606 / blast_aalookup scan): a
    // copy of the loop body of `word_finder_plain`, expanded in place in the
    // flush callbacks (local names resolve here, by macro hygiene).  `break 'hit`
    // is the `continue` of that loop; `seq` is the hit's position in the scan
    // stream, kept so that the saved HSPs can be put back into scan order.
    macro_rules! x_one_hit {
        ($qq:expr, $ss:expr, $sq:expr) => {{
            let q: i32 = $qq;
            let s: i32 = $ss;
            let seq: u32 = $sq;
            'hit: {
                // EXPERIMENT (LOSAT_X_BXLAZYCTX): the context of a seed is only
                // needed once the diagonal tests let it through; looking it up
                // has no side effect, so doing it later changes nothing.
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:560,783-784
                // ```c
                //                 curr_context = BSearchContextInfo(query_offset, query_info);
                // ...
                //                 Int4 curr_context = BSearchContextInfo(query_offset,
                //                                                        query_info);
                // ```
                // The two-hit path looks up the context (BSearchContextInfo) only after the
                // diagonal tests have let the pair through, and the one-hit path only when diff >=
                // 0. LOSAT_X_BXLAZYCTX does the same here; the lookup has no side effect.
                let context_of = |q: i32| {
                    batch
                        .contexts
                        .partition_point(|c| c.offset as i32 <= q)
                        .saturating_sub(1)
                };
                let eager_context = if lazy_context {
                    usize::MAX
                } else {
                    context_of(q)
                };
                let context;
                let index = if options.window_size == 0 {
                    (s - q) & mask
                } else {
                    (q - s) & mask
                } as usize;
                let (last, flag) = &mut diagonals.entries[index];
                let (u, end, extended) = if options.window_size == 0 {
                    if s - (*last - diagonals.offset) < 0 {
                        break 'hit;
                    }
                    context = if lazy_context {
                        context_of(q)
                    } else {
                        eager_context
                    };
                    let p = &parameters[context];
                    let Some(result) = extend_one_hit_blosum62(
                        query,
                        subject,
                        q as usize,
                        s as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        break 'hit;
                    };
                    (result.ungapped_data, result.s_last_off, true)
                } else {
                    if *flag {
                        if s + diagonals.offset < *last {
                            break 'hit;
                        }
                        *last = s + diagonals.offset;
                        *flag = false;
                        break 'hit;
                    }
                    let previous = *last - diagonals.offset;
                    let diff = s - previous;
                    if diff >= options.window_size {
                        *last = s + diagonals.offset;
                        break 'hit;
                    }
                    if diff < word {
                        break 'hit;
                    }
                    context = if lazy_context {
                        context_of(q)
                    } else {
                        eager_context
                    };
                    let p = &parameters[context];
                    if q - diff < batch.contexts[context].offset as i32 {
                        *last = s + diagonals.offset;
                        break 'hit;
                    }
                    let Some(result) = extend_two_hit_blosum62(
                        query,
                        subject,
                        (previous + word) as usize,
                        s as usize,
                        q as usize,
                        p.word_xdrop,
                        word as usize,
                    ) else {
                        break 'hit;
                    };
                    (result.ungapped_data, result.s_last_off, result.right_extend)
                };
                let p = &parameters[context];
                if u.score >= p.word_cutoff {
                    hits.push(InitHsp {
                        q_seed: q,
                        s_seed: s,
                        q_start: u.q_start,
                        s_start: u.s_start,
                        length: u.length,
                        score: u.score,
                    });
                    // No NCBI counterpart: records the scan position of a saved HSP so that scan
                    // order can be restored; it does not change any value NCBI computes.
                    seq_keys.push(seq);
                }
                if options.window_size == 0 {
                    *last = end - (word - 1) + diagonals.offset;
                } else if extended {
                    *flag = true;
                    *last = end - (word - 1) + diagonals.offset;
                } else {
                    *last = s + diagonals.offset;
                }
            }
        }};
    }
    lookup.scan(
        subject,
        subject.len().saturating_sub(1),
        options.window_size,
        |seeds| {
            trace(seeds);
            for &(q, s) in seeds {
                buckets.push(q as u32, s as u32, x_seq);
                x_seq += 1;
            }
            if buckets.flush_due(x_budget) {
                buckets.flush(|q, s, seq| x_one_hit!(q as i32, s as i32, seq));
            }
        },
    );
    buckets.flush(|q, s, seq| x_one_hit!(q as i32, s as i32, seq));
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:588-591
    // ```c
    //                 if (score >= cutoffs->cutoff_score)
    //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    //                                      query_offset, subject_offset, hsp_len,
    //                                      score);
    // ```
    // BlastSaveInitHsp appends HSPs in scan order. restore_order puts the HSPs saved in bucket
    // order back into that order, so the sort that follows sees the same input.
    crate::algorithm::tblastx::x_seed_bucket::restore_order(&mut hits, &seq_keys);
    diagonals.finish_subject(subject.len().saturating_sub(1))?;
    hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.s_start.cmp(&b.s_start))
            .then(b.length.cmp(&a.length))
            .then(a.q_start.cmp(&b.q_start))
    });
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:87-107
// ```c++
// static Int4 s_BlastDiagClear(BLAST_DiagTable * diag)
// {
//     Int4 i, n;
//     DiagStruct *diag_struct_array;
//
//     if (diag == NULL)
//         return 0;
//
//     n = diag->diag_array_length;
//
//     diag->offset = diag->window;
//
//     diag_struct_array = diag->hit_level_array;
//
//     for (i = 0; i < n; i++) {
//         diag_struct_array[i].flag = 0;
//         diag_struct_array[i].last_hit = -diag->window;
//         if (diag->hit_len_array) diag->hit_len_array[i] = 0;
//     }
//     return 0;
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
// ```c++
//     if (ewp->diag_table) {
//         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//             ewp->diag_table->offset = ewp->diag_table->window;
//             s_BlastDiagClear(ewp->diag_table);
//         } else {
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//         }
//     } else if (ewp->hash_table) {
// ```
#[cfg(test)]
mod tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    #[test]
    fn diagonal_overflow_threshold_clears_flags_and_negative_window() {
        let mut d = Diagonals::new(500, 40).unwrap();
        d.offset = i32::MAX / 4 - 1;
        d.entries[0] = (123, true);
        d.finish_subject(83).unwrap();
        assert_eq!(d.offset, i32::MAX / 4 - 1 + 83 + 40);
        assert_eq!(d.entries[0], (123, true));
        d.finish_subject(83).unwrap();
        assert_eq!(d.offset, 40);
        assert!(d.entries.iter().all(|x| *x == (-40, false)));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:49-59
    // ```c++
    //
    //         if (diag_table)
    //         {
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    // ```
    #[test]
    fn diagonal_size_power_and_int4_overflow_boundary() {
        assert_eq!(Diagonals::new(512, 0).unwrap().entries.len(), 512);
        assert_eq!(Diagonals::new(512, 1).unwrap().entries.len(), 1024);
        assert!(Diagonals::new(i32::MAX as usize, 1).is_err());
        assert!(Diagonals::new(usize::MAX, 0).is_err());
    }
}

#[cfg(test)]
mod session_f_diagonal_boundary_tests {
    use super::Diagonals;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:46-65
    // ```c++
    //         Int4 diag_array_length;
    //
    //         diag_table= (BLAST_DiagTable*) calloc(1, sizeof(BLAST_DiagTable));
    //
    //         if (diag_table)
    //         {
    //                 diag_array_length = 1;
    //                 /* What power of 2 is just longer than the query? */
    //                 while (diag_array_length < (qlen+window_size))
    //                 {
    //                         diag_array_length = diag_array_length << 1;
    //                 }
    //                 /* These are used in the word finders to shift and mask
    //                 rather than dividing and taking the remainder. */
    //                 diag_table->diag_array_length = diag_array_length;
    //                 diag_table->diag_mask = diag_array_length-1;
    //                 diag_table->multiple_hits = multiple_hits;
    //                 diag_table->offset = window_size;
    //                 diag_table->window = window_size;
    //         }
    // ```
    #[test]
    fn allocation_power_of_two_and_window_boundaries() {
        for (query, window, length) in [
            (63, 0, 64),
            (64, 0, 64),
            (65, 0, 128),
            (24, 40, 64),
            (25, 40, 128),
        ] {
            let d = Diagonals::new(query, window).unwrap();
            assert_eq!(d.entries.len(), length);
            assert_eq!(d.offset, window);
        }
        assert!(Diagonals::new(1, -1).is_err());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    #[test]
    fn reset_occurs_at_entry_threshold_and_clears_stale_flags() {
        for window in [0, 1, 39, 40, 41] {
            let mut d = Diagonals::new(65, window).unwrap();
            d.entries.fill((123, true));
            d.offset = 536_870_910; // INT4_MAX / 4 - 1, independent literal boundary.
            d.finish_subject(1).unwrap();
            assert_eq!(d.offset, 536_870_911 + window);
            assert!(d.entries.iter().all(|x| *x == (123, true)));
            d.finish_subject(5_000_000).unwrap();
            assert_eq!(d.offset, window);
            assert!(d.entries.iter().all(|x| *x == (-window, false)));
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // Independent worker slots advance only for their assigned subjects, including
    // a partial final wave; no query state is shared between slots or batches.
    #[test]
    fn uneven_wave_slots_and_new_batch_keep_separate_offsets() {
        let mut slots = (0..4)
            .map(|_| Diagonals::new(128, 40).unwrap())
            .collect::<Vec<_>>();
        for (oid, length) in [10, 20, 30, 40, 50, 60, 70, 80, 90].into_iter().enumerate() {
            slots[oid % 4].finish_subject(length).unwrap();
        }
        assert_eq!(
            slots.iter().map(|x| x.offset).collect::<Vec<_>>(),
            [310, 200, 220, 240]
        );
        assert_eq!(Diagonals::new(128, 40).unwrap().offset, 40);
    }
}

#[cfg(test)]
mod session_f_word_finder_reset_tests {
    use super::*;
    use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
    // ```c++
    //     if (ewp->diag_table) {
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    //     } else if (ewp->hash_table) {
    // ```
    // Expected SEED/INIT rows are frozen fresh pinned-NCBI trace bytes, not Rust output.
    #[test]
    fn fresh_retained_and_threshold_reset_word_finder_match_ncbi() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
        let q = parse_fasta(&std::fs::read(dir.join("query.fna")).unwrap(), false, false).unwrap();
        let s = parse_fasta(
            &std::fs::read(dir.join("subject.faa")).unwrap(),
            true,
            false,
        )
        .unwrap();
        let mut subject = s[0]
            .sequence
            .iter()
            .copied()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect::<Vec<_>>();
        subject.push(0);
        for window in [0, 1, 39, 40, 41] {
            let args = [
                "losat",
                "blastx",
                "-query",
                "query",
                "-subject",
                "subject",
                "-seg",
                "no",
                "-comp_based_stats",
                "0",
                "-window_size",
                &window.to_string(),
            ];
            let crate::cli::Commands::Blastx(args) =
                crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                    .unwrap()
                    .command
            else {
                unreachable!()
            };
            let options = args.resolve().unwrap();
            let batch = prepare_queries(&q, &options).unwrap();
            let mut params = parameters::score_block(&batch);
            parameters::effective_lengths(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                None,
            )
            .unwrap();
            let lookup = Lookup::new(&batch, &params, &options).unwrap();
            parameters::subject_parameters(
                &batch,
                &mut params,
                &options,
                s[0].sequence.len() as i64,
                1,
                s[0].sequence.len(),
                None,
            )
            .unwrap();
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            let expected =
                std::fs::read_to_string(dir.join(format!("window{window}.stage"))).unwrap();
            let seeds = expected
                .lines()
                .filter(|x| x.starts_with("SEED\t"))
                .map(|x| {
                    let f = x.split('\t').collect::<Vec<_>>();
                    (f[2].parse::<i32>().unwrap(), f[3].parse::<i32>().unwrap())
                })
                .collect::<Vec<_>>();
            let initial = expected
                .lines()
                .filter(|x| x.starts_with("INIT\t"))
                .map(|x| {
                    let f = x.split('\t').collect::<Vec<_>>();
                    InitHsp {
                        q_seed: f[3].parse().unwrap(),
                        s_seed: f[4].parse().unwrap(),
                        q_start: f[5].parse().unwrap(),
                        s_start: f[6].parse().unwrap(),
                        length: f[7].parse().unwrap(),
                        score: f[8].parse().unwrap(),
                    }
                })
                .collect::<Vec<_>>();
            let mut d = Diagonals::new(length, window).unwrap();
            for state in 0..4 {
                if state == 2 {
                    d.offset = 536_870_910;
                    d.entries.fill((123, true));
                }
                if state == 3 {
                    d.offset = 536_870_911;
                    d.entries.fill((536_870_900, true));
                    d.finish_subject(83).unwrap();
                }
                let mut actual = Vec::new();
                let hits = word_finder(&batch, &params, &options, &subject, &lookup, &mut d, |x| {
                    actual.extend_from_slice(x)
                })
                .unwrap();
                assert_eq!(actual, seeds, "window={window} state={state}");
                assert_eq!(hits, initial, "window={window} state={state}");
            }
        }
    }
}
