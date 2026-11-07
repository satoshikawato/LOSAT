//! Mask helpers shared by the programs: the lower-case runs of a record's letters, and the mask
//! coordinate conversions between a DNA query and its translated frames (`dna_to_protein_masks`,
//! `protein_to_dna_masks`) and the merge of overlapping masks (`combine_masks`).
//!
//! `lowercase_runs` moved from `tblastn/archive.rs` and the translated-mask functions from
//! `blastx/query_setup.rs` (SW6 items 2 and 3) so that `algorithm/common`
//! depends on no program.

/// The runs of lower-case letters of a record (inclusive 0-based positions).
pub(crate) fn lowercase_runs(letters: &[u8]) -> Vec<(i32, i32)> {
    let mut runs = Vec::new();
    let mut start = None;
    for (index, letter) in letters.iter().enumerate() {
        if letter.is_ascii_lowercase() {
            start.get_or_insert(index);
        } else if let Some(from) = start.take() {
            runs.push((from as i32, index as i32 - 1));
        }
    }
    if let Some(from) = start {
        runs.push((from as i32, letters.len() as i32 - 1));
    }
    runs
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:1604-1635
// ```c++
//             short frame = iter->first;
//             BlastSeqLoc * bsl = iter->second;
//
//             for (BlastSeqLoc* itr = bsl; itr; itr = itr->next) {
//                 int to(0), from(0);
//
//                 if (frame < 0) {
//                     from = ((int) dna_length + frame - itr->ssr->right) / CODON_LENGTH;
//                     to = ((int) dna_length + frame - itr->ssr->left) / CODON_LENGTH;
//                 } else {
//                     from = (itr->ssr->left - frame + 1) / CODON_LENGTH;
//                     to = (itr->ssr->right - frame + 1) / CODON_LENGTH;
//                 }
//                 if (from < 0)
//                     from = 0;
//                 if (to < 0)
//                     to = 0;
//                 const int kFrameLength = frame_lengths[(CSeqLocInfo::ETranslationFrame)frame];
//                 if (from >= kFrameLength)
//                     from = kFrameLength - 1;
//                 if (to >= kFrameLength)
//                     to = kFrameLength - 1;
//
//                 _ASSERT(from >= 0 && to >= 0);
//                 _ASSERT(from < kFrameLength && to < kFrameLength);
//                 itr->ssr->left  = from;
//                 itr->ssr->right = to;
//             }
//         }
//     }
// }
//
// ```
pub fn dna_to_protein_masks(masks: &[(i32, i32)], dna_length: usize, frame: i8) -> Vec<(i32, i32)> {
    let length = dna_length as i32;
    let frame = frame as i32;
    let aa_length = (length - (frame.abs() - 1)).max(0) / 3;
    masks
        .iter()
        .map(|&(left, right)| {
            let (from, to) = if frame < 0 {
                ((length + frame - right) / 3, (length + frame - left) / 3)
            } else {
                ((left - frame + 1) / 3, (right - frame + 1) / 3)
            };
            (from.max(0).min(aa_length - 1), to.max(0).min(aa_length - 1))
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:928-954
// ```c++
//                } else {
//                    from = CODON_LENGTH*seq_range->left + frame - 1;
//                    to = CODON_LENGTH*seq_range->right + frame - 1;
//                }
//
//                if (from < 0)
//                    from = 0;
//                if (to   < 0)
//                    to   = 0;
//                if (from >= dna_length)
//                    from = dna_length - 1;
//                if (to   >= dna_length)
//                    to   = dna_length - 1;
//
//                ASSERT(from >= 0);
//                ASSERT(to   >= 0);
//                ASSERT(from < dna_length);
//                ASSERT(to   < dna_length);
//
//                seq_range->left = from;
//                seq_range->right = to;
//            }
//        }
//    }
//    return status;
// }
//
// ```
pub fn protein_to_dna_masks(masks: &[(i32, i32)], dna_length: usize, frame: i8) -> Vec<(i32, i32)> {
    let length = dna_length as i32;
    let frame = frame as i32;
    masks
        .iter()
        .map(|&(left, right)| {
            let (from, to) = if frame < 0 {
                (length - 3 * right + frame + 1, length - 3 * left + frame)
            } else {
                (3 * left + frame - 1, 3 * right + frame - 1)
            };
            (from.max(0).min(length - 1), to.max(0).min(length - 1))
        })
        .collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:983-1004
// ```c++
//     qsort(ptrs, (size_t)num_elems, sizeof(*ptrs),
//           s_SeqRangeSortByStartPosition);
//
//     /* Merge the overlapping elements */
//     {
//         BlastSeqLoc* curr_tail = *mask_loc = ptrs[0];
//         for (i = 0; i < num_elems - 1; i++) {
//             const SSeqRange* next_ssr = ptrs[i+1]->ssr;
//             const Int4 stop = curr_tail->ssr->right;
//
//             if ((stop + link_value) > next_ssr->left) {
//                 curr_tail->ssr->right = MAX(stop, next_ssr->right);
//                 ptrs[i+1] = BlastSeqLocNodeFree(ptrs[i+1]);
//             } else {
//                 curr_tail = ptrs[i+1];
//             }
//         }
//     }
//
//     /* Rebuild the linked list */
//     {
//         BlastSeqLoc* tail = *mask_loc;
// ```
pub fn combine_masks(mut masks: Vec<(i32, i32)>) -> Vec<(i32, i32)> {
    masks.sort_by_key(|r| r.0);
    let mut result: Vec<(i32, i32)> = Vec::new();
    for mask in masks {
        if let Some(last) = result.last_mut().filter(|r| r.1 > mask.0) {
            last.1 = last.1.max(mask.1);
        } else {
            result.push(mask);
        }
    }
    result
}
