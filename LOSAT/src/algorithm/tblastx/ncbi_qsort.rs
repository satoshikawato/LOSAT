//! Native qsort wrappers for NCBI parity.

use super::chaining::UngappedHit;
use std::cmp::Ordering;

pub(crate) type UngappedHitCompare = fn(&UngappedHit, &UngappedHit) -> Ordering;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
// ```c
// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
//       ScoreCompareHSPs);
// ```
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-486
// ```c
// qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
//       s_RevCompareHSPsTbx);
// ```
// Sort an index array through the platform C qsort, then replay that order onto
// Rust HSP values. NCBI sorts pointer arrays; sorting usize indices gives the
// same element size and preserves qsort's native handling of comparator-equal
// entries without moving Rust structs through C.
pub(crate) fn qsort_ungapped_hits_by(hits: &mut [UngappedHit], compare: UngappedHitCompare) {
    if hits.len() <= 1 {
        return;
    }

    #[cfg(not(target_arch = "wasm32"))]
    {
        native_qsort_ungapped_hits_by(hits, compare);
    }

    #[cfg(target_arch = "wasm32")]
    {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
        // ```c
        // qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
        //       ScoreCompareHSPs);
        // ```
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-994
        // ```c
        // qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
        //       s_RevCompareHSPsTransl);
        // qsort(link_hsp_array, total_number_of_hsps,sizeof(LinkHSPStruct*),
        //       s_FwdCompareHSPsTransl);
        // ```
        // Wasm has no C qsort. Keep the native wrapper's index-array replay
        // shape so comparator-equal entries are carried forward in the current
        // HSPList order without adding LOSAT-only tie-break fields.
        index_replay_sort_ungapped_hits_by(hits, compare);
    }
}

#[cfg(any(target_arch = "wasm32", test))]
fn index_replay_sort_ungapped_hits_by(hits: &mut [UngappedHit], compare: UngappedHitCompare) {
    let original = hits.to_vec();
    let mut indices: Vec<usize> = (0..original.len()).collect();

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:476-486
    // ```c
    // for (index = 0; index < total_number_of_hsps; ++index) {
    //    link_hsp_array[index]->hsp = hsp_array[index];
    // }
    // qsort(link_hsp_array,total_number_of_hsps,sizeof(LinkHSPStruct*),
    //       s_RevCompareHSPsTbx);
    // ```
    indices.sort_by(|&lhs_index, &rhs_index| compare(&original[lhs_index], &original[rhs_index]));

    for (dst, src) in indices.into_iter().enumerate() {
        hits[dst] = original[src].clone();
    }
}

#[cfg(not(target_arch = "wasm32"))]
mod native {
    use super::{UngappedHit, UngappedHitCompare};
    use std::cell::Cell;
    use std::cmp::Ordering;
    use std::ffi::c_void;
    use std::mem;
    use std::os::raw::c_int;

    #[derive(Clone, Copy)]
    struct QsortContext {
        hits: *const UngappedHit,
        len: usize,
        compare: UngappedHitCompare,
    }

    thread_local! {
        static QSORT_CONTEXT: Cell<Option<QsortContext>> = Cell::new(None);
    }

    unsafe extern "C" {
        fn qsort(
            base: *mut c_void,
            nmemb: usize,
            size: usize,
            compar: extern "C" fn(*const c_void, *const c_void) -> c_int,
        );
    }

    extern "C" fn compare_indices(lhs: *const c_void, rhs: *const c_void) -> c_int {
        // Safety: qsort passes pointers to elements of the usize index array
        // supplied in native_qsort_ungapped_hits_by.
        let lhs_index = unsafe { *(lhs as *const usize) };
        let rhs_index = unsafe { *(rhs as *const usize) };

        QSORT_CONTEXT.with(|cell| {
            let Some(ctx) = cell.get() else {
                return 0;
            };
            if lhs_index >= ctx.len || rhs_index >= ctx.len {
                return 0;
            }

            // Safety: ctx.hits points at `original`, which is kept alive for the
            // whole qsort call; indices are checked against ctx.len above.
            let lhs_hit = unsafe { &*ctx.hits.add(lhs_index) };
            let rhs_hit = unsafe { &*ctx.hits.add(rhs_index) };

            match (ctx.compare)(lhs_hit, rhs_hit) {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            }
        })
    }

    pub(super) fn native_qsort_ungapped_hits_by(
        hits: &mut [UngappedHit],
        compare: UngappedHitCompare,
    ) {
        let original = hits.to_vec();
        let mut indices: Vec<usize> = (0..original.len()).collect();
        let ctx = QsortContext {
            hits: original.as_ptr(),
            len: original.len(),
            compare,
        };

        QSORT_CONTEXT.with(|cell| {
            let previous = cell.replace(Some(ctx));
            // Safety: indices is a valid mutable array of usize elements; the
            // comparator only reads from `original`, which outlives this call.
            unsafe {
                qsort(
                    indices.as_mut_ptr().cast::<c_void>(),
                    indices.len(),
                    mem::size_of::<usize>(),
                    compare_indices,
                );
            }
            cell.set(previous);
        });

        for (dst, src) in indices.into_iter().enumerate() {
            hits[dst] = original[src].clone();
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
use native::native_qsort_ungapped_hits_by;

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(raw_score: i32, order: usize) -> UngappedHit {
        UngappedHit {
            q_idx: 0,
            s_idx: 0,
            ctx_idx: 0,
            s_f_idx: 0,
            q_frame: 1,
            s_frame: 1,
            q_aa_start: 0,
            q_aa_end: 10,
            s_aa_start: 0,
            s_aa_end: 10,
            q_seed_off: 0,
            s_seed_off: 0,
            q_orig_len: 10,
            s_orig_len: 10,
            raw_score,
            e_value: 0.0,
            num_ident: 0,
            hsp_list_order: order,
            ordering_method: 0,
            linked_set: false,
            start_of_chain: false,
            link_id: order,
            chain_next_link_id: None,
        }
    }

    #[test]
    fn index_replay_sort_uses_ncbi_comparator_without_extra_ties() {
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1347-1355
        // ```c
        // if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
        //     0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
        //     0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
        //     0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
        //     result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
        // }
        // ```
        let mut hits = vec![hit(20, 0), hit(30, 1), hit(30, 2), hit(10, 3)];
        index_replay_sort_ungapped_hits_by(&mut hits, |a, b| b.raw_score.cmp(&a.raw_score));

        let orders: Vec<usize> = hits.iter().map(|h| h.hsp_list_order).collect();
        assert_eq!(orders, vec![1, 2, 0, 3]);
    }
}
