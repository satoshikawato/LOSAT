//! NCBI's merge of the preliminary hit lists of BLASTP query chunks (`BlastHSPStreamMerge`):
//! the chunks themselves are `blastp/protein_query_split.rs`, and each chunk's preliminary
//! stage is a search of its query parts (`blast_engine.rs`).

use std::cmp::Ordering;

use super::hsp::{score_compare_hsps, BlastpHitList, BlastpHspList, HitListHsp};

/// One query part of a chunk: the batch query it belongs to (the index of its preliminary
/// hit list among the batch's), its offset in that query and that query's length.
#[derive(Clone, Copy, Debug)]
pub(crate) struct BlastpChunkPart {
    pub batch_query: usize,
    pub q_idx: u32,
    pub offset: i32,
    pub query_length: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:476-534
/// The merge of one query chunk (RP-4): `chunk_lists[i]` is the preliminary hit list of the
/// chunk's part `parts[i]`; its HSPs move to the part's query (a blastp query has one
/// context, `query_context` its index, and frame 0) and the list merges with that query's hit
/// list of the chunks before, `merged[parts[i].batch_query]`. The split blocks of a gapped
/// search allow gaps (split_query_cxx.cpp:864-865, `GetGappedMode`).
pub(crate) fn merge_query_chunk(
    merged: &mut [Option<BlastpHitList>],
    chunk_lists: Vec<Option<BlastpHitList>>,
    parts: &[BlastpChunkPart],
    chunk_overlap_size: i32,
) {
    for (hitlist, part) in chunk_lists.into_iter().zip(parts) {
        let Some(mut hitlist) = hitlist else {
            continue;
        };
        rebase_hit_list(&mut hitlist, part.q_idx, part.offset, part.query_length);
        hit_list_merge(
            hitlist,
            &mut merged[part.batch_query],
            part.offset,
            chunk_overlap_size,
        );
    }
    for hitlist in merged.iter_mut().flatten() {
        let count = hitlist.hsplist_count;
        for list in hitlist.hsplist_array.iter_mut().take(count) {
            sort_by_score(list);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:508-516,518-526 (cited
// at `merge_query_chunk`)
/// The HSP lists of `hitlist` moved to query `q_idx`, whose letters start `offset` letters before
/// the HSPs' query coordinates, of `query_length` letters (`HitListHsp::move_to_query` on each HSP,
/// and each list's `query_index`): a query chunk's part merged into its query
/// (`merge_query_chunk`), or (LOSAT_X_BPQBLOCK) the compact list of a block's query moved to the
/// query's index among all queries (`blast_engine/prelim_run.rs`: `run_preliminary`'s loop over
/// the blocks of `x_prelim_block_plan`; offset 0: a blastp HSP's query coordinates are relative to its context).
pub(crate) fn rebase_hit_list<H: HitListHsp>(
    hitlist: &mut BlastpHitList<H>,
    q_idx: u32,
    offset: i32,
    query_length: usize,
) {
    for list in hitlist.hsplist_array.iter_mut().take(hitlist.hsplist_count) {
        for hsp in &mut list.hsps {
            hsp.move_to_query(q_idx, offset, query_length);
        }
        list.query_index = q_idx;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2132-2217
// The subjects of a hit list are distinct, so the OID sort is total. BLASTP has no
// `-max_hsps` merge limit before the traceback (`hsp_max` keeps every HSP).
fn hit_list_merge(
    mut hitlist1: BlastpHitList,
    combined: &mut Option<BlastpHitList>,
    split_offset: i32,
    chunk_overlap_size: i32,
) {
    let Some(mut hitlist2) = combined.take() else {
        *combined = Some(hitlist1);
        return;
    };
    let mut new_hitlist = BlastpHitList::new(hitlist1.hsplist_max);
    let mut lists1: Vec<BlastpHspList> = hitlist1
        .hsplist_array
        .drain(..)
        .take(hitlist1.hsplist_count)
        .collect();
    let mut lists2: Vec<BlastpHspList> = hitlist2
        .hsplist_array
        .drain(..)
        .take(hitlist2.hsplist_count)
        .collect();
    lists1.sort_by_key(|list| list.oid);
    lists2.sort_by_key(|list| list.oid);
    let query_is_split = split_offset > 0;
    let mut lists1 = lists1.into_iter().peekable();
    let mut lists2 = lists2.into_iter().peekable();
    while let (Some(list1), Some(list2)) = (lists1.peek(), lists2.peek()) {
        match list1.oid.cmp(&list2.oid) {
            Ordering::Less => new_hitlist.update(lists1.next().expect("peeked")),
            Ordering::Greater => new_hitlist.update(lists2.next().expect("peeked")),
            Ordering::Equal => {
                let list1 = lists1.next().expect("peeked");
                let mut list2 = lists2.next().expect("peeked");
                if query_is_split {
                    query_hsp_lists_merge(list1, &mut list2, split_offset, chunk_overlap_size);
                } else {
                    combine_by_score(list1, &mut list2);
                }
                new_hitlist.update(list2);
            }
        }
    }
    for list in lists1.chain(lists2) {
        new_hitlist.update(list);
    }
    *combined = Some(new_hitlist);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2921-3029
// The query split branch of `Blast_HSPListsMerge` for a protein query (one context, frame
// 0), with gaps allowed. `BlastpHsp` keeps `query.offset + 1` in `q_start` and `query.end`
// in `q_end` (likewise for the subject).
fn query_hsp_lists_merge(
    new: BlastpHspList,
    combined: &mut BlastpHspList,
    split_offset: i32,
    chunk_overlap_size: i32,
) {
    if new.hsps.is_empty() {
        return;
    }
    let mut incoming = new.hsps;
    let old = &mut combined.hsps;
    let mut hspcnt1 = 0;
    let mut hspcnt2 = 0;
    if split_offset >= 0 {
        for index1 in 0..old.len() {
            if old[index1].q_end as i64 > i64::from(split_offset) {
                old.swap(hspcnt1, index1);
                hspcnt1 += 1;
            }
        }
        for index2 in 0..incoming.len() {
            if (incoming[index2].q_start as i64 - 1)
                < i64::from(split_offset) + i64::from(chunk_overlap_size)
            {
                incoming.swap(hspcnt2, index2);
                hspcnt2 += 1;
            }
        }
    }
    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut slots: Vec<Option<super::hsp::BlastpHsp>> =
            incoming.into_iter().map(Some).collect();
        for hsp1 in old.iter_mut().take(hspcnt1) {
            for slot in slots.iter_mut().take(hspcnt2) {
                let Some(hsp2) = slot.as_ref() else {
                    continue;
                };
                if hsp1.query_context != hsp2.query_context {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1463-1478
                let end_diag = hsp1.q_end as i64 - hsp1.s_end as i64;
                let start_diag = (hsp2.q_start as i64 - 1) - (hsp2.s_start as i64 - 1);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1537
                if (end_diag - start_diag).abs() < 10
                    && crate::core::blast_hits::merge_two_chunk_hsps(hsp1, hsp2, true)
                        .expect("gapped split-query HSPs have positive length")
                {
                    *slot = None;
                }
            }
        }
        incoming = slots.into_iter().flatten().collect();
    }
    combined.hsps.extend(incoming);
    sort_by_score(combined);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2758-2766,2819-2829,2849
// `Blast_HSPListAppend` of a query that is not split here.
fn combine_by_score(new: BlastpHspList, combined: &mut BlastpHspList) {
    combined.hsps.extend(new.hsps);
    sort_by_score(combined);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1383
// Stable, as glibc's qsort under the pinned NCBI BLAST+ (TN-5): HSPs that two chunks found
// with other gapped starts keep their order.
fn sort_by_score(list: &mut BlastpHspList) {
    if list
        .hsps
        .windows(2)
        .any(|pair| score_compare_hsps(&pair[0], &pair[1]) == Ordering::Greater)
    {
        list.hsps.sort_by(score_compare_hsps);
    }
}
