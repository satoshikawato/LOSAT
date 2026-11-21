use super::*;

// No NCBI counterpart: checks the order bookkeeping of the bucketed two-hit stage; it does
// not compute any value NCBI computes.
#[test]
fn restore_flush_order_sorts_only_the_last_flush() {
    // Three items of earlier flushes, already in scan order, then four saved in bucket order.
    let mut items: Vec<u32> = vec![100, 101, 102, 13, 10, 12, 11];
    let mut keys: Vec<u32> = vec![7, 0, 5, 2];
    restore_flush_order(&mut items, 3, &mut keys);
    assert_eq!(items, vec![100, 101, 102, 10, 11, 12, 13]);
    assert!(keys.is_empty());

    // Keys up to the bound of a flush (2^32 - 3, see `budget`) order as unsigned numbers.
    let mut items: Vec<u32> = vec![1, 3, 2, 0];
    let mut keys: Vec<u32> = vec![3_000_000_000, u32::MAX - 2, 4_000_000_000, 5];
    restore_flush_order(&mut items, 0, &mut keys);
    assert_eq!(items, vec![0, 1, 2, 3]);

    // Nothing or one item saved: unchanged, keys emptied.
    let mut items: Vec<u32> = vec![5, 4];
    let mut keys: Vec<u32> = Vec::new();
    restore_flush_order(&mut items, 2, &mut keys);
    assert_eq!(items, vec![5, 4]);
    let mut keys: Vec<u32> = vec![9];
    restore_flush_order(&mut items, 1, &mut keys);
    assert_eq!(items, vec![5, 4]);
    assert!(keys.is_empty());
}

// No NCBI counterpart: drives the buckets the way the two callers do (push the hits of a
// scan call, flush when the budget is reached, restore after each flush) and checks that the
// saved items come out in scan order although the position restarts at every flush.
#[test]
fn saved_items_are_in_scan_order_across_flushes() {
    // 2^20 cells: several buckets (32 with the default bucket size), so a flush hands the
    // hits out in an order that is not the scan order.
    let diag_array_size: u32 = 1 << 20;
    let mut buckets = SeedBuckets::new(diag_array_size, diag_array_size - 1);
    assert!(
        buckets.buckets.len() > 1,
        "the test needs more than one bucket"
    );
    let budget = 1000usize;
    // Flushes whose saved items were out of scan order before the restore.
    let mut flushes_out_of_order = 0u32;
    // (scan position, per-flush position) of every saved hit.
    let mut saved: Vec<(u64, u32)> = Vec::new();
    let mut seq_keys: Vec<u32> = Vec::new();
    let mut flushes = 0u32;
    let mut largest_seq = 0u32;
    let mut scan_position = 0u64;
    // A deterministic stream: 300 scan calls of 0 to 96 hits on scattered diagonals.
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    // The hits buffered since the last flush, in scan order, with their scan position.
    let mut pending: Vec<(u32, u32, u64)> = Vec::new();
    for call in 0..300 {
        let hits = (next() % 97) as usize;
        for _ in 0..hits {
            let q_off = (next() % (1 << 20)) as u32;
            let s_off = (next() % (1 << 20)) as u32;
            buckets.push(q_off, s_off);
            pending.push((q_off, s_off, scan_position));
            scan_position += 1;
        }
        let last_call = call == 299;
        if buckets.flush_due(budget) || last_call {
            let start = saved.len();
            let flushed = buckets.total;
            let mut seen = 0usize;
            buckets.flush(|q, s, seq| {
                // The per-flush number is the index into this flush's hits.
                let (pq, ps, position) = pending[seq as usize];
                assert_eq!((pq, ps), (q, s));
                seen += 1;
                largest_seq = largest_seq.max(seq);
                // Save about a third of the hits.
                if position % 3 == 0 {
                    saved.push((position, seq));
                    seq_keys.push(seq);
                }
            });
            assert_eq!(seen, flushed);
            if saved[start..].windows(2).any(|w| w[0].0 > w[1].0) {
                flushes_out_of_order += 1;
            }
            restore_flush_order(&mut saved, start, &mut seq_keys);
            pending.clear();
            flushes += 1;
        }
    }
    assert!(flushes > 5, "the stream must span several flushes");
    assert!(
        flushes_out_of_order > 5,
        "the buckets must have changed the order within several flushes"
    );
    // The position restarted: it never reached the number of hits of the whole stream.
    assert!(u64::from(largest_seq) < scan_position / 4);
    assert!((largest_seq as usize) < budget + 96);
    // Every third hit was saved, in scan order.
    let positions: Vec<u64> = saved.iter().map(|&(p, _)| p).collect();
    let expected: Vec<u64> = (0..scan_position).filter(|p| p % 3 == 0).collect();
    assert_eq!(positions, expected);
}

// No NCBI counterpart: the cap that keeps the per-flush position within a u32.
#[test]
fn budget_cap_keeps_the_flush_position_in_u32() {
    assert_eq!(capped_budget(None), 1 << 21);
    assert_eq!(capped_budget(Some(4096)), 4096);
    assert_eq!(capped_budget(Some(MAX_BUDGET)), MAX_BUDGET);
    assert_eq!(capped_budget(Some(usize::MAX)), MAX_BUDGET);
    // budget - 1 buffered hits, then one scan call of at most i32::MAX hits.
    let most_buffered = (capped_budget(Some(usize::MAX)) as u64 - 1) + i32::MAX as u64;
    assert!(most_buffered < u64::from(u32::MAX));
}
