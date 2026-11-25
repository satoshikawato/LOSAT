use super::*;

fn contexts_of(lengths: &[usize]) -> (Vec<QueryContext>, usize) {
    // the layout of run.rs: a context after the one before and its sentinel, none for no letters
    let mut out = Vec::new();
    let mut offset = 0usize;
    for (i, &len) in lengths.iter().enumerate() {
        out.push(QueryContext {
            query_idx: (i / 2) as u32,
            frame: if i % 2 == 0 { 1 } else { -1 },
            query_offset: offset as i32,
            len,
            masks: Vec::new(),
        });
        offset += if len == 0 { 0 } else { len + 1 };
    }
    let total = out
        .iter()
        .map(|c| c.query_offset as usize + c.len)
        .max()
        .unwrap_or(0);
    (out, total)
}

// the table of QueryContextIndex::new
fn direct_map(contexts: &[QueryContext], total: usize) -> Vec<u32> {
    let mut map = vec![0u32; total];
    for (idx, ctx) in contexts.iter().enumerate() {
        let start = ctx.query_offset.max(0) as usize;
        let end = start + ctx.len;
        if start < end && end <= total {
            map[start..end].fill(idx as u32);
        }
    }
    map
}

// every offset: the table's value, or `None` in a block with two sentinels or more
fn check(lengths: &[usize], force_blocks: bool) -> &'static str {
    let (contexts, total) = contexts_of(lengths);
    let form = CtxO1::build_with(&contexts, total, force_blocks).expect("form");
    let (name, _, shift) = form.describe();
    let map = direct_map(&contexts, total);
    let sentinels: Vec<usize> = contexts
        .iter()
        .filter(|c| c.len > 0)
        .map(|c| c.query_offset as usize + c.len)
        .collect();
    let mut unanswered = 0usize;
    for (n, &want) in map.iter().enumerate() {
        match form.get(n) {
            Some(got) => assert_eq!(
                got, want as usize,
                "offset {n} of lengths {lengths:?} (force_blocks {force_blocks})"
            ),
            None => {
                let block = n >> shift;
                let in_block = sentinels.iter().filter(|&&e| e >> shift == block).count();
                assert!(
                    name == "blocks" && form.needs_runs() && in_block >= 2,
                    "offset {n}"
                );
                unanswered += 1;
            }
        }
    }
    assert_eq!(form.needs_runs(), unanswered > 0);
    assert_eq!(form.get(total), None);
    assert_eq!(form.get(total + 1000), None);
    // EXPERIMENT (LOSAT_X_DCEMITHOIST): the form as a value answers as `get` at every offset
    // (and past the end), in the variant of its form
    let view = form.x_view();
    let view_get = |n: usize| match &view {
        XView::Div(v) => v.get(n),
        XView::Blocks(v) => v.get(n),
        XView::BlocksRanks(v) => v.get(n),
    };
    match (&view, name) {
        (XView::Div(_), "div") | (XView::Blocks(_) | XView::BlocksRanks(_), "blocks") => {}
        _ => panic!("view of a {name} form"),
    }
    let empty_contexts = contexts.iter().any(|c| c.len == 0);
    if name == "blocks" {
        assert_eq!(matches!(view, XView::BlocksRanks(_)), empty_contexts);
    }
    for n in (0..total + 3).chain([total + 1000, 1 << 31, usize::MAX >> 1]) {
        assert_eq!(view_get(n), form.get(n), "view at offset {n} of {lengths:?}");
    }
    assert_eq!(XNone.get(0), None);
    name
}

fn lcg(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

#[test]
fn equal_lengths_divide_and_blocks_agree_with_the_table() {
    for len in [1usize, 2, 3, 7, 15, 16, 17, 24, 31, 32, 33, 150, 1000] {
        for queries in [1usize, 2, 3, 100] {
            let lengths = vec![len; 2 * queries];
            assert_eq!(check(&lengths, false), "div");
            assert_eq!(check(&lengths, true), "blocks");
        }
    }
}

#[test]
fn mixed_lengths_and_empty_contexts_agree_with_the_table() {
    let mut state = 0x5eed_u64;
    for round in 0..400 {
        let queries = 1 + (lcg(&mut state) % 40) as usize;
        let (lo, span) = match round % 4 {
            0 => (18usize, 13usize),
            1 => (1, 4),
            2 => (1, 300),
            _ => (60, 2000),
        };
        let mut lengths = Vec::new();
        for _ in 0..queries {
            let len = if lcg(&mut state).is_multiple_of(10) {
                0
            } else {
                lo + (lcg(&mut state) as usize % span)
            };
            lengths.push(len);
            lengths.push(len);
        }
        if lengths.iter().all(|&l| l == 0) {
            lengths[0] = 5;
            lengths[1] = 5;
        }
        let form = check(&lengths, false);
        let equal = lengths.iter().all(|&l| l == lengths[0]);
        assert_eq!(form, if equal { "div" } else { "blocks" }, "{lengths:?}");
    }
}

#[test]
fn tight_shortest_lengths_agree_with_the_table() {
    // the shortest context sets the block: lengths 2^k - 1 and 2^k - 2 are the tight cases
    for k in 1..12u32 {
        let short = (1usize << k) - 1;
        for lengths in [
            vec![short, short + 5, short, 3 * short + 1],
            vec![short + 1, short, short, short + 2, short + 9, short],
            vec![0, 0, short, short + 1, 0, 0, short, short],
        ] {
            assert_eq!(check(&lengths, false), "blocks");
        }
    }
}

#[test]
fn a_short_context_among_long_ones_takes_the_long_ones_block_size() {
    let mut state = 0x0b10c_u64;
    let mut lengths = Vec::new();
    for i in 0..200 {
        let len = if i == 77 {
            5
        } else {
            300 + (lcg(&mut state) as usize % 400)
        };
        lengths.push(len);
        lengths.push(len);
    }
    assert_eq!(check(&lengths, false), "blocks");
    let (contexts, total) = contexts_of(&lengths);
    let form = CtxO1::build_with(&contexts, total, false).unwrap();
    // the shortest context alone gives blocks of 4; 256 still holds every context but two
    assert_eq!(form.describe().2, 8);
    assert!(form.needs_runs());
    // without the short one, the same blocks and none crowded
    lengths[154] = 300;
    lengths[155] = 300;
    let (contexts, total) = contexts_of(&lengths);
    let form = CtxO1::build_with(&contexts, total, false).unwrap();
    assert_eq!(form.describe().2, 8);
    assert!(!form.needs_runs());
}

#[test]
fn a_short_last_context_leaves_the_last_block_uncrowded() {
    // the end of the last context is the table's length, not a sentinel inside it: the block
    // that holds the sentinel before it and that end has one sentinel
    let mut lengths = vec![400usize; 600];
    lengths.push(20);
    assert_eq!(check(&lengths, false), "blocks");
    let (contexts, total) = contexts_of(&lengths);
    let form = CtxO1::build_with(&contexts, total, false).unwrap();
    assert_eq!(form.describe().2, 8);
    assert!(!form.needs_runs());
    // with contexts without letters after it, that end is a sentinel (0) inside the table
    lengths.extend([0, 0]);
    assert_eq!(check(&lengths, false), "blocks");
}

#[test]
fn contexts_without_letters_at_the_end_add_the_last_sentinel() {
    for lengths in [
        vec![7usize, 7, 0, 0],
        vec![5, 9, 0, 0, 0, 0],
        vec![0, 0, 30, 30, 0, 0],
    ] {
        let (contexts, total) = contexts_of(&lengths);
        let ends = contexts
            .iter()
            .filter(|c| c.len > 0)
            .map(|c| c.query_offset as usize + c.len)
            .max()
            .unwrap();
        assert_eq!(total, ends + 1, "{lengths:?}");
        assert_eq!(check(&lengths, false), "blocks");
    }
}

#[test]
fn other_layouts_are_left_to_the_previous_form() {
    // a gap of two letters between contexts
    let (mut contexts, total) = contexts_of(&[10, 10, 12, 12]);
    contexts[2].query_offset += 1;
    contexts[3].query_offset += 1;
    assert!(CtxO1::build_with(&contexts, total + 1, false).is_none());
    // overlapping contexts
    let (mut contexts, total) = contexts_of(&[10, 10]);
    contexts[1].query_offset -= 3;
    assert!(CtxO1::build_with(&contexts, total - 3, false).is_none());
    // no letters at all
    let (contexts, _) = contexts_of(&[0, 0]);
    assert!(CtxO1::build_with(&contexts, 0, false).is_none());
    // a length that does not end at the last context (or one past it, see the next test)
    let (contexts, total) = contexts_of(&[10, 10]);
    assert!(CtxO1::build_with(&contexts, total + 2, false).is_none());
    assert!(CtxO1::build_with(&contexts, total - 1, false).is_none());
}
