use super::*;

// No NCBI counterpart: the grid limit of a tree made for a known number of HSPs
// (`x_grid_cells`).
#[test]
fn grid_follows_the_initial_hsps() {
    assert_eq!(x_grid_cells(0), X_GRID_CELLS_PER_HSP);
    assert_eq!(x_grid_cells(1), X_GRID_CELLS_PER_HSP);
    assert_eq!(x_grid_cells(10), 10 * X_GRID_CELLS_PER_HSP);
    assert_eq!(x_grid_cells(1 << 20), X_GRID_MAX_CELLS_FOR_HSPS);
    assert_eq!(x_grid_cells(usize::MAX), X_GRID_MAX_CELLS_FOR_HSPS);
}

fn make_tree_hsp(
    query_offset: i32,
    query_end: i32,
    subject_offset: i32,
    subject_end: i32,
    score: i32,
    subject_frame_sign: i32,
) -> TreeHsp {
    TreeHsp {
        query_offset,
        query_end,
        subject_offset,
        subject_end,
        score,
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:577-584,931-935,992-994
/// The test feeds the same random HSPs to a tree that uses the walk and a tree that uses the
/// side indexes. Every containment answer must agree, and the trees must end up with the
/// same nodes and the same number of unlinked HSPs. The HSPs are generated so that
/// containment, shared start points, shared end points and replacements all occur. The
/// indexed tree's grid limit changes between trials (`reset_with_bounds_cells`).
#[test]
fn x_side_indexes_agree_with_tree_walk_on_random_hsps() {
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }
    let trials: usize = std::env::var("LOSAT_FUZZ_CASES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(60);
    let mut rng = Rng(0xD6E8_FEB8_6659_FD93);
    let mut contained = 0usize;
    let mut queries = 0usize;
    let mut unlinked = 0usize;
    // One tree for all trials, as the scratch keeps it: its grid changes size between
    // trials and must come back empty each time.
    let mut indexed = BlastIntervalTree::new(0, 1, 0, 1);
    for trial in 0..trials {
        let q_len = [300i32, 5_000, 80_000, 3_000_000][trial % 4];
        let s_len = [400i32, 9_000, 120_000, 2_000_000][(trial / 4) % 4];
        let separation = [0i32, 6, 50][trial % 3];
        let mut walk = BlastIntervalTree::new(0, 2 * q_len + 1, 0, s_len + 1);
        // The walk answers while the index is dropped (as when tracing is on).
        walk.x.usable = false;
        // The grid at its own limit, at blastp's limit under LOSAT_X_BPTREECELLS, at the limit
        // of the trees made per subject or unit (`x_grid_cells`), and coarse.
        let cells = [X_MAX_CELLS, 1 << 14, X_GRID_MAX_CELLS_FOR_HSPS, 64][(trial / 16) % 4];
        if cells == X_MAX_CELLS {
            indexed.reset_with_bounds(0, 2 * q_len + 1, 0, s_len + 1);
        } else {
            indexed.reset_with_bounds_cells(0, 2 * q_len + 1, 0, s_len + 1, cells);
        }
        assert!(indexed.x.nq * indexed.x.ns <= cells, "trial {trial}");
        let mut added: Vec<(TreeHsp, i32)> = Vec::new();
        let count = 200 + rng.below(1500) as usize;
        for _ in 0..count {
            let context_offset = if rng.below(4) == 0 { q_len + 1 } else { 0 };
            let kind = rng.below(10);
            let (mut hsp, context_offset) = if kind < 3 && !added.is_empty() {
                // a piece of an earlier HSP (often contained in it)
                let (base, offset) = added[rng.below(added.len() as u64) as usize];
                let q_span = (base.query_end - base.query_offset).max(1);
                let s_span = (base.subject_end - base.subject_offset).max(1);
                let a = rng.below(q_span as u64) as i32;
                let b = a + rng.below((q_span - a) as u64 + 1) as i32;
                let shift = rng.below(5) as i32 - 2;
                let mut piece = base;
                piece.query_offset = base.query_offset + a;
                piece.query_end = base.query_offset + b;
                piece.subject_offset = (base.subject_offset + a.min(s_span) + shift).max(0);
                piece.subject_end =
                    (base.subject_offset + b.min(s_span) + shift).max(piece.subject_offset);
                piece.score = (base.score - rng.below(40) as i32 + 5).max(1);
                (piece, offset)
            } else if kind < 5 && !added.is_empty() {
                // shares the start or the end point of an earlier HSP
                let (base, offset) = added[rng.below(added.len() as u64) as usize];
                let mut other = base;
                let change = 1 + rng.below(200) as i32;
                if rng.below(2) == 0 {
                    other.query_end = (base.query_end + change - 100).max(base.query_offset);
                    other.subject_end = (base.subject_end + change - 100).max(base.subject_offset);
                } else {
                    other.query_offset =
                        (base.query_offset - change + 100).clamp(0, base.query_end);
                    other.subject_offset =
                        (base.subject_offset - change + 100).clamp(0, base.subject_end);
                }
                other.score = (base.score + rng.below(41) as i32 - 20).max(1);
                (other, offset)
            } else {
                let longest = [50u64, 800, 20_000][rng.below(3) as usize];
                let length = 5 + rng.below(longest) as i32;
                let length = length.min(q_len - 1).min(s_len - 1).max(1);
                let query_offset = rng.below((q_len - length) as u64) as i32;
                let subject_offset = rng.below((s_len - length) as u64) as i32;
                let skew = rng.below(7) as i32 - 3;
                (
                    TreeHsp {
                        query_offset,
                        query_end: query_offset + length,
                        subject_offset,
                        subject_end: (subject_offset + length + skew).max(subject_offset),
                        score: length + rng.below(30) as i32,
                        query_frame: 1,
                        query_length: q_len,
                        query_context_offset: context_offset,
                        subject_frame_sign: if rng.below(8) == 0 { -1 } else { 1 },
                    },
                    context_offset,
                )
            };
            hsp.query_context_offset = context_offset;
            // inside the ranges the tree was built for
            hsp.query_offset = hsp.query_offset.clamp(0, q_len);
            hsp.query_end = hsp.query_end.clamp(hsp.query_offset, q_len);
            hsp.subject_offset = hsp.subject_offset.clamp(0, s_len);
            hsp.subject_end = hsp.subject_end.clamp(hsp.subject_offset, s_len);

            let a = walk.contains_hsp(&hsp, context_offset, separation);
            let b = indexed.contains_hsp(&hsp, context_offset, separation);
            assert!(a == b, "trial {trial}: walk={a} indexed={b} for {hsp:?}");
            queries += 1;
            contained += a as usize;
            if !a || rng.below(4) == 0 {
                walk.add_hsp(hsp, context_offset, IndexMethod::QueryAndSubject);
                indexed.add_hsp(hsp, context_offset, IndexMethod::QueryAndSubject);
                added.push((hsp, context_offset));
            }
        }
        assert_eq!(walk.node_count(), indexed.node_count(), "trial {trial}");
        let dead = |tree: &BlastIntervalTree| tree.nodes.iter().filter(|n| n.x_dead).count();
        assert_eq!(dead(&walk), dead(&indexed), "trial {trial}");
        assert!(indexed.x.usable, "trial {trial}: index abandoned");
        unlinked += dead(&walk);
        // every HSP ever offered, asked again of the final trees
        for &(hsp, context_offset) in &added {
            let a = walk.contains_hsp(&hsp, context_offset, separation);
            let b = indexed.contains_hsp(&hsp, context_offset, separation);
            assert!(a == b, "trial {trial}: final walk={a} indexed={b}");
        }
    }
    // the generator must exercise both answers and the unlinking
    assert!(
        contained * 10 > queries,
        "{contained} of {queries} contained"
    );
    assert!(
        (queries - contained) * 10 > queries,
        "{contained} of {queries} contained"
    );
    assert!(
        unlinked > trials,
        "{unlinked} HSPs unlinked in {trials} trials"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:102-109,603-793 (cited at
// `try_add_hsp`)
// Ranges in order inside the tree's bounds never make an inverted node: `try_add_hsp` sends
// every one of them to `add_hsp` (`x_insertion_ends`) and builds the tree `add_hsp` builds. Two
// ranges past `s_max` in one subject subtree make one (NCBI's endless split); one alone is
// added. The walk itself is tested in `x_insert_inverts_agrees_with_ncbi_insertion_loop`.
#[test]
fn try_add_hsp_refuses_only_the_endless_split() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = |n: i32| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as i32
    };
    for round in 0..200 {
        let (q_max, s_max) = (20 + next(200), 20 + next(300));
        let mut plain = BlastIntervalTree::new(0, q_max + 1, 0, s_max + 1);
        let mut tried = BlastIntervalTree::new(0, q_max + 1, 0, s_max + 1);
        for _ in 0..(1 + next(60)) {
            let (q0, s0) = (next(q_max), next(s_max));
            let (q1, s1) = (q0 + next(q_max - q0 + 1), s0 + next(s_max - s0 + 1));
            let hsp = make_tree_hsp(q0, q1, s0, s1, 1 + next(50), 1);
            plain.add_hsp(hsp, 0, IndexMethod::QueryAndSubject);
            assert!(tried.try_add_hsp(hsp, 0), "round {round}");
        }
        assert_eq!(plain.node_count(), tried.node_count(), "round {round}");
    }
    // Query 0..20 (middle 10), subject 0..31: [32, 33] and [32, 32] past `s_max` 31.
    let mut tree = BlastIntervalTree::new(0, 21, 0, 31);
    assert!(tree.try_add_hsp(make_tree_hsp(9, 10, 32, 33, 22, 1), 0));
    let nodes = tree.node_count();
    assert!(!tree.try_add_hsp(make_tree_hsp(10, 10, 32, 32, 11, 1), 0));
    assert_eq!(tree.node_count(), nodes, "the tree is unchanged");
}

/// What NCBI's `BlastIntervalTreeAddHSP` (eQueryAndSubject) does with `hsp`: the end-point
/// checks drop it, or the insertion loop ends (the tree after it), or the loop does not end
/// (more splits than any tree of these sizes can hold).
enum NcbiInsertion {
    Dropped,
    Added(Box<BlastIntervalTree>),
    Endless,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:577-584,603-793 (quoted
// at `add_hsp`, `add_hsp_internal`, `split_and_insert`)
/// NCBI's insertion on a copy of `tree`, step by step with the tree's own pieces (the
/// end-point checks, the leaf, `split_and_insert`, the subject root), with a cap on the
/// splits: the reference that `x_insert_inverts` and `try_add_hsp` are tested against.
/// Also returns the tree after the end-point checks (the one the insertion walks).
fn ncbi_insertion_capped(
    tree: &BlastIntervalTree,
    hsp: TreeHsp,
    query_start: i32,
) -> (NcbiInsertion, BlastIntervalTree) {
    let mut t = tree.clone();
    let x_on = t.x.usable;
    if t.x_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Left, x_on)
        || t.x_has_hsp_endpoint(&hsp, query_start, IntervalDirection::Right, x_on)
    {
        return (NcbiInsertion::Dropped, t);
    }
    let checked = t.clone();
    let new_index = t.alloc_leaf_node(hsp, query_start);
    let (mut start, mut end) = (query_start + hsp.query_offset, query_start + hsp.query_end);
    let mut subject = false;
    let mut root = 0usize;
    for _ in 0..10_000 {
        let (l, r) = (t.nodes[root].leftend, t.nodes[root].rightend);
        let middle = ((i64::from(l) + i64::from(r)) / 2) as i32;
        let dir = if end < middle {
            IntervalDirection::Left
        } else if start > middle {
            IntervalDirection::Right
        } else {
            if subject {
                t.nodes[new_index].midptr = t.nodes[root].midptr;
                t.nodes[root].midptr = new_index as i32;
                return (NcbiInsertion::Added(Box::new(t)), checked);
            }
            subject = true;
            if t.nodes[root].midptr == 0 {
                let mid = t.alloc_subject_root_node();
                t.nodes[root].midptr = mid as i32;
            }
            root = t.nodes[root].midptr as usize;
            (start, end) = (hsp.subject_offset, hsp.subject_end);
            continue;
        };
        let child = if dir == IntervalDirection::Left {
            t.nodes[root].leftptr
        } else {
            t.nodes[root].rightptr
        };
        if child == 0 {
            if dir == IntervalDirection::Left {
                t.nodes[root].leftptr = new_index as i32;
            } else {
                t.nodes[root].rightptr = new_index as i32;
            }
            return (NcbiInsertion::Added(Box::new(t)), checked);
        }
        if !t.nodes[child as usize].is_leaf() {
            root = child as usize;
            continue;
        }
        root = t.split_and_insert(
            root,
            child as usize,
            dir,
            subject,
            IndexMethod::QueryAndSubject,
        );
    }
    (NcbiInsertion::Endless, checked)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:603-793 (cited at
// `try_add_hsp`)
// Review 3 of RG X2 (M3): the walk compared the half of a new subject root with the query
// node it came from and called NCBI's ordinary split endless when the two were equal (query
// root [0, 30] and subject root [0, 30]: both HSPs in the query node [0, 15] across its middle
// 7, both left of the subject middle 15; `tblastn -window_size 0` stopped where NCBI reports
// two HSPs). Review 4 (M4): a range whose start is one past its end (a gapped alignment of
// score 0) never straddles a middle, so two of them inside the bounds split a node `[p, p]`
// without end too ([11, 10] x [1, 0], then [11, 10] x [5, 4]). Against NCBI's loop run on the
// tree's own pieces: random trees whose query and subject roots are equal, halves or doubles
// of each other, or unrelated, with query context offsets, subject ranges inside the bounds
// and past `s_max`, and ranges of no letters; the walk, `try_add_hsp`'s answer and the nodes
// of the tree it builds agree with NCBI on every HSP. An endless case is checked to be refused
// before `try_add_hsp` runs, so that a slip fails here instead of filling the memory.
#[test]
fn x_insert_inverts_agrees_with_ncbi_insertion_loop() {
    let mut tree = BlastIntervalTree::new(0, 30, 0, 30);
    assert!(tree.try_add_hsp(make_tree_hsp(4, 11, 0, 7, 20, 1), 0));
    let second = make_tree_hsp(4, 11, 7, 14, 20, 1);
    assert!(!tree.x_insert_inverts(&second, 0));
    assert!(tree.try_add_hsp(second, 0));
    let mut tree = BlastIntervalTree::new(0, 30, 0, 30);
    assert!(tree.try_add_hsp(make_tree_hsp(11, 10, 1, 0, 0, 1), 0));
    let second = make_tree_hsp(11, 10, 5, 4, 0, 1);
    assert!(!tree.x_insertion_ends(&second, 0) && tree.x_insert_inverts(&second, 0));
    let nodes = tree.node_count();
    assert!(!tree.try_add_hsp(second, 0));
    assert_eq!(tree.node_count(), nodes, "the tree is unchanged");

    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = |n: i32| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as i32
    };
    let (mut endless, mut endless_inverted, mut added) = (0, 0, 0);
    for round in 0..3000 {
        let q_root = 8 + next(120);
        let s_root = match next(5) {
            0 => q_root,
            1 => q_root / 2,
            2 => q_root * 2,
            3 => q_root + next(3) - 1,
            _ => 8 + next(240),
        };
        let mut tried = BlastIntervalTree::new(0, q_root, 0, s_root);
        // The query's context offset in the batch (TBLASTN's queries after the first).
        let qoff = if next(2) == 0 { 0 } else { next(q_root / 2) };
        let q_len = q_root - qoff;
        for step in 0..(1 + next(40)) {
            let q0 = next(q_len + 1);
            let q1 = q0 + next(q_len - q0 + 1);
            // A third of the subject ranges are short (up to 4 letters), and those that
            // start near `s_max` end past it (TBLASTN's one-hit starts past their chunk).
            let s0 = next(s_root + 3);
            let s1 = if s0 > s_root || next(3) == 0 {
                s0 + next(4)
            } else {
                s0 + next(s_root - s0 + 1)
            };
            // One range in eight has no letters: start one past its end in both sequences
            // (a gapped alignment of score 0 whose right extension did not run).
            let inverted = next(8) == 0;
            let hsp = if inverted {
                make_tree_hsp(q0 + 1, q0, s0 + 1, s0, next(3), 1)
            } else {
                make_tree_hsp(q0, q1, s0, s1, 1 + next(50), 1)
            };
            let (ncbi, checked) = ncbi_insertion_capped(&tried, hsp, qoff);
            let at = format!("round {round} step {step}: q {q_root} +{qoff} s {s_root} {hsp:?}");
            if !matches!(ncbi, NcbiInsertion::Dropped) {
                assert_eq!(
                    checked.x_insert_inverts(&hsp, qoff),
                    matches!(ncbi, NcbiInsertion::Endless),
                    "walk, {at}"
                );
            }
            if matches!(ncbi, NcbiInsertion::Endless) {
                assert!(
                    !tried.x_insertion_ends(&hsp, qoff) && tried.x_insert_inverts(&hsp, qoff),
                    "refused before the insertion, {at}"
                );
            }
            let nodes = tried.node_count();
            let ok = tried.try_add_hsp(hsp, qoff);
            match ncbi {
                NcbiInsertion::Endless => {
                    endless += 1;
                    endless_inverted += usize::from(inverted);
                    assert!(!ok, "{at}");
                    assert_eq!(tried.node_count(), nodes, "unchanged, {at}");
                }
                NcbiInsertion::Added(t) => {
                    added += 1;
                    assert!(ok, "{at}");
                    assert_eq!(
                        format!("{:?}", tried.nodes),
                        format!("{:?}", t.nodes),
                        "{at}"
                    );
                }
                // The end-point checks dropped the HSP (and maybe leaves of lower scores).
                NcbiInsertion::Dropped => {
                    assert!(ok, "{at}");
                    let (t, nt) = (format!("{:?}", tried.nodes), format!("{:?}", checked.nodes));
                    assert_eq!(t, nt, "{at}");
                }
            }
        }
    }
    assert!(
        endless > 100 && endless_inverted > 20 && added > 10_000,
        "endless {endless} (no letters {endless_inverted}) added {added}"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:603-793 (cited at
// `try_add_hsp`) and c++/src/algo/blast/core/blast_gapalign.c:3830-3831
// RPS-BLAST's tree (line RG D1a-4): the root `[0, Lc + 1]` of a concatenated length `Lc` that
// an alias file's `LENGTH` and `NSEQ` make shorter than the profiles, also negative (a root
// whose end is below its start), with query ranges in order that start and end past the root's
// end or cross it, in profiles at several offsets, and subject ranges inside the bounds. A
// second range that starts past the end splits without end (a root end `R` of 0 or more:
// `[R + 1, R]`; below 0 every split), unless NCBI's end-point checks drop it or the leaf it
// meets; `try_add_hsp` answers as NCBI's loop and builds its trees.
#[test]
fn try_add_hsp_agrees_with_ncbi_on_query_ranges_past_the_root() {
    let mut state = 0x6a09_e667_f3bc_c908u64;
    let mut next = |n: i32| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as i32
    };
    // Root [0, 21] (`Lc` 20): [30, 40] goes to the right spine's end; [25, 26] meets it.
    let mut tree = BlastIntervalTree::new(0, 21, 0, 101);
    assert!(tree.try_add_hsp(make_tree_hsp(30, 40, 10, 20, 30, 1), 0));
    // One range crossing the end straddles a middle of the spine and is added.
    assert!(tree.try_add_hsp(make_tree_hsp(15, 35, 50, 70, 25, 1), 0));
    let nodes = tree.node_count();
    assert!(!tree.try_add_hsp(make_tree_hsp(25, 26, 80, 81, 20, 1), 0));
    assert_eq!(tree.node_count(), nodes, "the tree is unchanged");
    // Root [0, -1] (`Lc` -2): a range at 0 straddles the middle 0, one at 1 goes right.
    let mut tree = BlastIntervalTree::new(0, -1, 0, 101);
    assert!(tree.try_add_hsp(make_tree_hsp(0, 3, 10, 13, 30, 1), 0));
    assert!(tree.try_add_hsp(make_tree_hsp(1, 4, 20, 23, 30, 1), 0));
    assert!(!tree.try_add_hsp(make_tree_hsp(5, 9, 40, 44, 30, 1), 0));

    let (mut endless, mut added, mut dropped) = (0, 0, 0);
    for round in 0..3000 {
        // The root's end `Lc + 1`, from -12 to 47, and the subject's.
        let q_root = next(60) - 12;
        let s_root = 8 + next(120);
        let mut tried = BlastIntervalTree::new(0, q_root, 0, s_root);
        let mut last: Option<TreeHsp> = None;
        for step in 0..(1 + next(30)) {
            // A profile at its offset in the concatenation, a range inside it; one in four
            // starts where the last one did (NCBI's end-point checks drop it, or unlink the
            // leaf of the lower score).
            let offset = next(4) * 25;
            let (offset, q0, s0) = match last {
                Some(last) if next(4) == 0 => (
                    last.query_context_offset,
                    last.query_offset,
                    last.subject_offset,
                ),
                _ => (offset, next(24), next(s_root + 1)),
            };
            let q1 = q0 + next(24 - q0);
            let s1 = s0 + next(s_root - s0 + 1);
            let mut hsp = make_tree_hsp(q0, q1, s0, s1, 1 + next(50), 1);
            hsp.query_context_offset = offset;
            last = Some(hsp);
            let (ncbi, checked) = ncbi_insertion_capped(&tried, hsp, offset);
            let at = format!("round {round} step {step}: root {q_root} +{offset} {hsp:?}");
            let nodes = tried.node_count();
            let ok = tried.try_add_hsp(hsp, offset);
            match ncbi {
                NcbiInsertion::Endless => {
                    endless += 1;
                    assert!(offset + q0 > q_root, "a range past the end, {at}");
                    assert!(!ok, "{at}");
                    assert_eq!(tried.node_count(), nodes, "unchanged, {at}");
                }
                NcbiInsertion::Added(t) => {
                    added += 1;
                    assert!(ok, "{at}");
                    assert_eq!(
                        format!("{:?}", tried.nodes),
                        format!("{:?}", t.nodes),
                        "{at}"
                    );
                }
                NcbiInsertion::Dropped => {
                    dropped += 1;
                    assert!(ok, "{at}");
                    assert_eq!(
                        format!("{:?}", tried.nodes),
                        format!("{:?}", checked.nodes),
                        "{at}"
                    );
                }
            }
            if !ok {
                break;
            }
        }
    }
    assert!(
        endless > 500 && added > 2_000 && dropped > 50,
        "endless {endless} added {added} dropped {dropped}"
    );
}

#[test]
fn test_interval_tree_basic() {
    let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);

    // Add an HSP
    let hsp1 = TreeHsp {
        query_offset: 100,
        query_end: 200,
        subject_offset: 500,
        subject_end: 600,
        score: 100,
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign: 1,
    };
    tree.add_hsp(hsp1, 0, IndexMethod::QueryAndSubject);

    // Check containment of a smaller HSP within the first
    let hsp2 = TreeHsp {
        query_offset: 120,
        query_end: 180,
        subject_offset: 520,
        subject_end: 580,
        score: 50,
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign: 1,
    };

    // hsp2 is contained within hsp1 (lower score, both endpoints inside)
    assert!(tree.contains_hsp(&hsp2, 0, 0));
}

#[test]
fn test_interval_tree_not_contained() {
    let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);

    let hsp1 = TreeHsp {
        query_offset: 100,
        query_end: 200,
        subject_offset: 500,
        subject_end: 600,
        score: 100,
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign: 1,
    };
    tree.add_hsp(hsp1, 0, IndexMethod::QueryAndSubject);

    // HSP with higher score should not be contained
    let hsp2 = TreeHsp {
        query_offset: 120,
        query_end: 180,
        subject_offset: 520,
        subject_end: 580,
        score: 150, // Higher score
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign: 1,
    };

    assert!(!tree.contains_hsp(&hsp2, 0, 0));
}

#[test]
fn test_interval_tree_different_strand() {
    let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);

    let hsp1 = TreeHsp {
        query_offset: 100,
        query_end: 200,
        subject_offset: 500,
        subject_end: 600,
        score: 100,
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign: 1, // Forward strand
    };
    tree.add_hsp(hsp1, 0, IndexMethod::QueryAndSubject);

    // HSP on different strand should not be contained
    let hsp2 = TreeHsp {
        query_offset: 120,
        query_end: 180,
        subject_offset: 520,
        subject_end: 580,
        score: 50,
        query_frame: 1,
        query_length: 1000,
        query_context_offset: 0,
        subject_frame_sign: -1, // Reverse strand
    };

    assert!(!tree.contains_hsp(&hsp2, 0, 0));
}

#[test]
fn test_mb_hsp_close() {
    // Same diagonal - should be close
    assert!(BlastIntervalTree::mb_hsp_close(100, 200, 150, 250, 50));

    // Different diagonals - should not be close
    assert!(!BlastIntervalTree::mb_hsp_close(100, 200, 150, 300, 50));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:269-301
#[test]
fn test_interval_tree_common_endpoint_prefers_shorter_equal_score_hsp() {
    let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);
    let longer_hsp = make_tree_hsp(100, 200, 500, 600, 100, 1);
    let shorter_hsp = make_tree_hsp(100, 180, 500, 580, 100, 1);
    let long_only_candidate = make_tree_hsp(181, 190, 581, 590, 50, 1);

    tree.add_hsp(longer_hsp, 0, IndexMethod::QueryAndSubject);
    tree.add_hsp(shorter_hsp, 0, IndexMethod::QueryAndSubject);

    assert!(!tree.contains_hsp(&long_only_candidate, 0, 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819-823
#[test]
fn test_interval_tree_containment_uses_subject_frame_sign() {
    let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);
    let tree_hsp = make_tree_hsp(100, 200, 500, 600, 100, 3);
    let contained_hsp = make_tree_hsp(120, 180, 520, 580, 50, 1);

    tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);

    assert!(tree.contains_hsp(&contained_hsp, 0, 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819-847
#[test]
fn test_containing_hsp_reports_ncbi_score_and_endpoint_semantics() {
    let mut tree = BlastIntervalTree::new(0, 1000, 0, 2000);
    let tree_hsp = make_tree_hsp(100, 200, 500, 600, 100, 1);
    tree.add_hsp(tree_hsp, 0, IndexMethod::QueryAndSubject);

    let lower_score_inside = make_tree_hsp(120, 180, 520, 580, 50, 1);
    let equal_score_inside = make_tree_hsp(125, 175, 525, 575, 100, 1);
    let higher_score_inside = make_tree_hsp(120, 180, 520, 580, 101, 1);
    let endpoint_outside = make_tree_hsp(120, 220, 520, 620, 50, 1);
    let opposite_subject_sign = make_tree_hsp(120, 180, 520, 580, 50, -1);

    let container = tree
        .containing_hsp(&lower_score_inside, 0, 0)
        .expect("lower-score HSP with both endpoints inside must be contained");
    assert_eq!(container.query_offset, tree_hsp.query_offset);
    assert_eq!(container.subject_offset, tree_hsp.subject_offset);
    assert!(tree.contains_hsp(&equal_score_inside, 0, 0));
    assert!(!tree.contains_hsp(&higher_score_inside, 0, 0));
    assert!(!tree.contains_hsp(&endpoint_outside, 0, 0));
    assert!(!tree.contains_hsp(&opposite_subject_sign, 0, 0));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:819-843
#[test]
fn test_containment_respects_query_context_and_diag_separation() {
    let mut tree = BlastIntervalTree::new(0, 3000, 0, 2000);
    let mut tree_hsp = make_tree_hsp(100, 200, 500, 600, 100, 1);
    tree_hsp.query_context_offset = 1001;
    tree.add_hsp(tree_hsp, 1001, IndexMethod::QueryAndSubject);

    let same_context_far_diag = make_tree_hsp(150, 160, 500, 510, 50, 1);
    let same_context_close_diag = make_tree_hsp(149, 159, 500, 510, 50, 1);
    let different_context = make_tree_hsp(149, 159, 500, 510, 50, 1);

    assert!(tree.contains_hsp(&same_context_far_diag, 1001, 0));
    assert!(!tree.contains_hsp(&same_context_far_diag, 1001, 50));
    assert!(tree.contains_hsp(&same_context_close_diag, 1001, 50));
    assert!(!tree.contains_hsp(&different_context, 0, 0));
}
