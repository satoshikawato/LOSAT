use super::super::super::{advance_diag_table_offset, x_diag_lv::LvSeeds};
use super::super::copy_hash;
use super::*;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

type HashState = (
    Vec<u32>,
    Vec<(i32, i32, bool, i32, u32)>,
    u32,
    u32,
    i32,
    Vec<i32>,
    Vec<super::super::super::x_diag_lv::LvLast>,
);

/// Everything of a table a later seed can read: the heads, the cells up to
/// the occupancy, the occupancy, the capacity, the offset, the level bounds and the memo
/// (LOSAT_X_DIAGLV2).
fn hash_state(h: &DiagHashTable) -> HashState {
    assert_eq!(
        h.chain.len(),
        h.capacity as usize,
        "capacity and chain length"
    );
    (
        h.backbone.clone(),
        h.chain[..h.occupancy as usize]
            .iter()
            .map(|c| (c.diag, c.level, c.hit_saved, c.hit_len, c.next))
            .collect(),
        h.occupancy,
        h.capacity,
        h.offset,
        h.max_level.clone(),
        h.x_last.clone(),
    )
}

// A seed outcome of the one-hit fold, as `fold` takes it: the retrieve, the `last_hit`
// test, then an insert of a seed that is extended below the cutoff (`hit_len` =
// `s_end_pos - s_off_pos`, not saved) or above it (a hit; the ungapped end, saved).
#[derive(Clone, Copy, Debug)]
struct Op {
    diag: i32,
    s_off_pos: i32,
    s_end_pos: i32,
    pass: bool,
    ungapped_end: i32,
}

/// The seed of op `index`: its diagonal as `s - q`.
fn seed_of(index: usize, op: &Op) -> Seed {
    const Q0: i64 = 1 << 20;
    Seed {
        context: index as u32,
        q: Q0 as u32,
        s: (Q0 + i64::from(op.diag)) as u32,
        s_range: 0,
    }
}

fn last_hit<const MEMO: bool>(hash: &DiagHashTable, op: &Op) -> i32 {
    let found = if hash.x_lv_on() {
        hash.retrieve_below::<MEMO>(op.diag, &LvSeeds::new(false, 11, 0), op.s_off_pos as usize)
    } else {
        hash.retrieve(op.diag)
    };
    found.map_or(0, |(level, _, _)| level)
}

/// The reference: op after op on one table; the hits are op indices. With `MEMO` it reads
/// and keeps the table's memo (LOSAT_X_DIAGLV2) as the folds of a strips stage do, without it
/// neither, as the reference fold.
fn reference<const MEMO: bool>(
    hash: &mut DiagHashTable,
    ops: &[Op],
    first: usize,
    window: i32,
) -> Vec<usize> {
    let insert = if MEMO {
        DiagHashTable::insert_memo
    } else {
        DiagHashTable::insert_plain
    };
    let mut hits = Vec::new();
    for (k, op) in ops.iter().enumerate() {
        if op.s_off_pos < last_hit::<MEMO>(hash, op) {
            continue;
        }
        if op.pass {
            hits.push(first + k);
            insert(
                hash,
                op.diag,
                op.ungapped_end,
                0,
                true,
                op.s_off_pos,
                window,
            );
        } else {
            insert(
                hash,
                op.diag,
                op.s_end_pos,
                op.s_end_pos - op.s_off_pos,
                false,
                op.s_off_pos,
                window,
            );
        }
    }
    hits
}

fn hit_of(index: usize) -> UngappedHit {
    UngappedHit {
        context_idx: 0,
        query_idx: 0,
        query_frame: 1,
        query_context_offset: 0,
        seed_q_off: index as u32,
        seed_s_off: 0,
        qs: 0,
        qe: 0,
        ss: 0,
        se: 0,
        score: 0,
    }
}

// LOSAT_X_PAIRPARMERGE: `hits_by_strips_on_pool` on pools of 2, 3 and 8 threads, cut into 1,
// 2, 7, one per strip and more ranges than strips, and `hits_by_strips_here` off a pool (its
// `runs` left over from the stage before), put the hits after those already in the list in
// the order of `merge`'s sort of all of them (the same states): random stages of 1 to 9
// partitions whose keys ascend, 1 to 300 strips, strips without hits, an empty partition.
#[test]
fn hits_by_strips_is_the_sort_of_merge() {
    let mut rng = Lcg(0xA15_5EED);
    let mut runs: Vec<(u64, u32, u32)> = Vec::new();
    let pools: Vec<rayon::ThreadPool> = [2, 3, 8]
        .iter()
        .map(|&t| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build()
                .unwrap()
        })
        .collect();
    for round in 0..48u64 {
        let parts = 1 + (rng.next() % 9) as usize;
        let strips = 1 + (rng.next() % 300) as usize;
        let mut states: Vec<PartState> = (0..parts).map(|_| PartState::new()).collect();
        let mut index = 0usize;
        for strip in 0..strips {
            // a strip of `seeds` seeds in scan order, a third of them without a hit
            let seeds = rng.next() % 40;
            for pos in 0..seeds {
                if rng.next().is_multiple_of(3) {
                    continue;
                }
                let mut p = (rng.next() % parts as u64) as usize;
                if parts >= 3 && round % 2 == 0 && p == parts - 1 {
                    p = 0;
                }
                states[p].out.hit_keys.push(((strip as u64) << 32) | pos);
                states[p].hits.push(hit_of(index));
                index += 1;
            }
        }
        let before: Vec<UngappedHit> = (0..(round % 3) as usize)
            .map(|i| hit_of(1_000_000 + i))
            .collect();
        // `merge`'s sort
        let mut order: Vec<(u64, u32, u32)> = Vec::new();
        for (p, st) in states.iter().enumerate() {
            order.extend(
                st.out
                    .hit_keys
                    .iter()
                    .enumerate()
                    .map(|(j, &key)| (key, p as u32, j as u32)),
            );
        }
        order.sort_unstable();
        let mut want: Vec<usize> = before.iter().map(|h| h.seed_q_off as usize).collect();
        want.extend(
            order
                .iter()
                .map(|&(_, p, j)| states[p as usize].hits[j as usize].seed_q_off as usize),
        );
        for pool in &pools {
            for ranges in [1, 2, 7, strips, 4 * strips] {
                let got = pool.install(|| {
                    let mut hits = before.clone();
                    hits_by_strips_on_pool(&mut hits, &states, strips, ranges);
                    hits
                });
                let got: Vec<usize> = got.iter().map(|h| h.seed_q_off as usize).collect();
                assert_eq!(
                    got,
                    want,
                    "round {round}: {parts} partitions, {strips} strips, {ranges} ranges, {} threads",
                    pool.current_num_threads()
                );
            }
        }
        let mut hits = before.clone();
        hits_by_strips_here(&mut hits, &states, strips, &mut runs);
        let got: Vec<usize> = hits.iter().map(|h| h.seed_q_off as usize).collect();
        assert_eq!(
            got, want,
            "round {round}: {parts} partitions, {strips} strips, this thread"
        );
    }
}

/// The same ops as a partitioned stage: strips cut at `cuts`, `parts` partitions on
/// `threads` threads (in a pool of that many). `nested` runs a parallel iterator inside
/// every scan and fold call. With `owner_waits` a fold on the owner's thread first waits
/// until a helper has folded in this stage (100 us at a time, 200 ms at most in the stage),
/// so that the helpers get to fold on any machine (even on one CPU, or on a loaded one,
/// where the owner could otherwise fold the whole stage before a helper job starts). `home`: LOSAT_X_PAIRPARHOME's way of sharing the folds.
/// Returns the hits (op indices), the stage's counters and the number of threads that
/// folded.
#[allow(clippy::too_many_arguments)]
fn partitioned(
    pool: &rayon::ThreadPool,
    hash: &mut DiagHashTable,
    scratch: &mut PartScratch,
    ops: &[Op],
    first: usize,
    cuts: &[usize],
    parts: usize,
    ahead: usize,
    window: i32,
    nested: bool,
    owner_waits: bool,
    home: bool,
) -> (Vec<usize>, PartStats, usize) {
    let threads = pool.current_num_threads();
    partitioned_on(
        pool,
        threads,
        None,
        hash,
        scratch,
        ops,
        first,
        cuts,
        parts,
        ahead,
        window,
        nested,
        owner_waits,
        home,
    )
}

/// How `partitioned_on` slows a participant down.
#[derive(Clone, Copy, Debug)]
enum Slow {
    /// every scan call on this pool thread sleeps 1 ms first
    Scans(usize),
    /// the first scan on the owner's thread holds the stage (`hold_until_helpers_leave`)
    OwnerFirstScan,
}

/// (`Slow::OwnerFirstScan`, `Forced`) Sleeps in 1 ms steps until `done`, or until 60 s have
/// passed since `from` (then the test's assertions fail: a stall exit that never fires, or a
/// participant that never comes, fails the test instead of hanging it).
fn wait_until(from: std::time::Instant, done: impl Fn() -> bool) {
    while !done() && from.elapsed() < std::time::Duration::from_secs(60) {
        std::thread::sleep(std::time::Duration::from_millis(1));
    }
}

/// `Slow::OwnerFirstScan`, on the owner's thread in its first scan (no partition can fold
/// past that strip, so the stage makes no progress once the helpers have scanned the strips
/// the window lets them): waits until a helper has started a scan (it is in the stage; a
/// helper's scans wait for the owner's first one, so that the helpers cannot do the whole
/// stage before the owner holds a strip), then until every other thread of the pool has run
/// a job of `spawn_broadcast`. A helper's thread runs that job only once its participant has
/// returned (the participants of these stages make no rayon calls, the test's scans and
/// folds none without `nested`), and while the owner holds the strip a helper returns only
/// through the stall exit, so the hold ends after every helper that came has left; how soon
/// the OS runs the helpers does not matter. The owner's own job runs after the hold.
fn hold_until_helpers_leave(pool: &rayon::ThreadPool, helper_scanned: &AtomicBool) {
    let from = std::time::Instant::now();
    wait_until(from, || helper_scanned.load(Ordering::Acquire));
    let owner = rayon::current_thread_index();
    let others = Arc::new(AtomicUsize::new(0));
    let ran = Arc::clone(&others);
    pool.spawn_broadcast(move |ctx| {
        if Some(ctx.index()) != owner {
            ran.fetch_add(1, Ordering::AcqRel);
        }
    });
    wait_until(from, || {
        others.load(Ordering::Acquire) + 1 >= pool.current_num_threads()
    });
}

/// `partitioned` with the stage's thread count given (`threads` participants whatever the
/// pool's size) and a participant slowed down as `slow` says.
#[allow(clippy::too_many_arguments)]
fn partitioned_on(
    pool: &rayon::ThreadPool,
    threads: usize,
    slow: Option<Slow>,
    hash: &mut DiagHashTable,
    scratch: &mut PartScratch,
    ops: &[Op],
    first: usize,
    cuts: &[usize],
    parts: usize,
    ahead: usize,
    window: i32,
    nested: bool,
    owner_waits: bool,
    home: bool,
) -> (Vec<usize>, PartStats, usize) {
    partitioned_with(
        pool,
        threads,
        slow,
        hash,
        scratch,
        ops,
        first,
        cuts,
        parts,
        ahead,
        window,
        nested,
        owner_waits,
        home,
        // (a test that asks for pieces asks for the lists' stage they replace)
        scratch.x_pieces.is_some_and(|k| k > 0),
    )
}

/// `partitioned_on`; with `lists` the stage of LOSAT_X_DCSEEDLEAN (`seed_stage_part_lists`):
/// the folds do not record their buckets (`insert_t` without `TRACK`) and the home stage
/// scans the strips into the partitions' lists (`PartSink`), or with `home` and
/// `scratch.x_pieces` `Some(k)` (k > 0) into pieces of `k` seeds (LOSAT_X_SEEDPIECE,
/// `piece::seed_stage_pieces`; the scan stops at the first seed a piece has no room for and
/// goes on from it).
#[allow(clippy::too_many_arguments)]
fn partitioned_with(
    pool: &rayon::ThreadPool,
    threads: usize,
    slow: Option<Slow>,
    hash: &mut DiagHashTable,
    scratch: &mut PartScratch,
    ops: &[Op],
    first: usize,
    cuts: &[usize],
    parts: usize,
    ahead: usize,
    window: i32,
    nested: bool,
    owner_waits: bool,
    home: bool,
    lists: bool,
) -> (Vec<usize>, PartStats, usize) {
    let pieces = if lists && home {
        scratch.x_pieces.filter(|&k| k > 0)
    } else {
        None
    };
    partitioned_in(
        pool,
        threads,
        slow,
        hash,
        scratch,
        ops,
        first,
        cuts,
        parts,
        ahead,
        window,
        nested,
        owner_waits,
        home,
        lists,
        pieces,
    )
}

/// `partitioned_with` with the pieces given (`None`: the strips' lists).
#[allow(clippy::too_many_arguments)]
fn partitioned_in(
    pool: &rayon::ThreadPool,
    threads: usize,
    slow: Option<Slow>,
    hash: &mut DiagHashTable,
    scratch: &mut PartScratch,
    ops: &[Op],
    first: usize,
    cuts: &[usize],
    parts: usize,
    ahead: usize,
    window: i32,
    nested: bool,
    owner_waits: bool,
    home: bool,
    lists: bool,
    pieces: Option<usize>,
) -> (Vec<usize>, PartStats, usize) {
    let folders = Mutex::new(std::collections::BTreeSet::new());
    let helper_folded = AtomicBool::new(false);
    let owner_waited = AtomicUsize::new(0);
    let owner_scanned = AtomicBool::new(false);
    // (`Slow::OwnerFirstScan`: a thread other than the owner's has started a scan)
    let helper_scanned = AtomicBool::new(false);
    let started = std::time::Instant::now();
    pool.install(|| {
        let owner = rayon::current_thread_index();
        let slow_scans = match slow {
            Some(Slow::Scans(t)) => Some(t),
            _ => None,
        };
        let before_scan = || {
            if slow_scans.is_some() && rayon::current_thread_index() == slow_scans {
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
            if matches!(slow, Some(Slow::OwnerFirstScan)) {
                if rayon::current_thread_index() != owner {
                    wait_until(started, || owner_scanned.load(Ordering::Acquire));
                    helper_scanned.store(true, Ordering::Release);
                } else if !owner_scanned.swap(true, Ordering::AcqRel) {
                    hold_until_helpers_leave(pool, &helper_scanned);
                }
            }
            if nested {
                use rayon::prelude::*;
                let sum: u64 = (0..64u64).into_par_iter().map(|v| v * v).sum();
                assert_eq!(sum, 85344);
            }
        };
        let scan = |i: usize, out: &mut Vec<Seed>| {
            before_scan();
            for k in cuts[i]..cuts[i + 1] {
                out.push(seed_of(first + k, &ops[k]));
            }
        };
        let scan_lists = |i: usize, out: &mut PartSink<'_>| {
            before_scan();
            for k in cuts[i]..cuts[i + 1] {
                out.push(seed_of(first + k, &ops[k]));
            }
        };
        // (a piece's scan: `Resume::s` is the seed's index in its strip)
        let scan_piece = |i: usize, from: Option<Resume>, out: &mut PieceSink<'_>| {
            before_scan();
            let start = cuts[i] + from.map_or(0, |r| r.s);
            for k in start..cuts[i + 1] {
                if out.room() == 0 {
                    return Some(Resume::at(k - cuts[i]));
                }
                out.push(seed_of(first + k, &ops[k]));
            }
            None
        };
        #[allow(clippy::too_many_arguments)]
        fn put(
            lists: bool,
            hash: &mut DiagHashTable,
            out: &mut PartOut,
            key: u64,
            diag: i32,
            level: i32,
            hit_len: i32,
            hit_saved: bool,
            s_off: i32,
            window: i32,
        ) {
            if lists {
                insert_t::<true, false, true, false>(
                    hash, out, key, diag, level, hit_len, hit_saved, s_off, window,
                );
            } else {
                insert::<true, false, true>(
                    hash, out, key, diag, level, hit_len, hit_saved, s_off, window,
                );
            }
        }
        let fold = |st: &mut PartState, seeds: &[Seed], pos: &[u32]| {
            if nested {
                use rayon::prelude::*;
                let sum: u64 = (0..64u64).into_par_iter().map(|v| v + 1).sum();
                assert_eq!(sum, 2080);
            }
            let here = rayon::current_thread_index();
            lock(&folders).insert(here);
            if owner_waits {
                if here != owner {
                    helper_folded.store(true, Ordering::Relaxed);
                } else {
                    // (200 ms at most in a stage: on a loaded machine a helper may start late)
                    while !helper_folded.load(Ordering::Relaxed)
                        && owner_waited.load(Ordering::Relaxed) < 2000
                    {
                        std::thread::sleep(std::time::Duration::from_micros(100));
                        owner_waited.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
            for (k, seed) in seeds.iter().enumerate() {
                let index = seed.context as usize;
                let op = &ops[index - first];
                let key = st.out.key(pos, k);
                st.counts.seeds += 1;
                if op.s_off_pos < last_hit::<true>(&st.hash, op) {
                    st.counts.skipped += 1;
                    continue;
                }
                if op.pass {
                    st.hits.push(hit_of(index));
                    st.out.hit_keys.push(key);
                    put(
                        lists,
                        &mut st.hash,
                        &mut st.out,
                        key,
                        op.diag,
                        op.ungapped_end,
                        0,
                        true,
                        op.s_off_pos,
                        window,
                    );
                } else {
                    put(
                        lists,
                        &mut st.hash,
                        &mut st.out,
                        key,
                        op.diag,
                        op.s_end_pos,
                        op.s_end_pos - op.s_off_pos,
                        false,
                        op.s_off_pos,
                        window,
                    );
                }
            }
        };
        let mut hits = Vec::new();
        let stats = if let Some(piece_seeds) = pieces {
            rebuild(scratch, hash, parts, true, None);
            let occ0 = hash.occupancy;
            piece::seed_stage_pieces(
                cuts.len() - 1,
                threads,
                ahead.max(1),
                hash,
                &mut hits,
                scratch,
                parts,
                occ0,
                piece_seeds,
                &scan_piece,
                &fold,
                false,
            )
        } else {
            seed_stage_in(
                home,
                cuts.len() - 1,
                threads,
                ahead,
                hash,
                &mut hits,
                FoldPart {
                    scratch,
                    parts,
                    latest: None,
                },
                &scan,
                lists.then_some(&scan_lists),
                &fold,
                false,
            )
        };
        let folders = lock(&folders).len();
        (
            hits.iter().map(|h| h.seed_q_off as usize).collect(),
            stats,
            folders,
        )
    })
}

/// Random ops: diagonals of a few buckets' worth (d + 512 k share a bucket) and many
/// others, positions from `s` on, a few per position.
fn random_ops(rng: &mut Lcg, n: usize, s: &mut i32, diagonals: usize) -> Vec<Op> {
    let base = rng.below(1000) as i32 - 500;
    let diags: Vec<i32> = (0..diagonals)
        .map(|k| {
            if k < 16 {
                base + 512 * k as i32
            } else {
                rng.below(40_000) as i32 - 20_000
            }
        })
        .collect();
    (0..n)
        .map(|_| {
            *s += rng.below(4) as i32;
            let s_end = *s + 11 + rng.below(4) as i32;
            Op {
                diag: diags[rng.below(diags.len() as u64) as usize],
                s_off_pos: *s,
                s_end_pos: s_end,
                pass: rng.below(3) == 0,
                ungapped_end: s_end + rng.below(200) as i32,
            }
        })
        .collect()
}

fn random_cuts(rng: &mut Lcg, n: usize) -> Vec<usize> {
    let mut cuts = vec![0usize];
    while *cuts.last().unwrap() < n {
        // some empty strips
        let step = if rng.below(10) == 0 {
            0
        } else {
            1 + rng.below(80) as usize
        };
        cuts.push((cuts.last().unwrap() + step).min(n));
    }
    cuts
}

fn advance(hash: &mut DiagHashTable, subject_len: usize) {
    let mut diag_table_offset = 0;
    let (mut levels, mut lens) = (Vec::new(), Vec::new());
    let window = hash.window;
    advance_diag_table_offset(
        &mut diag_table_offset,
        window,
        subject_len,
        false,
        &mut levels,
        &mut lens,
        hash,
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
// `seed_stage_part` (and LOSAT_X_DCSEEDLEAN's `seed_stage_part_lists`, alternating with it
// and with LOSAT_X_PAIRPARHOME's two ways) against the reference on random op streams, for 2,
// 3, 7 and 16 partitions on 1..8 threads, with and without level bounds, all stages of a case
// on one scratch (its strip buffers passing from one form to the other): a stage on a new table that grows past 256 cells (many diagonals within the
// insert window), a second stage after the offset moved on (the next subject or chunk), a
// third after the clear of the offset (occupancy 1), a chunk the reference folds, a stage
// after it, a stage of no strips, and a stage with another partition count. Every stage
// gives the reference's hits and table cell for cell (heads, cells, occupancy, capacity,
// bounds). In every other case on 2+ threads the owner sleeps in its folds until a helper
// has folded (`partitioned`, at most 200 ms in a stage), so that these cases check the
// equality on folds by several threads on any machine. That a helper folded is not
// guaranteed for one case (a helper may start after the owner's wait on a loaded machine:
// the review of 2026-10 and three gate runs to 2026-11-11), so the test asserts it for some
// case, not for each.
#[test]
fn seed_stage_part_equals_the_reference_fold() {
    let mut rng = Lcg(0x9A27_F01D);
    let mut case = 0usize;
    let mut waited_cases = 0usize;
    let mut shared_cases = 0usize;
    for parts in [2usize, 3, 7, 16] {
        for threads in 1..=8usize {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for levels in [false, true] {
                case += 1;
                let window = [12, 40, 5000][case % 3];
                // the bounds (LOSAT_X_DIAGLV) and the memo (LOSAT_X_DIAGLV2) as `levels` says,
                // whatever the environment says
                let mut reference_hash = DiagHashTable::new(16);
                reference_hash.max_level = Vec::new();
                reference_hash.x_last = Vec::new();
                if levels {
                    reference_hash.x_lv_start();
                    reference_hash.x_lv2_start();
                }
                let mut hash = copy_hash(&reference_hash);
                let mut scratch = PartScratch::default();
                let mut s = 50i32;
                let mut first = 0usize;
                let mut most_folders = 0usize;
                let mut most_ahead = 0usize;
                // one of the two cases of each (parts, threads), alternating with threads
                let owner_waits = threads >= 2 && (threads + usize::from(levels)).is_multiple_of(2);
                for stage in 0..7 {
                    let at = format!(
                        "case {case} parts {parts} threads {threads} levels {levels} stage {stage}"
                    );
                    match stage {
                        1 | 4 => {
                            advance(&mut reference_hash, 3000);
                            advance(&mut hash, 3000);
                        }
                        2 => {
                            // the clear: offset >= INT4_MAX / 4
                            reference_hash.offset = i32::MAX / 4;
                            hash.offset = i32::MAX / 4;
                            advance(&mut reference_hash, 3000);
                            advance(&mut hash, 3000);
                            assert_eq!(hash.occupancy, 1);
                            s = 50;
                        }
                        _ => {}
                    }
                    let n = if stage == 5 {
                        0
                    } else {
                        500 + rng.below(2500) as usize
                    };
                    let diagonals = if stage == 0 {
                        1200
                    } else {
                        60 + rng.below(600) as usize
                    };
                    let ops = random_ops(&mut rng, n, &mut s, diagonals);
                    // as the search: a strips stage makes the memo (LOSAT_X_DIAGLV2) again at
                    // its start (`seed_stage`) and keeps it; the reference fold (stage 3)
                    // neither reads nor keeps it
                    let expected = if stage == 3 {
                        reference::<false>(&mut reference_hash, &ops, first, window)
                    } else {
                        if !reference_hash.x_last.is_empty() {
                            reference_hash.x_lv2_reset();
                        }
                        reference::<true>(&mut reference_hash, &ops, first, window)
                    };
                    if stage == 3 {
                        // a chunk folded by the reference (the partition tables do not see it)
                        let _ = reference::<false>(&mut hash, &ops, first, window);
                    } else {
                        if !hash.x_last.is_empty() {
                            hash.x_lv2_reset();
                        }
                        let stage_parts = if stage == 6 { parts + 1 } else { parts };
                        let cuts = if n == 0 {
                            vec![0]
                        } else {
                            random_cuts(&mut rng, n)
                        };
                        let ahead = 1 + rng.below(2 * threads as u64) as usize;
                        most_ahead = most_ahead.max(ahead);
                        // LOSAT_X_SEEDPIECE in the stages of the lists and home folds: the
                        // strips' lists, or pieces of 1 seed, a few, or hundreds to thousands
                        let piece_seeds = match (case * 7 + stage) % 4 {
                            0 => 0,
                            1 => 1,
                            2 => 2 + rng.below(40) as usize,
                            _ => 100 + rng.below(4000) as usize,
                        };
                        scratch.x_pieces = Some(piece_seeds);
                        let (got, stats, folders) = partitioned_with(
                            &pool,
                            threads,
                            None,
                            &mut hash,
                            &mut scratch,
                            &ops,
                            first,
                            &cuts,
                            stage_parts,
                            ahead,
                            window,
                            false,
                            owner_waits,
                            // both ways of sharing the folds (LOSAT_X_PAIRPARHOME), the
                            // stages of a case alternating
                            (case + stage) % 2 == 1,
                            // and with or without LOSAT_X_DCSEEDLEAN's lists and untracked
                            // folds: every pairing with the above, for each `levels`, over
                            // two neighbouring cases' stages
                            (case / 2 + stage / 2).is_multiple_of(2),
                        );
                        most_folders = most_folders.max(folders);
                        assert_eq!(got, expected, "{at}: hits");
                        assert_eq!(stats.counts.seeds, n as u64, "{at}");
                        // the strips' buffers stay within the window: the strips claimed
                        // past the slowest partition (the window and the claims that passed
                        // its test at the same time; a strip given back never to the spare
                        // buffers would make one buffer per strip)
                        assert!(
                            scratch.spare.len() <= most_ahead + threads,
                            "{at}: {} strip buffers, window {most_ahead}",
                            scratch.spare.len()
                        );
                        // the pieces: the pool's arenas at most, and every seed a piece of
                        // its own where a piece holds one seed
                        let used = scratch.arenas.iter().filter(|a| a.allocated()).count();
                        assert!(
                            used <= piece::PIECES_PER_THREAD * threads + piece::RESERVED,
                            "{at}: {used} arenas"
                        );
                        // a piece takes the slot of one every partition has folded before a
                        // slot no piece has used, and a strip's next piece in the same turn
                        // only a slot given back: one thread folds a turn's pieces before it
                        // makes the next turn's, so it holds at most two
                        if threads == 1 {
                            assert!(used <= 2, "{at}: {used} arenas on one thread");
                        }
                        let pieced = (case + stage) % 2 == 1
                            && (case / 2 + stage / 2).is_multiple_of(2)
                            && piece_seeds > 0;
                        if pieced && piece_seeds == 1 {
                            assert!(
                                stats.units >= n as u64 * stage_parts as u64,
                                "{at}: {} piece units for {n} seeds",
                                stats.units
                            );
                        }
                    }
                    assert_eq!(
                        hash_state(&hash),
                        hash_state(&reference_hash),
                        "{at}: table"
                    );
                    if stage == 0 {
                        assert!(
                            hash.occupancy > 256,
                            "{at}: the table did not grow ({})",
                            hash.occupancy
                        );
                    }
                    first += n;
                }
                if owner_waits {
                    waited_cases += 1;
                    if most_folders >= 2 {
                        shared_cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(waited_cases, 28);
    // the helpers folded too, in a stage of some case where the owner waited
    assert!(
        shared_cases >= 1,
        "the owner waited in its folds in {waited_cases} cases, but no stage of them was folded on more than one thread"
    );
}

// A parallel iterator inside every scan and fold call of every participant (the pattern
// that hung LOSAT_X_BNSTREAM, #218): a thread waiting inside it (the participant, or a thread
// that stole a piece of the call) may steal a helper job of the stage, which must not wait
// for the stage while the stage waits for the call below it (the helper returns at once on a
// participant's thread, and leaves a stage that has stalled for `HELPER_STALL`). The stages
// end and give the reference's hits and table.
#[test]
fn a_nested_parallel_call_in_every_participant_does_not_hang() {
    let mut rng = Lcg(0x4E57_ED00);
    for threads in 2..=8usize {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut reference_hash = DiagHashTable::new(16);
        let mut hash = copy_hash(&reference_hash);
        let mut scratch = PartScratch::default();
        let mut s = 50i32;
        let mut first = 0usize;
        for stage in 0..4 {
            let n = 2000 + rng.below(2000) as usize;
            let ops = random_ops(&mut rng, n, &mut s, 400);
            let expected = reference::<true>(&mut reference_hash, &ops, first, 40);
            let cuts = random_cuts(&mut rng, n);
            // stages 2 and 3 in pieces (LOSAT_X_SEEDPIECE) of a few seeds and of hundreds
            scratch.x_pieces = match stage {
                2 => Some(1 + rng.below(16) as usize),
                3 => Some(100 + rng.below(400) as usize),
                _ => None,
            };
            let (got, _, _) = partitioned(
                &pool,
                &mut hash,
                &mut scratch,
                &ops,
                first,
                &cuts,
                threads,
                2 * threads,
                40,
                true,
                false,
                // both ways of sharing the folds (LOSAT_X_PAIRPARHOME): stage 1 without
                stage != 1,
            );
            assert_eq!(got, expected, "threads {threads} stage {stage}");
            assert_eq!(
                hash_state(&hash),
                hash_state(&reference_hash),
                "threads {threads} stage {stage}"
            );
            first += n;
        }
    }
}

/// (`a_panicking_participant_ends_the_stage*`) `random_cuts` with strips 0 and 2 of `parts`
/// seeds or more, whose first `parts` ops are moved to diagonals of partitions 0, 1, ...,
/// `parts - 1` (`bucket % parts`; the next such diagonal from the op's own): the first unit
/// of strips 0 and 2 (the strip, or its first piece with pieces of `parts` seeds or more)
/// holds a seed of every partition, so that no partition passes it without a fold call.
fn cuts_with_every_partition(rng: &mut Lcg, ops: &mut [Op], parts: usize) -> Vec<usize> {
    let mut cuts = random_cuts(rng, ops.len());
    for i in [0, 2] {
        cuts[i + 1] = cuts[i + 1].max(cuts[i] + parts);
        for k in i + 2..cuts.len() {
            cuts[k] = cuts[k].max(cuts[k - 1]);
        }
        for (p, op) in ops[cuts[i]..cuts[i] + parts].iter_mut().enumerate() {
            while bucket_of_diag(op.diag) % parts != p {
                op.diag += 1;
            }
        }
    }
    cuts
}

/// (`a_panicking_participant_ends_the_stage*`) Makes a case's panicking fold happen on the
/// participant the case names, whatever the OS runs when. (A 1 ms sleep in every fold made it
/// only likely: on a loaded machine a helper that came while the owner slept in a fold of the
/// helper's home partition, which the owner folded as an orphan, `Run::Orphan`, found nothing
/// to do in the closed window, left the stalled stage after `HELPER_STALL` for good, and the
/// owner folded the whole stage: "helper fold, home true: no panic".) The stage's data is
/// `cuts_with_every_partition`'s, with a partition per thread.
/// - "owner fold": a helper's fold of strip 2 or later waits until the owner has entered one
///   (which panics). Until then no partition passes strip 2's first unit: a helper that folds
///   it waits there holding that partition, and the folds before it run freely. Every
///   partition reaches that unit (the owner alone can fold every unit before it), which is
///   then the slowest one, which the window and the pool of pieces always let the stage scan;
///   the helpers are fewer than the partitions, so the owner, which never leaves the stage,
///   folds a partition there that no helper holds: its home, else by `Run::Slowest` (without
///   `home`, any partition).
/// - "helper fold": the owner's folds wait until a helper has entered one (which panics), and
///   a helper's scans other than of strip 0's first unit wait until the owner has entered a
///   fold (or a helper has). The owner's first fold is in strip 0's first unit, where every
///   other partition then waits with a ready unit, so a helper in the stage folds one (its
///   home, else the slowest partition by `Run::Slowest`; without `home`, any). Before that
///   fold a helper finds nothing to do (and may leave a stalled stage) only while strip 0's
///   first unit is being scanned: by a helper, which then folds its home there, or by the
///   owner, and then the next strip (or piece) went to a helper, whose scan waits. Either way
///   a helper stays in the stage, and a helper that starts later finds the folds.
///
/// Each wait gives up 60 s after the case started (`wait_until`), so that a broken guarantee
/// fails the test (no panic) instead of hanging it.
struct Forced {
    place: &'static str,
    owner: Option<usize>,
    from: std::time::Instant,
    owner_folded: AtomicBool,
    helper_folded: AtomicBool,
}

impl Forced {
    /// In the owner's `pool.install`.
    fn new(place: &'static str) -> Self {
        Forced {
            place,
            owner: rayon::current_thread_index(),
            from: std::time::Instant::now(),
            owner_folded: AtomicBool::new(false),
            helper_folded: AtomicBool::new(false),
        }
    }

    /// Before a scan; `first`: of strip 0's first unit.
    fn scan(&self, first: bool) {
        if self.place == "helper fold" && !first && rayon::current_thread_index() != self.owner {
            wait_until(self.from, || {
                self.owner_folded.load(Ordering::Acquire)
                    || self.helper_folded.load(Ordering::Acquire)
            });
        }
    }

    /// Before a fold of strip `strip`: panics or waits as the place says.
    fn fold(&self, strip: u64) {
        let on_owner = rayon::current_thread_index() == self.owner;
        match self.place {
            "owner fold" if strip >= 2 => {
                if on_owner {
                    self.owner_folded.store(true, Ordering::Release);
                    panic!("a test fold panics on the owner");
                }
                wait_until(self.from, || self.owner_folded.load(Ordering::Acquire));
            }
            "helper fold" => {
                if !on_owner {
                    self.helper_folded.store(true, Ordering::Release);
                    panic!("a test fold panics on a helper");
                }
                self.owner_folded.store(true, Ordering::Release);
                wait_until(self.from, || self.helper_folded.load(Ordering::Acquire));
            }
            _ => {}
        }
    }
}

// A participant that panics (in a fold on the owner, in a fold on a helper, in a scan) makes
// the stage panic instead of hang and leaves the table untouched, and the next stage on the
// same scratch gives the reference's result. Every case panics whatever the OS runs when
// (`Forced`).
#[test]
fn a_panicking_participant_ends_the_stage() {
    let mut rng = Lcg(0x9A41_C000);
    for threads in [2usize, 4, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for (case, (place, home)) in ["owner fold", "helper fold", "scan"]
            .into_iter()
            .flat_map(|place| [(place, false), (place, true)])
            .enumerate()
        {
            let mut reference_hash = DiagHashTable::new(16);
            let mut hash = copy_hash(&reference_hash);
            let mut scratch = PartScratch::default();
            let mut s = 50i32;
            let n = 3000;
            let mut ops = random_ops(&mut rng, n, &mut s, 300);
            let cuts = cuts_with_every_partition(&mut rng, &mut ops, threads);
            let expected = reference::<true>(&mut reference_hash, &ops, 0, 40);
            let strips = cuts.len() - 1;
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                pool.install(|| {
                    let forced = Forced::new(place);
                    let scan = |i: usize, out: &mut Vec<Seed>| {
                        if place == "scan" && i == strips / 2 {
                            panic!("a test scan panics");
                        }
                        forced.scan(i == 0);
                        for k in cuts[i]..cuts[i + 1] {
                            out.push(seed_of(k, &ops[k]));
                        }
                    };
                    let fold = |st: &mut PartState, seeds: &[Seed], _pos: &[u32]| {
                        forced.fold(st.out.strip);
                        st.counts.seeds += seeds.len() as u64;
                    };
                    let mut hits = Vec::new();
                    seed_stage_in(
                        home,
                        strips,
                        threads,
                        2,
                        &mut hash,
                        &mut hits,
                        FoldPart {
                            scratch: &mut scratch,
                            parts: threads,
                            latest: None,
                        },
                        &scan,
                        None::<&NoLists>,
                        &fold,
                        false,
                    )
                    .counts
                    .seeds
                })
            }));
            assert!(
                result.is_err(),
                "threads {threads} case {case} ({place}, home {home}): no panic"
            );
            // the table is the one before the stage: a merge never ran
            let fresh = DiagHashTable::new(16);
            assert_eq!(hash_state(&hash), hash_state(&fresh), "{place}");
            let (got, _, _) = partitioned(
                &pool,
                &mut hash,
                &mut scratch,
                &ops,
                0,
                &cuts,
                threads,
                2,
                40,
                false,
                false,
                !home,
            );
            assert_eq!(got, expected, "threads {threads} {place} home {home}");
            assert_eq!(
                hash_state(&hash),
                hash_state(&reference_hash),
                "threads {threads} {place}"
            );
        }
    }
}

// EXPERIMENT (LOSAT_X_SEEDPIECE): a participant that panics in the stage in pieces (in a
// piece's scan, in a fold on the owner, in a fold on a helper) makes the stage panic instead
// of hang and leaves the table untouched, and the next stage on the same scratch gives the
// reference's result. Every case panics whatever the OS runs when (`Forced`; pieces of 8
// seeds, so that a piece can hold a seed of each of 8 partitions).
#[test]
fn a_panicking_participant_ends_the_stage_in_pieces() {
    // (the seeds of a piece)
    const PIECE: usize = 8;
    let mut rng = Lcg(0x9A41_C0DE);
    for threads in [2usize, 4, 8] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for place in ["scan", "owner fold", "helper fold"] {
            let mut reference_hash = DiagHashTable::new(16);
            let mut hash = copy_hash(&reference_hash);
            let mut scratch = PartScratch::default();
            let mut s = 50i32;
            let n = 3000;
            let mut ops = random_ops(&mut rng, n, &mut s, 300);
            let cuts = cuts_with_every_partition(&mut rng, &mut ops, threads);
            let expected = reference::<true>(&mut reference_hash, &ops, 0, 40);
            let strips = cuts.len() - 1;
            // ("scan": the first strip from the middle on with a second piece)
            let panics = (strips / 2..strips)
                .find(|&i| cuts[i + 1] - cuts[i] > PIECE)
                .unwrap();
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                pool.install(|| {
                    let forced = Forced::new(place);
                    let scan = |i: usize, from: Option<Resume>, out: &mut PieceSink<'_>| {
                        if place == "scan" && i == panics && from.is_some() {
                            panic!("a test piece scan panics");
                        }
                        forced.scan(i == 0 && from.is_none());
                        let start = cuts[i] + from.map_or(0, |r| r.s);
                        for k in start..cuts[i + 1] {
                            if out.room() == 0 {
                                return Some(Resume::at(k - cuts[i]));
                            }
                            out.push(seed_of(k, &ops[k]));
                        }
                        None
                    };
                    let fold = |st: &mut PartState, seeds: &[Seed], _pos: &[u32]| {
                        forced.fold(st.out.strip);
                        st.counts.seeds += seeds.len() as u64;
                    };
                    let mut hits = Vec::new();
                    rebuild(&mut scratch, &hash, threads, true, None);
                    let occ0 = hash.occupancy;
                    piece::seed_stage_pieces(
                        strips,
                        threads,
                        2,
                        &mut hash,
                        &mut hits,
                        &mut scratch,
                        threads,
                        occ0,
                        PIECE,
                        &scan,
                        &fold,
                        false,
                    )
                    .counts
                    .seeds
                })
            }));
            assert!(result.is_err(), "threads {threads} ({place}): no panic");
            let fresh = DiagHashTable::new(16);
            assert_eq!(hash_state(&hash), hash_state(&fresh), "{place}");
            // (the arenas of the failed stage are lost with it; the next stage makes its own)
            scratch.x_pieces = Some(29);
            let (got, _, _) = partitioned(
                &pool,
                &mut hash,
                &mut scratch,
                &ops,
                0,
                &cuts,
                threads,
                2,
                40,
                false,
                false,
                true,
            );
            assert_eq!(got, expected, "threads {threads} {place}");
            assert_eq!(
                hash_state(&hash),
                hash_state(&reference_hash),
                "threads {threads} {place}"
            );
        }
    }
}

// LOSAT_X_PAIRPARHOME: a home participant whose scans are slow (1 ms each on one pool
// thread; its partition lags while it scans, and it does not hold its partition's lock
// then) has its partition folded by the others (`Run::Slowest`, and at the tail
// `Run::Tail`), and a stage on more participants than the pool has threads (the helper jobs
// start only after the owner finished the stage alone) has its partitions folded without a
// participant at home (`Run::Orphan`). Both give the reference's hits and table. (Who folds
// what depends on the timing: the slow cases get up to five stages to fold a strip away
// from home.)
#[test]
fn home_partitions_are_folded_away_from_home() {
    let mut rng = Lcg(0x4053_A11D);
    for (pool_threads, threads, slow) in [(4usize, 4usize, true), (8, 8, true), (1, 4, false)] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(pool_threads)
            .build()
            .unwrap();
        let owner = pool.install(rayon::current_thread_index);
        let slow = slow.then(|| Slow::Scans((owner.unwrap_or(0) + 1) % pool_threads));
        let mut reference_hash = DiagHashTable::new(16);
        let mut hash = copy_hash(&reference_hash);
        let mut scratch = PartScratch::default();
        let mut s = 50i32;
        let mut first = 0usize;
        let mut away = false;
        for stage in 0..if slow.is_some() { 5 } else { 2 } {
            let at = format!("pool {pool_threads} threads {threads} slow {slow:?} stage {stage}");
            let n = 3000;
            let ops = random_ops(&mut rng, n, &mut s, 400);
            let expected = reference::<true>(&mut reference_hash, &ops, first, 40);
            let cuts: Vec<usize> = (0..=60).map(|k| k * n / 60).collect();
            let (got, stats, _) = partitioned_on(
                &pool,
                threads,
                slow,
                &mut hash,
                &mut scratch,
                &ops,
                first,
                &cuts,
                threads,
                2,
                40,
                false,
                false,
                true,
            );
            assert_eq!(got, expected, "{at}: hits");
            assert_eq!(
                hash_state(&hash),
                hash_state(&reference_hash),
                "{at}: table"
            );
            if pool_threads == 1 {
                // the owner alone, every partition but its own without a participant at home
                assert_eq!(stats.away, 60 * (threads as u64 - 1), "{at}");
            }
            first += n;
            away |= stats.away > 0;
            if away && slow.is_some() {
                break;
            }
        }
        assert!(
            away,
            "pool {pool_threads} threads {threads} slow {slow:?}: no strip folded away from home"
        );
    }
}

// A stage that stalls (the owner holds its first strip, `Slow::OwnerFirstScan`, and no
// partition can fold past that strip) makes the idle helpers leave once it has made no
// progress for `HELPER_STALL` (`helpers_left`); the owner finishes the stage alone, with
// LOSAT_X_PAIRPARHOME folding the helpers' home partitions away from home (they have no
// participant at home any more). Both ways give the reference's hits and table. The owner
// holds the strip from the time a helper is in the stage until every helper has left, and
// the helpers' scans wait for the owner's first one (`hold_until_helpers_leave`), so one
// stage of each way shows it on any machine, whatever its CPU count and load: this test
// needs no stage in which the OS happens to run a helper while the owner sleeps, nor one in
// which the owner takes a strip before the helpers have done the whole stage.
// (`a_stall_is_seen_on_the_first_turn_after_helper_stall` checks when a helper looks.)
#[test]
fn a_stalled_stage_lets_the_helpers_leave() {
    let mut rng = Lcg(0x57A1_1ED0);
    for threads in [2usize, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        // (and LOSAT_X_SEEDPIECE's stage, in pieces of 37 seeds: the owner holds its first
        // piece, strip 0 half-scanned)
        for (home, pieces) in [(false, None), (true, None), (true, Some(37))] {
            let at = format!("threads {threads} home {home} pieces {pieces:?}");
            let mut reference_hash = DiagHashTable::new(16);
            let mut hash = copy_hash(&reference_hash);
            let mut scratch = PartScratch {
                x_pieces: pieces,
                ..PartScratch::default()
            };
            let mut s = 50i32;
            let n = 2000;
            let ops = random_ops(&mut rng, n, &mut s, 300);
            let expected = reference::<true>(&mut reference_hash, &ops, 0, 40);
            let cuts: Vec<usize> = (0..=40).map(|k| k * n / 40).collect();
            let started = std::time::Instant::now();
            let (got, stats, _) = partitioned_on(
                &pool,
                threads,
                Some(Slow::OwnerFirstScan),
                &mut hash,
                &mut scratch,
                &ops,
                0,
                &cuts,
                threads,
                2 * threads,
                40,
                false,
                false,
                home,
            );
            assert_eq!(got, expected, "{at}: hits");
            assert_eq!(
                hash_state(&hash),
                hash_state(&reference_hash),
                "{at}: table"
            );
            assert!(
                stats.helpers_left > 0,
                "{at}: no helper left the stalled stage in {:?}",
                started.elapsed()
            );
            if home {
                assert!(stats.away > 0, "{at}: no strip folded away from home");
            }
        }
    }
}

// A helper looks at the stage on every turn without work after the spinning ones (`Stall`),
// so it leaves on the first look after the stage has made no progress for `HELPER_STALL`,
// however few turns that took (one turn that yields can take a scheduler quantum on a loaded
// machine); the spinning turns do not read the progress, and a change of the progress or work
// done starts the watch again.
#[test]
fn a_stall_is_seen_on_the_first_turn_after_helper_stall() {
    let mut stall = Stall::default();
    for idle in 1..SPIN_TURNS {
        assert!(!stall.leaves(idle, || unreachable!("a spinning turn reads the progress")));
    }
    assert!(!stall.leaves(SPIN_TURNS, || 7));
    std::thread::sleep(HELPER_STALL);
    // the stage moved: the watch starts again
    assert!(!stall.leaves(SPIN_TURNS + 1, || 8));
    assert!(!stall.leaves(SPIN_TURNS + 2, || 8));
    std::thread::sleep(HELPER_STALL);
    assert!(stall.leaves(SPIN_TURNS + 3, || 8));
    // work done: the watch starts again
    stall.reset();
    assert!(!stall.leaves(SPIN_TURNS, || 8));
    std::thread::sleep(HELPER_STALL);
    assert!(stall.leaves(SPIN_TURNS + 1, || 8));
}

// Every stage starts from a copy of the stage's table: a change of the table between two
// partitioned stages on one scratch that keeps the occupancy (here the level of every cell
// raised, which no fold does) is seen by the next stage, which gives the reference's result
// on the changed table and not the one on the table before the change.
#[test]
fn a_table_changed_between_stages_is_seen() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let mut rng = Lcg(0xBAD0_7AB1);
    let mut reference_hash = DiagHashTable::new(16);
    let mut hash = copy_hash(&reference_hash);
    let mut scratch = PartScratch::default();
    let mut s = 50i32;
    let ops = random_ops(&mut rng, 2000, &mut s, 300);
    let expected = reference::<true>(&mut reference_hash, &ops, 0, 40);
    let cuts = random_cuts(&mut rng, ops.len());
    let (got, _, _) = partitioned(
        &pool,
        &mut hash,
        &mut scratch,
        &ops,
        0,
        &cuts,
        2,
        2,
        40,
        false,
        false,
        true,
    );
    assert_eq!(got, expected, "first stage: hits");
    assert_eq!(
        hash_state(&hash),
        hash_state(&reference_hash),
        "first stage: table"
    );
    // the same diagonals again, further on
    let more: Vec<Op> = ops
        .iter()
        .map(|op| Op {
            s_off_pos: op.s_off_pos + 1000,
            s_end_pos: op.s_end_pos + 1000,
            ungapped_end: op.ungapped_end + 1000,
            ..*op
        })
        .collect();
    let unchanged = reference::<true>(&mut copy_hash(&reference_hash), &more, 2000, 40);
    let occupancy = hash.occupancy;
    for h in [&mut hash, &mut reference_hash] {
        for cell in &mut h.chain[1..occupancy as usize] {
            cell.level = i32::MAX / 2;
        }
        // the bounds (LOSAT_X_DIAGLV) and the memo (LOSAT_X_DIAGLV2) describe the cells the
        // inserts wrote; made again from the changed cells
        if h.x_lv_on() {
            let memo = !h.x_last.is_empty();
            h.x_lv_start();
            if memo {
                h.x_lv2_start();
            }
        }
    }
    let expected = reference::<true>(&mut reference_hash, &more, 2000, 40);
    assert_ne!(expected, unchanged, "the change does not matter");
    let cuts = random_cuts(&mut rng, more.len());
    let (got, _, _) = partitioned(
        &pool,
        &mut hash,
        &mut scratch,
        &more,
        2000,
        &cuts,
        2,
        2,
        40,
        false,
        false,
        false,
    );
    assert_eq!(got, expected, "second stage: hits");
    assert_eq!(
        hash_state(&hash),
        hash_state(&reference_hash),
        "second stage: table"
    );
}
