use super::super::{build_disc_mb_lookup, compute_mb_pv_params};
use super::*;
use crate::algorithm::blastn::disc_lookup::{get_disc_template_type, DiscWordType};
use crate::utils::dust::MaskedInterval;
use std::collections::BTreeSet;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:785-786,803-804
// Random queries (with runs of ambiguity codes, masks, queries shorter than the template and
// queries without a mask list) through the original loop and through the partitioned fill:
// every element of every table must be the same.
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

struct Case {
    queries: Vec<Vec<u8>>,
    offsets: Vec<i32>,
    masks: Vec<Vec<MaskedInterval>>,
}

/// One to four queries in BLASTNA: lengths 1..3,000 (one in four shorter than any template),
/// half of the letters from a short repeated motif (long chains), 2 % ambiguity codes (some in
/// runs), random masks or none; the last query may have no mask list at all.
fn random_case(rng: &mut Lcg) -> Case {
    let n_queries = 1 + rng.below(4) as usize;
    let mut queries = Vec::new();
    let mut offsets = Vec::new();
    let mut masks = Vec::new();
    let mut offset = 0i32;
    for _ in 0..n_queries {
        let len = if rng.below(4) == 0 {
            1 + rng.below(20) as usize
        } else {
            1 + rng.below(3000) as usize
        };
        let motif: Vec<u8> = (0..1 + rng.below(40)).map(|_| rng.below(4) as u8).collect();
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            let r = rng.below(100);
            if r < 1 {
                let run = 1 + rng.below(30) as usize;
                seq.extend(std::iter::repeat_n(14u8, run)); // N
            } else if r < 2 {
                seq.push(4 + rng.below(12) as u8); // another ambiguity code
            } else if r < 50 {
                seq.push(motif[seq.len() % motif.len()]);
            } else {
                seq.push(rng.below(4) as u8);
            }
        }
        seq.truncate(len);
        let mut m = Vec::new();
        let mut cursor = 0usize;
        while cursor < len && rng.below(3) != 0 {
            let start = cursor + rng.below(400) as usize;
            let end = start + 1 + rng.below(200) as usize;
            if start >= len {
                break;
            }
            m.push(MaskedInterval {
                start,
                end: end.min(len),
            });
            cursor = end;
        }
        offsets.push(offset);
        offset += len as i32 + 1;
        queries.push(seq);
        masks.push(m);
    }
    if rng.below(4) == 0 {
        masks.pop();
    }
    Case {
        queries,
        offsets,
        masks,
    }
}

/// Every (word size, template type, template length) that a discontiguous table has.
fn combos() -> Vec<(usize, DiscWordType, usize)> {
    let mut v = Vec::new();
    for word_size in [11usize, 12] {
        for word_type in [
            DiscWordType::Coding,
            DiscWordType::Optimal,
            DiscWordType::TwoTemplates,
        ] {
            for template_length in [16usize, 18, 21] {
                v.push((word_size, word_type, template_length));
            }
        }
    }
    v
}

/// The prefault rule of the tests' partitioned builds: a round whatever the first chunk's size,
/// and regions of one page (every marked page is touched).
const EVERY_PAGE: PrefaultRule = PrefaultRule {
    min_words_per_page: 0,
    region_shift: 0,
};

#[derive(Clone, Debug, PartialEq, Eq)]
struct Built {
    hashtable: Vec<u32>,
    next_pos: Vec<u32>,
    helper_array: Vec<u32>,
    hashtable2: Vec<u32>,
    next_pos2: Vec<u32>,
    helper_array2: Vec<u32>,
    pv_array: Vec<PvArrayType>,
    longest_chain: usize,
}

enum How {
    Serial,
    /// `fill_with` off the pool with that many parts, positions per chunk and prefault rule
    /// (the other partitioned variants use `EVERY_PAGE`).
    Rule(usize, usize, PrefaultRule),
    /// `fill` (the default chunk size), on the pool when `true`.
    Parts(usize, bool),
    /// `fill_with` with that many parts and positions per chunk, on the pool when `true`.
    Chunked(usize, usize, bool),
}

/// The tables of `case` allocated as `build_disc_mb_lookup` allocates them and filled `how`.
fn build(case: &Case, combo: (usize, DiscWordType, usize), how: How) -> Built {
    build_stats(case, combo, how).0
}

/// `build`, and what the fill did (`How::Chunked` only; else the default).
fn build_stats(case: &Case, combo: (usize, DiscWordType, usize), how: How) -> (Built, FillStats) {
    build_from(case, combo, how, None)
}

/// `build_stats`, with the counters set to `start(h)` (counter `h` of `helper_array`) and
/// `start(len + h)` (of `helper_array2`) before the fill when `start` is given; else zero, as
/// `build_disc_mb_lookup` allocates them.
fn build_from(
    case: &Case,
    combo: (usize, DiscWordType, usize),
    how: How,
    start: Option<&dyn Fn(usize) -> u32>,
) -> (Built, FillStats) {
    let (word_size, word_type, template_length) = combo;
    let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
    let hashsize = 1usize << (2 * word_size);
    let approx: usize = case.queries.iter().map(Vec::len).sum();
    let (pv_size, pv_array_bts) = compute_mb_pv_params(hashsize, approx, false, word_size);
    let template_type = get_disc_template_type(word_size as i32, template_length as u8, word_type);
    let two_templates = word_type == DiscWordType::TwoTemplates;
    let second_template_type = if two_templates {
        template_type.next()
    } else {
        crate::algorithm::blastn::disc_lookup::DiscTemplateType::Contiguous
    };
    let query_length = case
        .queries
        .iter()
        .zip(&case.offsets)
        .map(|(q, &o)| o as usize + q.len())
        .max()
        .unwrap_or(0);
    let second = |n: usize| {
        if two_templates {
            vec![0u32; n]
        } else {
            Vec::new()
        }
    };
    let mut b = Built {
        hashtable: vec![0; hashsize],
        next_pos: vec![0; query_length + 1],
        helper_array: vec![0; hashsize / 2048],
        hashtable2: second(hashsize),
        next_pos2: second(query_length + 1),
        helper_array2: second(hashsize / 2048),
        pv_array: vec![0; pv_size],
        longest_chain: 0,
    };
    if let Some(start) = start {
        let len = b.helper_array.len();
        for (h, c) in b.helper_array.iter_mut().enumerate() {
            *c = start(h);
        }
        for (h, c) in b.helper_array2.iter_mut().enumerate() {
            *c = start(len + h);
        }
    }
    let words = DiscWords {
        queries_blastna: &query_refs,
        query_offsets: &case.offsets,
        query_masks: &case.masks,
        template_length,
        template_type,
        second_template_type,
        two_templates,
    };
    let tables = DiscTables {
        hashtable: &mut b.hashtable,
        next_pos: &mut b.next_pos,
        helper_array: &mut b.helper_array,
        hashtable2: &mut b.hashtable2,
        next_pos2: &mut b.next_pos2,
        helper_array2: &mut b.helper_array2,
        pv_array: &mut b.pv_array,
        pv_array_bts,
    };
    let stats = match how {
        How::Serial => {
            fill_disc_mb_tables(&words, tables);
            FillStats::default()
        }
        // `fill` is `fill_with` in the default chunks without a probe; here the prefault round
        // runs whatever the first chunk's size and touches every marked page (`EVERY_PAGE`),
        // so that every partitioned build checks it (`fill`'s rule has tests of its own).
        How::Parts(parts, in_pool) => fill_with(
            parts,
            &words,
            tables,
            in_pool,
            CHUNK_POSITIONS,
            EVERY_PAGE,
            &|_| {},
        ),
        How::Chunked(parts, chunk, in_pool) => {
            fill_with(parts, &words, tables, in_pool, chunk, EVERY_PAGE, &|_| {})
        }
        How::Rule(parts, chunk, rule) => {
            fill_with(parts, &words, tables, false, chunk, rule, &|_| {})
        }
    };
    if !matches!(how, How::Serial) {
        check_prefault(&b, &stats);
    }
    b.longest_chain = disc_longest_chain(&b.helper_array, &b.helper_array2, two_templates);
    (b, stats)
}

/// The prefault round of a partitioned fill of `b` (the arrays it filled in place):
/// `check_prefault_arrays` on `hashtable`, `hashtable2` and `pv_array`.
fn check_prefault(b: &Built, stats: &FillStats) -> [BTreeSet<usize>; 3] {
    check_prefault_arrays([&b.hashtable, &b.hashtable2, &b.pv_array], stats)
}

/// The prefault round of a fill of `arrays` (`hashtable`, `hashtable2`, `pv_array`, as filled
/// in place), from its log of writes (`FillStats::touch_log`): each wrote 0 over 0; there was
/// at most one round; and the pages it wrote to are, array by array, exactly the pages the
/// first chunk wrote (`FillStats::first_pages`) in the regions of `2^region_shift` pages that lie
/// within the array's pages and that the first chunk wrote in full, when the round ran, none
/// otherwise. So it maps no page (and no region) the fill does not write, and misses none it
/// should touch; and those pages hold a value after the fill. Returns the pages touched per
/// array.
fn check_prefault_arrays(arrays: [&[u32]; 3], stats: &FillStats) -> [BTreeSet<usize>; 3] {
    let mut touched: [BTreeSet<usize>; 3] = Default::default();
    for t in &stats.touch_log {
        assert!(
            t.before == 0 && t.after == 0,
            "the prefault wrote {} over {} at element {} of array {}",
            t.after,
            t.before,
            t.index,
            t.array
        );
        let v = arrays[t.array];
        assert!(t.index < v.len(), "{t:?}");
        touched[t.array].insert((v.as_ptr() as usize + 4 * t.index) / PAGE_BYTES);
    }
    assert_eq!(stats.touches, stats.touch_log.len());
    assert!(
        stats.prefault_rounds <= 1,
        "{} rounds",
        stats.prefault_rounds
    );
    for k in 0..3 {
        // The first chunk's pages in the regions (`2^region_shift` pages by address) that lie
        // within the array's pages and that the first chunk wrote in full.
        let first = if stats.prefault_rounds == 1 && !arrays[k].is_empty() {
            let (a0, a1) = (
                arrays[k].as_ptr() as usize / PAGE_BYTES,
                (arrays[k].as_ptr() as usize + 4 * arrays[k].len() - 1) / PAGE_BYTES,
            );
            let shift = stats.region_shift;
            stats.first_pages[k]
                .iter()
                .copied()
                .filter(|&p| {
                    let region = p >> shift;
                    let (lo, hi) = (region << shift, ((region + 1) << shift) - 1);
                    lo >= a0 && hi <= a1 && (lo..=hi).all(|q| stats.first_pages[k].contains(&q))
                })
                .collect()
        } else {
            BTreeSet::new()
        };
        assert!(
            touched[k] == first,
            "array {k}: touched {} pages, the first chunk wrote {}; touched but not written {:?}, written but not touched {:?}",
            touched[k].len(),
            first.len(),
            touched[k].difference(&first).take(5).collect::<Vec<_>>(),
            first.difference(&touched[k]).take(5).collect::<Vec<_>>()
        );
        assert!(touched[k].is_subset(&nonzero_pages(arrays[k])), "array {k}");
    }
    touched
}

/// The serial tables have work in them: chains (counters above zero), and for two templates
/// both tables.
fn assert_not_trivial(b: &Built, two_templates: bool) {
    assert!(b.hashtable.iter().any(|&c| c != 0));
    assert!(b.pv_array.iter().any(|&w| w != 0));
    if two_templates {
        assert!(b.hashtable2.iter().any(|&c| c != 0));
    }
}

// The helper above allocates as `build_disc_mb_lookup` does: its serial tables are the
// function's (which runs the original loop on this test thread, which is no pool worker).
#[test]
fn helper_builds_what_build_disc_mb_lookup_builds() {
    let mut rng = Lcg(0xD15C_0001);
    for combo in combos() {
        let case = random_case(&mut rng);
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let approx: usize = case.queries.iter().map(Vec::len).sum();
        let lookup = build_disc_mb_lookup(
            &query_refs,
            &case.offsets,
            combo.0,
            combo.2,
            combo.1,
            &case.masks,
            approx,
        );
        let b = build(&case, combo, How::Serial);
        let disc = lookup.disc().unwrap();
        assert_eq!(lookup.mb_lookup.hashtable, b.hashtable, "{combo:?}");
        assert_eq!(lookup.mb_lookup.next_pos, b.next_pos, "{combo:?}");
        assert_eq!(lookup.mb_lookup.pv_array, b.pv_array, "{combo:?}");
        assert_eq!(disc.hashtable2, b.hashtable2, "{combo:?}");
        assert_eq!(disc.next_pos2, b.next_pos2, "{combo:?}");
        assert_eq!(lookup.longest_chain(), b.longest_chain, "{combo:?}");
    }
}

// Partitioned in P = 1, 2, 3, 7 and 16 parts one after another (every template type and
// length, one and two templates), with the default chunk and with chunks of 1 to 4,096
// positions (one per part count, rotating by round): the same tables as the original loop.
#[test]
fn partitioned_fill_builds_the_same_tables() {
    let mut rng = Lcg(0xD15C_5EED);
    let mut chains = 0;
    for round in 0..3 {
        for combo in combos() {
            let case = random_case(&mut rng);
            let reference = build(&case, combo, How::Serial);
            if round == 0 && case.queries.iter().map(Vec::len).sum::<usize>() > 200 {
                assert_not_trivial(&reference, combo.1 == DiscWordType::TwoTemplates);
            }
            chains += usize::from(reference.longest_chain > 3);
            for (k, parts) in [1usize, 2, 3, 7, 16].into_iter().enumerate() {
                let got = build(&case, combo, How::Parts(parts, false));
                assert!(
                    got == reference,
                    "round {round} {combo:?} parts {parts}: tables differ"
                );
                let chunk = [1usize, 17, 21, 333, 4096][(k + round) % 5];
                let got = build(&case, combo, How::Chunked(parts, chunk, false));
                assert!(
                    got == reference,
                    "round {round} {combo:?} parts {parts} chunk {chunk}: tables differ"
                );
            }
        }
    }
    assert!(chains > 10, "the cases must have long chains ({chains})");
}

// Chunk and span boundaries anywhere relative to the words (chunks of 1 to 1,000 positions, so
// spans start inside words, at ambiguity codes and at range starts, and many spans are empty):
// the same tables as the original loop, for 1, 3 and 16 parts.
#[test]
fn chunk_boundaries_anywhere() {
    let mut rng = Lcg(0xC4A2_B0DE);
    for case_no in 0..6 {
        let case = random_case(&mut rng);
        let combo = [
            (11usize, DiscWordType::Coding, 18usize),
            (11, DiscWordType::TwoTemplates, 21),
            (11, DiscWordType::Optimal, 16),
        ][case_no % 3];
        let reference = build(&case, combo, How::Serial);
        for chunk in [1usize, 2, 5, 17, 18, 19, 20, 21, 22, 64, 1000] {
            for parts in [1usize, 3, 16] {
                let got = build(&case, combo, How::Chunked(parts, chunk, false));
                assert!(
                    got == reference,
                    "case {case_no} {combo:?} chunk {chunk} parts {parts}: tables differ"
                );
            }
        }
    }
}

// Queries shorter than every template and empty mask lists (no mask list at all, and lists
// without intervals), and a query that is all ambiguity codes: no word, and the same (empty)
// tables either way.
#[test]
fn short_and_unmasked_queries() {
    let case = Case {
        queries: vec![vec![0, 1, 2, 3, 0, 1, 2], vec![14; 40], vec![3; 15]],
        offsets: vec![0, 8, 49],
        masks: vec![Vec::new()],
    };
    for combo in combos() {
        let reference = build(&case, combo, How::Serial);
        assert!(reference.hashtable.iter().all(|&c| c == 0), "{combo:?}");
        for parts in [2usize, 3, 7, 16] {
            assert!(build(&case, combo, How::Parts(parts, false)) == reference);
        }
    }
    // One long query without masks (`masks` empty), then the same with an empty list.
    let mut rng = Lcg(0xD15C_0E11);
    let seq: Vec<u8> = (0..5000).map(|_| rng.below(4) as u8).collect();
    for masks in [Vec::new(), vec![Vec::new()]] {
        let case = Case {
            queries: vec![seq.clone()],
            offsets: vec![0],
            masks,
        };
        for combo in combos() {
            let reference = build(&case, combo, How::Serial);
            assert_not_trivial(&reference, combo.1 == DiscWordType::TwoTemplates);
            for parts in [2usize, 3, 7, 16] {
                assert!(build(&case, combo, How::Parts(parts, false)) == reference);
            }
        }
    }
}

// On real pools (`with_search_pool`, the search's pool) of 1, 2, 3 and 8 threads: the
// dispatch gives the pool's size from 2 threads on and nothing on 1, and the fill on the pool
// (also with more parts than threads, and with chunks of 97 positions) builds the same tables
// as the original loop.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn fill_on_a_pool_builds_the_same_tables() {
    let mut rng = Lcg(0xFACE_D15C);
    let cases: Vec<Case> = (0..6).map(|_| random_case(&mut rng)).collect();
    for threads in [1usize, 2, 3, 8] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            let parts = parts_for(true);
            assert_eq!(parts, (threads > 1).then_some(threads), "threads {threads}");
            for (i, combo) in combos().into_iter().enumerate() {
                let case = &cases[i % cases.len()];
                let reference = build(case, combo, How::Serial);
                let got = match parts {
                    Some(parts) => build(case, combo, How::Parts(parts, true)),
                    None => build(case, combo, How::Serial),
                };
                assert!(got == reference, "threads {threads} {combo:?}");
                if threads > 1 {
                    let more = build(case, combo, How::Parts(16, true));
                    assert!(more == reference, "threads {threads} {combo:?} 16 parts");
                    let small = build(case, combo, How::Chunked(threads, 97, true));
                    assert!(
                        small == reference,
                        "threads {threads} {combo:?} chunks of 97"
                    );
                }
            }
            Ok(())
        })
        .unwrap();
    }
}

// The dispatch: nothing with the switch off, nothing off a pool (this test thread is no
// rayon worker), nothing on a search "pool" of one thread, the pool's size from two threads.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn dispatch_needs_a_search_pool_of_two_threads() {
    assert_eq!(rayon::current_thread_index(), None);
    assert_eq!(parts_for(false), None);
    assert_eq!(parts_for(true), None);
    for (threads, want) in [(1usize, None), (2, Some(2)), (5, Some(5))] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            assert_eq!(parts_for(true), want, "threads {threads}");
            assert_eq!(parts_for(false), None, "threads {threads}");
            Ok(())
        })
        .unwrap();
    }
}

/// The fill on a pool of three threads runs its parts on the workers of that pool, and on
/// more than one of them: the probe that each phase-B job (owner) calls first sees a rayon
/// worker of a pool of three threads, and the owners are spread over at least two workers. (Each probe call
/// waits until a second worker has called it (the first call that waits 10 s in vain ends
/// all waiting), so that the spread does not depend on the scheduler; a fill that ran its
/// parts one after another on one thread has no second worker and fails after those 10 s.)
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_fill_runs_on_several_workers_of_the_pool() {
    use std::collections::BTreeSet;
    use std::sync::atomic::{AtomicBool, AtomicUsize};
    use std::sync::Mutex;
    use std::time::{Duration, Instant};

    let mut rng = Lcg(0x0B5E_12FE);
    let case = Case {
        queries: vec![(0..20_000).map(|_| rng.below(4) as u8).collect()],
        offsets: vec![0],
        masks: Vec::new(),
    };
    let combo = (11usize, DiscWordType::TwoTemplates, 18usize);
    let reference = build(&case, combo, How::Serial);
    crate::utils::threading::with_search_pool(3, "blastn", |_pool| {
        let parts = parts_for(true).expect("a pool of three threads");
        assert_eq!(parts, 3);
        let calls = AtomicUsize::new(0);
        let in_pool_of_three = AtomicUsize::new(0);
        let workers = Mutex::new(BTreeSet::new());
        let gave_up = AtomicBool::new(false);
        let probe = |_part: u8| {
            calls.fetch_add(1, Ordering::SeqCst);
            if let Some(index) = rayon::current_thread_index() {
                workers.lock().unwrap().insert(index);
                if rayon::current_num_threads() == 3 {
                    in_pool_of_three.fetch_add(1, Ordering::SeqCst);
                }
            }
            let start = Instant::now();
            while !gave_up.load(Ordering::SeqCst) && workers.lock().unwrap().len() < 2 {
                if start.elapsed() > Duration::from_secs(10) {
                    gave_up.store(true, Ordering::SeqCst);
                }
                std::thread::yield_now();
            }
        };
        let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
        let mut got = build(&case, combo, How::Serial);
        got.hashtable.fill(0);
        got.next_pos.fill(0);
        got.helper_array.fill(0);
        got.hashtable2.fill(0);
        got.next_pos2.fill(0);
        got.helper_array2.fill(0);
        got.pv_array.fill(0);
        let (_, pv_array_bts) = compute_mb_pv_params(1 << 22, 20_000, false, 11);
        let template_type = get_disc_template_type(11, 18, DiscWordType::TwoTemplates);
        fill_with(
            parts,
            &DiscWords {
                queries_blastna: &query_refs,
                query_offsets: &case.offsets,
                query_masks: &case.masks,
                template_length: 18,
                template_type,
                second_template_type: template_type.next(),
                two_templates: true,
            },
            DiscTables {
                hashtable: &mut got.hashtable,
                next_pos: &mut got.next_pos,
                helper_array: &mut got.helper_array,
                hashtable2: &mut got.hashtable2,
                next_pos2: &mut got.next_pos2,
                helper_array2: &mut got.helper_array2,
                pv_array: &mut got.pv_array,
                pv_array_bts,
            },
            true,
            CHUNK_POSITIONS,
            EVERY_PAGE,
            &probe,
        );
        got.longest_chain = disc_longest_chain(&got.helper_array, &got.helper_array2, true);
        assert!(got == reference, "the tables differ");
        assert_eq!(calls.load(Ordering::SeqCst), 3);
        assert_eq!(in_pool_of_three.load(Ordering::SeqCst), 3);
        let workers = workers.lock().unwrap();
        assert!(
            workers.len() >= 2 && workers.iter().all(|&index| index < 3),
            "workers {workers:?}"
        );
        Ok(())
    })
    .unwrap();
}

/// `DiscWords` of `case` with the templates of `combo`.
fn words_of<'a>(
    case: &'a Case,
    query_refs: &'a [&'a [u8]],
    combo: (usize, DiscWordType, usize),
) -> DiscWords<'a> {
    let template_type = get_disc_template_type(combo.0 as i32, combo.2 as u8, combo.1);
    let two_templates = combo.1 == DiscWordType::TwoTemplates;
    DiscWords {
        queries_blastna: query_refs,
        query_offsets: &case.offsets,
        query_masks: &case.masks,
        template_length: combo.2,
        template_type,
        second_template_type: if two_templates {
            template_type.next()
        } else {
            crate::algorithm::blastn::disc_lookup::DiscTemplateType::Contiguous
        },
        two_templates,
    }
}

// The premise of the warm-up of phase A: every index function of the twelve templates (the
// second template of two included) reads only the last `template_length` bases of `accum`
// (its low `2 * template_length` bits), and the index is a cell of the table.
#[test]
fn template_indices_use_only_the_last_template_length_bases() {
    let mut rng = Lcg(0x1234_5678);
    let mut uses = 0;
    for word_size in [11i32, 12] {
        for template_length in [16u8, 18, 21] {
            for word_type in [
                DiscWordType::Coding,
                DiscWordType::Optimal,
                DiscWordType::TwoTemplates,
            ] {
                let first = get_disc_template_type(word_size, template_length, word_type);
                let mut templates = vec![first];
                if word_type == DiscWordType::TwoTemplates {
                    templates.push(first.next());
                }
                let mask = (1u64 << (2 * template_length as u64)) - 1;
                for t in templates {
                    assert!(
                        t != crate::algorithm::blastn::disc_lookup::DiscTemplateType::Contiguous
                    );
                    for _ in 0..20_000 {
                        let accum = (rng.next() << 33) ^ rng.next() ^ (rng.next() << 62);
                        let index = compute_discontiguous_index(accum, t);
                        assert_eq!(
                            index,
                            compute_discontiguous_index(accum & mask, t),
                            "{t:?} template_length {template_length}"
                        );
                        assert!((index as u64) < (1u64 << (2 * word_size as u64)), "{t:?}");
                    }
                    uses += 1;
                }
            }
        }
    }
    assert_eq!(uses, 24);
}

/// Random layouts: up to 12 queries with gaps of 0-2 letters between their offsets, ranges of
/// a few letters, masks that touch or overlap, short motifs, 2 % ambiguity codes.
fn random_layout(rng: &mut Lcg) -> Case {
    let n_queries = 1 + rng.below(12) as usize;
    let mut queries = Vec::new();
    let mut offsets = Vec::new();
    let mut masks = Vec::new();
    let mut offset = rng.below(3) as i32;
    for _ in 0..n_queries {
        let len = match rng.below(5) {
            0 => rng.below(25) as usize,
            1 => 1 + rng.below(200) as usize,
            _ => 1 + rng.below(4000) as usize,
        };
        let motif: Vec<u8> = (0..1 + rng.below(6)).map(|_| rng.below(4) as u8).collect();
        let mut seq = Vec::with_capacity(len);
        while seq.len() < len {
            let r = rng.below(100);
            if r < 2 {
                seq.push(4 + rng.below(12) as u8);
            } else if r < 60 {
                seq.push(motif[seq.len() % motif.len()]);
            } else {
                seq.push(rng.below(4) as u8);
            }
        }
        let mut m = Vec::new();
        let mut cursor = 0usize;
        while cursor < len && rng.below(4) != 0 {
            let start = cursor + rng.below(30) as usize;
            let end = start + rng.below(25) as usize;
            if start >= len {
                break;
            }
            m.push(MaskedInterval {
                start,
                end: end.min(len),
            });
            cursor = if rng.below(3) == 0 { start } else { end };
        }
        offsets.push(offset);
        offset += len as i32 + 1 + rng.below(3) as i32;
        queries.push(seq);
        masks.push(m);
    }
    Case {
        queries,
        offsets,
        masks,
    }
}

// Random layouts (gaps between query offsets, touching and overlapping masks, many short
// ranges) partitioned in 1-20 parts with chunks of 1-5,000 positions: the original loop's
// tables.
#[test]
fn random_layouts_partitioned_equal_the_loop() {
    let mut rng = Lcg(0xBEEF_0001);
    let combos = combos();
    let mut words = 0usize;
    for round in 0..200 {
        let case = random_layout(&mut rng);
        let combo = combos[rng.below(combos.len() as u64) as usize];
        let reference = build(&case, combo, How::Serial);
        words += reference.hashtable.iter().filter(|&&c| c != 0).count();
        for _ in 0..4 {
            let parts = 1 + rng.below(20) as usize;
            let chunk = 1 + rng.below(5000) as usize;
            let got = build(&case, combo, How::Chunked(parts, chunk, false));
            assert!(
                got == reference,
                "round {round} {combo:?} parts {parts} chunk {chunk}"
            );
        }
    }
    assert!(words > 10_000, "{words}");
}

// On search pools of 2, 5 and 16 threads: random layouts in 1-24 parts and chunks of 1-3,000
// positions, then a query of 1.3 million letters and its reverse complement in the default
// chunks (five chunk boundaries, one inside a mask), one and two templates: the original
// loop's tables.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn pool_fill_in_default_chunks() {
    let mut rng = Lcg(0xBEEF_0002);
    let combos = combos();
    for threads in [2usize, 5, 16] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            for round in 0..30 {
                let case = random_layout(&mut rng);
                let combo = combos[rng.below(combos.len() as u64) as usize];
                let reference = build(&case, combo, How::Serial);
                let chunk = 1 + rng.below(3000) as usize;
                let parts = 1 + rng.below(24) as usize;
                let got = build(&case, combo, How::Chunked(parts, chunk, true));
                assert!(
                    got == reference,
                    "threads {threads} round {round} {combo:?} parts {parts} chunk {chunk}"
                );
            }
            let seq: Vec<u8> = (0..1_300_000)
                .map(|i| if i % 9973 < 3 { 14 } else { rng.below(4) as u8 })
                .collect();
            let reverse: Vec<u8> = seq
                .iter()
                .rev()
                .map(|&b| if b < 4 { 3 - b } else { b })
                .collect();
            let case = Case {
                queries: vec![seq, reverse],
                offsets: vec![0, 1_300_001],
                masks: vec![
                    vec![MaskedInterval {
                        start: 524_280,
                        end: 524_300,
                    }],
                    Vec::new(),
                ],
            };
            for combo in [
                (11usize, DiscWordType::TwoTemplates, 21usize),
                (12, DiscWordType::Coding, 16),
            ] {
                let reference = build(&case, combo, How::Serial);
                let (got, stats) =
                    build_stats(&case, combo, How::Chunked(threads, CHUNK_POSITIONS, true));
                assert_eq!(stats.chunks, 5, "threads {threads}");
                assert!(got == reference, "threads {threads} long {combo:?}");
            }
            Ok(())
        })
        .unwrap();
    }
}

/// Blocks of `chunk` letters, each a different repeat of period 1 to 4 (unmasked): each
/// chunk's words fall on a few cells, of one owner or two, which change from block to block.
fn repeat_blocks(chunk: usize, blocks: usize) -> Case {
    let mut seq = Vec::with_capacity(chunk * blocks);
    for b in 0..blocks {
        let period = 1 + b % 4;
        let motif: Vec<u8> = (0..period)
            .map(|k| ((b * 7 + k * 3 + b / 4) % 4) as u8)
            .collect();
        seq.extend((0..chunk).map(|i| motif[i % period]));
    }
    Case {
        queries: vec![seq],
        offsets: vec![0],
        masks: Vec::new(),
    }
}

// The buckets' memory has a bound that does not grow with the parts (`bucket_keep`): at most
// `6 * chunk + 68 * parts^2` entries per template, also when the chunks' words fall on the
// cells of one owner after another (repeat blocks, where buckets that kept their room grew to
// `2 * parts` times a chunk); and about twice a chunk per template on random letters. In the
// default chunks, 4 and 16 parts, one and two templates.
#[test]
fn bucket_memory_is_bounded() {
    let chunk = CHUNK_POSITIONS;
    let blocks = 20;
    let mut rng = Lcg(0xBEEF_0003);
    let random = Case {
        queries: vec![(0..chunk * blocks).map(|_| rng.below(4) as u8).collect()],
        offsets: vec![0],
        masks: Vec::new(),
    };
    let repeats = repeat_blocks(chunk, blocks);
    for (name, case) in [("random", &random), ("repeats", &repeats)] {
        for word_type in [DiscWordType::Coding, DiscWordType::TwoTemplates] {
            let templates = if word_type == DiscWordType::TwoTemplates {
                2
            } else {
                1
            };
            for parts in [4usize, 16] {
                let (_, stats) =
                    build_stats(case, (11, word_type, 18), How::Chunked(parts, chunk, false));
                let bound = templates * (6 * chunk + 68 * parts * parts);
                eprintln!(
                    "bucket peak {name} {word_type:?} parts {parts}: {} entries = {:.1} MiB (bound {:.1} MiB)",
                    stats.peak_entries,
                    stats.peak_entries as f64 * 8.0 / 1048576.0,
                    bound as f64 * 8.0 / 1048576.0
                );
                assert_eq!(stats.chunks, blocks);
                assert!(
                    stats.peak_entries <= bound,
                    "{name} {word_type:?} parts {parts}"
                );
                if name == "random" {
                    assert!(
                        stats.peak_entries <= templates * (2 * chunk + 64 * parts * parts),
                        "{name} {word_type:?} parts {parts}: {}",
                        stats.peak_entries
                    );
                }
            }
        }
    }
}

/// One query of 60,000 letters: random letters with ten runs of 600 letters, each a repeat of
/// period 1 to 5 (unmasked), so a few cells get hundreds of words (counters far above zero)
/// and most blocks get no repeated cell (counters that stay where they start).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn counter_case() -> Case {
    let mut rng = Lcg(0xC0DE_0C01);
    let mut seq: Vec<u8> = (0..60_000).map(|_| rng.below(4) as u8).collect();
    for run in 0..10usize {
        let period = 1 + run % 5;
        let motif: Vec<u8> = (0..period).map(|_| rng.below(4) as u8).collect();
        let at = 3_000 + run * 5_600;
        for (i, letter) in seq[at..at + 600].iter_mut().enumerate() {
            *letter = motif[i % period];
        }
    }
    Case {
        queries: vec![seq],
        offsets: vec![0],
        masks: Vec::new(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:781-783,799-801
// The owners' private counters, added into the counters after the last chunk, give the loop's
// `helper_array` and `helper_array2` on search pools of 2, 3, 5 and 16 threads (the default
// chunk and chunks of 1,000 positions, one and two templates, word sizes 11 and 12), also
// when the counters start where the loop's increments wrap them past u32::MAX (start values
// u32::MAX, u32::MAX - 6 .. u32::MAX - 0, 0 and others): the case has blocks whose counters
// wrap, blocks whose counters stay zero, and blocks no increment reaches in every class.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn private_counters_give_the_loops_counters() {
    let case = counter_case();
    let start = |h: usize| -> u32 {
        match h % 4 {
            0 => u32::MAX - (h % 7) as u32,
            1 => 0,
            2 => u32::MAX,
            _ => (h as u32).wrapping_mul(2_654_435_761),
        }
    };
    let combos = [
        (11usize, DiscWordType::Coding, 18usize),
        (11, DiscWordType::TwoTemplates, 21),
        (12, DiscWordType::TwoTemplates, 16),
        (12, DiscWordType::Optimal, 18),
    ];
    // The loop's tables from these starts, and the classes of blocks they have.
    let references: Vec<Built> = combos
        .iter()
        .map(|&combo| {
            let reference = build_from(&case, combo, How::Serial, Some(&start)).0;
            let approx = case.queries[0].len();
            let (_, bts) = compute_mb_pv_params(1 << (2 * combo.0), approx, false, combo.0);
            let (block_shift, _) = plan(1 << (2 * combo.0), bts, 2);
            let per_block = 1usize << (block_shift - K_COMPRESSION_FACTOR_SHIFT);
            let arrays: Vec<&[u32]> = if combo.1 == DiscWordType::TwoTemplates {
                vec![&reference.helper_array, &reference.helper_array2]
            } else {
                vec![&reference.helper_array]
            };
            for (t, got) in arrays.into_iter().enumerate() {
                let base = t * reference.helper_array.len();
                let (mut wrapped, mut zero, mut untouched_high) = (0, 0, 0);
                for (b, block) in got.chunks(per_block).enumerate() {
                    let starts: Vec<u32> =
                        (0..block.len()).map(|i| start(base + b * per_block + i)).collect();
                    if block.iter().zip(&starts).any(|(&v, &s)| v < s) {
                        wrapped += 1;
                    }
                    if block.iter().zip(&starts).all(|(&v, &s)| v == s) {
                        if starts.iter().all(|&s| s == 0) {
                            zero += 1;
                        } else if starts.iter().any(|&s| s >= u32::MAX - 6) {
                            untouched_high += 1;
                        }
                    }
                }
                assert!(
                    wrapped > 0 && zero > 0 && untouched_high > 0,
                    "{combo:?} template {t}: blocks wrapped {wrapped}, zero {zero}, untouched near the wrap {untouched_high}"
                );
            }
            reference
        })
        .collect();
    for threads in [2usize, 3, 5, 16] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            assert_eq!(parts_for(true), Some(threads));
            for (&combo, reference) in combos.iter().zip(&references) {
                let got = build_from(&case, combo, How::Parts(threads, true), Some(&start)).0;
                assert!(got == *reference, "threads {threads} {combo:?}");
                let (got, stats) = build_from(
                    &case,
                    combo,
                    How::Chunked(threads, 1000, true),
                    Some(&start),
                );
                assert!(
                    got == *reference,
                    "threads {threads} {combo:?} chunks of 1,000"
                );
                // The owners count into copies of their own (not all into one shared copy,
                // which is just as exact but brings back the false sharing): the case has
                // repeats in several blocks, so at least two owners (of each template) have
                // counts before the flush.
                let templates = 1 + usize::from(combo.1 == DiscWordType::TwoTemplates);
                assert!(
                    stats.counter_owners >= 2 * templates,
                    "threads {threads} {combo:?}: owners with counts {} (templates {templates})",
                    stats.counter_owners
                );
            }
            Ok(())
        })
        .unwrap();
    }
    // From zero, as `build_disc_mb_lookup` allocates them: the loop's counters, not all zero.
    for &combo in &combos {
        let reference = build(&case, combo, How::Serial);
        assert!(reference.helper_array.iter().any(|&c| c > 100), "{combo:?}");
        for parts in [2usize, 3, 5, 16] {
            assert!(build(&case, combo, How::Parts(parts, false)) == reference);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:785-786
// The prefault round (`Fill::prefault`, here without the size rule) on search pools of 2, 5
// and 16 threads and off a pool in 1, 3 and 64 parts, word sizes 11 and 12, one and two
// templates: a query of 400 letters (chunks of 97 positions and the default chunk) and one of
// 1.2 million letters (the default chunk, and chunks of 5,000 on the pool of 5). The tables are
// the loop's, and `check_prefault` (which every partitioned build of these tests runs) finds 0
// written over 0, one round, and the pages touched equal to the pages the first chunk wrote,
// array by array (`next_pos` is never touched). The short query leaves most pages of the cells
// untouched by both; with the long one in the default chunks every page of the cells is
// written but perhaps the end pages (with the allocator's 16-byte offset the last page holds 4
// cells), and each page touched is touched exactly once (one toucher per page). Every reader
// job of the later chunks ran without the marking code (`read::<false>`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn prefault_touches_the_pages_the_fill_writes() {
    let mut rng = Lcg(0xFA17_0001);
    let small = Case {
        queries: vec![(0..400).map(|_| rng.below(4) as u8).collect()],
        offsets: vec![0],
        masks: Vec::new(),
    };
    let long = Case {
        queries: vec![(0..1_200_000)
            .map(|i| if i % 7919 < 2 { 14 } else { rng.below(4) as u8 })
            .collect()],
        offsets: vec![0],
        masks: vec![vec![MaskedInterval {
            start: 600_000,
            end: 610_000,
        }]],
    };
    let combos = [
        (11usize, DiscWordType::Coding, 18usize),
        (12, DiscWordType::TwoTemplates, 21),
    ];
    let cell_pages = |b: &Built, k: usize| {
        let v = if k == 0 { &b.hashtable } else { &b.hashtable2 };
        pages_of(page_offset(v.as_ptr()), v.len())
    };
    let check = |case: &Case, combo, how: How, reference: &Built, long: bool| {
        let (parts, chunk) = match how {
            How::Parts(p, _) => (p, CHUNK_POSITIONS),
            How::Chunked(p, c, _) => (p, c),
            _ => unreachable!(),
        };
        let what = format!("{combo:?} parts {parts} chunk {chunk}");
        let (got, stats) = build_stats(case, combo, how);
        assert!(got == *reference, "{what}: tables differ");
        let touched = check_prefault(&got, &stats);
        let templates = 1 + usize::from(combo.1 == DiscWordType::TwoTemplates);
        for t in 0..templates {
            let pages = cell_pages(&got, t);
            if long {
                let n = touched[HASHTABLE + t].len();
                if chunk == CHUNK_POSITIONS {
                    assert!(n + 2 >= pages, "{what} template {t}: {n} of {pages} pages");
                }
                let writes = stats.touch_log.iter().filter(|w| w.array == HASHTABLE + t);
                assert_eq!(writes.count(), n, "{what} template {t}: one touch per page");
            } else {
                let n = touched[HASHTABLE + t].len();
                assert!(
                    n > 0 && n < pages / 2,
                    "{what} template {t}: {n} of {pages} pages"
                );
            }
        }
        assert!(!touched[PV_ARRAY].is_empty(), "{what}");
        assert_eq!(stats.prefault_rounds, 1, "{what}");
        assert!(stats.chunks > 1 || !long, "{what}");
        assert_eq!(stats.unmarked_reads, (stats.chunks - 1) * parts, "{what}");
    };
    for &combo in &combos {
        let small_reference = build(&small, combo, How::Serial);
        let long_reference = build(&long, combo, How::Serial);
        for parts in [1usize, 3, 64] {
            for chunk in [97usize, CHUNK_POSITIONS] {
                check(
                    &small,
                    combo,
                    How::Chunked(parts, chunk, false),
                    &small_reference,
                    false,
                );
            }
        }
        check(&long, combo, How::Parts(3, false), &long_reference, true);
        for threads in [2usize, 5, 16] {
            crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
                check(
                    &small,
                    combo,
                    How::Chunked(threads, 97, true),
                    &small_reference,
                    false,
                );
                check(
                    &small,
                    combo,
                    How::Parts(threads, true),
                    &small_reference,
                    false,
                );
                check(
                    &long,
                    combo,
                    How::Parts(threads, true),
                    &long_reference,
                    true,
                );
                if threads == 5 {
                    check(
                        &long,
                        combo,
                        How::Chunked(threads, 5000, true),
                        &long_reference,
                        true,
                    );
                }
                Ok(())
            })
            .unwrap();
        }
    }
}

// Each page of the cells (one and two templates, word sizes 11 and 12) and of the presence
// vector has one toucher, the owner of the block of its middle element, for 1, 2, 3, 16, 64
// and 255 parts and page offsets 0, 4 and 1,023 elements; the pages come in rising order.
#[test]
fn every_page_has_one_toucher() {
    for word_size in [11usize, 12] {
        let table_size = 1usize << (2 * word_size);
        for approx in [100usize, 100_000] {
            let (pv_size, bts) = compute_mb_pv_params(table_size, approx, false, word_size);
            for parts in [1usize, 2, 3, 16, 64, 255] {
                let (block_shift, owner) = plan(table_size, bts, parts);
                for offset in [0usize, 4, 1023] {
                    let offsets = [offset, (offset + 7) % 1024, (offset + 512) % 1024];
                    let lens = [table_size, table_size, pv_size];
                    let mut jobs: Vec<FaultJob> = (0..parts).map(|_| FaultJob::default()).collect();
                    page_touchers(&mut jobs, &owner, block_shift, bts, offsets, lens);
                    for k in [HASHTABLE, HASHTABLE + 1, PV_ARRAY] {
                        let pages = pages_of(offsets[k], lens[k]);
                        let mut seen = vec![0u8; pages];
                        for (j, job) in jobs.iter().enumerate() {
                            let list = &job.pages[k];
                            assert!(list.windows(2).all(|w| w[0] < w[1]));
                            for &page in list {
                                seen[page as usize] += 1;
                                let m = middle(offsets[k], lens[k], page as usize);
                                let cell = if k == PV_ARRAY { m << bts } else { m };
                                assert_eq!(owner[cell >> block_shift] as usize, j);
                            }
                        }
                        assert!(
                            seen.iter().all(|&n| n == 1),
                            "word size {word_size} parts {parts} offset {offset} array {k}"
                        );
                    }
                }
            }
        }
    }
}

// `fill`'s size rule (`PREFAULT_MIN_WORDS_PER_PAGE`): the round runs only when the first chunk
// has at least `k` positions per page of the cells. Word size 11 (4,096 or 4,097 pages), one
// and two templates, 3 parts: a query of `k * 4,097 + 100` letters gets one round, one of
// `k * 4,096 - 100` none, and so does the longer query in chunks of `k * 4,096 - 100`; with no
// round nothing is touched and every reader job runs without the marking code. The tables
// are the loop's either way.
#[test]
fn the_round_runs_only_when_the_first_chunk_is_dense() {
    let k = PREFAULT_MIN_WORDS_PER_PAGE;
    assert!(k >= 1);
    let mut rng = Lcg(0xFA17_0002);
    let letters = |n: usize, rng: &mut Lcg| -> Case {
        Case {
            queries: vec![(0..n).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        }
    };
    let dense = letters(k * 4097 + 100, &mut rng);
    let sparse = letters(k * 4096 - 100, &mut rng);
    for combo in [
        (11usize, DiscWordType::Coding, 18usize),
        (11, DiscWordType::TwoTemplates, 16),
    ] {
        for (case, chunk, rounds) in [
            (&dense, CHUNK_POSITIONS, 1usize),
            (&sparse, CHUNK_POSITIONS, 0),
            (&dense, k * 4096 - 100, 0),
        ] {
            let reference = build(case, combo, How::Serial);
            let rule = PrefaultRule {
                min_words_per_page: k,
                region_shift: 0,
            };
            let (got, stats) = build_stats(case, combo, How::Rule(3, chunk, rule));
            let what = format!("{combo:?} chunk {chunk}");
            assert!(got == reference, "{what}: tables differ");
            assert_eq!(stats.prefault_rounds, rounds, "{what}");
            if rounds == 0 {
                assert_eq!(stats.touches, 0, "{what}");
                assert_eq!(stats.unmarked_reads, stats.chunks * 3, "{what}");
            }
        }
    }
}

// `fill`'s regions (`PREFAULT_RULE.region_shift`, 512 pages = 2 MiB) and smaller ones: the round
// touches the first chunk's pages only in the regions that lie within the array's pages and
// that the first chunk wrote in full (`check_prefault`, run by `build_from`, compares the
// exact set). 300,000 random letters write every page of a word size 11 table, so all of its
// whole regions are touched (whole multiples of 512 pages; the presence vector, 128 KiB,
// holds none); 400 letters complete no region, so nothing is touched; with regions of 8
// pages, 5,000 letters (word size 11, one template) and 40,000 (12, two templates) complete
// some regions and not others. The tables are the loop's.
#[test]
fn the_round_touches_whole_regions_only() {
    let mut rng = Lcg(0xFA17_0004);
    let letters = |n: usize, rng: &mut Lcg| -> Case {
        Case {
            queries: vec![(0..n).map(|_| rng.below(4) as u8).collect()],
            offsets: vec![0],
            masks: Vec::new(),
        }
    };
    let (big, small, medium) = (
        letters(300_000, &mut rng),
        letters(400, &mut rng),
        letters(5_000, &mut rng),
    );
    let rule = |shift: u32| PrefaultRule {
        min_words_per_page: 0,
        region_shift: shift,
    };
    let coding = (11usize, DiscWordType::Coding, 18usize);
    let shift = PREFAULT_RULE.region_shift;
    let (got, stats) = build_stats(&big, coding, How::Rule(3, CHUNK_POSITIONS, rule(shift)));
    assert!(got == build(&big, coding, How::Serial));
    let pages = pages_of(page_offset(got.hashtable.as_ptr()), got.hashtable.len());
    let touched = check_prefault(&got, &stats);
    let cells = touched[HASHTABLE].len();
    assert!(
        cells + 2 * 512 >= pages && cells.is_multiple_of(512) && cells > 0,
        "{cells} of {pages}"
    );
    assert!(touched[PV_ARRAY].is_empty());
    let (got, stats) = build_stats(&small, coding, How::Rule(3, CHUNK_POSITIONS, rule(shift)));
    assert!(got == build(&small, coding, How::Serial));
    assert_eq!(stats.prefault_rounds, 1);
    assert_eq!(stats.touches, 0);
    let medium12 = letters(40_000, &mut rng);
    for (combo, medium) in [
        (coding, &medium),
        ((12, DiscWordType::TwoTemplates, 21), &medium12),
    ] {
        let reference = build(medium, combo, How::Serial);
        let (got, stats) = build_stats(medium, combo, How::Rule(3, CHUNK_POSITIONS, rule(3)));
        assert!(got == reference, "{combo:?}");
        let touched = check_prefault(&got, &stats)[HASHTABLE].len();
        let first = stats.first_pages[HASHTABLE].len();
        assert!(
            touched > 0 && touched < first,
            "{combo:?}: touched {touched} of the first chunk's {first} pages"
        );
    }
}

// The round at other page offsets of the arrays than the allocator's: `hashtable`,
// `hashtable2` and `pv_array` as slices of page-aligned buffers starting 0, 4 and 1,023
// elements into a page (one and two templates, word size 11, 3 parts, a query of 30,000
// letters in chunks of 7,000): the tables are the loop's, and `check_prefault_arrays` finds 0
// over 0, one round and exactly the first chunk's pages.
#[test]
fn the_round_at_any_page_offset() {
    let mut rng = Lcg(0xFA17_0003);
    let case = Case {
        queries: vec![(0..30_000).map(|_| rng.below(4) as u8).collect()],
        offsets: vec![0],
        masks: Vec::new(),
    };
    let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
    /// A zeroed buffer with room for `len` elements from `offset` elements after a page start;
    /// returns it and that start.
    fn aligned(len: usize, offset: usize) -> (Vec<u32>, usize) {
        let v = vec![0u32; len + offset + 1024];
        let skip = (PAGE_BYTES - v.as_ptr() as usize % PAGE_BYTES) % PAGE_BYTES / 4;
        (v, skip + offset)
    }
    for combo in [
        (11usize, DiscWordType::Coding, 18usize),
        (11, DiscWordType::TwoTemplates, 21),
    ] {
        let reference = build(&case, combo, How::Serial);
        let two = combo.1 == DiscWordType::TwoTemplates;
        let (_, bts) = compute_mb_pv_params(1 << 22, 30_000, false, 11);
        for offset in [0usize, 4, 1023] {
            let mut got = reference.clone();
            let (mut cells, c0) = aligned(got.hashtable.len(), offset);
            let (mut cells2, c2) = aligned(got.hashtable2.len(), (offset + 7) % 1024);
            let (mut pv, p0) = aligned(got.pv_array.len(), (offset + 512) % 1024);
            got.next_pos.fill(0);
            got.helper_array.fill(0);
            got.next_pos2.fill(0);
            got.helper_array2.fill(0);
            let (n1, n2, np) = (
                got.hashtable.len(),
                got.hashtable2.len(),
                got.pv_array.len(),
            );
            let stats = fill_with(
                3,
                &words_of(&case, &query_refs, combo),
                DiscTables {
                    hashtable: &mut cells[c0..c0 + n1],
                    next_pos: &mut got.next_pos,
                    helper_array: &mut got.helper_array,
                    hashtable2: &mut cells2[c2..c2 + n2],
                    next_pos2: &mut got.next_pos2,
                    helper_array2: &mut got.helper_array2,
                    pv_array: &mut pv[p0..p0 + np],
                    pv_array_bts: bts,
                },
                false,
                7_000,
                EVERY_PAGE,
                &|_| {},
            );
            let what = format!("{combo:?} offset {offset}");
            assert_eq!(page_offset(cells[c0..].as_ptr()), offset, "{what}");
            assert_eq!(&cells[c0..c0 + n1], &reference.hashtable[..], "{what}");
            assert_eq!(&cells2[c2..c2 + n2], &reference.hashtable2[..], "{what}");
            assert_eq!(&pv[p0..p0 + np], &reference.pv_array[..], "{what}");
            assert_eq!(got.next_pos, reference.next_pos, "{what}");
            assert_eq!(got.next_pos2, reference.next_pos2, "{what}");
            got.longest_chain = disc_longest_chain(&got.helper_array, &got.helper_array2, two);
            assert_eq!(got.helper_array, reference.helper_array, "{what}");
            assert_eq!(got.longest_chain, reference.longest_chain, "{what}");
            assert_eq!(stats.prefault_rounds, 1, "{what}");
            check_prefault_arrays(
                [&cells[c0..c0 + n1], &cells2[c2..c2 + n2], &pv[p0..p0 + np]],
                &stats,
            );
        }
    }
}

// The bucket headers and the owners' copies of the counters lie on 64-byte lines of their own:
// a `Bucket` is one line (size and alignment 64), so the P x P headers of the grid are 64
// bytes apart on line boundaries; a `CounterLine` is one line, and each owner's copy starts on
// a line boundary and spans whole lines (also for a number of counters that is no multiple of
// 16).
#[test]
fn bucket_headers_and_counter_copies_have_lines_of_their_own() {
    assert_eq!(std::mem::align_of::<Bucket>(), 64);
    assert_eq!(std::mem::size_of::<Bucket>(), 64);
    assert_eq!(std::mem::align_of::<CounterLine>(), 64);
    assert_eq!(std::mem::size_of::<CounterLine>(), 64);
    // The prefault's per-job state and the readers' page marks too.
    assert_eq!(std::mem::align_of::<FaultJob>(), 64);
    assert_eq!(std::mem::size_of::<FaultJob>() % 64, 0);
    assert_eq!(std::mem::align_of::<BitLine>(), 64);
    assert_eq!(std::mem::size_of::<BitLine>(), 64);
    for (pages, parts) in [(1usize, 1usize), (4097, 16), (16385, 3), (33, 255)] {
        let mut marks = PageBits::new(pages, parts);
        let shared = marks.shared();
        assert!(shared.stride.is_multiple_of(8) && shared.stride * 64 >= pages);
        for r in 0..parts {
            assert_eq!(shared.ptr.wrapping_add(r * shared.stride) as usize % 64, 0);
        }
    }
    for parts in [1usize, 2, 3, 4, 16] {
        let grid: Vec<Bucket> = (0..parts * parts).map(|_| Bucket::default()).collect();
        for (i, bucket) in grid.iter().enumerate() {
            let at = bucket as *const Bucket as usize;
            assert_eq!(at % 64, 0, "parts {parts} bucket {i}");
            assert_eq!(
                at - grid.as_ptr() as usize,
                64 * i,
                "parts {parts} bucket {i}"
            );
        }
        for len in [2048usize, 8192, 100, 1] {
            let mut counters = OwnerCounters::new(len, parts);
            let shared = counters.shared();
            assert_eq!(shared.stride % 16, 0);
            assert!(shared.stride >= len && shared.stride < len + 16);
            for o in 0..parts {
                let copy = shared.of(o);
                assert_eq!(
                    copy.ptr as usize % 64,
                    0,
                    "parts {parts} len {len} owner {o}"
                );
                assert_eq!(copy.len, len);
                assert!(
                    copy.ptr as usize + 4 * shared.stride
                        <= shared.all.ptr as usize + 4 * shared.all.len
                );
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_DCCSR): the CSR built on the pool (`csr_on_pool`, and LOSAT_X_DCSEEDLEAN's
// `csr_on_pool_steps`), in chunks of any size, is the serial build's and lists the chains of
// the original loop (`x_dc_csr::x_assert_runs_are_the_chains` compares every cell's run with
// its chain, the presence vector, the counters, the longest chain and the two builds; it
// panics on a difference).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn check_csr(case: &Case, combo: (usize, DiscWordType, usize), parts: usize, chunk: usize) {
    let reference = build(case, combo, How::Serial);
    let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
    let words = words_of(case, &query_refs, combo);
    let cells = reference.hashtable.len();
    let approx: usize = case.queries.iter().map(Vec::len).sum();
    let (_, pv_array_bts) = compute_mb_pv_params(cells, approx, false, combo.0);
    // EXPERIMENT (LOSAT_X_DCSEEDLEAN): `csr_on_pool` and `csr_on_pool_steps`, the latter with and
    // without LOSAT_X_DCQSAHEAD's `qs` ahead.
    for (steps, qs_ahead) in [(false, false), (true, false), (true, true)] {
        let build = |pv: &mut [PvArrayType], helper: &mut [u32], x_qs: Option<Vec<u32>>| {
            if steps {
                csr_on_pool_steps(
                    parts,
                    &words,
                    cells,
                    pv,
                    pv_array_bts,
                    helper,
                    chunk,
                    Vec::new(),
                    x_qs,
                    qs_ahead,
                )
            } else {
                csr_on_pool(
                    parts,
                    &words,
                    cells,
                    pv,
                    pv_array_bts,
                    helper,
                    chunk,
                    Vec::new(),
                    x_qs,
                )
            }
        };
        let mut pv = vec![0; reference.pv_array.len()];
        let mut helper = vec![0; reference.helper_array.len()];
        let csr = build(&mut pv, &mut helper, None);
        // EXPERIMENT (LOSAT_X_CSRQSKEEP): on a kept `qs` of stale values, the same CSR.
        super::super::x_dc_csr::x_assert_same_on_a_kept_qs(
            &csr,
            |x_qs| {
                let mut pv = vec![0; reference.pv_array.len()];
                let mut helper = vec![0; reference.helper_array.len()];
                build(&mut pv, &mut helper, x_qs)
            },
            &format!("parts {parts} chunk {chunk} steps {steps} qs_ahead {qs_ahead}"),
        );
        super::super::x_dc_csr::x_assert_runs_are_the_chains(
            &words,
            &csr,
            true,
            &pv,
            &helper,
            pv_array_bts,
            &reference.hashtable,
            &reference.next_pos,
            &reference.pv_array,
            &reference.helper_array,
            reference.longest_chain,
        );
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn csr_on_a_pool_lists_the_chains() {
    let mut rng = Lcg(0xC5A0_0001);
    let one_template: Vec<_> = combos()
        .into_iter()
        .filter(|c| c.1 != DiscWordType::TwoTemplates)
        .collect();
    for threads in [2usize, 3, 5] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            for round in 0..16 {
                let case = if round % 2 == 0 {
                    random_case(&mut rng)
                } else {
                    random_layout(&mut rng)
                };
                let combo = one_template[rng.below(one_template.len() as u64) as usize];
                let chunk = [1usize, 97, 5000, CHUNK_POSITIONS][round % 4];
                check_csr(&case, combo, threads, chunk);
            }
            // `csr` takes the pool there, and builds on this thread off it.
            let case = random_case(&mut rng);
            let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
            let words = words_of(&case, &query_refs, (11, DiscWordType::Coding, 16));
            let mut pv = vec![0; 1 << 15];
            let mut helper = vec![0; (1 << 22) / 2048];
            let (_, on_pool) = csr(
                threads,
                &words,
                1 << 22,
                &mut pv,
                7,
                &mut helper,
                Vec::new(),
                None,
            );
            assert!(on_pool, "threads {threads}");
            Ok(())
        })
        .unwrap();
    }
    let case = random_case(&mut rng);
    let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
    let words = words_of(&case, &query_refs, (11, DiscWordType::Coding, 16));
    let mut pv = vec![0; 1 << 15];
    let mut helper = vec![0; (1 << 22) / 2048];
    let (_, on_pool) = csr(
        4,
        &words,
        1 << 22,
        &mut pv,
        7,
        &mut helper,
        Vec::new(),
        None,
    );
    assert!(!on_pool, "off a pool");
}

// EXPERIMENT (LOSAT_X_DCQSAHEAD): one query of 0 to 40 words, of random letters and of one letter
// repeated (every word in one cell, so in one owner's buckets), in one chunk and in chunks of 7
// positions: buckets of every length around `QS_AHEAD` and 16, built by `csr_on_pool` and by
// `csr_on_pool_steps` with and without `qs` ahead (`check_csr`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn csr_on_a_pool_at_every_short_length() {
    let mut rng = Lcg(0xC5A0_0003);
    let combo = (11, DiscWordType::Coding, 16);
    for threads in [2usize, 3] {
        crate::utils::threading::with_search_pool(threads, "blastn", |_pool| {
            for n in 0..=40usize {
                for repeated in [false, true] {
                    let queries = vec![(0..combo.2 - 1 + n)
                        .map(|_| if repeated { 0 } else { rng.below(4) as u8 })
                        .collect::<Vec<u8>>()];
                    let case = Case {
                        queries,
                        offsets: vec![0],
                        masks: vec![Vec::new()],
                    };
                    for chunk in [7usize, 5000] {
                        check_csr(&case, combo, threads, chunk);
                    }
                }
            }
            Ok(())
        })
        .unwrap();
    }
}

// EXPERIMENT (LOSAT_X_DCCSR): a table of more positions than a quarter of its cells, built on
// a pool, is laid out as CSR (with the switch on) and answers `count_hits`, `contains_hit`
// and `longest_chain` as the chains of the original loop; its runs are the chains.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn a_large_table_on_a_pool_answers_from_its_csr() {
    let mut rng = Lcg(0xC5A0_0002);
    let seq: Vec<u8> = (0..1_200_000)
        .map(|i| {
            if i % 7919 < 2 {
                14
            } else if i % 3 == 0 {
                (i / 5 % 4) as u8
            } else {
                rng.below(4) as u8
            }
        })
        .collect();
    let case = Case {
        queries: vec![seq.clone(), seq[..300_000].to_vec()],
        offsets: vec![0, 1_200_001],
        masks: vec![
            vec![MaskedInterval {
                start: 524_280,
                end: 524_300,
            }],
            Vec::new(),
        ],
    };
    let combo = (11usize, DiscWordType::Coding, 16usize);
    let reference = build(&case, combo, How::Serial);
    let query_refs: Vec<&[u8]> = case.queries.iter().map(Vec::as_slice).collect();
    let approx: usize = case.queries.iter().map(Vec::len).sum();
    crate::utils::threading::with_search_pool(4, "blastn", |_pool| {
        let lookup = build_disc_mb_lookup(
            &query_refs,
            &case.offsets,
            combo.0,
            combo.2,
            combo.1,
            &case.masks,
            approx,
        );
        let disc = lookup.disc().expect("a discontiguous table");
        assert_eq!(disc.x_csr().is_some(), super::super::x_dc_csr::on());
        assert_eq!(lookup.longest_chain(), reference.longest_chain);
        let mut chain = Vec::new();
        for cell in (0..reference.hashtable.len()).step_by(3) {
            chain.clear();
            let mut q = reference.hashtable[cell];
            while q != 0 {
                chain.push(q);
                q = reference.next_pos[q as usize];
            }
            assert_eq!(lookup.count_hits(cell as u64), chain.len(), "cell {cell}");
            for &q in chain.iter().take(3) {
                assert!(lookup.contains_hit(cell as u64, q), "cell {cell} {q}");
            }
            assert!(!lookup.contains_hit(cell as u64, 1_500_002), "cell {cell}");
            if let Some(csr) = disc.x_csr() {
                assert_eq!(csr.run(cell), chain.as_slice(), "cell {cell}");
            }
        }
        Ok(())
    })
    .unwrap();
}
