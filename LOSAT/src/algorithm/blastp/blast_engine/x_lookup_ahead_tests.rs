use super::*;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastp::x_query_blocks::XBlockSubjects;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_assert_same_flat_tables_as;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::XFlatCompressedLookup;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::tblastx::lookup::x_pair_lut::x_assert_same_contexts_as;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::tblastx::lookup::{x_assert_same_table_as, x_same_table};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use std::path::PathBuf;

// No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
/// The first difference between two query frame sets (every field), or `None`.
fn x_frames_difference(got: &[Vec<QueryFrame>], want: &[Vec<QueryFrame>]) -> Option<String> {
    if got.len() != want.len() {
        return Some(format!(
            "{} queries' frames built ahead, {} inline",
            got.len(),
            want.len()
        ));
    }
    for (query, (got, want)) in got.iter().zip(want).enumerate() {
        let same = got.len() == want.len()
            && got.iter().zip(want).all(|(a, b)| {
                a.frame == b.frame
                    && a.aa_seq == b.aa_seq
                    && a.aa_seq_nomask == b.aa_seq_nomask
                    && a.aa_len == b.aa_len
                    && a.orig_len == b.orig_len
                    && a.seg_masks == b.seg_masks
            });
        if !same {
            return Some(format!("the frames of query {query} differ"));
        }
    }
    None
}

// No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
/// Whether two flat tables are equal, field for field and slot for slot
/// (`x_assert_same_flat_tables_as` names the first difference).
fn x_same_flat_tables(a: &XFlatCompressedLookup, b: &XFlatCompressedLookup) -> bool {
    a.word_length == b.word_length
        && a.compressed_alphabet_size == b.compressed_alphabet_size
        && a.reciprocal_alphabet_size == b.reciprocal_alphabet_size
        && a.longest_chain == b.longest_chain
        && a.pv_array_bts == b.pv_array_bts
        && a.scaled_compress_table_u8 == b.scaled_compress_table_u8
        && a.pv == b.pv
        && a.starts == b.starts
        && a.hits == b.hits
}

// No NCBI counterpart: a check of the unit tests; it changes no value NCBI computes.
/// Whether two context sets are equal, every field, the Karlin blocks bit for bit
/// (`x_assert_same_contexts_as` names the first difference).
fn x_same_contexts(a: &[QueryContext], b: &[QueryContext]) -> bool {
    let bits = |context: &QueryContext| {
        let p = &context.karlin_params;
        [
            p.lambda.to_bits(),
            p.k.to_bits(),
            p.h.to_bits(),
            p.alpha.to_bits(),
            p.beta.to_bits(),
        ]
    };
    a.len() == b.len()
        && a.iter().zip(b).all(|(a, b)| {
            a.q_idx == b.q_idx
                && a.frame == b.frame
                && a.aa_seq == b.aa_seq
                && a.aa_seq_nomask == b.aa_seq_nomask
                && a.aa_len == b.aa_len
                && a.orig_len == b.orig_len
                && a.frame_base == b.frame_base
                && a.is_valid == b.is_valid
                && bits(a) == bits(b)
        })
}

// Nested parallel work (two levels of `par_iter`) of block `block`, as a block's search and a
// build run it on the pool.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn nested_work(block: usize) -> u64 {
    use rayon::prelude::*;
    (0..48u64)
        .into_par_iter()
        .map(|x| {
            (0..x + block as u64)
                .into_par_iter()
                .map(|y| std::hint::black_box(y * y))
                .sum::<u64>()
        })
        .sum()
}

/// How `fake_loop` builds a block's "table": inline, ahead in one task
/// (`x_run_with_build_ahead`), or ahead in two parts (`x_run_with_split_build_ahead`, item SC2d)
/// with the vectors given after the block's nested work (part 1 may end before or after), before
/// part 1 starts (the vectors come first), or once part 1 has ended (part 1 comes first, mostly).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[derive(Clone, Copy, Debug, PartialEq)]
enum XFakeAhead {
    Inline,
    Whole,
    Split,
    SplitVectorsFirst,
    SplitPartFirst,
}

/// Sets the flag when dropped (also on a panic).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XSetOnDrop<'a>(&'a std::sync::atomic::AtomicBool);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl Drop for XSetOnDrop<'_> {
    fn drop(&mut self) {
        self.0.store(true, Ordering::Release);
    }
}

/// A block loop with fake blocks, as `search` drives the blocks: block b's "table" is `build(b)`
/// (an `Err` for a failing or refused build; it may panic), made by its search inline, or ahead,
/// during block b - 1, as `ahead` says (in two parts: `build(b)` in part 1, nested work in part
/// 2); the search then starts the next block's build, runs nested parallel work (its "scan"),
/// gives the vectors to part 2 of a split build, fails where `fails` says, and writes its row to
/// `rows`. Returns the loop's outcome (the first error ends it); a panic leaves the loop.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn fake_loop(
    threads: usize,
    blocks: usize,
    ahead: XFakeAhead,
    build: &(dyn Fn(usize) -> Result<u64, String> + Sync),
    fails: Option<usize>,
    rows: &Mutex<Vec<String>>,
) -> Result<(), String> {
    use std::sync::atomic::AtomicBool;
    crate::utils::threading::with_search_pool(threads, "blastp", |_pool| {
        let mut next: Option<XBuilt<u64>> = None;
        for block in 0..blocks {
            let prebuilt = x_take_built(next.take());
            let took = prebuilt.is_some();
            let part_ended = AtomicBool::new(false);
            // `spawn`: the build's start, and for a split build the vectors' hand-over.
            let search = |spawn: Option<(&dyn Fn(), Option<&dyn Fn()>)>| -> Result<(), String> {
                let table = match prebuilt {
                    Some(table) => table,
                    None => build(block)?,
                };
                let (start, vectors) = spawn.unzip();
                let vectors = vectors.flatten();
                if let (XFakeAhead::SplitVectorsFirst, Some(vectors)) = (ahead, vectors) {
                    vectors();
                }
                if let Some(start) = start {
                    start();
                }
                if ahead == XFakeAhead::SplitPartFirst && vectors.is_some() {
                    while !part_ended.load(Ordering::Acquire) {
                        if !matches!(rayon::yield_now(), Some(rayon::Yield::Executed)) {
                            std::thread::yield_now();
                        }
                    }
                }
                let work = nested_work(block);
                if let (true, Some(vectors)) = (ahead != XFakeAhead::SplitVectorsFirst, vectors) {
                    vectors();
                }
                if fails == Some(block) {
                    return Err(format!("block {block} failed"));
                }
                rows.lock()
                    .unwrap()
                    .push(format!("block {block} table {table} work {work}"));
                Ok(())
            };
            let outcome = match ahead {
                _ if block + 1 == blocks => search(None),
                XFakeAhead::Inline => search(None),
                XFakeAhead::Whole => {
                    let (outcome, built) = x_run_with_build_ahead(
                        move || {
                            nested_work(block + 1);
                            build(block + 1).ok()
                        },
                        |spawn| search(Some((spawn, None))),
                    );
                    next = built;
                    outcome
                }
                XFakeAhead::Split | XFakeAhead::SplitVectorsFirst | XFakeAhead::SplitPartFirst => {
                    let part_ended = &part_ended;
                    let (outcome, built) = x_run_with_split_build_ahead(
                        move || {
                            let _ended = XSetOnDrop(part_ended);
                            nested_work(block + 1);
                            build(block + 1).ok()
                        },
                        |table| {
                            nested_work(block + 2);
                            Some(table)
                        },
                        |start, vectors| search(Some((start, Some(vectors)))),
                    );
                    next = built;
                    outcome
                }
            };
            assert!(!took || block > 0, "the first block is never built ahead");
            if let Err(error) = outcome {
                return Ok(Err(error));
            }
        }
        Ok(Ok(()))
    })
    .expect("pool")
}

/// The loop of `fake_loop` without and with the builds ahead (whole and split, item SC2d): the
/// rows written, the outcome, and a panic's message (`None` without one); all runs must agree.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn same_loops(
    threads: usize,
    blocks: usize,
    build: &(dyn Fn(usize) -> Result<u64, String> + Sync),
    fails: Option<usize>,
) -> (Vec<String>, Result<(), String>, Option<String>) {
    let run = |ahead: XFakeAhead| {
        let rows = Mutex::new(Vec::new());
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            fake_loop(threads, blocks, ahead, build, fails, &rows)
        }));
        let rows = rows.into_inner().unwrap();
        match outcome {
            Ok(outcome) => (rows, outcome, None),
            Err(payload) => {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| payload.downcast_ref::<&str>().map(|text| text.to_string()))
                    .unwrap_or_default();
                (rows, Ok(()), Some(message))
            }
        }
    };
    let off = run(XFakeAhead::Inline);
    for ahead in [
        XFakeAhead::Whole,
        XFakeAhead::Split,
        XFakeAhead::SplitVectorsFirst,
        XFakeAhead::SplitPartFirst,
    ] {
        let on = run(ahead);
        assert_eq!(
            on, off,
            "threads {threads}, {blocks} blocks, fails {fails:?}, {ahead:?}"
        );
    }
    off
}

// Every block after the first takes the values built ahead during the block before it; the rows
// are those of the loop without the build ahead, at 2, 3 and 7 threads, with one block (nothing
// built ahead), two and six; nested parallel work in the searches and the builds ends.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn blocks_take_the_values_built_ahead_and_write_the_same_rows() {
    let build = |block: usize| -> Result<u64, String> { Ok(10 * block as u64 + 7) };
    for threads in [2, 3, 7] {
        for blocks in [1, 2, 6] {
            let (rows, outcome, panic) = same_loops(threads, blocks, &build, None);
            assert_eq!(rows.len(), blocks);
            assert_eq!(outcome, Ok(()));
            assert_eq!(panic, None);
        }
    }
}

// A build that fails (an error, a refusal) gives no values ahead; its block builds inline and its
// error ends the loop where it ends without the build ahead: after the rows of the blocks before
// it, before its own. The first block's error, and the last block's, too.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn a_failing_build_ahead_fails_in_its_block_after_the_rows_before_it() {
    for failing in [0, 1, 3, 5] {
        let build = move |block: usize| -> Result<u64, String> {
            if block == failing {
                Err(format!("build {block} refused"))
            } else {
                Ok(block as u64)
            }
        };
        for threads in [2, 7] {
            let (rows, outcome, panic) = same_loops(threads, 6, &build, None);
            assert_eq!(rows.len(), failing);
            assert_eq!(outcome, Err(format!("build {failing} refused")));
            assert_eq!(panic, None);
        }
    }
}

// A build that panics: its panic is kept and resumed where its block's search would start, after
// the rows of the blocks before it, with the same payload as the inline build's panic.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn a_panic_of_a_build_ahead_surfaces_where_its_block_would_raise_it() {
    for panicking in [1, 4] {
        let build = move |block: usize| -> Result<u64, String> {
            if block == panicking {
                panic!("build {block} panicked");
            }
            Ok(block as u64)
        };
        for threads in [2, 3] {
            let (rows, outcome, panic) = same_loops(threads, 5, &build, None);
            assert_eq!(rows.len(), panicking);
            assert_eq!(outcome, Ok(()));
            assert_eq!(panic, Some(format!("build {panicking} panicked")));
        }
    }
}

// A block whose search fails ends the loop with its error; the values built ahead for the next
// block, an error of that build or its panic are dropped (no panic is resumed).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn a_failing_block_drops_what_was_built_ahead_for_the_next() {
    let fine = |block: usize| -> Result<u64, String> { Ok(block as u64) };
    let failing = |block: usize| -> Result<u64, String> {
        if block == 3 {
            Err("build 3 refused".to_string())
        } else {
            Ok(block as u64)
        }
    };
    let panicking = |block: usize| -> Result<u64, String> {
        if block == 3 {
            panic!("build 3 panicked");
        }
        Ok(block as u64)
    };
    let builds: [&(dyn Fn(usize) -> Result<u64, String> + Sync); 3] = [&fine, &failing, &panicking];
    for build in builds {
        for threads in [2, 7] {
            let (rows, outcome, panic) = same_loops(threads, 5, build, Some(2));
            assert_eq!(rows.len(), 2);
            assert_eq!(outcome, Err("block 2 failed".to_string()));
            assert_eq!(panic, None);
        }
    }
}

// `x_run_with_build_ahead` spawns the build once, on the first call of `spawn` only, and gives no
// outcome when `spawn` is not called.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_build_ahead_is_spawned_once_when_asked() {
    use std::sync::atomic::AtomicUsize;
    crate::utils::threading::with_search_pool(3, "blastp", |_pool| {
        let builds = AtomicUsize::new(0);
        let (result, built) = x_run_with_build_ahead(
            || Some(builds.fetch_add(1, Ordering::Relaxed)),
            |spawn| {
                spawn();
                spawn();
                5
            },
        );
        assert_eq!(result, 5);
        assert!(matches!(built, Some(XBuilt::Values(0))));
        assert_eq!(builds.load(Ordering::Relaxed), 1);
        let (_, built) = x_run_with_build_ahead(|| Some(1u8), |_spawn| ());
        assert!(built.is_none());
        let (_, built) = x_run_with_build_ahead(|| None::<u8>, |spawn| spawn());
        assert!(matches!(built, Some(XBuilt::Dropped)));
        Ok(())
    })
    .expect("pool");
}

// Item SC2d: the hand-off of a split build ahead starts part 2 exactly once, on the side that
// comes second, in either order: part 1's values first (the vectors take them), or the vectors
// first (part 1's task goes on); a second call of the vectors, and the vectors after part 1 gave
// no values, start nothing.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_hand_off_starts_part_two_once_on_the_side_that_comes_second() {
    let slot = Mutex::new(XHandoff::Empty);
    assert_eq!(XHandoff::part(&slot, 7u8), None, "part 1 first: kept");
    assert_eq!(
        XHandoff::vectors(&slot),
        Some(7),
        "then the vectors start part 2"
    );
    assert_eq!(
        XHandoff::vectors(&slot),
        None,
        "a second call starts nothing"
    );
    let slot = Mutex::new(XHandoff::Empty);
    assert_eq!(XHandoff::vectors(&slot), None, "the vectors first: kept");
    assert_eq!(XHandoff::vectors(&slot), None);
    assert_eq!(
        XHandoff::part(&slot, 9u8),
        Some(9),
        "then part 1's task starts part 2"
    );
    assert_eq!(XHandoff::vectors(&slot), None);
    let slot = Mutex::new(XHandoff::<u8>::Empty);
    XHandoff::close(&slot);
    assert_eq!(XHandoff::vectors(&slot), None, "part 1 gave no values");
    assert!(matches!(slot.into_inner().unwrap(), XHandoff::Closed));
}

// Item SC2d: `x_run_with_split_build_ahead` runs each part once, with part 1's values, in either
// order of the hand-over; gives no values where part 1 or part 2 gives none or the vectors never
// come (part 1's values dropped), keeps a panic of either part, and gives no outcome when part 1
// is not started.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_split_build_ahead_runs_each_part_once_and_keeps_its_outcome() {
    use std::sync::atomic::AtomicUsize;
    crate::utils::threading::with_search_pool(3, "blastp", |_pool| {
        for vectors_first in [false, true] {
            let (firsts, seconds) = (AtomicUsize::new(0), AtomicUsize::new(0));
            let (result, built) = x_run_with_split_build_ahead(
                || Some(10 + firsts.fetch_add(1, Ordering::Relaxed)),
                |part| Some(part * 100 + seconds.fetch_add(1, Ordering::Relaxed)),
                |start, vectors| {
                    if vectors_first {
                        vectors();
                    }
                    start();
                    start();
                    nested_work(3);
                    vectors();
                    vectors();
                    5
                },
            );
            assert_eq!(result, 5);
            assert!(
                matches!(built, Some(XBuilt::Values(1000))),
                "{vectors_first}"
            );
            assert_eq!(firsts.load(Ordering::Relaxed), 1);
            assert_eq!(seconds.load(Ordering::Relaxed), 1);
        }
        let both = |start: &dyn Fn(), vectors: &dyn Fn()| {
            start();
            vectors();
        };
        let (_, built) = x_run_with_split_build_ahead(|| None::<u8>, Some, both);
        assert!(matches!(built, Some(XBuilt::Dropped)), "part 1 gives none");
        let (_, built) = x_run_with_split_build_ahead(|| Some(1u8), |_| None::<u8>, both);
        assert!(matches!(built, Some(XBuilt::Dropped)), "part 2 gives none");
        let (_, built) = x_run_with_split_build_ahead(|| Some(1u8), Some, |start, _| start());
        assert!(matches!(built, Some(XBuilt::Dropped)), "no vectors");
        let (_, built) = x_run_with_split_build_ahead(|| Some(1u8), Some, |_, vectors| vectors());
        assert!(built.is_none(), "part 1 not started");
        let panics = |built: Option<XBuilt<u8>>, what: &str| match built {
            Some(XBuilt::Panicked(payload)) => {
                assert_eq!(payload.downcast_ref::<&str>(), Some(&what));
            }
            _ => panic!("{what}: no panic kept"),
        };
        let (_, built) = x_run_with_split_build_ahead(
            || -> Option<u8> { panic!("part 1 panicked") },
            Some,
            both,
        );
        panics(built, "part 1 panicked");
        for vectors_first in [false, true] {
            let (_, built) = x_run_with_split_build_ahead(
                || Some(1u8),
                |_| -> Option<u8> { panic!("part 2 panicked") },
                |start, vectors| {
                    if vectors_first {
                        vectors();
                    }
                    start();
                    nested_work(2);
                    vectors();
                },
            );
            panics(built, "part 2 panicked");
        }
        Ok(())
    })
    .expect("pool");
}

fn context(lambda: f64) -> QueryContext {
    QueryContext {
        q_idx: 1,
        frame: 0,
        aa_seq: vec![0, 4, 7, 0],
        aa_seq_nomask: None,
        aa_len: 2,
        orig_len: 2,
        frame_base: 3,
        is_valid: true,
        karlin_params: crate::stats::KarlinParams {
            lambda,
            ..Default::default()
        },
    }
}

fn table(hits: Vec<i32>) -> XFlatCompressedLookup {
    XFlatCompressedLookup {
        word_length: 5,
        compressed_alphabet_size: 15,
        reciprocal_alphabet_size: 7,
        longest_chain: 2,
        pv: vec![1, 0],
        pv_array_bts: 5,
        scaled_compress_table_u8: [-1; 256],
        starts: vec![0, 2, 2],
        hits,
    }
}

// The comparisons see a difference in one slot of the table's hits, in one bit of a Karlin
// block, and in one letter of a frame; equal values pass.
#[test]
fn the_comparisons_see_one_slot_one_bit_and_one_letter() {
    assert!(x_same_flat_tables(&table(vec![3, 9]), &table(vec![3, 9])));
    assert!(!x_same_flat_tables(&table(vec![3, 9]), &table(vec![3, 8])));
    let lambda = 0.3176;
    assert!(x_same_contexts(&[context(lambda)], &[context(lambda)]));
    let next = f64::from_bits(lambda.to_bits() + 1);
    assert!(!x_same_contexts(&[context(lambda)], &[context(next)]));
    assert!(!x_same_contexts(&[context(lambda)], &[]));
    let frame = crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg(
        b"MKVLAAGIVG",
        None,
    );
    let mut other = frame.clone();
    assert!(x_frames_difference(&[vec![frame.clone()]], &[vec![other.clone()]]).is_none());
    other.aa_seq[3] ^= 1;
    assert!(x_frames_difference(&[vec![frame.clone()]], &[vec![other]]).is_some());
    assert!(x_frames_difference(&[vec![frame]], &[]).is_some());
}

/// The arguments of `losat blastp -task <task>` (SEG on) with the search matrix `matrix`
/// (`-matrix`; `None`: the default, BLOSUM62) as the search resolves them.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn test_args_of(task: &str, matrix: Option<&str>) -> ResolvedBlastpArgs {
    BlastpArgs {
        query: PathBuf::from("query.faa"),
        subject: Some(PathBuf::from("subject.faa")),
        db: None,
        db_filter: Default::default(),
        db_mask: Default::default(),
        dbsize: None,
        show_gis: false,
        html: false,
        parse_deflines: false,
        id_scope: None,
        searchsp: None,
        task: task.to_string(),
        evalue: None,
        threshold: None,
        word_size: None,
        num_threads: 1,
        mt_mode: None,
        out: None,
        max_target_seqs: None,
        num_descriptions: None,
        num_alignments: None,
        max_hsps_per_subject: None,
        qcov_hsp_perc: None,
        culling_limit: None,
        best_hit_overhang: None,
        best_hit_score_edge: None,
        subject_besthit: false,
        ungapped: false,
        xdrop_ungap: None,
        xdrop_gap: None,
        xdrop_gap_final: None,
        window_size: None,
        matrix: matrix.map(str::to_string),
        gap_open: None,
        gap_extend: None,
        comp_based_stats: "2".to_string(),
        seg: Some("yes".to_string()),
        soft_masking: None,
        lcase_masking: false,
        use_sw_tback: false,
        query_loc: None,
        subject_loc: None,
        outfmt: "6".to_string(),
        archive: None,
        export_search_strategy: None,
        strategy: None,
        request_inputs: None,
        app: BlastpApp::Blastp,
        taxdb_setting: None,
        matrix_search: None,
        write_verification: None,
        read_verification: Default::default(),
        inclusion_ethresh: None,
        gap_trigger: None,
        sum_stats: false,
        first_query_ordinal: 0,
        scope_refusal: None,
        psiblast_pssm: None,
        formatting: Default::default(),
        rps: None,
        delta: None,
    }
    .resolve()
    .expect("resolved arguments")
}

/// `count` protein records of 40 to 239 residues from a xorshift generator; every third has a
/// run of one residue (a SEG mask).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn protein_records(seed: u64, count: usize) -> Vec<FastaRecord> {
    let letters = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state = seed;
    let mut below = |n: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as usize
    };
    (0..count)
        .map(|i| {
            let length = 40 + below(200);
            let mut sequence: Vec<u8> = (0..length).map(|_| letters[below(20)]).collect();
            if i % 3 == 0 {
                let at = below(length - 30);
                let run = letters[below(20)];
                sequence[at..at + 24].fill(run);
            }
            FastaRecord::new(format!("Query_{}", i + 1), b"q", &sequence)
        })
        .collect()
}

/// The table and contexts of a block as its own search builds them (`run_resolved_in_pool`):
/// `blastp_aa_lookup` for exact blastp, `blastp_compressed_lookup` with the block's batches and
/// lengths for blastp-fast.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn own_build(
    args: &ResolvedBlastpArgs,
    records: &[FastaRecord],
    frames: &[Vec<QueryFrame>],
    batch_size: u32,
    pool: &SearchPool<'_>,
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    match args.lookup_table_type {
        BlastpLookupTableType::AaLookupTable => {
            super::super::blastp_aa_lookup(args, frames, false, pool)
        }
        BlastpLookupTableType::CompressedAaLookupTable => {
            let lengths: Vec<usize> = records.iter().map(|record| record.seq().len()).collect();
            let batches = seq_range::QueryInput::whole(records).batches(batch_size);
            super::super::blastp_compressed_lookup(
                args,
                frames,
                &batches,
                &lengths,
                pool,
                XLookupFor::Search,
                0,
            )
            .expect("the block's own build")
        }
    }
}

// A block's values built ahead (as tasks of the pool, on the second spare's vectors made ready
// as the block loop makes them) are the values the block's own search builds: the query frames,
// the table and the contexts, with the Karlin blocks bit for bit, on a pool of three threads, SEG
// masks included. Exact blastp: `x_build_ahead`, once the block before it has its table.
// blastp-fast (item SC2d): part 1 (`x_flat_ahead_appends`) and part 2 (`x_flat_ahead_table`) on
// the vectors of the first block's table, given to the second spare as after its preliminary
// stage (`x_make_room_flat`; the table holds its `starts` buffer): part 1 before the vectors,
// the vectors before part 1, and both through `x_run_with_split_build_ahead`; for BLOSUM62 and
// (item SC2d-2) for other search matrices, whose contexts have their Karlin blocks before the
// build.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn values_built_ahead_are_the_values_the_block_builds_itself() {
    let batch_size = 3000;
    for (task, matrix, order) in [
        ("blastp-fast", None, 0),
        ("blastp-fast", None, 1),
        ("blastp-fast", None, 2),
        ("blastp-fast", Some("BLOSUM45"), 0),
        ("blastp-fast", Some("BLOSUM45"), 1),
        ("blastp-fast", Some("PAM250"), 2),
        ("blastp", None, 0),
    ] {
        let args = test_args_of(task, matrix);
        let task = format!("{task} {matrix:?}");
        let first = protein_records(0x9e37_79b9_7f4a_7c15, 24);
        let next = protein_records(0xd1b5_4a32_d192_ed03, 30);
        crate::utils::threading::with_search_pool(3, "blastp", |pool| {
            let subjects = XBlockSubjects::new(&[], Vec::new());
            let pool = pool.x_with_ext(&subjects);
            let no_spares = "no spares in the pool (LOSAT_X_BPCLUT or LOSAT_X_BPLKPFLAT is \
                                 set to 0)";
            // The block loop: the first block's own table is ready.
            let first_frames = blastp_query_frames(&args, &first);
            let (first_lookup, _) = own_build(&args, &first, &first_frames, batch_size, &pool);
            let mut first_starts = None;
            let built = match first_lookup {
                BlastpRuntimeLookup::CompressedFlat(first_table) => {
                    first_starts = Some(first_table.starts.as_ptr() as usize);
                    let largest = Cell::new(0);
                    let values = |built: Option<XAheadLookup>| {
                        Some(built.map_or(XBuilt::Dropped, XBuilt::Values))
                    };
                    match order {
                        // Part 1, then the first table's vectors, then part 2.
                        0 => {
                            let part = x_flat_ahead_appends(&args, &next, batch_size, &pool);
                            assert!(
                                x_make_room_flat(&pool, first_table, &largest),
                                "{no_spares}"
                            );
                            values(part.and_then(|part| x_flat_ahead_table(part, &pool)))
                        }
                        // The vectors, then part 1 and part 2.
                        1 => {
                            assert!(
                                x_make_room_flat(&pool, first_table, &largest),
                                "{no_spares}"
                            );
                            let part = x_flat_ahead_appends(&args, &next, batch_size, &pool);
                            values(part.and_then(|part| x_flat_ahead_table(part, &pool)))
                        }
                        _ => {
                            x_run_with_split_build_ahead(
                                || x_flat_ahead_appends(&args, &next, batch_size, &pool),
                                |part| x_flat_ahead_table(part, &pool),
                                |start, vectors| {
                                    start();
                                    assert!(
                                        x_make_room_flat(&pool, first_table, &largest),
                                        "{no_spares}"
                                    );
                                    vectors();
                                },
                            )
                            .1
                        }
                    }
                }
                // Exact blastp: the vectors for the next block's table are made ready, and its
                // build is spawned.
                first_lookup => {
                    assert!(
                        x_make_room(&pool, &first_lookup, &Cell::new(0)),
                        "{no_spares}"
                    );
                    x_run_with_build_ahead(|| x_build_ahead(&args, &next, &pool), |spawn| spawn()).1
                }
            };
            let Some(XBuilt::Values(ahead)) = built else {
                panic!(
                    "{task}: no values built ahead (LOSAT_X_BPLKPAHEAD is 0, LOSAT_X_BPQBLOCK \
                         is changed or LOSAT_DIAGNOSTICS is set)"
                );
            };
            // The next block's own search.
            let frames = blastp_query_frames(&args, &next);
            assert!(
                frames
                    .iter()
                    .flatten()
                    .any(|frame| !frame.seg_masks.is_empty()),
                "{task}: no SEG mask in the fixture"
            );
            if let Some(difference) = x_frames_difference(&ahead.query_frames, &frames) {
                panic!("{task}: {difference}");
            }
            let (table, contexts) = own_build(&args, &next, &frames, batch_size, &pool);
            assert!(!contexts.is_empty(), "{task}: no contexts");
            match (&ahead.table, &table) {
                (
                    BlastpRuntimeLookup::CompressedFlat(got),
                    BlastpRuntimeLookup::CompressedFlat(want),
                ) => {
                    assert!(!want.hits.is_empty(), "{task}: the flat table has no hits");
                    assert_eq!(
                        Some(got.starts.as_ptr() as usize),
                        first_starts,
                        "{task}, order {order}: the table is not on the first table's vectors"
                    );
                    if !x_same_flat_tables(got, want) {
                        x_assert_same_flat_tables_as("ahead:", got, want);
                        panic!("{task}: the tables differ");
                    }
                }
                (BlastpRuntimeLookup::Aa(got), BlastpRuntimeLookup::Aa(want)) => {
                    assert!(
                        want.backbone.iter().any(|cell| cell.num_used > 0)
                            || !want.overflow.is_empty(),
                        "{task}: the table has no hits"
                    );
                    if !x_same_table(got, want) {
                        x_assert_same_table_as("ahead:", got, want);
                        panic!("{task}: the tables differ");
                    }
                }
                // LOSAT_X_BPLK16: the tables of 16-bit hits, field by field as
                // `x_same_table`.
                (BlastpRuntimeLookup::Aa16(got), BlastpRuntimeLookup::Aa16(want)) => {
                    assert!(
                        want.backbone.iter().any(|cell| cell.num_used > 0)
                            || !want.overflow.is_empty(),
                        "{task}: the table has no hits"
                    );
                    assert!(
                        got.backbone.len() == want.backbone.len()
                            && got.backbone.iter().zip(&want.backbone).all(|(a, b)| {
                                a.num_used == b.num_used && a.entries == b.entries
                            })
                            && got.overflow == want.overflow
                            && got.pv == want.pv
                            && got.frame_bases == want.frame_bases
                            && got.num_contexts == want.num_contexts
                            && got.query_length == want.query_length
                            && got.word_length == want.word_length
                            && got.alphabet_size == want.alphabet_size
                            && got.charsize == want.charsize
                            && got.mask == want.mask
                            && got.longest_chain == want.longest_chain
                            && got.threshold == want.threshold
                            && got.row_max == want.row_max,
                        "{task}: the tables of 16-bit hits differ"
                    );
                }
                _ => panic!("{task}: the tables are of different kinds"),
            }
            if !x_same_contexts(&ahead.contexts, &contexts) {
                x_assert_same_contexts_as("ahead:", &ahead.contexts, &contexts);
                panic!("{task}: the contexts differ");
            }
            Ok(())
        })
        .expect("pool");
    }
}

// A block's search takes values built ahead only for its own records (the same slice); others are
// dropped and the search builds its own.
#[test]
fn values_built_ahead_are_taken_only_by_their_block() {
    let records = FastaRecord::new("Query_1", b"q1", b"MKVLAAGIVG");
    let records = vec![records.clone(), records];
    let values = |key| XAheadLookup {
        key,
        query_frames: Vec::new(),
        table: BlastpRuntimeLookup::CompressedFlat(table(Vec::new())),
        contexts: Vec::new(),
    };
    let mut block = XAheadBlock {
        prebuilt: Some(values(x_key(&records))),
        start_next: None,
        after_prelim_flat: None,
        table_taken_later: true,
    };
    assert!(block.take_prebuilt(&records[..1]).is_none());
    let mut block = XAheadBlock {
        prebuilt: Some(values(x_key(&records))),
        start_next: None,
        after_prelim_flat: None,
        table_taken_later: true,
    };
    assert!(block.take_prebuilt(&records).is_some());
    assert!(block.take_prebuilt(&records).is_none());
    assert!(XAheadBlock::none().take_prebuilt(&records).is_none());
}
