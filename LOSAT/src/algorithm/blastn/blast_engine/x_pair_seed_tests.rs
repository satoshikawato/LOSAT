#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::super::{
    build_blastna_matrix, build_nucl_score_table, build_query_blastna_buffer,
    diag_hash_insert_window, scan_subject_kmers_with_ranges, select_mb_scan_kind,
};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::*;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::algorithm::blastn::lookup::{build_two_stage_lookup, reverse_complement};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::core::blast_encoding::{encode_iupac_to_blastna, encode_subject_ncbi2na_packed};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::utils::dust::MaskedInterval;

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Lcg(u64);
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
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

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
type HashState = (Vec<u32>, Vec<(i32, i32, bool, i32, u32)>, u32, u32, i32);

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn hash_state(h: &DiagHashTable) -> HashState {
    (
        h.backbone.clone(),
        h.chain[..h.occupancy as usize]
            .iter()
            .map(|c| (c.diag, c.level, c.hit_saved, c.hit_len, c.next))
            .collect(),
        h.occupancy,
        h.capacity,
        h.offset,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-425
// A random stream of seed outcomes on the diagonal hash: retrieve, the `last_hit` test,
// then an insert as `fold` does (below the cutoff: `hit_len` = s_end_pos - s_off_pos,
// not saved; saved: the ungapped end, `hit_len` 0). Diagonals d + 512 k share a bucket,
// so cells of other diagonals are taken over (the stale rule); some cells start at
// `level = -window` (a table reset with a window, the round-1 sign-extension case) and
// the subject offsets start below zero around them. The stream folded from strips
// scanned on 1..8 threads must give the decisions and the table of the serial fold.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[derive(Clone, Copy, Debug)]
struct Op {
    diag: i32,
    s_off_pos: i32,
    s_end_pos: i32,
    pass: bool,
    ungapped_end: i32,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn apply(hash: &mut DiagHashTable, op: &Op, window: i32, decisions: &mut Vec<u8>) {
    let last_hit = hash.retrieve(op.diag).map_or(0, |(level, _, _)| level);
    if op.s_off_pos < last_hit {
        decisions.push(0);
        return;
    }
    if op.pass {
        decisions.push(2);
        hash.insert_plain(op.diag, op.ungapped_end, 0, true, op.s_off_pos, window);
    } else {
        decisions.push(1);
        hash.insert_plain(
            op.diag,
            op.s_end_pos,
            op.s_end_pos - op.s_off_pos,
            false,
            op.s_off_pos,
            window,
        );
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn diag_hash_fold_from_strips_equals_serial_fold() {
    let mut rng = Lcg(0xD1A6_0001);
    for case in 0..40 {
        let word_length = if case % 2 == 0 { 28 } else { 11 };
        let window = diag_hash_insert_window(0, 0, word_length);
        let table_window = 16;
        let mut start = DiagHashTable::new(table_window);
        let base = rng.below(1000) as i32 - 500;
        let diags: Vec<i32> = (0..12)
            .map(|k| {
                if k < 8 {
                    base + 512 * k
                } else {
                    rng.below(4000) as i32 - 2000
                }
            })
            .collect();
        // cells left at `last_hit = -window` by a reset with a window
        for &d in diags.iter().take(3) {
            start.insert_plain(d, -table_window, 0, false, -table_window, window);
        }
        let n = 200 + rng.below(3000) as usize;
        let mut s = -20i32;
        let ops: Vec<Op> = (0..n)
            .map(|_| {
                s += rng.below(9) as i32 - 2;
                let s_end = s + word_length as i32 + rng.below(4) as i32;
                Op {
                    diag: diags[rng.below(diags.len() as u64) as usize],
                    s_off_pos: s,
                    s_end_pos: s_end,
                    pass: rng.below(3) == 0,
                    ungapped_end: s_end + rng.below(200) as i32,
                }
            })
            .collect();
        let mut serial = copy_hash(&start);
        let mut serial_decisions = Vec::new();
        for op in &ops {
            apply(&mut serial, op, window, &mut serial_decisions);
        }
        for threads in 1..=8usize {
            let mut cuts = vec![0usize];
            while *cuts.last().unwrap() < n {
                let next = cuts.last().unwrap() + 1 + rng.below(64) as usize;
                cuts.push(next.min(n));
            }
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let (state, decisions) = pool.install(|| {
                let mut hash = copy_hash(&start);
                let mut decisions = Vec::new();
                let produce = |i: usize, out: &mut Vec<Op>| {
                    out.extend_from_slice(&ops[cuts[i]..cuts[i + 1]]);
                };
                let (mut a, mut b) = (0u64, 0u64);
                run_strips(
                    cuts.len() - 1,
                    threads - 1,
                    (case % 3) * threads,
                    &produce,
                    |_, strip| {
                        for op in strip {
                            apply(&mut hash, op, window, &mut decisions);
                        }
                    },
                    &mut a,
                    &mut b,
                    false,
                );
                (hash_state(&hash), decisions)
            });
            assert_eq!(decisions, serial_decisions, "case {case} threads {threads}");
            assert_eq!(state, hash_state(&serial), "case {case} threads {threads}");
        }
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn random_dna(rng: &mut Lcg, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
// The whole stage on random sequences (the subject made of mutated copies of query
// pieces on both strands, query masks on one context) for the scan kernels of every
// lookup width and step: the pairs of the reference wrapper
// `scan_subject_kmers_with_ranges` folded in one go, against `seed_stage` on 1..8
// threads (strips planned for each thread count). Hit lists and diagonal hashes must be
// the same.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn seed_stage_on_strips_equals_one_scan() {
    let mut rng = Lcg(0x5EED_5EED);
    let configs: [(usize, usize); 14] = [
        (28, 12),
        (19, 12),
        (11, 11),
        (16, 11),
        (18, 11),
        (21, 11),
        (9, 9),
        (10, 9),
        (10, 10),
        (11, 10),
        (12, 10),
        (16, 8),
        (15, 8),
        (24, 12),
    ];
    for (case, &(word_length, lut_word_length)) in configs.iter().enumerate() {
        let scan_step = word_length - lut_word_length + 1;
        let q_len = 3000 + rng.below(3000) as usize;
        let mut query = random_dna(&mut rng, q_len);
        // a repeat inside the query (long lookup chains)
        let motif = random_dna(&mut rng, 40);
        for k in 0..20 {
            let at = (k * 97) % (q_len - 40);
            query[at..at + 40].copy_from_slice(&motif);
        }
        let minus = reverse_complement(&query);
        let mut subject = Vec::new();
        while subject.len() < 40_000 {
            if rng.below(3) == 0 {
                let len = 50 + rng.below(400) as usize;
                subject.extend(random_dna(&mut rng, len));
            } else {
                let src = if rng.below(2) == 0 { &query } else { &minus };
                let len = 30 + rng.below(1500) as usize;
                let at = rng.below((q_len - len.min(q_len - 1)) as u64) as usize;
                let end = (at + len).min(q_len);
                for &b in &src[at..end] {
                    let r = rng.below(100);
                    if r < 3 {
                        subject.push(b"ACGT"[rng.below(4) as usize]);
                    } else if r < 4 {
                        // an indel
                    } else {
                        subject.push(b);
                    }
                }
            }
        }
        let s_len = subject.len();
        let packed = encode_subject_ncbi2na_packed(&subject);
        let ctx_seqs = [
            encode_iupac_to_blastna(&query),
            encode_iupac_to_blastna(&minus),
        ];
        let offsets = vec![0i32, q_len as i32 + 1];
        let masks = vec![
            vec![
                MaskedInterval {
                    start: 100,
                    end: 180,
                },
                MaskedInterval {
                    start: 1000,
                    end: 1013,
                },
            ],
            Vec::new(),
        ];
        let concat_len = 2 * q_len + 1;
        let concat_sentinels =
            build_query_blastna_buffer(&[query.as_slice()][..], &offsets, concat_len);
        let concat = concat_sentinels[1..concat_len + 1].to_vec();
        let ctx_refs: Vec<&[u8]> = ctx_seqs.iter().map(Vec::as_slice).collect();
        let lookup = build_two_stage_lookup(
            &ctx_refs,
            &offsets,
            word_length,
            lut_word_length,
            &masks,
            None,
            0,
            2 * q_len,
            false,
        );
        let contexts: Vec<QueryContext> = (0..2)
            .map(|c| QueryContext {
                query_idx: 0,
                frame: if c == 0 { 1 } else { -1 },
                query_offset: offsets[c],
                len: ctx_seqs[c].len(),
                masks: masks[c].clone(),
            })
            .collect();
        let index = QueryContextIndex::new(&contexts, false);
        let score_matrix = build_blastna_matrix(1, -2);
        let nucl_score_table = build_nucl_score_table(1, -2);
        // every other case a masked subject: scan ranges (NCBI's `seq_ranges`, right end
        // exclusive) over most of it with gaps of 1 to 120 bases, some shorter than a word
        // (the scan skips them), the last perhaps at the subject's end
        let masked = case % 2 == 1;
        let seq_ranges = if masked {
            let mut ranges = Vec::new();
            let mut at = rng.below(30) as usize;
            while at < s_len {
                let len = if rng.below(10) == 0 {
                    1 + rng.below(20) as usize
                } else {
                    30 + rng.below(2000) as usize
                };
                let end = (at + len).min(s_len);
                ranges.push((at as i32, end as i32));
                at = end + 1 + rng.below(120) as usize;
            }
            assert!(ranges.len() > 10);
            ranges
        } else {
            vec![(0i32, s_len as i32)]
        };
        let diag_offset = 37 + rng.below(1000) as isize;
        // the batch's forms (`ScanBatch::forms`), chosen as `ScanBatch::new` chooses them
        let forms = ScanForms::choose(Some(&lookup), lut_word_length, scan_step);
        let x = Inputs {
            two_stage: &lookup,
            forms: &forms,
            search_seq_packed: &packed,
            s_len,
            subject_seq_ranges: &seq_ranges,
            subject_masked: masked,
            scan_step,
            query_context_index: &index,
            query_contexts: &contexts,
            encoded_query_concat_blastna: &concat,
            encoded_query_concat_blastna_with_sentinels: &concat_sentinels,
            cutoff_scores: &[22, 22],
            x_dropoff_scores: &[16, 16],
            reduced_cutoff_scores: &[18, 18],
            score_matrix: &score_matrix,
            nucl_score_table: &nucl_score_table,
            exact_ungapped: word_length < 11,
            diag_offset,
            diag_hash_window: diag_hash_insert_window(0, 0, word_length),
            use_array_indexing: false,
            window_size: 0,
            scan_range: 0,
            small_na_word: false,
            diagnostics: false,
            one_subject_owner: true,
            recording: crate::utils::xreach::recording(),
            timing: None,
        };
        let new_hash = || {
            let mut h = DiagHashTable::new(0);
            h.offset = diag_offset as i32;
            h
        };
        // the reference wrapper visits the pairs in one scan
        let mut seeds = Vec::new();
        let kind = select_mb_scan_kind(lut_word_length, scan_step, masked);
        scan_subject_kmers_with_ranges(
            &packed,
            s_len,
            word_length,
            lut_word_length,
            scan_step,
            &seq_ranges,
            masked,
            kind,
            |kmer_start, s_range, kmer| {
                lookup.for_each_hit(kmer, |q_off_1| {
                    seed_of_pair(&x, q_off_1 as usize - 1, kmer_start, s_range, &mut seeds);
                });
            },
        );
        let mut serial_hash = new_hash();
        let mut serial_hits = Vec::new();
        let mut counts = FoldCounts::default();
        let mut pause = 0u64;
        fold::<false>(
            &x,
            &seeds,
            &mut serial_hash,
            &mut serial_hits,
            None,
            &mut pause,
            &mut counts,
            &mut part::PartOut::default(),
            &[],
        );
        assert!(
            serial_hits.len() > 5 && counts.skipped > 0,
            "case {case}: the test sequences give too few hits ({} hits, {} skipped)",
            serial_hits.len(),
            counts.skipped
        );
        let key = |h: &UngappedHit| {
            (
                h.context_idx,
                h.seed_q_off,
                h.seed_s_off,
                h.qs,
                h.qe,
                h.ss,
                h.se,
                h.score,
            )
        };
        let serial_keys: Vec<_> = serial_hits.iter().map(key).collect();
        for threads in 1..=8usize {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let (state, keys, stats_seeds) = pool.install(|| {
                let mut hash = new_hash();
                let mut hits = Vec::new();
                let stats = seed_stage(&x, threads, &mut hash, &mut hits, None, false, None);
                (
                    hash_state(&hash),
                    hits.iter().map(key).collect::<Vec<_>>(),
                    stats.counts.seeds,
                )
            });
            assert_eq!(
                stats_seeds,
                seeds.len() as u64,
                "case {case} threads {threads}"
            );
            assert_eq!(keys, serial_keys, "case {case} threads {threads}");
            assert_eq!(
                state,
                hash_state(&serial_hash),
                "case {case} threads {threads}"
            );
        }
        // LOSAT_X_DIAGLV (`x_diag_lv.rs`): on a table with level bounds the serial fold and the
        // strips on 1, 2, 3 and 7 threads skip walks and give the hits and the table of
        // the fold without bounds, and the same bounds.
        let new_lv_hash = || {
            let mut h = new_hash();
            h.x_lv_start();
            h
        };
        let mut lv_hash = new_lv_hash();
        let mut lv_hits = Vec::new();
        let before = super::super::x_diag_lv::test_skips();
        fold::<false>(
            &x,
            &seeds,
            &mut lv_hash,
            &mut lv_hits,
            None,
            &mut pause,
            &mut FoldCounts::default(),
            &mut part::PartOut::default(),
            &[],
        );
        let skipped = super::super::x_diag_lv::test_skips() - before;
        assert!(skipped > 0, "case {case}: no walk skipped");
        assert_eq!(lv_hits.iter().map(key).collect::<Vec<_>>(), serial_keys);
        assert_eq!(
            hash_state(&lv_hash),
            hash_state(&serial_hash),
            "case {case}"
        );
        for threads in [1usize, 2, 3, 7] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let (state, keys, max_level, skipped) = pool.install(|| {
                let mut hash = new_lv_hash();
                let mut hits = Vec::new();
                let before = super::super::x_diag_lv::test_skips();
                seed_stage(&x, threads, &mut hash, &mut hits, None, false, None);
                (
                    hash_state(&hash),
                    hits.iter().map(key).collect::<Vec<_>>(),
                    hash.max_level.clone(),
                    super::super::x_diag_lv::test_skips() - before,
                )
            });
            let at = format!("case {case} threads {threads}");
            assert!(skipped > 0, "{at}: no walk skipped");
            assert_eq!(keys, serial_keys, "{at}");
            assert_eq!(state, hash_state(&serial_hash), "{at}");
            assert_eq!(max_level, lv_hash.max_level, "{at}");
        }
        // LOSAT_X_PAIRPARFOLD (`x_pair_seed_part.rs`): the fold split into 2, 3, 7, 8 and 16
        // partitions on 1..8 threads, on tables with and without level bounds, for two
        // chunks in a row on one scratch (the second with the hash offset moved on), in both
        // forms of LOSAT_X_BNSEEDLISTS's stage (the scratch passing its strip buffers from one
        // form to the other): the hits and the table of the serial fold.
        let x2 = Inputs {
            diag_offset: diag_offset + s_len as isize + 16,
            ..x
        };
        for levels in [false, true] {
            let new_table = || {
                let mut h = new_hash();
                if levels {
                    h.x_lv_start();
                }
                h
            };
            let mut expected_hash = new_table();
            let mut expected = Vec::new();
            for (chunk, xc) in [&x, &x2].into_iter().enumerate() {
                expected_hash.offset = xc.diag_offset as i32;
                let mut hits = Vec::new();
                fold::<false>(
                    xc,
                    &seeds,
                    &mut expected_hash,
                    &mut hits,
                    None,
                    &mut pause,
                    &mut FoldCounts::default(),
                    &mut part::PartOut::default(),
                    &[],
                );
                expected.push((
                    hits.iter().map(key).collect::<Vec<_>>(),
                    hash_state(&expected_hash),
                    expected_hash.max_level.clone(),
                ));
                assert!(chunk == 0 || !hits.is_empty());
            }
            for (threads, parts) in [
                (1usize, 2usize),
                (2, 2),
                (3, 3),
                (3, 7),
                (4, 4),
                (5, 16),
                (7, 7),
                (8, 8),
                (8, 16),
            ] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                pool.install(|| {
                    let mut scratch = part::PartScratch::default();
                    // EXPERIMENT (LOSAT_X_SEEDPIECE): the lists' stage also in pieces of 1
                    // seed, a few and hundreds (the kernels stop inside words' runs and
                    // restart at many offsets), on the same scratch
                    for (lists, pieces) in [
                        (false, 0usize),
                        (true, 0),
                        (true, 1),
                        (true, 2 + (case + threads) % 9),
                        (true, 200 + 31 * parts),
                    ] {
                        scratch.x_pieces = Some(pieces);
                        let mut hash = new_table();
                        for (chunk, xc) in [&x, &x2].into_iter().enumerate() {
                            let at = format!(
                                "case {case} threads {threads} parts {parts} levels {levels} lists {lists} pieces {pieces} chunk {chunk}"
                            );
                            hash.offset = xc.diag_offset as i32;
                            let mut hits = Vec::new();
                            let stats = seed_stage_of(
                                xc,
                                threads,
                                &mut hash,
                                &mut hits,
                                None,
                                false,
                                Some(part::FoldPart {
                                    scratch: &mut scratch,
                                    parts,
                                    latest: None,
                                }),
                                lists,
                            );
                            assert_eq!(stats.counts.seeds, seeds.len() as u64, "{at}");
                            let (keys, state, max_level) = &expected[chunk];
                            assert_eq!(
                                &hits.iter().map(key).collect::<Vec<_>>(),
                                keys,
                                "{at}"
                            );
                            assert_eq!(&hash_state(&hash), state, "{at}");
                            assert_eq!(&hash.max_level, max_level, "{at}");
                        }
                    }
                });
            }
        }
    }
}
