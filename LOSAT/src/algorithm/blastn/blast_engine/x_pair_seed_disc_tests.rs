#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use super::super::super::scan_subject_disc_words_with_ranges;
use super::super::super::two_stage_scan::ScanForms;
use super::super::super::{
    build_blastna_matrix, build_nucl_score_table, build_query_blastna_buffer,
    diag_hash_insert_window, select_mb_scan_kind, QueryContext, QueryContextIndex,
};
use super::*;
use crate::algorithm::blastn::disc_lookup::DiscWordType;
use crate::algorithm::blastn::lookup::{
    build_disc_mb_lookup, build_two_stage_lookup, reverse_complement, TwoStageLookup,
};
use crate::core::blast_encoding::{encode_iupac_to_blastna, encode_subject_ncbi2na_packed};
use crate::utils::dust::MaskedInterval;

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

fn random_dna(rng: &mut Lcg, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
type HashState = (Vec<u32>, Vec<(i32, i32, bool, i32, u32)>, u32, u32, i32);

/// `seed_stage_of` in the instance of `lean` (LOSAT_X_FOLDLEAN) and of the table's memo
/// (LOSAT_X_DIAGLV2), with LOSAT_X_DCSEEDLEAN's stage or not (`lists`), whatever the
/// environment says (`seed_stage` reads the switches), after the memo is made again at the
/// stage's start as `seed_stage` does.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn stage_of(
    lean: bool,
    lists: bool,
    x: &Inputs<'_>,
    disc: &DiscTemplates,
    threads: usize,
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    fold_part: Option<part::FoldPart<'_>>,
) -> StageStats {
    if !hash.x_last.is_empty() {
        hash.x_lv2_reset();
    }
    match (lean, !hash.x_last.is_empty()) {
        (false, false) => seed_stage_of::<false, false, false>(
            x, disc, threads, hash, hits, None, false, fold_part, lists, None,
        ),
        (false, true) => seed_stage_of::<false, true, false>(
            x, disc, threads, hash, hits, None, false, fold_part, lists, None,
        ),
        (true, false) => seed_stage_of::<true, false, false>(
            x, disc, threads, hash, hits, None, false, fold_part, lists, None,
        ),
        (true, true) => seed_stage_of::<true, true, false>(
            x, disc, threads, hash, hits, None, false, fold_part, lists, None,
        ),
    }
}

/// `fold_two_hit` without partitions in the instance of `lean` (LOSAT_X_FOLDLEAN) and of
/// the table's memo (LOSAT_X_DIAGLV2), as `seed_stage` picks it.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[allow(clippy::too_many_arguments)]
fn fold_serial(
    lean: bool,
    x: &Inputs<'_>,
    seeds: &[Seed],
    hash: &mut DiagHashTable,
    hits: &mut Vec<UngappedHit>,
    pause: &mut u64,
    counts: &mut FoldCounts,
) {
    let out = &mut part::PartOut::default();
    match (lean, !hash.x_last.is_empty()) {
        (false, false) => {
            fold_two_hit::<false, false, false>(x, seeds, hash, hits, None, pause, counts, out, &[])
        }
        (false, true) => {
            fold_two_hit::<false, false, true>(x, seeds, hash, hits, None, pause, counts, out, &[])
        }
        (true, false) => {
            fold_two_hit::<false, true, false>(x, seeds, hash, hits, None, pause, counts, out, &[])
        }
        (true, true) => {
            fold_two_hit::<false, true, true>(x, seeds, hash, hits, None, pause, counts, out, &[])
        }
    }
}

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

/// A query of two contexts (plus and minus strand, masks on the plus strand) and a
/// subject made of mutated copies of query pieces on both strands and random pieces.
struct Case {
    contexts: Vec<QueryContext>,
    index: QueryContextIndex,
    ctx_seqs: Vec<Vec<u8>>,
    offsets: Vec<i32>,
    masks: Vec<Vec<MaskedInterval>>,
    concat: Vec<u8>,
    concat_sentinels: Vec<u8>,
    packed: Vec<u8>,
    s_len: usize,
    seq_ranges: Vec<(i32, i32)>,
    masked: bool,
    score_matrix: [i32; 256],
    nucl_score_table: [i32; 256],
}

fn make_case(rng: &mut Lcg, subject_len: usize) -> Case {
    let q_len = 4100 + rng.below(3000) as usize;
    make_case_with(rng, subject_len, q_len)
}

/// `make_case` with a query of `q_len` bases (more than 1013: the masks).
fn make_case_with(rng: &mut Lcg, subject_len: usize, q_len: usize) -> Case {
    let mut query = random_dna(rng, q_len);
    // a repeat inside the query (long lookup chains)
    let motif = random_dna(rng, 40);
    for k in 0..20 {
        let at = (k * 97) % (q_len - 40);
        query[at..at + 40].copy_from_slice(&motif);
    }
    let minus = reverse_complement(&query);
    let mut subject = Vec::new();
    while subject.len() < subject_len {
        if rng.below(3) == 0 {
            let len = 50 + rng.below(400) as usize;
            subject.extend(random_dna(rng, len));
        } else {
            let src = if rng.below(2) == 0 { &query } else { &minus };
            let len = 30 + rng.below(1500) as usize;
            let at = rng.below((q_len - len.min(q_len - 1)) as u64) as usize;
            let end = (at + len).min(q_len);
            for &b in &src[at..end] {
                // dc-megablast finds diverged copies: 10 % substitutions
                let r = rng.below(100);
                if r < 10 {
                    subject.push(b"ACGT"[rng.below(4) as usize]);
                } else if r < 11 {
                    // an indel
                } else {
                    subject.push(b);
                }
            }
        }
    }
    let s_len = subject.len();
    let packed = encode_subject_ncbi2na_packed(&subject);
    let ctx_seqs = vec![
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
    Case {
        contexts,
        index,
        ctx_seqs,
        offsets,
        masks,
        concat,
        concat_sentinels,
        packed,
        s_len,
        seq_ranges: vec![(0, s_len as i32)],
        masked: false,
        score_matrix: build_blastna_matrix(2, -3),
        nucl_score_table: build_nucl_score_table(2, -3),
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
/// `c` as a masked subject: its scan ranges (NCBI's `seq_ranges`, right end exclusive) cover
/// most of it with gaps of 1 to 120 bases between them; some are shorter than any template
/// (the scan skips them) and the last may end at the subject's end.
fn mask_case(rng: &mut Lcg, mut c: Case) -> Case {
    let mut ranges = Vec::new();
    let mut at = rng.below(30) as usize;
    while at < c.s_len {
        let len = if rng.below(10) == 0 {
            1 + rng.below(20) as usize
        } else {
            30 + rng.below(2000) as usize
        };
        let end = (at + len).min(c.s_len);
        ranges.push((at as i32, end as i32));
        at = end + 1 + rng.below(120) as usize;
    }
    assert!(ranges.len() > 10);
    c.seq_ranges = ranges;
    c.masked = true;
    c
}

fn ctx_refs(c: &Case) -> Vec<&[u8]> {
    c.ctx_seqs.iter().map(Vec::as_slice).collect()
}

/// The batch's forms (`ScanBatch::forms`) of a case's table `lookup`, chosen as `ScanBatch::new`
/// chooses them: from the table, its lookup word length and the scan step of `inputs` (1). The
/// discontiguous stages do not read them (a discontiguous table has no CSR head, and its scans
/// take their scanners from its templates); the parent's `run` checks them against the table in
/// debug builds (`ScanForms::chosen_from`).
fn forms_of(lookup: &TwoStageLookup) -> ScanForms<'_> {
    ScanForms::choose(Some(lookup), lookup.lut_word_length(), 1)
}

/// The inputs of a chunk of case `c` with the table `lookup` and its forms `forms` (`forms_of`).
fn inputs<'a>(
    c: &'a Case,
    lookup: &'a TwoStageLookup,
    forms: &'a ScanForms<'a>,
    diag_offset: isize,
    window_size: usize,
) -> Inputs<'a> {
    Inputs {
        two_stage: lookup,
        forms,
        search_seq_packed: &c.packed,
        s_len: c.s_len,
        subject_seq_ranges: &c.seq_ranges,
        subject_masked: c.masked,
        scan_step: 1,
        query_context_index: &c.index,
        query_contexts: &c.contexts,
        encoded_query_concat_blastna: &c.concat,
        encoded_query_concat_blastna_with_sentinels: &c.concat_sentinels,
        cutoff_scores: &[30, 30],
        x_dropoff_scores: &[20, 20],
        reduced_cutoff_scores: &[24, 24],
        score_matrix: &c.score_matrix,
        nucl_score_table: &c.nucl_score_table,
        exact_ungapped: lookup.word_length() < 11,
        diag_offset,
        diag_hash_window: diag_hash_insert_window(window_size, 0, lookup.word_length()),
        use_array_indexing: false,
        window_size,
        scan_range: 0,
        small_na_word: false,
        diagnostics: false,
        one_subject_owner: true,
        recording: crate::utils::xreach::recording(),
        timing: None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1651-1655
// The whole two-hit stage on random sequences for every discontiguous template (word
// size 11 and 12, coding, optimal and both, lengths 16, 18 and 21; each one-template
// type picks one of the scan kernels, two templates the two-template kernel): the pairs
// of the reference wrapper `scan_subject_disc_words_with_ranges` (one scan of the whole
// range) folded in one go, against `seed_stage` on 1..8 threads (strips planned for each
// thread count, so the kernels restart at many offsets). Every second case is a masked
// subject (`mask_case`): the reference scans its ranges one by one, each with its own
// `s_range`. Hit lists and diagonal hashes must be the same, and the serial fold must have
// taken every branch: the `last_hit` rejection, the first hit of a window (insert with
// `hit_ready` 0), the second hit (ungapped extension) both below and above the cutoff.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn disc_two_hit_stage_on_strips_equals_one_scan() {
    let mut rng = Lcg(0xD15C_5EED);
    let mut case_no = 0;
    for word_size in [11usize, 12] {
        for word_type in [
            DiscWordType::Coding,
            DiscWordType::Optimal,
            DiscWordType::TwoTemplates,
        ] {
            for template_length in [16usize, 18, 21] {
                case_no += 1;
                let c = make_case(&mut rng, 30_000);
                // every second case: a masked subject, scanned in its ranges
                let c = if case_no % 2 == 0 {
                    mask_case(&mut rng, c)
                } else {
                    c
                };
                let lookup = build_disc_mb_lookup(
                    &ctx_refs(&c),
                    &c.offsets,
                    word_size,
                    template_length,
                    word_type,
                    &c.masks,
                    c.concat.len(),
                );
                let forms = forms_of(&lookup);
                let disc = lookup.disc().expect("a discontiguous table");
                assert_eq!(lookup.word_length(), template_length);
                let diag_offset = 37 + rng.below(1000) as isize;
                let window_size = 40;
                let x = inputs(&c, &lookup, &forms, diag_offset, window_size);
                assert!(in_scope(&x));
                let new_hash = || {
                    let mut h = DiagHashTable::new(window_size as i32);
                    h.offset = diag_offset as i32;
                    h
                };
                let mut seeds = Vec::new();
                scan_subject_disc_words_with_ranges(
                    &c.packed,
                    c.s_len,
                    disc,
                    &c.seq_ranges,
                    c.masked,
                    |s_off, s_range, index, index2| {
                        lookup.for_each_hit(index, |q_off_1| {
                            if q_off_1 != 0 {
                                seed_of_pair(&x, q_off_1 as usize - 1, s_off, s_range, &mut seeds);
                            }
                        });
                        if let Some(index2) = index2 {
                            lookup.for_each_hit2(index2, |q_off_1| {
                                seed_of_pair(&x, q_off_1 as usize - 1, s_off, s_range, &mut seeds);
                            });
                        }
                    },
                );
                let mut serial_hash = new_hash();
                let mut serial_hits = Vec::new();
                let mut counts = FoldCounts::default();
                let mut pause = 0u64;
                fold_serial(
                    false,
                    &x,
                    &seeds,
                    &mut serial_hash,
                    &mut serial_hits,
                    &mut pause,
                    &mut counts,
                );
                let saved = serial_hash.chain[1..serial_hash.occupancy as usize]
                    .iter()
                    .filter(|cell| cell.hit_saved)
                    .count();
                let first_hits = serial_hash.chain[1..serial_hash.occupancy as usize]
                    .iter()
                    .filter(|cell| !cell.hit_saved && cell.hit_len > 0)
                    .count();
                assert!(
                    serial_hits.len() > 5 && counts.skipped > 0 && saved > 0 && first_hits > 0,
                    "case {case_no}: the test sequences give too few hits ({} hits, {} skipped, {saved} saved cells, {first_hits} first-hit cells)",
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
                    let (state, keys, stats_seeds, strips) = pool.install(|| {
                        let mut hash = new_hash();
                        let mut hits = Vec::new();
                        let stats =
                            seed_stage(&x, disc, threads, &mut hash, &mut hits, None, false, None);
                        (
                            hash_state(&hash),
                            hits.iter().map(key).collect::<Vec<_>>(),
                            stats.counts.seeds,
                            stats.strips,
                        )
                    });
                    assert!(strips > 1, "case {case_no} threads {threads}");
                    assert_eq!(
                        stats_seeds,
                        seeds.len() as u64,
                        "case {case_no} threads {threads}"
                    );
                    assert_eq!(keys, serial_keys, "case {case_no} threads {threads}");
                    assert_eq!(
                        state,
                        hash_state(&serial_hash),
                        "case {case_no} threads {threads}"
                    );
                }
                // LOSAT_X_DIAGLV (`x_diag_lv.rs`): on a table with level bounds the serial fold
                // and the strips on 1, 2, 3 and 7 threads skip walks and give the hits and
                // the table of the fold without bounds, and the same bounds.
                let new_lv_hash = || {
                    let mut h = new_hash();
                    h.x_lv_start();
                    h
                };
                let mut lv_hash = new_lv_hash();
                let mut lv_hits = Vec::new();
                let before = super::super::super::x_diag_lv::test_skips();
                fold_serial(
                    false,
                    &x,
                    &seeds,
                    &mut lv_hash,
                    &mut lv_hits,
                    &mut pause,
                    &mut FoldCounts::default(),
                );
                let skipped = super::super::super::x_diag_lv::test_skips() - before;
                assert!(skipped > 0, "case {case_no}: no walk skipped");
                assert_eq!(lv_hits.iter().map(key).collect::<Vec<_>>(), serial_keys);
                assert_eq!(
                    hash_state(&lv_hash),
                    hash_state(&serial_hash),
                    "case {case_no}"
                );
                for threads in [1usize, 2, 3, 7] {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()
                        .unwrap();
                    let (state, keys, max_level, skipped) = pool.install(|| {
                        let mut hash = new_lv_hash();
                        let mut hits = Vec::new();
                        let before = super::super::super::x_diag_lv::test_skips();
                        seed_stage(&x, disc, threads, &mut hash, &mut hits, None, false, None);
                        (
                            hash_state(&hash),
                            hits.iter().map(key).collect::<Vec<_>>(),
                            hash.max_level.clone(),
                            super::super::super::x_diag_lv::test_skips() - before,
                        )
                    });
                    let at = format!("case {case_no} threads {threads}");
                    assert!(skipped > 0, "{at}: no walk skipped");
                    assert_eq!(keys, serial_keys, "{at}");
                    assert_eq!(state, hash_state(&serial_hash), "{at}");
                    assert_eq!(max_level, lv_hash.max_level, "{at}");
                }
                // LOSAT_X_DIAGLV2 (`x_diag_lv.rs`) and LOSAT_X_FOLDLEAN: on a table with
                // level bounds, with and without the memo, the serial fold's four instances
                // (`LEAN` x `MEMO`) give the hits and the table of the fold on a table with
                // neither; with the memo it answers retrieves, without it never.
                let mut plain_hash = new_hash();
                plain_hash.max_level = Vec::new();
                plain_hash.x_last = Vec::new();
                let mut plain_hits = Vec::new();
                fold_serial(
                    false,
                    &x,
                    &seeds,
                    &mut plain_hash,
                    &mut plain_hits,
                    &mut pause,
                    &mut FoldCounts::default(),
                );
                assert_eq!(plain_hits.iter().map(key).collect::<Vec<_>>(), serial_keys);
                assert_eq!(hash_state(&plain_hash), hash_state(&serial_hash));
                for (lean, memo) in (0..4).map(|k| (k & 2 != 0, k & 1 != 0)) {
                    let mut h = new_hash();
                    h.x_lv_start();
                    if memo {
                        h.x_lv2_start();
                    } else {
                        h.x_last = Vec::new();
                    }
                    let mut hits = Vec::new();
                    let before = super::super::super::x_diag_lv::test_memos();
                    fold_serial(
                        lean,
                        &x,
                        &seeds,
                        &mut h,
                        &mut hits,
                        &mut pause,
                        &mut FoldCounts::default(),
                    );
                    let memos = super::super::super::x_diag_lv::test_memos() - before;
                    let at = format!("case {case_no} lean {lean} memo {memo}");
                    assert_eq!(memos > 0, memo, "{at}: memo answers {memos}");
                    assert_eq!(
                        hits.iter().map(key).collect::<Vec<_>>(),
                        serial_keys,
                        "{at}"
                    );
                    assert_eq!(hash_state(&h), hash_state(&plain_hash), "{at}");
                }
                // LOSAT_X_PAIRPARFOLD (`x_pair_seed_part.rs`): the two-hit fold split into
                // 2, 3, 7, 8 and 16 partitions on 1..8 threads (and two configurations
                // without partitions), on tables without level bounds, with bounds, and with
                // bounds and the memo (LOSAT_X_DIAGLV2), in the stage's instances with and
                // without LOSAT_X_FOLDLEAN (alternating over the configurations), for two
                // stages on one scratch with a chunk that the fold without the memo folds
                // between them (as the reference fold does where a subject takes it: the memo
                // is not kept there, and the next stage makes it again), each with the hash
                // offset moved on: the hits and the table of the serial fold, and with every
                // configuration those of the table without bounds. Where the environment's
                // LOSAT_X_FOLDLEAN is the configuration's, the stage is `seed_stage` itself, so
                // that its reset of the memo after the fold without it is the one tested (a
                // debug build panics on a stale memo); elsewhere `stage_of`.
                let x2 = Inputs {
                    diag_offset: diag_offset + c.s_len as isize + window_size as isize,
                    ..x
                };
                let x3 = Inputs {
                    diag_offset: x2.diag_offset + c.s_len as isize + window_size as isize,
                    ..x
                };
                // the chunks: (inputs, folded by a stage, else by the fold without the memo)
                let chunks = [(&x, true), (&x2, false), (&x3, true)];
                let fold_memo_free =
                    |xc: &Inputs<'_>, hash: &mut DiagHashTable, hits: &mut Vec<UngappedHit>| {
                        let mut pause = 0u64;
                        fold_two_hit::<false, false, false>(
                            xc,
                            &seeds,
                            hash,
                            hits,
                            None,
                            &mut pause,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                    };
                let mut plain_expected: Vec<(Vec<_>, HashState)> = Vec::new();
                for (levels, memo) in [(false, false), (true, false), (true, true)] {
                    // whatever the environment says
                    let new_table = || {
                        let mut h = new_hash();
                        h.max_level = Vec::new();
                        h.x_last = Vec::new();
                        if levels {
                            h.x_lv_start();
                        }
                        if memo {
                            h.x_lv2_start();
                        }
                        h
                    };
                    let mut expected_hash = new_table();
                    let mut expected = Vec::new();
                    for (k, &(xc, stage)) in chunks.iter().enumerate() {
                        expected_hash.offset = xc.diag_offset as i32;
                        let mut hits = Vec::new();
                        if stage {
                            // as `seed_stage`: the memo made again at the stage's start
                            if memo {
                                expected_hash.x_lv2_reset();
                            }
                            fold_serial(
                                false,
                                xc,
                                &seeds,
                                &mut expected_hash,
                                &mut hits,
                                &mut pause,
                                &mut FoldCounts::default(),
                            );
                        } else {
                            fold_memo_free(xc, &mut expected_hash, &mut hits);
                        }
                        let keys = hits.iter().map(key).collect::<Vec<_>>();
                        if levels {
                            let (plain_keys, plain_state) = &plain_expected[k];
                            let at =
                                format!("case {case_no} levels {levels} memo {memo} chunk {k}");
                            assert_eq!(&keys, plain_keys, "{at}");
                            assert_eq!(&hash_state(&expected_hash), plain_state, "{at}");
                        } else {
                            plain_expected.push((keys.clone(), hash_state(&expected_hash)));
                        }
                        expected.push((
                            hits.iter().map(key).collect::<Vec<_>>(),
                            hash_state(&expected_hash),
                            expected_hash.max_level.clone(),
                            expected_hash.x_last.clone(),
                        ));
                    }
                    for (k, (threads, parts)) in [
                        (1usize, 2usize),
                        (2, 2),
                        (3, 3),
                        (3, 7),
                        (4, 4),
                        (5, 16),
                        (8, 8),
                        (8, 16),
                        (4, 0),
                        (6, 0),
                    ]
                    .into_iter()
                    .enumerate()
                    {
                        let lean = k % 2 == 1;
                        // LOSAT_X_DCSEEDLEAN's stage or not, in both instances of `lean`
                        let lists = (k / 2) % 2 == 0;
                        let pool = rayon::ThreadPoolBuilder::new()
                            .num_threads(threads)
                            .build()
                            .unwrap();
                        pool.install(|| {
                            let mut hash = new_table();
                            let mut scratch = PartScratch::default();
                            for (chunk, &(xc, stage)) in chunks.iter().enumerate() {
                                // EXPERIMENT (LOSAT_X_SEEDPIECE): the lists' stage in pieces
                                // of 1 seed, a few, hundreds, or the strips' lists, the
                                // stages of a configuration alternating (the real kernels
                                // stop inside words' runs and restart at many offsets)
                                let pieces = [1usize, 3 + case_no % 11, 0, 250 + 37 * k]
                                    [(k + chunk + 2 * usize::from(levels) + usize::from(memo))
                                        % 4];
                                scratch.x_pieces = Some(pieces);
                                let at = format!(
                                    "case {case_no} threads {threads} parts {parts} levels {levels} memo {memo} lean {lean} lists {lists} pieces {pieces} chunk {chunk}"
                                );
                                hash.offset = xc.diag_offset as i32;
                                let mut hits = Vec::new();
                                if stage {
                                    let fold_part = (parts > 0).then_some(part::FoldPart {
                                        scratch: &mut scratch,
                                        parts,
                                        latest: None,
                                    });
                                    let stats =
                                        if lean == lean_on() && lists == seed_lean_on() {
                                            seed_stage(
                                                xc, disc, threads, &mut hash, &mut hits, None,
                                                false, fold_part,
                                            )
                                        } else {
                                            stage_of(
                                                lean, lists, xc, disc, threads, &mut hash,
                                                &mut hits, fold_part,
                                            )
                                        };
                                    assert_eq!(stats.counts.seeds, seeds.len() as u64, "{at}");
                                } else {
                                    fold_memo_free(xc, &mut hash, &mut hits);
                                }
                                let (keys, state, max_level, x_last) = &expected[chunk];
                                assert_eq!(
                                    &hits.iter().map(key).collect::<Vec<_>>(),
                                    keys,
                                    "{at}"
                                );
                                assert_eq!(&hash_state(&hash), state, "{at}");
                                assert_eq!(&hash.max_level, max_level, "{at}");
                                assert_eq!(&hash.x_last, x_last, "{at}");
                            }
                        });
                    }
                }
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_PAIRPARDISCARR): the two-hit stage with the diagonal array (a query
// of 1500 to 3899 bases, so a block of at most 7799 letters; the array has NCBI's length,
// the power of two just above the block plus the window) for every discontiguous template:
// the pairs of the reference wrapper folded in one go by `fold_two_hit_array`, against
// `seed_stage_array` on 1..8 threads. Hit lists and both arrays must be the same, and the
// serial fold must have taken every branch: the `last_hit` rejection, the first hit of a
// window (`hit_len_array` set), the second hit both below and above the cutoff (a saved
// hit leaves `flag` 1).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn disc_two_hit_array_stage_on_strips_equals_one_scan() {
    let mut rng = Lcg(0xA22A_5EED);
    let mut case_no = 0;
    for word_size in [11usize, 12] {
        for word_type in [
            DiscWordType::Coding,
            DiscWordType::Optimal,
            DiscWordType::TwoTemplates,
        ] {
            for template_length in [16usize, 18, 21] {
                case_no += 1;
                let q_len = 1500 + rng.below(2400) as usize;
                let c = make_case_with(&mut rng, 30_000, q_len);
                // every second case: a masked subject, scanned in its ranges
                let c = if case_no % 2 == 0 {
                    mask_case(&mut rng, c)
                } else {
                    c
                };
                let lookup = build_disc_mb_lookup(
                    &ctx_refs(&c),
                    &c.offsets,
                    word_size,
                    template_length,
                    word_type,
                    &c.masks,
                    c.concat.len(),
                );
                let forms = forms_of(&lookup);
                let disc = lookup.disc().expect("a discontiguous table");
                let diag_offset = 37 + rng.below(1000) as isize;
                let window_size = 40;
                let x = Inputs {
                    use_array_indexing: true,
                    ..inputs(&c, &lookup, &forms, diag_offset, window_size)
                };
                assert!(in_scope_array(&x) && !in_scope(&x));
                let mut length = 1usize;
                while length < c.concat.len() + window_size {
                    length <<= 1;
                }
                let mut seeds = Vec::new();
                scan_subject_disc_words_with_ranges(
                    &c.packed,
                    c.s_len,
                    disc,
                    &c.seq_ranges,
                    c.masked,
                    |s_off, s_range, index, index2| {
                        lookup.for_each_hit(index, |q_off_1| {
                            if q_off_1 != 0 {
                                seed_of_pair(&x, q_off_1 as usize - 1, s_off, s_range, &mut seeds);
                            }
                        });
                        if let Some(index2) = index2 {
                            lookup.for_each_hit2(index2, |q_off_1| {
                                seed_of_pair(&x, q_off_1 as usize - 1, s_off, s_range, &mut seeds);
                            });
                        }
                    },
                );
                let mut serial_levels = vec![DiagStruct::default(); length];
                let mut serial_lens = vec![0u8; length];
                let mut serial_hits = Vec::new();
                let mut counts = FoldCounts::default();
                let mut pause = 0u64;
                fold_two_hit_array::<false>(
                    &x,
                    &seeds,
                    &mut DiagArray {
                        levels: &mut serial_levels,
                        lens: &mut serial_lens,
                        length: length as isize,
                        mask: length as isize - 1,
                        size: length,
                    },
                    &mut serial_hits,
                    None,
                    &mut pause,
                    &mut counts,
                );
                let saved = serial_levels.iter().filter(|d| d.flag != 0).count();
                let first_hits = serial_lens.iter().filter(|&&l| l > 0).count();
                assert!(
                    serial_hits.len() > 5 && counts.skipped > 0 && saved > 0 && first_hits > 0,
                    "case {case_no}: the test sequences give too few hits ({} hits, {} skipped, {saved} saved diagonals, {first_hits} first-hit diagonals)",
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
                let state = |levels: &[DiagStruct], lens: &[u8]| {
                    (
                        levels
                            .iter()
                            .map(|d| (d.last_hit, d.flag))
                            .collect::<Vec<_>>(),
                        lens.to_vec(),
                    )
                };
                let serial_keys: Vec<_> = serial_hits.iter().map(key).collect();
                let serial_state = state(&serial_levels, &serial_lens);
                for threads in 1..=8usize {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()
                        .unwrap();
                    let (got_state, keys, stats_seeds, strips) = pool.install(|| {
                        let mut levels = vec![DiagStruct::default(); length];
                        let mut lens = vec![0u8; length];
                        let mut hits = Vec::new();
                        let stats = seed_stage_array(
                            &x,
                            disc,
                            threads,
                            &mut DiagArray {
                                levels: &mut levels,
                                lens: &mut lens,
                                length: length as isize,
                                mask: length as isize - 1,
                                size: length,
                            },
                            &mut hits,
                            None,
                            false,
                        );
                        (
                            state(&levels, &lens),
                            hits.iter().map(key).collect::<Vec<_>>(),
                            stats.counts.seeds,
                            stats.strips,
                        )
                    });
                    assert!(strips > 1, "case {case_no} threads {threads}");
                    assert_eq!(
                        stats_seeds,
                        seeds.len() as u64,
                        "case {case_no} threads {threads}"
                    );
                    assert_eq!(keys, serial_keys, "case {case_no} threads {threads}");
                    assert!(
                        got_state == serial_state,
                        "case {case_no} threads {threads}: the diagonal arrays differ"
                    );
                }
            }
        }
    }
}

// The scope: a discontiguous table in two-hit mode with the diagonal hash only. A
// contiguous table in two-hit mode (megablast or blastn with `-window_size`, which the
// CLI rejects) is in neither this scope nor the parent's: the parent's `run` refuses it
// (window_size != 0) and the reference runs.
#[test]
fn scope_is_the_two_hit_discontiguous_stage() {
    let mut rng = Lcg(0x5C09E);
    let c = make_case(&mut rng, 5_000);
    let disc = build_disc_mb_lookup(
        &ctx_refs(&c),
        &c.offsets,
        11,
        18,
        DiscWordType::Coding,
        &c.masks,
        c.concat.len(),
    );
    let disc_forms = forms_of(&disc);
    let x = inputs(&c, &disc, &disc_forms, 40, 40);
    assert!(in_scope(&x));
    assert!(!in_scope(&Inputs {
        window_size: 0,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    }));
    assert!(!in_scope(&Inputs {
        use_array_indexing: true,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    }));
    assert!(!in_scope(&Inputs {
        scan_range: 1,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    }));
    // a masked subject is scanned in its ranges by the stage (the stage tests compare both)
    assert!(in_scope(&Inputs {
        subject_masked: true,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    }));
    assert!(!in_scope(&Inputs {
        diagnostics: true,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    }));
    assert!(!in_scope(&Inputs {
        one_subject_owner: false,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    }));
    let contiguous = build_two_stage_lookup(
        &ctx_refs(&c),
        &c.offsets,
        28,
        12,
        &c.masks,
        None,
        0,
        c.concat.len(),
        false,
    );
    let contiguous_forms = forms_of(&contiguous);
    assert!(!in_scope(&inputs(
        &c,
        &contiguous,
        &contiguous_forms,
        40,
        40
    )));
    // EXPERIMENT (LOSAT_X_PAIRPARDISCARR): the same scope with the diagonal array.
    fn array(x: Inputs<'_>) -> Inputs<'_> {
        Inputs {
            use_array_indexing: true,
            ..x
        }
    }
    assert!(in_scope_array(&array(inputs(
        &c,
        &disc,
        &disc_forms,
        40,
        40
    ))));
    assert!(!in_scope_array(&inputs(&c, &disc, &disc_forms, 40, 40)));
    assert!(!in_scope_array(&array(Inputs {
        window_size: 0,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    })));
    assert!(!in_scope_array(&array(Inputs {
        scan_range: 1,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    })));
    assert!(in_scope_array(&array(Inputs {
        subject_masked: true,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    })));
    assert!(!in_scope_array(&array(Inputs {
        small_na_word: true,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    })));
    assert!(!in_scope_array(&array(Inputs {
        diagnostics: true,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    })));
    assert!(!in_scope_array(&array(Inputs {
        one_subject_owner: false,
        ..inputs(&c, &disc, &disc_forms, 40, 40)
    })));
    assert!(!in_scope_array(&array(inputs(
        &c,
        &contiguous,
        &contiguous_forms,
        40,
        40
    ))));
    // The parent's `run` in a pool of two threads (where its one-hit stage would run):
    // a contiguous table with a window is refused and nothing is touched.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        pool.install(|| {
            let mut hash = DiagHashTable::new(40);
            let mut hits = Vec::new();
            let x = inputs(&c, &contiguous, &contiguous_forms, 40, 40);
            let mut fold_scratch = PartScratch::default();
            let (mut levels, mut lens) = (Vec::new(), Vec::new());
            assert!(!super::super::run(
                &x,
                &mut hash,
                DiagArray {
                    levels: &mut levels,
                    lens: &mut lens,
                    length: 0,
                    mask: 0,
                    size: 0,
                },
                &mut hits,
                &mut fold_scratch
            ));
            assert!(hits.is_empty());
            assert_eq!(hash.occupancy, 1);
        });
    }
}

// `ScanForms::choose`, the one place that chooses a batch's forms (`ScanBatch::forms`): a
// discontiguous table gets no CSR head and no ranked runs (its scans take their forms from its
// templates), and the megablast scanners of an unmasked and a masked chunk are those of
// `select_mb_scan_kind` for the lookup word length and the scan step, with a table or without.
// `chosen_from` holds for what `choose` was given, and not for another step whose scanners differ
// or another table with a CSR head, even one of the same content.
#[test]
fn choose_gives_the_forms_of_the_table_width_and_step() {
    let mut rng = Lcg(0xF04A5);
    let c = make_case(&mut rng, 2_000);
    for (word_size, template_length, word_type) in [
        (11, 16, DiscWordType::Coding),
        (12, 21, DiscWordType::Optimal),
        (11, 18, DiscWordType::TwoTemplates),
    ] {
        let disc = build_disc_mb_lookup(
            &ctx_refs(&c),
            &c.offsets,
            word_size,
            template_length,
            word_type,
            &c.masks,
            c.concat.len(),
        );
        let lut = disc.lut_word_length();
        for step in 1..=4 {
            let at = format!("disc {word_size}/{template_length} step {step}");
            let forms = ScanForms::choose(Some(&disc), lut, step);
            assert!(forms.head.is_none(), "{at}");
            assert!(!forms.rank_once, "{at}");
            for masked in [false, true] {
                assert_eq!(
                    forms.scan_kind(masked),
                    select_mb_scan_kind(lut, step, masked),
                    "{at} masked {masked}"
                );
            }
            assert!(forms.chosen_from(Some(&disc), lut, step), "{at}");
        }
    }
    for lut in 7..=13 {
        for step in 1..=5 {
            let at = format!("no table, lookup width {lut} step {step}");
            let forms = ScanForms::choose(None, lut, step);
            assert!(forms.head.is_none(), "{at}");
            assert!(!forms.rank_once, "{at}");
            for masked in [false, true] {
                assert_eq!(
                    forms.scan_kind(masked),
                    select_mb_scan_kind(lut, step, masked),
                    "{at} masked {masked}"
                );
            }
            assert!(forms.chosen_from(None, lut, step), "{at}");
        }
    }
    // lookup width 11: steps 1 and 2 take different scanners
    let build = || {
        build_two_stage_lookup(
            &ctx_refs(&c),
            &c.offsets,
            28,
            11,
            &c.masks,
            None,
            0,
            c.concat.len(),
            false,
        )
    };
    let (contiguous, same_content) = (build(), build());
    let forms = forms_of(&contiguous);
    assert!(forms.chosen_from(Some(&contiguous), 11, 1));
    assert!(!forms.chosen_from(Some(&contiguous), 11, 2));
    assert_eq!(forms.head.is_some(), contiguous.x_mb_csr().is_some());
    if forms.head.is_some() {
        assert!(!forms.chosen_from(Some(&same_content), 11, 1));
    }
}

// The parent's `run` checks in debug builds that its inputs' forms were chosen from their table
// and scan step (`ScanForms::chosen_from`): forms of another step stop it before any scan.
#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "forms not chosen from the inputs' table and scan step")]
fn run_stops_at_forms_of_another_step() {
    let mut rng = Lcg(0x5C09E);
    let c = make_case(&mut rng, 2_000);
    let contiguous = build_two_stage_lookup(
        &ctx_refs(&c),
        &c.offsets,
        28,
        11,
        &c.masks,
        None,
        0,
        c.concat.len(),
        false,
    );
    let forms = forms_of(&contiguous);
    let x = Inputs {
        scan_step: 2,
        ..inputs(&c, &contiguous, &forms, 0, 0)
    };
    let mut hash = DiagHashTable::new(0);
    let (mut levels, mut lens) = (Vec::new(), Vec::new());
    super::super::run(
        &x,
        &mut hash,
        DiagArray {
            levels: &mut levels,
            lens: &mut lens,
            length: 0,
            mask: 0,
            size: 0,
        },
        &mut Vec::new(),
        &mut PartScratch::default(),
    );
}

/// EXPERIMENT (LOSAT_X_DCLATEST): a fold's table that logs what the fold reads and writes: per
/// retrieve the diagonal, the position and the decisions it makes (`LvSeeds::decisions`:
/// rejected, two-hit branch), per insert the diagonal, level, saved flag and position.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct Traced<'t, T> {
    inner: T,
    seeds: super::super::super::x_diag_lv::LvSeeds,
    reads: &'t mut Vec<(i32, usize, (bool, bool))>,
    inserts: &'t mut Vec<(i32, i32, bool, i32)>,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl<T: FoldTable> FoldTable for Traced<'_, T> {
    const AHEAD: bool = T::AHEAD;

    fn ahead(&self, seed: &Seed) {
        self.inner.ahead(seed);
    }

    fn retrieve(&mut self, diag: i32, s_off_pos: usize) -> (i32, bool) {
        let (level, saved) = self.inner.retrieve(diag, s_off_pos);
        let decisions = self.seeds.decisions(Some((level, 0, saved)), s_off_pos);
        self.reads.push((diag, s_off_pos, decisions));
        (level, saved)
    }

    fn retrieve_off(&self, diag: i32) -> Option<(i32, i32, bool)> {
        self.inner.retrieve_off(diag)
    }

    fn insert(
        &mut self,
        out: &mut part::PartOut,
        key: u64,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) {
        self.inserts.push((diag, level, hit_saved, s_off));
        self.inner.insert(
            out,
            key,
            diag,
            level,
            hit_len,
            hit_saved,
            s_off,
            window_size,
        );
    }
}

/// The fields of a hit the stage tests compare.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn hit_keys(hits: &[UngappedHit]) -> Vec<(u32, u32, u32, i32, i32, u32, u32, i32)> {
    hits.iter()
        .map(|h| {
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
        })
        .collect()
}

/// The pairs of the reference wrapper `scan_subject_disc_words_with_ranges` as seeds (one scan of
/// the whole range; a masked subject's ranges one by one).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn all_seeds(c: &Case, lookup: &TwoStageLookup, x: &Inputs<'_>) -> Vec<Seed> {
    let disc = lookup.disc().expect("a discontiguous table");
    let mut seeds = Vec::new();
    scan_subject_disc_words_with_ranges(
        &c.packed,
        c.s_len,
        disc,
        &c.seq_ranges,
        c.masked,
        |s_off, s_range, index, index2| {
            lookup.for_each_hit(index, |q_off_1| {
                if q_off_1 != 0 {
                    seed_of_pair(x, q_off_1 as usize - 1, s_off, s_range, &mut seeds);
                }
            });
            if let Some(index2) = index2 {
                lookup.for_each_hit2(index2, |q_off_1| {
                    seed_of_pair(x, q_off_1 as usize - 1, s_off, s_range, &mut seeds);
                });
            }
        },
    );
    seeds
}

/// The windows of the map's tests for a word length `wl`: the band 2w + 2 < wl and around it,
/// around the word length, and the default 40 and above.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn latest_windows(wl: usize) -> Vec<usize> {
    let mut w: Vec<usize> = (1..=10).collect();
    w.extend([wl - 1, wl, wl + 1, 40, 64]);
    w.sort_unstable();
    w.dedup();
    w
}

/// A diagonal hash without level bounds or memo (whatever the environment says).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn plain_table(window: usize, offset: i32) -> DiagHashTable {
    let mut h = DiagHashTable::new(window as i32);
    h.max_level = Vec::new();
    h.x_last = Vec::new();
    h.offset = offset;
    h
}

/// `fold_two_hit_t::<PT, LEAN, MEMO, TRACK>` (the chains' own fold) against `fold_two_hit_on` on
/// `Chains` of the same constants, from equal tables: the hits, the table (with its level bounds
/// and memo), the counts, and what a partition's fold writes to `out`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn copies_agree_in<const PT: bool, const LEAN: bool, const MEMO: bool, const TRACK: bool>(
    x: &Inputs<'_>,
    seeds: &[Seed],
    w: usize,
    offset: i32,
    levels: bool,
    at: &str,
) {
    let table = || {
        let mut h = plain_table(w, offset);
        if levels {
            h.x_lv_start();
        }
        if MEMO {
            h.x_lv2_start();
        }
        h
    };
    let pos: Vec<u32> = (0..seeds.len() as u32).collect();
    let at = format!("{at}: PT {PT} LEAN {LEAN} MEMO {MEMO} TRACK {TRACK} levels {levels}");
    let (mut own, mut own_hits, mut own_out) =
        (table(), Vec::new(), part::PartOut::of_test_strip(TRACK));
    let mut own_counts = FoldCounts::default();
    fold_two_hit_t::<PT, LEAN, MEMO, TRACK>(
        x,
        seeds,
        &mut own,
        &mut own_hits,
        None,
        &mut 0,
        &mut own_counts,
        &mut own_out,
        if PT { &pos } else { &[] },
    );
    let (mut hash, mut hits, mut out) = (table(), Vec::new(), part::PartOut::of_test_strip(TRACK));
    let mut counts = FoldCounts::default();
    fold_two_hit_on::<PT, LEAN, _>(
        x,
        seeds,
        &mut Chains::<PT, LEAN, MEMO, TRACK>::new(x, &mut hash),
        &mut hits,
        None,
        &mut 0,
        &mut counts,
        &mut out,
        if PT { &pos } else { &[] },
    );
    assert_eq!(hit_keys(&own_hits), hit_keys(&hits), "{at}: hits");
    assert!(hash_state(&own) == hash_state(&hash), "{at}: table");
    assert_eq!(own.max_level, hash.max_level, "{at}: level bounds");
    assert_eq!(own.x_last, hash.x_last, "{at}: memo");
    assert_eq!(
        (own_counts.seeds, own_counts.skipped),
        (counts.seeds, counts.skipped),
        "{at}: counts"
    );
    assert_eq!(own_out.recorded(), out.recorded(), "{at}: out");
}

/// `copies_agree_in` for every instance the production folds use: `fold_two_hit` (PT either,
/// LEAN either, MEMO either, TRACK on) and the partitions' lists (`seed_stage_lists`: PT, LEAN
/// either, MEMO either, TRACK off); a table with level bounds where it has a memo, with and
/// without them where it has none.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn copies_agree_in_every_instance(x: &Inputs<'_>, seeds: &[Seed], w: usize, offset: i32, at: &str) {
    macro_rules! each {
        ($pt:literal, $track:literal) => {
            for levels in [false, true] {
                copies_agree_in::<$pt, false, false, $track>(x, seeds, w, offset, levels, at);
                copies_agree_in::<$pt, true, false, $track>(x, seeds, w, offset, levels, at);
            }
            copies_agree_in::<$pt, false, true, $track>(x, seeds, w, offset, true, at);
            copies_agree_in::<$pt, true, true, $track>(x, seeds, w, offset, true, at);
        };
    }
    each!(false, true);
    each!(true, true);
    each!(true, false);
}

// EXPERIMENT (LOSAT_X_DCLATEST): seed by seed, the fold on the map of the latest inserts makes the
// decisions (rejected, two-hit branch) and the inserts of the fold on the chains, and gives its
// hits (and that fold on the chains gives the hits and the table of `fold_two_hit_t`, the chains'
// own fold), for every template and windows from 1 to 64: with NCBI's take-over window in the scope (w
// + 1, module comment (iii) of `x_dc_latest.rs`) at every window, with LOSAT's
// (`diag_hash_insert_window`) outside the band 2w + 2 < wl; inside the band, and wherever w < wl,
// neither fold extends a seed (module comment (v)).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_fold_on_the_map_makes_the_decisions_of_the_chains() {
    let mut rng = Lcg(0x1A7E_57D1);
    let mut case_no = 0;
    let mut band = 0;
    for word_size in [11usize, 12] {
        for word_type in [
            DiscWordType::Coding,
            DiscWordType::Optimal,
            DiscWordType::TwoTemplates,
        ] {
            for template_length in [16usize, 18, 21] {
                case_no += 1;
                let c = make_case(&mut rng, 20_000);
                let c = if case_no % 2 == 0 {
                    mask_case(&mut rng, c)
                } else {
                    c
                };
                let lookup = build_disc_mb_lookup(
                    &ctx_refs(&c),
                    &c.offsets,
                    word_size,
                    template_length,
                    word_type,
                    &c.masks,
                    c.concat.len(),
                );
                let forms = forms_of(&lookup);
                let wl = template_length;
                let diag_offset = 37 + rng.below(1000) as isize;
                let seeds = all_seeds(&c, &lookup, &inputs(&c, &lookup, &forms, diag_offset, 40));
                for w in latest_windows(wl) {
                    let x = inputs(&c, &lookup, &forms, diag_offset, w);
                    let lv = super::super::super::x_diag_lv::LvSeeds::new(true, wl, w);
                    for ncbi in [true, false] {
                        let at = format!("case {case_no} wl {wl} w {w} ncbi {ncbi}");
                        let xw = Inputs {
                            diag_hash_window: if ncbi {
                                w as i32 + 1
                            } else {
                                x.diag_hash_window
                            },
                            ..x
                        };
                        let mut hash = plain_table(w, diag_offset as i32);
                        let (mut reads, mut inserts, mut hits) =
                            (Vec::new(), Vec::new(), Vec::new());
                        fold_two_hit_on::<false, false, _>(
                            &xw,
                            &seeds,
                            &mut Traced {
                                inner: Chains::<false, false, false, true>::new(&xw, &mut hash),
                                seeds: lv,
                                reads: &mut reads,
                                inserts: &mut inserts,
                            },
                            &mut hits,
                            None,
                            &mut 0,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                        // the chains' own fold (`fold_two_hit_t`, written out) is that fold
                        let mut own = plain_table(w, diag_offset as i32);
                        let mut own_hits = Vec::new();
                        fold_two_hit_t::<false, false, false, true>(
                            &xw,
                            &seeds,
                            &mut own,
                            &mut own_hits,
                            None,
                            &mut 0,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                        assert_eq!(hit_keys(&own_hits), hit_keys(&hits), "{at}: fold_two_hit_t");
                        assert!(
                            hash_state(&own) == hash_state(&hash),
                            "{at}: fold_two_hit_t's table"
                        );
                        // and so in every instance of the production folds
                        copies_agree_in_every_instance(&xw, &seeds, w, diag_offset as i32, &at);
                        let mut map = DiagLatest::default();
                        map.begin(lv.bound(diag_offset as usize));
                        let (mut map_reads, mut map_inserts, mut map_hits) =
                            (Vec::new(), Vec::new(), Vec::new());
                        fold_two_hit_on::<false, false, _>(
                            &xw,
                            &seeds,
                            &mut Traced {
                                inner: Latest::new(&xw, &mut map),
                                seeds: lv,
                                reads: &mut map_reads,
                                inserts: &mut map_inserts,
                            },
                            &mut map_hits,
                            None,
                            &mut 0,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                        if ncbi || 2 * w + 2 >= wl {
                            assert_eq!(map_reads, reads, "{at}");
                            assert_eq!(map_inserts, inserts, "{at}");
                            assert_eq!(hit_keys(&map_hits), hit_keys(&hits), "{at}");
                        } else {
                            band += 1;
                        }
                        if w < wl {
                            assert!(hits.is_empty() && map_hits.is_empty(), "{at}");
                        } else if w == 40 {
                            assert!(hits.len() > 5, "{at}: too few hits");
                        }
                    }
                }
            }
        }
    }
    assert!(band > 0);
}

// EXPERIMENT (LOSAT_X_DCLATEST): the stage on the maps (`seed_stage_latest`) gives the hits of the
// reference fold on the chains, chunk after chunk on one scratch: the owner's map (1 and 4 threads
// without partitions) and 2, 3, 7, 8 and 16 partitions on 1..8 threads, with and without
// LOSAT_X_FOLDLEAN and LOSAT_X_DCSEEDLEAN's lists (in pieces of 1 seed, a few, hundreds, or the
// strips' lists), at windows from the band to 64, on plain and masked subjects. The chunks move the
// offset as NCBI does (`advance_diag_table_offset`), from below a quarter of the largest Int4 past it, where it is
// cleared (and the maps are emptied at the next stage's start). Outside the band (module comment
// (v) of `x_dc_latest.rs`) the rejected seeds are the chains' too; the chunk's table is never read
// or changed.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn disc_two_hit_stage_on_maps_gives_the_hits_of_the_chains() {
    let mut rng = Lcg(0x0D1A_6111);
    let mut case_no = 0;
    for word_size in [11usize, 12] {
        for word_type in [
            DiscWordType::Coding,
            DiscWordType::Optimal,
            DiscWordType::TwoTemplates,
        ] {
            for template_length in [16usize, 18, 21] {
                case_no += 1;
                let c = make_case(&mut rng, 20_000);
                let c = if case_no % 2 == 0 {
                    mask_case(&mut rng, c)
                } else {
                    c
                };
                let lookup = build_disc_mb_lookup(
                    &ctx_refs(&c),
                    &c.offsets,
                    word_size,
                    template_length,
                    word_type,
                    &c.masks,
                    c.concat.len(),
                );
                let forms = forms_of(&lookup);
                let disc = lookup.disc().expect("a discontiguous table");
                let wl = template_length;
                let seeds = all_seeds(&c, &lookup, &inputs(&c, &lookup, &forms, 40, 40));
                let mut windows = vec![1 + case_no % 6, wl - 1, wl + 1, 40, 64];
                windows.dedup();
                for w in windows {
                    let x0 = inputs(&c, &lookup, &forms, w as isize, w);
                    // the reference: the chunks' offsets, hits and rejected seeds
                    let mut hash = plain_table(w, w as i32);
                    let mut offsets = Vec::new();
                    let mut expected = Vec::new();
                    for chunk in 0..5 {
                        if chunk == 2 {
                            // as after many chunks
                            hash.offset = i32::MAX / 4 - 3;
                        }
                        let xc = Inputs {
                            diag_offset: hash.offset as isize,
                            ..x0
                        };
                        let (mut hits, mut counts) = (Vec::new(), FoldCounts::default());
                        fold_serial(
                            false,
                            &xc,
                            &seeds,
                            &mut hash,
                            &mut hits,
                            &mut 0,
                            &mut counts,
                        );
                        offsets.push(hash.offset);
                        expected.push((hit_keys(&hits), counts.skipped));
                        super::super::super::advance_diag_table_offset(
                            &mut 0,
                            w as i32,
                            c.s_len,
                            false,
                            &mut [],
                            &mut [],
                            &mut hash,
                        );
                    }
                    assert_eq!(offsets[4], w as i32, "the offset was not cleared");
                    if w >= wl {
                        assert!(expected.iter().all(|(hits, _)| hits.len() > 5));
                    }
                    for (k, (threads, parts)) in [
                        (1usize, None),
                        (4, None),
                        (1, Some(2usize)),
                        (2, Some(2)),
                        (3, Some(3)),
                        (3, Some(7)),
                        (5, Some(16)),
                        (8, Some(8)),
                    ]
                    .into_iter()
                    .enumerate()
                    {
                        let lean = k % 2 == 1;
                        let lists = (k / 2) % 2 == 0;
                        let pool = rayon::ThreadPoolBuilder::new()
                            .num_threads(threads)
                            .build()
                            .unwrap();
                        pool.install(|| {
                            let mut scratch = PartScratch::default();
                            let mut table = plain_table(w, w as i32);
                            let before = hash_state(&table);
                            for (chunk, &offset) in offsets.iter().enumerate() {
                                let pieces = [1usize, 3 + case_no % 11, 0, 250 + 37 * k]
                                    [(k + chunk) % 4];
                                scratch.x_pieces = Some(pieces);
                                let at = format!(
                                    "case {case_no} wl {wl} w {w} threads {threads} parts {parts:?} lean {lean} lists {lists} pieces {pieces} chunk {chunk}"
                                );
                                let xc = Inputs {
                                    diag_offset: offset as isize,
                                    ..x0
                                };
                                let bound0 = latest_bound(&xc).expect("a chunk with bounds");
                                let mut hits = Vec::new();
                                let stats = seed_stage_latest(
                                    &xc,
                                    disc,
                                    threads,
                                    &mut table,
                                    &mut hits,
                                    None,
                                    false,
                                    &mut scratch,
                                    parts,
                                    bound0,
                                    lean,
                                    lists,
                                );
                                let (keys, skipped) = &expected[chunk];
                                assert_eq!(&hit_keys(&hits), keys, "{at}");
                                assert_eq!(stats.counts.seeds, seeds.len() as u64, "{at}");
                                if 2 * w + 2 >= wl {
                                    assert_eq!(stats.counts.skipped, *skipped, "{at}");
                                }
                            }
                            assert_eq!(hash_state(&table), before, "the chunk's table changed");
                        });
                    }
                }
            }
        }
    }
}

// EXPERIMENT (LOSAT_X_DCLATEST): masked subjects whose scan ranges are out of order, as a subject
// of `-db` with `-db_soft_mask` keeps its stored masks' order (`subject_split.rs`): `mask_case`'s
// ranges rotated, reversed, two swapped, or one starting a position before the last position of
// the one before. The seeds' S then steps back, and the decision of `run` (`latest_bound`, then
// `latest_ranges_ascend`) sends the chunk to the chains. At the edges of the rule the ranges in
// order, unmasked, touching (a range starting at the last position of the one before) or with a
// range too short to scan placed out of order keep S from stepping back and take the maps. Seed by
// seed the fold the decision picks makes the reads, inserts and hits of the reference fold on the
// chains over the seeds in scan order (NCBI's); this holds for the first kind only because the
// decision picks the chains (`put` asserts the order of S in debug builds). The parent's `run`
// on 2 and 4 threads then gives the reference's hits, and the reference's table where the chunk
// went to the chains (the maps' stage leaves the table as it was).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn a_chunk_whose_ranges_step_back_goes_to_the_chains() {
    let mut rng = Lcg(0x0BAC_C5EE);
    let (mut to_chains, mut to_maps) = (0, 0);
    for (case_no, (word_size, word_type, template_length)) in [
        (11usize, DiscWordType::Coding, 16usize),
        (12, DiscWordType::Optimal, 18),
        (11, DiscWordType::TwoTemplates, 21),
        (12, DiscWordType::Coding, 21),
    ]
    .into_iter()
    .enumerate()
    {
        let plain = make_case(&mut rng, 20_000);
        let mut c = mask_case(&mut rng, plain);
        let lookup = build_disc_mb_lookup(
            &ctx_refs(&c),
            &c.offsets,
            word_size,
            template_length,
            word_type,
            &c.masks,
            c.concat.len(),
        );
        let forms = forms_of(&lookup);
        let tl = template_length as i32;
        let sorted = c.seq_ranges.clone();
        // two neighbouring ranges with seeds, both long enough to scan
        let sorted_seeds = all_seeds(&c, &lookup, &inputs(&c, &lookup, &forms, 40, 40));
        let seeded = |(left, right): (i32, i32)| {
            right - left > 2 * tl
                && sorted_seeds
                    .iter()
                    .any(|seed| (left..=right - tl).contains(&(seed.s as i32)))
        };
        let k = (1..sorted.len() - 2)
            .find(|&k| seeded(sorted[k]) && seeded(sorted[k + 1]))
            .expect("two ranges with seeds");
        fn with(r: &[(i32, i32)], f: impl FnOnce(&mut Vec<(i32, i32)>)) -> Vec<(i32, i32)> {
            let mut r = r.to_vec();
            f(&mut r);
            r
        }
        // (the variant, its ranges, masked, whether its ranges step back)
        let variants: Vec<(&str, Vec<(i32, i32)>, bool, bool)> = vec![
            ("in order", sorted.clone(), true, false),
            ("unmasked", vec![(0, c.s_len as i32)], false, false),
            (
                "rotated",
                with(&sorted, |r| {
                    let half = r.len() / 2;
                    r.rotate_left(half)
                }),
                true,
                true,
            ),
            ("reversed", with(&sorted, |r| r.reverse()), true, true),
            ("swapped", with(&sorted, |r| r.swap(k, k + 1)), true, true),
            (
                "one back",
                with(&sorted, |r| r[k + 1].0 = r[k].1 - tl - 1),
                true,
                true,
            ),
            (
                "touching",
                with(&sorted, |r| r[k + 1].0 = r[k].1 - tl),
                true,
                false,
            ),
            (
                "short one back",
                with(&sorted, |r| r.insert(k + 2, (r[k].0, r[k].0 + tl - 1))),
                true,
                false,
            ),
        ];
        for (name, ranges, masked, steps_back) in variants {
            c.seq_ranges = ranges;
            c.masked = masked;
            let diag_offset = 37 + rng.below(1000) as isize;
            let seeds = all_seeds(&c, &lookup, &inputs(&c, &lookup, &forms, diag_offset, 40));
            // S steps back only where the ranges do (and does, with the seeds of these sequences,
            // where whole ranges come out of order)
            let seeds_back = seeds.windows(2).any(|p| p[1].s < p[0].s);
            assert!(
                !seeds_back || steps_back,
                "case {case_no} {name}: the seeds' order"
            );
            if matches!(name, "rotated" | "reversed" | "swapped") {
                assert!(seeds_back, "case {case_no} {name}: the seeds' order");
            }
            for w in [template_length + 1, 40] {
                let at = format!("case {case_no} {name} w {w}");
                let x = inputs(&c, &lookup, &forms, diag_offset, w);
                let lv = super::super::super::x_diag_lv::LvSeeds::new(true, template_length, w);
                assert_eq!(
                    latest_ranges_ascend(&x, template_length),
                    !steps_back,
                    "{at}"
                );
                let decided =
                    latest_bound(&x).filter(|_| latest_ranges_ascend(&x, template_length));
                assert_eq!(decided.is_none(), steps_back, "{at}: the decision");
                // the reference: the chains, seed by seed
                let mut hash = plain_table(w, diag_offset as i32);
                let (mut reads, mut inserts, mut hits) = (Vec::new(), Vec::new(), Vec::new());
                fold_two_hit_on::<false, false, _>(
                    &x,
                    &seeds,
                    &mut Traced {
                        inner: Chains::<false, false, false, true>::new(&x, &mut hash),
                        seeds: lv,
                        reads: &mut reads,
                        inserts: &mut inserts,
                    },
                    &mut hits,
                    None,
                    &mut 0,
                    &mut FoldCounts::default(),
                    &mut part::PartOut::default(),
                    &[],
                );
                let reference = hash_state(&hash);
                if w == 40 {
                    assert!(hits.len() > 5, "{at}: too few hits");
                }
                // the fold the decision picks
                let (mut got_reads, mut got_inserts, mut got_hits) =
                    (Vec::new(), Vec::new(), Vec::new());
                let mut got_hash = plain_table(w, diag_offset as i32);
                let mut map = DiagLatest::default();
                match decided {
                    Some(bound0) => {
                        to_maps += 1;
                        map.begin(bound0);
                        fold_two_hit_on::<false, false, _>(
                            &x,
                            &seeds,
                            &mut Traced {
                                inner: Latest::new(&x, &mut map),
                                seeds: lv,
                                reads: &mut got_reads,
                                inserts: &mut got_inserts,
                            },
                            &mut got_hits,
                            None,
                            &mut 0,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                    }
                    None => {
                        to_chains += 1;
                        fold_two_hit_on::<false, false, _>(
                            &x,
                            &seeds,
                            &mut Traced {
                                inner: Chains::<false, false, false, true>::new(&x, &mut got_hash),
                                seeds: lv,
                                reads: &mut got_reads,
                                inserts: &mut got_inserts,
                            },
                            &mut got_hits,
                            None,
                            &mut 0,
                            &mut FoldCounts::default(),
                            &mut part::PartOut::default(),
                            &[],
                        );
                    }
                }
                assert_eq!(got_reads, reads, "{at}: reads");
                assert_eq!(got_inserts, inserts, "{at}: inserts");
                assert_eq!(hit_keys(&got_hits), hit_keys(&hits), "{at}: hits");
                // the decision point itself
                for threads in [2usize, 4] {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()
                        .unwrap();
                    pool.install(|| {
                        let mut table = plain_table(w, diag_offset as i32);
                        let before = hash_state(&table);
                        let mut stage_hits = Vec::new();
                        let mut fold_scratch = PartScratch::default();
                        let (mut levels, mut lens) = (Vec::new(), Vec::new());
                        let ran = run(
                            &x,
                            &mut table,
                            DiagArray {
                                levels: &mut levels,
                                lens: &mut lens,
                                length: 0,
                                mask: 0,
                                size: 0,
                            },
                            &mut stage_hits,
                            &mut fold_scratch,
                        );
                        assert_eq!(ran, on(), "{at} threads {threads}");
                        if ran {
                            assert_eq!(
                                hit_keys(&stage_hits),
                                hit_keys(&hits),
                                "{at} threads {threads}: the stage's hits"
                            );
                            let on_maps = latest_on() && !steps_back;
                            assert_eq!(
                                hash_state(&table),
                                if on_maps { before } else { reference.clone() },
                                "{at} threads {threads}: the stage's table"
                            );
                        }
                    });
                }
            }
        }
    }
    assert!(to_chains > 0 && to_maps > 0);
}

/// EXPERIMENT (LOSAT_X_DCFOLDREG): what one fold call on a map did: its trace (`latest::Event`:
/// each seed's retrieve and insert), the map after it (`DiagLatest::state`), its hits, its counts
/// (seeds, skipped), what a partition's fold wrote to `out`, and its ungapped extensions
/// (LOSAT_TIMING's count, without partitions).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[derive(PartialEq)]
struct MapFold {
    trace: Vec<latest::Event>,
    state: latest::MapState,
    hits: Vec<(u32, u32, u32, i32, i32, u32, u32, i32)>,
    counts: (u64, u64),
    out: (Vec<u64>, Vec<(u64, u32)>, Vec<u16>, u64),
    extensions: u64,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl MapFold {
    /// Where `self` and `other` first differ (nothing when they are equal).
    fn difference(&self, other: &MapFold) -> Option<String> {
        if let Some(k) = (0..self.trace.len().max(other.trace.len()))
            .find(|&k| self.trace.get(k) != other.trace.get(k))
        {
            return Some(format!(
                "trace step {k}: {:?} against {:?}",
                self.trace.get(k),
                other.trace.get(k)
            ));
        }
        let (a, b) = (&self.state, &other.state);
        if let Some(k) = (0..a.0.len().max(b.0.len())).find(|&k| a.0.get(k) != b.0.get(k)) {
            return Some(format!(
                "bucket {k}: {:?} against {:?}",
                a.0.get(k),
                b.0.get(k)
            ));
        }
        if a != b {
            return Some(format!(
                "map: stash {:?} bits {} shift {} entries {} top {} against stash {:?} bits {} shift {} entries {} top {}",
                a.1, a.2, a.3, a.4, a.5, b.1, b.2, b.3, b.4, b.5
            ));
        }
        if self.hits != other.hits {
            return Some(format!("hits {:?} against {:?}", self.hits, other.hits));
        }
        if self != other {
            return Some(format!(
                "counts {:?} out {:?} extensions {} against counts {:?} out {:?} extensions {}",
                self.counts, self.out, self.extensions, other.counts, other.out, other.extensions
            ));
        }
        None
    }
}

/// One fold call on `map` (with its trace) over `seeds`, seed `k` at place `pos[k]` of the stage
/// (`PT`): `fold_reg::fold_two_hit_reg` with `TOW` = `LEAN` (`reg`; the word length of a
/// discontiguous table is its lookup width), or the fold of LOSAT_X_DCFOLDREG=0, `fold_two_hit_on`
/// on `Latest` with `LEAN`. Timed, so that the extensions are counted.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn map_fold<const PT: bool, const LEAN: bool>(
    reg: bool,
    x: &Inputs<'_>,
    seeds: &[Seed],
    pos: &[u32],
    map: &mut DiagLatest,
) -> MapFold {
    assert_eq!(x.two_stage.word_length(), x.two_stage.lut_word_length());
    map.trace = Some(Vec::new());
    let timing = super::super::super::BlastnTiming::new();
    let (mut hits, mut counts) = (Vec::new(), FoldCounts::default());
    let mut out = part::PartOut::of_test_strip(false);
    let pos = if PT { pos } else { &[] };
    if reg {
        fold_reg::fold_two_hit_reg::<PT, LEAN>(
            x,
            seeds,
            map,
            &mut hits,
            Some(&timing),
            &mut 0,
            &mut counts,
            &mut out,
            pos,
        );
    } else {
        fold_two_hit_on::<PT, LEAN, _>(
            x,
            seeds,
            &mut Latest::new(x, map),
            &mut hits,
            Some(&timing),
            &mut 0,
            &mut counts,
            &mut out,
            pos,
        );
    }
    MapFold {
        trace: map.trace.take().expect("the trace"),
        state: map.state(),
        hits: hit_keys(&hits),
        counts: (counts.seeds, counts.skipped),
        out: out.recorded(),
        extensions: timing.ungapped_calls.load(Ordering::Relaxed),
    }
}

/// What the folds of `reg_agrees` met: seeds that found an entry (in a slot or the stash), found
/// one in the stash, were rejected, were extended, and made a hit; the most buckets and stash
/// entries a map had, and the stages that emptied it.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[derive(Clone, Copy, Debug, Default)]
struct RegCover {
    found: u64,
    stash_found: u64,
    stash_put: u64,
    rejected: u64,
    extended: u64,
    hits: u64,
    most_buckets: usize,
    most_stash: usize,
    calls: u64,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl RegCover {
    fn add(&mut self, f: &MapFold) {
        let mut stashed = None;
        for e in &f.trace {
            match *e {
                latest::Event::Read {
                    diag,
                    s,
                    level,
                    stash,
                    ..
                } => {
                    self.found += u64::from(level != 0 || stash);
                    self.stash_found += u64::from(stash);
                    stashed = stash.then_some((diag, s));
                }
                latest::Event::Put { diag, s, .. } => {
                    self.stash_put += u64::from(stashed == Some((diag, s)));
                    stashed = None;
                }
            }
        }
        self.rejected += f.counts.1;
        self.extended += f.extensions + f.out.3;
        self.hits += f.hits.len() as u64;
        self.most_buckets = self.most_buckets.max(f.state.0.len());
        self.most_stash = self.most_stash.max(f.state.1.len());
        self.calls += 1;
    }

    fn sum(mut self, o: RegCover) -> RegCover {
        self.found += o.found;
        self.stash_found += o.stash_found;
        self.stash_put += o.stash_put;
        self.rejected += o.rejected;
        self.extended += o.extended;
        self.hits += o.hits;
        self.most_buckets = self.most_buckets.max(o.most_buckets);
        self.most_stash = self.most_stash.max(o.most_stash);
        self.calls += o.calls;
        self
    }
}

/// The hash offsets of `stages` chunks of `s_len` letters at window `w` from `first`, as NCBI moves
/// the offset after each chunk (`advance_diag_table_offset`), with the second chunk's offset set
/// just below a quarter of the largest Int4, so that the offset is cleared after it (and the maps
/// are emptied at the next stage's start).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn stage_offsets(first: i32, w: usize, s_len: usize, stages: usize) -> Vec<i32> {
    let mut hash = plain_table(w, first);
    let mut offsets = Vec::new();
    for stage in 0..stages {
        if stage == 1 {
            hash.offset = i32::MAX / 4 - 3;
        }
        offsets.push(hash.offset);
        super::super::super::advance_diag_table_offset(
            &mut 0,
            w as i32,
            s_len,
            false,
            &mut [],
            &mut [],
            &mut hash,
        );
    }
    offsets
}

/// `fold_two_hit_reg` (`PT`, `TOW` = `LEAN`) against the fold of LOSAT_X_DCFOLDREG=0 on maps of
/// their own: the stages at `offsets` (each map begun at the stage's first bound, `latest_bound`),
/// each stage's `seeds` folded in calls cut at random places (one call, a few, or one seed per call
/// for a stretch), and after every call the same trace (every retrieve and insert, seed by seed),
/// map, hits, counts and `out`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn reg_agrees<const PT: bool, const LEAN: bool>(
    x0: &Inputs<'_>,
    seeds: &[Seed],
    offsets: &[i32],
    rng: &mut Lcg,
    at: &str,
) -> RegCover {
    let mut cover = RegCover::default();
    let (mut map, mut reg_map) = (DiagLatest::default(), DiagLatest::default());
    let pos: Vec<u32> = (0..seeds.len() as u32).collect();
    for (stage, &offset) in offsets.iter().enumerate() {
        let x = Inputs {
            diag_offset: offset as isize,
            ..*x0
        };
        let bound0 = latest_bound(&x).expect("a chunk with bounds");
        map.begin(bound0);
        reg_map.begin(bound0);
        let n = seeds.len();
        let mut cuts: Vec<usize> = match rng.below(3) {
            0 => vec![],
            1 => (0..1 + rng.below(8))
                .map(|_| rng.below(n as u64) as usize)
                .collect(),
            _ => {
                let from = rng.below(n as u64) as usize;
                (from..(from + 40).min(n)).collect()
            }
        };
        cuts.extend([0, n]);
        cuts.sort_unstable();
        cuts.dedup();
        for (call, part) in cuts.windows(2).enumerate() {
            let (from, to) = (part[0], part[1]);
            let at =
                format!("{at} PT {PT} LEAN {LEAN} stage {stage} call {call} seeds {from}..{to}");
            let expected =
                map_fold::<PT, LEAN>(false, &x, &seeds[from..to], &pos[from..to], &mut map);
            let got =
                map_fold::<PT, LEAN>(true, &x, &seeds[from..to], &pos[from..to], &mut reg_map);
            if let Some(difference) = expected.difference(&got) {
                panic!("{at}: {difference}");
            }
            cover.add(&expected);
        }
    }
    cover
}

/// `reg_agrees` in the four instances the stage runs (`PT` either, `LEAN` either).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn reg_agrees_in_every_instance(
    x0: &Inputs<'_>,
    seeds: &[Seed],
    offsets: &[i32],
    rng: &mut Lcg,
    at: &str,
) -> RegCover {
    reg_agrees::<false, true>(x0, seeds, offsets, rng, at)
        .sum(reg_agrees::<false, false>(x0, seeds, offsets, rng, at))
        .sum(reg_agrees::<true, true>(x0, seeds, offsets, rng, at))
        .sum(reg_agrees::<true, false>(x0, seeds, offsets, rng, at))
}

// EXPERIMENT (LOSAT_X_DCFOLDREG): the fold with the map in registers makes every retrieve, decision
// and insert of the fold on the map of LOSAT_X_DCFOLDREG=0, seed by seed, and leaves the same map,
// hits and counts after every call, on the seeds of the stage (the reference scan's pairs) for
// every template, at windows in the band 2w + 2 < wl, around the word length and above it, over
// stages through the offset's clear, cut into calls at random places, in the four instances of the
// stage (with and without partitions, with and without LOSAT_X_FOLDLEAN).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_fold_in_registers_makes_the_steps_of_the_fold_on_the_map() {
    let mut rng = Lcg(0xF01D_4E65);
    let mut case_no = 0;
    let mut cover = RegCover::default();
    for word_size in [11usize, 12] {
        for word_type in [
            DiscWordType::Coding,
            DiscWordType::Optimal,
            DiscWordType::TwoTemplates,
        ] {
            for template_length in [16usize, 18, 21] {
                case_no += 1;
                let c = make_case(&mut rng, 20_000);
                let c = if case_no % 2 == 0 {
                    mask_case(&mut rng, c)
                } else {
                    c
                };
                let lookup = build_disc_mb_lookup(
                    &ctx_refs(&c),
                    &c.offsets,
                    word_size,
                    template_length,
                    word_type,
                    &c.masks,
                    c.concat.len(),
                );
                let forms = forms_of(&lookup);
                let wl = template_length;
                let seeds = all_seeds(&c, &lookup, &inputs(&c, &lookup, &forms, 40, 40));
                let mut windows = vec![1 + case_no % 6, wl + 1, 40];
                windows.dedup();
                for w in windows {
                    let at = format!("case {case_no} wl {wl} w {w}");
                    let x0 = inputs(&c, &lookup, &forms, 0, w);
                    let offsets = stage_offsets(37 + rng.below(1000) as i32, w, c.s_len, 3);
                    let got = reg_agrees_in_every_instance(&x0, &seeds, &offsets, &mut rng, &at);
                    if w == 40 {
                        assert!(got.hits > 20, "{at}: too few hits {got:?}");
                        if case_no == 1 {
                            // also cover below-cutoff extensions on real scan seeds
                            let high_cutoffs = vec![10_000; c.contexts.len()];
                            let x_hi = Inputs {
                                cutoff_scores: &high_cutoffs,
                                ..x0
                            };
                            cover = cover.sum(reg_agrees_in_every_instance(
                                &x_hi,
                                &seeds,
                                &offsets,
                                &mut rng,
                                &format!("{at} high-cutoff"),
                            ));
                        }
                    }
                    cover = cover.sum(got);
                }
            }
        }
    }
    assert!(
        cover.found > 0 && cover.rejected > 0 && cover.extended > cover.hits && cover.hits > 0,
        "{cover:?}"
    );
}

/// Seeds of `c` (template length `wl`) for the map's edge cases: positions that mostly stay or step
/// by 1 or 2, sometimes jump past every window, with bursts of up to 150 seeds at one position;
/// diagonals from `crowd` (diagonals that share their home bucket up to a few doublings), seeds
/// repeated on a recent diagonal (found, rejected, extended), and random ones. Each seed is a
/// pair inside a context of the query and the subject (any pair can be extended).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn edge_seeds(rng: &mut Lcg, c: &Case, wl: usize, crowd: &[i32], n: usize) -> Vec<Seed> {
    let context_of = |q: i64| -> Option<u32> {
        c.contexts
            .iter()
            .position(|ctx| {
                let start = i64::from(ctx.query_offset);
                q >= start && q + wl as i64 <= start + ctx.len as i64
            })
            .map(|k| k as u32)
    };
    let mut seeds: Vec<Seed> = Vec::new();
    let mut s = wl as i64;
    let mut burst = 0;
    while seeds.len() < n {
        if burst > 0 {
            burst -= 1;
        } else {
            match rng.below(100) {
                0..=1 => burst = 50 + rng.below(100),
                2..=4 => s += 50 + rng.below(150) as i64,
                5..=49 => {}
                _ => s += 1 + rng.below(2) as i64,
            }
        }
        if s + wl as i64 > c.s_len as i64 {
            break;
        }
        let diag = match rng.below(10) {
            // a crowded diagonal that a pair at s can have (q = s - diag in the query)
            0..=5 => {
                let lo = crowd.partition_point(|&d| i64::from(d) <= s - c.concat.len() as i64);
                let hi = crowd.partition_point(|&d| i64::from(d) <= s);
                if lo == hi {
                    continue;
                }
                i64::from(crowd[lo + rng.below((hi - lo) as u64) as usize])
            }
            // a diagonal of the last seeds
            6..=7 if !seeds.is_empty() => {
                let back = 1 + rng.below(seeds.len().min(30) as u64) as usize;
                let seed = &seeds[seeds.len() - back];
                i64::from(seed.s) - i64::from(seed.q)
            }
            _ => s - rng.below(c.concat.len() as u64) as i64,
        };
        let q = s - diag;
        if let Some(context) = context_of(q) {
            seeds.push(Seed {
                context,
                q: q as u32,
                s: s as u32,
                s_range: c.s_len as u32,
            });
        }
    }
    seeds
}

// EXPERIMENT (LOSAT_X_DCFOLDREG): as `the_fold_in_registers_makes_the_steps_of_the_fold_on_the_map`
// on seeds made for the map's edges (`edge_seeds`): diagonals crowding one home bucket (full
// buckets, the stash, puts over stash entries, the buckets doubled several times), bursts at one
// position, jumps past the window, and seeds repeated on a diagonal (found, rejected, two hits).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_fold_in_registers_makes_the_steps_of_the_fold_on_the_map_at_its_edges() {
    let mut rng = Lcg(0x0ED6_E5EE);
    let mut cover = RegCover::default();
    for (case_no, (word_size, word_type, template_length)) in [
        (11usize, DiscWordType::Coding, 16usize),
        (12, DiscWordType::Optimal, 18),
        (11, DiscWordType::TwoTemplates, 21),
    ]
    .into_iter()
    .enumerate()
    {
        let c = make_case_with(&mut rng, 30_000, 20_000);
        let lookup = build_disc_mb_lookup(
            &ctx_refs(&c),
            &c.offsets,
            word_size,
            template_length,
            word_type,
            &c.masks,
            c.concat.len(),
        );
        let forms = forms_of(&lookup);
        let wl = template_length;
        // the diagonals whose home is one bucket while the map has up to 2^7 buckets
        let top = |d: i32| (d as u32).wrapping_mul(0x9E37_79B1) >> 25;
        let want = top(5_000);
        let crowd: Vec<i32> = (-(c.concat.len() as i32)..c.s_len as i32)
            .filter(|&d| top(d) == want)
            .collect();
        assert!(crowd.len() > 400, "{}", crowd.len());
        let seeds = edge_seeds(&mut rng, &c, wl, &crowd, 60_000);
        for w in [1 + case_no % 6, wl - 1, wl + 1, 40, 64] {
            let at = format!("edge case {case_no} wl {wl} w {w}");
            let x0 = inputs(&c, &lookup, &forms, 0, w);
            // the hash offset starts at the window and returns to it at the clear (`latest_bound`),
            // so a stage's offset is never below it (below it, at w 64, the chunk has no bounds)
            let first = (37 + rng.below(1000) as i32).max(w as i32);
            let offsets = stage_offsets(first, w, c.s_len, 3);
            cover = cover.sum(reg_agrees_in_every_instance(
                &x0, &seeds, &offsets, &mut rng, &at,
            ));
        }
    }
    assert!(
        cover.found > 0
            && cover.rejected > 0
            && cover.extended > 0
            && cover.stash_found > 0
            && cover.stash_put > 0
            && cover.most_stash > 0
            && cover.most_buckets >= 1 << 9,
        "{cover:?}"
    );
}

// EXPERIMENT (LOSAT_X_DCEMITHOIST): `DcPairs` and `scan_strip_disc_csr`'s `HOIST` instance against
// `seed_of_pair` and the instance without `HOIST`, the code `LOSAT_X_DCEMITHOIST=0` runs.

/// Contexts of `lengths` (0: a context without letters) as the batch lays them out: each context
/// with letters right after the one before and its sentinel, one without letters at the place the
/// next one starts. `gap` letters more after each sentinel leave no constant-time form (the runs
/// answer); a `back` above 1 makes contexts overlap (the table answers).
fn hoist_contexts(lengths: &[usize], gap: usize, back: usize) -> Vec<QueryContext> {
    let mut at = 0usize;
    lengths
        .iter()
        .enumerate()
        .map(|(k, &len)| {
            let query_offset = at as i32;
            if len > 0 {
                at = (at + len + 1 + gap).saturating_sub(back);
            }
            QueryContext {
                query_idx: (k / 2) as u32,
                frame: if k % 2 == 0 { 1 } else { -1 },
                query_offset,
                len,
                masks: Vec::new(),
            }
        })
        .collect()
}

/// The form of a context index that `DcPairs` asks (`QueryContextIndex::x_view`, and for offsets it
/// leaves, `x_context_for_offset_ref`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum HoistForm {
    Div,
    Blocks,
    BlocksRanks,
    /// `Blocks` with crowded blocks (their offsets go to the runs).
    Crowded,
    /// No constant-time form: the runs.
    Runs,
    /// No constant-time form, overlapping contexts: `direct_map`.
    Table,
}

fn hoist_form(index: &QueryContextIndex) -> HoistForm {
    match (&index.x_o1, index.x_view()) {
        (Some(o1), _) if o1.needs_runs() => HoistForm::Crowded,
        (_, Some(XView::Div(_))) => HoistForm::Div,
        (_, Some(XView::Blocks(_))) => HoistForm::Blocks,
        (_, Some(XView::BlocksRanks(_))) => HoistForm::BlocksRanks,
        (_, None) if !index.direct_map.is_empty() => HoistForm::Table,
        (_, None) => HoistForm::Runs,
    }
}

/// Random contexts of the form `form` (asserted by the caller) over about `cover` offsets:
/// lengths, contexts without letters, gaps and overlaps.
fn hoist_layout(rng: &mut Lcg, form: HoistForm, cover: usize) -> Vec<QueryContext> {
    // lengths from `len` until they cover `cover` offsets (at least two)
    let lengths = |rng: &mut Lcg, len: &mut dyn FnMut(&mut Lcg, usize) -> usize| {
        let mut lengths = Vec::new();
        let mut sum = 0usize;
        while sum < cover || lengths.len() < 2 {
            let l = len(rng, lengths.len());
            sum += l + 1;
            lengths.push(l);
        }
        lengths
    };
    let long = |rng: &mut Lcg| 20 + rng.below(3000) as usize;
    // lengths within 1024..1536 keep 2^shift = 1024 with zero crowded blocks (no shift step can
    // stay below 1 / CROWDED_SHARE)
    let even = |rng: &mut Lcg| 1024 + rng.below(512) as usize;
    match form {
        HoistForm::Div => {
            let len = [1, 17, 18, 21, 64, 1000][rng.below(6) as usize];
            hoist_contexts(&vec![len; (cover / (len + 1)).clamp(1, 4000)], 0, 0)
        }
        HoistForm::Blocks => {
            let mut lengths = lengths(rng, &mut |rng, _| even(rng));
            lengths[0] = lengths[1] + 1;
            hoist_contexts(&lengths, 0, 0)
        }
        HoistForm::BlocksRanks => {
            // a plus-strand search (every minus context without letters), or zeros anywhere
            let alternate = rng.below(2) == 0;
            let mut lengths = lengths(rng, &mut |rng, k| {
                if (alternate && k % 2 == 1) || (!alternate && rng.below(3) == 0) {
                    0
                } else {
                    even(rng)
                }
            });
            lengths[0] = even(rng);
            lengths[1] = 0;
            hoist_contexts(&lengths, 0, 0)
        }
        HoistForm::Crowded => {
            // long contexts with three consecutive 1-letter contexts: at shift 3 (block 8) at
            // least one adjacent pair of sentinels shares a block while crowded_blocks <= blocks/256
            let mut lengths = Vec::new();
            let mut sum = 0usize;
            let need = cover.max(20_000);
            while sum < need || lengths.len() < 10 {
                let k = lengths.len();
                let l = if (3..=5).contains(&k) { 1 } else { even(rng) };
                sum += l + 1;
                lengths.push(l);
            }
            hoist_contexts(&lengths, 0, 0)
        }
        HoistForm::Runs => {
            let lengths = lengths(rng, &mut |rng, _| long(rng));
            hoist_contexts(&lengths, 1 + rng.below(5) as usize, 0)
        }
        HoistForm::Table => {
            let lengths = lengths(rng, &mut |rng, _| long(rng) + 20);
            hoist_contexts(&lengths, 0, 2 + rng.below(15) as usize)
        }
    }
}

/// The offsets and lengths of `contexts`, for messages.
fn hoist_spans(contexts: &[QueryContext]) -> Vec<(i32, usize)> {
    contexts.iter().map(|c| (c.query_offset, c.len)).collect()
}

/// A sink with `push` alone (`SeedOut::x_push_run`'s default, which `part::PartSink` and
/// `part::PieceSink` take).
struct PushOnly(Vec<Seed>);

impl SeedOut for PushOnly {
    fn push(&mut self, seed: Seed) {
        self.0.push(seed);
    }
}

/// A sink of at most `limit` seeds, closed as a seed piece is (`part::PieceSink`'s `room` and
/// `close`).
struct Capped {
    seeds: Vec<Seed>,
    limit: usize,
}

impl SeedOut for Capped {
    fn push(&mut self, seed: Seed) {
        assert!(self.seeds.len() < self.limit, "a push past the room");
        self.seeds.push(seed);
    }

    fn room(&self) -> usize {
        self.limit - self.seeds.len()
    }

    fn close(&mut self) {
        self.limit = self.seeds.len();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1106-1144
// `DcPairs::run` against `seed_of_pair` called for each pair of a run, as `scan_strip_disc_csr`
// calls it: random runs of query offsets (random, and at the edges: the first and the last offset
// of a context where a word fits, the first where it does not, sentinels, offsets at and past the
// contexts' end, which the binary search answers) at random subject offsets (and at the subject's
// end: the last word that fits and the first that does not), for every form of the context index
// (`Div`, `Blocks`, `Blocks` with ranks, crowded `Blocks`, the runs, the table; both strands' and
// one strand's bounds of the binary search) and the templates of lengths 16, 18 and 21. Every
// field of every seed, in order, after seeds already in the list (a list with no spare room, and
// one with room), and through `push` alone (`x_push_run`'s default).
#[test]
fn the_hoisted_emit_gives_the_seeds_of_seed_of_pair() {
    let mut rng = Lcg(0xE317_0001);
    let c = make_case_with(&mut rng, 3000, 1100);
    let mut forms_seen = std::collections::BTreeSet::new();
    let (mut kept, mut ctx_dropped, mut s_dropped) = (0usize, 0usize, 0usize);
    let forms = [
        HoistForm::Div,
        HoistForm::Blocks,
        HoistForm::BlocksRanks,
        HoistForm::Crowded,
        HoistForm::Runs,
        HoistForm::Table,
    ];
    for (word_type, template_length) in [
        (DiscWordType::Coding, 16usize),
        (DiscWordType::Optimal, 18),
        (DiscWordType::Coding, 21),
    ] {
        let lookup = build_disc_mb_lookup(
            &ctx_refs(&c),
            &c.offsets,
            11,
            template_length,
            word_type,
            &c.masks,
            c.concat.len(),
        );
        let scan_forms = forms_of(&lookup);
        let word = lookup.word_length();
        assert_eq!(word, template_length);
        for round in 0..60 {
            let form = forms[round % forms.len()];
            let cover = 2000 + rng.below(60_000) as usize;
            let contexts = hoist_layout(&mut rng, form, cover);
            let one_strand = rng.below(2) == 0;
            let index = QueryContextIndex::new(&contexts, one_strand);
            assert_eq!(hoist_form(&index), form, "{:?}", hoist_spans(&contexts));
            forms_seen.insert(form);
            let x = Inputs {
                query_context_index: &index,
                query_contexts: &contexts,
                ..inputs(&c, &lookup, &scan_forms, 0, 40)
            };
            let view = index.x_view();
            let total = contexts
                .iter()
                .map(|ctx| ctx.query_offset as usize + ctx.len)
                .max()
                .unwrap();
            // the edges of every context (1-based offsets of the run)
            let mut edges = vec![1usize, total, total + 1, total + 40];
            for ctx in &contexts {
                let start = ctx.query_offset as usize;
                let end = start + ctx.len;
                edges.extend([start + 1, end + 1, end + 2]);
                if ctx.len >= word {
                    edges.extend([end - word + 1, end - word + 2]);
                }
            }
            for _ in 0..400 {
                let len = match rng.below(10) {
                    0 => 0,
                    1 => 40 + rng.below(60) as usize,
                    _ => 1 + rng.below(12) as usize,
                };
                let run: Vec<u32> = (0..len)
                    .map(|_| {
                        if rng.below(3) == 0 {
                            edges[rng.below(edges.len() as u64) as usize] as u32
                        } else {
                            1 + rng.below(total as u64 + 60) as u32
                        }
                    })
                    .collect();
                let kmer_start = match rng.below(8) {
                    0 => c.s_len - word,
                    1 => c.s_len - word + 1,
                    2 => c.s_len,
                    3 => 0,
                    _ => rng.below((c.s_len + 1) as u64) as usize,
                };
                let s_range = kmer_start + word + rng.below(5000) as usize;
                let prefix: Vec<Seed> = (0..rng.below(4))
                    .map(|k| Seed {
                        context: k as u32,
                        q: 7,
                        s: 9,
                        s_range: 11,
                    })
                    .collect();
                let mut want = prefix.clone();
                for &q_off_1 in &run {
                    seed_of_pair(&x, q_off_1 as usize - 1, kmer_start, s_range, &mut want);
                }
                let pairs = DcPairs::new(&x, s_range);
                let at = format!("{form:?} word {word} run {run:?} at {kmer_start}");
                // a list without spare room, then one with room
                let mut tight = prefix.clone();
                tight.shrink_to_fit();
                pairs.run(view, &run, kmer_start, &mut tight);
                assert_eq!(tight, want, "{at}");
                let mut roomy = Vec::with_capacity(prefix.len() + run.len() + 5);
                roomy.extend_from_slice(&prefix);
                pairs.run(view, &run, kmer_start, &mut roomy);
                assert_eq!(roomy, want, "{at}");
                let mut push_only = PushOnly(prefix.clone());
                pairs.run(view, &run, kmer_start, &mut push_only);
                assert_eq!(push_only.0, want, "{at}");
                let added = want.len() - prefix.len();
                kept += added;
                if kmer_start + word > c.s_len {
                    s_dropped += run.len();
                } else {
                    ctx_dropped += run.len() - added;
                }
            }
        }
    }
    assert_eq!(forms_seen.len(), forms.len());
    assert!(
        kept > 10_000 && ctx_dropped > 1000 && s_dropped > 1000,
        "{kept} {ctx_dropped} {s_dropped}"
    );
}

/// The pieces of the scan `scan` of a strip from `start` into sinks of `limit` seeds, as
/// LOSAT_X_SEEDPIECE cuts a strip: each piece's seeds and where it stopped.
fn hoist_pieces(
    limit: usize,
    start: usize,
    mut scan: impl FnMut(Resume, &mut Capped) -> Option<Resume>,
) -> Vec<(Vec<Seed>, Option<Resume>)> {
    let mut from = Resume::at(start);
    let mut pieces = Vec::new();
    loop {
        let mut sink = Capped {
            seeds: Vec::new(),
            limit,
        };
        let stop = scan(from, &mut sink);
        pieces.push((sink.seeds, stop));
        match stop {
            Some(at) => from = at,
            None => return pieces,
        }
        assert!(pieces.len() < 1 << 22, "no progress");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1406-1418
// `scan_strip_disc_csr` with `HOIST` against it without on CSR tables (a query of 2 x 560,000
// letters: above LOSAT_X_DCCSR's size rule for 4^11 cells; repeats give runs of 20 and more) of
// the templates of lengths 18 and 21, the strips of 1, 5 and 33 per subject of the strips stage
// (`plan_strips`), the subject whole and masked (scanned in its ranges), the contexts of the
// search and contexts of every form of the index laid over the same query: the same seeds in the
// same order, field by field. Then in pieces of 1, 2, 3 and 7 seeds (LOSAT_X_SEEDPIECE's cut:
// `Cut::fit` and its `Err` path `Cut::run_cold`, which stops inside a run): the same seeds in each
// piece and the same stop.
#[test]
fn the_hoisted_csr_scan_gives_the_seeds_and_stops_of_the_csr_scan() {
    let mut rng = Lcg(0xE317_0002);
    let q_len = 560_000;
    let c = make_case_with(&mut rng, 1000, q_len);
    // a subject of diverged copies of the query's repeat (`make_case_with`: at offset 0, twenty
    // times in its first 2,000 letters) and of other pieces of both strands
    let plus: Vec<u8> = c.ctx_seqs[0].iter().map(|&b| b"ACGT"[b as usize]).collect();
    let minus = reverse_complement(&plus);
    let mut subject = Vec::new();
    while subject.len() < 24_000 {
        let src = if rng.below(2) == 0 { &plus } else { &minus };
        let len = 30 + rng.below(800) as usize;
        let at = match rng.below(3) {
            0 => rng.below(2000) as usize,
            1 => q_len - 2000 + rng.below(1000) as usize,
            _ => rng.below((q_len - len) as u64) as usize,
        };
        for &b in &src[at..(at + len).min(q_len)] {
            if rng.below(100) < 6 {
                subject.push(b"ACGT"[rng.below(4) as usize]);
            } else {
                subject.push(b);
            }
        }
    }
    let s_len = subject.len();
    let packed = encode_subject_ncbi2na_packed(&subject);
    let whole = vec![(0i32, s_len as i32)];
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
    // the search's contexts (`Div`: both strands of one query), and others over the same offsets
    let mut layouts = vec![c.contexts.clone()];
    for form in [
        HoistForm::Blocks,
        HoistForm::BlocksRanks,
        HoistForm::Crowded,
        HoistForm::Runs,
    ] {
        let contexts = hoist_layout(&mut rng, form, c.concat.len());
        let index = QueryContextIndex::new(&contexts, false);
        assert_eq!(hoist_form(&index), form, "{:?}", hoist_spans(&contexts));
        layouts.push(contexts);
    }
    let (mut seeds_n, mut stops_in_runs) = (0usize, 0usize);
    for (word_type, template_length) in
        [(DiscWordType::Optimal, 18usize), (DiscWordType::Coding, 21)]
    {
        let lookup = build_disc_mb_lookup(
            &ctx_refs(&c),
            &c.offsets,
            11,
            template_length,
            word_type,
            &c.masks,
            c.concat.len(),
        );
        let forms = forms_of(&lookup);
        let disc = lookup.disc().expect("a discontiguous table");
        let Some(csr) = disc.x_csr() else {
            assert_eq!(
                std::env::var("LOSAT_X_DCCSR").as_deref(),
                Ok("0"),
                "no CSR table"
            );
            return;
        };
        assert!(
            (1..1u32 << 22).any(|cell| csr.run(cell as usize).len() >= 20),
            "no long run"
        );
        for (masked, seq_ranges) in [(false, &whole), (true, &ranges)] {
            for (li, contexts) in layouts.iter().enumerate() {
                let index = QueryContextIndex::new(contexts, false);
                let x = Inputs {
                    search_seq_packed: &packed,
                    s_len,
                    subject_seq_ranges: seq_ranges,
                    subject_masked: masked,
                    query_context_index: &index,
                    query_contexts: contexts,
                    ..inputs(&c, &lookup, &forms, 0, 40)
                };
                for strips in [1usize, 5, 33] {
                    let jobs = plan_strips(
                        seq_ranges,
                        s_len,
                        template_length,
                        template_length,
                        1,
                        masked,
                        strips,
                    );
                    assert!(!jobs.is_empty());
                    for (j, job) in jobs.iter().enumerate() {
                        let at = format!(
                            "template {template_length} masked {masked} layout {li} strips {strips} job {j} {job:?}"
                        );
                        let from = Resume::at(job.start);
                        let mut want: Vec<Seed> = Vec::new();
                        let r = scan_strip_disc_csr::<_, false>(
                            &x, disc, csr, job, from, &mut want, false,
                        );
                        assert_eq!(r, None, "{at}");
                        let mut got: Vec<Seed> = Vec::new();
                        let r = scan_strip_disc_csr::<_, false>(
                            &x, disc, csr, job, from, &mut got, true,
                        );
                        assert_eq!(r, None, "{at}");
                        assert_eq!(got, want, "{at}");
                        seeds_n += got.len();
                        // pieces (one layout, the subject whole and masked, one strip count)
                        if li > 1 || strips != 5 {
                            continue;
                        }
                        for limit in [1usize, 2, 3, 7] {
                            let want = hoist_pieces(limit, job.start, |from, out| {
                                scan_strip_disc_csr::<_, true>(&x, disc, csr, job, from, out, false)
                            });
                            let got = hoist_pieces(limit, job.start, |from, out| {
                                scan_strip_disc_csr::<_, true>(&x, disc, csr, job, from, out, true)
                            });
                            assert_eq!(got, want, "{at} pieces of {limit}");
                            stops_in_runs += got
                                .iter()
                                .filter(|(_, stop)| stop.is_some_and(|r| r.skip > 0))
                                .count();
                        }
                    }
                }
            }
        }
    }
    assert!(
        seeds_n > 100_000 && stops_in_runs > 1000,
        "{seeds_n} seeds, {stops_in_runs} stops inside runs"
    );
}
