use super::super::super::{
    build_blastna_matrix, build_nucl_score_table, diag_hash_insert_window,
    scan_subject_kmers_with_ranges, select_mb_scan_kind, DIAGHASH_NUM_BUCKETS,
};
use super::*;
use crate::algorithm::blastn::lookup::{build_pv_direct_lookup, reverse_complement};
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
        self.next() % n
    }
}

fn random_dna(rng: &mut Lcg, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

type HashState = (Vec<u32>, Vec<(i32, i32, bool, i32, u32)>, u32, u32, i32);

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684
// The whole stage on random sequences (the subject made of mutated copies of query pieces
// on both strands, query masks on one context, so `s_TypeOfWord` checks words there) for
// the direct table at word sizes 6, 7 and 8: the pairs of the reference wrapper
// `scan_subject_kmers_with_ranges` folded in one go, against `seed_stage` on 1..6 threads
// with the owner's fold and (2 threads or more) the partitioned fold in both forms of
// LOSAT_X_BNSEEDLISTS's stage, for the strip target
// of `run` and for two small targets whose strips start at offsets that are not a
// multiple of 4 (`plan_strips`: per = max(ceil(grid points / target), 512)). Hit lists and
// diagonal hashes must be the same.
#[test]
fn seed_stage_on_strips_equals_one_scan() {
    let mut rng = Lcg(0x0B13_5EED);
    // each word size on an unmasked and a masked subject
    for (word_length, masked) in [6usize, 7, 8]
        .into_iter()
        .flat_map(|w| [(w, false), (w, true)])
    {
        let q_len = 3000 + rng.below(3000) as usize;
        let mut query = random_dna(&mut rng, q_len);
        // a repeat inside the query (long lookup lists)
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
        let ctx_refs: Vec<&[u8]> = ctx_seqs.iter().map(Vec::as_slice).collect();
        // The concatenated query (the contexts at `offsets`, a sentinel between them).
        let concat: Vec<u8> = [ctx_seqs[0].as_slice(), &[15], ctx_seqs[1].as_slice()].concat();
        let lookup = build_pv_direct_lookup(
            &ctx_refs,
            &offsets,
            word_length,
            word_length,
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
        let score_matrix = build_blastna_matrix(1, -3);
        let nucl_score_table = build_nucl_score_table(1, -3);
        // a masked subject's scan ranges (NCBI's `seq_ranges`, right end exclusive) cover most
        // of it with gaps of 1 to 120 bases, some shorter than a word (the scan skips them)
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
        let x = NaInputs {
            pv_direct_lookup: Some(&lookup),
            na_lookup: None,
            use_direct_lookup: true,
            word_length,
            search_seq_packed: &packed,
            s_len,
            subject_seq_ranges: &seq_ranges,
            subject_masked: masked,
            scan_step: 1,
            query_context_index: &index,
            query_contexts: &contexts,
            encoded_query_concat_blastna: &concat,
            cutoff_scores: &[16, 16],
            x_dropoff_scores: &[20, 20],
            reduced_cutoff_scores: &[16, 16],
            score_matrix: &score_matrix,
            nucl_score_table: &nucl_score_table,
            exact_ungapped: word_length < 11,
            diag_offset,
            diag_hash_window: diag_hash_insert_window(0, 0, word_length),
            use_array_indexing: false,
            window_size: 0,
            scan_range: 0,
            diagnostics: false,
            one_subject_owner: true,
            timing: None,
        };
        let new_hash = || {
            let mut h = DiagHashTable::new(0);
            h.offset = diag_offset as i32;
            h
        };
        // the reference wrapper visits the pairs in one scan
        let mut seeds = Vec::new();
        scan_subject_kmers_with_ranges(
            &packed,
            s_len,
            word_length,
            word_length,
            1,
            &seq_ranges,
            masked,
            select_mb_scan_kind(word_length, 1, masked),
            |kmer_start, s_range, kmer| {
                for &q_off_1 in x.hits(kmer) {
                    if q_off_1 == 0 {
                        continue;
                    }
                    let q_off0 = q_off_1 as usize - 1;
                    seeds.push(Seed {
                        context: index.context_for_offset(q_off0) as u32,
                        q: q_off0 as u32,
                        s: kmer_start as u32,
                        s_range: s_range as u32,
                    });
                }
            },
        );
        let mut serial_hash = new_hash();
        let mut serial_hits = Vec::new();
        let mut counts = FoldCounts::default();
        let mut pause = 0u64;
        fold_na::<false>(
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
        let below_cutoff = seeds.len() as u64 - counts.skipped - serial_hits.len() as u64;
        assert!(
            serial_hits.len() > 5 && counts.skipped > 0 && below_cutoff > 0,
            "word {word_length} masked {masked}: the test sequences give too few cases ({} hits, {} skipped, {} other)",
            serial_hits.len(),
            counts.skipped,
            below_cutoff
        );
        let key = |h: &UngappedHit| {
            (
                h.context_idx,
                h.query_idx,
                h.query_frame,
                h.query_context_offset,
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
        let grid = s_len - word_length + 1;
        // two small targets whose strips (but the first) start off a multiple of 4
        let odd_targets: Vec<usize> = (2..64)
            .filter(|&t| grid.div_ceil(t) > 512 && grid.div_ceil(t) % 4 != 0)
            .take(2)
            .collect();
        assert_eq!(
            odd_targets.len(),
            2,
            "word {word_length} masked {masked}: grid {grid}"
        );
        for threads in 1..=6usize {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            for target in [
                STRIPS_PER_THREAD * threads.max(4),
                odd_targets[0],
                odd_targets[1],
            ] {
                let per = grid.div_ceil(target).max(512);
                // unpartitioned, and partitioned in both forms of LOSAT_X_BNSEEDLISTS's stage
                // on one scratch (its strip buffers passing from one form to the other)
                let mut scratch = PartScratch::default();
                // (and the lists' stage in pieces of 1 seed, a few and hundreds,
                // EXPERIMENT LOSAT_X_SEEDPIECE: the kernels stop inside words' runs)
                for (partitioned, lists, pieces) in [
                    (false, false, 0usize),
                    (true, false, 0),
                    (true, true, 0),
                    (true, true, 1),
                    (true, true, 2 + (threads + target) % 7),
                    (true, true, 300 + 17 * threads),
                    (true, false, 0),
                ] {
                    if partitioned && threads < 2 {
                        continue;
                    }
                    scratch.x_pieces = Some(pieces);
                    let at =
                        format!("word {word_length} masked {masked} threads {threads} target {target} partitioned {partitioned} lists {lists} pieces {pieces}");
                    let (state, keys, stats_seeds, strips) = pool.install(|| {
                        let mut hash = new_hash();
                        let mut hits = Vec::new();
                        let fold_part = partitioned.then(|| part::FoldPart {
                            scratch: &mut scratch,
                            parts: threads.min(DIAGHASH_NUM_BUCKETS),
                        });
                        let stats = seed_stage_of(
                            &x, threads, target, &mut hash, &mut hits, None, false, fold_part,
                            lists,
                        );
                        (
                            hash_state(&hash),
                            hits.iter().map(key).collect::<Vec<_>>(),
                            stats.counts.seeds,
                            stats.strips,
                        )
                    });
                    assert_eq!(stats_seeds, seeds.len() as u64, "{at}");
                    // (a masked subject's plan has a strip or more per scan range)
                    if !masked {
                        assert_eq!(strips, grid.div_ceil(per), "{at}");
                    }
                    assert_eq!(keys, serial_keys, "{at}");
                    assert_eq!(state, hash_state(&serial_hash), "{at}");
                }
            }
        }
    }
}
