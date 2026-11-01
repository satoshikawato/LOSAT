//! No NCBI counterpart (`LOSAT_X_TBNLINKCTX`, strict set; `=0` links every query context in one
//! pass, as before): the HSPs of a subject are linked one query context at a time.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1048-1060
//! ```c
//!       /* Find the first HSP in the context of the head */
//!       for (index = 0; index < hsp_list->hspcnt &&
//!                       hsp_list->hsp_array[index]->context < head_hsp->context;
//! ```
//! NCBI links a subject's list in one pass over every query context: the heads come in sum-score
//! order whatever their context, and a head's candidates and the set it joins are HSPs of its own
//! context. Nothing a head does reads or changes another context's HSPs, so each context's HSPs
//! linked on their own, in that order restricted to the context, give the same chains, numbers and
//! e-values; the values go back to the order of the list before the final sort by score. A subject
//! whose units ran on the pool (the single subject of a search, or one longer than a chunk's
//! budget) links its contexts there side by side, the ones with the most HSPs first; any other
//! subject links them one after another on its thread. A list of one context is linked as before.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

use super::super::search_gapped::GappedHsp;
use super::LinkedHsp;
use crate::utils::threading::SearchPool;

/// `LOSAT_X_TBNLINKCTX` (strict set; `=0` off).
pub(super) fn x_link_ctx_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNLINKCTX"))
}

// No NCBI counterpart: LOSAT_X_STATS counters; they do not change any value NCBI computes.
/// Lists linked one context at a time, their contexts, and those of them linked on a pool.
static STATS: [AtomicU64; 3] = [const { AtomicU64::new(0) }; 3];

fn x_stats_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_STATS").is_some())
}

/// The `[X_STATS] TBNLINKCTX` line (printed by `main` at exit) with `LOSAT_X_STATS`, when a list
/// of two or more contexts was linked.
pub fn print_summary() {
    let [lists, contexts, pooled] = STATS.each_ref().map(|v| v.load(Ordering::Relaxed));
    if lists > 0 && x_stats_on() {
        eprintln!("[X_STATS] TBNLINKCTX lists={lists} contexts={contexts} pooled={pooled}");
    }
}

/// The linked values of `input` in its order (`source_index` = the position in `input`), each
/// query context linked on its own by `link_values` (which links a list of one or more contexts
/// and gives its values in the order of that list): on `pool` when it is enabled, else one after
/// another.
pub(super) fn x_link_by_context(
    pool: Option<&SearchPool<'_>>,
    input: &[(usize, GappedHsp)],
    link_values: impl Fn(&[(usize, GappedHsp)]) -> Vec<LinkedHsp> + Sync,
) -> Vec<LinkedHsp> {
    if input.windows(2).all(|pair| pair[0].0 == pair[1].0) {
        return link_values(input);
    }
    let contexts = input
        .iter()
        .map(|&(context, _)| context + 1)
        .max()
        .unwrap_or(0);
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); contexts];
    for (index, &(context, _)) in input.iter().enumerate() {
        groups[context].push(index);
    }
    groups.retain(|group| !group.is_empty());
    let link_group = |group: &Vec<usize>| {
        let part: Vec<(usize, GappedHsp)> = group.iter().map(|&index| input[index]).collect();
        let mut values = link_values(&part);
        for value in &mut values {
            value.source_index = group[value.source_index];
        }
        values
    };
    let mut linked: Vec<Vec<LinkedHsp>> = Vec::new();
    let stats = x_stats_on();
    if stats {
        STATS[0].fetch_add(1, Ordering::Relaxed);
        STATS[1].fetch_add(groups.len() as u64, Ordering::Relaxed);
    }
    #[cfg(feature = "parallel")]
    if let Some(pool) = pool.filter(|pool| pool.enabled()) {
        use rayon::prelude::*;
        let mut order: Vec<&Vec<usize>> = groups.iter().collect();
        order.sort_by_key(|group| std::cmp::Reverse(group.len()));
        if stats {
            STATS[2].fetch_add(1, Ordering::Relaxed);
        }
        linked = pool.install(|| {
            order
                .par_iter()
                .with_max_len(1)
                .map(|group| link_group(group))
                .collect()
        });
    }
    #[cfg(not(feature = "parallel"))]
    let _ = pool;
    if linked.is_empty() {
        linked = groups.iter().map(link_group).collect();
    }
    let mut slots: Vec<Option<LinkedHsp>> = vec![None; input.len()];
    for value in linked.into_iter().flatten() {
        slots[value.source_index] = Some(value);
    }
    slots
        .into_iter()
        .map(|value| value.expect("each HSP lies in one context"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::super::{x_link_list, x_link_values, LinkedHspList};
    use super::*;
    use crate::algorithm::tblastn::stage_d_stats::{LocalContextLength, LocalLinkParameters};
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    use crate::stats::spouge::lookup_protein_gumbel_params;
    use crate::stats::tables::lookup_protein_params_gapped;

    /// xorshift64*: the test's own random numbers.
    struct Random(u64);
    impl Random {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }
        fn range(&mut self, low: i32, high: i32) -> i32 {
            low + (self.next() % (high - low + 1) as u64) as i32
        }
    }

    /// A subject list of `count` HSPs on `contexts` query contexts (the contexts drawn with
    /// unequal weights, some never), both strands: chains of pieces that keep or just break the
    /// linking gaps, equal scores, repeated offsets.
    fn hsps(
        random: &mut Random,
        count: usize,
        contexts: usize,
        subject_aa: i32,
        link: &LocalLinkParameters,
    ) -> Vec<(usize, GappedHsp)> {
        let mut out = Vec::with_capacity(count);
        while out.len() < count {
            let a = random.range(0, contexts as i32 - 1);
            let context = a.min(random.range(0, contexts as i32 - 1)) as usize;
            let sign = if random.next() % 2 == 0 { 1 } else { -1 };
            let mut q = random.range(0, 300);
            let mut s = random.range(0, subject_aa - 2_000);
            for _ in 0..random.range(1, 6) {
                if out.len() == count {
                    break;
                }
                let length = random.range(8, 60);
                let hsp = GappedHsp {
                    frame: sign * random.range(1, 3) as i8,
                    score: 22 + 4 * random.range(0, 30),
                    q_start: q,
                    q_end: q + length.min(400 - q).max(1),
                    q_gapped_start: q,
                    s_start: s,
                    s_end: s + length + random.range(-3, 3),
                    s_gapped_start: s,
                };
                out.push((context, hsp));
                q = (hsp.q_end + random.range(-link.overlap_size - 2, link.gap_size + 2))
                    .clamp(0, 390);
                s = hsp.s_end + random.range(-link.overlap_size - 2, link.longest_intron + 5);
                if random.next() % 8 == 0 {
                    s = random.range(0, subject_aa - 2_000);
                }
            }
        }
        out
    }

    fn assert_same(actual: &LinkedHspList, expected: &LinkedHspList, label: &str) {
        assert_eq!(
            actual.best_evalue.to_bits(),
            expected.best_evalue.to_bits(),
            "{label}"
        );
        assert_eq!(actual.hsps.len(), expected.hsps.len(), "{label}");
        for (a, e) in actual.hsps.iter().zip(&expected.hsps) {
            assert_eq!(
                (a.context, a.hsp, a.num, a.evalue.to_bits(), a.source_index),
                (e.context, e.hsp, e.num, e.evalue.to_bits(), e.source_index),
                "{label}"
            );
        }
    }

    // No NCBI counterpart (LOSAT_X_TBNLINKCTX): each context linked on its own, on one thread
    // or on a pool of 4, gives the list of the one pass over every context (`=0`): order,
    // context, HSP, number, e-value bits, input position and best e-value; on lists of 1 to
    // 3,000 HSPs over 1 to 35 contexts, with sets of two or more in most lists, and on lists whose
    // HSPs are each in the first and the last context (a query given twice), where the final
    // sort meets ties in every field it reads and keeps the order of the list.
    #[test]
    fn contexts_linked_apart_give_the_list_of_one_pass() {
        let link = LocalLinkParameters {
            gap_decay_rate: 0.1,
            gap_size: 40,
            overlap_size: 9,
            longest_intron: 122,
            cutoff_small_gap: 0,
        };
        let subject_aa = 60_000;
        let subject_nt = subject_aa * 3 + 2;
        let spec = ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        };
        let gumbel = lookup_protein_gumbel_params(&spec, i64::from(subject_aa)).unwrap();
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let mut random = Random(0x9E37_79B9_7F4A_7C15);
        let mut linked_sets = 0;
        crate::utils::threading::with_search_pool(4, "tblastn", |pool| {
            for contexts in [1usize, 2, 3, 7, 35] {
                let lengths = vec![
                    LocalContextLength {
                        length_adjustment: 20,
                        eff_searchsp: 380 * (i64::from(subject_aa) - 20),
                        subject_stat_length: i64::from(subject_aa),
                    };
                    contexts
                ];
                let query_lengths = vec![400; contexts];
                let params = vec![gapped; contexts];
                let link_values = |input: &[(usize, GappedHsp)]| {
                    x_link_values(
                        input,
                        &query_lengths,
                        &lengths,
                        subject_nt,
                        &params,
                        Some(&gumbel),
                        &link,
                    )
                };
                for count in [1usize, 2, 5, 60, 400, 3_000] {
                    let input = hsps(&mut random, count, contexts, subject_aa, &link);
                    let expected = x_link_list(link_values(&input));
                    linked_sets += expected.hsps.iter().filter(|h| h.num > 1).count();
                    let label = format!("{contexts} contexts, {count} HSPs");
                    let serial = x_link_list(x_link_by_context(None, &input, link_values));
                    assert_same(&serial, &expected, &format!("{label}, one thread"));
                    let pooled = x_link_list(x_link_by_context(Some(pool), &input, link_values));
                    assert_same(&pooled, &expected, &format!("{label}, pool of 4"));
                    if contexts == 1 {
                        continue;
                    }
                    let twice: Vec<(usize, GappedHsp)> =
                        hsps(&mut random, count, 1, subject_aa, &link)
                            .into_iter()
                            .flat_map(|(_, hsp)| [(contexts - 1, hsp), (0, hsp)])
                            .collect();
                    let expected = x_link_list(link_values(&twice));
                    let label = format!("{label}, each HSP in contexts 0 and {}", contexts - 1);
                    let serial = x_link_list(x_link_by_context(None, &twice, link_values));
                    assert_same(&serial, &expected, &format!("{label}, one thread"));
                    let pooled = x_link_list(x_link_by_context(Some(pool), &twice, link_values));
                    assert_same(&pooled, &expected, &format!("{label}, pool of 4"));
                }
            }
            Ok(())
        })
        .unwrap();
        assert!(linked_sets > 1_000, "sets of two or more: {linked_sets}");
    }
}
