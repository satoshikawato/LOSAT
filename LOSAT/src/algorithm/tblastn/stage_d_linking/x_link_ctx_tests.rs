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
        let sign = if random.next().is_multiple_of(2) {
            1
        } else {
            -1
        };
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
            q = (hsp.q_end + random.range(-link.overlap_size - 2, link.gap_size + 2)).clamp(0, 390);
            s = hsp.s_end + random.range(-link.overlap_size - 2, link.longest_intron + 5);
            if random.next().is_multiple_of(8) {
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
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
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
                let twice: Vec<(usize, GappedHsp)> = hsps(&mut random, count, 1, subject_aa, &link)
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
