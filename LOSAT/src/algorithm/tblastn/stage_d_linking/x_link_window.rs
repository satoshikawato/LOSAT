//! EXPERIMENT (LOSAT_X_TBNLINKWIN): the candidates of a head in tblastn's uneven-gap linking
//! chosen by subject position.
//!
//! NCBI's loop (c++/src/algo/blast/core/link_hsps.c:1602-1757, `s_BlastUnevenGapLinkHSPs`) tests
//! every HSP whose query offsets fall in the window of the head set
//! (`s_HSPOffsetEndBinarySearch` .. `s_HSPOffsetBinarySearch`) with
//! `s_LinkedHSPSetsAdmissible`, which merges the two sets and checks each pair of neighbours. A
//! short query puts nearly every HSP of the subject list in that window, so a genome searched
//! with one protein (thousands of HSPs per list, a handful of links) costs a number of
//! admissibility tests that grows with the square of the list (perf-si-20261012 session 21: a
//! 411-residue protein x the human genome, 85,908 HSPs, 159 million candidates, 13 admissible).
//!
//! Every admissible merged set is ordered by query offset (the merge) and, by the checks, has
//! each subject start before the next one and each next start at most `longest_intron` after the
//! end before it. So a candidate set whose subject range [lowest start, highest end] lies more
//! than `longest_intron` after the head set's range, or ends more than `longest_intron` before
//! it, is not admissible: all HSPs of the earlier set come first in the merged order, and the
//! pair where the sets meet breaks the gap check. Here the HSPs that are alone in their set are
//! looked up by subject start (same context and frame sign as the head, start within
//! `longest_intron` plus the longest HSP of the list), the HSPs of sets of two or more are all
//! kept, and the candidates are tested in NCBI's order (query offset order, inside the same
//! query-offset window) with NCBI's tests. When the query-offset window has no more HSPs than the
//! subject starts looked up, the whole window is tested, as NCBI does. Only candidates that cannot
//! be admissible are left out, and the loop keeps the first best one in that order, so every
//! link, E-value and number is NCBI's.
//!
//! `LOSAT_X_TBNLINKWIN` is in the strict set (`utils/xswitch.rs`): on unless it is `0`.
//! With `LOSAT_X_TBNLINKWINSHADOW` set as well (any
//! value), the whole window of each head whose candidates were chosen is tested too (a head that
//! takes its whole window needs no check), and the search stops (panic) if an admissible
//! candidate is not among the chosen ones. With `LOSAT_X_STATS`
//! set, `main` prints `[X_STATS] TBNLINKWIN heads= window= chosen= shadow_heads=` at exit (the
//! heads, the candidates of their query-offset windows, the candidates tested), when a list was
//! linked with the switch. blastx keeps its own linking (`blastx/linking.rs`).

use super::{sets_admissible, LocalLinkParameters, WorkHsp};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

// No NCBI counterpart: reads LOSAT_X_TBNLINKWIN once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_TBNLINKWIN): whether the linking loop tests the chosen candidates only; a
/// switch of the strict set, on unless it is `0`.
pub(super) fn x_tbn_link_window() -> bool {
    #[cfg(test)]
    if let Some((on, _)) = tests::X_FORCE.get() {
        return on;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_TBNLINKWIN"))
}

// No NCBI counterpart: reads LOSAT_X_TBNLINKWINSHADOW once.
/// EXPERIMENT (LOSAT_X_TBNLINKWIN): whether each head's whole window is tested as well.
fn x_shadow() -> bool {
    #[cfg(test)]
    if let Some((_, shadow)) = tests::X_FORCE.get() {
        return shadow;
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_X_TBNLINKWINSHADOW").is_some())
}

// No NCBI counterpart: LOSAT_X_STATS counters; they do not change any value NCBI computes.
/// Heads, candidates of their query-offset windows, candidates tested, heads checked by the
/// shadow.
static STATS: [AtomicU64; 4] = [const { AtomicU64::new(0) }; 4];

/// The `[X_STATS] TBNLINKWIN` line (printed by `main` at exit) with `LOSAT_X_STATS`, when a list
/// was linked with the switch.
// No NCBI counterpart: prints the counters above; it does not change any value NCBI computes.
pub fn print_summary() {
    let [heads, window, chosen, shadow] = STATS.each_ref().map(|v| v.load(Ordering::Relaxed));
    if heads > 0 && std::env::var_os("LOSAT_X_STATS").is_some() {
        eprintln!(
            "[X_STATS] TBNLINKWIN heads={heads} window={window} chosen={chosen} shadow_heads={shadow}"
        );
    }
}

/// EXPERIMENT (LOSAT_X_TBNLINKWIN): the lookup of one subject list's linking call.
pub(super) struct XLinkWindow {
    /// (context, frame sign, subject start, rank in the query offset order) of every HSP, sorted.
    by_subject: Vec<(usize, i8, i32, usize)>,
    /// The rank in the query offset order of every HSP.
    rank: Vec<usize>,
    /// The largest `s_end - s_start` of the list (0 at least).
    span: i64,
    /// `longest_intron`.
    intron: i64,
    /// The ranks of the HSPs of sets of two or more (sorted), and whether an HSP is among them.
    linked: Vec<usize>,
    is_linked: Vec<bool>,
    /// The chosen ranks, then the chosen HSPs in the query offset order.
    ranks: Vec<usize>,
    chosen: Vec<usize>,
    /// heads, window, chosen, shadow heads (flushed to `STATS` when dropped).
    counts: [u64; 4],
}

impl XLinkWindow {
    // No NCBI counterpart: an index of the list; it computes no value NCBI computes.
    /// The lookup of the list `work`, whose query offset order is `offset_order`.
    pub(super) fn new(
        work: &[WorkHsp],
        offset_order: &[usize],
        link: &LocalLinkParameters,
    ) -> Self {
        let mut rank = vec![0; work.len()];
        for (position, &id) in offset_order.iter().enumerate() {
            rank[id] = position;
        }
        let mut by_subject: Vec<_> = work
            .iter()
            .enumerate()
            .map(|(id, item)| {
                let hsp = &item.value.hsp;
                (
                    item.value.context,
                    hsp.frame.signum(),
                    hsp.s_start,
                    rank[id],
                )
            })
            .collect();
        by_subject.sort_unstable();
        let span = work
            .iter()
            .map(|item| item.value.hsp.s_end as i64 - item.value.hsp.s_start as i64)
            .fold(0, i64::max);
        XLinkWindow {
            by_subject,
            rank,
            span,
            intron: link.longest_intron as i64,
            linked: Vec::new(),
            is_linked: vec![false; work.len()],
            ranks: Vec::new(),
            chosen: Vec::new(),
            counts: [0; 4],
        }
    }

    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:1713-1723
    // ```c
    //        for (index1 = hsp_index_left; index1 < hsp_index_right; ++index1) {
    //            BlastLinkedHSPSet* lhsp = offset_hsp_array[index1];
    //            ...
    //            if (s_LinkedHSPSetsAdmissible(head_hsp, lhsp,
    //                                          link_hsp_params, program)) {
    // ```
    // The same candidates in the same order, less those `s_LinkedHSPSetsAdmissible` refuses
    // whatever their query offsets (see the module text).
    /// The candidates of the head set `head` among `offset_order[left..right]`, in that order:
    /// all but those that cannot be admissible (see the module text), or the whole window when
    /// it is no longer than the subject starts looked up.
    pub(super) fn candidates<'a>(
        &'a mut self,
        work: &[WorkHsp],
        offset_order: &'a [usize],
        head: usize,
        left: usize,
        right: usize,
        link: &LocalLinkParameters,
    ) -> &'a [usize] {
        let context = work[head].value.context;
        let sign = work[head].value.hsp.frame.signum();
        let (mut low, mut high) = (i64::MAX, i64::MIN);
        let mut member = Some(head);
        while let Some(id) = member {
            low = low.min(work[id].value.hsp.s_start as i64);
            high = high.max(work[id].value.hsp.s_end as i64);
            member = work[id].next;
        }
        // An HSP alone in its set can be admissible only with a start in [from, to]: one that
        // starts after `to` begins more than `intron` after the head set's highest end, and one
        // that starts before `from` ends more than `intron` before its lowest start.
        let (from, to) = (low - self.intron - self.span, high + self.intron);
        let first = self
            .by_subject
            .partition_point(|&(c, s, start, _)| (c, s, start as i64) < (context, sign, from));
        let last = self
            .by_subject
            .partition_point(|&(c, s, start, _)| (c, s, start as i64) <= (context, sign, to))
            .max(first);
        self.counts[0] += 1;
        self.counts[1] += (right - left) as u64;
        // A window no longer than the looked-up starts: NCBI's candidates as they are (a
        // superset of the chosen ones, in the same order).
        let plain = right - left <= last - first;
        #[cfg(test)]
        tests::X_PATHS.set({
            let mut paths = tests::X_PATHS.get();
            paths[usize::from(plain)] += 1;
            paths
        });
        if plain {
            self.counts[2] += (right - left) as u64;
            return &offset_order[left..right];
        }
        let window = left..right;
        self.ranks.clear();
        self.ranks.extend(
            self.by_subject[first..last]
                .iter()
                .map(|&(_, _, _, rank)| rank)
                .filter(|rank| window.contains(rank)),
        );
        let linked_left = self.linked.partition_point(|&rank| rank < left);
        let linked_right = self.linked.partition_point(|&rank| rank < right);
        self.ranks
            .extend_from_slice(&self.linked[linked_left..linked_right]);
        self.ranks.sort_unstable();
        self.ranks.dedup();
        self.chosen.clear();
        self.chosen
            .extend(self.ranks.iter().map(|&rank| offset_order[rank]));
        self.counts[2] += self.chosen.len() as u64;
        if x_shadow() {
            self.counts[3] += 1;
            for &candidate in &offset_order[left..right] {
                if sets_admissible(work, head, candidate, link)
                    && self.ranks.binary_search(&self.rank[candidate]).is_err()
                {
                    panic!(
                        "LOSAT_X_TBNLINKWINSHADOW: HSP {candidate} is admissible with the set of HSP {head} but was not chosen"
                    );
                }
            }
        }
        &self.chosen
    }

    // No NCBI counterpart: bookkeeping of the lookup after `s_CombineLinkedHSPSets`.
    /// Notes the set whose head is `head` after a link: its HSPs are kept as candidates from
    /// then on.
    pub(super) fn linked(&mut self, work: &[WorkHsp], head: usize) {
        let mut member = Some(head);
        while let Some(id) = member {
            if !self.is_linked[id] {
                self.is_linked[id] = true;
                let rank = self.rank[id];
                let at = self.linked.partition_point(|&linked| linked < rank);
                self.linked.insert(at, rank);
            }
            member = work[id].next;
        }
    }
}

impl Drop for XLinkWindow {
    fn drop(&mut self) {
        for (total, count) in STATS.iter().zip(self.counts) {
            total.fetch_add(count, Ordering::Relaxed);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::{link_preliminary_hsps, LinkedHspList};
    use crate::algorithm::tblastn::search_gapped::GappedHsp;
    use crate::algorithm::tblastn::stage_d_stats::{LocalContextLength, LocalLinkParameters};
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    use crate::stats::spouge::lookup_protein_gumbel_params;
    use crate::stats::tables::lookup_protein_params_gapped;
    use std::cell::Cell;

    thread_local! {
        /// (switch, shadow) for the linking calls of this test thread.
        pub(super) static X_FORCE: Cell<Option<(bool, bool)>> = const { Cell::new(None) };
        /// Heads of this test thread whose candidates were chosen (`[0]`) and whose whole window
        /// was taken (`[1]`).
        pub(super) static X_PATHS: Cell<[u64; 2]> = const { Cell::new([0; 2]) };
    }

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

    /// A subject list of `count` HSPs: chains of pieces that keep or just break the linking
    /// gaps (`link`), some long ones, repeated offsets, and others anywhere, on `contexts`
    /// query contexts and both strands.
    fn hsps(
        random: &mut Random,
        count: usize,
        contexts: usize,
        subject_aa: i32,
        link: &LocalLinkParameters,
    ) -> Vec<(usize, GappedHsp)> {
        let mut out = Vec::with_capacity(count);
        while out.len() < count {
            let context = random.range(0, contexts as i32 - 1) as usize;
            let sign = if random.next() % 2 == 0 { 1 } else { -1 };
            let mut q = random.range(0, 300);
            let mut s = random.range(0, subject_aa - 2_000);
            let pieces = random.range(1, 6);
            for _ in 0..pieces {
                if out.len() == count {
                    break;
                }
                let length = if random.next() % 50 == 0 {
                    random.range(300, 1_500)
                } else {
                    random.range(8, 60)
                };
                let hsp = GappedHsp {
                    frame: sign * random.range(1, 3) as i8,
                    score: random.range(22, 160),
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
                    // the same query offset again, elsewhere on the subject
                    s = random.range(0, subject_aa - 2_000);
                }
            }
        }
        out
    }

    fn link_with(
        force: (bool, bool),
        input: &[(usize, GappedHsp)],
        contexts: usize,
        subject_nt: i32,
        link: &LocalLinkParameters,
    ) -> LinkedHspList {
        let spec = ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        };
        let gumbel = lookup_protein_gumbel_params(&spec, (subject_nt / 3) as i64).unwrap();
        let gapped = lookup_protein_params_gapped(ScoringMatrix::Blosum62);
        let lengths = vec![
            LocalContextLength {
                length_adjustment: 20,
                eff_searchsp: 380 * (subject_nt as i64 / 3 - 20),
                subject_stat_length: subject_nt as i64 / 3,
            };
            contexts
        ];
        X_FORCE.set(Some(force));
        let linked = link_preliminary_hsps(
            input,
            &vec![400; contexts],
            &lengths,
            subject_nt,
            &vec![gapped; contexts],
            &gumbel,
            link,
        )
        .unwrap();
        X_FORCE.set(None);
        linked
    }

    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:1414-1491 (s_LinkedHSPSetsAdmissible)
    // ```c
    // if (left_hsp->hsp->subject.end < right_hsp->hsp->subject.offset - gap_s) break;
    // if (left_hsp->hsp->subject.end > right_hsp->hsp->subject.offset + overlap || ...
    //     left_hsp->hsp->subject.offset >= right_hsp->hsp->subject.offset || ...) break;
    // ```
    // The chosen candidates give NCBI's links: on random lists the linked list (order, sets,
    // numbers, E-value bits, best E-value) is the one of the whole window, with the shadow
    // checking every head as well, for several sizes, contexts and intron lengths.
    #[test]
    fn x_chosen_candidates_give_the_links_of_the_whole_window() {
        let mut random = Random(0x9E37_79B9_7F4A_7C15);
        let mut links = 0;
        X_PATHS.set([0; 2]);
        for (count, contexts, intron) in [
            (2, 1, 122),
            (7, 1, 40),
            (40, 2, 122),
            (300, 1, 122),
            (300, 3, 900),
            (1_500, 1, 122),
            (1_500, 2, 41),
        ] {
            for round in 0..4 {
                let link = LocalLinkParameters {
                    gap_decay_rate: 0.1,
                    gap_size: 40,
                    overlap_size: 9,
                    longest_intron: intron,
                    cutoff_small_gap: 0,
                };
                let subject_aa = [3_000, 20_000, 200_000, 2_000_000][round];
                let input = hsps(&mut random, count, contexts, subject_aa, &link);
                let subject_nt = subject_aa * 3 + 2;
                let expected = link_with((false, false), &input, contexts, subject_nt, &link);
                for force in [(true, false), (true, true)] {
                    let actual = link_with(force, &input, contexts, subject_nt, &link);
                    let label = format!("{count} HSPs, {contexts} contexts, intron {intron}, subject {subject_aa}, {force:?}");
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
                links += expected.hsps.iter().filter(|hsp| hsp.num > 1).count();
            }
        }
        assert!(links > 100, "the lists must link: {links} HSPs in sets");
        let paths = X_PATHS.get();
        assert!(
            paths[0] > 0 && paths[1] > 0,
            "both kinds of heads must occur: {paths:?}"
        );
    }
}
