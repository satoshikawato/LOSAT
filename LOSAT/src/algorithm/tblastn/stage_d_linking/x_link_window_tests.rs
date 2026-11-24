use super::super::{link_preliminary_hsps, LinkedHspList};
use crate::algorithm::tblastn::search_gapped::GappedHsp;
use crate::algorithm::tblastn::stage_d_stats::{LocalContextLength, LocalLinkParameters};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::spouge::lookup_protein_gumbel_params;
use crate::stats::tables::lookup_protein_params_gapped;
use std::cell::Cell;

thread_local! {
    /// `Some(false)`: the linking calls of this test thread test NCBI's whole window.
    pub(super) static X_FORCE: Cell<Option<bool>> = const { Cell::new(None) };
    /// `Some(on)`: `LOSAT_X_TBNLINKSETS` of the linking calls of this test thread.
    pub(super) static X_SETS: Cell<Option<bool>> = const { Cell::new(None) };
    /// Heads of this test thread whose candidates were chosen (`[0]`) and whose whole window
    /// was taken (`[1]`), sets brought in by the lookup (`[2]`), and those of them whose first
    /// HSP lies in the window, not met by the lookup (`[3]`).
    pub(super) static X_PATHS: Cell<[u64; 4]> = const { Cell::new([0; 4]) };
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
        let sign = if random.next().is_multiple_of(2) {
            1
        } else {
            -1
        };
        let mut q = random.range(0, 300);
        let mut s = random.range(0, subject_aa - 2_000);
        let pieces = random.range(1, 6);
        for _ in 0..pieces {
            if out.len() == count {
                break;
            }
            let length = if random.next().is_multiple_of(50) {
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
            q = (hsp.q_end + random.range(-link.overlap_size - 2, link.gap_size + 2)).clamp(0, 390);
            s = hsp.s_end + random.range(-link.overlap_size - 2, link.longest_intron + 5);
            if random.next().is_multiple_of(8) {
                // the same query offset again, elsewhere on the subject
                s = random.range(0, subject_aa - 2_000);
            }
        }
    }
    out
}

/// A subject list of repeat copies: `families` patterns of two to five pieces that keep the
/// linking gaps in the query, each copied to `count / families` places on the subject (on
/// `contexts` query contexts and both strands) with the subject gaps drawn per copy (most
/// within `longest_intron`, some just past it) and scores that fall along the copy, so that
/// the lists form many small sets spread over the subject, whose HSPs reach past the
/// query-offset windows of the heads, and later sets of a copy link with its earlier ones.
fn repeats(
    random: &mut Random,
    count: usize,
    families: usize,
    contexts: usize,
    subject_aa: i32,
    link: &LocalLinkParameters,
) -> Vec<(usize, GappedHsp)> {
    let mut out = Vec::with_capacity(count);
    let patterns: Vec<(usize, Vec<(i32, i32)>)> = (0..families)
        .map(|_| {
            let context = random.range(0, contexts as i32 - 1) as usize;
            let mut q = random.range(0, 150);
            let pieces = (0..random.range(2, 5))
                .map(|_| {
                    let piece = (q, q + random.range(8, 30));
                    q = piece.1 + random.range(-link.overlap_size + 1, link.gap_size / 2);
                    piece
                })
                .collect();
            (context, pieces)
        })
        .collect();
    while out.len() < count {
        let (context, pieces) = &patterns[out.len() % families];
        let sign = if random.next().is_multiple_of(2) {
            1
        } else {
            -1
        };
        let frame = sign * random.range(1, 3) as i8;
        let mut s = random.range(0, subject_aa - 2_000);
        for (piece, &(q_start, q_end)) in pieces.iter().enumerate() {
            if out.len() == count {
                break;
            }
            let fall = 12 * piece as i32;
            let s_end = s + (q_end - q_start) + random.range(-2, 2);
            out.push((
                *context,
                GappedHsp {
                    frame,
                    score: random.range(60 - fall, 110 - fall),
                    q_start,
                    q_end,
                    q_gapped_start: q_start,
                    s_start: s,
                    s_end,
                    s_gapped_start: s,
                },
            ));
            s = s_end
                + if random.next().is_multiple_of(6) {
                    link.longest_intron + random.range(1, 5)
                } else {
                    random.range(-link.overlap_size + 1, link.longest_intron)
                };
        }
    }
    out
}

/// The linked list of `input`: NCBI's whole window (`sets` = `None`), or the chosen
/// candidates with `LOSAT_X_TBNLINKSETS` on or off.
fn link_with(
    sets: Option<bool>,
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
    X_FORCE.set(Some(sets.is_some()));
    X_SETS.set(sets);
    let linked = link_preliminary_hsps(
        input,
        &vec![400; contexts],
        &lengths,
        subject_nt,
        &vec![gapped; contexts],
        Some(&gumbel),
        link,
    )
    .unwrap();
    X_FORCE.set(None);
    X_SETS.set(None);
    linked
}

/// Links `input` the three ways and compares the chosen candidates' lists (switch on and off)
/// with the whole window's: order, sets, numbers, E-value bits, best E-value. Returns the
/// HSPs in sets.
fn check(
    input: &[(usize, GappedHsp)],
    contexts: usize,
    subject_aa: i32,
    link: &LocalLinkParameters,
    label: &str,
) -> usize {
    let subject_nt = subject_aa * 3 + 2;
    let expected = link_with(None, input, contexts, subject_nt, link);
    for sets in [true, false] {
        let actual = link_with(Some(sets), input, contexts, subject_nt, link);
        let label = format!("{label}, LOSAT_X_TBNLINKSETS {sets}");
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
    expected.hsps.iter().filter(|hsp| hsp.num > 1).count()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1482-1488 (s_LinkedHSPSetsAdmissible)
// The chosen candidates give NCBI's links: on random lists and on lists of repeat copies
// (many small sets over the whole subject, several contexts) the linked list (order, sets,
// numbers, E-value bits, best E-value) is the one of the whole window, with
// `LOSAT_X_TBNLINKSETS` on and off, for several sizes, contexts and intron lengths; heads take
// both paths and the lookup brings sets in.
#[test]
fn x_chosen_candidates_give_the_links_of_the_whole_window() {
    let mut random = Random(0x9E37_79B9_7F4A_7C15);
    let mut links = 0;
    X_PATHS.set([0; 4]);
    for (count, contexts, intron) in [
        (2, 1, 122),
        (7, 1, 40),
        (40, 2, 122),
        (300, 1, 122),
        (300, 3, 900),
        (1_500, 1, 122),
        (1_500, 2, 41),
        // `-max_intron_length` 5 to 7 and 8 to 10 (gaps of 1 and 2 letters) and a gap longer
        // than every subject here (`core::link_hsps::translated_longest_intron`; the lists
        // place up to six pieces a gap apart, which stays within `i32` for this value).
        (300, 2, 1),
        (1_500, 1, 2),
        (300, 1, 100_000_000),
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
            let label =
                format!("{count} HSPs, {contexts} contexts, intron {intron}, subject {subject_aa}");
            links += check(&input, contexts, subject_aa, &link, &label);
        }
    }
    assert!(links > 100, "the lists must link: {links} HSPs in sets");
    let random_paths = X_PATHS.get();
    let mut repeat_links = 0;
    for (count, families, contexts, intron) in [
        (60, 2, 1, 122),
        (600, 3, 1, 122),
        (600, 5, 3, 40),
        (2_000, 4, 2, 122),
        (2_000, 1, 1, 900),
    ] {
        for subject_aa in [5_000, 20_000, 200_000, 2_000_000] {
            let link = LocalLinkParameters {
                gap_decay_rate: 0.1,
                gap_size: 40,
                overlap_size: 9,
                longest_intron: intron,
                cutoff_small_gap: 0,
            };
            let input = repeats(&mut random, count, families, contexts, subject_aa, &link);
            let label = format!(
                "repeats: {count} HSPs, {families} families, {contexts} contexts, \
                     intron {intron}, subject {subject_aa}"
            );
            repeat_links += check(&input, contexts, subject_aa, &link, &label);
        }
    }
    let paths = X_PATHS.get();
    assert!(
        repeat_links > 1_000,
        "the repeats must link: {repeat_links} HSPs in sets"
    );
    assert!(
        paths[0] > 0 && paths[1] > 0,
        "both kinds of heads must occur: {paths:?}"
    );
    assert!(
        random_paths[2] > 0 && paths[2] - random_paths[2] > 1_000,
        "the lookup must bring sets in, on both kinds of lists: {random_paths:?} {paths:?}"
    );
    assert!(
        paths[3] - random_paths[3] > 100,
        "the repeats must bring in sets whose first HSP the lookup did not meet: {paths:?}"
    );
}
