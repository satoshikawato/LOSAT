use super::*;

const ALL_TEMPLATES: [DiscTemplateType; 12] = [
    DiscTemplateType::T11_16Coding,
    DiscTemplateType::T11_16Optimal,
    DiscTemplateType::T12_16Coding,
    DiscTemplateType::T12_16Optimal,
    DiscTemplateType::T11_18Coding,
    DiscTemplateType::T11_18Optimal,
    DiscTemplateType::T12_18Coding,
    DiscTemplateType::T12_18Optimal,
    DiscTemplateType::T11_21Coding,
    DiscTemplateType::T11_21Optimal,
    DiscTemplateType::T12_21Coding,
    DiscTemplateType::T12_21Optimal,
];

fn template_length(t: DiscTemplateType) -> usize {
    match t as i32 {
        1..=4 => 16,
        5..=8 => 18,
        _ => 21,
    }
}

fn pack(bases: &[u8]) -> Vec<u8> {
    let mut packed = vec![0u8; bases.len().div_ceil(4)];
    for (i, &b) in bases.iter().enumerate() {
        packed[i / 4] |= b << (2 * (3 - (i % 4)));
    }
    packed
}

fn bases(len: usize, mut seed: u64) -> Vec<u8> {
    (0..len)
        .map(|_| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 33) & 3) as u8
        })
        .collect()
}

/// The index of the word at `start` read base by base, as `s_FillDiscMBTable` builds
/// the accumulator (blast_nalookup.c:752-775).
fn index_from_bases(seq: &[u8], start: usize, t: DiscTemplateType) -> u32 {
    let mut accum = 0u64;
    for &b in &seq[..start + template_length(t)] {
        accum = (accum << 2) | b as u64;
    }
    compute_discontiguous_index(accum, t)
}

// The NCBI template strings (blast_nalookup.h:195-210) read as masks over the template
// positions: an index uses exactly the weight's bases, in order.
#[test]
fn templates_select_the_documented_positions() {
    let patterns: [(DiscTemplateType, &str); 12] = [
        (DiscTemplateType::T12_16Optimal, "1110110110110111"),
        (DiscTemplateType::T11_16Optimal, "1110010110110111"),
        (DiscTemplateType::T12_18Optimal, "111010110010110111"),
        (DiscTemplateType::T11_18Optimal, "111010010110010111"),
        (DiscTemplateType::T12_21Optimal, "111010010110010010111"),
        (DiscTemplateType::T11_21Optimal, "111010010100010010111"),
        (DiscTemplateType::T12_16Coding, "1111101101101101"),
        (DiscTemplateType::T11_16Coding, "1101101101101101"),
        (DiscTemplateType::T12_18Coding, "101101101101101101"),
        (DiscTemplateType::T11_18Coding, "101101100101101101"),
        (DiscTemplateType::T12_21Coding, "100101101101100101101"),
        (DiscTemplateType::T11_21Coding, "100101100101100101101"),
    ];
    for (t, pattern) in patterns {
        let len = pattern.len();
        assert_eq!(len, template_length(t));
        // Each weighted position alone, as a 3 base, gives the index bits of that
        // position; unweighted positions give no bits.
        let weighted: Vec<usize> = pattern
            .bytes()
            .enumerate()
            .filter(|&(_, c)| c == b'1')
            .map(|(i, _)| i)
            .collect();
        for i in 0..len {
            let mut seq = vec![0u8; len];
            seq[i] = 3;
            let index = index_from_bases(&seq, 0, t);
            match weighted.iter().position(|&w| w == i) {
                Some(rank) => {
                    let shift = 2 * (weighted.len() - 1 - rank);
                    assert_eq!(index, 3 << shift, "{t:?} position {i}");
                }
                None => assert_eq!(index, 0, "{t:?} position {i}"),
            }
        }
    }
}

#[test]
fn scans_give_the_index_of_every_word_from_any_start() {
    let seq = bases(300, 7);
    let packed = pack(&seq);
    for t in ALL_TEMPLATES {
        let tl = template_length(t) as i32;
        for kind in [
            DiscScanSubjectKind::General,
            DiscScanSubjectKind::TwoTemplates,
            DiscScanSubjectKind::Scan11_18,
            DiscScanSubjectKind::Scan11_21,
        ] {
            if (kind == DiscScanSubjectKind::Scan11_18 && t != DiscTemplateType::T11_18Coding)
                || (kind == DiscScanSubjectKind::Scan11_21 && t != DiscTemplateType::T11_21Coding)
            {
                continue;
            }
            let second = if kind == DiscScanSubjectKind::TwoTemplates {
                if t as i32 % 2 == 0 {
                    continue;
                }
                t.next()
            } else {
                DiscTemplateType::Contiguous
            };
            for start in 0..9 {
                let end = seq.len() as i32 - tl;
                let mut seen = Vec::new();
                disc_word_scan_subject(
                    kind,
                    &packed,
                    t,
                    second,
                    tl,
                    [start, end],
                    &mut |s_off, index, index2| seen.push((s_off, index, index2)),
                );
                assert_eq!(seen.len() as i32, end - start + 1);
                for (s_off, index, index2) in seen {
                    assert_eq!(index, index_from_bases(&seq, s_off as usize, t));
                    if kind == DiscScanSubjectKind::TwoTemplates {
                        assert_eq!(index2, Some(index_from_bases(&seq, s_off as usize, second)));
                    } else {
                        assert_eq!(index2, None);
                    }
                }
            }
        }
    }
}

// No NCBI counterpart (a check of the scans' rounds, `x_disc_rounds`): every scan (two
// templates' too) gives the index of every word of its range for subjects of every length near
// a template's, ranges
// ending at each of the last positions and from each place of a byte, and a scan that its
// callback stops (`disc_word_scan_subject_until`) gives exactly the words up to that one.
#[test]
fn scans_give_every_word_of_any_range_and_stop_where_told() {
    for t in ALL_TEMPLATES {
        let tl = template_length(t);
        for kind in [
            DiscScanSubjectKind::General,
            DiscScanSubjectKind::TwoTemplates,
            DiscScanSubjectKind::Scan11_18,
            DiscScanSubjectKind::Scan11_21,
        ] {
            if (kind == DiscScanSubjectKind::Scan11_18 && t != DiscTemplateType::T11_18Coding)
                || (kind == DiscScanSubjectKind::Scan11_21 && t != DiscTemplateType::T11_21Coding)
            {
                continue;
            }
            let second = if kind == DiscScanSubjectKind::TwoTemplates {
                if t as i32 % 2 == 0 {
                    continue;
                }
                Some(t.next())
            } else {
                None
            };
            for len in (tl..tl + 41).chain([257, 258, 259, 260]) {
                let seq = bases(len, len as u64 + 11);
                let packed = pack(&seq);
                let last = (len - tl) as i32;
                for start in 0..8.min(last + 1) {
                    for end in (start..=last).rev().take(6) {
                        let words = |stop_after: usize| {
                            let mut seen = Vec::new();
                            disc_word_scan_subject_until(
                                kind,
                                &packed,
                                t,
                                second.unwrap_or(DiscTemplateType::Contiguous),
                                tl as i32,
                                [start, end],
                                &mut |s_off, index, index2| {
                                    seen.push((s_off, index, index2));
                                    seen.len() < stop_after
                                },
                            );
                            seen
                        };
                        let all = words(usize::MAX);
                        let want: Vec<_> = (start..=end)
                            .map(|s| {
                                let at = |t| index_from_bases(&seq, s as usize, t);
                                (s, at(t), second.map(at))
                            })
                            .collect();
                        assert_eq!(all, want, "{t:?} {kind:?} len {len} [{start}, {end}]");
                        for stop_after in [1, 2, 3, 4, 5, 7, 9] {
                            let n = stop_after.min(want.len());
                            assert_eq!(words(stop_after), want[..n], "stop after {stop_after}");
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn template_choice_follows_weight_length_and_type() {
    assert_eq!(
        get_disc_template_type(11, 18, DiscWordType::Coding),
        DiscTemplateType::T11_18Coding
    );
    assert_eq!(
        get_disc_template_type(11, 18, DiscWordType::TwoTemplates),
        DiscTemplateType::T11_18Coding
    );
    assert_eq!(
        get_disc_template_type(12, 21, DiscWordType::Optimal),
        DiscTemplateType::T12_21Optimal
    );
    assert_eq!(
        get_disc_template_type(13, 18, DiscWordType::Coding),
        DiscTemplateType::Contiguous
    );
    assert_eq!(
        get_disc_template_type(11, 17, DiscWordType::Coding),
        DiscTemplateType::Contiguous
    );
    assert_eq!(
        DiscTemplateType::T11_18Coding.next(),
        DiscTemplateType::T11_18Optimal
    );
    assert_eq!(
        DiscTemplateType::T12_21Coding.next(),
        DiscTemplateType::T12_21Optimal
    );
}
