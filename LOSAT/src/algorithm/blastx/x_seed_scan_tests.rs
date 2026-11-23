use super::*;
use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
// The frozen pinned-NCBI INIT rows of the session F fixture (one-hit and two-hit windows,
// fresh, retained and reset diagonal tables): the direct loop gives the NCBI hits and the
// same diagonal table as the reference loop.
/// An entry with `last_hit` and the flag set.
fn entry(last_hit: i32) -> DiagStruct {
    let mut entry = DiagStruct::default();
    entry.set_last_hit(last_hit);
    entry.set_flag(1);
    entry
}

/// The direct loop's arms: (the scan on presence bytes, the branch-free two-hit test).
/// EXPERIMENT (LOSAT_X_BXSCANBYTE): its plain loops for every table and for sparse ones, and
/// its AVX2 forms where the CPU has AVX2; the arms with it off scan with the branch-free
/// candidate loop.
fn arms() -> Vec<(XScanByte, bool)> {
    let plain = XScanByte::On { avx2: false };
    #[cfg_attr(
        not(target_arch = "x86_64"),
        expect(unused_mut, reason = "the AVX2 arms are pushed on x86_64 only")
    )]
    let mut arms = vec![
        (XScanByte::Off, false),
        (XScanByte::Off, true),
        (plain, false),
        (plain, true),
        (XScanByte::Sparse, true),
    ];
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        arms.push((XScanByte::On { avx2: true }, false));
        arms.push((XScanByte::On { avx2: true }, true));
    }
    arms
}

#[test]
fn direct_scan_loop_matches_ncbi_hits_and_reference_diagonal_table() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
    let q = parse_fasta(&std::fs::read(dir.join("query.fna")).unwrap(), false, false).unwrap();
    let s = parse_fasta(
        &std::fs::read(dir.join("subject.faa")).unwrap(),
        true,
        false,
    )
    .unwrap();
    let mut subject = s[0]
        .sequence
        .iter()
        .copied()
        .map(crate::utils::matrix::aa_char_to_ncbistdaa)
        .collect::<Vec<_>>();
    subject.push(0);
    for window in [0, 1, 39, 40, 41] {
        let args = [
            "losat",
            "blastx",
            "-query",
            "query",
            "-subject",
            "subject",
            "-seg",
            "no",
            "-comp_based_stats",
            "0",
            "-window_size",
            &window.to_string(),
        ];
        let crate::cli::Commands::Blastx(args) =
            crate::cli::try_parse_from(args)
                .unwrap()
                .command
        else {
            unreachable!()
        };
        let options = args.resolve().unwrap();
        let batch = prepare_queries(&q, &options).unwrap();
        let mut params = parameters::score_block(&batch);
        parameters::effective_lengths(
            &batch,
            &mut params,
            &options,
            s[0].sequence.len() as i64,
            1,
            None,
        )
        .unwrap();
        let lookup = Lookup::new(&batch, &params, &options).unwrap();
        parameters::subject_parameters(
            &batch,
            &mut params,
            &options,
            parameters::SubjectNumbers::one(s[0].sequence.len()),
            None,
        )
        .unwrap();
        let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
        let expected = std::fs::read_to_string(dir.join(format!("window{window}.stage"))).unwrap();
        let initial = expected
            .lines()
            .filter(|x| x.starts_with("INIT\t"))
            .map(|x| {
                let f = x.split('\t').collect::<Vec<_>>();
                InitHsp {
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect::<Vec<_>>();
        let mut reference = Diagonals::new(length, window).unwrap();
        // EXPERIMENT (LOSAT_X_BXSCANBYTE, LOSAT_X_BXHITBF): the direct loop with and without
        // the scan on presence bytes and the branch-free two-hit test (the other arms scan
        // with the branch-free candidate loop), each on its own table, compared alike. (This fixture
        // never fills the offset array, and its `Sparse` arm scans on presence bytes only if
        // its table is sparse: `x_scan_byte_gives_ncbi_scan_calls_on_blastx_tables` takes the
        // capacity stops, `the_byte_scan_takes_sparse_tables_...` the rule.)
        let arms = arms();
        let mut tables: Vec<_> = arms
            .iter()
            .map(|_| Diagonals::new(length, window).unwrap())
            .collect();
        // EXPERIMENT (LOSAT_X_BXUNGAPSENT): the same arms with the sentinel walks, without and
        // with the shadow's comparison.
        let mut sent_tables: Vec<_> = arms
            .iter()
            .map(|_| {
                [
                    Diagonals::new(length, window).unwrap(),
                    Diagonals::new(length, window).unwrap(),
                ]
            })
            .collect();
        let lead: Vec<u8> = std::iter::once(0).chain(subject.iter().copied()).collect();
        for state in 0..4 {
            for d in tables
                .iter_mut()
                .chain(sent_tables.iter_mut().flatten())
                .chain([&mut reference])
            {
                if state == 2 {
                    d.offset = 536_870_910;
                    d.entries.fill(entry(123));
                }
                if state == 3 {
                    d.offset = 536_870_911;
                    d.entries.fill(entry(536_870_900));
                    d.finish_subject(83).unwrap();
                }
            }
            // The branch-free test acts on every two-hit window above the word size here.
            assert_eq!(
                x_hitbf_acts(true, &options, &subject, &reference),
                window > options.word_size,
                "window={window} state={state}"
            );
            let plain = word_finder(
                &batch,
                &params,
                &options,
                &subject,
                &lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            assert_eq!(plain, initial, "window={window} state={state}");
            for (&(scan_byte, hit_bf), d) in arms.iter().zip(&mut tables) {
                let hits = x_word_finder_direct(
                    &batch, &params, &options, &subject, &lookup, d, scan_byte, hit_bf,
                )
                .unwrap();
                assert_eq!(
                    hits, plain,
                    "window={window} state={state} byte={scan_byte:?} hit_bf={hit_bf}"
                );
                assert!(
                    d.offset == reference.offset && d.entries == reference.entries,
                    "window={window} state={state} byte={scan_byte:?} hit_bf={hit_bf}"
                );
            }
            for (&(scan_byte, hit_bf), pair) in arms.iter().zip(&mut sent_tables) {
                for (d, shadow) in pair.iter_mut().zip([false, true]) {
                    let hits = by_word_size!(
                        options.word_size,
                        W => x_word_finder_direct_sent::<_, W>(
                            &batch,
                            &params,
                            &options,
                            &subject,
                            &lookup,
                            d,
                            scan_byte,
                            hit_bf,
                            XSent::for_test(Blosum62Bx, &batch.sequence_start, &lead, shadow),
                        )
                    )
                    .unwrap();
                    let arm = format!(
                        "window={window} state={state} \
                             byte={scan_byte:?} hit_bf={hit_bf} sent shadow={shadow}"
                    );
                    assert_eq!(hits, plain, "{arm}");
                    assert!(
                        d.offset == reference.offset && d.entries == reference.entries,
                        "{arm}"
                    );
                }
            }
        }
        // EXPERIMENT (LOSAT_X_BXHITBF): the edge of the branch-free test's range. At the most
        // offset `XHitBf::fits` allows for this subject (cleared cells) the test acts, and one
        // more it does not; the loops agree at both.
        if window > options.word_size {
            let most = (1 << 30) - 1 - (subject.len() as i32 - 1) - 2 * window;
            for (offset, acts) in [(most, true), (most + 1, false)] {
                for d in tables.iter_mut().chain([&mut reference]) {
                    d.offset = offset;
                    d.entries.fill(DiagStruct::clear(window));
                }
                assert_eq!(
                    x_hitbf_acts(true, &options, &subject, &reference),
                    acts,
                    "window={window} offset={offset}"
                );
                let plain = word_finder(
                    &batch,
                    &params,
                    &options,
                    &subject,
                    &lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap();
                // Only differences of offsets count below 2^30: the NCBI rows of a new table.
                assert_eq!(plain, initial, "window={window} offset={offset}");
                for (&(scan_byte, hit_bf), d) in arms.iter().zip(&mut tables) {
                    let hits = x_word_finder_direct(
                        &batch, &params, &options, &subject, &lookup, d, scan_byte, hit_bf,
                    )
                    .unwrap();
                    let what = format!(
                        "window={window} offset={offset} byte={scan_byte:?} hit_bf={hit_bf}"
                    );
                    assert_eq!(hits, plain, "{what}");
                    assert!(
                        d.offset == reference.offset && d.entries == reference.entries,
                        "{what}"
                    );
                }
            }
        }
        // A subject of one residue (shorter than a word less one letter) and one of two: with a
        // two-hit window the one scan position reads a word from position 0, which runs past
        // the sentinel of the one-residue subject. The loops scan it from a padded copy
        // (seed.rs `x_scan_bytes`; the slices here hold exactly the residues and the
        // sentinel, so an unchecked read past them trips the debug checks), and they agree.
        for n in [1, 2] {
            let short = subject[..n].iter().copied().chain([0]).collect::<Vec<u8>>();
            let mut reference = Diagonals::new(length, window).unwrap();
            let plain = word_finder(
                &batch,
                &params,
                &options,
                &short,
                &lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            assert!(plain.is_empty(), "window={window} {n} residues: {plain:?}");
            for &(scan_byte, hit_bf) in &arms {
                let mut d = Diagonals::new(length, window).unwrap();
                let hits = x_word_finder_direct(
                    &batch, &params, &options, &short, &lookup, &mut d, scan_byte, hit_bf,
                )
                .unwrap();
                assert_eq!(
                    hits, plain,
                    "window={window} {n} residues, byte={scan_byte:?} hit_bf={hit_bf}"
                );
                assert!(
                    d.offset == reference.offset && d.entries == reference.entries,
                    "window={window} {n} residues, byte={scan_byte:?} hit_bf={hit_bf}"
                );
            }
        }
        // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:57-60
        // An offset just below 2^30 (NCBI reaches such offsets at a -window_size near 2^29,
        // from the second subject on): the stores past 2^30 wrap in the 31-bit field, and the
        // loops read the wrapped values back alike. tests/blastx_diag_bitfield.rs compares the
        // rows with NCBI 2.17.0 at such a window. (The offset is past INT4_MAX / 4, so the
        // tables are cleared at the end of the subject.) The offset is outside the range of
        // the branch-free test, so its arms run the loop without it.
        for d in tables.iter_mut().chain([&mut reference]) {
            d.offset = (1 << 30) - 64;
            d.entries.fill(DiagStruct::clear(window));
        }
        assert!(!x_hitbf_acts(true, &options, &subject, &reference));
        let plain = word_finder(
            &batch,
            &params,
            &options,
            &subject,
            &lookup,
            &mut reference,
            |_| {},
        )
        .unwrap();
        for (&(scan_byte, hit_bf), d) in arms.iter().zip(&mut tables) {
            let hits = x_word_finder_direct(
                &batch, &params, &options, &subject, &lookup, d, scan_byte, hit_bf,
            )
            .unwrap();
            assert_eq!(
                hits, plain,
                "window={window} wrapped, byte={scan_byte:?} hit_bf={hit_bf}"
            );
            assert!(
                d.offset == reference.offset && d.entries == reference.entries,
                "window={window} wrapped, byte={scan_byte:?} hit_bf={hit_bf}"
            );
        }
        if !initial.is_empty() {
            assert_ne!(
                plain, initial,
                "window={window}: the wrapped stores change the seeds"
            );
        }
    }
}

/// The statements of `word_finder` up to the two `continue`s of a flag-0 cell
/// (aa_ungapped.c:519-551): the cell after a hit that is too far or overlaps, or `None` for a
/// flag-1 cell or a candidate for the extension.
fn two_hit_reference(
    cell: DiagStruct,
    s: i32,
    offset: i32,
    window: i32,
    word: i32,
) -> Option<DiagStruct> {
    if cell.flag() != 0 {
        return None;
    }
    let diff = s.wrapping_sub(cell.last_hit().wrapping_sub(offset));
    if diff >= window {
        let mut next = cell;
        next.set_last_hit(s.wrapping_add(offset));
        return Some(next);
    }
    (diff < word).then_some(cell)
}

/// `XHitBf::test` on a copy of `cell`: the cell after it and whether it took the hit.
fn bf_outcome(cell: DiagStruct, s: i32, offset: i32, window: i32, word: i32) -> (DiagStruct, bool) {
    let mut next = cell;
    let bf = XHitBf::in_range(offset, window, word);
    let taken = bf.test(&mut next, bf.here(s));
    (next, taken)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:519-551 (module header).
// EXPERIMENT (LOSAT_X_BXHITBF): the branch-free test against the reference statements at the
// edges of its range: offsets from 0 to the most `XHitBf::fits` allows for the subject's length
// (and that it allows no more), word sizes from 1 to window - 1, windows from 2 to 2^28 - 1,
// subject positions at both ends, and cells around the outcomes' edges (diff = wordsize - 1,
// wordsize, window - 1, window), the cleared cell (`-window`), a new table's 0, the lowest and
// highest `last_hit` a flag-0 cell can hold, and flag-1 cells. Past the range a cleared cell
// can be read wrong, which is why the range is checked.
#[test]
fn the_branch_free_test_gives_the_reference_outcomes_in_its_range() {
    // too far, overlap, candidate, flag 1
    let mut outcomes = [0u64; 4];
    for window in [2, 3, 4, 40, 41, 1000, (1 << 28) - 1] {
        let mut words = vec![1, 2, 3, window / 2, window - 1];
        words.retain(|&w| 1 <= w && w < window);
        words.sort_unstable();
        words.dedup();
        for word in words {
            for length in [1usize, 2, 100, 1 << 20] {
                let most = (1i64 << 30) - 1 - length as i64 - 2 * i64::from(window);
                let most = i32::try_from(most).unwrap();
                assert!(XHitBf::fits(most, length, window, word));
                assert!(!XHitBf::fits(most + 1, length, window, word));
                assert!(!XHitBf::fits(-1, length, window, word));
                for offset in [0, most / 2, most] {
                    let high = offset + length as i32;
                    for s in [0, 1, length as i32 / 2, length as i32 - 1] {
                        if !(0..length as i32).contains(&s) {
                            continue;
                        }
                        let here = s + offset;
                        let mut last_hits = vec![-window, 0, high - 1, offset, here + 1];
                        for d in [word - 1, word, window - 1, window, window + 1, 0] {
                            last_hits.push(here - d);
                        }
                        for last_hit in last_hits {
                            for flag in [0, 1] {
                                // A flag-0 cell holds a last_hit in [-window, high).
                                if flag == 0 && !(-window..high).contains(&last_hit) {
                                    continue;
                                }
                                let mut cell = DiagStruct::default();
                                cell.set_last_hit(last_hit);
                                cell.set_flag(flag);
                                let (next, taken) = bf_outcome(cell, s, offset, window, word);
                                let what = format!(
                                    "window={window} word={word} length={length} offset={offset} s={s} last_hit={last_hit} flag={flag}"
                                );
                                match two_hit_reference(cell, s, offset, window, word) {
                                    Some(expected) => {
                                        assert!(taken, "{what}");
                                        assert_eq!(next, expected, "{what}");
                                        outcomes[usize::from(next == cell)] += 1;
                                    }
                                    None => {
                                        assert!(!taken, "{what}");
                                        assert_eq!(next, cell, "{what}");
                                        outcomes[2 + flag as usize] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(outcomes.iter().all(|&n| n > 0), "{outcomes:?}");
    assert!(!XHitBf::fits(0, 10, 3, 3));
    assert!(!XHitBf::fits(0, 10, 3, 0));
    assert!(!XHitBf::fits(0, 10, 0, 3));
    assert!(XHitBf::fits(0, 10, 3, 2));
    // Just past the range: a cleared cell at the end of a subject is too far for the
    // reference (diff = 2^30 + 499), but twice that wraps in an i32 and the branch-free test
    // reads it as an overlap.
    let (window, word, length) = (1000, 3, 1000usize);
    let offset = (1 << 30) - 1500;
    assert!(!XHitBf::fits(offset, length, window, word));
    let cleared = DiagStruct::clear(window);
    let expected = two_hit_reference(cleared, 999, offset, window, word).unwrap();
    assert_ne!(expected, cleared);
    assert_ne!(
        bf_outcome(cleared, 999, offset, window, word),
        (expected, true)
    );
}

/// A blastx lookup table as `Lookup::new` builds it for one nucleotide query of `letters`, with
/// `extra` arguments (the word size), and its offset-array size.
fn blastx_lookup(letters: &[u8], extra: &[&str]) -> Lookup {
    let fasta = [b">q\n".as_slice(), letters, b"\n"].concat();
    let q = parse_fasta(&fasta, false, false).unwrap();
    let mut args = vec![
        "losat", "blastx", "-query", "query", "-subject", "subject", "-seg", "no",
    ];
    args.extend_from_slice(extra);
    let crate::cli::Commands::Blastx(args) =
        crate::cli::try_parse_from(args)
            .unwrap()
            .command
    else {
        unreachable!()
    };
    let options = args.resolve().unwrap();
    let batch = prepare_queries(&q, &options).unwrap();
    let mut params = parameters::score_block(&batch);
    parameters::effective_lengths(&batch, &mut params, &options, 1_000_000, 1, None).unwrap();
    Lookup::new(&batch, &params, &options).unwrap()
}

/// The scan calls of the direct loop over `subject` (residues and the sentinel) with an offset
/// array of `array_size`: (pairs, `scan_range[1]`) after each call, by NCBI's loop or, with
/// `mode`, on presence bytes (`x_scan_one_range_byte`). The buffer has the slack the direct
/// loop gives it and guard entries after it, which no call may write.
fn direct_calls(
    lookup: &Lookup,
    subject: &[u8],
    array_size: usize,
    mode: Option<XScanByte>,
) -> Vec<(Vec<(u32, u32)>, i32)> {
    let Lookup::Standard(table, _) = lookup else {
        unreachable!("a standard table")
    };
    let guard = OffsetPair {
        q_off: 0xdead_beef,
        s_off: 0xfeed_f00d,
    };
    let end = array_size + X_SCAN_BYTE_SLACK;
    let mut pairs = vec![OffsetPair::default(); end];
    pairs.resize(end + 2 * X_SCAN_BYTE_SLACK, guard);
    let length = subject.len() - 1;
    let mut range = [0, 0, length as i32 - table.word_length];
    if range[2] < range[1] {
        range[2] = range[1];
    }
    let mut pad = Vec::new();
    let scanned = x_scan_bytes(subject, range[2], table.word_length, &mut pad);
    let mut out = Vec::new();
    while range[1] <= range[2] {
        let n = match mode {
            Some(mode) => x_scan_one_range_byte(
                mode,
                table,
                lookup.x_scan_presence(mode).expect("presence bytes"),
                scanned,
                &mut pairs,
                array_size as i32,
                &mut range,
            ),
            None => s_blast_aa_scan_subject_one_range(
                table,
                scanned,
                &mut pairs,
                array_size as i32,
                &mut range,
            ),
        };
        out.push((
            pairs[..n as usize]
                .iter()
                .map(|pair| (pair.q_off, pair.s_off))
                .collect(),
            range[1],
        ));
        assert!(
            pairs[end..]
                .iter()
                .all(|pair| (pair.q_off, pair.s_off) == (guard.q_off, guard.s_off)),
            "write past the slack (array_size {array_size}, {mode:?})"
        );
        assert!(out.len() < 1_000_000, "scan does not advance");
    }
    out
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:75-127 (the function compared
// against, as `s_blast_aa_scan_subject_one_range`). EXPERIMENT (LOSAT_X_BXSCANBYTE): on
// blastx's tables (`Lookup::new`: six frames of a nucleotide query, threshold 12) and the
// direct loop's subjects (residues, then the sentinel) and offset array (`Lookup::capacity`,
// and smaller ones that stop the scan after a few words), the scan on presence bytes gives
// NCBI's pairs and `scan_range[1]` after every call, in both forms.
#[test]
fn x_scan_byte_gives_ncbi_scan_calls_on_blastx_tables() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut below = |n: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % n as u64) as usize
    };
    #[cfg_attr(
        not(target_arch = "x86_64"),
        expect(unused_mut, reason = "the AVX2 form is pushed on x86_64 only")
    )]
    let mut forms = vec![XScanByte::On { avx2: false }];
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        forms.push(XScanByte::On { avx2: true });
    }
    let (mut compared, mut stops) = (0usize, 0usize);
    for (case, query_len) in [300usize, 1_200, 6_000].into_iter().enumerate() {
        let query: Vec<u8> = (0..query_len).map(|_| b"ACGT"[below(4)]).collect();
        let lookup = blastx_lookup(&query, &[]);
        let longest = lookup.capacity() - 4096;
        for len in [1, 2, 3, 4, 255, 256, 257, 259, 1_000, 5_000] {
            // ncbistdaa letters 0..=27 (a database protein has 1..=27), the sentinel after.
            let mut subject: Vec<u8> = (0..len).map(|_| below(28) as u8).collect();
            subject.push(0);
            for array_size in [
                longest.max(1),
                longest + 1,
                3 * longest + 2,
                lookup.capacity(),
            ] {
                let want = direct_calls(&lookup, &subject, array_size, None);
                for &mode in &forms {
                    let got = direct_calls(&lookup, &subject, array_size, Some(mode));
                    assert_eq!(
                        got, want,
                        "case {case} subject {len} array_size {array_size} {mode:?}"
                    );
                }
                compared += want.iter().map(|(pairs, _)| pairs.len()).sum::<usize>();
                stops += want.len() - 1;
            }
        }
    }
    assert!(compared > 100_000, "{compared} pairs compared");
    assert!(stops > 100, "{stops} capacity stops");
}

// No NCBI counterpart: LOSAT_X_BXSCANBYTE's choice of the scan per table (speed only; both
// scans give NCBI's pairs). The presence bytes are those of `pv` within the backbone, built
// once and kept with the table; the sparse rule takes the table of a 300-nt query (few
// neighbourhood words) and not that of a 6-kb query; the compressed table keeps its scan.
#[test]
fn the_byte_scan_takes_sparse_tables_where_asked_and_only_standard_ones() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut base = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        b"ACGT"[(state % 4) as usize]
    };
    let short: Vec<u8> = (0..300).map(|_| base()).collect();
    let long: Vec<u8> = (0..6_000).map(|_| base()).collect();
    let plain = XScanByte::On { avx2: false };
    for (letters, sparse) in [(&short, true), (&long, false)] {
        let lookup = blastx_lookup(letters, &[]);
        let Lookup::Standard(table, kept) = &lookup else {
            unreachable!()
        };
        assert!(lookup.x_scan_presence(XScanByte::Off).is_none());
        assert!(kept.presence.get().is_none(), "built before a scan asked");
        let bytes = lookup
            .x_scan_presence(plain)
            .expect("every blastx table fits");
        let same = std::ptr::eq(bytes, lookup.x_scan_presence(plain).unwrap());
        assert!(same, "built once per table");
        for (index, &byte) in bytes.iter().enumerate() {
            let bit = table
                .pv
                .get(index >> 5)
                .map_or(0, |word| (word >> (index & 31)) & 1);
            let inside = index < table.backbone.len();
            assert_eq!(byte as u32, bit * inside as u32, "index {index}");
        }
        assert_eq!(
            lookup.x_scan_presence(XScanByte::Sparse).is_some(),
            sparse,
            "{} letters",
            letters.len()
        );
    }
    let compressed = blastx_lookup(&short, &["-word_size", "5"]);
    assert!(matches!(compressed, Lookup::Compressed(_)));
    assert!(compressed.x_scan_presence(plain).is_none());
}
