use super::super::input::parse_fasta;
use super::super::{parameters, query_setup::prepare_queries};
use super::*;

/// A blastx search's lookup table, contexts and options for one nucleotide query of `letters`
/// with `extra` arguments (the task), as `seed.rs`'s tests build them.
fn search(
    letters: &[u8],
    extra: &[&str],
) -> (
    PreparedQueryBatch,
    Vec<ContextParameters>,
    ResolvedOptions,
    Lookup,
) {
    let fasta = [b">q\n".as_slice(), letters, b"\n"].concat();
    let q = parse_fasta(&fasta, false, false).unwrap();
    let mut args = vec![
        "losat", "blastx", "-query", "query", "-subject", "subject", "-seg", "no",
    ];
    args.extend_from_slice(extra);
    let crate::cli::Commands::Blastx(args) = crate::cli::try_parse_from(args).unwrap().command
    else {
        unreachable!()
    };
    let options = args.resolve().unwrap();
    let batch = prepare_queries(&q, &options).unwrap();
    let mut params = parameters::score_block(&batch);
    parameters::effective_lengths(&batch, &mut params, &options, 1_000_000, 1, None).unwrap();
    let lookup = Lookup::new(&batch, &params, &options).unwrap();
    (batch, params, options, lookup)
}

/// `length` letters of a fixed pseudo-random protein, encoded, with the sentinel after them.
fn subject(length: usize) -> Vec<u8> {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut letters: Vec<u8> = (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            crate::utils::matrix::aa_char_to_ncbistdaa(
                b"ARNDCQEGHILKMFPSTWYV"[(state % 20) as usize],
            )
        })
        .collect();
    letters.push(0);
    letters
}

/// The word hits of the word finders' scan of `subject` over `ranges` (the whole subject
/// without).
fn seeds(
    lookup: &Lookup,
    subject: &[u8],
    ranges: Option<&[(i32, i32)]>,
    window: i32,
) -> Vec<(i32, i32)> {
    let mut seeds = Vec::new();
    match ranges {
        None => lookup.scan(subject, subject.len() - 1, window, |pairs| {
            seeds.extend_from_slice(pairs)
        }),
        Some(ranges) => x_scan(lookup, subject, ranges, window, &mut |pairs| {
            seeds.extend_from_slice(pairs)
        }),
    }
    seeds
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-59 and
// aa_ungapped.c:495-505, :762-769 (cited at `x_scan`): a range `(left, right)` holds the words
// that start at `left..=right - word`; a first range too short for a word holds its first word
// in the two-hit search and stops the one-hit search's scan.
#[test]
fn the_scan_over_ranges_takes_the_words_of_the_ranges() {
    let query = b"ATGGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTG";
    for extra in [&[][..], &["-task", "blastx-fast"][..]] {
        let (_, _, options, lookup) = search(query, extra);
        let word = options.word_size;
        for length in [400usize, 2] {
            let subject = subject(length);
            let length = length as i32;
            for window in [0, 40] {
                let whole = seeds(&lookup, &subject, None, window);
                assert_eq!(
                    seeds(&lookup, &subject, Some(&[(0, length)]), window),
                    whole
                );
                let runs: &[&[(i32, i32)]] = if length == 2 {
                    &[&[(1, 1)], &[(0, 1)]]
                } else {
                    &[
                        &[(10, 20)],
                        &[(0, 0), (50, 51)],
                        &[(2, 30), (100, 200)],
                        &[(length - 5, length - 1)],
                        &[(0, length - 1)],
                    ]
                };
                for &runs in runs {
                    let ranges = scan_ranges(&record(&vec![b'A'; length as usize], runs));
                    let short = ranges[0].1 - word < ranges[0].0;
                    let expected: Vec<(i32, i32)> = whole
                        .iter()
                        .copied()
                        .filter(|&(_, s): &(i32, i32)| {
                            (short && window > 0 && s == 0)
                                || ranges
                                    .iter()
                                    .any(|&(left, right)| left <= s && s <= right - word)
                        })
                        .filter(|_| !(short && window == 0))
                        .collect();
                    assert_eq!(
                        seeds(&lookup, &subject, Some(&ranges), window),
                        expected,
                        "{extra:?} length {length} window {window} runs {runs:?}"
                    );
                }
            }
        }
    }
}

// The reference word finder over one range that covers each subject is the word finder of the
// whole subjects: the same hits, here and on the subjects after it (the same diagonal table).
#[test]
fn one_range_over_the_subject_is_the_whole_scan() {
    let query = b"ATGGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTGAAAATGTTCCCAAGCACCTGGTACGTCGCTCGTAACGATTGCCAAGAAGGACATATCCTG";
    for window in ["0", "40"] {
        let (batch, params, options, lookup) = search(query, &["-window_size", window]);
        let bx = BxLookup::Wide(lookup);
        let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
        let mut whole = Diagonals::new(length, options.window_size).unwrap();
        let mut ranged = Diagonals::new(length, options.window_size).unwrap();
        for letters in [300, 7, 300] {
            let subject = subject(letters);
            let expected = super::super::seed::word_finder_bx(
                &batch,
                &params,
                &options,
                &subject,
                &bx,
                &mut whole,
                |_| {},
            )
            .unwrap();
            let mut seeds = Vec::new();
            let hits = x_word_finder(
                super::super::seed::Blosum62Bx,
                &batch,
                &params,
                &options,
                &subject,
                &bx,
                &mut ranged,
                &[(0, letters as i32)],
                Some(&mut seeds),
            )
            .unwrap();
            assert!(!seeds.is_empty());
            assert_eq!(hits, expected, "window {window} letters {letters}");
        }
    }
}

fn record(letters: &[u8], runs: &[(i32, i32)]) -> FastaRecord {
    FastaRecord {
        internal_id: String::new(),
        title: Vec::new(),
        sequence: letters.to_vec(),
        lowercase_masks: runs.to_vec(),
        source: None,
        around: None,
        skipped: None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213 (cited at
// `scan_ranges`): the left end is the last letter of the run before, included.
#[test]
fn scan_ranges_run_from_the_last_lower_case_letter() {
    assert_eq!(scan_ranges(&record(b"MKVVA", &[(2, 3)])), [(0, 2), (3, 5)]);
    assert_eq!(
        scan_ranges(&record(b"MKVVA", &[(0, 0), (3, 3)])),
        [(0, 0), (0, 3), (3, 5)]
    );
    assert_eq!(scan_ranges(&record(b"MKV", &[(2, 2)])), [(0, 2), (2, 3)]);
    assert_eq!(scan_ranges(&record(b"MKV", &[(0, 2)])), [(0, 0), (2, 3)]);
}

#[test]
fn only_fasta_records_with_runs_have_ranges() {
    let records = [record(b"MKVVA", &[]), record(b"MKVVA", &[(1, 2)])];
    let subjects = Subjects::Fasta(&records);
    assert_eq!(x_subject_scan_ranges(subjects, 0), None);
    assert_eq!(
        x_subject_scan_ranges(subjects, 1),
        Some(vec![(0, 1), (2, 5)])
    );
    assert_eq!(x_subject_scan_ranges(Subjects::Unreadable, 0), None);
}
