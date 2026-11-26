use super::*;
use crate::dbindex::{test_letters, test_volume, Options};

/// BLASTNA codes of `ACGT` letters.
fn blastna(letters: &str) -> Vec<u8> {
    letters
        .bytes()
        .map(|b| match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        })
        .collect()
}

fn search_of(fasta: &str, options: &Options, query: &str, word_size: u64) -> VolumeSeeds {
    search_with(fasta, options, query, word_size, 0)
}

/// The search of `query` with the two-hit window `two_hits` (0: the one-hit search).
fn search_with(
    fasta: &str,
    options: &Options,
    query: &str,
    word_size: u64,
    two_hits: u64,
) -> VolumeSeeds {
    search_masked(fasta, options, query, word_size, two_hits, &[])
}

/// `search_with` with the masked letters `masks` of `-db_hard_mask`.
fn search_masked(
    fasta: &str,
    options: &Options,
    query: &str,
    word_size: u64,
    two_hits: u64,
    masks: &[(u32, Vec<(u32, u32)>)],
) -> VolumeSeeds {
    let (_dir, path) = test_volume(&format!("search{word_size}-{two_hits}"), fasta, options);
    let volume = IndexVolume::open(&path).unwrap();
    let query = blastna(query);
    search(
        &volume,
        &query,
        &[(0, query.len() as u32)],
        SearchOptions {
            word_size,
            two_hits,
        },
        masks,
        &SearchPool::x_serial(),
    )
    .unwrap()
}

/// Each letter of `letters` replaced by another one.
fn other_letters(letters: &str) -> String {
    letters
        .chars()
        .map(|c| match c {
            'A' => 'C',
            'C' => 'G',
            'G' => 'T',
            _ => 'A',
        })
        .collect()
}

/// Two runs on one diagonal, each of at least the word size and shorter than the two-hit
/// search's contig length (twice the word size), `gap` letters apart: the one-hit search saves
/// both; the two-hit search saves the second where its end is at most the window past the first's
/// end and the second's length, and neither otherwise.
#[test]
fn the_two_hit_search_saves_a_seed_with_an_earlier_hit_within_the_window() {
    let s1 = test_letters(400, 23);
    let fasta = format!(">s1\n{s1}\n");
    let gap = 5;
    let second = 114 + gap;
    let query = format!(
        "{}{}{}",
        &s1[100..114],
        other_letters(&s1[114..second]),
        &s1[second..second + 14]
    );
    let both = [
        Seed {
            q_off: 0,
            s_off: 100,
        },
        Seed {
            q_off: 14 + gap as u32,
            s_off: second as u32,
        },
    ];
    let one_hit = search_with(&fasta, &small_options(), &query, 12, 0);
    assert_eq!(one_hit.seeds(0, 0), Some(&both[..]));
    for window in [gap as u64, 40] {
        let two_hits = search_with(&fasta, &small_options(), &query, 12, window);
        assert_eq!(two_hits.seeds(0, 0), Some(&both[1..]), "window {window}");
    }
    let two_hits = search_with(&fasta, &small_options(), &query, 12, gap as u64 - 1);
    assert_eq!(two_hits.seeds(0, 0), None);
    // A run of the contig length is saved without an earlier hit.
    let long = search_with(&fasta, &small_options(), &s1[100..124], 12, 1);
    assert_eq!(
        long.seeds(0, 0),
        Some(
            &[Seed {
                q_off: 0,
                s_off: 100
            }][..]
        )
    );
}

/// With `-db_hard_mask` a masked letter matches no query letter: an exact run of the query over
/// masked letters is the seeds of its unmasked stretches (the one-hit search), and the two-hit
/// search pairs no two hits of a diagonal across a masked letter; a run over masked letters alone
/// is no seed.
#[test]
fn a_masked_letter_ends_a_run_and_pairs_no_hits() {
    let s1 = test_letters(400, 23);
    let fasta = format!(">s1\n{s1}\n");
    let gap = 5;
    let second = 114 + gap;
    // The query is the subject's letters 100..133: one run of 33 letters.
    let query = &s1[100..second + 14];
    let masks = [(0, vec![(114, second as u32)])];
    let whole = [Seed {
        q_off: 0,
        s_off: 100,
    }];
    let parts = [
        Seed {
            q_off: 0,
            s_off: 100,
        },
        Seed {
            q_off: 14 + gap as u32,
            s_off: second as u32,
        },
    ];
    let one_hit = search_with(&fasta, &small_options(), query, 12, 0);
    assert_eq!(one_hit.seeds(0, 0), Some(&whole[..]));
    let one_hit = search_masked(&fasta, &small_options(), query, 12, 0, &masks);
    assert_eq!(one_hit.seeds(0, 0), Some(&parts[..]));
    // The run is of the contig length (24) without the masks; the parts are shorter, and the
    // second's earlier hit is across the masked letters.
    let two_hits = search_with(&fasta, &small_options(), query, 12, 40);
    assert_eq!(two_hits.seeds(0, 0), Some(&whole[..]));
    let two_hits = search_masked(&fasta, &small_options(), query, 12, 40, &masks);
    assert_eq!(two_hits.seeds(0, 0), None);
    // Masked letters 90..140 hold the whole run.
    let masks = [(0, vec![(90, 140)])];
    let none = search_masked(&fasta, &small_options(), query, 12, 0, &masks);
    assert!(!none.has_seeds(0));
}

/// A sequence without masked letters keeps the seeds of the search without `-db_hard_mask` where
/// another sequence of its local id has some, and a mask of a sequence is of its own letters.
#[test]
fn the_masks_are_of_their_own_sequence() {
    let s1 = test_letters(300, 29);
    let s2 = test_letters(300, 31);
    let fasta = format!(">s1\n{s1}\n>s2\n{s2}\n");
    let query = format!("{}{}", &s2[100..160], &s1[100..160]);
    let plain = search_with(&fasta, &small_options(), &query, 28, 0);
    let masks = [(0, vec![(10, 20), (200, 220)])];
    let masked = search_masked(&fasta, &small_options(), &query, 28, 0, &masks);
    assert_eq!(masked.seeds(1, 0), plain.seeds(1, 0));
    assert_eq!(masked.seeds(0, 0), plain.seeds(0, 0));
    assert_eq!(
        masked.seeds(1, 0),
        Some(
            &[Seed {
                q_off: 0,
                s_off: 100
            }][..]
        )
    );
    // Sequence 1's letters 120..130 masked: its run is cut there (s1's run stays whole).
    let masks = [(1, vec![(120, 130)])];
    let masked = search_masked(&fasta, &small_options(), &query, 28, 0, &masks);
    assert_eq!(masked.seeds(0, 0), plain.seeds(0, 0));
    assert_eq!(
        masked.seeds(1, 0),
        Some(
            &[Seed {
                q_off: 30,
                s_off: 130
            }][..]
        )
    );
}

fn small_options() -> Options {
    let mut options = Options::default_options();
    options.hkey_width = 8;
    options
}

#[test]
fn an_exact_run_is_one_seed_at_its_left_end() {
    let s1 = test_letters(400, 3);
    let s2 = test_letters(200, 5);
    let fasta = format!(">s1\n{s1}\n>s2\n{s2}\n");
    // The query is letters 100..200 of the second sequence.
    let seeds = search_of(&fasta, &small_options(), &s2[100..200], 28);
    assert_eq!(seeds.word_size(), 28);
    assert!(!seeds.has_seeds(0));
    assert!(seeds.has_seeds(1));
    assert_eq!(seeds.seeds(0, 0), None);
    assert_eq!(
        seeds.seeds(1, 0),
        Some(
            &[Seed {
                q_off: 0,
                s_off: 100
            }][..]
        )
    );
    // No sequence past the volume's.
    assert!(!seeds.has_seeds(2));
}

#[test]
fn a_seed_goes_to_the_chunk_whose_data_hold_its_left_end() {
    let mut options = small_options();
    options.chunk_size = 100;
    options.chunk_overlap = 20;
    let s1 = test_letters(300, 9);
    let fasta = format!(">s1\n{s1}\n");
    // Chunks start at 0, 80, 160 and 240; the run starts at 170.
    let seeds = search_of(&fasta, &options, &s1[170..230], 28);
    assert_eq!(
        seeds.seeds(0, 2),
        Some(
            &[Seed {
                q_off: 0,
                s_off: 10
            }][..]
        )
    );
    assert_eq!(seeds.seeds(0, 1), None);
    assert_eq!(seeds.seeds(0, 3), None);
}

#[test]
fn a_run_shorter_than_the_word_is_no_seed() {
    let s1 = test_letters(300, 13);
    let fasta = format!(">s1\n{s1}\n");
    let mut query = s1[100..120].to_string();
    query.push_str(&test_letters(40, 21));
    let seeds = search_of(&fasta, &small_options(), &query, 28);
    assert!(!seeds.has_seeds(0));
}

#[test]
fn a_run_stops_at_an_ambiguity_of_the_subject() {
    let left = test_letters(80, 17);
    let right = test_letters(80, 19);
    let fasta = format!(">s1\n{left}NNNN{right}\n");
    // The query runs across the ambiguity with other letters there: the seed is the left part.
    let query = format!("{left}ACGT{right}");
    let seeds = search_of(&fasta, &small_options(), &query, 28);
    let found = seeds.seeds(0, 0).unwrap();
    assert!(found.contains(&Seed { q_off: 0, s_off: 0 }));
    assert!(found.iter().all(|seed| seed.q_off == seed.s_off));
}

#[test]
fn an_empty_sequence_checks_every_chunk_as_ncbi() {
    // `CheckResults` of a sequence whose first chunk is 0 scans from chunk id 0 to the end.
    let seeds = VolumeSeeds {
        word_size: 28,
        map: vec![1, 2, 0],
        results: vec![None, Some(vec![Seed { q_off: 1, s_off: 2 }]), None],
    };
    assert!(!seeds.has_seeds(0));
    assert!(seeds.has_seeds(1));
    assert!(seeds.has_seeds(2));
    assert_eq!(seeds.seeds(2, 0), None);
    assert_eq!(seeds.seeds(1, 0), Some(&[Seed { q_off: 1, s_off: 2 }][..]));
}

#[test]
fn an_oid_is_located_in_its_volume() {
    let volume = |start_oid, n_oids| IndexedVolume {
        path: Vec::new(),
        start_oid,
        n_oids,
        has_index: true,
    };
    let new = LoadedIndex::new(false, vec![volume(0, 3), volume(3, 2), volume(5, 4)]);
    assert_eq!(new.locate(0), (0, 0));
    assert_eq!(new.locate(4), (1, 1));
    assert_eq!(new.locate(5), (2, 0));
    assert_eq!(new.locate(12), (2, 7));
    let old = LoadedIndex::new(true, new.volumes.clone());
    assert_eq!(old.locate(4), (1, 1));
    assert_eq!(old.locate(8), (2, 3));
    // Past every volume: the first volume, with the OID itself.
    assert_eq!(old.locate(9), (0, 9));
}

/// A chunk of a sequence past the sequence's chunk ids lands in a later sequence's (refused);
/// past the last sequence's, or for a sequence past the volume's, NCBI finds no list.
#[test]
fn a_chunk_past_the_sequences_chunks_is_a_later_sequences() {
    let seeds = VolumeSeeds {
        word_size: 28,
        map: vec![1, 3, 4],
        results: vec![None; 4],
    };
    assert!(!seeds.chunk_of_a_later_sequence(0, 1));
    assert!(seeds.chunk_of_a_later_sequence(0, 2));
    assert!(!seeds.chunk_of_a_later_sequence(1, 0));
    assert!(seeds.chunk_of_a_later_sequence(1, 1));
    assert!(!seeds.chunk_of_a_later_sequence(2, 7));
    assert!(!seeds.chunk_of_a_later_sequence(3, 5));
}

/// `Pass::divides` is `n % d == 0` (NCBI's offset iterator) for the passes it can hold.
#[test]
fn a_pass_divides_as_the_remainder() {
    let mut ns: Vec<u32> = (0..3000).collect();
    ns.extend([
        u32::MAX,
        u32::MAX - 1,
        1 << 31,
        (1 << 31) - 1,
        4_294_967_291,
        123_456_789,
    ]);
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    for _ in 0..2000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ns.push(state as u32);
    }
    let mut ds: Vec<u64> = (1..=600).collect();
    ds.extend([
        u64::from(u32::MAX) - 1,
        u64::from(u32::MAX),
        u64::from(u32::MAX) + 1,
        u64::MAX,
        1 << 16,
        65_537,
        2_147_483_648,
    ]);
    for &d in &ds {
        let pass = Pass::new(d);
        for &n in &ns {
            assert_eq!(pass.divides(n), u64::from(n) % d == 0, "{n} % {d}");
        }
    }
}

/// `query_byte` packs four letters as NCBI's byte loops do one by one, and names the letters that
/// are not bases: the first from the left (`extend_right`) and from the right (`extend_left`).
#[test]
fn four_query_letters_pack_as_one_byte() {
    for code in 0u32..1 << 16 {
        let letters: Vec<u8> = (0..4).map(|k| ((code >> (4 * k)) & 0xF) as u8).collect();
        let (byte, ambiguous) = query_byte(&letters, 0).unwrap();
        let first = letters.iter().position(|&letter| letter > 3);
        let last = letters.iter().rposition(|&letter| letter > 3);
        assert_eq!(ambiguous != 0, first.is_some(), "{letters:?}");
        if let (Some(first), Some(last)) = (first, last) {
            assert_eq!(ambiguous.trailing_zeros() / 8, first as u32, "{letters:?}");
            assert_eq!(
                ambiguous.leading_zeros() / 8,
                3 - last as u32,
                "{letters:?}"
            );
        } else {
            let scalar = letters
                .iter()
                .fold(0u8, |byte, &letter| (byte << 2) + letter);
            assert_eq!(byte, scalar, "{letters:?}");
        }
    }
    assert_eq!(query_byte(&[0, 1, 2], 0), Err(Corrupt));
    assert_eq!(query_byte(&[0, 1, 2, 3], -1), Err(Corrupt));
}

/// The seeds do not depend on when the roots are turned into seeds (NCBI's `LIMIT_ROOTS`) nor on
/// the threads: a volume of several local ids, lookup segments with gaps, an empty one and an
/// ambiguity, Nmers in many blocks, the one-hit and the two-hit search, with and without masks.
#[test]
fn any_round_limit_and_thread_count_give_the_same_seeds() {
    let mut options = small_options();
    options.chunk_size = 10_000;
    options.chunk_overlap = 100;
    let seqs: Vec<String> = (0..12)
        .map(|k| test_letters(9_000 + 700 * k, 31 + k as u64))
        .collect();
    let fasta: String = seqs
        .iter()
        .enumerate()
        .map(|(k, s)| format!(">s{k}\n{s}\n"))
        .collect();
    let (_dir, path) = test_volume("rounds", &fasta, &options);
    let volume = IndexVolume::open(&path).unwrap();
    assert!(volume.lid_map().len() / 4 >= 3);
    let mut letters = String::new();
    for k in [1usize, 4, 7, 10] {
        letters.push_str(&seqs[k][1000..2600]);
    }
    let mut query = blastna(&letters);
    for at in (50..query.len()).step_by(97) {
        query[at] = (query[at] + 1) % 4;
    }
    query[3000] = 14;
    let segments = vec![
        (0, 1200),
        (1200, 1200),
        (1300, 3500),
        (3600, query.len() as u32),
    ];
    let masks = [(4u32, vec![(1500, 1550), (2000, 2003)])];
    for masks in [&[][..], &masks[..]] {
        for two_hits in [0, 40] {
            let options = SearchOptions {
                word_size: 16,
                two_hits,
            };
            let serial = SearchPool::x_serial();
            let of = |pool: &SearchPool<'_>, limit: u64| {
                search_with_limit(&volume, &query, &segments, options, masks, pool, limit).unwrap()
            };
            let reference = of(&serial, LIMIT_ROOTS);
            let seeds: usize = reference.results.iter().flatten().map(Vec::len).sum();
            assert!(seeds > 10, "{seeds} seeds");
            for limit in [1, 2, 3, 7, 64, 1000] {
                assert_eq!(of(&serial, limit), reference, "limit {limit}");
            }
            crate::utils::threading::with_search_pool(4, "blastn", |pool| {
                for limit in [1, 5, 300, LIMIT_ROOTS] {
                    assert_eq!(of(pool, limit), reference, "limit {limit}, 4 threads");
                }
                Ok(())
            })
            .unwrap();
        }
    }
}

/// `query_bytes16` is four `query_byte`s, and none where a letter is not a base.
#[test]
fn sixteen_query_letters_pack_as_four_bytes() {
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    for round in 0..20_000 {
        let letters: Vec<u8> = (0..16)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                if round % 3 == 0 {
                    (state % 16) as u8
                } else {
                    (state % 4) as u8
                }
            })
            .collect();
        let four: Vec<(u8, u32)> = (0..4)
            .map(|k| query_byte(&letters, 4 * k).unwrap())
            .collect();
        let expected = four.iter().all(|&(_, ambiguous)| ambiguous == 0).then(|| {
            four.iter().enumerate().fold(0u32, |word, (k, &(byte, _))| {
                word | u32::from(byte) << (8 * k)
            })
        });
        assert_eq!(query_bytes16(&letters, 0), expected, "{letters:?}");
    }
    assert_eq!(query_bytes16(&[0; 15], 0), None);
    assert_eq!(query_bytes16(&[0; 16], -1), None);
}

#[test]
fn the_seed_list_moves_as_a_list_with_a_cursor() {
    // NCBI's std::list and its cursor, as a Vec and an index, against SeedList under the same
    // random steps, insertions, erasures at and before the cursor, changes of the seed before it
    // and returns to the start.
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    for _ in 0..200 {
        let (mut model, mut it) = (Vec::<TrackedSeed>::new(), 0usize);
        let mut list = SeedList::default();
        for step in 0..400u32 {
            match next(6) {
                0 if it < model.len() => {
                    it += 1;
                    list.step();
                }
                1 if it < model.len() => {
                    model.remove(it);
                    list.erase();
                }
                2 => {
                    let seed = TrackedSeed::new(step, next(1000) as u32, 12, step + 11);
                    model.insert(it, seed);
                    it += 1;
                    list.insert(seed);
                }
                3 if it > 0 => {
                    if next(2) == 0 {
                        model[it - 1].qright += 1;
                        list.before_cursor().unwrap().qright += 1;
                    } else {
                        model.remove(it - 1);
                        it -= 1;
                        list.erase_before_cursor();
                    }
                }
                4 => {
                    it = 0;
                    list.rewind();
                }
                _ => {}
            }
            assert_eq!(list.at(), model.get(it).copied());
            let gap = list.gap_end - list.gap_start;
            assert!(list.at < list.gap_start || list.at >= list.gap_end || gap == 0);
            assert_eq!(list.buf.len() - gap, model.len());
            let before = if list.at >= list.gap_end {
                list.at - gap
            } else {
                list.at
            };
            assert_eq!(before, it);
            assert_eq!(
                list.before_cursor().copied(),
                it.checked_sub(1).map(|at| model[at])
            );
        }
        assert_eq!(list.into_ordered(), model);
    }
}
