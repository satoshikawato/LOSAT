use std::sync::atomic::AtomicU32;

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
    let path = test_volume(&format!("search{word_size}"), fasta, options);
    let volume = IndexVolume::open(&path).unwrap();
    let query = blastna(query);
    search(
        &volume,
        &query,
        &[(0, query.len() as u32)],
        SearchOptions {
            word_size,
            two_hits: 0,
        },
    )
    .unwrap()
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
        path: String::new(),
        start_oid,
        n_oids,
    };
    let new = LoadedIndex {
        old_style: false,
        volumes: vec![volume(0, 3), volume(3, 2), volume(5, 4)],
        searches: AtomicU32::new(0),
    };
    assert_eq!(new.locate(0), (0, 0));
    assert_eq!(new.locate(4), (1, 1));
    assert_eq!(new.locate(5), (2, 0));
    assert_eq!(new.locate(12), (2, 7));
    let old = LoadedIndex {
        old_style: true,
        volumes: new.volumes.clone(),
        searches: AtomicU32::new(0),
    };
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
