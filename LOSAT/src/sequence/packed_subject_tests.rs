use super::*;

/// Every byte the reader can give a nucleotide record: IUPAC letters and `-`, both cases.
const LETTERS: &[u8] = b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn";

fn random_letters(length: usize, seed: u64, ambiguous_share: u64) -> Vec<u8> {
    let mut state = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        state >> 33
    };
    let mut out = Vec::with_capacity(length);
    while out.len() < length {
        let letter = if next() % 100 < ambiguous_share {
            LETTERS[(next() % LETTERS.len() as u64) as usize]
        } else {
            b"ACGTacgt"[(next() % 8) as usize]
        };
        // runs of one letter, so that the runs and intervals are longer than one base
        let repeat = 1 + (next() % 5) as usize;
        for _ in 0..repeat.min(length - out.len()) {
            out.push(letter);
        }
    }
    out
}

/// `length` upper-case bases A, C, G, T.
fn random_bases(length: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn cases() -> Vec<Vec<u8>> {
    let mut cases: Vec<Vec<u8>> = vec![
        Vec::new(),
        b"A".to_vec(),
        b"n".to_vec(),
        b"AC".to_vec(),
        b"ACG".to_vec(),
        b"ACGT".to_vec(),
        b"ACGTN".to_vec(),
        b"NNNNACGTNNNN".to_vec(),
        b"acgtACGTnnNN--RYry".to_vec(),
        LETTERS.to_vec(),
    ];
    // every letter at both ends and alone, with lengths across the remainder byte
    for &letter in LETTERS {
        for length in 1..=9 {
            let mut seq = vec![b'A'; length];
            seq[0] = letter;
            cases.push(seq.clone());
            seq[length - 1] = letter;
            cases.push(seq);
        }
    }
    for (seed, share) in [(1, 0), (2, 5), (3, 50), (4, 100)] {
        for length in [1, 2, 3, 4, 5, 63, 64, 65, 1000, 4099] {
            cases.push(random_letters(length, seed * 10_000 + length as u64, share));
        }
    }
    // blocks of `BLOCK` letters of one case (the fast path), with case changes and other
    // letters inside a block, at its first or last letter and between blocks, a run of N
    // across blocks, and every length around a block
    for length in [31, 32, 33, 35, 36, 63, 64, 65, 67, 68, 96, 97, 4099] {
        for (seed, every) in [(5u64, 0usize), (6, 29), (7, 32), (8, 61)] {
            let mut seq = random_bases(length, seed * 10_000 + length as u64);
            if every > 0 {
                for (n, at) in (0..length).step_by(every).enumerate() {
                    match n % 4 {
                        0 => seq[at] = b'N',
                        1 => seq[at..(at + 40).min(length)].make_ascii_lowercase(),
                        2 => seq[at] = seq[at].to_ascii_lowercase(),
                        _ => seq[at..(at + 70).min(length)].fill(b'N'),
                    }
                }
            }
            cases.push(seq.clone());
            cases.push(seq.to_ascii_lowercase());
        }
    }
    cases
}

#[test]
fn every_form_is_the_one_made_from_the_letters() {
    for seq in cases() {
        let store = PackedNucleotide::from_subject_letters(&seq).expect("IUPAC letters");
        assert_eq!(store.len(), seq.len());
        assert_eq!(store.is_empty(), seq.is_empty());
        let mut out = Vec::new();
        store.ascii_into(&mut out);
        assert_eq!(out, seq, "letters of {:?}", String::from_utf8_lossy(&seq));
        assert_eq!(
            store.packed(),
            crate::core::blast_encoding::encode_subject_ncbi2na_packed(&seq).as_slice()
        );
        store.ncbi4na_into(&mut out);
        let ncbi4na: Vec<u8> = seq
            .iter()
            .map(|&letter| iupacna_to_ncbi4na(letter).unwrap())
            .collect();
        assert_eq!(out, ncbi4na);
        store.blastna_into(&mut out);
        assert_eq!(out, encode_iupac_to_blastna(&seq));
        // The letters' ncbi4na codes resolved by `CRandom(length)`, as letters, the way the
        // live TBLASTN subject translation (`x_subject_ncbi2na`) resolves a record's letters.
        store.resolved_letters_into(&mut out);
        let resolved: Vec<u8> = crate::core::blast_encoding::resolve_ncbi4na_to_ncbi2na(&ncbi4na)
            .iter()
            .map(|&code| b"ACGT"[usize::from(code)])
            .collect();
        assert_eq!(out, resolved);
        let upper = seq.to_ascii_uppercase();
        let n = seq.len();
        // windows at every offset within a byte, of every length across a byte, and
        // windows that start or end inside a run
        let mut windows = vec![(0, n), (0, n / 2), (n / 3, n), (n / 2, n / 2)];
        for from in 0..n.min(9) {
            for to in from..n.min(from + 11) {
                windows.push((from, to));
                windows.push((n - to, n - from));
            }
        }
        for run in store.ambiguity_runs().iter().take(20) {
            let (start, end) = (run.offset as usize, (run.offset + run.length) as usize);
            windows.push((start + 1, end.min(n)));
            windows.push((start.saturating_sub(2), (end + 3).min(n)));
            windows.push((start, (start + 1).min(n)));
        }
        let mut codes = Vec::new();
        for (from, to) in windows {
            store.letters_into(from, to, &mut out);
            assert_eq!(out, &upper[from..to], "letters [{from}, {to})");
            let letters = SubjectNt::Letters(&seq);
            let packed = SubjectNt::Packed(&store);
            assert_eq!((letters.len(), packed.len()), (n, n));
            let mut buffer = Vec::new();
            assert_eq!(
                packed.fetch(from, to, &mut buffer),
                letters
                    .fetch(from, to, &mut Vec::new())
                    .to_ascii_uppercase()
            );
            codes.resize(to - from, 0xFF);
            store.codes_into(from, &mut codes);
            let as_letters: Vec<u8> = codes
                .iter()
                .map(|&code| b"ACGT"[usize::from(code)])
                .collect();
            assert_eq!(as_letters, &resolved[from..to], "codes [{from}, {to})");
            let mut bases = vec![0u8; to - from];
            store.bases_into(from, &mut bases);
            assert_eq!(bases, as_letters, "bases [{from}, {to})");
        }
        let lower: Vec<_> = crate::algorithm::blastn::coordination::collect_lowercase_masks(&seq)
            .iter()
            .map(|interval| (interval.start, interval.end))
            .collect();
        let held: Vec<_> = store
            .lowercase()
            .iter()
            .map(|interval| (interval.start, interval.end))
            .collect();
        assert_eq!(held, lower);
    }
}

/// `blastna_into_par` on a pool writes the bytes of `blastna_into`: every case above, and
/// subjects of several pieces with runs and intervals across the pieces' ends, every length
/// remainder, into a buffer that held other bytes.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[test]
fn the_pool_decode_writes_the_bytes_of_blastna_into() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .expect("a pool of 4 threads");
    let mut cases = cases();
    let bases = 4 * DECODE_PIECE;
    for (seed, share) in [(9u64, 0u64), (10, 3), (11, 60)] {
        for length in [bases - 1, bases, bases + 1, 3 * bases + 2, 3 * bases + 3] {
            let mut seq = random_letters(length, seed * 1_000_003 + length as u64, share);
            // a run of N and a lower-case interval across each piece's end
            for end in (bases..length).step_by(bases) {
                seq[end - 3..(end + 3).min(length)].fill(b'N');
                seq[end - 9..end - 4].make_ascii_lowercase();
            }
            cases.push(seq);
        }
    }
    let mut serial = Vec::new();
    let mut out = vec![0xAB; 17];
    for seq in cases {
        let store = PackedNucleotide::from_subject_letters(&seq).expect("IUPAC letters");
        store.blastna_into(&mut serial);
        pool.install(|| store.blastna_into_par(&mut out));
        assert_eq!(out, serial, "length {}", seq.len());
        assert_eq!(out, encode_iupac_to_blastna(&seq), "length {}", seq.len());
    }
}

#[test]
fn runs_are_maximal_runs_of_one_code() {
    let store = PackedNucleotide::from_subject_letters(b"ANNnNRRACGT--N").unwrap();
    let runs: Vec<_> = store
        .ambiguity_runs()
        .iter()
        .map(|run| (run.offset, run.length, run.ncbi4na))
        .collect();
    assert_eq!(runs, [(1, 4, 15), (5, 2, 5), (11, 2, 0), (13, 1, 15)]);
    assert_eq!(
        store.heap_bytes(),
        14 / 4 + 1 + 4 * 12 + std::mem::size_of::<MaskedInterval>()
    );
}

#[test]
fn a_byte_without_a_code_keeps_the_letters() {
    for seq in [&b"ACGU"[..], b"ACG*", b"AC GT", b"ACGX", &[b'A', 0xFF]] {
        assert!(
            PackedNucleotide::from_subject_letters(seq).is_none(),
            "{seq:?}"
        );
    }
}

/// NCBI's reader applying `runs` to the letters of `packed` (c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:834-842), as
/// ncbi4na codes.
fn ncbi_letters(packed: &[u8], length: usize, runs: &[AmbiguityRun]) -> Vec<u8> {
    let mut seq: Vec<u8> = (0..length)
        .map(|index| 1 << ((packed[index >> 2] >> (6 - 2 * (index & 3))) & 3))
        .collect();
    for run in runs {
        let mut position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= length as i64 {
            break;
        }
        for _ in 0..row_len {
            if position >= 0 && position < length as i64 {
                seq[position as usize] = run.ncbi4na;
            }
            position += 1;
        }
    }
    seq
}

/// Every form of `store` against the ncbi4na codes `codes` and the stored `packed` bytes.
fn assert_stored_forms(store: &PackedNucleotide, packed: &[u8], codes: &[u8], what: &str) {
    let length = codes.len();
    assert_eq!(store.len(), length, "{what}");
    assert_eq!(store.packed(), packed, "{what}");
    let mut out = Vec::new();
    store.ncbi4na_into(&mut out);
    assert_eq!(out, codes, "{what}: ncbi4na");
    store.ascii_into(&mut out);
    let letters: Vec<u8> = codes
        .iter()
        .map(|&code| NCBI4NA_TO_IUPACNA[usize::from(code)])
        .collect();
    assert_eq!(out, letters, "{what}: letters");
    store.letters_into(length / 3, length, &mut out);
    assert_eq!(out, &letters[length / 3..], "{what}: window");
    store.resolved_letters_into(&mut out);
    let resolved: Vec<u8> = (0..length)
        .map(|index| b"ACGT"[usize::from((packed[index >> 2] >> (6 - 2 * (index & 3))) & 3)])
        .collect();
    assert_eq!(out, resolved, "{what}: resolved");
    let runs = store.ambiguity_runs();
    assert!(
        runs.windows(2)
            .all(|pair| pair[0].offset + pair[0].length <= pair[1].offset),
        "{what}: runs in order and apart"
    );
    assert!(store.lowercase().is_empty(), "{what}");
}

#[test]
fn a_stored_sequence_reads_as_ncbis_reader_reads_it() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut checked = 0;
    for name in [
        "seqdb/db/nucl_v4",
        "seqdb/db/nucl_v5",
        "seqdb/db/nucl_mv",
        "db_search/nslice",
    ] {
        let path = root.join("tests").join(name);
        let db =
            crate::seqdb::SeqDb::open(path.to_str().unwrap(), crate::seqdb::SeqType::Nucleotide)
                .unwrap_or_else(|error| panic!("{name}: {error}"));
        let mut next = 0;
        while let Some(range) = db.next_block(next, u64::MAX) {
            next = range.end_ordinal();
            let block = db.read_block(range).unwrap();
            for i in 0..block.len() {
                let packed = block.nucleotide_packed(i);
                let length = block.nucleotide_len(i).unwrap();
                let runs = block.ambiguity_runs(i).unwrap();
                let store = PackedNucleotide::from_stored(packed, length, &runs);
                let blastna = block.nucleotide_blastna(i).unwrap();
                let mut out = Vec::new();
                store.blastna_into(&mut out);
                assert_eq!(out, &blastna[1..=length], "{name} {i}: blastna");
                // The letters of BLASTN's reported sequences (`reported_store`): those of
                // `nucleotide_iupac`, whole and in any range.
                assert_eq!(
                    block.nucleotide_stored(i).unwrap(),
                    (length, runs.clone()),
                    "{name} {i}"
                );
                let iupac = block.nucleotide_iupac(i).unwrap();
                store.letters_into(0, length, &mut out);
                assert_eq!(out, iupac, "{name} {i}: IUPAC");
                for (from, to) in [
                    (0, 1),
                    (length / 3, length - length / 4),
                    (length - 1, length),
                ] {
                    store.letters_into(from, to, &mut out);
                    assert_eq!(out, &iupac[from..to], "{name} {i}: IUPAC {from}..{to}");
                }
                let codes = ncbi_letters(packed, length, &runs);
                assert_stored_forms(&store, packed, &codes, &format!("{name} {i}"));
                // The store of the sequence read alone, in
                // the buffer of its read, is the same.
                let alone = db
                    .x_read_sequence_mapped(range.first_ordinal + i as u32)
                    .unwrap();
                assert_eq!(
                    PackedNucleotide::x_from_stored_buffer(
                        alone.x_into_nucleotide_packed(),
                        length,
                        &runs
                    ),
                    store,
                    "{name} {i}: the store of the read buffer"
                );
                checked += 1;
            }
        }
        assert_eq!(next, db.num_seqs(), "{name}");
    }
    assert!(checked > 20, "{checked} sequences");
}

#[test]
fn stored_runs_out_of_order_are_applied_as_ncbi_applies_them() {
    let letters = random_letters(203, 7, 0);
    let plain = PackedNucleotide::from_subject_letters(&letters).unwrap();
    let packed = plain.packed().to_vec();
    let run = |offset: u32, length: u32, ncbi4na: u8| AmbiguityRun {
        offset,
        length,
        ncbi4na,
    };
    let cases: Vec<Vec<AmbiguityRun>> = vec![
        vec![],
        vec![run(0, 3, 15), run(3, 1, 5), run(200, 3, 0)],
        // overlapping, out of order, a plain base's code, and a run that passes the end
        vec![
            run(10, 8, 15),
            run(12, 2, 10),
            run(5, 1, 1),
            run(198, 40, 6),
        ],
        // a run past the end stops the ones after it
        vec![run(20, 2, 15), run(203, 1, 3), run(30, 2, 15)],
        // a negative offset (NCBI's signed position) starts before the sequence
        vec![run(u32::MAX - 1, 4, 9), run(1, 1, 15)],
        vec![run(50, 4, 15), run(50, 4, 15), run(52, 6, 3)],
    ];
    for (number, runs) in cases.iter().enumerate() {
        let store = PackedNucleotide::from_stored(&packed, letters.len(), runs);
        let codes = ncbi_letters(&packed, letters.len(), runs);
        assert_stored_forms(&store, &packed, &codes, &format!("case {number}"));
    }
}

/// A view of stored bytes with the runs of `x_stored_runs` (a sequence of a database block,
/// read in place) reads the letters and bases of the store of the same bytes and runs, which
/// are NCBI's (`ncbi_letters`): random bytes, with random runs in order and apart (some of the
/// codes of A, C, G and T), across either end and at negative offsets, with bits above the four
/// of a code, out of order and overlapping, and out of order with the stored bases' codes only
/// (no run is left).
#[test]
fn a_view_of_stored_bytes_reads_as_their_store() {
    let mut state = 5u64;
    let mut next = |below: u64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) % below
    };
    let mut left_none = 0;
    for case in 0..600 {
        let length = 1 + next(300) as usize;
        let packed: Vec<u8> = (0..length.div_ceil(4)).map(|_| next(256) as u8).collect();
        let base = |index: usize| (packed[index >> 2] >> (6 - 2 * (index & 3))) & 3;
        let mut runs = Vec::new();
        let mut offset = next(8) as i64 - 4;
        while offset < length as i64 + 4 {
            let run_length = 1 + next(6) as u32;
            let ncbi4na = match case % 3 {
                0 => 1 << next(4),
                1 => next(16) as u8,
                _ => next(256) as u8,
            };
            runs.push(AmbiguityRun {
                offset: offset as i32 as u32,
                length: run_length,
                ncbi4na,
            });
            offset += i64::from(run_length) + next(12) as i64;
        }
        match case % 4 {
            0 => {}
            1 => runs.reverse(),
            2 => {
                for index in (1..runs.len()).rev() {
                    runs.swap(index, next(index as u64 + 1) as usize);
                }
            }
            _ => {
                // out of order, each run a stored base's own code
                runs = (0..length)
                    .rev()
                    .step_by(3)
                    .map(|index| AmbiguityRun {
                        offset: index as u32,
                        length: 1,
                        ncbi4na: 1 << base(index),
                    })
                    .collect();
            }
        }
        let what = format!("case {case}");
        let store = PackedNucleotide::from_stored(&packed, length, &runs);
        let mut kept = runs.clone();
        x_stored_runs(&packed, length, &mut kept);
        let view = PackedRef::x_stored(&packed, length, &kept);
        assert_eq!(view.len(), length, "{what}");
        assert_eq!(view.ambiguity_runs(), store.ambiguity_runs(), "{what}");
        if case % 4 == 3 && length > 3 {
            assert!(kept.is_empty(), "{what}");
            left_none += 1;
        }
        let masked: Vec<AmbiguityRun> = runs
            .iter()
            .map(|run| AmbiguityRun {
                ncbi4na: run.ncbi4na & 0x0F,
                ..*run
            })
            .collect();
        let letters: Vec<u8> = ncbi_letters(&packed, length, &masked)
            .iter()
            .map(|&code| NCBI4NA_TO_IUPACNA[usize::from(code)])
            .collect();
        let bases: Vec<u8> = (0..length)
            .map(|index| b"ACGT"[usize::from(base(index))])
            .collect();
        for _ in 0..8 {
            let from = next(length as u64 + 1) as usize;
            let to = from + next((length - from) as u64 + 1) as usize;
            let mut out = Vec::new();
            for reader in [view, store.view()] {
                reader.letters_into(from, to, &mut out);
                assert_eq!(out, &letters[from..to], "{what} letters {from}..{to}");
                reader.letters_to(from, &mut out);
                assert_eq!(out, &letters[from..to], "{what} letters_to {from}..{to}");
                reader.bases_into(from, &mut out);
                assert_eq!(out, &bases[from..to], "{what} bases {from}..{to}");
            }
            store.letters_into(from, to, &mut out);
            assert_eq!(out, &letters[from..to], "{what} the store's letters");
        }
    }
    assert!(left_none > 100, "{left_none}");
}
