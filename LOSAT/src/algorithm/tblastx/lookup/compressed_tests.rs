use super::*;
use crate::utils::matrix::ncbistdaa;

// A small deterministic generator (xorshift64) for the tests of the query batches' tables.
struct XorShift(u64);

impl XorShift {
    fn below(&mut self, n: usize) -> usize {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x % n as u64) as usize
    }
}

// Queries of `lengths` letters drawn from `letters`, concatenated as BLASTP's lookup query
// (each query followed by a sentinel), with one lookup segment per query.
fn concatenated_queries(
    rng: &mut XorShift,
    lengths: &[usize],
    letters: &[u8],
) -> (Vec<u8>, Vec<(i32, i32)>) {
    let mut query = Vec::new();
    let mut locations = Vec::new();
    for &length in lengths {
        let left = query.len() as i32;
        for _ in 0..length {
            query.push(letters[rng.below(letters.len())]);
        }
        locations.push((left, left + length as i32 - 1));
        query.push(0);
    }
    (query, locations)
}

// A table of `locations` with no limit on its overflow cells (not finalized).
fn unlimited_table(
    threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
) -> BlastCompressedAaLookupTable {
    let mut table = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
    table.threshold = threshold;
    table.overflow_cell_capacity = usize::MAX;
    let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
    table.add_neighboring_words(&matrix, query, locations);
    table
}

fn assert_same_table(lhs: &BlastCompressedAaLookupTable, rhs: &BlastCompressedAaLookupTable) {
    assert!(lhs.backbone == rhs.backbone, "backbone cells differ");
    assert!(
        lhs.overflow_cells == rhs.overflow_cells,
        "overflow cells differ"
    );
    assert_eq!(lhs.pv, rhs.pv);
    assert_eq!(lhs.pv_array_bts, rhs.pv_array_bts);
    assert_eq!(lhs.longest_chain, rhs.longest_chain);
    assert_eq!(lhs.neighbor_matches, rhs.neighbor_matches);
    assert_eq!(lhs.exact_matches, rhs.exact_matches);
    for index in 0..lhs.backbone.len() {
        if lhs.backbone[index].num_used > 5 {
            assert_eq!(
                lhs.query_offsets_in_scan_order(index),
                rhs.query_offsets_in_scan_order(index)
            );
        }
    }
}

// The queries of the tests below: few letters and a low threshold (11, in NCBI's matrix
// units), so that words repeat and cells overflow.
const FEW_LETTERS: [u8; 4] = [ncbistdaa::A, ncbistdaa::L, ncbistdaa::S, ncbistdaa::K];
const LOW_THRESHOLD: i32 = 1100;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-213
// (a) The table of all batches takes more cells than the capacity, no batch's own table
// does: the table is built, and it is the table without a limit.
#[test]
fn compressed_batches_within_capacity_give_the_unlimited_table() {
    let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
    let (query, locations) =
        concatenated_queries(&mut rng, &[300, 250, 400, 150, 350, 300], &FEW_LETTERS);
    let tables = [0..2, 2..3, 3..6];
    let batch_cells: Vec<usize> = tables
        .iter()
        .map(|table| {
            unlimited_table(LOW_THRESHOLD, &query, &locations[table.clone()])
                .overflow_cells
                .len()
        })
        .collect();
    let capacity = *batch_cells.iter().max().unwrap();
    let mut unlimited = unlimited_table(LOW_THRESHOLD, &query, &locations);
    assert!(capacity > 0);
    assert!(unlimited.overflow_cells.len() > capacity);
    unlimited.finalize();

    let table = try_build_compressed_lookup_with_capacity(
        Blosum62Compressed,
        5,
        LOW_THRESHOLD,
        &query,
        &locations,
        &tables,
        capacity,
    )
    .expect("word size 5")
    .expect("no batch fills the capacity");
    assert_same_table(&table, &unlimited);

    // Locations outside NCBI's tables (a split batch) are in the table, untested.
    let table = try_build_compressed_lookup_with_capacity(
        Blosum62Compressed,
        5,
        LOW_THRESHOLD,
        &query,
        &locations,
        &[0..2, 3..6],
        capacity,
    )
    .expect("word size 5")
    .expect("no batch fills the capacity");
    assert_same_table(&table, &unlimited);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4888-4898
// The compressed score matrix of each search matrix: BLOSUM62's is the default builder's,
// IDENTITY (no frequency ratios) has none, and each other matrix has its own. The table of
// BLOSUM62 as a `ScoringMatrix` is the default builder's table.
#[test]
fn compressed_score_matrix_of_each_search_matrix() {
    for size in [10, 15] {
        let blosum62 = build_blosum62_compressed_score_matrix(size).expect("BLOSUM62");
        assert_eq!(
            build_compressed_score_matrix(ScoringMatrix::Blosum62, size),
            Some(blosum62)
        );
        assert_eq!(
            build_compressed_score_matrix(ScoringMatrix::Identity, size),
            None
        );
        for matrix in ScoringMatrix::ALL {
            if !matches!(matrix, ScoringMatrix::Blosum62 | ScoringMatrix::Identity) {
                let scores = build_compressed_score_matrix(matrix, size).expect("ratios");
                assert_ne!(scores, blosum62, "{matrix:?}");
            }
        }
    }
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    let (query, locations) = concatenated_queries(&mut rng, &[300, 250], &FEW_LETTERS);
    let tables = [0..2];
    let want = try_build_blosum62_compressed_lookup(5, LOW_THRESHOLD, &query, &locations, &tables)
        .expect("word size 5")
        .expect("no batch fills the banks");
    let got = try_build_compressed_lookup(
        ScoringMatrix::Blosum62,
        5,
        LOW_THRESHOLD,
        &query,
        &locations,
        &tables,
    )
    .expect("word size 5")
    .expect("no batch fills the banks");
    assert_same_table(&got, &want);
    // Another matrix's words fill other cells.
    let pam30 = try_build_compressed_lookup(
        ScoringMatrix::Pam30,
        5,
        LOW_THRESHOLD,
        &query,
        &locations,
        &tables,
    )
    .expect("word size 5")
    .expect("no batch fills the banks");
    assert!(
        pam30.backbone != want.backbone,
        "PAM30's cells are BLOSUM62's"
    );
    assert!(try_build_compressed_lookup(
        ScoringMatrix::Identity,
        5,
        LOW_THRESHOLD,
        &query,
        &locations,
        &tables
    )
    .is_none());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
// (b) One batch's own table needs more cells than the capacity: refused, whether the
// batch is the first (its table so far is the whole table) or a later one.
#[test]
fn compressed_batch_past_capacity_is_refused() {
    let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
    let (query, locations) = concatenated_queries(&mut rng, &[100, 400, 60, 80], &FEW_LETTERS);
    let small = [0..1, 2..3, 3..4]
        .iter()
        .map(|table| {
            unlimited_table(LOW_THRESHOLD, &query, &locations[table.clone()])
                .overflow_cells
                .len()
        })
        .max()
        .unwrap();
    let large = unlimited_table(LOW_THRESHOLD, &query, &locations[1..2])
        .overflow_cells
        .len();
    assert!(small > 0 && large > small + 1);
    let build = |tables: &[std::ops::Range<usize>], capacity: usize| {
        try_build_compressed_lookup_with_capacity(
            Blosum62Compressed,
            5,
            LOW_THRESHOLD,
            &query,
            &locations,
            tables,
            capacity,
        )
        .expect("word size 5")
    };
    // A later batch: its own table is built for the test, and the refusal names it.
    assert_eq!(
        build(&[0..1, 1..2, 2..3, 3..4], large - 1).err(),
        Some(BanksFilled(1))
    );
    assert!(build(&[0..1, 1..2, 2..3, 3..4], large).is_ok());
    assert_eq!(build(&[1..2, 2..4], large - 1).err(), Some(BanksFilled(0)));
    // The first batch: the table so far is its own.
    let first = unlimited_table(LOW_THRESHOLD, &query, &locations[0..2])
        .overflow_cells
        .len();
    let rest = unlimited_table(LOW_THRESHOLD, &query, &locations[2..4])
        .overflow_cells
        .len();
    assert!(first >= large && rest < first - 1);
    assert_eq!(build(&[0..2, 2..4], first - 1).err(), Some(BanksFilled(0)));
    assert!(build(&[0..2, 2..4], first).is_ok());
    // The large query in a split batch (no table of its own) is not tested.
    assert!(build(&[0..1, 2..3, 3..4], small).is_ok());
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:825-841
// (c) On random queries: a batch's own table takes cells(n) = 0 (n <= 5) or floor(n / 4)
// cells per backbone cell, the batch takes at least as many cells of the table of all
// batches (the shortcut is sound), and the build refuses exactly when some batch's own
// table needs more cells than the capacity.
#[test]
fn compressed_batch_shortcut_agrees_with_exact_test() {
    let mut rng = XorShift(0x1234_5678_9abc_def1);
    let mut shortcut_undecided = 0usize;
    for round in 0..12 {
        let letters: &[u8] = if round % 3 == 2 {
            &[
                ncbistdaa::A,
                ncbistdaa::L,
                ncbistdaa::S,
                ncbistdaa::K,
                ncbistdaa::E,
                ncbistdaa::G,
                ncbistdaa::X,
            ]
        } else {
            &FEW_LETTERS
        };
        let threshold = [1100, 1600, 2000, 2400][rng.below(4)];
        let lengths: Vec<usize> = (0..2 + rng.below(6)).map(|_| 5 + rng.below(400)).collect();
        let (query, locations) = concatenated_queries(&mut rng, &lengths, letters);
        // Batches of consecutive queries; some are left out (split batches).
        let mut tables = Vec::new();
        let mut start = 0;
        while start < locations.len() {
            let end = (start + 1 + rng.below(3)).min(locations.len());
            if rng.below(4) != 0 {
                tables.push(start..end);
            }
            start = end;
        }
        let mut whole = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
        whole.threshold = threshold;
        whole.overflow_cell_capacity = usize::MAX;
        let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
        let mut next = 0;
        let mut exact = Vec::new();
        for table in &tables {
            whole.add_neighboring_words(&matrix, &query, &locations[next..table.start]);
            let before = whole.overflow_cells.len();
            whole.add_neighboring_words(&matrix, &query, &locations[table.clone()]);
            let shortcut = whole.overflow_cells.len() - before;
            let own = unlimited_table(threshold, &query, &locations[table.clone()]);
            let cells: usize = own
                .backbone
                .iter()
                .map(|cell| match cell.num_used {
                    n if n <= 5 => 0,
                    n => n as usize / 4,
                })
                .sum();
            assert_eq!(own.overflow_cells.len(), cells);
            assert!(shortcut >= cells, "round {round}: {shortcut} < {cells}");
            exact.push((cells, shortcut));
            next = table.end;
        }
        whole.add_neighboring_words(&matrix, &query, &locations[next..]);
        whole.finalize();

        let mut capacities = vec![1usize, whole.overflow_cells.len() + 1];
        for &(cells, shortcut) in &exact {
            capacities.extend([cells, cells + 1, shortcut, shortcut + 1]);
            if cells > 1 {
                capacities.push(cells - 1);
            }
        }
        capacities.retain(|&capacity| capacity > 0);
        for capacity in capacities {
            // The first table whose own words fill the banks.
            let refused = exact.iter().position(|&(cells, _)| cells > capacity);
            shortcut_undecided += exact
                .iter()
                .filter(|&&(cells, shortcut)| cells <= capacity && shortcut > capacity)
                .count();
            let built = try_build_compressed_lookup_with_capacity(
                Blosum62Compressed,
                5,
                threshold,
                &query,
                &locations,
                &tables,
                capacity,
            )
            .expect("word size 5");
            match built {
                Err(filled) => assert_eq!(
                    Some(filled),
                    refused.map(BanksFilled),
                    "round {round}, capacity {capacity}"
                ),
                Ok(table) => {
                    assert_eq!(refused, None, "round {round}, capacity {capacity}");
                    assert_same_table(&table, &whole);
                }
            }
        }
    }
    // Batches that the shortcut leaves to the exact test, and that pass it, were seen.
    assert!(shortcut_undecided > 0);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1071-1110
// (d) `CompressedWordCount` gives the number of words that `add_word_hits` adds for each
// offset, on random words of all NCBISTDAA letters (some without a compressed letter), with
// and without a threshold, and the count of a word does not depend on the order of its
// letters (the counts are kept per set of letters); at each word size of the table (6 and 7
// with shorter queries and higher thresholds, whose tables are larger), for BLOSUM62 and two
// other matrices (`CompressedWordCount::new_of`, BLASTP's block search with them, SC2b).
#[test]
fn compressed_word_count_is_the_words_of_each_offset() {
    for search_matrix in [
        ScoringMatrix::Blosum62,
        ScoringMatrix::Pam30,
        ScoringMatrix::Blosum45,
    ] {
        compressed_word_count_is_the_words_of_each_offset_of(search_matrix);
    }
}

fn compressed_word_count_is_the_words_of_each_offset_of(search_matrix: ScoringMatrix) {
    let mut rng = XorShift(0x0bad_cafe_1234_5678);
    let standard = (*b"ACDEFGHIKLMNPQRSTVWY").map(aa_char_to_ncbistdaa);
    let mut letters: Vec<u8> = (1..BLASTAA_SIZE as u8).collect();
    for _ in 0..3 {
        letters.extend_from_slice(&standard);
    }
    for (word_size, query_length, thresholds) in [
        (5, 40, &[0, 1100, 1600, 1930, 2000, 2400][..]),
        (6, 16, &[1600, 2100, 2400][..]),
        (7, 16, &[1600, 2025, 2400][..]),
    ] {
        let (alphabet, _) = compressed_lookup_parameters(word_size).expect("word size");
        let matrix =
            build_compressed_score_matrix(search_matrix, alphabet).expect("compressed matrix");
        let info = CompressedNeighborInfo::new(&matrix, alphabet);
        for &threshold in thresholds {
            let counter = CompressedWordCount::new_of(search_matrix, word_size, threshold)
                .expect("word size");
            let mut table =
                BlastCompressedAaLookupTable::new(word_size, 0.0).expect("compressed lookup");
            table.threshold = threshold;
            table.overflow_cell_capacity = usize::MAX;
            let query: Vec<u8> = (0..query_length)
                .map(|_| letters[rng.below(letters.len())])
                .collect();
            let mut words_before = 0u64;
            for offset in 0..=query.len() - word_size {
                table.add_word_hits(&info, &matrix, &query, offset as i32);
                let words_after: u64 = table.backbone.iter().map(|cell| cell.num_used as u64).sum();
                let added = words_after - words_before;
                words_before = words_after;
                let word = &query[offset..offset + word_size];
                assert_eq!(
                    counter.words(word),
                    added,
                    "{search_matrix:?} word size {word_size}, threshold {threshold}, offset {offset}"
                );
                if word.iter().all(|&residue| {
                    usize::from(table.compress_table[usize::from(residue)]) < alphabet
                }) {
                    let mut reversed = word.to_vec();
                    reversed.reverse();
                    assert_eq!(
                        counter.count(&reversed),
                        added,
                        "word size {word_size}, reversed, offset {offset}"
                    );
                }
            }
        }
    }
}

// (d2) Each set of `word_size` letters that have a compressed letter has its own rank
// (`CompressedWordCount::rank`), below the kept counts' length, whatever the order of the
// word's letters; a word with a letter that has none has no rank; at each word size of the
// table.
#[test]
fn compressed_word_count_ranks_each_set_of_letters_once() {
    for word_size in 5..=7 {
        let counter = CompressedWordCount::new(word_size, 1100).expect("word size");
        let letters: Vec<u8> = (0..BLASTAA_SIZE as u8)
            .filter(|&residue| counter.letter_index[usize::from(residue)] != u8::MAX)
            .collect();
        assert_eq!(letters.len(), 23);
        let mut ranked = vec![false; counter.counts.len()];
        // Each multiset of letters as its sorted indices, in lexicographic order.
        let mut set = vec![0usize; word_size];
        loop {
            let word: Vec<u8> = set.iter().map(|&i| letters[i]).collect();
            let rank = counter.rank(&word).expect("a rank");
            assert!(!ranked[rank], "word size {word_size}: {word:?}");
            ranked[rank] = true;
            let mut reordered = word.clone();
            reordered.rotate_left(word_size / 2);
            reordered.swap(0, word_size - 1);
            assert_eq!(counter.rank(&reordered), Some(rank));
            let Some(last) = set.iter().rposition(|&i| i + 1 < letters.len()) else {
                break;
            };
            let next = set[last] + 1;
            set[last..].fill(next);
        }
        assert!(ranked.iter().all(|&r| r), "word size {word_size}");
        for residue in 0..BLASTAA_SIZE as u8 {
            if counter.letter_index[usize::from(residue)] == u8::MAX {
                let mut word = letters[..word_size].to_vec();
                word[1] = residue;
                assert_eq!(counter.rank(&word), None);
                assert_eq!(counter.words(&word), 0);
            }
        }
    }
}

// (e) The count clears a table exactly when a quarter of its words fits the capacity, and a
// table it clears is not refused (`try_build_compressed_lookup_with_capacity`).
#[test]
fn compressed_word_count_clears_only_tables_within_the_capacity() {
    let mut rng = XorShift(0x5851_f42d_4c95_7f2d);
    let (query, locations) = concatenated_queries(&mut rng, &[120, 300, 80, 200], &FEW_LETTERS);
    let queries: Vec<&[u8]> = locations
        .iter()
        .map(|&(left, right)| &query[left as usize..=right as usize])
        .collect();
    let counter = CompressedWordCount::new(5, LOW_THRESHOLD).expect("word size 5");
    let words: u64 = queries
        .iter()
        .flat_map(|query| query.windows(5))
        .map(|word| counter.words(word))
        .sum();
    let table = unlimited_table(LOW_THRESHOLD, &query, &locations);
    let table_words: u64 = table.backbone.iter().map(|cell| cell.num_used as u64).sum();
    assert_eq!(words, table_words);
    let cells = table.overflow_cells.len();
    let cleared = words.div_ceil(4) as usize;
    assert!(cells > 0 && cells <= cleared);
    assert_eq!(
        counter.words_within(&mut queries.iter().copied(), cleared),
        Some(words)
    );
    assert_eq!(
        counter.words_within(&mut queries.iter().copied(), cleared - 1),
        None
    );
    let table = try_build_compressed_lookup_with_capacity(
        Blosum62Compressed,
        5,
        LOW_THRESHOLD,
        &query,
        &locations,
        &[0..locations.len()],
        cleared,
    )
    .expect("word size 5");
    assert!(table.is_ok());
    // The public count uses NCBI's capacity.
    assert_eq!(
        counter.batch_words(&mut queries.iter().copied()),
        Some(words)
    );
    // NCBI's backbone of word size 5 (15 letters) and 7 (10 letters).
    assert_eq!(counter.backbone_cells(), 15u64.pow(5));
    assert_eq!(
        CompressedWordCount::new(7, LOW_THRESHOLD)
            .expect("word size 7")
            .backbone_cells(),
        10u64.pow(7)
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:400-405
// (f) BLASTP's block search (EXPERIMENT LOSAT_X_BPQBLOCK, `blastp/x_query_blocks.rs`): the
// table of some queries (a block) gives, at each subject offset, the hits of those queries
// that the table of all queries gives, not always in the same order (cells of more than 5
// offsets); and a scan gives the same hits in the same order with any offset-pair array
// that takes its longest cell, wherever it stops and resumes.
#[test]
fn compressed_block_table_scans_the_whole_tables_hits_per_subject_offset() {
    let mut rng = XorShift(0x6a09_e667_f3bc_c908);
    let (query, locations) = concatenated_queries(&mut rng, &[150, 90, 200, 60, 120], &FEW_LETTERS);
    let mut whole = unlimited_table(LOW_THRESHOLD, &query, &locations);
    whole.finalize();
    let block = 1..3;
    let mut part = unlimited_table(LOW_THRESHOLD, &query, &locations[block.clone()]);
    part.finalize();
    let block_offsets = locations[block.start].0..=locations[block.end - 1].1;
    let mut subject: Vec<u8> = (0..400)
        .map(|_| FEW_LETTERS[rng.below(FEW_LETTERS.len())])
        .collect();
    subject.push(0);
    let scan = |table: &BlastCompressedAaLookupTable, array_size: i32| {
        let seq_ranges = [(0, subject.len() as i32 - 1)];
        let mut s_range = [0, 0, seq_ranges[0].1 - 5];
        let mut offset_pairs = vec![OffsetPair::default(); array_size as usize];
        let mut hits = Vec::new();
        let mut stops = 0;
        while s_range[1] <= s_range[2] {
            let count = table.scan_subject(
                &subject,
                &seq_ranges,
                &mut offset_pairs,
                array_size,
                &mut s_range,
            );
            hits.extend(
                offset_pairs[..count as usize]
                    .iter()
                    .map(|pair| (pair.s_off, pair.q_off as i32)),
            );
            stops += 1;
        }
        (hits, stops)
    };
    let (whole_hits, _) = scan(&whole, 4096 + whole.longest_chain);
    let (part_hits, _) = scan(&part, 4096 + part.longest_chain);
    let whole_in_block: Vec<(u32, i32)> = whole_hits
        .iter()
        .copied()
        .filter(|(_, q_off)| block_offsets.contains(q_off))
        .collect();
    assert_ne!(
        part_hits, whole_in_block,
        "some subject offset's order differs"
    );
    let mut sorted_part = part_hits.clone();
    sorted_part.sort_unstable();
    let mut sorted_whole = whole_in_block.clone();
    sorted_whole.sort_unstable();
    assert_eq!(sorted_part, sorted_whole);
    assert!(part_hits.windows(2).all(|pair| pair[0].0 <= pair[1].0));
    for array_size in [
        part.longest_chain,
        part.longest_chain + 7,
        3 * part.longest_chain,
    ] {
        let (hits, stops) = scan(&part, array_size);
        assert!(stops > 1, "array {array_size}: the scan stops and resumes");
        assert_eq!(hits, part_hits, "array {array_size}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:782-849
// The cast to `Int4` of `opt->threshold` times `kMatrixScale`, when past `Int4`, is `INT_MIN`
// (every word is a neighbor); words past the overflow banks make the table `Err` (NCBI writes past
// its banks).
#[test]
fn a_threshold_past_int4_takes_every_word_up_to_the_banks() {
    assert_eq!(
        BlastCompressedAaLookupTable::new(5, 19.3)
            .unwrap()
            .threshold,
        1930
    );
    for threshold in [21_474_836.48, 1e10, f64::INFINITY] {
        let lookup = BlastCompressedAaLookupTable::new(5, threshold).unwrap();
        assert_eq!(lookup.threshold, i32::MIN, "{threshold}");
    }
    let query = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0];
    let words = 15usize.pow(5);
    let lookup = build_blosum62_compressed_lookup(5, 1e10, &query, &[(1, 6)])
        .unwrap()
        .expect("two offsets take no overflow cell");
    let offsets: i32 = lookup.backbone.iter().map(|cell| cell.num_used).sum();
    assert_eq!(offsets as usize, 2 * words);
    assert!(lookup.overflow_cells.is_empty());
    // Seven offsets take one overflow cell in each backbone cell.
    assert!(
        build_blosum62_compressed_lookup_with_capacity(5, 1e10, &query, &[(1, 11)], 1000)
            .unwrap()
            .is_err()
    );
}

#[test]
fn compressed_word_hit_overflow_matches_ncbi_scan_order() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    for offset in 10..18 {
        lookup.add_word_hit(42, offset);
    }

    assert_eq!(lookup.backbone[42].num_used, 8);
    assert_eq!(lookup.overflow_cells.len(), 2);
    assert_eq!(
        lookup.query_offsets_in_scan_order(42),
        vec![10, 11, 12, 17, 13, 14, 15, 16]
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1182-1242
#[test]
fn compressed_finalize_sets_pv_and_longest_chain() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    lookup.add_word_hit(7, 100);
    lookup.add_word_hit(33, 200);
    lookup.add_word_hit(33, 201);
    lookup.finalize();

    assert_eq!(lookup.longest_chain, 2);
    assert!(compressed_pv_test(
        &lookup.pv,
        7,
        lookup.pv_array_bts as usize
    ));
    assert!(compressed_pv_test(
        &lookup.pv,
        33,
        lookup.pv_array_bts as usize
    ));
    assert!(!compressed_pv_test(
        &lookup.pv,
        34,
        lookup.pv_array_bts as usize
    ));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:859-929
#[test]
fn compressed_encoded_and_unencoded_indices_match() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    let word = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
    ];
    let encoded_word: Vec<u8> = word
        .iter()
        .map(|&residue| lookup.compress_table[residue as usize])
        .collect();
    let encoded_index = compute_encoded_compressed_index(5, &encoded_word, 15);
    let (unencoded_index, skip) =
        compute_compressed_index(5, &word, 15, &lookup.scaled_compress_table);

    assert_eq!(skip, 0);
    assert_eq!(encoded_index as i32, unencoded_index);
    lookup.add_unencoded(&word, 321);
    assert_eq!(lookup.query_offsets_in_scan_order(encoded_index), vec![321]);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:322-384
#[test]
fn compressed_scan_subject_emits_ncbi_subject_major_order() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    let word = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
    ];
    lookup.add_unencoded(&word, 10);
    lookup.add_unencoded(&word, 11);
    lookup.finalize();

    let subject = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
        0,
    ];
    let seq_ranges = [(0, 10)];
    let mut s_range = [0, 0, 5];
    let mut offset_pairs = vec![OffsetPair::default(); 8];
    let hits = lookup.scan_subject(&subject, &seq_ranges, &mut offset_pairs, 8, &mut s_range);

    assert_eq!(hits, 4);
    assert_eq!(
        offset_pairs[..hits as usize]
            .iter()
            .map(|pair| (pair.q_off, pair.s_off))
            .collect::<Vec<_>>(),
        vec![(10, 0), (11, 0), (10, 5), (11, 5)]
    );
    assert_eq!(s_range[1], 6);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:373-378
#[test]
fn compressed_scan_subject_small_buffer_resumes_without_splitting_subject_offset() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    let word = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
    ];
    for query_offset in 10..16 {
        lookup.add_unencoded(&word, query_offset);
    }
    lookup.finalize();

    let subject = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
        ncbistdaa::X,
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
        0,
    ];
    let seq_ranges = [(0, 11)];
    let mut s_range = [0, 0, 6];
    let mut offset_pairs = vec![OffsetPair::default(); 6];
    let mut observed = Vec::new();

    while s_range[1] <= s_range[2] {
        let hits = lookup.scan_subject(&subject, &seq_ranges, &mut offset_pairs, 6, &mut s_range);
        assert_eq!(hits, 6);
        observed.extend(
            offset_pairs[..hits as usize]
                .iter()
                .map(|pair| (pair.q_off, pair.s_off)),
        );
    }

    assert_eq!(
        observed,
        vec![
            (10, 0),
            (11, 0),
            (12, 0),
            (13, 0),
            (14, 0),
            (15, 0),
            (10, 6),
            (11, 6),
            (12, 6),
            (13, 6),
            (14, 6),
            (15, 6),
        ]
    );
}
