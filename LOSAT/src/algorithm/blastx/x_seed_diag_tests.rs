use super::*;
use crate::algorithm::blastp::x_seed_fused::XDiagMode;
use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

/// The cell tests this machine can run (the masked AVX2 batch, with LOSAT_X_BXCROSS its
/// batches, their calls deferred with LOSAT_X_BXCROSSDEFER and not, and LOSAT_X_BXCROSS128's
/// 128-bit batches, which `2` runs where the AVX2 ones would), each with the branch-free scan
/// (BXSCANBF, graduated).
fn modes() -> Vec<XCellTest> {
    let mut modes = vec![
        XCellTest::Diag(XDiagMode::PerHit),
        XCellTest::Diag(XDiagMode::Scalar),
    ];
    for (cross, defer, cross128) in [
        (false, false, 0),
        (true, false, 1),
        (true, true, 1),
        (true, false, 2),
        (true, true, 2),
    ] {
        let auto = x_cell_test("", cross, defer, cross128);
        if auto != XCellTest::Diag(XDiagMode::PerHit) && !modes.contains(&auto) {
            modes.push(auto);
        }
    }
    modes
}

/// A search of `query` (nucleotide FASTA) with `-window_size window` (default options
/// otherwise), its statistics set for subjects of `subject_length` residues; its standard table
/// has `Int4` hits, or with `narrow` 16-bit hits.
fn search(
    query: &[u8],
    window: i32,
    subject_length: usize,
    narrow: bool,
) -> (
    PreparedQueryBatch,
    Vec<ContextParameters>,
    ResolvedOptions,
    BxLookup,
) {
    search_with(query, window, subject_length, narrow, "BLOSUM62", &[])
}

/// [`search`] with `-matrix matrix` (its default gap costs and threshold, its score block) and
/// the options `extra`; another matrix than BLOSUM62 takes its build (`Lookup::x_new_of`, item
/// SC6b).
fn search_with(
    query: &[u8],
    window: i32,
    subject_length: usize,
    narrow: bool,
    matrix: &str,
    extra: &[&str],
) -> (
    PreparedQueryBatch,
    Vec<ContextParameters>,
    ResolvedOptions,
    BxLookup,
) {
    let q = parse_fasta(query, false, false).unwrap();
    let window = window.to_string();
    let mut args = vec![
        "losat",
        "blastx",
        "-query",
        "query",
        "-subject",
        "subject",
        "-window_size",
        &window,
        "-matrix",
        matrix,
    ];
    args.extend_from_slice(extra);
    let crate::cli::Commands::Blastx(args) = crate::cli::try_parse_from(args).unwrap().command
    else {
        unreachable!()
    };
    let options = args.resolve().unwrap();
    let batch = prepare_queries(&q, &options).unwrap();
    let mut params = parameters::score_block_for(&batch, options.scoring);
    parameters::effective_lengths(
        &batch,
        &mut params,
        &options,
        subject_length as i64,
        1,
        None,
    )
    .unwrap();
    let lookup = Lookup::x_new_for(&batch, &params, &options, 1, narrow.then_some(0)).unwrap();
    assert_eq!(matches!(lookup, BxLookup::Narrow(_)), narrow);
    parameters::subject_parameters(
        &batch,
        &mut params,
        &options,
        parameters::SubjectNumbers::one(subject_length),
        None,
    )
    .unwrap();
    (batch, params, options, lookup)
}

/// The subject records of a protein FASTA, encoded with the sentinel after the last residue.
fn subjects(fasta: &[u8]) -> Vec<Vec<u8>> {
    parse_fasta(fasta, true, false)
        .unwrap()
        .iter()
        .map(|record| {
            let mut subject = record
                .sequence
                .iter()
                .copied()
                .map(crate::utils::matrix::aa_char_to_ncbistdaa)
                .collect::<Vec<_>>();
            subject.push(0);
            subject
        })
        .collect()
}

/// Runs every subject through the reference and through the fused loop (cell test `mode`;
/// the branch-free scan), each on its own table starting from
/// `start(table)`, and compares the hit lists and the tables after each subject; returns the
/// hits and the counts of the fused loop.
fn compare(
    search: &(
        PreparedQueryBatch,
        Vec<ContextParameters>,
        ResolvedOptions,
        BxLookup,
    ),
    subjects: &[Vec<u8>],
    mode: XCellTest,
    start: impl Fn(&mut Diagonals),
) -> (usize, [u64; 4]) {
    match &search.3 {
        BxLookup::Wide(Lookup::Standard(table, kept)) => {
            compare_on(search, table, kept, subjects, mode, start)
        }
        BxLookup::Narrow(NarrowLookup { table, kept }) => {
            compare_on(search, table, kept, subjects, mode, start)
        }
        BxLookup::Wide(Lookup::Compressed(_)) => panic!("the standard lookup table"),
    }
}

/// [`compare`] for the standard table `table` of the search's lookup, of either width.
fn compare_on<K: XCell>(
    (batch, params, options, lookup): &(
        PreparedQueryBatch,
        Vec<ContextParameters>,
        ResolvedOptions,
        BxLookup,
    ),
    table: &BlastAaLookupTable<K>,
    kept: &XKeptPresence,
    subjects: &[Vec<u8>],
    mode: XCellTest,
    start: impl Fn(&mut Diagonals),
) -> (usize, [u64; 4]) {
    let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
    let mut reference = Diagonals::new(length, options.window_size).unwrap();
    let mut fused = Diagonals::new(length, options.window_size).unwrap();
    // EXPERIMENT (LOSAT_X_BXUNGAPSENT): the instances with the sentinel walks, without and
    // with the shadow's comparison, each on a table of its own.
    let mut sent = [
        Diagonals::new(length, options.window_size).unwrap(),
        Diagonals::new(length, options.window_size).unwrap(),
    ];
    start(&mut reference);
    start(&mut fused);
    sent.iter_mut().for_each(&start);
    let (mut hits, mut counts) = (0, [0u64; 4]);
    let matrix = options.scoring;
    for (i, subject) in subjects.iter().enumerate() {
        // BLOSUM62's reference and fused loop, or the instances of another matrix (item SC6b).
        let (expected, (actual, c)) = if matrix == ScoringMatrix::Blosum62 {
            (
                word_finder_bx(
                    batch,
                    params,
                    options,
                    subject,
                    lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap(),
                x_word_finder_fused(
                    Blosum62Bx, mode, batch, params, options, subject, table, kept, &mut fused,
                )
                .unwrap(),
            )
        } else {
            let scoring = MatrixBx::of(matrix);
            (
                word_finder_bx_of(
                    scoring,
                    batch,
                    params,
                    options,
                    subject,
                    lookup,
                    &mut reference,
                    |_| {},
                )
                .unwrap(),
                x_word_finder_fused(
                    scoring, mode, batch, params, options, subject, table, kept, &mut fused,
                )
                .unwrap(),
            )
        };
        let run = format!("{matrix:?} {mode:?} subject {i}");
        assert_eq!(actual, expected, "{run}");
        assert_eq!(fused.offset, reference.offset, "{run}");
        assert!(fused.entries == reference.entries, "{run}");
        let lead: Vec<u8> = std::iter::once(0).chain(subject.iter().copied()).collect();
        for (table_sent, shadow) in sent.iter_mut().zip([false, true]) {
            let (got, c2) = if matrix == ScoringMatrix::Blosum62 {
                let x = XSent::for_test(Blosum62Bx, &batch.sequence_start, &lead, shadow);
                by_word_size!(
                    options.word_size,
                    W => x_word_finder_fused_sent::<_, _, W>(
                        mode, batch, params, options, subject, table, kept, table_sent, x,
                    )
                )
            } else {
                let x = XSent::for_test(MatrixBx::of(matrix), &batch.sequence_start, &lead, shadow);
                by_word_size!(
                    options.word_size,
                    W => x_word_finder_fused_sent::<_, _, W>(
                        mode, batch, params, options, subject, table, kept, table_sent, x,
                    )
                )
            }
            .unwrap();
            assert_eq!(got, expected, "{run} sent shadow={shadow}");
            assert_eq!(c2, c, "{run} sent shadow={shadow}");
            assert_eq!(table_sent.offset, reference.offset, "{run} sent");
            assert!(table_sent.entries == reference.entries, "{run} sent");
        }
        hits += actual.len();
        for (total, n) in counts.iter_mut().zip(c) {
            *total += n;
        }
    }
    (hits, counts)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
// The session F fixture (the frozen NCBI INIT rows of tests/blastx_f_diagonal) at its two-hit
// windows, from a fresh table, a table near INT4_MAX / 4, a reset one and one whose offset is
// just below 2^30 (the 31-bit `last_hit` wraps): the fused loop gives the reference's hits and
// table with every cell test.
#[test]
fn the_fused_loop_matches_the_reference_on_the_session_f_fixture() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
    let query = std::fs::read(dir.join("query.fna")).unwrap();
    let subject = subjects(&std::fs::read(dir.join("subject.faa")).unwrap());
    let length = subject[0].len() - 1;
    for (window, narrow) in [1, 39, 40, 41]
        .into_iter()
        .flat_map(|window| [(window, false), (window, true)])
    {
        let s = search(&query, window, length, narrow);
        let flagged = |last_hit: i32| {
            let mut entry = DiagStruct::default();
            entry.set_last_hit(last_hit);
            entry.set_flag(1);
            entry
        };
        for mode in modes() {
            let repeated = [subject[0].clone(), subject[0].clone(), subject[0].clone()];
            compare(&s, &repeated, mode, |_| {});
            compare(&s, &repeated, mode, |d| {
                d.offset = 536_870_910;
                d.entries.fill(flagged(123));
            });
            compare(&s, &repeated, mode, |d| {
                d.offset = 536_870_911;
                d.entries.fill(flagged(536_870_900));
                d.finish_subject(83).unwrap();
            });
            compare(&s, &repeated, mode, |d| {
                d.offset = (1 << 30) - 64;
                d.entries.fill(DiagStruct::clear(window));
            });
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:211-242
// The context table of a batch gives the binary search's index at every query offset, and a
// second batch (other queries) gets its own table.
#[test]
fn the_context_table_gives_the_binary_search_index_for_every_offset() {
    let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let genome = parse_fasta(
        &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
        false,
        false,
    )
    .unwrap();
    let letters = &genome[0].sequence;
    for pieces in [[1_500, 2_200, 900, 3_100], [4_000, 333, 2_500, 1_200]] {
        let mut query = Vec::new();
        for (i, &n) in pieces.iter().enumerate() {
            query.extend_from_slice(format!(">p{i}\n").as_bytes());
            query.extend_from_slice(&letters[i * 50_000..i * 50_000 + n]);
            query.push(b'\n');
        }
        let (batch, _, _, lookup) = search(&query, 40, 500, false);
        let BxLookup::Wide(Lookup::Standard(_, kept)) = &lookup else {
            panic!("a standard table");
        };
        let span = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
        assert!(batch.contexts.len() >= 24);
        x_with_context_table(&batch, kept, |table| {
            for q in 0..span as i32 {
                let expected = batch
                    .contexts
                    .partition_point(|c| c.offset as i32 <= q)
                    .saturating_sub(1);
                assert_eq!(table.get(q), Some(expected), "q={q}");
            }
            assert_eq!(table.get(span as i32), None);
        });
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-552
// 6 kb of the AP027131 genome against the first 60 of its proteins, one table for all of
// them: many lookup cells with eight hits or more (so the batches run, not only the scalar
// rest of a cell), extensions with and without a right extension, and the entries a subject
// leaves for the next. Windows 40 (the default) and 400, and a run starting just below 2^30.
#[test]
fn the_fused_loop_matches_the_reference_on_many_cells_and_subjects() {
    let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let genome = parse_fasta(
        &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
        false,
        false,
    )
    .unwrap();
    let piece = String::from_utf8(genome[0].sequence[200_000..206_000].to_vec()).unwrap();
    let query = format!(">piece\n{piece}\n").into_bytes();
    let all = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
    let subjects = &all[..60];
    let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
    for (window, narrow) in [(40, false), (40, true), (400, false), (400, true)] {
        let s = search(&query, window, mean, narrow);
        // The lookup has cells of a full batch (eight hits) and more.
        let cells: Vec<usize> = match &s.3 {
            BxLookup::Wide(Lookup::Standard(table, _)) => {
                table.backbone.iter().map(|c| c.count()).collect()
            }
            BxLookup::Narrow(NarrowLookup { table, .. }) => {
                table.backbone.iter().map(|c| c.count()).collect()
            }
            BxLookup::Wide(Lookup::Compressed(_)) => panic!("the standard lookup table"),
        };
        assert!(cells.iter().any(|&n| n >= 16));
        for mode in modes() {
            let (hits, counts) = compare(&s, subjects, mode, |_| {});
            assert!(hits > 0, "window={window} {mode:?}: no hit");
            if mode == XCellTest::Diag(XDiagMode::PerHit) {
                assert_eq!(counts, [0, 0, 0, 0]);
            } else {
                assert!(counts[0] > 0 && counts[2] > 0, "{mode:?} {counts:?}");
                assert!(counts[1] > counts[0], "{mode:?} {counts:?}");
            }
            compare(&s, subjects, mode, |d| {
                d.offset = (1 << 30) - 3_000;
                d.entries.fill(DiagStruct::clear(window));
            });
        }
    }
}

/// The 6 kb of AP027131 as a query and its first 40 proteins as subjects (with their mean
/// length) of the tests of the other matrices.
fn other_matrix_inputs() -> (Vec<u8>, Vec<Vec<u8>>, usize) {
    let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let genome = parse_fasta(
        &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
        false,
        false,
    )
    .unwrap();
    let piece = String::from_utf8(genome[0].sequence[200_000..206_000].to_vec()).unwrap();
    let query = format!(">piece\n{piece}\n").into_bytes();
    let mut subjects = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
    subjects.truncate(40);
    let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
    (query, subjects, mean)
}

/// BLASTX's other matrices, gapped, and IDENTITY, which an ungapped search takes with word size
/// 3 (its byte table has keys for the sentinel walks).
const OTHER_MATRICES: [(&str, &[&str]); 8] = [
    ("BLOSUM45", &[]),
    ("BLOSUM50", &[]),
    ("BLOSUM80", &[]),
    ("BLOSUM90", &[]),
    ("PAM30", &[]),
    ("PAM70", &[]),
    ("PAM250", &[]),
    ("IDENTITY", &["-ungapped", "-comp_based_stats", "0"]),
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-606 (as above), with
// `matrix->data` of the search's matrix (`s_BlastAaExtendTwoHit`).
// The test above with every other matrix of BLASTX's table and IDENTITY (item SC6b): their
// tables of either width, the fused loop with their walks (`bx_two_hit`) and the sentinel walks
// with their keys give the reference's hits and diagonal table.
#[test]
fn the_fused_loop_matches_the_reference_with_other_matrices() {
    let (query, subjects, mean) = other_matrix_inputs();
    for (matrix, extra) in OTHER_MATRICES {
        for narrow in [false, true] {
            let s = search_with(&query, 40, mean, narrow, matrix, extra);
            assert_ne!(s.2.scoring, ScoringMatrix::Blosum62);
            for mode in modes() {
                let (hits, counts) = compare(&s, &subjects, mode, |_| {});
                assert!(hits > 0, "{matrix} {mode:?}: no hit");
                if mode != XCellTest::Diag(XDiagMode::PerHit) {
                    assert!(counts[2] > 0, "{matrix} {mode:?} {counts:?}");
                }
            }
        }
    }
}

// No NCBI counterpart (speed only): the table of 16-bit hits is the table of 32-bit hits of the
// same query in every property that chooses the loop (its hits per used cell are the build's
// choice of the width, `x_narrow_cutoff`), and the entry point of the seed stage scans it, with
// the fused loop or (where that does not act) the reference scan, to the reference's hits and
// diagonal table.
#[test]
fn a_table_of_16_bit_hits_is_scanned_to_the_reference_hits_by_the_seed_stage() {
    let fasta = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let genome = parse_fasta(
        &std::fs::read(fasta.join("AP027131.fasta")).unwrap(),
        false,
        false,
    )
    .unwrap();
    let piece = String::from_utf8(genome[0].sequence[200_000..206_000].to_vec()).unwrap();
    let query = format!(">piece\n{piece}\n").into_bytes();
    let all = subjects(&std::fs::read(fasta.join("AP027131.faa")).unwrap());
    let subjects = &all[..40];
    let mean = subjects.iter().map(|s| s.len() - 1).sum::<usize>() / subjects.len();
    let wide = search(&query, 40, mean, false);
    let narrow = search(&query, 40, mean, true);
    let (
        BxLookup::Wide(Lookup::Standard(wide_table, _)),
        BxLookup::Narrow(NarrowLookup {
            table: narrow_table,
            ..
        }),
    ) = (&wide.3, &narrow.3)
    else {
        panic!("the tables of the two widths");
    };
    let ratio = x_hits_per_cell(wide_table);
    assert!(ratio > 0);
    assert_eq!(x_hits_per_cell(narrow_table), ratio);
    // The build's choice at the table's own ratio and one hundredth above it.
    for (cutoff, narrower) in [(ratio, true), (ratio + 1, false)] {
        let (batch, params, options, _) = &narrow;
        let query_batch = Lookup::x_new_as(batch, params, options, 1, Some(cutoff)).unwrap();
        assert_eq!(matches!(query_batch, BxLookup::Narrow(_)), narrower);
    }
    let (batch, params, options, lookup) = &narrow;
    let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
    let mut reference = Diagonals::new(length, options.window_size).unwrap();
    let mut scanned = Diagonals::new(length, options.window_size).unwrap();
    let mut hits = 0;
    for subject in subjects {
        let expected = word_finder_bx(
            batch,
            params,
            options,
            subject,
            lookup,
            &mut reference,
            |_| {},
        )
        .unwrap();
        let actual = crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
            Blosum62Bx,
            batch,
            params,
            options,
            subject,
            lookup,
            &mut scanned,
            None,
        )
        .unwrap();
        assert_eq!(actual, expected);
        assert!(scanned.offset == reference.offset && scanned.entries == reference.entries);
        hits += actual.len();
    }
    assert!(hits > 0);
}

// No NCBI counterpart (speed only): the entry point of another matrix's seed stage
// (`x_word_finder_scan` with `MatrixBx`, item SC6b) scans its tables as the search builds them
// (16-bit hits where `x_narrow_cutoff` lets the fused loop take them), with and without the
// sentinel walks, and with `-window_size 0` (the one-hit stage, the direct loop), to the
// reference's hits and diagonal table.
#[test]
fn another_matrix_s_seed_stage_gives_the_reference_hits() {
    let (query, subjects, mean) = other_matrix_inputs();
    let mut widths = Vec::new();
    for (matrix, extra) in OTHER_MATRICES {
        for window in [40, 0] {
            let (batch, params, options, _) =
                search_with(&query, window, mean, false, matrix, extra);
            let lookup =
                Lookup::x_new_of(&batch, &params, &options, 1, x_narrow_cutoff(&options)).unwrap();
            if window > 0 {
                widths.push((matrix, matches!(lookup, BxLookup::Narrow(_))));
            } else {
                assert!(
                    matches!(lookup, BxLookup::Wide(_)),
                    "{matrix}: one-hit table"
                );
            }
            let scoring = MatrixBx::of(options.scoring);
            let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
            // The sentinel walks only with a two-hit window (`query_fits_of`).
            for with_sent in [false, window > 0] {
                let mut reference = Diagonals::new(length, options.window_size).unwrap();
                let mut scanned = Diagonals::new(length, options.window_size).unwrap();
                let mut hits = 0;
                for subject in &subjects {
                    let expected = word_finder_bx_of(
                        scoring,
                        &batch,
                        &params,
                        &options,
                        subject,
                        &lookup,
                        &mut reference,
                        |_| {},
                    )
                    .unwrap();
                    let lead: Vec<u8> = std::iter::once(0).chain(subject.iter().copied()).collect();
                    let sent = with_sent
                        .then(|| XSent::for_test(scoring, &batch.sequence_start, &lead, false));
                    let actual = crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
                        scoring,
                        &batch,
                        &params,
                        &options,
                        subject,
                        &lookup,
                        &mut scanned,
                        sent,
                    )
                    .unwrap();
                    let run = format!("{matrix} window={window} sent={with_sent}");
                    assert_eq!(actual, expected, "{run}");
                    assert!(
                        scanned.offset == reference.offset && scanned.entries == reference.entries,
                        "{run}"
                    );
                    hits += actual.len();
                }
                assert!(hits > 0, "{matrix} window={window}");
            }
        }
    }
    // Where the fused loop acts (`x_bxdiag_mode`: AVX2 and POPCNT, or a value of
    // LOSAT_X_BXDIAG), some of these tables take 16-bit hits and the others the direct loop.
    if x_bxdiag_mode() == 1 {
        assert!(widths.iter().any(|w| w.1), "{widths:?}");
        assert!(widths.iter().any(|w| !w.1), "{widths:?}");
    }
    eprintln!("16-bit tables: {widths:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:513-606 (as above), with the
// word length of the search's lookup table (`wordsize` is `lookup->word_length`).
// No NCBI counterpart (speed only): the seed stage of the word sizes other than 3 and 5 (item SC9)
// scans the tables the search builds (`Lookup::x_new_for` with `x_narrow_cutoff`: for 2 and 4
// the lookup table, of 16-bit hits where the fused loop takes it, and for 6 and 7 the
// compressed table), with BLOSUM62 and another matrix, two-hit and one-hit, to the reference's
// hits and diagonal table. Item SC9b-1: the word sizes 2 to 4 of a two-hit search also take the
// sentinel walks (`query_fits` takes them, `XSent` with the instance of the word size, chosen
// by the loops once per subject), without and with the shadow's comparison, to the same hits
// and table.
#[test]
fn every_word_size_s_seed_stage_gives_the_reference_hits() {
    fn seed_stage<S: BxScoring>(
        scoring: S,
        search: &(PreparedQueryBatch, Vec<ContextParameters>, ResolvedOptions),
        lookup: &BxLookup,
        subjects: &[Vec<u8>],
        run: &str,
    ) -> usize {
        let (batch, params, options) = search;
        let length = batch.contexts.last().map(|c| c.offset + c.length).unwrap();
        let mut reference = Diagonals::new(length, options.window_size).unwrap();
        let mut scanned = Diagonals::new(length, options.window_size).unwrap();
        // The sentinel walks (`query_fits`: word size 2 to 4, a two-hit window), each on its
        // own table.
        let with_sent = (2..=4).contains(&options.word_size) && options.window_size > 0;
        let mut sent_tables = [
            Diagonals::new(length, options.window_size).unwrap(),
            Diagonals::new(length, options.window_size).unwrap(),
        ];
        let mut sent_runs = 0;
        let mut hits = 0;
        for subject in subjects {
            let expected = word_finder_bx_of(
                scoring,
                batch,
                params,
                options,
                subject,
                lookup,
                &mut reference,
                |_| {},
            )
            .unwrap();
            let actual = crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
                scoring,
                batch,
                params,
                options,
                subject,
                lookup,
                &mut scanned,
                None,
            )
            .unwrap();
            assert_eq!(actual, expected, "{run}");
            assert!(
                scanned.offset == reference.offset && scanned.entries == reference.entries,
                "{run}"
            );
            if with_sent {
                let lead: Vec<u8> = std::iter::once(0).chain(subject.iter().copied()).collect();
                for (d, shadow) in sent_tables.iter_mut().zip([false, true]) {
                    let sent = XSent::for_test(scoring, &batch.sequence_start, &lead, shadow);
                    let got = crate::algorithm::blastx::seed::x_seed_scan::x_word_finder_scan(
                        scoring,
                        batch,
                        params,
                        options,
                        subject,
                        lookup,
                        d,
                        Some(sent),
                    )
                    .unwrap();
                    assert_eq!(got, expected, "{run} sent shadow={shadow}");
                    assert!(
                        d.offset == reference.offset && d.entries == reference.entries,
                        "{run} sent shadow={shadow}"
                    );
                    sent_runs += 1;
                }
            }
            hits += actual.len();
        }
        assert_eq!(sent_runs > 0, with_sent, "{run}");
        hits
    }
    let (query, subjects, mean) = other_matrix_inputs();
    let mut tables = Vec::new();
    // Each word size with BLOSUM62 and PAM30, two-hit or one-hit (a table of 4 letters
    // takes seconds to build in a debug build).
    for (matrix, word, window) in [
        ("BLOSUM62", "2", 40),
        ("BLOSUM62", "2", 0),
        ("PAM30", "2", 40),
        ("BLOSUM62", "4", 40),
        ("PAM30", "4", 0),
        ("BLOSUM62", "6", 40),
        ("PAM30", "6", 0),
        ("BLOSUM62", "7", 0),
        ("PAM30", "7", 40),
    ] {
        let (batch, params, options, _) =
            search_with(&query, window, mean, false, matrix, &["-word_size", word]);
        let lookup =
            Lookup::x_new_for(&batch, &params, &options, 1, x_narrow_cutoff(&options)).unwrap();
        let kind = match &lookup {
            BxLookup::Narrow(l) => ("16-bit", l.table.word_length),
            BxLookup::Wide(Lookup::Standard(l, _)) => ("standard", l.word_length),
            BxLookup::Wide(Lookup::Compressed(l)) => ("compressed", x_compressed_word_length(l)),
        };
        assert_eq!(kind.1.to_string(), word, "{matrix} {word}");
        assert_eq!(kind.0 == "compressed", options.compressed_lookup());
        let run = format!("{matrix} -word_size {word} window={window} {}", kind.0);
        let search = (batch, params, options);
        let hits = if matrix == "BLOSUM62" {
            seed_stage(Blosum62Bx, &search, &lookup, &subjects, &run)
        } else {
            seed_stage(
                MatrixBx::of(search.2.scoring),
                &search,
                &lookup,
                &subjects,
                &run,
            )
        };
        assert!(hits > 0, "{run}");
        tables.push(run);
    }
    // Where the fused loop acts, the two-hit tables of 2 letters take 16-bit hits.
    if x_bxdiag_mode() == 1 {
        assert!(tables.iter().any(|t| t.contains("16-bit")), "{tables:?}");
    }
    eprintln!("tables: {tables:?}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:446-469 (cited at
// `Lookup::new`).
// Item SC9: BLASTX's table of 4-letter words, as the search builds it (`Lookup::x_new_for`: the
// joint builder on BLASTX's own lookup query, on one thread and in 4 parts on a pool of 4
// threads), is the reference builder's table of the same lookup query (the batch's buffer,
// lookup segments and contexts), field by field: BLOSUM62 two-hit and one-hit, PAM30, IDENTITY
// ungapped, and a -query_loc range.
#[test]
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn the_table_of_four_letters_is_the_reference_table() {
    let (query, _, mean) = other_matrix_inputs();
    for (matrix, window, extra) in [
        ("BLOSUM62", 40, &[][..]),
        ("BLOSUM62", 0, &[][..]),
        ("PAM30", 15, &[][..]),
        ("IDENTITY", 40, &["-ungapped", "-comp_based_stats", "0"][..]),
        ("BLOSUM62", 40, &["-query_loc", "1001-4000"][..]),
    ] {
        let mut args = vec!["-word_size", "4"];
        args.extend_from_slice(extra);
        let (batch, params, options, _) = search_with(&query, window, mean, false, matrix, &args);
        let concat = batch.sequence_start[1..].to_vec();
        let reference = crate::algorithm::tblastx::lookup::x_reference_table_from_prepared(
            options.scoring,
            4,
            concat.clone(),
            batch.lookup_segments.clone(),
            lookup_contexts(&batch, &params, &concat),
            crate::core::blast_util::ncbi_int4_from_double(options.threshold),
        );
        let check = |parts: usize| {
            let lookup = Lookup::x_new_for(&batch, &params, &options, parts, None).unwrap();
            let BxLookup::Wide(Lookup::Standard(table, _)) = &lookup else {
                panic!("{matrix}: the standard table of Int4 hits");
            };
            crate::algorithm::tblastx::lookup::x_assert_same_table_as(
                &format!("{matrix} window {window} {extra:?} parts {parts}:"),
                table,
                &reference,
            );
        };
        check(1);
        crate::utils::threading::with_search_pool(4, "blastx", |pool| {
            pool.install(|| check(4));
            Ok(())
        })
        .unwrap();
        assert!(reference.longest_chain > 0, "{matrix}: an empty table");
    }
}

// No NCBI counterpart (memory only): a search of two threads or more builds its tables of 2 and 3
// letters pooled and gives their vectors back when they are dropped; a table of 4 letters (the
// joint builder, which takes no kept vectors, item SC9) is not marked pooled, so that the pool
// does not keep every dropped table of the search.
#[test]
fn a_table_of_four_letters_gives_no_vectors_back() {
    let (query, _, mean) = other_matrix_inputs();
    for (matrix, word) in [
        ("BLOSUM62", "2"),
        ("BLOSUM62", "4"),
        ("PAM30", "3"),
        ("PAM30", "4"),
    ] {
        let extra = ["-word_size", word, "-num_threads", "8"];
        let (batch, params, options, _) = search_with(&query, 40, mean, false, matrix, &extra);
        let lookup = Lookup::x_new_for(&batch, &params, &options, 1, None).unwrap();
        let BxLookup::Wide(Lookup::Standard(table, kept)) = &lookup else {
            panic!("{matrix} {word}: the standard table of Int4 hits");
        };
        assert_eq!(table.word_length.to_string(), word);
        let pooled = crate::algorithm::tblastx::lookup::x_bx_pooled_for(8);
        assert_eq!(kept.pooled, pooled && word != "4", "{matrix} {word}");
    }
}

// No NCBI counterpart (speed only): the choice of the loop names a table by its build number,
// so a table at the address and with the vectors of the one before it is worked out again.
#[test]
fn the_loop_choice_is_worked_out_again_for_each_table_build() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/blastx_f_diagonal");
    let query = std::fs::read(dir.join("query.fna")).unwrap();
    let (_, _, _, lookup) = search(&query, 40, 500, false);
    let BxLookup::Wide(Lookup::Standard(table, kept)) = &lookup else {
        panic!("a standard table");
    };
    let next = crate::algorithm::blastx::seed::x_seed_scan::XKeptPresence::new();
    assert!(kept.build > 0 && next.build > kept.build);
    let (ratio, first) = x_table_ratio(table, kept.build);
    assert!(first);
    assert_eq!(x_table_ratio(table, kept.build), (ratio, false));
    // The same table under the next build number, as a pooled rebuild would give it.
    assert_eq!(x_table_ratio(table, next.build), (ratio, true));
    assert_eq!(x_table_ratio(table, next.build), (ratio, false));
}
