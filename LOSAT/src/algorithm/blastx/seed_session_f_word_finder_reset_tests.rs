use super::session_f_diagonal_boundary_tests::flagged;
use super::*;
use crate::algorithm::blastx::{input::parse_fasta, parameters, query_setup::prepare_queries};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-174
// Expected SEED/INIT rows are frozen fresh pinned-NCBI trace bytes, not Rust output.
#[test]
fn fresh_retained_and_threshold_reset_word_finder_match_ncbi() {
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
            crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
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
        let seeds = expected
            .lines()
            .filter(|x| x.starts_with("SEED\t"))
            .map(|x| {
                let f = x.split('\t').collect::<Vec<_>>();
                (f[2].parse::<i32>().unwrap(), f[3].parse::<i32>().unwrap())
            })
            .collect::<Vec<_>>();
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
        let mut d = Diagonals::new(length, window).unwrap();
        for state in 0..4 {
            if state == 2 {
                d.offset = 536_870_910;
                d.entries.fill(flagged(123));
            }
            if state == 3 {
                d.offset = 536_870_911;
                d.entries.fill(flagged(536_870_900));
                d.finish_subject(83).unwrap();
            }
            let mut actual = Vec::new();
            let hits = word_finder(&batch, &params, &options, &subject, &lookup, &mut d, |x| {
                actual.extend_from_slice(x)
            })
            .unwrap();
            assert_eq!(actual, seeds, "window={window} state={state}");
            assert_eq!(hits, initial, "window={window} state={state}");
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:742-798
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:137
// A negative window (a search strategy's) runs the one-hit word finder on a table whose entries
// start at `last_hit` 0 and whose offset starts at the window: a word of a diagonal that no
// extension has touched is extended only from `s_off >= -window`. The scan is the one-hit scan's,
// and the table's state (offset and entries) carries into the next call, so the same words on a
// table that is old enough are all extended.
#[test]
fn negative_window_one_hit_word_finder_reads_the_table_at_the_window_and_carries_it() {
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
    // One search of the fixture's subject with `window` on a table that `table` gives (the table
    // `Diagonals::new` makes, by default); the seeds, the hits and the table after it.
    let run = |window: i32,
               table: &dyn Fn(usize) -> Diagonals|
     -> (Vec<(i32, i32)>, Vec<InitHsp>, Diagonals) {
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
            "0",
        ];
        let crate::cli::Commands::Blastx(args) =
            crate::cli::try_parse_from::<crate::cli::Cli, _, _>(args)
                .unwrap()
                .command
        else {
            unreachable!()
        };
        let mut options = args.resolve().unwrap();
        options.window_size = window;
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
        let mut d = table(length);
        let mut seeds = Vec::new();
        let hits = word_finder(&batch, &params, &options, &subject, &lookup, &mut d, |x| {
            seeds.extend_from_slice(x)
        })
        .unwrap();
        (seeds, hits, d)
    };
    let fresh = |window: i32| move |length: usize| Diagonals::new(length, window).unwrap();
    let (seeds, hits, after) = run(-40, &fresh(-40));
    // The table's own length for the one-hit reference below.
    let entries = after.entries.len();
    let shaped = |offset: i64, window: i32| {
        move |length: usize| {
            let mut d = Diagonals::new(length, window).unwrap();
            d.entries = vec![DiagStruct::default(); entries];
            d.wide_offset = offset;
            d.offset = offset as i32;
            d
        }
    };
    // The same scan as `-window_size 0`'s, on a table of the same length.
    let (seeds_0, hits_0, _) = run(0, &shaped(0, 0));
    assert_eq!(seeds, seeds_0);
    // A fresh table extends no word of an untouched diagonal before `s_off` 40.
    assert!(hits.iter().all(|hit| hit.s_seed >= 40), "{hits:?}");
    // The offset moved by the subject's length plus the window.
    let residues = subject.len() as i64 - 1;
    assert_eq!(after.wide_offset, -40 + residues - 40);
    assert_eq!(i64::from(after.offset), after.wide_offset);
    // A table whose entries are all older than the window (`last_hit` 0 at offset 1000) extends
    // every word, as `-window_size 0` on a fresh table does.
    let (seeds_old, hits_old, _) = run(-40, &shaped(1000, -40));
    assert_eq!(seeds_old, seeds_0);
    assert_eq!(hits_old, hits_0);
}
