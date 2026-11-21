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
