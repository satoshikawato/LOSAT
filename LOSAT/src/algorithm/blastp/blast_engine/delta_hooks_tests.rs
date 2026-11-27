//! Tests of `delta_hooks` (series PR, item PS13a, fix round 1).

use crate::algorithm::pssm_search::in_pssm::query_ungapped_block;
use crate::config::ScoringMatrix;

/// `DeltaSide::cdd_pssm` refuses a query whose letters give a Karlin-Altschul block with the domain
/// search's matrix and none with `-matrix`, the case where NCBI's engine would read a null block.
/// With the matrices LOSAT knows (NCBI's), whether a query has the block depends on its set of
/// letters only: each letter gives a block under every matrix or under none, and every letter's
/// mean score against the standard composition is negative under every matrix, so a mix of letters
/// gives a block wherever one of its letters does. So the refusal is not reached; this test (every
/// letter, and every pair of letters) keeps it so if a matrix is added.
#[test]
fn every_letter_set_gives_a_karlin_block_under_every_matrix_or_none() {
    let matrices: Vec<ScoringMatrix> = ScoringMatrix::ALL
        .into_iter()
        .filter(|&matrix| matrix != ScoringMatrix::Identity)
        .collect();
    for first in 1u8..28 {
        for second in first..28 {
            let mut query = vec![first; 20];
            query.extend(std::iter::repeat_n(second, 20));
            let blocks: Vec<bool> = matrices
                .iter()
                .map(|&matrix| query_ungapped_block(matrix, &query).is_some())
                .collect();
            assert!(
                blocks.iter().all(|&block| block == blocks[0]),
                "letters {first} and {second}: {blocks:?} under {matrices:?}"
            );
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
// NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:195-215
/// deltablast's `-searchsp` sets the process's OLD_FSC, which the domain search reads too (its
/// options are a new RPS-BLAST options handle, which takes neither `-searchsp`'s value nor
/// deltablast's other options): the run's `DeltaRun` carries the run's `old_fsc` (row
/// `db.dl.fsc.6`), and the domain search's arguments keep no search space.
#[test]
fn the_domain_search_reads_the_runs_old_fsc() {
    let parse = |extra: &[&str]| {
        let mut argv = vec![
            "losat",
            "deltablast",
            "-query",
            "q.faa",
            "-subject",
            "s.faa",
        ];
        argv.extend_from_slice(extra);
        match crate::cli::try_parse_from(argv) {
            Ok(crate::cli::Cli {
                command: crate::cli::Commands::Deltablast(args),
            }) => args,
            Ok(other) => panic!("not deltablast: {other:?}"),
            Err(error) => panic!("{error}"),
        }
    };
    let environment = crate::blastinput::app::old_fsc_environment();
    for (extra, searchsp) in [(&[][..], false), (&["-searchsp", "1000000"][..], true)] {
        let resolved = parse(extra).into_blastp_args().unwrap().resolve().unwrap();
        assert_eq!(resolved.old_fsc, searchsp || environment, "{extra:?}");
        let delta = resolved.delta.as_deref().expect("deltablast's run");
        assert_eq!(delta.old_fsc, resolved.old_fsc, "{extra:?}");
        assert_eq!(delta.options.domain.searchsp, None, "{extra:?}");
    }
}
