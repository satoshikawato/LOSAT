//! The composition redo's two forms of the matches it makes ahead (`LOSAT_X_BXREDOFLOW`), the flow
//! and the fork-join batches of `=0`, chosen per search by the observer
//! (`RuntimeObserver::x_redo_flow`) in one process: the same hit lists, every field of every HSP
//! list, and the same early-termination decisions of the loop in the same order. The inputs (one
//! query of 9,000 letters of AP027131 against its 595 proteins, a short hit list) make the loop's
//! early-termination test skip most matches while the flow's ring of matches ahead wraps; that the
//! report equals one thread's is `tests/blastx_unit_cap.rs`.

use super::super::{search_internal_observed, RuntimeObserver};
use crate::algorithm::blastx::input::{read_fasta, FastaRecord};
use crate::algorithm::blastx::results::ResultsObserver;
use crate::cli::Commands;
use std::path::Path;

/// An observer that lets the redo be made ahead (`x_passive`), chooses its form, and keeps the
/// loop's early-termination decisions (`early` is called by the loop in stream order, also when
/// the redo is made ahead).
struct Choose {
    flow: bool,
    decisions: Vec<bool>,
}

impl ResultsObserver for Choose {}

impl RuntimeObserver for Choose {
    fn early(&mut self, _evalue: f64, _queries: usize, decision: bool) {
        self.decisions.push(decision);
    }
    fn x_passive(&self) -> bool {
        true
    }
    fn x_redo_flow(&self) -> Option<bool> {
        Some(self.flow)
    }
}

fn fixture(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The first 9,000 letters of AP027131 and its 595 proteins.
fn inputs() -> (Vec<FastaRecord>, Vec<FastaRecord>) {
    let mut genome = read_fasta(&fixture("AP027131.fasta"), false, false).unwrap();
    genome.truncate(1);
    genome[0].sequence.truncate(9_000);
    let proteins = read_fasta(&fixture("AP027131.faa"), true, false).unwrap();
    assert_eq!(proteins.len(), 595);
    (genome, proteins)
}

fn options(threads: usize, extra: &[&str]) -> crate::algorithm::blastx::args::ResolvedOptions {
    let threads = threads.to_string();
    let mut argv = vec![
        "losat",
        "blastx",
        "-query",
        "unused.fna",
        "-subject",
        "unused.faa",
        "-num_threads",
        threads.as_str(),
    ];
    argv.extend_from_slice(extra);
    let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
    let Commands::Blastx(args) = cli.command else {
        unreachable!()
    };
    args.resolve().unwrap()
}

/// One search with the form `flow`: the hit lists of each batch (`Debug`: every field) and the
/// loop's early-termination decisions.
fn search(
    queries: &[FastaRecord],
    proteins: &[FastaRecord],
    threads: usize,
    extra: &[&str],
    flow: bool,
) -> (String, Vec<bool>) {
    let mut observer = Choose {
        flow,
        decisions: Vec::new(),
    };
    let results =
        search_internal_observed(queries, proteins, &options(threads, extra), &mut observer)
            .unwrap();
    let lists: Vec<_> = results
        .iter()
        .map(|batch| (batch.query_ordinal, batch.search_skipped, &batch.queries))
        .collect();
    (format!("{lists:?}"), observer.decisions)
}

#[test]
fn the_flow_and_the_batches_give_the_same_hit_lists() {
    let (queries, proteins) = inputs();
    for extra in [
        &["-comp_based_stats", "1", "-max_target_seqs", "2"][..],
        &["-use_sw_tback", "-max_target_seqs", "3"][..],
    ] {
        for threads in [2, 4] {
            let (flow, flow_decisions) = search(&queries, &proteins, threads, extra, true);
            let (batches, batch_decisions) = search(&queries, &proteins, threads, extra, false);
            assert!(
                flow == batches,
                "{extra:?}, {threads} threads: hit lists differ"
            );
            assert_eq!(
                flow_decisions, batch_decisions,
                "{extra:?}, {threads} threads"
            );
            // The loop's early-termination test skipped some matches and not others.
            assert!(
                flow_decisions.contains(&true) && flow_decisions.contains(&false),
                "{extra:?}, {threads} threads: decisions {flow_decisions:?}"
            );
        }
    }
}
