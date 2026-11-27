//! The parts of a search's results, one per search call of the engine: the whole search, or with
//! psiblast's rounds one part per round of each query (`FormatObserver::round_begin`). Each part
//! has its own hit list, HSP numbers, statistics and database ranks; this is where `run_collect`
//! and `run_batches` turn them into the records of the whole search.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:250-285 (cited at
//! `local_blast::FormatObserver::round_begin`: each round's `CSearchResultSet`, formatted with its
//! iteration number)
#![warn(missing_docs, unused)]

use std::sync::Arc;

use super::input::Inputs;
use super::local_blast::ReportedSubject;
use super::program::Program;
use super::run::{query_ordinal, records, HspRecord, QueryRounds, QueryStats};
use crate::blastinput::seq_range::SequenceRange;
use crate::record::FastaRecord;
use crate::report::blast_xml::XmlStats;
use crate::report::PairwiseHit;

/// One search call's part of the results.
struct Part {
    /// NCBI's round number: 1 in a search without rounds.
    round: u32,
    /// The searched query of a round, which the call numbers 0; 0 in a search without rounds,
    /// whose calls number the searched queries themselves.
    query: usize,
    /// By the call's query number: the statistics.
    query_stats: Vec<Option<QueryStats>>,
    /// The sequences that a database search's hit list names, by the call's ranks.
    database_subjects: Option<Vec<ReportedSubject>>,
    /// PHI-BLAST's pattern occurrence of each HSP of the call's hit list
    /// (`FormatObserver::phi_occurrences`); `None` for a call without a pattern.
    occurrences: Option<Vec<u32>>,
}

impl Part {
    fn new(query: usize, round: u32) -> Self {
        Self {
            round,
            query,
            query_stats: Vec::new(),
            database_subjects: None,
            occurrences: None,
        }
    }
}

/// The parts of one search, from the observer's events (`FormatObserver`), and the queries whose
/// rounds converged.
pub(super) struct Parts {
    parts: Vec<Part>,
    /// Whether the search has rounds (a `round_begin` came).
    rounds: bool,
    /// The searched queries whose rounds converged.
    converged: Vec<usize>,
    /// The query of a PSSM search (`FormatObserver::pssm_query`), which the inputs do not hold.
    pssm_query: Option<FastaRecord>,
}

/// The records of a search's parts (`Parts::finish`).
pub(super) struct Finished {
    /// The HSPs of every part, in the order of the parts: `index` counts over them from the
    /// `first_index` given.
    pub(super) hsps: Vec<HspRecord>,
    /// By part: how many of `hsps` it holds.
    pub(super) lengths: Vec<usize>,
    /// By query record of the inputs: the statistics of its first round.
    pub(super) query_stats: Vec<Option<QueryStats>>,
    /// By query record of the inputs: its rounds.
    pub(super) rounds: Vec<QueryRounds>,
}

impl Parts {
    /// A search's parts before any event: one part, the whole search.
    pub(super) fn new() -> Self {
        Self {
            parts: vec![Part::new(0, 1)],
            rounds: false,
            converged: Vec::new(),
            pssm_query: None,
        }
    }

    /// The part that the call now running fills.
    pub(super) fn current(&self) -> usize {
        self.parts.len() - 1
    }

    /// `FormatObserver::round_begin`: the next call is round `round` of the searched query
    /// `query`. The first round of the search takes the part that no call has filled yet.
    pub(super) fn round_begin(&mut self, query: usize, round: u32) {
        if self.rounds {
            self.parts.push(Part::new(query, round));
        } else {
            self.rounds = true;
            self.parts[0] = Part::new(query, round);
        }
    }

    /// `FormatObserver::query_stats`.
    pub(super) fn query_stats(&mut self, q_idx: usize, stats: XmlStats) {
        let part = self.parts.last_mut().expect("a search has a part");
        if part.query_stats.len() <= q_idx {
            part.query_stats.resize(q_idx + 1, None);
        }
        part.query_stats[q_idx] = QueryStats::of_context(stats);
    }

    /// `FormatObserver::database_subjects`.
    pub(super) fn database_subjects(&mut self, subjects: Vec<ReportedSubject>) {
        self.parts
            .last_mut()
            .expect("a search has a part")
            .database_subjects = Some(subjects);
    }

    /// `FormatObserver::rounds_converged`.
    pub(super) fn rounds_converged(&mut self, query: usize) {
        self.converged.push(query);
    }

    /// `FormatObserver::pssm_query`.
    pub(super) fn pssm_query(&mut self, query: &FastaRecord) {
        self.pssm_query = Some(query.clone());
    }

    /// `FormatObserver::phi_occurrences`.
    pub(super) fn phi_occurrences(&mut self, occurrences: Vec<u32>) {
        self.parts
            .last_mut()
            .expect("a search has a part")
            .occurrences = Some(occurrences);
    }

    /// The records of the parts, from the hit lists that the calls handed over (`hits`, one per
    /// call, in the order of the calls), with `index` counted from `first_index`; the statistics
    /// and rounds of each query of `inputs`, whose searched queries are `ordinals`
    /// (`run::searched_queries`); the query of a PSSM search is the one the search reported
    /// (`pssm_query`). A search that handed over no hit list has no records. The error is a fault
    /// of the engine: a search of rounds whose hit lists are not one per round, a PSSM search with
    /// hits that reported no query.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish(
        self,
        program: Program,
        inputs: Inputs<'_>,
        hits: &[Arc<Vec<PairwiseHit>>],
        query_range: Option<SequenceRange>,
        ordinals: Option<&[usize]>,
        first_index: usize,
    ) -> Result<Finished, String> {
        if !hits.is_empty() && hits.len() != self.parts.len() {
            return Err(format!(
                "{} hit lists for {} parts of the search",
                hits.len(),
                self.parts.len()
            ));
        }
        let inputs = match (&self.pssm_query, inputs.query_file) {
            (Some(query), Some(_)) => Inputs {
                query: std::slice::from_ref(query),
                ..inputs
            },
            (None, Some(_)) if hits.iter().any(|list| !list.is_empty()) => {
                return Err("a PSSM search reported hits but not its query".to_string());
            }
            _ => inputs,
        };
        let mut hsps = Vec::new();
        let mut lengths = Vec::with_capacity(self.parts.len());
        for (part, list) in self.parts.iter().zip(hits) {
            let mut records = records(
                program,
                inputs,
                list,
                query_range,
                ordinals,
                part.database_subjects.as_deref(),
                part.query,
            );
            for (at, record) in records.iter_mut().enumerate() {
                record.index += first_index + hsps.len();
                record.round = part.round;
                record.occurrence = part
                    .occurrences
                    .as_ref()
                    .and_then(|occurrences| occurrences.get(at).copied());
            }
            lengths.push(records.len());
            hsps.append(&mut records);
        }
        let queries = inputs.query_count();
        let mut query_stats = vec![None; queries];
        let mut round_stats: Vec<Vec<Option<QueryStats>>> = vec![Vec::new(); queries];
        for part in &self.parts {
            for (q_idx, stats) in part.query_stats.iter().enumerate() {
                let ordinal = query_ordinal(ordinals, part.query + q_idx);
                if ordinal >= queries {
                    continue;
                }
                if part.round == 1 {
                    query_stats[ordinal] = *stats;
                }
                round_stats[ordinal].push(*stats);
            }
        }
        let rounds = round_stats
            .into_iter()
            .enumerate()
            .map(|(ordinal, mut stats)| {
                // A query without a round has the statistics of `query_stats` (`None`).
                if stats.is_empty() {
                    stats.push(query_stats[ordinal]);
                }
                QueryRounds {
                    stats,
                    converged: self
                        .converged
                        .iter()
                        .any(|&query| query_ordinal(ordinals, query) == ordinal),
                }
            })
            .collect();
        Ok(Finished {
            hsps,
            lengths,
            query_stats,
            rounds,
        })
    }
}
