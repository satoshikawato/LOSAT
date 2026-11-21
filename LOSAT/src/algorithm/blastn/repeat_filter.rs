//! blastn's `-filtering_db`: the queries are searched against a database of repeats, and the
//! query ranges of the hits mask the queries as their lower case does.
//!
//! NCBI runs the search where it makes the masks of each query batch, after DUST and before
//! WindowMasker (the option excludes WindowMasker's), and again for each chunk of a split batch,
//! whose query source makes its masks anew (`coordination::query_masks`,
//! `coordination::chunk_query_masks`). The search is BLASTN's own engine
//! (`blast_engine::repeat_search_hits`) with the options of `SEARCH_OPTIONS`, whatever the main
//! search's options are; only the query strand (`-strand`) and the query locations are those of
//! the main search.
//!
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:176-224
// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:236-390

use std::sync::Arc;

use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::record::FastaRecord;
use crate::report::PairwiseHit;
use crate::utils::dust::MaskedInterval;
use crate::utils::threading::SearchPool;

use super::args::BlastnArgs;
use super::coordination::QueryLetters;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_filter.h:72
/// `REPEAT_MASK_LINK_VALUE`: a hit's interval joins the masks before it when it starts less
/// than 5 letters after the end of one (`common::masks::combine_linked`).
const REPEAT_MASK_LINK_VALUE: usize = 5;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:306-322
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_filter.h:57-66
/// The options of the repeat search (`s_CreateRepeatsSearchOptions`): blastn's traditional
/// defaults (`-task blastn`), then the mismatch penalty, the match reward, the cutoff score
/// (`-min_raw_gapped_score`), the final gapped X-drop, the ungapped X-drop, the gap costs, no
/// DUST and the word size of `REPEATS_SEARCH_*`. Every other option has the task's default (an
/// E-value of 10, 500 subjects): NCBI makes a new options handle for the search.
const SEARCH_OPTIONS: [&str; 20] = [
    "-task",
    "blastn",
    "-penalty",
    "-1",
    "-reward",
    "1",
    "-min_raw_gapped_score",
    "26",
    "-xdrop_gap_final",
    "90",
    "-xdrop_ungap",
    "40",
    "-gapopen",
    "2",
    "-gapextend",
    "1",
    "-dust",
    "no",
    "-word_size",
    "11",
];

/// The repeat search of `-filtering_db`: the name of the database, as the option gives it, and the
/// settings that look it up (`BLASTDB`, the registry), the caller's (`BlastnArgs::resolve_repeat_filter`).
#[derive(Debug)]
pub struct RepeatFilter {
    database: String,
    settings: ApplicationSettings,
}

/// Why the repeat search of a batch or chunk failed.
#[derive(Debug)]
pub(crate) enum RepeatError {
    /// An error of NCBI's database reader that NCBI reports for each query of the batch
    /// (`batch_error`): its text.
    Database(Vec<u8>),
    /// Any other error: LOSAT's own rejection of the search, or one that NCBI does not report
    /// as LOSAT does (`rejection`).
    Search(anyhow::Error),
}

impl RepeatFilter {
    /// The repeat search of the database `database`, looked up with `settings`.
    pub(crate) fn new(database: String, settings: ApplicationSettings) -> Self {
        Self { database, settings }
    }

    /// The settings that look the database up.
    pub(crate) fn settings(&self) -> &ApplicationSettings {
        &self.settings
    }

    /// The arguments of the search: `SEARCH_OPTIONS` against the database, the main search's
    /// `strand`, `threads` threads (the pool's, `blast_engine::repeat_search_hits`), every query
    /// in one batch, as NCBI's search of the repeats takes the queries of its caller's batch.
    pub(crate) fn search_args(&self, strand: &str, threads: usize) -> anyhow::Result<BlastnArgs> {
        let mut argv: Vec<&str> = vec!["losat", "blastn", "-db", &self.database, "-strand", strand];
        argv.extend(SEARCH_OPTIONS);
        let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).map_err(|error| {
            anyhow::anyhow!("the repeat search of -filtering_db has arguments that LOSAT's BLASTN does not read ({error})")
        })?;
        let crate::cli::Commands::Blastn(mut args) = cli.command else {
            unreachable!("the argv names blastn")
        };
        args.num_threads = threads;
        args.one_query_batch = true;
        Ok(args)
    }

    /// The hits of the repeat search of `queries` on the strands of `strand`: for each query, the
    /// query range of each HSP (0-based, inclusive, in the letters' coordinates; the range of a
    /// minus-strand HSP on the plus strand), in the order of the hit lists. The search reads the
    /// letters as they are (NCBI searches the queries without their masks).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:143-176
    pub(crate) fn hits<L: QueryLetters + ?Sized>(
        &self,
        strand: &str,
        queries: &L,
        pool: &SearchPool<'_>,
    ) -> Result<Vec<Vec<(usize, usize)>>, RepeatError> {
        let records: Vec<FastaRecord> = (0..queries.count())
            .map(|index| {
                queries.read(index, |letters| {
                    FastaRecord::new(format!("Query_{}", index + 1), b"", letters)
                })
            })
            .collect();
        let args = self
            .search_args(strand, pool.threads())
            .map_err(RepeatError::Search)?;
        let lists = super::blast_engine::repeat_search_hits(args, &records, self, pool)?;
        let mut found: Vec<Vec<(usize, usize)>> = vec![Vec::new(); records.len()];
        for hit in lists.iter().flat_map(|list: &Arc<Vec<PairwiseHit>>| list.iter()) {
            let (from, to) = (
                hit.hit.q_start.min(hit.hit.q_end),
                hit.hit.q_start.max(hit.hit.q_end),
            );
            found[hit.hit.q_idx as usize].push((from - 1, to - 1));
        }
        Ok(found)
    }

    /// Adds the masks of the repeats to `masks`, the masks of `queries` so far (DUST's and the
    /// lower case; for a chunk, the batch's masks cut to the chunk too): each query's hits
    /// (`hits`) join its masks (`combine`); a query without a hit keeps its masks as they are.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:236-306
    pub(crate) fn add_masks<L: QueryLetters + ?Sized>(
        &self,
        strand: &str,
        queries: &L,
        masks: &mut [Vec<MaskedInterval>],
        pool: &SearchPool<'_>,
    ) -> Result<(), RepeatError> {
        let found = self.hits(strand, queries, pool)?;
        for (query_masks, mut intervals) in masks.iter_mut().zip(found) {
            if intervals.is_empty() {
                continue;
            }
            intervals.extend(query_masks.iter().map(|mask| (mask.start, mask.end - 1)));
            *query_masks = combine(intervals)
                .into_iter()
                .map(|(left, right)| MaskedInterval::new(left, right + 1))
                .collect();
        }
        Ok(())
    }

    /// The database's name, as `-filtering_db` gives it.
    pub(crate) fn database(&self) -> &str {
        &self.database
    }
}

/// A query's masks with its repeat hits (0-based, inclusive): `BlastSeqLocCombine` with
/// `REPEAT_MASK_LINK_VALUE`, sorted by their starts (`common::masks::combine_linked`).
pub(crate) fn combine(intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    crate::algorithm::common::masks::combine_linked(intervals, REPEAT_MASK_LINK_VALUE)
}

/// NCBI's stop at the set-up of a query batch whose repeat search failed with `error`, with
/// `queries` queries in it: each query's set-up runs the search again (the masks are made
/// only once they are made) and posts its error as a warning, and the batch, left without a valid
/// query, stops with them (`TSearchMessages::ToString`: the severity, the text and a space for
/// each query) after `BLAST engine error: ` (exit 3), as WindowMasker's error of a batch does
/// (`window_masker::batch_error`). `None` for LOSAT's own errors (`rejection`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:176-228
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
pub(crate) fn batch_error(error: &RepeatError, queries: usize) -> Option<crate::cli::NativeError> {
    let RepeatError::Database(text) = error else {
        return None;
    };
    let mut message = b"BLAST engine error: ".to_vec();
    for _ in 0..queries {
        message.extend_from_slice(b"Warning: ");
        message.extend_from_slice(text);
        message.push(b' ');
    }
    message.push(b'\n');
    Some(crate::cli::NativeError { exit: 3, message })
}

/// LOSAT's rejection of a repeat search's error that NCBI does not report as `batch_error` does:
/// LOSAT's own refusals of the search, and an error of the masks of a chunk of a split batch
/// whose whole queries were masked (NCBI masks each chunk when it searches it, after the batch's
/// set-up).
pub(crate) fn rejection(filter: &RepeatFilter, error: RepeatError) -> anyhow::Error {
    match error {
        RepeatError::Search(error) => error.context(format!(
            "the repeat search of -filtering_db {} failed",
            filter.database
        )),
        RepeatError::Database(text) => anyhow::anyhow!(
            "the database of -filtering_db {} gives an error where NCBI BLAST+ reports it otherwise ({}); this is not supported by LOSAT's BLASTN",
            filter.database,
            String::from_utf8_lossy(&text)
        ),
    }
}

impl BlastnArgs {
    /// `-filtering_db`'s repeat search (`repeat_filter`), with the caller's `settings`; nothing
    /// without the option.
    pub(crate) fn resolve_repeat_filter(&mut self, settings: &ApplicationSettings) {
        self.repeat_filter = self
            .filtering_db
            .clone()
            .map(|database| Arc::new(RepeatFilter::new(database, settings.clone())));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn search_args_are_the_repeat_options() {
        let filter = RepeatFilter::new(
            "rpt".to_string(),
            crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        );
        let args = filter.search_args("plus", 4).expect("valid arguments");
        assert_eq!(args.db.as_deref(), Some("rpt"));
        assert_eq!(args.task, "blastn");
        assert_eq!((args.reward, args.penalty), (Some(1), Some(-1)));
        assert_eq!((args.gap_open, args.gap_extend), (Some(2), Some(1)));
        assert_eq!(args.word_size, Some(11));
        assert_eq!(args.min_raw_gapped_score, Some(26));
        assert_eq!(args.xdrop_ungap, Some(40.0));
        assert_eq!(args.xdrop_gap_final, Some(90.0));
        assert_eq!(args.strand, "plus");
        assert_eq!(args.num_threads, 4);
        assert!(args.one_query_batch);
        assert!(args.filtering_db.is_none());
    }

    #[test]
    fn batch_error_repeats_the_text_for_each_query() {
        let error = RepeatError::Database(b"Database name is required.".to_vec());
        let native = batch_error(&error, 2).expect("NCBI's error");
        assert_eq!(native.exit, 3);
        assert_eq!(
            native.message,
            b"BLAST engine error: Warning: Database name is required. Warning: Database name is required. \n"
        );
        assert!(batch_error(&RepeatError::Search(anyhow::anyhow!("x")), 2).is_none());
    }
}
