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

use std::ffi::{OsStr, OsString};
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
const REPEAT_MASK_LINK_VALUE: u8 = 5;

/// The letters next to a query's searched letters whose masks can join theirs (`combine`): a mask
/// that ends `REPEAT_MASK_LINK_VALUE` or more letters before an interval does not join it.
const REACH: usize = REPEAT_MASK_LINK_VALUE as usize;

/// The masks of a query's record next to its `-query_loc` interval, which NCBI's masks of the
/// query hold (the lower case of the whole record, `-lcase_masking`) and which join the
/// interval's own masks when its repeat hits are combined with them (`RepeatFilter::add_masks`).
/// The masks of a query without a range, or without `-lcase_masking`, are the default (none).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct MasksOutside {
    /// Bit `k`: the letter `k + 1` letters before the interval is masked (`k < REACH`): lower
    /// case, or a letter that the combination of the batch's masks joined to the interval's
    /// (`add_masks`), which a chunk of a split batch restricts (`query_split::masked_before_part`).
    pub(crate) before: u8,
    /// The first masked letter after the interval, in letters after its last (up to `REACH`).
    pub(crate) after: Option<usize>,
}

impl MasksOutside {
    /// The lower case of `record` next to its letters `offset..offset + length`.
    pub(crate) fn lower_case(record: &[u8], offset: usize, length: usize) -> Self {
        let lower = |index: Option<usize>| {
            index
                .and_then(|index| record.get(index))
                .is_some_and(u8::is_ascii_lowercase)
        };
        let before = (0..REACH)
            .filter(|&k| lower(offset.checked_sub(k + 1)))
            .fold(0u8, |bits, k| bits | 1 << k);
        let after = (1..=REACH).find(|&a| lower(Some(offset + length + a - 1)));
        Self { before, after }
    }

    /// The masks of a part of a split batch's chunk: those the chunk keeps before the part's
    /// letters end `before` letters before its first (`query_split::masked_before_part`).
    pub(crate) fn before_part(before: Option<usize>) -> Self {
        let before = before
            .filter(|&letters| letters <= REACH)
            .map_or(0, |letters| 1 << (letters - 1));
        Self {
            before,
            after: None,
        }
    }

    /// The nearest masked letter before the interval, in letters before its first.
    fn nearest_before(&self) -> Option<usize> {
        (self.before != 0).then(|| self.before.trailing_zeros() as usize + 1)
    }
}

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
    database: OsString,
    settings: ApplicationSettings,
}

/// Why the repeat search of a batch or chunk failed.
#[derive(Debug)]
pub(crate) enum RepeatError {
    /// An exception of NCBI's search of the repeats, which the set-up of each query of the batch
    /// posts as a warning (`batch_error`): its message (`GetMsg`). The errors of the database
    /// reader where the database is opened, and the stop of a search whose own set-up found no
    /// valid query (`CBlastException` eSetup, the text after `BLAST engine error: `).
    Posted(Vec<u8>),
    /// Any other error: LOSAT's own rejection of the search, or one that NCBI does not report
    /// as LOSAT does (`rejection`).
    Search(anyhow::Error),
}

impl RepeatFilter {
    /// The repeat search of the database `database`, looked up with `settings`.
    pub(crate) fn new(database: OsString, settings: ApplicationSettings) -> Self {
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
        let mut argv: Vec<OsString> = vec![
            OsString::from("losat"),
            OsString::from("blastn"),
            OsString::from("-db"),
            self.database.clone(),
            OsString::from("-strand"),
            OsString::from(strand),
        ];
        argv.extend(SEARCH_OPTIONS.into_iter().map(OsString::from));
        let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).map_err(|error| {
            anyhow::anyhow!("the repeat search of -filtering_db has arguments that LOSAT's BLASTN does not read ({error})")
        })?;
        let crate::cli::Commands::Blastn(mut args) = cli.command else {
            unreachable!("the argv names blastn")
        };
        args.num_threads = threads;
        args.one_query_batch = true;
        args.hits_only = true;
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
        for hit in lists
            .iter()
            .flat_map(|list: &Arc<Vec<PairwiseHit>>| list.iter())
        {
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
    /// Returns the hits, which the BLAST archive of a query batch reads (`archive::CallBatches`'s
    /// `repeat_hits`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:236-306
    ///
    /// NCBI combines the masks in the record's letters, so with `-query_loc` the masks outside
    /// the interval join too (`outside`, one for each query, or empty without a range); the
    /// search then takes the masks' part in the interval. The letters before the interval that
    /// a query's combination joins are added to its `outside` (`MasksOutside::before`).
    pub(crate) fn add_masks<L: QueryLetters + ?Sized>(
        &self,
        strand: &str,
        queries: &L,
        masks: &mut [Vec<MaskedInterval>],
        outside: &mut [MasksOutside],
        pool: &SearchPool<'_>,
    ) -> Result<Vec<Vec<(usize, usize)>>, RepeatError> {
        let found = self.hits(strand, queries, pool)?;
        for (index, (query_masks, hits)) in masks.iter_mut().zip(&found).enumerate() {
            if !hits.is_empty() {
                *query_masks = join(
                    hits,
                    query_masks,
                    queries.length(index),
                    outside.get_mut(index),
                );
            }
        }
        Ok(found)
    }

    /// The database's name, as `-filtering_db` gives it.
    pub(crate) fn database(&self) -> &OsStr {
        &self.database
    }
}

/// The masks of a query of `length` letters with its repeat `hits` (0-based, inclusive) joined to
/// its `masks` (`combine`), in the record's letters with the masks `outside` its interval, then
/// cut to the interval; the letters before the interval that the combination joins are added to
/// `outside` (`RepeatFilter::add_masks`).
fn join(
    hits: &[(usize, usize)],
    masks: &[MaskedInterval],
    length: usize,
    outside: Option<&mut MasksOutside>,
) -> Vec<MaskedInterval> {
    let last = length as i64 - 1;
    let mut intervals: Vec<(i64, i64)> = hits
        .iter()
        .copied()
        .chain(masks.iter().map(|mask| (mask.start, mask.end - 1)))
        .map(|(left, right)| (left as i64, right as i64))
        .collect();
    if let Some(outside) = outside.as_deref() {
        if let Some(before) = outside.nearest_before() {
            intervals.push((-(before as i64), -(before as i64)));
        }
        if let Some(after) = outside.after {
            intervals.push((last + after as i64, last + after as i64));
        }
    }
    let combined = combine(intervals);
    if let Some(outside) = outside {
        for &(left, right) in combined.iter().filter(|interval| interval.0 < 0) {
            for k in 0..REACH {
                let letter = -(k as i64) - 1;
                if left <= letter && letter <= right {
                    outside.before |= 1 << k;
                }
            }
        }
    }
    combined
        .into_iter()
        .filter_map(|(left, right)| {
            let (left, right) = (left.max(0), right.min(last));
            (left <= right).then(|| MaskedInterval::new(left as usize, right as usize + 1))
        })
        .collect()
}

/// A query's masks with its repeat hits (0-based, inclusive): `BlastSeqLocCombine` with
/// `REPEAT_MASK_LINK_VALUE`, sorted by their starts (`common::masks::combine_linked`).
pub(crate) fn combine<T>(intervals: Vec<(T, T)>) -> Vec<(T, T)>
where
    T: Copy + Ord + std::ops::Add<Output = T> + From<u8>,
{
    crate::algorithm::common::masks::combine_linked(intervals, T::from(REPEAT_MASK_LINK_VALUE))
}

/// NCBI's stop at the set-up of a query batch whose repeat search failed with `error`, with
/// `queries` queries in it: each query's set-up runs the search again (the masks are made
/// only once they are made) and posts its error as a warning, and the batch, left without a valid
/// query, stops with them (`report::query_warnings::each_query_error`; exit 3), as WindowMasker's
/// error of a batch does (`window_masker::batch_error`). `None` for LOSAT's own errors
/// (`rejection`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:176-228
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
pub(crate) fn batch_error(error: &RepeatError, queries: usize) -> Option<crate::cli::NativeError> {
    let RepeatError::Posted(text) = error else {
        return None;
    };
    Some(crate::report::query_warnings::each_query_error(
        text, queries,
    ))
}

/// LOSAT's rejection of a repeat search's error that NCBI does not report as `batch_error` does:
/// LOSAT's own refusals of the search, and an error of the masks of a chunk of a split batch
/// whose whole queries were masked (NCBI masks each chunk when it searches it, after the batch's
/// set-up).
pub(crate) fn rejection(filter: &RepeatFilter, error: RepeatError) -> anyhow::Error {
    match error {
        // The text alone: the search's own error (a `NativeError` that the command line would
        // write as NCBI's) is not this run's.
        RepeatError::Search(error) => anyhow::anyhow!(
            "the repeat search of -filtering_db {} failed: {error:#}",
            filter.database.to_string_lossy()
        ),
        RepeatError::Posted(text) => anyhow::anyhow!(
            "the repeat search of -filtering_db {} stops with an error ({}) where NCBI BLAST+ reports it otherwise; this is not supported by LOSAT's BLASTN",
            filter.database.to_string_lossy(),
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
#[path = "repeat_filter_tests.rs"]
mod tests;
