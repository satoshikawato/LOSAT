//! The reader batches of BLASTX and RPSTBLASTN as their loops search and report them
//! (`BatchTake::deliver_views`): under `-parse_deflines`, a query whose set-up NCBI reads past
//! the end of its earlier record in a batch that NCBI searches (`query_batch::PastFirstQuery`;
//! line RG-KB item 7d) is NCBI's view of it, an invalid query of the earlier record's title.
//!
//! This module is a child of `native` for the parts of `BatchTake` and `add_batch_scope` that the
//! delivery shares with the loop's take of a batch.

use super::{add_batch_scope, BatchTake};
use crate::algorithm::blastx::args::BlastxArgs;
use crate::algorithm::blastx::input::{scope_records, FastaRecord};
use crate::algorithm::blastx::input_source::QueryStreamSettings;
use crate::blastinput::query_batch::{PastFirstQuery, QueryScopes};
use anyhow::Result;

/// What `BatchTake::deliver_views` keeps from one reader batch to the next.
#[derive(Default)]
pub(super) struct Delivered<'a> {
    /// The run's searched queries and its input records (every record read, those that
    /// `-query_loc` skips among them: `input::scope_records`) of the batches delivered so far.
    counts: std::cell::Cell<(usize, usize)>,
    /// The titles of the input records of delivered batches that a later batch's query read past
    /// them shows, by input record, where the scopes knew that query when the record was delivered
    /// (read ahead, `QueryScopes::read_ahead`).
    titles: std::cell::RefCell<Vec<(usize, Vec<u8>)>>,
    /// The title of any other input record that the batch being delivered does not hold
    /// (`BatchTake::with_earlier_titles`): `None` where it cannot be read.
    pub(super) earlier: Option<&'a dyn Fn(usize) -> Result<Option<Vec<u8>>>>,
}

impl BatchTake<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:629-652 (cited at
    // `query_batch::PastFirstQuery`)
    /// Delivers the reader batch `queries` (after the queries of the batches delivered before it)
    /// to the loop (`deliver`) as the loop searches and reports it. Under `-parse_deflines` the
    /// batch joins the scopes of the queries first (`add_batch_scope`, which the loop's `take` of
    /// the batch then finds done): BLASTX's loop takes each batch at its delivery, and
    /// RPSTBLASTN's window takes it after the batches read behind it, so LOSAT's refusal of a
    /// batch comes there, after the reports of the batches before it, as it did. A searched query
    /// whose set-up NCBI reads past the end of its earlier record (`PastFirstQuery`) is delivered
    /// as NCBI's set-up holds it (`PastFirstQuery::set_up_parts`, nucleotide: BLASTX's and
    /// RPSTBLASTN's queries are): the earlier record's title, whose Bioseq NCBI's scope finds for
    /// the id, and a filler of its interval's length, whose contexts are invalid as NCBI's are.
    /// Every other batch is delivered as it was read, without a copy.
    pub(in crate::algorithm::blastx) fn deliver_views(
        &self,
        queries: &[FastaRecord],
        deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    ) -> Result<()> {
        // An empty batch is NCBI's `Empty CBlastQueryVector`, before it is taken.
        if queries.is_empty() || self.query_scopes.borrow().is_none() {
            return deliver(queries);
        }
        let (queries_before, input_before) = self.delivered.counts.get();
        add_batch_scope(
            self.query_scopes,
            &self.scope_taken,
            queries,
            queries_before,
            self.strategy.is_some(),
            self.options,
            self.subjects,
            self.search_error,
        )?;
        let inputs = queries.len() + queries[0].skipped.as_deref().map_or(0, <[_]>::len);
        self.delivered
            .counts
            .set((queries_before + queries.len(), input_before + inputs));
        let views = self.views(queries, queries_before, input_before)?;
        deliver(views.as_deref().unwrap_or(queries))
    }

    /// NCBI's view of the batch `queries` (the run's searched queries from `queries_before`, its
    /// input records from `input_before`): `None` where no query of the batch reads past its
    /// earlier record. The titles of the batch's records that the scopes' later queries show are
    /// kept first (`Delivered::titles`).
    fn views(
        &self,
        queries: &[FastaRecord],
        queries_before: usize,
        input_before: usize,
    ) -> Result<Option<Vec<FastaRecord>>> {
        let scopes = self.query_scopes.borrow();
        let end = queries_before + queries.len();
        self.delivered.keep_titles(
            queries,
            input_before,
            QueryScopes::past_first_in(scopes.as_ref(), end..usize::MAX),
        );
        let past_first = QueryScopes::past_first_in(scopes.as_ref(), queries_before..end);
        if past_first.is_empty() {
            return Ok(None);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:154-156
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:170-184
        // NCBI's error of a batch whose score block cannot load the `-matrix` of an ungapped
        // search holds each query's set-up message, this read's too, which LOSAT's error of the
        // batch, made of the records searched (`query_setup::unknown_matrix_message`), lacks.
        let program = self.options.program().upper;
        if self.options.unknown_matrix && self.options.matrix_file.is_none() {
            anyhow::bail!(
                "a query record under -parse_deflines that reads past the end of an earlier \
                 record with its Seq-id, in a search whose -matrix NCBI BLAST+ does not know, is \
                 not supported by LOSAT's {program}"
            );
        }
        let mut views = queries.to_vec();
        for &(query, past) in past_first {
            let title = self.delivered.title(queries, input_before, past, program)?;
            let view = &mut views[query - queries_before];
            past.set_up_parts(&mut view.title, &mut view.sequence, &title, true);
        }
        Ok(Some(views))
    }
}

impl Delivered<'_> {
    /// Keeps the titles of the input records of the batch `queries` (from `input_before`) that the
    /// queries `later` of later batches show (`PastFirstQuery::first_input`); nothing where there
    /// is none.
    fn keep_titles(
        &self,
        queries: &[FastaRecord],
        input_before: usize,
        later: &[(usize, PastFirstQuery)],
    ) {
        if later.is_empty() {
            return;
        }
        let mut titles = self.titles.borrow_mut();
        for &(_, past) in later {
            let kept = titles.iter().any(|&(input, _)| input == past.first_input);
            let title = past
                .first_input
                .checked_sub(input_before)
                .and_then(|at| title_in_batch(queries, at));
            if let (false, Some(title)) = (kept, title) {
                titles.push((past.first_input, title.to_vec()));
            }
        }
    }

    /// The title of the earlier record of the query `past` of the batch `queries` (whose input
    /// records start at `input_before`): a record of the batch, a kept one (`keep_titles`), or
    /// one read again (`earlier`); else LOSAT's refusal (`program`).
    fn title(
        &self,
        queries: &[FastaRecord],
        input_before: usize,
        past: PastFirstQuery,
        program: &str,
    ) -> Result<Vec<u8>> {
        let input = past.first_input;
        if let Some(title) = input
            .checked_sub(input_before)
            .and_then(|at| title_in_batch(queries, at))
        {
            return Ok(title.to_vec());
        }
        if let Some((_, title)) = self
            .titles
            .borrow()
            .iter()
            .find(|&&(kept, _)| kept == input)
        {
            return Ok(title.clone());
        }
        match self
            .earlier
            .map(|earlier| earlier(input))
            .transpose()?
            .flatten()
        {
            Some(title) => Ok(title),
            None => anyhow::bail!(
                "a query record under -parse_deflines reads past the end of an earlier record with \
                 its Seq-id, whose title NCBI BLAST+ shows for it, and LOSAT cannot read that \
                 record's title again from this query input (the standard input, a pipe or an \
                 imported search strategy's queries); this is not supported by LOSAT's {program}"
            ),
        }
    }
}

/// The title of the batch's input record `at` (`input::scope_records`), where the batch holds it:
/// a record that `-query_loc` skips keeps none (`input::SkippedRecord`).
fn title_in_batch(queries: &[FastaRecord], at: usize) -> Option<&[u8]> {
    let place = scope_records(queries).get(at)?.2?;
    Some(&queries[place].title)
}

// No NCBI counterpart: NCBI's scope holds every query record of the run, LOSAT's reader one
// batch at a time.
/// The title of the input record `input` (`PastFirstQuery::first_input`: the records of the query
/// input in the order read, those that `-query_loc` skips among them) of the command line's query
/// file `args.query`, read again from its start as the read before the search reads it
/// (`QueryStreamSettings::stream_before_search` of `settings`, in batches of `batch_size`
/// letters): the title of an earlier record that the batches no longer hold, which a query read
/// past it shows (`BatchTake::deliver_views`). `None` where the input cannot be read again: the
/// standard input, a file that is not a regular one (a pipe), an imported strategy's queries.
pub(super) fn reread_query_title(
    args: &BlastxArgs,
    settings: &QueryStreamSettings,
    batch_size: usize,
    input: usize,
) -> Result<Option<Vec<u8>>> {
    let imported = args
        .import
        .as_ref()
        .is_some_and(|run| run.queries.is_some());
    if imported
        || args.query.as_os_str() == "-"
        || !std::fs::metadata(&args.query).is_ok_and(|metadata| metadata.is_file())
    {
        return Ok(None);
    }
    let Ok(file) =
        crate::blastinput::input_files::open_input(&args.query, "query", settings.program.upper)
    else {
        return Ok(None);
    };
    let mut stream = settings.stream_before_search(file);
    let (mut read, mut title) = (0usize, None);
    // The read ends at the record, with an error of its own that is dropped.
    let _ = stream.read_batches(
        args.lcase_masking,
        batch_size.max(1),
        None,
        &mut |batch| {
            if let Some(record) = input.checked_sub(read).and_then(|at| batch.get(at)) {
                title = Some(record.title.clone());
                anyhow::bail!("the record is read");
            }
            read += batch.len();
            Ok(())
        },
        &mut |_| Ok(()),
    );
    Ok(title)
}
