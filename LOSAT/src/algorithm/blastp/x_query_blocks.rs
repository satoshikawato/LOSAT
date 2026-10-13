//! EXPERIMENT (LOSAT_X_BPQBLOCK): blastp searches its queries in blocks of whole query batches,
//! each block a search of its own (its lookup table, its diagonal table, its hit lists, its
//! rows), in query order.
//!
//! The reference path searches every query of the batches that NCBI does not split at once:
//! one lookup table and one diagonal table for all queries, then the hit lists of all queries
//! until the rows are written. The lookup table of a whole proteome holds about 35 entries of
//! 4 bytes per query letter (7.5 million letters: about 1 GB, 2.7 GB while it is built), and the
//! diagonal table of 2^23 entries does not fit in a cache. NCBI reads one batch of about 10,000
//! letters at a time and searches it (`blastp_app.cpp` below), so its tables stay small.
//!
//! Why the rows are the reference's: a blastp query's hit list does not depend on the other
//! queries of the search (`blast_engine.rs` `search_split_query_batches`), so a block of whole
//! queries finds the hit lists that the search of all queries finds for them:
//! - the two-hit test of a hit reads and writes only the entry of its diagonal. Two queries
//!   share an entry when a diagonal runs from the end of one query into the next one, or when
//!   two diagonals are `diag_mask + 1` apart. In the first case an extension stops at the
//!   sentinel between the queries (the sentinel is NCBISTDAA 0, `-`, whose scores are
//!   `BLAST_SCORE_MIN`, so the walk breaks there), and the next query's first hit on the
//!   diagonal finds either a flag whose `last_hit` lies before the sentinel (a new hit), a
//!   `last_hit` at least `window` letters back (a new hit), or a `last_hit` of the other query
//!   at least one word and the sentinel back, which the context check turns into a new hit: in
//!   every case the entry becomes (this hit, no flag), as on a table that never saw the other
//!   query. In the second case the later hit is more than `window` letters further along the
//!   subject than the earlier one, and an extension on the earlier diagonal ends before the
//!   query set does, so the later hit again starts a new hit.
//! - the gapped stage keeps its state per query (`found_high_score`, `restricted_align_array`,
//!   the redo of one query) and the interval tree holds HSPs in disjoint query ranges.
//! - the three things that do depend on NCBI's batches are the batches' own: the split of a
//!   long batch, the heaps that the composition-based redo of a query reads
//!   (BlastCompo_EarlyTermination: the queries of its batch), and a batch of invalid queries
//!   only. A block is made of whole batches, and the batches inside a block are those of the
//!   whole input (they are counted greedily from the block's first record, which starts a
//!   batch of the whole input).
//! The rows of tabular format 6 depend only on each query's hit list and come in query order,
//! so the blocks' rows written one after the other are the rows of the whole search. The block
//! search is used only where nothing else is written around the rows (see
//! `blast_engine.rs` `search`); everything else takes the reference path.
//!
//! Off unless `LOSAT_X_BPQBLOCK` is set to a positive number of query letters (the block closes
//! at the first batch end at or past it). Not in the strict set.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:562-572
//! ```c
//!                 /* Check if the last hit hits current query. Because last_hit
//!                    is never reset, it may contain a hit to the previous
//!                    concatenated query */
//!
//!                 if (query_offset - diff <
//!                     query_info->contexts[curr_context].query_offset) {
//!
//!                     /* there was no last hit for this diagnol; start a new hit */
//!                     diag_array[diag_coord].last_hit =
//!                         subject_offset + diag_offset;
//!                     continue;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1575
//! ```c
//!     /* Initialize with BLAST_SCORE_MIN */
//!     for (i = 0; i < sbp->alphabet_size; i++) {
//!         for (j = 0; j < sbp->alphabet_size; j++) {
//!             matrix[i][j] = BLAST_SCORE_MIN;
//!     ...
//!                 i == AMINOACID_TO_NCBISTDAA['-'] ||
//!     ...
//!                 continue;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:859-860
//! ```c
//!         if (score <= 0 || (maxscore - score) >= dropoff)
//!             break;
//! ```

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::ops::Range;
use std::path::Path;

use crate::api::local_blast::OutputSink;
use crate::blastinput::seq_range::RangedBatch;

// No NCBI counterpart: reads the LOSAT_X_BPQBLOCK switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the number of query letters of a block, or `None` (off) unless
/// the variable is a positive integer.
pub(crate) fn x_bp_qblock() -> Option<usize> {
    use std::sync::OnceLock;
    static N: OnceLock<Option<usize>> = OnceLock::new();
    *N.get_or_init(|| {
        std::env::var("LOSAT_X_BPQBLOCK")
            .ok()
            .and_then(|value| value.trim().parse::<usize>().ok())
            .filter(|&letters| letters > 0)
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-287
// ```c++
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//             ...
//                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
//                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
//                 results = lcl_blast.Run();
//             ...
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
// NCBI searches and reports one batch at a time; a block is one or more whole batches.
/// The blocks of searched queries: consecutive batches (`batches`, of the whole input), each
/// block closed at the first batch end where its letters (`lengths` of the searched queries)
/// reach `block_letters`. `None` when a batch has no searched query or the batches do not cover
/// the searched queries in order.
pub(crate) fn x_query_blocks(
    batches: &[RangedBatch],
    lengths: &[usize],
    block_letters: usize,
) -> Option<Vec<Range<usize>>> {
    let mut blocks = Vec::new();
    let mut start = 0usize;
    let mut end = 0usize;
    let mut letters = 0usize;
    for batch in batches {
        if batch.searched.is_empty() || batch.searched.start != end {
            return None;
        }
        end = batch.searched.end;
        letters += lengths.get(batch.searched.clone())?.iter().sum::<usize>();
        if letters >= block_letters {
            blocks.push(start..end);
            start = end;
            letters = 0;
        }
    }
    if end != lengths.len() {
        return None;
    }
    if start < end {
        blocks.push(start..end);
    }
    Some(blocks)
}

/// The output of the block search: the one sink of the format, opened once for all blocks.
/// A file is created at its first write or flush, which the first block's report makes, as
/// the reference path creates it when its report is written.
pub(crate) enum XBlockSink<'a> {
    Stdout(BufWriter<io::Stdout>),
    File {
        path: &'a Path,
        file: Option<BufWriter<File>>,
    },
    Writer(&'a mut (dyn Write + Send)),
}

impl<'a> XBlockSink<'a> {
    // No NCBI counterpart: the destination of `OutputSink::open`, kept open across the blocks.
    pub(crate) fn new(sink: &'a mut OutputSink<'_>) -> Self {
        match sink {
            OutputSink::Stdout => Self::Stdout(BufWriter::new(io::stdout())),
            OutputSink::File(path) => Self::File { path, file: None },
            OutputSink::Writer(writer) => Self::Writer(&mut **writer),
        }
    }

    fn writer(&mut self) -> io::Result<&mut dyn Write> {
        Ok(match self {
            Self::Stdout(stdout) => stdout,
            Self::File { path, file } => {
                if file.is_none() {
                    *file = Some(BufWriter::new(File::create(*path)?));
                }
                file.as_mut().expect("file created above")
            }
            Self::Writer(writer) => &mut **writer,
        })
    }
}

impl Write for XBlockSink<'_> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.writer()?.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer()?.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn batch(input: Range<usize>, searched: Range<usize>) -> RangedBatch {
        RangedBatch { input, searched }
    }

    #[test]
    fn blocks_close_at_the_first_batch_end_past_the_size() {
        let lengths = [6, 6, 3, 9, 2, 2, 7];
        let batches = [
            batch(0..2, 0..2),
            batch(2..4, 2..4),
            batch(4..6, 4..6),
            batch(6..7, 6..7),
        ];
        // Letters per batch: 12, 12, 4, 7.
        assert_eq!(
            x_query_blocks(&batches, &lengths, 20),
            Some(vec![0..4, 4..7])
        );
        assert_eq!(
            x_query_blocks(&batches, &lengths, 12),
            Some(vec![0..2, 2..4, 4..7])
        );
        assert_eq!(
            x_query_blocks(&batches, &lengths, 1),
            Some(vec![0..2, 2..4, 4..6, 6..7])
        );
        assert_eq!(x_query_blocks(&batches, &lengths, 1000), Some(vec![0..7]));
    }

    #[test]
    fn batches_that_do_not_cover_the_queries_give_no_blocks() {
        let lengths = [5, 5, 5];
        assert_eq!(
            x_query_blocks(&[batch(0..2, 0..2), batch(2..3, 2..2)], &lengths, 1),
            None
        );
        assert_eq!(x_query_blocks(&[batch(0..2, 0..2)], &lengths, 1), None);
        assert_eq!(
            x_query_blocks(&[batch(0..1, 0..1), batch(1..3, 2..3)], &lengths, 1),
            None
        );
        assert_eq!(x_query_blocks(&[], &[], 1), Some(vec![]));
    }

    #[test]
    fn a_file_is_created_at_the_first_write_or_flush() {
        let dir = std::env::temp_dir().join(format!("losat-x-qblock-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("out.tsv");
        let _ = std::fs::remove_file(&path);
        let mut output = OutputSink::File(&path);
        {
            let mut sink = XBlockSink::new(&mut output);
            assert!(!path.exists());
            sink.flush().unwrap();
            assert!(path.exists());
            sink.write_all(b"a\t1\n").unwrap();
            sink.write_all(b"b\t2\n").unwrap();
            sink.flush().unwrap();
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"a\t1\nb\t2\n");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
