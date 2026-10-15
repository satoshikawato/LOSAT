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
//! With `-db` the kernel runs over the database's subjects one after the other on one diagonal
//! table (`blast_engine/database.rs` `search_blocks`, whose blocks start at the one-thread
//! offset). From one subject to the next NCBI's offset grows by the subject's length plus
//! `window`, and the table is cleared where the offset reaches INT4_MAX / 4
//! (`Blast_ExtendWordExit` below), so an entry left by an earlier subject lies at least `window`
//! letters behind every hit of a later one: without its flag it gives `diff >= window`, with its
//! flag the hit is not behind it; in both cases the entry becomes (this hit, no flag), whatever
//! query wrote it, as a cleared entry would. A subject's hits therefore depend on that subject
//! only, and the argument above holds subject by subject. The batches of `-db` are those of
//! `-subject`. A database search can also fail after the rows of some batches where it reads the
//! definition lines of the reported subjects; the block search runs only when no definition line
//! of the database can fail that way (`database::x_titles_allow_blocks`).
//!
//! With the compressed lookup table of blastp-fast (`CompressedAaLookupTable`), the hits of a
//! block's scan reach the two-hit test in another order within a subject offset, and that order
//! changes nothing:
//! - a cell of the block's table holds the offsets of the block's queries, added in the order of
//!   the whole table (the order of the queries' offsets, shifted by the block's start), but the
//!   scan reads a cell of more than 5 offsets as NCBI does: the first three, then the others in
//!   groups of four, the group added last first (`query_offsets_in_scan_order`; blast_aalookup.c
//!   and blast_aascan.c below). Where the groups fall depends on the cell's count, so a block's
//!   hits at one subject offset do not always come in the whole search's order.
//! - the scan stops where the offset-pair array cannot take the next cell's hits and resumes at
//!   that subject offset. The array's size (`OFFSET_ARRAY_SIZE` plus the table's longest chain)
//!   and the offsets where the scan stops depend on the other queries' hits, but a cell is never
//!   split (the array takes the longest cell), and the resumed scan reads the same word there (a
//!   word that hit has no letter outside the compressed alphabet, so the priming stops at once,
//!   and its index is the one the rolling index gives). So both searches give the block's
//!   queries the same hits per subject offset, the offsets in increasing order.
//! - the hits of one subject offset lie on different diagonal entries: their query offsets differ
//!   by less than the concatenated query's length, and the diagonal table has at least that many
//!   entries (blast_extend.c below). The test of a hit reads and writes only its entry and reads
//!   no other state that a hit changes, so their order changes no test and no extension, only the
//!   order in which their ungapped alignments are saved.
//! - those alignments are sorted before anything reads them (`ncbi_qsort_init_hsps_by_score`,
//!   stable as glibc's qsort, by score, subject start, length and query start). Two alignments
//!   that tie have the same subject start and query start, so the same diagonal, and come from
//!   hits at different subject offsets, saved in the same order in both searches; the sorted list
//!   is the reference's (its query starts shifted by the block's start).
//! The same argument covers the reference path's one table of all queries against NCBI's table
//! per query batch. A block's table is tested for NCBI's overflow banks with the block's own
//! batches (`blast_engine.rs` `ncbi_compressed_lookup_tables`), and the block search is taken only
//! where the words of each batch's offsets show that no batch's table fills them
//! (`x_block_words_fit`, tested before the first block), so that the reference path's refusal,
//! which comes before any row, is never made in a later block.
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
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:167-173
//! ```c
//!     if (ewp->diag_table) {
//!         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//!             ewp->diag_table->offset = ewp->diag_table->window;
//!             s_BlastDiagClear(ewp->diag_table);
//!         } else {
//!             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//!         }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:831-840
//! ```c
//!         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
//!
//!         if (cell_index == 0 ) { /* can't be empty => it's full  */
//!
//!             /* fetch next vacant cell */
//!         	CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
//!
//!             /* shuffle the pointers */
//!             new_cell->next = backbone_cell->payload.overflow_list.head;
//!             backbone_cell->payload.overflow_list.head = new_cell;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:274-284
//! ```c
//!     /* prime the index */
//!     for(s = s_first; s <= s_last; s++){
//!         index = s_ComputeCompressedIndex(word_length - 1, s,
//!                                          compressed_alphabet_size,
//!                                          &skip, lookup);
//!         if(!skip)
//!           break;
//!     }
//!
//!     next_char = ((s <= s_last)? s[word_length-1] : 0);
//!     preshift = (Int4)((((Int8)index) * recip) >> 32);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:384-405
//! ```c
//!                 /* handle the rest of the list */
//!                 if (first_cell_entries) {
//!                    curr_cell = curr_cell->next;
//!                 }
//!                    while (curr_cell != NULL) {
//!     ...
//!              totalhits += numhits;
//!           }
//!           else
//!               /* not enough space in the destination array */
//!           {
//!               s_range[1] = (Int4)(s - subject->sequence);
//!               return totalhits;
//!           }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-56
//! ```c
//!                 diag_array_length = 1;
//!                 /* What power of 2 is just longer than the query? */
//!                 while (diag_array_length < (qlen+window_size))
//!                 {
//!                         diag_array_length = diag_array_length << 1;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:234
//! ```c
//!     Blast_InitHitListSortByScore(init_hitlist);
//! ```

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::ops::Range;
use std::path::Path;

use super::encoding::EncodedProtein;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::OutputSink;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::utils::threading::SearchPool;

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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
// ```c
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
//        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
//            continue;
//        }
// ```
// NCBI reads each subject again for each query batch; the encoded subject is a function of the
// subject's record only, so one encoding serves every block. Reuse only.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the subjects encoded once for all blocks (attached to the search's
/// pool, `SearchPool::x_with_ext`), with the records they encode.
pub(crate) struct XBlockSubjects {
    /// The address and count of the records encoded.
    records: (usize, usize),
    encoded: Vec<EncodedProtein>,
}

impl XBlockSubjects {
    // No NCBI counterpart: keeps the encoding of `records` (module comment above).
    pub(crate) fn new(records: &[FastaRecord], encoded: Vec<EncodedProtein>) -> Self {
        assert_eq!(records.len(), encoded.len());
        Self {
            records: (records.as_ptr() as usize, records.len()),
            encoded,
        }
    }
}

// No NCBI counterpart: reads back the encoding attached by the block search (reuse only).
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the encoded subjects attached to `pool`, when they encode
/// `records` (the same slice).
pub(crate) fn x_block_subjects<'p>(
    pool: &SearchPool<'p>,
    records: &[FastaRecord],
) -> Option<&'p [EncodedProtein]> {
    pool.x_ext::<XBlockSubjects>()
        .filter(|subjects| subjects.records == (records.as_ptr() as usize, records.len()))
        .map(|subjects| &subjects.encoded[..])
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803 (the contexts' Karlin
// blocks, `x_prepare_blosum62_lookup_query_in_pool`)
// ```c
//       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ```
/// EXPERIMENT (LOSAT_X_BPQBLOCK): `prepare_blosum62_lookup_query_for_word_size(query_frames,
/// word_size, false)`, with the contexts' Karlin blocks computed on `pool`'s threads when the switch
/// is set and the pool has more than one thread (the same values).
pub(crate) fn x_bp_lookup_query(
    query_frames: &[Vec<QueryFrame>],
    word_size: usize,
    pool: &SearchPool<'_>,
) -> (Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>) {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_bp_qblock().is_some() && pool.enabled() && pool.threads() > 1 {
        return crate::algorithm::tblastx::lookup::x_prepare_blosum62_lookup_query_in_pool(
            query_frames,
            word_size,
            false,
            pool,
        );
    }
    let _ = pool;
    crate::algorithm::tblastx::lookup::prepare_blosum62_lookup_query_for_word_size(
        query_frames,
        word_size,
        false,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803
    // ```c
    //       if ( !contexts[context].is_valid )
    //           continue;
    //       ...
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // ```
    // The lookup query with the contexts' Karlin blocks computed on a pool of 1, 3 and 7 threads is
    // the reference's (`prepare_blosum62_lookup_query_for_word_size`): the query, the segments and
    // every context (Karlin blocks bit for bit, validity). Random queries with low-complexity runs
    // (SEG-masked, several segments), a query of X only (no Karlin block: invalid), a short one.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn lookup_query_with_karlin_blocks_on_a_pool_is_the_references() {
        use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
        use crate::algorithm::tblastx::lookup::x_pair_lut::x_assert_same_contexts;
        use crate::algorithm::tblastx::lookup::{
            prepare_blosum62_lookup_query_for_word_size, x_prepare_blosum62_lookup_query_in_pool,
        };
        let letters = b"ACDEFGHIKLMNPQRSTVWYBZXUO*";
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut below = |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        for _ in 0..40 {
            let length = 1 + below(400);
            let mut sequence: Vec<u8> = (0..length).map(|_| letters[below(20)]).collect();
            if length > 60 && below(2) == 0 {
                let at = below(length - 40);
                let run = letters[below(letters.len())];
                sequence[at..at + 30].fill(run);
            }
            if below(10) == 0 {
                sequence[below(length)] = letters[20 + below(6)];
            }
            sequences.push(sequence);
        }
        sequences.push(b"XXXXXXXXXXXX".to_vec());
        sequences.push(b"MK".to_vec());
        let seg = crate::utils::seg::SegParams::default();
        for seg in [None, Some(&seg)] {
            let frames: Vec<Vec<QueryFrame>> = sequences
                .iter()
                .map(|sequence| vec![encode_protein_query_frame_with_seg(sequence, seg)])
                .collect();
            let (query, locations, contexts) =
                prepare_blosum62_lookup_query_for_word_size(&frames, 5, false);
            assert!(contexts.iter().any(|context| !context.is_valid));
            for threads in [1, 3, 7] {
                let (pool_query, pool_locations, pool_contexts) =
                    crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                        Ok(x_prepare_blosum62_lookup_query_in_pool(
                            &frames, 5, false, pool,
                        ))
                    })
                    .expect("pool");
                assert_eq!(pool_query, query);
                assert_eq!(pool_locations, locations);
                x_assert_same_contexts(&pool_contexts, &contexts);
            }
        }
    }

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
