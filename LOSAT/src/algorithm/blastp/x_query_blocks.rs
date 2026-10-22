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
//! search is used only where the blocks' reports written one after the other are the report of
//! the whole search (see `blast_engine.rs` `search` and `x_query_block_plan`); everything else
//! takes the reference path.
//!
//! LOSAT_X_BPQBLOCKFMT (strict set, `x_bp_qblock_fmt`) also takes outfmt 0 and 7, whose reports
//! hold more than the rows (NCBI prints the prolog, then each batch's result sets in its batch
//! loop, then the epilog: `blastp_app.cpp` below):
//! - each query's part depends only on that query's final list and on fixed inputs: outfmt 7's
//!   header (program, version, the query's title, the `# Database:` name, the query's hit count)
//!   and outfmt 0's section (the query's title and length, its description table and alignments,
//!   the subjects' titles, the footer with the Karlin-Altschul block of the query's context and its
//!   search space, its validity and whether its batch was searched: a batch property, and blocks
//!   are whole batches). With `-db` a reported sequence is numbered by its rank among the
//!   sequences that the call's lists name; the rank only indexes the call's records and titles,
//!   so a rank among one block's sequences writes the same bytes as one among all. A query's
//!   warnings come before its part in both (`QueryWarnings`, the block's reader messages).
//! - what is not per query: the outfmt 0 prolog (the program, the `Database:` header of the
//!   subject set or the database: fixed inputs) is written by the first block only
//!   (`BlastpReports::x_prolog`); the outfmt 0 epilog (the database lines, the gapped
//!   Karlin-Altschul and Gumbel blocks of the scoring and the database length, the options: fixed
//!   inputs) by the last block only, as before; outfmt 7's last line counts the queries formatted
//!   (NCBI's `m_QueriesFormatted`), so the last block adds the queries of the blocks before
//!   (`BlastpReports::x_queries_before`).
//! - the refusals of the reference path that come before it writes anything and depend on the
//!   hits are, for outfmt 0, those of a shown subject's title (a title of punctuation that NCBI's
//!   `x_CleanAndCompress` reads past, ABI v1's titles): the blocks are searched only when no
//!   subject's title is refused (`blast_engine.rs` `x_block_shown_titles_allowed`, with `-db`
//!   `database::x_titles_allow_blocks`), so that none comes after the reports of the blocks
//!   before. The other refusals are those of outfmt 6 above.
//! - nothing else is written to standard error while a block is searched (only the debug trace and
//!   `LOSAT_TIMING`'s lines, which each block prints for itself).
//! The shadow (LOSAT_X_BPQBLOCKFMTSHADOW, implies the switch) keeps a copy of the bytes that the
//! blocks write to the output and to standard error (`XCopyWriter`), then searches all queries
//! again, writes their report and warnings into memory and compares both streams
//! (`x_block_shadow_compare`), for every format the block search takes.
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
//! LOSAT_X_BPQBLOCKSA (strict set, `x_bp_qblock_sa`, outfmt 8, 9 and 12) and LOSAT_X_BPQBLOCKXML
//! (strict set, `x_bp_qblock_xml`, outfmt 5) take the other formats of one blastp report; both act
//! only with `LOSAT_X_BPQBLOCK`, and the shadow implies both:
//! - outfmt 8, 9 and 12 (Seq-annot of Seq-aligns as ASN.1 text, BER and JSON): NCBI writes one
//!   self-contained Seq-annot per query that has alignments, at once, a new object stream for each,
//!   and writes no prolog and no epilog for them. LOSAT's writer keeps nothing from one query to the
//!   next (new encoder per query); what it reads is the query's final list and fixed inputs (the
//!   program, the database title, the ids `Query_N` and `Subject_N` of the records, the oid of a
//!   `-db` subject). So the blocks' documents written one after the other are the documents of
//!   the whole search. The refusals that depend on the hits are those of a `-db` definition line
//!   (`database::x_titles_read_allow_blocks`: every sequence has one definition line that reads and
//!   holds no byte outside ' ' ..= '~'); a `-subject` search reads no title, and an encoder's real
//!   (e-value, bit score) is finite, so the reference path fails nowhere that the blocks do not.
//! - outfmt 5 (BLAST XML): NCBI makes one formatter for the run; its first formatted query writes
//!   the header (from that query), each query appends an Iteration numbered over the run, and the
//!   epilog writes the end once. The writer takes the count of the queries before its block
//!   (`BlastpReports::x_queries_before`, which starts the header and `iter-num`) and the epilog
//!   flag (the last block), as outfmt 7's count does; the run values (program, version,
//!   `BlastOutput_db`, matrix, gap costs, filter, `db-num`, `db-len`, kappa, lambda, entropy) are
//!   fixed inputs; `hsp-len`, `eff-space`, the masks and the message are the query's own. The
//!   writer refuses what NCBI cannot write (a byte >= 0x80 or a control byte in a name, a title it
//!   cannot guess the encoding of or reads past the end of) before its first byte, which in a later
//!   block would come after the earlier blocks' bytes. So the block search is taken only when
//!   none is refused, checked with the writer's own functions before the first block: the database
//!   name and every query's definition line and message id (`x_block_xml_run_allowed`), every
//!   `-subject` record's names (`x_block_xml_subjects_allowed`; conservative, the writer refuses
//!   only a subject with hits), and with `-db` the definition lines as for outfmt 8
//!   (`database::x_titles_read_allow_blocks`). Otherwise the search of all queries refuses where
//!   it does today; the file sink is not created before then.
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
//! ```c++
//!         if (results.HasAlignments()) {
//!     ...
//!             if(m_FormatType == CFormattingArgs::eAsnText)
//!             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:360-368,2232-2246
//! ```c++
//!     // no header for some output types
//!     if (m_FormatType >= CFormattingArgs::eXml) {
//!     ...
//!     } else if (m_FormatType >= CFormattingArgs::eTabular)
//!         return;  // No footer for these.
//!     ...
//!         m_Outfile << m_BlastXMLIncremental->m_SerialXmlEnd << endl;
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:725-744
//! ```c++
//!     if( incremental_output && incremental_struct->m_IterationNum == 0) {
//!     ...
//!         if (incremental_struct)
//!         	incremental_struct->m_IterationNum++;
//! ```
//!
//! LOSAT_X_BPQBLOCKWEB (strict set, `x_bp_qblock_web`; acts only with `LOSAT_X_BPQBLOCK`, and the
//! shadow implies it) takes a `run_local` call with several outputs, the hit records or an observer
//! (ABI v2 of the Web build asks for outfmt 0, 6 and 7, the records and the byte ranges of each HSP
//! in each output), when every output of several is a writer of the caller, no `-db` is searched and
//! an observer comes with the records. Each block's call writes all the outputs, in the order asked:
//! - each output's bytes are those of the one-output case above: what is not per query (the prolog,
//!   the epilog, the counts of outfmt 7 and 5) is set per block for every output alike
//!   (`BlastpReports`), each output's writer writes only to its own sink, and with several outputs
//!   the sinks are the caller's writers, unbuffered (a stream that two outputs share would get each
//!   block's formats in turn, so several outputs need writers; one output may also be standard
//!   output or a file, as above). The checks of the refusals that depend on the hits are made for
//!   every kind of report among the outputs before the first block. The warnings go with the first
//!   output, query by query, in every block.
//! - the final list of all queries is the blocks' lists one after the other: it is ordered by query,
//!   then by subject, then by HSP (`hsp.rs`'s final sort), and a query's list does not depend on the
//!   other queries (above). A record's `q_idx` is its block's first query plus its index in the
//!   block; `s_idx` and the subject's length come from the call's one subject set (`-subject` only:
//!   with `-db` they would rank the subjects that the call's lists name); no other field holds a
//!   value of the whole search. The block search hands the records of all blocks to the caller's
//!   sink once, after the last block, as the search of all queries hands its one list once (before
//!   its reports; no caller reads the records before the reports).
//! - an observer event's index is that of the HSP in the final list: the block's index plus the
//!   number of HSPs of the blocks before (`XOffsetObserver`, counted from the block's records, which is
//!   why an observer needs them); its position is read by the caller from its own writer, whose bytes
//!   are those of the search of all queries. Each format's events come in the same order as in the
//!   search of all queries; between formats they alternate block by block, and an observer reads
//!   each format's events apart (the ABI orders the events within a stream only).
//! The shadow also compares, with these, the records (`Debug` text) and each format's events with
//! their positions in that format's bytes (`XShadowObserver`).
//!
//! `LOSAT_X_BPQBLOCK` is a value switch of the strict set (`utils/xswitch.rs`): a positive number
//! of query letters (the block closes at the first batch end at or past it), 50,000 when it is
//! not set or is not a positive number. `0`, or `LOSAT_X_DEFAULT=0` with the variable not set,
//! turns the block path off.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-295
//! ```c++
//!         formatter.PrintProlog();
//!     ...
//!                     formatter.PrintOneResultSet(**result, query_batch);
//!     ...
//!         formatter.PrintEpilog(opt);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1430,2235
//! ```c++
//!     m_QueriesFormatted++;
//!     ...
//!         tabinfo.PrintNumProcessed(m_QueriesFormatted);
//! ```
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
use std::sync::Mutex;

use super::encoding::EncodedProtein;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_bpclut_mode, x_bplkpflat_mode, x_pool_builds, XFlatArgs, XFlatSpare, XLkpFlat,
};
use crate::algorithm::tblastx::lookup::{QueryContext, XAaSpares};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::{FormatObserver, HspIndex, OutputSink};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::utils::threading::SearchPool;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
        x_lkpflat_count_contexts, XExtraTasks,
    },
    crate::algorithm::tblastx::lookup::x_pair_lut::x_assert_same_contexts_as,
    crate::algorithm::tblastx::lookup::XBlosum62LookupQuery,
    crate::stats::KarlinParams,
    std::sync::OnceLock,
};

// No NCBI counterpart: reads the LOSAT_X_BPQBLOCK switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the number of query letters of a block, or `None` (off). A
/// strict-set value switch: the positive integer given, else 50,000; `None` for `0`, or for
/// `LOSAT_X_DEFAULT=0` with the variable not set.
pub(crate) fn x_bp_qblock() -> Option<usize> {
    use std::sync::OnceLock;
    static N: OnceLock<Option<usize>> = OnceLock::new();
    *N.get_or_init(|| crate::utils::xswitch::value("LOSAT_X_BPQBLOCK", 1..=usize::MAX, 50_000))
}

// No NCBI counterpart: reads LOSAT_X_BPQBLOCKFMT (a strict-set switch: on unless it is `0`, or
// `LOSAT_X_DEFAULT=0` and it is not set) and its shadow LOSAT_X_BPQBLOCKFMTSHADOW (off unless set to
// any value; it implies the switch) once. They choose which reports the block search writes, not
// any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCKFMT): 0 = the block search only for one outfmt 6 report, 1 = for one
/// report of outfmt 0, 6 or 7 (with LOSAT_X_BPQBLOCKWEB also for several reports of these and the
/// others the block search takes), 2 = the shadow (`x_block_shadow_compare`).
pub(crate) fn x_bp_qblock_fmt() -> u8 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| {
        if std::env::var_os("LOSAT_X_BPQBLOCKFMTSHADOW").is_some() {
            2
        } else if crate::utils::xswitch::on("LOSAT_X_BPQBLOCKFMT") {
            1
        } else {
            0
        }
    })
}

// No NCBI counterpart: reads LOSAT_X_BPQBLOCKSA and LOSAT_X_BPQBLOCKXML (strict-set switches: on
// unless `0`, or `LOSAT_X_DEFAULT=0` and not set) once; the shadow (LOSAT_X_BPQBLOCKFMTSHADOW) implies
// both, as it implies LOSAT_X_BPQBLOCKFMT. They choose which reports the block search writes, not any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCKSA): whether the block search also takes one Seq-annot report
/// (outfmt 8, 9 or 12). Acts only with `LOSAT_X_BPQBLOCK`.
pub(crate) fn x_bp_qblock_sa() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var_os("LOSAT_X_BPQBLOCKFMTSHADOW").is_some()
            || crate::utils::xswitch::on("LOSAT_X_BPQBLOCKSA")
    })
}

/// EXPERIMENT (LOSAT_X_BPQBLOCKXML): whether the block search also takes one BLAST XML report
/// (outfmt 5). Acts only with `LOSAT_X_BPQBLOCK`.
pub(crate) fn x_bp_qblock_xml() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var_os("LOSAT_X_BPQBLOCKFMTSHADOW").is_some()
            || crate::utils::xswitch::on("LOSAT_X_BPQBLOCKXML")
    })
}

// No NCBI counterpart: reads LOSAT_X_BPQBLOCKWEB (a strict-set switch: on unless `0`, or
// `LOSAT_X_DEFAULT=0` and not set) once; the shadow (LOSAT_X_BPQBLOCKFMTSHADOW) implies it. It
// chooses which calls the block search takes, not any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCKWEB): whether the block search also takes a `run_local` call with
/// more than one output, the hit records or an observer (ABI v2 of the Web build: outfmt 0, 6 and
/// 7, the records and the byte ranges), each output a writer of the caller and no `-db`. Acts only
/// with `LOSAT_X_BPQBLOCK`.
pub(crate) fn x_bp_qblock_web() -> bool {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| {
        std::env::var_os("LOSAT_X_BPQBLOCKFMTSHADOW").is_some()
            || crate::utils::xswitch::on("LOSAT_X_BPQBLOCKWEB")
    })
}

// No NCBI counterpart: the reports the block search writes; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): which report the block search writes, which the driver's checks
/// before the first block depend on.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum XBlockReport {
    /// Outfmt 6 and 7.
    Tabular,
    /// Outfmt 0 (LOSAT_X_BPQBLOCKFMT).
    Pairwise,
    /// Outfmt 8, 9 and 12 (LOSAT_X_BPQBLOCKSA).
    Seqalign,
    /// Outfmt 5 (LOSAT_X_BPQBLOCKXML).
    Xml,
}

// No NCBI counterpart: names which switches decide the block search of a call for `x_reach!`; it
// does not change any value NCBI computes.
/// A switch besides LOSAT_X_BPQBLOCK that governs a search whose block plan is made
/// (`LOSAT_X_REACH`): the format switch of each output (outfmt 6 has none), and for a call of
/// several outputs, the hit records or an observer the switch of `run_local`. The plan's refusals
/// before it looks at a call's outputs record none.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum XBlockFormat {
    /// Outfmt 0 and 7 (LOSAT_X_BPQBLOCKFMT).
    Fmt,
    /// Outfmt 8, 9 and 12 (LOSAT_X_BPQBLOCKSA).
    Seqalign,
    /// Outfmt 5 (LOSAT_X_BPQBLOCKXML).
    Xml,
    /// Several outputs, the hit records or an observer (LOSAT_X_BPQBLOCKWEB).
    Web,
}

// No NCBI counterpart: records which side the block-search switches took; it does not change any
// value NCBI computes.
/// `LOSAT_X_REACH`: the search took the block search (`own`) or the search of all queries, for
/// LOSAT_X_BPQBLOCK and each switch of `formats`, each while it is on.
pub(crate) fn x_reach_block_search(formats: &[XBlockFormat], own: bool) {
    use crate::utils::xreach::x_reach;
    if x_bp_qblock().is_none() {
        return;
    }
    x_reach!("LOSAT_X_BPQBLOCK", own);
    for format in formats {
        match format {
            XBlockFormat::Fmt if x_bp_qblock_fmt() > 0 => x_reach!("LOSAT_X_BPQBLOCKFMT", own),
            XBlockFormat::Seqalign if x_bp_qblock_sa() => x_reach!("LOSAT_X_BPQBLOCKSA", own),
            XBlockFormat::Xml if x_bp_qblock_xml() => x_reach!("LOSAT_X_BPQBLOCKXML", own),
            XBlockFormat::Web if x_bp_qblock_web() => x_reach!("LOSAT_X_BPQBLOCKWEB", own),
            _ => {}
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:725-733
// ```c++
//     if( incremental_output && incremental_struct->m_IterationNum == 0) {
//     ...
//         *out_stream << serial_xml_start << "\n<BlastOutput_iterations>" ;
// ```
// The refusals of the XML report (`write_blastp_xml`) come before its first byte; a block would
// make them after the bytes of the blocks before. These are the writer's own checks, made for the
// whole run before the first block (`names::check_text`, `names::xml_query_title`,
// `query_warnings::warning_query_id`, `names::xml_subject_names`); when one fails the search of
// all queries runs and refuses where it does today.
/// EXPERIMENT (LOSAT_X_BPQBLOCKXML): whether the database name and every query's definition line and
/// message id can be written by the XML report. The message of a query holds only its warning id
/// (`warning_query_id`) besides ASCII text. The database name check is defensive: the CLI refuses
/// the name in block 0 before any byte with the same output, and `run_local` takes no `-db`; it
/// would bind if `run_local` ever searched a database (its file is created at `open()`).
pub(crate) fn x_block_xml_run_allowed(
    db_name: &str,
    queries: &[FastaRecord],
    pool: &SearchPool<'_>,
) -> bool {
    use crate::report::blast_xml::names;
    names::check_text(db_name.as_bytes()).is_ok()
        && x_all_in_groups(queries.len(), pool, |range| {
            queries[range].iter().all(|record| {
                names::xml_query_title(&record.title, true).is_ok()
                    && names::check_text(&crate::report::query_warnings::warning_query_id(
                        record.local_id.as_bytes(),
                        &record.title,
                    ))
                    .is_ok()
            })
        })
}

/// EXPERIMENT (LOSAT_X_BPQBLOCKXML): whether the names of every `-subject` record can be written
/// by the XML report (`names::xml_subject_names`; the report refuses only a subject with hits,
/// so this is conservative).
pub(crate) fn x_block_xml_subjects_allowed(
    subjects: &[FastaRecord],
    pool: &SearchPool<'_>,
) -> bool {
    x_all_in_groups(subjects.len(), pool, |range| {
        subjects[range].iter().all(|record| {
            crate::report::blast_xml::names::xml_subject_names(
                record.local_id.as_bytes(),
                &record.title,
                true,
            )
            .is_ok()
        })
    })
}

// No NCBI counterpart: runs a check over groups of consecutive records, on the pool's threads
// when it has more than one; the outcome is whether all pass.
fn x_all_in_groups(
    count: usize,
    pool: &SearchPool<'_>,
    allowed: impl Fn(Range<usize>) -> bool + Sync,
) -> bool {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() && pool.threads() > 1 {
        use rayon::prelude::*;
        let group = count.div_ceil(pool.threads()).max(1);
        return pool.install(|| {
            (0..count)
                .into_par_iter()
                .step_by(group)
                .all(|start| allowed(start..(start + group).min(count)))
        });
    }
    let _ = pool;
    allowed(0..count)
}

// No NCBI counterpart: a copy of what LOSAT writes, for the shadow (LOSAT_X_BPQBLOCKFMTSHADOW); it
// does not change any byte written.
/// EXPERIMENT (LOSAT_X_BPQBLOCKFMTSHADOW): a writer that keeps a copy of the bytes it passes on
/// when `copy` is set (the shadow), and only passes them on otherwise. The copy is shared with the
/// shadow's observer (`XShadowObserver`), which reads its length as the position in the stream.
pub(crate) struct XCopyWriter<'a, W: Write + ?Sized> {
    pub(crate) inner: &'a mut W,
    pub(crate) copy: Option<&'a Mutex<Vec<u8>>>,
}

impl<W: Write + ?Sized> Write for XCopyWriter<'_, W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.inner.write(buf)?;
        if let Some(copy) = self.copy {
            copy.lock()
                .expect("the shadow's copy")
                .extend_from_slice(&buf[..written]);
        }
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
// ```c
// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
// ```
// Every format prints from the final result set of its call; a block's call numbers the HSPs of
// its own list from 0. No NCBI counterpart for the offset: it maps a block's numbers to those of
// the list of all queries, the blocks' lists one after the other.
/// EXPERIMENT (LOSAT_X_BPQBLOCKWEB): the caller's observer for one block, whose HSP indices are
/// those of the block's list plus `base`, the number of HSPs of the blocks before; with the shadow
/// (LOSAT_X_BPQBLOCKFMTSHADOW) each event is also recorded in `shadow`, with the index passed on.
pub(crate) struct XOffsetObserver<'a, 's> {
    pub(crate) inner: &'a mut (dyn FormatObserver + Send),
    pub(crate) base: HspIndex,
    pub(crate) shadow: Option<&'a mut XShadowObserver<'s>>,
}

impl FormatObserver for XOffsetObserver<'_, '_> {
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex) {
        if let Some(shadow) = self.shadow.as_deref_mut() {
            shadow.hsp_begin(format, self.base + hsp);
        }
        self.inner.hsp_begin(format, self.base + hsp);
    }

    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        if let Some(shadow) = self.shadow.as_deref_mut() {
            shadow.hsp_end(format, self.base + hsp);
        }
        self.inner.hsp_end(format, self.base + hsp);
    }

    fn subject_begin(&mut self, format: usize, first_hsp: HspIndex) {
        if let Some(shadow) = self.shadow.as_deref_mut() {
            shadow.subject_begin(format, self.base + first_hsp);
        }
        self.inner.subject_begin(format, self.base + first_hsp);
    }

    fn subject_end(&mut self, format: usize, first_hsp: HspIndex) {
        if let Some(shadow) = self.shadow.as_deref_mut() {
            shadow.subject_end(format, self.base + first_hsp);
        }
        self.inner.subject_end(format, self.base + first_hsp);
    }
}

// No NCBI counterpart: the record of the shadow (LOSAT_X_BPQBLOCKFMTSHADOW); it writes nothing.
/// EXPERIMENT (LOSAT_X_BPQBLOCKFMTSHADOW): an observer that records, per format, each event with its
/// HSP index and its position in the format's copy (`XCopyWriter`).
pub(crate) struct XShadowObserver<'a> {
    copies: &'a [Mutex<Vec<u8>>],
    events: Vec<Vec<(u8, HspIndex, usize)>>,
}

impl<'a> XShadowObserver<'a> {
    pub(crate) fn new(copies: &'a [Mutex<Vec<u8>>]) -> Self {
        Self {
            copies,
            events: vec![Vec::new(); copies.len()],
        }
    }

    fn record(&mut self, format: usize, event: u8, hsp: HspIndex) {
        let position = self.copies[format].lock().expect("the shadow's copy").len();
        self.events[format].push((event, hsp, position));
    }

    /// The number of events of all formats.
    pub(crate) fn count(&self) -> usize {
        self.events.iter().map(Vec::len).sum()
    }

    /// The events of each format as text, one line each: `begin`, `end`, `subject_begin` or
    /// `subject_end`, the HSP index and the position.
    pub(crate) fn text(&self) -> Vec<Vec<u8>> {
        const NAMES: [&str; 4] = ["begin", "end", "subject_begin", "subject_end"];
        self.events
            .iter()
            .map(|events| {
                events
                    .iter()
                    .map(|&(event, hsp, position)| {
                        format!("{} {hsp} {position}\n", NAMES[usize::from(event)])
                    })
                    .collect::<String>()
                    .into_bytes()
            })
            .collect()
    }
}

impl FormatObserver for XShadowObserver<'_> {
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex) {
        self.record(format, 0, hsp);
    }

    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        self.record(format, 1, hsp);
    }

    fn subject_begin(&mut self, format: usize, first_hsp: HspIndex) {
        self.record(format, 2, first_hsp);
    }

    fn subject_end(&mut self, format: usize, first_hsp: HspIndex) {
        self.record(format, 3, first_hsp);
    }
}

// No NCBI counterpart: the check of the shadow (LOSAT_X_BPQBLOCKFMTSHADOW); it does not change any
// value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCKFMTSHADOW): compares what the block search wrote with what the search
/// of all queries wrote into memory, stream by stream (`(name, blocks, reference)`: each output, the
/// warnings, and with LOSAT_X_BPQBLOCKWEB the hit records (`Debug`, one per line) and each format's
/// observer events (`XShadowObserver::text`)). On a difference it prints
/// `[X_SHADOW] BPQBLOCKFMT <stream> DIFFERENT ...` and panics; otherwise
/// `[X_SHADOW] BPQBLOCKFMT blocks= formats= bytes= stderr_bytes= hits= events= identical` (`bytes` of
/// all outputs, `hits` the records, `events` of all formats).
pub(crate) fn x_block_shadow_compare(
    blocks: usize,
    formats: usize,
    streams: &[(String, &[u8], &[u8])],
    hits: usize,
    events: usize,
) {
    let mut bytes = 0;
    let mut stderr_bytes = 0;
    for (stream, block_bytes, reference_bytes) in streams {
        if block_bytes != reference_bytes {
            let at = block_bytes
                .iter()
                .zip(reference_bytes.iter())
                .position(|(a, b)| a != b)
                .unwrap_or(block_bytes.len().min(reference_bytes.len()));
            eprintln!(
                "[X_SHADOW] BPQBLOCKFMT {stream} DIFFERENT at byte {at} (blocks {} bytes, search of all queries {} bytes)",
                block_bytes.len(),
                reference_bytes.len()
            );
            panic!("LOSAT_X_BPQBLOCKFMTSHADOW: the block search's {stream} differs from the search of all queries");
        }
        if stream.starts_with("output") {
            bytes += block_bytes.len();
        } else if stream == "stderr" {
            stderr_bytes += block_bytes.len();
        }
    }
    eprintln!(
        "[X_SHADOW] BPQBLOCKFMT blocks={blocks} formats={formats} bytes={bytes} stderr_bytes={stderr_bytes} hits={hits} events={events} identical"
    );
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
/// pool, `SearchPool::x_with_ext`), with the records they encode; and the vectors of the last
/// block's flat lookup table, for the next block's build (LOSAT_X_BPLKPFLAT, `x_bp_lookup_spare`),
/// and those of exact blastp's tables built ahead (LOSAT_X_BPLKPAHEAD, `x_bp_aa_spares`).
pub(crate) struct XBlockSubjects {
    /// The address and count of the records encoded.
    records: (usize, usize),
    /// EXPERIMENT (LOSAT_X_BPSUBJMETA): the records' tabular IDs and titles, made with the switch.
    x_meta: Option<super::blast_engine::x_subject_meta::XSubjectMeta>,
    encoded: Vec<EncodedProtein>,
    x_lookup_spare: XFlatSpare,
    /// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the second spare: the vectors of the next block's table,
    /// allocated or moved here by the block search's thread for the build of that table ahead of its
    /// block (`blast_engine/x_lookup_ahead.rs`, `x_bp_lookup_spares`).
    x_lookup_ahead: XFlatSpare,
    /// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the vectors and neighbour lists of exact blastp's tables
    /// built ahead of their blocks (`x_bp_aa_spares`).
    x_aa_spares: XAaSpares,
}

impl XBlockSubjects {
    // No NCBI counterpart: keeps the encoding of `records` (module comment above).
    pub(crate) fn new(records: &[FastaRecord], encoded: Vec<EncodedProtein>) -> Self {
        assert_eq!(records.len(), encoded.len());
        Self {
            records: (records.as_ptr() as usize, records.len()),
            x_meta: super::blast_engine::x_subject_meta::XSubjectMeta::x_for_block_search(records),
            encoded,
            x_lookup_spare: XFlatSpare::default(),
            x_lookup_ahead: XFlatSpare::default(),
            x_aa_spares: XAaSpares::default(),
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

// No NCBI counterpart: reads back the IDs and titles attached by the block search (reuse only).
/// EXPERIMENT (LOSAT_X_BPSUBJMETA): the subjects' tabular IDs and titles attached to `pool`
/// (`blast_engine/x_subject_meta.rs`), when the switch made them and they are those of `records`
/// (the same slice).
pub(crate) fn x_block_subject_meta<'p>(
    pool: &SearchPool<'p>,
    records: &[FastaRecord],
) -> Option<&'p super::blast_engine::x_subject_meta::XSubjectMeta> {
    pool.x_ext::<XBlockSubjects>()
        .filter(|subjects| subjects.records == (records.as_ptr() as usize, records.len()))
        .and_then(|subjects| subjects.x_meta.as_ref())
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

// No NCBI counterpart: where LOSAT_X_BPLKPFLAT acts; it reads switches and changes no value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the switch's mode (`x_bplkpflat_mode`) where `pool` builds
/// blastp-fast's flat table (`LOSAT_X_BPCLUT` or `LOSAT_X_BPCLUTSHADOW`, and `x_pool_builds`: a pool
/// of more than one thread, `lookup/x_clut_par.rs`), else `Off`.
fn x_lkpflat_mode_for_pool(pool: &SearchPool<'_>) -> XLkpFlat {
    if x_bpclut_mode() != 0 && x_pool_builds(pool) {
        x_bplkpflat_mode()
    } else {
        XLkpFlat::Off
    }
}

// No NCBI counterpart: buffer reuse between the blocks' lookup tables (`XFlatSpare`).
/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the block search's spare vectors of the flat tables
/// (`XBlockSubjects`), where the switch acts: the search of a block gives its table's vectors back
/// here, and the next block's build takes them.
pub(crate) fn x_bp_lookup_spare<'p>(pool: &SearchPool<'p>) -> Option<&'p XFlatSpare> {
    if x_lkpflat_mode_for_pool(pool) == XLkpFlat::Off {
        return None;
    }
    pool.x_ext::<XBlockSubjects>()
        .map(|subjects| &subjects.x_lookup_spare)
}

// No NCBI counterpart: buffer reuse between the blocks' lookup tables (`XFlatSpare`).
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the block search's two spares, where LOSAT_X_BPLKPFLAT acts: the
/// block spare (`x_bp_lookup_spare`) and the second spare, which the build of the next block's table
/// ahead of that block takes (`blast_engine/x_lookup_ahead.rs`).
pub(crate) fn x_bp_lookup_spares<'p>(
    pool: &SearchPool<'p>,
) -> Option<(&'p XFlatSpare, &'p XFlatSpare)> {
    x_bp_lookup_spare(pool)?;
    pool.x_ext::<XBlockSubjects>()
        .map(|subjects| (&subjects.x_lookup_spare, &subjects.x_lookup_ahead))
}

// No NCBI counterpart: where LOSAT_X_BPLKPAHEAD acts for exact blastp, and buffer reuse between its
// tables (`XAaSpares`); it reads switches and changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the block search's spares of exact blastp's tables
/// (`XAaSpares`), where the switch builds them ahead of their blocks
/// (`blast_engine/x_lookup_ahead.rs`): the switch is on, the table is built in place
/// (`XAaSpares::x_applies`), and the pool has more than one thread (`x_pool_builds`).
pub(crate) fn x_bp_aa_spares<'p>(pool: &SearchPool<'p>) -> Option<&'p XAaSpares> {
    use super::blast_engine::x_lookup_ahead::{x_bp_lkp_ahead, XLkpAhead};
    if x_bp_lkp_ahead() == XLkpAhead::Off || !XAaSpares::x_applies() || !x_pool_builds(pool) {
        return None;
    }
    pool.x_ext::<XBlockSubjects>()
        .map(|subjects| &subjects.x_aa_spares)
}

/// EXPERIMENT (LOSAT_X_BPLKPFLAT, LOSAT_X_BPLKPAHEAD): which build of a block search the lookup query
/// is for, which decides the spare vectors the pool's build of the flat table takes (`XLookupFlat`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XLookupFor {
    /// A block's (or a search's) own build: the block spare.
    Search,
    /// A query chunk of a split batch (`x_query_chunk`): none.
    QueryChunk,
    /// The next block's build ahead of it, on a pool task (LOSAT_X_BPLKPAHEAD): the second spare,
    /// and no vector of the table allocated by the build (`XFlatArgs::x_ahead`).
    Ahead,
}

/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the lookup query of blastp-fast's compressed table and what the
/// pool's build of its flat table gets with the switch (`lookup/x_clut_par.rs`): the block search's
/// spare vectors, and with `LOSAT_X_BPQBLOCK` the contexts' Karlin blocks as extra tasks of the
/// build's appends pass, in place of their own pass on the pool (`x_bp_lookup_query`).
pub(crate) struct XLookupFlat<'a> {
    mode: XLkpFlat,
    spare: Option<&'a XFlatSpare>,
    /// LOSAT_X_BPLKPAHEAD: the build ahead of its block (`XLookupFor::Ahead`).
    x_ahead: bool,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    karlin: Option<XKarlinTasks<'a>>,
}

impl<'a> XLookupFlat<'a> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803 (the contexts'
    // Karlin blocks, `x_bp_lookup_query`)
    // ```c
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // ```
    /// The values of `x_bp_lookup_query(query_frames, word_size, pool)`; where the switch acts with
    /// `LOSAT_X_BPQBLOCK`, the contexts come without their Karlin blocks, which `finish` sets after
    /// the build (`XKarlinTasks`). A query chunk of a split batch (`XLookupFor::QueryChunk`, a search
    /// that ends after its preliminary stage, inside its block's search, and drops its table there)
    /// takes no spare vectors: the block's table keeps them. The build of the next block's table
    /// ahead of it (`XLookupFor::Ahead`, LOSAT_X_BPLKPAHEAD) takes the second spare's.
    pub(crate) fn lookup_query<'p: 'a>(
        query_frames: &'a [Vec<QueryFrame>],
        word_size: usize,
        pool: &SearchPool<'p>,
        x_for: XLookupFor,
    ) -> (Self, Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>) {
        let mode = x_lkpflat_mode_for_pool(pool);
        if x_bplkpflat_mode() != XLkpFlat::Off {
            crate::utils::xreach::x_reach!("LOSAT_X_BPLKPFLAT", mode != XLkpFlat::Off);
        }
        let spare = match x_for {
            XLookupFor::Search => x_bp_lookup_spare(pool),
            XLookupFor::QueryChunk => None,
            XLookupFor::Ahead => x_bp_lookup_spares(pool).map(|(_, ahead)| ahead),
        };
        let x_ahead = x_for == XLookupFor::Ahead;
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if mode != XLkpFlat::Off && x_bp_qblock().is_some() {
            let (mut flat, query, locations, contexts) =
                Self::with_karlin_tasks(query_frames, word_size, mode, spare);
            flat.x_ahead = x_ahead;
            return (flat, query, locations, contexts);
        }
        let (query, locations, contexts) = x_bp_lookup_query(query_frames, word_size, pool);
        let flat = Self {
            mode,
            spare,
            x_ahead,
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            karlin: None,
        };
        (flat, query, locations, contexts)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2781
    // ```c
    //       if ( !contexts[context].is_valid )
    //           continue;
    //
    //       query_length = contexts[context].query_length;
    //       context_offset = contexts[context].query_offset;
    //       buffer = &query[context_offset];
    //
    //       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
    //       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    //       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
    //       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
    // ```
    // NCBI makes the lookup segments (c++/src/algo/blast/core/blast_setup.c:635) before these
    // blocks (BlastSetup_ScoreBlkInit, blast_setup.c:652); the blocks are left to the build's tasks.
    /// The lookup query without the contexts' Karlin blocks, which become the build's extra tasks.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn with_karlin_tasks(
        query_frames: &'a [Vec<QueryFrame>],
        word_size: usize,
        mode: XLkpFlat,
        spare: Option<&'a XFlatSpare>,
    ) -> (Self, Vec<u8>, Vec<(i32, i32)>, Vec<QueryContext>) {
        let lookup_query = XBlosum62LookupQuery::new(false);
        let (query, locations, contexts) = lookup_query.without_karlin(query_frames, word_size);
        let frames: Vec<&QueryFrame> = query_frames.iter().flatten().collect();
        let karlin = XKarlinTasks {
            query_frames,
            word_size,
            done: frames.iter().map(|_| OnceLock::new()).collect(),
            frames,
            lookup_query,
        };
        let flat = Self {
            mode,
            spare,
            x_ahead: false,
            karlin: Some(karlin),
        };
        (flat, query, locations, contexts)
    }

    // No NCBI counterpart: the arguments of the build (`x_bp_build_compressed_lookup`).
    /// The switch's mode, the spare vectors and the extra tasks, for the build.
    pub(crate) fn build_args(&self) -> XFlatArgs<'_> {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let extra = self
            .karlin
            .as_ref()
            .map(|karlin| karlin as &dyn XExtraTasks);
        #[cfg(any(
            not(feature = "parallel"),
            all(target_arch = "wasm32", not(feature = "wasm-threads"))
        ))]
        let extra = None;
        XFlatArgs {
            mode: self.mode,
            spare: self.spare,
            extra,
            x_ahead: self.x_ahead,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803
    // ```c
    //       if ( !contexts[context].is_valid )
    //           continue;
    //       ...
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // ```
    /// After the build: with Karlin tasks (`lookup_query`), sets the contexts' validity and Karlin
    /// blocks that the build's tasks computed and computes the others here; the tasks do not run
    /// where the pool's build falls back to the serial builder (a backbone too large for its
    /// stripes, `x_try_build_flat_parallel_with`); a build that fails before its passes ends the search before
    /// this call. Without Karlin tasks it does nothing. With the shadow mode, compares the query,
    /// the segments and the contexts with `x_bp_lookup_query`'s, the Karlin blocks bit for bit
    /// (panics on a difference). Returns the number of contexts the tasks computed (for the tests;
    /// the shadow counts them as well).
    pub(crate) fn finish(
        self,
        contexts: &mut [QueryContext],
        query: &[u8],
        locations: &[(i32, i32)],
        pool: &SearchPool<'_>,
    ) -> usize {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if let Some(karlin) = self.karlin {
            let from_tasks = karlin.fill(contexts);
            if self.mode == XLkpFlat::Shadow {
                let (reference_query, reference_locations, reference_contexts) =
                    x_bp_lookup_query(karlin.query_frames, karlin.word_size, pool);
                assert!(
                    reference_query == query && reference_locations == locations,
                    "LOSAT_X_BPLKPFLATSHADOW: lookup query differs"
                );
                x_assert_same_contexts_as(
                    "LOSAT_X_BPLKPFLATSHADOW:",
                    contexts,
                    &reference_contexts,
                );
                x_lkpflat_count_contexts(contexts.len(), from_tasks);
            }
            return from_tasks;
        }
        let _ = (contexts, query, locations, pool);
        0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803
// ```c
//    for (context = query_info->first_context;
//         context <= query_info->last_context; ++context) {
//       ...
//       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
// ```
// NCBI computes the contexts' Karlin blocks one after the other; each reads only its own context's
// letters, so computing them as tasks of the lookup build (on any thread, in any order) gives each
// context the same block (`XBlosum62LookupQuery`). The table's build reads only the query, the
// segments and the contexts' `q_idx` and `frame_base` (`ncbi_compressed_lookup_tables`), never a
// Karlin block, and the computation has no error path (an invalid context is a flag), so no error
// moves. Scheduling only.
/// EXPERIMENT (LOSAT_X_BPLKPFLAT): the contexts' Karlin blocks as extra tasks of the build
/// (`XExtraTasks`), one per context.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XKarlinTasks<'a> {
    query_frames: &'a [Vec<QueryFrame>],
    word_size: usize,
    /// The frames in the order of the contexts.
    frames: Vec<&'a QueryFrame>,
    lookup_query: XBlosum62LookupQuery,
    /// Each context's validity and Karlin block, once its task has run.
    done: Vec<OnceLock<(bool, KarlinParams)>>,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl XExtraTasks for XKarlinTasks<'_> {
    fn count(&self) -> usize {
        self.frames.len()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2778-2782
    // ```c
    //       Blast_ResFreqString(sbp, rfp, (char*)buffer, query_length);
    //       sbp->sfp[context] = Blast_ScoreFreqNew(sbp->loscore, sbp->hiscore);
    //       BlastScoreFreqCalc(sbp, sbp->sfp[context], rfp, stdrfp);
    //       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    // ```
    fn run(&self, index: usize) {
        let karlin = self.lookup_query.context_karlin(self.frames[index]);
        let stored = self.done[index].set(karlin).is_ok();
        assert!(stored, "the Karlin block of context {index} computed twice");
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl XKarlinTasks<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2782-2784
    // ```c
    //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
    //       if (loop_status) {
    //           contexts[context].is_valid = FALSE;
    // ```
    /// Sets each context's validity and Karlin block from its task, or computes it here when its
    /// task has not run. Returns the number of contexts set from the tasks.
    fn fill(&self, contexts: &mut [QueryContext]) -> usize {
        assert_eq!(contexts.len(), self.frames.len());
        let mut from_tasks = 0usize;
        for ((context, frame), done) in contexts.iter_mut().zip(&self.frames).zip(&self.done) {
            let (is_valid, karlin_params) = match done.get() {
                Some(&computed) => {
                    from_tasks += 1;
                    computed
                }
                None => self.lookup_query.context_karlin(frame),
            };
            context.is_valid = is_valid;
            context.karlin_params = karlin_params;
        }
        from_tasks
    }
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
    // The same with the blocks computed as extra tasks of the flat table's build
    // (LOSAT_X_BPLKPFLAT, `XLookupFlat`).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn lookup_query_with_karlin_blocks_on_a_pool_is_the_references() {
        use crate::algorithm::blastp::encoding::encode_protein_query_frame_with_seg;
        use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
            x_build_flat_in_pool, XFlatBuild,
        };
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
                // LOSAT_X_BPLKPFLAT: the Karlin blocks as extra tasks of the pool's build of the
                // flat table (every one of them at 3 and 7 threads; at one thread the serial builder
                // runs none and `finish` computes them all), with the shadow's checks too.
                for mode in [XLkpFlat::On, XLkpFlat::Shadow] {
                    let from_tasks =
                        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                            let (flat, flat_query, flat_locations, mut flat_contexts) =
                                XLookupFlat::with_karlin_tasks(&frames, 5, mode, None);
                            assert_eq!(flat_query, query);
                            assert_eq!(flat_locations, locations);
                            let built = x_build_flat_in_pool(
                                5,
                                1100,
                                &flat_query,
                                &flat_locations,
                                &[0..flat_locations.len()],
                                pool,
                                flat.build_args(),
                            );
                            assert!(matches!(built, Some(XFlatBuild::Table(_))));
                            let from_tasks =
                                flat.finish(&mut flat_contexts, &flat_query, &flat_locations, pool);
                            x_assert_same_contexts(&flat_contexts, &contexts);
                            Ok(from_tasks)
                        })
                        .expect("pool");
                    let expected = if threads > 1 { contexts.len() } else { 0 };
                    assert_eq!(from_tasks, expected, "threads {threads}, {mode:?}");
                }
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

    // The XML checks of the block search agree with the writer's own functions on the probe titles
    // of `report/blast_xml/names.rs`: names that pass (periods, `TPA:`, doubled spaces, an empty
    // title), a byte >= 0x80 (UTF-8 and Latin-1), a control byte, and punctuation that
    // `x_CleanAndCompress` reads past the end of. A one-record check and the same record among
    // others, on a pool of 1 and 3 threads.
    #[test]
    fn the_xml_checks_agree_with_the_writers_functions() {
        use crate::report::blast_xml::{names, XmlSubjects};
        let titles: [&[u8]; 12] = [
            b"q1 title ending with periods...",
            b"TPA: q3 prefixed",
            b"q4  double  space ,  comma ;",
            b"",
            b"s2 &amp; <y> \"q\"",
            b"tab\there",
            "caf\u{e9}".as_bytes(),
            b"caf\xe9",
            b"a\x01b",
            b", ,",
            b",",
            b"x \xc3\xa9 y",
        ];
        let record = |index: usize, title: &[u8]| {
            FastaRecord::new(
                format!("Query_{}", index + 1),
                title,
                b"MKVLAAGIVALLLAAGCSS",
            )
        };
        let ok_query = |record: &FastaRecord| {
            names::query_texts(
                "BLASTP",
                std::slice::from_ref(record),
                true,
                &[Vec::new()],
                |_| false,
            )
            .is_ok()
        };
        let ok_subject = |record: &FastaRecord| {
            XmlSubjects::Records {
                records: std::slice::from_ref(record),
                length: &|_| 0,
                protein: true,
            }
            .names(0, "BLASTP")
            .is_ok()
        };
        let records: Vec<FastaRecord> = titles
            .iter()
            .enumerate()
            .map(|(index, title)| record(index, title))
            .collect();
        let passing = records.iter().filter(|record| ok_query(record)).count();
        assert!(
            passing >= 5 && passing < records.len(),
            "the probe titles pass and fail: {passing}"
        );
        for threads in [1, 3] {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                for record in &records {
                    let one = std::slice::from_ref(record);
                    assert_eq!(
                        x_block_xml_run_allowed("db", one, pool),
                        ok_query(record),
                        "query title {:?}, threads {threads}",
                        String::from_utf8_lossy(&record.title)
                    );
                    assert_eq!(
                        x_block_xml_subjects_allowed(one, pool),
                        ok_subject(record),
                        "subject title {:?}, threads {threads}",
                        String::from_utf8_lossy(&record.title)
                    );
                }
                // One refusable record among many fails the whole set, wherever it lies.
                let many = |bad: Option<usize>| -> Vec<FastaRecord> {
                    (0..50)
                        .map(|index| {
                            let title: &[u8] = if Some(index) == bad {
                                "caf\u{e9}".as_bytes()
                            } else {
                                b"ok title"
                            };
                            record(index, title)
                        })
                        .collect()
                };
                assert!(x_block_xml_run_allowed("db", &many(None), pool));
                assert!(x_block_xml_subjects_allowed(&many(None), pool));
                for bad in [0, 17, 49] {
                    assert!(!x_block_xml_run_allowed("db", &many(Some(bad)), pool));
                    assert!(!x_block_xml_subjects_allowed(&many(Some(bad)), pool));
                }
                // The database name is checked as the writer checks it.
                assert!(x_block_xml_run_allowed("", &many(None), pool));
                assert!(!x_block_xml_run_allowed("d\u{e9}", &many(None), pool));
                assert!(!x_block_xml_run_allowed("d\x01", &many(None), pool));
                Ok(())
            })
            .expect("pool");
        }
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
