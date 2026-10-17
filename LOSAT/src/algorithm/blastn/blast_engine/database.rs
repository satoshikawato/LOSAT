//! The database path of BLASTN (`-db NAME`; tasks megablast, blastn, dc-megablast and
//! blastn-short): the search of an NCBI BLAST nucleotide database read in blocks through
//! `crate::seqdb`, as NCBI BLAST+ 2.17.0 `blastn -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, stage DB-4).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the set-level numbers come from the index files (`BlastnDatabase::metadata`): the
//!   number of sequences and the total length (an `Int8`, not the `int` of the subject
//!   set's report);
//! - the preliminary stage reads the stored 2-bit data in blocks and searches them in the
//!   order of their ordinal ids (`search_blocks`); a sequence's preliminary HSPs go into
//!   the preliminary hit lists as soon as its block's round is searched, so nothing but the
//!   hit lists is kept for the whole database;
//! - the traceback reads again, by ordinal id, the sequences of the preliminary hit lists
//!   (`traceback_listed`), with their ambiguity codes (BLASTNA);
//! - the definition lines of the final hit lists are read after each query batch's
//!   traceback (`read_titles`), and the reports show the database
//!   (`BlastnDatabase::pairwise_database`, `name`) and the stored titles and letters of the
//!   reported sequences only (`reported_subjects`).
//!
//! The per-subject kernel, the hit lists and the report writers are those of the
//! `-subject` path (`run.rs`); the opening of the database, the reader's errors and the
//! rules for the stored titles are shared with BLASTP (`algorithm::common::database`).

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::{Context, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::algorithm::blastn::coordination::SubjectMetadata;
use crate::algorithm::blastn::hsp::BlastnHitList;
use crate::algorithm::common::database::{self as shared, StoredTitle};
pub(super) use crate::algorithm::common::database::{HeaderReadError, PreliminaryReadError};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{BlockRange, SeqDb, SeqDbError, SeqType, SequenceBlock};
use crate::utils::threading::SearchPool;

use super::BlastnSubjectRef;

/// The program that LOSAT's refusals name.
const PROGRAM: &str = "BLASTN";

/// Bytes of the sequence file that one block of the preliminary stage holds at most
/// (`SeqDb::next_block`; a longer sequence is a block of its own). No NCBI counterpart:
/// NCBI maps the volume files and hands its threads chunks of ordinal ids
/// (seqdbimpl.cpp:328-350). The size does not change which sequences are searched or
/// their order.
const BLOCK_BYTES: u64 = 1 << 20;

/// Blocks that one parallel round of the preliminary stage reads and searches per
/// thread before their results are merged in ordinal order. No NCBI counterpart (as
/// `BLOCK_BYTES`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const BLOCKS_PER_THREAD: usize = 8;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:75-78
// ```c
// const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
//     'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K',
//     'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
// };
// ```
/// The IUPAC letter of a BLASTNA code.
const BLASTNA_TO_IUPACNA: [u8; 16] = *b"ACGTRYMKWSBDHVN-";

/// A reader error outside the preliminary search as the command line reports it
/// (`algorithm::common::database::database_error`).
fn database_error(error: SeqDbError) -> anyhow::Error {
    shared::database_error(error, PROGRAM)
}

/// A reader error while the preliminary search reads the sequences
/// (`algorithm::common::database::search_read_error`).
fn search_read_error(error: SeqDbError) -> anyhow::Error {
    shared::search_read_error(error, PROGRAM)
}

/// The database of a `-db` search: the opened database, its name as NCBI keeps it
/// (`BlastnArgs::database_name`: the `-db` value without white space at its ends) and its
/// numbers as the search uses them.
pub(super) struct BlastnDatabase {
    db: SeqDb,
    name: String,
    metadata: SubjectMetadata,
}

impl BlastnDatabase {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    //         _ASSERT(db_adapter && scope);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:182-187
    // ```c++
    // static Int8
    // s_SeqDbGetTotLen(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetTotalLength();
    // }
    // ```
    /// Opens the database of `-db name` (`algorithm::common::database::open_database`)
    /// with the command line's errors. The total length is NCBI's `Int8`; one that the
    /// engine cannot hold in a `usize` (a 32-bit target) is LOSAT's refusal.
    pub(super) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        let db = shared::open_database(name, SeqType::Nucleotide, settings, PROGRAM)?;
        let total = db.total_length();
        let db_len_total = usize::try_from(total)
            .ok()
            .filter(|&letters| i64::try_from(letters).is_ok())
            .with_context(|| {
                format!(
                    "a database of {total} letters is not supported by LOSAT's BLASTN on this platform"
                )
            })?;
        let metadata = SubjectMetadata {
            db_len_total,
            db_num_seqs: db.num_seqs() as usize,
        };
        Ok(Self {
            db,
            name: name.to_string(),
            metadata,
        })
    }

    /// The `-db` value as NCBI keeps it (outfmt 7's `# Database:` line).
    pub(super) fn name(&self) -> &str {
        &self.name
    }

    /// The number of sequences (`CSeqDB::GetNumSeqs`) and the total length
    /// (`CSeqDB::GetTotalLength`): the real numbers, also with `-dbsize`.
    pub(super) fn metadata(&self) -> &SubjectMetadata {
        &self.metadata
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`).
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db, &self.name)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1405
// ```c
//     memset((void*) &seq_arg, 0, sizeof(seq_arg));
//
//     /* Encoding is set so there are no sentinel bytes, and protein/nucleotide
//        sequences are retieved in ncbistdaa/ncbi2na encodings respectively. */
//     seq_arg.encoding = eBlastEncodingProtein;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1805-1838
// ```c++
//     } else if ('n' == seqtype) {
//         // The last byte is partially full; the last two bits of the
//         // last byte store the number of nucleotides in the last byte
//         // (0 to 3).
//     ...
//         *buffer = m_Seq->GetFileDataPtr(start_offset);
//     ...
//         int whole_bytes = int(end_offset - start_offset - 1);
//
//         char last_char = (*buffer)[whole_bytes];
//
//         int remainder = last_char & 3;
//         length = (whole_bytes * 4) + remainder;
//     }
// ```
// The preliminary stage searches the stored 2-bit data as makeblastdb wrote it (its
// ambiguous bases already replaced by makeblastdb's random bases), in the layout of the
// subject set's packed letters (four bases per byte from the high bits down, the count of
// the last byte's bases in its two low bits; `encode_subject_ncbi2na_packed`).
/// Sequence `i` of a block as the preliminary stage reads it: its length and its stored
/// 2-bit data.
fn block_subject(block: &SequenceBlock, i: usize) -> Result<BlastnSubjectRef<'_>> {
    let length = block.nucleotide_len(i).map_err(search_read_error)?;
    Ok(BlastnSubjectRef {
        length,
        packed: block.nucleotide_packed(i),
        letters: None,
        blastna: None,
        store: None,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
// ```c
//     db_length = BlastSeqSrcGetTotLen(seq_src);
//
//     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
//        Int4 stat_length;
//        if (seq_arg.oid == BLAST_SEQSRC_ERROR)
//            break;
//        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0)
//            continue;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
// ```c
//           /* Save the results. */
//           status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
/// The preliminary stage of a database search: every sequence, in the order of its
/// ordinal id, through `search` (the per-subject kernel), and what it keeps of a sequence
/// (`Some`) through `merge` in that order (NCBI's loop with one thread; `merge` puts the
/// sequence's preliminary HSPs into the preliminary hit lists).
///
/// Without `parallel`: the blocks in order, one scratch (`new_scratch`) for the whole
/// database, as the one-thread loop of the `-subject` path.
///
/// With `parallel`: rounds of blocks. The blocks of a round are read and searched in
/// parallel with a scratch per worker (`map_init`), as the parallel `-subject` path
/// searches its subjects (the result of a subject does not depend on what the worker's
/// scratch searched before it, which that path relies on). The round's results are merged
/// in ordinal order and its blocks dropped before the next round is read: a worker holds
/// one block at a time, so the sequence memory is bounded by the threads and the block
/// size, not by the database. A block that cannot be read does not stop the other blocks
/// of its round; the round's first error in ordinal order is returned, the error of the
/// one-thread order (the blocks are cut alike).
pub(super) fn search_blocks<T, R, I, K, M>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R),
{
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let blocks_per_thread = BLOCKS_PER_THREAD;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let blocks_per_thread = 1;
    search_blocks_of(
        database,
        pool,
        parallel,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        search,
        merge,
    )
}

/// `search_blocks` with blocks of at most `block_bytes` and rounds of `blocks_per_thread`
/// blocks per thread.
#[allow(clippy::too_many_arguments)]
fn search_blocks_of<T, R, I, K, M>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &I,
    search: &K,
    merge: &mut M,
) -> Result<()>
where
    R: Send,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>) -> Option<R> + Sync + Send,
    M: FnMut(u32, R),
{
    let db = &database.db;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parallel {
        let round_blocks = pool.threads().saturating_mul(blocks_per_thread).max(1);
        let mut next_ordinal = 0u32;
        loop {
            let mut round: Vec<BlockRange> = Vec::with_capacity(round_blocks);
            while round.len() < round_blocks {
                let Some(range) = db.next_block(next_ordinal, block_bytes) else {
                    break;
                };
                round.push(range);
                next_ordinal = range.end_ordinal();
            }
            if round.is_empty() {
                break;
            }
            let found: Vec<Result<Vec<(u32, R)>>> = pool.install(|| {
                round
                    .par_iter()
                    .map_init(new_scratch, |scratch, &range| -> Result<Vec<(u32, R)>> {
                        let block = db.read_block(range).map_err(search_read_error)?;
                        let mut kept = Vec::new();
                        for i in 0..block.len() {
                            let ordinal = range.first_ordinal + i as u32;
                            if let Some(result) =
                                search(scratch, ordinal, block_subject(&block, i)?)
                            {
                                kept.push((ordinal, result));
                            }
                        }
                        Ok(kept)
                    })
                    .collect()
            });
            // The sequences that keep nothing change no hit list, so the order of the
            // merges is the one-thread order.
            for block_results in found {
                for (ordinal, result) in block_results? {
                    merge(ordinal, result);
                }
            }
        }
        return Ok(());
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, parallel, blocks_per_thread);

    let mut scratch = new_scratch();
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block = db.read_block(range).map_err(search_read_error)?;
        for i in 0..block.len() {
            let ordinal = range.first_ordinal + i as u32;
            if let Some(result) = search(&mut scratch, ordinal, block_subject(&block, i)?) {
                merge(ordinal, result);
            }
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:196-203
// ```c
//    /* sort in order of decreasing subject OID. HSPLists will be
//       read out from the end of hsplist_array later */
//
//    hsp_stream->num_hsplists = num_hsplists;
//    if (num_hsplists > 1) {
//       qsort(hsp_stream->sorted_hsplists, num_hsplists,
//                     sizeof(BlastHSPList *), s_SortHSPListByOid);
//    }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:819-828
// ```c
// EBlastEncoding Blast_TracebackGetEncoding(EBlastProgramType program_number)
// {
//     if (Blast_SubjectIsProtein(program_number)) {
//         return eBlastEncodingProtein;
//     } else if (Blast_SubjectIsTranslated(program_number)) {
//         return eBlastEncodingNcbi4na;
//     } else {
//         return eBlastEncodingNucleotide;
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1567-1573
// ```c
//                 seq_arg.oid = batch->hsplist_array[0]->oid;
//                 seq_arg.encoding = encoding;
//                 seq_arg.check_oid_exclusion = TRUE;
//                 seq_arg.reset_ranges = FALSE;
//                 seq_arg.ranges = ranges;
//
//                 if (BlastSeqSrcGetSequence(seqsrc, &seq_arg) < 0) {
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
// ```c++
//     m_TbackSearch.Reset(new CBlastTracebackSearch(m_QueryFactory,
//                                                   m_InternalData,
//                                                   m_Opts,
//                                                   seqinfo_src,
//                                                   search_msgs));
//     ...
//     m_TbackSearch->SetNumberOfThreads(GetNumberOfThreads());
//     CRef<CSearchResultSet> retval = m_TbackSearch->Run();
// ```
// The sequence comes in BLASTNA with its ambiguity codes, as the whole sequence (NCBI's
// partial fetch of the HSP ranges with fences gives the same alignments; decision 87 of
// the campaign: the whole sequence first).
// The traceback is outside the handler of the reader's errors of the preliminary search
// (local_blast.cpp:240-251), so a sequence that cannot be read there is the application's
// `BLAST Database error` (`database_error`).
/// The traceback of a database search: the sequences of `listed` (their ordinal ids,
/// ascending, with what the preliminary stage kept of each) read again one at a time in
/// the worker that searches them (the stored 2-bit data and the BLASTNA letters with the
/// ambiguity codes, `SequenceBlock::nucleotide_blastna`, dropped after the sequence), and
/// what `traceback` returns of each, in the order of `listed`. With `parallel` the
/// sequences are searched in parallel with a scratch per worker (`map_init`), as the
/// parallel `-subject` path searches its subjects.
pub(super) fn traceback_listed<T, P, R, I, K>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    listed: &[(u32, P)],
    new_scratch: &I,
    traceback: &K,
) -> Result<Vec<R>>
where
    P: Sync,
    R: Send,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>, &P) -> Option<R> + Sync + Send,
{
    let db = &database.db;
    let one = |scratch: &mut T, &(ordinal, ref kept): &(u32, P)| -> Result<Option<R>> {
        let stored = db.read_sequence(ordinal).map_err(database_error)?;
        let length = stored.nucleotide_len(0).map_err(database_error)?;
        let blastna = stored.nucleotide_blastna(0).map_err(database_error)?;
        let subject = BlastnSubjectRef {
            length,
            packed: stored.nucleotide_packed(0),
            letters: None,
            // The kernel's BLASTNA has no sentinel bytes.
            blastna: Some(&blastna[1..=length]),
            store: None,
        };
        Ok(traceback(scratch, ordinal, subject, kept))
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parallel {
        let found: Vec<Result<Option<R>>> = pool.install(|| {
            listed
                .par_iter()
                .map_init(new_scratch, |scratch, entry| one(scratch, entry))
                .collect()
        });
        let mut results = Vec::with_capacity(found.len());
        for result in found {
            results.extend(result?);
        }
        return Ok(results);
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, parallel);

    let mut scratch = new_scratch();
    let mut results = Vec::new();
    for entry in listed {
        results.extend(one(&mut scratch, entry)?);
    }
    Ok(results)
}

/// The stored titles of the database sequences that the final hit lists of the query
/// batches searched so far name, by ordinal id (`read_titles`).
pub(super) type DatabaseTitles = BTreeMap<u32, StoredTitle>;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
// ```c++
//     CRef<ILocalQueryData> qdata = m_QueryFactory->MakeLocalQueryData(m_Options);
//
//     vector<TSeqLocInfoVector> subj_masks;
//     TSeqAlignVector aligns =
//         LocalBlastResults2SeqAlign(hsp_results,
//                                    *qdata,
//                                    *m_SeqInfoSrc,
//                                    m_OptsMemento->m_ProgramType,
//                                    m_Options->GetGappedMode(),
//                                    m_Options->GetOutOfFrameMode(),
//                                    subj_masks,
//                                    m_ResultType);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1582
// ```c++
//     for (int index = 0; index < hit_list->hsplist_count; index++) {
//         BlastHSPList* hsp_list = hit_list->hsplist_array[index];
//         if (!hsp_list)
//             continue;
//
//         // Sort HSPs with e-values as first priority and scores as
//         // tie-breakers, since that is the order we want to see them in
//         // in Seq-aligns.
//         Blast_HSPListSortByEvalue(hsp_list);
//
//         const Uint4 kOid = hsp_list->oid;
//         TSeqPos subj_length = 0;
//         CRef<CSeq_id> subject_id;
//         GetSequenceLengthAndId(seqinfo_src, kOid, subject_id, &subj_length);
// ```
// NCBI first reads a header file (and opens it, seqdbvol.cpp:127-133) when the traceback's
// results become Seq-aligns: after the batch's traceback, for the final hit lists only
// (after their sort by e-value and their cut to the hit list size, before the HSPs of
// score 0 are dropped), query by query and in the order of each hit list. A header file
// that is missing or cannot be read fails the batch there (`BLAST Database error`, exit
// status 2, blast_app_util.hpp:185-203, outside the search stage's handler of
// local_blast.cpp:240-251), after the outfmt 0 prolog and the reports of the batches
// before; a batch without hits reads no header. A title read in an earlier batch is not
// read again: the header file stays open, and the title is the same.
/// Reads the definition lines of the sequences that the final hit lists `hit_lists` of a
/// query batch name, in NCBI's order, into `titles`. A reader error of NCBI's own is a
/// `HeaderReadError` (it fails the batch); LOSAT's refusals are errors of their own.
pub(super) fn read_titles(
    database: &BlastnDatabase,
    hit_lists: &[Option<BlastnHitList>],
    titles: &mut DatabaseTitles,
) -> Result<()> {
    for (query, hit_list) in hit_lists.iter().enumerate() {
        let Some(hit_list) = hit_list else {
            continue;
        };
        let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
        for hsp_list in &hit_list.hsplist_array[..count] {
            if titles.contains_key(&hsp_list.oid) {
                continue;
            }
            let title = shared::read_title(&database.db, hsp_list.oid, query, PROGRAM)?;
            titles.insert(hsp_list.oid, title);
        }
    }
    Ok(())
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(super) struct ReportedSubjects {
    /// The ordinal id of each rank, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The record of each rank: the id `gnl|BL_ORD_ID|n`, the stored title, and the
    /// stored letters (empty when the reports do not show them).
    pub(super) records: Vec<FastaRecord>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1406-1414
// ```c++
//             x_GetAmbChar(oid, ambchars);
//
//             if (ambchars.empty()) {
//                 // keep as 2 bit
//                 s_SeqDBWriteSeqDataNucl(seqinst, seq_buffer, length);
//             } else {
//                 // translate to 4 bit
//                 s_SeqDBWriteSeqDataNucl(seqinst, seq_buffer, length, ambchars);
//             }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:80-83
// ```c
// const char NCBI4NA_TO_IUPACNA[BLASTNA_SIZE] = {
//     '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
//     'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
// };
// ```
// The formatter reads a database sequence from its Bioseq: the stored 2-bit data, or with
// ambiguity codes their NCBI4na values (`x_GetAmbigSeq`'s letters before their mapping to
// BLASTNA; `SequenceBlock::nucleotide_blastna` maps them, and `BLASTNA_TO_IUPACNA` gives
// the same letter for a BLASTNA code as `NCBI4NA_TO_IUPACNA` for its NCBI4na value).
// Upper case: a database has no lowercase letters.
/// The reported database sequences of `hit_lists` (the final hit lists of all query
/// batches, after their HSPs of score 0 were dropped): the hit lists name them by rank
/// from here on (the rank grows with the ordinal id, so every comparison that breaks a tie
/// by the ordinal id gives the same order), and their records are made for the reports
/// from `titles`, with their stored letters when `letters` (outfmt 0 and the hit records
/// read them; the tabular formats do not). A reported sequence whose stored title LOSAT
/// does not report as NCBI does is LOSAT's refusal.
pub(super) fn reported_subjects(
    database: &BlastnDatabase,
    titles: &DatabaseTitles,
    hit_lists: &mut [Option<BlastnHitList>],
    letters: bool,
) -> Result<ReportedSubjects> {
    let mut ordinals: Vec<u32> = hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| hit_list.hsplist_array.iter())
        .filter(|hsp_list| !hsp_list.hsps.is_empty())
        .map(|hsp_list| hsp_list.oid)
        .collect();
    ordinals.sort_unstable();
    ordinals.dedup();
    let rank = |ordinal: u32| -> u32 {
        ordinals
            .binary_search(&ordinal)
            .expect("every reported subject was listed") as u32
    };
    for hit_list in hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            if hsp_list.hsps.is_empty() {
                continue;
            }
            hsp_list.oid = rank(hsp_list.oid);
            for hsp in &mut hsp_list.hsps {
                hsp.s_idx = rank(hsp.s_idx);
            }
        }
    }
    let mut records = Vec::with_capacity(ordinals.len());
    for &ordinal in &ordinals {
        let stored = titles
            .get(&ordinal)
            .with_context(|| format!("database sequence {ordinal} has no title read"))?;
        shared::check_reported_title(stored.replaced_bytes, ordinal, PROGRAM)?;
        let sequence = if letters {
            let block = database.db.read_sequence(ordinal).map_err(database_error)?;
            let length = block.nucleotide_len(0).map_err(database_error)?;
            let blastna = block.nucleotide_blastna(0).map_err(database_error)?;
            blastna[1..=length]
                .iter()
                .map(|&code| BLASTNA_TO_IUPACNA[usize::from(code & 0xF)])
                .collect()
        } else {
            Vec::new()
        };
        records.push(FastaRecord {
            local_id: format!("gnl|BL_ORD_ID|{ordinal}"),
            title: stored.title.clone(),
            sequence,
            warnings: Vec::new(),
        });
    }
    Ok(ReportedSubjects { ordinals, records })
}

impl ReportedSubjects {
    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every
    /// reported sequence, by rank (`algorithm::common::database::tabular_subject_id`). A
    /// sequence whose id LOSAT cannot write as NCBI is an error before any format is
    /// written.
    pub(super) fn tabular_subject_ids(&self) -> Result<Vec<Arc<[u8]>>> {
        self.records
            .iter()
            .zip(&self.ordinals)
            .map(|(record, &ordinal)| {
                shared::tabular_subject_id(&record.title, ordinal, false, PROGRAM).map(Arc::from)
            })
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:310-312
// ```c++
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
/// EXPERIMENT (LOSAT_X_RPTSTREAM): whether every definition line of the database reads as the
/// reports read those of the reported sequences (`read_title`: NCBI's reader, one definition
/// line), holds no byte outside ' ' ..= '~' (`check_reported_title`) and gives a tabular subject
/// id (`ReportedSubjects::tabular_subject_ids`). Then no batch after the first can meet one of
/// LOSAT's refusals of a title, which come before any of the search's reports without the
/// switch, so each batch's reports can be written when the loop takes it (`x_rpt_stream`; the
/// test of blastp's block search, `blastp::blast_engine::database::x_titles_allow_blocks`). Read
/// in `pool` over the ordinal ids; a reader error is `false` (the reports then stay after the
/// loop, and the batch that reads the title fails where NCBI's does).
pub(super) fn x_titles_allow_streaming(database: &BlastnDatabase, pool: &SearchPool<'_>) -> bool {
    let db = &database.db;
    let allowed = |ordinal: u32| -> bool {
        match shared::read_title(db, ordinal, 0, PROGRAM) {
            Ok(stored) => {
                !stored.replaced_bytes
                    && shared::tabular_subject_id(&stored.title, ordinal, false, PROGRAM).is_ok()
            }
            Err(_) => false,
        }
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() {
        return pool.install(|| {
            (0..db.num_seqs())
                .into_par_iter()
                .with_min_len(4096)
                .all(allowed)
        });
    }
    let _ = pool;
    (0..db.num_seqs()).all(allowed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastn::hsp::{BlastnHsp, BlastnHspList};
    use crate::cli::NativeError;

    fn fixture(name: &str) -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/db")
            .join(name)
            .to_str()
            .expect("fixture path is UTF-8")
            .to_string()
    }

    fn open(name: &str) -> BlastnDatabase {
        BlastnDatabase::open(&fixture(name), &ApplicationSettings::default())
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    /// LOSAT/tests/seqdb/expected/nucl_mv.entries.tsv (blastdbcmd): the lengths, letters and
    /// titles in ordinal order.
    fn entries() -> Vec<(usize, String, String)> {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/expected/nucl_mv.entries.tsv");
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                (
                    fields[1].parse().unwrap(),
                    fields[3].to_string(),
                    fields[4].to_string(),
                )
            })
            .collect()
    }

    /// A temporary directory of this test (removed by the test).
    fn temp_dir(tag: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "losat-blastn-db-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn test_hsp(s_idx: u32, raw_score: i32) -> BlastnHsp {
        BlastnHsp {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            e_value: 1e-5,
            bit_score: 30.0,
            num_ident: 10,
            query_frame: 1,
            query_length: 20,
            q_idx: 0,
            s_idx,
            raw_score,
            internal_q_offset_0: 0,
            internal_q_end_0: 10,
            internal_s_offset_0: 0,
            internal_s_end_0: 10,
            internal_query_context_offset: 0,
            gap_info: None,
            num_positives: 10,
        }
    }

    fn hit_list(lists: Vec<(u32, f64)>, count: usize) -> Option<BlastnHitList> {
        let mut hit_list = BlastnHitList::new(10);
        hit_list.hsplist_array = lists
            .into_iter()
            .map(|(oid, best_evalue)| BlastnHspList {
                oid,
                query_index: 0,
                hsps: vec![test_hsp(oid, 60)],
                best_evalue,
            })
            .collect();
        hit_list.hsplist_count = count;
        Some(hit_list)
    }

    // NCBI BLAST+ 2.17.0 `blastdbcmd -info` of the fixture
    // (LOSAT/tests/seqdb/expected/nucl_mv.info.tsv), and align_format_util.cpp:385-391.
    #[test]
    fn database_numbers_come_from_the_index_files() {
        let database = open("nucl_mv");
        assert_eq!(database.metadata().db_num_seqs, 12);
        assert_eq!(database.metadata().db_len_total, 850);
        let shown = database.pairwise_database();
        assert_eq!(shown.definition, b"nucleotide in volumes");
        assert_eq!(shown.date, "Oct 12, 2026  4:12 PM");
        assert_eq!((shown.num_sequences, shown.total_letters), (12, 850));
        assert_eq!(database.name(), fixture("nucl_mv"));
    }

    // blast_app_util.hpp:185-203 and seqdb.cpp:155-165: NCBI's errors with exit status 2,
    // LOSAT's refusals in LOSAT's wording, naming BLASTN.
    // (A path in an error text: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn open_errors_are_ncbis_database_errors_or_losats_refusals() {
        let native = |error: anyhow::Error| {
            let native = error
                .downcast_ref::<NativeError>()
                .unwrap_or_else(|| panic!("not a native error: {error:#}"));
            (native.exit, native.message.clone())
        };
        let settings = ApplicationSettings::default();
        let Err(error) = BlastnDatabase::open("", &settings) else {
            panic!("an empty name opened")
        };
        assert_eq!(
            native(error),
            (
                2,
                "BLAST Database error: Database name is required.\n".to_string()
            )
        );
        // A protein database is not found as a nucleotide one.
        for name in ["no_such_database", "prot_v4"] {
            let Err(error) = BlastnDatabase::open(&fixture(name), &settings) else {
                panic!("{name} opened")
            };
            let (exit, message) = native(error);
            assert_eq!(exit, 2);
            assert!(
                message.starts_with(
                    "BLAST Database error: No alias or index file found for nucleotide database ["
                ),
                "{message}"
            );
        }
        // The reader's refusals name BLASTN.
        let error = database_error(SeqDbError::Unsupported {
            path: "/data/db.nal".into(),
            message: "the alias file filters sequences".to_string(),
        });
        assert!(error.downcast_ref::<NativeError>().is_none());
        assert!(
            format!("{error:#}").ends_with("; this database is not supported by LOSAT's BLASTN"),
            "{error:#}"
        );
        let read = search_read_error(SeqDbError::Open {
            path: "/data/db.nsq".into(),
            source: std::io::Error::from_raw_os_error(13),
        })
        .downcast::<PreliminaryReadError>()
        .expect("NCBI's error of the search stage");
        assert_eq!(
            read.to_string(),
            "BLAST engine error: Database memory map file error"
        );
    }

    /// A kernel that records what it is given and keeps every third sequence.
    fn recording_kernel(
        seen: &std::sync::Mutex<Vec<(u32, usize, Vec<u8>)>>,
    ) -> impl Fn(&mut u32, u32, BlastnSubjectRef<'_>) -> Option<u32> + Sync + Send + '_ {
        move |calls, ordinal, subject| {
            *calls += 1;
            assert!(subject.letters.is_none() && subject.blastna.is_none());
            seen.lock()
                .unwrap()
                .push((ordinal, subject.length, subject.packed.to_vec()));
            (ordinal % 3 == 0).then_some(ordinal)
        }
    }

    fn run_blocks(
        database: &BlastnDatabase,
        threads: usize,
        block_bytes: u64,
    ) -> (Vec<(u32, usize, Vec<u8>)>, Vec<u32>) {
        let seen = std::sync::Mutex::new(Vec::new());
        let mut merged = Vec::new();
        crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
            search_blocks_of(
                database,
                pool,
                threads > 1,
                block_bytes,
                2,
                &|| 0u32,
                &recording_kernel(&seen),
                &mut |ordinal, kept| {
                    assert_eq!(ordinal, kept);
                    merged.push(ordinal);
                },
            )
        })
        .unwrap();
        let mut seen = seen.into_inner().unwrap();
        seen.sort_by_key(|entry| entry.0);
        (seen, merged)
    }

    // blast_engine.c:1407-1426 (the sequences in the order of their ordinal ids) and
    // blast_engine.c:1401-1405 (the stored 2-bit data): the block driver gives the kernel
    // every sequence once with its stored data, and merges what it keeps in ordinal order,
    // whatever the block size and the thread count, across the six volumes of the fixture.
    #[test]
    fn blocks_give_every_sequence_once_and_merge_in_ordinal_order() {
        let database = open("nucl_mv");
        let single = open("nucl_v4");
        let entries = entries();
        let (reference, merged) = run_blocks(&database, 1, u64::MAX);
        assert_eq!(reference.len(), 12);
        for (ordinal, (seen_ordinal, length, packed)) in reference.iter().enumerate() {
            assert_eq!(*seen_ordinal as usize, ordinal);
            assert_eq!(*length, entries[ordinal].0, "sequence {ordinal}");
            // The stored layout of the subject set's packed letters: length / 4 + 1 bytes,
            // the count of the last byte's bases in its two low bits.
            assert_eq!(packed.len(), length / 4 + 1, "sequence {ordinal}");
            assert_eq!(usize::from(packed[length / 4] & 3), length % 4);
        }
        assert_eq!(merged, [0, 3, 6, 9]);
        assert_eq!(run_blocks(&single, 1, u64::MAX).0, reference);
        for threads in [1usize, 2, 4] {
            for block_bytes in [1u64, 16, 64, 1 << 20] {
                for db in [&database, &single] {
                    let (seen, order) = run_blocks(db, threads, block_bytes);
                    assert_eq!(
                        seen, reference,
                        "{threads} threads, blocks of {block_bytes}"
                    );
                    assert_eq!(order, merged, "{threads} threads, blocks of {block_bytes}");
                }
            }
        }
    }

    // blast_traceback.c:819-828 and 1567-1573: the traceback reads each listed sequence
    // again by ordinal id, in BLASTNA with its ambiguity codes; the results come in the
    // order of the listed ordinal ids at every thread count.
    #[test]
    fn listed_sequences_are_read_again_with_their_ambiguity_codes() {
        let database = open("nucl_mv");
        let entries = entries();
        // n02, n06 (every ambiguity code), n08 (all N), n12 (ambiguity at both ends).
        let listed: Vec<(u32, u32)> = vec![(1, 10), (5, 50), (7, 70), (11, 110)];
        for threads in [1usize, 2, 4] {
            let found = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                traceback_listed(
                    &database,
                    pool,
                    threads > 1,
                    &listed,
                    &|| (),
                    &|_: &mut (), ordinal: u32, subject: BlastnSubjectRef<'_>, kept: &u32| {
                        assert_eq!(*kept, ordinal * 10);
                        let blastna = subject.blastna.expect("the traceback's letters");
                        assert_eq!(blastna.len(), subject.length);
                        assert_eq!(subject.packed.len(), subject.length / 4 + 1);
                        let letters: String = blastna
                            .iter()
                            .map(|&code| BLASTNA_TO_IUPACNA[usize::from(code)] as char)
                            .collect();
                        (ordinal != 7).then_some((ordinal, letters))
                    },
                )
            })
            .unwrap();
            let expected: Vec<(u32, String)> = [1u32, 5, 11]
                .iter()
                .map(|&ordinal| (ordinal, entries[ordinal as usize].1.clone()))
                .collect();
            assert_eq!(found, expected, "{threads} threads");
        }
    }

    // blast_seqalign.cpp:1569-1582 and seqdbvol.cpp:2302-2310: NCBI opens a volume's
    // header file when the Seq-aligns of a batch's final hit lists are made, query by query
    // and in the order of each hit list; a list past the hit list's count is not read, and
    // a title read before is not read again.
    // (A path in the system's temporary directory: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn titles_are_read_for_the_final_hit_lists_in_ncbis_order() {
        let source = open("nucl_mv");
        // The volume of each ordinal id.
        let mut volume_of = Vec::new();
        for (volume, vol) in source.db.volumes().iter().enumerate() {
            volume_of.extend(std::iter::repeat_n(volume, vol.num_oids() as usize));
        }
        assert_eq!(volume_of.len(), 12);
        let last = 11u32;
        let first_volume_ordinal = 0u32;
        assert_ne!(volume_of[last as usize], volume_of[0]);
        // A sequence of a third volume, read in the first batch.
        let middle = (0..12u32)
            .find(|&ordinal| {
                volume_of[ordinal as usize] != volume_of[0]
                    && volume_of[ordinal as usize] != volume_of[last as usize]
            })
            .unwrap();
        let dir = temp_dir("nohdr");
        let missing = [
            format!("nucl_mv.{:02}.nhr", volume_of[0]),
            format!("nucl_mv.{:02}.nhr", volume_of[last as usize]),
        ];
        for entry in
            std::fs::read_dir(fixture("nucl_mv.nal").trim_end_matches("nucl_mv.nal")).unwrap()
        {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_string();
            if name.starts_with("nucl_mv.") && !missing.contains(&name) {
                std::fs::copy(&path, dir.join(&name)).unwrap();
            }
        }
        let name = dir.join("nucl_mv").to_str().unwrap().to_string();
        let database = BlastnDatabase::open(&name, &ApplicationSettings::default()).unwrap();
        let mut titles = DatabaseTitles::new();
        // Query 1 reports `middle`; ordinal id 0 (no header file) is past the count.
        let batch = vec![
            None,
            hit_list(vec![(middle, 1e-9), (first_volume_ordinal, 1e-3)], 1),
        ];
        read_titles(&database, &batch, &mut titles).unwrap();
        assert_eq!(titles.keys().copied().collect::<Vec<_>>(), [middle]);
        assert_eq!(
            titles[&middle].title,
            entries()[middle as usize].2.as_bytes()
        );
        // The next batch: query 0 names `middle` again (not read again), query 1 names the
        // last sequence before ordinal id 0: the first header file read is the last
        // volume's.
        let batch = vec![
            hit_list(vec![(middle, 1e-9)], 1),
            hit_list(vec![(last, 1e-9), (first_volume_ordinal, 1e-5)], 2),
        ];
        let error = read_titles(&database, &batch, &mut titles)
            .err()
            .expect("a missing header file")
            .downcast::<HeaderReadError>()
            .expect("NCBI's error, with the query that reached it");
        assert_eq!(error.query, 1);
        let native = error
            .error
            .downcast_ref::<NativeError>()
            .expect("NCBI's error");
        assert_eq!(native.exit, 2);
        // The volume is recorded under its directory normalised with symbolic links
        // followed (seqdbalias.cpp:851-862), as in BLASTP's test.
        assert_eq!(
            native.message,
            format!(
                "BLAST Database error: Error: File ({}) not found.\n",
                dir.canonicalize().unwrap().join(&missing[1]).display()
            )
        );
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // The reports name the reported sequences by rank (the order of the ordinal ids), with
    // the stored title, the id `gnl|BL_ORD_ID|n` and the stored letters in upper case with
    // their ambiguity codes (seqdbvol.cpp:1406-1414, blast_encoding.c:80-83), and the
    // tabular ids are the first words of the titles.
    #[test]
    fn reported_subjects_are_renamed_to_ranks_with_their_stored_records() {
        let database = open("nucl_mv");
        let entries = entries();
        let mut hit_lists = vec![
            hit_list(vec![(9, 1e-9), (5, 1e-5)], 2),
            None,
            hit_list(vec![(5, 1e-9), (11, 1e-5)], 2),
        ];
        let mut titles = DatabaseTitles::new();
        read_titles(&database, &hit_lists, &mut titles).unwrap();
        // A list left without HSPs (its HSPs had score 0) is not reported.
        hit_lists[2].as_mut().unwrap().hsplist_array[1].hsps.clear();
        let reported = reported_subjects(&database, &titles, &mut hit_lists, true).unwrap();
        assert_eq!(reported.ordinals, [5, 9]);
        let ranks = |hit_list: &Option<BlastnHitList>| -> Vec<(u32, Vec<u32>)> {
            hit_list
                .as_ref()
                .unwrap()
                .hsplist_array
                .iter()
                .map(|list| (list.oid, list.hsps.iter().map(|hsp| hsp.s_idx).collect()))
                .collect()
        };
        assert_eq!(ranks(&hit_lists[0]), [(1, vec![1]), (0, vec![0])]);
        assert_eq!(ranks(&hit_lists[2]), [(0, vec![0]), (11, vec![])]);
        for (rank, &ordinal) in reported.ordinals.iter().enumerate() {
            let record = &reported.records[rank];
            assert_eq!(record.local_id, format!("gnl|BL_ORD_ID|{ordinal}"));
            assert_eq!(record.title, entries[ordinal as usize].2.as_bytes());
            assert_eq!(record.seq(), entries[ordinal as usize].1.as_bytes());
        }
        // n10 was stored from lower case letters: the database has upper case only.
        assert!(reported.records[1]
            .seq()
            .iter()
            .all(|letter| !letter.is_ascii_lowercase()));
        assert_eq!(
            reported
                .tabular_subject_ids()
                .unwrap()
                .iter()
                .map(|id| id.as_ref())
                .collect::<Vec<_>>(),
            [&b"n06"[..], b"n10"]
        );
        // Without the letters (the tabular formats only): the records are empty.
        let mut hit_lists = vec![hit_list(vec![(9, 1e-9)], 1)];
        let reported = reported_subjects(&database, &titles, &mut hit_lists, false).unwrap();
        assert!(reported.records[0].seq().is_empty());
    }

    // NCBI writes an error message for a stored title byte outside ' ' ..= '~' each time
    // its formatter reads the definition line (`check_reported_title`).
    // (A path in the system's temporary directory: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn a_reported_title_with_a_byte_outside_the_visible_range_is_rejected() {
        let dir = temp_dir("badtitle");
        for ext in ["nin", "nsq"] {
            std::fs::copy(fixture("nucl_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
        }
        let mut header = std::fs::read(fixture("nucl_v4") + ".nhr").unwrap();
        let at = header
            .windows(13)
            .position(|window| window == b"n02 length 41")
            .expect("the title of n02");
        header[at + 5] = 0xE9;
        std::fs::write(dir.join("b.nhr"), header).unwrap();
        let name = dir.join("b").to_str().unwrap().to_string();
        let database = BlastnDatabase::open(&name, &ApplicationSettings::default()).unwrap();
        let mut titles = DatabaseTitles::new();
        let mut clean = vec![hit_list(vec![(0, 1e-9)], 1)];
        read_titles(&database, &clean, &mut titles).unwrap();
        reported_subjects(&database, &titles, &mut clean, false).unwrap();
        let mut bad = vec![hit_list(vec![(1, 1e-9)], 1)];
        read_titles(&database, &bad, &mut titles).unwrap();
        assert_eq!(titles[&1].title, b"n02 l#ngth 41");
        let error = reported_subjects(&database, &titles, &mut bad, false)
            .err()
            .expect("the bad byte");
        assert!(
            format!("{error:#}").starts_with(
                "database sequence 1 has a stored title with a byte outside ' ' to '~'"
            ),
            "{error:#}"
        );
        assert!(format!("{error:#}").ends_with("not supported by LOSAT's BLASTN"));
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // EXPERIMENT (LOSAT_X_RPTSTREAM): each batch's reports are written when the loop takes it
    // only when every definition line of the database reads and none is refused
    // (`x_titles_allow_streaming`), at one thread and in a pool. (A path in the system's
    // temporary directory: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn streaming_is_allowed_only_when_every_title_reads_and_is_reported() {
        let allowed = |database: &BlastnDatabase| -> [bool; 2] {
            [1, 4].map(|threads| {
                crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                    Ok(x_titles_allow_streaming(database, pool))
                })
                .unwrap()
            })
        };
        assert_eq!(allowed(&open("nucl_v4")), [true, true]);
        assert_eq!(allowed(&open("nucl_mv")), [true, true]);
        let dir = temp_dir("streamtitles");
        // A header file that cannot be read (the fifth volume's): the batch that reads one of
        // its titles fails where NCBI's does, after the reports of the batches before.
        for entry in
            std::fs::read_dir(fixture("nucl_mv.nal").trim_end_matches("nucl_mv.nal")).unwrap()
        {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_string();
            if name.starts_with("nucl_mv.") && name != "nucl_mv.04.nhr" {
                std::fs::copy(&path, dir.join(&name)).unwrap();
            }
        }
        let missing =
            BlastnDatabase::open(dir.join("nucl_mv").to_str().unwrap(), &Default::default())
                .unwrap();
        assert_eq!(allowed(&missing), [false, false]);
        // A stored title with a byte outside ' ' ..= '~' (LOSAT's refusal when it is reported).
        for ext in ["nin", "nsq"] {
            std::fs::copy(fixture("nucl_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
        }
        let mut header = std::fs::read(fixture("nucl_v4") + ".nhr").unwrap();
        let at = header
            .windows(13)
            .position(|window| window == b"n02 length 41")
            .expect("the title of n02");
        header[at + 5] = 0xE9;
        std::fs::write(dir.join("b.nhr"), header).unwrap();
        let bad =
            BlastnDatabase::open(dir.join("b").to_str().unwrap(), &Default::default()).unwrap();
        assert_eq!(allowed(&bad), [false, false]);
        drop((missing, bad));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
