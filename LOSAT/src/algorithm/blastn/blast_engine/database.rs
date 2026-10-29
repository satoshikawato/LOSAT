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
//!   the preliminary hit lists as soon as its block is handed over (LOSAT_X_BNDBSTREAM; without
//!   it, when its block's round is searched), so nothing but the hit lists is kept for the
//!   whole database;
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

use std::cell::RefCell;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

use anyhow::{Context, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::algorithm::blastn::coordination::SubjectMetadata;
use crate::algorithm::blastn::hsp::BlastnHitList;
use crate::algorithm::common::database::{self as shared, ParsedShown, ReportedHeader};
pub(super) use crate::algorithm::common::database::{HeaderReadError, PreliminaryReadError};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::db_subject::ParsedSubject;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{BlockRange, NucleotideAmbiguity, SeqDb, SeqDbError, SeqType, SequenceBlock};
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

// EXPERIMENT (LOSAT_X_BNDBCUT): the blocks of the preliminary stage cut
// by the database size and the thread count (a child module, so that it uses the private
// driver of this file).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_db_cut.rs"]
mod x_db_cut;

// EXPERIMENT (LOSAT_X_BNDBSTREAM): the parallel preliminary stage without rounds, the blocks'
// results handed to the hit lists in ordinal order while the threads search the next blocks (a
// child module, so that it uses the private block reader of this file).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_db_stream.rs"]
mod x_db_stream;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:75-78
// ```c
// const char BLASTNA_TO_IUPACNA[BLASTNA_SIZE] = {
//     'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K',
//     'W', 'S', 'B', 'D', 'H', 'V', 'N', '-'
// };
// ```
/// The IUPAC letter of a BLASTNA code.
const BLASTNA_TO_IUPACNA: [u8; 16] = *b"ACGTRYMKWSBDHVN-";

/// `SeqDb::read_block`, mapped with LOSAT_X_BNMAPPED.
fn x_read_block(db: &SeqDb, range: BlockRange) -> Result<SequenceBlock, SeqDbError> {
    if super::x_bn_mapped() {
        db.x_read_block_mapped(range)
    } else {
        db.read_block(range)
    }
}

/// `SeqDb::read_sequence`, mapped with LOSAT_X_BNMAPPED.
fn x_read_sequence(db: &SeqDb, ordinal: u32) -> Result<SequenceBlock, SeqDbError> {
    if super::x_bn_mapped() {
        db.x_read_sequence_mapped(ordinal)
    } else {
        db.read_sequence(ordinal)
    }
}

/// With LOSAT_X_BNMAPPED, gives a block read by `x_read_block` or `x_read_sequence` back to
/// `utils::x_mapped::X_BYTES`.
fn x_done_block(block: SequenceBlock) {
    if super::x_bn_mapped() {
        block.x_recycle();
    }
}

/// With LOSAT_X_BNMAPPED, gives letters of `x_blastna` back to `utils::x_mapped::X_BYTES`.
fn x_done_bytes(mut bytes: Vec<u8>) {
    if super::x_bn_mapped() {
        super::x_bn_mapped_reach(bytes.len());
        crate::utils::x_mapped::X_BYTES.give(&mut bytes);
    }
}

/// `SequenceBlock::nucleotide_blastna(0)`, mapped with LOSAT_X_BNMAPPED.
fn x_blastna(block: &SequenceBlock) -> Result<Vec<u8>, SeqDbError> {
    if super::x_bn_mapped() {
        block.x_nucleotide_blastna_mapped(0)
    } else {
        block.nucleotide_blastna(0)
    }
}

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
    /// The taxonomy database of this run, opened at the first name asked for (outfmt 15 and 16,
    /// and the batches' prefetch of outfmt 0, 5, 6, 7 and 17).
    taxdb: shared::TaxdbRun,
}

/// The length and the ambiguity runs of the database sequences that one report call wrote
/// letters of (`BlastnDatabase::sam_letters`), by ordinal id.
#[derive(Default)]
pub(super) struct SamLetterCache(RefCell<HashMap<u32, NucleotideAmbiguity>>);

impl BlastnDatabase {
    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:185
    // ```c++
    //             "\tLN:" + NStr::UInt8ToString(refseq.GetBioseqLength()));
    // ```
    /// The length in bases of each of the database sequences `ordinals` (outfmt 17: the
    /// Bioseq's length for `@SQ` and the clips of a record), from the index and the last byte
    /// of each sequence's packed bases; no letters are read.
    pub(super) fn sam_lengths(&self, ordinals: &[u32]) -> Result<Vec<usize>> {
        ordinals
            .iter()
            .map(|&ordinal| self.db.nucleotide_length(ordinal).map_err(database_error))
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-851
    // ```c++
    //             	CRef<CBlast_def_line_set> bdlRef =
    //             			CSeqDB::ExtractBlastDefline(subject_bh);
    //             	if(setSubjectIds) {
    //             		x_SetSubjectIds(subject_bh, bdlRef);
    //             	}
    //             	if(setSubjectTaxInfoAll) {
    //             		x_SetTaxInfoAll(subject_bh, bdlRef);
    //             	}
    // ```
    /// The taxonomy of the tabular fields (`staxid`, `staxids`;
    /// `report::tabular::TabularTaxonomy::of_deflines`) of the reported sequences `ordinals`, by
    /// rank: from the definition lines of a sequence with parsed Seq-ids (`parsed[rank]`), else
    /// from its stored ones (`SeqDb::defline_set`, which keeps whether a line has a `taxid`), for
    /// the ranks that `hit_lists` name; the others are left empty. A header that does not decode
    /// is refused before the report.
    pub(super) fn tabular_taxonomy(
        &self,
        ordinals: &[u32],
        parsed: &[Option<Arc<ParsedSubject>>],
        hit_lists: &[Option<BlastnHitList>],
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        use crate::report::tabular::TabularTaxonomy;
        let mut taxonomy = vec![TabularTaxonomy::default(); ordinals.len()];
        let mut read = vec![false; ordinals.len()];
        for hsp_list in hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| &hit_list.hsplist_array)
        {
            let Some(rank) = hsp_list.hsps.first().map(|hsp| hsp.s_idx as usize) else {
                continue;
            };
            if rank >= ordinals.len() || std::mem::replace(&mut read[rank], true) {
                continue;
            }
            taxonomy[rank] = match parsed.get(rank).and_then(Option::as_deref) {
                Some(subject) => TabularTaxonomy::of_deflines(subject.deflines()),
                None => {
                    let set = self
                        .db
                        .defline_set(ordinals[rank])
                        .map_err(database_error)?;
                    TabularTaxonomy::of_deflines(&set.deflines)
                }
            };
        }
        Ok(taxonomy)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:222-234
    // ```c++
    //             vect.GetSeqData(tgt_rg.GetFrom(), tgt_rg.GetToOpen(), seq_data);
    // ```
    /// The IUPAC letters `[from, to)` of database sequence `ordinal` (outfmt 17 with `SQ`):
    /// only the bytes of the range are read, and the sequence's length and ambiguity runs are
    /// read once per call and kept in `cache` (NCBI's `CSeqVector` also keeps them).
    pub(super) fn sam_letters(
        &self,
        cache: &SamLetterCache,
        ordinal: u32,
        from: usize,
        to: usize,
    ) -> Result<Vec<u8>> {
        let mut entries = cache.0.borrow_mut();
        let ambiguity = match entries.entry(ordinal) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => entry.insert(
                self.db
                    .nucleotide_ambiguity(ordinal)
                    .map_err(database_error)?,
            ),
        };
        self.db
            .nucleotide_iupac_range_with(ordinal, from, to, ambiguity)
            .map_err(database_error)
    }

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
        // A database with parsed Seq-ids is read (line DBID, ID2b), unless a variable changes
        // NCBI's choice of its Seq-ids.
        shared::check_seqids(&db, settings, PROGRAM)?;
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
            taxdb: shared::TaxdbRun::default(),
        })
    }

    /// The Blast-def-line taxid of each reported sequence by rank (`ordinals`: the ordinal id of
    /// each rank, `ReportedSubjects::ordinals`; `parsed`: the sequences with parsed Seq-ids,
    /// `ReportedSubjects::parsed`) and what the taxonomy database says about it, for the XML2
    /// and JSON formats (`algorithm::common::database::read_reported_taxonomy`); the sequences
    /// of `hit_lists` with HSPs are the reported ones.
    #[allow(clippy::type_complexity)]
    pub(super) fn read_taxids(
        &self,
        ordinals: &[u32],
        parsed: &[Option<Arc<ParsedSubject>>],
        hit_lists: &[Option<BlastnHitList>],
    ) -> Result<(
        Vec<i32>,
        Vec<crate::report::blast_xml::HitTaxonomy>,
        Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
        shared::TaxdbNotes,
    )> {
        let reported: std::collections::BTreeSet<usize> = hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| &hit_list.hsplist_array)
            .filter_map(|hsp_list| hsp_list.hsps.first())
            .map(|hsp| hsp.s_idx as usize)
            .collect();
        shared::read_reported_taxonomy(&self.db, &self.taxdb, ordinals, &reported, parsed, PROGRAM)
    }

    /// Whether `taxdb.bti` is on the search path, asked once for the run
    /// (`shared::TaxdbRun::found`): the batches' prefetch reads the taxonomy database only then
    /// (`prefetch_messages`).
    pub(super) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
    }

    /// NCBI's messages about the taxonomy database when the query batch whose final hit lists
    /// are `hit_lists` (after their HSPs of score 0 were dropped, still naming the sequences by
    /// ordinal id) prefetches the subjects it reports (`shared::prefetch_messages`): the
    /// sequences with HSPs in them, the ones `read_taxids` reads, with the definition lines of
    /// those with parsed Seq-ids from `titles` (`read_titles`). Empty, and nothing read, without
    /// `taxdb.bti` on the search path. A taxonomy database that LOSAT cannot read as NCBI does,
    /// and a header whose taxid does not decode, are refused.
    pub(super) fn prefetch_messages(
        &self,
        titles: &DatabaseTitles,
        hit_lists: &[Option<BlastnHitList>],
    ) -> Result<Vec<u8>> {
        let subjects = hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| {
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                &hit_list.hsplist_array[..count]
            })
            .filter(|hsp_list| !hsp_list.hsps.is_empty())
            .map(|hsp_list| {
                let parsed = match titles.get(&hsp_list.oid) {
                    Some(ReportedHeader::Parsed(subject)) => Some(&**subject),
                    _ => None,
                };
                (hsp_list.oid, parsed)
            });
        shared::prefetch_messages(&self.db, &self.taxdb, subjects, PROGRAM)
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
        x_fenced: None,
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
///
/// EXPERIMENT (LOSAT_X_BNDBCUT): with more than one thread the blocks are cut by the
/// database size and the thread count, and a failed search returns the one-thread order's
/// error (`x_db_cut.rs`, exactness argument in its header).
///
/// EXPERIMENT (LOSAT_X_BNDBSTREAM): with more than one thread there are no rounds: the
/// threads claim the blocks in ordinal order and one of them hands the results over in that
/// order while the others search the next blocks, at most a window of blocks ahead
/// (`x_db_stream.rs`, exactness argument in its header).
#[allow(clippy::too_many_arguments)]
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
    M: FnMut(u32, R) + Send,
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
    // EXPERIMENT (LOSAT_X_BNDBCUT): `x_db_cut.rs`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_db_cut::x_bn_db_cut() {
        return x_db_cut::x_search_blocks_cut(
            database,
            pool,
            parallel,
            blocks_per_thread,
            new_scratch,
            search,
            merge,
        );
    }
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

// No NCBI counterpart (scheduling only): the rounds of LOSAT's parallel stage; the merges in
// ordinal order are NCBI's one-thread order (blast_engine.c:1553-1554, quoted above
// `search_blocks`).
/// The rounds of `search_blocks_of` with `parallel` (without LOSAT_X_BNDBSTREAM):
/// `blocks_per_thread` blocks of at most `block_bytes` per thread
/// searched with a scratch per worker (`map_init`), and the round's results merged in ordinal
/// order when the whole round is searched.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn search_rounds<T, R, I, K, M>(
    db: &SeqDb,
    pool: &SearchPool<'_>,
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
                    let block = x_read_block(db, range).map_err(search_read_error)?;
                    let mut kept = Vec::new();
                    for i in 0..block.len() {
                        let ordinal = range.first_ordinal + i as u32;
                        if let Some(result) = search(scratch, ordinal, block_subject(&block, i)?) {
                            kept.push((ordinal, result));
                        }
                    }
                    x_done_block(block);
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
    Ok(())
}

/// `search_blocks` with blocks of at most `block_bytes` and rounds of `blocks_per_thread`
/// blocks per thread (with LOSAT_X_BNDBSTREAM, the window of `x_db_stream.rs` instead of
/// rounds).
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
    M: FnMut(u32, R) + Send,
{
    let db = &database.db;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if parallel {
        // EXPERIMENT (LOSAT_X_BNDBSTREAM): `x_db_stream.rs`.
        if x_db_stream::x_bn_db_stream() {
            crate::utils::xreach::x_reach!("LOSAT_X_BNDBSTREAM", true);
            return x_db_stream::x_search_blocks_stream(
                db,
                pool,
                block_bytes,
                new_scratch,
                search,
                merge,
            );
        }
        return search_rounds(
            db,
            pool,
            block_bytes,
            blocks_per_thread,
            new_scratch,
            search,
            merge,
        );
    }
    // One thread: the blocks are searched one after another (the reference path).
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_db_stream::x_bn_db_stream() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNDBSTREAM", false);
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, parallel, blocks_per_thread);

    let mut scratch = new_scratch();
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block = x_read_block(db, range).map_err(search_read_error)?;
        for i in 0..block.len() {
            let ordinal = range.first_ordinal + i as u32;
            if let Some(result) = search(&mut scratch, ordinal, block_subject(&block, i)?) {
                merge(ordinal, result);
            }
        }
        x_done_block(block);
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
///
/// EXPERIMENT (LOSAT_X_BNTBRANGE): a sequence of at least `x_tb_range_min` bases is decoded in
/// the ranges of its kept HSPs only (`spans`: `[from, to)` on the sequence; `x_partial_ranges`),
/// with fences; when a traceback alignment meets one, the sequence is decoded whole and traced
/// again from `kept`.
pub(super) fn traceback_listed<T, P, R, I, K, S>(
    database: &BlastnDatabase,
    pool: &SearchPool<'_>,
    parallel: bool,
    listed: &[(u32, P)],
    new_scratch: &I,
    traceback: &K,
    spans: &S,
) -> Result<Vec<R>>
where
    P: Sync,
    R: Send + std::fmt::Debug,
    I: Fn() -> T + Sync + Send,
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>, &P) -> Option<R> + Sync + Send,
    S: Fn(&P) -> Vec<(usize, usize)> + Sync + Send,
{
    let db = &database.db;
    let x_ranged = super::x_bn_tb_range();
    let x_min = super::x_tb_range_min();
    let one = |scratch: &mut T, &(ordinal, ref kept): &(u32, P)| -> Result<Option<R>> {
        let stored = x_read_sequence(db, ordinal).map_err(database_error)?;
        let length = stored.nucleotide_len(0).map_err(database_error)?;
        let x_in_ranges = x_ranged && length >= x_min;
        if x_ranged {
            crate::utils::xreach::x_reach!("LOSAT_X_BNTBRANGE", x_in_ranges);
        }
        if x_in_ranges {
            let ranges = x_partial_ranges(spans(kept), length);
            let blastna = stored
                .x_nucleotide_blastna_ranges(0, &ranges)
                .map_err(database_error)?;
            let fenced = std::sync::atomic::AtomicBool::new(false);
            let subject = BlastnSubjectRef {
                length,
                packed: stored.nucleotide_packed(0),
                letters: None,
                blastna: Some(&blastna[1..=length]),
                store: None,
                x_fenced: Some(&fenced),
            };
            let result = traceback(scratch, ordinal, subject, kept);
            drop(blastna);
            let fenced = fenced.into_inner();
            x_tb_range_count(length, &ranges, fenced);
            if !fenced {
                x_done_block(stored);
                return Ok(result);
            }
        }
        let result = x_traceback_whole(&stored, length, scratch, ordinal, kept, traceback)?;
        x_done_block(stored);
        Ok(result)
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

/// The traceback of `traceback_listed` on the whole sequence of `stored`.
fn x_traceback_whole<T, P, R, K>(
    stored: &SequenceBlock,
    length: usize,
    scratch: &mut T,
    ordinal: u32,
    kept: &P,
    traceback: &K,
) -> Result<Option<R>>
where
    K: Fn(&mut T, u32, BlastnSubjectRef<'_>, &P) -> Option<R>,
{
    let blastna = x_blastna(stored).map_err(database_error)?;
    let subject = BlastnSubjectRef {
        length,
        packed: stored.nucleotide_packed(0),
        letters: None,
        // The kernel's BLASTNA has no sentinel bytes.
        blastna: Some(&blastna[1..=length]),
        store: None,
        x_fenced: None,
    };
    let result = traceback(scratch, ordinal, subject, kept);
    x_done_bytes(blastna);
    Ok(result)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_seqsrc.c:467-519
// ```c
//     begin = MAX(0, begin - BLAST_SEQSRC_OVERHANG);
//     end = MIN(end + BLAST_SEQSRC_OVERHANG, len);
//     arg->ranges[arg->num_ranges++] = begin;
//     arg->ranges[arg->num_ranges++] = end;
// ...
//     arg->num_ranges /= 2;
//     if (arg->num_ranges <= 1) return;
//     qsort(arg->ranges, arg->num_ranges, 2*sizeof(Int4), BeginCompareHSPs);
//     i=0;
//     for (j=1; j<arg->num_ranges; ++j) {
//         Int4 begin = arg->ranges[j*2];
//         Int4 end = arg->ranges[j*2+1];
//         ASSERT(begin >= arg->ranges[2*i]);
//         if (begin > arg->ranges[2*i+1] + BLAST_SEQSRC_MINGAP) {
//             /* insert as a new range */
//             ++i;
//             arg->ranges[i*2] = begin;
//             arg->ranges[i*2+1] = end;
//         } else if (end > arg->ranges[2*i+1]) {
//             /* merge into the previous range */
//             arg->ranges[i*2+1] = end;
//         }
//     }
//     arg->num_ranges = i+1;
// ```
// (blast_seqsrc.h:203-204: BLAST_SEQSRC_MINGAP 1024, BLAST_SEQSRC_OVERHANG 1024.) The ranges
// come from every HSP that the batch keeps of the sequence (`BLAST_SetupPartialFetching`,
// blast_traceback.c:1320-1370); the merged ranges do not depend on the order of equal starts.
/// EXPERIMENT (LOSAT_X_BNTBRANGE): NCBI's partial-fetch ranges of a sequence of `length` bases
/// whose HSPs cover `spans` (`[from, to)`): each widened by 1024 bases on both sides within the
/// sequence, sorted by start, and merged into the range before unless it starts more than 1024
/// bases after that one's end. Sorted and disjoint.
pub(super) fn x_partial_ranges(
    spans: impl IntoIterator<Item = (usize, usize)>,
    length: usize,
) -> Vec<(usize, usize)> {
    x_partial_ranges_padded(spans, length, x_tb_range_pad())
}

/// EXPERIMENT (LOSAT_X_BNTBRANGEPAD=<bases>, a test knob): NCBI's 1024 of `x_partial_ranges`
/// (both the widening and the gap), 16 to 1,048,576; a small one makes the traceback meet
/// fences, so that the tests reach the whole sequence's second traceback. Not below 16: the
/// start of a gapped alignment is found by reading up to 12 bases on each side of the seed
/// (`blast_get_start_for_gapped_alignment_nucl`), outside the kernels' fence tests.
fn x_tb_range_pad() -> usize {
    static PAD: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *PAD.get_or_init(|| {
        std::env::var("LOSAT_X_BNTBRANGEPAD")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|pad| (16..=1 << 20).contains(pad))
            .unwrap_or(1024)
    })
}

fn x_partial_ranges_padded(
    spans: impl IntoIterator<Item = (usize, usize)>,
    length: usize,
    pad: usize,
) -> Vec<(usize, usize)> {
    #[allow(non_snake_case)]
    let (OVERHANG, MINGAP) = (pad, pad);
    let mut ranges: Vec<(usize, usize)> = spans
        .into_iter()
        .map(|(from, to)| {
            (
                from.saturating_sub(OVERHANG).min(length),
                to.saturating_add(OVERHANG).min(length),
            )
        })
        .collect();
    ranges.sort_unstable_by_key(|range| range.0);
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
    for (from, to) in ranges {
        match merged.last_mut() {
            Some(last) if from <= last.1.saturating_add(MINGAP) => last.1 = last.1.max(to),
            _ => merged.push((from, to)),
        }
    }
    merged
}

/// EXPERIMENT (LOSAT_X_BNTBRANGE): sequences decoded in ranges, those traced again whole, their
/// bases and the bases of their ranges (`x_tb_range_print_summary`).
static X_TB_RANGE_COUNTS: [std::sync::atomic::AtomicU64; 4] =
    [const { std::sync::atomic::AtomicU64::new(0) }; 4];

fn x_tb_range_count(length: usize, ranges: &[(usize, usize)], fenced: bool) {
    use std::sync::atomic::Ordering::Relaxed;
    let in_ranges: usize = ranges.iter().map(|range| range.1 - range.0).sum();
    X_TB_RANGE_COUNTS[0].fetch_add(1, Relaxed);
    X_TB_RANGE_COUNTS[1].fetch_add(u64::from(fenced), Relaxed);
    X_TB_RANGE_COUNTS[2].fetch_add(length as u64, Relaxed);
    X_TB_RANGE_COUNTS[3].fetch_add(in_ranges as u64, Relaxed);
}

/// EXPERIMENT (LOSAT_X_BNTBRANGELOG=1): one stderr line with the counts of LOSAT_X_BNTBRANGE.
pub(super) fn x_tb_range_print_summary() {
    use std::sync::atomic::Ordering::Relaxed;
    if std::env::var_os("LOSAT_X_BNTBRANGELOG").is_none_or(|v| v != "1") {
        return;
    }
    let [ranged, fenced, bases, in_ranges] = X_TB_RANGE_COUNTS.each_ref().map(|c| c.load(Relaxed));
    eprintln!(
        "[X_BNTBRANGE] sequences in ranges {ranged}, traced again whole {fenced}, bases {bases}, in ranges {in_ranges}"
    );
}

/// The headers of the database sequences that the final hit lists of the query batches searched
/// so far name, by ordinal id (`read_titles`): the stored title, or the definition lines of a
/// sequence with parsed Seq-ids.
pub(super) type DatabaseTitles = BTreeMap<u32, ReportedHeader>;

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
/// query batch name, in NCBI's order, into `titles` (`common::database::read_header`). A reader
/// error of NCBI's own is a `HeaderReadError` (it fails the batch); LOSAT's refusals are errors
/// of their own.
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
            let title = shared::read_header(&database.db, hsp_list.oid, query, PROGRAM)?;
            titles.insert(hsp_list.oid, title);
        }
    }
    Ok(())
}

/// The database sequences that the reports show, by rank (`reported_subjects`).
pub(super) struct ReportedSubjects {
    /// The ordinal id of each rank, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The record of each rank: the id `gnl|BL_ORD_ID|n`, the stored title (a sequence with
    /// parsed Seq-ids: the title of its first definition line), and the stored letters (empty
    /// when the reports do not show them).
    pub(super) records: Vec<FastaRecord>,
    /// The sequences with parsed Seq-ids by rank (empty when none is reported): their ids and
    /// titles come from their definition lines (`report::db_subject`).
    pub(super) parsed: Vec<Option<Arc<ParsedSubject>>>,
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
/// does not report as NCBI does is LOSAT's refusal, and so is one with parsed Seq-ids that the
/// reports `shown` do not write as NCBI (`common::database::check_reported_parsed`, and with
/// `shown.joined` `refuse_joined_hits`; `first_query`: the index of the query of `hit_lists[0]`).
pub(super) fn reported_subjects(
    database: &BlastnDatabase,
    titles: &DatabaseTitles,
    hit_lists: &mut [Option<BlastnHitList>],
    letters: bool,
    shown: ParsedShown,
    first_query: usize,
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
    let mut parsed: Vec<Option<Arc<ParsedSubject>>> = Vec::new();
    // Whether `taxdb.bti` is on the search path, looked up at the first sequence with parsed
    // Seq-ids that outfmt 0 shows.
    let mut taxdb_found: Option<bool> = None;
    for (rank, &ordinal) in ordinals.iter().enumerate() {
        let header = titles
            .get(&ordinal)
            .with_context(|| format!("database sequence {ordinal} has no title read"))?;
        let title = match header {
            ReportedHeader::Stored(stored) => {
                shared::check_reported_title(stored.replaced_bytes, ordinal, PROGRAM)?;
                stored.title.clone()
            }
            ReportedHeader::Parsed(subject) => {
                let taxdb =
                    shown.pairwise && *taxdb_found.get_or_insert_with(|| database.db.taxdb_found());
                shared::check_reported_parsed(subject, shown, false, taxdb, ordinal, PROGRAM)?;
                if parsed.is_empty() {
                    parsed = vec![None; ordinals.len()];
                }
                parsed[rank] = Some(Arc::clone(subject));
                subject.first_title().to_vec()
            }
        };
        let sequence = if letters {
            let block = x_read_sequence(&database.db, ordinal).map_err(database_error)?;
            let length = block.nucleotide_len(0).map_err(database_error)?;
            let blastna = x_blastna(&block).map_err(database_error)?;
            let letters = blastna[1..=length]
                .iter()
                .map(|&code| BLASTNA_TO_IUPACNA[usize::from(code & 0xF)])
                .collect();
            x_done_bytes(blastna);
            x_done_block(block);
            letters
        } else {
            Vec::new()
        };
        records.push(FastaRecord {
            local_id: format!("gnl|BL_ORD_ID|{ordinal}"),
            title,
            sequence,
            warnings: Vec::new(),
            source: Default::default(),
        });
    }
    if shown.joined && !parsed.is_empty() {
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            shared::refuse_joined_hits(
                first_query + query,
                hit_list
                    .hsplist_array
                    .iter()
                    .filter(|hsp_list| !hsp_list.hsps.is_empty())
                    .map(|hsp_list| {
                        let rank = hsp_list.oid as usize;
                        Some((ordinals[rank], parsed[rank].as_deref()?))
                    }),
                PROGRAM,
            )?;
        }
    }
    Ok(ReportedSubjects {
        ordinals,
        records,
        parsed,
    })
}

impl ReportedSubjects {
    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every
    /// reported sequence, by rank (`algorithm::common::database::tabular_subject_ids`). A
    /// sequence whose id LOSAT cannot write as NCBI is an error before any format is
    /// written.
    pub(super) fn tabular_subject_ids(&self) -> Result<Vec<Arc<[u8]>>> {
        shared::tabular_subject_ids(&self.records, &self.ordinals, &self.parsed, false, PROGRAM)
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:310-312
// ```c++
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
// ```
/// EXPERIMENT (LOSAT_X_RPTSTREAM): which of LOSAT's refusals of the database's stored data could
/// come in a batch after the first, where it comes before any of the search's reports without
/// the switch (`x_rpt_stream`; the same test as blastp's block search,
/// `blastp::blast_engine::database::x_titles_allow_blocks`), or `None`:
/// - `Off::Titles`: a definition line that does not read as the reports read those of the
///   reported sequences (`read_title`: NCBI's reader, one definition line), holds a byte outside
///   ' ' ..= '~' (`check_reported_title`) or gives no tabular subject id
///   (`ReportedSubjects::tabular_subject_ids`);
/// - `Off::Ambiguity`: a sequence whose ambiguity data the traceback cannot decode
///   (`traceback_listed` and `reported_subjects` read them with
///   `SequenceBlock::nucleotide_blastna`, first for the batch whose lists name the sequence; the
///   preliminary stage reads the 2-bit letters only);
/// - `Off::Taxonomy`, with `taxids` (the batches' prefetch reads the taxonomy database: outfmt 0,
///   6, 7 and 17 with `taxdb.bti` on the search path, `BlastnDatabase::prefetch_messages`, DBID
///   item T): a lookup of the taxonomy database that could be refused
///   (`shared::taxdb_may_refuse` without the scientific names, which only XML2 and JSON write
///   and XML is never streamed). The prefetch's refusal of NCBI's taxonomy cache
///   (`shared::TaxdbRun::cache_collision`) comes at the run's first prefetch, in the first batch
///   before its reports. A definition line whose taxid does not decode is `Off::Titles` then:
///   the titles' pass decodes the taxids too (`SeqDb::visit_deflines_taxids`; one sequence at a
///   time, `SeqDb::defline_taxids`), as blastp's block search does.
/// Every ordinal id is read in `pool` (the titles, then the ambiguity data of the sequences that
/// have any; `SeqDb::check_ambiguity` reads only those bytes); a reader error is a refusal too
/// (the reports then stay after the loop, and the batch that meets the error fails where it does
/// without the switch). EXPERIMENT (LOSAT_X_DBIDXREAD): both are read in ranges of
/// `HEADER_RANGE_SEQS` sequences (`SeqDb::visit_deflines`, `SeqDb::visit_ambiguity`, the same
/// results as `read_title` and `SeqDb::check_ambiguity`) instead of one sequence at a time.
/// A database with parsed Seq-ids gives `Off::Titles` at its first definition line (`read_title`
/// reads ordinal ids alone), so its reports come after the whole search, as BLASTP's do without
/// its block search (line DBID, DECISIONS 9 (6)), except SAM (`Ready::sam_despite`).
pub(super) fn x_streaming_refusal(
    database: &BlastnDatabase,
    taxids: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<super::x_rpt_stream::Off> {
    let db = &database.db;
    let title_allowed = |ordinal: u32, stored: Result<shared::StoredTitle>| -> bool {
        match stored {
            Ok(stored) => {
                !stored.replaced_bytes
                    && shared::tabular_subject_id(&stored.title, ordinal, false, PROGRAM).is_ok()
            }
            Err(_) => false,
        }
    };
    if crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
        let all = |allowed: &(dyn Fn(&BlockRange) -> bool + Sync)| -> bool {
            #[cfg(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            ))]
            if pool.enabled() {
                return pool.install(|| ranges.par_iter().all(allowed));
            }
            let _ = pool;
            ranges.iter().all(allowed)
        };
        // With `taxids` a record whose taxid does not decode is an `Err`
        // (`visit_deflines_taxids`).
        let titles = |range: &BlockRange| -> bool {
            let visit = |ordinal: u32, deflines| {
                title_allowed(ordinal, shared::stored_title(deflines, ordinal, 0, PROGRAM))
            };
            if taxids {
                db.visit_deflines_taxids(*range, visit)
            } else {
                db.visit_deflines(*range, visit)
            }
        };
        if !all(&titles) {
            return Some(super::x_rpt_stream::Off::Titles);
        }
        let ambiguity =
            |range: &BlockRange| -> bool { db.visit_ambiguity(*range, |_, result| result.is_ok()) };
        if !all(&ambiguity) {
            return Some(super::x_rpt_stream::Off::Ambiguity);
        }
        return x_taxonomy_refusal(database, taxids);
    }
    // One sequence at a time (`LOSAT_X_DBIDXREAD=0`): the taxid after the title.
    let title = |ordinal: u32| -> bool {
        title_allowed(ordinal, shared::read_title(db, ordinal, 0, PROGRAM))
            && (!taxids || db.defline_taxids(ordinal).is_ok())
    };
    let ambiguity = |ordinal: u32| -> bool { db.check_ambiguity(ordinal).is_ok() };
    let all = |allowed: &(dyn Fn(u32) -> bool + Sync)| -> bool {
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
    };
    if !all(&title) {
        return Some(super::x_rpt_stream::Off::Titles);
    }
    if !all(&ambiguity) {
        return Some(super::x_rpt_stream::Off::Ambiguity);
    }
    x_taxonomy_refusal(database, taxids)
}

/// EXPERIMENT (LOSAT_X_RPTSTREAM): `Off::Taxonomy` of `x_streaming_refusal`, asked with `taxids`
/// only (`shared::taxdb_may_refuse` without the names, `false` without `taxdb.bti` on the search
/// path).
fn x_taxonomy_refusal(database: &BlastnDatabase, taxids: bool) -> Option<super::x_rpt_stream::Off> {
    (taxids && shared::taxdb_may_refuse(&database.db, false))
        .then_some(super::x_rpt_stream::Off::Taxonomy)
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

    // DECISIONS 3 Q8 of line DBID: NCBI's choice of the Seq-ids of a database with parsed
    // Seq-ids depends on SEQ_ID_PREFER_ACCESSION_OVER_GI and SEQ_ID_AVOID_GI.
    #[test]
    fn a_seq_id_variable_refuses_a_database_with_parsed_seqids_only() {
        let path = |name: &str| {
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join(name)
                .to_str()
                .expect("fixture path is UTF-8")
                .to_string()
        };
        let settings = ApplicationSettings {
            seq_id_rank_variable: Some("SEQ_ID_PREFER_ACCESSION_OVER_GI"),
            ..ApplicationSettings::default()
        };
        let parsed = path("tests/db_search/ids/nucl_v5");
        BlastnDatabase::open(&parsed, &ApplicationSettings::default())
            .expect("a database with parsed Seq-ids opens");
        let Err(error) = BlastnDatabase::open(&parsed, &settings) else {
            panic!("a database with parsed Seq-ids opened with SEQ_ID_PREFER_ACCESSION_OVER_GI")
        };
        assert_eq!(
            format!("{error:#}"),
            "the environment variable SEQ_ID_PREFER_ACCESSION_OVER_GI, which changes the Seq-ids that NCBI BLAST+ chooses for a database with parsed Seq-ids, is not supported by LOSAT's BLASTN"
        );
        BlastnDatabase::open(&path("tests/seqdb/db/nucl_v4"), &settings)
            .expect("a database without parsed Seq-ids opens with SEQ_ID_PREFER_ACCESSION_OVER_GI");
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
                    &|_: &u32| Vec::new(),
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

    // blast_seqsrc.c:467-519: each HSP widened by 1024 bases within the sequence, sorted by
    // start, merged unless it starts more than 1024 bases after the range before.
    #[test]
    fn partial_ranges_are_ncbis() {
        assert_eq!(x_partial_ranges_padded([], 100, 1024), vec![]);
        assert_eq!(
            x_partial_ranges_padded([(10, 20)], 100, 1024),
            vec![(0, 100)]
        );
        // 3976..7024, 5476..8024, 8976..11124: each within 1024 of the one before.
        assert_eq!(
            x_partial_ranges_padded([(6500, 7000), (5000, 6000), (10000, 10100)], 20000, 1024),
            vec![(3976, 11124)]
        );
        // 0..1124 and 3172..6124: 3172 > 1124 + 1024.
        assert_eq!(
            x_partial_ranges_padded([(4196, 5100), (0, 100)], 5500, 1024),
            vec![(0, 1124), (3172, 5500)]
        );
        // 0..1124 and 2148..: exactly 1024 apart merges.
        assert_eq!(
            x_partial_ranges_padded([(0, 100), (3172, 3200)], 10000, 1024),
            vec![(0, 4224)]
        );
        // A range inside the one before.
        assert_eq!(
            x_partial_ranges_padded([(5000, 9000), (6000, 6100)], 20000, 1024),
            vec![(3976, 10024)]
        );
    }

    // seqdbvol.cpp:1565-1595 (`GetAmbigPartialSeq`): the bases of the ranges as the whole
    // decode gives them, a fence before and after each range inside the sequence, the two
    // sentinels, and nothing written elsewhere.
    #[test]
    fn blastna_in_ranges_is_the_whole_decode_inside_them() {
        let database = open("nucl_mv");
        let mut checked = 0;
        for ordinal in 0..12u32 {
            let block = database.db.read_sequence(ordinal).unwrap();
            let length = block.nucleotide_len(0).unwrap();
            let whole = block.nucleotide_blastna(0).unwrap();
            let mut lists: Vec<Vec<(usize, usize)>> = vec![
                vec![(0, length)],
                vec![(1, length.saturating_sub(1).max(1))],
            ];
            if length >= 3 {
                lists.push(vec![(0, 1), (length - 1, length)]);
            }
            if length >= 12 {
                lists.push(vec![(2, 5), (7, 9), (11, length)]);
                lists.push(vec![(0, 3), (length / 2, length / 2 + 2)]);
            }
            for ranges in lists {
                let ranges: Vec<(usize, usize)> =
                    ranges.into_iter().filter(|r| r.0 < r.1).collect();
                let ranged = block.x_nucleotide_blastna_ranges(0, &ranges).unwrap();
                assert_eq!(ranged.len(), length + 2);
                let mut expected = vec![0u8; length + 2];
                for &(from, to) in &ranges {
                    if from > 0 {
                        expected[from] = crate::seqdb::FENCE_SENTRY;
                    }
                    if to < length {
                        expected[to + 1] = crate::seqdb::FENCE_SENTRY;
                    }
                }
                for &(from, to) in &ranges {
                    expected[from + 1..=to].copy_from_slice(&whole[from + 1..=to]);
                }
                expected[0] = whole[0];
                expected[length + 1] = whole[length + 1];
                assert_eq!(ranged, expected, "sequence {ordinal}, ranges {ranges:?}");
                checked += 1;
            }
        }
        assert!(checked >= 36);
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
            titles[&middle].stored().title,
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
        let reported = reported_subjects(
            &database,
            &titles,
            &mut hit_lists,
            true,
            ParsedShown::default(),
            0,
        )
        .unwrap();
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
        let reported = reported_subjects(
            &database,
            &titles,
            &mut hit_lists,
            false,
            ParsedShown::default(),
            0,
        )
        .unwrap();
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
        reported_subjects(
            &database,
            &titles,
            &mut clean,
            false,
            ParsedShown::default(),
            0,
        )
        .unwrap();
        let mut bad = vec![hit_list(vec![(1, 1e-9)], 1)];
        read_titles(&database, &bad, &mut titles).unwrap();
        assert_eq!(titles[&1].stored().title, b"n02 l#ngth 41");
        let error = reported_subjects(
            &database,
            &titles,
            &mut bad,
            false,
            ParsedShown::default(),
            0,
        )
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
    // only when every definition line of the database reads and none is refused, and the
    // ambiguity data of every sequence decode (`x_streaming_refusal`), at one thread and in a
    // pool. (A path in the system's temporary directory: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn streaming_needs_every_title_and_every_ambiguity_record_to_read() {
        use super::super::x_rpt_stream::Off;
        let refusal = |database: &BlastnDatabase| -> [Option<Off>; 2] {
            [1, 4].map(|threads| {
                crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                    Ok(x_streaming_refusal(database, false, pool))
                })
                .unwrap()
            })
        };
        assert_eq!(refusal(&open("nucl_v4")), [None, None]);
        assert_eq!(refusal(&open("nucl_mv")), [None, None]);
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
        assert_eq!(refusal(&missing), [Some(Off::Titles); 2]);
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
        assert_eq!(refusal(&bad), [Some(Off::Titles); 2]);
        // The ambiguity count of each sequence with ambiguity data raised by one entry (the
        // old format of n05 and the new format of n06): the traceback's decoding refuses it.
        let clean = open("nucl_v4");
        let volume = &clean.db.volumes()[0];
        let num_oids = volume.num_oids();
        let (seq_offsets, amb_offsets) = (
            volume.index().sequence_offsets(0, num_oids).unwrap(),
            volume.index().ambiguity_offsets(0, num_oids).unwrap(),
        );
        drop(clean);
        for (ordinal, tag) in [(4u32, "c"), (5, "d")] {
            let at = amb_offsets[ordinal as usize] as usize;
            assert!(
                at < seq_offsets[ordinal as usize + 1] as usize,
                "n0{}",
                ordinal + 1
            );
            for ext in ["nin", "nhr"] {
                std::fs::copy(
                    fixture("nucl_v4") + "." + ext,
                    dir.join(format!("{tag}.{ext}")),
                )
                .unwrap();
            }
            let mut sequences = std::fs::read(fixture("nucl_v4") + ".nsq").unwrap();
            let count = u32::from_be_bytes(sequences[at..at + 4].try_into().unwrap());
            sequences[at..at + 4].copy_from_slice(&(count + 1).to_be_bytes());
            std::fs::write(dir.join(format!("{tag}.nsq")), sequences).unwrap();
            let corrupt =
                BlastnDatabase::open(dir.join(tag).to_str().unwrap(), &Default::default()).unwrap();
            let error = corrupt
                .db
                .check_ambiguity(ordinal)
                .expect_err("the raised count");
            assert!(
                error
                    .to_string()
                    .contains("the ambiguity data are shorter than their count"),
                "{error}"
            );
            // The traceback's read of the same sequence fails with the same error.
            let block = corrupt.db.read_sequence(ordinal).unwrap();
            assert_eq!(
                block.nucleotide_blastna(0).unwrap_err().to_string(),
                error.to_string()
            );
            assert_eq!(
                refusal(&corrupt),
                [Some(Off::Ambiguity); 2],
                "n0{}",
                ordinal + 1
            );
        }
        // A sequence without ambiguity data reads nothing more (n01).
        open("nucl_v4").db.check_ambiguity(0).unwrap();
        drop((missing, bad));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// The fixture database `name` of `tests/db_search`, with the taxonomy database variant
    /// `taxdb` of `tests/seqdb/taxdb` (`None`: none) as `BLASTDB`, from the directory `cwd`.
    fn open_with_taxdb(name: &str, taxdb: Option<&str>, cwd: &std::path::Path) -> BlastnDatabase {
        use crate::seqdb::alias::SearchPath;
        let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
        let taxdb = taxdb.map(|variant| tests.join("seqdb/taxdb").join(variant));
        let search =
            SearchPath::new(cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None).unwrap();
        let path = tests.join("db_search").join(name);
        let db = SeqDb::open_in(path.to_str().unwrap(), SeqType::Nucleotide, &search).unwrap();
        BlastnDatabase {
            metadata: SubjectMetadata {
                db_len_total: db.total_length() as usize,
                db_num_seqs: db.num_seqs() as usize,
            },
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
        }
    }

    // DBID item T: a batch's prefetch looks up the sequences with HSPs in its final hit lists
    // (within the hit list size), once per volume and taxid over the batches of a run, and opens
    // the taxonomy database at the run's first lookup; nothing without `taxdb.bti`.
    // `x16_ntaxid_99990` (two sequences of taxid 99990) with `offset_last_past_end` (the last
    // entry, 99990, is past the end of `taxdb.btd`), `x16_ntaxid` (taxid 9606) with `bad_magic`.
    #[test]
    fn the_prefetch_looks_up_the_sequences_with_hsps_once_per_run() {
        let cwd = temp_dir("prefetch");
        let text = |database: &BlastnDatabase, hit_lists: &[Option<BlastnHitList>]| -> String {
            let titles = DatabaseTitles::new();
            String::from_utf8(database.prefetch_messages(&titles, hit_lists).unwrap()).unwrap()
        };
        let database = open_with_taxdb("x16_ntaxid_99990", Some("offset_last_past_end"), &cwd);
        assert!(database.taxdb_found());
        // Sequence 0 without HSPs (its HSPs of score 0 were dropped), sequence 1 past the hit
        // list size: nothing is fetched, the taxonomy database is not opened.
        let mut none = hit_list(vec![(0, 1e-9), (1, 1e-8)], 1);
        none.as_mut().unwrap().hsplist_array[0].hsps.clear();
        assert_eq!(text(&database, &[None, none]), "");
        let offset = "Error: [blastn] Error: Offset error at end of taxdb file.\n";
        assert_eq!(text(&database, &[hit_list(vec![(1, 1e-9)], 1)]), offset);
        // The same volume and taxid in a later batch.
        assert_eq!(
            text(&database, &[hit_list(vec![(0, 1e-9), (1, 1e-8)], 2)]),
            ""
        );
        let database = open_with_taxdb("x16_ntaxid", Some("bad_magic"), &cwd);
        assert_eq!(
            text(&database, &[hit_list(vec![(0, 1e-9), (1, 1e-8)], 2)]),
            "Error: [blastn] Error: Tax database file has wrong magic number.\n"
        );
        let database = open_with_taxdb("x16_ntaxid", None, &cwd);
        assert!(!database.taxdb_found());
        assert_eq!(
            text(&database, &[hit_list(vec![(0, 1e-9), (1, 1e-8)], 2)]),
            ""
        );
        std::fs::remove_dir_all(&cwd).unwrap();
    }

    // EXPERIMENT (LOSAT_X_RPTSTREAM), DBID item T: where the batches' prefetch reads the taxonomy
    // database (`taxids`), each batch's reports are written when the loop takes it only when no
    // lookup can be refused (`Off::Taxonomy`). `x16_ntaxid` (taxid 9606) with no taxonomy
    // database, ones whose lookups cannot be refused (a wrong magic number is a message; a
    // scientific name of 9606 that XML2 cannot write is not asked for by the prefetch) and one
    // whose record of 9606 is (an empty `taxdb.btd`, `shared::taxdb_may_refuse`), at one thread
    // and in a pool.
    #[test]
    fn streaming_with_the_prefetch_needs_lookups_that_cannot_be_refused() {
        use super::super::x_rpt_stream::Off;
        let cwd = temp_dir("streamtaxdb");
        for (taxdb, refused) in [
            (None, None),
            (Some("good"), None),
            (Some("bad_magic"), None),
            (Some("high_bytes_utf8"), None),
            (Some("btd_empty"), Some(Off::Taxonomy)),
        ] {
            let database = open_with_taxdb("x16_ntaxid", taxdb, &cwd);
            for threads in [1, 4] {
                crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                    assert_eq!(x_streaming_refusal(&database, false, pool), None);
                    assert_eq!(
                        x_streaming_refusal(&database, true, pool),
                        refused,
                        "{taxdb:?} threads {threads}"
                    );
                    Ok(())
                })
                .unwrap();
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }
}
