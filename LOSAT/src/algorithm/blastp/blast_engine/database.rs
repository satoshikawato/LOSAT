//! The database path of BLASTP (`-db NAME`): the search of an NCBI BLAST database read
//! in blocks through `crate::seqdb`, as NCBI BLAST+ 2.17.0 `blastp -db NAME
//! -num_threads 1` searches it (PD-LOSAT-DATABASE-SEARCH, P1 and P1.1).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the set-level numbers come from the index files (`BlastpDatabase`): the number of
//!   sequences, the total length, and 10 for the shortest sequence behind the hit-saving
//!   cutoff (`BLAST_SEQSRC_MINLENGTH`);
//! - the preliminary stage reads the sequences in blocks and searches them in the order
//!   of their ordinal ids (`search_blocks`); nothing but the hit lists is kept;
//! - the sequences and definition lines of the hit lists are read again by ordinal id
//!   after the preliminary stage (`ListedSubjects`) for the redo and the reports;
//! - the reports show the database (`BlastpDatabase::pairwise_database`, `name`) and
//!   the subject ids and titles of sequences without parsed Seq-ids
//!   (`ListedSubjects::tabular_subject_ids`, `stitles`).
//!
//! The per-subject kernel, the composition redo and the report writers are those of the
//! `-subject` path.

use std::sync::Arc;

use anyhow::{bail, Context, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::NativeError;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{BlockRange, SearchPath, SeqDb, SeqDbError, SeqType, SequenceBlock};
use crate::stats::search_space::SearchSpace;
use crate::stats::{compute_length_adjustment_ncbi, KarlinParams};
use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE};

use super::{
    ncbistdaa_to_ascii, update_blastp_preliminary_hit_lists, BlastpHitList, BlastpSubjectRef,
    BlastpSubjectResult, BlastpSubjectScratch, EncodedProtein,
};

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
// ```c
// #define BLAST_SEQSRC_MINLENGTH  10    /**< Default minimal sequence length */
// ```
/// The length that a database answers for its shortest sequence
/// (`BlastSeqSrcGetMinSeqLen`, seqdbimpl.cpp:128-136).
pub(super) const BLAST_SEQSRC_MINLENGTH: i32 = 10;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-119
// ```c
// const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
// '-','A','B','C','D','E','F','G','H','I','K','L','M',
// 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
// 'O', 'J'};
// ```
/// NCBIstdaa code of pyrrolysine (O).
const NCBISTDAA_O: u8 = 26;
/// NCBIstdaa code of the gap character (`-`), which is also the sentinel byte.
const NCBISTDAA_GAP: u8 = 0;

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

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`ResolvedBlastpArgs::db`: the `-db` value without white space at its ends).
pub(super) struct BlastpDatabase {
    db: SeqDb,
    name: String,
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbatlas.hpp:509-535
// ```c++
//     static const string GenerateSearchPath() {
//         string splitter;
//         string path;
// #if defined(NCBI_OS_UNIX)
//         splitter = ":";
// #else
//         splitter = ";";
// #endif
//         // Local directory first;
//         path  = CDirEntry::NormalizePath(CDir::GetCwd(),eFollowLinks);
//         path += splitter;
//         // Then, BLASTDB;
//         CNcbiEnvironment env;
//         path += CDirEntry::NormalizePath(env.Get("BLASTDB"),eFollowLinks);
//         path += splitter;
//         // Finally, the config file.
//         CNcbiApplicationAPI* app = CNcbiApplicationAPI::Instance();
//         if (app) {
//             const CNcbiRegistry& registry = app->GetConfig();
//             if (registry.HasEntry("BLAST", "BLASTDB")) {
//                 path += CDirEntry::NormalizePath(registry.Get("BLAST", "BLASTDB"),eFollowLinks);
//                 path += splitter;
//             }
//         }
//         return path;
//     }
// ```
/// Opens the protein database `name` (the `-db` value as NCBI keeps it,
/// `seqdb::db_argument_name`) for a BLASTP search: the one function through which BLASTP
/// opens a database. The name is looked for as NCBI looks for it: in the current
/// directory, then in the directories of the `BLASTDB` variable, then in those of the
/// registry entry `[BLAST] BLASTDB` (`settings.blastdb`, read from
/// `NCBI_CONFIG__BLAST__BLASTDB`, `blastp.ini` and `.ncbirc`).
pub(super) fn open_search_database(
    name: &str,
    settings: &ApplicationSettings,
) -> Result<SeqDb, DatabaseOpenError> {
    let registry = match &settings.blastdb {
        Ok(value) => value.as_deref(),
        Err(reason) => return Err(DatabaseOpenError::Registry(reason.clone())),
    };
    let path = SearchPath::from_environment(registry).map_err(DatabaseOpenError::Reader)?;
    SeqDb::open_in(name, SeqType::Protein, &path).map_err(DatabaseOpenError::Reader)
}

/// Why `open_search_database` could not open a database.
#[derive(Debug)]
pub(super) enum DatabaseOpenError {
    /// The registry entry `[BLAST] BLASTDB` has a value that LOSAT does not read as NCBI
    /// does (`ApplicationSettings::blastdb`).
    Registry(String),
    /// An error of the reader (NCBI's own, or LOSAT's refusal: `database_error`).
    Reader(SeqDbError),
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203
// ```c++
//     catch (const CSeqDBException& e) {                                      \
// 		if (e.GetErrCode() == CSeqDBException::eTooManyOpenFiles) {         \
// 			string err_msg =                                                \
//                 "Too many open files, please raise the open file limit";    \
//         	LOG_POST(Error << "BLAST Database error: " << err_msg);         \
//             exit_code = BLAST_DATABASE_ERROR;                               \
// 		}                                                                   \
// 		else if (e.GetErrCode() == CSeqDBException::eMemoryMappingFailure){ \
// 			string err_msg = "Database memory map file error";              \
//             if (!e.GetMsg().empty())                                        \
//                 err_msg += " - " + e.GetMsg();                              \
//         	LOG_POST(Error << "BLAST Database error: " << err_msg);         \
//             exit_code = BLAST_DATABASE_ERROR;                               \
// 		}                                                                   \
//         else {                                                              \
//         	LOG_POST(Error << "BLAST Database error: " << e.GetMsg());      \
//         	exit_code = BLAST_DATABASE_ERROR;                               \
//         }                                                                   \
//     }                                                                       \
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:153
// ```c++
// #define BLAST_DATABASE_ERROR                    2
// ```
/// A reader error outside the preliminary search (when the database is opened, and when
/// the sequences and titles of the hit lists are read) as the command line reports it:
/// NCBI's own errors (`SeqDbError::ncbi_message`: a database that is not found, a file
/// that is no database file or cannot be opened) as `BLAST Database error: <text>` with
/// exit status 2; everything else is LOSAT's refusal, in LOSAT's wording.
pub(super) fn database_error(error: SeqDbError) -> anyhow::Error {
    match error.ncbi_message() {
        Some(text) => NativeError {
            exit: 2,
            message: format!("BLAST Database error: {text}\n"),
        }
        .into(),
        None => refusal(&error.to_string()),
    }
}

/// LOSAT's refusal of a database that it does not search as NCBI does.
fn refusal(reason: &str) -> anyhow::Error {
    anyhow::anyhow!("{reason}; this database is not supported by LOSAT's BLASTP")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251
// ```c++
// 	catch (CSeqDBException & e) {
// 		if (e.GetErrCode() == CSeqDBException::eTooManyOpenFiles) {
//             _TRACE(e.GetMsg());
// 			string err_msg = BlastErrorCode2String(BLASTERR_DB_TOO_MANY_OPEN_FILES);
// 			NCBI_THROW(CBlastException, eCoreBlastError, err_msg);
// 		}
// 		else {
//             _TRACE(e.GetMsg());
// 			string err_msg = BlastErrorCode2String(BLASTERR_DB_MEMORY_MAP);
// 			NCBI_THROW(CBlastException, eCoreBlastError, err_msg);
// 		}
// 	}
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
// ```c++
//         } else {                                                            \
//             LOG_POST(Error << "BLAST engine error: " << e.GetMsg());        \
//             exit_code = BLAST_ENGINE_ERROR;                                 \
// ```
/// A reader error while the preliminary search reads the sequences: NCBI's search stage
/// replaces the reader's text (`SeqDbError::ncbi_search_message`), and the application
/// reports `BLAST engine error: <text>` with exit status 3. The error is kept apart from
/// the others (`PreliminaryReadError`), because NCBI has written the outfmt 0 prolog and
/// read the first query batch before it (`blast_engine::run_resolved_in_pool`).
fn search_read_error(error: SeqDbError) -> anyhow::Error {
    match error.ncbi_search_message() {
        Some(text) => PreliminaryReadError(text).into(),
        None => refusal(&error.to_string()),
    }
}

/// NCBI's error of a database file that cannot be read during the preliminary search
/// (`search_read_error`), with the text after `BLAST engine error: `.
#[derive(Debug)]
pub(super) struct PreliminaryReadError(pub(super) &'static str);

impl std::fmt::Display for PreliminaryReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BLAST engine error: {}", self.0)
    }
}

impl std::error::Error for PreliminaryReadError {}

impl BlastpDatabase {
    /// Opens the database of `-db name` (`open_search_database`) with the command line's
    /// errors (`database_error`).
    pub(super) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        let db = open_search_database(name, settings).map_err(|error| match error {
            DatabaseOpenError::Registry(reason) => refusal(&reason),
            DatabaseOpenError::Reader(error) => database_error(error),
        })?;
        Ok(Self {
            db,
            name: name.to_string(),
        })
    }

    /// The `-db` value as NCBI keeps it.
    pub(super) fn name(&self) -> &str {
        &self.name
    }

    /// The number of sequences (`CSeqDB::GetNumSeqs`).
    pub(super) fn num_sequences(&self) -> u32 {
        self.db.num_seqs()
    }

    /// The total length (`CSeqDB::GetTotalLength`).
    pub(super) fn total_letters(&self) -> u64 {
        self.db.total_length()
    }

    /// The total length as the engine holds it.
    pub(super) fn total_letters_usize(&self) -> Result<usize> {
        let total = self.db.total_length();
        usize::try_from(total)
            .ok()
            .filter(|&letters| i64::try_from(letters).is_ok())
            .with_context(|| {
                format!(
                    "a database of {total} letters is not supported by LOSAT's BLASTP on this platform"
                )
            })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
    // ```c++
    //     info.name = seqdb->GetDBNameList();
    //     info.definition = seqdb->GetTitle();
    //     if (info.definition.empty())
    //         info.definition = info.name;
    //     info.date = seqdb->GetDate();
    //     info.total_length = seqdb->GetTotalLength();
    //     info.number_seqs = seqdb->GetNumSeqs();
    // ```
    /// What the `Database:` lines of the pairwise report show: the title (the name as
    /// given for an empty title), the date and the real totals, also with `-dbsize`.
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        let title = self.db.title();
        PairwiseDatabase {
            definition: if title.is_empty() {
                self.name.as_bytes().to_vec()
            } else {
                title.to_vec()
            },
            date: self.db.date().to_string(),
            num_sequences: u64::from(self.db.num_seqs()),
            total_letters: self.db.total_length(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-848
// ```c
//          BLAST_ComputeLengthAdjustment(kbp->K, kbp->logK,
//                                        alpha/kbp->Lambda, beta,
//                                        query_length, db_length,
//                                        db_num_seqs, &length_adjustment);
//
//          if (effective_search_space == 0) {
//
//              /* if the database length was specified, do not
//                 adjust it when calculating the search space;
//                 it's counter-intuitive to specify a value and
//                 not have that value be used */
//         	 /* Changing this rule for now sicne cutoff score depends
//         	  * on the effective seach space length. SB-902
//         	  */
//
//         	 Int8 effective_db_length = db_length - ((Int8)db_num_seqs * length_adjustment);
//
//         	 // Just in case effective_db_length < 0
//         	 if (effective_db_length <= 0)
//         		 effective_db_length = 1;
//
//              effective_search_space = effective_db_length *
//                              (query_length - length_adjustment);
//          }
//       }
//       query_info->contexts[index].eff_searchsp = effective_search_space;
//       query_info->contexts[index].length_adjustment = length_adjustment;
// ```
/// The length adjustment and the effective search space of a query of `query_len`
/// letters against a database of `db_length` letters (the real total, or `-dbsize`) and
/// `db_num_seqs` sequences, in NCBI's `Int8`: the search space as the reports print it,
/// and as the statistics read it (a `double`). A product that does not fit `Int8`, where
/// NCBI's arithmetic overflows, is rejected.
pub(super) fn blastp_database_search_space(
    query_len: usize,
    db_length: i64,
    db_num_seqs: i64,
    params: &KarlinParams,
) -> Result<(SearchSpace, i64)> {
    let overflow = || {
        format!(
            "a database length of {db_length} letters with a query of {query_len} letters overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's BLASTP"
        )
    };
    let query_length = i64::try_from(query_len).ok().with_context(overflow)?;
    let length_adjustment =
        compute_length_adjustment_ncbi(query_length, db_length, db_num_seqs, params)
            .length_adjustment;
    let mut effective_db_length = db_num_seqs
        .checked_mul(length_adjustment)
        .and_then(|adjusted| db_length.checked_sub(adjusted))
        .with_context(overflow)?;
    if effective_db_length <= 0 {
        effective_db_length = 1;
    }
    let effective_query_length = query_length - length_adjustment;
    let effective_search_space = effective_db_length
        .checked_mul(effective_query_length)
        .with_context(overflow)?;
    Ok((
        SearchSpace {
            effective_query_len: effective_query_length as f64,
            effective_db_len: effective_db_length as f64,
            effective_space: effective_search_space as f64,
            length_adjustment,
        },
        effective_search_space,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593
// ```c
//     /* Initialize with BLAST_SCORE_MIN */
//     for (i = 0; i < sbp->alphabet_size; i++) {
//         for (j = 0; j < sbp->alphabet_size; j++) {
//             matrix[i][j] = BLAST_SCORE_MIN;
//         }
//     }
//
//     for (i = 0; i < sbp->alphabet_size; i++) {
//         for (j = 0; j < sbp->alphabet_size; j++) {
//             /* skip special characters */
//             if (i == AMINOACID_TO_NCBISTDAA['U'] ||
//                 i == AMINOACID_TO_NCBISTDAA['O'] ||
//                 i == AMINOACID_TO_NCBISTDAA['-'] ||
//                 j == AMINOACID_TO_NCBISTDAA['U'] ||
//                 j == AMINOACID_TO_NCBISTDAA['O'] ||
//                 j == AMINOACID_TO_NCBISTDAA['-']) {
//                 continue;
//             }
//             matrix[i][j] = NCBISM_GetScore((const SNCBIPackedScoreMatrix*) psm,
//                                            i, j);
//         }
//     }
//     /* Use the C scores for U and X scores for the O characters;
//        if this is not done then they will never align to non-gap residues */
//     x_index = AMINOACID_TO_NCBISTDAA['X'];
//     u_index = AMINOACID_TO_NCBISTDAA['U'];
//     o_index = AMINOACID_TO_NCBISTDAA['O'];
//     c_index = AMINOACID_TO_NCBISTDAA['C'];
//     for (i = 0; i < sbp->alphabet_size; i++) {
//         matrix[u_index][i] = matrix[c_index][i];
//         matrix[i][u_index] = matrix[i][c_index];
//         matrix[o_index][i] = matrix[x_index][i];
//         matrix[i][o_index] = matrix[i][x_index];
//     }
// ```
// NCBI searches the stored bytes in place: a stored O (26) scores as X in every matrix
// row and column, so the preliminary stage, which reads scores only, reads it as X (21),
// as the `-subject` encoder does (`encode_protein_sequence`, after
// blast_setup_cxx.cpp:894-899); the redo and the reports read the stored letter
// (`ListedSubjects`). A stored gap character (0, which
// makeblastdb writes for `-` in a protein FASTA) never aligns in NCBI (row and column 0
// stay `BLAST_SCORE_MIN`) and NCBI's formatter then counts it as a gap of the alignment;
// LOSAT's kernels read byte 0 as the end of a sequence, so a sequence that holds one is
// rejected. A byte that is no NCBIstdaa residue (28 or more) indexes NCBI's matrix past
// its rows and is rejected too.
/// Prepares the block `block` (the sequences from the ordinal id `range.first_ordinal`)
/// for the search: O is read as X, and a sequence that LOSAT cannot search as NCBI does
/// is an error that names its ordinal id.
fn prepare_search_block(block: &mut SequenceBlock) -> Result<()> {
    let sequences = block.len();
    let mut gaps = 0usize;
    let mut invalid = false;
    for residue in block.protein_bytes_mut() {
        if *residue == NCBISTDAA_O {
            *residue = ncbistdaa::X;
        }
        gaps += usize::from(*residue == NCBISTDAA_GAP);
        invalid |= usize::from(*residue) >= BLASTAA_SIZE;
    }
    // A well-formed block holds one sentinel byte before its first sequence and one
    // after every sequence, and no other byte 0.
    if gaps == sequences + 1 && !invalid {
        return Ok(());
    }
    for i in 0..sequences {
        let ordinal = u64::from(block.first_ordinal()) + i as u64;
        let with_sentinels = block.protein_with_sentinels(i);
        let (first, last) = (with_sentinels[0], with_sentinels[with_sentinels.len() - 1]);
        if first != NCBISTDAA_GAP || last != NCBISTDAA_GAP {
            bail!(
                "database sequence {ordinal} is not between two NUL bytes in its sequence file (the volume is damaged), which is not supported by LOSAT's BLASTP"
            );
        }
        let residues = block.protein(i);
        if let Some(position) = residues
            .iter()
            .position(|&residue| residue == NCBISTDAA_GAP)
        {
            bail!(
                "database sequence {ordinal} holds a gap character ('-') at position {}; a protein database sequence with gap characters is not supported by LOSAT's BLASTP (NCBI BLAST+ never aligns the character and then counts it as a gap of the alignment)",
                position + 1
            );
        }
        if let Some(position) = residues
            .iter()
            .position(|&residue| usize::from(residue) >= BLASTAA_SIZE)
        {
            bail!(
                "database sequence {ordinal} holds a byte that is no NCBIstdaa residue at position {}, which is not supported by LOSAT's BLASTP",
                position + 1
            );
        }
    }
    bail!("a database block is not well-formed, which is not supported by LOSAT's BLASTP")
}

/// Reads the block `range` and prepares it for the search (`prepare_search_block`).
fn read_search_block(db: &SeqDb, range: BlockRange) -> Result<SequenceBlock> {
    let mut block = db.read_block(range).map_err(search_read_error)?;
    prepare_search_block(&mut block)?;
    Ok(block)
}

/// Sequence `i` of a prepared block as the kernel reads it.
fn block_subject(block: &SequenceBlock, i: usize) -> BlastpSubjectRef<'_> {
    let aa_seq = block.protein_with_sentinels(i);
    BlastpSubjectRef {
        aa_seq,
        aa_len: aa_seq.len() - 2,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:160-172
// ```c
// Int2
// Blast_ExtendWordExit(Blast_ExtendWord * ewp, Int4 subject_length)
// {
//     if (!ewp)
//         return -1;
//
//     if (ewp->diag_table) {
//         if (ewp->diag_table->offset >= INT4_MAX / 4) {
//             ewp->diag_table->offset = ewp->diag_table->window;
//             s_BlastDiagClear(ewp->diag_table);
//         } else {
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
//         }
// ```
/// The offset of the diagonal table after the sequences of `range` when it is
/// `diag_offset` before them: the recurrence of `BlastpSubjectScratch::finish_subject`
/// over the sequence lengths of the index. A sequence without letters is not searched
/// (blast_engine.c:1429-1432), so the offset does not move for it.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn diag_offset_after(db: &SeqDb, range: BlockRange, mut diag_offset: i32, window: i32) -> i32 {
    let offsets = db.volumes()[range.volume].index().sequence_offsets();
    for bounds in offsets[range.first as usize..=range.end as usize].windows(2) {
        // The sequence and the sentinel byte after it.
        let subject_length = (bounds[1] - bounds[0]).saturating_sub(1);
        if subject_length == 0 {
            continue;
        }
        if diag_offset >= i32::MAX / 4 {
            diag_offset = window;
        } else {
            diag_offset += subject_length as i32 + window;
        }
    }
    diag_offset
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
/// ordinal id, through the per-subject kernel `process_subject`, with its HSP lists put
/// into `preliminary_hit_lists` in that order (NCBI's loop with one thread).
///
/// One thread: the blocks in order, one scratch (one diagonal table) for the whole
/// database, the hit lists updated after every subject, as the one-thread loop of the
/// `-subject` path.
///
/// More threads: rounds of blocks. The blocks of a round are read and searched in
/// parallel; each worker starts its block from a cleared diagonal table at the offset
/// that the one-thread order gives the block's first sequence
/// (`prepare_independent_subject`, `diag_offset_after`) and advances it inside the block
/// as the one-thread path does. A table entry left by an earlier subject is below the
/// offset of every later subject by at least the window, so it reads as a cleared entry
/// does (the property the parallel `-subject` path uses per subject). The round's
/// results are put into the hit lists in ordinal order, and its blocks are dropped
/// before the next round is read: the memory held is bounded by the threads and the
/// block size, not by the database.
pub(super) fn search_blocks<S, K>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    new_scratch: &S,
    process_subject: &K,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    prelim_hitlist_size: usize,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
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
        window,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        process_subject,
        preliminary_hit_lists,
        prelim_hitlist_size,
    )
}

/// `search_blocks` with blocks of at most `block_bytes` and rounds of `blocks_per_thread`
/// blocks per thread.
#[allow(clippy::too_many_arguments)]
fn search_blocks_of<S, K>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    prelim_hitlist_size: usize,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
{
    let db = &database.db;
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if super::blastp_should_parallelize_stage(
        "prelim",
        pool.threads(),
        pool.threads(),
        db.num_seqs() as usize,
        db.num_seqs() as usize,
        0,
    ) && pool.enabled()
    {
        let round_blocks = pool.threads().saturating_mul(blocks_per_thread).max(1);
        let mut next_ordinal = 0u32;
        let mut diag_offset = window;
        loop {
            // The blocks of the round and the diagonal offset of each one's first
            // sequence in the one-thread order.
            let mut round: Vec<(BlockRange, i32)> = Vec::with_capacity(round_blocks);
            while round.len() < round_blocks {
                let Some(range) = db.next_block(next_ordinal, block_bytes) else {
                    break;
                };
                round.push((range, diag_offset));
                diag_offset = diag_offset_after(db, range, diag_offset, window);
                next_ordinal = range.end_ordinal();
            }
            if round.is_empty() {
                break;
            }
            let found: Vec<Result<Vec<BlastpSubjectResult>>> = pool.install(|| {
                round
                    .par_iter()
                    .map_init(
                        new_scratch,
                        |scratch,
                         &(range, first_diag_offset)|
                         -> Result<Vec<BlastpSubjectResult>> {
                            let block = read_search_block(db, range)?;
                            scratch.prepare_independent_subject(first_diag_offset);
                            let mut with_hits = Vec::new();
                            let mut subject_result = BlastpSubjectResult::new(0);
                            for i in 0..block.len() {
                                let s_idx = range.first_ordinal as usize + i;
                                process_subject(
                                    scratch,
                                    s_idx,
                                    block_subject(&block, i),
                                    &mut subject_result,
                                );
                                if !subject_result.query_hsp_lists.is_empty() {
                                    with_hits.push(std::mem::replace(
                                        &mut subject_result,
                                        BlastpSubjectResult::new(0),
                                    ));
                                }
                            }
                            Ok(with_hits)
                        },
                    )
                    .collect()
            });
            // The subjects without HSP lists change no hit list, so the order of the
            // updates is the one-thread order.
            for block_results in found {
                for mut subject_result in block_results? {
                    update_blastp_preliminary_hit_lists(
                        preliminary_hit_lists,
                        &mut subject_result,
                        prelim_hitlist_size,
                    );
                }
            }
        }
        return Ok(());
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = (pool, window, blocks_per_thread);

    let mut scratch = new_scratch();
    let mut subject_result = BlastpSubjectResult::new(0);
    let mut next_ordinal = 0u32;
    while let Some(range) = db.next_block(next_ordinal, block_bytes) {
        let block = read_search_block(db, range)?;
        for i in 0..block.len() {
            process_subject(
                &mut scratch,
                range.first_ordinal as usize + i,
                block_subject(&block, i),
                &mut subject_result,
            );
            update_blastp_preliminary_hit_lists(
                preliminary_hit_lists,
                &mut subject_result,
                prelim_hitlist_size,
            );
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

/// The database sequences that the preliminary hit lists name, read by ordinal id after
/// the preliminary stage: they stand for the subject set in the rest of the search (the
/// composition redo and the reports).
///
/// From here on a hit list names a sequence by its rank among the listed ordinal ids
/// (`ordinals[rank]`), so that the code of the `-subject` path indexes these vectors as
/// it indexes the encoded subject set. The rank grows with the ordinal id, so every
/// comparison that breaks a tie by the ordinal id (`s_EvalueCompareHSPLists`,
/// `ScoreCompareHSPs`; `hsp.rs`) gives the same order.
pub(super) struct ListedSubjects {
    /// The listed ordinal ids, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The sequences as stored (O stays O), between two sentinel bytes, by rank: what
    /// the composition redo reads.
    pub(super) encoded: Vec<EncodedProtein>,
    /// The stored titles and the stored letters as the reports read them, by rank.
    pub(super) records: Vec<FastaRecord>,
    /// Whether the stored title of a subject holds a byte outside ' ' ..= '~' (read as
    /// '#'), by rank (`check_reported_titles`).
    replaced_title_bytes: Vec<bool>,
}

impl ListedSubjects {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1381-1385
    // ```c
    //         if( program_number == eBlastTypeTblastn ) {
    //             seq_info->seq_arg.encoding = eBlastEncodingNcbi4na;
    //         } else {
    //             seq_info->seq_arg.encoding = eBlastEncodingProtein;
    //         }
    // ```
    // NCBI's redo fetches each sequence of the hit lists again from the database, in the
    // encoding of the search: the stored bytes.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1620-1636
    // ```c
    //     /* Copy the sequence data */
    //     for (idx = 0;  idx < seqData->length;  idx++) {
    //         seqData->data[idx] = origData[idx];
    //     }
    //
    //     if ( !(KAPPA_BLASTP_NO_SEG_SEQUENCE) ) {
    //         if (compo_adjust_mode
    //             && (!subject_maybe_biased || *subject_maybe_biased)) {
    //
    //             if ( (!shouldTestIdentical)
    //                  || (shouldTestIdentical
    //                      && (!s_TestNearIdentical(seqData, 0, queryData,
    //                                               q_range->begin, query_words,
    //                                               align)))) {
    //
    //                 status = s_DoSegSequenceData(seqData, eBlastTypeBlastp,
    //                                              subject_maybe_biased);
    // ```
    // Here the stored letter matters: a stored O (26) is not the X (21) that the query
    // holds for its own O (blast_setup_cxx.cpp:894-899), so an O aligned to an O counts
    // against "near identical" and the subject is then SEG-filtered, where the same
    // records given as `-subject` (O read as X on both sides) are near identical. Shown
    // with the oracle (campaign dbsearch-20261012, `$BUILD_ROOT/dbsearch/db3/exp_o_near`;
    // fixture rows `db.o_stored.*`): NCBI `-db` gives 266 bits where `-subject` gives
    // 307. So the redo reads the stored bytes; its matrices score O as X
    // (`adjust_scores.rs`, composition_adjustment.c:681-703). The preliminary stage uses
    // scores only, which are the same for O and X (`prepare_search_block`).
    //
    // NCBI reference for the titles:
    // c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:970-1017 (`s_GetBioseqTitle`: the
    // title of the first definition line, the others appended as ` >id title`).
    /// Reads the sequences and the definition lines of the ordinal ids in `hit_lists`
    /// and renames the subjects of the hit lists to their ranks.
    pub(super) fn read(
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
    ) -> Result<Self> {
        let db = &database.db;
        let mut ordinals: Vec<u32> = hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| hit_list.hsplist_array.iter().map(|hsp_list| hsp_list.oid))
            .collect();
        ordinals.sort_unstable();
        ordinals.dedup();
        let mut encoded = Vec::with_capacity(ordinals.len());
        let mut records = Vec::with_capacity(ordinals.len());
        let mut replaced_title_bytes = Vec::with_capacity(ordinals.len());
        for &ordinal in &ordinals {
            let stored = db.read_sequence(ordinal).map_err(database_error)?;
            let residues = stored.protein(0);
            // (A sequence with a gap residue or a foreign byte never reaches a hit
            // list: `prepare_search_block` rejects its block.)
            let mut aa_seq = Vec::with_capacity(residues.len() + 2);
            aa_seq.push(NCBISTDAA_GAP);
            aa_seq.extend_from_slice(residues);
            aa_seq.push(NCBISTDAA_GAP);
            encoded.push(EncodedProtein {
                aa_seq,
                aa_len: residues.len(),
            });
            let mut deflines = db.deflines(ordinal).map_err(database_error)?;
            // makeblastdb without -parse_seqids writes one definition line per sequence;
            // NCBI shows the others of a sequence as ` >id title` after the first.
            if deflines.len() != 1 {
                bail!(
                    "database sequence {ordinal} has {} definition lines; a sequence with more or fewer than one is not supported by LOSAT's BLASTP",
                    deflines.len()
                );
            }
            let defline = deflines.remove(0);
            replaced_title_bytes.push(defline.replaced_title_bytes > 0);
            let letters: Vec<u8> = residues
                .iter()
                .map(|&residue| ncbistdaa_to_ascii(residue) as u8)
                .collect();
            records.push(FastaRecord::new(
                format!("gnl|BL_ORD_ID|{}", defline.bl_ord_id),
                &defline.title,
                &letters,
            ));
        }
        let rank = |ordinal: u32| -> u32 {
            ordinals
                .binary_search(&ordinal)
                .expect("every subject of the hit lists was listed") as u32
        };
        for hit_list in hit_lists.iter_mut().flatten() {
            for hsp_list in &mut hit_list.hsplist_array {
                hsp_list.oid = rank(hsp_list.oid);
                for hsp in &mut hsp_list.hsps {
                    hsp.s_idx = rank(hsp.s_idx);
                }
            }
        }
        Ok(Self {
            ordinals,
            encoded,
            records,
            replaced_title_bytes,
        })
    }

    // NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1865-1886
    // ```c++
    //         string message;
    //         if (io != NULL) {
    //             message += io->GetStackTrace() + "\n";
    //         }
    //         message += "Bad char [0x" +
    //                    NStr::NumericToString((unsigned char)c,0,16)+
    //                    "] in string";
    //         if (io != NULL) {
    //             message += " at " + io->GetPosition();
    //         }
    //         if (!str.empty()) {
    //             message += "\n" + str;
    //         }
    //         switch (fix_method) {
    //         case eFNP_ReplaceAndWarn:
    // #if 0
    //             CNcbiDiag(eDiag_Error, eDPF_Default)
    //                 << ErrCode(NCBI_ERRCODE_X, 7) << message << Endm;
    // #else
    //             ERR_POST_X(7, message);
    // #endif
    //             break;
    // ```
    // NCBI BLAST+ 2.17.0 `blastp -db` on a version 4 database whose header file holds
    // the bytes 0xE9 and 0x7F in two titles (campaign dbsearch-20261012,
    // `$BUILD_ROOT/dbsearch/db3/exp_badbyte`): the report shows `#` for the byte, and
    // standard error gets, for each reported subject and bad byte,
    //   Error: (803.7) [blastp] Blast-def-line-set.E.title
    //   Bad char [0xE9] in string at byte 21
    //   <the stored title>
    // as often as the formatter reads the definition line: twice per subject for the
    // default outfmt 6, three times for `6 sseqid stitle`, five times for outfmt 0. The
    // number of readings is the formatter's, not ported; makeblastdb and `LOSAT makedb`
    // store `#` for such a byte, so only a header file written otherwise holds one.
    /// LOSAT's rejection of a reported subject (one with rows in the final `hit_lists`)
    /// whose stored title holds a byte outside ' ' ..= '~': NCBI writes an error message
    /// for the byte each time its formatter reads the definition line.
    pub(super) fn check_reported_titles(&self, hit_lists: &[Option<BlastpHitList>]) -> Result<()> {
        for rank in Self::reported(hit_lists) {
            if self.replaced_title_bytes[rank] {
                bail!(
                    "database sequence {} has a stored title with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line); this is not supported by LOSAT's BLASTP",
                    self.ordinals[rank]
                );
            }
        }
        Ok(())
    }

    /// The ranks of the subjects with rows in `hit_lists` (the final hit lists).
    fn reported(hit_lists: &[Option<BlastpHitList>]) -> std::collections::BTreeSet<usize> {
        hit_lists
            .iter()
            .flatten()
            .flat_map(|hit_list| {
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                hit_list.hsplist_array[..count]
                    .iter()
                    .flat_map(|hsp_list| hsp_list.hsps.iter().map(|hsp| hsp.s_idx as usize))
            })
            .collect()
    }

    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every
    /// listed subject, by rank (`report::defline::database_subject_id`); a subject
    /// without rows gets an empty id. A subject with rows whose id LOSAT cannot write as
    /// NCBI is an error before any format is written.
    pub(super) fn tabular_subject_ids(
        &self,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<Vec<Arc<[u8]>>> {
        let reported = Self::reported(hit_lists);
        let mut ids: Vec<Arc<[u8]>> = Vec::with_capacity(self.records.len());
        for (rank, record) in self.records.iter().enumerate() {
            if !reported.contains(&rank) {
                ids.push(Arc::from(&b""[..]));
                continue;
            }
            let ordinal = self.ordinals[rank];
            let id = crate::report::defline::database_subject_id(&record.title, true)
                .map_err(|problem| {
                    anyhow::anyhow!(
                        "database sequence {ordinal} has a title {problem}, which is not supported by LOSAT's BLASTP"
                    )
                })?;
            ids.push(Arc::from(id));
        }
        Ok(ids)
    }

    /// The `stitle` text of every listed subject, by rank: the stored title (an empty
    /// one is written `N/A`, `blastp_subject_stitle`).
    pub(super) fn stitles(&self) -> Vec<Arc<[u8]>> {
        self.records
            .iter()
            .map(|record| Arc::from(record.title()))
            .collect()
    }

    /// LOSAT's rejection of the outfmt 0 titles that it does not write as NCBI, of the
    /// subjects that the report shows (`shown`, by rank): a title of punctuation that
    /// NCBI's `x_CleanAndCompress` reads past the end of when it makes the title of the
    /// description table (NCBI BLAST+ 2.17.0 crashes there), as for the records of a
    /// subject set (`check_shown_subject_titles`).
    pub(super) fn check_shown_titles(
        &self,
        shown: &std::collections::BTreeSet<usize>,
    ) -> Result<()> {
        for &rank in shown {
            let title = &self.records[rank].title;
            let reads_past_end = crate::report::defline::generate_defline(title, true, true)
                .is_ok_and(|defline| defline.reads_past_end);
            if reads_past_end {
                bail!(
                    "database sequence {} has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0); this is not supported by LOSAT's BLASTP",
                    self.ordinals[rank]
                );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::blastp::hsp::{BlastpHsp, BlastpHspList};
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    use crate::stats::lookup_protein_params;

    fn fixture(name: &str) -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/seqdb/db")
            .join(name)
            .to_str()
            .expect("fixture path is UTF-8")
            .to_string()
    }

    fn open(name: &str) -> BlastpDatabase {
        BlastpDatabase::open(&fixture(name), &ApplicationSettings::default())
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    fn blosum62_gapped() -> KarlinParams {
        lookup_protein_params(&ProteinScoringSpec {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        })
    }

    // NCBI BLAST+ 2.17.0, `blastp -db sp_1k -query q1.faa -num_threads 1` (1,000 sequences,
    // 357,308 letters; a query of 379 letters): `Effective search space used: 84653016`,
    // and with `-dbsize 1000000000`: `245967282000` (the task's
    // logs/db3_probe_cli1.log, cases dbsize.db.fmt0_none and dbsize.db.fmt0).
    #[test]
    fn search_space_is_the_one_of_ncbis_footer() {
        let params = blosum62_gapped();
        let (space, exact) = blastp_database_search_space(379, 357_308, 1_000, &params).unwrap();
        assert_eq!(exact, 84_653_016);
        assert_eq!(space.effective_space, 84_653_016.0);
        let (space, exact) =
            blastp_database_search_space(379, 1_000_000_000, 1_000, &params).unwrap();
        assert_eq!(exact, 245_967_282_000);
        assert_eq!(space.effective_space, 245_967_282_000.0);
        assert_eq!(
            space.effective_query_len * space.effective_db_len,
            space.effective_space
        );
    }

    // blast_setup.c:821-844: the adjustment is taken off the query once and off the
    // database once per sequence.
    #[test]
    fn search_space_of_a_small_database_length() {
        let params = blosum62_gapped();
        // `-dbsize 1`: no length adjustment is possible (blast_stat.c: c < 0), so the
        // search space is the query length.
        let (space, exact) = blastp_database_search_space(379, 1, 1_000, &params).unwrap();
        assert_eq!((space.length_adjustment, exact), (0, 379));
        // A total of five letters per sequence.
        let (space, exact) = blastp_database_search_space(379, 5_000, 1_000, &params).unwrap();
        assert_eq!(space.length_adjustment, 4);
        assert_eq!(space.effective_db_len, 1_000.0);
        assert_eq!(exact, 1_000 * 375);
    }

    // NCBI's `Int8` product overflows for such a `-dbsize`; LOSAT rejects it.
    #[test]
    fn search_space_beyond_int8_is_rejected() {
        let error =
            blastp_database_search_space(379, 9_000_000_000_000_000_000, 1_000, &blosum62_gapped())
                .unwrap_err();
        assert!(
            format!("{error:#}").contains("overflows NCBI BLAST+'s 64-bit effective search space"),
            "{error:#}"
        );
    }

    #[test]
    fn database_numbers_come_from_the_index_files() {
        // LOSAT/tests/seqdb/expected/prot_mv.info.tsv (blastdbcmd -info).
        let database = open("prot_mv");
        assert_eq!(database.num_sequences(), 10);
        assert_eq!(database.total_letters(), 801);
        assert_eq!(database.total_letters_usize().unwrap(), 801);
        let shown = database.pairwise_database();
        assert_eq!(shown.definition, b"protein in volumes");
        assert_eq!(shown.date, "Oct 12, 2026  4:12 PM");
        assert_eq!((shown.num_sequences, shown.total_letters), (10, 801));
        // The name is the `-db` value as given.
        assert_eq!(database.name(), fixture("prot_mv"));
    }

    // align_format_util.cpp:385-387: an empty title is shown as the name.
    // (A path in the system's temporary directory or in an error text: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn a_database_without_a_title_shows_its_name() {
        let dir = std::env::temp_dir().join(format!(
            "losat-db-notitle-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        // The fixture's index with an empty title: version, type, title length 0, then
        // the rest of the file after the stored title.
        let index = std::fs::read(fixture("prot_v4") + ".pin").unwrap();
        let title_len = u32::from_be_bytes([index[8], index[9], index[10], index[11]]) as usize;
        let mut untitled = index[..8].to_vec();
        untitled.extend_from_slice(&0u32.to_be_bytes());
        untitled.extend_from_slice(&index[12 + title_len..]);
        // The offsets of the index point into the header and sequence files, not into
        // the index itself, so they stay valid.
        std::fs::write(dir.join("u.pin"), untitled).unwrap();
        std::fs::copy(fixture("prot_v4") + ".psq", dir.join("u.psq")).unwrap();
        std::fs::copy(fixture("prot_v4") + ".phr", dir.join("u.phr")).unwrap();
        let name = dir.join("u").to_str().unwrap().to_string();
        let database = BlastpDatabase::open(&name, &ApplicationSettings::default()).unwrap();
        assert_eq!(database.pairwise_database().definition, name.as_bytes());
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // blast_app_util.hpp:185-203 and seqdb.cpp:155-165.
    // (A path in the system's temporary directory or in an error text: Unix spelling.)
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
        let Err(error) = BlastpDatabase::open("", &settings) else {
            panic!("an empty name opened")
        };
        assert_eq!(
            native(error),
            (
                2,
                "BLAST Database error: Database name is required.\n".to_string()
            )
        );
        let Err(error) = BlastpDatabase::open(&fixture("no_such_database"), &settings) else {
            panic!("a missing database opened")
        };
        let (exit, message) = native(error);
        assert_eq!(exit, 2);
        assert!(
            message.starts_with(
                "BLAST Database error: No alias or index file found for protein database ["
            ),
            "{message}"
        );
        // A database with parsed Seq-ids: LOSAT's refusal, not NCBI's error.
        let Err(error) = BlastpDatabase::open(&fixture("prot_ids_v4"), &settings) else {
            panic!("a database with parsed Seq-ids opened")
        };
        assert!(error.downcast_ref::<NativeError>().is_none());
        assert!(
            format!("{error:#}").ends_with("this database is not supported by LOSAT's BLASTP"),
            "{error:#}"
        );
    }

    // blast_stat.c:1583-1593: O scores as X. p03 of the fixture holds every NCBIstdaa
    // letter but the gap (`ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT`).
    #[test]
    fn a_search_block_reads_o_as_x_and_keeps_the_other_letters() {
        let database = open("prot_v4");
        let range = database.db.next_block(0, u64::MAX).unwrap();
        assert_eq!((range.first_ordinal, range.len()), (0, 10));
        let stored = database.db.read_block(range).unwrap();
        let block = read_search_block(&database.db, range).unwrap();
        for i in 0..block.len() {
            let expected: Vec<u8> = stored
                .protein(i)
                .iter()
                .map(|&residue| if residue == 26 { 21 } else { residue })
                .collect();
            assert_eq!(block.protein(i), expected.as_slice(), "sequence {i}");
            let subject = block_subject(&block, i);
            assert_eq!(subject.aa_len, expected.len());
            assert_eq!(subject.aa_seq.len(), expected.len() + 2);
            assert_eq!(
                (subject.aa_seq[0], subject.aa_seq[expected.len() + 1]),
                (0, 0)
            );
            assert_eq!(&subject.aa_seq[1..=expected.len()], expected.as_slice());
        }
        let letters: String = block
            .protein(2)
            .iter()
            .map(|&residue| ncbistdaa_to_ascii(residue))
            .collect();
        assert_eq!(letters, "ABCDEFGHIKLMNPQRSTVWXYZUXJMHCENQT");
        assert!(stored.protein(2).contains(&26) && !block.protein(2).contains(&26));
    }

    #[test]
    fn a_search_block_with_a_gap_residue_or_a_foreign_byte_is_rejected() {
        let database = open("prot_v4");
        let range = database.db.next_block(0, u64::MAX).unwrap();
        // The gap character in sequence 3 (its third residue).
        let mut block = database.db.read_block(range).unwrap();
        let at = block.protein_bytes_mut().len() - block.protein_bytes_mut().len(); // 0: the leading sentinel
        assert_eq!(block.protein_bytes_mut()[at], 0);
        let offset_of = |block: &SequenceBlock, i: usize| -> usize {
            (0..i).map(|j| block.protein(j).len() + 1).sum::<usize>() + 1
        };
        let position = offset_of(&block, 3) + 2;
        block.protein_bytes_mut()[position] = 0;
        let error = prepare_search_block(&mut block).unwrap_err();
        assert!(
            format!("{error:#}")
                .starts_with("database sequence 3 holds a gap character ('-') at position 3;"),
            "{error:#}"
        );
        // A byte that is no NCBIstdaa residue in sequence 9.
        let mut block = database.db.read_block(range).unwrap();
        let position = offset_of(&block, 9);
        block.protein_bytes_mut()[position] = 28;
        let error = prepare_search_block(&mut block).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(
                "database sequence 9 holds a byte that is no NCBIstdaa residue at position 1,"
            ),
            "{error:#}"
        );
        // A sentinel that is not NUL (the byte after sequence 0) with a NUL elsewhere,
        // so that the count of NUL bytes alone would not show it.
        let mut block = database.db.read_block(range).unwrap();
        let after_first = offset_of(&block, 1) - 1;
        block.protein_bytes_mut()[after_first] = 1;
        let error = prepare_search_block(&mut block).unwrap_err();
        assert!(
            format!("{error:#}").contains("is not between two NUL bytes"),
            "{error:#}"
        );
    }

    /// A kernel that records what it is given and "finds" every third sequence.
    fn recording_kernel(
        window: i32,
        seen: &std::sync::Mutex<Vec<(usize, i32, Vec<u8>)>>,
    ) -> impl Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult)
           + Sync
           + '_ {
        move |scratch, s_idx, subject, result| {
            result.reset(s_idx);
            assert_eq!(subject.aa_seq.len(), subject.aa_len + 2);
            seen.lock()
                .unwrap()
                .push((s_idx, scratch.diag_offset, subject.aa_seq.to_vec()));
            scratch.finish_subject(subject.aa_len, window);
            if s_idx % 3 == 0 {
                result.query_hsp_lists.push((
                    0,
                    BlastpHspList {
                        oid: s_idx as u32,
                        query_index: 0,
                        hsps: vec![test_hsp(s_idx as u32, 50 + (s_idx as i32 % 2))],
                        best_evalue: 1e-5,
                    },
                ));
            }
        }
    }

    fn test_hsp(s_idx: u32, raw_score: i32) -> BlastpHsp {
        BlastpHsp {
            identity: 100.0,
            length: 10,
            mismatch: 0,
            gapopen: 0,
            q_start: 1,
            q_end: 10,
            s_start: 1,
            s_end: 10,
            gapped_q_start: 0,
            gapped_s_start: 0,
            e_value: 1e-5,
            bit_score: 30.0,
            num_ident: 10,
            query_context: 0,
            query_frame: 0,
            subject_frame: 0,
            query_length: 20,
            q_idx: 0,
            s_idx,
            raw_score,
            gap_info: None,
            num_positives: 10,
            comp_adjustment_method: 0,
        }
    }

    fn run_blocks(
        database: &BlastpDatabase,
        threads: usize,
        block_bytes: u64,
        window: i32,
    ) -> (Vec<(usize, i32, Vec<u8>)>, Vec<u32>) {
        let seen = std::sync::Mutex::new(Vec::new());
        let mut hit_lists: Vec<Option<BlastpHitList>> = vec![None];
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            search_blocks_of(
                database,
                pool,
                window,
                block_bytes,
                2,
                &|| BlastpSubjectScratch::new(16, 64, window, 21),
                &recording_kernel(window, &seen),
                &mut hit_lists,
                500,
            )
        })
        .unwrap();
        let mut seen = seen.into_inner().unwrap();
        seen.sort_by_key(|entry| entry.0);
        let listed = hit_lists[0]
            .as_ref()
            .map(|list| {
                list.hsplist_array[..list.hsplist_count]
                    .iter()
                    .map(|l| l.oid)
                    .collect()
            })
            .unwrap_or_default();
        (seen, listed)
    }

    // blast_extend.c:160-172 (the diagonal offset moves by the subject length and the
    // window after every subject) and blast_engine.c:1411-1554 (the subjects in the order
    // of their ordinal ids): the block driver gives the kernel every sequence once, with
    // the diagonal offset of the one-thread order, whatever the block size and the thread
    // count, across the six volumes of the fixture.
    #[test]
    fn blocks_give_every_sequence_the_diagonal_offset_of_the_one_thread_order() {
        let database = open("prot_mv");
        let single = open("prot_v4");
        let window = 40;
        // LOSAT/tests/seqdb/expected/prot_v4.entries.tsv: the lengths in ordinal order.
        let lengths = [61, 1, 33, 46, 130, 88, 45, 46, 51, 300];
        let mut expected_offsets = Vec::new();
        let mut offset = window;
        for length in lengths {
            expected_offsets.push(offset);
            offset += length + window;
        }
        let (reference, reference_listed) = run_blocks(&database, 1, u64::MAX, window);
        assert_eq!(reference.len(), 10);
        for (ordinal, (s_idx, diag_offset, aa_seq)) in reference.iter().enumerate() {
            assert_eq!(*s_idx, ordinal);
            assert_eq!(
                *diag_offset, expected_offsets[ordinal],
                "sequence {ordinal}"
            );
            assert_eq!(aa_seq.len() as i32, lengths[ordinal] + 2);
        }
        assert_eq!(reference_listed.len(), 4);
        // The same sequences from the single-volume fixture.
        assert_eq!(run_blocks(&single, 1, u64::MAX, window).0, reference);
        for threads in [1usize, 2, 4] {
            for block_bytes in [1u64, 64, 200, 1 << 20] {
                for db in [&database, &single] {
                    let (seen, listed) = run_blocks(db, threads, block_bytes, window);
                    assert_eq!(
                        seen, reference,
                        "{threads} threads, blocks of {block_bytes}"
                    );
                    assert_eq!(
                        listed, reference_listed,
                        "{threads} threads, blocks of {block_bytes}"
                    );
                }
            }
        }
    }

    // blast_extend.c:168-170: at INT4_MAX / 4 the offset returns to the window.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    #[test]
    fn the_offset_recurrence_of_a_block_is_the_one_of_the_scratch() {
        let database = open("prot_v4");
        let range = database.db.next_block(0, u64::MAX).unwrap();
        let lengths = [61usize, 1, 33, 46, 130, 88, 45, 46, 51, 300];
        for (start, window) in [
            (40, 40),
            (0, 0),
            (i32::MAX / 4 - 100, 40),
            (i32::MAX / 4, 11),
        ] {
            let mut scratch = BlastpSubjectScratch::new(16, 64, window, 21);
            scratch.diag_offset = start;
            for length in lengths {
                scratch.finish_subject(length, window);
            }
            assert_eq!(
                diag_offset_after(&database.db, range, start, window),
                scratch.diag_offset,
                "from {start} with window {window}"
            );
        }
    }

    // The listed subjects stand for the subject set after the preliminary stage: ranks
    // in ordinal order, with the stored letters (O is O; p03 and p08 of the fixture).
    #[test]
    fn listed_subjects_are_read_by_ordinal_and_renamed_to_ranks() {
        let database = open("prot_mv");
        let list = |oid: u32, raw_score: i32| BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![test_hsp(oid, raw_score), test_hsp(oid, raw_score - 1)],
            best_evalue: 1e-5,
        };
        let mut first = BlastpHitList::new(10);
        for (oid, score) in [(7, 60), (2, 55), (9, 50)] {
            first.update(list(oid, score));
        }
        let mut second = BlastpHitList::new(10);
        second.update(list(2, 70));
        second.update(list(0, 65));
        let mut hit_lists = vec![Some(first), None, Some(second)];
        let listed = ListedSubjects::read(&database, &mut hit_lists).unwrap();
        assert_eq!(listed.ordinals, [0, 2, 7, 9]);
        let oids = |hit_list: &Option<BlastpHitList>| -> Vec<(u32, Vec<u32>)> {
            let hit_list = hit_list.as_ref().unwrap();
            hit_list.hsplist_array[..hit_list.hsplist_count]
                .iter()
                .map(|l| (l.oid, l.hsps.iter().map(|hsp| hsp.s_idx).collect()))
                .collect()
        };
        let mut renamed = oids(&hit_lists[0]);
        renamed.sort();
        assert_eq!(renamed, [(1, vec![1, 1]), (2, vec![2, 2]), (3, vec![3, 3])]);
        let mut renamed = oids(&hit_lists[2]);
        renamed.sort();
        assert_eq!(renamed, [(0, vec![0, 0]), (1, vec![1, 1])]);
        // p03: every letter; the redo and the report read the stored letters (O is O).
        assert_eq!(
            listed.records[1].title(),
            b"p03 every letter of NCBIstdaa but the gap"
        );
        assert_eq!(
            listed.records[1].seq(),
            b"ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT"
        );
        let redone: String = listed.encoded[1].aa_seq[1..=listed.encoded[1].aa_len]
            .iter()
            .map(|&residue| ncbistdaa_to_ascii(residue))
            .collect();
        assert_eq!(redone, "ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT");
        assert_eq!(listed.encoded[1].aa_seq.len(), 35);
        assert_eq!(
            listed.records[2].seq(),
            b"MUGYDNEDMDNEPCNDCRACIGFRYKYMMVKWOFASVIHWKHVWDJ"
        );
        assert_eq!(listed.encoded[3].aa_len, 300);
        // The tabular ids of the subjects with rows, the stored titles for `stitle`.
        let ids = listed.tabular_subject_ids(&hit_lists).unwrap();
        assert_eq!(
            ids.iter().map(|id| id.as_ref()).collect::<Vec<_>>(),
            [&b"p01"[..], b"p03", b"p08", b"p10"]
        );
        assert_eq!(
            listed.stitles()[0].as_ref(),
            b"p01 first protein [Test organism]"
        );
        listed
            .check_shown_titles(&[0usize, 1, 2, 3].into_iter().collect())
            .unwrap();
        listed.check_reported_titles(&hit_lists).unwrap();
    }

    // NCBI BLAST+ 2.17.0 writes `Error: (803.7) [blastp] Blast-def-line-set.E.title`,
    // `Bad char [0xE9] in string at byte 21` and the stored title to standard error
    // each time its formatter reads such a definition line (`check_reported_titles`).
    // (A path in the system's temporary directory or in an error text: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn a_reported_title_with_a_byte_outside_the_visible_range_is_rejected() {
        let dir = std::env::temp_dir().join(format!(
            "losat-db-badtitle-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        for ext in ["pin", "psq"] {
            std::fs::copy(fixture("prot_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
        }
        // The title of p10 (`p10 last protein`) with a byte of Latin-1 in it.
        let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
        let at = header
            .windows(16)
            .position(|window| window == b"p10 last protein")
            .expect("the title of p10");
        header[at + 5] = 0xE9;
        std::fs::write(dir.join("b.phr"), header).unwrap();
        let name = dir.join("b").to_str().unwrap().to_string();
        let database = BlastpDatabase::open(&name, &ApplicationSettings::default()).unwrap();
        let list = |oid: u32| BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![test_hsp(oid, 60)],
            best_evalue: 1e-5,
        };
        // A subject whose title is clean: nothing to reject.
        let mut clean = BlastpHitList::new(10);
        clean.update(list(0));
        let mut hit_lists = vec![Some(clean)];
        let listed = ListedSubjects::read(&database, &mut hit_lists).unwrap();
        listed.check_reported_titles(&hit_lists).unwrap();
        // The subject with the bad byte: read as `#`, and rejected when it is reported.
        let mut bad = BlastpHitList::new(10);
        bad.update(list(9));
        bad.update(list(0));
        let mut hit_lists = vec![Some(bad)];
        let listed = ListedSubjects::read(&database, &mut hit_lists).unwrap();
        assert_eq!(listed.records[1].title(), b"p10 l#st protein");
        let error = listed.check_reported_titles(&hit_lists).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(
                "database sequence 9 has a stored title with a byte outside ' ' to '~'"
            ),
            "{error:#}"
        );
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
