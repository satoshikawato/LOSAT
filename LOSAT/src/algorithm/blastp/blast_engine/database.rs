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
//! - the sequences of the hit lists are read again by ordinal id after the preliminary
//!   stage (`ListedSubjects`) for the redo, and the definition lines of the final hit
//!   lists after the redo (`ListedSubjects::read_titles`) for the reports;
//! - the reports show the database (`BlastpDatabase::pairwise_database`, `name`) and
//!   the subject ids and titles of sequences without parsed Seq-ids
//!   (`ListedTitles::tabular_subject_ids`, `stitles`).
//!
//! The per-subject kernel, the composition redo and the report writers are those of the
//! `-subject` path. The opening of the database, the reader's errors and the rules for the
//! stored titles are shared with BLASTN (`algorithm::common::database`).

use std::sync::Arc;

use anyhow::{bail, Context, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

use crate::algorithm::common::database as shared;
pub(super) use crate::algorithm::common::database::{HeaderReadError, PreliminaryReadError};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
#[cfg(test)]
use crate::cli::NativeError;
use crate::report::blast_xml::HitTaxonomy;
use crate::report::db_subject::ParsedSubject;
use crate::report::pairwise::PairwiseDatabase;
use crate::seqdb::{BlockRange, SeqDb, SeqDbError, SeqType, SequenceBlock};
use crate::stats::search_space::SearchSpace;
use crate::stats::{compute_length_adjustment_ncbi, KarlinParams};
use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE};

use super::{
    ncbistdaa_to_ascii, BlastpHitList, BlastpSubjectRef, BlastpSubjectResult, BlastpSubjectScratch,
    EncodedProtein, PrelimSink, XPrelimHsp,
};

// EXPERIMENT (LOSAT_X_BPDBCUT): the blocks of the preliminary stage cut by the database size and
// the thread count.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_db_cut;
// EXPERIMENT (LOSAT_X_BPDBKEEP): one search scratch and one block buffer per worker for the whole
// preliminary stage.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_db_keep;
mod x_listed_lean;
mod x_redo_db_subjects;
#[cfg(test)]
pub(super) use x_listed_lean::x_assert_same_titles_and_letters;
pub(super) use x_listed_lean::x_bp_listed_lean;
pub(super) use x_redo_db_subjects::XShownResidues;

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
    /// The taxonomy database of this run, opened at the first name asked for (outfmt 15 and 16,
    /// and the batches' prefetch of outfmt 0, 5, 6 and 7).
    taxdb: shared::TaxdbRun,
    /// The program that LOSAT's refusals name: BLASTP or PSIBLAST (`BlastpApp::name`).
    program: &'static str,
}

impl BlastpDatabase {
    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(super) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    /// The ids of the queries of the search (`SeqDb::id_scope`).
    pub(super) fn in_id_scope(
        mut self,
        id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    ) -> Self {
        self.db.set_id_scope(id_scope);
        self
    }

    /// Opens the database of `-db name` (`algorithm::common::database::open_database`)
    /// with the command line's errors, which name `program`.
    pub(super) fn open(
        name: &str,
        settings: &ApplicationSettings,
        program: &'static str,
    ) -> Result<Self> {
        Self::open_restricted(name, settings, program, None, &mut std::io::sink())
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line (`shared::restrict_database`, whose warnings go to `diagnostics`).
    pub(super) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        program: &'static str,
        list: Option<&crate::seqdb::DbList>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Protein, settings, program)?;
        // A database with parsed Seq-ids (line DBID): BLASTP (ID2a) and psiblast read it.
        shared::check_seqids(&db, settings, program)?;
        shared::restrict_database(&mut db, list, diagnostics, program)?;
        Ok(Self {
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
            program,
        })
    }

    /// The `-db` value as NCBI keeps it.
    pub(super) fn name(&self) -> &str {
        &self.name
    }

    /// Whether `taxdb.bti` is on the search path, asked once for the run
    /// (`shared::TaxdbRun::found`): the batches' prefetch reads the taxonomy database only then
    /// (`ListedTitles::prefetch_messages`).
    pub(super) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
    }

    /// The number of sequences (`CSeqDB::GetNumSeqs`).
    pub(super) fn num_sequences(&self) -> u32 {
        self.db.num_seqs()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:145-150
    // ```c++
    // static Int4
    // s_SeqDbGetMinLength(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetMinLength();
    // }
    // ```
    /// The shortest sequence length of the preliminary hit cutoff (`BlastSeqSrcGetMinSeqLen`,
    /// `min_subject_length` of `BLAST_GapAlignSetUp`): `BLAST_SEQSRC_MINLENGTH` for a whole
    /// database, the scanned shortest length with a list (line DI-8: `SeqDb::min_length`).
    pub(super) fn cutoff_subject_length(&self) -> i32 {
        match self.db.min_length() {
            Some(shortest) => i32::try_from(shortest).unwrap_or(i32::MAX),
            None => BLAST_SEQSRC_MINLENGTH,
        }
    }

    /// The total length the statistics read (`SeqDb::statistics_numbers`).
    pub(super) fn total_letters(&self) -> u64 {
        self.db.statistics_numbers().0
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(super) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
    }

    /// The taxids of the search's list (`SeqDb::list_taxids`).
    pub(super) fn list_taxids(&self) -> &[i32] {
        self.db.list_taxids()
    }

    /// The number of sequences the statistics read (`SeqDb::statistics_numbers`).
    pub(super) fn statistics_num_sequences(&self) -> u32 {
        self.db.statistics_numbers().1
    }

    /// The number of sequences the reports print: the whole database's, also with a taxid
    /// list (`SeqDb::printed_num_seqs`).
    pub(super) fn printed_num_sequences(&self) -> u32 {
        self.db.printed_num_seqs()
    }

    /// The total length the reports print (`SeqDb::printed_total_length`).
    pub(super) fn printed_total_letters(&self) -> u64 {
        self.db.printed_total_length()
    }

    /// The total length as the engine holds it.
    pub(super) fn total_letters_usize(&self) -> Result<usize> {
        let total = self.total_letters();
        usize::try_from(total)
            .ok()
            .filter(|&letters| i64::try_from(letters).is_ok())
            .with_context(|| {
                format!(
                    "a database of {total} letters is not supported by LOSAT's {} on this platform",
                    self.program
                )
            })
    }

    /// What the `Database:` lines of the pairwise report show
    /// (`algorithm::common::database::pairwise_database`): the title (the name as given
    /// for an empty title), the date and the real totals, also with `-dbsize`.
    pub(super) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// The names of the taxonomy database that the rows of the custom tabular report `spec`
    /// look up for the subjects `taxonomy` (`shared::tabular_names`, with the run's `TaxdbRun`);
    /// `program` names NCBI's messages and LOSAT's refusals (`BLASTP`, `PSIBLAST`).
    pub(super) fn tabular_names(
        &self,
        spec: &crate::report::tabular::TabularSpec,
        taxonomy: &[crate::report::tabular::TabularTaxonomy],
        program: &str,
    ) -> Result<Option<crate::report::tabular::TabularNames>> {
        shared::tabular_names(&self.db, &self.taxdb, spec, taxonomy, program)
    }

    /// The open of the taxonomy database by the rows of a query batch (`shared::TaxdbRun::open_for_rows`).
    pub(super) fn open_for_tabular_rows(&self, program: &str) -> Result<()> {
        self.taxdb.open_for_rows(&self.db, program)
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
/// NCBI's arithmetic overflows, is rejected (the refusal names `program`).
pub(super) fn blastp_database_search_space(
    query_len: usize,
    db_length: i64,
    db_num_seqs: i64,
    params: &KarlinParams,
    program: &str,
) -> Result<(SearchSpace, i64)> {
    let overflow = || {
        format!(
            "a database length of {db_length} letters with a query of {query_len} letters overflows NCBI BLAST+'s 64-bit effective search space, which is not supported by LOSAT's {program}"
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
// (`ListedSubjects`). A stored gap character (0, which makeblastdb writes for `-` in a
// protein FASTA) is searched as the residue it is in NCBI: row and column 0 stay
// `BLAST_SCORE_MIN`, so an ungapped extension stops at it and a gapped alignment goes
// round it with a gap, and NCBI's reports show it as `-` (`recount_shown_identities`).
// A byte that is no NCBIstdaa residue (28 or more) indexes NCBI's matrix past its rows
// and is rejected.
/// Prepares the block `block` (the sequences from the ordinal id `range.first_ordinal`)
/// for the search: O is read as X, and a sequence that LOSAT cannot search as NCBI does
/// is an error that names its ordinal id (and `program`).
fn prepare_search_block(block: &mut SequenceBlock, program: &str) -> Result<()> {
    let sequences = block.len();
    let invalid = read_o_as_x(block.protein_bytes_mut());
    // A well-formed block holds one sentinel byte before its first sequence and one
    // after every sequence.
    let sentinels_in_place = (0..sequences).all(|i| {
        let with_sentinels = block.protein_with_sentinels(i);
        (i > 0 || with_sentinels[0] == NCBISTDAA_GAP)
            && with_sentinels[with_sentinels.len() - 1] == NCBISTDAA_GAP
    });
    if !invalid && sentinels_in_place {
        return Ok(());
    }
    for i in 0..sequences {
        let ordinal = u64::from(block.first_ordinal()) + i as u64;
        shared::check_protein_sequence(block.protein_with_sentinels(i), ordinal, program)?;
    }
    bail!("a database block is not well-formed, which is not supported by LOSAT's {program}")
}

// No NCBI counterpart: LOSAT's check of a block's bytes (`prepare_search_block`; O read as X is
// NCBI's rule, quoted there). Every byte is written back (O as X, any other byte unchanged) and
// the flag is made with `u8` lanes, so the compiler works on 16 bytes (32 with AVX2) at a time
// with no branch per byte. With a branch for the write of each byte and a `usize` count of the
// gap bytes, the loop took 38 % of the cycles of a one-thread blastp-fast `-db` search (one
// 357-residue query against UniRef50 volume 00). The values are those of the byte-by-byte loop.
/// Reads O as X in `bytes`; returns whether a byte is no NCBIstdaa residue (28 or more).
fn read_o_as_x(bytes: &mut [u8]) -> bool {
    const _: () = assert!(BLASTAA_SIZE <= u8::MAX as usize);
    const BLASTAA_SIZE_U8: u8 = BLASTAA_SIZE as u8;
    let mut invalid = 0u8;
    for residue in bytes.iter_mut() {
        let byte = if *residue == NCBISTDAA_O {
            ncbistdaa::X
        } else {
            *residue
        };
        *residue = byte;
        invalid |= u8::from(byte >= BLASTAA_SIZE_U8);
    }
    invalid != 0
}

/// Reads the block `range` and prepares it for the search (`prepare_search_block`); the errors
/// name `program`.
fn read_search_block(db: &SeqDb, range: BlockRange, program: &str) -> Result<SequenceBlock> {
    let mut block = db
        .read_block(range)
        .map_err(|error| shared::search_read_error(error, program))?;
    prepare_search_block(&mut block, program)?;
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
///
/// No search result depends on this value: from a cleared table, every starting offset
/// from `window` up to INT4_MAX / 4 gives the same decisions (a cleared entry reads as
/// `diff >= window`, or `diff >= 0` for one hit, and its flag is 0). A block starts at
/// the one-thread offset so that its table holds the values of the one-thread table and
/// the INT4_MAX / 4 reset (which keeps the 31-bit `last_hit` in range) falls on the
/// same sequences.
///
/// The lengths come from the index entries that the block's read reads (`Volume::read_block`):
/// when they cannot be read here, the offset is returned unchanged, and the block's read meets
/// the same error, which the round returns.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn diag_offset_after(db: &SeqDb, range: BlockRange, mut diag_offset: i32, window: i32) -> i32 {
    let Ok(offsets) = db.volumes()[range.volume]
        .index()
        .sequence_offsets(range.first, range.end)
    else {
        return diag_offset;
    };
    for bounds in offsets.windows(2) {
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
/// into `x_compact_hit_lists` in that order (NCBI's loop with one thread).
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
/// before the next round is read: during this stage a worker holds one block at a
/// time, so the sequence memory is bounded by the threads and the block size, not by
/// the database (the results of a round wait for the merge; after the stage,
/// `ListedSubjects` holds the sequences of the listed subjects). A block that cannot be
/// read or searched does not stop the other blocks of its round; the round's first
/// error in ordinal order is returned, the error of the one-thread order (the blocks
/// are cut alike). The hit lists hold the compact element (`x_prelim_compact.rs`).
/// EXPERIMENT (LOSAT_X_BPDBCUT): with more than one thread the blocks are cut by the
/// database size and the thread count, and a failed search returns the one-thread order's
/// error (`database/x_db_cut.rs`, exactness argument in its header).
#[allow(clippy::too_many_arguments)]
pub(super) fn search_blocks<S, K, P>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
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
    // EXPERIMENT (LOSAT_X_BPDBCUT): `database/x_db_cut.rs`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        if x_db_cut::x_bp_db_cut() {
            return x_db_cut::x_search_blocks_cut(
                database,
                pool,
                window,
                blocks_per_thread,
                new_scratch,
                process_subject,
                sink,
            );
        }
    }
    search_blocks_of(
        database,
        pool,
        window,
        BLOCK_BYTES,
        blocks_per_thread,
        new_scratch,
        process_subject,
        sink,
    )
}

/// `search_blocks` with blocks of at most `block_bytes` and rounds of `blocks_per_thread`
/// blocks per thread.
#[allow(clippy::too_many_arguments)]
fn search_blocks_of<S, K, P>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink + Send,
{
    // EXPERIMENT (LOSAT_X_BPDBKEEP): `database/x_db_keep.rs`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if x_db_keep::x_bp_db_keep() {
        return x_db_keep::x_search_blocks_kept(
            database,
            pool,
            window,
            block_bytes,
            blocks_per_thread,
            new_scratch,
            process_subject,
            sink,
        );
    }
    search_blocks_made(
        database,
        pool,
        window,
        block_bytes,
        blocks_per_thread,
        new_scratch,
        process_subject,
        sink,
    )
}

/// `search_blocks_of` without LOSAT_X_BPDBKEEP: each task makes its scratch and reads its block
/// into a new buffer.
#[allow(clippy::too_many_arguments)]
fn search_blocks_made<S, K, P>(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
    window: i32,
    block_bytes: u64,
    blocks_per_thread: usize,
    new_scratch: &S,
    process_subject: &K,
    sink: &mut P,
) -> Result<()>
where
    S: Fn() -> BlastpSubjectScratch + Sync,
    K: Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult) + Sync,
    P: PrelimSink,
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
                            let block = read_search_block(db, range, database.program)?;
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
                                if !subject_result.x_compact_hsp_lists.is_empty() {
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
                    sink.add(&mut subject_result);
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
        let block = read_search_block(db, range, database.program)?;
        for i in 0..block.len() {
            process_subject(
                &mut scratch,
                range.first_ordinal as usize + i,
                block_subject(&block, i),
                &mut subject_result,
            );
            sink.add(&mut subject_result);
        }
        next_ordinal = range.end_ordinal();
    }
    Ok(())
}

/// The database sequences that the preliminary hit lists name, read by ordinal id after
/// the preliminary stage: they stand for the subject set in the composition redo. Their
/// definition lines are read after the redo, for the final hit lists only
/// (`read_titles`).
///
/// From here on a hit list names a sequence by its rank among the listed ordinal ids
/// (`ordinals[rank]`), so that the code of the `-subject` path indexes these vectors as
/// it indexes the encoded subject set. The rank grows with the ordinal id, so every
/// comparison that breaks a tie by the ordinal id (`s_EvalueCompareHSPLists`,
/// `ScoreCompareHSPs`; `hsp.rs`) gives the same order.
pub(super) struct ListedSubjects {
    /// The listed ordinal ids, ascending.
    pub(super) ordinals: Vec<u32>,
    /// The sequences as stored (O stays O), between two sentinel bytes, by rank: empty
    /// until the redo has run (it reads each match's subject from the database, `x_list`), and
    /// then what the reports read (`x_read_shown`): a subject whose residues no report reads holds
    /// its length only (`aa_seq` empty).
    pub(super) encoded: Vec<EncodedProtein>,
    /// The number of residues of each subject, by rank (`x_list`).
    x_lengths: Vec<u32>,
    /// Whether a read of the redo or of the first Newton problems found a pyrrolysine (O) or a
    /// gap character, by rank (`x_list`).
    x_recounted: Vec<std::sync::atomic::AtomicBool>,
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
    /// The listed subjects of a search with every residue read, as the reports of a search
    /// that makes records read them after the redo (`x_list`, then `x_read_shown` with
    /// `XShownResidues::Listed`); the unit tests of the titles start here. The searches
    /// themselves read each match's subject from the database in the redo (`x_list`).
    #[cfg(test)]
    pub(super) fn read(
        database: &BlastpDatabase,
        hit_lists: &mut [Option<BlastpHitList>],
        x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    ) -> Result<Self> {
        let mut listed = Self::x_list(database, hit_lists, x_compact_hit_lists)?;
        listed.x_read_shown(database, &[], XShownResidues::Listed, &[])?;
        Ok(listed)
    }

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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2302-2310
    // ```c++
    // CTempString
    // CSeqDBVol::x_GetHdrAsn1Binary(int oid) const
    // {
    //     TIndx hdr_start = 0;
    //     TIndx hdr_end   = 0;
    //
    //     //m_Atlas.Lock(locked);
    //
    //     if (!m_HdrFileOpened) x_OpenHdrFile();
    // ```
    // NCBI first reads a header file (and opens it, seqdbvol.cpp:127-133) when the
    // traceback's results become Seq-aligns: after the batch's traceback, for the final
    // hit lists only, query by query and in the order of each hit list. A header file
    // that is missing or cannot be read fails the batch there (`BLAST Database error`,
    // exit status 2, blast_app_util.hpp:185-203, outside the search stage's handler of
    // local_blast.cpp:240-251), after the outfmt 0 prolog and the reports of the batches
    // before; a batch without hits reads no header. Shown with the oracle (campaign
    // dbsearch-20261012, `scripts/db3_nophr_probe.py`; comparison cases `cli.path.no_phr*`).
    //
    // NCBI reference for the titles:
    // c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:970-1017 (`s_GetBioseqTitle`: the
    // title of the first definition line, the others appended as ` >id title`).
    /// Reads the definition lines of the subjects of the final hit lists (`hit_lists`,
    /// after the redo and before anything is dropped from them), in NCBI's order. A reader
    /// error of NCBI's own is a `HeaderReadError` that names the query whose hit list
    /// reached it; LOSAT's refusals are errors of their own.
    pub(super) fn read_titles(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<ListedTitles<'_>> {
        let mut titles: Vec<Option<(Vec<u8>, bool)>> = vec![None; self.ordinals.len()];
        let mut parsed = ParsedSubjects::default();
        let headers = shared::HeaderReader::new(&database.db);
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            for hsp_list in &hit_list.hsplist_array[..count] {
                let rank = hsp_list.oid as usize;
                if titles[rank].is_some() {
                    continue;
                }
                let ordinal = self.ordinals[rank];
                titles[rank] = Some(parsed.read(
                    database,
                    &headers,
                    rank,
                    ordinal,
                    query,
                    self.ordinals.len(),
                )?);
            }
        }
        let mut records = Vec::with_capacity(self.ordinals.len());
        let mut replaced_title_bytes = Vec::with_capacity(self.ordinals.len());
        for ((&ordinal, encoded), title) in self.ordinals.iter().zip(&self.encoded).zip(titles) {
            let letters: Vec<u8> = encoded.aa_seq[1..=encoded.aa_len]
                .iter()
                .map(|&residue| ncbistdaa_to_ascii(residue) as u8)
                .collect();
            // A subject that is in no final hit list keeps an empty title: no report
            // shows it.
            let (title, replaced) = title.unwrap_or_default();
            replaced_title_bytes.push(replaced);
            records.push(FastaRecord::new(
                // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the
                // OID (line DI-8).
                format!("gnl|BL_ORD_ID|{}", database.db.oid(ordinal)),
                &title,
                &letters,
            ));
        }
        Ok(ListedTitles {
            listed: self,
            shown: ListedShown::Records(records),
            replaced_title_bytes,
            parsed: parsed.by_rank,
            program: database.program,
        })
    }
}

/// The ordinal ids of the subjects of `hit_lists` and `x_compact_hit_lists`, ascending, each once.
fn listed_ordinals(
    hit_lists: &[Option<BlastpHitList>],
    x_compact_hit_lists: &[Option<BlastpHitList<XPrelimHsp>>],
) -> Vec<u32> {
    let mut ordinals: Vec<u32> = hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| hit_list.hsplist_array.iter().map(|hsp_list| hsp_list.oid))
        .chain(
            x_compact_hit_lists
                .iter()
                .flatten()
                .flat_map(|hit_list| hit_list.hsplist_array.iter().map(|hsp_list| hsp_list.oid)),
        )
        .collect();
    ordinals.sort_unstable();
    ordinals.dedup();
    ordinals
}

/// Stored residues between two sentinel bytes.
fn encoded_stored(residues: &[u8]) -> EncodedProtein {
    // (A sequence with a gap residue or a foreign byte never reaches a hit
    // list: `prepare_search_block` rejects its block.)
    let mut aa_seq = Vec::with_capacity(residues.len() + 2);
    aa_seq.push(NCBISTDAA_GAP);
    aa_seq.extend_from_slice(residues);
    aa_seq.push(NCBISTDAA_GAP);
    EncodedProtein {
        aa_seq,
        aa_len: residues.len(),
    }
}

/// Renames the subjects of `hit_lists` and `x_compact_hit_lists` to their ranks among `ordinals`.
fn rename_to_ranks(
    ordinals: &[u32],
    hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
) {
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
    for hit_list in x_compact_hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            hsp_list.oid = rank(hsp_list.oid);
            for hsp in &mut hsp_list.hsps {
                hsp.hsp.s_idx = rank(hsp.hsp.s_idx);
            }
        }
    }
}

/// The stored titles and letters of the listed subjects as the reports read them, by rank
/// (`ListedSubjects::read_titles`).
pub(super) struct ListedTitles<'a> {
    listed: &'a ListedSubjects,
    /// The titles (and letters) as the reports read them, by rank.
    pub(super) shown: ListedShown,
    /// Whether the stored title of a subject holds a byte outside ' ' ..= '~' (read as
    /// '#'), by rank (`check_reported_titles`).
    replaced_title_bytes: Vec<bool>,
    /// The subjects with parsed Seq-ids, by rank (empty when no reported subject has them;
    /// their entry of `shown` is the title of their first definition line).
    pub(super) parsed: Vec<Option<Arc<ParsedSubject>>>,
    /// The program that LOSAT's refusals name (`BlastpDatabase::program`).
    program: &'static str,
}

/// The reported subjects with parsed Seq-ids of `ListedSubjects::read_titles`, by rank.
#[derive(Default)]
pub(super) struct ParsedSubjects {
    pub(super) by_rank: Vec<Option<Arc<ParsedSubject>>>,
}

impl ParsedSubjects {
    /// Reads the header of the subject `ordinal` of rank `rank` (of `ranks`) of `database` with
    /// `headers` (its reader), where the hit list of query `query` reaches it
    /// (`shared::read_header`): its stored title and whether it held a replaced byte, or for a
    /// sequence with parsed Seq-ids, kept here, the title of its first definition line.
    pub(super) fn read(
        &mut self,
        database: &BlastpDatabase,
        headers: &shared::HeaderReader<'_>,
        rank: usize,
        ordinal: u32,
        query: usize,
        ranks: usize,
    ) -> Result<(Vec<u8>, bool)> {
        match shared::read_header(headers, ordinal, query, database.program)? {
            shared::ReportedHeader::Stored(stored) => Ok((stored.title, stored.replaced_bytes)),
            shared::ReportedHeader::Parsed(subject) => {
                let title = subject.first_title().to_vec();
                if self.by_rank.is_empty() {
                    self.by_rank = vec![None; ranks];
                }
                self.by_rank[rank] = Some(subject);
                Ok((title, false))
            }
        }
    }
}

/// The titles (and letters) of the listed subjects as the reports read them, by rank.
pub(super) enum ListedShown {
    /// The stored titles (empty for a subject in no final hit list) and the stored letters.
    Records(Vec<FastaRecord>),
    /// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the stored titles, where the subjects in no final hit
    /// list share one empty title; the letters are the stored residues
    /// (`ListedSubjects::encoded`, `x_listed_lean.rs`).
    Lean(Vec<Arc<[u8]>>),
}

impl ListedTitles<'_> {
    /// The stored title of the subject of rank `rank` (empty for a subject in no final hit
    /// list).
    fn title(&self, rank: usize) -> &[u8] {
        match &self.shown {
            ListedShown::Records(records) => &records[rank].title,
            ListedShown::Lean(titles) => &titles[rank],
        }
    }

    /// The records of `read_titles`.
    #[cfg(test)]
    fn records(&self) -> &[FastaRecord] {
        match &self.shown {
            ListedShown::Records(records) => records,
            ListedShown::Lean(_) => panic!("lean titles have no records"),
        }
    }

    /// LOSAT's rejection of a reported subject (one with rows in the final `hit_lists`)
    /// whose stored title holds a byte outside ' ' ..= '~'
    /// (`algorithm::common::database::check_reported_title`): NCBI writes an error message
    /// for the byte each time its formatter reads the definition line.
    pub(super) fn check_reported_titles(&self, hit_lists: &[Option<BlastpHitList>]) -> Result<()> {
        for rank in Self::reported(hit_lists) {
            shared::check_reported_title(
                self.replaced_title_bytes[rank],
                self.listed.ordinals[rank],
                self.program,
            )?;
        }
        Ok(())
    }

    /// The Blast-def-line taxid of each reported subject by rank (0 for none and for a
    /// subject in no final hit list), which the XML2 and JSON formats write in the `Hit`
    /// description, and what the taxonomy database says about it (for a subject with parsed
    /// Seq-ids, about the taxid of each definition line: `shared::read_reported_taxonomy`).
    /// They are read apart from the titles, for those formats only.
    ///
    /// NCBI writes the scientific name of a taxid when it finds `taxdb.bti` on the search path
    /// (`data4xml2format.cpp x_InitCommon`, `blastxml2_format.cpp:380-386`; validity is not
    /// checked): `taxonomy[rank]` is then `Named`, or `NotFound` when the taxonomy database does
    /// not know the taxid (the description is the id alone). Without `taxdb.bti`, or without a
    /// taxid, `taxonomy` is empty. A scientific name that NCBI's XML stream cannot write, and a
    /// taxonomy database that LOSAT cannot read as NCBI does, are refusals made here, before any
    /// byte is written.
    ///
    /// NCBI's messages about the taxonomy database go to standard error: those of opening and of
    /// the lookup that the volume makes for each taxid of a fetched subject (`CSeqDBVol::
    /// x_GetTaxonomy`, once per volume and taxid) are returned as `TaxdbNotes`, which
    /// `write_blastp_xml` writes at the first `Hit` that needs them, and those of the lookup of
    /// each `Hit` are written there too.
    pub(super) fn read_taxids(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<(
        Vec<i32>,
        Vec<HitTaxonomy>,
        Vec<Vec<HitTaxonomy>>,
        shared::TaxdbNotes,
    )> {
        shared::read_reported_taxonomy(
            &database.db,
            &database.taxdb,
            &self.listed.ordinals,
            &Self::reported(hit_lists),
            &self.parsed,
            self.program,
        )
    }

    /// NCBI's messages about the taxonomy database when the query batch whose final hit lists
    /// are `hit_lists` prefetches the subjects it reports (`shared::prefetch_messages`): the
    /// subjects with rows in them, the ones `read_taxids` reads. Empty, and nothing read,
    /// without `taxdb.bti` on the search path. A taxonomy database that LOSAT cannot read as
    /// NCBI does, and a header whose taxid does not decode, are refused.
    pub(super) fn prefetch_messages(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<Vec<u8>> {
        shared::prefetch_messages(
            &database.db,
            &database.taxdb,
            Self::reported_ranks(hit_lists)
                .map(|rank| (self.listed.ordinals[rank], self.parsed_subject(rank))),
            self.program,
        )
    }

    /// The ranks of the subjects with rows in `hit_lists` (the final hit lists).
    fn reported(hit_lists: &[Option<BlastpHitList>]) -> std::collections::BTreeSet<usize> {
        Self::reported_ranks(hit_lists).collect()
    }

    /// The rank of the subject of each row of `hit_lists`, in the order of the lists.
    fn reported_ranks(hit_lists: &[Option<BlastpHitList>]) -> impl Iterator<Item = usize> + '_ {
        hit_lists.iter().flatten().flat_map(|hit_list| {
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            hit_list.hsplist_array[..count]
                .iter()
                .flat_map(|hsp_list| hsp_list.hsps.iter().map(|hsp| hsp.s_idx as usize))
        })
    }

    /// The subject id of the tabular formats (`sseqid`, `sacc`, `saccver`) of every
    /// listed subject, by rank (`report::defline::database_subject_id`); a subject
    /// without rows gets an empty id. A subject with rows whose id LOSAT cannot write as
    /// NCBI is an error before any format is written.
    pub(super) fn tabular_subject_ids(
        &self,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<(Vec<Arc<[u8]>>, Vec<Option<[Arc<[u8]>; 2]>>)> {
        let reported = Self::reported(hit_lists);
        let ranks = self.listed.ordinals.len();
        let mut ids: Vec<Arc<[u8]>> = Vec::with_capacity(ranks);
        // `sacc` and `saccver` of the subjects with parsed Seq-ids, by rank (empty without
        // them; the others' are their `sseqid`).
        let mut accessions: Vec<Option<[Arc<[u8]>; 2]>> = Vec::new();
        // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): one empty id for every subject without rows.
        let empty: Arc<[u8]> = Arc::from(&b""[..]);
        for rank in 0..ranks {
            if !reported.contains(&rank) {
                ids.push(match self.shown {
                    ListedShown::Lean(_) => Arc::clone(&empty),
                    ListedShown::Records(_) => Arc::from(&b""[..]),
                });
                continue;
            }
            let ordinal = self.listed.ordinals[rank];
            if let Some(subject) = self.parsed_subject(rank) {
                let [id, accession, version] = subject.tabular_ids(true).map_err(|problem| {
                    shared::parsed_refusal(ordinal, &format!("has {problem}"), self.program)
                })?;
                ids.push(Arc::from(id));
                if accessions.is_empty() {
                    accessions = vec![None; ranks];
                }
                accessions[rank] = Some([Arc::from(accession), Arc::from(version)]);
                continue;
            }
            let id = shared::tabular_subject_id(self.title(rank), ordinal, true, self.program)?;
            ids.push(Arc::from(id));
        }
        Ok((ids, accessions))
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
    /// The taxonomy of the custom tabular fields (`staxid`, `staxids`;
    /// `report::tabular::TabularTaxonomy::of_deflines`) of the reported subjects (those with
    /// rows in `hit_lists`), by rank: from the definition lines of a subject with parsed Seq-ids,
    /// else from its stored ones (`SeqDb::defline_set`, which keeps whether a line has a
    /// `taxid`); the other ranks are left empty. A header that does not decode is refused before
    /// the report, as BLASTN's (`BlastnDatabase::tabular_taxonomy`).
    pub(super) fn tabular_taxonomy(
        &self,
        database: &BlastpDatabase,
        hit_lists: &[Option<BlastpHitList>],
    ) -> Result<Vec<crate::report::tabular::TabularTaxonomy>> {
        use crate::report::tabular::TabularTaxonomy;
        let ranks = self.listed.ordinals.len();
        let mut taxonomy = vec![TabularTaxonomy::default(); ranks];
        for rank in Self::reported(hit_lists) {
            if rank >= ranks {
                continue;
            }
            taxonomy[rank] = match self.parsed_subject(rank) {
                Some(subject) => TabularTaxonomy::of_deflines(subject.deflines()),
                None => {
                    let set = database
                        .db
                        .defline_set(self.listed.ordinals[rank])
                        .map_err(|error| shared::database_error(error, self.program))?;
                    TabularTaxonomy::of_deflines(&set.deflines)
                }
            };
        }
        Ok(taxonomy)
    }

    /// The subject of rank `rank` when it has parsed Seq-ids.
    pub(super) fn parsed_subject(&self, rank: usize) -> Option<&ParsedSubject> {
        self.parsed.get(rank).and_then(Option::as_deref)
    }

    /// LOSAT's refusals of the subjects with parsed Seq-ids that the reports of the kinds asked
    /// for do not write as NCBI, before any format is written: a subject that NCBI's formatters
    /// cannot find by its Seq-align id where a format reads its Bioseq (`fetched`: outfmt 0, 5,
    /// 6, 7, 13 to 16; `ParsedSubject::fetched`), two of them that follow each
    /// other in a hit list with matching Seq-align ids where outfmt 0, 5 or 13 to 16 is written
    /// (`joined`, `shared::refuse_joined_hits`). `hit_lists` are in the order of the reports.
    pub(super) fn check_parsed_reports(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        fetched: bool,
        joined: bool,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        if fetched {
            for rank in Self::reported(hit_lists) {
                if let Some(subject) = self.parsed_subject(rank) {
                    subject.fetched().map_err(|reason| {
                        shared::parsed_refusal(self.listed.ordinals[rank], &reason, self.program)
                    })?;
                }
            }
        }
        if joined {
            for (query, hit_list) in hit_lists.iter().enumerate() {
                let Some(hit_list) = hit_list else {
                    continue;
                };
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                shared::refuse_joined_hits(
                    query,
                    hit_list.hsplist_array[..count].iter().map(|hsp_list| {
                        let rank = hsp_list.oid as usize;
                        self.parsed_subject(rank)
                            .map(|subject| (self.listed.ordinals[rank], subject))
                    }),
                    self.program,
                )?;
            }
        }
        Ok(())
    }

    /// LOSAT's refusal of psiblast's cut of its traceback hit lists `hit_lists` (in the order of
    /// the reports) to `max_target_seqs` subjects where NCBI's formatter keeps one more, for
    /// every format: a list longer than `max_target_seqs` whose first `max_target_seqs + 1` HSP
    /// lists hold two parsed subjects that follow each other with matching Seq-align ids
    /// (`shared::refuse_joined_cut`).
    pub(super) fn check_joined_cut(
        &self,
        hit_lists: &[Option<BlastpHitList>],
        max_target_seqs: usize,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        for (query, hit_list) in hit_lists.iter().enumerate() {
            let Some(hit_list) = hit_list else {
                continue;
            };
            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
            if count <= max_target_seqs {
                continue;
            }
            shared::refuse_joined_cut(
                query,
                max_target_seqs,
                hit_list.hsplist_array[..count].iter().map(|hsp_list| {
                    let rank = hsp_list.oid as usize;
                    self.parsed_subject(rank)
                        .map(|subject| (self.listed.ordinals[rank], subject))
                }),
                self.program,
            )?;
        }
        Ok(())
    }

    /// The `stitle` text of every listed subject, by rank: the stored title (an empty
    /// one is written `N/A` by the shared tabular writer, `report::tabular`).
    pub(super) fn stitles(&self) -> Vec<Arc<[u8]>> {
        match &self.shown {
            ListedShown::Records(records) => records
                .iter()
                .map(|record| Arc::from(record.title()))
                .collect(),
            ListedShown::Lean(titles) => titles.clone(),
        }
    }

    /// LOSAT's rejection of the outfmt 0 titles that it does not write as NCBI, of the
    /// subjects that the report shows (`shown`, by rank): a title of punctuation that
    /// NCBI's `x_CleanAndCompress` reads past the end of when it makes the title of the
    /// description table (NCBI BLAST+ 2.17.0 crashes there), as for the records of a
    /// subject set (`check_shown_subject_titles`).
    pub(super) fn check_shown_titles(
        &self,
        database: &BlastpDatabase,
        shown: &std::collections::BTreeSet<usize>,
    ) -> Result<()> {
        for &rank in shown {
            let ordinal = self.listed.ordinals[rank];
            match self.parsed_subject(rank) {
                Some(subject) => shared::check_shown_parsed_subject(
                    subject,
                    true,
                    database.db.taxdb_found(),
                    ordinal,
                    self.program,
                )?,
                None => shared::check_shown_protein_title(self.title(rank), ordinal, self.program)?,
            }
        }
        Ok(())
    }
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search writes a block's rows before the next block is
// searched (`x_query_blocks.rs`), so with `-db` it may run only where the search of all queries
// would not fail after some rows. Past the preliminary stage (whose reader errors come in the
// first block, which reads the whole database as every batch does), the failures of a database
// search that depend on its hits are those of the definition lines of the reported subjects:
// NCBI's reader error (`HeaderReadError`, the batch of the query that reaches it fails after the
// reports of the batches before) and LOSAT's refusals (`ListedSubjects::read_titles` and
// `ListedSubjects::x_read_titles_lean`,
// `ListedTitles::check_reported_titles`, `ListedTitles::tabular_subject_ids`, and for outfmt 0
// `ListedTitles::check_shown_titles`). When every definition line of the database reads and none
// of them is refused, none of these failures can happen in any batch. With `taxdb.bti` on the
// search path, the prefetch of each batch (`ListedTitles::prefetch_messages`) also reads the
// reported subjects' taxids and looks them up in the taxonomy database, which can refuse as
// XML2's lookups can (DBID item T): `taxids` then also asks that every sequence's taxid decodes
// and that no lookup can be refused (`shared::taxdb_may_refuse` without the names, which the
// prefetch does not ask for), as `x_titles_taxids_allow_blocks` asks for outfmt 15 and 16.
// The prefetch's refusal of NCBI's taxonomy cache (`shared::TaxdbRun::cache_collision`) comes at
// the run's first prefetch, in block 0 before its rows.

// No NCBI counterpart: the order of LOSAT's own checks; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPDBKEEP): whether `check` holds for every range of
/// `db.header_ranges(HEADER_RANGE_SEQS)`, read in `pool` as the title checks read them, with each
/// range made from its index (`SeqDb::x_header_range`): no list of the ranges is allocated, whose
/// free (1.5 MiB for UniRef50) would raise glibc's mmap threshold before the search.
fn x_all_header_ranges(
    db: &crate::seqdb::SeqDb,
    pool: &crate::utils::threading::SearchPool<'_>,
    check: impl Fn(&BlockRange) -> bool + Sync,
) -> bool {
    let seqs = crate::seqdb::HEADER_RANGE_SEQS;
    let count = db.x_header_range_count(seqs);
    let check_index = |index: usize| check(&db.x_header_range(seqs, index));
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() {
        return pool.install(|| (0..count).into_par_iter().all(check_index));
    }
    let _ = pool;
    (0..count).all(check_index)
}

// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): whether every sequence of the database has one definition line
/// that reads, holds no byte outside ' ' ..= '~' and gives a tabular subject id as NCBI writes
/// it, or with `pairwise` (an outfmt 0 report, which writes no tabular id) a
/// title that the outfmt 0 report writes as NCBI, and with `taxids` (the batches' prefetch reads
/// the taxonomy database) a taxid that decodes, where no lookup of the taxonomy database can be
/// refused. Read in `pool` over the ordinal ids; a reader error is `false` (the search of all
/// queries then reports it where NCBI does).
/// EXPERIMENT (LOSAT_X_DBIDXREAD): the definition lines are read in ranges of
/// `HEADER_RANGE_SEQS` sequences (`SeqDb::visit_deflines`, the same results as
/// `SeqDb::deflines`; with `taxids`, `SeqDb::visit_deflines_taxids`, the titles and the taxids in
/// one pass) instead of one sequence at a time.
pub(super) fn x_titles_allow_blocks(
    database: &BlastpDatabase,
    pairwise: bool,
    taxids: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    let allowed_deflines = |ordinal: u32,
                            deflines: std::result::Result<
        Vec<crate::seqdb::OrdinalDefline>,
        crate::seqdb::SeqDbError,
    >|
     -> bool {
        match deflines {
            Ok(deflines) => {
                deflines.len() == 1
                    && deflines[0].replaced_title_bytes == 0
                    && if pairwise {
                        shared::check_shown_protein_title(
                            &deflines[0].title,
                            ordinal,
                            database.program,
                        )
                        .is_ok()
                    } else {
                        crate::report::defline::database_subject_id(&deflines[0].title, true)
                            .is_ok()
                    }
            }
            Err(_) => false,
        }
    };
    if crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        // A record whose taxid does not decode is an `Err` (`visit_deflines_taxids`).
        let allowed_range = |range: &BlockRange| -> bool {
            if taxids {
                db.visit_deflines_taxids(*range, |ordinal, deflines| {
                    allowed_deflines(ordinal, deflines)
                })
            } else {
                db.visit_deflines(*range, |ordinal, deflines| {
                    allowed_deflines(ordinal, deflines)
                })
            }
        };
        if crate::seqdb::x_bp_db_keep() {
            return x_all_header_ranges(db, pool, allowed_range)
                && !(taxids && shared::taxdb_may_refuse(db, false));
        }
        let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let titles = if pool.enabled() {
            pool.install(|| ranges.par_iter().all(allowed_range))
        } else {
            ranges.iter().all(allowed_range)
        };
        #[cfg(not(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        )))]
        let titles = {
            let _ = pool;
            ranges.iter().all(allowed_range)
        };
        return titles && !(taxids && shared::taxdb_may_refuse(db, false));
    }
    let allowed = |ordinal: u32| -> bool { allowed_deflines(ordinal, db.deflines(ordinal)) };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let titles = if pool.enabled() {
        pool.install(|| {
            (0..db.num_seqs())
                .into_par_iter()
                .with_min_len(4096)
                .all(allowed)
        })
    } else {
        (0..db.num_seqs()).all(allowed)
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let titles = (0..db.num_seqs()).all(allowed);
    // One sequence at a time (`LOSAT_X_DBIDXREAD=0`): the taxids after the titles.
    titles && (!taxids || x_taxids_allow_blocks(database, false, pool))
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search of outfmt 8, 9, 12 and 5.
// These reports show no tabular id and no outfmt 0 title, so of the failures that
// `x_titles_allow_blocks` rules out only those of the definition lines' reading remain: NCBI's
// reader error (`HeaderReadError`), more than one definition line, and a stored title with a byte
// outside ' ' ..= '~' (`ListedSubjects::read_titles`, `ListedTitles::check_reported_titles`; for
// outfmt 5 that range is what `names::check_text` passes). A new function, so that
// `x_titles_allow_blocks` and the per-definition-line code of the outfmt 0, 6 and 7 path stay as
// they are.
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): whether every sequence of the database has
/// one definition line that reads and holds no byte outside ' ' ..= '~'. Read in `pool` as
/// `x_titles_allow_blocks` reads it; a reader error is `false`.
pub(super) fn x_titles_read_allow_blocks(
    database: &BlastpDatabase,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    let readable = |deflines: std::result::Result<
        Vec<crate::seqdb::OrdinalDefline>,
        crate::seqdb::SeqDbError,
    >|
     -> bool {
        match deflines {
            Ok(deflines) => deflines.len() == 1 && deflines[0].replaced_title_bytes == 0,
            Err(_) => false,
        }
    };
    if crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
        let readable_range = |range: &BlockRange| -> bool {
            db.visit_deflines(*range, |_ordinal, deflines| readable(deflines))
        };
        if crate::seqdb::x_bp_db_keep() {
            return x_all_header_ranges(db, pool, readable_range);
        }
        let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        if pool.enabled() {
            return pool.install(|| ranges.par_iter().all(readable_range));
        }
        let _ = pool;
        return ranges.iter().all(readable_range);
    }
    let allowed = |ordinal: u32| -> bool { readable(db.deflines(ordinal)) };
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

// EXPERIMENT (LOSAT_X_BPQBLOCK): the block search of outfmt 15 and 16. Besides what outfmt 5
// reads (`x_titles_read_allow_blocks`) they read, for each subject with hits, the Blast-def-line
// taxid (`SeqDb::defline_taxids`, a decode apart from the titles') and, when `taxdb.bti` is on the
// search path, the taxonomy database (`shared::read_reported_taxids`): its opening, the lookup of
// the taxid, and a scientific name that NCBI's stream can write (`names::check_text`). The search
// of all queries makes those refusals before its first byte; a block would make them after the
// earlier blocks' bytes. So the block search is taken only where none can fail for any subject
// of the database: every sequence's taxid reads, and `shared::taxdb_may_refuse` (OF5b, the one
// implementation of "can a lookup of this database refuse?", DECISIONS 34 and 38) says no. That
// check answers from the smaller side: one pass over the taxonomy database when it is the smaller
// one (NCBI's real taxdb has no refusable record, so no taxid lookup and no pass over the
// database's taxids), else the database's taxids, sorted and deduplicated, looked up. It has its
// own taxonomy handle, so the run's open messages stay for the report. This is conservative,
// since the report reads the taxids of the subjects with hits only.
// The titles and the taxid decodes are read in one pass over the headers (`SeqDb::visit_deflines_taxids`,
// one read per `HEADER_RANGE_SEQS` sequences, as `x_titles_read_allow_blocks` reads them); no set of
// taxids is built here.
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): `x_titles_read_allow_blocks` and `x_taxids_allow_blocks`
/// in one pass over the headers: whether every sequence of the database has one definition line
/// that reads and holds no byte outside ' ' ..= '~', its taxid reads, and no lookup of the
/// taxonomy database can be refused (`shared::taxdb_may_refuse`; with `names`, the outfmt 15 and
/// 16 that write the scientific names, also no name that their writer refuses; without, outfmt
/// 5's prefetch). Read in `pool`; a failure is `false`. It writes nothing and changes no state of
/// the run. With `LOSAT_X_DBIDXREAD=0` the two checks run one after the other, one sequence at a
/// time.
pub(super) fn x_titles_taxids_allow_blocks(
    database: &BlastpDatabase,
    names: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    if crate::seqdb::x_dbidxread() == crate::seqdb::XDbIdxRead::Off {
        return x_titles_read_allow_blocks(database, pool)
            && x_taxids_allow_blocks(database, names, pool);
    }
    let readable_range = |range: &BlockRange| -> bool {
        db.visit_deflines_taxids(*range, |_ordinal, deflines| match deflines {
            // The taxid `read_reported_taxids` takes of a sequence: its first definition
            // line's.
            Ok(deflines) => deflines.len() == 1 && deflines[0].replaced_title_bytes == 0,
            Err(_) => false,
        })
    };
    if crate::seqdb::x_bp_db_keep() {
        return x_all_header_ranges(db, pool, readable_range)
            && !shared::taxdb_may_refuse(db, names);
    }
    let ranges = db.header_ranges(crate::seqdb::HEADER_RANGE_SEQS);
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let readable = if pool.enabled() {
        pool.install(|| ranges.par_iter().all(readable_range))
    } else {
        ranges.iter().all(readable_range)
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let readable = {
        let _ = pool;
        ranges.iter().all(readable_range)
    };
    readable && !shared::taxdb_may_refuse(db, names)
}

/// EXPERIMENT (LOSAT_X_BPQBLOCK): whether every sequence's taxid reads and no lookup of the
/// taxonomy database can be refused (`shared::taxdb_may_refuse` with `names`, which is `false`
/// where `taxdb.bti` is not found). One sequence at a time (`LOSAT_X_DBIDXREAD=0`); see
/// `x_titles_taxids_allow_blocks`.
fn x_taxids_allow_blocks(
    database: &BlastpDatabase,
    names: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let db = &database.db;
    let reads = |ordinal: u32| db.defline_taxids(ordinal).is_ok();
    let count = db.num_seqs();
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let all_read = if pool.enabled() {
        pool.install(|| (0..count).into_par_iter().with_min_len(4096).all(reads))
    } else {
        (0..count).all(reads)
    };
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let all_read = {
        let _ = pool;
        (0..count).all(reads)
    };
    all_read && !shared::taxdb_may_refuse(db, names)
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
        BlastpDatabase::open(&fixture(name), &ApplicationSettings::default(), "BLASTP")
            .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    }

    // The byte-by-byte loop that `read_o_as_x` replaces.
    fn read_o_as_x_by_byte(bytes: &mut [u8]) -> bool {
        let mut invalid = false;
        for residue in bytes {
            if *residue == NCBISTDAA_O {
                *residue = ncbistdaa::X;
            }
            invalid |= usize::from(*residue) >= BLASTAA_SIZE;
        }
        invalid
    }

    // Lengths across vectors of 16 and 32 bytes, gap bytes from none to all, O, the last residue
    // (27), the first byte that is none (28) and 255: the bytes and the flag are those of the
    // byte-by-byte loop.
    #[test]
    fn read_o_as_x_gives_the_values_of_the_byte_by_byte_loop() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut flagged = 0usize;
        for length in [0usize, 1, 2, 223, 224, 225, 447, 448, 449, 1_000, 4_096] {
            for gap_every in [0u64, 1, 2, 7, 300] {
                for wrong in [None, Some(27u8), Some(28), Some(255)] {
                    let bytes: Vec<u8> = (0..length)
                        .map(|_| {
                            let r = next();
                            if gap_every != 0 && r % gap_every == 0 {
                                NCBISTDAA_GAP
                            } else if r % 11 == 0 {
                                NCBISTDAA_O
                            } else if r % 97 == 0 {
                                wrong.unwrap_or(1)
                            } else {
                                1 + (r % 25) as u8
                            }
                        })
                        .collect();
                    let mut runs = bytes.clone();
                    let mut by_byte = bytes;
                    let expected = read_o_as_x_by_byte(&mut by_byte);
                    flagged += usize::from(expected);
                    assert_eq!(
                        read_o_as_x(&mut runs),
                        expected,
                        "length {length}, a gap every {gap_every}, {wrong:?}"
                    );
                    assert_eq!(runs, by_byte);
                }
            }
        }
        assert!(flagged > 20, "cases with a byte of 28 or more: {flagged}");
        // Every byte a gap.
        let mut gaps = vec![NCBISTDAA_GAP; 1_000];
        assert!(!read_o_as_x(&mut gaps));
        // A byte of 28 alone: second, after the first vectors, last of the input; O next to it.
        for (length, at) in [(2usize, 1usize), (300, 224), (511, 510)] {
            let mut bytes = vec![1u8; length];
            bytes[at] = 28;
            assert!(read_o_as_x(&mut bytes), "28 at {at} of {length}");
        }
        let mut bytes = vec![1u8; 300];
        bytes[223] = NCBISTDAA_O;
        bytes[224] = NCBISTDAA_O;
        assert!(!read_o_as_x(&mut bytes));
        assert_eq!((bytes[223], bytes[224]), (ncbistdaa::X, ncbistdaa::X));
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
        let (space, exact) =
            blastp_database_search_space(379, 357_308, 1_000, &params, "BLASTP").unwrap();
        assert_eq!(exact, 84_653_016);
        assert_eq!(space.effective_space, 84_653_016.0);
        let (space, exact) =
            blastp_database_search_space(379, 1_000_000_000, 1_000, &params, "BLASTP").unwrap();
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
        let (space, exact) =
            blastp_database_search_space(379, 1, 1_000, &params, "BLASTP").unwrap();
        assert_eq!((space.length_adjustment, exact), (0, 379));
        // A total of five letters per sequence.
        let (space, exact) =
            blastp_database_search_space(379, 5_000, 1_000, &params, "BLASTP").unwrap();
        assert_eq!(space.length_adjustment, 4);
        assert_eq!(space.effective_db_len, 1_000.0);
        assert_eq!(exact, 1_000 * 375);
    }

    // NCBI's `Int8` product overflows for such a `-dbsize`; LOSAT rejects it.
    #[test]
    fn search_space_beyond_int8_is_rejected() {
        let error = blastp_database_search_space(
            379,
            9_000_000_000_000_000_000,
            1_000,
            &blosum62_gapped(),
            "BLASTP",
        )
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
        assert_eq!(shown.databases[0].date, "Oct 12, 2026  4:12 PM");
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
        let database =
            BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
        assert_eq!(database.pairwise_database().definition, name.as_bytes());
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // The first header of each volume is read when the database is opened (`shared::
    // check_seqids`), as before line DBID ID2a: a header that LOSAT does not read is refused
    // there, and only one with parsed Seq-ids is let through (review 1 of ID2a, M2).
    #[cfg(unix)]
    #[test]
    fn the_first_header_of_each_volume_is_read_at_open() {
        let dir = std::env::temp_dir().join(format!(
            "losat-db-badheader-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
        header[0] = 0x31;
        std::fs::write(dir.join("b.phr"), header).unwrap();
        std::fs::copy(fixture("prot_v4") + ".pin", dir.join("b.pin")).unwrap();
        std::fs::copy(fixture("prot_v4") + ".psq", dir.join("b.psq")).unwrap();
        let name = dir.join("b").to_str().unwrap().to_string();
        let settings = ApplicationSettings::default();
        for program in ["BLASTP", "PSIBLAST"] {
            let Err(error) = BlastpDatabase::open(&name, &settings, program) else {
                panic!("a damaged first header opened ({program})")
            };
            assert!(
                format!("{error:#}").ends_with(&format!(
                    "b.phr: header record, byte 1: tag 0x31 where Blast-def-line-set (0x30) is expected; this database is not supported by LOSAT's {program}"
                )),
                "{error:#}"
            );
        }
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
        let Err(error) = BlastpDatabase::open("", &settings, "BLASTP") else {
            panic!("an empty name opened")
        };
        assert_eq!(
            native(error),
            (
                2,
                "BLAST Database error: Database name is required.\n".to_string()
            )
        );
        let Err(error) = BlastpDatabase::open(&fixture("no_such_database"), &settings, "BLASTP")
        else {
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
        // A database with parsed Seq-ids opens for BLASTP and psiblast.
        for program in ["BLASTP", "PSIBLAST"] {
            BlastpDatabase::open(&fixture("prot_ids_v4"), &settings, program).unwrap_or_else(|e| {
                panic!("{program} opens a database with parsed Seq-ids: {e:#}")
            });
        }
    }

    // local_blast.cpp:240-251 and blast_app_util.hpp:185-203: a reader error of NCBI's own
    // is the search stage's text during the preliminary search (`BLAST engine error`,
    // exit 3) and the reader's text elsewhere (`BLAST Database error`, exit 2); what the
    // reader does not read is LOSAT's refusal in both places.
    #[cfg(unix)]
    #[test]
    fn reader_errors_of_the_search_and_of_the_other_reads() {
        let unreadable = || SeqDbError::Open {
            path: "/data/db.psq".into(),
            source: std::io::Error::from_raw_os_error(13),
        };
        let read = shared::search_read_error(unreadable(), "BLASTP")
            .downcast::<PreliminaryReadError>()
            .expect("NCBI's error of the search stage");
        assert_eq!(read.0, "Database memory map file error");
        assert_eq!(
            read.to_string(),
            "BLAST engine error: Database memory map file error"
        );
        let native = shared::database_error(unreadable(), "BLASTP");
        let native = native.downcast_ref::<NativeError>().expect("NCBI's error");
        assert_eq!(native.exit, 2);
        assert_eq!(
            native.message,
            "BLAST Database error: Database memory map file error - CMemoryFile: Cannot memory map file '/data/db.psq':Permission denied\n"
        );
        let too_many = || SeqDbError::Open {
            path: "/data/db.psq".into(),
            source: std::io::Error::from_raw_os_error(24),
        };
        let read = shared::search_read_error(too_many(), "BLASTP")
            .downcast::<PreliminaryReadError>()
            .expect("NCBI's error of the search stage");
        assert_eq!(
            read.0,
            "Too many open files, please raise the open file limit"
        );
        let missing = SeqDbError::NotFound("Error: File (/data/db.psq) not found.".to_string());
        let read = shared::search_read_error(missing, "BLASTP")
            .downcast::<PreliminaryReadError>()
            .expect("NCBI's error of the search stage");
        assert_eq!(read.0, "Database memory map file error");
        for refused in [shared::search_read_error, shared::database_error] {
            let error = refused(
                SeqDbError::Unsupported {
                    path: "/data/db.pal".into(),
                    message: "the alias file filters sequences".to_string(),
                },
                "BLASTP",
            );
            assert!(error.downcast_ref::<NativeError>().is_none());
            assert!(error.downcast_ref::<PreliminaryReadError>().is_none());
            assert!(
                format!("{error:#}")
                    .ends_with("; this database is not supported by LOSAT's BLASTP"),
                "{error:#}"
            );
        }
    }

    // blast_stat.c:1583-1593: O scores as X. p03 of the fixture holds every NCBIstdaa
    // letter but the gap (`ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT`).
    #[test]
    fn a_search_block_reads_o_as_x_and_keeps_the_other_letters() {
        let database = open("prot_v4");
        let range = database.db.next_block(0, u64::MAX).unwrap();
        assert_eq!((range.first_ordinal, range.len()), (0, 10));
        let stored = database.db.read_block(range).unwrap();
        let block = read_search_block(&database.db, range, "BLASTP").unwrap();
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
    fn a_search_block_with_a_foreign_byte_is_rejected_and_one_with_a_gap_residue_is_not() {
        let database = open("prot_v4");
        let range = database.db.next_block(0, u64::MAX).unwrap();
        // The gap character in sequence 3 (its third residue), and as the whole of sequence 5:
        // residues that NCBI searches (row and column 0 of the matrix).
        let mut block = database.db.read_block(range).unwrap();
        let at = block.protein_bytes_mut().len() - block.protein_bytes_mut().len(); // 0: the leading sentinel
        assert_eq!(block.protein_bytes_mut()[at], 0);
        let offset_of = |block: &SequenceBlock, i: usize| -> usize {
            (0..i).map(|j| block.protein(j).len() + 1).sum::<usize>() + 1
        };
        let position = offset_of(&block, 3) + 2;
        block.protein_bytes_mut()[position] = 0;
        let fifth = offset_of(&block, 5);
        let fifth_len = block.protein(5).len();
        block.protein_bytes_mut()[fifth..fifth + fifth_len].fill(0);
        prepare_search_block(&mut block, "BLASTP").unwrap();
        assert_eq!(block.protein(3)[2], 0);
        assert!(block.protein(5).iter().all(|&residue| residue == 0));
        // A byte that is no NCBIstdaa residue in sequence 9.
        let mut block = database.db.read_block(range).unwrap();
        let position = offset_of(&block, 9);
        block.protein_bytes_mut()[position] = 28;
        let error = prepare_search_block(&mut block, "BLASTP").unwrap_err();
        assert!(
            format!("{error:#}").starts_with(
                "database sequence 9 holds a byte that is no NCBIstdaa residue at position 1,"
            ),
            "{error:#}"
        );
        // A sentinel that is not NUL (the byte after sequence 0), alone and with a NUL
        // inside sequence 9, so that the count of NUL bytes alone would not show it.
        for compensated in [false, true] {
            let mut block = database.db.read_block(range).unwrap();
            let after_first = offset_of(&block, 1) - 1;
            block.protein_bytes_mut()[after_first] = 1;
            if compensated {
                let inside = offset_of(&block, 9) + 5;
                block.protein_bytes_mut()[inside] = 0;
            }
            let error = prepare_search_block(&mut block, "BLASTP").unwrap_err();
            assert!(
                format!("{error:#}")
                    .starts_with("database sequence 0 is not between two NUL bytes"),
                "{error:#}"
            );
        }
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
                result.x_compact_hsp_lists.push((
                    0,
                    BlastpHspList {
                        oid: s_idx as u32,
                        query_index: 0,
                        hsps: vec![test_prelim_hsp(s_idx as u32, 50 + (s_idx as i32 % 2))],
                        best_evalue: 1e-5,
                    },
                ));
            }
        }
    }

    fn test_prelim_hsp(s_idx: u32, raw_score: i32) -> XPrelimHsp {
        use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
        XPrelimHsp {
            hsp: BlastpPreliminaryHsp {
                query_start: 0,
                query_end: 10,
                subject_start: 0,
                subject_end: 10,
                gapped_query_start: 1,
                gapped_subject_start: 1,
                raw_score,
                query_context: 0,
                query_frame: 0,
                subject_frame: 0,
                query_length: 20,
                q_idx: 0,
                s_idx,
            },
            e_value: 1e-5,
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
        let mut hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>> = vec![None];
        crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            search_blocks_of(
                database,
                pool,
                window,
                block_bytes,
                2,
                &|| BlastpSubjectScratch::new(16, 64, window, 21),
                &recording_kernel(window, &seen),
                &mut crate::algorithm::blastp::blast_engine::PrelimCollector {
                    x_compact_hit_lists: &mut hit_lists,
                    prelim_hitlist_size: 500,
                },
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

    /// EXPERIMENT (LOSAT_X_BPDBCUT): the blocks cut by the database size and the thread count
    /// (`database/x_db_cut.rs`) against the blocks of `BLOCK_BYTES`.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    mod x_db_cut_blocks {
        use super::super::x_db_cut::{
            x_count_blocks, x_cut_block_bytes, x_db_cut_block_bytes, x_search_blocks_cut,
            x_sequence_bytes,
        };
        use super::*;
        use std::path::{Path, PathBuf};

        type Seen = Vec<(usize, i32, Vec<u8>)>;

        /// What the recording kernel saw (sorted by ordinal id) and the listed ordinal ids of a
        /// search of `database` at `threads` threads in rounds of `blocks_per_thread` blocks per
        /// thread: `search_blocks_of` with blocks of `block_bytes` (`Some`), or the cut's driver
        /// (`None`).
        fn search(
            database: &BlastpDatabase,
            threads: usize,
            blocks_per_thread: usize,
            block_bytes: Option<u64>,
        ) -> Result<(Seen, Vec<u32>)> {
            let driver = match block_bytes {
                None => Driver::Cut,
                Some(block_bytes) => Driver::Of(block_bytes),
            };
            search_with(database, threads, blocks_per_thread, driver)
        }

        /// The drivers of the preliminary stage that the tests call.
        #[derive(Clone, Copy, Debug)]
        enum Driver {
            /// `x_search_blocks_cut`.
            Cut,
            /// `search_blocks_of` with blocks of this many bytes.
            Of(u64),
            /// `x_db_keep::x_search_blocks_kept` with blocks of this many bytes.
            Kept(u64),
            /// `search_blocks_made` (without LOSAT_X_BPDBKEEP) with blocks of this many bytes.
            Made(u64),
            /// `Kept` with the kernel's hook.
            KeptWith(u64, XHook),
            /// `x_db_keep::x_search_window` with blocks of this many bytes, the window and the hook.
            Window(u64, super::super::x_db_keep::XWindowSize, Option<XHook>),
        }

        /// What the kernel does before it records a sequence (`Driver::KeptWith`).
        #[derive(Clone, Copy, Debug)]
        enum XHook {
            /// Sleeps 50 ms at this ordinal id.
            SleepAt(usize),
            /// Panics at this ordinal id.
            PanicAt(usize),
        }

        /// `kernel` after the hook `hook`.
        fn hooked_kernel<'a>(
            kernel: impl Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult)
                + Sync
                + 'a,
            hook: Option<XHook>,
        ) -> impl Fn(&mut BlastpSubjectScratch, usize, BlastpSubjectRef<'_>, &mut BlastpSubjectResult)
               + Sync
               + 'a {
            move |scratch, s_idx, subject, result| {
                match hook {
                    Some(XHook::SleepAt(ordinal)) if ordinal == s_idx => {
                        std::thread::sleep(std::time::Duration::from_millis(50));
                    }
                    Some(XHook::PanicAt(ordinal)) if ordinal == s_idx => {
                        panic!("the kernel's hook at {ordinal}");
                    }
                    _ => {}
                }
                kernel(scratch, s_idx, subject, result);
            }
        }

        /// `search` with the driver `driver`.
        fn search_with(
            database: &BlastpDatabase,
            threads: usize,
            blocks_per_thread: usize,
            driver: Driver,
        ) -> Result<(Seen, Vec<u32>)> {
            let window = 40;
            let seen = std::sync::Mutex::new(Vec::new());
            let mut hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>> = vec![None];
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                let new_scratch = || BlastpSubjectScratch::new(16, 64, window, 21);
                let hook = match driver {
                    Driver::KeptWith(_, hook) => Some(hook),
                    Driver::Window(_, _, hook) => hook,
                    _ => None,
                };
                let kernel = hooked_kernel(recording_kernel(window, &seen), hook);
                let mut sink = crate::algorithm::blastp::blast_engine::PrelimCollector {
                    x_compact_hit_lists: &mut hit_lists,
                    prelim_hitlist_size: 500,
                };
                match driver {
                    Driver::Cut => x_search_blocks_cut(
                        database,
                        pool,
                        window,
                        blocks_per_thread,
                        &new_scratch,
                        &kernel,
                        &mut sink,
                    ),
                    Driver::Of(block_bytes) => search_blocks_of(
                        database,
                        pool,
                        window,
                        block_bytes,
                        blocks_per_thread,
                        &new_scratch,
                        &kernel,
                        &mut sink,
                    ),
                    Driver::Kept(block_bytes) | Driver::KeptWith(block_bytes, _) => {
                        super::super::x_db_keep::x_search_blocks_kept(
                            database,
                            pool,
                            window,
                            block_bytes,
                            blocks_per_thread,
                            &new_scratch,
                            &kernel,
                            &mut sink,
                        )
                    }
                    Driver::Made(block_bytes) => search_blocks_made(
                        database,
                        pool,
                        window,
                        block_bytes,
                        blocks_per_thread,
                        &new_scratch,
                        &kernel,
                        &mut sink,
                    ),
                    Driver::Window(block_bytes, size, _) => {
                        super::super::x_db_keep::x_search_window(
                            &database,
                            pool,
                            window,
                            block_bytes,
                            size,
                            &new_scratch,
                            &kernel,
                            &mut sink,
                        )
                    }
                }
            })?;
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
            Ok((seen, listed))
        }

        /// The cut of `database` on a pool of `threads` threads.
        fn cut_of(database: &BlastpDatabase, threads: usize) -> u64 {
            crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                Ok(x_db_cut_block_bytes(database, pool))
            })
            .unwrap()
        }

        /// The letters of the generated database: 300 sequences of 1 to 1,500 letters from a fixed
        /// generator, with sequences of 40,000 and 70,000 letters at ordinal ids 100 and 250.
        fn generated_lengths() -> Vec<usize> {
            let mut state = 0x2545_f491_4f6c_dd1du64;
            (0..300)
                .map(|ordinal| match ordinal {
                    100 => 40_000,
                    250 => 70_000,
                    _ => {
                        state = state
                            .wrapping_mul(6_364_136_223_846_793_005)
                            .wrapping_add(1_442_695_040_888_963_407);
                        1 + (state >> 33) as usize % 1_500
                    }
                })
                .collect()
        }

        /// The generated database written by `losat makeblastdb` into a new directory of the test
        /// (`label`), with volumes of at most `max_file_sz` bytes; returns the directory and the
        /// database name.
        fn generated_database(label: &str, max_file_sz: &str) -> (PathBuf, String) {
            const LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
            let mut fasta = Vec::new();
            let mut state = 7u64;
            for (ordinal, length) in generated_lengths().into_iter().enumerate() {
                fasta.extend_from_slice(format!(">s{ordinal} generated\n").as_bytes());
                for i in 0..length {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    fasta.push(LETTERS[(state >> 33) as usize % LETTERS.len()]);
                    if i % 80 == 79 || i + 1 == length {
                        fasta.push(b'\n');
                    }
                }
            }
            let dir =
                std::env::temp_dir().join(format!("losat-bpdbcut-{label}-{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            let name = dir.join("gen").to_str().unwrap().to_string();
            let args = crate::writedb::makeblastdb::MakedbArgs {
                input: "-".to_string(),
                input_type: "fasta".to_string(),
                dbtype: "prot".to_string(),
                title: Some("gen".to_string()),
                out: Some(name.clone()),
                blastdb_version: 4,
                max_file_sz: max_file_sz.to_string(),
                logfile: None,
                parse_seqids: false,
                hash_index: false,
                taxid: None,
                taxid_map: None,
            };
            let now = crate::writedb::time::LocalTime {
                year: 2026,
                month: 10,
                day: 16,
                hour: 13,
                minute: 30,
                second: 0,
            };
            let (mut log, mut messages) = (Vec::new(), Vec::new());
            crate::writedb::makeblastdb::build_database(
                &args,
                Some(&now),
                crate::writedb::Platform::host(),
                &mut fasta.as_slice(),
                &mut log,
                &mut messages,
            )
            .unwrap_or_else(|error| panic!("makeblastdb: {error}"));
            (dir, name)
        }

        /// EXPERIMENT (LOSAT_X_BPDBKEEP): `cut` with the switch's largest block.
        fn x_keep_cut(cut: u64) -> u64 {
            if crate::seqdb::x_bp_db_keep() {
                cut.min(super::super::x_db_keep::X_KEEP_BLOCK_BYTES)
            } else {
                cut
            }
        }

        fn open_path(name: &str) -> BlastpDatabase {
            BlastpDatabase::open(name, &ApplicationSettings::default(), "BLASTP")
                .unwrap_or_else(|error| panic!("{name}: {error:#}"))
        }

        // The cut of the fixture databases (812 bytes of sequence, one volume or six): one thread
        // keeps 1 MiB, more threads the 4 KiB floor, and every volume is one block either way. The
        // sequence bytes are those of the one-volume sequence file.
        #[test]
        fn the_cut_of_tiny_and_multi_volume_databases() {
            let single = open("prot_v4");
            let volumes = open("prot_mv");
            let psq_bytes = |names: &[&str]| -> u64 {
                names
                    .iter()
                    .map(|name| std::fs::metadata(fixture(name)).unwrap().len())
                    .sum()
            };
            assert_eq!(x_sequence_bytes(&single), psq_bytes(&["prot_v4.psq"]));
            assert_eq!(x_sequence_bytes(&single), 812);
            let mv: Vec<String> = (0..6).map(|v| format!("prot_mv.{v:02}.psq")).collect();
            let mv: Vec<&str> = mv.iter().map(String::as_str).collect();
            // Each further volume has its own leading sentinel byte.
            assert_eq!(x_sequence_bytes(&volumes) + 5, psq_bytes(&mv));
            assert_eq!(x_sequence_bytes(&volumes), x_sequence_bytes(&single));
            assert_eq!(cut_of(&single, 1), x_keep_cut(BLOCK_BYTES));
            assert_eq!(cut_of(&volumes, 1), x_keep_cut(BLOCK_BYTES));
            for threads in [2, 3, 4, 8] {
                assert_eq!(cut_of(&single, threads), x_cut_block_bytes(812, threads));
                assert_eq!(cut_of(&single, threads), 4096);
                assert_eq!(cut_of(&volumes, threads), 4096);
            }
            for block_bytes in [BLOCK_BYTES, 4096] {
                assert_eq!(x_count_blocks(&single.db, block_bytes), 1);
                assert_eq!(x_count_blocks(&volumes.db, block_bytes), 6);
            }
        }

        // blast_engine.c:1411-1554 and blast_extend.c:160-172: with the cut, the kernel gets every
        // sequence once with the one-thread order's diagonal offset, and the hit list gets the
        // one-thread order's updates, at 1, 2, 3 and 7 threads, in rounds of 2 and of
        // `BLOCKS_PER_THREAD` blocks per thread, on one volume and on several; the sequences of
        // 40,000 and 70,000 letters are larger than the cut at 2 threads and more (blocks of their
        // own).
        #[test]
        fn blocks_of_the_cut_give_the_one_thread_search() {
            let (single_dir, single_name) = generated_database("one", "3GB");
            let (multi_dir, multi_name) = generated_database("multi", "100000");
            let single = open_path(&single_name);
            let multi = open_path(&multi_name);
            let lengths = generated_lengths();
            let letters: usize = lengths.iter().sum();
            assert_eq!(single.total_letters(), letters as u64);
            assert_eq!(multi.total_letters(), letters as u64);
            assert_eq!(single.db.volumes().len(), 1);
            assert!(multi.db.volumes().len() > 2, "{}", multi.db.volumes().len());
            let (reference, reference_listed) = search(&single, 1, 2, Some(BLOCK_BYTES)).unwrap();
            assert_eq!(reference.len(), lengths.len());
            assert_eq!(reference_listed.len(), 100);
            assert_eq!(x_count_blocks(&single.db, BLOCK_BYTES), 1);
            for threads in [1usize, 2, 3, 7] {
                let cut = cut_of(&single, threads);
                if threads == 1 {
                    assert_eq!(cut, x_keep_cut(BLOCK_BYTES));
                } else {
                    assert_eq!(
                        cut,
                        x_cut_block_bytes(letters as u64 + lengths.len() as u64 + 1, threads)
                    );
                    // Both long sequences at 3 threads and more, the longer one at 2.
                    assert!(cut < 70_000 + 2, "{threads} threads: {cut}");
                    assert!(
                        x_count_blocks(&single.db, cut) > 2 * threads as u64,
                        "{threads} threads: {} blocks",
                        x_count_blocks(&single.db, cut)
                    );
                }
                for blocks_per_thread in [2, BLOCKS_PER_THREAD] {
                    for (what, db) in [("one volume", &single), ("volumes", &multi)] {
                        let (seen, listed) = search(db, threads, blocks_per_thread, None).unwrap();
                        let context = format!(
                            "{what}, {threads} threads, {blocks_per_thread} blocks per thread"
                        );
                        assert!(seen == reference, "{context}");
                        assert_eq!(listed, reference_listed, "{context}");
                    }
                }
            }
            let _ = std::fs::remove_dir_all(single_dir);
            let _ = std::fs::remove_dir_all(multi_dir);
        }

        // No NCBI counterpart: LOSAT_X_BPDBKEEP keeps a scratch and a block buffer per worker
        // (`x_db_keep.rs`); the kernel sees the same sequences from the same diagonal offsets, and
        // the hit lists receive the same updates, as with a scratch and a buffer per task.
        // At 1, 2, 3 and 7 threads, blocks of 1 MiB and of 4 KiB (the 40,000 and 70,000-letter
        // sequences in blocks of their own, so buffers grow and are refilled shorter), in rounds of 2 and of
        // `BLOCKS_PER_THREAD` blocks per thread, on one volume and on several.
        #[test]
        fn kept_scratches_and_buffers_give_the_search_without_them() {
            let (single_dir, single_name) = generated_database("keep1", "3GB");
            let (multi_dir, multi_name) = generated_database("keepn", "100000");
            let single = open_path(&single_name);
            let multi = open_path(&multi_name);
            let (reference, reference_listed) =
                search_with(&single, 1, 2, Driver::Made(BLOCK_BYTES)).unwrap();
            assert_eq!(reference.len(), generated_lengths().len());
            for threads in [1usize, 2, 3, 7] {
                for block_bytes in [BLOCK_BYTES, 4096] {
                    for blocks_per_thread in [2, BLOCKS_PER_THREAD] {
                        for (what, db) in [("one volume", &single), ("volumes", &multi)] {
                            let context = format!(
                                "{what}, {threads} threads, {block_bytes}-byte blocks, {blocks_per_thread} blocks per thread"
                            );
                            let made = search_with(
                                db,
                                threads,
                                blocks_per_thread,
                                Driver::Made(block_bytes),
                            )
                            .unwrap();
                            let kept = search_with(
                                db,
                                threads,
                                blocks_per_thread,
                                Driver::Kept(block_bytes),
                            )
                            .unwrap();
                            assert!(made.0 == reference, "{context}");
                            assert_eq!(made.1, reference_listed, "{context}");
                            assert!(kept.0 == made.0, "{context}");
                            assert_eq!(kept.1, made.1, "{context}");
                        }
                    }
                }
            }
            let _ = std::fs::remove_dir_all(single_dir);
            let _ = std::fs::remove_dir_all(multi_dir);
        }

        /// The windows that wait on the database of the tests: of 2 blocks per thread, and wider
        /// with bytes that narrow it to those 2 per thread as soon as any result waits.
        fn waiting_windows(threads: usize) -> [super::super::x_db_keep::XWindowSize; 2] {
            let floor = 2 * threads;
            [
                super::super::x_db_keep::XWindowSize {
                    floor,
                    width: floor,
                    max_waiting_bytes: usize::MAX,
                },
                super::super::x_db_keep::XWindowSize {
                    floor,
                    width: 64 * threads,
                    max_waiting_bytes: 1,
                },
            ]
        }

        // No NCBI counterpart: the window of LOSAT_X_BPDBKEEP and LOSAT_X_BPDBWIDE (`x_db_keep.rs`)
        // with a slow block, the oldest one (the later blocks fill the window and their takers wait
        // for its merge) or a later one: the kernel sees the same sequences from the same diagonal
        // offsets and the hit lists receive the same updates as with a scratch and a buffer per
        // task. A panic of the kernel in either block ends the search with the panic (the waiting
        // workers stop). At 2, 3 and 7 threads, blocks of 4 KiB, windows that wait
        // (`waiting_windows`; the default window is wider than this database).
        #[test]
        fn the_window_waits_for_its_oldest_block_and_ends_on_a_panic() {
            let (dir, name) = generated_database("keepw", "3GB");
            let database = open_path(&name);
            let block_bytes = 4096;
            assert!(x_count_blocks(&database.db, block_bytes) > 40);
            for threads in [2usize, 3, 7] {
                let made = search_with(&database, threads, 2, Driver::Made(block_bytes)).unwrap();
                for size in waiting_windows(threads) {
                    for ordinal in [0usize, 150] {
                        let context = format!("{threads} threads, ordinal {ordinal}, {size:?}");
                        let slow = search_with(
                            &database,
                            threads,
                            2,
                            Driver::Window(block_bytes, size, Some(XHook::SleepAt(ordinal))),
                        )
                        .unwrap();
                        assert!(slow.0 == made.0, "{context}");
                        assert_eq!(slow.1, made.1, "{context}");
                        let panicked =
                            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                search_with(
                                    &database,
                                    threads,
                                    2,
                                    Driver::Window(
                                        block_bytes,
                                        size,
                                        Some(XHook::PanicAt(ordinal)),
                                    ),
                                )
                            }));
                        assert!(panicked.is_err(), "{context}");
                    }
                }
            }
            let _ = std::fs::remove_dir_all(dir);
        }

        fn error_text(result: Result<(Seen, Vec<u32>)>) -> String {
            match result {
                Ok(_) => panic!("the search of the damaged database succeeded"),
                Err(error) => format!("{error:#}"),
            }
        }

        // A damaged database fails as with one thread at every thread count: a sequence file 500
        // bytes shorter than its index says (the byte range of the read names the block), and
        // then also a byte 28 in sequence 100 (a refusal of the search block that a smaller block
        // meets before the read that fails, while the one-thread order's 1 MiB block fails to read
        // first).
        #[test]
        fn a_damaged_database_gives_the_one_thread_error_with_the_cut() {
            let (dir, name) = generated_database("damaged", "3GB");
            let psq = Path::new(&name).with_extension("psq");
            let sequence_100 = open_path(&name).db.volumes()[0]
                .index()
                .sequence_start(100)
                .unwrap() as usize;
            let mut bytes = std::fs::read(&psq).unwrap();
            bytes.truncate(bytes.len() - 500);
            std::fs::write(&psq, &bytes).unwrap();
            for foreign in [false, true] {
                if foreign {
                    bytes[sequence_100 + 1_000] = BLASTAA_SIZE as u8;
                    std::fs::write(&psq, &bytes).unwrap();
                }
                let database = open_path(&name);
                let expected = error_text(search(&database, 1, 2, Some(BLOCK_BYTES)));
                assert!(
                    expected.contains("passes the end of the file"),
                    "{expected}"
                );
                assert_eq!(
                    error_text(search(&database, 4, 2, Some(BLOCK_BYTES))),
                    expected
                );
                for threads in [2usize, 4, 7] {
                    let cut = cut_of(&database, threads);
                    // Without the one-thread order's error, the cut's own error differs.
                    let own = error_text(search(&database, threads, 2, Some(cut)));
                    assert_ne!(own, expected, "{threads} threads");
                    // EXPERIMENT (LOSAT_X_BPDBKEEP): the kept buffers give the same errors.
                    for block_bytes in [BLOCK_BYTES, cut] {
                        assert_eq!(
                            error_text(search_with(
                                &database,
                                threads,
                                2,
                                Driver::Kept(block_bytes)
                            )),
                            error_text(search_with(
                                &database,
                                threads,
                                2,
                                Driver::Made(block_bytes)
                            )),
                            "{threads} threads, {block_bytes}-byte blocks"
                        );
                        // Windows whose workers wait when the failed block stops the search.
                        for size in waiting_windows(threads) {
                            assert_eq!(
                                error_text(search_with(
                                    &database,
                                    threads,
                                    2,
                                    Driver::Window(block_bytes, size, None)
                                )),
                                error_text(search_with(
                                    &database,
                                    threads,
                                    2,
                                    Driver::Made(block_bytes)
                                )),
                                "{threads} threads, {block_bytes}-byte blocks, {size:?}"
                            );
                        }
                    }
                    if foreign {
                        assert!(own.contains("is no NCBIstdaa residue"), "{own}");
                    }
                    assert_eq!(
                        error_text(search(&database, threads, 2, None)),
                        expected,
                        "{threads} threads"
                    );
                }
            }
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    // `x_list` lists the subjects and renames them to their ranks with the stored lengths and no
    // copy; the redo reads the stored residues by rank; after the redo `x_read_shown` keeps the
    // residues that the reports read (p03 and p08 hold an O).
    #[test]
    fn listed_subjects_read_the_stored_residues_by_rank() {
        let database = open("prot_mv");
        let list = |oid: u32, raw_score: i32| BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![test_hsp(oid, raw_score), test_hsp(oid, raw_score - 1)],
            best_evalue: 1e-5,
        };
        let make = || {
            let mut first = BlastpHitList::new(10);
            for (oid, score) in [(7, 60), (2, 55), (9, 50)] {
                first.update(list(oid, score));
            }
            let mut second = BlastpHitList::new(10);
            second.update(list(2, 70));
            second.update(list(0, 65));
            vec![Some(first), None, Some(second)]
        };
        let mut lists = make();
        let listed = ListedSubjects::x_list(&database, &mut lists, &mut []).unwrap();
        assert_eq!(listed.ordinals, [0, 2, 7, 9]);
        assert!(listed.encoded.is_empty());
        // The subjects are renamed to their ranks: ordinals 7, 2, 9 and 2, 0 are ranks 2, 1, 3
        // and 1, 0.
        let ranks = |lists: &[Option<BlastpHitList>]| -> Vec<Vec<u32>> {
            lists
                .iter()
                .flatten()
                .map(|hit_list| hit_list.hsplist_array.iter().map(|l| l.oid).collect())
                .collect()
        };
        assert_eq!(ranks(&lists), [vec![2, 1, 3], vec![1, 0]]);
        // The stored residues of each rank (between the sentinels), straight from the database.
        let stored: Vec<Vec<u8>> = listed
            .ordinals
            .iter()
            .map(|&ordinal| {
                database
                    .db
                    .read_sequence(ordinal)
                    .unwrap()
                    .protein(0)
                    .to_vec()
            })
            .collect();
        let lengths: Vec<usize> = stored.iter().map(Vec::len).collect();
        assert_eq!(lengths[3], 300);
        let subjects = listed.redo_subjects(&database);
        assert_eq!(subjects.len(), 4);
        for (rank, &length) in lengths.iter().enumerate() {
            assert_eq!(subjects.aa_len(rank), Some(length));
            let read = subjects.read(rank as u32).unwrap();
            assert_eq!(read.residues(), &stored[rank][..]);
        }
        assert_eq!(subjects.aa_len(4), None);
        // A final hit list of ranks 1 (p03, with an O) and 3 (without); rank 2 (p08) holds an O but
        // is in no final list.
        let mut final_list = BlastpHitList::new(10);
        final_list.update(list(1, 60));
        final_list.update(list(3, 50));
        let final_lists = vec![None, Some(final_list)];
        // With `WithO`, the HSPs of a query with an O (query 0 here) keep every subject's residues.
        for (residues, query_has_o, kept) in [
            (XShownResidues::WithO, &[][..], [false, true, false, false]),
            (
                XShownResidues::WithO,
                &[false, true],
                [false, true, false, false],
            ),
            (XShownResidues::WithO, &[true], [false, true, false, true]),
            (XShownResidues::Final, &[][..], [false, true, false, true]),
            (XShownResidues::Listed, &[][..], [true; 4]),
        ] {
            let mut shown = ListedSubjects::x_list(&database, &mut make(), &mut []).unwrap();
            // The redo reads the subjects of the final lists (and finds the O of p03).
            for rank in [1, 3] {
                shown.redo_subjects(&database).read(rank).unwrap();
            }
            shown
                .x_read_shown(&database, &final_lists, residues, query_has_o)
                .unwrap();
            assert_eq!(shown.encoded.len(), 4, "{residues:?}");
            for (rank, kept) in kept.into_iter().enumerate() {
                assert_eq!(
                    shown.encoded[rank].aa_len, lengths[rank],
                    "{residues:?} {rank}"
                );
                if kept {
                    // Between two sentinel bytes.
                    assert_eq!(
                        shown.encoded[rank].aa_seq.len(),
                        lengths[rank] + 2,
                        "{residues:?} {rank}"
                    );
                    assert_eq!(
                        &shown.encoded[rank].aa_seq[1..=lengths[rank]],
                        &stored[rank][..],
                        "{residues:?} {rank}"
                    );
                } else {
                    assert!(shown.encoded[rank].aa_seq.is_empty(), "{residues:?} {rank}");
                }
            }
        }
        // `read` (the unit tests' listed subjects) is `Listed`.
        let read = ListedSubjects::read(&database, &mut make(), &mut []).unwrap();
        assert_eq!(read.encoded[2].aa_len, lengths[2]);
        assert_eq!(read.encoded[3].aa_seq.len(), 302);
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
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
        assert_eq!(listed.ordinals, [0, 2, 7, 9]);
        let titles = listed.read_titles(&database, &hit_lists).unwrap();
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
            titles.records()[1].title(),
            b"p03 every letter of NCBIstdaa but the gap"
        );
        assert_eq!(
            titles.records()[1].seq(),
            b"ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT"
        );
        let redone: String = listed.encoded[1].aa_seq[1..=listed.encoded[1].aa_len]
            .iter()
            .map(|&residue| ncbistdaa_to_ascii(residue))
            .collect();
        assert_eq!(redone, "ABCDEFGHIKLMNPQRSTVWXYZUOJMHCENQT");
        assert_eq!(listed.encoded[1].aa_seq.len(), 35);
        assert_eq!(
            titles.records()[2].seq(),
            b"MUGYDNEDMDNEPCNDCRACIGFRYKYMMVKWOFASVIHWKHVWDJ"
        );
        assert_eq!(listed.encoded[3].aa_len, 300);
        // The tabular ids of the subjects with rows, the stored titles for `stitle`.
        let (ids, accessions) = titles.tabular_subject_ids(&hit_lists).unwrap();
        assert!(accessions.is_empty());
        assert_eq!(
            ids.iter().map(|id| id.as_ref()).collect::<Vec<_>>(),
            [&b"p01"[..], b"p03", b"p08", b"p10"]
        );
        assert_eq!(
            titles.stitles()[0].as_ref(),
            b"p01 first protein [Test organism]"
        );
        titles
            .check_shown_titles(&database, &[0usize, 1, 2, 3].into_iter().collect())
            .unwrap();
        titles.check_reported_titles(&hit_lists).unwrap();
    }

    // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles are the records' titles, rank by rank,
    // with one shared empty title for the subjects in no final hit list, and the stored residues
    // give the records' letters (`x_listed_lean.rs`).
    #[test]
    fn lean_titles_are_the_records_titles_and_letters() {
        let database = open("prot_mv");
        let list = |oid: u32, raw_score: i32| BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![test_hsp(oid, raw_score)],
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
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
        // Ordinal ids 7 and 9 (ranks 2 and 3) stay listed but leave the final hit lists.
        let first = hit_lists[0].as_mut().unwrap();
        first.hsplist_array.retain(|hsp_list| hsp_list.oid == 1);
        first.hsplist_count = first.hsplist_array.len();
        let records = listed.read_titles(&database, &hit_lists).unwrap();
        let lean = listed.x_read_titles_lean(&database, &hit_lists).unwrap();
        let ListedShown::Lean(lean_titles) = &lean.shown else {
            panic!("lean titles");
        };
        assert!(matches!(records.shown, ListedShown::Records(_)));
        assert_eq!(lean_titles.len(), 4);
        for rank in 0..4 {
            assert_eq!(
                lean.title(rank),
                records.records()[rank].title(),
                "rank {rank}"
            );
            assert_eq!(
                lean.replaced_title_bytes[rank],
                records.replaced_title_bytes[rank]
            );
        }
        assert_eq!(lean.title(1), b"p03 every letter of NCBIstdaa but the gap");
        assert!(lean.title(2).is_empty() && lean.title(3).is_empty());
        assert!(Arc::ptr_eq(&lean_titles[2], &lean_titles[3]));
        assert_eq!(
            lean.tabular_subject_ids(&hit_lists).unwrap(),
            records.tabular_subject_ids(&hit_lists).unwrap()
        );
        lean.check_shown_titles(&database, &[0usize, 1].into_iter().collect())
            .unwrap();
        lean.check_reported_titles(&hit_lists).unwrap();
        assert_eq!(lean.stitles(), records.stitles());
        // Titles, flags and every letter against the stored residues.
        x_assert_same_titles_and_letters(&lean, &records);
    }

    // blast_seqalign.cpp:1569-1582 and seqdbvol.cpp:2302-2310: NCBI opens a volume's
    // header file when the Seq-aligns of the final hit lists are made, query by query and
    // in the order of each hit list; a subject that is in no final hit list is never read.
    // prot_mv: volume 00 holds the ordinal ids 0-2, 01 holds 3, 05 holds 9.
    // (A path in the system's temporary directory: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn titles_are_read_for_the_final_hit_lists_in_ncbis_order() {
        let dir = std::env::temp_dir().join(format!(
            "losat-db-nohdr-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        for entry in
            std::fs::read_dir(fixture("prot_mv.pal").trim_end_matches("prot_mv.pal")).unwrap()
        {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_string();
            if name.starts_with("prot_mv.") && name != "prot_mv.00.phr" && name != "prot_mv.05.phr"
            {
                std::fs::copy(&path, dir.join(&name)).unwrap();
            }
        }
        let name = dir.join("prot_mv").to_str().unwrap().to_string();
        let database =
            BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
        let list = |oid: u32, best_evalue: f64| BlastpHspList {
            oid,
            query_index: 0,
            hsps: vec![test_hsp(oid, 60)],
            best_evalue,
        };
        let hit_list = |lists: Vec<BlastpHspList>, count: usize| {
            let mut hit_list = BlastpHitList::new(10);
            hit_list.hsplist_array = lists;
            hit_list.hsplist_count = count;
            Some(hit_list)
        };
        // Query 1 reports ordinal id 3 (volume 01); ordinal id 1 (volume 00, no header
        // file) is in its list past the count, and is not read.
        let mut hit_lists = vec![None, hit_list(vec![list(3, 1e-9), list(1, 1e-3)], 1)];
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
        assert_eq!(listed.ordinals, [1, 3]);
        let titles = listed.read_titles(&database, &hit_lists).unwrap();
        assert_eq!(
            titles.records()[1].title(),
            b"p04 lower case and ambiguity codes"
        );
        assert_eq!(titles.records()[0].title(), b"");
        assert_eq!(titles.records()[0].seq(), b"I");
        // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles do not read it either.
        let lean = listed.x_read_titles_lean(&database, &hit_lists).unwrap();
        assert_eq!(lean.title(1), b"p04 lower case and ambiguity codes");
        assert_eq!(lean.title(0), b"");
        // Query 2 reports ordinal id 9 (volume 05) before ordinal id 0 (volume 00): the
        // first header file read is that of volume 05, in query 2's batch.
        let mut hit_lists = vec![
            None,
            hit_list(vec![list(3, 1e-9), list(1, 1e-3)], 1),
            hit_list(vec![list(9, 1e-9), list(0, 1e-5)], 2),
        ];
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
        assert_eq!(listed.ordinals, [0, 1, 3, 9]);
        let error = listed
            .read_titles(&database, &hit_lists)
            .err()
            .expect("a missing header file")
            .downcast::<HeaderReadError>()
            .expect("NCBI's error, with the query that reached it");
        assert_eq!(error.query, 2);
        let native = error
            .error
            .downcast_ref::<NativeError>()
            .expect("NCBI's error");
        assert_eq!(native.exit, 2);
        // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles are read in the same order, with
        // the same error.
        let lean_error = listed
            .x_read_titles_lean(&database, &hit_lists)
            .err()
            .expect("a missing header file")
            .downcast::<HeaderReadError>()
            .expect("NCBI's error, with the query that reached it");
        assert_eq!(lean_error.query, 2);
        assert_eq!(
            lean_error
                .error
                .downcast_ref::<NativeError>()
                .unwrap()
                .message,
            native.message
        );
        // The volume is recorded under its directory normalised with symbolic links
        // followed (seqdbalias.cpp:851-862, `NormalizePath(dir_name, eFollowLinks)`), so
        // the message names that directory: on macOS the temporary directory under `/var`
        // is shown under `/private/var`.
        assert_eq!(
            native.message,
            format!(
                "BLAST Database error: Error: File ({}) not found.\n",
                dir.canonicalize().unwrap().join("prot_mv.05.phr").display()
            )
        );
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
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
        let database =
            BlastpDatabase::open(&name, &ApplicationSettings::default(), "BLASTP").unwrap();
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
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
        let titles = listed.read_titles(&database, &hit_lists).unwrap();
        titles.check_reported_titles(&hit_lists).unwrap();
        // The subject with the bad byte: read as `#`, and rejected when it is reported.
        let mut bad = BlastpHitList::new(10);
        bad.update(list(9));
        bad.update(list(0));
        let mut hit_lists = vec![Some(bad)];
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut []).unwrap();
        let titles = listed.read_titles(&database, &hit_lists).unwrap();
        assert_eq!(titles.records()[1].title(), b"p10 l#st protein");
        let error = titles.check_reported_titles(&hit_lists).unwrap_err();
        assert!(
            format!("{error:#}").starts_with(
                "database sequence 9 has a stored title with a byte outside ' ' to '~'"
            ),
            "{error:#}"
        );
        // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the lean titles hold the same title and flag, and
        // give the same refusal.
        let lean = listed.x_read_titles_lean(&database, &hit_lists).unwrap();
        assert_eq!(lean.title(1), b"p10 l#st protein");
        let lean_error = lean.check_reported_titles(&hit_lists).unwrap_err();
        assert_eq!(format!("{lean_error:#}"), format!("{error:#}"));
        drop(database);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // The compact preliminary lists (`x_prelim_compact.rs`) are listed and renamed to ranks with
    // the lists of full HSPs, before they are expanded.
    #[test]
    fn compact_lists_are_listed_and_renamed_with_the_others() {
        use crate::algorithm::blastp::gapalign::BlastpPreliminaryHsp;
        let database = open("prot_mv");
        let compact = |oid: u32| BlastpHspList {
            oid,
            query_index: 1,
            hsps: vec![XPrelimHsp {
                hsp: BlastpPreliminaryHsp {
                    query_start: 0,
                    query_end: 10,
                    subject_start: 0,
                    subject_end: 10,
                    gapped_query_start: 1,
                    gapped_subject_start: 1,
                    raw_score: 60,
                    query_context: 1,
                    query_frame: 0,
                    subject_frame: 0,
                    query_length: 20,
                    q_idx: 1,
                    s_idx: oid,
                },
                e_value: 1e-5,
            }],
            best_evalue: 1e-5,
        };
        let mut off = BlastpHitList::new(10);
        off.update(BlastpHspList {
            oid: 7,
            query_index: 0,
            hsps: vec![test_hsp(7, 60)],
            best_evalue: 1e-5,
        });
        let mut compact_list = BlastpHitList::new(10);
        compact_list.update(compact(9));
        compact_list.update(compact(2));
        let mut hit_lists = vec![Some(off), None];
        let mut compact_lists = vec![None, Some(compact_list)];
        let listed = ListedSubjects::read(&database, &mut hit_lists, &mut compact_lists).unwrap();
        assert_eq!(listed.ordinals, [2, 7, 9]);
        assert_eq!(listed.encoded.len(), 3);
        let off = &hit_lists[0].as_ref().unwrap().hsplist_array[0];
        assert_eq!((off.oid, off.hsps[0].s_idx), (1, 1));
        let mut renamed: Vec<(u32, u32)> = compact_lists[1]
            .as_ref()
            .unwrap()
            .hsplist_array
            .iter()
            .map(|list| (list.oid, list.hsps[0].hsp.s_idx))
            .collect();
        renamed.sort();
        assert_eq!(renamed, [(0, 0), (2, 2)]);
    }

    // EXPERIMENT (LOSAT_X_BPQBLOCK): the one pass over the headers (`SeqDb::visit_deflines_taxids`)
    // gives the answer of the two checks it replaces (`x_titles_read_allow_blocks`, then
    // `x_taxids_allow_blocks` one sequence at a time), for a database whose headers read, one with
    // taxids, an alias of volumes, and one whose header file is missing, at one thread and in a pool.
    #[test]
    fn titles_and_taxids_in_one_pass_give_the_answers_of_the_two_checks() {
        let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let paja = manifest_dir.join("tests/db_search/paja_taxid");
        let dir = std::env::temp_dir().join(format!(
            "losat-db-qblock-onepass-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        for ext in ["pin", "psq"] {
            std::fs::copy(paja.with_extension(ext), dir.join(format!("nohdr.{ext}"))).unwrap();
        }
        let mut databases = vec![
            (open("prot_v4"), true),
            (open("prot_mv"), true),
            (
                BlastpDatabase::open(
                    paja.to_str().unwrap(),
                    &ApplicationSettings::default(),
                    "BLASTP",
                )
                .unwrap(),
                true,
            ),
        ];
        // A database whose header file is missing: the search of all queries reports NCBI's error.
        if let Ok(missing) = BlastpDatabase::open(
            dir.join("nohdr").to_str().unwrap(),
            &ApplicationSettings::default(),
            "BLASTP",
        ) {
            databases.push((missing, false));
        }
        for (database, readable) in &databases {
            for threads in [1, 4] {
                crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                    for names in [true, false] {
                        let two = x_titles_read_allow_blocks(database, pool)
                            && x_taxids_allow_blocks(database, names, pool);
                        assert_eq!(two, *readable);
                        assert_eq!(
                            x_titles_taxids_allow_blocks(database, names, pool),
                            two,
                            "threads {threads}"
                        );
                    }
                    Ok(())
                })
                .unwrap();
            }
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // EXPERIMENT (LOSAT_X_BPQBLOCK): the block search of `-db` runs only when every definition
    // line of the database reads and none is refused (`x_titles_allow_blocks`), at one thread and
    // in a pool. (A path in the system's temporary directory: Unix spelling.)
    #[cfg(unix)]
    #[test]
    fn blocks_are_allowed_only_when_every_title_reads_and_is_reported() {
        let allowed = |database: &BlastpDatabase| -> [bool; 2] {
            [1, 4].map(|threads| {
                crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                    Ok(x_titles_allow_blocks(database, false, false, pool))
                })
                .unwrap()
            })
        };
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the same databases for the
        // read-only check of outfmt 8, 9, 12 and 5.
        let read_allowed = |database: &BlastpDatabase| -> [bool; 2] {
            [1, 4].map(|threads| {
                crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                    Ok(x_titles_read_allow_blocks(database, pool))
                })
                .unwrap()
            })
        };
        assert_eq!(allowed(&open("prot_v4")), [true, true]);
        assert_eq!(read_allowed(&open("prot_v4")), [true, true]);
        let dir = std::env::temp_dir().join(format!(
            "losat-db-qblock-titles-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        // A header file that cannot be read: the search of all queries reports NCBI's error.
        for entry in
            std::fs::read_dir(fixture("prot_mv.pal").trim_end_matches("prot_mv.pal")).unwrap()
        {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_str().unwrap().to_string();
            if name.starts_with("prot_mv.") && name != "prot_mv.05.phr" {
                std::fs::copy(&path, dir.join(&name)).unwrap();
            }
        }
        let missing = BlastpDatabase::open(
            dir.join("prot_mv").to_str().unwrap(),
            &Default::default(),
            "BLASTP",
        )
        .unwrap();
        assert_eq!(allowed(&missing), [false, false]);
        assert_eq!(read_allowed(&missing), [false, false]);
        // A stored title with a byte outside ' ' ..= '~' (LOSAT's refusal when it is reported).
        for ext in ["pin", "psq"] {
            std::fs::copy(fixture("prot_v4") + "." + ext, dir.join(format!("b.{ext}"))).unwrap();
        }
        let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
        let at = header
            .windows(16)
            .position(|window| window == b"p10 last protein")
            .expect("the title of p10");
        header[at + 5] = 0xE9;
        std::fs::write(dir.join("b.phr"), header).unwrap();
        let bad = BlastpDatabase::open(
            dir.join("b").to_str().unwrap(),
            &Default::default(),
            "BLASTP",
        )
        .unwrap();
        assert_eq!(allowed(&bad), [false, false]);
        assert_eq!(read_allowed(&bad), [false, false]);
        // A title of punctuation only (`, ,`): the tabular id NCBI writes cannot be made from it
        // (outfmt 6 and 7 refuse it, `x_titles_allow_blocks`), but outfmt 8, 9, 12 and 5 write no
        // id, so `x_titles_read_allow_blocks` takes the blocks (the stricter function must not
        // replace it).
        for ext in ["pin", "psq"] {
            std::fs::copy(fixture("prot_v4") + "." + ext, dir.join(format!("c.{ext}"))).unwrap();
        }
        let mut header = std::fs::read(fixture("prot_v4") + ".phr").unwrap();
        let at = header
            .windows(16)
            .position(|window| window == b"p10 last protein")
            .expect("the title of p10");
        header[at..at + 16].copy_from_slice(b", , , , , , , , ");
        std::fs::write(dir.join("c.phr"), header).unwrap();
        let punctuation = BlastpDatabase::open(
            dir.join("c").to_str().unwrap(),
            &Default::default(),
            "BLASTP",
        )
        .unwrap();
        assert_eq!(allowed(&punctuation), [false, false]);
        assert_eq!(read_allowed(&punctuation), [true, true]);
        drop((missing, bad, punctuation));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // EXPERIMENT (LOSAT_X_BPQBLOCK): where the batches' prefetch reads the taxonomy database
    // (`taxids`, DBID item T), the block search of `-db` runs only when no lookup can be refused:
    // the one pass over the headers gives the answer of the titles' check followed by
    // `x_taxids_allow_blocks`. `sprot_tax` of `tests/fixtures/outfmt16/db` (two sequences, taxid
    // 9606) with no taxonomy database, ones whose lookups cannot be refused (a wrong magic number
    // is a message; a scientific name of 9606 that XML2 cannot write is not asked for by the
    // prefetch) and one whose record of 9606 is (an empty `taxdb.btd`, `shared::taxdb_may_refuse`),
    // at one thread and in a pool.
    #[test]
    fn blocks_with_the_prefetch_need_lookups_that_cannot_be_refused() {
        use crate::seqdb::alias::SearchPath;
        let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
        let cwd = std::env::temp_dir().join(format!(
            "losat-db-qblock-prefetch-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&cwd).unwrap();
        let path = tests.join("fixtures/outfmt16/db/sprot_tax");
        for (taxdb, refused) in [
            (None, false),
            (Some("good"), false),
            (Some("bad_magic"), false),
            (Some("high_bytes_utf8"), false),
            (Some("btd_empty"), true),
        ] {
            let taxdb = taxdb.map(|name| tests.join("seqdb/taxdb").join(name));
            let search =
                SearchPath::new(&cwd, taxdb.as_deref().map(std::path::Path::as_os_str), None)
                    .unwrap();
            let database = BlastpDatabase {
                db: SeqDb::open_in(path.to_str().unwrap(), SeqType::Protein, &search).unwrap(),
                name: "sprot_tax".to_string(),
                taxdb: shared::TaxdbRun::default(),
                program: "BLASTP",
            };
            assert_eq!(database.taxdb_found(), taxdb.is_some());
            for threads in [1, 4] {
                crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
                    for pairwise in [false, true] {
                        let titles = x_titles_allow_blocks(&database, pairwise, false, pool);
                        assert!(titles, "{taxdb:?} threads {threads}");
                        let two = titles && x_taxids_allow_blocks(&database, false, pool);
                        assert_eq!(two, !refused, "{taxdb:?} threads {threads}");
                        assert_eq!(
                            x_titles_allow_blocks(&database, pairwise, true, pool),
                            two,
                            "{taxdb:?} threads {threads}"
                        );
                    }
                    Ok(())
                })
                .unwrap();
            }
        }
        std::fs::remove_dir_all(&cwd).unwrap();
    }
}
