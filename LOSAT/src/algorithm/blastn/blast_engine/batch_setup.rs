//! The batch-wide setup that `search_query_batch` does before its subjects: the query block (the
//! contexts of the queries and their BLASTNA buffer), the lookup table with the offset array
//! size, the preliminary stage's search spaces, and the cutoffs of every context (of the batch
//! and of one subject).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:152-273
//! (`SetupQueryInfo_OMF`) and :484-658 (`SetupQueries_OMF`);
//! c++/src/algo/blast/core/blast_nalookup.c:1227-1368 (`BlastMBLookupTableNew`);
//! c++/src/algo/blast/core/blast_parameters.c:279-419 (`BlastInitialWordParametersUpdate`) and
//! :830-994 (`BlastHitSavingParametersUpdate`);
//! c++/src/algo/blast/api/prelim_stage.cpp:219-221 (the preliminary stage's option spaces) and
//! c++/src/algo/blast/core/blast_setup.c:1001-1027 (`BLAST_OneSubjectUpdateParameters`).
//!
//! This module is a child of `run` for the types it keeps private (`QueryContext`,
//! `QueryContextIndex`, `SubjectSource`). The functions are statements of
//! `search_query_batch` moved here in their order; what they read is a struct of references or a
//! parameter, and what they make is what `search_query_batch` then reads under the same names.

use super::{
    diag_container_query_length, fill_strand_query_blastna_buffer, mask_query_letters,
    reverse_mask_intervals, ungapped, ungapped_cutoff_limit_error, x_mb_lut_keep, x_query_drop,
    x_query_prep_par, BatchStage, BlastnArgs, ContextKarlin, LookupTables, QueryContext,
    QueryContextIndex, QueryLetters, QueryStrand, SequenceData, SubjectSource, OFFSET_ARRAY_SIZE,
};
use crate::algorithm::blastn::coordination::{build_lookup_tables, TaskConfig};
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:342-344,370-374
// new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
use crate::algorithm::blastn::ncbi_cutoffs::{
    cutoff_score_for_ungapped_extension, cutoff_score_max_from_evalue, gap_trigger_raw_score,
    GAP_TRIGGER_BIT_SCORE_NUCL,
};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::utils::dust::MaskedInterval;
use anyhow::Result;
use std::sync::Arc;

/// The query block of a batch: one context per strand of each query, with the masks, the
/// offsets in the concatenated query and the lengths NCBI's `BlastQueryInfo` holds
/// (`query_contexts` and what is made of it).
pub(super) struct QueryBlock {
    pub(super) query_contexts: Vec<QueryContext>,
    pub(super) query_context_masks: Vec<Vec<MaskedInterval>>,
    pub(super) query_concat_length: usize,
    pub(super) diag_container_query_length: usize,
    pub(super) query_context_offsets: Vec<i32>,
    pub(super) query_context_index: QueryContextIndex,
    /// The strands of `-strand`, which the encoding reads.
    pub(super) strand: QueryStrand,
}

/// The contexts and offsets of the batch's queries (the strand that is not searched has a context
/// of length 0), before the letters are encoded.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:152-273
/// (`SetupQueryInfo_OMF`, which calls `s_QueryInfo_SetContext`, :64-88, for each context).
pub(super) fn query_block(
    args: &BlastnArgs,
    letters: &x_query_drop::Letters<'_>,
    seq_data: &SequenceData<'_>,
) -> QueryBlock {
    // Build per-context query sequences (plus/minus) for blastn.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:839-847
    // NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    let mut query_contexts: Vec<QueryContext> = Vec::new();
    let mut query_context_masks: Vec<Vec<MaskedInterval>> = Vec::new();
    // NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-81
    // Each context starts after the context before and its sentinel; the contexts of a
    // query without letters take no room (and are invalid: `context_ungapped` in
    // `search_query_batch`).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:229-243
    // `-strand` sets the strand of the query's Seq-loc: the strand that is not searched has a
    // context of length 0 (invalid, without masks and without room in the buffer). Chosen once
    // here; `both` keeps the plain two-strand code below.
    let strand = QueryStrand::from_option(&args.strand);
    let (plus_on, minus_on) = strand.searched();
    let mut query_concat_offset: usize = 0;
    let context_room = |length: usize| if length == 0 { 0 } else { length + 1 };
    for q_idx in 0..seq_data.queries.len() {
        let q_len = letters.length(q_idx);
        let plus_len = if plus_on { q_len } else { 0 };
        let minus_len = if minus_on { q_len } else { 0 };
        let plus_masks = if plus_on {
            seq_data.query_masks[q_idx].clone()
        } else {
            Vec::new()
        };
        let minus_masks = if minus_on {
            reverse_mask_intervals(&seq_data.query_masks[q_idx], q_len)
        } else {
            Vec::new()
        };
        let plus_offset = query_concat_offset;
        let minus_offset = plus_offset + context_room(plus_len);

        query_contexts.push(QueryContext {
            query_idx: q_idx as u32,
            frame: 1,
            query_offset: plus_offset as i32,
            len: plus_len,
            masks: plus_masks.clone(),
        });
        query_contexts.push(QueryContext {
            query_idx: q_idx as u32,
            frame: -1,
            query_offset: minus_offset as i32,
            len: minus_len,
            masks: minus_masks.clone(),
        });

        query_context_masks.push(plus_masks);
        query_context_masks.push(minus_masks);
        query_concat_offset = minus_offset + context_room(minus_len);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:245-250
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:654
    // The searched block leaves out the first and the last sentinel of the buffer. A batch
    // has a query with letters (`run_in_pool` stops on a batch without one).
    let query_concat_length = query_contexts.last().map_or(0, |last| {
        let length = last.len;
        (last.query_offset.max(0) as usize + length + if length == 0 { 1 } else { 2 })
            .saturating_sub(2)
    });
    let diag_container_query_length = diag_container_query_length(&query_contexts);
    // NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:166-174
    let query_context_offsets: Vec<i32> =
        query_contexts.iter().map(|ctx| ctx.query_offset).collect();
    let query_context_index = QueryContextIndex::new(&query_contexts, strand != QueryStrand::Both);
    QueryBlock {
        query_contexts,
        query_context_masks,
        query_concat_length,
        diag_container_query_length,
        query_context_offsets,
        query_context_index,
        strand,
    }
}

/// Encodes the batch's queries into BLASTNA in `encoded_query_buffer` (the caller's, kept from
/// batch to batch), and with `-soft_masking false` writes the masks into the letters of the
/// contexts: returns the copy before the masks (NCBI's `sequence_nomask`), `None` otherwise.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:484-658
/// (`SetupQueries_OMF`) and c++/src/algo/blast/core/blast_filter.c:1349-1406
/// (`BlastSetUp_MaskQuery`, which `BLAST_MainSetUp` calls without mask-at-hash).
pub(super) fn encode_query_block(
    args: &BlastnArgs,
    letters: &x_query_drop::Letters<'_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    block: &mut QueryBlock,
    encoded_query_buffer: &mut Vec<u8>,
) -> Option<Vec<u8>> {
    let strand = block.strand;
    let query_concat_length = block.query_concat_length;
    let query_context_offsets = &block.query_context_offsets;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:498-604
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:120-121
    // The buffer is the one BLASTNA copy of the query contexts; the searched block is the
    // buffer without its first and last sentinel, and a context is a slice of it.
    // EXPERIMENT (LOSAT_X_QPREPPAR): with the switch on and a search pool of two or more threads
    // the queries are encoded by the pool, each into its own part of the buffer; the buffer is
    // the same (`x_query_prep_par.rs`).
    // Diagnostic (LOSAT_TIMING): the batch's query encoding for the common line
    // (utils/xtiming.rs). No NCBI counterpart.
    let t_encode = crate::utils::xtiming::clock();
    // `-strand plus` or `minus` builds the buffer serially from the contexts that are searched
    // (`fill_strand_query_blastna_buffer`); `both` is the code of the search without the option.
    if strand == QueryStrand::Both {
        x_query_prep_par::query_buffer(
            letters,
            query_context_offsets,
            query_concat_length,
            parallel_pool,
            encoded_query_buffer,
        );
    } else {
        fill_strand_query_blastna_buffer(
            letters,
            query_context_offsets,
            query_concat_length,
            strand,
            encoded_query_buffer,
        );
    }
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:609-623
    // With `-soft_masking true` (the default) the DUST, lower-case and WindowMasker masks only
    // restrict the lookup table through query_context_masks, and there is one copy of the
    // query. `false` searches a copy with the masks written into it.
    if args.soft_masking() {
        None
    } else {
        mask_query_letters(
            encoded_query_buffer,
            &mut block.query_contexts,
            &block.query_context_masks,
        )
    }
}

/// What `batch_lookup_tables` reads.
pub(super) struct LookupInputs<'a> {
    pub(super) config: &'a TaskConfig,
    pub(super) args: &'a BlastnArgs,
    /// The batch's query block (`QueryBlock`) as the lookup table sees it.
    pub(super) encoded_queries: &'a [&'a [u8]],
    pub(super) query_context_masks: &'a [Vec<MaskedInterval>],
    pub(super) query_context_offsets: &'a [i32],
    pub(super) subject_records: Option<&'a [FastaRecord]>,
    pub(super) subject_packed_cache: Option<&'a [Vec<u8>]>,
    pub(super) approx_table_entries: usize,
    pub(super) x_keep: Option<&'a Arc<x_mb_lut_keep::Keep>>,
    pub(super) subjects: SubjectSource<'a>,
    /// The batch is split into query chunks, which have their own tables.
    pub(super) split: bool,
}

/// The lookup table of the batch (empty for a batch split into query chunks), the scan step and
/// the size of the offset array of the subjects' word finder.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1227-1368
/// (`BlastMBLookupTableNew`), c++/src/algo/blast/core/blast_engine.c:970-1069
/// (`s_BlastSetUpAuxStructures`, which sizes the offset pairs) and
/// c++/src/algo/blast/core/lookup_wrap.c:255-289 (`GetOffsetArraySize`).
pub(super) fn batch_lookup_tables(inputs: &LookupInputs<'_>) -> (LookupTables, usize, usize) {
    let &LookupInputs {
        config,
        args,
        encoded_queries,
        query_context_masks,
        query_context_offsets,
        subject_records,
        subject_packed_cache: subject_packed_cache_ref,
        approx_table_entries,
        x_keep,
        subjects,
        split,
    } = inputs;
    // Build lookup tables
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1270-1306
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    let subjects_for_lookup: &[FastaRecord] = subject_records.unwrap_or(&[]);
    // A split batch has no lookup table of its own; its chunks have theirs.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-207
    let (lookup_tables, scan_step) = if split {
        (
            LookupTables {
                two_stage_lookup: None,
                pv_direct_lookup: None,
                na_lookup: None,
            },
            config.scan_step,
        )
    } else {
        // Diagnostic (LOSAT_TIMING): the batch's lookup table for the common line
        // (utils/xtiming.rs). No NCBI counterpart.
        let t_build_lookup = crate::utils::xtiming::clock();
        let tables = build_lookup_tables(
            config,
            args,
            encoded_queries,
            query_context_masks,
            query_context_offsets,
            subjects_for_lookup,
            subject_packed_cache_ref,
            approx_table_entries,
            x_keep,
            Some(subjects.metadata().db_len_total as u64),
        );
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::BuildLookup, t_build_lookup);
        tables
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:991-1041
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
    let offset_array_size = if let Some(two_stage) = lookup_tables.two_stage_lookup.as_ref() {
        OFFSET_ARRAY_SIZE + two_stage.longest_chain()
    } else if let Some(pv_direct) = lookup_tables.pv_direct_lookup.as_ref() {
        OFFSET_ARRAY_SIZE + pv_direct.longest_chain()
    } else if let Some(na_lookup) = lookup_tables.na_lookup.as_ref() {
        OFFSET_ARRAY_SIZE + na_lookup.longest_chain()
    } else {
        OFFSET_ARRAY_SIZE
    };
    (lookup_tables, scan_step, offset_array_size)
}

/// What `batch_spaces` reads.
pub(super) struct SpaceInputs<'a> {
    pub(super) args: &'a BlastnArgs,
    pub(super) stage: BatchStage<'a>,
    pub(super) subjects: SubjectSource<'a>,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) context_valid: &'a [bool],
    pub(super) search_karlin: &'a [ContextKarlin],
    /// The setup's effective search space of each context (`context_statistics`).
    pub(super) query_eff_searchsp: &'a [i64],
}

/// The preliminary stage's search spaces of a batch.
pub(super) struct BatchSpaces {
    /// The database has the statistics of each subject (`DbStatistics::per_subject`): its subjects
    /// read their own cutoffs and spaces.
    pub(super) per_subject: bool,
    /// The preliminary stage's option spaces of the batch's contexts: the setup's reads them in
    /// place of its own and the statistics of each subject read them. `None`: the setup's own.
    pub(super) option_spaces: Option<Vec<i64>>,
    /// The preliminary stage's effective search spaces where they differ from the setup's
    /// `query_eff_searchsp` (which the traceback and the reports read): a context's option space
    /// when it is not 0. `None`: the setup's.
    pub(super) prelim_eff_searchsp: Option<Vec<i64>>,
}

/// The preliminary stage's option spaces of the batch's contexts, which its setup reads in place
/// of its own where the statistics' total before the fallback to OID 0 is not positive, and which
/// the statistics of each subject read (`DbStatistics::per_subject`). A query chunk has the
/// batch's (`BatchStage::ChunkPrelim`): its own setup's spaces are already them.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-221 (the memento and
/// `SplitQuery_SetEffectiveSearchSpace`, whose spaces the preliminary stage reads; the spaces of
/// one option are `stats::setup::preliminary_option_space`)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1444 (the statistics of
/// each subject, `BLAST_OneSubjectUpdateParameters`, where the database length is 0)
pub(super) fn batch_spaces(inputs: &SpaceInputs<'_>) -> Result<BatchSpaces> {
    let &SpaceInputs {
        args,
        stage,
        subjects,
        query_contexts,
        context_valid,
        search_karlin,
        query_eff_searchsp,
    } = inputs;
    let x_statistics = match subjects {
        SubjectSource::Database(database) => Some(database.statistics()),
        SubjectSource::Prepared(_) => None,
    };
    let per_subject = x_statistics.is_some_and(|statistics| statistics.per_subject());
    let option_spaces: Option<Vec<i64>> = match (stage, x_statistics) {
        (BatchStage::Search { .. }, Some(statistics))
            if statistics.gumbel_db_length <= 0 || per_subject =>
        {
            Some(
                query_contexts
                    .iter()
                    .zip(search_karlin)
                    .zip(context_valid)
                    .map(|((ctx, blocks), &valid)| {
                        crate::stats::setup::preliminary_option_space(
                            if valid { ctx.len as i64 } else { 0 },
                            statistics.gumbel_db_length,
                            statistics.request_count,
                            args.db_length(),
                            args.options_searchsp(),
                            false,
                            &blocks.gapped,
                            "BLASTN",
                        )
                    })
                    .collect::<Result<_, _>>()?,
            )
        }
        (
            BatchStage::ChunkPrelim {
                batch_eff_searchsp, ..
            },
            _,
        ) if per_subject => Some(batch_eff_searchsp.to_vec()),
        _ => None,
    };
    // The preliminary stage's spaces where they differ from the setup's (`query_eff_searchsp`,
    // which the traceback and the reports read): a context's option space when it is not 0.
    let prelim_eff_searchsp: Option<Vec<i64>> = match (stage, &option_spaces, x_statistics) {
        (BatchStage::Search { .. }, Some(spaces), Some(statistics))
            if statistics.gumbel_db_length <= 0 =>
        {
            Some(
                query_eff_searchsp
                    .iter()
                    .zip(spaces)
                    .map(|(&setup, &option)| if option != 0 { option } else { setup })
                    .collect(),
            )
        }
        _ => None,
    };
    Ok(BatchSpaces {
        per_subject,
        option_spaces,
        prelim_eff_searchsp,
    })
}

/// What `batch_cutoffs` reads. The sets of the batch, of its preliminary stage and of one subject
/// differ in `query_eff_searchsp` and `subject_length` only
/// (`CutoffInputs { query_eff_searchsp: spaces, ..inputs }`).
#[derive(Clone, Copy)]
pub(super) struct CutoffInputs<'a> {
    pub(super) subjects: SubjectSource<'a>,
    /// NCBI's `Int8` (`SubjectMetadata::db_len_total`).
    pub(super) db_len_total: i64,
    /// NCBI's `Int4` (`SubjectMetadata::db_num_seqs`).
    pub(super) db_num_seqs: i32,
    pub(super) query_contexts: &'a [QueryContext],
    pub(super) context_valid: &'a [bool],
    pub(super) query_eff_searchsp: &'a [i64],
    pub(super) search_karlin: &'a [ContextKarlin],
    pub(super) x_dropoff_init: &'a [i32],
    /// `-min_raw_gapped_score` (`TaskConfig::min_raw_gapped_score`).
    pub(super) min_raw_gapped_score: Option<i32>,
    /// NCBI's matrix-only scoring (`TaskConfig::matrix_only`): the word cutoff of a gapped search
    /// is the ungapped search's formula with the gapped blocks, `-min_raw_gapped_score` caps it at
    /// half its value, and the X-dropoffs keep their initial values.
    pub(super) matrix_only: bool,
    pub(super) evalue_threshold: f64,
    /// The subject length of the word cutoffs of an ungapped search: the subject's own where NCBI
    /// computes the statistics of each subject (`BlastInitialWordParametersUpdate`'s `subj_length`,
    /// blast_setup.c:1018-1021), `None` for the average of the subject set (blast_engine.c:1372-1374).
    pub(super) subject_length: Option<u64>,
}

/// The cutoff scores, the hit-saving cutoff scores, the X-dropoffs and the reduced cutoff scores
/// of each context, in that order.
pub(super) type BatchCutoffs = (Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>);

/// The cutoffs of every context of the batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:279-419
/// (`BlastInitialWordParametersUpdate`) and c++/src/algo/blast/core/blast_parameters.c:830-994
/// (`BlastHitSavingParametersUpdate`).
pub(super) fn batch_cutoffs<const UNGAPPED: bool>(
    inputs: &CutoffInputs<'_>,
) -> Result<BatchCutoffs> {
    let &CutoffInputs {
        subjects,
        db_len_total,
        db_num_seqs,
        query_contexts,
        context_valid: context_valid_ref,
        query_eff_searchsp,
        search_karlin: search_karlin_ref,
        x_dropoff_init: x_dropoff_init_ref,
        min_raw_gapped_score,
        matrix_only,
        evalue_threshold,
        subject_length,
    } = inputs;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348
    // The word cutoff from the e-value of `BLAST_Cutoffs` (`ungapped::word_cutoff_score`): every
    // ungapped search, and a gapped one under matrix-only scoring.
    let word_cutoff_formula = UNGAPPED || matrix_only;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:336-377
    // NCBI updates these per subject (`BlastInitialWordParametersUpdate`); LOSAT's values do not
    // depend on the subject, so they are computed once for the batch instead of in every subject
    // (`search_subject`): the cutoff scores, the hit-saving cutoff scores, the X-dropoffs and the
    // reduced cutoff scores of each context.
    // The subject length that NCBI's word cutoff of an ungapped search takes: the average
    // length of the subjects (a `-subject` set's `(Uint4)(total / count)`, a database's
    // `(Int4)(total / MAX(1, count))`); the dbscan mode of the command line holds the total
    // of a `-subject` set, so `BlastInitialWordParametersUpdate` runs once for the batch.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1445,
    // c++/src/algo/blast/api/seqsrc_multiseq.cpp:245-266, seqsrc_seqdb.cpp:203-209
    let x_avg_subject_length: u64 =
        if let (true, Some(length)) = (word_cutoff_formula, subject_length) {
            length
        } else if word_cutoff_formula {
            match subjects {
                // `BlastSeqSrcGetAvgSeqLen` reads `GetTotalLength` and `GetNumSeqs`, not the
                // statistics' `STATS_*` keys (`SeqDb::average_length`).
                SubjectSource::Database(database) => database.average_length() as u32 as u64,
                SubjectSource::Prepared(_) => ungapped::average_subject_length_subject_set(
                    db_len_total as u64,
                    db_num_seqs as u64,
                ),
            }
        } else {
            0
        };
    let mut cutoff_scores: Vec<i32> = Vec::with_capacity(query_contexts.len());
    let mut hit_saving_cutoff_scores: Vec<i32> = Vec::with_capacity(query_contexts.len());
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:206-211;
    // c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181;
    // c++/src/algo/blast/core/blast_engine.c:1434-1445
    // if (dbscan_mode) { ... m_iTotalLength += (Int8) (*iter)->length; }
    // if (db_length == 0) { BLAST_OneSubjectUpdateParameters(...); }
    // The CLI subject set has a nonzero total length. Its context search
    // space is retained across subjects, just as for output statistics.
    for context_idx in 0..query_contexts.len() {
        // An invalid context has an infinite cutoff and, unset, a zero X-drop and
        // reduced cutoff (blast_parameters.c:324-332, in `context_statistics::context_statistics`).
        if !context_valid_ref[context_idx] {
            cutoff_scores.push(i32::MAX);
            hit_saving_cutoff_scores.push(i32::MAX);
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:925-946
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:876-882
        // `-min_raw_gapped_score` replaces the e-value's cutoff (the scale factor of
        // BLASTN is 1); the word cutoff, the reduced cutoff and megablast's ungapped
        // X-drop below follow it.
        let hit_saving_cutoff = match min_raw_gapped_score {
            Some(score) => score,
            None => cutoff_score_max_from_evalue(
                evalue_threshold,
                query_eff_searchsp[context_idx],
                &search_karlin_ref[context_idx].gapped,
            ),
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:882-893
        // The cap of the word cutoff (`cutoff_score_max`): the hit-saving cutoff, or half of
        // `-min_raw_gapped_score` (C's integer division) under matrix-only scoring.
        let cutoff_score_max = match min_raw_gapped_score {
            Some(score) if matrix_only => score / 2,
            _ => hit_saving_cutoff,
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:342-374
        // new_cutoff = MIN(new_cutoff, hit_params->cutoffs[context].cutoff_score_max);
        let gap_trigger = gap_trigger_raw_score(
            GAP_TRIGGER_BIT_SCORE_NUCL,
            &search_karlin_ref[context_idx].ungapped,
        );
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:346-374
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:299-303
        let cutoff = if word_cutoff_formula {
            let blocks = &search_karlin_ref[context_idx];
            ungapped::word_cutoff_score(
                query_contexts[context_idx].len,
                x_avg_subject_length,
                if UNGAPPED {
                    &blocks.ungapped
                } else {
                    &blocks.gapped
                },
                cutoff_score_max,
            )
            .map_err(ungapped_cutoff_limit_error)?
        } else {
            cutoff_score_for_ungapped_extension(gap_trigger, cutoff_score_max, 1.0)
        };
        cutoff_scores.push(cutoff);
        hit_saving_cutoff_scores.push(hit_saving_cutoff);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:380-397
    let mut x_dropoff_scores: Vec<i32> = Vec::with_capacity(cutoff_scores.len());
    for ((cutoff, &x_dropoff_init), &valid) in cutoff_scores
        .iter()
        .zip(x_dropoff_init_ref)
        .zip(context_valid_ref)
    {
        let x_dropoff = if !valid {
            0
        } else if x_dropoff_init == 0 && !matrix_only {
            *cutoff
        } else {
            x_dropoff_init
        };
        x_dropoff_scores.push(x_dropoff);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:408-412
    let mut reduced_cutoff_scores: Vec<i32> = Vec::with_capacity(cutoff_scores.len());
    for (cutoff, &valid) in cutoff_scores.iter().zip(context_valid_ref) {
        reduced_cutoff_scores.push(if valid {
            (0.8 * (*cutoff as f64)) as i32
        } else {
            0
        });
    }
    Ok((
        cutoff_scores,
        hit_saving_cutoff_scores,
        x_dropoff_scores,
        reduced_cutoff_scores,
    ))
}

/// The cutoffs of the batch and of its preliminary stage.
pub(super) struct BatchCutoffSet {
    /// The cutoffs of the batch's contexts with the setup's spaces
    /// (`CutoffInputs::query_eff_searchsp`).
    pub(super) batch: BatchCutoffs,
    /// The cutoffs of the preliminary stage's own spaces where they differ (`prelim_eff_searchsp`),
    /// `None` where they are the batch's.
    pub(super) prelim: Option<BatchCutoffs>,
}

/// The cutoffs of the batch and of the preliminary stage's own spaces.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:348-374
///
/// NCBI uses UNGAPPED params (kbp_std) for gap_trigger calculation
/// and GAPPED params (kbp_gap) for cutoff_score_max calculation.
/// The cutoffs of every context of the batch (`batch_cutoffs`): they depend on
/// the context only (the gapped calculation's branch above: gap_trigger capped by
/// cutoff_score_max), so the batch computes them once, and both stages read them from their
/// batch (`PrelimBatch`, `TraceBatch`).
///
/// The preliminary stage's own spaces `prelim_eff_searchsp` give the second set. The batch's
/// cutoffs are computed first; both can fail with the same refusal (`ungapped_cutoff_limit_error`).
pub(super) fn batch_cutoff_set<const UNGAPPED: bool>(
    inputs: &CutoffInputs<'_>,
    prelim_eff_searchsp: Option<&[i64]>,
) -> Result<BatchCutoffSet> {
    let batch = batch_cutoffs::<UNGAPPED>(inputs)?;
    let prelim = match prelim_eff_searchsp {
        Some(spaces) => Some(batch_cutoffs::<UNGAPPED>(&CutoffInputs {
            query_eff_searchsp: spaces,
            ..*inputs
        })?),
        None => None,
    };
    Ok(BatchCutoffSet { batch, prelim })
}

/// The statistics of one subject (`per_subject`): each context's effective search space for the
/// subject's length with the option spaces `spaces` (the preliminary stage's, or the traceback's
/// own: `-searchsp`, which the memento gives it back, prelim_stage.cpp:219), and the cutoffs of
/// `BlastHitSavingParametersUpdate` and `BlastInitialWordParametersUpdate` (the X-drops keep their
/// initial values). The traceback of a list reads the subject's hit-saving cutoffs
/// (`Blast_HSPReevaluateWithAmbiguitiesGapped`, its `hit_params` being thread 0's in NCBI's one
/// thread, blast_traceback.c:1598-1650) and the e-values of the setup's spaces. `inputs` are the
/// batch's (`query_eff_searchsp` and `subject_length` are replaced for each subject);
/// `db_length` is `-dbsize`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1011-1022 (the function's first
/// lines, 1001-1014, are cited at `stats::setup::subject_eff_lengths`)
pub(super) fn subject_statistics<'a, const UNGAPPED: bool>(
    inputs: CutoffInputs<'a>,
    db_length: i64,
) -> impl Fn(usize, &[i64]) -> Result<(BatchCutoffs, Vec<i64>)> + Sync + 'a {
    move |subject_length: usize, spaces: &[i64]| -> Result<(BatchCutoffs, Vec<i64>)> {
        let subject_spaces: Vec<i64> = inputs
            .query_contexts
            .iter()
            .zip(inputs.search_karlin)
            .zip(inputs.context_valid)
            .zip(spaces)
            .map(|(((ctx, blocks), &valid), &space)| {
                if !valid {
                    return Ok(space);
                }
                Ok(crate::stats::setup::subject_eff_lengths(
                    ctx.len as i64,
                    subject_length as i64,
                    inputs.db_num_seqs,
                    db_length,
                    space,
                    false,
                    &blocks.gapped,
                    "BLASTN",
                )?
                .eff_searchsp)
            })
            .collect::<Result<_, crate::stats::setup::SearchSpaceOverflow>>()?;
        let cutoffs = batch_cutoffs::<UNGAPPED>(&CutoffInputs {
            query_eff_searchsp: &subject_spaces,
            subject_length: Some(subject_length as u64),
            ..inputs
        })?;
        Ok((cutoffs, subject_spaces))
    }
}
