//! BLASTX score-block and preliminary parameter setup, using Rust calculations.
use super::{
    args::ResolvedOptions,
    query_setup::{PreparedQueryBatch, QueryContext},
};
use crate::algorithm::tblastx::ncbi_cutoffs::{gap_trigger_raw_score, x_drop_raw_score};
use crate::config::{ProteinScoringSpec, ScoringMatrix};
use crate::stats::cutoffs::{cutoff_score, CutoffConversion, SumStatsCutoffs};
use crate::stats::karlin_calc::{
    apply_check_ideal, compute_aa_composition, compute_karlin_params_ungapped,
    compute_score_freq_profile_for_matrix, compute_std_aa_composition,
};
use crate::stats::spouge::{blast_spouge_etos, lookup_protein_gumbel_params, BlastGumbelBlk};
use crate::stats::{lookup_protein_params, KarlinParams};
use anyhow::{bail, ensure, Result};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2743-2754
#[derive(Clone, Debug)]
pub struct ContextParameters {
    pub valid: bool,
    pub ungapped: KarlinParams,
    pub search_space: i64,
    pub length_adjustment: i32,
    pub word_cutoff: i32,
    pub word_xdrop_init: i32,
    pub word_xdrop: i32,
    pub hit_cutoff: i32,
    pub hit_cutoff_max: i32,
}

/// No NCBI counterpart: whether `LOSAT_X_BXKAONCE` is on (a strict-set switch, on unless `=0`):
/// a query batch's Karlin blocks computed once when the batch is not split (`search.rs`), and the
/// BLOSUM62 ideal block once per process (`ideal_block`).
pub(crate) fn ka_once_on() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BXKAONCE"))
}

/// The BLOSUM62 ideal Karlin block of the standard composition against itself, as `score_block`
/// computes it for each call without `LOSAT_X_BXKAONCE`.
fn ideal_block_computed(std: &[f64; crate::utils::matrix::BLASTAA_SIZE]) -> KarlinParams {
    let ideal_profile =
        compute_score_freq_profile_for_matrix(std, std, -4, 11, ScoringMatrix::Blosum62);
    compute_karlin_params_ungapped(&ideal_profile).expect("BLOSUM62 ideal KA block")
}

/// No NCBI counterpart (`LOSAT_X_BXKAONCE`): the ideal block is a constant of the program (no
/// input), so it is computed on first use and shared; the value `ideal_block_computed` gives.
/// BLOSUM62's only: another matrix's searches compute theirs in `score_block_of`.
fn ideal_block() -> KarlinParams {
    static IDEAL: std::sync::OnceLock<KarlinParams> = std::sync::OnceLock::new();
    *IDEAL.get_or_init(|| ideal_block_computed(&compute_std_aa_composition()))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2776-2801
/// `score_block` with the search's matrix: BLOSUM62's is `score_block`, another matrix's
/// `score_block_of`.
pub fn score_block_for(
    batch: &PreparedQueryBatch,
    matrix: ScoringMatrix,
) -> Vec<ContextParameters> {
    score_block_for_in(batch, matrix, false)
}

/// `score_block_for`; with `on_pool` (a reader batch searched whole, which is left the threads of
/// its places, `search::x_bx_mem::whole_batch_threads`) the contexts' blocks are computed on the
/// pool this thread runs on. Each block reads its own context alone, and the blocks are collected
/// in context order, so the values and their order are those of the loop. No NCBI counterpart of
/// its own (scheduling only): NCBI fills the contexts' blocks one after another in this loop.
pub(crate) fn score_block_for_in(
    batch: &PreparedQueryBatch,
    matrix: ScoringMatrix,
    on_pool: bool,
) -> Vec<ContextParameters> {
    if matrix == ScoringMatrix::Blosum62 {
        score_block_in(batch, on_pool)
    } else {
        score_block_of(batch, matrix, on_pool)
    }
}

/// No NCBI counterpart (scheduling only): `block` of each context of `batch`, in context order;
/// on the pool this thread runs on where `on_pool` (`score_block_for_in`).
fn context_blocks(
    batch: &PreparedQueryBatch,
    on_pool: bool,
    block: impl Fn(&QueryContext) -> ContextParameters + Sync + Send,
) -> Vec<ContextParameters> {
    #[cfg(feature = "parallel")]
    if on_pool {
        use rayon::prelude::*;
        return batch.contexts.par_iter().map(block).collect();
    }
    let _ = on_pool;
    batch.contexts.iter().map(block).collect()
}

/// `score_block` with another matrix than BLOSUM62: the contexts' score profiles over the
/// matrix's score range (`sbp->loscore`, `sbp->hiscore`) and its ideal block
/// (`ideal_karlin_params_of`, the block TBLASTN's and BLASTP's searches with that matrix take).
/// Not inlined, so that BLOSUM62's `score_block` keeps its code.
#[cold]
#[inline(never)]
fn score_block_of(
    batch: &PreparedQueryBatch,
    matrix: ScoringMatrix,
    on_pool: bool,
) -> Vec<ContextParameters> {
    use crate::algorithm::tblastx::lookup::{ideal_karlin_params_of, matrix_score_bounds};
    let std = compute_std_aa_composition();
    let ideal = ideal_karlin_params_of(matrix);
    let (low, high) = matrix_score_bounds(matrix);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3127-3130 (the alpha and
    // beta of an ungapped search, which the length adjustment reads, blast_setup.c:815)
    // The matrix's ungapped row (BLOSUM62's are the computed block's defaults), as TBLASTX's
    // contexts take them (`tblastx/blast_engine/run_impl.rs`).
    let row = crate::stats::tables::lookup_protein_params_ungapped(matrix);
    context_blocks(batch, on_pool, |c| {
        let mut seq = vec![0];
        seq.extend_from_slice(&batch.sequence_start[1 + c.offset..1 + c.offset + c.length]);
        seq.push(0);
        let comp = compute_aa_composition(&seq, c.length);
        let profile = compute_score_freq_profile_for_matrix(&comp, &std, low, high, matrix);
        let ka = compute_karlin_params_ungapped(&profile).ok().map(|p| {
            let mut ka = apply_check_ideal(p, ideal);
            ka.alpha = row.alpha;
            ka.beta = row.beta;
            ka
        });
        ContextParameters {
            valid: c.is_valid && ka.is_some(),
            ungapped: ka.unwrap_or_default(),
            search_space: 0,
            length_adjustment: 0,
            word_cutoff: i32::MAX,
            word_xdrop_init: 0,
            word_xdrop: 0,
            hit_cutoff: i32::MAX,
            hit_cutoff_max: 0,
        }
    })
}

pub fn score_block(batch: &PreparedQueryBatch) -> Vec<ContextParameters> {
    score_block_in(batch, false)
}

/// `score_block`, on the pool where `on_pool` (`score_block_for_in`).
fn score_block_in(batch: &PreparedQueryBatch, on_pool: bool) -> Vec<ContextParameters> {
    let std = compute_std_aa_composition();
    let ideal = if ka_once_on() {
        ideal_block()
    } else {
        ideal_block_computed(&std)
    };
    context_blocks(batch, on_pool, |c| {
        let mut seq = vec![0];
        seq.extend_from_slice(&batch.sequence_start[1 + c.offset..1 + c.offset + c.length]);
        seq.push(0);
        let comp = compute_aa_composition(&seq, c.length);
        let profile =
            compute_score_freq_profile_for_matrix(&comp, &std, -4, 11, ScoringMatrix::Blosum62);
        let ka = compute_karlin_params_ungapped(&profile)
            .ok()
            .map(|p| apply_check_ideal(p, ideal));
        ContextParameters {
            valid: c.is_valid && ka.is_some(),
            ungapped: ka.unwrap_or_default(),
            search_space: 0,
            length_adjustment: 0,
            word_cutoff: i32::MAX,
            word_xdrop_init: 0,
            word_xdrop: 0,
            hit_cutoff: i32::MAX,
            hit_cutoff_max: 0,
        }
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:913-920
/// The database length of the Gumbel block: `-dbsize` when it is not 0 (also below 0, as
/// NCBI), else the real total `db_length`; not divided by 3, since BLASTX's subject is not
/// translated (blast_setup.c:921-922).
pub fn gumbel_db_length(options: &ResolvedOptions, db_length: i64) -> i64 {
    if options.db_length != 0 {
        options.db_length
    } else {
        db_length
    }
}

/// The gapped Karlin block of the search's matrix and gap costs (`kbp_gap`, from NCBI's tables);
/// an ungapped search reads none (its gap costs need not be in the matrix's table, `args.rs`), and
/// gets the default block, which it never reads.
pub(crate) fn gapped_karlin(options: &ResolvedOptions) -> KarlinParams {
    if options.gapped {
        lookup_protein_params(&ProteinScoringSpec {
            matrix: options.scoring,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        })
    } else {
        KarlinParams::default()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
/// The Gumbel block of a gapped search, none with `-searchsp` (at any value) or `OLD_FSC`, nor in
/// an ungapped search (`ResolvedOptions::has_gumbel_block`).
pub fn gumbel_block(options: &ResolvedOptions, db_length: i64) -> Option<BlastGumbelBlk> {
    if !options.has_gumbel_block() {
        return None;
    }
    lookup_protein_gumbel_params(
        &ProteinScoringSpec {
            matrix: options.scoring,
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
        },
        gumbel_db_length(options, db_length),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:784-824
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:836-847
/// LOSAT's refusal of a search whose effective search space overflows NCBI's `Int8` (the shared
/// one of `stats::setup::calc_eff_lengths`).
pub(crate) use crate::stats::setup::SearchSpaceOverflow;

pub fn effective_lengths(
    batch: &PreparedQueryBatch,
    parameters: &mut [ContextParameters],
    options: &ResolvedOptions,
    db_length: i64,
    db_count: i32,
    fixed_spaces: Option<&[i64]>,
) -> Result<()> {
    ensure!(
        parameters.len() == batch.contexts.len(),
        "one statistical block per BLASTX context"
    );
    let options_space = options.effective_search_space;
    // `stats::setup::eff_lengths_database`: BLASTX's subject is not translated.
    let Some((db_length, db_count)) = crate::stats::setup::eff_lengths_database(
        db_length,
        db_count,
        options.db_length,
        fixed_spaces.is_some() || options_space != 0,
        false,
    ) else {
        return Ok(());
    };
    let gapped = gapped_karlin(options);
    for (i, (c, p)) in batch.contexts.iter().zip(parameters).enumerate() {
        // NCBI's loop starts at the batch's first context (3 with `-strand minus`): the contexts
        // before it keep the zeros of `BlastQueryInfoNew`'s `calloc`, which RPSTBLASTN's redo
        // reads for the first query's frame +1 (`one_query_info->contexts[0]`).
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:770-772
        if i < batch.first_context {
            p.length_adjustment = 0;
            p.search_space = 0;
            continue;
        }
        let space = fixed_spaces.map_or(options_space, |spaces| spaces[i]);
        if !p.valid {
            p.length_adjustment = 0;
            p.search_space = space;
            continue;
        }
        let ka = if options.gapped { gapped } else { p.ungapped };
        let lengths = crate::stats::setup::calc_eff_lengths(
            c.length as i64,
            db_length,
            db_count,
            &ka,
            space,
            options.program().upper,
        )?;
        p.length_adjustment = lengths.length_adjustment;
        p.search_space = lengths.eff_searchsp;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:908-944
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-987
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:935-975
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:340-383
/// The numbers of the subjects that the parameters of a search read (`Subjects::numbers`), as
/// NCBI's sequence source answers them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SubjectNumbers {
    /// The total length of the statistics (blast_setup.c:908-911): `STATS_TOTLEN` or
    /// `BlastSeqSrcGetTotLen` (`SeqDb::statistics_numbers`), the sum of a FASTA set's records.
    pub total: i64,
    /// The number of sequences of the statistics (blast_setup.c:925-928), NCBI's `Int4`:
    /// `STATS_NSEQ` or `GetNumSeqs`, 0 for an alias file's `NSEQ 0` and negative for a negative
    /// one; a FASTA set's records.
    pub count: i32,
    /// The Gumbel block's database length before `-dbsize` (blast_setup.c:913-920): a
    /// database's total before the fallback to OID 0 (`DbStatistics::gumbel_db_length`, 0 or
    /// negative for an alias file's `LENGTH` of 0 or of 2^63 or more), a FASTA set's total.
    pub gumbel_total: i64,
    /// `BlastSeqSrcGetMinSeqLen` (blast_setup.c:969-973): a database's `GetMinLength`
    /// (`SeqDb::min_length`, signed 32-bit `Int4`), a FASTA set's shortest record raised to
    /// `BLAST_SEQSRC_MINLENGTH`.
    pub shortest: i32,
    /// `BlastSeqSrcGetAvgSeqLen` (blast_engine.c:1372-1374, the word parameters): a database's
    /// `GetTotalLength` divided by the larger of 1 and `GetNumSeqs` (`SeqDb::average_length`, no `STATS_*` key), a
    /// FASTA set's total over its records.
    pub average: i32,
    /// The number of sequences of the preliminary stage's option spaces (`STATS_NSEQ` or
    /// `GetNumSeqs` without the fallback: `DbStatistics::request_count`), a FASTA set's records.
    pub request_count: i32,
}

impl SubjectNumbers {
    /// The numbers of one FASTA subject of `length` residues.
    #[cfg(test)]
    pub fn one(length: usize) -> Self {
        Self {
            total: length as i64,
            count: 1,
            gumbel_total: length as i64,
            shortest: i32::try_from(length.max(super::subjects::BLAST_SEQSRC_MINLENGTH))
                .unwrap_or(i32::MAX),
            average: length as i32,
            request_count: 1,
        }
    }
}

pub fn subject_parameters(
    batch: &PreparedQueryBatch,
    parameters: &mut [ContextParameters],
    options: &ResolvedOptions,
    numbers: SubjectNumbers,
    fixed_spaces: Option<&[i64]>,
) -> Result<i32> {
    let SubjectNumbers {
        total: db_length,
        count: db_count,
        gumbel_total,
        shortest,
        average,
        request_count: _,
    } = numbers;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
    // BLASTX's subject is not translated, so the shortest length is not divided by 3. Without the
    // Gumbel block (`ResolvedOptions::has_gumbel_block`: `-ungapped`, `-searchsp`, `OLD_FSC`) the
    // statistics' average divides by their number of sequences, which an alias file's `NSEQ 0`
    // makes 0: NCBI crashes there (SIGFPE), LOSAT refuses (`common::database::check_setup_count`,
    // DECISIONS 14 of line RG). A set of subjects whose average is below one letter fails the
    // batch's setup: the engine error after the start of the outfmt 0 report, exit status 3, as
    // every other error of this setup (`native::database_search_error`; NCBI 2.17.0 blastx
    // `-subject` of three empty records, or of `M`, `` and ``, with `-searchsp`, `OLD_FSC` or
    // `-ungapped -comp_based_stats 0`, env parity probe f_probe of 2026-11-12).
    let gumbel_block_used = options.has_gumbel_block();
    crate::algorithm::common::database::check_setup_count(
        db_length,
        db_count,
        gumbel_block_used,
        "BLASTX",
    )?;
    let Some(n) = crate::stats::setup::setup_subject_length(
        db_length,
        db_count,
        gumbel_block_used.then_some((shortest, false)),
    ) else {
        bail!(crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT);
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:759-760 (cited at
    // `stats::setup::hit_saving_parameters_fail`)
    if crate::stats::setup::hit_saving_parameters_fail(options.sum_stats, options.gapped, n) {
        bail!(crate::stats::setup::HIT_SAVING_PARAMETERS_ERROR);
    }
    effective_lengths(
        batch,
        parameters,
        options,
        db_length,
        db_count,
        fixed_spaces,
    )?;
    let gumbel = gumbel_block(options, gumbel_total);
    context_cutoffs(
        batch,
        parameters,
        options,
        gumbel.as_ref(),
        n,
        average,
        options.composition > 1,
    )?;
    Ok(if options.gapped {
        gapped_x_drops(options, &gapped_karlin(options)).0
    } else {
        0
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-220 and
// c++/src/algo/blast/api/split_query_aux_priv.cpp:155-181 (cited at
// `stats::setup::preliminary_option_space`)
/// The option spaces that NCBI's preliminary stage sets for the contexts of the query batch
/// `batch` (`SplitQuery_SetEffectiveSearchSpace`), where they differ from the setup's own spaces
/// (`effective_lengths` of the batch): a database whose statistics' total before the fallback to
/// OID 0 is not positive (`SubjectNumbers::gumbel_total`), without `-searchsp`. `None` otherwise,
/// where the setup's spaces are the options' (`x_search_core`'s `spaces`).
pub(crate) fn preliminary_spaces(
    batch: &PreparedQueryBatch,
    parameters: &[ContextParameters],
    options: &ResolvedOptions,
    numbers: SubjectNumbers,
    database: bool,
) -> Result<Option<Vec<i64>>> {
    if !database || numbers.gumbel_total > 0 || options.effective_search_space != 0 {
        return Ok(None);
    }
    let gapped = gapped_karlin(options);
    let spaces = batch
        .contexts
        .iter()
        .zip(parameters)
        .map(|(c, p)| {
            let ka = if options.gapped { gapped } else { p.ungapped };
            crate::stats::setup::preliminary_option_space(
                if p.valid { c.length as i64 } else { 0 },
                numbers.gumbel_total,
                numbers.request_count,
                options.db_length,
                0,
                false,
                &ka,
                options.program().upper,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(spaces))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1001-1024 (cited at
// `stats::setup::subject_eff_lengths`)
/// The parameters of the contexts of a query chunk `batch` for the subject of `subject_length`
/// residues, where NCBI's preliminary search computes the statistics of each subject
/// (`Subjects::per_subject`, a database length of 0): the chunk's preliminary parameters
/// `preliminary` with the effective lengths of the subject's length (`numbers.count` sequences,
/// `-dbsize`, the options' spaces of the preliminary stage), the hit-saving and word cutoffs of
/// the subject's length without the stretch of a composition-based statistics. `option_spaces` are
/// the whole query's options' spaces of the preliminary stage (`PreparedQueryBatch::
/// subject_option_spaces`), which a chunk reads by its own context index.
pub(crate) fn subject_update(
    batch: &PreparedQueryBatch,
    preliminary: &[ContextParameters],
    options: &ResolvedOptions,
    numbers: SubjectNumbers,
    subject_length: usize,
    option_spaces: &[i64],
) -> Result<Vec<ContextParameters>> {
    let mut parameters = preliminary.to_vec();
    effective_lengths(
        batch,
        &mut parameters,
        options,
        subject_length as i64,
        numbers.count,
        Some(option_spaces),
    )?;
    let gumbel = gumbel_block(options, numbers.gumbel_total);
    let length = i32::try_from(subject_length)?;
    context_cutoffs(
        batch,
        &mut parameters,
        options,
        gumbel.as_ref(),
        length,
        length,
        false,
    )?;
    Ok(parameters)
}

/// The gap decay rate of the hit-saving and word cutoffs (`BlastLinkHSPParametersNew`, cited at
/// `linking::link_parameters`, the one rule of whether BLASTX links), 0 without linking: no sum
/// statistics, or a gapped search whose intron length links nothing (blast_parameters.c:793-805).
fn gap_decay_rate(options: &ResolvedOptions) -> f64 {
    super::linking::link_parameters(options).map_or(0.0, |link| link.gap_decay_rate)
}

/// `CUTOFF_E_BLASTX`, the e-value of an ungapped search's word cutoff (`s_GetCutoffEvalue`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:78
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:143-144
const CUTOFF_E_BLASTX: f64 = 1.0;

/// The word cutoff of a valid context of `length` residues (`BlastInitialWordParametersUpdate`,
/// cited above `SubjectNumbers`): the gap trigger of its ungapped block, below it an ungapped
/// search's `BLAST_Cutoffs` for `CUTOFF_E_BLASTX` with the gap decay `decay` in the search space
/// of the shorter of the context and the subject length `average` times `average`, at most the
/// hit-saving cutoff `p.hit_cutoff_max`. NCBI casts the two lengths to `Uint8`; both are positive
/// here (a context's length and a subject length that fit `Int4`), so their `Int8` product is
/// the same.
fn word_cutoff(
    options: &ResolvedOptions,
    length: usize,
    average: i32,
    decay: f64,
    p: &ContextParameters,
) -> i32 {
    let trigger = gap_trigger_raw_score(options.gap_trigger, &p.ungapped);
    if options.gapped {
        trigger
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:356-363
        let searchsp = i64::from((length as i32).min(average)) * i64::from(average);
        cutoff_score(
            CUTOFF_E_BLASTX,
            &p.ungapped,
            searchsp,
            Some(decay),
            CutoffConversion::FiniteQuotient,
        )
        .min(trigger)
    }
    .min(p.hit_cutoff_max)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:396-416
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1448-1455 (cited at
// `search::compute_subject_with`): the link cutoffs read it (`linking::link_cutoffs`).
// NCBI's query batch is the reader batch (`reader_batch.rs`): the statistics of each subject update
// the word cutoffs of all its contexts, also those of the units LOSAT searches apart.
/// The smallest word cutoff over the contexts `contexts` of a reader batch (`ReaderBatchShape::
/// batch_contexts`: length, validity, ungapped block) for the subject of `subject_length` residues,
/// where NCBI computes the statistics of each subject: each context's effective lengths with its
/// preliminary option space (`stats::setup::preliminary_option_space`: the spaces `x_search_core`
/// gives the batch), its hit-saving cutoff (an ungapped search has no Gumbel block) and its word
/// cutoff, as `subject_update` computes them. `None` without a valid context. Only an ungapped
/// search with sum statistics reads it.
pub(crate) fn subject_word_cutoff_min(
    contexts: &[(usize, bool, KarlinParams)],
    options: &ResolvedOptions,
    numbers: SubjectNumbers,
    subject_length: usize,
) -> Result<Option<i32>> {
    let decay = gap_decay_rate(options);
    let average = i32::try_from(subject_length)?;
    let mut minimum: Option<i32> = None;
    for &(query_length, valid, ungapped) in contexts {
        if !valid || query_length == 0 {
            continue;
        }
        let option_space = crate::stats::setup::preliminary_option_space(
            query_length as i64,
            numbers.gumbel_total,
            numbers.request_count,
            options.db_length,
            options.effective_search_space,
            false,
            &ungapped,
            options.program().upper,
        )?;
        let lengths = crate::stats::setup::subject_eff_lengths(
            query_length as i64,
            subject_length as i64,
            numbers.count,
            options.db_length,
            option_space,
            false,
            &ungapped,
            options.program().upper,
        )?;
        let p = ContextParameters {
            valid,
            ungapped,
            search_space: lengths.eff_searchsp,
            length_adjustment: lengths.length_adjustment,
            word_cutoff: 0,
            word_xdrop_init: 0,
            word_xdrop: 0,
            hit_cutoff: 0,
            hit_cutoff_max: cutoff_score(
                options.evalue,
                &ungapped,
                lengths.eff_searchsp,
                None,
                CutoffConversion::FiniteQuotient,
            ),
        };
        let cutoff = word_cutoff(options, query_length, average, decay, &p);
        minimum = Some(minimum.map_or(cutoff, |m| m.min(cutoff)));
    }
    Ok(minimum)
}

/// The cutoffs of each valid context (`BlastHitSavingParametersNew`/`Update`,
/// `BlastInitialWordParametersNew`/`Update`, cited above `SubjectNumbers`) after its effective
/// lengths: the hit-saving cutoff from `evalue` (stretched 5 times for the Gumbel block with
/// `stretch`, a composition-based statistics above 1) at the subject length `n` (Spouge's, and the
/// sum statistics' search space), and the word cutoff, which an ungapped search computes at the
/// subject length `average`. The setup passes the shortest or average subject length, the
/// statistics of each subject the subject's own length and no stretch (`subject_update`).
fn context_cutoffs(
    batch: &PreparedQueryBatch,
    parameters: &mut [ContextParameters],
    options: &ResolvedOptions,
    gumbel: Option<&BlastGumbelBlk>,
    n: i32,
    average: i32,
    stretch: bool,
) -> Result<()> {
    let gapped = gapped_karlin(options);
    let decay = gap_decay_rate(options);
    // The average context length is that of the query batch (a reader batch searched in several
    // units has one: `reader_batch.rs`, c++/src/algo/blast/core/blast_parameters.c:954-957).
    let avg_q = i32::try_from(batch.batch_concatenated_length() / batch.batch_context_count())?;
    // The sum statistics' cutoffs carry their e-value from one valid context to the next
    // (`stats::cutoffs::SumStatsCutoffs`, blast_parameters.c:949-976), over the query batch: a
    // unit after a valid context of its reader batch starts where those contexts left the e-value
    // (`PreparedQueryBatch::batch_valid_before`). They all read the gapped block and the space
    // of `avg_q` and `n`, so one step gives the state after any number of them
    // (`SumStatsCutoffs::next_valid_context`).
    let mut sum_stats = SumStatsCutoffs::new(avg_q, n, decay, CutoffConversion::FiniteQuotient);
    if options.gapped && decay != 0.0 && batch.batch_valid_before {
        sum_stats.next_valid_context(&gapped);
    }
    for (c, p) in batch.contexts.iter().zip(parameters) {
        if !p.valid {
            continue;
        }
        let ka = if options.gapped { gapped } else { p.ungapped };
        p.hit_cutoff_max = if let Some(gbp) = gumbel {
            blast_spouge_etos(
                if stretch {
                    5.0 * options.evalue
                } else {
                    options.evalue
                },
                &ka,
                gbp,
                c.length as i32,
                n,
            )
        } else {
            cutoff_score(
                options.evalue,
                &ka,
                p.search_space,
                None,
                CutoffConversion::FiniteQuotient,
            )
        };
        p.hit_cutoff = p.hit_cutoff_max;
        if options.gapped && decay != 0.0 {
            p.hit_cutoff = p.hit_cutoff.min(sum_stats.next_valid_context(&ka));
        }
        p.word_cutoff = word_cutoff(options, c.length, average, decay, p);
        p.word_xdrop_init = x_drop_raw_score(options.x_dropoff, &p.ungapped, 1.0);
        p.word_xdrop = if p.word_xdrop_init == 0 {
            p.word_cutoff
        } else {
            p.word_xdrop_init
        };
    }
    Ok(())
}

/// The preliminary and final gapped X-drops in raw scores (`gap_x_dropoff`,
/// `gap_x_dropoff_final`) of the search's X-drops in bits (`-xdrop_gap`, `-xdrop_gap_final` or
/// BLASTX's defaults) at the gapped Lambda `gapped` (`gapped_karlin`; the contexts share it, the
/// smallest one), converted as the protein programs convert them (`core::x_drop`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:454-463 (cited at
/// `core::x_drop::gapped_x_dropoffs`)
pub(crate) fn gapped_x_drops(options: &ResolvedOptions, gapped: &KarlinParams) -> (i32, i32) {
    crate::core::x_drop::gapped_x_dropoffs(
        options.gap_x_dropoff,
        options.gap_x_dropoff_final,
        gapped.lambda,
    )
}

/// The refusal, before the search, of a gapped X-drop that NCBI's kernels cannot size
/// (`core::x_drop::check_gapped_x_dropoffs`); an ungapped search has none.
pub(crate) fn check_gapped_x_drops(options: &ResolvedOptions) -> Result<()> {
    if !options.gapped {
        return Ok(());
    }
    // The final X-drop reaches the unscaled `ALIGN_EX` only on the ordinary traceback
    // (composition 0, `runtime.rs`).
    crate::core::x_drop::check_gapped_x_dropoffs(
        "BLASTX",
        (options.gap_x_dropoff, options.gap_x_dropoff_final),
        gapped_x_drops(options, &gapped_karlin(options)),
        (options.gap_open, options.gap_extend),
        options.composition == 0,
    )
}
