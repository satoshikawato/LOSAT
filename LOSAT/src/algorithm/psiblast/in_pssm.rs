//! psiblast `-in_pssm F` (series PS, items PS4a and PS4b): the warnings of `PsiBlastSetupScoreBlock` for psiblast's
//! options (the frequency ratios all zero, the composition-based statistics reset to mode 1) and LOSAT's refusals of
//! a psiblast `-in_pssm` search that the PSSM search does not run yet (DECISIONS 11 of line PS). The PSSM file is
//! read by `algorithm::pssm_search::in_pssm` (psiblast's and tblastn's reader).

use std::sync::Arc;

use anyhow::{bail, Result};

use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::pssm_search::PsiSearch;
use crate::blastinput::app::CompositionMode;

/// The warnings of `PsiBlastSetupScoreBlock` for a PSSM search in the composition-based mode
/// `mode` with the matrix `matrix_name` (the frequency ratios all zero when the statistics are on;
/// modes 2 and 3 reset to 1), which the query's warning line carries. None for a PSSM that
/// `CPsiBlastImpl`'s constructor stops before (`CPsiBlastValidate::Pssm`, the PSSM engine's
/// scores).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:179-207
/// ```c++
///     _ASSERT(options->GetCompositionBasedStats() < eNumCompoAdjustModes);
///     // the message below is meaningless for deltablast
///     if (options->GetProgram() != eDeltaBlast &&
///         (options->GetCompositionBasedStats() != eNoCompositionBasedStats) &&
///         freq_ratios_all_zeros) {
///         ostringstream os;
///         os << "Frequency ratios for PSSM are all zeros, frequency ratios for ";
///         os << options->GetMatrixName() << " will be used during traceback ";
///         os << "in composition based statistics";
///         ...
///     }
///
///     if (options->GetCompositionBasedStats() > eCompositionBasedStats) {
///         // ugly, but necessary
///         const_cast<CBlastOptions*>(&*options)
///             ->SetCompositionBasedStats(eCompositionBasedStats);
///         ostringstream os;
///         os << "Composition-based score adjustment conditioned on "
///            << "sequence properties and unconditional composition-based score "
///            << "adjustment is not supported with PSSMs, resetting to default "
///            << "value of standard composition-based statistics";
///         ...
///     }
/// ```
pub(crate) fn setup_messages(
    search: &PsiSearch,
    mode: CompositionMode,
    matrix_name: &str,
) -> Vec<String> {
    let mut messages = Vec::new();
    if search.validation_error.is_some() || search.scores_error.is_some() {
        return messages;
    }
    if mode != CompositionMode::NoCompositionBasedStats
        && crate::core::composition_adjustment::compo_scale::freq_ratios_all_zeros(
            &search.freq_ratios,
        )
    {
        messages.push(format!(
            "Frequency ratios for PSSM are all zeros, frequency ratios for {matrix_name} will be used during traceback in composition based statistics"
        ));
    }
    if position_based_mode(mode) != mode {
        messages.push(
            "Composition-based score adjustment conditioned on sequence properties and unconditional composition-based score adjustment is not supported with PSSMs, resetting to default value of standard composition-based statistics"
                .to_string(),
        );
    }
    messages
}

/// The composition-based mode of a PSSM search for the mode `mode` of the options: 2 and 3 are
/// 1 (`PsiBlastSetupScoreBlock` sets it in the options, so the preliminary cutoff's stretch and
/// the redo read 1; `setup_messages`).
pub(crate) fn position_based_mode(mode: CompositionMode) -> CompositionMode {
    match mode {
        CompositionMode::CompositionMatrixAdjust | CompositionMode::CompoForceFullMatrixAdjust => {
            CompositionMode::CompositionBasedStats
        }
        other => other,
    }
}

/// LOSAT's refusals of a psiblast `-in_pssm` search that the PSSM search does not run yet
/// (DECISIONS 11 of line PS), checked when the options have been processed, before the search.
/// A run that NCBI stops before the search (`PsiSearch::stop_before_search`: the search strategy's
/// "Remote PSI-BL2SEQ" with `-subject`, NCBI's check of the PSSM, the PSSM engine's scores) writes
/// the prolog and NCBI's error, whatever its format or query (line RG F1b-2, review 4 R5-6): these
/// refusals of the search's results do not apply to it.
pub(crate) fn check_in_pssm_options(args: &ResolvedBlastpArgs, psi: &Arc<PsiSearch>) -> Result<()> {
    if psi
        .stop_before_search(args.strategy.is_some(), args.db.is_some())
        .is_some()
    {
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:191-198 (the SEG mask
    // of the query, then `PsiBlastSetupScoreBlock`)
    // The query's SEG mask changes the PSSM search's hits (DESIGN_ps4 §1.6).
    if args.seg.params().is_some() || args.soft_masking {
        bail!("-seg and -soft_masking true with -in_pssm are not supported by LOSAT's PSIBLAST");
    }
    // The outfmt 0, 6 and 7 reports (the others are item PS15b's).
    let choice = crate::blastinput::app::parse_formatting_string(&args.outfmt)?;
    if !matches!(choice.number, 0 | 6 | 7) {
        bail!(
            "output format {} with -in_pssm is not supported by LOSAT's PSIBLAST yet",
            choice.number
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:91-94
    // ```c++
    //     else if (Blast_ProgramIsPsiBlast(program))
    //         retval = 10000;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:191-196 (the
    // chunks of a split query are set up without the PSSM)
    // A PSSM query that NCBI splits is searched with BLOSUM62 in its preliminary stage. Sizes
    // that NCBI cannot convert split nothing: they stop the search's query splitter, before any
    // split (`query_batch::read_psiblast_query_batches`).
    let Ok(sizes) = crate::algorithm::blastp::protein_query_split::protein_query_split_sizes(10000)
    else {
        return Ok(());
    };
    if crate::algorithm::blastp::protein_query_split::split_protein_batch(
        &[psi.query.seq().len()],
        sizes,
    )
    .is_some()
    {
        bail!(
            "a PSSM query of {} letters, which NCBI BLAST+ splits into chunks searched without the PSSM, is not supported by LOSAT's PSIBLAST",
            psi.query.seq().len()
        );
    }
    Ok(())
}
