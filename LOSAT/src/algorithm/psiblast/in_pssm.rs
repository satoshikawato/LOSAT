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
/// scores), and for one whose query set-up fails (`PsiSearch::setup_error`): that set-up is the
/// first step of the search's preparation, before `PsiBlastSetupScoreBlock`
/// (blast_aux_priv.cpp:149-151, 197).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:179-207
pub(crate) fn setup_messages(
    search: &PsiSearch,
    mode: CompositionMode,
    matrix_name: &str,
) -> Vec<String> {
    if search.validation_error.is_some()
        || search.scores_error.is_some()
        || search.setup_error.is_some()
    {
        return Vec::new();
    }
    pssm_setup_messages(&search.freq_ratios, mode, matrix_name)
}

/// `setup_messages` of a PSSM whose frequency ratios are `freq_ratios`: every column of the
/// PSSM, also those that the search does not read (a later round of a `-query_loc` range reads
/// the record's PSSM only up to the range's length, `PsiSearch::from_engine`; NCBI's check loops
/// over `freq_ratios->GetCols()`, the record's columns).
pub(crate) fn pssm_setup_messages(
    freq_ratios: &[[f64; crate::core::blast_util::BLASTAA_SIZE]],
    mode: CompositionMode,
    matrix_name: &str,
) -> Vec<String> {
    let mut messages = Vec::new();
    if mode != CompositionMode::NoCompositionBasedStats
        && crate::core::composition_adjustment::compo_scale::freq_ratios_all_zeros(freq_ratios)
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
/// refusals of the search's results do not apply to it. `outfmts` are the numbers of the requested
/// formats: the command line's one, or each of the records' entries (the library API), checked as
/// the command line checks it alone.
pub(crate) fn check_in_pssm_options(
    args: &ResolvedBlastpArgs,
    psi: &Arc<PsiSearch>,
    outfmts: &[i32],
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:506, 537 (cited at
    // `PsiSearch::write_annot_posts`: the prolog, then the query's Bioseq joins the scope)
    // An exception of NCBI's object manager when the query joins the scope (`PsiSearch::
    // annot_abort`) ends the run after the prolog of its format, before any report: NCBI BLAST+
    // 2.17.0 `om/` cases t04 and k02 with -outfmt 10 (nothing on the standard output) and 20 (its
    // header line). No refusal of the reports below applies to it.
    if psi.annot.abort.is_some() {
        return Ok(());
    }
    // The XML and JSON query definition of a query whose annotations NCBI reads for it
    // (`in_pssm::refuse_annot_defline`), in every run that writes a report and for each
    // requested format.
    for &number in outfmts {
        crate::algorithm::pssm_search::in_pssm::refuse_annot_defline(psi, number, "PSIBLAST")?;
    }
    if psi
        .stop_before_search(args.strategy.is_some(), args.db.is_some())
        .is_some()
    {
        return Ok(());
    }
    // The pairwise report, the tabular reports 6, 7, 10 and 20, BLAST XML, the Seq-annots of
    // outfmt 8 and 9 and the archive (the others are item PS15b-2's).
    if let Some(number) = outfmts
        .iter()
        .find(|number| !matches!(number, 0 | 5..=11 | 20))
    {
        let option = if psi.from_msa { "-in_msa" } else { "-in_pssm" };
        bail!("output format {number} with {option} is not supported by LOSAT's PSIBLAST yet");
    }
    // The archive's request holds the PSSM that the search holds (`PsiSearch::archive_pssm`):
    // `-in_pssm`'s object, whose members the XML and ASN.1 binary streams do not all write yet
    // (`report::archive::pssm_unwritten`; the later rounds' PSSMs keep its query). Lifted where
    // those streams write them (line RG item F3's encoder, for XML).
    if outfmts.contains(&crate::blastinput::app::output_format::ARCHIVE_FORMAT) {
        if let (Some(handle), Some(pssm)) = (args.archive.as_ref(), psi.archive_pssm.as_deref()) {
            // NCBI's stream stops at the matrix name that a PSSM of frequency ratios only lacks,
            // in every encoding that writes (`InPssmFile::matrix_name_unwritable`, which
            // `-out_pssm` refuses too; a value of `ARCHIVE_FORMAT` that NCBI writes nothing for
            // has no stream to stop). With `-subject` NCBI's archive stops earlier, at the
            // PSSM round's empty database reference (rc 3), which LOSAT answers as it does.
            if psi.archive_matrix_name_unwritable
                && args.db.is_some()
                && handle.encoding() != crate::report::archive::ArchiveEncoding::Nothing
            {
                bail!(
                    "{}",
                    super::in_pssm_files::matrix_name_refusal("in the archive of -outfmt 11")
                );
            }
            if let Some(what) = crate::report::archive::pssm_unwritten(handle.encoding(), pssm) {
                let stream = match handle.encoding() {
                    crate::report::archive::ArchiveEncoding::Xml => "XML (ARCHIVE_FORMAT=xml)",
                    crate::report::archive::ArchiveEncoding::Text => {
                        "ASN.1 text (no ARCHIVE_FORMAT)"
                    }
                    _ => "ASN.1 binary (ARCHIVE_FORMAT=bin)",
                };
                bail!(
                    "-outfmt 11 in {stream} with a PSSM whose members hold {what} is not supported by LOSAT's PSIBLAST yet"
                );
            }
        }
    }
    Ok(())
}
