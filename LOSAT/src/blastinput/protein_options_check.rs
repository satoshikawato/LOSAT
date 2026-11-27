//! NCBI's `BLAST_ValidateOptions` for blastp, tblastn and psiblast (`validate_protein_options`):
//! the option checks over the protein matrix tables of `stats::protein_options`.

use crate::stats::protein_options::{
    karlin_blk_gapped_load_from_tables, print_allowed_values, print_matrix_message, TableStatus,
};

/// The options that NCBI's `BLAST_ValidateOptions` checks for blastp and tblastn, as their
/// option handlers set them.
pub struct ProteinOptionsCheck<'a> {
    /// Not `-ungapped`.
    pub gapped: bool,
    pub gap_open: i32,
    pub gap_extend: i32,
    /// The matrix name as typed.
    pub matrix_name: &'a [u8],
    pub threshold: f64,
    pub word_size: i32,
    /// Whether the compressed-alphabet lookup table is used (`eCompressedAaLookupTable`).
    pub compressed_lookup: bool,
    /// The hit list size (`hit_options->hitlist_size`).
    pub hitlist_size: usize,
    pub evalue: f64,
    /// The program is psiblast (`eBlastTypePsiBlast`) rather than blastp or tblastn.
    pub psiblast: bool,
    /// The program is tblastn with `-in_pssm` (PSI-TBLASTN, `eBlastTypePsiTblastn`): neither blastp
    /// nor tblastn nor psiblast.
    pub psi_tblastn: bool,
    /// The program is PHI-BLAST (psiblast `-phi_pattern`, `eBlastTypePhiBlastp`): its scoring is
    /// checked as psiblast's (`psiblast` is set too), and it has no lookup table or ungapped
    /// extension to check.
    pub phi: bool,
    /// The ungapped X-drop in bits (`word_options->x_dropoff`) where the command line gives it
    /// (`-xdrop_ungap` of blastp, tblastn and rpsblast); `None`: the program's default, which is
    /// above 0.
    pub x_dropoff: Option<f64>,
    /// `-culling_limit`, `-best_hit_overhang` and `-best_hit_score_edge` as the options hold them
    /// (an imported search strategy's values unchecked: `app::check_hsp_filtering`).
    pub culling_limit: Option<i32>,
    pub best_hit_overhang: Option<f64>,
    pub best_hit_score_edge: Option<f64>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:409,3687
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3266
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4059-4061
/// LOSAT's limit of `-evalue` for blastp, psiblast, rpsblast and tblastn (`program`, in
/// capitals), checked where NCBI starts the search (decision D12 of
/// docs/evidence/losat_web_e2e/AUTHORITY.md):
///
/// - An infinite value (`+inf`, `1e999`): NCBI's score cutoff of an e-value is the logarithm of
///   0 converted to an `Int4`, a conversion out of range. Refused in every mode.
/// - Exactly DBL_MAX (1.7976931348623157e308) with the composition-based redo
///   (`composition_redo`: modes 1 to 3 after IDENTITY's reset): a list whose HSPs all fail the
///   redo keeps the best e-value DBL_MAX, which passes `-evalue`, so NCBI puts the empty list
///   into the hit list and crashes reading its first HSP (SIGSEGV) when the search keeps enough
///   HSPs; LOSAT cannot tell beforehand. Without the redo (`-comp_based_stats 0`, IDENTITY,
///   `-ungapped`) NCBI never inserts an empty list and runs: NCBI BLAST+ 2.17.0 gives the rows
///   of 1.7976931348623156e308 there (blastp, psiblast and tblastn, gapped and ungapped, line
///   O4c's review D1), so LOSAT searches it.
pub(crate) fn check_evalue_limit(
    evalue: f64,
    composition_redo: bool,
    program: &str,
) -> anyhow::Result<()> {
    let ncbi_program = program.to_ascii_lowercase();
    if evalue > f64::MAX {
        anyhow::bail!(
            "an infinite -evalue ({evalue}), whose score cutoff NCBI BLAST+'s {ncbi_program} cannot hold in an integer, is not supported by LOSAT's {program}"
        );
    }
    if evalue == f64::MAX && composition_redo {
        anyhow::bail!(
            "an -evalue of DBL_MAX ({evalue}) with composition-based statistics (-comp_based_stats 1 to 3), with which NCBI BLAST+'s {ncbi_program} crashes on some inputs, is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// NCBI's `BLAST_ValidateOptions` for blastp and tblastn, in its order, with NCBI's
/// messages (`CInputException`): the scoring options, the lookup table options, the
/// ungapped X-drop (`-xdrop_ungap`), the hit saving options and the word size of the
/// identity matrix. LOSAT has no greedy or Smith-Waterman extension option or cutoff score
/// option for these programs, so NCBI's other checks cannot fail.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1759-1776
pub fn validate_protein_options(check: &ProteinOptionsCheck<'_>) -> anyhow::Result<()> {
    use crate::blastinput::app::{options_error, options_error_bytes};
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:910-936 (the scoring
    // check of BLAST_ValidateOptions; blastp and tblastn may use IDENTITY)
    // psiblast and PSI-TBLASTN are neither blastp nor tblastn: their matrices are the standard
    // ones (`std_matrix_only`).
    let std_matrix_only = check.psiblast || check.psi_tblastn;
    if check.gapped {
        match karlin_blk_gapped_load_from_tables(
            check.gap_open,
            check.gap_extend,
            check.matrix_name,
            std_matrix_only,
        ) {
            Ok(_) => {}
            Err(TableStatus::UnknownMatrix) => {
                return Err(options_error_bytes(&print_matrix_message(
                    check.matrix_name,
                    std_matrix_only,
                )))
            }
            Err(TableStatus::UnknownGapCosts) => {
                return Err(options_error_bytes(&print_allowed_values(
                    check.matrix_name,
                    check.gap_open,
                    check.gap_extend,
                )))
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1298-1300
    // PHI-BLAST skips the lookup table's checks and the initial word's (blast_options.c:609-612):
    // the hit saving checks follow.
    if check.phi {
        return check_hit_saving(check);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1395 (the lookup
    // check of BLAST_ValidateOptions for blastp and tblastn)
    // (`threshold <= 0` is false for NaN, as in C.)
    if check.threshold <= 0.0 {
        return Err(options_error("Non-zero threshold required"));
    }
    if check.word_size <= 0 {
        return Err(options_error("Word-size must be greater than zero"));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1350-1357
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1358-1363
    // The checks of the compressed lookup table below are blastp's, tblastn's and blastx's
    // only (blast_options.c:1376-1378); PSI-TBLASTN takes the `else` branch.
    if check.psiblast {
        if check.word_size > 4 {
            return Err(options_error("Word-size must be less than 5 for psiblast"));
        }
    } else if check.psi_tblastn {
        if check.word_size > 4 {
            return Err(options_error(
                "Word-size must be less than 6 for protein comparison",
            ));
        }
    } else if check.word_size > 7 {
        return Err(options_error(
            "Word-size must be less than 8 for a tblastn, blastp or blastx search",
        ));
    }
    let blastp_tblastn = !check.psiblast && !check.psi_tblastn;
    if blastp_tblastn && check.word_size > 5 && !check.compressed_lookup {
        return Err(options_error(
            "Blastp, Blastx or Tblastn with word size > 5 requires a compressed alphabet lookup table",
        ));
    }
    if blastp_tblastn && check.compressed_lookup && !(5..=7).contains(&check.word_size) {
        return Err(options_error(
            "Compressed alphabet lookup table requires word size 5, 6 or 7",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:608-618 (the initial
    // word check of BLAST_ValidateOptions, after the lookup table's and before the hit saving
    // options')
    if check.x_dropoff.is_some_and(|x_dropoff| x_dropoff <= 0.0) {
        return Err(options_error("x_dropoff must be greater than zero"));
    }
    check_hit_saving(check)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1783-1795
    if check.word_size > 5 && check.matrix_name.eq_ignore_ascii_case(b"IDENTITY") {
        return Err(options_error(
            "Word size larger than 5 is not supported for the identity scoring matrix",
        ));
    }
    Ok(())
}

/// The hit saving checks of `BLAST_ValidateOptions` (`validate_protein_options`).
fn check_hit_saving(check: &ProteinOptionsCheck<'_>) -> anyhow::Result<()> {
    use crate::blastinput::app::options_error;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1511-1517 (the hit saving
    // check of the hit list size, before the e-value's)
    if check.hitlist_size < 1 {
        return Err(options_error("No hits are being saved"));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1521 (the hit
    // saving check; LOSAT's blastp and tblastn have no cutoff score option)
    if check.evalue <= 0.0 {
        return Err(options_error(
            "expect value or cutoff score must be greater than zero",
        ));
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1543-1548 (the hit
    // saving check of the HSP filtering options, after the culling limit's, which the options'
    // setter keeps at 0 or above)
    crate::blastinput::app::check_hsp_filtering(
        check.culling_limit,
        check.best_hit_overhang,
        check.best_hit_score_edge,
    )?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4886-4896
/// LOSAT's refusal of IDENTITY with the compressed lookup table: word size 5, which `-task
/// blastp-fast` and `tblastn-fast` set by default (larger sizes are NCBI's error,
/// `validate_protein_options`). NCBI's compressed lookup table scores its words with a matrix
/// made from the search matrix's frequency ratios, as LOSAT's does
/// (`tblastx/lookup/compressed.rs`, `build_compressed_score_matrix`). IDENTITY has none, and NCBI
/// BLAST+ 2.17.0's blastp and tblastn die of SIGSEGV (`s_BlastCompressedAaScanSubject`) when they
/// search a query batch with it, subjects without letters too. A run that NCBI stops at its first
/// query batch (a batch of queries without letters: NCBI's error) does not reach the search, so
/// the programs check this where they search their first batch; the options-only check of the Web
/// (`check_options`) checks it without the inputs.
pub fn check_identity_compressed_lookup(
    matrix_name: impl AsRef<[u8]>,
    compressed_lookup: bool,
    word_size: i32,
    ncbi_program: &str,
    program: &str,
) -> anyhow::Result<()> {
    let matrix_name = matrix_name.as_ref();
    if compressed_lookup && matrix_name.eq_ignore_ascii_case(b"IDENTITY") {
        anyhow::bail!(
            "-matrix {} with the compressed lookup table (-word_size {word_size}), with which NCBI BLAST+'s {ncbi_program} crashes, is not supported by LOSAT's {program}",
            String::from_utf8_lossy(matrix_name)
        );
    }
    Ok(())
}

#[cfg(test)]
#[path = "protein_options_check_tests.rs"]
mod tests;
