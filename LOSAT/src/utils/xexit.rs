//! The reports that the end of a command-line run writes (`write`, called by `cli/main_entry.rs`
//! once the program has run): the `LOSAT_X_STATS` counters, the `LOSAT_TIMING` lines, the
//! `LOSAT_X_REACH` file, the `[X_SHADOW]` summaries of the shadow switches and the
//! `[X_BNTBRANGE]` line.
//!
//! Each of them writes only under a variable of `NAMES` or, in a diagnostic build, of
//! `xdiag::NAMES`. `write` tests those variables first and returns at once when none is set, so a
//! run that asks for no report runs none of their code: each report's own test of its variable,
//! spread over the binary, was the last code such a run touched, about 1.1 MiB of the binary's
//! pages, which in a small run raised the most memory the run held (line RG-KB item 9a-4). A
//! report under a new variable adds the variable to `NAMES`; `LOSAT/tests/exit_reports.rs` reads
//! the functions `write` calls, and the functions they call, and fails where one names a `LOSAT_*`
//! variable outside `NAMES`, `xdiag::NAMES` and the switches (a switch read passes).
//!
//! No NCBI counterpart: LOSAT's own diagnostics. NCBI BLAST+ has none of them, and this module
//! changes nothing that NCBI's programs write.

/// No NCBI counterpart: the variables under which a report of `write` writes, sorted. A report
/// read through `xdiag` (a diagnostic build's) is asked for by its variable in `xdiag::NAMES`.
pub(crate) const NAMES: &[&str] = &[
    "LOSAT_TIMING",
    "LOSAT_X_BNTBRANGELOG",
    "LOSAT_X_BXDIAGSHADOW",
    "LOSAT_X_BXSUBJFLOWSHADOW",
    "LOSAT_X_BXUNGAPSENTSHADOW",
    "LOSAT_X_REACH",
    "LOSAT_X_STATS",
];

/// No NCBI counterpart: whether the run asks for a report: a variable of `NAMES` is set, or one of
/// `xdiag::NAMES` in a diagnostic build (which reads them). Each name is one `getenv`, as the
/// reports ask, so the match rules of the platform (case-insensitive names on Windows) are the
/// same; the value does not matter here (a report that wants a given value, or an absolute path,
/// still tests it).
fn asked() -> bool {
    use crate::utils::xdiag;
    NAMES.iter().any(|name| std::env::var_os(name).is_some())
        || xdiag::NAMES
            .iter()
            .any(|name| xdiag::var_os(name).is_some())
}

/// No NCBI counterpart: writes the reports the run asks for, in this order, to the standard error
/// (`LOSAT_X_REACH`'s to its directory); nothing, and none of their code runs, when it asks for
/// none (`asked`). Every LOSAT program that ends here (the search programs, `makeblastdb`,
/// `rpsbproc`) runs it; the others end before.
pub(crate) fn write() {
    if !asked() {
        return;
    }
    // Diagnostic (LOSAT_X_REACH): the sides the conditional switches took (utils/xreach.rs).
    crate::utils::xreach::write_report();
    crate::utils::xstats::print();
    crate::core::composition_adjustment::adjust_scores::x_adjmemo_print_stats();
    crate::core::composition_adjustment::adjust_scores::x_newtonexact_print_stats();
    crate::algorithm::tblastx::x_seed_bucket::print_summary();
    crate::algorithm::tblastx::blast_engine::x_pair_par_print_summary();
    crate::algorithm::tblastx::blast_engine::x_seed_release_print_summary();
    crate::algorithm::tblastx::sum_stats_linking::x_hsp_inplace_print_summary();
    crate::algorithm::tblastx::sum_stats_linking::x_link_sort_print_summary();
    crate::algorithm::tblastx::sum_stats_linking::x_link_slices_print_summary();
    crate::algorithm::tblastn::x_subject_translate_print_summary();
    crate::algorithm::tblastn::x_subject_store_print_summary();
    crate::algorithm::tblastn::x_scan_bf_print_summary();
    crate::algorithm::tblastn::x_scan_byte_print_summary();
    crate::algorithm::tblastn::x_lut_lean_print_summary();
    crate::algorithm::tblastn::x_link_window_print_summary();
    crate::algorithm::tblastn::x_link_ctx_print_summary();
    crate::algorithm::blastn::x_subject_store_print_summary();
    crate::algorithm::blastn::x_query_dust_par_print_summary();
    crate::algorithm::blastn::x_greedy_row_print_summary();
    crate::algorithm::blastn::blast_engine::x_query_prep_par_print_summary();
    crate::algorithm::blastn::blast_engine::x_qgroup_par_print_summary();
    crate::algorithm::blastn::blast_engine::x_tb_range_print_summary();
    crate::algorithm::blastn::blast_engine::x_lut_drop_print_summary();
    crate::algorithm::blastn::blast_engine::x_query_drop_print_summary();
    crate::algorithm::blastp::x_bp_hit_bf_print_summary();
    crate::algorithm::blastp::x_bp_flat_print_summary();
    crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bplkpflat_print_summary();
    crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_clutpairs_print_summary();
    crate::algorithm::blastp::x_bp_lkp_ahead_print_summary();
    crate::algorithm::blastp::x_restricted_gap_print_summary();
    crate::algorithm::blastn::lookup::x_disc_lookup_par::print_summary();
    crate::algorithm::blastn::lookup::x_dc_csr::print_summary();
    crate::algorithm::blastn::lookup::x_mb_lookup_print_summary();
    crate::algorithm::blastx::seed::x_seed_scan::print_shadow_stats();
    crate::utils::x_logclone::print_shadow_stats();
    // Diagnostic (LOSAT_TIMING): the common line, last (utils/xtiming.rs).
    crate::utils::xtiming::print();
}
