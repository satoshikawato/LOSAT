//! tblastn `-in_pssm F` (PSI-TBLASTN; line PS, item PS7, DECISIONS 35): the PSSM read as psiblast
//! reads it (`algorithm::pssm_search::in_pssm`, with tblastn's program name and matrix), its query
//! record, the lookup table of its rows, the program name of its reports, and LOSAT's refusals of
//! what NCBI BLAST+ 2.17.0's tblastn cannot search with a PSSM.
//!
//! The search is tblastn's (the translation of the subjects, the units, chunks and threads, the
//! containment trees, the fence retry, the `-db` stream, the hit lists and the reports): only the
//! source of the scores (`score::PssmScore`), the Karlin-Altschul blocks and the lookup table are
//! the PSSM's, as NCBI's `eBlastTypePsiTblastn` (`PSSM_QUERY_MASK | eBlastTypeTblastn`) takes
//! tblastn's code everywhere except where `Blast_QueryIsPssm` selects the PSSM.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:82-83

use std::path::Path;

use anyhow::{bail, Result};

use super::args::{ResolvedTblastnArgs, TblastnArgs};
use crate::algorithm::pssm_search::in_pssm::{read_in_pssm, PssmReader};
use crate::algorithm::pssm_search::rows::QueryPssm;
use crate::algorithm::pssm_search::PsiSearch;
use crate::algorithm::tblastx::lookup::BlastAaLookupTable;

/// The program line of PSI-TBLASTN's reports (`Blast_ProgramNameFromType` of
/// `eBlastTypePsiTblastn`, upper-cased by the formatter).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:690
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:94
pub(super) const PROGRAM: &str = "PSITBLASTN";

/// Reads `-in_pssm path` for tblastn: psiblast's reader (NCBI's `CPsiBlastArgs`, the same handler
/// for both programs), whose refusals name LOSAT's TBLASTN and whose scores of a PSSM of
/// frequency ratios only are made with the options' matrix (`matrix_name`); `strategy` and
/// `database` are whether the run writes a search strategy and searches `-db`
/// (`PssmReader::strategy`); the messages of NCBI's reader go to `diagnostics`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:125-130
pub(super) fn read(
    path: &Path,
    matrix_name: &str,
    strategy: bool,
    database: bool,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    diagnostics: &mut dyn std::io::Write,
) -> Result<PsiSearch> {
    let (search, _file) = read_in_pssm(
        path,
        PssmReader {
            program: "TBLASTN",
            matrix_name,
            strategy,
            database,
            archive,
        },
        diagnostics,
    )?;
    Ok(search)
}

/// The check of `-db_gencode` with `-in_pssm`: the value is not used (PSI-TBLASTN translates with
/// code 1, `check_options`), and NCBI's command line takes the codes of its list only, without
/// the code 32 that LOSAT's TBLASTN adds for its own translations (`PD-TLOSAN-LOCAL-GENCODE-32`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1011-1015
pub(super) fn check_db_gencode(code: u8) -> Result<()> {
    if code == 32 {
        bail!(
            "-db_gencode 32 with -in_pssm is not supported by LOSAT's TBLASTN (NCBI BLAST+ 2.17.0 does not take the code 32, and its PSI-TBLASTN translates the subjects with code 1)"
        );
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1338-1340
/// Whether `error` is NCBI's error of a `-in_pssm` file that cannot be opened (`cli::inaccessible`),
/// which NCBI raises inside `Run`: the options' handler opens the file there, so its `CATCH_ALL`
/// writes the error archive of outfmt 11 for it (NCBI BLAST+ 2.17.0, `-in_pssm` of a missing file
/// with `-subject` and with `-db`), as for the files that `archive::error_lines` takes as the
/// subjects' (`subject_file`).
pub(super) fn file_error_in_run(error: &anyhow::Error) -> bool {
    error
        .downcast_ref::<crate::cli::NativeError>()
        .is_some_and(|native| {
            native
                .message
                .starts_with(b"Command line argument error: Argument \"in_pssm\".")
        })
}

/// LOSAT's refusals of a PSI-TBLASTN search, after NCBI's checks of the options and the opening of
/// the database, before anything is written: NCBI's crash of `-mt_mode 1` (DECISIONS 35 Q4), the
/// output formats other than 0 to 4, 6, 7 and 10 (item PS15b of line PS; 11 crashes NCBI).
pub(super) fn check_limits(
    args: &TblastnArgs,
    resolved: &ResolvedTblastnArgs,
    psi: &PsiSearch,
    outfmt: i32,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:222-226, 241-270 (cited at
    // `tblastn_app::search_pssm_query`: the query's Bioseq joins the scope before the formatter
    // is made)
    // An exception of NCBI's object manager when the query joins the scope (`PsiSearch::
    // annot_abort`) ends the run before any report: its posts and its error, with outfmt 11 the
    // error archive (NCBI BLAST+ 2.17.0, line RG D3-1 `om/` cases t04, k02 and the other cases of
    // exit status 255 with -outfmt 11). No refusal of the reports below applies to it.
    let index_stops = psi.annot.abort.is_some();
    // The XML and JSON query definition of a query whose annotations NCBI reads for it
    // (`in_pssm::refuse_annot_defline`), in every run that writes a report.
    if !index_stops {
        crate::algorithm::pssm_search::in_pssm::refuse_annot_defline(psi, outfmt, "TBLASTN")?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:159-164
    // NCBI's split by queries reads the query file, which a PSSM search has not opened: NCBI
    // BLAST+ 2.17.0 aborts (probe `ts1.mt1.n2.6`). With `-subject` it runs on one thread.
    if args.mt_mode == Some(1) && resolved.num_threads > 1 && args.db.is_some() {
        bail!(
            "-mt_mode 1 with -num_threads above 1, -db and -in_pssm, with which NCBI BLAST+ 2.17.0's tblastn aborts, is not supported by LOSAT's TBLASTN"
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:791-796 (the tabular
    // formats of tblastn; the pairwise report, `PrintOneResultSet`)
    // The pairwise report and the tabular formats 6, 7 and 10 (their field lists too). A run that
    // NCBI stops before the search (`PsiSearch::stop_before_search`: the search strategy's "Remote
    // PSI-BL2SEQ" with `-subject`, NCBI's check of the PSSM, the PSSM engine's scores) writes the
    // prolog of its format and NCBI's error, with outfmt 11 the error archive (NCBI 2.17.0, line
    // RG F1b-2 fix 1: `t.sub.f*`, `t.scaled.db.f11`, `t.ratios45.f11`; review 4 R5-6): these
    // refusals of the search's reports do not apply to it.
    let searched = !index_stops
        && psi
            .stop_before_search(args.strategy.is_some(), args.db.is_some())
            .is_none();
    // The query-anchored formats 1 to 4 are outfmt 0's report with another alignment section
    // (line RG T3b: NCBI 2.17.0's `-in_pssm` sections, probes `tn.t3b.pssm.*`, are those of the
    // PSSM's query).
    if searched && !matches!(outfmt, 0..=10 | 20) {
        bail!(
            "output format {outfmt} with -in_pssm is not supported by LOSAT's TBLASTN yet{}",
            if outfmt == 11 {
                " (NCBI BLAST+ 2.17.0's tblastn stops with a C++ exception)"
            } else {
                ""
            }
        );
    }
    Ok(())
}

/// The lookup table of a PSSM query: the neighbouring words of the PSSM's rows over the lookup
/// segments, the complement of the query's filtering locations (`-seg`'s masks, also with
/// `-soft_masking true`, which keeps them out of the table only). `query` is the query record's
/// letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:631-637
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:86-99 (cited at
/// `pssm_search::lookup`: a score block with a PSSM indexes the PSSM's rows)
pub(super) fn lookup_table(
    rows: &QueryPssm,
    query: &[u8],
    seg: Option<&crate::utils::seg::SegParams>,
    mask_lowercase: bool,
    word_length: usize,
    threshold: i32,
) -> BlastAaLookupTable {
    let masks =
        super::search_seed::encode_tblastn_lookup_query(query, seg, mask_lowercase).seg_masks;
    let segments: Vec<(i64, i64)> =
        crate::algorithm::tblastx::lookup::build_lookup_segments_for_context(
            0,
            query.len(),
            &masks,
        )
        .into_iter()
        .map(|(left, right)| (i64::from(left), i64::from(right)))
        .collect();
    crate::algorithm::pssm_search::lookup::build_pssm_lookup_in(
        rows,
        word_length,
        threshold,
        &segments,
    )
}

#[cfg(test)]
#[path = "in_pssm_tests.rs"]
mod tests;
