//! The values of the options that blast_formatter's reports print, one record for the run
//! (`FormatterOptions`): NCBI's `CBlastOptions` of the first search object (`GetSearchOptions`
//! builds them once, `BuiltOptions`) as the formatters read them. The pairwise report's epilog
//! and tables read the scores, gap costs and modes; every report reads the matrix name, which an
//! `rpsblast` or `rpstblastn` object's `GetResultSet` sets for the objects after it
//! (`archive::Run::print`).

use crate::report::archive::read::options::{BuiltOptions, EProgram, Setter};

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:60-178
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:136-219
// NCBI reference (598d8ae6): c++/src/algo/blast/api/disc_nucl_options.cpp:59-89
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:61-99
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:788-813
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:58-157
/// The values of the options that the reports print (`CBlastFormat::PrintEpilog`,
/// blast_format.cpp:2249-2292, and the formatters' matrix): those of the factory options of the
/// program (`CBlastOptionsFactory::Create`: `SetTraditionalBlastnDefaults` for blastn,
/// megablast's, dc-megablast's, rpsblast's handle, which sets no word threshold, and the protein
/// handles of blastp, psiblast, deltablast and PHI-BLAST), then the archive's calls on them in
/// order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct FormatterOptions {
    pub(super) reward: i32,
    pub(super) penalty: i32,
    pub(super) gap_open: i32,
    pub(super) gap_extend: i32,
    pub(super) word_threshold: i32,
    pub(super) window_size: i32,
    pub(super) gapped: bool,
    pub(super) sum_statistics: bool,
    /// `GetMatrixName()` (`matrix_name`), `None` for NULL; `GetResultSet` of an `rpsblast` or
    /// `rpstblastn` object sets it to its search statistics' matrix (`SetMatrixName`).
    pub(super) matrix_name: Option<Vec<u8>>,
}

impl FormatterOptions {
    /// The values of the options `options` of the first search object.
    pub(super) fn of(options: &BuiltOptions) -> Self {
        let (reward, penalty, gap_open, gap_extend, word_threshold, window_size) =
            match options.program {
                EProgram::Blastn => (2, -3, 5, 2, 0, 0),
                EProgram::Megablast => (1, -2, 0, 0, 0, 0),
                EProgram::DiscMegablast => (2, -3, 5, 2, 0, 40),
                EProgram::RpsBlast => (1, -3, 11, 1, 0, 40),
                _ => (1, -3, 11, 1, 11, 40),
            };
        let mut values = FormatterOptions {
            reward,
            penalty,
            gap_open,
            gap_extend,
            word_threshold,
            window_size,
            gapped: true,
            sum_statistics: false,
            matrix_name: matrix_name(options),
        };
        for setter in &options.setters {
            match *setter {
                Setter::MatchReward(value) => values.reward = value,
                Setter::MismatchPenalty(value) => values.penalty = value,
                Setter::GapOpeningCost(value) => values.gap_open = value,
                Setter::GapExtensionCost(value) => values.gap_extend = value,
                Setter::WordThreshold(value) => values.word_threshold = value,
                Setter::WindowSize(value) => values.window_size = value,
                Setter::GappedMode(gapped) => values.gapped = gapped,
                Setter::SumStatisticsMode(on) => values.sum_statistics = on,
                _ => {}
            }
        }
        values
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:121
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_nucl_options.cpp:200
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.cpp:58
/// `GetMatrixName()` of the options `options`: the last `MatrixName` setter, else the factory's:
/// BLOSUM62 for the protein handles; none for the nucleotide handles (`SetMatrixName` with a null
/// name leaves the scoring options' NULL) and for RPS-BLAST's (whose defaults set no matrix).
fn matrix_name(options: &BuiltOptions) -> Option<Vec<u8>> {
    let set = options
        .setters
        .iter()
        .rev()
        .find_map(|setter| match setter {
            Setter::MatrixName(name) => Some(name.as_bytes().to_vec()),
            _ => None,
        });
    set.or_else(|| {
        matches!(
            options.program,
            EProgram::Blastp
                | EProgram::Blastx
                | EProgram::Tblastn
                | EProgram::Tblastx
                | EProgram::PsiBlast
                | EProgram::PsiTblastn
                | EProgram::PhiBlastp
                | EProgram::DeltaBlast
        )
        .then(|| b"BLOSUM62".to_vec())
    })
}

#[cfg(test)]
#[path = "formatter_options_tests.rs"]
mod tests;
