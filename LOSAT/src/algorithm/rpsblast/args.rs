//! Command-line arguments for RPS-BLAST (`losat rpsblast`).

use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
use crate::algorithm::rps_search::engine::RpsSearch;
use crate::blastinput::report_counts::FormattingArgs;
use crate::blastinput::value_parsers::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:70-131
// ```c++
// CRPSBlastAppArgs::CRPSBlastAppArgs()
// {
//     const bool kQueryIsProtein = true;
//     const bool kIsRpsBlast = true;
//     const bool kIsCBS2and3Supported = false;
//     const bool kFilterByDefault = false;
//     ...
//     static const char kDefaultTask[] = "rpsblast";
//     SetTask(kDefaultTask);
//
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs(false, kIsRpsBlast));
//     ...
//     m_StdCmdLineArgs.Reset(new CStdCmdLineArgs);
//     ...
//     arg.Reset(new CGenericSearchArgs(kQueryIsProtein, kIsRpsBlast));
//     ...
//     arg.Reset(new CFilteringArgs(kQueryIsProtein, kFilterByDefault));
//     ...
//     arg.Reset(new CWindowSizeArg);
//     ...
//     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
//     ...
//     m_FormattingArgs.Reset(new CFormattingArgs);
//     ...
//     m_MTArgs.Reset(new CRPSBlastMTArgs());
//     ...
//     m_RemoteArgs.Reset(new CRemoteArgs);
//     ...
//     m_DebugArgs.Reset(new CDebugArgs);
//     ...
//     arg.Reset(new CCompositionBasedStatsArgs(kIsCBS2and3Supported,
//                                              kDfltArgCompBasedStatsRPS,
//                                              cbs_opt_zero_descr));
// ```
// rpsblast's arguments: blastp's search, filtering, query, formatting and composition-based
// statistics arguments without the matrix, the gap costs, the word size and threshold (the
// profile database holds them), with RPS-BLAST's thread arguments. LOSAT runs the search with
// blastp's driver (`BlastpApp::Rpsblast`) and the RPS engine (`algorithm::rps_search`).
// Composition-based statistics (modes 1 to 3, the default 1) run the redo of `rps_search::redo`
// with the database's frequency ratios, or without them where the database has none (NCBI then
// searches with mode 0); the arguments that LOSAT does not run yet (`-searchsp`, item PS11c-2;
// `-use_sw_tback`, PS11c-3; `-html`, PS15a) are refused by the parser (`cli.rs`
// `is_unported_rpsblast_arg`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct RpsblastArgs {
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // ```c
    //     arg_desc.AddDefaultKey(kArgQuery, "input_file",
    //                      "Input file name",
    //                      CArgDescriptions::eInputFile, kDfltArgQuery);
    // ```
    #[arg(long, default_value = "-", value_parser = ncbi_input_path(), value_name = "PATH")]
    pub query: PathBuf,
    #[arg(
        long,
        value_name = "NAME",
        help = "BLAST database name (a profile database)"
    )]
    pub db: Option<String>,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    #[arg(long, value_parser = rpsblast_real, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2036-2042
    // ```c++
    //         arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
    //                                 "Effective length of the database ",
    //                                 CArgDescriptions::eInt8);
    // ```
    /// `-dbsize`: replaces the database's length in every statistic (`x_AdjustDbSize` then keeps
    /// the volumes' own sequence counts, `rps_hooks::rps_search`).
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // ```c++
    //     arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
    //                            CArgDescriptions::eString, m_DefaultOpt);
    // ```
    #[arg(
        long,
        default_value = "1",
        help = "Use composition-based statistics: D or d default (as 0 for RPS-BLAST), 0 or F or f simplified, 1 or T or t composition-based statistics, 2 and 3 as 1 with a larger preliminary e-value"
    )]
    pub comp_based_stats: String,
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,
    #[arg(
        long = "show_gis",
        default_value_t = false,
        help = "Show NCBI GIs in deflines?"
    )]
    pub show_gis: bool,
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// blast_args.cpp:2683-2721, quoted at `FormattingArgs`): the width of the outfmt 0
    /// alignment rows; the sorts are read with NCBI's warnings and sort nothing: NCBI's
    /// rpsblast does not hand them to its formatter (`app/rpsblast_app.cpp:187`,
    /// `SetLineLength` only).
    #[command(flatten)]
    pub formatting: FormattingArgs,
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]", conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = rpsblast_query_coverage, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229
    // ```c++
    //         // ungapped X-drop
    //         // Default values: blastn=20, megablast=10, others=7
    //         arg_desc.AddOptionalKey(kArgUngappedXDropoff, "float_value",
    //                             "X-dropoff value (in bits) for ungapped extensions",
    //                             CArgDescriptions::eDouble);
    //         ...
    //             arg_desc.AddOptionalKey(kArgGappedXDropoff, "float_value",
    //                  "X-dropoff value (in bits) for preliminary gapped extensions",
    //                  CArgDescriptions::eDouble);
    //         ...
    //             arg_desc.AddOptionalKey(kArgFinalGappedXDropoff, "float_value",
    //                          "X-dropoff value (in bits) for final gapped alignment",
    //                          CArgDescriptions::eDouble);
    // ```
    /// `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`: the RPS engine's X-drops in bits
    /// (`RpsSearch`); NCBI's checks of their values are made with the options' others.
    #[arg(long = "xdrop_ungap", value_name = "float_value", value_parser = rpsblast_real, help = "X-dropoff value (in bits) for ungapped extensions")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long = "xdrop_gap", value_name = "float_value", value_parser = rpsblast_real, help = "X-dropoff value (in bits) for preliminary gapped extensions")]
    pub xdrop_gap: Option<f64>,
    #[arg(long = "xdrop_gap_final", value_name = "float_value", value_parser = rpsblast_real, help = "X-dropoff value (in bits) for final gapped alignment")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-244
    // ```c++
    //     if (!m_SuppressSumStats) {
    //         arg_desc.AddOptionalKey(kArgSumStats, "bool_value",
    //                      	 	"Use sum statistics",
    //                      	 	CArgDescriptions::eBoolean);
    //     }
    // ```
    /// `-sum_stats`: `false` is rpsblast's default; `true` is refused where the search starts
    /// (`rps_hooks::rps_search`), after NCBI's errors of the options, the database and the files.
    #[arg(long = "sum_stats", value_name = "bool_value", value_parser = ncbi_boolean, help = "Use sum statistics")]
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-342
    // ```c++
    //         arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options",
    //                         "Filter query sequence with SEG "
    //                         "(Format: '" + kDfltArgApplyFiltering + "', " +
    //                         "'window locut hicut', or '" + kDfltArgNoFiltering +
    //                         "' to disable)",
    //                         CArgDescriptions::eString, m_FilterByDefault
    //                         ? kDfltArgSegFiltering : kDfltArgNoFiltering);
    //         arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                         "Apply filtering locations as soft masks",
    //                         CArgDescriptions::eBoolean,
    //                         kDfltArgLookupTableMaskingOnlyProt);
    // ```
    /// `-seg`, `-soft_masking` and `-lcase_masking` are blastp's query encoding: a hard mask
    /// turns the letters into X for the whole search, the reports show the masks in lower case
    /// and count identities on the letters as read; soft masking keeps the letters (it only
    /// leaves words out of a query lookup table, which an RPS search does not build), so it
    /// changes the outfmt 0 display alone.
    #[arg(
        long,
        value_name = "SEG_options",
        help = "Filter query sequence with SEG (Format: 'yes', 'window locut hicut', or 'no' to disable) [default: no]"
    )]
    pub seg: Option<String>,
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40]")]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
    // ```c++
    //     // lowercase masking
    //     arg_desc.AddFlag(kArgUseLCaseMasking,
    //          "Use lower case filtering in query and subject sequence(s)?", true);
    // ```
    #[arg(
        long = "lcase_masking",
        default_value_t = false,
        help = "Use lower case filtering in query and subject sequence(s)?"
    )]
    pub lcase_masking: bool,
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    #[arg(
        long = "parse_deflines",
        default_value_t = false,
        help = "Should the query and subject defline(s) be parsed?"
    )]
    pub parse_deflines: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:48-66
    // ```c++
    //     const int kDfltRpsThreadingMode = 1;
    //     arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
    //                            "Number of threads to use in RPS BLAST search:\n "
    //                            "0 (auto = num of databases)\n "
    //                            "1 (disable)\n max number of threads = num of databases",
    //                            CArgDescriptions::eInteger,
    //                            NStr::IntToString(kDfltRpsThreadingMode));
    //     arg_desc.SetConstraint(kArgNumThreads,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    //     arg_desc.AddDefaultKey(kArgMTMode, "int_value",
    //                            "Multi-thread mode to use in RPS BLAST search:\n "
    //                            "0 (auto) split by database vols\n "
    //                            "1 split by queries",
    //                             CArgDescriptions::eInteger,
    //                             NStr::IntToString(0));
    //        arg_desc.SetConstraint(kArgMTMode,
    //                               new CArgAllowValuesBetween(0, 1, true));
    // ```
    /// `-num_threads`: LOSAT's pool size (0: the number of volumes, NCBI's "auto"); the output
    /// is the same at every value.
    #[arg(long, default_value_t = 1, value_parser = nonnegative_ncbi_integer)]
    pub num_threads: i32,
    /// `-mt_mode` is read and has no effect: LOSAT searches the queries in parallel and writes
    /// NCBI's `-num_threads 1` output at every value.
    #[arg(long, value_name = "int_value", value_parser = rpsblast_mt_mode, help = "Multi-thread mode: 0 split by database volumes, 1 split by queries (accepted; the output does not depend on it)")]
    pub mt_mode: Option<u8>,
}

/// An `-mt_mode` value of rpsblast: an integer from 0 to 1.
fn rpsblast_mt_mode(value: &str) -> Result<u8, String> {
    let mode = ncbi_integer_at_least(value, 0)?;
    u8::try_from(mode)
        .ok()
        .filter(|mode| *mode <= 1)
        .ok_or_else(|| "Illegal value: expected an integer from 0 to 1".to_string())
}

/// The composition-based statistics mode of a `-comp_based_stats` value as RPS-BLAST reads it:
/// its first letter only, `D`/`d` and any unknown letter as 0, `T`/`t` as 1.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:839-865
/// ```c++
///         ECompoAdjustModes compo_mode = eNoCompositionBasedStats;
///
///         switch (comp_stat_string[0]) {
///             case '0': case 'F': case 'f':
///                 compo_mode = eNoCompositionBasedStats;
///                 break;
///             case '1':
///                 compo_mode = eCompositionBasedStats;
///                 break;
///             case 'D': case 'd':
///                 if ((program == eRPSBlast) || (program == eRPSTblastn)) {
///                     compo_mode = eNoCompositionBasedStats;
///                 }
///                 ...
///             case '2':
///                 compo_mode = eCompositionMatrixAdjust;
///                 break;
///             case '3':
///                 compo_mode = eCompoForceFullMatrixAdjust;
///                 break;
///             case 'T': case 't':
///                 compo_mode = (program == eRPSBlast || program == eRPSTblastn || program == eDeltaBlast) ?
///                     eCompositionBasedStats : eCompositionMatrixAdjust;
///                 break;
///         }
/// ```
pub(crate) fn rps_composition_mode(value: &str) -> u8 {
    match value.as_bytes().first() {
        Some(b'1' | b'T' | b't') => 1,
        Some(b'2') => 2,
        Some(b'3') => 3,
        _ => 0,
    }
}

impl RpsblastArgs {
    /// blastp's arguments of the RPS-BLAST search, after LOSAT's refusal of the output formats of
    /// a later item.
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        Self::check_outfmt(&self.outfmt)?;
        let requested_cbs = rps_composition_mode(&self.comp_based_stats);
        Ok(BlastpArgs {
            query: self.query,
            subject: None,
            db: self.db,
            db_filter: Default::default(),
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            parse_deflines: self.parse_deflines,
            id_scope: None,
            searchsp: None,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:60-66
            // ```c++
            // CBlastRPSOptionsHandle::SetLookupTableDefaults()
            // {
            //     m_Opts->SetLookupTableType(eRPSLookupTable);
            //     // N.B.: the word threshold is not set because for RPS-BLAST it CANNOT be
            //     // done, i.e.: the RPS-BLAST databases already have this value encoded in
            //     // them and therefore this value cannot be changed on an RPS-BLAST run
            // }
            // ```
            // blastp's task gives BLOSUM62 11/1 (the profile database's `.aux` must name these,
            // `rps_search`), window 40, word size 3 and no SEG; the threshold is 0 for the
            // reports (`ResolvedBlastpArgs::threshold`).
            task: "blastp".to_string(),
            evalue: self.evalue,
            threshold: None,
            word_size: None,
            // `-num_threads 0`: one thread per volume of the database (NCBI's "auto"), set once
            // the database is open (`RpsSearch::auto_threads`).
            num_threads: usize::try_from(self.num_threads.max(1)).unwrap_or(1),
            mt_mode: None,
            out: self.out,
            max_target_seqs: self.max_target_seqs,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            // `-max_hsps` and `-qcov_hsp_perc` filter in the RPS traceback (`RpsSearch`).
            max_hsps_per_subject: None,
            qcov_hsp_perc: None,
            culling_limit: None,
            ungapped: false,
            window_size: self.window_size,
            matrix: None,
            gap_open: None,
            gap_extend: None,
            // The requested mode; `rps_hooks::check_freq_ratio_file` sets it to 0 where the
            // database has no frequency ratios.
            comp_based_stats: requested_cbs.to_string(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:68-73
            // ```c++
            // CBlastRPSOptionsHandle::SetQueryOptionDefaults()
            // {
            //     SetSegFiltering(false);
            //     m_Opts->SetStrandOption(objects::eNa_strand_unknown);
            // }
            // ```
            // The query masks are blastp's (`blastp_query_frames`, the SEG of the filtering
            // handler and the lower-case letters of the query reader), applied before the score
            // block (`blast_setup.c:601-624`), so the scan, the extensions, the CBS 0 rescale
            // and the CBS > 0 redo read the masked letters as NCBI's do.
            seg: self.seg,
            soft_masking: self.soft_masking,
            lcase_masking: self.lcase_masking,
            use_sw_tback: false,
            query_loc: self.query_loc,
            subject_loc: None,
            outfmt: self.outfmt,
            archive: None,
            app: BlastpApp::Rpsblast,
            inclusion_ethresh: None,
            first_query_ordinal: 0,
            psiblast_pssm: None,
            formatting: self.formatting,
            rps: Some(RpsSearch {
                max_hsps_per_subject: self.max_hsps_per_subject.unwrap_or(0),
                query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
                requested_cbs,
                auto_threads: self.num_threads == 0,
                x_drop_ungapped_bits: self.xdrop_ungap,
                x_drop_gapped_bits: self.xdrop_gap,
                x_drop_final_bits: self.xdrop_gap_final,
                sum_stats: self.sum_stats == Some(true),
            }),
        })
    }

    /// NCBI parses `-outfmt` before the option handlers (`app::parse_formatting_string`, its
    /// errors are NCBI's); the formats that NCBI runs and LOSAT's RPSBLAST does not write yet are
    /// refused (item PS15a). Outfmt 17 goes on to NCBI's own error (`formatting_handler_check`).
    pub fn check_outfmt(outfmt: &str) -> Result<()> {
        let choice = crate::blastinput::app::parse_formatting_string(outfmt)?;
        if let Some(format) = crate::blastinput::app::ReportFormat::from_number(choice.number) {
            use crate::blastinput::app::ReportFormat;
            if !matches!(
                format,
                ReportFormat::Pairwise
                    | ReportFormat::Tabular
                    | ReportFormat::TabularWithComments
                    | ReportFormat::CommaSeparatedValues
                    | ReportFormat::Sam
            ) {
                bail!(
                    "output format {} is not supported by LOSAT's RPSBLAST",
                    choice.number
                );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};

    fn parse(extra: &[&str]) -> Result<RpsblastArgs, String> {
        let mut argv = vec!["losat", "rpsblast", "-query", "q.faa", "-db", "cdd"];
        argv.extend_from_slice(extra);
        match crate::cli::try_parse_from::<Cli, _, _>(argv) {
            Ok(Cli {
                command: Commands::Rpsblast(args),
                ..
            }) => Ok(args),
            Ok(_) => Err("not rpsblast".to_string()),
            Err(error) => Err(error.to_string()),
        }
    }

    #[test]
    fn the_composition_modes_read_rps_blasts_letters() {
        for (value, mode) in [
            ("0", 0),
            ("F", 0),
            ("f", 0),
            ("D", 0),
            ("d", 0),
            ("1", 1),
            ("T", 1),
            ("t", 1),
            ("2", 2),
            ("3", 3),
            ("x", 0),
            ("", 0),
            ("1u", 1),
        ] {
            assert_eq!(rps_composition_mode(value), mode, "{value:?}");
        }
    }

    #[test]
    fn the_search_is_blastps_with_rpsblasts_program() {
        let args = parse(&[
            "-comp_based_stats",
            "F",
            "-num_threads",
            "0",
            "-mt_mode",
            "1",
        ])
        .unwrap()
        .into_blastp_args()
        .unwrap();
        assert_eq!(args.app, BlastpApp::Rpsblast);
        assert_eq!(args.comp_based_stats, "0");
        assert_eq!(args.num_threads, 1);
        assert_eq!(args.rps.as_ref().unwrap().requested_cbs, 0);
        let args = parse(&["-max_hsps", "2", "-qcov_hsp_perc", "50"])
            .unwrap()
            .into_blastp_args()
            .unwrap();
        // RPS-BLAST's defaults (blastinput_unit_test.cpp:3296-3303 `RpsBlastCBS`): composition-based
        // statistics 1 and no SEG (`blast_rps_options.cpp:68-73`).
        assert_eq!(args.seg, None);
        assert_eq!(
            args.resolve().unwrap().seg,
            crate::algorithm::blastp::args::BlastpSegSpec::No
        );
        let rps = args.rps.unwrap();
        assert_eq!(rps.requested_cbs, 1);
        assert_eq!(
            (rps.max_hsps_per_subject, rps.query_cov_hsp_perc),
            (2, 50.0)
        );
        assert_eq!(args.max_hsps_per_subject, None);
    }

    #[test]
    fn options_of_later_items_are_refused() {
        for extra in [&["-outfmt", "5"][..], &["-outfmt", "1"]] {
            let error = parse(extra).unwrap().into_blastp_args().unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains("not supported by LOSAT's RPSBLAST"),
                "{extra:?}: {error}"
            );
        }
        for name in ["searchsp", "use_sw_tback", "html"] {
            let error = parse(&[&format!("-{name}"), "1"]).unwrap_err();
            assert!(
                error.contains(&format!("-{name} is not supported by LOSAT's RPSBLAST")),
                "{name}: {error}"
            );
        }
        assert!(parse(&["-subject", "s.faa"]).is_err());
        assert!(parse(&["-mt_mode", "2"]).is_err());
    }

    // blast_args.cpp:209-244, 332-342, 1939-1942, 2683-2721 (quoted at the fields): the masks, the
    // X-drops, the formatting options and `-sum_stats false` reach blastp's arguments.
    #[test]
    fn the_options_of_ps11c_1_reach_blastps_arguments() {
        let args = parse(&[
            "-seg",
            "10 1.8 2.1",
            "-soft_masking",
            "true",
            "-lcase_masking",
            "-xdrop_ungap",
            "3",
            "-xdrop_gap",
            "-1",
            "-xdrop_gap_final",
            "0",
            "-line_length",
            "1",
            "-sorthits",
            "4",
            "-sorthsps",
            "0",
            "-sum_stats",
            "F",
        ])
        .unwrap()
        .into_blastp_args()
        .unwrap();
        assert_eq!(args.seg.as_deref(), Some("10 1.8 2.1"));
        assert_eq!(args.soft_masking, Some(true));
        assert!(args.lcase_masking);
        assert_eq!(
            args.formatting,
            FormattingArgs {
                line_length: Some(1),
                sort_hits: Some(4),
                sort_hsps: Some(0),
            }
        );
        let rps = args.rps.as_ref().unwrap();
        assert_eq!(
            (
                rps.x_drop_ungapped_bits,
                rps.x_drop_gapped_bits,
                rps.x_drop_final_bits
            ),
            (Some(3.0), Some(-1.0), Some(0.0))
        );
        assert!(!rps.sum_stats);
        let resolved = args.resolve().unwrap();
        assert_eq!(resolved.line_length, 1);
        assert!(resolved.lcase_masking && resolved.soft_masking);
        // Nothing given: NCBI's defaults.
        let args = parse(&[]).unwrap().into_blastp_args().unwrap();
        assert_eq!(args.formatting, FormattingArgs::default());
        assert_eq!(args.resolve().unwrap().line_length, 60);
        let rps = args.rps.unwrap();
        assert_eq!(
            (
                rps.x_drop_ungapped_bits,
                rps.x_drop_gapped_bits,
                rps.x_drop_final_bits
            ),
            (None, None, None)
        );
    }

    // blast_args.cpp:2683-2721: `-line_length` is 1 or more, `-sorthits` and `-sorthsps` 0 to 4
    // (LOSAT's parser text and rc 2, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 1); the X-drops are
    // NCBI's doubles (`nan` and `inf` are its conversion errors).
    #[test]
    fn the_formatting_and_x_drop_values_are_ncbis() {
        for extra in [
            &["-line_length", "0"][..],
            &["-line_length", "-3"],
            &["-sorthits", "5"],
            &["-sorthsps", "-1"],
            &["-sorthsps", "x"],
            &["-xdrop_gap", "nan"],
            &["-xdrop_ungap", "inf"],
            &["-sum_stats", "maybe"],
        ] {
            assert!(parse(extra).is_err(), "{extra:?}");
        }
        assert_eq!(
            parse(&["-line_length", "2147483647"])
                .unwrap()
                .formatting
                .line_length,
            Some(2_147_483_647)
        );
        // `-sum_stats true` reaches the search, which refuses it (`rps_hooks::rps_search`).
        assert!(
            parse(&["-sum_stats", "true"])
                .unwrap()
                .into_blastp_args()
                .unwrap()
                .rps
                .unwrap()
                .sum_stats
        );
        assert_eq!(
            parse(&["-xdrop_gap", "1e300"]).unwrap().xdrop_gap,
            Some(1e300)
        );
    }

    // blast_options.c:608-618 (quoted at `validate_protein_options`): an ungapped X-drop of 0 or
    // below is NCBI's options error, made when the options are validated (after the SEG
    // option's own error, which its handler raises).
    #[test]
    fn an_ungapped_x_drop_of_zero_or_below_is_ncbis_options_error() {
        for value in ["0", "-1", "-0.0"] {
            let error = parse(&["-xdrop_ungap", value])
                .unwrap()
                .into_blastp_args()
                .unwrap()
                .resolve()
                .unwrap_err();
            assert_eq!(
                format!("{error:#}"),
                "BLAST query/options error: x_dropoff must be greater than zero\nPlease refer to the BLAST+ user manual.\n",
                "{value}"
            );
        }
        let error = parse(&["-xdrop_ungap", "0", "-seg", "maybe"])
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .resolve()
            .unwrap_err();
        assert!(
            format!("{error:#}").contains("Invalid number of arguments to filtering option"),
            "{error:#}"
        );
        // The hit saving options are checked after the X-drop.
        let error = parse(&[
            "-xdrop_ungap",
            "0",
            "-num_descriptions",
            "0",
            "-num_alignments",
            "0",
        ])
        .unwrap()
        .into_blastp_args()
        .unwrap()
        .resolve()
        .unwrap_err();
        assert!(
            format!("{error:#}").contains("x_dropoff must be greater than zero"),
            "{error:#}"
        );
        assert!(parse(&["-xdrop_ungap", "1e-300"])
            .unwrap()
            .into_blastp_args()
            .unwrap()
            .resolve()
            .is_ok());
    }
}
