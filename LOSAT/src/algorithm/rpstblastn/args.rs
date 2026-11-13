//! Command-line arguments for RPSTBLASTN (`losat rpstblastn`).

use anyhow::Result;
use clap::Args;
use std::path::PathBuf;
use std::sync::Arc;

use super::search::TranslatedRpsSearch;
use crate::algorithm::blastx::BlastxArgs;
use crate::algorithm::rps_search::engine::RpsSearch;
use crate::algorithm::rpsblast::args::{check_rps_outfmt, rps_composition_mode, rpsblast_mt_mode};
use crate::blastinput::report_counts::FormattingArgs;
use crate::blastinput::value_parsers::*;

/// The text of `rpstblastn -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage`), frozen as NCBI
/// writes it. `losat rpstblastn -help` writes it (`RpstblastnArgs::ncbi_help`).
pub const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3829-3834
// ```c++
//     if (detailed) {
//         x.AddDetails(arr);
//     } else {
//         arr.push_back(kEmptyStr);
//         arr.push_back("Use '-help' to print detailed descriptions of command line arguments");
//     }
// ```
/// The text of `rpstblastn -h`: the head of `HELP` before `OPTIONAL ARGUMENTS`, then NCBI's
/// pointer to `-help` (frozen as NCBI writes it).
pub const USAGE: &str = include_str!("usage.txt");

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:44-110
// ```c++
// CRPSTBlastnAppArgs::CRPSTBlastnAppArgs()
// {
//     CRef<IBlastCmdLineArgs> arg;
//     bool const kFilterByDefault = false;
//     static const string kProgram("rpstblastn");
//     arg.Reset(new CProgramDescriptionArgs(kProgram,
//                                "Translated Reverse Position Specific BLAST"));
//     const bool kQueryIsProtein = false;
//     const bool kIsRpsBlast = true;
//     const bool kIsCBS2and3Supported = false;
//     ...
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs(false, kIsRpsBlast));
//     ...
//     arg.Reset(new CGenericSearchArgs( !kQueryIsProtein, kIsRpsBlast ));
//     ...
//     arg.Reset(new CGeneticCodeArgs(CGeneticCodeArgs::eQuery));
//     ...
//     arg.Reset(new CFilteringArgs( !kQueryIsProtein, kFilterByDefault));
//     ...
//     arg.Reset(new CWindowSizeArg);
//     ...
//     arg.Reset(new CGappedArgs);
//     ...
//     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
//     ...
//     m_FormattingArgs.Reset(new CFormattingArgs);
//     ...
//     m_MTArgs.Reset(new CRPSBlastMTArgs());
//     ...
//     arg.Reset(new CCompositionBasedStatsArgs(kIsCBS2and3Supported, kDfltArgCompBasedStatsRPS ));
// ```
// rpstblastn's arguments: rpsblast's, with the genetic code of the query (`-query_gencode`), the
// query strands (`-strand`) and `-ungapped`, in the order of NCBI's `-help` (the struct is
// rpsblast's copied, not flattened: a flattened group would change the order of the options).
// The search runs with BLASTX's driver and the RPS engine (`into_blastx_args`). Composition-based
// statistics 1 to 3 (the default 1), `-searchsp`, `-dbsize` and `OLD_FSC` are refused until the
// later items of series PR run them (PS12b, PS12a phase B); `-html`, the remote and the search
// strategy options by the parser (`cli.rs` `is_unported_rpstblastn_arg`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct RpstblastnArgs {
    // *** Input query options
    #[arg(long, default_value = "-", value_parser = ncbi_input_path(), value_name = "PATH")]
    pub query: PathBuf,
    #[arg(long, default_value_t = 1, value_parser = genetic_code)]
    pub query_gencode: u8,
    #[arg(long = "query_loc", value_name = "RANGE")]
    pub query_loc: Option<String>,
    #[arg(long, default_value = "both", value_parser = ["both", "plus", "minus"])]
    pub strand: String,
    // *** General search options
    #[arg(long, value_name = "NAME")]
    pub db: Option<String>,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    #[arg(long, value_parser = rpstblastn_real)]
    pub evalue: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // ```c++
    //     arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
    //                            CArgDescriptions::eString, m_DefaultOpt);
    // ```
    /// Read as RPS-BLAST reads it (`rpsblast::args::rps_composition_mode`).
    #[arg(long, default_value = "1")]
    pub comp_based_stats: String,
    // *** Formatting options
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,
    #[arg(long = "show_gis", default_value_t = false)]
    pub show_gis: bool,
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = report_count)]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = report_count)]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler, `FormattingArgs`),
    /// handed to BLASTX's arguments: the width of the outfmt 0 alignment rows; the sorts are read
    /// with NCBI's warnings and sort nothing (`ResolvedOptions::report_sort`: NCBI's
    /// rpstblastn hands its formatter the line length only, as rpsblast's).
    #[command(flatten)]
    pub formatting: FormattingArgs,
    // *** Query filtering options
    #[arg(long, value_name = "SEG_options")]
    pub seg: Option<String>,
    #[arg(long, value_parser = ncbi_boolean)]
    pub soft_masking: Option<bool>,
    #[arg(long = "lcase_masking", default_value_t = false)]
    pub lcase_masking: bool,
    // *** Restrict search or results
    #[arg(long = "qcov_hsp_perc", value_parser = rpstblastn_query_coverage)]
    pub qcov_hsp_perc: Option<f64>,
    #[arg(long = "max_hsps", value_parser = blastn_count)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long, value_parser = blastn_count, conflicts_with_all = ["num_descriptions", "num_alignments"])]
    pub max_target_seqs: Option<usize>,
    // *** Statistical options
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS")]
    pub dbsize: Option<i64>,
    #[arg(long, value_parser = nonnegative_ncbi_int8, value_name = "INT_VALUE")]
    pub searchsp: Option<i64>,
    /// `-sum_stats`: `false` is rpstblastn's default; `true` is refused where the search starts
    /// (`rps_search::translated`), after NCBI's errors of the options, the database and the
    /// queries.
    #[arg(long = "sum_stats", value_name = "bool_value", value_parser = ncbi_boolean)]
    pub sum_stats: Option<bool>,
    // *** Extension options
    #[arg(long = "xdrop_ungap", value_name = "float_value", value_parser = rpstblastn_real)]
    pub xdrop_ungap: Option<f64>,
    #[arg(long = "xdrop_gap", value_name = "float_value", value_parser = rpstblastn_real)]
    pub xdrop_gap: Option<f64>,
    #[arg(long = "xdrop_gap_final", value_name = "float_value", value_parser = rpstblastn_real)]
    pub xdrop_gap_final: Option<f64>,
    #[arg(long, value_parser = nonnegative_ncbi_integer)]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917
    // ```c++
    //     arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    // ```
    /// `-ungapped`: NCBI's options error with composition-based statistics (BLASTX's
    /// `resolve`); without them NCBI BLAST+ 2.17.0 crashes after the search, which LOSAT refuses
    /// (`rps_search::translated`, DECISIONS 38 D3).
    #[arg(long, default_value_t = false)]
    pub ungapped: bool,
    // *** Miscellaneous options
    #[arg(long = "parse_deflines", default_value_t = false)]
    pub parse_deflines: bool,
    /// `-num_threads`: LOSAT's pool size (0: the number of volumes, NCBI's "auto"); the output
    /// is the same at every value.
    #[arg(long, default_value_t = 1, value_parser = nonnegative_ncbi_integer)]
    pub num_threads: i32,
    /// `-mt_mode` is read and has no effect, as rpsblast's.
    #[arg(long, value_name = "int_value", value_parser = rpsblast_mt_mode)]
    pub mt_mode: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    // ```c++
    //     arg_desc.AddFlag(kArgUseSWTraceback,
    //                      "Compute locally optimal Smith-Waterman alignments?",
    //                      true);
    // ```
    /// `-use_sw_tback`: without composition-based statistics NCBI's search does not read it
    /// (accepted); with them it is refused until item PS12c.
    #[arg(long = "use_sw_tback", default_value_t = false)]
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3829-3834 (quoted at `USAGE`)
    /// `-help` or `-help-full` (`full`: NCBI 2.17.0 writes one text for both) or `-h` (`short`):
    /// NCBI's reader stops at them, the words before them read first. The command line's reader sets this hidden argument there
    /// (`cli::try_parse_from_command`; no word of the user names it, `NCBI_HELP_ARG`), and the
    /// run writes NCBI's text (`ncbi_help_text`).
    #[arg(long = NCBI_HELP_ARG, hide = true, value_parser = ["full", "short"])]
    pub(crate) ncbi_help: Option<String>,
}

/// The long name of the reader's own argument `RpstblastnArgs::ncbi_help`, which no word of the
/// command line names (the toolkit's hidden `-logfile` and `-conffile` are named by theirs).
pub(crate) const NCBI_HELP_ARG: &str = "ncbi_help";

impl RpstblastnArgs {
    /// NCBI's text of `-help` or `-h` when the command line asks for it (`ncbi_help`).
    pub fn ncbi_help_text(&self) -> Option<&'static str> {
        match self.ncbi_help.as_deref() {
            Some("full") => Some(HELP),
            Some(_) => Some(USAGE),
            None => None,
        }
    }

    /// The RPS search of the arguments (what BLASTX's arguments do not hold), after LOSAT's
    /// refusal of the output formats that a later item writes (NCBI parses `-outfmt` before its
    /// option handlers). The refusals of the composition-based statistics that a later item
    /// runs come once per run before the first byte, after NCBI's own errors of the options and
    /// the database (`TranslatedRpsSearch::check_options`; review A of PS12a, N1).
    pub(crate) fn rps_search(&self) -> Result<TranslatedRpsSearch> {
        // NCBI parses `-outfmt` before the option handlers (`app::parse_formatting_string`);
        // outfmt 13 and 14 to the standard output get NCBI's own error (BLASTX's).
        let to_standard_output = self
            .out
            .as_deref()
            .is_none_or(|path| path.as_os_str() == "-");
        check_rps_outfmt(&self.outfmt, "RPSTBLASTN", to_standard_output)?;
        let requested_cbs = rps_composition_mode(&self.comp_based_stats);
        Ok(TranslatedRpsSearch {
            search: RpsSearch {
                max_hsps_per_subject: self.max_hsps_per_subject.unwrap_or(0),
                query_cov_hsp_perc: self.qcov_hsp_perc.unwrap_or(0.0),
                requested_cbs,
                auto_threads: self.num_threads == 0,
                x_drop_ungapped_bits: self.xdrop_ungap,
                x_drop_gapped_bits: self.xdrop_gap,
                x_drop_final_bits: self.xdrop_gap_final,
                sum_stats: self.sum_stats == Some(true),
            },
            ungapped: self.ungapped,
            use_sw_tback: self.use_sw_tback,
            comp_based_stats: self.comp_based_stats.clone(),
            several_volumes: false,
        })
    }

    /// BLASTX's arguments of the RPSTBLASTN search, with the RPS search that takes the place of
    /// BLASTX's (`BlastxArgs::rpstblastn`, `rps_search`).
    pub fn into_blastx_args(self) -> Result<BlastxArgs> {
        let search = self.rps_search()?;
        let requested_cbs = search.search.requested_cbs;
        Ok(BlastxArgs {
            query: self.query,
            subject: None,
            subject_loc: None,
            db: self.db,
            db_filter: Default::default(),
            // `-dbsize`, `-searchsp` and `OLD_FSC` (BLASTX's resolve reads the environment):
            // BLASTX's effective lengths of the frames, the search space over the frames in the
            // hit cutoff (`rps_search::engine`), per volume with `-dbsize` on several volumes.
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            parse_deflines: self.parse_deflines,
            searchsp: self.searchsp,
            // `-xdrop_ungap` reaches NCBI's options check in BLASTX's resolve ("x_dropoff must be
            // greater than zero"); the search's X-drops are the RPS engine's, at the database's
            // scale (`RpsSearch`), so BLASTX's gapped ones keep their defaults.
            xdrop_ungap: self.xdrop_ungap,
            xdrop_gap: None,
            xdrop_gap_final: None,
            out: self.out,
            task: "blastx".to_string(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:60-67
            // ```c++
            // CBlastRPSOptionsHandle::SetLookupTableDefaults()
            // {
            //     m_Opts->SetLookupTableType(eRPSLookupTable);
            //     // N.B.: the word threshold is not set because for RPS-BLAST it CANNOT be
            //     // done, i.e.: the RPS-BLAST databases already have this value encoded in
            //     // them and therefore this value cannot be changed on an RPS-BLAST run
            // }
            // ```
            // The word size, the threshold, the matrix and the gap costs are the profile
            // database's (its `.aux`, which must name BLOSUM62 11/1, `rps_search`): BLASTX's
            // defaults are those values for the reports, whose threshold is 0
            // (`ResolvedOptions::threshold`).
            evalue: self.evalue.unwrap_or(10.0),
            word_size: None,
            // `-max_hsps` and `-qcov_hsp_perc` filter in the RPS traceback (`RpsSearch`).
            max_hsps: None,
            qcov_hsp_perc: None,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1460-1465
            // ```c
            //    /* By default, sum statistics is used for all translated searches
            //     * (except RPS BLAST), and for all ungapped searches.
            //     */
            //    if (program_number == eBlastTypeRpsTblastn) {
            // 	   (*options)->do_sum_stats = FALSE;
            //    } else if (!gapped_calculation ||
            // ```
            sum_stats: Some(false),
            gap_open: None,
            gap_extend: None,
            max_intron_length: 0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/rpstblastn_options.cpp:60-64
            // ```c++
            // void
            // CRPSTBlastnOptionsHandle::SetQueryOptionDefaults()
            // {
            //     m_Opts->SetQueryGeneticCode(BLAST_GENETIC_CODE);
            // }
            // ```
            // No SEG by default (`kFilterByDefault`); the query masks are BLASTX's, on the
            // translated frames.
            seg: self.seg.unwrap_or_else(|| "no".to_string()),
            soft_masking: self.soft_masking.unwrap_or(false),
            matrix: None,
            threshold: None,
            window_size: self.window_size,
            culling_limit: None,
            subject_besthit: false,
            ungapped: self.ungapped,
            lcase_masking: self.lcase_masking,
            strand: self.strand,
            query_loc: self.query_loc,
            query_gencode: self.query_gencode,
            outfmt: self.outfmt,
            max_target_seqs: self
                .max_target_seqs
                .map(|n| i32::try_from(n).unwrap_or(i32::MAX)),
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            formatting: self.formatting,
            // `-num_threads 0`: one thread per volume of the database (NCBI's "auto"), set once
            // the database is open (`RpsSearch::auto_threads`).
            num_threads: self.num_threads.max(1),
            mt_mode: None,
            // The mode as BLASTX's handler reads a digit, so that NCBI's error of an ungapped
            // search with composition-based statistics is BLASTX's (`resolve`).
            comp_based_stats: requested_cbs.to_string(),
            rpstblastn: Some(Arc::new(search)),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};

    fn parse(extra: &[&str]) -> Result<RpstblastnArgs, String> {
        let mut argv = vec!["losat", "rpstblastn", "-query", "q.fna", "-db", "cdd"];
        argv.extend_from_slice(extra);
        match crate::cli::try_parse_from::<Cli, _, _>(argv) {
            Ok(Cli {
                command: Commands::Rpstblastn(args),
                ..
            }) => Ok(args),
            Ok(_) => Err("not rpstblastn".to_string()),
            Err(error) => Err(crate::cli::render_message(&error)),
        }
    }

    fn refusal(extra: &[&str]) -> String {
        parse(extra)
            .unwrap()
            .into_blastx_args()
            .unwrap_err()
            .to_string()
    }

    /// The refusal of `check_options`, once per run.
    fn late_refusal(extra: &[&str]) -> String {
        search(extra).check_options().unwrap_err().to_string()
    }

    fn search(extra: &[&str]) -> TranslatedRpsSearch {
        parse(extra).unwrap().rps_search().unwrap()
    }

    // ncbiargs.cpp:3829-3834 (quoted at `USAGE`): `-help` and `-h` write NCBI's texts, whatever
    // follows them; the words before them are read first.
    #[test]
    fn help_and_usage_are_ncbis_texts() {
        assert_eq!(HELP.len(), 9266);
        assert!(HELP.ends_with("   Compute locally optimal Smith-Waterman alignments?\n\n"));
        let head = &HELP[..HELP.find("OPTIONAL ARGUMENTS").unwrap()];
        assert_eq!(
            USAGE,
            format!("{head}Use '-help' to print detailed descriptions of command line arguments\n")
        );
        assert!(USAGE.contains(&format!(
            "   Translated Reverse Position Specific BLAST {}\n",
            crate::NCBI_BLAST_VERSION
        )));
        for (words, text) in [
            (&["-help"][..], HELP),
            (&["-help", "-no_such_option"], HELP),
            (&["--help"], HELP),
            (&["-h"], USAGE),
            (&["-evalue", "1", "-h", "-xyz", "2"], USAGE),
            // NCBI checks the dependencies of the words after its reader (`x_PostCheck`).
            (
                &["-max_target_seqs", "1", "-num_alignments", "1", "-help"],
                HELP,
            ),
        ] {
            assert_eq!(
                parse(words).unwrap().ncbi_help_text(),
                Some(text),
                "{words:?}"
            );
        }
        assert_eq!(parse(&[]).unwrap().ncbi_help_text(), None);
        // The hidden argument is the reader's: no word of the user names it.
        for words in [&["-ncbi_help", "full"][..], &["--ncbi_help=full"]] {
            assert!(
                parse(words).unwrap_err().contains("unknown option"),
                "{words:?}"
            );
        }
        // A word before the help that NCBI's reader refuses is refused first.
        assert!(parse(&["-no_such_option", "-help"])
            .unwrap_err()
            .contains("unknown option"));
        assert!(parse(&["-evalue", "x", "-help"])
            .unwrap_err()
            .contains("-evalue"));
    }

    // blast_args.cpp:839-865 (quoted at `rps_composition_mode`): RPS-BLAST's letters; mode 0 runs,
    // the others are refused until item PS12b; with -ungapped BLASTX's resolve gives NCBI's error.
    #[test]
    fn composition_modes_run_without_statistics_only() {
        for value in ["0", "F", "f", "D", "d"] {
            let args = parse(&["-comp_based_stats", value])
                .unwrap()
                .into_blastx_args()
                .unwrap();
            assert_eq!(args.comp_based_stats, "0", "{value}");
            assert!(args.rpstblastn.is_some());
            assert_eq!(
                search(&["-comp_based_stats", value]).search.requested_cbs,
                0
            );
        }
        for value in ["1", "T", "t", "2", "3"] {
            assert!(
                late_refusal(&["-comp_based_stats", value]).contains(&format!(
                    "-comp_based_stats {value} is not supported by LOSAT's RPSTBLASTN yet"
                )),
                "{value}"
            );
        }
        // The default (1).
        assert_eq!(
            late_refusal(&[]),
            "-comp_based_stats 1 is not supported by LOSAT's RPSTBLASTN yet"
        );
        // ungapped with the default: NCBI's options error (BLASTX's resolve).
        let error = parse(&["-ungapped"])
            .unwrap()
            .into_blastx_args()
            .unwrap()
            .resolve()
            .unwrap_err();
        assert!(
            error.to_string().contains(
                "Composition-adjusted searched are not supported with an ungapped search"
            ),
            "{error}"
        );
        // ungapped without statistics reaches the search, which refuses it.
        assert!(search(&["-ungapped", "-comp_based_stats", "0"]).ungapped);
    }

    // The refusals of the options that later items run, before the search.
    #[test]
    fn options_of_later_items_are_refused() {
        for (extra, text) in [
            (
                &["-outfmt", "5"][..],
                "output format 5 is not supported by LOSAT's RPSTBLASTN",
            ),
            (
                &["-outfmt", "1"],
                "output format 1 is not supported by LOSAT's RPSTBLASTN",
            ),
            (
                &["-outfmt", "11"],
                "output format 11 is not supported by LOSAT's RPSTBLASTN",
            ),
            (
                &["-outfmt", "20"],
                "output format 20 is not supported by LOSAT's RPSTBLASTN",
            ),
            (
                &["-outfmt", "13", "-out", "x.json"],
                "output format 13 is not supported by LOSAT's RPSTBLASTN",
            ),
        ] {
            assert!(
                refusal(extra).contains(text),
                "{extra:?}: {}",
                refusal(extra)
            );
        }
        // The composition-based statistics of later items, once per run (`check_options`).
        assert!(late_refusal(&["-use_sw_tback"])
            .contains("-use_sw_tback with composition-based statistics is not supported by LOSAT's RPSTBLASTN yet"));
        // -searchsp and -dbsize reach BLASTX's options.
        let args = parse(&[
            "-comp_based_stats",
            "0",
            "-searchsp",
            "1000",
            "-dbsize",
            "5000",
        ])
        .unwrap()
        .into_blastx_args()
        .unwrap();
        assert_eq!((args.searchsp, args.dbsize), (Some(1000), Some(5000)));
        // The formats LOSAT writes, a custom field list and NCBI's own error of 17 (BLASTX's).
        // 13 and 14 to the standard output go on to NCBI's error too (BLASTX's).
        for outfmt in [
            "0",
            "6",
            "7",
            "10",
            "6 qseqid sseqid qframe",
            "17",
            "13",
            "14",
        ] {
            assert!(
                parse(&["-comp_based_stats", "0", "-outfmt", outfmt])
                    .unwrap()
                    .into_blastx_args()
                    .is_ok(),
                "{outfmt}"
            );
        }
        // -use_sw_tback without statistics: accepted, NCBI's search does not read it.
        assert!(parse(&["-use_sw_tback", "-comp_based_stats", "0"])
            .unwrap()
            .into_blastx_args()
            .is_ok());
        for name in ["html", "remote", "entrez_query", "import_search_strategy"] {
            let error = parse(&[&format!("-{name}"), "x"]).unwrap_err();
            assert!(
                error.contains(&format!("-{name} is not supported by LOSAT's RPSTBLASTN")),
                "{name}: {error}"
            );
        }
        // NCBI's options of the other programs: LOSAT's parser (PD exception 1).
        for name in [
            "subject",
            "word_size",
            "matrix",
            "threshold",
            "gapopen",
            "task",
        ] {
            assert!(
                parse(&[&format!("-{name}"), "1"])
                    .unwrap_err()
                    .contains("unknown option"),
                "{name}"
            );
        }
    }

    // NCBI's values: the genetic codes, the strands, the threads and modes, the formatting values.
    #[test]
    fn values_are_ncbis() {
        for extra in [
            &["-query_gencode", "7"][..],
            &["-query_gencode", "0"],
            &["-strand", "foo"],
            &["-num_threads", "-1"],
            &["-mt_mode", "2"],
            &["-line_length", "0"],
            &["-sorthits", "5"],
            &["-max_hsps", "0"],
            &["-max_target_seqs", "0"],
            &["-qcov_hsp_perc", "101"],
            &["-sum_stats", "maybe"],
            &["-xdrop_gap", "nan"],
            &["-window_size", "-1"],
            &["-searchsp", "-1"],
            &["-max_target_seqs", "5", "-num_alignments", "2"],
        ] {
            assert!(parse(extra).is_err(), "{extra:?}");
        }
        for code in [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 25, 31, 33,
        ] {
            assert_eq!(
                parse(&["-query_gencode", &code.to_string()])
                    .unwrap()
                    .query_gencode,
                code
            );
        }
    }

    // What reaches BLASTX's arguments and the RPS search: the query handling, the masks (no SEG
    // by default), the counts, the X-drops, `-max_hsps` and `-qcov_hsp_perc` (the RPS
    // traceback's), the formatting options, `-num_threads 0`, no sum statistics.
    #[test]
    fn the_search_is_blastxs_with_the_rps_search() {
        let words = [
            "-comp_based_stats",
            "F",
            "-strand",
            "minus",
            "-query_gencode",
            "4",
            "-query_loc",
            "10-500",
            "-lcase_masking",
            "-max_hsps",
            "2",
            "-qcov_hsp_perc",
            "50",
            "-xdrop_ungap",
            "3",
            "-xdrop_gap_final",
            "0",
            "-line_length",
            "70",
            "-sorthits",
            "1",
            "-num_threads",
            "0",
            "-mt_mode",
            "1",
            "-max_target_seqs",
            "7",
            "-evalue",
            "0.5",
            "-window_size",
            "0",
        ];
        let args = parse(&words).unwrap().into_blastx_args().unwrap();
        assert_eq!(
            (
                args.strand.as_str(),
                args.query_gencode,
                args.query_loc.as_deref()
            ),
            ("minus", 4, Some("10-500"))
        );
        assert!(args.lcase_masking);
        assert_eq!((args.seg.as_str(), args.soft_masking), ("no", false));
        assert_eq!((args.max_hsps, args.qcov_hsp_perc), (None, None));
        assert_eq!((args.max_target_seqs, args.evalue), (Some(7), 0.5));
        assert_eq!((args.num_threads, args.window_size), (1, Some(0)));
        assert_eq!(args.sum_stats, Some(false));
        let rps = search(&words);
        assert_eq!(
            (
                rps.search.max_hsps_per_subject,
                rps.search.query_cov_hsp_perc,
                rps.search.auto_threads
            ),
            (2, 50.0, true)
        );
        assert_eq!(
            (
                rps.search.x_drop_ungapped_bits,
                rps.search.x_drop_gapped_bits,
                rps.search.x_drop_final_bits
            ),
            (Some(3.0), None, Some(0.0))
        );
        // The formatting options are BLASTX's (one owner): the line length reaches the report, the
        // sorts sort nothing (NCBI's rpstblastn hands its formatter the line length only).
        assert_eq!(
            args.formatting,
            FormattingArgs {
                line_length: Some(70),
                sort_hits: Some(1),
                sort_hsps: None,
            }
        );
        let resolved = args.resolve().unwrap();
        assert_eq!(resolved.formatting.line_length(), 70);
        assert_eq!(resolved.report_sort(), Default::default());
        // RPS-BLAST prints no word threshold (`blast_rps_options.cpp:60-67`, quoted above).
        assert_eq!(resolved.threshold, 0.0);
        assert!(!resolved.sum_stats);
        assert_eq!(resolved.composition, 0);
        assert!(resolved.rpstblastn.is_some());
        // `-sum_stats true` reaches the search, which refuses it.
        let words = [
            "-comp_based_stats",
            "0",
            "-sum_stats",
            "true",
            "-seg",
            "yes",
        ];
        assert!(search(&words).search.sum_stats);
        assert_eq!(
            parse(&words).unwrap().into_blastx_args().unwrap().seg,
            "yes"
        );
    }
}
