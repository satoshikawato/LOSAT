//! Command-line arguments for RPSTBLASTN (`losat rpstblastn`).

use anyhow::Result;
use std::path::PathBuf;
use std::sync::Arc;

use super::search::TranslatedRpsSearch;
use crate::algorithm::blastx::BlastxArgs;
use crate::algorithm::rps_search::engine::RpsSearch;
use crate::algorithm::rpsblast::args::{check_rps_outfmt, rps_composition_mode};
use crate::blastinput::report_counts::FormattingArgs;

/// The text of `rpstblastn -help` (NCBI 2.17.0's, `CArgDescriptions::PrintUsage`), frozen as NCBI
/// writes it: its `Tool`'s help (`cli::reader`), which `-h` shortens as NCBI does.
pub const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:44-110
// rpstblastn's arguments: rpsblast's, with the genetic code of the query (`-query_gencode`), the
// query strands (`-strand`) and `-ungapped`, in the order of NCBI's `-help` (the struct is
// rpsblast's copied, not flattened: a flattened group would change the order of the options).
// The search runs with BLASTX's driver and the RPS engine (`into_blastx_args`), composition-based
// statistics 0 to 3 included (the default 1, item PS12b), with `-use_sw_tback` NCBI's
// Smith-Waterman redo of the frames (item PS12c); `-html`, the remote options and
// `-import_search_strategy` are refused by the parser (their rows of `ncbi_app::args::refusal`).
#[derive(Debug, Clone)]
pub struct RpstblastnArgs {
    // *** Input query options
    pub query: PathBuf,
    pub query_gencode: u8,
    pub query_loc: Option<String>,
    pub strand: String,
    // *** General search options
    pub db: Option<String>,
    pub out: Option<PathBuf>,
    /// The search strategy file (`-export_search_strategy`, cited at
    /// `BlastpArgs::export_search_strategy`), created by the search's `run` before the other options
    /// are read.
    pub export_search_strategy: Option<PathBuf>,
    pub evalue: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    /// Read as RPS-BLAST reads it (`rpsblast::args::rps_composition_mode`).
    pub comp_based_stats: String,
    // *** Formatting options
    pub outfmt: String,
    pub show_gis: bool,
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler, `FormattingArgs`),
    /// handed to BLASTX's arguments: the width of the outfmt 0 alignment rows; the sorts are read
    /// with NCBI's warnings and sort nothing (`ResolvedOptions::report_sort`: NCBI's
    /// rpstblastn hands its formatter the line length only, as rpsblast's).
    pub formatting: FormattingArgs,
    // *** Query filtering options
    pub seg: Option<String>,
    pub soft_masking: Option<bool>,
    pub lcase_masking: bool,
    // *** Restrict search or results
    pub qcov_hsp_perc: Option<f64>,
    pub max_hsps_per_subject: Option<usize>,
    pub max_target_seqs: Option<usize>,
    // *** Statistical options
    pub dbsize: Option<i64>,
    pub searchsp: Option<i64>,
    /// `-sum_stats`: `false` is rpstblastn's default; `true` is refused where the search starts
    /// (`rps_search::translated`), after NCBI's errors of the options, the database and the
    /// queries.
    pub sum_stats: Option<bool>,
    // *** Extension options
    pub xdrop_ungap: Option<f64>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917
    /// `-ungapped`: NCBI's options error with composition-based statistics (BLASTX's
    /// `resolve`); without them NCBI BLAST+ 2.17.0 crashes after the search, which LOSAT refuses
    /// (`rps_search::translated`, DECISIONS 38 D3).
    pub ungapped: bool,
    // *** Miscellaneous options
    pub parse_deflines: bool,
    /// `-num_threads`: LOSAT's pool size (0: the number of volumes, NCBI's "auto"); the output
    /// is the same at every value.
    pub num_threads: i32,
    /// `-mt_mode` is read and has no effect, as rpsblast's.
    pub mt_mode: Option<u8>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    /// `-use_sw_tback`: the composition-based redo (modes 1 to 3) is NCBI's Smith-Waterman redo
    /// of the frames (`RpsSearch::use_sw_tback`); without composition-based statistics it changes
    /// nothing.
    pub use_sw_tback: bool,
}

impl RpstblastnArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, rpstblastn's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(RpstblastnArgs {
            query: args.option_or_default("query")?,
            query_gencode: args.option_or_default("query_gencode")?,
            query_loc: args.given_option("query_loc")?,
            strand: args.option_or_default("strand")?,
            db: args.given_option("db")?,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
            evalue: args.given_option("evalue")?,
            comp_based_stats: args.option_or_default("comp_based_stats")?,
            outfmt: args.option_or_default("outfmt")?,
            show_gis: args.flag("show_gis"),
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: FormattingArgs::from_args(args)?,
            seg: args.given_option("seg")?,
            soft_masking: args.given_option("soft_masking")?,
            lcase_masking: args.flag("lcase_masking"),
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            max_hsps_per_subject: args.given_option("max_hsps")?,
            max_target_seqs: args.given_option("max_target_seqs")?,
            dbsize: args.given_option("dbsize")?,
            searchsp: args.given_option("searchsp")?,
            sum_stats: args.given_option("sum_stats")?,
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            xdrop_gap: args.given_option("xdrop_gap")?,
            xdrop_gap_final: args.given_option("xdrop_gap_final")?,
            window_size: args.given_option("window_size")?,
            ungapped: args.flag("ungapped"),
            parse_deflines: args.flag("parse_deflines"),
            num_threads: args.option_or_default("num_threads")?,
            mt_mode: args.given_option("mt_mode")?,
            use_sw_tback: args.flag("use_sw_tback"),
        })
    }
}

impl RpstblastnArgs {
    /// The RPS search of the arguments (what BLASTX's arguments do not hold), after NCBI's reading
    /// of `-outfmt`, which comes before its option handlers. LOSAT's refusals of options come once
    /// per run before the first byte, after NCBI's own errors of the options and the database
    /// (`TranslatedRpsSearch::check_options`; review A of PS12a, N1).
    pub(crate) fn rps_search(&self) -> Result<TranslatedRpsSearch> {
        // NCBI reads `-outfmt` before the option handlers (`rpsblast::args::check_rps_outfmt`);
        // outfmt 13 and 14 to the standard output get NCBI's own error (BLASTX's).
        check_rps_outfmt(&self.outfmt, self.export_search_strategy.is_some())?;
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
                use_sw_tback: self.use_sw_tback,
            },
            ungapped: self.ungapped,
            several_volumes: false,
            costs: None,
            strategy_db_length: None,
            units: None,
            strategy: false,
            setup: None,
        })
    }

    /// BLASTX's arguments of the RPSTBLASTN search, with the RPS search that takes the place of
    /// BLASTX's (`BlastxArgs::rpstblastn`, `rps_search`).
    pub fn into_blastx_args(self) -> Result<BlastxArgs> {
        let search = self.rps_search()?;
        let requested_cbs = search.search.requested_cbs;
        let request_inputs = Some(Box::new(self.option_inputs()));
        Ok(BlastxArgs {
            query: self.query,
            subject: None,
            subject_loc: None,
            db: self.db,
            db_filter: Default::default(),
            // rpstblastn has no database masks (NCBI's unknown arguments).
            db_mask: Default::default(),
            // `-dbsize`, `-searchsp` and `OLD_FSC` (BLASTX's resolve reads the environment):
            // BLASTX's effective lengths of the frames, the search space over the frames in the
            // hit cutoff (`rps_search::engine`), per volume with `-dbsize` on several volumes.
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            // `-html` is not a psiblast, rpsblast or rpstblastn option (`cli::reader::named_refusal`).
            html: false,
            parse_deflines: self.parse_deflines,
            searchsp: self.searchsp,
            // `-xdrop_ungap` reaches NCBI's options check in BLASTX's resolve ("x_dropoff must be
            // greater than zero"); the search's X-drops are the RPS engine's, at the database's
            // scale (`RpsSearch`), so BLASTX's gapped ones keep their defaults.
            xdrop_ungap: self.xdrop_ungap,
            xdrop_gap: None,
            xdrop_gap_final: None,
            out: self.out,
            export_search_strategy: self.export_search_strategy,
            import_search_strategy: None,
            import: None,
            given: Default::default(),
            task: "blastx".to_string(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:60-67
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
            sum_stats: Some(false),
            gap_open: None,
            gap_extend: None,
            max_intron_length: 0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/rpstblastn_options.cpp:60-64
            // No SEG by default (`kFilterByDefault`); the query masks are BLASTX's, on the
            // translated frames.
            seg: self.seg.unwrap_or_else(|| "no".to_string()),
            soft_masking: self.soft_masking.unwrap_or(false),
            matrix: None,
            threshold: None,
            window_size: self.window_size,
            culling_limit: None,
            // rpstblastn has no best hit options (NCBI's unknown arguments).
            best_hit_overhang: None,
            best_hit_score_edge: None,
            subject_besthit: false,
            ungapped: self.ungapped,
            lcase_masking: self.lcase_masking,
            // RPSTBLASTN's `-use_sw_tback` is its RPS search's (`rps_search`).
            use_sw_tback: false,
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
            request_inputs,
            taxdb_setting: None,
            matrix_search: None,
            first_query_ordinal: 0,
            // The command line's `run` makes the file (`RunSetup::open`); only LOSAT Web's
            // records searches take it here, and they do not run RPSTBLASTN.
            strategy: None,
            scope_refusal: None,
        })
    }

    /// The options that rpstblastn's search strategy requests record
    /// (`options::algorithm_options` of `RequestKind::Rpstblastn`: NCBI's argument classes of
    /// `rpstblastn_args.cpp`), from the arguments as given (`-strand` is never recorded).
    /// `HitlistSize` is the request's (`StrategyExport::for_run`).
    pub(crate) fn option_inputs(&self) -> crate::report::archive::options::OptionInputs {
        use crate::report::archive::options::{seg_filter, OptionInputs, RequestKind};
        let mut inputs = OptionInputs::new(RequestKind::Rpstblastn, "rpstblastn");
        inputs.dbsize = self.dbsize;
        inputs.searchsp = self.searchsp;
        inputs.evalue = self.evalue;
        // `-xdrop_ungap` is never recorded (`SetXDropoff` has no remote option).
        inputs.xdrop_gap = self.xdrop_gap;
        inputs.xdrop_gap_final = self.xdrop_gap_final;
        inputs.qcov_hsp_perc = self.qcov_hsp_perc;
        inputs.max_hsps = self.max_hsps_per_subject.map(|count| count as i32);
        inputs.sum_stats = self.sum_stats;
        inputs.query_gencode = i32::from(self.query_gencode);
        inputs.window_size = self.window_size;
        inputs.soft_masking = self.soft_masking;
        inputs.seg = self.seg.as_deref().and_then(seg_filter);
        inputs.ungapped = self.ungapped;
        inputs.comp_based_stats = i32::from(
            crate::algorithm::rpsblast::args::rps_composition_mode(&self.comp_based_stats),
        );
        inputs.use_sw_tback = self.use_sw_tback;
        inputs
    }
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
