//! Command-line arguments for RPS-BLAST (`losat rpsblast`).

use anyhow::Result;
use clap::Args;
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
use crate::algorithm::rps_search::engine::RpsSearch;
use crate::blastinput::report_counts::FormattingArgs;
use crate::blastinput::value_parsers::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:70-131
// rpsblast's arguments: blastp's search, filtering, query, formatting and composition-based
// statistics arguments without the matrix, the gap costs, the word size and threshold (the
// profile database holds them), with RPS-BLAST's thread arguments. LOSAT runs the search with
// blastp's driver (`BlastpApp::Rpsblast`) and the RPS engine (`algorithm::rps_search`).
// Composition-based statistics (modes 1 to 3, the default 1) run the redo of `rps_search::redo`
// with the database's frequency ratios, or without them where the database has none (NCBI then
// searches with mode 0), and with `-use_sw_tback` NCBI's Smith-Waterman redo (PS11c-3); the
// argument that LOSAT does not run yet (`-html`, PS15a) is refused by the parser (its row of
// `ncbi_app::args::refusal`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct RpsblastArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    #[arg(long, default_value = "-", value_parser = ncbi_path(), value_name = "PATH")]
    pub query: PathBuf,
    #[arg(
        long,
        value_name = "NAME",
        help = "BLAST database name (a profile database)"
    )]
    pub db: Option<String>,
    #[arg(long, value_name = "PATH", value_parser = ncbi_path())]
    pub out: Option<PathBuf>,
    /// The search strategy file (`-export_search_strategy`, cited at
    /// `BlastpArgs::export_search_strategy`), created by the search's `run` before the other options
    /// are read.
    #[arg(long, value_name = "filename", value_parser = ncbi_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    #[arg(long, value_parser = ncbi_double, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2234-2236
    /// `-dbsize`: replaces the database's length in every statistic (`x_AdjustDbSize` then keeps
    /// the volumes' own sequence counts, `rps_hooks::rps_search`).
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
    /// `-searchsp`: at any value (0 too) the search runs without the Gumbel block
    /// (`ResolvedBlastpArgs::old_fsc`); above 0 it is every query's and volume's effective
    /// search space (`rps_hooks::rps_search`).
    #[arg(long, value_parser = ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
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
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = ncbi_integer_as::<usize>, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = ncbi_integer_as::<usize>, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    // `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    // blast_args.cpp:2683-2721, cited at `FormattingArgs`): the width of the outfmt 0
    // alignment rows; the sorts are read with NCBI's warnings and sort nothing: NCBI's
    // rpsblast does not hand them to its formatter (`app/rpsblast_app.cpp:187`,
    // `SetLineLength` only).
    #[command(flatten)]
    pub formatting: FormattingArgs,
    #[arg(long, value_parser = ncbi_integer_as::<usize>, help = "Maximum number of aligned sequences to keep [default: 500]")]
    pub max_target_seqs: Option<usize>,
    #[arg(long = "max_hsps", value_parser = ncbi_integer_as::<usize>)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = ncbi_double, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229
    /// `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`: the RPS engine's X-drops in bits
    /// (`RpsSearch`); NCBI's checks of their values are made with the options' others.
    #[arg(long = "xdrop_ungap", value_name = "float_value", value_parser = ncbi_double, help = "X-dropoff value (in bits) for ungapped extensions")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long = "xdrop_gap", value_name = "float_value", value_parser = ncbi_double, help = "X-dropoff value (in bits) for preliminary gapped extensions")]
    pub xdrop_gap: Option<f64>,
    #[arg(long = "xdrop_gap_final", value_name = "float_value", value_parser = ncbi_double, help = "X-dropoff value (in bits) for final gapped alignment")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-244
    /// `-sum_stats`: `false` is rpsblast's default; `true` is refused where the search starts
    /// (`rps_hooks::rps_search`), after NCBI's errors of the options, the database and the files.
    #[arg(long = "sum_stats", value_name = "bool_value", value_parser = ncbi_boolean, help = "Use sum statistics")]
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:879-899 and 892-899
    /// `-use_sw_tback`: the composition-based redo (modes 1 to 3) is NCBI's Smith-Waterman redo
    /// (`RpsSearch::use_sw_tback`); without composition-based statistics it changes nothing.
    #[arg(
        long = "use_sw_tback",
        default_value_t = false,
        help = "Compute locally optimal Smith-Waterman alignments?"
    )]
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-342
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
    #[arg(long, value_parser = ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40]")]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
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
    /// `-num_threads`: LOSAT's pool size (0: the number of volumes, NCBI's "auto"); the output
    /// is the same at every value.
    #[arg(long, default_value_t = 1, value_parser = ncbi_integer)]
    pub num_threads: i32,
    /// `-mt_mode` is read and has no effect: LOSAT searches the queries in parallel and writes
    /// NCBI's `-num_threads 1` output at every value.
    #[arg(long, value_name = "int_value", value_parser = ncbi_integer_as::<u8>, help = "Multi-thread mode: 0 split by database volumes, 1 split by queries (accepted; the output does not depend on it)")]
    pub mt_mode: Option<u8>,
}

/// The composition-based statistics mode of a `-comp_based_stats` value as RPS-BLAST reads it:
/// its first letter only, `D`/`d` and any unknown letter as 0, `T`/`t` as 1.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:836-865
pub(crate) fn rps_composition_mode(value: &str) -> u8 {
    match value.as_bytes().first() {
        Some(b'1' | b'T' | b't') => 1,
        Some(b'2') => 2,
        Some(b'3') => 3,
        _ => 0,
    }
}

impl RpsblastArgs {
    /// blastp's arguments of the RPS-BLAST search, after NCBI's reading of `-outfmt`.
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        Self::check_outfmt(&self.outfmt, self.export_search_strategy.is_some())?;
        let requested_cbs = rps_composition_mode(&self.comp_based_stats);
        let request_inputs = Some(Box::new(self.option_inputs()));
        Ok(BlastpArgs {
            query: self.query,
            subject: None,
            db: self.db,
            db_filter: Default::default(),
            // rpsblast has no database masks (NCBI's unknown arguments).
            db_mask: Default::default(),
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            // `-html` is not a psiblast, rpsblast or rpstblastn option (`cli::reader::named_refusal`).
            html: false,
            parse_deflines: self.parse_deflines,
            id_scope: None,
            searchsp: self.searchsp,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:60-66
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
            // rpsblast has no `-best_hit_*` (NCBI's parser neither).
            best_hit_overhang: None,
            best_hit_score_edge: None,
            subject_besthit: false,
            ungapped: false,
            // rpsblast's X-drops are its RPS engine's (`rps`).
            xdrop_ungap: None,
            xdrop_gap: None,
            xdrop_gap_final: None,
            window_size: self.window_size,
            matrix: None,
            gap_open: None,
            gap_extend: None,
            // The requested mode; `rps_hooks::check_freq_ratio_file` sets it to 0 where the
            // database has no frequency ratios.
            comp_based_stats: requested_cbs.to_string(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:68-73
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
            export_search_strategy: self.export_search_strategy,
            import_search_strategy: None,
            import: None,
            given: Default::default(),
            strategy: None,
            request_inputs,
            app: BlastpApp::Rpsblast,
            taxdb_setting: None,
            matrix_search: None,
            inclusion_ethresh: None,
            gap_trigger: None,
            sum_stats: false,
            first_query_ordinal: 0,
            scope_refusal: None,
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
                use_sw_tback: self.use_sw_tback,
            }),
        })
    }

    /// The options that rpsblast's search strategy requests record (`options::algorithm_options`
    /// of `RequestKind::Rpsblast`: NCBI's argument classes of `rpsblast_args.cpp`), from the
    /// arguments as given. `HitlistSize` is the request's (`StrategyExport::for_run`).
    pub(crate) fn option_inputs(&self) -> crate::report::archive::options::OptionInputs {
        use crate::report::archive::options::{seg_filter, OptionInputs, RequestKind};
        let mut inputs = OptionInputs::new(RequestKind::Rpsblast, "blastp");
        inputs.dbsize = self.dbsize;
        inputs.searchsp = self.searchsp;
        inputs.evalue = self.evalue;
        // `-xdrop_ungap` is never recorded (`SetXDropoff` has no remote option).
        inputs.xdrop_gap = self.xdrop_gap;
        inputs.xdrop_gap_final = self.xdrop_gap_final;
        inputs.qcov_hsp_perc = self.qcov_hsp_perc;
        inputs.max_hsps = self.max_hsps_per_subject.map(|count| count as i32);
        inputs.sum_stats = self.sum_stats;
        inputs.window_size = self.window_size;
        inputs.soft_masking = self.soft_masking;
        inputs.seg = self.seg.as_deref().and_then(seg_filter);
        inputs.comp_based_stats = i32::from(rps_composition_mode(&self.comp_based_stats));
        inputs.use_sw_tback = self.use_sw_tback;
        inputs
    }

    /// NCBI reads `-outfmt` before the option handlers, with its own errors (`check_rps_outfmt`;
    /// `strategy`: the run writes a search strategy file).
    pub fn check_outfmt(outfmt: &str, strategy: bool) -> Result<()> {
        check_rps_outfmt(outfmt, strategy)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
/// NCBI's early reading of `-outfmt` before the option handlers of its RPS programs (rpsblast,
/// rpstblastn), whose errors are NCBI's (`app::early_format_choice`; `strategy`: the run writes a
/// search strategy file, before which NCBI reads no `-outfmt`, so a value that does not parse
/// goes on to the formatting handler's error, line RG F/d36). Every format that they run is
/// written (items PS15a-1 and PS15a-2): RPSBLAST's through BLASTP's reports, RPSTBLASTN's through
/// BLASTX's. A value that NCBI refuses later (outfmt 17, a custom delimiter that the field list
/// uses, outfmt 13 and 14 without a file name in `-out`) goes on to NCBI's error at its place,
/// after the arguments' own errors and before the database is opened (line RG F/d36 probes
/// rps.conf20, rps.conf20.q, rps.conf20.nodb, rpst.conf20, rps.o14, rpst.o14).
pub(crate) fn check_rps_outfmt(outfmt: &str, strategy: bool) -> Result<()> {
    crate::blastinput::app::early_format_choice(outfmt, strategy)?;
    Ok(())
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
