//! Command-line arguments for RPS-BLAST (`losat rpsblast`).

use anyhow::Result;
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
use crate::algorithm::rps_search::engine::RpsSearch;
use crate::blastinput::report_counts::FormattingArgs;

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
#[derive(Debug, Clone)]
pub struct RpsblastArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    pub query: PathBuf,
    pub db: Option<String>,
    pub out: Option<PathBuf>,
    /// The search strategy file (`-export_search_strategy`, cited at
    /// `BlastpArgs::export_search_strategy`), created by the search's `run` before the other options
    /// are read.
    pub export_search_strategy: Option<PathBuf>,
    pub evalue: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2234-2236
    /// `-dbsize`: replaces the database's length in every statistic (`x_AdjustDbSize` then keeps
    /// the volumes' own sequence counts, `rps_hooks::rps_search`).
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
    /// `-searchsp`: at any value (0 too) the search runs without the Gumbel block
    /// (`ResolvedBlastpArgs::old_fsc`); above 0 it is every query's and volume's effective
    /// search space (`rps_hooks::rps_search`).
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    pub comp_based_stats: String,
    pub outfmt: String,
    pub show_gis: bool,
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    // `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    // blast_args.cpp:2683-2721, cited at `FormattingArgs`): the width of the outfmt 0
    // alignment rows; the sorts are read with NCBI's warnings and sort nothing: NCBI's
    // rpsblast does not hand them to its formatter (`app/rpsblast_app.cpp:187`,
    // `SetLineLength` only).
    pub formatting: FormattingArgs,
    pub max_target_seqs: Option<usize>,
    pub max_hsps_per_subject: Option<usize>,
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229
    /// `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`: the RPS engine's X-drops in bits
    /// (`RpsSearch`); NCBI's checks of their values are made with the options' others.
    pub xdrop_ungap: Option<f64>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-244
    /// `-sum_stats`: `false` is rpsblast's default; `true` is refused where the search starts
    /// (`rps_hooks::rps_search`), after NCBI's errors of the options, the database and the files.
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:879-899 and 892-899
    /// `-use_sw_tback`: the composition-based redo (modes 1 to 3) is NCBI's Smith-Waterman redo
    /// (`RpsSearch::use_sw_tback`); without composition-based statistics it changes nothing.
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-342
    /// `-seg`, `-soft_masking` and `-lcase_masking` are blastp's query encoding: a hard mask
    /// turns the letters into X for the whole search, the reports show the masks in lower case
    /// and count identities on the letters as read; soft masking keeps the letters (it only
    /// leaves words out of a query lookup table, which an RPS search does not build), so it
    /// changes the outfmt 0 display alone.
    pub seg: Option<String>,
    pub soft_masking: Option<bool>,
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
    pub lcase_masking: bool,
    pub query_loc: Option<String>,
    pub parse_deflines: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:48-66
    /// `-num_threads`: LOSAT's pool size (0: the number of volumes, NCBI's "auto"); the output
    /// is the same at every value.
    pub num_threads: i32,
    /// `-mt_mode` is read and has no effect: LOSAT searches the queries in parallel and writes
    /// NCBI's `-num_threads 1` output at every value.
    pub mt_mode: Option<u8>,
}

impl RpsblastArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, rpsblast's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(RpsblastArgs {
            query: args.option_or_default("query")?,
            db: args.given_option("db")?,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
            evalue: args.given_option("evalue")?,
            dbsize: args.given_option("dbsize")?,
            searchsp: args.given_option("searchsp")?,
            comp_based_stats: args.option_or_default("comp_based_stats")?,
            outfmt: args.option_or_default("outfmt")?,
            show_gis: args.flag("show_gis"),
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: FormattingArgs::from_args(args)?,
            max_target_seqs: args.given_option("max_target_seqs")?,
            max_hsps_per_subject: args.given_option("max_hsps")?,
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            xdrop_gap: args.given_option("xdrop_gap")?,
            xdrop_gap_final: args.given_option("xdrop_gap_final")?,
            sum_stats: args.given_option("sum_stats")?,
            use_sw_tback: args.flag("use_sw_tback"),
            seg: args.given_option("seg")?,
            soft_masking: args.given_option("soft_masking")?,
            window_size: args.given_option("window_size")?,
            lcase_masking: args.flag("lcase_masking"),
            query_loc: args.given_option("query_loc")?,
            parse_deflines: args.flag("parse_deflines"),
            num_threads: args.option_or_default("num_threads")?,
            mt_mode: args.given_option("mt_mode")?,
        })
    }
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
            write_verification: None,
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
