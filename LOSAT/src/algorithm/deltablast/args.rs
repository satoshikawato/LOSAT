//! Command-line arguments for DELTABLAST (`losat deltablast`)

use anyhow::Result;
use clap::{Args, Parser};
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs, DeltaOptions};
use crate::algorithm::psiblast::PsiblastArgs;
use crate::algorithm::rpsblast::args::RpsblastArgs;
use crate::blastinput::value_parsers::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/deltablast_args.cpp:44-120
// deltablast's arguments are psiblast's without `-ipglist`, `-negative_ipglist`, `-in_pssm`,
// `-in_msa`, `-msa_master_idx`, `-ignore_msa_master` and `-phi_pattern`, with `-rpsdb`,
// `-domain_inclusion_ethresh` and `-show_domain_hits`, and with composition-based statistics 1 by
// default. The fields are psiblast's (`PsiblastArgs`: one definition of each type and parser);
// `into_blastp_args` hands them to psiblast's conversion, whose search deltablast runs from round 1
// on (`BlastpApp::Deltablast`, `blastp::blast_engine::delta_hooks`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct DeltablastArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // `-` (the default) is standard input.
    #[arg(long, default_value = "-", value_parser = ncbi_path(), value_name = "PATH")]
    pub query: PathBuf,
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastinput/app.rs` `missing_subject_error`).
    #[arg(long, value_parser = ncbi_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // NCBI's argument reader reports the exclusion (`cli::reader`).
    #[arg(
        long,
        value_name = "NAME",
        help = "BLAST database name (a protein database, with or without parsed Seq-ids)"
    )]
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist`,
    /// `-negative_seqidlist` and `-no_taxid_expansion`: blastp's (line DI-8), which hold in every
    /// round (NCBI's `CLocalDbAdapter` of the first round searches every later one too).
    #[command(flatten)]
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    #[arg(long, value_parser = ncbi_int8, value_name = "LETTERS", help = "Effective length of the database")]
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    #[arg(
        long = "show_gis",
        default_value_t = false,
        help = "Show NCBI GIs in deflines?"
    )]
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1961-1962
    // (`CQueryOptionsArgs`, psiblast_args.cpp:86)
    #[arg(
        long = "parse_deflines",
        default_value_t = false,
        help = "Should the query and subject defline(s) be parsed?"
    )]
    pub parse_deflines: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
    // psiblast's arguments include `CGenericSearchArgs` (psiblast_args.cpp:66-67). The value is
    // blastp's (`BlastpArgs::searchsp`): every round, those of `-in_pssm` too, then runs without
    // the Gumbel block (`ResolvedBlastpArgs::old_fsc`).
    #[arg(long, value_parser = ncbi_int8, value_name = "INT_VALUE", help = "Effective length of the search space")]
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-245, 321-323
    // psiblast's `CGenericSearchArgs` keeps it (psiblast_args.cpp:67; blastp suppresses it). A
    // gapped protein search does not link HSPs by default (blast_options.c:1460-1471), so `false`
    // is psiblast's search; `true` links the HSPs of each subject in every round
    // (`BlastpArgs::sum_stats`, `core::link_hsps`).
    #[arg(long, value_parser = ncbi_boolean, help = "Use sum statistics")]
    pub sum_stats: Option<bool>,
    #[arg(long, value_parser = ncbi_double, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = ncbi_double, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 11, or the matrix's suggestion]")]
    pub threshold: Option<f64>,
    #[arg(long, value_parser = ncbi_integer, help = "Word size for wordfinder algorithm [default: 3]")]
    pub word_size: Option<i32>,
    #[arg(long, default_value_t = 1, value_parser = ncbi_integer_as::<usize>)]
    pub num_threads: usize,
    #[arg(long, value_name = "PATH", value_parser = ncbi_path())]
    pub out: Option<PathBuf>,
    /// The search strategy file (`-export_search_strategy`, cited at
    /// `BlastpArgs::export_search_strategy`), created by the search's `run` before the other options
    /// are read (`-in_pssm` included: the request holds the PSSM, `PsiSearch::request`).
    #[arg(long, value_name = "filename", value_parser = ncbi_path(), help = "File name to record the search strategy used")]
    pub export_search_strategy: Option<PathBuf>,
    #[arg(long, value_parser = ncbi_integer_as::<usize>, help = "Maximum number of aligned sequences to keep [default: 500]")]
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    #[arg(long = "num_descriptions", value_name = "int_value", value_parser = ncbi_integer_as::<usize>, help = "Number of database sequences to show one-line descriptions for. Not applicable for outfmt > 4 [default: 500]")]
    pub num_descriptions: Option<usize>,
    #[arg(long = "num_alignments", value_name = "int_value", value_parser = ncbi_integer_as::<usize>, help = "Number of database sequences to show alignments for [default: 250]")]
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    #[command(flatten)]
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    #[arg(long = "max_hsps", value_parser = ncbi_integer_as::<usize>)]
    pub max_hsps_per_subject: Option<usize>,
    #[arg(long = "qcov_hsp_perc", value_parser = ncbi_double, help = "Percent query coverage per hsp")]
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3340
    // psiblast's `CHspFilteringArgs` (psiblast_args.cpp:79), blastp's (`BlastpArgs::culling_limit`):
    // every round culls, the PSSM rounds and `-in_pssm` too (`pssm_hooks::pssm_prelim`).
    #[arg(long, value_parser = ncbi_integer, help = "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit")]
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // psiblast's `CHspFilteringArgs`, blastp's (`BlastpArgs::best_hit_overhang`): every round,
    // the PSSM rounds and `-in_pssm` too (`pssm_hooks::pssm_prelim`).
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    #[arg(long, value_parser = ncbi_double, value_name = "float_value")]
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    #[arg(long, value_parser = ncbi_double, value_name = "float_value")]
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329,3347-3349
    // psiblast's `CHspFilteringArgs`, blastp's (`BlastpArgs::subject_besthit`): every round.
    /// Return only the best HSP for each non overlapping query region
    #[arg(long)]
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229, 276-286
    // (cited at `BlastpArgs::xdrop_ungap`)
    // psiblast's `CGenericSearchArgs` (psiblast_args.cpp:67): blastp's X-drops in bits, which
    // every round converts with its own Karlin-Altschul blocks (`pssm_hooks`).
    #[arg(long = "xdrop_ungap", value_name = "float_value", value_parser = ncbi_double, help = "X-dropoff value (in bits) for ungapped extensions [default: 7]")]
    pub xdrop_ungap: Option<f64>,
    #[arg(long = "xdrop_gap", value_name = "float_value", value_parser = ncbi_double, help = "X-dropoff value (in bits) for preliminary gapped extensions [default: 15]")]
    pub xdrop_gap: Option<f64>,
    #[arg(long = "xdrop_gap_final", value_name = "float_value", value_parser = ncbi_double, help = "X-dropoff value (in bits) for final gapped alignment [default: 25]")]
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1058-1079
    // psiblast's `CGapTriggerArgs` (psiblast_args.cpp:105; blastp has none): the cutoff of the
    // ungapped alignments that a gapped extension starts from, in every round
    // (`ResolvedBlastpArgs::gap_trigger_bits`).
    #[arg(long = "gap_trigger", value_name = "float_value", default_value = "22", value_parser = ncbi_double, help = "Number of bits to trigger gapping")]
    pub gap_trigger: f64,
    #[arg(long, value_parser = ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11, or the matrix's best value]")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1, or the matrix's best value]")]
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    #[arg(
        long,
        default_value = "1",
        help = "Use composition-based statistics: D or d (default, as 1), 0 or F or f no composition-based statistics, 1 or T or t composition-based statistics; 2 and 3 are read and run as 1"
    )]
    pub comp_based_stats: String,
    // NCBI's filtering handler (`blastinput/app.rs` `parse_seg_option`); no filtering when
    // omitted (`CPSIBlastOptionsHandle::SetQueryOptionDefaults`).
    #[arg(long, help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: no]")]
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // psiblast's filtering handler is blastp's (`CFilteringArgs`, psiblast_args.cpp:70).
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1940-1942
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    #[arg(long)]
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // NCBI's argument is a string (c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660), parsed before the option
    // handlers (`blastinput/app.rs` `parse_formatting_string`).
    #[arg(long, default_value = "0", value_name = "SPEC")]
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1147-1152
    #[arg(
        long = "num_iterations",
        default_value = "1",
        value_name = "int_value",
        value_parser = ncbi_integer,
        help = "Number of iterations to perform (0 means run until convergence)"
    )]
    pub num_iterations: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1086-1090
    // Only the PSSM engine reads it (`PSIBlastOptions::pseudo_count`), which one round runs only
    // for `-save_pssm_after_last_round` (`pssm.rs`); without it the value changes nothing (NCBI
    // 2.17.0: the same reports with 5 and -1).
    #[arg(
        long,
        default_value = "0",
        value_parser = ncbi_integer,
        help = "Pseudo-count value used when constructing PSSM"
    )]
    pub pseudocount: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1100-1104
    // The first round reads it too: it is the cutoff of the composition-based redo's heaps
    // (`ResolvedBlastpArgs::inclusion_ethresh`).
    #[arg(
        long,
        default_value = "0.002",
        value_parser = ncbi_double,
        help = "E-value inclusion threshold for pairwise alignments"
    )]
    pub inclusion_ethresh: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1156-1160
    // The file is opened when a PSSM is written (`pssm.rs` `AsnPssmWriter`, as
    // `CAutoOutputFileReset`), not by the parser.
    #[arg(
        long = "out_pssm",
        value_name = "checkpoint_file",
        help = "File name to store checkpoint file"
    )]
    pub out_pssm: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1161-1170
    // The file is opened when a PSSM is written (`pssm.rs` `write_ascii_pssm_files`, as
    // `CAutoOutputFileReset`), not by the parser.
    #[arg(
        long = "out_ascii_pssm",
        value_name = "ascii_mtx_file",
        help = "File name to store ASCII version of PSSM"
    )]
    pub out_ascii_pssm: Option<String>,
    #[arg(
        long = "save_pssm_after_last_round",
        help = "Save PSSM after the last database search"
    )]
    pub save_pssm_after_last_round: bool,
    #[arg(
        long = "save_each_pssm",
        help = "Save PSSM after each iteration (file name is given in -save_pssm or -save_ascii_pssm options)"
    )]
    pub save_each_pssm: bool,

    /// `BlastpArgs::first_query_ordinal`: the number of the first query of the records of a call
    /// that starts there (`api::run_batches`). 0 on the command line.
    #[arg(skip)]
    pub first_query_ordinal: usize,
    /// LOSAT's refusal of a record of these records with the id of a query of an earlier call of
    /// the run (`query_batch::ScopeRefusal`), where the library searches each query batch as a call
    /// of its own (`api::run_batches`): the call's add of its batch gives it, so that the batch
    /// stops in the command line's form (line WQA D10 round 5, 2026-11-20). `None` on the command
    /// line. Not a public option.
    #[arg(skip)]
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
    /// The run's numbers of the files of `-save_each_pssm` (`pssm::PssmVersions`, NCBI's
    /// `CAutoOutputFileReset` of the run): a command line's run makes them; the library's query
    /// batches, each a call of its own, share the run's (`api::batches`).
    #[arg(skip)]
    pub pssm_versions: crate::algorithm::psiblast::pssm::PssmVersions,
    /// The queries that the checks of the PSSM files before the search cover
    /// (`pssm::PssmCheckScope`): the call's on the command line; the run's in the first of the
    /// library's query batches, none in the later ones (`api::batches`).
    #[arg(skip)]
    pub pssm_check: crate::algorithm::psiblast::pssm::PssmCheckScope,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/deltablast_args.cpp:44-120
    #[arg(
        long = "domain_inclusion_ethresh",
        default_value = "0.05",
        value_parser = ncbi_double,
        value_name = "ethresh",
        help = "E-value inclusion threshold for alignments with conserved domains"
    )]
    pub domain_inclusion_ethresh: f64,
    #[arg(
        long = "rpsdb",
        default_value = "cdd_delta",
        value_name = "database_name",
        help = "BLAST domain database name"
    )]
    pub rpsdb: String,
    #[arg(long = "show_domain_hits", help = "Show domain hits")]
    pub show_domain_hits: bool,
}

impl DeltablastArgs {
    /// blastp's arguments for deltablast's search, after LOSAT's refusal of the output formats
    /// that it does not write yet.
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        Self::check_outfmt(
            &self.outfmt,
            self.out.as_deref(),
            self.export_search_strategy.is_some(),
        )?;
        let delta = DeltaOptions {
            domain_inclusion_ethresh: self.domain_inclusion_ethresh,
            show_domain_hits: self.show_domain_hits,
            domain: Box::new(Self::domain_search_args(
                &self.rpsdb,
                self.domain_inclusion_ethresh,
                self.lcase_masking,
            )?),
        };
        let num_iterations = self.num_iterations;
        let psi = PsiblastArgs {
            query: self.query,
            subject: self.subject,
            db: self.db,
            db_filter: self.db_filter,
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            parse_deflines: self.parse_deflines,
            searchsp: self.searchsp,
            sum_stats: self.sum_stats,
            evalue: self.evalue,
            threshold: self.threshold,
            word_size: self.word_size,
            num_threads: self.num_threads,
            out: self.out,
            export_search_strategy: self.export_search_strategy,
            max_target_seqs: self.max_target_seqs,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            formatting: self.formatting,
            max_hsps_per_subject: self.max_hsps_per_subject,
            qcov_hsp_perc: self.qcov_hsp_perc,
            culling_limit: self.culling_limit,
            best_hit_overhang: self.best_hit_overhang,
            best_hit_score_edge: self.best_hit_score_edge,
            subject_besthit: self.subject_besthit,
            xdrop_ungap: self.xdrop_ungap,
            xdrop_gap: self.xdrop_gap,
            xdrop_gap_final: self.xdrop_gap_final,
            gap_trigger: self.gap_trigger,
            window_size: self.window_size,
            matrix: self.matrix,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            comp_based_stats: Self::comp_based_stats_as_psiblast(&self.comp_based_stats),
            seg: self.seg,
            soft_masking: self.soft_masking,
            lcase_masking: self.lcase_masking,
            use_sw_tback: self.use_sw_tback,
            query_loc: self.query_loc,
            subject_loc: self.subject_loc,
            outfmt: self.outfmt,
            num_iterations: self.num_iterations,
            pseudocount: self.pseudocount,
            inclusion_ethresh: self.inclusion_ethresh,
            out_pssm: self.out_pssm,
            out_ascii_pssm: self.out_ascii_pssm,
            save_pssm_after_last_round: self.save_pssm_after_last_round,
            save_each_pssm: self.save_each_pssm,
            in_pssm: None,
            in_msa: None,
            msa_master_idx: None,
            ignore_msa_master: false,
            first_query_ordinal: self.first_query_ordinal,
            scope_refusal: self.scope_refusal,
            pssm_versions: self.pssm_versions,
            pssm_check: self.pssm_check,
        };
        let mut args = psi.into_blastp_args()?;
        args.app = BlastpApp::Deltablast;
        args.delta = Some(delta);
        // The options of the search strategy and the archive are psiblast's: LOSAT refuses both
        // for deltablast (`-export_search_strategy`, outfmt 11).
        args.request_inputs = None;
        // NCBI reference (598d8ae6): c++/src/app/blast/deltablast_app.cpp:466-580
        // Every run goes through the round driver: round 1 is the search of the PSSM of the
        // conserved domains, also with `-num_iterations 1`.
        if let Some(pssm) = args.psiblast_pssm.as_mut() {
            pssm.rounds = Some(u32::try_from(num_iterations).unwrap_or(0));
        }
        Ok(args)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:195-215
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:53-78
    /// The domain search's arguments: rpsblast's defaults (NCBI's RPS-BLAST options handle) on the
    /// database `rpsdb`, whose E-value is the domain inclusion threshold, and the query's
    /// lower-case masking (the domain search reads the query that round 1 reads). A threshold of 0
    /// or less keeps the E-value 10 here: NCBI's options error comes after the domain database is
    /// opened (`blastp::blast_engine::delta_hooks`).
    fn domain_search_args(rpsdb: &str, ethresh: f64, lcase_masking: bool) -> Result<BlastpArgs> {
        /// rpsblast's arguments with the database alone.
        #[derive(Parser)]
        struct DomainCommand {
            #[command(flatten)]
            args: RpsblastArgs,
        }
        let command = DomainCommand::try_parse_from(["domain", "--db", rpsdb]).map_err(|error| {
            anyhow::anyhow!("LOSAT internal error: the domain search's options: {error}")
        })?;
        let mut args = command.args.into_blastp_args()?;
        args.evalue = Some(if ethresh > 0.0 { ethresh } else { 10.0 });
        args.lcase_masking = lcase_masking;
        Ok(args)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:837-870
    /// deltablast's reading of `-comp_based_stats`: `D`, `d`, `T` and `t` are mode 1 (psiblast's
    /// are 2), so they are given to psiblast's reader as `1`; every other value is read as
    /// psiblast reads it (2 and 3 keep the warning of PSSM searches).
    fn comp_based_stats_as_psiblast(value: &str) -> String {
        match value.as_bytes().first() {
            Some(b'D' | b'd' | b'T' | b't') => format!("1{}", &value[1..]),
            _ => value.to_string(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2744-2748
    /// NCBI's own reading of `-outfmt` comes first (psiblast's, `PsiblastArgs::check_outfmt`); the
    /// formats that NCBI runs and LOSAT's deltablast does not write yet (all but 0, 6 and 7, item
    /// PS15b-3) are refused, except those that NCBI refuses itself later (`ncbi_refuses_format`).
    pub fn check_outfmt(outfmt: &str, out: Option<&std::path::Path>, strategy: bool) -> Result<()> {
        use crate::blastinput::app::{early_format_choice, ncbi_refuses_format, Program};
        PsiblastArgs::check_outfmt(outfmt, out, strategy)?;
        let Some(choice) = early_format_choice(outfmt, strategy)?.parsed().cloned() else {
            return Ok(());
        };
        if ncbi_refuses_format(&choice, Program::Blastp, out) {
            return Ok(());
        }
        if !matches!(choice.number, 0 | 6 | 7) {
            anyhow::bail!(
                "output format {} is not supported by LOSAT's DELTABLAST yet",
                choice.number
            );
        }
        Ok(())
    }
}
