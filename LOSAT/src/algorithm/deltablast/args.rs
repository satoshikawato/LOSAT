//! Command-line arguments for DELTABLAST (`losat deltablast`)

use anyhow::Result;
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs, DeltaOptions};
use crate::algorithm::psiblast::PsiblastArgs;
use crate::algorithm::rpsblast::args::RpsblastArgs;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/deltablast_args.cpp:44-120
// deltablast's arguments are psiblast's without `-ipglist`, `-negative_ipglist`, `-in_pssm`,
// `-in_msa`, `-msa_master_idx`, `-ignore_msa_master` and `-phi_pattern`, with `-rpsdb`,
// `-domain_inclusion_ethresh` and `-show_domain_hits`, and with composition-based statistics 1 by
// default. The fields are psiblast's (`PsiblastArgs`: one definition of each type);
// `into_blastp_args` hands them to psiblast's conversion, whose search deltablast runs from round 1
// on (`BlastpApp::Deltablast`, `blastp::blast_engine::delta_hooks`).
#[derive(Debug, Clone)]
pub struct DeltablastArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // `-` (the default) is standard input.
    pub query: PathBuf,
    // Optional for the parser: a missing `-subject` is NCBI's error after the parsing
    // (`blastinput/app.rs` `missing_subject_error`).
    pub subject: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2371
    // NCBI's argument reader reports the exclusion (`cli::reader`).
    pub db: Option<String>,
    /// `-taxids`, `-taxidlist`, `-negative_taxids`, `-negative_taxidlist`, `-seqidlist`,
    /// `-negative_seqidlist` and `-no_taxid_expansion`: blastp's (line DI-8), which hold in every
    /// round (NCBI's `CLocalDbAdapter` of the first round searches every later one too).
    pub db_filter: crate::blastinput::db_filter::DbFilterArgs,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2231-2237
    pub dbsize: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2663
    pub show_gis: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1961-1962
    // (`CQueryOptionsArgs`, psiblast_args.cpp:86)
    pub parse_deflines: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:234-238
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:303-307
    // psiblast's arguments include `CGenericSearchArgs` (psiblast_args.cpp:66-67). The value is
    // blastp's (`BlastpArgs::searchsp`): every round, those of `-in_pssm` too, then runs without
    // the Gumbel block (`ResolvedBlastpArgs::old_fsc`).
    pub searchsp: Option<i64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:240-245, 321-323
    // psiblast's `CGenericSearchArgs` keeps it (psiblast_args.cpp:67; blastp suppresses it). A
    // gapped protein search does not link HSPs by default (blast_options.c:1460-1471), so `false`
    // is psiblast's search; `true` links the HSPs of each subject in every round
    // (`BlastpArgs::sum_stats`, `core::link_hsps`).
    pub sum_stats: Option<bool>,
    pub evalue: Option<f64>,
    pub threshold: Option<f64>,
    pub word_size: Option<i32>,
    pub num_threads: usize,
    pub out: Option<PathBuf>,
    /// The search strategy file (`-export_search_strategy`, cited at
    /// `BlastpArgs::export_search_strategy`), created by the search's `run` before the other options
    /// are read (`-in_pssm` included: the request holds the PSSM, `PsiSearch::request`).
    pub export_search_strategy: Option<PathBuf>,
    pub max_target_seqs: Option<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2665-2685,2732-2738: 0 or more; -max_target_seqs excludes both.
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
    /// `-line_length`, `-sorthits` and `-sorthsps` (NCBI's formatting handler,
    /// `blastinput::report_counts::FormattingArgs`).
    pub formatting: crate::blastinput::report_counts::FormattingArgs,
    pub max_hsps_per_subject: Option<usize>,
    pub qcov_hsp_perc: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3340
    // psiblast's `CHspFilteringArgs` (psiblast_args.cpp:79), blastp's (`BlastpArgs::culling_limit`):
    // every round culls, the PSSM rounds and `-in_pssm` too (`pssm_hooks::pssm_prelim`).
    pub culling_limit: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3304-3328
    // psiblast's `CHspFilteringArgs`, blastp's (`BlastpArgs::best_hit_overhang`): every round,
    // the PSSM rounds and `-in_pssm` too (`pssm_hooks::pssm_prelim`).
    /// Best Hit algorithm overhang value (recommended value: 0.1)
    pub best_hit_overhang: Option<f64>,
    /// Best Hit algorithm score edge value (recommended value: 0.1)
    pub best_hit_score_edge: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329,3347-3349
    // psiblast's `CHspFilteringArgs`, blastp's (`BlastpArgs::subject_besthit`): every round.
    /// Return only the best HSP for each non overlapping query region
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:209-229, 276-286
    // (cited at `BlastpArgs::xdrop_ungap`)
    // psiblast's `CGenericSearchArgs` (psiblast_args.cpp:67): blastp's X-drops in bits, which
    // every round converts with its own Karlin-Altschul blocks (`pssm_hooks`).
    pub xdrop_ungap: Option<f64>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1058-1079
    // psiblast's `CGapTriggerArgs` (psiblast_args.cpp:105; blastp has none): the cutoff of the
    // ungapped alignments that a gapped extension starts from, in every round
    // (`ResolvedBlastpArgs::gap_trigger_bits`).
    pub gap_trigger: f64,
    pub window_size: Option<i32>,
    pub matrix: Option<String>,
    pub gap_open: Option<i32>,
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    pub comp_based_stats: String,
    // NCBI's filtering handler (`blastinput/app.rs` `parse_seg_option`); no filtering when
    // omitted (`CPSIBlastOptionsHandle::SetQueryOptionDefaults`).
    pub seg: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // psiblast's filtering handler is blastp's (`CFilteringArgs`, psiblast_args.cpp:70).
    pub soft_masking: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1940-1942
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:804-806
    pub use_sw_tback: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    pub query_loc: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    pub subject_loc: Option<String>,
    // NCBI's argument is a string (c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660), parsed before the option
    // handlers (`blastinput/app.rs` `parse_formatting_string`).
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1147-1152
    pub num_iterations: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1086-1090
    // Only the PSSM engine reads it (`PSIBlastOptions::pseudo_count`), which one round runs only
    // for `-save_pssm_after_last_round` (`pssm.rs`); without it the value changes nothing (NCBI
    // 2.17.0: the same reports with 5 and -1).
    pub pseudocount: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1100-1104
    // The first round reads it too: it is the cutoff of the composition-based redo's heaps
    // (`ResolvedBlastpArgs::inclusion_ethresh`).
    pub inclusion_ethresh: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1156-1160
    // The file is opened when a PSSM is written (`pssm.rs` `AsnPssmWriter`, as
    // `CAutoOutputFileReset`), not by the parser.
    pub out_pssm: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1161-1170
    // The file is opened when a PSSM is written (`pssm.rs` `write_ascii_pssm_files`, as
    // `CAutoOutputFileReset`), not by the parser.
    pub out_ascii_pssm: Option<String>,
    pub save_pssm_after_last_round: bool,
    pub save_each_pssm: bool,

    /// `BlastpArgs::first_query_ordinal`: the number of the first query of the records of a call
    /// that starts there (`api::run_batches`). 0 on the command line.
    pub first_query_ordinal: usize,
    /// LOSAT's refusal of a record of these records with the id of a query of an earlier call of
    /// the run (`query_batch::ScopeRefusal`), where the library searches each query batch as a call
    /// of its own (`api::run_batches`): the call's add of its batch gives it, so that the batch
    /// stops in the command line's form (line WQA D10 round 5, 2026-11-20). `None` on the command
    /// line. Not a public option.
    pub(crate) scope_refusal: Option<crate::blastinput::query_batch::ScopeRefusal>,
    /// The run's numbers of the files of `-save_each_pssm` (`pssm::PssmVersions`, NCBI's
    /// `CAutoOutputFileReset` of the run): a command line's run makes them; the library's query
    /// batches, each a call of its own, share the run's (`api::batches`).
    pub pssm_versions: crate::algorithm::psiblast::pssm::PssmVersions,
    /// The queries that the checks of the PSSM files before the search cover
    /// (`pssm::PssmCheckScope`): the call's on the command line; the run's in the first of the
    /// library's query batches, none in the later ones (`api::batches`).
    pub pssm_check: crate::algorithm::psiblast::pssm::PssmCheckScope,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/deltablast_args.cpp:44-120
    pub domain_inclusion_ethresh: f64,
    pub rpsdb: String,
    pub show_domain_hits: bool,
}

impl DeltablastArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, deltablast's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(DeltablastArgs {
            query: args.option_or_default("query")?,
            subject: args.given_option("subject")?,
            db: args.given_option("db")?,
            db_filter: crate::blastinput::db_filter::DbFilterArgs::from_args(args)?,
            dbsize: args.given_option("dbsize")?,
            show_gis: args.flag("show_gis"),
            parse_deflines: args.flag("parse_deflines"),
            searchsp: args.given_option("searchsp")?,
            sum_stats: args.given_option("sum_stats")?,
            evalue: args.given_option("evalue")?,
            threshold: args.given_option("threshold")?,
            word_size: args.given_option("word_size")?,
            num_threads: args.option_or_default("num_threads")?,
            out: args.given_option("out")?,
            // LOSAT refuses `-export_search_strategy` by name (`cli::reader::named_refusal`,
            // item PS15b-3).
            export_search_strategy: None,
            max_target_seqs: args.given_option("max_target_seqs")?,
            num_descriptions: args.given_option("num_descriptions")?,
            num_alignments: args.given_option("num_alignments")?,
            formatting: crate::blastinput::report_counts::FormattingArgs::from_args(args)?,
            max_hsps_per_subject: args.given_option("max_hsps")?,
            qcov_hsp_perc: args.given_option("qcov_hsp_perc")?,
            culling_limit: args.given_option("culling_limit")?,
            best_hit_overhang: args.given_option("best_hit_overhang")?,
            best_hit_score_edge: args.given_option("best_hit_score_edge")?,
            subject_besthit: args.flag("subject_besthit"),
            xdrop_ungap: args.given_option("xdrop_ungap")?,
            xdrop_gap: args.given_option("xdrop_gap")?,
            xdrop_gap_final: args.given_option("xdrop_gap_final")?,
            gap_trigger: args.option_or_default("gap_trigger")?,
            window_size: args.given_option("window_size")?,
            matrix: args.given_option("matrix")?,
            gap_open: args.given_option("gapopen")?,
            gap_extend: args.given_option("gapextend")?,
            comp_based_stats: args.option_or_default("comp_based_stats")?,
            seg: args.given_option("seg")?,
            soft_masking: args.given_option("soft_masking")?,
            lcase_masking: args.flag("lcase_masking"),
            use_sw_tback: args.flag("use_sw_tback"),
            // LOSAT refuses `-query_loc` by name (`cli::reader::named_refusal`, item PS13b).
            query_loc: None,
            subject_loc: args.given_option("subject_loc")?,
            outfmt: args.option_or_default("outfmt")?,
            num_iterations: args.option_or_default("num_iterations")?,
            pseudocount: args.option_or_default("pseudocount")?,
            inclusion_ethresh: args.option_or_default("inclusion_ethresh")?,
            out_pssm: args.given_option("out_pssm")?,
            out_ascii_pssm: args.given_option("out_ascii_pssm")?,
            save_pssm_after_last_round: args.flag("save_pssm_after_last_round"),
            save_each_pssm: args.flag("save_each_pssm"),
            first_query_ordinal: Default::default(),
            scope_refusal: Default::default(),
            pssm_versions: Default::default(),
            pssm_check: Default::default(),
            domain_inclusion_ethresh: args.option_or_default("domain_inclusion_ethresh")?,
            rpsdb: args.option_or_default("rpsdb")?,
            // LOSAT refuses `-show_domain_hits` by name (`cli::reader::named_refusal`, item PS13b).
            show_domain_hits: false,
        })
    }

    /// blastp's arguments for deltablast's search, after LOSAT's refusal of the output formats
    /// that it does not write yet (all but 0, 6 and 7, item PS15b-3; psiblast's reading of
    /// `-outfmt`, `PsiblastArgs::check_outfmt`, names deltablast in it).
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        PsiblastArgs::check_outfmt(
            &self.outfmt,
            self.out.as_deref(),
            self.export_search_strategy.is_some(),
            BlastpApp::Deltablast,
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
        let rps = RpsblastArgs {
            query: PathBuf::from("-"),
            db: Some(rpsdb.to_string()),
            out: None,
            export_search_strategy: None,
            evalue: Some(if ethresh > 0.0 { ethresh } else { 10.0 }),
            dbsize: None,
            searchsp: None,
            comp_based_stats: "1".to_string(),
            outfmt: "0".to_string(),
            show_gis: false,
            html: false,
            num_descriptions: None,
            num_alignments: None,
            formatting: Default::default(),
            max_target_seqs: None,
            max_hsps_per_subject: None,
            qcov_hsp_perc: None,
            xdrop_ungap: None,
            xdrop_gap: None,
            xdrop_gap_final: None,
            sum_stats: None,
            use_sw_tback: false,
            seg: None,
            soft_masking: None,
            window_size: None,
            lcase_masking,
            query_loc: None,
            parse_deflines: false,
            num_threads: 1,
            mt_mode: None,
        };
        let mut blastp = rps.into_blastp_args()?;
        if let Some(rps) = blastp.rps.as_mut() {
            rps.domain_search = true;
        }
        Ok(blastp)
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
}
