//! Command-line arguments for PSIBLAST

use anyhow::{bail, Result};
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:58-116
// psiblast's arguments are blastp's without `-task` (its task is "psiblast"), `-ungapped`
// (no `CGappedArgs`), the database masks and `-mt_mode`, with the PSI-BLAST, PSSM engine and
// PHI-BLAST arguments. LOSAT runs psiblast's first round as blastp's search (`run`), its later
// rounds (`-num_iterations`) with the PSSM made from each round (`rounds.rs`), and writes the
// ASN.1 and ASCII PSSMs (`pssm.rs`); with `-in_pssm` the
// query is a PSSM, searched with the PSSM's stages (`in_pssm.rs`, `algorithm::pssm_search`), in
// rounds and with the PSSM files too (`in_pssm_files.rs`);
// `-in_msa` makes the PSSM from an alignment (`in_msa.rs`) and searches it as `-in_pssm` does;
// `-phi_pattern` runs round 1 as a PHI-BLAST search (`algorithm::phi_search`; rounds after it are a
// later item's refusal, `phi_rounds_refusal`).
// The values are those that NCBI's argument reader read (`ncbi_app::args`, `from_args`); the
// checks that NCBI makes after reading them are blastp's (`BlastpArgs::check_options`, with
// psiblast's program).
#[derive(Debug, Clone)]
pub struct PsiblastArgs {
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2692-2693
    // NCBI's `CBlastFormat::m_IsHTML` (blast_format.cpp:90,224).
    pub html: bool,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1231-1247
    // NCBI's argument reader reports the exclusions (`cli::reader`), `-phi_pattern`'s too. The file is
    // opened and read with the options (`algorithm::psiblast::in_pssm`).
    pub in_pssm: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1184-1229
    // `msa_exclusions` are `-in_pssm`, `-query`, `-query_loc` and `-phi_pattern`; NCBI's argument
    // reader reports the exclusions and requirements (`cli::reader`). The file is read with the
    // options (`algorithm::psiblast::in_msa`).
    pub in_msa: Option<PathBuf>,
    pub msa_master_idx: Option<i32>,
    pub ignore_msa_master: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1366-1380
    // NCBI's argument reader reports its exclusions (`cli::reader`). The file is read with the
    // options (`algorithm::phi_search::file`) and its pattern compiled before the search
    // (`ResolvedBlastpArgs::phi`).
    pub phi_pattern: Option<PathBuf>,

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
}

impl PsiblastArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, psiblast's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(PsiblastArgs {
            query: args.option_or_default("query")?,
            subject: args.given_option("subject")?,
            db: args.given_option("db")?,
            db_filter: crate::blastinput::db_filter::DbFilterArgs::from_args_with_ipg(args)?,
            dbsize: args.given_option("dbsize")?,
            show_gis: args.flag("show_gis"),
            html: args.flag("html"),
            parse_deflines: args.flag("parse_deflines"),
            searchsp: args.given_option("searchsp")?,
            sum_stats: args.given_option("sum_stats")?,
            evalue: args.given_option("evalue")?,
            threshold: args.given_option("threshold")?,
            word_size: args.given_option("word_size")?,
            num_threads: args.option_or_default("num_threads")?,
            out: args.given_option("out")?,
            export_search_strategy: args.given_option("export_search_strategy")?,
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
            query_loc: args.given_option("query_loc")?,
            subject_loc: args.given_option("subject_loc")?,
            outfmt: args.option_or_default("outfmt")?,
            num_iterations: args.option_or_default("num_iterations")?,
            pseudocount: args.option_or_default("pseudocount")?,
            inclusion_ethresh: args.option_or_default("inclusion_ethresh")?,
            out_pssm: args.given_option("out_pssm")?,
            out_ascii_pssm: args.given_option("out_ascii_pssm")?,
            save_pssm_after_last_round: args.flag("save_pssm_after_last_round"),
            save_each_pssm: args.flag("save_each_pssm"),
            in_pssm: args.given_option("in_pssm")?,
            in_msa: args.given_option("in_msa")?,
            msa_master_idx: args.given_option("msa_master_idx")?,
            ignore_msa_master: args.flag("ignore_msa_master"),
            phi_pattern: args.given_option("phi_pattern")?,
            first_query_ordinal: Default::default(),
            scope_refusal: Default::default(),
            pssm_versions: Default::default(),
            pssm_check: Default::default(),
        })
    }
}

impl PsiblastArgs {
    /// blastp's arguments for psiblast's search, after LOSAT's refusals of what it does not run
    /// (the output formats of the later items of the series).
    pub fn into_blastp_args(self) -> Result<BlastpArgs> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1277-1281
        // A query given as a PSSM (`-in_pssm`) runs its rounds as a FASTA query's, round 1 with the
        // input PSSM (item PS5b), and a range of the query (`-query_loc`) runs them as NCBI does:
        // the later rounds search the range with the PSSM of the whole record
        // (`PsiSearch::from_engine`, series Q, Q4a).
        Self::check_outfmt(
            &self.outfmt,
            self.out.as_deref(),
            self.export_search_strategy.is_some(),
            self.phi_pattern.is_some(),
        )?;
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248 (cited at
        // `in_pssm_files`: the loop writes the PSSM of each round at its start, the input PSSM of
        // `-in_pssm` in round 1)
        // The round driver runs `-num_iterations` other than 1, and with `-in_pssm` one round too
        // when a PSSM file is written (the input PSSM before the search).
        let rounds = (self.num_iterations != 1
            || ((self.in_pssm.is_some() || self.in_msa.is_some())
                && (self.out_pssm.is_some() || self.out_ascii_pssm.is_some())))
        .then(|| u32::try_from(self.num_iterations).unwrap_or(0));
        // The options of the search strategy's request and of the archive's (outfmt 11), which
        // record psiblast's own (`service "psi"`, the PSSM engine's options).
        let request_inputs = Some(Box::new(self.option_inputs()));
        Ok(BlastpArgs {
            query: self.query,
            subject: self.subject,
            db: self.db,
            // blastp's database filters (line DI-8), which every round searches with.
            db_filter: self.db_filter,
            // psiblast has no database masks (NCBI's unknown arguments).
            db_mask: Default::default(),
            dbsize: self.dbsize,
            show_gis: self.show_gis,
            html: self.html,
            parse_deflines: self.parse_deflines,
            id_scope: None,
            searchsp: self.searchsp,
            lcase_masking: self.lcase_masking,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:46-50
            // psiblast's options are those of blastp's task (`CBlastAdvancedProteinOptionsHandle`)
            // with no filtering (`SetQueryOptionDefaults`, as blastp's).
            task: "blastp".to_string(),
            evalue: self.evalue,
            threshold: self.threshold,
            word_size: self.word_size,
            num_threads: self.num_threads,
            // psiblast has no `-mt_mode` (its rows of `ncbi_app::args::refusal`).
            mt_mode: None,
            out: self.out,
            max_target_seqs: self.max_target_seqs,
            num_descriptions: self.num_descriptions,
            num_alignments: self.num_alignments,
            max_hsps_per_subject: self.max_hsps_per_subject,
            qcov_hsp_perc: self.qcov_hsp_perc,
            culling_limit: self.culling_limit,
            best_hit_overhang: self.best_hit_overhang,
            best_hit_score_edge: self.best_hit_score_edge,
            subject_besthit: self.subject_besthit,
            ungapped: false,
            xdrop_ungap: self.xdrop_ungap,
            xdrop_gap: self.xdrop_gap,
            xdrop_gap_final: self.xdrop_gap_final,
            gap_trigger: Some(self.gap_trigger),
            sum_stats: self.sum_stats == Some(true),
            window_size: self.window_size,
            matrix: self.matrix,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            comp_based_stats: self.comp_based_stats,
            seg: self.seg,
            soft_masking: self.soft_masking,
            use_sw_tback: self.use_sw_tback,
            query_loc: self.query_loc,
            subject_loc: self.subject_loc,
            outfmt: self.outfmt,
            archive: None,
            export_search_strategy: self.export_search_strategy,
            import_search_strategy: None,
            import: None,
            given: Default::default(),
            strategy: None,
            request_inputs,
            app: BlastpApp::Psiblast,
            taxdb_setting: None,
            matrix_search: None,
            write_verification: None,
            read_verification: Default::default(),
            inclusion_ethresh: Some(self.inclusion_ethresh),
            first_query_ordinal: self.first_query_ordinal,
            scope_refusal: self.scope_refusal.clone(),
            psiblast_pssm: Some(crate::algorithm::psiblast::pssm::PsiblastPssmOptions {
                save_pssm_after_last_round: self.save_pssm_after_last_round,
                save_each_pssm: self.save_each_pssm,
                out_pssm: self.out_pssm,
                out_ascii_pssm: self.out_ascii_pssm,
                pseudocount: self.pseudocount,
                rounds,
                in_pssm: self.in_pssm,
                in_msa: self
                    .in_msa
                    .map(|path| crate::algorithm::psiblast::in_msa::MsaInput {
                        path,
                        master_index: self.msa_master_idx,
                        ignore_master: self.ignore_msa_master,
                    }),
                versions: self.pssm_versions,
                check_scope: self.pssm_check,
                phi_pattern: self.phi_pattern,
            }),
            // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:502
            // The line length only: `-sorthits` and `-sorthsps` are read with their warnings
            // (`BlastpArgs::check_options`) and sort nothing (`ResolvedBlastpArgs::sort`).
            formatting: self.formatting,
            rps: None,
        })
    }

    /// The options that psiblast's search strategy requests record (`options::algorithm_options`
    /// of `RequestKind::Psiblast`: NCBI's argument classes of `psiblast_args.cpp` and their
    /// conditions), from the arguments as given. `HitlistSize` is the request's
    /// (`StrategyExport::for_run`); the database lists are the database's
    /// (`OptionsHandle::set_database_filter`).
    pub(crate) fn option_inputs(&self) -> crate::report::archive::options::OptionInputs {
        use crate::report::archive::options::{seg_filter, OptionInputs, RequestKind};
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:123-131
        // `-phi_pattern` makes the handle PHI-BLAST's (`RequestKind::Phiblastp`), whose request
        // records the pattern as read; a file that does not read stops the run at the PHI
        // handler with NCBI's error (`BlastpArgs::check_options`), after the strategy file is
        // created and before it is written.
        let kind = match self.phi_pattern {
            Some(_) => RequestKind::Phiblastp,
            None => RequestKind::Psiblast,
        };
        let mut inputs = OptionInputs::new(kind, "blastp");
        inputs.phi_pattern = self
            .phi_pattern
            .as_deref()
            .and_then(crate::algorithm::phi_search::file::pattern_text);
        inputs.dbsize = self.dbsize;
        inputs.searchsp = self.searchsp;
        inputs.sum_stats = self.sum_stats;
        inputs.evalue = self.evalue;
        inputs.gapopen = self.gap_open;
        inputs.gapextend = self.gap_extend;
        inputs.matrix = self.matrix.clone();
        // `-xdrop_ungap` is never recorded (`SetXDropoff` has no remote option).
        inputs.xdrop_gap = self.xdrop_gap;
        inputs.xdrop_gap_final = self.xdrop_gap_final;
        inputs.word_size = self.word_size;
        inputs.qcov_hsp_perc = self.qcov_hsp_perc;
        inputs.max_hsps = self.max_hsps_per_subject.map(|count| count as i32);
        inputs.culling_limit = self.culling_limit;
        inputs.best_hit_overhang = self.best_hit_overhang;
        inputs.best_hit_score_edge = self.best_hit_score_edge;
        inputs.subject_besthit = self.subject_besthit;
        // `static_cast<int>` of the double, as BLASTP's (`blastp::archive::option_inputs`).
        inputs.threshold = self
            .threshold
            .map(crate::core::blast_util::ncbi_int4_from_double);
        inputs.window_size = self.window_size;
        inputs.soft_masking = self.soft_masking;
        inputs.seg = self.seg.as_deref().and_then(seg_filter);
        if let Ok((mode, unified_p)) =
            crate::blastinput::app::parse_comp_based_stats(&self.comp_based_stats, true, false)
        {
            inputs.comp_based_stats = mode as i32;
            inputs.unified_p = unified_p;
        }
        inputs.use_sw_tback = self.use_sw_tback;
        inputs.gap_trigger = self.gap_trigger;
        inputs.pseudocount = self.pseudocount;
        inputs.inclusion_ethresh = self.inclusion_ethresh;
        inputs.num_iterations = Some(self.num_iterations);
        inputs.ignore_msa_master = self.ignore_msa_master;
        inputs
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2744-2748
    /// NCBI parses `-outfmt` before the option handlers; its own errors come from that parse
    /// and from blastp's (`run`). The formats that NCBI runs and LOSAT's psiblast does not
    /// write yet are refused here, with `-out` `out`, but a value that NCBI refuses itself later
    /// (`app::ncbi_refuses_format`: outfmt 17, a custom delimiter that the field list uses,
    /// outfmt 13 and 14 without a file name) goes on to NCBI's error at its place (blastp's).
    /// The tabular reports with a field list, and outfmt 10 and 20, are blastp's
    /// (`blastp::tabular`), without a custom delimiter, which psiblast's formatter does not get
    /// (psiblast_app.cpp:482-505, `report::tabular::TabularSpec::of_app`). `strategy`: the run
    /// writes a search strategy file (`-export_search_strategy`), before which NCBI reads no
    /// `-outfmt`. `phi`: the run has `-phi_pattern`, whose formats 8, 9 and 12 end with NCBI's
    /// engine error (`blast_engine/phi_hooks.rs`); its refusal of formats 5 and 11 comes after
    /// NCBI's errors of the pattern's search (`phi_hooks::phi_limits`).
    pub fn check_outfmt(
        outfmt: &str,
        out: Option<&std::path::Path>,
        strategy: bool,
        phi: bool,
    ) -> Result<()> {
        use crate::blastinput::app::{Program, ReportFormat};
        // A value not read yet (`EarlyFormat::Deferred`) is the formatting handler's.
        let Some(choice) = crate::blastinput::app::early_format_choice(outfmt, strategy)?
            .parsed()
            .cloned()
        else {
            return Ok(());
        };
        if crate::blastinput::app::ncbi_refuses_format(&choice, Program::Blastp, out) {
            return Ok(());
        }
        if ReportFormat::from_number(choice.number).is_some() {
            // With a pattern: PHI-BLAST's reports (`report_write::write_phi_report`), formats
            // 8, 9 and 12 (which reach NCBI's engine error), and 5 and 11 (`phi_limits`).
            // Without a pattern: `writes_outfmt` (0..=12, 18, 20).
            let allowed = if phi {
                matches!(choice.number, 0 | 5..=12 | 20)
            } else {
                writes_outfmt(choice.number)
            };
            if !allowed {
                bail!(
                    "output format {} is not supported by LOSAT's PSIBLAST yet",
                    choice.number
                );
            }
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:265-287
/// Whether LOSAT's psiblast writes the output format `number` for its queries (FASTA, `-in_pssm`
/// and `-in_msa`, at any `-num_iterations`, `-db` and `-subject`): the one rule of the command
/// line's refusal (`PsiblastArgs::check_outfmt`). NCBI writes every round through the formatter of
/// the format (`WriteArchive`, or `PrintOneResultSet` with the round's number); LOSAT writes all
/// but the XML2 and JSON documents of outfmt 13 to 16 (item PS15b-2 of line PS, its second part).
pub(crate) fn writes_outfmt(number: i32) -> bool {
    matches!(number, 0..=12 | 18 | 20)
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
