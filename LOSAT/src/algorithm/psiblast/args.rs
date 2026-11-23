//! Command-line arguments for PSIBLAST

use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

use crate::algorithm::blastp::args::{BlastpApp, BlastpArgs};
use crate::blastinput::value_parsers::*;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:58-116
// psiblast's arguments are blastp's without `-task` (its task is "psiblast"), `-ungapped`
// (no `CGappedArgs`), the database masks and `-mt_mode`, with the PSI-BLAST, PSSM engine and
// PHI-BLAST arguments. LOSAT runs psiblast's first round as blastp's search (`run`), its later
// rounds (`-num_iterations`) with the PSSM made from each round (`rounds.rs`), and writes the
// ASN.1 and ASCII PSSMs (`pssm.rs`); with `-in_pssm` the
// query is a PSSM, searched with the PSSM's stages (`in_pssm.rs`, `algorithm::pssm_search`), in
// rounds and with the PSSM files too (`in_pssm_files.rs`);
// `-in_msa` makes the PSSM from an alignment (`in_msa.rs`) and searches it as `-in_pssm` does;
// `-phi_pattern` is refused by name after the reader (its rows of `ncbi_app::args::refusal`).
// The values are
// read as NCBI's argument parser reads them (`value_parsers.rs`); the checks that NCBI makes
// after reading them are blastp's (`BlastpArgs::check_options`, with psiblast's program).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct PsiblastArgs {
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
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1231-1247
    // NCBI's argument reader reports the exclusions (`cli::reader`). The file is opened and read with
    // the options (`algorithm::psiblast::in_pssm`); `-phi_pattern` is refused by name after the reader (its rows of
    // `ncbi_app::args::refusal`).
    #[arg(
        long = "in_pssm",
        value_name = "psi_chkpt_file",
        value_parser = ncbi_path(),
        help = "PSI-BLAST checkpoint file"
    )]
    pub in_pssm: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1184-1229
    // `msa_exclusions` are `-in_pssm`, `-query`, `-query_loc` and `-phi_pattern`; the last is
    // refused by name (its rows of `ncbi_app::args::refusal`) after NCBI's argument reader, which reports the
    // exclusions and requirements (`cli::reader`). The file is read with the options (`algorithm::psiblast::in_msa`).
    #[arg(
        long = "in_msa",
        value_name = "align_restart",
        value_parser = ncbi_path(),
        help = "File name of multiple sequence alignment to restart PSI-BLAST"
    )]
    pub in_msa: Option<PathBuf>,
    #[arg(
        long = "msa_master_idx",
        value_name = "index",
        value_parser = ncbi_integer,
        help = "Ordinal number (1-based index) of the sequence to use as a master in the multiple sequence alignment. If not provided, the first sequence in the multiple sequence alignment will be used"
    )]
    pub msa_master_idx: Option<i32>,
    #[arg(
        long = "ignore_msa_master",
        help = "Ignore the master sequence when creating PSSM"
    )]
    pub ignore_msa_master: bool,

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
        let request_inputs = self
            .export_search_strategy
            .is_some()
            .then(|| Box::new(self.option_inputs()));
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
            strategy: None,
            request_inputs,
            app: BlastpApp::Psiblast,
            taxdb_setting: None,
            matrix_search: None,
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
        let mut inputs = OptionInputs::new(RequestKind::Psiblast, "blastp");
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
    /// `-outfmt`.
    pub fn check_outfmt(outfmt: &str, out: Option<&std::path::Path>, strategy: bool) -> Result<()> {
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
        if let Some(format) = ReportFormat::from_number(choice.number) {
            if !matches!(
                format,
                ReportFormat::Pairwise
                    | ReportFormat::Tabular
                    | ReportFormat::TabularWithComments
                    | ReportFormat::CommaSeparatedValues
                    | ReportFormat::CommaSeparatedValuesWithHeader
            ) {
                bail!(
                    "output format {} is not supported by LOSAT's PSIBLAST",
                    choice.number
                );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "args_tests.rs"]
mod tests;
