//! Pinned NCBI TBLASTN local-subject search and Stage E report options.
use anyhow::{Context, Result};
use clap::Args;
use std::io::Write;
use std::path::PathBuf;

use crate::blastinput::query_batch::query_batches;
use crate::blastinput::value_parsers::*;
use crate::utils::genetic_code::GeneticCode;

use super::stage_d_pipeline::{
    run_local_for_report_threads, LocalStageDProfile, LocalStageDScoring,
};
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::render;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::config::ScoringMatrix;

// NCBI reference: c++/src/algo/blast/blastinput/tblastn_args.cpp:54-61
// ```c++
//     static const char kDefaultTask[] = "tblastn";
//     SetTask(kDefaultTask);
//     set<string> tasks;
//     tasks.insert(kDefaultTask);
//     tasks.insert("tblastn-fast");
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
// ```
// NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-85
// ```c++
// SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
// m_Opts->SetSumStatisticsMode();
// m_Opts->SetCompositionBasedStats(eCompositionMatrixAdjust);
// SetDbGeneticCode(BLAST_GENETIC_CODE);
// ```
// The arguments are read with NCBI's grammar (`blastinput/value_parsers.rs`); the values
// that NCBI checks after parsing are checked by `check_options`, in NCBI's order. The
// NCBI options that LOSAT's TBLASTN does not implement are rejected by the parser
// (`cli.rs`).
#[derive(Args, Debug, Clone)]
#[command(rename_all = "snake_case")]
pub struct TblastnArgs {
    // NCBI blast_args.cpp:3425-3427:
    // arg_desc.AddDefaultKey(kArgQuery, "input_file", "Input file name",
    //                        CArgDescriptions::eInputFile, kDfltArgQuery);
    // The default `-` is standard input.
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH", default_value = "-")]
    pub query: PathBuf,
    #[arg(long, value_parser = ncbi_input_path(), value_name = "PATH")]
    pub subject: Option<PathBuf>,
    #[arg(long, default_value = "tblastn", value_parser = ["tblastn", "tblastn-fast"])]
    pub task: String,
    #[arg(long, value_name = "PATH", value_parser = ncbi_output_path())]
    pub out: Option<PathBuf>,
    // An omitted value keeps the default of the task (10).
    #[arg(long, value_parser = tblastn_real, help = "Expectation value (E) threshold for saving hits [default: 10]")]
    pub evalue: Option<f64>,
    #[arg(long, value_parser = protein_word_size, help = "Word size for wordfinder algorithm [default: 3; 5 for tblastn-fast]")]
    pub word_size: Option<i32>,
    // NCBI blast_args.cpp:258-273: when omitted, each cost keeps the selected matrix's
    // BLAST_MATRIX_BEST value.
    #[arg(long = "gapopen", value_parser = ncbi_integer, help = "Cost to open a gap [default: 11, or the matrix's best value]")]
    pub gap_open: Option<i32>,
    #[arg(long = "gapextend", value_parser = ncbi_integer, help = "Cost to extend a gap [default: 1, or the matrix's best value]")]
    pub gap_extend: Option<i32>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1029-1056
    // ```c++
    // arg_desc.AddDefaultKey(kArgDbGeneticCode, "int_value",
    //                        "Genetic code to use to translate "
    //                        "database/subjects (see user manual for details)\n",
    //                        CArgDescriptions::eInteger,
    //                        NStr::IntToString(BLAST_GENETIC_CODE));
    // opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
    // ```
    // PD-TLOSAN-LOCAL-GENCODE-32 adds gc.prt ID 32 only for TBLASTN.
    #[arg(long, default_value_t = 1, value_parser = tblastn_genetic_code)]
    pub db_gencode: u8,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments [default: 0]")]
    pub max_intron_length: Option<i32>,
    #[arg(long, help = "Scoring matrix name [default: BLOSUM62]")]
    pub matrix: Option<String>,
    #[arg(long, value_parser = tblastn_threshold_value, help = "Minimum word score such that the word is added to the BLAST lookup table [default: 13, or the matrix's suggestion; 19.3, 21 or 20.25 for word sizes 5, 6 and 7; 20 for tblastn-fast]")]
    pub threshold: Option<f64>,
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:220-229,280-286:
    // arg_desc.AddOptionalKey(kArgGappedXDropoff, "float_value", ...);
    // arg_desc.AddOptionalKey(kArgFinalGappedXDropoff, "float_value", ...);
    // opt.SetGapXDropoff(args[kArgGappedXDropoff].AsDouble());
    // opt.SetGapXDropoffFinal(args[kArgFinalGappedXDropoff].AsDouble());
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for preliminary gapped extensions [default: 15]")]
    pub xdrop_gap: Option<f64>,
    #[arg(long, value_parser = tblastn_real, help = "X-dropoff value (in bits) for final gapped alignment [default: 25]")]
    pub xdrop_gap_final: Option<f64>,
    // Read when NCBI's handler reads it (`blastinput/app.rs` `parse_comp_based_stats`).
    #[arg(
        long,
        default_value = "2",
        help = "Use composition-based statistics: 0 (F, f) none, 1 composition-based statistics, 2 (D, d, T, t) conditional compositional score matrix adjustment, 3 unconditional compositional score matrix adjustment"
    )]
    pub comp_based_stats: String,
    #[arg(long, default_value = "0")]
    pub outfmt: String,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1945-1949
    // ```c++
    //     // query location
    //     arg_desc.AddOptionalKey(kArgQueryLocation, "range",
    //                             "Location on the query sequence in 1-based offsets "
    //                             "(Format: start-stop)",
    //                             CArgDescriptions::eString);
    // ```
    // Read by the query options handler (`check_options`).
    #[arg(
        long = "query_loc",
        value_name = "RANGE",
        help = "Location on the query sequence in 1-based offsets (Format: start-stop)"
    )]
    pub query_loc: Option<String>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2372-2376
    // ```c++
    //         // subject location
    //         arg_desc.AddOptionalKey(kArgSubjectLocation, "range",
    //                         "Location on the subject sequence in 1-based offsets "
    //                         "(Format: start-stop)",
    //                         CArgDescriptions::eString);
    // ```
    // Read by the database arguments handler when it reads the subjects (`run`).
    #[arg(
        long = "subject_loc",
        value_name = "RANGE",
        help = "Location on the subject sequence in 1-based offsets (Format: start-stop)"
    )]
    pub subject_loc: Option<String>,
    // Read when NCBI's filtering handler reads it (`blastinput/app.rs` `parse_seg_option`);
    // an omitted value is NCBI's default for tblastn, "12 2.2 2.5".
    #[arg(
        long,
        help = "SEG: no, yes, or \"WINDOW LOCUT HICUT\" [default: 12 2.2 2.5]"
    )]
    pub seg: Option<String>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942,2549-2557
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //     "Use lower case filtering in query and subject sequence(s)?", true);
    // ReadSequencesToBlast(..., use_lcase_masks, subjects, ...);
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // ```c++
    //     arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                            "Apply filtering locations as soft masks",
    //                            CArgDescriptions::eBoolean,
    //                            kDfltArgLookupTableMaskingOnlyProt);
    // ```
    #[arg(long, value_parser = ncbi_boolean, help = "Apply filtering locations as soft masks [default: false]")]
    pub soft_masking: Option<bool>,
    #[arg(long, value_parser = ncbi_boolean, help = "Use sum statistics [default: true]")]
    pub sum_stats: Option<bool>,
    #[arg(long, value_parser = nonnegative_ncbi_integer, help = "Multiple hits window size, use 0 to specify 1-hit algorithm [default: 40, or the matrix's suggestion]")]
    pub window_size: Option<i32>,
    // An omitted value keeps the default hit list size (500).
    #[arg(long, value_parser = blastn_count, help = "Maximum number of aligned sequences to keep [default: 500]")]
    pub max_target_seqs: Option<usize>,
    #[arg(long, default_value_t = 1, value_parser = blastn_count)]
    pub num_threads: usize,
    #[arg(long)]
    pub ungapped: bool,
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:997-1056
// ```c++
// static int gcs[] = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,
//                     21,22,23,24,25,26,27,28,29,30,31,33};
// return (genetic_codes.find(val) != genetic_codes.end());
// ```
// Stage A product decision admits 32 from gc.prt:340-347 for TBLASTN only.
pub fn tblastn_genetic_code(value: &str) -> std::result::Result<u8, String> {
    let id = value.parse::<u8>().map_err(|_| "invalid genetic code ID")?;
    GeneticCode::try_from_id(id).map(|_| id)
}

/// The options of a TBLASTN search as NCBI's option handlers set them.
#[derive(Debug, Clone)]
pub struct ResolvedTblastnArgs {
    pub task: String,
    pub evalue: f64,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub xdrop_gap: f64,
    pub xdrop_gap_final: f64,
    pub word_size: i32,
    /// Whether the compressed-alphabet lookup table is used (word sizes 5-7, tblastn-fast).
    pub compressed_lookup: bool,
    pub sum_stats: bool,
    pub db_gencode: u8,
    pub gapped: bool,
    pub max_intron_length: i32,
    pub seg: SegSpec,
    pub soft_masking: bool,
    pub lcase_masking: bool,
    /// The matrix name as typed (the report prints it so).
    pub matrix_name: String,
    pub threshold: f64,
    pub window_size: i32,
    pub hitlist_size: usize,
    /// The -max_target_seqs value, which also sets the numbers of descriptions and
    /// alignments of the pairwise report (500 and 250 when it is not given).
    pub max_target_seqs_given: Option<usize>,
    pub composition_mode: crate::blastinput::app::CompositionMode,
    /// The `-comp_based_stats` value as typed.
    pub comp_based_stats: String,
    pub num_threads: usize,
    /// The `-query_loc` range, as given.
    pub query_range: Option<crate::blastinput::seq_range::SequenceRange>,
}

impl TblastnArgs {
    /// The options as `check_options` sets them, without its warnings.
    pub fn resolve(&self) -> Result<ResolvedTblastnArgs> {
        self.check_options(&mut std::io::sink())
    }

    /// NCBI's `CBlastAppArgs::SetOptions` for tblastn after the files are opened: the task's
    /// options, each argument handler in the order of `CTblastnAppArgs` (with the warnings
    /// that they post to `diagnostics`), and `Validate`, with NCBI's errors.
    ///
    /// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
    /// ```c
    ///     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    ///         (*arg)->ExtractAlgorithmOptions(args, opts);
    ///     }
    ///
    ///     m_IsUngapped = !opts.GetGappedMode();
    ///     try { retval->Validate(); }
    ///     catch (const CBlastException& e) {
    ///         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    ///     }
    /// ```
    pub fn check_options(&self, diagnostics: &mut dyn Write) -> Result<ResolvedTblastnArgs> {
        use crate::blastinput::app::{
            formatting_handler_check, parse_comp_based_stats, parse_formatting_string,
            parse_seg_option,
        };
        use crate::stats::protein_options::{
            protein_gap_existence_extend_params, suggested_threshold, suggested_window_size,
            validate_protein_options, ProteinOptionsCheck, SuggestionProgram,
        };
        // NCBI reference: c++/src/algo/blast/api/blast_options_handle.cpp:436-447
        // ```c
        //     else if (!NStr::CompareNocase(task, "tblastn") ||
        //     		 !NStr::CompareNocase(task, "tblastn-fast"))
        //     {
        //     	CTBlastnOptionsHandle* opts =
        //     	            dynamic_cast<CTBlastnOptionsHandle*>
        //          (CBlastOptionsFactory::Create(eTblastn, locality));
        //     	if(task == "tblastn-fast") {
        //     		opts->SetWordSize(5);
        //             opts->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
        //             opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTP_FAST);
        //     	}
        //     	retval = opts;
        // ```
        // NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-85 (the handle's
        // defaults: threshold 13, word size 3, window 40, X-drops 15 and 25, BLOSUM62
        // 11/1, e-value 10, sum statistics, composition mode 2, genetic code 1)
        // ```c++
        // SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
        // m_Opts->SetSumStatisticsMode();
        // m_Opts->SetCompositionBasedStats(eCompositionMatrixAdjust);
        // SetDbGeneticCode(BLAST_GENETIC_CODE);
        // ```
        let fast = self.task == "tblastn-fast";
        let mut word_size: i32 = if fast { 5 } else { 3 };
        let mut compressed_lookup = fast;
        let mut threshold: f64 = if fast { 20.0 } else { 13.0 };
        let mut evalue = 10.0;
        let (mut gap_open, mut gap_extend) = (11, 1);
        let (mut xdrop_gap, mut xdrop_gap_final) = (15.0, 25.0);
        let mut sum_stats = true;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:251-322 (the
        // generic search handler)
        // ```c
        //     if (args.Exist(kArgEvalue) && args[kArgEvalue]) {
        //         opt.SetEvalueThreshold(args[kArgEvalue].AsDouble());
        //     }
        //
        //     int gap_open=0, gap_extend=0;
        //     if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
        //          BLAST_GetProteinGapExistenceExtendParams
        //              (args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);
        // ...
        //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
        //         if (m_QueryIsProtein && args[kArgWordSize].AsInteger() > 4){
        //            opt.SetLookupTableType(eCompressedAaLookupTable);
        //            opt.SetWordThreshold(19.3);
        // ...
        //     if (args.Exist(kArgSumStats) && args[kArgSumStats]) {
        //         opt.SetSumStatisticsMode(args[kArgSumStats].AsBoolean());
        //     }
        // ```
        if let Some(value) = self.evalue {
            evalue = value;
        }
        let (matrix_gap_open, matrix_gap_extend) = self
            .matrix
            .as_deref()
            .and_then(protein_gap_existence_extend_params)
            .unwrap_or((0, 0));
        match (self.gap_open, &self.matrix) {
            (Some(value), _) => gap_open = value,
            (None, Some(_)) => gap_open = matrix_gap_open,
            (None, None) => {}
        }
        match (self.gap_extend, &self.matrix) {
            (Some(value), _) => gap_extend = value,
            (None, Some(_)) => gap_extend = matrix_gap_extend,
            (None, None) => {}
        }
        if let Some(value) = self.xdrop_gap {
            xdrop_gap = value;
        }
        if let Some(value) = self.xdrop_gap_final {
            xdrop_gap_final = value;
        }
        if let Some(value) = self.word_size {
            if value > 4 {
                compressed_lookup = true;
                threshold = 19.3;
                if value > 5 {
                    threshold = 21.0;
                }
                if value > 6 {
                    threshold = 20.25;
                }
            }
            word_size = value;
            // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:612-619
            // ```c
            // CBlastOptionsLocal::SetWordSize(int ws)
            // {
            //     m_LutOpts->word_size = ws;
            //     if (m_LutOpts->lut_type == eCompressedAaLookupTable && ws <= 4)
            // 	m_LutOpts->lut_type = eAaLookupTable;
            //     else if (m_LutOpts->lut_type == eAaLookupTable && ws > 4)
            // 	m_LutOpts->lut_type = eCompressedAaLookupTable;
            // }
            // ```
            // (tblastn-fast with -word_size 3 or 4 uses the plain lookup table.)
            if compressed_lookup && value <= 4 {
                compressed_lookup = false;
            }
        }
        if let Some(value) = self.sum_stats {
            sum_stats = value;
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:387-408 (the
        // filtering handler; the default -seg of tblastn is "12 2.2 2.5")
        // ```c
        //     if (args[kArgLookupTableMaskingOnly]) {
        //         opt.SetMaskAtHash(args[kArgLookupTableMaskingOnly].AsBoolean());
        //     }
        //
        //     vector<string> tokens;
        //
        //     try {
        //         if (m_QueryIsProtein && args[kArgSegFiltering]) {
        //             const string& seg_opts = args[kArgSegFiltering].AsString();
        // ```
        let soft_masking = self.soft_masking.unwrap_or(false);
        let seg = parse_seg_option(self.seg.as_deref().unwrap_or("12 2.2 2.5"), "TBLASTN")?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:586-623 (the word
        // threshold handler: the matrix's suggestion replaces only tblastn's default 13)
        // ```c
        //     } else if (program == eTblastn &&
        //                word_threshold != BLAST_WORD_THRESHOLD_TBLASTN) {
        //         retval = false;
        // ...
        //     if (args[kArgWordScoreThreshold]) {
        //         opt.SetWordThreshold(args[kArgWordScoreThreshold].AsDouble());
        //     } else if (s_IsDefaultWordThreshold(opt.GetProgram(),
        //                                         opt.GetWordThreshold())) {
        // ```
        let matrix_name = self
            .matrix
            .clone()
            .unwrap_or_else(|| "BLOSUM62".to_string());
        if let Some(value) = self.threshold {
            threshold = value;
        } else if threshold as i32 == 13 {
            threshold = suggested_threshold(SuggestionProgram::TranslatedSubject, &matrix_name);
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:485-497
        // ```c
        // CWindowSizeArg::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& opt)
        // {
        //     if (args[kArgWindowSize]) {
        //         opt.SetWindowSize(args[kArgWindowSize].AsInteger());
        //     } else {
        //         int window = -1;
        //         BLAST_GetSuggestedWindowSize(opt.GetProgramType(),
        //                                      opt.GetMatrixName(),
        //                                      &window);
        // ```
        let window_size = self
            .window_size
            .unwrap_or_else(|| suggested_window_size(&matrix_name));
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
        // ```c++
        //     // set the sequence range
        //     if (args.Exist(kArgQueryLocation) && args[kArgQueryLocation]) {
        //         m_Range = ParseSequenceRange(args[kArgQueryLocation].AsString(),
        //                                      "Invalid specification of query location");
        //     }
        // ```
        // The query options handler comes after the window size handler and before the
        // formatting handler (tblastn_args.cpp:44-135).
        let query_range = crate::blastinput::seq_range::parse_optional_range(
            self.query_loc.as_deref(),
            crate::blastinput::seq_range::RangeRole::Query,
            "TBLASTN",
        )?;
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886 and
        // 2960-2977 (the formatting handler: other programs' formats, the hit list size
        // and its warning)
        // ```c
        //     if(hitlist_size < 5){
        //    		ERR_POST(Warning << "Examining 5 or more matches is recommended");
        //     }
        // ```
        let choice = parse_formatting_string(&self.outfmt)?;
        formatting_handler_check(&choice, false)?;
        if self.max_target_seqs.is_some_and(|size| size < 5) {
            diagnostics.write_all(&crate::report::query_warnings::few_matches_warning(
                "tblastn",
            ))?;
        }
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:893-903 (the
        // composition-based statistics handler; `blastinput/app.rs`)
        // ```c
        //     if (args[kArgCompBasedStats]) {
        //         unique_ptr<bool> ungapped(args.Exist(kArgUngapped)
        //             ? new bool(args[kArgUngapped]) : 0);
        //         s_SetCompositionBasedStats(opt,
        //                                    args[kArgCompBasedStats].AsString(),
        //                                    args[kArgUseSWTraceback],
        //                                    ungapped.get());
        //     }
        // ```
        let (composition_mode, _) =
            parse_comp_based_stats(&self.comp_based_stats, false, self.ungapped)?;
        validate_protein_options(&ProteinOptionsCheck {
            gapped: !self.ungapped,
            gap_open,
            gap_extend,
            matrix_name: &matrix_name,
            threshold,
            word_size,
            compressed_lookup,
            evalue,
        })?;
        Ok(ResolvedTblastnArgs {
            task: self.task.clone(),
            evalue,
            gap_open,
            gap_extend,
            xdrop_gap,
            xdrop_gap_final,
            word_size,
            compressed_lookup,
            sum_stats,
            db_gencode: self.db_gencode,
            gapped: !self.ungapped,
            max_intron_length: self.max_intron_length.unwrap_or(0),
            seg,
            soft_masking,
            lcase_masking: self.lcase_masking,
            matrix_name,
            threshold,
            window_size,
            hitlist_size: self.max_target_seqs.unwrap_or(500),
            max_target_seqs_given: self.max_target_seqs,
            composition_mode,
            comp_based_stats: self.comp_based_stats.clone(),
            num_threads: self.num_threads,
            query_range,
        })
    }
}

/// LOSAT's limits, checked where NCBI starts the search: the options that NCBI runs and
/// LOSAT's TBLASTN does not implement are rejected, and the others select the search.
fn search_settings(args: &ResolvedTblastnArgs) -> Result<SearchSettings> {
    use crate::blastinput::app::CompositionMode;
    // NCBI's compressed lookup table (blast_aalookup.c, word size 5 and more or
    // tblastn-fast) is not implemented by LOSAT's TBLASTN.
    if args.compressed_lookup {
        if args.task == "tblastn-fast" {
            anyhow::bail!(
                "-task tblastn-fast (word size {} with the compressed lookup table) is not supported by LOSAT's TBLASTN",
                args.word_size
            );
        }
        anyhow::bail!(
            "-word_size {} (the compressed lookup table) is not supported by LOSAT's TBLASTN",
            args.word_size
        );
    }
    if !args.gapped {
        anyhow::bail!("-ungapped (an ungapped search) is not supported by LOSAT's TBLASTN");
    }
    // NCBI reference: c++/src/algo/blast/core/blast_parameters.c:774-815 (uneven gap
    // linking of the HSPs, `longest_intron`)
    if args.max_intron_length != 0 {
        anyhow::bail!(
            "-max_intron_length {} (uneven gap linking) is not supported by LOSAT's TBLASTN",
            args.max_intron_length
        );
    }
    // NCBI reference: c++/src/algo/blast/core/blast_traceback.c:1486-1499
    // ```c
    //     } else if (ext_params->options->compositionBasedStats > 0 ||
    //                ext_params->options->eTbackExt == eSmithWatermanTbck) {
    // ...
    //         retval =
    //                 Blast_RedoAlignmentCore_MT(program_number,
    // ...
    //     } else {
    // ```
    // The composition mode selects the ordinary traceback (mode 0) or the composition
    // redo; LOSAT's TBLASTN implements modes 0 and 2.
    let composition_mode2 = match args.composition_mode {
        CompositionMode::NoCompositionBasedStats => false,
        CompositionMode::CompositionMatrixAdjust => true,
        CompositionMode::CompositionBasedStats | CompositionMode::CompoForceFullMatrixAdjust => {
            anyhow::bail!(
                "-comp_based_stats {} (composition mode {}) is not supported by LOSAT's TBLASTN",
                args.comp_based_stats,
                args.composition_mode as i32
            )
        }
    };
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187:
    // arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
    // NCBI c++/src/algo/blast/api/prelim_stage.cpp:145-147:
    // TBlastThreads the_threads(GetNumberOfThreads());
    crate::utils::threading::validate_threads(args.num_threads)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:409,3687
    // ```c
    //     *pbestEvalue = DBL_MAX;
    //                 if (best_evalue <= hitParams->options->expect_value) {
    // ```
    // A list whose HSPs all fail keeps the best e-value DBL_MAX, so an -evalue of DBL_MAX
    // or more (+inf, 1e999, 1.7976931348623157e308) puts an empty list into the hit list,
    // and NCBI's tblastn dies of SIGSEGV reading its first HSP (blast_hits.c:3266) when
    // the search keeps enough HSPs, and runs otherwise; LOSAT cannot tell beforehand, so
    // it rejects the value (decision D12 of docs/evidence/losat_web_e2e/AUTHORITY.md). A
    // smaller -evalue such as 1e308 runs.
    if args.evalue >= f64::MAX {
        anyhow::bail!(
            "an -evalue of DBL_MAX or more ({}), with which NCBI BLAST+'s tblastn crashes on some inputs, is not supported by LOSAT's TBLASTN",
            args.evalue
        );
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:43-70 (the
    // preliminary hit list size through NCBI's `int`; `blastn/hsp.rs`)
    // A hit list size whose preliminary size is not positive crashes NCBI.
    if crate::algorithm::blastn::hsp::get_prelim_hitlist_size(
        args.hitlist_size,
        composition_mode2,
        true,
    ) < 1
    {
        anyhow::bail!(
            "a -max_target_seqs of {}, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's TBLASTN",
            args.hitlist_size
        );
    }
    let matrix: Option<ScoringMatrix> = args.matrix_name.parse().ok();
    let threshold = args.threshold;
    let window = args.window_size;
    let profile_supported = (matrix == Some(ScoringMatrix::Blosum62)
        && args.word_size == 3
        && args.gap_open == 11
        && args.gap_extend == 1
        && threshold == 13.0
        && window == 40)
        || (matrix == Some(ScoringMatrix::Blosum45)
            && !composition_mode2
            && args.word_size == 2
            && args.gap_open == 14
            && args.gap_extend == 2
            && threshold == 16.0
            && window == 60);
    let Some(matrix) = matrix.filter(|_| profile_supported) else {
        anyhow::bail!(
            "the matrix {} with gap costs {}/{}, word size {}, threshold {} and window size {} (composition mode {}) is not supported by LOSAT's TBLASTN; it implements BLOSUM62 11/1 with word size 3, threshold 13 and window size 40, and BLOSUM45 14/2 with word size 2, threshold 16, window size 60 and -comp_based_stats 0",
            args.matrix_name,
            args.gap_open,
            args.gap_extend,
            args.word_size,
            threshold,
            window,
            args.composition_mode as i32
        );
    };
    Ok(SearchSettings {
        composition_mode2,
        scoring: LocalStageDScoring {
            matrix,
            gap_open: args.gap_open,
            gap_extend: args.gap_extend,
            word_size: args.word_size as usize,
            threshold: threshold as i32,
            window,
            gap_xdrop_bits: args.xdrop_gap,
            final_xdrop_bits: args.xdrop_gap_final,
        },
    })
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options (`TblastnArgs::check_options`) and LOSAT's limits.
pub fn check_options(args: &TblastnArgs) -> Result<()> {
    search_settings(&args.resolve()?).map(|_| ())
}

impl TblastnArgs {
    // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
    // results = lcl_blast.Run();
    // formatter.PrintOneResultSet(**result, query);
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
    // PrintOneResultSet dispatches the complete result to outfmt 0/6/7.
    /// `settings` are the NCBI application settings that LOSAT reproduces
    /// (`ncbi_environment::check_ncbi_application_settings`): the subject reader
    /// (`fasta_reader`) uses them; the query keeps today's reading until step S8 of the port
    /// plan.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
    /// ```c++
    ///     CNcbiApplication* app = CNcbiApplication::Instance();
    ///     if (app) {
    ///         const CNcbiRegistry& registry = app->GetConfig();
    ///         x_LoadDataLoadersConfig(registry);
    ///         x_LoadBlastDbDataLoaderConfig(registry);
    ///     }
    /// ```
    pub fn run(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    ) -> Result<()> {
        use crate::algorithm::blastn::input as fasta_input;
        use crate::blastinput::app;
        use crate::blastinput::fasta_reader::{read_subjects, FastaInputSource, ReaderConfig};
        use crate::blastinput::seq_range;
        let args = self;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
        // ```c
        //     if (GetExportSearchStrategyStream(args) ||
        //            m_FormattingArgs->ArchiveFormatRequested(args)) {
        //         locality = CBlastOptions::eBoth;
        //     }
        // ```
        // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
        // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
        // there too.
        let choice = app::parse_formatting_string(&args.outfmt)?;
        let format = app::report_format(&choice, "TBLASTN", false, args.out.as_deref(), false)?;
        // LOSAT's thread capability is checked before any input or output, as for BLASTN.
        crate::utils::threading::validate_threads(args.num_threads)?;
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2562
        // ```c
        //         CRef<blast::CBlastQueryVector> subjects;
        //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
        //                                        subj_range, parse_deflines,
        //                                        use_lcase_masks, subjects, m_IsMapper);
        //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
        //
        //     } else if (!m_IsIgBlast){
        //         // IgBlast permits use of germline database
        //         NCBI_THROW(CInputException, eInvalidInput,
        //            "Either a BLAST database or subject sequence(s) must be specified");
        //     }
        // ```
        // The handler of the database arguments reads the subjects (an empty subject set
        // fails there), before the query and the output are opened.
        let Some(subject_path) = args.subject.clone() else {
            return Err(app::missing_subject_error());
        };
        fasta_input::check_utf8_file_name(&subject_path, "subject", "TBLASTN")?;
        let subject_file = fasta_input::open_input(&subject_path, "subject", "TBLASTN")?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
        // ```c++
        //             subj_input_stream = &args[kArgSubject].AsInputFile();
        //         }
        //
        //         TSeqRange subj_range;
        //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
        //             subj_range =
        //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
        //                             "Invalid specification of subject location");
        //         }
        // ```
        // The subject range is read after the subject file is opened, before it is read.
        let subject_range = seq_range::parse_optional_range(
            args.subject_loc.as_deref(),
            seq_range::RangeRole::Subject,
            "TBLASTN",
        )?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-246
        // ```c++
        //     SDataLoaderConfig dlconfig(read_proteins);
        //     dlconfig.OptimizeForWholeLargeSequenceRetrieval();
        //
        //     CBlastInputSourceConfig iconfig(dlconfig);
        //     iconfig.SetRange(range);
        //     iconfig.SetBelieveDeflines(parse_deflines);
        //     iconfig.SetLowercaseMask(use_lcase_masking);
        //     iconfig.SetSubjectLocalIdMode();
        //     if (!read_proteins && gaps_to_Ns) {
        //         iconfig.SetConvertGapsToNs(true);
        //     }
        //
        //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
        //     CRef<CBlastInput> input(new CBlastInput(fasta));
        //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
        //     sequences = input->GetAllSeqs(*scope);
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2466-2470
        // ```c++
        //     EMoleculeType mol_type = Blast_SubjectIsNucleotide(opts.GetProgramType())
        //         ? CSearchDatabase::eBlastDbIsNucleotide
        //         : CSearchDatabase::eBlastDbIsProtein;
        //     m_IsProtein = (mol_type == CSearchDatabase::eBlastDbIsProtein);
        // ```
        // TBLASTN's subjects are nucleotide (`read_proteins` = `IsProtein()` is false, so the
        // reader assumes nucleotides, `fAssumeNuc`). NCBI reads the records one at a time, writes
        // each one's messages as it reads it, and checks each one's range after reading it: a
        // range that starts past the end of a record stops the reading there
        // (`read_subjects`). Records without residues (and intervals without letters) are
        // read without a message and reported when the formatter sets up the subjects
        // (`before_prolog`).
        let mut subject_source = FastaInputSource::from_argument(
            &subject_path,
            subject_file,
            ReaderConfig::subject("TBLASTN", false, settings.data_loaders),
        );
        let read_subject_records = {
            let mut stderr = std::io::stderr();
            read_subjects(
                &mut subject_source,
                subject_range.as_ref(),
                &mut |message: &[u8]| std::io::Write::write_all(&mut stderr, message),
            )?
        };
        drop(subject_source);
        let (subjects, subject_placements) =
            match seq_range::cut_subjects(&read_subject_records, subject_range.as_ref())? {
                Some((cut, placements)) => (cut, placements),
                None => (read_subject_records, seq_range::Placements::default()),
            };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // ```c
        // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
        //     : m_QueryVector(& queries)
        // {
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // `m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects))` (blast_args.cpp:2557)
        // fails for a subject input without records.
        if subjects.is_empty() {
            return Err(app::empty_subjects_error());
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
        // ```c
        //     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
        //         m_InputStream == NULL) {
        // ...
        //         else {
        //             m_InputStream = &args[kArgQuery].AsInputFile();
        //         }
        //     }
        // ...
        //     else {
        //         m_OutputStream = &args[kArgOutput].AsOutputFile();
        //     }
        // ```
        // The query is opened and the output file created before the options are checked;
        // `-` is standard input and standard output.
        fasta_input::check_utf8_file_name(&args.query, "query", "TBLASTN")?;
        let query_file = fasta_input::open_input(&args.query, "query", "TBLASTN")?;
        if let Some(path) = args.out.as_deref() {
            fasta_input::check_utf8_file_name(path, "out", "TBLASTN")?;
        }
        let out_file = match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
            Some(path) => Some(std::io::BufWriter::new(
                std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
            )),
            None => None,
        };
        let pairwise = format == Some(app::ReportFormat::Pairwise);
        let outfmt = choice.normalized();
        let mut stderr = std::io::stderr();
        let mut stream = crate::cli::ReportStream {
            inner: match out_file {
                Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
                None => crate::cli::report_standard_output(),
            },
            failed: false,
        };
        let result = {
            let mut outputs =
                ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut stderr);
            search_cli(
                args,
                query_file,
                &subjects,
                &subject_placements,
                &mut outputs,
                &choice,
            )
        };
        let flushed = std::io::Write::flush(&mut stream);
        // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
        // ```c
        //     catch (const std::ios::failure&) {                                      \
        //         LOG_POST(Error << "BLAST failed to write output");                  \
        //         exit_code = BLAST_OUTPUT_ERROR;                                     \
        //     }                                                                       \
        // ```
        // The outfmt 0 formatter's stream throws when a write fails; for outfmt 6/7 NCBI
        // aborts and LOSAT reports the error (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 3).
        if stream.failed && pairwise {
            return Err(crate::cli::NativeError {
                exit: 6,
                message: "BLAST failed to write output\n".to_string(),
            }
            .into());
        }
        result?;
        flushed.context("failed to write the output")
    }
}

/// The part of `run` after the files are opened: NCBI's processing of the options,
/// `Query is Empty!`, the query batch size and the subjects without letters, LOSAT's checks
/// of the query and the search.
fn search_cli(
    args: TblastnArgs,
    mut query_file: std::fs::File,
    subjects: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    outputs: &mut ReportOutputs<'_>,
    choice: &crate::blastinput::app::FormatChoice,
) -> Result<()> {
    use crate::algorithm::blastn::input as fasta_input;
    let resolved = args.check_options(outputs.diagnostics)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:213-216
    // ```c
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // When the subjects were read from standard input too, `cin` has reached its end, so
    // NCBI gets no position for the query either.
    let seekable = !(args.query.as_os_str() == "-"
        && args
            .subject
            .as_deref()
            .is_some_and(|path| path.as_os_str() == "-"))
        && std::io::Seek::stream_position(&mut query_file).is_ok();
    let query_bytes = fasta_input::read_fasta_bytes(&mut query_file, &args.query, "query")?;
    drop(query_file);
    if fasta_input::is_blank(&query_bytes) {
        if !seekable {
            anyhow::bail!(
                "an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTN"
            );
        }
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    let batch_size = before_prolog(subjects, outputs.diagnostics)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:267-270
    // ```c
    // 	if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made.
    crate::blastinput::app::xinclude_check(choice)?;
    // The query keeps LOSAT's checks of its `bio` records until step S8 of the port plan;
    // the records checked are searched as NCBI's reader's records (`FastaRecord::from_bio`).
    fasta_input::check_protein_sequence_lines_of(&query_bytes, "query", "TBLASTN")?;
    let queries = fasta_input::bio_records_of(&query_bytes, &args.query, "query", "TBLASTN")?;
    fasta_input::check_protein_input_of(&query_bytes, &queries, "query", "TBLASTN")?;
    // NCBI warns about a query without residues ("Sequence contains no data") and fails
    // when no query has any; LOSAT rejects such a query, as BLASTN does.
    fasta_input::check_records_have_residues_of(&queries, "query", "TBLASTN")?;
    drop(query_bytes);
    let queries: Vec<FastaRecord> = queries
        .iter()
        .enumerate()
        .map(|(index, record)| FastaRecord::from_bio(record, index + 1, "Query_", true))
        .collect();
    search(
        &args,
        resolved,
        batch_size,
        &queries,
        subjects,
        subject_placements,
        outputs,
    )
}

/// NCBI's steps of tblastn between `Query is Empty!` and the outfmt 0 prolog that concern
/// the subjects: the query batch size, then the formatter, which sets up the subjects and
/// warns about each one without letters (a record without residues, or an interval that
/// starts just past its record's end). Such a subject stays in the database statistics
/// with no letters and is never searched (`scan_protein_words_by_chunk` skips subjects
/// shorter than a word). Returns the query batch size.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:217-221
/// ```c
///             fasta.Reset(new CBlastFastaInputSource(
///                                          m_CmdLineArgs->GetInputStream(),
///                                          iconfig));
///             input.Reset(new CBlastInput(&*fasta,
///                                         m_CmdLineArgs->GetQueryBatchSize()));
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
/// ```c
///     if(m_IsDbScan) {
/// 	int num_seqs=0;
///         int total_length=0;
/// 	if (!is_remote_search)
///         {
///                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
///                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
///                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
/// ```c
///         catch(CBlastException & e ) {
///         	// Skip bad subject sequence
///         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
///         		seqblk_vec->push_back(subj);
///         ...
///         		warning += "Subject sequence contains no data";
///         		ERR_POST(Warning << warning);
///         		continue;
/// ```
/// NCBI's gapped tblastn never stops for the average subject length, also when every
/// subject is without letters: its shortest subject is raised to `BLAST_SEQSRC_MINLENGTH`
/// before the division by 3 (`run_local_search_with_pool` ports blast_setup.c:969-973 and
/// seqsrc_multiseq.cpp:226-240; LOSAT rejects `-ungapped` and `OLD_FSC`, the two ways to a
/// search without `sbp->gbp`).
fn before_prolog(subjects: &[FastaRecord], diagnostics: &mut dyn Write) -> Result<u32> {
    let batch_size = crate::blastinput::app::query_batch_size("TBLASTN", 20000)?;
    crate::algorithm::blastn::input::write_empty_subject_warnings(
        subjects,
        "tblastn",
        diagnostics,
    )?;
    Ok(batch_size)
}

// NCBI reference: c++/src/algo/blast/api/tblastn_options.cpp:54-58
// ```c++
// CTBlastnOptionsHandle::SetLookupTableDefaults()
// {
//     CBlastProteinOptionsHandle::SetLookupTableDefaults();
//     SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
// }
// ```
// The options handle keeps the composition mode and the scoring choices that the
// search reads; these are the ones the Stage D port uses.
/// The search settings that the options select.
struct SearchSettings {
    composition_mode2: bool,
    scoring: LocalStageDScoring,
}

// NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
// results = lcl_blast.Run();
// formatter.PrintOneResultSet(**result, query);
// NCBI c++/src/app/blast/blast_formatter.cpp:429-467:
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one TBLASTN search over already parsed records and writes every requested
/// output format from the same result (the shared entry of the CLI and web ABI v2).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the
/// options and of the subjects come before the search. The `-query` and `-subject` values
/// of `args` are used only as display names.
///
/// The records are those of NCBI's reader (`fasta_reader`): the subjects as NCBI's
/// nucleotide reader reads them, the queries made from `bio` records with
/// `FastaRecord::from_bio` (LOSAT's checks of the query keep `bio`'s reading until step S8
/// of the port plan).
pub fn run_local(
    args: TblastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        if app::report_format(&choice, "TBLASTN", false, None, false)?.is_none() {
            app::formatting_handler_check(&choice, false)?;
            app::xinclude_check(&choice)?;
        }
    }
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        seq_range::RangeRole::Subject,
        "TBLASTN",
    )?;
    for record in
        &subject_records[..seq_range::subjects_read(subject_records, subject_range.as_ref())]
    {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = seq_range::cut_subjects(subject_records, subject_range.as_ref())?;
    if subject_records.is_empty() {
        return Err(app::empty_subjects_error());
    }
    let resolved = args.check_options(outputs.diagnostics)?;
    if query_records.is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    let whole_subjects = seq_range::Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    let batch_size = before_prolog(searched_subjects, outputs.diagnostics)?;
    search(
        &args,
        resolved,
        batch_size,
        query_records,
        searched_subjects,
        subject_placements,
        outputs,
    )
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options, the
/// reading of the subjects and the steps before the prolog (`before_prolog`): LOSAT's
/// checks of the queries and limits, the environment, and the search.
fn search(
    args: &TblastnArgs,
    resolved: ResolvedTblastnArgs,
    batch_size: u32,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::algorithm::blastn::input as fasta_input;
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    // The subjects are NCBI's reader's records: their residues are read as NCBI reads
    // them (`U` as `T`), and a subject without letters stays in the database (its warning
    // was written before the prolog, `before_prolog`). The queries keep LOSAT's checks.
    fasta_input::check_protein_residues_of(query_records, "query", "TBLASTN")?;
    fasta_input::check_records_have_residues_of(query_records, "query", "TBLASTN")?;
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`seq_range::cut_queries`): the records were checked whole, and are searched cut.
    let input_records = query_records;
    let ranged_queries = resolved
        .query_range
        .as_ref()
        .map(|range| seq_range::cut_queries(query_records, range));
    let (query_records, query_input) = match ranged_queries {
        Some(ranged) => (std::borrow::Cow::Owned(ranged.records), ranged.input),
        None => (
            std::borrow::Cow::Borrowed(query_records),
            seq_range::QueryInput::whole(query_records),
        ),
    };
    let query_records: &[FastaRecord] = &query_records;
    seq_range::check_no_empty_interval(
        query_records,
        &query_input.placements,
        &query_input.ordinals,
        seq_range::RangeRole::Query,
        "TBLASTN",
    )?;
    let SearchSettings {
        composition_mode2,
        scoring,
    } = search_settings(&resolved)?;
    app::check_unsupported_environment("TBLASTN")?;
    app::check_old_fsc("TBLASTN")?;
    app::check_query_split_environment("TBLASTN", false)??;
    let subjects = subject_records;
    let subject_path = args
        .subject
        .as_ref()
        .context("TBLASTN subject file missing")?;
    // The batches count whole records (`seq_range::QueryInput::batches`); with
    // `-query_loc`, a batch whose records are all skipped stops NCBI.
    let query_batches = query_input.batches(batch_size);
    if let Some(first_batch) = query_batches
        .first()
        .filter(|batch| batch.searched.is_empty())
    {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // ```c
        // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
        //     : m_QueryVector(& queries)
        // {
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // A batch size of 0 reads no query into the first batch, after the prolog; so does
        // a first batch whose records are all skipped (`-query_loc`), after their title
        // warnings.
        for format in outputs.formats.iter_mut() {
            if app::parse_formatting_string(format.outfmt)?.number == 0 {
                let mut writer = format.sink.open()?;
                crate::report::pairwise::write_tblastn_pairwise_prolog(
                    &mut writer,
                    "2.17.0+",
                    &format!(
                        "User specified sequence set (Input: {})",
                        subject_path.display()
                    ),
                    subjects.len(),
                    subjects.iter().map(|record| record.seq().len()).sum(),
                )?;
                writer.flush()?;
            }
        }
        for record in &input_records[first_batch.input.clone()] {
            outputs.diagnostics.write_all(&record.warnings)?;
        }
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    let last_batch_empty = query_batches
        .last()
        .is_some_and(|batch| batch.searched.is_empty());
    let query_offsets: Vec<usize> = (0..query_records.len())
        .map(|index| query_input.placements.offset(index))
        .collect();
    let queries = query_records;
    let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
    let subject_seqs: Vec<_> = subjects
        .iter()
        .map(|record| record.seq().to_vec())
        .collect();
    let seg = resolved.seg.params();
    // NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
    // (*thread)->Run(); (*thread)->Join(&result);
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1417:
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    //                         CConstRef<blast::CBlastQueryVector> queries,
    // NCBI c++/src/app/blast/tblastn_app.cpp:251-252,275-301:
    // input.Reset(new CBlastInput(&*fasta, GetQueryBatchSize()));
    // for (; !input->End(); ...) {
    //     query = input->GetNextSeqBatch(*scope);
    //     CLocalBlast lcl_blast(query_factory, ..., db_adapter);
    //     results = lcl_blast.Run();
    // }
    // NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:70-127:
    // case eTblastn: retval = 20000;
    // Search each accumulated batch independently; preliminary linking
    // and sum-statistics depend on its concatenated query contexts.
    let mut results = Vec::with_capacity(queries.len());
    let mut lengths: Option<LocalSubjectParameters> = None;
    let mut ungapped_karlin = Vec::with_capacity(queries.len());
    let mut query_validity = Vec::with_capacity(queries.len());
    let mut query_batch_skipped = Vec::with_capacity(queries.len());
    for range in query_batches
        .iter()
        .map(|batch| batch.searched.clone())
        .filter(|searched| !searched.is_empty())
    {
        let (mut batch_results, batch_lengths, batch_karlin, batch_validity) =
            run_local_for_report_threads(
                &query_seqs[range.clone()],
                &subject_seqs,
                LocalStageDProfile {
                    seg: seg.as_ref(),
                    soft_masking: resolved.soft_masking,
                    mask_lowercase: resolved.lcase_masking,
                    genetic_code: resolved.db_gencode,
                    expect_value: resolved.evalue,
                    max_target_seqs: resolved.hitlist_size,
                    query_offsets: &query_offsets[range],
                },
                composition_mode2,
                resolved.sum_stats,
                scoring,
                resolved.num_threads,
            )?;
        results.append(&mut batch_results);
        if let Some(all_lengths) = &mut lengths {
            // NCBI tblastn_app.cpp:288-301 formats each batch's own
            // query-context lengths; retain those values in input order.
            all_lengths.lengths.extend(batch_lengths.lengths);
            all_lengths.cutoffs.extend(batch_lengths.cutoffs);
        } else {
            lengths = Some(batch_lengths);
        }
        ungapped_karlin.extend(batch_karlin);
        // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
        // if (CheckInternalData() != 0) each result in this Run() batch
        // receives an unsearched ancillary state and null align set.
        let skipped = batch_validity.iter().all(|&valid| !valid);
        query_batch_skipped.extend(std::iter::repeat(skipped).take(batch_validity.len()));
        query_validity.extend(batch_validity);
    }
    let lengths = lengths.context("TBLASTN query file is empty")?;
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1443-1451:
    // if (results.HasWarnings()) ERR_POST(Warning << results.GetWarningStrings());
    // The warnings belong to the result: each query's are written before its report, once
    // (`QueryWarnings`), on one line in NCBI's order (`query_warning`). The reader's
    // warnings about the titles of a query batch come before the batch's first report.
    let mut query_warning_lines: Vec<Vec<u8>> = queries
        .iter()
        .zip(&query_validity)
        .enumerate()
        .map(|(index, (query, valid))| {
            let mut messages = Vec::new();
            messages.extend(crate::report::query_warnings::replaced_o_message(
                query.seq(),
            ));
            if !valid {
                messages.push(crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string());
            }
            crate::report::query_warnings::query_warning(
                "tblastn",
                query_input.ordinal(index),
                query,
                &messages,
            )
        })
        .collect();
    let unsearched_title_warnings =
        crate::report::query_warnings::prepend_ranged_batch_title_warnings(
            &mut query_warning_lines,
            input_records,
            &query_batches,
            |record: &FastaRecord| record.warnings.clone(),
        );
    render(
        outputs,
        queries,
        subjects,
        subject_path,
        &mut results,
        &lengths,
        &ungapped_karlin,
        &query_validity,
        &query_batch_skipped,
        scoring,
        &resolved.matrix_name,
        resolved.max_target_seqs_given,
        resolved.db_gencode,
        seg.as_ref(),
        resolved.lcase_masking,
        &query_warning_lines,
        &super::stage_e_report::ReportRanges {
            queries: &query_input.placements,
            subjects: subject_placements,
            epilog: !last_batch_empty,
        },
    )?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // ```c
    // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
    //     : m_QueryVector(& queries)
    // {
    //     if (queries.Empty()) {
    //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
    //     }
    // ```
    // With `-query_loc`, a last batch whose records are all skipped stops NCBI after the
    // reports of the batches before (written without the epilog), when it has read that
    // batch (its title warnings).
    if last_batch_empty {
        outputs.diagnostics.write_all(&unsearched_title_warnings)?;
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: Empty CBlastQueryVector\n".to_string(),
        }
        .into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::cli::{try_parse_from, Cli, Commands};

    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:1938-1942
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //     "Use lower case filtering in query and subject sequence(s)?", true);
    #[test]
    fn lcase_masking_flag_is_explicit() {
        assert!(!parse(&[]).unwrap().lcase_masking);
        assert!(parse(&["-lcase_masking"]).unwrap().lcase_masking);
    }

    fn parse(extra: &[&str]) -> std::result::Result<TblastnArgs, clap::Error> {
        let mut argv = vec!["losat", "tblastn", "-query", "q.faa", "-subject", "s.fna"];
        argv.extend_from_slice(extra);
        let cli: Cli = try_parse_from(argv)?;
        let Commands::Tblastn(args) = cli.command else {
            unreachable!()
        };
        Ok(args)
    }

    fn error(extra: &[&str]) -> String {
        parse(extra).unwrap().resolve().unwrap_err().to_string()
    }

    fn limit(extra: &[&str]) -> String {
        check_options(&parse(extra).unwrap())
            .unwrap_err()
            .to_string()
    }

    // NCBI tblastn_options.cpp:54-85 and blast_prot_options.cpp:85-141:
    // SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
    // m_Opts->SetSumStatisticsMode(); SetDbGeneticCode(BLAST_GENETIC_CODE);
    #[test]
    fn ordinary_defaults_and_code32_are_accepted() {
        let args = parse(&[]).unwrap();
        assert_eq!(args.outfmt, "0");
        assert_eq!(args.comp_based_stats, "2");
        let resolved = args.resolve().unwrap();
        assert_eq!(resolved.task, "tblastn");
        assert_eq!(resolved.db_gencode, 1);
        assert_eq!(resolved.word_size, 3);
        assert_eq!(resolved.threshold, 13.0);
        assert_eq!(resolved.window_size, 40);
        assert_eq!((resolved.gap_open, resolved.gap_extend), (11, 1));
        assert_eq!(resolved.evalue, 10.0);
        assert_eq!(resolved.matrix_name, "BLOSUM62");
        assert!(resolved.sum_stats && !resolved.soft_masking);
        assert_eq!(resolved.max_intron_length, 0);
        assert_eq!(resolved.hitlist_size, 500);
        check_options(&args).unwrap();
        let resolved = parse(&["-soft_masking", "true", "-sum_stats", "off"])
            .unwrap()
            .resolve()
            .unwrap();
        assert!(resolved.soft_masking && !resolved.sum_stats);
        assert!(parse(&["-sum_stats", "maybe"]).is_err());
        // NCBI blast_args.cpp:586-623: the matrix's suggestion (+2 for a translated
        // subject) replaces only the default 13.
        let resolved = parse(&["-matrix", "BLOSUM45"]).unwrap().resolve().unwrap();
        assert_eq!(
            (
                resolved.threshold,
                resolved.window_size,
                resolved.gap_open,
                resolved.gap_extend
            ),
            (16.0, 60, 14, 2)
        );
        for id in [
            1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            30, 31, 32, 33,
        ] {
            let id_string = id.to_string();
            assert_eq!(parse(&["-db_gencode", &id_string]).unwrap().db_gencode, id);
        }
    }

    // NCBI blast_args.cpp:997-1056 has 26 CLI values; gc.prt:340-347
    // defines ID 32 for the TBLASTN-only product decision.
    #[test]
    fn invalid_codes_are_explicit_errors() {
        for id in ["0", "7", "8", "17", "20", "34", "255", "256", "abc", "-1"] {
            let error = parse(&["-db_gencode", id]).unwrap_err().to_string();
            assert!(error.contains("genetic code"), "{id}: {error}");
        }
    }

    // NCBI's checks of the options (blast_options.c BLAST_ValidateOptions, the option
    // handlers of blast_args.cpp) and LOSAT's limits, where NCBI starts the search.
    #[test]
    fn ncbi_checks_and_losat_limits_are_explicit() {
        assert!(error(&["-matrix", "BAD"]).contains("BAD is not a supported matrix"));
        assert!(error(&["-gapopen", "1", "-gapextend", "1"])
            .contains("Gap existence and extension values of 1 and 1 not supported for BLOSUM62"));
        assert!(error(&["-matrix", "IDENTITY", "-word_size", "6"])
            .contains("Word size larger than 5 is not supported for the identity scoring matrix"));
        assert!(error(&["-evalue", "0"])
            .contains("expect value or cutoff score must be greater than zero"));
        assert!(error(&["-word_size", "8"]).contains("Word-size must be less than 8"));
        assert!(error(&["-threshold", "0"]).contains("Non-zero threshold required"));
        assert!(error(&["-ungapped"]).contains("Composition-adjusted searched"));
        assert!(error(&["-seg", "12 2.2"]).contains("Invalid number of arguments"));
        assert!(parse(&["-word_size", "1"]).is_err());
        // PAM30 with its best gap costs is an NCBI search that LOSAT does not implement.
        let resolved = parse(&["-matrix", "PAM30"]).unwrap().resolve().unwrap();
        assert_eq!((resolved.gap_open, resolved.gap_extend), (9, 1));
        for extra in [
            &["-matrix", "PAM30"][..],
            &["-word_size", "5"],
            &["-task", "tblastn-fast"],
            &["-ungapped", "-comp_based_stats", "F"],
            &["-max_intron_length", "100"],
            &["-comp_based_stats", "1"],
            &["-comp_based_stats", "3"],
            &["-window_size", "0"],
            &["-threshold", "13.5"],
            &["-max_target_seqs", "1073741799"],
            &["-comp_based_stats", "0", "-max_target_seqs", "2147483598"],
        ] {
            assert!(
                limit(extra).contains("not supported by LOSAT's TBLASTN"),
                "{extra:?}"
            );
        }
        let resolved = parse(&["-task", "tblastn-fast"])
            .unwrap()
            .resolve()
            .unwrap();
        assert_eq!((resolved.word_size, resolved.threshold), (5, 20.0));
        // NCBI blast_args.cpp:834-892: the first character chooses the mode (any other is
        // mode 0); a `u` asks for unified P-values only for blastp.
        for value in ["bogus", "0", "F", "2u", "D", "t"] {
            check_options(&parse(&["-comp_based_stats", value]).unwrap()).unwrap();
        }
        // NCBI blast_hits.c:43-70 through an int: 2^30 without composition statistics
        // keeps a preliminary hit list of 10 subjects.
        check_options(
            &parse(&["-comp_based_stats", "0", "-max_target_seqs", "1073741824"]).unwrap(),
        )
        .unwrap();
        for extra in [
            &["-max_hsps", "2"][..],
            &["-db", "db"],
            &["-remote"],
            &["-in_pssm", "p.asn"],
            &["-use_sw_tback"],
            &["-xdrop_ungap", "7"],
        ] {
            let error = crate::cli::render_message(&parse(extra).unwrap_err());
            assert!(
                error.contains("is not supported by LOSAT's TBLASTN"),
                "{extra:?}: {error}"
            );
        }
    }
}
