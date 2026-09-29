use crate::blastinput::value_parsers::*;
use clap::Args;
use std::path::PathBuf;

// CLI v2 names/defaults: NCBI c++/src/algo/blast/blastinput/blast_args.cpp:
// 166-170,203-207,332-349,2657-2660,3158-3163; cmdline_flags.cpp:46-94.
// arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value", ..., eInteger);
// arg_desc.SetConstraint(kArgMaxHSPsPerSubject, new CArgAllowValuesGreaterThanOrEqual(1));
// arg_desc.AddDefaultKey(kArgNumThreads, "int_value", ..., NStr::IntToString(kDfltValue));
// The single-dash lexical translation is owned by crate::cli.
#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
pub struct BlastnArgs {
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub query: PathBuf,
    #[arg(long, value_parser = file_path(), value_name = "PATH")]
    pub subject: PathBuf,
    #[arg(long, default_value = "megablast", long_help = "Implemented tasks: megablast and blastn. Task defaults: megablast uses word size 28, reward 1, penalty -2, gaps 0/0; blastn uses word size 11, reward 2, penalty -3, gaps 5/2. An omitted option takes the default of the task.", value_parser = ["megablast", "blastn"])]
    pub task: String,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:166-170,288-300
    // ```c
    //         arg_desc.AddOptionalKey(kArgWordSize, "int_value", description,
    //                                 CArgDescriptions::eInteger);
    // ...
    //     if ( args.Exist(kArgWordSize) && args[kArgWordSize]) {
    // ...
    //         opt.SetWordSize(args[kArgWordSize].AsInteger());
    // ```
    // An omitted value keeps the default of the task (`coordination.rs`).
    #[arg(long, value_parser = blastn_word_size, help = "Word size for wordfinder algorithm (default: 28 for megablast, 11 for blastn)")]
    pub word_size: Option<usize>,
    #[arg(long, default_value_t = 1, value_parser = positive_usize)]
    pub num_threads: usize,
    #[arg(long, default_value_t = 10.0, value_parser = blastn_evalue)]
    pub evalue: f64,
    /// Percent identity threshold for filtering HSPs (Blast_HSPTest).
    ///
    /// Default: 0.0 (disabled, matches NCBI BLAST default from calloc initialization).
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    /// ```c
    /// return ((hsp->num_ident * 100.0 <
    ///         align_length * hit_options->percent_identity) ||
    ///         align_length < hit_options->min_hit_length) ;
    /// ```
    #[arg(long = "perc_identity", default_value_t = 0.0, value_parser = percentage)]
    pub percent_identity: f64,
    /// Minimum hit length for filtering HSPs (Blast_HSPTest).
    ///
    /// Default: 0 (disabled, matches NCBI BLAST default from calloc initialization).
    ///
    /// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:993-1001 (s_HSPTest)
    /// ```c
    /// return ((hsp->num_ident * 100.0 <
    ///         align_length * hit_options->percent_identity) ||
    ///         align_length < hit_options->min_hit_length) ;
    /// ```
    #[arg(long, default_value_t = 0, help = "LOSAT-specific engine parameter")]
    pub min_hit_length: usize,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2960-2968
    // ```c
    // if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //    m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //    m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    //    hitlist_size = m_NumAlignments;
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:219,221
    // ```c
    // const size_t kDfltArgNumDescriptions = 500;
    // const size_t kDfltArgNumAlignments = 250;
    // ```
    // An omitted value keeps the default hit list size (500, `hitlist_size`) but the
    // pairwise report then shows 250 alignments, so the option has no clap default.
    #[arg(long, value_parser = positive_usize, help = "Maximum number of aligned sequences to keep (default: 500)")]
    pub max_target_seqs: Option<usize>,
    /// Maximum number of hits to save (NCBI BLAST hitlist_size)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    // NCBI blastinput/blast_args.cpp:2960-2968:
    // hitlist_size = m_NumAlignments; // resolved from max_target_seqs above
    // This internal fallback is not a separate public option.
    #[arg(skip = 500usize)]
    pub hitlist_size: usize,
    /// Remove word seeds with high frequency in the searched database.
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:257 (limit_lookup)
    #[arg(long = "limit_lookup", default_value_t = false)]
    pub limit_lookup: bool,
    /// Maximum database word count for lookup filtering.
    /// Reference: ncbi-blast/c++/include/algo/blast/core/blast_options.h:172-174
    /// #define MAX_DB_WORD_COUNT_MAPPER 30
    #[arg(
        long = "max_db_word_count",
        default_value_t = 30,
        help = "LOSAT-specific engine parameter"
    )]
    pub max_db_word_count: u8,
    /// Maximum number of HSPs per subject (unlimited when omitted)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    #[arg(long = "max_hsps", value_parser = positive_usize)]
    pub max_hsps_per_subject: Option<usize>,
    /// Minimum diagonal separation between HSPs on the same subject (0 = auto, task-specific)
    /// Reference: ncbi-blast/c++/src/algo/blast/api/blast_nucl_options.cpp:231-270
    /// Default: 50 for blastn, 6 for megablast
    // NCBI api/blast_nucl_options.cpp:239,259:
    // SetMinDiagSeparation(50); SetMinDiagSeparation(6);
    // Task resolution owns this internal value; no CLI override.
    #[arg(skip = 0usize)]
    pub min_diag_separation: usize,
    #[arg(long, value_name = "PATH")]
    pub out: Option<PathBuf>,
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:648-660,674-679
    // ```c
    //     arg_desc.AddOptionalKey(kArgMismatch, "penalty",
    //                            "Penalty for a nucleotide mismatch",
    //                            CArgDescriptions::eInteger);
    //     arg_desc.SetConstraint(kArgMismatch,
    //                            new CArgAllowValuesLessThanOrEqual(0));
    // ...
    //     arg_desc.SetConstraint(kArgMatch,
    //                            new CArgAllowValuesGreaterThanOrEqual(0));
    // ...
    //     if (cmd_line_args.Exist(kArgMismatch) && cmd_line_args[kArgMismatch]) {
    //         options.SetMismatchPenalty(cmd_line_args[kArgMismatch].AsInteger());
    //     }
    //     if (cmd_line_args.Exist(kArgMatch) && cmd_line_args[kArgMatch]) {
    //         options.SetMatchReward(cmd_line_args[kArgMatch].AsInteger());
    //     }
    // ```
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:175-182,262-271
    // ```c
    //         arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
    //                                 "Cost to open a gap",
    //                                 CArgDescriptions::eInteger);
    // ...
    //     if (args.Exist(kArgGapOpen) && args[kArgGapOpen]) {
    //         opt.SetGapOpeningCost(args[kArgGapOpen].AsInteger());
    //     }
    // ```
    // An omitted option keeps the default of the task (`coordination.rs`); a given one is
    // checked as NCBI checks it (`scoring.rs`). NCBI keeps the reward and the penalty in 16
    // bits; a reward of 0 or less (NCBI's rmblastn matrix scoring, or no valid query) is
    // not implemented and is rejected after NCBI's checks (`scoring.rs`).
    #[arg(long, value_parser = blastn_reward, help = "Reward for a nucleotide match (default: 1 for megablast, 2 for blastn; LOSAT does not support 0)")]
    pub reward: Option<i32>,
    #[arg(long, value_parser = nonpositive_i32, help = "Penalty for a nucleotide mismatch (default: -2 for megablast, -3 for blastn)")]
    pub penalty: Option<i32>,
    #[arg(
        long = "gapopen",
        help = "Cost to open a gap (default: 0 for megablast, 5 for blastn)"
    )]
    pub gap_open: Option<i32>,
    #[arg(
        long = "gapextend",
        help = "Cost to extend a gap (default: 0 for megablast, 2 for blastn)"
    )]
    pub gap_extend: Option<i32>,
    // NCBI blast_args.cpp:410-420: opt.SetDustFiltering(false/true);
    // opt.SetDustFilteringLevel(...); opt.SetDustFilteringWindow(...);
    // opt.SetDustFilteringLinker(...);
    #[arg(long, default_value = "20 64 1", value_parser = parse_dust_filtering, help = "DUST: no, yes, or LEVEL WINDOW LINKER")]
    pub dust: DustSpec,
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:1939-1942
    // ```c
    // arg_desc.AddFlag(kArgUseLCaseMasking,
    //      "Use lower case filtering in query and subject sequence(s)?", true);
    // ```
    #[arg(long = "lcase_masking", default_value_t = false)]
    pub lcase_masking: bool,
    /// Apply subject best hit filtering (disabled by default).
    /// Reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:135
    #[arg(long = "subject_besthit", default_value_t = false)]
    pub subject_besthit: bool,
    #[arg(
        long,
        default_value_t = false,
        help = "LOSAT-specific engine parameter"
    )]
    pub verbose: bool,
    /// Scan stride for subject sequence scanning (NCBI BLAST optimization).
    /// Higher values skip more positions, reducing k-mer lookups but potentially missing some seeds.
    /// Default: 0 (auto-calculate based on word_size: 1 for word_size < 16, 4 for word_size >= 16).
    /// For megablast (word_size=28), scan_step=4 reduces lookups by ~4x with minimal sensitivity loss.
    // NCBI core/blast_nalookup.c:397,1334:
    // lookup->scan_step = lookup->word_length - lookup->lut_word_length + 1;
    // Lookup construction owns this internal value; no CLI override.
    #[arg(skip = 0usize)]
    pub scan_step: usize,

    /// Output format (NCBI BLAST compatible).
    ///
    /// Supported formats:
    ///   0 = Pairwise (the NCBI report; the default)
    ///   6 = Tabular (tab-separated values)
    ///   7 = Tabular with comment lines (headers)
    ///
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782
    // ```c
    // case 6: case 7:
    //     formatter = new CBlastTabularFormatter(...);
    //     break;
    // ```
    // NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
    // ```c
    // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //     x_PrintField(*iter);
    // }
    // ```
    /// Custom field specifications are not yet ported and fail explicitly.
    ///
    /// Default fields: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
    #[arg(long, default_value = "0", value_name = "SPEC", value_parser = blastn_outfmt)]
    pub outfmt: String,
}
