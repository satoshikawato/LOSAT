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
    #[arg(long, default_value = "megablast", long_help = "Implemented tasks: megablast and blastn. Task-dependent engine defaults: megablast uses word size 28, reward 1, penalty -2, gaps 0/0; blastn uses word size 11, reward 2, penalty -3, gaps 5/2. The existing engine resolves sentinel/default-valued scoring fields by task.", value_parser = ["megablast", "blastn"])]
    pub task: String,
    #[arg(long, default_value_t = 28, value_parser = blastn_word_size)]
    pub word_size: usize,
    #[arg(long, default_value_t = 1, value_parser = positive_usize)]
    pub num_threads: usize,
    #[arg(long, default_value_t = 10.0, value_parser = nonnegative_f64)]
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
    #[arg(long, default_value = "500", value_parser = positive_usize)]
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
    // Scoring parameters - defaults are for megablast task
    // For blastn task, these are overridden in run() based on --task
    #[arg(long, default_value_t = 1, value_parser = positive_i32)]
    pub reward: i32,
    #[arg(long, default_value_t = -2, value_parser = negative_i32)]
    pub penalty: i32,
    #[arg(long = "gapopen", default_value_t = 0, value_parser = nonnegative_i32)]
    pub gap_open: i32,
    #[arg(long = "gapextend", default_value_t = 0, value_parser = nonnegative_i32)]
    pub gap_extend: i32,
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
    ///   0 = Pairwise (not yet implemented; fails explicitly)
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
