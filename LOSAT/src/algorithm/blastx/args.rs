//! BLASTX options and native serial local-subject search/report dispatch.
use crate::blastinput::value_parsers::{genetic_code, nonnegative_f64};
use anyhow::{bail, Result};
use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:49-62
// ```c++
//                                   "Translated Query-Protein Subject BLAST"));
//     const bool kQueryIsProtein = false;
//     m_Args.push_back(arg);
//     m_ClientId = kProgram + " " + CBlastVersion().Print();
//
//     static const char kDefaultTask[] = "blastx";
//     SetTask(kDefaultTask);
//     set<string> tasks;
//     tasks.insert(kDefaultTask);
//     tasks.insert("blastx-fast");
//     arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
//     m_Args.push_back(arg);
//
//     m_BlastDbArgs.Reset(new CBlastDatabaseArgs);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:72-78
// ```c++
//     // N.B.: query is not protein because the options are applied on the
//     // translated query
//     arg.Reset(new CGenericSearchArgs( !kQueryIsProtein ));
//     m_Args.push_back(arg);
//
//     //Disable until OOF is supported in align manager
//     //SB-1043
// ```
pub struct BlastxArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3425-3427
    // ```c++
    //     arg_desc.AddDefaultKey(kArgQuery, "input_file",
    //                      "Input file name",
    //                      CArgDescriptions::eInputFile, kDfltArgQuery);
    // ```
    #[arg(long, value_name = "PATH")]
    pub query: PathBuf,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2364-2366
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSubject, "subject_input_file",
    //                                 "Subject sequence(s) to search",
    //                                 CArgDescriptions::eInputFile);
    // ```
    #[arg(long, value_name = "PATH")]
    pub subject: PathBuf,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3441-3443
    // ```c++
    //     arg_desc.AddDefaultKey(kArgOutput, "output_file",
    //                    "Output file name",
    //                    CArgDescriptions::eOutputFile, "-");
    // ```
    #[arg(long, value_name = "PATH")]
    pub out: Option<PathBuf>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:99-100
    // ```c++
    //         arg_desc.AddDefaultKey(kTask, "task_name", "Task to execute",
    //                                CArgDescriptions::eString, m_DefaultTask);
    // ```
    #[arg(long, default_value = "blastx", value_parser = ["blastx"])]
    pub task: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:146-146
    // ```c++
    //         arg_desc.AddOptionalKey(kArgEvalue, "evalue", des, CArgDescriptions::eDouble);
    // ```
    #[arg(long, default_value_t = 10.0, value_parser = finite_f64)]
    pub evalue: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:166-167
    // ```c++
    //         arg_desc.AddOptionalKey(kArgWordSize, "int_value", description,
    //                                 CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = word_size)]
    pub word_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:203-205
    // ```c++
    //         arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value",
    //                            "Set maximum number of HSPs per subject sequence to save for each query",
    //                            CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = positive)]
    pub max_hsps: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:241-243
    // ```c++
    //         arg_desc.AddOptionalKey(kArgSumStats, "bool_value",
    //                      	 	"Use sum statistics",
    //                      	 	CArgDescriptions::eBoolean);
    // ```
    #[arg(long, value_parser = boolean, action = clap::ArgAction::Set)]
    pub sum_stats: Option<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:175-177
    // ```c++
    //         arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
    //                                 "Cost to open a gap",
    //                                 CArgDescriptions::eInteger);
    // ```
    #[arg(long = "gapopen")]
    pub gap_open: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:180-182
    // ```c++
    //         arg_desc.AddOptionalKey(kArgGapExtend, "extend_penalty",
    //                                "Cost to extend a gap",
    //                                CArgDescriptions::eInteger);
    // ```
    #[arg(long = "gapextend")]
    pub gap_extend: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:937-942
    // ```c++
    //     arg_desc.AddDefaultKey(kArgMaxIntronLength, "length",
    //                     "Length of the largest intron allowed in a translated "
    //                     "nucleotide sequence when linking multiple distinct "
    //                     "alignments",
    //                     CArgDescriptions::eInteger,
    //                     NStr::IntToString(kDfltArgMaxIntronLength));
    // ```
    #[arg(long, default_value_t = 0, value_parser = nonnegative)]
    pub max_intron_length: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:332-338
    // ```c++
    //         arg_desc.AddDefaultKey(kArgSegFiltering, "SEG_options",
    //                         "Filter query sequence with SEG "
    //                         "(Format: '" + kDfltArgApplyFiltering + "', " +
    //                         "'window locut hicut', or '" + kDfltArgNoFiltering +
    //                         "' to disable)",
    //                         CArgDescriptions::eString, m_FilterByDefault
    //                         ? kDfltArgSegFiltering : kDfltArgNoFiltering);
    // ```
    #[arg(long, default_value = "12 2.2 2.5", value_parser = seg)]
    pub seg: BlastxSeg,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:339-342
    // ```c++
    //         arg_desc.AddDefaultKey(kArgLookupTableMaskingOnly, "soft_masking",
    //                         "Apply filtering locations as soft masks",
    //                         CArgDescriptions::eBoolean,
    //                         kDfltArgLookupTableMaskingOnlyProt);
    // ```
    #[arg(long, default_value = "false", value_parser = boolean, action = clap::ArgAction::Set)]
    pub soft_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:531-534
    // ```c++
    //     arg_desc.AddDefaultKey(kArgMatrixName, "matrix_name",
    //                            "Scoring matrix name",
    //                            CArgDescriptions::eString,
    //                            string(""));
    // ```
    #[arg(long)]
    pub matrix: Option<String>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:578-581
    // ```c++
    //     arg_desc.AddOptionalKey(kArgWordScoreThreshold, "float_value",
    //                  "Minimum word score such that the word is added to the "
    //                  "BLAST lookup table",
    //                  CArgDescriptions::eDouble);
    // ```
    #[arg(long, value_parser = nonnegative_f64)]
    pub threshold: Option<f64>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:475-478
    // ```c++
    //     arg_desc.AddOptionalKey(kArgWindowSize, "int_value",
    //                             "Multiple hits window size, use 0 to specify "
    //                             "1-hit algorithm",
    //                             CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = nonnegative)]
    pub window_size: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3296-3299
    // ```c++
    //     arg_desc.AddOptionalKey(kArgCullingLimit, "int_value",
    //                      "If the query range of a hit is enveloped by that of at "
    //                      "least this many higher-scoring hits, delete the hit",
    //                      CArgDescriptions::eInteger);
    // ```
    #[arg(long, default_value_t = 0, value_parser = nonnegative)]
    pub culling_limit: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3329-3329
    // ```c++
    //     arg_desc.AddFlag(kArgSubjectBestHit, "Return only the best HSP for each non overlapping query region", true);
    // ```
    #[arg(long)]
    pub subject_besthit: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-917
    // ```c++
    //     arg_desc.AddFlag(kArgUngapped, "Perform ungapped alignment only?", true);
    // ```
    #[arg(long)]
    pub ungapped: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1941-1942
    // ```c++
    //     arg_desc.AddFlag(kArgUseLCaseMasking,
    //          "Use lower case filtering in query and subject sequence(s)?", true);
    // ```
    #[arg(long)]
    pub lcase_masking: bool,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1953-1955
    // ```c++
    //         arg_desc.AddDefaultKey(kArgStrand, "strand",
    //                      "Query strand(s) to search against database/subject",
    //                                CArgDescriptions::eString, kDfltArgStrand);
    // ```
    #[arg(long, default_value = "both", value_parser = ["both", "plus", "minus"])]
    pub strand: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1021-1024
    // ```c++
    //         arg_desc.AddDefaultKey(kArgQueryGeneticCode, "int_value",
    //                                "Genetic code to use to translate query (see https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes for details)\n",
    //                                CArgDescriptions::eInteger,
    //                                NStr::IntToString(BLAST_GENETIC_CODE));
    // ```
    #[arg(long, default_value_t = 1, value_parser = genetic_code)]
    pub query_gencode: u8,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
    // ```c++
    //     arg_desc.AddDefaultKey(kArgOutputFormat, "format",
    //                            kOutputFormatDescription,
    //                            CArgDescriptions::eString,
    //                            NStr::IntToString(dft_outfmt));
    // ```
    #[arg(long, default_value = "0", value_parser = outfmt)]
    pub outfmt: String,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2726-2730
    // ```c++
    //         arg_desc.AddOptionalKey(kArgMaxTargetSequences, "num_sequences",
    //                             "Maximum number of aligned sequences to keep \n"
    //     						"(value of 5 or more is recommended)\n"
    //     						"Default = `" + NStr::IntToString(BLAST_HITLIST_SIZE) + "'",
    //                             CArgDescriptions::eInteger);
    // ```
    #[arg(long, value_parser = positive)]
    pub max_target_seqs: Option<i32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3158-3161
    // ```c++
    //     arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
    //                            "Number of threads (CPUs) to use in the BLAST search",
    //                            CArgDescriptions::eInteger,
    //                            NStr::IntToString(kDfltValue));
    // ```
    #[arg(long, default_value_t = 1, value_parser = positive, help = "Number of threads (CPUs); serial WASI accepts only 1")]
    pub num_threads: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:796-797
    // ```c++
    //     arg_desc.AddDefaultKey(kArgCompBasedStats, "compo", legend,
    //                            CArgDescriptions::eString, m_DefaultOpt);
    // ```
    #[arg(long, default_value = "2", value_parser = composition)]
    pub comp_based_stats: u8,
}

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-407
// ```c++
//         if (m_QueryIsProtein && args[kArgSegFiltering]) {
//             const string& seg_opts = args[kArgSegFiltering].AsString();
//             if (seg_opts == kDfltArgNoFiltering) {
//                 opt.SetSegFiltering(false);
//             } else if (seg_opts == kDfltArgApplyFiltering) {
//                 opt.SetSegFiltering(true);
//             } else {
//                 x_TokenizeFilteringArgs(seg_opts, tokens);
//                 opt.SetSegFilteringWindow(NStr::StringToInt(tokens[0]));
//                 opt.SetSegFilteringLocut(NStr::StringToDouble(tokens[1]));
//                 opt.SetSegFilteringHicut(NStr::StringToDouble(tokens[2]));
//             }
// ```
pub struct BlastxSeg {
    pub enabled: bool,
    pub window: i32,
    pub locut: f64,
    pub hicut: f64,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:375-384
// ```c++
// CFilteringArgs::x_TokenizeFilteringArgs(const string& filtering_args,
//                                         vector<string>& output) const
// {
//     output.clear();
//     NStr::Split(filtering_args, " ", output);
//     if (output.size() != 3) {
//         NCBI_THROW(CInputException, eInvalidInput,
//                    "Invalid number of arguments to filtering option");
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-407
// ```c++
//         if (m_QueryIsProtein && args[kArgSegFiltering]) {
//             const string& seg_opts = args[kArgSegFiltering].AsString();
//             if (seg_opts == kDfltArgNoFiltering) {
//                 opt.SetSegFiltering(false);
//             } else if (seg_opts == kDfltArgApplyFiltering) {
//                 opt.SetSegFiltering(true);
//             } else {
//                 x_TokenizeFilteringArgs(seg_opts, tokens);
//                 opt.SetSegFilteringWindow(NStr::StringToInt(tokens[0]));
//                 opt.SetSegFilteringLocut(NStr::StringToDouble(tokens[1]));
//                 opt.SetSegFilteringHicut(NStr::StringToDouble(tokens[2]));
//             }
// ```
fn seg(value: &str) -> std::result::Result<BlastxSeg, String> {
    let mut s = BlastxSeg {
        enabled: value != "no",
        window: 12,
        locut: 2.2,
        hicut: 2.5,
    };
    if matches!(value, "yes" | "no") {
        return Ok(s);
    }
    let t: Vec<_> = value.split(' ').collect();
    if t.len() != 3 {
        return Err("invalid number of arguments to filtering option".into());
    }
    s.window = t[0].parse::<i32>().map_err(|_| "invalid SEG window")?;
    s.locut = t[1].parse::<f64>().map_err(|_| "invalid SEG locut")?;
    s.hicut = t[2].parse::<f64>().map_err(|_| "invalid SEG hicut")?;
    if !s.locut.is_finite() || !s.hicut.is_finite() {
        return Err("SEG requires finite cutoffs".into());
    }
    Ok(s)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:146-154
// ```c++
//         arg_desc.AddOptionalKey(kArgEvalue, "evalue", des, CArgDescriptions::eDouble);
//     } else if (m_QueryIsProtein) {
//         arg_desc.AddDefaultKey(kArgEvalue, "evalue",
//                      "Expectation value (E) threshold for saving hits ",
//                      CArgDescriptions::eDouble,
//                      NStr::DoubleToString(1.0));
//     } else {
//         //igblastn
//         arg_desc.AddDefaultKey(kArgEvalue, "evalue",
// ```
fn finite_f64(value: &str) -> std::result::Result<f64, String> {
    let x = value.parse::<f64>().map_err(|_| "requires a real value")?;
    if !x.is_finite() {
        Err("requires a finite real value".into())
    } else {
        Ok(x)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:937-945
// ```c++
//     arg_desc.AddDefaultKey(kArgMaxIntronLength, "length",
//                     "Length of the largest intron allowed in a translated "
//                     "nucleotide sequence when linking multiple distinct "
//                     "alignments",
//                     CArgDescriptions::eInteger,
//                     NStr::IntToString(kDfltArgMaxIntronLength));
//     arg_desc.SetConstraint(kArgMaxIntronLength,
//                            new CArgAllowValuesGreaterThanOrEqual(0));
//     arg_desc.SetCurrentGroup("");
// ```
fn nonnegative(value: &str) -> std::result::Result<i32, String> {
    let n = value
        .parse::<i32>()
        .map_err(|_| "requires an Int4 integer")?;
    if n < 0 {
        Err("requires an integer >= 0".into())
    } else {
        Ok(n)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2750-2754
// ```c++
// }
//
//
// static void s_ValidateCustomDelim(const string& customFmtSpec,const string& customDelim)
// {
// ```
fn positive(value: &str) -> std::result::Result<i32, String> {
    let n = nonnegative(value)?;
    if n == 0 {
        Err("requires an integer >= 1".into())
    } else {
        Ok(n)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:291-305
// ```c++
//            opt.SetWordThreshold(19.3);
//            if (args[kArgWordSize].AsInteger() > 5) {
//                opt.SetWordThreshold(21.0);
//            }
//            if (args[kArgWordSize].AsInteger() > 6) {
//                opt.SetWordThreshold(20.25);
//            }
//         }
//         opt.SetWordSize(args[kArgWordSize].AsInteger());
//
//     }
//
//     if (args.Exist(kArgEffSearchSpace) && args[kArgEffSearchSpace]) {
//         CNcbiEnvironment env;
//         env.Set("OLD_FSC", "true");
// ```
fn word_size(value: &str) -> std::result::Result<i32, String> {
    let n = nonnegative(value)?;
    if matches!(n, 3 | 5) {
        Ok(n)
    } else {
        Err("unsupported BLASTX word_size: use 3 or 5".into())
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:838-871
// ```c++
//         switch (comp_stat_string[0]) {
//             case '0': case 'F': case 'f':
//                 compo_mode = eNoCompositionBasedStats;
//                 break;
//             case '1':
//                 compo_mode = eCompositionBasedStats;
//                 break;
//             case 'D': case 'd':
//                 if ((program == eRPSBlast) || (program == eRPSTblastn)) {
//                     compo_mode = eNoCompositionBasedStats;
//                 }
//                 else if (program == eDeltaBlast) {
//                     compo_mode = eCompositionBasedStats;
//                 }
//                 else {
//                     compo_mode = eCompositionMatrixAdjust;
//                 }
//                 break;
//             case '2':
//                 compo_mode = eCompositionMatrixAdjust;
//                 break;
//             case '3':
//                 compo_mode = eCompoForceFullMatrixAdjust;
//                 break;
//             case 'T': case 't':
//                 compo_mode = (program == eRPSBlast || program == eRPSTblastn || program == eDeltaBlast) ?
//                     eCompositionBasedStats : eCompositionMatrixAdjust;
//                 break;
//         }
//
//         if(program == ePSITblastn) {
//             compo_mode = eNoCompositionBasedStats;
//         }
//
// ```
fn composition(value: &str) -> std::result::Result<u8, String> {
    match value {
        "0" | "F" | "f" => Ok(0),
        "2" | "D" | "d" | "T" | "t" => Ok(2),
        _ => Err("unsupported BLASTX composition token: use 0/F/f or 2/D/d/T/t".into()),
    }
}
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:487-498
// ```c++
//
//
// inline CArg_Boolean::CArg_Boolean(const string& name, const string& value)
//     : CArg_String(name, value)
// {
//     try {
//         m_Boolean = NStr::StringToBool(value);
//     } catch (const CException& e) {
//         NCBI_RETHROW(e,CArgException,eConvert, s_ArgExptMsg(GetName(),
//             "Argument cannot be converted",value));
//     }
// }
// ```
fn boolean(value: &str) -> std::result::Result<bool, String> {
    match value.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" | "1" => Ok(true),
        "false" | "f" | "no" | "n" | "0" => Ok(false),
        _ => Err("invalid boolean".into()),
    }
}
// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-41
// ```c++
// const char* kDfltArgTabularOutputFmt =
//     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
//     "evalue bitscore";
// const char* kDfltArgTabularOutputFmtTag("std");
// ```
pub const STD_FIELDS: &str =
    "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore";
// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:43-49
// ```c++
// const size_t kNumTabularOutputFormatSpecifiers = 50;
// const SFormatSpec sc_FormatSpecifiers[kNumTabularOutputFormatSpecifiers] = {
//     SFormatSpec("qseqid",
//                 "Query Seq-id",
//                 eQuerySeqId),
//     SFormatSpec("qgi",
//                 "Query GI",
// ```
pub const FIELDS: &[&str] = &[
    "qseqid", "qacc", "qaccver", "qlen", "sseqid", "sacc", "saccver", "slen", "qstart", "qend",
    "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident", "nident",
    "mismatch", "positive", "gapopen", "gaps", "ppos", "qframe", "sframe", "frames", "btop",
    "stitle",
];
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2827-2845
// ```c++
//         }
//         int val = 0;
//         try { val = NStr::StringToInt(fmt_choice); }
//         catch (const CStringException&) {   // probably a conversion error
//             CNcbiOstrstream os;
//             os << "'" << fmt_choice << "' is not a valid output format";
//             string msg = CNcbiOstrstreamToString(os);
//             NCBI_THROW(CInputException, eInvalidInput, msg);
//         }
//         if (val < 0 || val >= static_cast<int>(eEndValue)) {
//             string msg("Formatting choice is out of range");
//             throw std::out_of_range(msg);
//         }
//         if (m_IsIgBlast && (val != 3 && val != 4 && val != 7 && val != eAirrRearrangement)) {
//             string msg("Formatting choice is not valid");
//             throw std::out_of_range(msg);
//         }
//         fmt_type = static_cast<EOutputFormat>(val);
//         if ( !(fmt_type == eTabular ||
// ```
fn outfmt(value: &str) -> std::result::Result<String, String> {
    let mut tokens = value.split_whitespace();
    let f = tokens.next().ok_or("empty outfmt")?;
    if !matches!(f, "0" | "6" | "7") {
        return Err("unsupported BLASTX outfmt: use 0, 6 or 7".into());
    }
    let fields: Vec<_> = tokens.collect();
    if f == "0" && !fields.is_empty() {
        return Err("custom fields require outfmt 6/7".into());
    }
    for field in fields {
        if field != "std" && !FIELDS.contains(&field) {
            return Err(format!("unsupported BLASTX outfmt field '{field}'"));
        }
    }
    Ok(value.to_string())
}

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blastx_options.cpp:57-70
// ```c++
//     CBlastProteinOptionsHandle::SetLookupTableDefaults();
//     m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTX);
// }
//
// void
// CBlastxOptionsHandle::SetQueryOptionDefaults()
// {
//     CBlastProteinOptionsHandle::SetQueryOptionDefaults();
//     m_Opts->SetStrandOption(objects::eNa_strand_both);
//     m_Opts->SetQueryGeneticCode(BLAST_GENETIC_CODE);
//     SetSegFiltering(false); // disable SEG filtering because of eCompositionMatrixAdjust mode
// }
//
// void
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:131-146
// ```c++
// }
//
// void
// CBlastProteinOptionsHandle::SetHitSavingOptionsDefaults()
// {
//     SetHitlistSize(500);
//     SetEvalueThreshold(BLAST_EXPECT_VALUE);
//     SetMinDiagSeparation(0);
//     SetPercentIdentity(0);
//     // set some default here, allow INT4MAX to mean infinity
//     SetMaxNumHspPerSequence(0);
//     SetMaxHspsPerSubject(0);
//
//     SetCutoffScore(0); // will be calculated based on evalue threshold,
//     // effective lengths and Karlin-Altschul params in BLAST_Cutoffs_simple
//     // and passed to the engine in the params structure
// ```
pub struct ResolvedOptions {
    pub word_size: i32,
    pub threshold: f64,
    pub lookup_type: i32,
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gapped: bool,
    pub composition: u8,
    pub seg: BlastxSeg,
    pub soft_masking: bool,
    pub strand: String,
    pub query_gencode: u8,
    pub window_size: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // ```c++
    //
    // void
    // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    // }
    //
    // void
    // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
    //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:123-145
    // ```c++
    //                                              searches except blastn */
    // #define BLAST_UNGAPPED_X_DROPOFF_NUCL 20 /**< ungapped dropoff score for
    //                                               blastn (and megablast) */
    //
    // /** default dropoff for preliminary gapped extensions */
    // #define BLAST_GAP_X_DROPOFF_PROT 15 /**< default dropoff (all protein-
    //                                          based gapped extensions) */
    // #define BLAST_GAP_X_DROPOFF_NUCL 30 /**< default dropoff for non-greedy
    //                                          nucleotide gapped extensions */
    // #define BLAST_GAP_X_DROPOFF_GREEDY 25 /**< default dropoff for greedy
    //                                          nucleotide gapped extensions */
    // #define BLAST_GAP_X_DROPOFF_TBLASTX 0 /**< default dropoff for tblastx */
    //
    // /** default bit score that will trigger gapped extension */
    // #define BLAST_GAP_TRIGGER_PROT 22.0 /**< default bit score that will trigger
    //                                          a gapped extension for all protein-
    //                                          based searches */
    // #define BLAST_GAP_TRIGGER_NUCL 27.0  /**< default bit score that will trigger
    //                                          a gapped extension for blastn */
    //
    // /** default dropoff for the final gapped extension with traceback */
    // #define BLAST_GAP_X_DROPOFF_FINAL_PROT 25 /**< default dropoff (all protein-
    //                                                based gapped extensions) */
    // ```
    pub x_dropoff: f64,
    pub gap_x_dropoff: f64,
    pub gap_x_dropoff_final: f64,
    pub gap_trigger: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
    // ```c++
    // void
    // CBlastProteinOptionsHandle::SetEffectiveLengthsOptionsDefaults()
    // {
    //     SetDbLength(0);
    //     SetDbSeqNum(0);
    //     SetEffectiveSearchSpace(0);
    // }
    // ```
    pub db_length: i64,
    pub effective_search_space: i64,

    pub evalue: f64,
    pub sum_stats: bool,
    pub max_intron_length: i32,
    pub max_hsps: i32,
    pub culling_limit: i32,
    pub subject_besthit: bool,
    pub hitlist_size: i32,
    pub max_target_seqs_explicit: bool,
    pub num_descriptions: i32,
    pub num_alignments: i32,
    pub outfmt: u8,
    pub fields: Vec<String>,
    pub num_threads: i32,
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3642
// ```c++
//     CRef<CBlastOptionsHandle> retval(x_CreateOptionsHandle(locality, args));
//     CBlastOptions& opts = retval->SetOptions();
//     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
//         (*arg)->ExtractAlgorithmOptions(args, opts);
//     }
//
//     m_IsUngapped = !opts.GetGappedMode();
//     try { retval->Validate(); }
//     catch (const CBlastException& e) {
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
//     }
//     return retval;
// }
//
// ```
impl BlastxArgs {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3629-3642
    // ```c++
    //     CRef<CBlastOptionsHandle> retval(x_CreateOptionsHandle(locality, args));
    //     CBlastOptions& opts = retval->SetOptions();
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    //     catch (const CBlastException& e) {
    //         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    //     }
    //     return retval;
    // }
    //
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2940-2947
    // ```c++
    //     	    m_LineLength = args[kArgLineLength].AsInteger();
    //     	}
    //         if(args.Exist(kArgSortHits) && args[kArgSortHits])
    //         {
    //        	    m_HitsSortOption = args[kArgSortHits].AsInteger();
    //         }
    //     }
    //     else
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310
    // ```c++
    //     if (program_number != eBlastTypeBlastn &&
    //         program_number != eBlastTypeMapping &&
    //         (!Blast_ProgramIsRpsBlast(program_number)) &&
    //         options->threshold <= 0)
    //     {
    //         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //                          "Non-zero threshold required");
    //         return BLASTERR_OPTION_VALUE_INVALID;
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
    // ```c++
    // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
    // 	{
    // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //          "expect value or cutoff score must be greater than zero");
    // 		return BLASTERR_OPTION_VALUE_INVALID;
    // 	}
    // ```
    pub fn resolve(&self) -> Result<ResolvedOptions> {
        if self.query.as_os_str() == "-" || self.subject.as_os_str() == "-" {
            bail!("unsupported BLASTX stdin query/subject");
        }
        if self.matrix.as_deref().is_some_and(|m| m != "BLOSUM62") {
            bail!("unsupported BLASTX matrix: use BLOSUM62");
        }
        if self.gap_open.is_some_and(|g| g != 11) || self.gap_extend.is_some_and(|g| g != 1) {
            bail!("unsupported BLASTX gap costs: use 11/1");
        }
        if self
            .out
            .as_ref()
            .is_some_and(|p| p.file_name().is_some_and(|n| n.len() >= 256))
        {
            bail!("BLASTX -out filename length must be < 256");
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:872-877
        // ```c++
        //         if (ungapped && *ungapped && compo_mode != eNoCompositionBasedStats) {
        //             NCBI_THROW(CInputException, eInvalidInput,
        //                        "Composition-adjusted searched are not supported with "
        //                        "an ungapped search, please add -comp_based_stats F or "
        //                        "do a gapped search");
        //         }
        // ```
        if self.ungapped && self.comp_based_stats != 0 {
            return Err(SemanticOptionsError("Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1765-1773
        // ```c++
        //    if ((status = LookupTableOptionsValidate(program_number,
        //                     lookup_options, blast_msg)) != 0)
        //        return status;
        //    if ((status = BlastInitialWordOptionsValidate(program_number,
        //                     word_options, blast_msg)) != 0)
        //        return status;
        //    if ((status = BlastHitSavingOptionsValidate(program_number, hit_options,
        //                                                blast_msg)) != 0)
        //        return status;
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1303-1310
        // ```c++
        //     if (program_number != eBlastTypeBlastn &&
        //         program_number != eBlastTypeMapping &&
        //         (!Blast_ProgramIsRpsBlast(program_number)) &&
        //         options->threshold <= 0)
        //     {
        //         Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //                          "Non-zero threshold required");
        //         return BLASTERR_OPTION_VALUE_INVALID;
        // ```
        if self.threshold.is_some_and(|x| x <= 0.0) {
            return Err(SemanticOptionsError("Non-zero threshold required").into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
        // ```c++
        // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
        // 	{
        // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
        //          "expect value or cutoff score must be greater than zero");
        // 		return BLASTERR_OPTION_VALUE_INVALID;
        // 	}
        // ```
        if self.evalue <= 0.0 {
            return Err(SemanticOptionsError(
                "expect value or cutoff score must be greater than zero",
            )
            .into());
        }
        let word_size = self.word_size.unwrap_or(3);
        let outfmt = self
            .outfmt
            .split_whitespace()
            .next()
            .unwrap()
            .parse::<u8>()?;
        let hitlist_size = self.max_target_seqs.unwrap_or(500);
        let mut fields = Vec::new();
        let spec: Vec<_> = self.outfmt.split_whitespace().skip(1).collect();
        if outfmt != 0 {
            let spec = if spec.is_empty() { vec!["std"] } else { spec };
            for token in spec {
                let expanded: Vec<_> = if token == "std" {
                    STD_FIELDS.split_whitespace().collect()
                } else {
                    vec![token]
                };
                for field in expanded {
                    if !fields.iter().any(|f| f == field) {
                        fields.push(field.to_string());
                    }
                }
            }
        }
        Ok(ResolvedOptions {
            word_size,
            threshold: self
                .threshold
                .unwrap_or(if word_size == 5 { 19.3 } else { 12.0 }),
            lookup_type: if word_size == 5 { 4 } else { 3 },
            matrix: "BLOSUM62".into(),
            gap_open: 11,
            gap_extend: 1,
            gapped: !self.ungapped,
            composition: self.comp_based_stats,
            seg: self.seg.clone(),
            soft_masking: self.soft_masking,
            strand: self.strand.clone(),
            query_gencode: self.query_gencode,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
            // ```c++
            //
            // void
            // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
            // {
            //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
            //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
            // }
            //
            // void
            // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
            // {
            //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
            //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
            //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
            //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
            //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
            // }
            // ```
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:123-145
            // ```c++
            //                                              searches except blastn */
            // #define BLAST_UNGAPPED_X_DROPOFF_NUCL 20 /**< ungapped dropoff score for
            //                                               blastn (and megablast) */
            //
            // /** default dropoff for preliminary gapped extensions */
            // #define BLAST_GAP_X_DROPOFF_PROT 15 /**< default dropoff (all protein-
            //                                          based gapped extensions) */
            // #define BLAST_GAP_X_DROPOFF_NUCL 30 /**< default dropoff for non-greedy
            //                                          nucleotide gapped extensions */
            // #define BLAST_GAP_X_DROPOFF_GREEDY 25 /**< default dropoff for greedy
            //                                          nucleotide gapped extensions */
            // #define BLAST_GAP_X_DROPOFF_TBLASTX 0 /**< default dropoff for tblastx */
            //
            // /** default bit score that will trigger gapped extension */
            // #define BLAST_GAP_TRIGGER_PROT 22.0 /**< default bit score that will trigger
            //                                          a gapped extension for all protein-
            //                                          based searches */
            // #define BLAST_GAP_TRIGGER_NUCL 27.0  /**< default bit score that will trigger
            //                                          a gapped extension for blastn */
            //
            // /** default dropoff for the final gapped extension with traceback */
            // #define BLAST_GAP_X_DROPOFF_FINAL_PROT 25 /**< default dropoff (all protein-
            //                                                based gapped extensions) */
            // ```
            x_dropoff: 7.0,
            gap_x_dropoff: 15.0,
            gap_x_dropoff_final: 25.0,
            gap_trigger: 22.0,
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
            // ```c++
            // void
            // CBlastProteinOptionsHandle::SetEffectiveLengthsOptionsDefaults()
            // {
            //     SetDbLength(0);
            //     SetDbSeqNum(0);
            //     SetEffectiveSearchSpace(0);
            // }
            // ```
            db_length: 0,
            effective_search_space: 0,
            window_size: self.window_size.unwrap_or(40),
            evalue: self.evalue,
            sum_stats: self.sum_stats.unwrap_or(true),
            max_intron_length: self.max_intron_length,
            max_hsps: self.max_hsps.unwrap_or(0),
            culling_limit: self.culling_limit,
            subject_besthit: self.subject_besthit,
            hitlist_size,
            max_target_seqs_explicit: self.max_target_seqs.is_some(),
            num_descriptions: hitlist_size,
            num_alignments: if outfmt == 0 && self.max_target_seqs.is_none() {
                250
            } else {
                hitlist_size
            },
            outfmt,
            fields,
            num_threads: self.num_threads,
        })
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:195-220
    // ```c++
    //         /*** Initialize the database/subject ***/
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope;
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    //         _ASSERT(db_adapter && scope);
    //
    //         /*** Get the query sequence(s) ***/
    //         CRef<CQueryOptionsArgs> query_opts =
    //             m_CmdLineArgs->GetQueryOptionsArgs();
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    //         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                      query_opts->UseLowercaseMasks(),
    //                                      query_opts->GetParseDeflines(),
    //                                      query_opts->GetRange());
    //         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //            	ERR_POST(Warning << "Query is Empty!");
    //            	return BLAST_EXIT_SUCCESS;
    //         }
    //         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
    //         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
    //
    //         /*** Get the formatting options ***/
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:276-290
    // ```c++
    // 	        BLAST_PROF_START( APP.LOOP.BLAST );
    //                 CLocalBlast lcl_blast(queries, m_OptsHndl, db_adapter);
    //                 lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
    //                 results = lcl_blast.Run();
    // 	        BLAST_PROF_STOP( APP.LOOP.BLAST );
    //             }
    // 	    BLAST_PROF_START( APP.LOOP.FMT );
    //             if (fmt_args->ArchiveFormatRequested(args)) {
    //                 formatter.WriteArchive(*queries, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
    //                 m_Bah.ResetMessages();
    //             } else {
    //                 BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		                            fmt_args->GetFormattedOutputChoice());
    //             	ITERATE(CSearchResultSet, result, *results) {
    //                	    formatter.PrintOneResultSet(**result, query_batch);
    // ```
    pub fn run(self) -> Result<()> {
        super::native::run(&self)
    }
}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3630-3642
// ```c++
//     CBlastOptions& opts = retval->SetOptions();
//     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
//         (*arg)->ExtractAlgorithmOptions(args, opts);
//     }
//
//     m_IsUngapped = !opts.GetGappedMode();
//     try { retval->Validate(); }
//     catch (const CBlastException& e) {
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
//     }
//     return retval;
// }
//
// ```
mod tests {
    use super::*;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3630-3642
    // ```c++
    //     CBlastOptions& opts = retval->SetOptions();
    //     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
    //         (*arg)->ExtractAlgorithmOptions(args, opts);
    //     }
    //
    //     m_IsUngapped = !opts.GetGappedMode();
    //     try { retval->Validate(); }
    //     catch (const CBlastException& e) {
    //         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
    //     }
    //     return retval;
    // }
    //
    // ```
    fn parse(options: &[&str]) -> BlastxArgs {
        let mut argv = vec!["LOSAT", "blastx", "-query", "q", "-subject", "s"];
        argv.extend_from_slice(options);
        let cli: crate::cli::Cli = crate::cli::try_parse_from(argv).unwrap();
        match cli.command {
            crate::cli::Commands::Blastx(a) => a,
            _ => unreachable!(),
        }
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1300-1306
    // ```c++
    //     if (kPhiBlast)
    //         return 0;
    //
    //     if (program_number != eBlastTypeBlastn &&
    //         program_number != eBlastTypeMapping &&
    //         (!Blast_ProgramIsRpsBlast(program_number)) &&
    //         options->threshold <= 0)
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
    // ```c++
    // 	if (options->expect_value <= 0.0 && options->cutoff_score <= 0)
    // 	{
    // 		Blast_MessageWrite(blast_msg, eBlastSevError, kBlastMessageNoContext,
    //          "expect value or cutoff score must be greater than zero");
    // 		return BLASTERR_OPTION_VALUE_INVALID;
    // 	}
    // ```
    fn cli_acceptance_precedes_core_validation() {
        assert_eq!(parse(&["-threshold", "0"]).threshold, Some(0.0));
        assert!(parse(&["-threshold", "0"]).resolve().is_err());
        assert_eq!(parse(&["-evalue", "0"]).evalue, 0.0);
        assert!(parse(&["-evalue", "0"]).resolve().is_err());
    }
    #[test]
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2910-2928
    // ```c++
    //     if(m_OutputFormat <= eFlatQueryAnchoredNoIdentities) {
    //
    //
    //     	 m_NumDescriptions = m_DfltNumDescriptions;
    //     	 m_NumAlignments = m_DfltNumAlignments;
    //
    //     	 if (args.Exist(kArgNumDescriptions) && args[kArgNumDescriptions]) {
    //     	    m_NumDescriptions = args[kArgNumDescriptions].AsInteger();
    //     	 }
    //
    //          if (args.Exist(kArgNumAlignments) && args[kArgNumAlignments]) {
    //     		m_NumAlignments = args[kArgNumAlignments].AsInteger();
    //     	}
    //
    //          if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //     	    m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //     		m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    //     		hitlist_size = m_NumAlignments;
    //          }
    // ```
    fn omission_and_explicit_target_count_remain_distinct() {
        let omitted = parse(&[]).resolve().unwrap();
        let explicit = parse(&["-max_target_seqs", "500"]).resolve().unwrap();
        assert_eq!(
            (
                omitted.hitlist_size,
                omitted.num_descriptions,
                omitted.num_alignments
            ),
            (500, 500, 250)
        );
        assert_eq!(
            (
                explicit.hitlist_size,
                explicit.num_descriptions,
                explicit.num_alignments
            ),
            (500, 500, 500)
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3636-3639
// ```c++
//     try { retval->Validate(); }
//     catch (const CBlastException& e) {
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
//     }
// ```
// The native app catches this category after stream extraction. Product scope
// refusals retain the separate registered CLI-v2 boundary.
#[derive(Debug)]
pub(crate) struct SemanticOptionsError(pub &'static str);
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-175
// ```c++
//     catch (const blast::CInputException& e) {                               \
//         LOG_POST(Error << "BLAST query/options error: " << e.GetMsg());     \
//         LOG_POST(Error << "Please refer to the BLAST+ user manual.");       \
//         exit_code = BLAST_INPUT_ERROR;                                      \
// ```
impl std::fmt::Display for SemanticOptionsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.0)
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3638-3638
// ```c++
//         NCBI_THROW(CInputException, eInvalidInput, e.GetMsg());
// ```
impl std::error::Error for SemanticOptionsError {}
