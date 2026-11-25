use super::*;
use crate::cli::{try_parse_from, Cli, Commands};
use crate::report::archive::options::{algorithm_options, ParamValue};
use crate::report::serial::real::double_to_string_posix;

// The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for blastp (S/of8/a2/out, the OF8a2
// probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
// archive adds: the name, the options given besides the inputs, and the list as
// `Name=variant:value`.
const CASES: &[(&str, &[&str], &str)] = &[
    ("def", &[], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("fast", &["-task", "blastp-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("short", &["-task", "blastp-short"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbs0", &["-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("cbs1", &["-comp_based_stats", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
    ("cbs3", &["-comp_based_stats", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
    ("cbs3sw", &["-comp_based_stats", "3", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:3 SmithWatermanMode=boolean:True"),
    ("cbsF", &["-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("cbst", &["-comp_based_stats", "t"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbsD", &["-comp_based_stats", "D"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbs2u", &["-comp_based_stats", "2u"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 UnifiedP=integer:1 SmithWatermanMode=boolean:False"),
    ("cbs1u", &["-comp_based_stats", "1U"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:1 UnifiedP=integer:1 SmithWatermanMode=boolean:False"),
    ("cbs0u", &["-comp_based_stats", "0u"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("cbs3u_sw", &["-comp_based_stats", "3u", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:3 UnifiedP=integer:1 SmithWatermanMode=boolean:True"),
    ("segyes", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("segno", &["-seg", "no"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("segpar", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    // SM-1 probes (ncbi-parity-20261019/SM/probe/arch): `-soft_masking` sets `MaskAtHash` after
    // the task's `ClearFilterOptions` (blast_args.cpp:387-391).
    ("opt_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_softf", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_soft_seg", &["-soft_masking", "true", "-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True SegFiltering=boolean:True WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_soft_fast", &["-soft_masking", "true", "-task", "blastp-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 MaskAtHash=boolean:True SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_soft_short", &["-soft_masking", "true", "-task", "blastp-short"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("segyes_short", &["-seg", "yes", "-task", "blastp-short"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:True RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt1", &["-evalue", "1e-5", "-threshold", "12", "-word_size", "2", "-max_hsps", "1", "-qcov_hsp_perc", "50", "-window_size", "30"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:2 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:12 WindowSize=integer:30 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws5", &["-word_size", "5"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:19 WordSize=integer:5 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws6", &["-word_size", "6"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:21 WordSize=integer:6 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws7", &["-word_size", "7", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:15 WordSize=integer:7 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_matrix80", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM80 WordThreshold=integer:12 WindowSize=integer:25 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_matrix45", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM45 WordThreshold=integer:14 WindowSize=integer:60 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_pam30", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:PAM30 WordThreshold=integer:16 WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_pam70", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:PAM70 WordThreshold=integer:14 WindowSize=integer:20 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_gaps", &["-gapopen", "9", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_dbsize", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_w0", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:0 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_mts1", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:1 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_mts4", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:4 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_naln3", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:3 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ndesc", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ndescaln", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:4 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_thr_frac", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ev", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_xdrop", &["-matrix", "BLOSUM62", "-gapopen", "11", "-gapextend", "1"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM62 WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_fast_ws3", &["-task", "blastp-fast", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 WordThreshold=integer:20 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_fast_thr", &["-task", "blastp-fast", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:15 MaskAtHash=boolean:False SegFiltering=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_short_matrix", &["-task", "blastp-short", "-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:BLOSUM62 GapOpeningCost=integer:11 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:11 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_short_ev", &["-task", "blastp-short", "-evalue", "100", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:100 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:3 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_short_cbs", &["-task", "blastp-short", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:16 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("opt_short_ws6", &["-task", "blastp-short", "-word_size", "6"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:6 WordThreshold=integer:21 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_short_thr", &["-task", "blastp-short", "-threshold", "12"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:PAM30 GapOpeningCost=integer:9 GapExtensionCost=integer:1 WordSize=integer:2 WordThreshold=integer:12 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_short_ws4", &["-task", "blastp-short", "-word_size", "4", "-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:20000 MatrixName=string:BLOSUM62 GapOpeningCost=integer:11 GapExtensionCost=integer:1 WordSize=integer:4 WordThreshold=integer:11 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:False WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws4", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws2_matrix", &["-word_size", "2", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 WordSize=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM80 WordThreshold=integer:12 WindowSize=integer:25 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_gapopen_only", &["-gapopen", "10"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_gapext_only", &["-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapExtensionCost=integer:2 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_fast_matrix", &["-task", "blastp-fast", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 GapOpeningCost=integer:10 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:BLOSUM80 WindowSize=integer:25 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_pam30_thr", &["-matrix", "PAM30", "-threshold", "11"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 MaskAtHash=boolean:False SegFiltering=boolean:False MatrixName=string:PAM30 WordThreshold=integer:11 WindowSize=integer:15 HitlistSize=integer:500 UngappedMode=boolean:False CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ungapped_cbs0", &["-ungapped", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:True CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("opt_ungapped_cbsF", &["-ungapped", "-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 UngappedMode=boolean:True CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
];

fn posix(value: f64) -> String {
    let mut digits = Vec::new();
    double_to_string_posix(value, 15, &mut digits);
    String::from_utf8(digits).unwrap()
}

fn render(params: &[crate::report::archive::options::Param]) -> String {
    params
        .iter()
        .map(|param| {
            let value = match &param.value {
                ParamValue::Evalue(value) => format!("cutoff:e-value:{}", posix(*value)),
                ParamValue::Integer(value) => format!("integer:{value}"),
                ParamValue::BigInteger(value) => format!("big-integer:{value}"),
                ParamValue::Boolean(true) => "boolean:True".to_string(),
                ParamValue::Boolean(false) => "boolean:False".to_string(),
                ParamValue::Real(value) => format!("real:{}", posix(*value)),
                ParamValue::Str(value) => format!("string:{value}"),
                ParamValue::IntegerList(values) => format!("integer-list:{values:?}"),
                ParamValue::BigIntegerList(values) => format!("big-integer-list:{values:?}"),
            };
            format!("{}={value}", param.name)
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// The arguments of `extra` after `-query q.fa -subject s.fa`, with `-outfmt 11`.
fn parse(extra: &[&str]) -> Option<BlastpArgs> {
    let mut argv = vec!["losat", "blastp", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend_from_slice(extra);
    argv.extend_from_slice(&["-outfmt", "11"]);
    let cli: Cli = try_parse_from(argv).ok()?;
    match cli.command {
        Commands::Blastp(args) => Some(args),
        _ => None,
    }
}

// A `-threshold` past `Int4` (2147483648, 1e10, +inf, 1e999) is recorded as `INT_MIN`, as
// NCBI BLAST+ 2.17.0 converts it (the lookup table takes the same value). NCBI's build writes
// it as `-8` in the archive, from a signed overflow that LOSAT does not reproduce
// (`report/serial/asn_text.rs` `itoa_i32`).
#[test]
fn a_threshold_past_int4_is_recorded_as_int_min() {
    let inputs = |value: &str| {
        option_inputs(&parse(&["-threshold", value]).unwrap_or_else(|| panic!("{value}")))
    };
    let listed = |threshold: &str| {
        let case = CASES
            .iter()
            .find(|case| case.0 == "opt_thr_frac")
            .expect("opt_thr_frac");
        case.2.replace(
            "WordThreshold=integer:11",
            &format!("WordThreshold=integer:{threshold}"),
        )
    };
    for value in ["2147483648", "1e10", "+inf", "1e999"] {
        assert_eq!(
            render(&algorithm_options(&inputs(value))),
            listed("-2147483648"),
            "{value}"
        );
    }
    assert_eq!(
        render(&algorithm_options(&inputs("2147483647.5"))),
        listed("2147483647")
    );
}

#[test]
fn the_option_lists_equal_those_of_ncbis_probes() {
    assert!(CASES.len() >= 50);
    for &(name, extra, expected) in CASES {
        let args = parse(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
        assert_eq!(
            render(&algorithm_options(&option_inputs(&args))),
            expected,
            "{name} {extra:?}"
        );
    }
}

#[test]
fn the_service_is_plain() {
    assert_eq!(
        crate::report::archive::options::service(Program::Blastp, "blastp-short"),
        "plain"
    );
}
