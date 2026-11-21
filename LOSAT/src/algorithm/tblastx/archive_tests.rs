use super::*;
use crate::cli::{try_parse_from, Cli, Commands};
use crate::report::archive::options::{algorithm_options, ParamValue};
use crate::report::serial::real::double_to_string_posix;

// The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for tblastx (S/of8/b3/out, the OF8b3
// probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
// archive adds: the name, the options given besides the inputs, and the list as
// `Name=variant:value`. LOSAT's TBLASTX refuses some of these options at the command line
// (`-parse_deflines`; `-xdrop_ungap`, which NCBI does not record, is read since line RG X1,
// `-sum_stats` since line RG X3, `-strand`, which NCBI does not record either, since line NP
// O2b, and `-subject_besthit` since line NP O6b): the list that NCBI would record is what the
// archive's model knows about them, pinned here.
const CASES: &[(&str, &[&str], &str)] = &[
    ("def_sub", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("seg_no_sub", &["-seg", "no"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:False WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("seg_yes_sub", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("seg_par_sub", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("seg_def_sub", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("batch100_sub", &[], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("dgc11_sub", &["-db_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:11 HitlistSize=integer:500"),
    ("dgc2_db", &["-db_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:2 HitlistSize=integer:500"),
    ("dgc4_db", &["-db_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:4 HitlistSize=integer:500"),
    ("gc44_sub", &["-query_gencode", "4", "-db_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:4 DbGeneticCode=integer:4 HitlistSize=integer:500"),
    ("mts_warn_sub", &["-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:2"),
    ("opt_besthit_sub", &["-subject_besthit"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 UseSubjectBestHit=boolean:True WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_cull0_sub", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:0 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_cull2_db", &["-culling_limit", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:2 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_cull5_sub", &["-culling_limit", "5"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:5 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_cull_mts_sub", &["-culling_limit", "3", "-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:3 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:2"),
    ("opt_dbsize_db", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_ev10_sub", &["-evalue", "10"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_ev_sub", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_evbig_sub", &["-evalue", "1000"], "EvalueThreshold=cutoff:e-value:1000 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_hsps_sub", &["-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:10 MaxHspsPerSubject=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_m45_sub", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:16 WindowSize=integer:60 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_m62_sub", &["-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM62 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_m62lc_sub", &["-matrix", "blosum62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:blosum62 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_m80_sub", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WordThreshold=integer:14 WindowSize=integer:25 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_many_db", &["-evalue", "1e-5", "-threshold", "12", "-window_size", "30", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-query_gencode", "4", "-db_gencode", "4", "-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:1e-05 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:12 WindowSize=integer:30 QueryGeneticCode=integer:4 DbGeneticCode=integer:4 HitlistSize=integer:500"),
    ("opt_many_sub", &["-evalue", "1e-5", "-word_size", "3", "-threshold", "12", "-window_size", "30", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-query_gencode", "4", "-db_gencode", "4", "-seg", "12 2.2 2.5", "-strand", "plus", "-culling_limit", "4"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:3 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 Culling=integer:4 WindowSize=integer:30 QueryGeneticCode=integer:4 DbGeneticCode=integer:4 HitlistSize=integer:500"),
    ("opt_mts1_db", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:1"),
    ("opt_mts4_sub", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:4"),
    ("opt_naln_sub", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:3"),
    ("opt_ndesc_sub", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_ndescaln_sub", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:4"),
    ("opt_pam30_sub", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM30 WordThreshold=integer:18 WindowSize=integer:15 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_pam70_sub", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM70 WordThreshold=integer:16 WindowSize=integer:20 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_parse_sub", &["-parse_deflines"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_qcov_sub", &["-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:10 QueryCovHspPerc=real:50 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_soft_sub", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:True SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_softf_sub", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_sumstats_sub", &["-sum_stats", "true"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_sumstatsf_sub", &["-sum_stats", "false"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_thr13_sub", &["-threshold", "13"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_thr_frac_sub", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:11 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_thr_sub", &["-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:15 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_win0_sub", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:0 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_win40_sub", &["-window_size", "40"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_win_sub", &["-window_size", "30"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:30 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_ws2_sub", &["-word_size", "2"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:2 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_ws3_db", &["-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_ws4_sub", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("opt_xdropu_sub", &["-xdrop_ungap", "20"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qgc11_sub", &["-query_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:11 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qgc2_sub", &["-query_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:2 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qgc4_db", &["-query_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:4 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qloc_db", &["-query_loc", "10-200"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qloc_full_sub", &["-query_loc", "1-366"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qloc_minus_sub", &["-query_loc", "10-200", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("qsloc_sub", &["-query_loc", "10-200", "-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("sloc_sub", &["-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("strand_both_sub", &["-strand", "both"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("strand_minus_db", &["-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
    ("strand_plus_db", &["-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 QueryGeneticCode=integer:1 DbGeneticCode=integer:1 HitlistSize=integer:500"),
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
            };
            format!("{}={value}", param.name())
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// The recorded inputs of a command line. The options that LOSAT's TBLASTX does not accept
/// (they stop at the command line) are applied to the inputs by hand.
fn inputs_of(extra: &[&str]) -> Option<OptionInputs> {
    let mut parsed: Vec<&str> = Vec::new();
    let mut index = 0;
    while index < extra.len() {
        match extra[index] {
            "-parse_deflines" => index += 1,
            other => {
                parsed.push(other);
                index += 1;
            }
        }
    }
    let mut argv = vec!["losat", "tblastx", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend_from_slice(&parsed);
    argv.extend_from_slice(&["-outfmt", "11"]);
    let cli: Cli = try_parse_from(argv).ok()?;
    let Commands::Tblastx(args) = cli.command else {
        return None;
    };
    Some(option_inputs(&args))
}

// A `-threshold` past `Int4` (2147483648, 1e10, +inf, 1e999) is recorded as `INT_MIN`, as
// NCBI BLAST+ 2.17.0 converts it (the lookup table takes the same value). NCBI's build writes
// it as `-8` in the archive, from a signed overflow that LOSAT does not reproduce
// (`report/serial/asn_text.rs` `itoa_i32`).
#[test]
fn a_threshold_past_int4_is_recorded_as_int_min() {
    let inputs =
        |value: &str| inputs_of(&["-threshold", value]).unwrap_or_else(|| panic!("{value}"));
    let listed = |threshold: &str| {
        let case = CASES
            .iter()
            .find(|case| case.0 == "opt_thr_frac_sub")
            .expect("opt_thr_frac_sub");
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
    assert!(CASES.len() >= 55);
    for &(name, extra, expected) in CASES {
        let inputs = inputs_of(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
        assert_eq!(
            render(&algorithm_options(&inputs)),
            expected,
            "{name} {extra:?}"
        );
    }
}

#[test]
fn the_service_is_plain() {
    assert_eq!(
        crate::report::archive::options::service(Program::Tblastx, "tblastx"),
        "plain"
    );
}

#[test]
fn results_masks_are_per_frame_in_ncbis_order_without_a_whole_query_mask() {
    // A 30-nt query: frame lengths 10, 9, 9 (BLAST_GetTranslatedProteinLength), so the
    // DNA length is 30; frames come in the order of the contexts, not of the stats.
    let record = FastaRecord::new("Query_1", b"q", &[b'A'; 30]);
    let stats = TblastxQueryStats {
        seg_masks: vec![(-1, vec![(2, 5)]), (1, vec![(2, 5)]), (2, vec![(0, 9)])],
        ..TblastxQueryStats::default()
    };
    let masks = query_results_masks(&record, &stats, 0);
    let frames: Vec<Frame> = masks.iter().map(|mask| mask.frame).collect();
    assert_eq!(frames, [Frame::Plus1, Frame::Plus2, Frame::Minus1]);
    // Residues 2..=4 of frame +1: NCBI's end is 3 * 4 (blast_filter.c:925-931), the formula
    // `query_masks_convert_with_ncbis_frame_formulas` pins.
    assert_eq!(masks[0].intervals, [(6, 12)]);
    // A mask that spans the whole query is dropped (blast_aux.cpp:936-951), and a
    // `-query_loc` start moves the rest into the record.
    let moved = query_results_masks(&record, &stats, 100);
    assert_eq!(moved[0].intervals, [(106, 112)]);
}

#[test]
fn frames_rank_as_ncbis_contexts() {
    let ranks: Vec<u8> = [1i8, 2, 3, -1, -2, -3]
        .into_iter()
        .map(frame_rank)
        .collect();
    assert_eq!(ranks, [1, 2, 3, 4, 5, 6]);
}
