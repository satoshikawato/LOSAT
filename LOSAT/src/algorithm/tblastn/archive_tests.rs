use super::*;
use crate::cli::{try_parse_from, Cli, Commands};
use crate::report::archive::options::{algorithm_options, ParamValue};
use crate::report::serial::real::double_to_string_posix;

// The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote for tblastn (S/of8/b2/out, the OF8b2
// probes; `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that each
// archive adds: the name, the options given besides the inputs, and the list as
// `Name=variant:value`. LOSAT's TBLASTN refuses many of these options at the command line or
// in its limits (word sizes 2 and 4, matrices, gap costs,
// `-comp_based_stats 1` and `3`, `-ungapped`, `-max_intron_length`): the list that NCBI would
// record is what the archive's model knows about them, pinned here.
const CASES: &[(&str, &[&str], &str)] = &[
    ("def", &[], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbs0", &["-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("cbs0u", &["-comp_based_stats", "0u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("cbs1", &["-comp_based_stats", "1"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
    ("cbs1u", &["-comp_based_stats", "1u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:1 SmithWatermanMode=boolean:False"),
    ("cbs2", &["-comp_based_stats", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbs2u", &["-comp_based_stats", "2u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbs3", &["-comp_based_stats", "3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
    ("cbs3sw", &["-comp_based_stats", "3", "-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:True"),
    ("cbs3u", &["-comp_based_stats", "3u"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:3 SmithWatermanMode=boolean:False"),
    ("cbsD", &["-comp_based_stats", "D"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("cbsF", &["-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("cbst", &["-comp_based_stats", "t"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("fast", &["-task", "tblastn-fast"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("gc11", &["-db_gencode", "11"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:11 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("gc2", &["-db_gencode", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:2 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("gc4", &["-db_gencode", "4"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:4 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("mts_warn", &["-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:2 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_besthit", &["-subject_besthit"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 UseSubjectBestHit=boolean:True WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_cull0", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:0 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_cull5", &["-culling_limit", "5"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_cull_mts", &["-culling_limit", "3", "-max_target_seqs", "2"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 Culling=integer:3 WindowSize=integer:40 HitlistSize=integer:2 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    // `-best_hit_overhang` / `-best_hit_score_edge`, each when given, where `Culling` is (line NP
    // O6a-3's strategy probes `a3/tn_ss/tn_qdb.ss`, `tn_q_se.ss`).
    ("opt_bh_both", &["-best_hit_overhang", "0.2", "-best_hit_score_edge", "0.3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 BestHitOverhang=real:0.2 BestHitScoreEdge=real:0.3 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_bh_se", &["-best_hit_score_edge", "0.3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 BestHitScoreEdge=real:0.3 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_dbsize", &["-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ev", &["-evalue", "0.01"], "EvalueThreshold=cutoff:e-value:0.01 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_evbig", &["-evalue", "1000"], "EvalueThreshold=cutoff:e-value:1000 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_fast_m80", &["-task", "tblastn-fast", "-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:20 GapOpeningCost=integer:10 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WindowSize=integer:25 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_fast_thr", &["-task", "tblastn-fast", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:5 WordThreshold=integer:15 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_fast_ws", &["-task", "tblastn-fast", "-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 WordThreshold=integer:20 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_gapext", &["-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_gapopen", &["-gapopen", "10"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_gaps", &["-gapopen", "9", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_hsps", &["-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:10 MaxHspsPerSubject=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_intron0", &["-max_intron_length", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_intron", &["-max_intron_length", "100"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:100 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_m45", &["-matrix", "BLOSUM45"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:16 WindowSize=integer:60 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_m45cbs0", &["-matrix", "BLOSUM45", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:14 GapExtensionCost=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM45 WordThreshold=integer:16 WindowSize=integer:60 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("opt_m62", &["-matrix", "BLOSUM62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM62 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_m62lc", &["-matrix", "blosum62"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:11 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:blosum62 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_m80", &["-matrix", "BLOSUM80"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:BLOSUM80 WordThreshold=integer:14 WindowSize=integer:25 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_many", &["-evalue", "1e-5", "-threshold", "12", "-window_size", "30", "-gapopen", "10", "-gapextend", "1", "-max_hsps", "2", "-qcov_hsp_perc", "30", "-db_gencode", "4", "-seg", "12 2.2 2.5"], "EvalueThreshold=cutoff:e-value:1e-05 GapOpeningCost=integer:10 GapExtensionCost=integer:1 QueryCovHspPerc=real:30 MaxHspsPerSubject=integer:2 DbGeneticCode=integer:4 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:12 WindowSize=integer:30 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_mts1", &["-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:1 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_mts4", &["-max_target_seqs", "4"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:4 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_naln", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:3 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ndesc", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ndescaln", &["-num_descriptions", "3", "-num_alignments", "4"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:4 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_pam30", &["-matrix", "PAM30"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:9 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM30 WordThreshold=integer:18 WindowSize=integer:15 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_pam70", &["-matrix", "PAM70"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:10 GapExtensionCost=integer:1 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 MatrixName=string:PAM70 WordThreshold=integer:16 WindowSize=integer:20 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_parse", &["-parse_deflines"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_qcov", &["-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:10 QueryCovHspPerc=real:50 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:True SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_softf", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_sumstats", &["-sum_stats", "true"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:True DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_sumstatsf", &["-sum_stats", "false"], "EvalueThreshold=cutoff:e-value:10 SumStatistics=boolean:False DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_thr13", &["-threshold", "13"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_thr_frac", &["-threshold", "11.7"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:11 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_thr", &["-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:15 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ungapped_cbs0", &["-ungapped", "-comp_based_stats", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("opt_ungapped_cbsF", &["-ungapped", "-comp_based_stats", "F"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:True LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:0 SmithWatermanMode=boolean:False"),
    ("opt_win0", &["-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:0 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_win", &["-window_size", "30"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:30 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws2", &["-word_size", "2"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:2 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws3", &["-word_size", "3"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:3 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws4", &["-word_size", "4"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:4 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws5", &["-word_size", "5"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:19 WordSize=integer:5 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws6", &["-word_size", "6"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:21 WordSize=integer:6 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_ws7", &["-word_size", "7", "-threshold", "15"], "EvalueThreshold=cutoff:e-value:10 WordThreshold=integer:15 WordSize=integer:7 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_xdrop", &["-xdrop_gap", "10", "-xdrop_gap_final", "30"], "EvalueThreshold=cutoff:e-value:10 GapXDropoff=real:10 GapXDropoffFinal=real:30 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("opt_xdropu", &["-xdrop_ungap", "20"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("qloc", &["-query_loc", "10-100"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("qloc_full", &["-query_loc", "1-140"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("qsloc", &["-query_loc", "10-100", "-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("sloc", &["-subject_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("sw", &["-use_sw_tback"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:12 SegFilteringLocut=real:2.2 SegFilteringHicut=real:2.5 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:True"),
    // `-seg yes` and `-seg "10 1.8 2.1"` (frozen by the dropped rows `sa.arch.tbn.segyes` / `.segpar`, S1).
    ("segyes", &["-seg", "yes"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFiltering=boolean:True WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
    ("segpar", &["-seg", "10 1.8 2.1"], "EvalueThreshold=cutoff:e-value:10 DbGeneticCode=integer:1 UngappedMode=boolean:False LongestIntronLength=integer:0 MaskAtHash=boolean:False SegFilteringWindow=integer:10 SegFilteringLocut=real:1.8 SegFilteringHicut=real:2.1 WordThreshold=integer:13 WindowSize=integer:40 HitlistSize=integer:500 CompositionBasedStats=integer:2 SmithWatermanMode=boolean:False"),
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
                ParamValue::StrandType(value) => format!("strand-type:{value}"),
                ParamValue::BigIntegerList(values) => format!("big-integer-list:{values:?}"),
            };
            format!("{}={value}", param.name())
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// The recorded inputs of a command line, as LOSAT's TBLASTN parses it (`-parse_deflines`, which
/// the archive does not record, left out).
fn inputs_of(extra: &[&str]) -> Option<OptionInputs> {
    let parsed: Vec<&str> = extra
        .iter()
        .copied()
        .filter(|&option| option != "-parse_deflines")
        .collect();
    let mut argv = vec!["losat", "tblastn", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend_from_slice(&parsed);
    argv.extend_from_slice(&["-outfmt", "11"]);
    let cli: Cli = try_parse_from(argv).ok()?;
    let Commands::Tblastn(args) = cli.command else {
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
    assert!(CASES.len() >= 70);
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
        crate::report::archive::options::service(Program::Tblastn, "tblastn-fast"),
        "plain"
    );
}

#[test]
fn lower_case_runs_are_inclusive_record_positions() {
    use crate::algorithm::common::masks::lowercase_runs;
    assert_eq!(lowercase_runs(b"MWmehSVlk"), [(2, 4), (7, 8)]);
    assert_eq!(lowercase_runs(b"mwME"), [(0, 1)]);
    assert_eq!(lowercase_runs(b"MWME"), []);
    assert_eq!(lowercase_runs(b""), []);
}
