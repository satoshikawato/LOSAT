use super::*;
use crate::cli::{try_parse_from, Cli, Commands};
use crate::report::archive::options::{self, algorithm_options, ParamValue};
use crate::report::serial::real::double_to_string_posix;

// The `algorithm-options` that NCBI BLAST+ 2.17.0 wrote in the OF8 design's probes
// (S/probe8/out, `check_archive_layout.py opts`), without the `EffectiveSearchSpace` that
// each archive adds: the probe's name, the options given besides the inputs, and the
// list as `Name=variant:value`.
const CASES: &[(&str, &[&str], &str)] = &[
("b_n_allN_db_fmt11", &[], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
// DI-6 (the frozen archives xdrop.*.11 of tests/outfmt0_manifest.tsv): the X-drops and the window
// as given, -no_greedy's algorithms (eDynProgScoreOnly, eDynProgTbck = 0), the window of
// dc-megablast replaced in its place, pyani ANIb's options.
("di6_nogreedy", &["-no_greedy", "-gapopen", "5", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:5 GapExtensionCost=integer:2 GapExtnAlgorithm=integer:0 GapTracebackAlgorithm=integer:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("di6_xdrops", &["-xdrop_ungap", "5", "-xdrop_gap", "12.5", "-xdrop_gap_final", "150", "-window_size", "40"], "EvalueThreshold=cutoff:e-value:10 GapXDropoff=real:12.5 GapXDropoffFinal=real:150 MaskAtHash=boolean:True UngappedMode=boolean:False WindowSize=integer:40 HitlistSize=integer:500"),
("di6_dc_win0", &["-task", "dc-megablast", "-window_size", "0"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:0 MBTemplateLength=integer:18 WordSize=integer:11 WindowSize=integer:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("di6_anib", &["-task", "blastn", "-xdrop_gap_final", "150", "-dust", "no", "-evalue", "1e-15", "-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:1e-15 GapXDropoffFinal=real:150 MaskAtHash=boolean:True DustFiltering=boolean:False UngappedMode=boolean:False HitlistSize=integer:1"),
// DI-6b (the frozen archives minraw.*.11): the value is not recorded (the remote options'
// eBlastOpt_CutoffScore sits under a condition that is always false, api/blast_options_cxx.cpp:435-442).
("di6b_minraw", &["-min_raw_gapped_score", "15"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("e_n_badpair_fmt11", &["-reward", "2", "-penalty", "-7"], "EvalueThreshold=cutoff:e-value:10 MismatchPenalty=integer:-7 MatchReward=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("in_n_dust_minus_fmt11", &["-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("in_n_lcase_mask_fmt11", &["-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
// SM-3 (the frozen archives sm3.*.11 of tests/outfmt0_manifest.tsv): the value as given.
("sm3_hard", &["-soft_masking", "false"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:False UngappedMode=boolean:False HitlistSize=integer:500"),
("sm3_soft", &["-soft_masking", "true"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("in_n_mts1_fmt11", &["-task", "blastn", "-max_target_seqs", "1"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:1"),
("in_n_qloc_fmt11", &["-query_loc", "10-300"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("in_n_qloc_full_fmt11", &["-query_loc", "1-480"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("m_n_dust_plus_fmt11", &["-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("m_n_lc_minus_fmt11", &["-lcase_masking", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g1", &["-word_size", "20", "-evalue", "1e-5", "-perc_identity", "90", "-max_hsps", "1", "-qcov_hsp_perc", "50"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:20 PercentIdentity=real:90 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g10", &["-gapopen", "2", "-gapextend", "2"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:2 GapExtensionCost=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g1rev", &["-qcov_hsp_perc", "50", "-max_hsps", "1", "-perc_identity", "90", "-evalue", "1e-5", "-word_size", "20"], "EvalueThreshold=cutoff:e-value:1e-05 WordSize=integer:20 PercentIdentity=real:90 QueryCovHspPerc=real:50 MaxHspsPerSubject=integer:1 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g2", &["-penalty", "-3", "-culling_limit", "1", "-strand", "plus", "-dust", "no"], "EvalueThreshold=cutoff:e-value:10 MismatchPenalty=integer:-3 MaskAtHash=boolean:True DustFiltering=boolean:False UngappedMode=boolean:False Culling=integer:1 HitlistSize=integer:500"),
// NP O6a (NCBI BLAST+ 2.17.0, O/O6/probe_a): each best hit option when given, between
// `Culling` (which they exclude) and `UseSubjectBestHit`.
("o6a_besthit_both", &["-subject_besthit", "-best_hit_score_edge", "0.05", "-best_hit_overhang", "0.2"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False BestHitOverhang=real:0.2 BestHitScoreEdge=real:0.05 UseSubjectBestHit=boolean:True HitlistSize=integer:500"),
("o6a_besthit_edge", &["-best_hit_score_edge", "0.3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False BestHitScoreEdge=real:0.3 HitlistSize=integer:500"),
("opt_blastn_g3", &["-reward", "2", "-penalty", "-3", "-gapopen", "5", "-gapextend", "2", "-subject_besthit", "-lcase_masking"], "EvalueThreshold=cutoff:e-value:10 GapOpeningCost=integer:5 GapExtensionCost=integer:2 MismatchPenalty=integer:-3 MatchReward=integer:2 MaskAtHash=boolean:True UngappedMode=boolean:False UseSubjectBestHit=boolean:True HitlistSize=integer:500"),
("opt_blastn_g4", &["-max_target_seqs", "1", "-dust", "20 64 1", "-query_loc", "10-300", "-subject_loc", "10-300", "-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:True DustFilteringLevel=integer:20 DustFilteringWindow=integer:64 DustFilteringLinker=integer:1 UngappedMode=boolean:False HitlistSize=integer:1"),
("opt_blastn_g5", &["-num_descriptions", "5", "-num_alignments", "5", "-strand", "minus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:5"),
("opt_blastn_g6", &["-task", "dc-megablast", "-template_type", "coding", "-template_length", "16", "-word_size", "12"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:0 MBTemplateLength=integer:16 WordSize=integer:12 WindowSize=integer:40 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g7", &["-task", "blastn-short"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("short_max_hsps1", &["-task", "blastn-short", "-max_hsps", "1"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True MaxHspsPerSubject=integer:1 UngappedMode=boolean:False HitlistSize=integer:500"),
("short_perc_identity90", &["-task", "blastn-short", "-perc_identity", "90"], "EvalueThreshold=cutoff:e-value:1000 MatchReward=integer:1 MismatchPenalty=integer:-3 WordSize=integer:7 DustFiltering=boolean:False SegFiltering=boolean:False RepeatFiltering=boolean:False MaskAtHash=boolean:True PercentIdentity=real:90 UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g8", &["-task", "blastn", "-word_size", "7", "-dust", "yes"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:7 MaskAtHash=boolean:True DustFiltering=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("opt_blastn_g9", &["-task", "megablast", "-word_size", "28", "-evalue", "10", "-strand", "both", "-dust", "yes"], "EvalueThreshold=cutoff:e-value:10 WordSize=integer:28 MaskAtHash=boolean:True DustFiltering=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("optdb_blastn_g4", &["-max_target_seqs", "1", "-dust", "20 64 1", "-query_loc", "10-300", "-dbsize", "1000000"], "EvalueThreshold=cutoff:e-value:10 DbLength=big-integer:1000000 MaskAtHash=boolean:True DustFilteringLevel=integer:20 DustFilteringWindow=integer:64 DustFilteringLinker=integer:1 UngappedMode=boolean:False HitlistSize=integer:1"),
("t_n_cull0_fmt11", &["-culling_limit", "0"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False Culling=integer:0 HitlistSize=integer:500"),
("t_n_dmb_opt_fmt11", &["-task", "dc-megablast", "-template_type", "optimal", "-template_length", "21"], "EvalueThreshold=cutoff:e-value:10 MBTemplateType=integer:1 MBTemplateLength=integer:21 WordSize=integer:11 WindowSize=integer:40 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("t_n_naln_fmt11", &["-num_alignments", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:3"),
("t_n_ndesc_fmt11", &["-num_descriptions", "3"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("t_n_pid0_fmt11", &["-perc_identity", "0"], "EvalueThreshold=cutoff:e-value:10 PercentIdentity=real:0 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
("t_n_strandplus_lc_fmt11", &["-lcase_masking", "-strand", "plus"], "EvalueThreshold=cutoff:e-value:10 MaskAtHash=boolean:True UngappedMode=boolean:False HitlistSize=integer:500"),
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
            format!("{}={value}", param.name)
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn parse(extra: &[&str]) -> Option<BlastnArgs> {
    let mut argv = vec!["losat", "blastn", "-query", "q.fa", "-subject", "s.fa"];
    argv.extend_from_slice(extra);
    argv.extend_from_slice(&["-outfmt", "11"]);
    let cli: Cli = try_parse_from(argv).ok()?;
    match cli.command {
        Commands::Blastn(args) => Some(args),
        _ => None,
    }
}

#[test]
fn the_option_lists_equal_those_of_ncbis_probes() {
    let mut compared = 0;
    for &(name, extra, expected) in CASES {
        let args = parse(extra).unwrap_or_else(|| panic!("{name} {extra:?} does not parse"));
        let dust = args.dust_filtering.clone();
        let inputs = option_inputs(&args, dust.as_deref());
        assert_eq!(
            render(&algorithm_options(&inputs)),
            expected,
            "{name} {extra:?}"
        );
        compared += 1;
    }
    assert!(compared >= 27, "{compared}");
}

fn dust(ranges: &[(usize, usize)]) -> Vec<MaskedInterval> {
    ranges
        .iter()
        .map(|&(start, end)| MaskedInterval::new(start, end))
        .collect()
}

#[test]
fn lcase_mask_is_the_query_vector_mask() {
    use crate::blastinput::fasta_reader::RecordSource;
    // Lower case at 0..20 of the record, the range starts at 29: the lower case outside
    // the range is kept, in record coordinates, once per strand.
    let mut whole = vec![b'a'; 20];
    whole.extend(std::iter::repeat_n(b'C', 280));
    let mut record = FastaRecord::new("Query_1", b"q", &whole[29..]);
    record.source = Some(Box::new(RecordSource {
        defline: b">q".to_vec(),
        uncut: Some(whole),
        ..RecordSource::default()
    }));
    assert_eq!(
        lcase_mask_intervals(&record, &[], 29, true, false, true, false),
        vec![(0, 19), (0, 19)]
    );
    // Without -lcase_masking there is no user mask.
    assert!(lcase_mask_intervals(&record, &[], 29, false, false, true, false).is_empty());
    // DUST intervals of the range are mapped by the range's start, merged with the lower
    // case when they overlap or abut, and sorted.
    assert_eq!(
        lcase_mask_intervals(
            &record,
            &dust(&[(0, 5), (30, 40)]),
            20,
            true,
            false,
            true,
            false
        ),
        vec![(0, 24), (0, 24), (50, 59), (50, 59)]
    );
    // A DUST-only mask: ascending once on a single strand, descending with -strand minus,
    // twice with both strands.
    let plain = FastaRecord::new("Query_2", b"q", &[b'A'; 100]);
    let intervals = dust(&[(5, 10), (20, 30)]);
    assert_eq!(
        lcase_mask_intervals(&plain, &intervals, 0, false, false, false, false),
        vec![(5, 9), (20, 29)]
    );
    assert_eq!(
        lcase_mask_intervals(&plain, &intervals, 0, false, false, false, true),
        vec![(20, 29), (5, 9)]
    );
    assert_eq!(
        lcase_mask_intervals(&plain, &intervals, 0, false, false, true, false),
        vec![(5, 9), (5, 9), (20, 29), (20, 29)]
    );
    // A query with lower case is ascending even with -strand minus.
    let mut lower = vec![b'A'; 100];
    lower[50] = b'c';
    let lower = FastaRecord::new("Query_3", b"q", &lower);
    assert_eq!(
        lcase_mask_intervals(&lower, &intervals, 0, true, false, false, true),
        vec![(5, 9), (5, 9), (20, 29), (20, 29), (50, 50), (50, 50)]
    );
    // WindowMasker's masks replace the others: each interval twice (the plus and minus
    // frames), in the masker's order, with every strand option (NCBI 2.17.0's archives of
    // `-window_masker_db`, line NP M2).
    for (both, minus) in [(true, false), (false, false), (false, true)] {
        assert_eq!(
            lcase_mask_intervals(&lower, &intervals, 3, true, true, both, minus),
            vec![(8, 12), (8, 12), (23, 32), (23, 32)]
        );
    }
}

#[test]
fn the_service_follows_the_task() {
    for (task, service) in [
        ("megablast", "megablast"),
        ("dc-megablast", "megablast"),
        ("blastn", "plain"),
        ("blastn-short", "plain"),
    ] {
        assert_eq!(options::service(Program::Blastn, task), service);
    }
}

#[test]
fn error_lines_are_ncbis_error_texts_only() {
    let native = |message: &str| {
        anyhow::Error::from(crate::cli::NativeError {
            exit: 1,
            message: message.to_string().into(),
        })
    };
    let lines = error_lines(
        &native("LOSAT query/options error: bad  option\nPlease refer to the LOSAT user manual.\n"),
        false,
    )
    .unwrap();
    assert_eq!(
        lines,
        vec![
            (ERROR, b"LOSAT query/options error: bad option".to_vec()),
            (ERROR, b"Please refer to the LOSAT user manual.".to_vec()),
        ]
    );
    assert!(error_lines(&native("Command line argument error: x\n"), false).is_none());
    let subject = error_lines(
        &native(
            "Command line argument error: Argument \"subject\". File is not accessible:  `a.fa'\n",
        ),
        true,
    )
    .unwrap();
    assert_eq!(
        subject,
        vec![(
            ERROR,
            b"Command line argument error: Argument \"subject\". File is not accessible: `a.fa'"
                .to_vec()
        )]
    );
    assert!(error_lines(&anyhow::anyhow!("unsupported by LOSAT"), false).is_none());
    // `CATCH_ALL`'s report of a `std::exception` (exit status 255, blast_app_util.hpp:260-263):
    // one error message (line RG D3-1, `Error: vector::reserve` of the `-in_pssm` reader).
    let std_exception =
        anyhow::Error::from(crate::cli::NativeError::std_exception("vector::reserve"));
    assert_eq!(
        error_lines(&std_exception, false).unwrap(),
        vec![(ERROR, b"Error: vector::reserve".to_vec())]
    );
    // A text with `Error: ` and another exit status is not one of `CATCH_ALL`'s reports.
    assert!(error_lines(&native("Error: something\n"), false).is_none());
}
