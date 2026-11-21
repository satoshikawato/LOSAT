use super::*;

fn names(inputs: &OptionInputs) -> Vec<String> {
    algorithm_options(inputs)
        .iter()
        .map(|param| param.name().to_string())
        .collect()
}

// NCBI BLAST+ 2.17.0 blastn -window_masker_db (line NP M2): the parameter follows the
// filtering handler's MaskAtHash and DUST calls, also after blastn-short's cleared filters.
#[test]
fn the_window_masker_database_follows_the_dust_options() {
    let mut inputs = OptionInputs::new(Program::Blastn, "blastn");
    inputs.window_masker_db = Some("counts.obinary".to_string());
    inputs.dust = Some(Filter::No);
    let params = algorithm_options(&inputs);
    let listed: Vec<&str> = params.iter().map(|param| param.name()).collect();
    let at = listed
        .iter()
        .position(|&name| name == "WindowMaskerDatabase")
        .unwrap();
    assert_eq!(&listed[at - 2..at], ["MaskAtHash", "DustFiltering"]);
    assert_eq!(
        params[at].value,
        ParamValue::Str("counts.obinary".to_string())
    );
    let mut short = OptionInputs::new(Program::Blastn, "blastn-short");
    short.window_masker_db = Some("counts.obinary".to_string());
    let listed = names(&short);
    let at = listed
        .iter()
        .position(|name| name == "WindowMaskerDatabase")
        .unwrap();
    assert_eq!(listed[at - 1], "MaskAtHash");
    assert_eq!(listed[at + 1], "UngappedMode");
}

// NCBI BLAST+ 2.17.0 blastn -window_masker_taxid (line NP M3, ws npm/m3 c01 and t00): the
// taxid takes the place of the database's parameter; the taxid 0 writes none
// (`SetWindowMaskerTaxId` resets it).
#[test]
fn the_window_masker_taxid_is_a_parameter_unless_it_is_0() {
    let mut inputs = OptionInputs::new(Program::Blastn, "blastn");
    inputs.window_masker_taxid = Some(562);
    let params = algorithm_options(&inputs);
    let listed: Vec<&str> = params.iter().map(|param| param.name()).collect();
    let at = listed
        .iter()
        .position(|&name| name == "WindowMaskerTaxId")
        .unwrap();
    assert_eq!(listed[at - 1], "MaskAtHash");
    assert_eq!(params[at].value, ParamValue::Integer(562));
    inputs.window_masker_taxid = Some(0);
    assert!(!names(&inputs)
        .iter()
        .any(|name| name == "WindowMaskerTaxId"));
}

#[test]
fn blastn_defaults_are_the_four_parameters_of_the_sample() {
    let inputs = OptionInputs::new(Program::Blastn, "megablast");
    assert_eq!(
        names(&inputs),
        [
            "EvalueThreshold",
            "MaskAtHash",
            "UngappedMode",
            "HitlistSize"
        ]
    );
}

#[test]
fn a_parameter_set_again_keeps_its_place() {
    let mut inputs = OptionInputs::new(Program::Blastn, "blastn-short");
    inputs.evalue = Some(1e-5);
    let list = algorithm_options(&inputs);
    assert_eq!(list[0].name(), "EvalueThreshold");
    assert_eq!(list[0].value, ParamValue::Evalue(1e-5));
    assert_eq!(list[1].name(), "MatchReward");
    assert_eq!(
        list.iter()
            .filter(|p| p.name() == "EvalueThreshold")
            .count(),
        1
    );
}

#[test]
fn matrix_gaps_are_the_default_rows_of_the_matrices() {
    // `BLAST_GetProteinGapExistenceExtendParams`: the row that the matrix's table prefers.
    for (matrix, want) in [
        ("BLOSUM62", (11, 1)),
        ("BLOSUM45", (14, 2)),
        ("BLOSUM50", (13, 2)),
        ("BLOSUM90", (10, 1)),
        ("PAM250", (15, 2)),
        ("pam250", (15, 2)),
        // A matrix that NCBI does not know: both costs stay 0.
        ("BLOSUM99", (0, 0)),
    ] {
        assert_eq!(matrix_gaps(matrix), want, "{matrix}");
    }
}

#[test]
fn the_matrix_suggestions_ignore_the_case_of_the_name() {
    // NCBI compares the names with `strcasecmp` (c++/src/algo/blast/core/blast_options.c:1174-1236).
    for (matrix, threshold, window) in [
        ("BLOSUM62", 11, 40),
        ("blosum62", 11, 40),
        ("BLOSUM45", 14, 60),
        ("Blosum45", 14, 60),
        ("BLOSUM80", 12, 25),
        ("blosum80", 12, 25),
        ("PAM30", 16, 15),
        ("pam30", 16, 15),
        ("Pam70", 14, 20),
    ] {
        assert_eq!(suggested_threshold(matrix), threshold, "{matrix}");
        assert_eq!(suggested_window(matrix), window, "{matrix}");
    }
    assert_eq!(suggested_threshold("Identity"), 27);
    assert_eq!(suggested_threshold("blosum62_20"), 100);
}

#[test]
fn protein_programs_follow_the_class_order() {
    let mut inputs = OptionInputs::new(Program::Blastp, "blastp");
    inputs.seg = None;
    assert_eq!(
        names(&inputs),
        [
            "EvalueThreshold",
            "MaskAtHash",
            "SegFiltering",
            "WordThreshold",
            "WindowSize",
            "HitlistSize",
            "UngappedMode",
            "CompositionBasedStats",
            "SmithWatermanMode"
        ]
    );
}

// NCBI BLAST+ 2.17.0 (line RG F1b, probe 9): the default lists of psiblast
// (`p9.ps.q.db.itdef.f6`), rpsblast (`p9.rps.tk.f6`) and rpstblastn (`p9.rtn.nq4.f6`), with
// their programs and services.
#[test]
fn the_applications_lists_are_their_handles_and_classes() {
    let psiblast = OptionInputs::new(RequestKind::Psiblast, "blastp");
    assert_eq!(
        names(&psiblast),
        [
            "EvalueThreshold",
            "InclusionThreshold",
            "PseudoCountWeight",
            "MaskAtHash",
            "SegFiltering",
            "WordThreshold",
            "WindowSize",
            "HitlistSize",
            "CompositionBasedStats",
            "SmithWatermanMode",
            "GapTrigger",
            "IgnoreMsaMaster"
        ]
    );
    let params = algorithm_options(&psiblast);
    assert_eq!(params[1].value, ParamValue::Real(0.002));
    assert_eq!(params[2].value, ParamValue::Integer(0));
    assert_eq!(params[10].value, ParamValue::Real(22.0));
    assert_eq!(
        names(&OptionInputs::new(RequestKind::Rpsblast, "blastp")),
        [
            "EvalueThreshold",
            "MaskAtHash",
            "SegFiltering",
            "HitlistSize",
            "CompositionBasedStats",
            "SmithWatermanMode"
        ]
    );
    assert_eq!(
        names(&OptionInputs::new(RequestKind::Rpstblastn, "blastx")),
        [
            "EvalueThreshold",
            "QueryGeneticCode",
            "MaskAtHash",
            "SegFiltering",
            "UngappedMode",
            "HitlistSize",
            "CompositionBasedStats",
            "SmithWatermanMode"
        ]
    );
    // tblastn `-in_pssm` (`F/f1d/p7ncbi/tbn.pssm.db.ss`): psiblast's constructor, then
    // `SetPSITblastnDefaults`, then tblastn's classes; mode 0, the code 1 of the constructor.
    let mut psi_tblastn = OptionInputs::new(RequestKind::PsiTblastn, "tblastn-fast");
    psi_tblastn.comp_based_stats = 2;
    psi_tblastn.db_gencode = 4;
    let params = algorithm_options(&psi_tblastn);
    let listed: Vec<String> = params
        .iter()
        .map(|param| format!("{}={:?}", param.name(), param.value))
        .collect();
    assert_eq!(
        listed,
        [
            "EvalueThreshold=Evalue(10.0)",
            "InclusionThreshold=Real(0.002)",
            "PseudoCountWeight=Integer(0)",
            "SmithWatermanMode=Boolean(false)",
            "DbGeneticCode=Integer(1)",
            "WordThreshold=Integer(13)",
            "UngappedMode=Boolean(false)",
            "LongestIntronLength=Integer(0)",
            "MaskAtHash=Boolean(false)",
            "SegFilteringWindow=Integer(12)",
            "SegFilteringLocut=Real(2.2)",
            "SegFilteringHicut=Real(2.5)",
            "WindowSize=Integer(40)",
            "HitlistSize=Integer(500)",
            "CompositionBasedStats=Integer(0)",
        ]
    );
    // The request's program and service, and the application's prefix of its posts.
    for (kind, program, named, application) in [
        (RequestKind::Psiblast, "blastp", "psi", "psiblast"),
        (RequestKind::PsiTblastn, "blastp", "psi", "tblastn"),
        (RequestKind::Rpsblast, "blastp", "rpsblast", "rpsblast"),
        (RequestKind::Rpstblastn, "tblastn", "rpsblast", "rpstblastn"),
        (RequestKind::Blastp, "blastp", "plain", "blastp"),
    ] {
        assert_eq!(kind.program(), program);
        assert_eq!(service(kind, "blastp"), named, "{kind:?}");
        assert_eq!(kind.application(), application, "{kind:?}");
    }
}

// NCBI BLAST+ 2.17.0 (probe 9): psiblast's given options in place (`p9.ps.q.db.opts2.f6`) or
// where their class sets them (`p9.ps.q.db.opts1.f6`, `opts3`, `p9.ps.three.db.it2.ssp.f6`);
// rpsblast's window only when given (`p9.rps.q40.opts.f6`).
#[test]
fn the_applications_options_take_their_places() {
    let mut inputs = OptionInputs::new(RequestKind::Psiblast, "blastp");
    inputs.inclusion_ethresh = 0.01;
    inputs.pseudocount = 2;
    inputs.gap_trigger = 20.0;
    inputs.xdrop_gap = Some(16.0);
    inputs.searchsp = Some(1_000_000);
    inputs.matrix = Some("BLOSUM62".to_string());
    inputs.culling_limit = Some(2);
    inputs.sum_stats = Some(true);
    let params = algorithm_options(&inputs);
    let listed: Vec<&str> = params.iter().map(|param| param.name()).collect();
    assert_eq!(
        listed[..9],
        [
            "EvalueThreshold",
            "InclusionThreshold",
            "PseudoCountWeight",
            "GapOpeningCost",
            "GapExtensionCost",
            "GapXDropoff",
            "EffectiveSearchSpace",
            "SumStatistics",
            "MaskAtHash"
        ]
    );
    assert_eq!(params[1].value, ParamValue::Real(0.01));
    assert_eq!(params[2].value, ParamValue::Integer(2));
    let at = |name: &str| listed.iter().position(|&listed| listed == name).unwrap();
    assert_eq!(at("MatrixName"), at("SegFiltering") + 1);
    assert_eq!(at("Culling"), at("WordThreshold") + 1);
    assert_eq!(params[at("GapTrigger")].value, ParamValue::Real(20.0));
    // `UnifiedP` is BLASTP's only (blast_args.cpp:880-884; review 3 R3-1).
    let unified = |kind: RequestKind| {
        let mut inputs = OptionInputs::new(kind, "blastp");
        inputs.unified_p = true;
        names(&inputs).iter().any(|name| name == "UnifiedP")
    };
    assert!(unified(RequestKind::Blastp));
    assert!(!unified(RequestKind::Psiblast));
    assert!(!unified(RequestKind::Rpsblast));
    let mut rps = OptionInputs::new(RequestKind::Rpsblast, "blastp");
    rps.window_size = Some(30);
    let listed = names(&rps);
    assert_eq!(listed[3], "WindowSize");
    assert_eq!(listed[2], "SegFiltering");
}
