use super::*;
use crate::report::archive::read::options::{BuiltOptions, EProgram, Setter};

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
    inputs.window_masker_db = Some(b"counts.obinary".to_vec());
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
        ParamValue::Str(b"counts.obinary".to_vec())
    );
    let mut short = OptionInputs::new(Program::Blastn, "blastn-short");
    short.window_masker_db = Some(b"counts.obinary".to_vec());
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
        assert_eq!(matrix_gaps(matrix.as_bytes()), want, "{matrix}");
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
        assert_eq!(
            suggested_threshold(matrix.as_bytes()),
            threshold,
            "{matrix}"
        );
        assert_eq!(suggested_window(matrix.as_bytes()), window, "{matrix}");
    }
    assert_eq!(suggested_threshold(b"Identity"), 27);
    assert_eq!(suggested_threshold(b"blosum62_20"), 100);
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

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blastsetup_unit_test.cpp:225-238 SubjectIsNucleotide
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_program.c:53-54
#[test]
fn subject_is_nucleotide_follows_ncbi_program_types() {
    for (kind, nucleotide) in [
        (RequestKind::Blastn, true),
        (RequestKind::Blastp, false),
        (RequestKind::Blastx, false),
        (RequestKind::Tblastn, true),
        (RequestKind::Tblastx, true),
        (RequestKind::Psiblast, false),
        (RequestKind::PsiTblastn, true),
        (RequestKind::Rpsblast, false),
        (RequestKind::Rpstblastn, false),
    ] {
        assert_eq!(kind.subject_is_nucleotide(), nucleotide, "{kind:?}");
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
    inputs.matrix = Some(b"BLOSUM62".to_vec());
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

/// The options that `CBlastOptionsBuilder` builds of `setters` for `program`, without the
/// request's lists and masks.
fn built(program: EProgram, setters: Vec<Setter>) -> BuiltOptions {
    BuiltOptions {
        program,
        setters,
        entrez_query: None,
        first_db_seq: None,
        final_db_seq: None,
        gi_list: None,
        negative_gi_list: None,
        taxid_list: None,
        negative_taxid_list: None,
        db_filtering_algorithm_id: None,
        db_filtering_algorithm_key: None,
        subject_masking_type: None,
        query_masks: None,
        query_range: crate::report::archive::read::options::SeqRange::EMPTY,
    }
}

// The remote options of a built handle (`RemoteOptions`: blast_formatter `-outfmt 11` and the
// import of a search strategy): the factory's, then the setters'.
#[test]
fn filter_string_expansion_matches_ncbi_rules() {
    // Nucleotide `L;m;` -> MaskAtHash TRUE, DustFiltering TRUE + defaults, RepeatFiltering FALSE.
    let opts = RemoteOptions::built(&built(
        EProgram::Blastn,
        vec![Setter::FilterString("L;m;".to_string())],
    ));
    let names: Vec<(&str, &ParamValue)> =
        opts.params().iter().map(|p| (p.name(), &p.value)).collect();
    assert_eq!(
        names,
        vec![
            ("EvalueThreshold", &ParamValue::Evalue(10.0)),
            ("DustFiltering", &ParamValue::Boolean(true)),
            ("RepeatFiltering", &ParamValue::Boolean(false)),
            ("MaskAtHash", &ParamValue::Boolean(true)),
            ("DustFilteringLevel", &ParamValue::Integer(20)),
            ("DustFilteringWindow", &ParamValue::Integer(64)),
            ("DustFilteringLinker", &ParamValue::Integer(1)),
        ]
    );

    // Custom DUST + Repeat DB on blastn.
    let opts = RemoteOptions::built(&built(
        EProgram::Blastn,
        vec![Setter::FilterString(
            "D 15 32 2; R -d custom.lib".to_string(),
        )],
    ));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "DustFilteringLevel" && p.value == ParamValue::Integer(15)));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "RepeatFiltering" && p.value == ParamValue::Boolean(true)));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "RepeatFilteringDB"
            && p.value == ParamValue::Str(b"custom.lib".to_vec())));

    // Protein `S` -> SegFiltering TRUE + defaults, no Dust/Repeat.
    let opts = RemoteOptions::built(&built(
        EProgram::Blastp,
        vec![Setter::FilterString("S".to_string())],
    ));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "SegFiltering" && p.value == ParamValue::Boolean(true)));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "SegFilteringWindow" && p.value == ParamValue::Integer(12)));
    assert!(!opts.params().iter().any(|p| p.name() == "DustFiltering"));

    // Invalid filter string falls back to cleared options.
    let opts = RemoteOptions::built(&built(
        EProgram::Blastn,
        vec![Setter::FilterString("D 15 32".to_string())],
    ));
    assert!(opts
        .params()
        .iter()
        .any(|p| p.name() == "DustFiltering" && p.value == ParamValue::Boolean(false)));
    assert!(!opts
        .params()
        .iter()
        .any(|p| p.name() == "DustFilteringLevel"));
}

#[test]
fn vecscreen_keeps_defaults_mode_stuck_except_evalue() {
    let mut opts = RemoteOptions::built(&built(
        EProgram::VecScreen,
        vec![
            Setter::WordSize(11),
            Setter::EvalueThreshold(500.0),
            Setter::HitlistSize(10),
        ],
    ));
    opts.set(fields::HITLIST_SIZE, ParamValue::Integer(25));
    assert_eq!(opts.program_and_service(), ("blastn", "vecscreen"));
    assert_eq!(
        opts.params(),
        &[Param::new(
            fields::EVALUE_THRESHOLD,
            ParamValue::Evalue(500.0)
        )]
    );
}

#[test]
fn factory_defaults_for_disc_psi_psitblastn_and_delta() {
    let disc = RemoteOptions::built(&built(EProgram::DiscMegablast, vec![]));
    let disc_names: Vec<&str> = disc.params().iter().map(Param::name).collect();
    assert_eq!(
        disc_names,
        vec![
            "EvalueThreshold",
            "MBTemplateType",
            "MBTemplateLength",
            "WordSize",
            "WindowSize",
        ]
    );

    let psi_tbn = RemoteOptions::built(&built(EProgram::PsiTblastn, vec![]));
    let psi_tbn_names: Vec<&str> = psi_tbn.params().iter().map(Param::name).collect();
    assert_eq!(
        psi_tbn_names,
        vec![
            "EvalueThreshold",
            "InclusionThreshold",
            "PseudoCountWeight",
            "SmithWatermanMode",
            "DbGeneticCode",
            "WordThreshold",
        ]
    );

    let delta = RemoteOptions::built(&built(EProgram::DeltaBlast, vec![]));
    assert_eq!(delta.program_and_service(), ("blastp", "delta_blast"));
    assert!(delta
        .params()
        .iter()
        .any(|p| p.name() == "InclusionThreshold" && p.value == ParamValue::Real(0.05)));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1352-1364,2619-2634
// `SetWindowMaskerTaxId` of 0 removes the parameter; `SetUseIndex` records `MbIndexName` only when
// the name is not empty (NCBI 2.17.0's own strategies: `F/f3/tmp/r6/f2/wm.asn`, `ui.asn`).
#[test]
fn a_taxid_of_0_and_an_index_without_a_name_record_what_ncbi_s_setters_do() {
    let options = RemoteOptions::built(&built(
        EProgram::Megablast,
        vec![
            Setter::WindowMaskerTaxId(9606),
            Setter::WindowMaskerTaxId(0),
            Setter::UseIndex {
                index_name: String::new(),
            },
        ],
    ));
    assert_eq!(
        options.params(),
        &[
            Param::new(fields::EVALUE_THRESHOLD, ParamValue::Evalue(10.0)),
            Param::new(fields::FORCE_MB_INDEX, ParamValue::Boolean(true)),
        ]
    );
    let named = RemoteOptions::built(&built(
        EProgram::Megablast,
        vec![Setter::UseIndex {
            index_name: "idx".to_string(),
        }],
    ));
    assert_eq!(
        named.params()[1..],
        [
            Param::new(fields::FORCE_MB_INDEX, ParamValue::Boolean(true)),
            Param::new(fields::MB_INDEX_NAME, ParamValue::Str(b"idx".to_vec())),
        ]
    );
}

// A `FilterString` number that `sscanf` does not store is undecided where the remote options
// record it (review R6 M1: `S 10 abc 2.5`, `D abc 64 1`), and decided where they do not.
#[test]
fn an_unread_filter_number_is_undecided_where_it_is_recorded() {
    let filter = |program, text: &str| {
        RemoteOptions::built(&built(
            program,
            vec![Setter::FilterString(text.to_string())],
        ))
    };
    assert_eq!(
        filter(EProgram::Blastp, "S 10 abc 2.5;").undecided(),
        Some(UNDECIDED_FILTER_VALUE)
    );
    assert_eq!(
        filter(EProgram::Blastn, "D abc 64 1;").undecided(),
        Some(UNDECIDED_FILTER_VALUE)
    );
    // SEG's values are not recorded for a nucleotide query, DUST's not for a protein one.
    assert_eq!(filter(EProgram::Blastn, "S 10 abc 2.5;").undecided(), None);
    assert_eq!(filter(EProgram::Blastp, "D abc 64 1;").undecided(), None);
    let level = filter(EProgram::Blastn, "D 99999999999999999999 64 1;");
    assert!(level.params().contains(&Param::new(
        fields::DUST_FILTERING_LEVEL,
        ParamValue::Integer(-1)
    )));
}

// The search programs' handles are the factory's (`RequestKind::factory_program`): dc-megablast's
// four parameters, PSI-TBLASTN's, and the request's program and service.
#[test]
fn the_search_programs_handles_are_the_factory_s() {
    let dc = OptionInputs::new(Program::Blastn, "dc-megablast");
    assert_eq!(
        names(&dc)[..5],
        [
            "EvalueThreshold",
            "MBTemplateType",
            "MBTemplateLength",
            "WordSize",
            "WindowSize"
        ]
    );
    assert_eq!(service(Program::Blastn, "dc-megablast"), "megablast");
    assert_eq!(service(Program::Blastn, "blastn-short"), "plain");
    assert_eq!(service(RequestKind::PsiTblastn, "tblastn"), "psi");
    assert_eq!(service(RequestKind::Rpstblastn, "rpstblastn"), "rpsblast");
    assert_eq!(RequestKind::PsiTblastn.program(), "blastp");
    assert_eq!(RequestKind::Rpstblastn.program(), "tblastn");
}
