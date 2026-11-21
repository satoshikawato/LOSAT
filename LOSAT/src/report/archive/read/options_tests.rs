//! The builder on the crafted requests of line RG F2's probes (`E.*` of `F/f2/SURVEY_F2.md`
//! section 2.2: NCBI BLAST+ 2.17.0's own blastn and blastp exports with one change each), with
//! the values the generic reader gives. The texts and exit statuses are NCBI's: rc 3 is the
//! `CBlastException` clause of `CATCH_ALL` (`BLAST engine error: <GetMsg()>`), rc 255 the
//! `CException` one (`Error: <what()>`).

use super::*;
use crate::ncbi_exception::NCBI_BUILD_ROOT;

fn blast4(value_type: ValueType, value: Value) -> Value {
    Value::Choice(u16::from(value_type.tag()), Box::new(value))
}

fn int(n: i64) -> Value {
    blast4(ValueType::Integer, Value::Int(n))
}

fn big(n: i64) -> Value {
    blast4(ValueType::BigInteger, Value::Int(n))
}

fn boolean(b: bool) -> Value {
    blast4(ValueType::Boolean, Value::Bool(b))
}

fn real(r: f64) -> Value {
    blast4(ValueType::Real, Value::Real(r))
}

fn string(s: &str) -> Value {
    blast4(ValueType::String, Value::Str(s.as_bytes().to_vec()))
}

fn evalue(e: f64) -> Value {
    blast4(
        ValueType::Cutoff,
        Value::Choice(0, Box::new(Value::Real(e))),
    )
}

fn raw_score(score: i64) -> Value {
    blast4(
        ValueType::Cutoff,
        Value::Choice(1, Box::new(Value::Int(score))),
    )
}

fn int_list(numbers: &[i32]) -> Value {
    blast4(ValueType::IntegerList, Value::Int4s(numbers.to_vec()))
}

fn strand(n: i32) -> Value {
    blast4(ValueType::StrandType, Value::Enum(n))
}

fn mask(tag: i64) -> Value {
    blast4(
        ValueType::QueryMask,
        Value::Sequence(vec![(1, Value::Int(tag))]),
    )
}

type List = Vec<(&'static str, Value)>;

fn raw<'a>(list: &'a [(&'static str, Value)]) -> Vec<RawParam<'a>> {
    list.iter()
        .map(|(name, value)| RawParam {
            name: name.as_bytes(),
            value,
        })
        .collect()
}

/// The algorithm options of NCBI's blastn export (`megablast` service) after `head`.
fn blastn_options(head: List) -> List {
    let mut list = head;
    list.extend([
        ("EvalueThreshold", evalue(10.0)),
        ("MaskAtHash", boolean(true)),
        ("UngappedMode", boolean(false)),
        ("HitlistSize", int(500)),
    ]);
    list
}

/// The algorithm options of NCBI's blastp export (`plain` service) after `head`.
fn blastp_options(head: List) -> List {
    let mut list = head;
    list.extend([
        ("EvalueThreshold", evalue(10.0)),
        ("MaskAtHash", boolean(false)),
        ("SegFiltering", boolean(false)),
        ("WordThreshold", int(11)),
        ("WindowSize", int(40)),
        ("HitlistSize", int(500)),
        ("UngappedMode", boolean(false)),
        ("CompositionBasedStats", int(2)),
        ("SmithWatermanMode", boolean(false)),
    ]);
    list
}

struct Request {
    program: &'static str,
    service: &'static str,
    algorithm: List,
    program_options: Option<List>,
    format: Option<List>,
}

impl Request {
    fn blastn(algorithm: List) -> Request {
        Request {
            program: "blastn",
            service: "megablast",
            algorithm,
            program_options: None,
            format: None,
        }
    }

    fn blastp(algorithm: List) -> Request {
        Request {
            program: "blastp",
            service: "plain",
            algorithm,
            program_options: None,
            format: None,
        }
    }

    fn with_program(mut self, list: List) -> Request {
        self.program_options = Some(list);
        self
    }

    fn with_format(mut self, list: List) -> Request {
        self.format = Some(list);
        self
    }

    fn import(&self) -> Result<ImportedOptions, NcbiException> {
        let algorithm = raw(&self.algorithm);
        let program = self.program_options.as_deref().map(raw);
        let format = self.format.as_deref().map(raw);
        import_options(
            self.program.as_bytes(),
            self.service.as_bytes(),
            OptionLists {
                algorithm: Some(&algorithm),
                program: program.as_deref(),
                format: format.as_deref(),
            },
        )
    }
}

fn engine_error(result: Result<ImportedOptions, NcbiException>) -> (&'static str, String) {
    let error = result.expect_err("NCBI's builder fails");
    (
        error.kind(),
        String::from_utf8_lossy(error.get_msg()).into_owned(),
    )
}

#[test]
fn ncbis_own_exports_build_their_programs_options() {
    let blastn = Request::blastn(blastn_options(Vec::new()))
        .import()
        .unwrap();
    assert_eq!(blastn.options.program, EProgram::Megablast);
    assert_eq!(blastn.options.program.task_name(), "megablast");
    assert_eq!(
        blastn.options.setters,
        [
            Setter::EvalueThreshold(10.0),
            Setter::MaskAtHash(true),
            Setter::GappedMode(true),
            Setter::HitlistSize(500),
        ]
    );
    assert_eq!(blastn.options.query_range, SeqRange::EMPTY);
    assert_eq!(blastn.psi_num_of_iterations, 0);
    assert_eq!(blastn.filtering_id, -1);
    assert_eq!(blastn.filtering_key, "");
    assert_eq!(blastn.subject_masking_type, 0);

    let blastp = Request::blastp(blastp_options(Vec::new()))
        .import()
        .unwrap();
    assert_eq!(blastp.options.program, EProgram::Blastp);
    assert_eq!(blastp.options.program.task_name(), "blastp");
    assert_eq!(
        blastp.options.setters,
        [
            Setter::EvalueThreshold(10.0),
            Setter::MaskAtHash(false),
            Setter::SegFiltering(false),
            Setter::WordThreshold(11),
            Setter::WindowSize(40),
            Setter::HitlistSize(500),
            Setter::GappedMode(true),
            Setter::CompositionBasedStats(2),
            Setter::SmithWatermanMode(false),
        ]
    );
}

/// `E.blastn.badparam`, `E.blastp.badparam`: rc 3 `BLAST engine error: Internal: Error processing
/// option [FooBar] type [5].`
#[test]
fn an_unknown_parameter_is_refused_with_its_type() {
    for request in [
        Request::blastn(blastn_options(vec![("FooBar", int(1))])),
        Request::blastp(blastp_options(vec![("FooBar", int(1))])),
    ] {
        assert_eq!(
            engine_error(request.import()),
            (
                "CRemoteBlastException::eServiceNotAvailable",
                "Internal: Error processing option [FooBar] type [5].".to_string()
            )
        );
    }
    // A known name of another type, a known name the switch has no arm for (NCBI's own export of
    // `-window_masker_db`, `I.blastn.wmdb.alone`), a name of no first letter the switch knows.
    for (name, value, which) in [
        ("WordSize", real(11.0), 7),
        ("Culling", string("x"), 12),
        ("WindowMaskerDatabase", string("wm.stats"), 12),
        ("Web_JobTitle", string("x"), 12),
        ("wordsize", int(11), 5),
    ] {
        let request = Request::blastn(blastn_options(vec![(name, value)]));
        assert_eq!(
            engine_error(request.import()).1,
            format!("Internal: Error processing option [{name}] type [{which}].")
        );
    }
    let request = Request::blastn(blastn_options(vec![("", int(1))]));
    assert_eq!(
        engine_error(request.import()),
        (
            "CBlastException::eInvalidArgument",
            "Option has no name.".to_string()
        )
    );
}

/// `E.blastn.badtype_evalue`: rc 3 `BLAST engine error: EvalueThreshold has unsupported type.`;
/// `E.*.evalue_cutoff`, `E.*.dupparam`: a real is read, the last value wins.
#[test]
fn evalue_threshold_is_taken_by_its_name() {
    let request = Request::blastn(blastn_options(vec![("EvalueThreshold", int(5))]));
    assert_eq!(
        engine_error(request.import()),
        (
            "CBlastException::eInvalidArgument",
            "EvalueThreshold has unsupported type.".to_string()
        )
    );
    let request = Request::blastn(blastn_options(vec![("EvalueThreshold", raw_score(5))]));
    assert_eq!(
        engine_error(request.import()).1,
        "EvalueThreshold has unsupported type."
    );
    let as_real = Request::blastn(blastn_options(vec![("EvalueThreshold", real(1e-30))]))
        .import()
        .unwrap();
    assert_eq!(
        as_real.options.setters[..2],
        [
            Setter::EvalueThreshold(1e-30),
            Setter::EvalueThreshold(10.0)
        ]
    );
    // The list as NCBI's export writes it with `-evalue 1e-30` and a second one before it.
    let mut list = vec![("EvalueThreshold", evalue(1e-30))];
    list.extend(blastn_options(Vec::new()));
    list.swap(0, 1);
    let duplicate = Request::blastn(list).import().unwrap();
    assert_eq!(
        duplicate.options.setters[..2],
        [
            Setter::EvalueThreshold(10.0),
            Setter::EvalueThreshold(1e-30)
        ]
    );
}

/// `E.blastn.algo_in_prog`: the program options come after the algorithm options.
#[test]
fn program_options_follow_the_algorithm_options() {
    let built = Request::blastn(blastn_options(Vec::new()))
        .with_program(vec![("EvalueThreshold", evalue(1e-30))])
        .import()
        .unwrap();
    assert_eq!(
        built.options.setters.last(),
        Some(&Setter::EvalueThreshold(1e-30))
    );
}

/// `E.blastn.cutoff_score`, `E.blastp.cutoff_score`: rc 255, `Error: ` and the exception's report
/// (NCBI's official 2.17.0 Linux build: its build tree).
#[test]
fn cutoff_score_fails_in_ncbis_getter() {
    for request in [
        Request::blastn(blastn_options(vec![("CutoffScore", raw_score(30))])),
        Request::blastp(blastp_options(vec![("CutoffScore", raw_score(30))])),
    ] {
        let error = request.import().unwrap_err();
        assert_eq!(
            String::from_utf8_lossy(&error.what()),
            crate::wording::losat_text(&format!(
                "NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/objects/blast/Blast4_value_.cpp\", \
                 line 252: Error: (CInvalidChoiceSelection::eFail) CBlast4_value::GetInteger(): \
                 Invalid choice selection: NCBI-Blast4::Blast4-value.cutoff\n"
            ))
        );
    }
    // Of another type, the name is not the field's.
    let request = Request::blastn(blastn_options(vec![("CutoffScore", int(30))]));
    assert_eq!(
        engine_error(request.import()).1,
        "Internal: Error processing option [CutoffScore] type [5]."
    );
}

/// `E.blastn.badprog`, `E.blastp.badprog`, `E.*.badservice`.
#[test]
fn the_program_and_service_are_ncbis_combinations() {
    let mut request = Request::blastn(blastn_options(Vec::new()));
    request.program = "blastz";
    assert_eq!(
        engine_error(request.import()),
        (
            "CBlastException::eInvalidArgument",
            "Unsupported combination of program (blastz) and service (megablast).".to_string()
        )
    );
    let mut request = Request::blastp(blastp_options(Vec::new()));
    request.program = "blastz";
    assert_eq!(
        engine_error(request.import()),
        (
            "CBlastException::eNotSupported",
            "Program type 'blastz' not supported".to_string()
        )
    );
    for (program, service) in [("blastn", "foo"), ("blastp", "foo")] {
        let mut request = Request::blastn(blastn_options(Vec::new()));
        request.program = program;
        request.service = service;
        assert_eq!(
            engine_error(request.import()).1,
            format!("Unsupported combination of program ({program}) and service ({service}).")
        );
    }
    // The texts keep the request's case; the program is found in lower case.
    assert_eq!(
        compute_program(b"BlastN", b"MegaBlast").unwrap(),
        EProgram::Megablast
    );
    assert_eq!(
        compute_program(b"BLASTZ", b"Foo").unwrap_err().get_msg(),
        b"Unsupported combination of program (BLASTZ) and service (Foo)."
    );
    assert_eq!(
        compute_program(b"BLASTZ", b"plain").unwrap_err().get_msg(),
        b"Program type 'blastz' not supported"
    );
    for (program, service, expected) in [
        ("blastp", "rpsblast", EProgram::RpsBlast),
        ("blastp", "psi", EProgram::PsiBlast),
        ("blastp", "phi", EProgram::Blastp),
        ("blastp", "delta_blast", EProgram::DeltaBlast),
        ("blastn", "vecscreen", EProgram::VecScreen),
        ("blastn", "plain", EProgram::Blastn),
        ("tblastn", "rpsblast", EProgram::RpsTblastn),
        ("tblastn", "psi", EProgram::PsiTblastn),
        ("tblastx", "anything", EProgram::Tblastx),
        ("blastx", "rpsblast", EProgram::RpsTblastn),
        ("blastx", "plain", EProgram::Blastx),
        ("blastn-short", "multi_blast", EProgram::Blastn),
        ("rmblastn", "plain", EProgram::Blastn),
        ("dc-megablast", "plain", EProgram::DiscMegablast),
        ("mapr2r", "plain", EProgram::Mapper),
    ] {
        assert_eq!(
            compute_program(program.as_bytes(), service.as_bytes()).unwrap(),
            expected,
            "{program} {service}"
        );
    }
}

/// `E.blastp.phi`: phiblastp with the pattern; `E.blastn.phi`: rc 3 `BLAST engine error: Incorrect
/// combination of option (PHIPattern) and program (blastn)` (the export's service, megablast, makes
/// a program that has no PHI-BLAST).
#[test]
fn a_phi_pattern_makes_phi_blast() {
    let built = Request::blastp(blastp_options(vec![("PHIPattern", string("C-x(2)-C"))]))
        .import()
        .unwrap();
    assert_eq!(built.options.program, EProgram::PhiBlastp);
    assert_eq!(built.options.program.task_name(), "phiblastp");
    assert_eq!(
        built.options.setters[0],
        Setter::PhiPattern {
            pattern: "C-x(2)-C".to_string(),
            is_na: false
        }
    );
    // The program of the message is the request's.
    let request = Request::blastn(blastn_options(vec![("PHIPattern", string("C-x(2)-C"))]));
    assert_eq!(
        engine_error(request.import()),
        (
            "CRemoteBlastException::eServiceNotAvailable",
            "Incorrect combination of option (PHIPattern) and program (blastn)".to_string()
        )
    );
    let mut request = Request::blastn(blastn_options(vec![("PHIPattern", string("ACGT"))]));
    request.service = "plain";
    let built = request.import().unwrap();
    assert_eq!(built.options.program, EProgram::PhiBlastn);
    assert_eq!(
        built.options.setters[0],
        Setter::PhiPattern {
            pattern: "ACGT".to_string(),
            is_na: true
        }
    );
    // An empty pattern makes PHI-BLAST and sets nothing.
    let built = Request::blastp(blastp_options(vec![("PHIPattern", string(""))]))
        .import()
        .unwrap();
    assert_eq!(built.options.program, EProgram::PhiBlastp);
    assert_eq!(built.options.setters[0], Setter::EvalueThreshold(10.0));
}

/// `E.blastp.web_step3`, `E.blastn.web_step3`: psiblast from the format options, 3 iterations.
#[test]
fn web_step_number_makes_psiblast_and_its_iterations() {
    for request in [
        Request::blastp(blastp_options(Vec::new())),
        Request::blastn(blastn_options(Vec::new())),
    ] {
        let built = request
            .with_format(vec![("Web_StepNumber", int(3))])
            .import()
            .unwrap();
        assert_eq!(built.options.program, EProgram::PsiBlast);
        assert_eq!(built.options.program.task_name(), "psiblast");
        assert_eq!(built.psi_num_of_iterations, 3);
    }
    // The web's string form, a string that is no number (0), a negative count (`unsigned int`),
    // and the first of two.
    for (value, expected) in [(string("4"), 4), (string("x4"), 0), (int(-1), u32::MAX)] {
        let built = Request::blastp(blastp_options(Vec::new()))
            .with_format(vec![("Web_StepNumber", value), ("Web_StepNumber", int(9))])
            .import()
            .unwrap();
        assert_eq!(built.psi_num_of_iterations, expected);
    }
    // A string makes no psiblast (the field is an integer); another type fails in `GetString()`.
    let built = Request::blastp(blastp_options(Vec::new()))
        .with_format(vec![("Web_StepNumber", string("2"))])
        .import()
        .unwrap();
    assert_eq!(built.options.program, EProgram::Blastp);
    let error = Request::blastp(blastp_options(Vec::new()))
        .with_format(vec![("Web_StepNumber", boolean(true))])
        .import()
        .unwrap_err();
    assert_eq!(error.kind(), "CInvalidChoiceSelection::eFail");
    assert_eq!(
        error.get_msg(),
        b"CBlast4_value::GetString(): Invalid choice selection: NCBI-Blast4::Blast4-value.boolean"
    );
    // The other PSI-BLAST parameters, and the first that tells: dc-megablast before psiblast.
    for name in ["Web_RunPsiBlast", "IgnoreMsaMaster"] {
        let built = Request::blastn(blastn_options(Vec::new()))
            .with_format(vec![(name, boolean(true))])
            .import()
            .unwrap();
        assert_eq!(built.options.program, EProgram::PsiBlast, "{name}");
    }
    let built = Request::blastn(blastn_options(vec![
        ("MBTemplateLength", int(0)),
        ("PseudoCountWeight", int(0)),
        ("MBTemplateLength", int(18)),
    ]))
    .import()
    .unwrap();
    assert_eq!(built.options.program, EProgram::PsiBlast);
    let built = Request::blastn(blastn_options(vec![
        ("MBTemplateLength", int(274)),
        ("PseudoCountWeight", int(0)),
    ]))
    .import()
    .unwrap();
    assert_eq!(built.options.program, EProgram::DiscMegablast);
    assert_eq!(built.options.program.task_name(), "dc-megablast");
    assert_eq!(built.options.setters[0], Setter::MbTemplateLength(18));
}

/// `E.blastn.qrange_out`, `E.blastn.required_start_only`: the range as set, its end open.
#[test]
fn the_query_range_is_kept_as_set() {
    let built = Request::blastn(blastn_options(Vec::new()))
        .with_program(vec![
            ("RequiredStart", int(5000)),
            ("RequiredEnd", int(6000)),
        ])
        .import()
        .unwrap();
    assert_eq!(
        built.options.query_range,
        SeqRange {
            from: 5000,
            to_open: 6000
        }
    );
    let built = Request::blastn(blastn_options(Vec::new()))
        .with_program(vec![("RequiredStart", int(10))])
        .import()
        .unwrap();
    assert_eq!(
        built.options.query_range,
        SeqRange {
            from: 10,
            to_open: u32::MAX
        }
    );
}

/// `E.*.taxid_empty`, `E.*.neg_taxid`, `E.*.gilist`, `E.*.entrez`, `E.*.firstdbseq`: kept for the
/// database (what NCBI's import then does with them is the import's).
#[test]
fn database_parameters_are_kept() {
    let built = Request::blastn(blastn_options(Vec::new()))
        .with_program(vec![
            ("TaxidList", int_list(&[])),
            ("NegativeTaxidList", int_list(&[9606])),
            ("GiList", int_list(&[1, 2, -1])),
            ("NegativeGiList", int_list(&[3])),
            ("EntrezQuery", string("human[orgn]")),
            ("FirstDbSeq", int(2)),
            ("FinalDbSeq", int(3)),
            ("DbFilteringAlgorithmKey", string("dust")),
            ("DbFilteringAlgorithmId", int(11)),
            ("SubjectMaskingType", int(2)),
        ])
        .import()
        .unwrap();
    let options = &built.options;
    assert_eq!(options.taxid_list.as_deref(), Some(&[][..]));
    assert_eq!(options.negative_taxid_list.as_deref(), Some(&[9606][..]));
    assert_eq!(options.gi_list.as_deref(), Some(&[1, 2, 4_294_967_295][..]));
    assert_eq!(options.negative_gi_list.as_deref(), Some(&[3][..]));
    assert_eq!(options.entrez_query.as_deref(), Some("human[orgn]"));
    assert_eq!(
        (options.first_db_seq, options.final_db_seq),
        (Some(2), Some(3))
    );
    assert_eq!(
        (built.filtering_id, built.filtering_key.as_str()),
        (11, "dust")
    );
    assert_eq!(built.subject_masking_type, 2);
    assert_eq!(options.setters.len(), 4);
}

/// `E.*.hsprange_culling`: `Culling TRUE` with `HspRangeMax 1` is the culling limit 1, after the
/// other calls; the megablast index of `ForceMbIndex`.
#[test]
fn interactions_come_last() {
    let built = Request::blastn(blastn_options(vec![
        ("Culling", boolean(true)),
        ("HspRangeMax", int(1)),
        ("Culling", int(7)),
        ("MbIndexName", string("idx")),
        ("ForceMbIndex", boolean(true)),
    ]))
    .import()
    .unwrap();
    assert_eq!(built.options.setters[0], Setter::CullingLimit(7));
    assert_eq!(
        built.options.setters[built.options.setters.len() - 2..],
        [
            Setter::CullingLimit(1),
            Setter::UseIndex {
                index_name: "idx".to_string()
            }
        ]
    );
    let built = Request::blastn(blastn_options(vec![
        ("Culling", boolean(true)),
        ("Culling", boolean(false)),
        ("ForceMbIndex", boolean(false)),
    ]))
    .import()
    .unwrap();
    assert_eq!(built.options.setters.len(), 4);
}

/// The query masks of the program options are dropped when the algorithm options had one.
#[test]
fn query_masks_come_from_one_list() {
    let built = Request::blastn(blastn_options(vec![("LCaseMask", mask(1))]))
        .with_program(vec![("LCaseMask", mask(2))])
        .import()
        .unwrap();
    assert_eq!(
        built.options.query_masks,
        Some(vec![Value::Sequence(vec![(1, Value::Int(1))])])
    );
    let built = Request::blastn(blastn_options(Vec::new()))
        .with_program(vec![("LCaseMask", mask(2)), ("LCaseMask", mask(3))])
        .import()
        .unwrap();
    assert_eq!(built.options.query_masks.map(|masks| masks.len()), Some(2));
}

/// Every arm's setter, with NCBI's conversions.
#[test]
fn every_setter_of_the_switch() {
    let built = Request::blastn(vec![
        ("BestHitScoreEdge", real(0.1)),
        ("BestHitOverhang", real(0.2)),
        ("CompositionBasedStats", int(3)),
        ("DbGeneticCode", int(4)),
        ("DbLength", big(5_000_000_000)),
        ("DustFiltering", boolean(true)),
        ("DustFilteringLevel", int(20)),
        ("DustFilteringWindow", int(64)),
        ("DustFilteringLinker", int(1)),
        ("DomainInclusionThreshold", real(0.05)),
        ("EffectiveSearchSpace", big(7)),
        ("FilterString", string("L;m;")),
        ("GapExtensionCost", int(2)),
        ("GapOpeningCost", int(5)),
        ("GapTracebackAlgorithm", int(1)),
        ("GapTrigger", real(22.0)),
        ("GapXDropoff", real(30.0)),
        ("GapXDropoffFinal", real(100.0)),
        ("GapExtnAlgorithm", int(0)),
        ("HitlistSize", int(10)),
        ("InclusionThreshold", real(0.002)),
        ("LongestIntronLength", int(0)),
        ("MBTemplateType", int(257)),
        ("MatchReward", int(1)),
        ("MatrixName", string("BLOSUM80")),
        (
            "MatrixTable",
            blast4(ValueType::Matrix, Value::Sequence(Vec::new())),
        ),
        ("MismatchPenalty", int(-2)),
        ("MaxHspsPerSubject", int(3)),
        ("OutOfFrameMode", boolean(false)),
        ("PercentIdentity", real(90.0)),
        ("QueryGeneticCode", int(2)),
        ("QueryCovHspPerc", real(50.0)),
        ("RepeatFiltering", boolean(true)),
        ("RepeatFilteringDB", string("repeats")),
        ("StrandOption", strand(2)),
        ("SegFiltering", boolean(true)),
        ("SegFilteringWindow", int(12)),
        ("SegFilteringLocut", real(2.2)),
        ("SegFilteringHicut", real(2.5)),
        ("SumStatistics", boolean(true)),
        ("SmithWatermanMode", boolean(true)),
        ("UngappedMode", boolean(true)),
        ("UnifiedP", int(1)),
        ("UseSubjectBestHit", boolean(false)),
        ("WindowMaskerTaxId", int(9606)),
        ("WindowSize", int(0)),
        ("WordSize", int(28)),
        ("WordThreshold", int(0)),
    ])
    .import()
    .unwrap();
    assert_eq!(
        built.options.setters,
        [
            Setter::BestHitScoreEdge(0.1),
            Setter::BestHitOverhang(0.2),
            Setter::CompositionBasedStats(3),
            Setter::DbGeneticCode(4),
            Setter::DbLength(5_000_000_000),
            Setter::DustFiltering(true),
            Setter::DustFilteringLevel(20),
            Setter::DustFilteringWindow(64),
            Setter::DustFilteringLinker(1),
            Setter::DomainInclusionThreshold(0.05),
            Setter::EffectiveSearchSpace(7),
            Setter::FilterString("L;m;".to_string()),
            Setter::GapExtensionCost(2),
            Setter::GapOpeningCost(5),
            Setter::GapTracebackAlgorithm(1),
            Setter::GapTrigger(22.0),
            Setter::GapXDropoff(30.0),
            Setter::GapXDropoffFinal(100.0),
            Setter::GapExtnAlgorithm(0),
            Setter::HitlistSize(10),
            Setter::InclusionThreshold(0.002),
            Setter::LongestIntronLength(0),
            Setter::MbTemplateType(1),
            Setter::MatchReward(1),
            Setter::MatrixName("BLOSUM80".to_string()),
            Setter::MismatchPenalty(-2),
            Setter::MaxHspsPerSubject(3),
            Setter::OutOfFrameMode(false),
            Setter::PercentIdentity(90.0),
            Setter::QueryGeneticCode(2),
            Setter::QueryCovHspPerc(50.0),
            Setter::RepeatFiltering(true),
            Setter::RepeatFilteringDb("repeats".to_string()),
            Setter::StrandOption(2),
            Setter::SegFiltering(true),
            Setter::SegFilteringWindow(12),
            Setter::SegFilteringLocut(2.2),
            Setter::SegFilteringHicut(2.5),
            Setter::SumStatisticsMode(true),
            Setter::SmithWatermanMode(true),
            Setter::GappedMode(false),
            Setter::UnifiedP(1),
            Setter::SubjectBestHit,
            Setter::WindowMaskerTaxId(9606),
            Setter::WindowSize(0),
            Setter::WordSize(28),
            Setter::WordThreshold(0),
        ]
    );
    // psiblast's own parameters (`PseudoCountWeight` and `IgnoreMsaMaster` make psiblast).
    let built = Request::blastp(vec![
        ("PseudoCountWeight", int(7)),
        ("IgnoreMsaMaster", boolean(true)),
    ])
    .import()
    .unwrap();
    assert_eq!(built.options.program, EProgram::PsiBlast);
    assert_eq!(
        built.options.setters,
        [Setter::PseudoCount(7), Setter::IgnoreMsaMaster(true)]
    );
}

/// `Web_StepNumber` is `FetchData`'s, not the builder's.
#[test]
fn the_step_number_is_the_fetchs() {
    let step = vec![("Web_StepNumber", string("x"))];
    let format = raw(&step);
    let built = build(
        b"blastp",
        b"plain",
        OptionLists {
            algorithm: None,
            program: None,
            format: Some(&format),
        },
    )
    .unwrap();
    assert_eq!(built.program, EProgram::Blastp);
}
