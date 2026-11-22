//! The options of a stored request (a search strategy, `-import_search_strategy`, or the request
//! of a BLAST archive, `blast_formatter`): NCBI's `CBlastOptionsBuilder` and the options part of
//! `CImportStrategy::FetchData`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:79-129
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:667-695
//! [`build`] is that function: the program, its task, and the `Set*` calls that the builder makes
//! on the program's factory options (`CBlastOptionsFactory::Create(program)`), in NCBI's order
//! ([`Setter`]); with what the builder keeps for the database and the query
//! ([`BuiltOptions`]). The caller applies the calls in order to the program's own options, over
//! the defaults of the task: the last call of an option wins, as NCBI's do. [`import_options`] is
//! `FetchData`'s part (the PSI-BLAST iteration count from the format options first).
//!
//! The builder sees each parameter as the reader gives it ([`RawParam`]): its name and its
//! `Blast4-value` (`Value::Choice`, the alternative's tag and its value). A parameter is applied
//! only when its name and value type are a field's (`fields::Field::matches`).

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::archive::fields::{self, Field, ValueType};
use crate::report::serial::asn_types::Value;

use super::request::{ParamSet, Request};

/// NCBI's `EProgram` (blast_types.hpp), the values that the builder can give.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EProgram {
    Blastn,
    Megablast,
    DiscMegablast,
    Blastp,
    Blastx,
    Tblastn,
    Tblastx,
    RpsBlast,
    RpsTblastn,
    PsiBlast,
    PsiTblastn,
    PhiBlastp,
    PhiBlastn,
    DeltaBlast,
    VecScreen,
    Mapper,
}

impl EProgram {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:676-699
    /// `EProgramToTaskName`: the task of the options (the strategy's task).
    pub fn task_name(self) -> &'static str {
        match self {
            EProgram::Blastn => "blastn",
            EProgram::Megablast => "megablast",
            EProgram::DiscMegablast => "dc-megablast",
            EProgram::Blastp => "blastp",
            EProgram::Blastx => "blastx",
            EProgram::Tblastn => "tblastn",
            EProgram::Tblastx => "tblastx",
            EProgram::RpsBlast => "rpsblast",
            EProgram::RpsTblastn => "rpstblastn",
            EProgram::PsiBlast => "psiblast",
            EProgram::PsiTblastn => "psitblastn",
            EProgram::PhiBlastp => "phiblastp",
            EProgram::PhiBlastn => "phiblastn",
            EProgram::DeltaBlast => "deltablast",
            EProgram::VecScreen => "vecscreen",
            EProgram::Mapper => "mapr2g",
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_program.c:43-44
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:708-755
    /// `Blast_QueryIsNucleotide` of the options' program type (`EProgramToEBlastProgramType`): the
    /// nucleotide searches and the translated queries.
    pub fn query_is_nucleotide(self) -> bool {
        matches!(
            self,
            EProgram::Blastn
                | EProgram::Megablast
                | EProgram::DiscMegablast
                | EProgram::VecScreen
                | EProgram::Mapper
                | EProgram::Blastx
                | EProgram::Tblastx
                | EProgram::RpsTblastn
                | EProgram::PhiBlastn
        )
    }
}

/// A `Blast4-parameter` as the reader gives it: the name (a VisibleString: the reader keeps bytes
/// 0x20..0x7E only, `asn_types::read`) and the `Blast4-value` (`Value::Choice`).
#[derive(Clone, Copy, Debug)]
pub struct RawParam<'a> {
    pub name: &'a [u8],
    pub value: &'a Value,
}

impl RawParam<'_> {
    /// `CBlast4_value::Which()`.
    fn value_type(&self) -> ValueType {
        blast4_value(self.value).0
    }

    /// `CBlast4Field::Match`.
    fn is(&self, field: Field) -> bool {
        field.matches(self.name, self.value_type())
    }

    /// The name as NCBI's messages print it.
    fn name_text(&self) -> String {
        String::from_utf8_lossy(self.name).into_owned()
    }
}

/// The alternative and value of a `Blast4-value`.
fn blast4_value(value: &Value) -> (ValueType, &Value) {
    match value {
        Value::Choice(tag, value) => match ValueType::of_tag(*tag) {
            Some(value_type) => (value_type, value),
            None => panic!("Blast4-value has no alternative {tag}"),
        },
        _ => panic!("a Blast4-value is read as a CHOICE"),
    }
}

// The report has no function (NCBI BLAST+ 2.17.0 `E.blastn.cutoff_score`, line RG F2 probes: the
// generated Blast4-value accessor fails with `CInvalidChoiceSelection::eFail`, an invalid choice
// selection, `NCBI-Blast4::Blast4-value.cutoff`).
const BLAST4_VALUE_SELECTION: ThrowSite = ThrowSite {
    file: "src/objects/blast/Blast4_value_.cpp",
    line: 252,
    kind: "CInvalidChoiceSelection::eFail",
    function: None,
};

// NCBI reference (598d8ae6): c++/src/serial/exception.cpp:134-141
/// `CBlast4_value::Get<requested>()` of a value of another alternative: NCBI's
/// `CInvalidChoiceSelection` (`SPrintIdentifier` capitalizes the name and writes `-` as `_`:
/// `GetBig_integer`).
fn invalid_selection(current: ValueType, requested: ValueType) -> NcbiException {
    let mut getter = String::new();
    for (i, c) in requested.variant().chars().enumerate() {
        getter.push(match c {
            '-' => '_',
            c if i == 0 => c.to_ascii_uppercase(),
            c => c,
        });
    }
    NcbiException::new(
        BLAST4_VALUE_SELECTION,
        format!(
            "CBlast4_value::Get{getter}(): Invalid choice selection: NCBI-Blast4::Blast4-value.{}",
            current.variant()
        ),
    )
}

/// The `Blast4-value` getters (`CBlast4_value::Get*`): the value of the alternative asked for, or
/// NCBI's `CInvalidChoiceSelection`. The value's own shape is the table's (an `integer` is read
/// as `Value::Int` of `int` storage, ...).
struct Getter<'a>(ValueType, &'a Value);

impl<'a> Getter<'a> {
    fn of(value: &'a Value) -> Self {
        let (value_type, value) = blast4_value(value);
        Getter(value_type, value)
    }

    fn check(&self, requested: ValueType) -> Result<&'a Value, NcbiException> {
        if self.0 == requested {
            Ok(self.1)
        } else {
            Err(invalid_selection(self.0, requested))
        }
    }

    fn integer(&self) -> Result<i32, NcbiException> {
        match self.check(ValueType::Integer)? {
            // `int` storage: the reader refuses what does not fit (`ReadInt4`).
            Value::Int(n) => Ok(i32::try_from(*n).expect("an integer of int storage")),
            _ => panic!("a Blast4-value integer is read as an INTEGER"),
        }
    }

    fn big_integer(&self) -> Result<i64, NcbiException> {
        match self.check(ValueType::BigInteger)? {
            Value::Int(n) => Ok(*n),
            _ => panic!("a Blast4-value big-integer is read as an INTEGER"),
        }
    }

    fn boolean(&self) -> Result<bool, NcbiException> {
        match self.check(ValueType::Boolean)? {
            Value::Bool(b) => Ok(*b),
            _ => panic!("a Blast4-value boolean is read as a BOOLEAN"),
        }
    }

    fn real(&self) -> Result<f64, NcbiException> {
        match self.check(ValueType::Real)? {
            Value::Real(r) => Ok(*r),
            _ => panic!("a Blast4-value real is read as a REAL"),
        }
    }

    fn string(&self) -> Result<&'a [u8], NcbiException> {
        match self.check(ValueType::String)? {
            Value::Str(text) => Ok(text),
            _ => panic!("a Blast4-value string is read as a VisibleString"),
        }
    }

    fn strand_type(&self) -> Result<i32, NcbiException> {
        match self.check(ValueType::StrandType)? {
            Value::Enum(n) => Ok(*n),
            _ => panic!("a Blast4-value strand-type is read as an ENUMERATED"),
        }
    }

    fn integer_list(&self) -> Result<&'a [i32], NcbiException> {
        match self.check(ValueType::IntegerList)? {
            Value::Int4s(numbers) => Ok(numbers),
            _ => panic!("a Blast4-value integer-list is read as INTEGERs of int storage"),
        }
    }

    /// `GetCutoff()`: `Blast4-cutoff` (`e-value` REAL, `raw-score` INTEGER).
    fn cutoff(&self) -> Result<Cutoff, NcbiException> {
        match self.check(ValueType::Cutoff)? {
            Value::Choice(0, value) => match **value {
                Value::Real(e_value) => Ok(Cutoff::EValue(e_value)),
                _ => panic!("a Blast4-cutoff e-value is read as a REAL"),
            },
            Value::Choice(1, value) => match **value {
                Value::Int(_) => Ok(Cutoff::RawScore),
                _ => panic!("a Blast4-cutoff raw-score is read as an INTEGER"),
            },
            _ => panic!("a Blast4-cutoff is read as a CHOICE of two"),
        }
    }
}

/// `Blast4-cutoff`: an e-value, or a raw score (whose number no reader of the builder takes).
enum Cutoff {
    EValue(f64),
    RawScore,
}

/// The text of a string value as NCBI's `const char*` and `std::string` setters take it (the
/// reader keeps bytes 0x20..0x7E only, so no NUL ends a `c_str()` early).
fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// One `Set*` call of the builder on the options (`CBlastOptions`, or the handle's, which calls
/// it), with the value as NCBI's setter takes it.
#[derive(Clone, Debug, PartialEq)]
pub enum Setter {
    BestHitScoreEdge(f64),
    BestHitOverhang(f64),
    /// `SetCompositionBasedStats((ECompoAdjustModes) v)`: the number as read.
    CompositionBasedStats(i32),
    CullingLimit(i32),
    DbGeneticCode(i32),
    DbLength(i64),
    DustFiltering(bool),
    DustFilteringLevel(i32),
    DustFilteringWindow(i32),
    DustFilteringLinker(i32),
    DomainInclusionThreshold(f64),
    EffectiveSearchSpace(i64),
    EvalueThreshold(f64),
    /// `SetFilterString(v, true)`.
    FilterString(String),
    GapExtensionCost(i32),
    GapOpeningCost(i32),
    /// `SetGapTracebackAlgorithm((EBlastTbackExt) v)`: the number as read.
    GapTracebackAlgorithm(i32),
    GapTrigger(f64),
    GapXDropoff(f64),
    GapXDropoffFinal(f64),
    /// `SetGapExtnAlgorithm(static_cast<EBlastPrelimGapExt>(v))`: the number as read.
    GapExtnAlgorithm(i32),
    HitlistSize(i32),
    InclusionThreshold(f64),
    IgnoreMsaMaster(bool),
    LongestIntronLength(i32),
    /// `SetMBTemplateLength(unsigned char)`: the `int` converted (modulo 256).
    MbTemplateLength(u8),
    /// `SetMBTemplateType(unsigned char)`: the `int` converted (modulo 256).
    MbTemplateType(u8),
    MatchReward(i32),
    MatrixName(String),
    MismatchPenalty(i32),
    MaskAtHash(bool),
    MaxHspsPerSubject(i32),
    OutOfFrameMode(bool),
    /// `SetPHIPattern(pattern, is_na)`, `is_na` = the options' query is nucleotide.
    PhiPattern {
        pattern: String,
        is_na: bool,
    },
    PercentIdentity(f64),
    PseudoCount(i32),
    QueryGeneticCode(i32),
    QueryCovHspPerc(f64),
    RepeatFiltering(bool),
    RepeatFilteringDb(String),
    /// `SetStrandOption((ENa_strand) v)`: the `Blast4-strand-type` number (forward 1, reverse 2,
    /// both 3, the numbers of `ENa_strand`'s plus, minus and both).
    StrandOption(i32),
    SegFiltering(bool),
    SegFilteringWindow(i32),
    SegFilteringLocut(f64),
    SegFilteringHicut(f64),
    SumStatisticsMode(bool),
    SmithWatermanMode(bool),
    /// `SetGappedMode(!v)` of `UngappedMode v`.
    GappedMode(bool),
    UnifiedP(i32),
    /// `SetSubjectBestHit()`, whatever the value.
    SubjectBestHit,
    WindowMaskerTaxId(i32),
    WindowSize(i32),
    WordSize(i32),
    /// `SetWordThreshold(double)` of an `int`.
    WordThreshold(i32),
    /// `SetUseIndex(true, m_MbIndexName, m_ForceMbIndex)` (`x_ApplyInteractions`): use the
    /// megablast index `index_name`, forced; `old_style_index` keeps its default (false).
    UseIndex {
        index_name: String,
    },
}

// NCBI reference (598d8ae6): c++/include/util/range.hpp:290-308
/// `TSeqRange` (`CRange<TSeqPos>`): `from` and `to_open`, the empty range `GetEmpty()` both
/// `GetPositionMax()`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeqRange {
    pub from: u32,
    pub to_open: u32,
}

impl SeqRange {
    /// `TSeqRange::GetEmpty()`.
    pub const EMPTY: SeqRange = SeqRange {
        from: u32::MAX,
        to_open: u32::MAX,
    };
}

/// What the builder makes of a request's options.
#[derive(Clone, Debug, PartialEq)]
pub struct BuiltOptions {
    /// The program of the options (`ComputeProgram`, then `AdjustProgram`).
    pub program: EProgram,
    /// The calls on the factory options of `program`, in NCBI's order.
    pub setters: Vec<Setter>,
    /// `EntrezQuery` (`HaveEntrezQuery`, `GetEntrezQuery`).
    pub entrez_query: Option<String>,
    /// `FirstDbSeq`, `FinalDbSeq`.
    pub first_db_seq: Option<i32>,
    pub final_db_seq: Option<i32>,
    /// `GiList`, `NegativeGiList`: each `int` as `(Uint4)` (`GI_FROM(Uint4, gi_v)`).
    pub gi_list: Option<Vec<i64>>,
    pub negative_gi_list: Option<Vec<i64>>,
    /// `TaxidList`, `NegativeTaxidList`.
    pub taxid_list: Option<Vec<i32>>,
    pub negative_taxid_list: Option<Vec<i32>>,
    /// `DbFilteringAlgorithmId`, `DbFilteringAlgorithmKey`.
    pub db_filtering_algorithm_id: Option<i32>,
    pub db_filtering_algorithm_key: Option<String>,
    /// `SubjectMaskingType` (`(ESubjectMaskingType) v`: the number as read).
    pub subject_masking_type: Option<i32>,
    /// `LCaseMask` values (`Blast4-mask`), those of the algorithm options, else of the program
    /// options (`HaveQueryMasks`, `GetQueryMasks`).
    pub query_masks: Option<Vec<Value>>,
    /// `GetRestrictedQueryRange()`: `RequiredStart` sets `from`, `RequiredEnd` `to_open`.
    pub query_range: SeqRange,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-229
// The throws of the builder. NCBI's applications print only their message (`GetMsg()`, in the
// `CBlastException` clause of `CATCH_ALL`, whose `eInvalidOptions` none of these is), so the
// function of the report is not known and not given.
const COMPUTE_PROGRAM: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/blast_options_builder.cpp",
    line: 134,
    kind: "CBlastException::eInvalidArgument",
    function: None,
};
const OPTION_HAS_NO_NAME: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/blast_options_builder.cpp",
    line: 154,
    kind: "CBlastException::eInvalidArgument",
    function: None,
};
const EVALUE_TYPE: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/blast_options_builder.cpp",
    line: 232,
    kind: "CBlastException::eInvalidArgument",
    function: None,
};
const UNKNOWN_OPTION: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/blast_options_builder.cpp",
    line: 541,
    kind: "CRemoteBlastException::eServiceNotAvailable",
    function: None,
};
const PHI_PROGRAM: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/blast_options_builder.cpp",
    line: 616,
    kind: "CRemoteBlastException::eServiceNotAvailable",
    function: None,
};
const PROGRAM_NAME: ThrowSite = ThrowSite {
    file: "src/algo/blast/api/blast_aux.cpp",
    line: 808,
    kind: "CBlastException::eNotSupported",
    function: None,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:63-138
/// `ComputeProgram`: the program of a request's program and service.
pub fn compute_program(program: &[u8], service: &[u8]) -> Result<EProgram, NcbiException> {
    // `NStr::ToLower` lowers the ASCII letters of the bytes (`tolower` of the "C" locale that
    // NCBI's applications keep).
    let p = program.to_ascii_lowercase();
    let s = service.to_ascii_lowercase();
    let mut name: &[u8] = &p;
    let mut found = false;
    match &p[..] {
        b"blastp" => match &s[..] {
            b"rpsblast" => {
                name = b"rpsblast";
                found = true;
            }
            b"psi" => {
                name = b"psiblast";
                found = true;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:85-87
            // The services "phi" and "indexed" keep the program BLASTP.
            b"phi" | b"indexed" => found = true,
            b"delta_blast" => {
                name = b"deltablast";
                found = true;
            }
            _ => {}
        },
        b"blastn" => {
            if s == b"megablast" {
                name = b"megablast";
                found = true;
            }
            if s == b"vecscreen" {
                name = b"vecscreen";
                found = true;
            }
            if s == b"sra" || s == b"wgs" || s == b"indexed" {
                found = true;
            }
        }
        b"tblastn" => {
            if s == b"rpsblast" {
                name = b"rpstblastn";
                found = true;
            } else if s == b"psi" {
                name = b"psitblastn";
                found = true;
            }
            if s == b"sra" || s == b"wgs" {
                found = true;
            }
        }
        b"tblastx" => found = true,
        b"blastx" => {
            if s == b"rpsblast" {
                name = b"rpstblastn";
                found = true;
            }
        }
        _ => {}
    }
    if !(found || s == b"plain" || s == b"multi_blast") {
        return Err(NcbiException::new(
            COMPUTE_PROGRAM,
            format!(
                "Unsupported combination of program ({}) and service ({}).",
                text(program),
                text(service)
            ),
        ));
    }
    program_name_to_enum(name)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:757-810
/// `ProgramNameToEnum` of the name that `ComputeProgram` made (lower case already).
fn program_name_to_enum(name: &[u8]) -> Result<EProgram, NcbiException> {
    let lower = name.to_ascii_lowercase();
    Ok(
        if lower.starts_with(b"blastn") || lower.starts_with(b"rmblastn") {
            EProgram::Blastn
        } else if lower.starts_with(b"blastp") {
            EProgram::Blastp
        } else {
            match &lower[..] {
                b"blastx" => EProgram::Blastx,
                b"tblastn" => EProgram::Tblastn,
                b"tblastx" => EProgram::Tblastx,
                b"rpsblast" => EProgram::RpsBlast,
                b"rpstblastn" => EProgram::RpsTblastn,
                b"megablast" => EProgram::Megablast,
                b"psiblast" => EProgram::PsiBlast,
                b"psitblastn" => EProgram::PsiTblastn,
                b"dc-megablast" => EProgram::DiscMegablast,
                b"deltablast" => EProgram::DeltaBlast,
                b"vecscreen" => EProgram::VecScreen,
                b"mapper" | b"mapr2g" | b"mapr2r" | b"mapg2g" => EProgram::Mapper,
                _ => {
                    return Err(NcbiException::new(
                        PROGRAM_NAME,
                        format!("Program type '{}' not supported", text(name)),
                    ))
                }
            }
        },
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:585-645
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:649-665
/// `AdjustProgram` over the merged lists (algorithm, program, format options in this order,
/// `s_MergeCBlast4_parameters`): a `PHIPattern` makes a
/// BLASTN or BLASTP program PHI-BLAST's, and refuses any other; then the first parameter that
/// tells, in the lists' order: a non-zero `MBTemplateLength` makes dc-megablast, a
/// `Web_StepNumber`, `Web_RunPsiBlast`, `PseudoCountWeight` or `IgnoreMsaMaster` psiblast, of any
/// program. `program_string` is the request's program as written.
pub fn adjust_program(
    all: &[RawParam<'_>],
    program: EProgram,
    program_string: &[u8],
) -> Result<EProgram, NcbiException> {
    if all.iter().any(|p| p.is(fields::PHI_PATTERN)) {
        return match program {
            EProgram::PhiBlastn | EProgram::Blastn => Ok(EProgram::PhiBlastn),
            EProgram::PhiBlastp | EProgram::Blastp => Ok(EProgram::PhiBlastp),
            _ => Err(NcbiException::new(
                PHI_PROGRAM,
                format!(
                    "Incorrect combination of option ({}) and program ({})",
                    fields::PHI_PATTERN.name(),
                    text(program_string)
                ),
            )),
        };
    }
    for p in all {
        if p.is(fields::MB_TEMPLATE_LENGTH) {
            if Getter::of(p.value).integer()? != 0 {
                return Ok(EProgram::DiscMegablast);
            }
        } else if p.is(fields::WEB_STEP_NUMBER)
            || p.is(fields::WEB_RUN_PSI_BLAST)
            || p.is(fields::PSEUDO_COUNT)
            || p.is(fields::IGNORE_MSA_MASTER)
        {
            return Ok(EProgram::PsiBlast);
        }
    }
    Ok(program)
}

/// The builder's state while it reads the lists (the members of `CBlastOptionsBuilder` that the
/// options do not hold).
struct Builder {
    program: EProgram,
    setters: Vec<Setter>,
    perform_culling: bool,
    hsp_range_max: i32,
    force_mb_index: bool,
    mb_index_name: String,
    ignore_query_masks: bool,
    out: BuiltOptions,
}

impl Builder {
    fn new(program: EProgram) -> Builder {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:47-61
        // The constructor starts with culling off, the HSP range maximum 0, the query range
        // empty and the forced MegaBLAST index off; the optional members are unset and the
        // index name is empty.
        Builder {
            program,
            setters: Vec::new(),
            perform_culling: false,
            hsp_range_max: 0,
            force_mb_index: false,
            mb_index_name: String::new(),
            ignore_query_masks: false,
            out: BuiltOptions {
                program,
                setters: Vec::new(),
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
                query_range: SeqRange::EMPTY,
            },
        }
    }

    fn set(&mut self, setter: Setter) {
        self.setters.push(setter);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:140-169
    /// `x_ProcessOneOption`: the arm of the parameter's first letter, then the field that matches.
    /// Each arm cites its lines of blast_options_builder.cpp.
    fn process_one_option(&mut self, p: &RawParam<'_>) -> Result<(), NcbiException> {
        if p.name.is_empty() {
            return Err(NcbiException::new(
                OPTION_HAS_NO_NAME,
                "Option has no name.",
            ));
        }
        let v = Getter::of(p.value);
        let found = match p.name[0] {
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:170-178
            b'B' => {
                if p.is(fields::BEST_HIT_SCORE_EDGE) {
                    self.set(Setter::BestHitScoreEdge(v.real()?));
                } else if p.is(fields::BEST_HIT_OVERHANG) {
                    self.set(Setter::BestHitOverhang(v.real()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:180-193
            // `CutoffScore`'s field is a cutoff, read with `GetInteger()`, which throws
            // for every value the field matches.
            b'C' => {
                if p.is(fields::COMPOSITION_BASED_STATS) {
                    self.set(Setter::CompositionBasedStats(v.integer()?));
                } else if p.is(fields::CULLING) {
                    self.perform_culling = v.boolean()?;
                } else if p.is(fields::CULLING_LIMIT) {
                    self.set(Setter::CullingLimit(v.integer()?));
                } else if p.is(fields::CUTOFF_SCORE) {
                    v.integer()?;
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:195-217
            b'D' => {
                if p.is(fields::DB_GENETIC_CODE) {
                    self.set(Setter::DbGeneticCode(v.integer()?));
                } else if p.is(fields::DB_LENGTH) {
                    self.set(Setter::DbLength(v.big_integer()?));
                } else if p.is(fields::DUST_FILTERING) {
                    self.set(Setter::DustFiltering(v.boolean()?));
                } else if p.is(fields::DUST_FILTERING_LEVEL) {
                    self.set(Setter::DustFilteringLevel(v.integer()?));
                } else if p.is(fields::DUST_FILTERING_WINDOW) {
                    self.set(Setter::DustFilteringWindow(v.integer()?));
                } else if p.is(fields::DUST_FILTERING_LINKER) {
                    self.set(Setter::DustFilteringLinker(v.integer()?));
                } else if p.is(fields::DB_FILTERING_ALGORITHM_ID) {
                    self.out.db_filtering_algorithm_id = Some(v.integer()?);
                } else if p.is(fields::DB_FILTERING_ALGORITHM_KEY) {
                    self.out.db_filtering_algorithm_key = Some(text(v.string()?));
                } else if p.is(fields::DOMAIN_INCLUSION_THRESHOLD) {
                    self.set(Setter::DomainInclusionThreshold(v.real()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:219-237
            // `EvalueThreshold` is taken by its name too, as a real or a cutoff e-value.
            b'E' => {
                if p.is(fields::EFFECTIVE_SEARCH_SPACE) {
                    self.set(Setter::EffectiveSearchSpace(v.big_integer()?));
                } else if p.is(fields::ENTREZ_QUERY) {
                    self.out.entrez_query = Some(text(v.string()?));
                } else if p.is(fields::EVALUE_THRESHOLD)
                    || p.name == fields::EVALUE_THRESHOLD.name().as_bytes()
                {
                    let evalue = match v.0 {
                        ValueType::Real => v.real()?,
                        ValueType::Cutoff => match v.cutoff()? {
                            Cutoff::EValue(evalue) => evalue,
                            Cutoff::RawScore => return Err(evalue_type()),
                        },
                        _ => return Err(evalue_type()),
                    };
                    self.set(Setter::EvalueThreshold(evalue));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:239-251
            b'F' => {
                if p.is(fields::FILTER_STRING) {
                    self.set(Setter::FilterString(text(v.string()?)));
                } else if p.is(fields::FINAL_DB_SEQ) {
                    self.out.final_db_seq = Some(v.integer()?);
                } else if p.is(fields::FIRST_DB_SEQ) {
                    self.out.first_db_seq = Some(v.integer()?);
                } else if p.is(fields::FORCE_MB_INDEX) {
                    self.force_mb_index = v.boolean()?;
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:253-297
            // The GI list's numbers are each taken as an unsigned 32-bit number (in the
            // `NCBI_STRICT_GI` and the `NCBI_INT8_GI` branch alike: the field's value is an integer
            // list, never a big-integer list).
            b'G' => {
                if p.is(fields::GAP_EXTENSION_COST) {
                    self.set(Setter::GapExtensionCost(v.integer()?));
                } else if p.is(fields::GAP_OPENING_COST) {
                    self.set(Setter::GapOpeningCost(v.integer()?));
                } else if p.is(fields::GI_LIST) {
                    self.out.gi_list = Some(gis(v.integer_list()?));
                } else if p.is(fields::GAP_TRACEBACK_ALGORITHM) {
                    self.set(Setter::GapTracebackAlgorithm(v.integer()?));
                } else if p.is(fields::GAP_TRIGGER) {
                    self.set(Setter::GapTrigger(v.real()?));
                } else if p.is(fields::GAP_X_DROPOFF) {
                    self.set(Setter::GapXDropoff(v.real()?));
                } else if p.is(fields::GAP_X_DROPOFF_FINAL) {
                    self.set(Setter::GapXDropoffFinal(v.real()?));
                } else if p.is(fields::GAP_EXTN_ALGORITHM) {
                    self.set(Setter::GapExtnAlgorithm(v.integer()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:299-307
            b'H' => {
                if p.is(fields::HITLIST_SIZE) {
                    self.set(Setter::HitlistSize(v.integer()?));
                } else if p.is(fields::HSP_RANGE_MAX) {
                    self.hsp_range_max = v.integer()?;
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:309-317
            b'I' => {
                if p.is(fields::INCLUSION_THRESHOLD) {
                    self.set(Setter::InclusionThreshold(v.real()?));
                } else if p.is(fields::IGNORE_MSA_MASTER) {
                    self.set(Setter::IgnoreMsaMaster(v.boolean()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:319-340
            // The masks of the program options are skipped when the algorithm options
            // had one (`m_IgnoreQueryMasks`, set by `GetSearchOptions`).
            b'L' => {
                if p.is(fields::L_CASE_MASK) {
                    if !self.ignore_query_masks {
                        self.out
                            .query_masks
                            .get_or_insert_with(Vec::new)
                            .push(v.1.clone());
                    }
                } else if p.is(fields::LONGEST_INTRON_LENGTH) {
                    self.set(Setter::LongestIntronLength(v.integer()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:342-365
            // The `MatrixTable` parameter is accepted and has no effect.
            b'M' => {
                if p.is(fields::MB_TEMPLATE_LENGTH) {
                    self.set(Setter::MbTemplateLength(v.integer()? as u8));
                } else if p.is(fields::MB_TEMPLATE_TYPE) {
                    self.set(Setter::MbTemplateType(v.integer()? as u8));
                } else if p.is(fields::MATCH_REWARD) {
                    self.set(Setter::MatchReward(v.integer()?));
                } else if p.is(fields::MATRIX_NAME) {
                    self.set(Setter::MatrixName(text(v.string()?)));
                } else if p.is(fields::MATRIX_TABLE) {
                } else if p.is(fields::MISMATCH_PENALTY) {
                    self.set(Setter::MismatchPenalty(v.integer()?));
                } else if p.is(fields::MASK_AT_HASH) {
                    self.set(Setter::MaskAtHash(v.boolean()?));
                } else if p.is(fields::MB_INDEX_NAME) {
                    self.mb_index_name = text(v.string()?);
                } else if p.is(fields::MAX_HSPS_PER_SUBJECT) {
                    self.set(Setter::MaxHspsPerSubject(v.integer()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:367-373
            b'O' => {
                if p.is(fields::OUT_OF_FRAME_MODE) {
                    self.set(Setter::OutOfFrameMode(v.boolean()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:375-416
            // `TTaxId` is an `int`: each list number is converted to the tax-id type as it is.
            b'N' => {
                if p.is(fields::NEGATIVE_GI_LIST) {
                    self.out.negative_gi_list = Some(gis(v.integer_list()?));
                } else if p.is(fields::NEGATIVE_TAXID_LIST) {
                    self.out.negative_taxid_list = Some(v.integer_list()?.to_vec());
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:418-431
            // An empty pattern sets nothing.
            b'P' => {
                if p.is(fields::PHI_PATTERN) {
                    let pattern = v.string()?;
                    if !pattern.is_empty() {
                        self.set(Setter::PhiPattern {
                            pattern: text(pattern),
                            is_na: self.program.query_is_nucleotide(),
                        });
                    }
                } else if p.is(fields::PERCENT_IDENTITY) {
                    self.set(Setter::PercentIdentity(v.real()?));
                } else if p.is(fields::PSEUDO_COUNT) {
                    self.set(Setter::PseudoCount(v.integer()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:433-441
            b'Q' => {
                if p.is(fields::QUERY_GENETIC_CODE) {
                    self.set(Setter::QueryGeneticCode(v.integer()?));
                } else if p.is(fields::QUERY_COV_HSP_PERC) {
                    self.set(Setter::QueryCovHspPerc(v.real()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:443-455
            // `SetFrom` and `SetToOpen` take a `TSeqPos` (`unsigned int`).
            b'R' => {
                if p.is(fields::REPEAT_FILTERING) {
                    self.set(Setter::RepeatFiltering(v.boolean()?));
                } else if p.is(fields::REPEAT_FILTERING_DB) {
                    self.set(Setter::RepeatFilteringDb(text(v.string()?)));
                } else if p.is(fields::REQUIRED_START) {
                    self.out.query_range.from = v.integer()? as u32;
                } else if p.is(fields::REQUIRED_END) {
                    self.out.query_range.to_open = v.integer()? as u32;
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:457-479
            b'S' => {
                if p.is(fields::STRAND_OPTION) {
                    self.set(Setter::StrandOption(v.strand_type()?));
                } else if p.is(fields::SEG_FILTERING) {
                    self.set(Setter::SegFiltering(v.boolean()?));
                } else if p.is(fields::SEG_FILTERING_WINDOW) {
                    self.set(Setter::SegFilteringWindow(v.integer()?));
                } else if p.is(fields::SEG_FILTERING_LOCUT) {
                    self.set(Setter::SegFilteringLocut(v.real()?));
                } else if p.is(fields::SEG_FILTERING_HICUT) {
                    self.set(Setter::SegFilteringHicut(v.real()?));
                } else if p.is(fields::SUM_STATISTICS_MODE) {
                    self.set(Setter::SumStatisticsMode(v.boolean()?));
                } else if p.is(fields::SMITH_WATERMAN_MODE) {
                    self.set(Setter::SmithWatermanMode(v.boolean()?));
                } else if p.is(fields::SUBJECT_MASKING_TYPE) {
                    self.out.subject_masking_type = Some(v.integer()?);
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:481-496
            b'T' => {
                if p.is(fields::TAXID_LIST) {
                    self.out.taxid_list = Some(v.integer_list()?.to_vec());
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:497-511
            // `GappedMode` is stored inverted: the request's value is the negation of the
            // option's.
            b'U' => {
                if p.is(fields::GAPPED_MODE) {
                    self.set(Setter::GappedMode(!v.boolean()?));
                } else if p.is(fields::UNIFIED_P) {
                    self.set(Setter::UnifiedP(v.integer()?));
                } else if p.is(fields::SUBJECT_BEST_HIT) {
                    self.set(Setter::SubjectBestHit);
                } else {
                    return self.not_found(p);
                }
                true
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:513-525
            // No arm for `WindowMaskerDatabase`, which NCBI's own export writes
            // (`I.blastn.wmdb.alone`: "Internal: Error processing option [WindowMaskerDatabase]
            // type [12].").
            b'W' => {
                if p.is(fields::WINDOW_MASKER_TAX_ID) {
                    self.set(Setter::WindowMaskerTaxId(v.integer()?));
                } else if p.is(fields::WINDOW_SIZE) {
                    self.set(Setter::WindowSize(v.integer()?));
                } else if p.is(fields::WORD_SIZE) {
                    self.set(Setter::WordSize(v.integer()?));
                } else if p.is(fields::WORD_THRESHOLD) {
                    self.set(Setter::WordThreshold(v.integer()?));
                } else {
                    return self.not_found(p);
                }
                true
            }
            _ => false,
        };
        if found {
            Ok(())
        } else {
            self.not_found(p)
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:531-544
    // NCBI reference (598d8ae6): c++/include/algo/blast/api/search_strategy.hpp:103-105
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2018-2027
    /// A parameter that no arm takes. Both NCBI callers build with `ignore_unsupported_options`
    /// false.
    fn not_found(&self, p: &RawParam<'_>) -> Result<(), NcbiException> {
        Err(NcbiException::new(
            UNKNOWN_OPTION,
            format!(
                "Internal: Error processing option [{}] type [{}].",
                p.name_text(),
                p.value_type().which()
            ),
        ))
    }

    /// `x_ProcessOptions`: the list in its order.
    fn process_options(&mut self, list: &[RawParam<'_>]) -> Result<(), NcbiException> {
        for p in list {
            self.process_one_option(p)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:561-571
    fn apply_interactions(&mut self) {
        if self.perform_culling {
            self.set(Setter::CullingLimit(self.hsp_range_max));
        }
        if self.force_mb_index {
            let index_name = self.mb_index_name.clone();
            self.set(Setter::UseIndex { index_name });
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:230-232
/// The exception for an `EvalueThreshold` parameter that is neither a real nor an e-value cutoff.
fn evalue_type() -> NcbiException {
    NcbiException::new(EVALUE_TYPE, "EvalueThreshold has unsupported type.")
}

/// A GI list: each number as `(Uint4)`, then a `TGi` (`Int8`).
fn gis(numbers: &[i32]) -> Vec<i64> {
    numbers.iter().map(|&n| i64::from(n as u32)).collect()
}

/// The three option lists of a request (`algorithm-options`, `program-options`,
/// `format-options`), each `None` when the request has none, in file order with duplicates.
#[derive(Clone, Copy, Debug, Default)]
pub struct OptionLists<'a> {
    pub algorithm: Option<&'a [RawParam<'a>]>,
    pub program: Option<&'a [RawParam<'a>]>,
    pub format: Option<&'a [RawParam<'a>]>,
}

/// `CBlastOptionsBuilder(program, service, locality, false).GetSearchOptions(algorithm,
/// program options, format options)`: the options of a request.
pub fn build(
    program: &[u8],
    service: &[u8],
    lists: OptionLists<'_>,
) -> Result<BuiltOptions, NcbiException> {
    let computed = compute_program(program, service)?;
    let all: Vec<RawParam<'_>> = [lists.algorithm, lists.program, lists.format]
        .into_iter()
        .flatten()
        .flatten()
        .copied()
        .collect();
    let adjusted = adjust_program(&all, computed, program)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:85-88
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:299-313
    // `CBlastOptionsFactory::Create(program, eBoth)`: the options have the local and the remote
    // part, so the Blast4 parameter list that `-outfmt 11` and `-export_search_strategy` write is
    // kept.
    let mut builder = Builder::new(adjusted);
    builder.ignore_query_masks = false;
    builder.process_options(lists.algorithm.unwrap_or_default())?;
    builder.ignore_query_masks = builder.out.query_masks.is_some();
    builder.process_options(lists.program.unwrap_or_default())?;
    builder.apply_interactions();
    let mut out = builder.out;
    out.setters = builder.setters;
    Ok(out)
}

// NCBI reference (598d8ae6): c++/include/algo/blast/api/search_strategy.hpp:60-95
/// What `CImportStrategy::FetchData` keeps (`CImportStrategyData`)
/// beside the builder's options (its `m_QueryRange` is the options' `query_range`).
#[derive(Clone, Debug, PartialEq)]
pub struct ImportedOptions {
    pub options: BuiltOptions,
    /// `m_PsiNumOfIterations` (`unsigned int`, 0 when the request has no `Web_StepNumber`).
    pub psi_num_of_iterations: u32,
    /// `m_FilteringID`: -1, or the builder's `DbFilteringAlgorithmId`.
    pub filtering_id: i32,
    /// `m_FilteringKey`: empty, or the builder's `DbFilteringAlgorithmKey`.
    pub filtering_key: String,
    /// `m_SubjectMaskingType`: `eNoSubjMasking` (0), or the builder's non-zero type.
    pub subject_masking_type: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:98-128
/// `CImportStrategy::FetchData`'s options (`-import_search_strategy`): the iteration count of the
/// first parameter named `Web_StepNumber` in the format options (`GetParamByName`; an integer, or
/// the number of a string, 0 when it is none: `NStr::StringToInt` with `fConvErr_NoThrow`; a value
/// of another type fails in `GetString()`), then the options built with `eBoth`.
pub fn import_options(
    program: &[u8],
    service: &[u8],
    lists: OptionLists<'_>,
) -> Result<ImportedOptions, NcbiException> {
    let mut psi_num_of_iterations = 0;
    if let Some(format) = lists.format {
        let step = format
            .iter()
            .find(|p| p.name == fields::WEB_STEP_NUMBER.name().as_bytes());
        if let Some(p) = step {
            let v = Getter::of(p.value);
            psi_num_of_iterations = match v.integer() {
                Ok(n) => n as u32,
                Err(_) => crate::ncbi_str::string_to_int(v.string()?, 0).unwrap_or(0) as u32,
            };
        }
    }
    let options = build(program, service, lists)?;
    Ok(ImportedOptions {
        filtering_id: options.db_filtering_algorithm_id.unwrap_or(-1),
        filtering_key: options
            .db_filtering_algorithm_key
            .clone()
            .unwrap_or_default(),
        subject_masking_type: options.subject_masking_type.unwrap_or(0),
        psi_num_of_iterations,
        options,
    })
}

/// [`import_options`] of a strategy's request as the archive model read it ([`Request::params`]:
/// each list in the file's order, a name as often as it is stored). An absent list is given as
/// an empty one: NCBI's builder and `FetchData` take nothing from either.
pub(crate) fn import_request(request: &Request) -> Result<ImportedOptions, NcbiException> {
    let list = |set: ParamSet| -> Vec<RawParam<'_>> {
        request
            .params()
            .iter()
            .filter(|param| param.set == set)
            .map(|param| RawParam {
                name: &param.name,
                value: param.value().value,
            })
            .collect()
    };
    let (algorithm, program, format) = (
        list(ParamSet::Algorithm),
        list(ParamSet::Program),
        list(ParamSet::Format),
    );
    import_options(
        &request.program,
        &request.service,
        OptionLists {
            algorithm: Some(&algorithm),
            program: Some(&program),
            format: Some(&format),
        },
    )
}

#[cfg(test)]
#[path = "options_tests.rs"]
mod tests;
