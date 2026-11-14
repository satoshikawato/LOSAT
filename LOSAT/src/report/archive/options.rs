//! The request options of a BLAST archive (`algorithm-options` of the `queue-search`
//! request): the parameters that the `Set*` calls of NCBI's application record.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:257-269
//! ```c++
//! void CBlastOptions::x_AttachValue(CRef<CBlast4_parameter> p) ...
//! ```
//! With `-outfmt 11` the options handle is local and remote (blast_args.cpp:3624-3627), so
//! every `Set*` call made outside defaults mode records a Blast4 parameter: appended when it
//! is first set and replaced in place when it is set again. The list is therefore a fold of
//! the call sequence: (1) the handle (`EvalueThreshold` 10); (2) the task
//! (`CBlastOptionsFactory`); (3) the argument classes in the order of `<program>_args.cpp`
//! (blast_args.cpp:3631-3633), each with the conditions below ("given" is on the command
//! line); (4) what the search sets in the handle after the options: `EffectiveSearchSpace` of
//! each searched batch, and RPS's database length and `.aux` costs (`super::HandleCarry`, applied
//! onto this list by the writer with the same rule).
//! `S/DESIGN_OF8.md` section 2.2 lists the classes; the unit tests of the programs' adapters
//! replay every option list that NCBI wrote in the design's probes.

use crate::blastinput::app::Program;

/// The value of a recorded parameter (`Blast4-value`).
#[derive(Clone, Debug, PartialEq)]
pub enum ParamValue {
    /// `cutoff e-value`.
    Evalue(f64),
    /// `integer` (an `int` setter keeps its type, e.g. `WordThreshold`).
    Integer(i32),
    /// `big-integer`.
    BigInteger(i64),
    /// `boolean`.
    Boolean(bool),
    /// `real`.
    Real(f64),
    /// `string`.
    Str(String),
    /// `integer-list`.
    IntegerList(Vec<i32>),
}

/// A recorded parameter (`Blast4-parameter`).
#[derive(Clone, Debug, PartialEq)]
pub struct Param {
    pub name: &'static str,
    pub value: ParamValue,
}

/// A filtering option as the command line gives it (`-seg`, `-dust`).
#[derive(Clone, Debug, PartialEq)]
pub enum Filter {
    /// `no`.
    No,
    /// `yes`.
    Yes,
    /// `-seg WINDOW LOCUT HICUT`.
    Seg { window: i32, locut: f64, hicut: f64 },
    /// `-dust LEVEL WINDOW LINKER`.
    Dust {
        level: i32,
        window: i32,
        linker: i32,
    },
}

impl From<crate::blastinput::value_parsers::SegSpec> for Filter {
    fn from(spec: crate::blastinput::value_parsers::SegSpec) -> Self {
        use crate::blastinput::value_parsers::SegSpec;
        match spec {
            SegSpec::No => Filter::No,
            SegSpec::Yes => Filter::Yes,
            SegSpec::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => Filter::Seg {
                window,
                locut,
                hicut,
            },
        }
    }
}

/// The filter of a `-seg` value as given (`app::parse_seg_option`; `label` names the program in
/// its error), `None` for a value the program refuses (no request is written then).
pub fn seg_filter(value: &str, label: &str) -> Option<Filter> {
    crate::blastinput::app::parse_seg_option(value, label)
        .ok()
        .map(Filter::from)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:73-88
// ```c++
// void
// CBlastOptionsHandle::SetDefaults()
// {
//     if (m_Opts->GetLocality() != CBlastOptions::eRemote) {
//         m_Opts->SetDefaultsMode(true);
//         SetLookupTableDefaults();
//         ...
//         m_Opts->SetDefaultsMode(false);
//     }
//     SetRemoteProgramAndService_Blast3();
// }
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/psiblast_options.hpp:82-84
// ```c++
//     virtual void SetRemoteProgramAndService_Blast3()
//     {
//         m_Opts->SetRemoteProgramAndService_Blast3("blastp", "psi");
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_rps_options.hpp:122-124
// ```c++
//     virtual void SetRemoteProgramAndService_Blast3()
//     {
//         m_Opts->SetRemoteProgramAndService_Blast3("blastp", "rpsblast");
// ```
// NCBI reference (598d8ae6): c++/include/algo/blast/api/rpstblastn_options.hpp:74-76
// ```c++
//     virtual void SetRemoteProgramAndService_Blast3()
//     {
//         m_Opts->SetRemoteProgramAndService_Blast3("tblastn", "rpsblast");
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:138-145
// ```c++
//     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
//        CPSIBlastOptionsHandle * rv(new CPSIBlastOptionsHandle(locality));
//        rv->SetPSITblastnDefaults();
//        if (args[kTask].AsString() == "blastx-fast") {
//            rv->SetWordSize(6);
//            rv->SetOptions().SetLookupTableType(eCompressedAaLookupTable);
//            rv->SetWordThreshold(21.0);
// ```
/// The kind of a run's request: its options handle (the `program` and `service` it names, the
/// setters of its constructor) and the argument classes of its application, in their order. The
/// five search programs, and the applications whose handle is another one's: psiblast with a
/// FASTA or a PSSM query (`CPSIBlastOptionsHandle`), tblastn with a PSSM query (PSI-TBLASTN:
/// `CPSIBlastOptionsHandle` with `SetPSITblastnDefaults`, whatever `-task`: the test names
/// `blastx-fast`, which is no tblastn task), rpsblast (`CBlastRPSOptionsHandle`) and rpstblastn
/// (`CRPSTBlastnOptionsHandle`). `Program` stays the reports' name of the search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RequestKind {
    Blastn,
    Blastp,
    Blastx,
    Tblastn,
    Tblastx,
    Psiblast,
    PsiTblastn,
    Rpsblast,
    Rpstblastn,
}

impl From<Program> for RequestKind {
    fn from(program: Program) -> Self {
        match program {
            Program::Blastn => RequestKind::Blastn,
            Program::Blastp => RequestKind::Blastp,
            Program::Blastx => RequestKind::Blastx,
            Program::Tblastn => RequestKind::Tblastn,
            Program::Tblastx => RequestKind::Tblastx,
        }
    }
}

impl RequestKind {
    /// The request's `program` (`GetRemoteProgramAndService_Blast3`).
    pub fn program(self) -> &'static str {
        match self {
            RequestKind::Blastn => "blastn",
            RequestKind::Blastp
            | RequestKind::Psiblast
            | RequestKind::PsiTblastn
            | RequestKind::Rpsblast => "blastp",
            RequestKind::Blastx => "blastx",
            RequestKind::Tblastn | RequestKind::Rpstblastn => "tblastn",
            RequestKind::Tblastx => "tblastx",
        }
    }

    /// Whether the handle is an RPS search's.
    fn is_rps(self) -> bool {
        matches!(self, RequestKind::Rpsblast | RequestKind::Rpstblastn)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:46-67
    // ```c++
    // CPSIBlastOptionsHandle::CPSIBlastOptionsHandle(EAPILocality locality)
    //     : CBlastAdvancedProteinOptionsHandle(locality)
    // {
    //     SetDefaults();
    //     m_Opts->SetProgram(ePSIBlast);
    //     if (m_Opts->GetLocality() == CBlastOptions::eRemote) {
    //         return;
    //     }
    //     SetPSIBlastDefaults();
    // }
    // ...
    // void CPSIBlastOptionsHandle::SetPSIBlastDefaults(void)
    // {
    //     m_Opts->SetInclusionThreshold(PSI_INCLUSION_ETHRESH);
    //     m_Opts->SetPseudoCount(PSI_PSEUDO_COUNT_CONST);
    // }
    //
    // void CPSIBlastOptionsHandle::SetPSITblastnDefaults(void)
    // {
    //     m_Opts->SetProgram(ePSITblastn);
    //     m_Opts->SetSmithWatermanMode(false);
    //     m_Opts->SetDbGeneticCode(BLAST_GENETIC_CODE);
    //     m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_TBLASTN);
    // }
    // ```
    /// The setters that the handle's constructor calls after defaults mode, after the
    /// `EvalueThreshold` that every handle records: psiblast's PSSM engine defaults
    /// (`PSI_INCLUSION_ETHRESH` 0.002, `PSI_PSEUDO_COUNT_CONST` 0), and for PSI-TBLASTN then
    /// `SetPSITblastnDefaults` (`SmithWatermanMode` FALSE, `DbGeneticCode` 1, `WordThreshold`
    /// `BLAST_WORD_THRESHOLD_TBLASTN` 13; `SetProgram` records nothing; NCBI 2.17.0
    /// `tbn.pssm.db.ss` of line RG F1d probe 7).
    fn constructor_calls(self) -> Vec<Param> {
        let psi = || {
            vec![
                Param {
                    name: "InclusionThreshold",
                    value: ParamValue::Real(0.002),
                },
                Param {
                    name: "PseudoCountWeight",
                    value: ParamValue::Integer(0),
                },
            ]
        };
        match self {
            RequestKind::Psiblast => psi(),
            RequestKind::PsiTblastn => {
                let mut calls = psi();
                calls.extend([
                    Param {
                        name: "SmithWatermanMode",
                        value: ParamValue::Boolean(false),
                    },
                    Param {
                        name: "DbGeneticCode",
                        value: ParamValue::Integer(1),
                    },
                    Param {
                        name: "WordThreshold",
                        value: ParamValue::Integer(13),
                    },
                ]);
                calls
            }
            _ => Vec::new(),
        }
    }
}

/// Every option whose value or presence an archive records, as `Option`s that each program
/// fills from its resolved arguments (`None` = not given).
#[derive(Clone, Debug)]
pub struct OptionInputs {
    pub kind: RequestKind,
    /// The task, its default included (`megablast`, `blastp`, ...).
    pub task: String,
    pub dbsize: Option<i64>,
    /// `-searchsp`: recorded where `CGenericSearchArgs::ExtractAlgorithmOptions` sets it
    /// (blast_args.cpp:303-307), else after the options (`writer`).
    pub searchsp: Option<i64>,
    pub evalue: Option<f64>,
    pub gapopen: Option<i32>,
    pub gapextend: Option<i32>,
    /// The `-matrix` name as typed when given (NCBI records it so, `blosum62` included; the
    /// suggestions below compare it without case, as NCBI's `strcasecmp` does).
    pub matrix: Option<String>,
    pub xdrop_gap: Option<f64>,
    pub xdrop_gap_final: Option<f64>,
    pub word_size: Option<i32>,
    pub perc_identity: Option<f64>,
    pub qcov_hsp_perc: Option<f64>,
    pub max_hsps: Option<i32>,
    pub sum_stats: Option<bool>,
    pub penalty: Option<i32>,
    pub reward: Option<i32>,
    /// `-template_type`: coding 0, optimal 1, coding_and_optimal 2.
    pub template_type: Option<i32>,
    pub template_length: Option<i32>,
    pub soft_masking: Option<bool>,
    /// `-seg` (protein queries): the argument's default when not given.
    pub seg: Option<Filter>,
    /// `-dust` (BLASTN), only when given.
    pub dust: Option<Filter>,
    /// `-window_masker_db` (BLASTN) as given.
    pub window_masker_db: Option<String>,
    /// `-window_masker_taxid` (BLASTN) as given.
    pub window_masker_taxid: Option<i32>,
    pub threshold: Option<i32>,
    /// `-culling_limit` (also 0), when given.
    pub culling_limit: Option<i32>,
    pub subject_besthit: bool,
    pub window_size: Option<i32>,
    pub ungapped: bool,
    /// `-no_greedy` (BLASTN): the dynamic programming algorithms are recorded.
    pub no_greedy: bool,
    /// The database whose MegaBLAST index a BLASTN run asked for and did not use (`-use_index
    /// true`, `-db`, task megablast; `blastn::blast_engine::run` `megablast_index`).
    pub mb_index_name: Option<String>,
    pub max_intron_length: i32,
    pub query_gencode: i32,
    pub db_gencode: i32,
    /// `-max_target_seqs`, else `-num_alignments`, else 500.
    pub hitlist_size: i32,
    pub comp_based_stats: i32,
    /// The `u` of `-comp_based_stats 2u` as given: `SetUnifiedP(1)` for BLASTP in an adjusting mode
    /// (`algorithm_options`).
    pub unified_p: bool,
    pub use_sw_tback: bool,
    /// psiblast's `-inclusion_ethresh` (`CPssmEngineArgs`: recorded in place of the constructor's).
    pub inclusion_ethresh: f64,
    /// psiblast's `-pseudocount` (`CPssmEngineArgs`).
    pub pseudocount: i32,
    /// psiblast's `-gap_trigger` (`CGapTriggerArgs`, always recorded).
    pub gap_trigger: f64,
    /// psiblast's `-num_iterations`: the request's `Web_StepNumber` with `-db`
    /// (`ArchiveRequest::step_number`).
    pub num_iterations: Option<i32>,
}

impl OptionInputs {
    /// The inputs of a run of `kind` (a program's or an application's, `RequestKind`) that gives
    /// no option: every field at its default.
    pub fn new(kind: impl Into<RequestKind>, task: &str) -> Self {
        Self {
            kind: kind.into(),
            task: task.to_string(),
            dbsize: None,
            searchsp: None,
            evalue: None,
            gapopen: None,
            gapextend: None,
            matrix: None,
            xdrop_gap: None,
            xdrop_gap_final: None,
            word_size: None,
            perc_identity: None,
            qcov_hsp_perc: None,
            max_hsps: None,
            sum_stats: None,
            penalty: None,
            reward: None,
            template_type: None,
            template_length: None,
            soft_masking: None,
            seg: None,
            dust: None,
            window_masker_db: None,
            window_masker_taxid: None,
            threshold: None,
            culling_limit: None,
            subject_besthit: false,
            window_size: None,
            ungapped: false,
            no_greedy: false,
            mb_index_name: None,
            max_intron_length: 0,
            query_gencode: 1,
            db_gencode: 1,
            hitlist_size: 500,
            comp_based_stats: 2,
            unified_p: false,
            use_sw_tback: false,
            inclusion_ethresh: 0.002,
            pseudocount: 0,
            gap_trigger: 22.0,
            num_iterations: None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:257-260
// ```c++
//     int gap_open=0, gap_extend=0;
//     if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
//          BLAST_GetProteinGapExistenceExtendParams
//              (args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);
// ```
/// The `(gap open, gap extend)` that `-matrix` implies (`BLAST_GetProteinGapExistenceExtendParams`:
/// the matrix's default row, which BLOSUM45 puts at 14 / 2). A matrix that NCBI does not know
/// leaves both at 0 (blast_stat.c:3383-3384), which NCBI BLAST+ 2.17.0 records in the search
/// strategy of a tblastx run whose search then fails on that matrix (RG-F1 probe 5 `mat`).
fn matrix_gaps(matrix: &str) -> (i32, i32) {
    crate::stats::protein_options::protein_gap_existence_extend_params(matrix).unwrap_or((0, 0))
}

/// The word threshold the matrix suggests before the program's addition
/// (`BLAST_GetSuggestedThreshold`, blast_options.c:1174-1209; `protein_options`).
fn suggested_threshold(matrix: &str) -> i32 {
    crate::stats::protein_options::suggested_threshold(
        crate::stats::protein_options::SuggestionProgram::Protein,
        matrix,
    ) as i32
}

/// The window size the matrix suggests (`BLAST_GetSuggestedWindowSize`,
/// blast_options.c:1211-1236; `protein_options`).
fn suggested_window(matrix: &str) -> i32 {
    crate::stats::protein_options::suggested_window_size(matrix)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:586-602
// ```c++
// static bool
// s_IsDefaultWordThreshold(EProgram program, double threshold)
// {
//     int word_threshold = static_cast<int>(threshold);
//     bool retval = true;
//     if (program == eBlastp &&
//         word_threshold != BLAST_WORD_THRESHOLD_BLASTP) {
//         retval = false;
//     } else if (program == eBlastx &&
//         word_threshold != BLAST_WORD_THRESHOLD_BLASTX) {
//         retval = false;
//     } else if (program == eTblastn &&
//                word_threshold != BLAST_WORD_THRESHOLD_TBLASTN) {
//         retval = false;
//     }
//     return retval;
// }
// ```
/// Whether the word threshold set before the threshold class (`value`) is the default that the
/// matrix's suggestion replaces: BLASTP, BLASTX and TBLASTN keep any other value; every other
/// program (psiblast's `ePSIBlast`, TBLASTX, ...) takes the suggestion whatever was set.
fn is_default_word_threshold(kind: RequestKind, value: i32) -> bool {
    match kind {
        RequestKind::Blastp => value == 11,
        RequestKind::Blastx => value == 12,
        RequestKind::Tblastn => value == 13,
        _ => true,
    }
}

/// The setters of the task (`CBlastOptionsFactory::Create(task)` outside defaults mode).
fn task_calls(task: &str) -> Vec<Param> {
    use ParamValue::{Boolean, Evalue, Integer, Str};
    let p = |name, value| Param { name, value };
    match task {
        "blastn-short" => vec![
            p("EvalueThreshold", Evalue(1000.0)),
            p("MatchReward", Integer(1)),
            p("MismatchPenalty", Integer(-3)),
            p("WordSize", Integer(7)),
            p("DustFiltering", Boolean(false)),
            p("SegFiltering", Boolean(false)),
            p("RepeatFiltering", Boolean(false)),
            // `ClearFilterOptions` ends with `SetMaskAtHash(false)` (api/blast_options_cxx.cpp:1022-1031).
            p("MaskAtHash", Boolean(false)),
        ],
        "dc-megablast" => vec![
            p("MBTemplateType", Integer(0)),
            p("MBTemplateLength", Integer(18)),
            p("WordSize", Integer(11)),
            p("WindowSize", Integer(40)),
        ],
        "blastp-fast" | "tblastn-fast" | "blastx-fast" => {
            vec![p("WordSize", Integer(5)), p("WordThreshold", Integer(20))]
        }
        "blastp-short" => vec![
            p("EvalueThreshold", Evalue(20000.0)),
            p("MatrixName", Str("PAM30".to_string())),
            p("GapOpeningCost", Integer(9)),
            p("GapExtensionCost", Integer(1)),
            p("WordSize", Integer(2)),
            p("WordThreshold", Integer(16)),
            p("DustFiltering", Boolean(false)),
            p("SegFiltering", Boolean(false)),
            p("RepeatFiltering", Boolean(false)),
            // `ClearFilterOptions` ends with `SetMaskAtHash(false)` (api/blast_options_cxx.cpp:1022-1031).
            p("MaskAtHash", Boolean(false)),
        ],
        _ => Vec::new(),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Class {
    Db,
    Generic,
    Nucl,
    Dmb,
    Filtering,
    Gapped,
    Hsp,
    Window,
    MbIndex,
    Formatting,
    Matrix,
    Threshold,
    Intron,
    QueryGeneticCode,
    DbGeneticCode,
    Cbs,
    /// psiblast's `CGapTriggerArgs`.
    GapTrigger,
    /// psiblast's `CPsiBlastArgs` (`SetIgnoreMsaMaster`).
    PsiBlast,
    /// psiblast's `CPssmEngineArgs`.
    PssmEngine,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:102-113 (the end of
// psiblast's `m_Args` pushes; rpsblast_args.cpp:70-127 and rpstblastn_args.cpp:45-118 for the RPS
// programs)
// ```c++
//     arg.Reset(new CCompositionBasedStatsArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CGapTriggerArgs(kQueryIsProtein));
//     m_Args.push_back(arg);
//
//     m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eProteinDb));
//     arg.Reset(m_PsiBlastArgs);
//     m_Args.push_back(arg);
//
//     arg.Reset(new CPssmEngineArgs);
//     m_Args.push_back(arg);
// ```
/// The argument classes of a program in the order of `<program>_args.cpp` (the order of
/// the `m_Args` pushes).
fn class_order(kind: RequestKind) -> &'static [Class] {
    use Class::*;
    match kind {
        RequestKind::Blastn => &[
            Db, Generic, Nucl, Dmb, Filtering, Gapped, Hsp, Window, MbIndex, Formatting,
        ],
        RequestKind::Blastp => &[
            Db, Generic, Filtering, Matrix, Threshold, Hsp, Window, Formatting, Gapped, Cbs,
        ],
        // psiblast has no gapped class; the PSI-BLAST classes come after the composition-based
        // statistics (`p9.ps.q.db.itdef.f6`, NCBI 2.17.0).
        RequestKind::Psiblast => &[
            Db, Generic, Filtering, Matrix, Threshold, Hsp, Window, Formatting, Cbs, GapTrigger,
            PsiBlast, PssmEngine,
        ],
        RequestKind::Rpsblast => &[Db, Generic, Filtering, Window, Formatting, Cbs],
        RequestKind::Rpstblastn => &[
            Db,
            Generic,
            QueryGeneticCode,
            Filtering,
            Window,
            Gapped,
            Formatting,
            Cbs,
        ],
        RequestKind::Blastx => &[
            Db,
            Generic,
            Intron,
            Filtering,
            Matrix,
            Threshold,
            Hsp,
            Window,
            Gapped,
            QueryGeneticCode,
            Formatting,
            Cbs,
        ],
        // PSI-TBLASTN's classes are tblastn's (`CTblastnAppArgs`, tblastn_args.cpp:44-131;
        // its `CPsiBlastArgs` of a nucleotide database records nothing).
        RequestKind::Tblastn | RequestKind::PsiTblastn => &[
            Db,
            Generic,
            DbGeneticCode,
            Gapped,
            Intron,
            Filtering,
            Matrix,
            Threshold,
            Hsp,
            Window,
            Formatting,
            Cbs,
        ],
        RequestKind::Tblastx => &[
            Db,
            Generic,
            Intron,
            Filtering,
            Matrix,
            Threshold,
            Hsp,
            Window,
            QueryGeneticCode,
            DbGeneticCode,
            Formatting,
        ],
    }
}

/// The folded parameter list: appended when first set, replaced in place when set again
/// (`x_AttachValue`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParamList {
    params: Vec<Param>,
}

impl ParamList {
    /// A `Set*` call.
    pub fn set(&mut self, name: &'static str, value: ParamValue) {
        match self.params.iter_mut().find(|param| param.name == name) {
            Some(param) => param.value = value,
            None => self.params.push(Param { name, value }),
        }
    }

    fn get(&self, name: &str) -> Option<&ParamValue> {
        self.params
            .iter()
            .find(|param| param.name == name)
            .map(|param| &param.value)
    }

    /// The parameters, in their order.
    pub fn params(&self) -> &[Param] {
        &self.params
    }

    pub fn into_params(self) -> Vec<Param> {
        self.params
    }
}

/// The recorded `algorithm-options` of a run, without `EffectiveSearchSpace` (which each
/// archive appends).
pub fn algorithm_options(inputs: &OptionInputs) -> Vec<Param> {
    use ParamValue::{BigInteger, Boolean, Evalue, Integer, Real, Str};
    let kind = inputs.kind;
    let mut list = ParamList::default();
    list.set("EvalueThreshold", Evalue(10.0));
    for param in kind.constructor_calls() {
        list.set(param.name, param.value);
    }
    // PSI-TBLASTN's handle is not the task's (`RequestKind`: tblastn_args.cpp:138-149).
    if kind != RequestKind::PsiTblastn {
        for param in task_calls(&inputs.task) {
            list.set(param.name, param.value);
        }
    }
    let matrix = inputs.matrix.clone().unwrap_or_else(|| {
        if inputs.task == "blastp-short" {
            "PAM30"
        } else {
            "BLOSUM62"
        }
        .to_string()
    });
    // NCBI's `CGenericSearchArgs` of BLASTX is made with the flag that turns the word threshold
    // rule on, as BLASTP's and TBLASTN's are (probes `opt_ws5` .. `opt_ws7`).
    let protein_query = matches!(
        kind,
        RequestKind::Blastp
            | RequestKind::Tblastn
            | RequestKind::Blastx
            | RequestKind::Psiblast
            | RequestKind::PsiTblastn
    );
    let mut threshold_set_by_word_size = false;
    for &class in class_order(kind) {
        match class {
            Class::Db => {
                if let Some(size) = inputs.dbsize {
                    list.set("DbLength", BigInteger(size));
                }
            }
            Class::Generic => {
                if let Some(evalue) = inputs.evalue {
                    list.set("EvalueThreshold", Evalue(evalue));
                }
                if let Some(open) = inputs.gapopen {
                    list.set("GapOpeningCost", Integer(open));
                } else if inputs.matrix.is_some() {
                    list.set("GapOpeningCost", Integer(matrix_gaps(&matrix).0));
                }
                if let Some(extend) = inputs.gapextend {
                    list.set("GapExtensionCost", Integer(extend));
                } else if inputs.matrix.is_some() {
                    list.set("GapExtensionCost", Integer(matrix_gaps(&matrix).1));
                }
                if let Some(xdrop) = inputs.xdrop_gap {
                    list.set("GapXDropoff", Real(xdrop));
                }
                if let Some(xdrop) = inputs.xdrop_gap_final {
                    list.set("GapXDropoffFinal", Real(xdrop));
                }
                if let Some(word_size) = inputs.word_size {
                    if protein_query && word_size > 4 {
                        threshold_set_by_word_size = true;
                        list.set("WordThreshold", Integer(19));
                        if word_size > 5 {
                            list.set("WordThreshold", Integer(21));
                        }
                        if word_size > 6 {
                            list.set("WordThreshold", Integer(20));
                        }
                    }
                    list.set("WordSize", Integer(word_size));
                }
                if let Some(searchsp) = inputs.searchsp {
                    list.set("EffectiveSearchSpace", BigInteger(searchsp));
                }
                if let Some(identity) = inputs.perc_identity {
                    list.set("PercentIdentity", Real(identity));
                }
                if let Some(coverage) = inputs.qcov_hsp_perc {
                    list.set("QueryCovHspPerc", Real(coverage));
                }
                if let Some(max_hsps) = inputs.max_hsps {
                    list.set("MaxHspsPerSubject", Integer(max_hsps));
                }
                if let Some(sum_stats) = inputs.sum_stats {
                    list.set("SumStatistics", Boolean(sum_stats));
                }
            }
            Class::Nucl => {
                if let Some(penalty) = inputs.penalty {
                    list.set("MismatchPenalty", Integer(penalty));
                }
                if let Some(reward) = inputs.reward {
                    list.set("MatchReward", Integer(reward));
                }
                // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:681-685
                // ```c++
                //     if (cmd_line_args.Exist(kArgNoGreedyExtension) &&
                //         cmd_line_args[kArgNoGreedyExtension]) {
                //         options.SetGapExtnAlgorithm(eDynProgScoreOnly);
                //         options.SetGapTracebackAlgorithm(eDynProgTbck);
                //     }
                // ```
                // Both enumerators are 0 (blast_options.h: eDynProgScoreOnly, eDynProgTbck).
                if inputs.no_greedy {
                    list.set("GapExtnAlgorithm", Integer(0));
                    list.set("GapTracebackAlgorithm", Integer(0));
                }
            }
            Class::Dmb => {
                if let Some(kind) = inputs.template_type {
                    list.set("MBTemplateType", Integer(kind));
                }
                if let Some(length) = inputs.template_length {
                    list.set("MBTemplateLength", Integer(length));
                }
            }
            Class::Filtering => {
                list.set(
                    "MaskAtHash",
                    Boolean(inputs.soft_masking.unwrap_or(kind == RequestKind::Blastn)),
                );
                if kind != RequestKind::Blastn {
                    // SEG is on by default for the translated searches only (`kFilterByDefault`).
                    let seg = inputs.seg.clone().unwrap_or(
                        if matches!(
                            kind,
                            RequestKind::Blastp
                                | RequestKind::Psiblast
                                | RequestKind::Rpsblast
                                | RequestKind::Rpstblastn
                        ) {
                            Filter::No
                        } else {
                            Filter::Seg {
                                window: 12,
                                locut: 2.2,
                                hicut: 2.5,
                            }
                        },
                    );
                    match seg {
                        Filter::No => list.set("SegFiltering", Boolean(false)),
                        Filter::Yes => list.set("SegFiltering", Boolean(true)),
                        Filter::Seg {
                            window,
                            locut,
                            hicut,
                        } => {
                            list.set("SegFilteringWindow", Integer(window));
                            list.set("SegFilteringLocut", Real(locut));
                            list.set("SegFilteringHicut", Real(hicut));
                        }
                        Filter::Dust { .. } => {}
                    }
                } else if let Some(dust) = inputs.dust.clone() {
                    match dust {
                        Filter::No => list.set("DustFiltering", Boolean(false)),
                        Filter::Yes => list.set("DustFiltering", Boolean(true)),
                        Filter::Dust {
                            level,
                            window,
                            linker,
                        } => {
                            list.set("DustFilteringLevel", Integer(level));
                            list.set("DustFilteringWindow", Integer(window));
                            list.set("DustFilteringLinker", Integer(linker));
                        }
                        Filter::Seg { .. } => {}
                    }
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1357-1363
                // ```c++
                //     if (m_Remote) {
                //         if (value) {
                //             m_Remote->SetValue(eBlastOpt_WindowMaskerTaxId, value);
                //         } else {
                //             m_Remote->ResetValue(CBlast4Field::Get(eBlastOpt_WindowMaskerTaxId));
                //         }
                // ```
                // (`SetWindowMaskerTaxId`, called by the filtering handler before it reads
                // `-window_masker_db`, blast_args.cpp:437-444: the taxid 0 is no parameter.)
                if let Some(taxid) = inputs.window_masker_taxid.filter(|&taxid| taxid != 0) {
                    list.set("WindowMaskerTaxId", Integer(taxid));
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-459
                // ```c++
                //     if (args.Exist(kArgWindowMaskerDatabase) &&
                //         args[kArgWindowMaskerDatabase]) {
                //         const string& stat_file = args[kArgWindowMaskerDatabase].AsString();
                //         ...
                //         opt.SetWindowMaskerDatabase(stat_file.c_str());
                // ```
                if let Some(database) = &inputs.window_masker_db {
                    list.set("WindowMaskerDatabase", Str(database.clone()));
                }
            }
            Class::Matrix => {
                if let Some(name) = &inputs.matrix {
                    list.set("MatrixName", Str(name.clone()));
                }
            }
            Class::Threshold => {
                if let Some(threshold) = inputs.threshold {
                    list.set("WordThreshold", Integer(threshold));
                } else {
                    let current = list.get("WordThreshold").cloned();
                    // The threshold that `blastp-short` records (16) is not the options'
                    // value (the handle's 11) until a word size above 4 sets one: the class
                    // then replaces it with the matrix's suggestion (probe `opt_short_matrix`).
                    let is_default = match current {
                        None => true,
                        Some(_) if inputs.task == "blastp-short" && !threshold_set_by_word_size => {
                            true
                        }
                        Some(ParamValue::Integer(value)) => is_default_word_threshold(kind, value),
                        Some(_) => false,
                    };
                    if is_default {
                        let mut value = suggested_threshold(&matrix);
                        // `BLAST_GetSuggestedThreshold` adds 2 for a translated subject
                        // (blast_options.c:1202-1205: "Covers tblastn, tblastx, psi-tblastn
                        // rpstblastn"), 1 for a translated query.
                        value += match kind {
                            RequestKind::Tblastn
                            | RequestKind::Tblastx
                            | RequestKind::PsiTblastn => 2,
                            RequestKind::Blastx => 1,
                            _ => 0,
                        };
                        list.set("WordThreshold", Integer(value));
                    }
                }
            }
            Class::Hsp => {
                if let Some(limit) = inputs.culling_limit {
                    list.set("Culling", Integer(limit));
                }
                if inputs.subject_besthit {
                    list.set("UseSubjectBestHit", Boolean(true));
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-498
            // ```c++
            //     if (args[kArgWindowSize]) {
            //         opt.SetWindowSize(args[kArgWindowSize].AsInteger());
            //     } else {
            //         int window = -1;
            //         BLAST_GetSuggestedWindowSize(opt.GetProgramType(),
            //                                      opt.GetMatrixName(),
            //                                      &window);
            //         if (window != -1) {
            //             opt.SetWindowSize(window);
            //         }
            // ```
            // The RPS searches have no suggestion (`p9.rps.tk.f6`, `p9.rtn.nq4.f6`).
            Class::Window => {
                if let Some(size) = inputs.window_size {
                    list.set("WindowSize", Integer(size));
                } else if kind != RequestKind::Blastn && !kind.is_rps() {
                    list.set("WindowSize", Integer(suggested_window(&matrix)));
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:2619-2632
            // ```c++
            // void CBlastOptions::SetUseIndex(
            //         bool use_index, const string & index_name,
            //         bool force_index, bool old_style_index )
            // {
            //     ...
            //     if (m_Remote) {
            //         m_Remote->SetValue(eBlastOpt_ForceMbIndex, force_index);
            //         if ( !index_name.empty() ) {
            //             m_Remote->SetValue(eBlastOpt_MbIndexName, index_name.c_str());
            // ```
            // `CMbIndexArgs` sets both (`force_index` true, the `-db` name); the index that
            // does not load sets `SetUseIndex(false)` (api/setup_factory.cpp:455-459), whose
            // defaults make the first FALSE in its place and leave the name.
            Class::MbIndex => {
                if let Some(name) = &inputs.mb_index_name {
                    list.set("ForceMbIndex", Boolean(false));
                    list.set("MbIndexName", Str(name.clone()));
                }
            }
            Class::Gapped => list.set("UngappedMode", Boolean(inputs.ungapped)),
            Class::Intron => list.set("LongestIntronLength", Integer(inputs.max_intron_length)),
            Class::QueryGeneticCode => list.set("QueryGeneticCode", Integer(inputs.query_gencode)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1050-1055
            // ```c++
            //     if (m_Target == eDatabase && args[kArgDbGeneticCode] &&
            //         (program == eTblastn || program == eTblastx) ) {
            //         opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
            //     }
            // ```
            // PSI-TBLASTN (`ePSITblastn`) keeps the constructor's code 1 in its place.
            Class::DbGeneticCode => {
                if kind != RequestKind::PsiTblastn {
                    list.set("DbGeneticCode", Integer(inputs.db_gencode));
                }
            }
            Class::Formatting => list.set("HitlistSize", Integer(inputs.hitlist_size)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:879-885
            // ```c++
            //         opt.SetCompositionBasedStats(compo_mode);
            //         if (program == eBlastp &&
            //             compo_mode != eNoCompositionBasedStats &&
            //             tolower(comp_stat_string[1]) == 'u') {
            //             opt.SetUnifiedP(1);
            //         }
            //         opt.SetSmithWatermanMode(smith_waterman_value);
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:866-868
            // ```c++
            //         if(program == ePSITblastn) {
            //             compo_mode = eNoCompositionBasedStats;
            //         }
            // ```
            // `UnifiedP` is BLASTP's only (psiblast's handle is `ePSIBlast`); `-use_sw_tback` is
            // recorded in every mode, 0 included; PSI-TBLASTN records mode 0 whatever is given.
            Class::Cbs => {
                let mode = if kind == RequestKind::PsiTblastn {
                    0
                } else {
                    inputs.comp_based_stats
                };
                list.set("CompositionBasedStats", Integer(mode));
                if kind == RequestKind::Blastp && inputs.unified_p && inputs.comp_based_stats != 0 {
                    list.set("UnifiedP", Integer(1));
                }
                list.set("SmithWatermanMode", Boolean(inputs.use_sw_tback));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1073-1078
            // ```c++
            // CGapTriggerArgs::ExtractAlgorithmOptions(const CArgs& args,
            //                                          CBlastOptions& opt)
            // {
            //     if (args[kArgGapTrigger]) {
            //         opt.SetGapTrigger(args[kArgGapTrigger].AsDouble());
            //     }
            // ```
            // The argument has a default (22), so it is always set.
            Class::GapTrigger => list.set("GapTrigger", Real(inputs.gap_trigger)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1329-1331
            // ```c++
            //         if (!m_IsDeltaBlast) {
            //             opt.SetIgnoreMsaMaster(args[kArgIgnoreMsaMaster]);
            //         }
            // ```
            // (`-ignore_msa_master` is refused by LOSAT's parser: always FALSE.)
            Class::PsiBlast => list.set("IgnoreMsaMaster", Boolean(false)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1113-1119
            // ```c++
            //     if (args[kArgPSIPseudocount]) {
            //         opt.SetPseudoCount(args[kArgPSIPseudocount].AsInteger());
            //     }
            //
            //     if (args[kArgPSIInclusionEThreshold]) {
            //         opt.SetInclusionThreshold(args[kArgPSIInclusionEThreshold].AsDouble());
            //     }
            // ```
            // Both have defaults: set in the places of the constructor's.
            Class::PssmEngine => {
                list.set("PseudoCountWeight", Integer(inputs.pseudocount));
                list.set("InclusionThreshold", Real(inputs.inclusion_ethresh));
            }
        }
    }
    list.into_params()
}

/// The `Blast4-queue-search-request` service of a run: `megablast` for the megablast and
/// dc-megablast tasks of BLASTN, `psi` for psiblast, `rpsblast` for the RPS searches (the handles'
/// `SetRemoteProgramAndService_Blast3`, `RequestKind`), `plain` otherwise.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp (`x_Process_Query`
/// callers set the service from the options: `blastn` task `megablast` and `dc-megablast`
/// give "megablast", every other task "plain").
pub fn service(kind: impl Into<RequestKind>, task: &str) -> &'static str {
    match kind.into() {
        RequestKind::Blastn if matches!(task, "megablast" | "dc-megablast") => "megablast",
        RequestKind::Psiblast | RequestKind::PsiTblastn => "psi",
        RequestKind::Rpsblast | RequestKind::Rpstblastn => "rpsblast",
        _ => "plain",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(inputs: &OptionInputs) -> Vec<String> {
        algorithm_options(inputs)
            .iter()
            .map(|param| param.name.to_string())
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
        let listed: Vec<&str> = params.iter().map(|param| param.name).collect();
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
        let listed: Vec<&str> = params.iter().map(|param| param.name).collect();
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
        assert_eq!(list[0].name, "EvalueThreshold");
        assert_eq!(list[0].value, ParamValue::Evalue(1e-5));
        assert_eq!(list[1].name, "MatchReward");
        assert_eq!(
            list.iter().filter(|p| p.name == "EvalueThreshold").count(),
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
        // NCBI compares the names with `strcasecmp` (blast_options.c:1174-1236).
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
            .map(|param| format!("{}={:?}", param.name, param.value))
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
        for (kind, program, named) in [
            (RequestKind::Psiblast, "blastp", "psi"),
            (RequestKind::PsiTblastn, "blastp", "psi"),
            (RequestKind::Rpsblast, "blastp", "rpsblast"),
            (RequestKind::Rpstblastn, "tblastn", "rpsblast"),
            (RequestKind::Blastp, "blastp", "plain"),
        ] {
            assert_eq!(kind.program(), program);
            assert_eq!(service(kind, "blastp"), named, "{kind:?}");
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
        let listed: Vec<&str> = params.iter().map(|param| param.name).collect();
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
}
