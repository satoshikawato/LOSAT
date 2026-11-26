//! The request options of a BLAST archive (`algorithm-options` of the `queue-search`
//! request): the parameters that the `Set*` calls of NCBI's application record.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:257-269
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

use super::fields::{self, Field, ValueType};
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
    /// `strand-type` (`Blast4-strand-type`: forward 1, reverse 2, both 3).
    StrandType(i32),
    // NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:512,529
    /// `big-integer-list` (`Blast4-value` choice 12, `SEQUENCE OF BigInt`).
    BigIntegerList(Vec<i64>),
}

impl ParamValue {
    /// The `Blast4-value` alternative of the value.
    pub fn value_type(&self) -> ValueType {
        match self {
            ParamValue::Evalue(_) => ValueType::Cutoff,
            ParamValue::Integer(_) => ValueType::Integer,
            ParamValue::BigInteger(_) => ValueType::BigInteger,
            ParamValue::Boolean(_) => ValueType::Boolean,
            ParamValue::Real(_) => ValueType::Real,
            ParamValue::Str(_) => ValueType::String,
            ParamValue::IntegerList(_) => ValueType::IntegerList,
            ParamValue::StrandType(_) => ValueType::StrandType,
            ParamValue::BigIntegerList(_) => ValueType::BigIntegerList,
        }
    }
}

/// A recorded parameter (`Blast4-parameter`): a field of the registry (`fields`) with a value of
/// the field's type, so that NCBI's `CBlast4Field::Match` takes every parameter LOSAT writes.
#[derive(Clone, Debug, PartialEq)]
pub struct Param {
    field: Field,
    pub value: ParamValue,
}

impl Param {
    pub fn new(field: Field, value: ParamValue) -> Param {
        debug_assert_eq!(
            value.value_type(),
            field.value_type(),
            "the value of {} is not of its field's type",
            field.name()
        );
        Param { field, value }
    }

    /// The parameter's field.
    pub fn field(&self) -> Field {
        self.field
    }

    /// The parameter's name.
    pub fn name(&self) -> &'static str {
        self.field.name()
    }
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

/// The filter of a `-seg` value as given (`app::parse_seg_option`), `None` for a value that NCBI's
/// filtering handler refuses (no request is written then).
pub fn seg_filter(value: &str) -> Option<Filter> {
    crate::blastinput::app::parse_seg_option(value)
        .ok()
        .map(Filter::from)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:73-88
// NCBI reference (598d8ae6): c++/include/algo/blast/api/psiblast_options.hpp:82-84
// NCBI reference (598d8ae6): c++/include/algo/blast/api/blast_rps_options.hpp:122-124
// NCBI reference (598d8ae6): c++/include/algo/blast/api/rpstblastn_options.hpp:74-76
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:138-146
/// The kind of a run's request: its options handle (the `program` and `service` it names, the
/// setters of its constructor) and the argument classes of its application, in their order. The
/// five search programs, and the applications whose handle is another one's: psiblast with a
/// FASTA or a PSSM query (`CPSIBlastOptionsHandle`), tblastn with a PSSM query (PSI-TBLASTN:
/// `CPSIBlastOptionsHandle` with `SetPSITblastnDefaults`, whatever `-task`: the test names
/// `blastx-fast`, which is no tblastn task), psiblast with `-phi_pattern`
/// (`CPHIBlastProtOptionsHandle`, `Phiblastp`), rpsblast (`CBlastRPSOptionsHandle`) and rpstblastn
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
    /// psiblast with `-phi_pattern` (PHI-BLAST, item PS8a).
    Phiblastp,
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
            | RequestKind::Phiblastp
            | RequestKind::Rpsblast => "blastp",
            RequestKind::Blastx => "blastx",
            RequestKind::Tblastn | RequestKind::Rpstblastn => "tblastn",
            RequestKind::Tblastx => "tblastx",
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:418-421
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:120-122
    /// The application's name, which NCBI's diagnostics write as the prefix of its posts
    /// (`[psiblast]`; rpsblast_app.cpp:119, rpstblastn_app.cpp:119, and the program of the five
    /// search programs). rpstblastn's `-mt_mode 1` with several threads posts as `rpstblastn_mt`
    /// (rpstblastn_app.cpp:254), a run whose output LOSAT writes as `-num_threads 1`'s
    /// (DECISIONS F1b-2 of line RG). tblastn's `-in_pssm` (PSI-TBLASTN), whose request names
    /// `blastp`, posts as `tblastn` (NCBI 2.17.0, `tblastn -in_pssm` with a `-db` name holding
    /// `é`: `Error: (803.7) [tblastn] Blast4-request.body.queue-search.subject.database`, probe RG
    /// K/k4b/agent7/probe4 `tbn_dbe_f6`).
    pub fn application(self) -> &'static str {
        match self {
            RequestKind::Psiblast | RequestKind::Phiblastp => "psiblast",
            RequestKind::PsiTblastn => "tblastn",
            RequestKind::Rpsblast => "rpsblast",
            RequestKind::Rpstblastn => "rpstblastn",
            kind => kind.program(),
        }
    }

    /// Whether the handle is an RPS search's.
    fn is_rps(self) -> bool {
        matches!(self, RequestKind::Rpsblast | RequestKind::Rpstblastn)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_program.c:53-54
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_program.h:72-90
    /// Whether the subject sequences are nucleotide (`Blast_SubjectIsNucleotide`):
    /// BLASTN, TBLASTN (including PSI-TBLASTN with `-in_pssm`) and TBLASTX search nucleotide
    /// subjects; BLASTP, PSIBLAST, BLASTX, RPSBLAST and RPSTBLASTN search protein or profile
    /// databases.
    pub fn subject_is_nucleotide(self) -> bool {
        matches!(
            self,
            RequestKind::Blastn
                | RequestKind::Tblastn
                | RequestKind::PsiTblastn
                | RequestKind::Tblastx
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:46-74
    /// The setters that the handle's constructor calls after defaults mode, after the
    /// `EvalueThreshold` that every handle records: psiblast's PSSM engine defaults
    /// (`PSI_INCLUSION_ETHRESH` 0.002, `PSI_PSEUDO_COUNT_CONST` 0), and for PSI-TBLASTN then
    /// `SetPSITblastnDefaults` (`SmithWatermanMode` FALSE, `DbGeneticCode` 1, `WordThreshold`
    /// `BLAST_WORD_THRESHOLD_TBLASTN` 13; `SetProgram` records nothing; NCBI 2.17.0
    /// `tbn.pssm.db.ss` of line RG F1d probe 7).
    fn constructor_calls(self) -> Vec<Param> {
        let psi = || {
            vec![
                Param::new(fields::INCLUSION_THRESHOLD, ParamValue::Real(0.002)),
                Param::new(fields::PSEUDO_COUNT, ParamValue::Integer(0)),
            ]
        };
        match self {
            RequestKind::Psiblast => psi(),
            RequestKind::PsiTblastn => {
                let mut calls = psi();
                calls.extend([
                    Param::new(fields::SMITH_WATERMAN_MODE, ParamValue::Boolean(false)),
                    Param::new(fields::DB_GENETIC_CODE, ParamValue::Integer(1)),
                    Param::new(fields::WORD_THRESHOLD, ParamValue::Integer(13)),
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
    /// `-filtering_db` (BLASTN) as given.
    pub filtering_db: Option<String>,
    /// `-window_masker_db` (BLASTN) as given.
    pub window_masker_db: Option<String>,
    /// `-window_masker_taxid` (BLASTN) as given.
    pub window_masker_taxid: Option<i32>,
    pub threshold: Option<i32>,
    /// `-culling_limit` (also 0), when given.
    pub culling_limit: Option<i32>,
    /// `-best_hit_overhang` and `-best_hit_score_edge`, each when given.
    pub best_hit_overhang: Option<f64>,
    pub best_hit_score_edge: Option<f64>,
    pub subject_besthit: bool,
    pub window_size: Option<i32>,
    pub ungapped: bool,
    /// `-no_greedy` (BLASTN): the dynamic programming algorithms are recorded.
    pub no_greedy: bool,
    /// The name of the MegaBLAST index that a BLASTN run asked for (`-use_index true`, task
    /// megablast: the `-index_name` value, else `-db`; `blastn::blast_engine::mb_index`).
    pub mb_index_name: Option<String>,
    pub max_intron_length: i32,
    pub query_gencode: i32,
    pub db_gencode: i32,
    /// `-max_target_seqs`, else `-num_alignments`, else 500.
    pub hitlist_size: i32,
    pub comp_based_stats: i32,
    /// The `u` of `-comp_based_stats 2u` as given: `SetUnifiedP` set to 1 for BLASTP in an adjusting mode
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
    /// An imported strategy's parameters as its setters record them, in the file's order
    /// (`CBlastOptionsBuilder::GetSearchOptions`): they replace the command line's classes.
    pub imported: Option<Vec<Param>>,
    /// NCBI's name of the running application where it is not the kind's: the application of
    /// a command line that imports another program's strategy, whose posts all carry its prefix
    /// (`ArchiveRequest::application`); `None`: the kind's own (`RequestKind::application`).
    pub application: Option<&'static str>,
    /// psiblast's `-ignore_msa_master` (`CPsiBlastArgs`, always recorded).
    pub ignore_msa_master: bool,
    /// psiblast's `-phi_pattern`: the pattern as the file gives it (`phi_search::file`), of
    /// `RequestKind::Phiblastp`.
    pub phi_pattern: Option<String>,
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
            filtering_db: None,
            window_masker_db: None,
            window_masker_taxid: None,
            threshold: None,
            culling_limit: None,
            best_hit_overhang: None,
            best_hit_score_edge: None,
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
            imported: None,
            application: None,
            ignore_msa_master: false,
            phi_pattern: None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:257-260
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

/// The setters of the task (`CBlastOptionsFactory::Create` of `task` outside defaults mode).
fn task_calls(task: &str) -> Vec<Param> {
    use ParamValue::{Boolean, Evalue, Integer, Str};
    let p = Param::new;
    match task {
        "blastn-short" => vec![
            p(fields::EVALUE_THRESHOLD, Evalue(1000.0)),
            p(fields::MATCH_REWARD, Integer(1)),
            p(fields::MISMATCH_PENALTY, Integer(-3)),
            p(fields::WORD_SIZE, Integer(7)),
            p(fields::DUST_FILTERING, Boolean(false)),
            p(fields::SEG_FILTERING, Boolean(false)),
            p(fields::REPEAT_FILTERING, Boolean(false)),
            // `ClearFilterOptions` ends with `SetMaskAtHash` set to false (api/blast_options_cxx.cpp:1022-1031).
            p(fields::MASK_AT_HASH, Boolean(false)),
        ],
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:365-370
        "rmblastn" => vec![
            p(fields::MATCH_REWARD, Integer(0)),
            p(fields::MISMATCH_PENALTY, Integer(0)),
        ],
        "dc-megablast" => vec![
            p(fields::MB_TEMPLATE_TYPE, Integer(0)),
            p(fields::MB_TEMPLATE_LENGTH, Integer(18)),
            p(fields::WORD_SIZE, Integer(11)),
            p(fields::WINDOW_SIZE, Integer(40)),
        ],
        "blastp-fast" | "tblastn-fast" | "blastx-fast" => {
            vec![
                p(fields::WORD_SIZE, Integer(5)),
                p(fields::WORD_THRESHOLD, Integer(20)),
            ]
        }
        "blastp-short" => vec![
            p(fields::EVALUE_THRESHOLD, Evalue(20000.0)),
            p(fields::MATRIX_NAME, Str("PAM30".to_string())),
            p(fields::GAP_OPENING_COST, Integer(9)),
            p(fields::GAP_EXTENSION_COST, Integer(1)),
            p(fields::WORD_SIZE, Integer(2)),
            p(fields::WORD_THRESHOLD, Integer(16)),
            p(fields::DUST_FILTERING, Boolean(false)),
            p(fields::SEG_FILTERING, Boolean(false)),
            p(fields::REPEAT_FILTERING, Boolean(false)),
            // `ClearFilterOptions` ends with `SetMaskAtHash` set to false (api/blast_options_cxx.cpp:1022-1031).
            p(fields::MASK_AT_HASH, Boolean(false)),
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
    /// psiblast's `CPhiBlastArgs`.
    Phi,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:102-113 (the end of
// psiblast's `m_Args` pushes; rpsblast_args.cpp:70-127 and rpstblastn_args.cpp:45-118 for the RPS
// programs)
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:102-116
        // psiblast's classes with the PHI-BLAST class last; the composition-based statistics
        // class records nothing for `ePHIBlastp` (`Class::Cbs`; `fi.export.it1`, NCBI 2.17.0).
        RequestKind::Phiblastp => &[
            Db, Generic, Filtering, Matrix, Threshold, Hsp, Window, Formatting, GapTrigger,
            PsiBlast, PssmEngine, Phi,
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
    /// A `Set*` call: the parameter of the same name is replaced (`x_AttachValue` compares names).
    pub fn set(&mut self, field: Field, value: ParamValue) {
        let param = Param::new(field, value);
        match self
            .params
            .iter_mut()
            .find(|recorded| recorded.name() == field.name())
        {
            Some(recorded) => *recorded = param,
            None => self.params.push(param),
        }
    }

    fn get(&self, field: Field) -> Option<&ParamValue> {
        self.params
            .iter()
            .find(|param| param.name() == field.name())
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
    list.set(fields::EVALUE_THRESHOLD, Evalue(10.0));
    for param in kind.constructor_calls() {
        list.set(param.field(), param.value);
    }
    // PSI-TBLASTN's handle is not the task's (`RequestKind`: tblastn_args.cpp:138-149).
    if kind != RequestKind::PsiTblastn {
        for param in task_calls(&inputs.task) {
            list.set(param.field(), param.value);
        }
    }
    if let Some(imported) = &inputs.imported {
        return replayed_options(list, imported, inputs);
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
            | RequestKind::Phiblastp
            | RequestKind::PsiTblastn
    );
    let mut threshold_set_by_word_size = false;
    for &class in class_order(kind) {
        match class {
            Class::Db => {
                if let Some(size) = inputs.dbsize {
                    list.set(fields::DB_LENGTH, BigInteger(size));
                }
            }
            Class::Generic => {
                if let Some(evalue) = inputs.evalue {
                    list.set(fields::EVALUE_THRESHOLD, Evalue(evalue));
                }
                if let Some(open) = inputs.gapopen {
                    list.set(fields::GAP_OPENING_COST, Integer(open));
                } else if inputs.matrix.is_some() {
                    list.set(fields::GAP_OPENING_COST, Integer(matrix_gaps(&matrix).0));
                }
                if let Some(extend) = inputs.gapextend {
                    list.set(fields::GAP_EXTENSION_COST, Integer(extend));
                } else if inputs.matrix.is_some() {
                    list.set(fields::GAP_EXTENSION_COST, Integer(matrix_gaps(&matrix).1));
                }
                if let Some(xdrop) = inputs.xdrop_gap {
                    list.set(fields::GAP_X_DROPOFF, Real(xdrop));
                }
                if let Some(xdrop) = inputs.xdrop_gap_final {
                    list.set(fields::GAP_X_DROPOFF_FINAL, Real(xdrop));
                }
                if let Some(word_size) = inputs.word_size {
                    if protein_query && word_size > 4 {
                        threshold_set_by_word_size = true;
                        list.set(fields::WORD_THRESHOLD, Integer(19));
                        if word_size > 5 {
                            list.set(fields::WORD_THRESHOLD, Integer(21));
                        }
                        if word_size > 6 {
                            list.set(fields::WORD_THRESHOLD, Integer(20));
                        }
                    }
                    list.set(fields::WORD_SIZE, Integer(word_size));
                }
                if let Some(searchsp) = inputs.searchsp {
                    list.set(fields::EFFECTIVE_SEARCH_SPACE, BigInteger(searchsp));
                }
                if let Some(identity) = inputs.perc_identity {
                    list.set(fields::PERCENT_IDENTITY, Real(identity));
                }
                if let Some(coverage) = inputs.qcov_hsp_perc {
                    list.set(fields::QUERY_COV_HSP_PERC, Real(coverage));
                }
                if let Some(max_hsps) = inputs.max_hsps {
                    list.set(fields::MAX_HSPS_PER_SUBJECT, Integer(max_hsps));
                }
                if let Some(sum_stats) = inputs.sum_stats {
                    list.set(fields::SUM_STATISTICS_MODE, Boolean(sum_stats));
                }
            }
            Class::Nucl => {
                if let Some(penalty) = inputs.penalty {
                    list.set(fields::MISMATCH_PENALTY, Integer(penalty));
                }
                if let Some(reward) = inputs.reward {
                    list.set(fields::MATCH_REWARD, Integer(reward));
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:681-685
                // Both enumerators are 0 (blast_options.h: eDynProgScoreOnly, eDynProgTbck).
                if inputs.no_greedy {
                    list.set(fields::GAP_EXTN_ALGORITHM, Integer(0));
                    list.set(fields::GAP_TRACEBACK_ALGORITHM, Integer(0));
                }
            }
            Class::Dmb => {
                if let Some(kind) = inputs.template_type {
                    list.set(fields::MB_TEMPLATE_TYPE, Integer(kind));
                }
                if let Some(length) = inputs.template_length {
                    list.set(fields::MB_TEMPLATE_LENGTH, Integer(length));
                }
            }
            Class::Filtering => {
                list.set(
                    fields::MASK_AT_HASH,
                    Boolean(inputs.soft_masking.unwrap_or(kind == RequestKind::Blastn)),
                );
                if kind != RequestKind::Blastn {
                    // SEG is on by default for the translated searches only (`kFilterByDefault`).
                    let seg = inputs.seg.clone().unwrap_or(
                        if matches!(
                            kind,
                            RequestKind::Blastp
                                | RequestKind::Psiblast
                                | RequestKind::Phiblastp
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
                        Filter::No => list.set(fields::SEG_FILTERING, Boolean(false)),
                        Filter::Yes => list.set(fields::SEG_FILTERING, Boolean(true)),
                        Filter::Seg {
                            window,
                            locut,
                            hicut,
                        } => {
                            list.set(fields::SEG_FILTERING_WINDOW, Integer(window));
                            list.set(fields::SEG_FILTERING_LOCUT, Real(locut));
                            list.set(fields::SEG_FILTERING_HICUT, Real(hicut));
                        }
                        Filter::Dust { .. } => {}
                    }
                } else if let Some(dust) = inputs.dust.clone() {
                    match dust {
                        Filter::No => list.set(fields::DUST_FILTERING, Boolean(false)),
                        Filter::Yes => list.set(fields::DUST_FILTERING, Boolean(true)),
                        Filter::Dust {
                            level,
                            window,
                            linker,
                        } => {
                            list.set(fields::DUST_FILTERING_LEVEL, Integer(level));
                            list.set(fields::DUST_FILTERING_WINDOW, Integer(window));
                            list.set(fields::DUST_FILTERING_LINKER, Integer(linker));
                        }
                        Filter::Seg { .. } => {}
                    }
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:432-435
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1332-1339
                // (`SetRepeatFilteringDB`: the name as given, an empty one too; the handler sets
                // no `RepeatFiltering` boolean.)
                if let Some(database) = &inputs.filtering_db {
                    list.set(fields::REPEAT_FILTERING_DB, Str(database.clone()));
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:1357-1363
                // (`SetWindowMaskerTaxId`, called by the filtering handler before it reads
                // `-window_masker_db`, blast_args.cpp:437-444: the taxid 0 is no parameter.)
                if let Some(taxid) = inputs.window_masker_taxid.filter(|&taxid| taxid != 0) {
                    list.set(fields::WINDOW_MASKER_TAX_ID, Integer(taxid));
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-459
                if let Some(database) = &inputs.window_masker_db {
                    list.set(fields::WINDOW_MASKER_DATABASE, Str(database.clone()));
                }
            }
            Class::Matrix => {
                if let Some(name) = &inputs.matrix {
                    list.set(fields::MATRIX_NAME, Str(name.clone()));
                }
            }
            Class::Threshold => {
                if let Some(threshold) = inputs.threshold {
                    list.set(fields::WORD_THRESHOLD, Integer(threshold));
                } else {
                    let current = list.get(fields::WORD_THRESHOLD).cloned();
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
                        list.set(fields::WORD_THRESHOLD, Integer(value));
                    }
                }
            }
            Class::Hsp => {
                if let Some(limit) = inputs.culling_limit {
                    list.set(fields::CULLING_LIMIT, Integer(limit));
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3338-3349
                // Each when given, after `Culling` (they exclude each other) and before
                // `UseSubjectBestHit` (NCBI BLAST+ 2.17.0, line NP O6a's probe).
                if let Some(overhang) = inputs.best_hit_overhang {
                    list.set(fields::BEST_HIT_OVERHANG, Real(overhang));
                }
                if let Some(score_edge) = inputs.best_hit_score_edge {
                    list.set(fields::BEST_HIT_SCORE_EDGE, Real(score_edge));
                }
                if inputs.subject_besthit {
                    list.set(fields::SUBJECT_BEST_HIT, Boolean(true));
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:485-498
            // The RPS searches have no suggestion (`p9.rps.tk.f6`, `p9.rtn.nq4.f6`).
            Class::Window => {
                if let Some(size) = inputs.window_size {
                    list.set(fields::WINDOW_SIZE, Integer(size));
                } else if kind != RequestKind::Blastn && !kind.is_rps() {
                    list.set(fields::WINDOW_SIZE, Integer(suggested_window(&matrix)));
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:2619-2632
            // `CMbIndexArgs` sets both (`force_index` true, the `-db` name); the index that
            // does not load calls `SetUseIndex` with false (api/setup_factory.cpp:455-459), whose
            // defaults make the first FALSE in its place and leave the name. An empty name
            // (`-index_name ''`) is not set. A search with the index sets it TRUE again in the
            // run's options handle (`blastn::blast_engine::mb_index::note_index_use`).
            Class::MbIndex => {
                if let Some(name) = &inputs.mb_index_name {
                    list.set(fields::FORCE_MB_INDEX, Boolean(false));
                    if !name.is_empty() {
                        list.set(fields::MB_INDEX_NAME, Str(name.clone()));
                    }
                }
            }
            Class::Gapped => list.set(fields::GAPPED_MODE, Boolean(inputs.ungapped)),
            Class::Intron => list.set(
                fields::LONGEST_INTRON_LENGTH,
                Integer(inputs.max_intron_length),
            ),
            Class::QueryGeneticCode => {
                list.set(fields::QUERY_GENETIC_CODE, Integer(inputs.query_gencode))
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1050-1055
            // PSI-TBLASTN (`ePSITblastn`) keeps the constructor's code 1 in its place.
            Class::DbGeneticCode => {
                if kind != RequestKind::PsiTblastn {
                    list.set(fields::DB_GENETIC_CODE, Integer(inputs.db_gencode));
                }
            }
            Class::Formatting => list.set(fields::HITLIST_SIZE, Integer(inputs.hitlist_size)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:879-885
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:866-869
            // `UnifiedP` is BLASTP's only (psiblast's handle is `ePSIBlast`); `-use_sw_tback` is
            // recorded in every mode, 0 included; PSI-TBLASTN records mode 0 whatever is given.
            Class::Cbs => {
                let mode = if kind == RequestKind::PsiTblastn {
                    0
                } else {
                    inputs.comp_based_stats
                };
                list.set(fields::COMPOSITION_BASED_STATS, Integer(mode));
                if kind == RequestKind::Blastp && inputs.unified_p && inputs.comp_based_stats != 0 {
                    list.set(fields::UNIFIED_P, Integer(1));
                }
                list.set(fields::SMITH_WATERMAN_MODE, Boolean(inputs.use_sw_tback));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1073-1078
            // The argument has a default (22), so it is always set.
            Class::GapTrigger => list.set(fields::GAP_TRIGGER, Real(inputs.gap_trigger)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1329-1331
            Class::PsiBlast => {
                list.set(fields::IGNORE_MSA_MASTER, Boolean(inputs.ignore_msa_master))
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1113-1119
            // Both have defaults: set in the places of the constructor's.
            Class::PssmEngine => {
                list.set(fields::PSEUDO_COUNT, Integer(inputs.pseudocount));
                list.set(fields::INCLUSION_THRESHOLD, Real(inputs.inclusion_ethresh));
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1382-1411
            // The pattern as read (`SetPHIPattern`).
            Class::Phi => {
                if let Some(pattern) = &inputs.phi_pattern {
                    list.set(fields::PHI_PATTERN, Str(pattern.clone()));
                }
            }
        }
    }
    list.into_params()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3811-3828
/// The options of an imported strategy: the handle `CBlastOptionsFactory::Create` makes for the
/// strategy's program (`list`), each setter of the file replayed over it in the file's order,
/// then the classes `SetOptionsForSavedStrategy` extracts again: `CFormattingArgs` sets the
/// hitlist size, `CBlastDatabaseArgs` the command line's `-dbsize` (`dbsize` holds it, or the
/// strategy's value, which then keeps its place) and `CMbIndexArgs` the index.
fn replayed_options(mut list: ParamList, imported: &[Param], inputs: &OptionInputs) -> Vec<Param> {
    for param in imported {
        list.set(param.field(), param.value.clone());
    }
    list.set(
        fields::HITLIST_SIZE,
        ParamValue::Integer(inputs.hitlist_size),
    );
    if let Some(size) = inputs.dbsize {
        list.set(fields::DB_LENGTH, ParamValue::BigInteger(size));
    }
    if let Some(name) = &inputs.mb_index_name {
        list.set(fields::FORCE_MB_INDEX, ParamValue::Boolean(false));
        if !name.is_empty() {
            list.set(fields::MB_INDEX_NAME, ParamValue::Str(name.clone()));
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
#[path = "options_tests.rs"]
mod tests;
