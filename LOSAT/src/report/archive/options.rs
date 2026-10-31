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
//! line); (4) after each batch's search `EffectiveSearchSpace` (appended last by the writer).
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

/// Every option whose value or presence an archive records, as `Option`s that each program
/// fills from its resolved arguments (`None` = not given).
#[derive(Clone, Debug)]
pub struct OptionInputs {
    pub program: Program,
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
    pub threshold: Option<i32>,
    /// `-culling_limit` (also 0), when given.
    pub culling_limit: Option<i32>,
    pub subject_besthit: bool,
    pub window_size: Option<i32>,
    pub ungapped: bool,
    /// `-no_greedy` (BLASTN): the dynamic programming algorithms are recorded.
    pub no_greedy: bool,
    pub max_intron_length: i32,
    pub query_gencode: i32,
    pub db_gencode: i32,
    /// `-max_target_seqs`, else `-num_alignments`, else 500.
    pub hitlist_size: i32,
    pub comp_based_stats: i32,
    /// `-comp_based_stats 2u` (BLASTP, an adjusting mode): `SetUnifiedP(1)`.
    pub unified_p: bool,
    pub use_sw_tback: bool,
}

impl OptionInputs {
    /// The inputs of a run that gives no option: every field at its default.
    pub fn new(program: Program, task: &str) -> Self {
        Self {
            program,
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
            threshold: None,
            culling_limit: None,
            subject_besthit: false,
            window_size: None,
            ungapped: false,
            no_greedy: false,
            max_intron_length: 0,
            query_gencode: 1,
            db_gencode: 1,
            hitlist_size: 500,
            comp_based_stats: 2,
            unified_p: false,
            use_sw_tback: false,
        }
    }
}

/// The `(gap open, gap extend)` that `-matrix` implies (`BLAST_GetProteinGapExistenceExtendParams`:
/// the matrix's default row, which BLOSUM45 puts at 14 / 2).
fn matrix_gaps(matrix: &str) -> (i32, i32) {
    crate::stats::protein_options::protein_gap_existence_extend_params(matrix).unwrap_or((11, 1))
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

fn default_threshold(program: Program) -> i32 {
    match program {
        Program::Blastp => 11,
        Program::Blastx => 12,
        _ => 13,
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
    Formatting,
    Matrix,
    Threshold,
    Intron,
    QueryGeneticCode,
    DbGeneticCode,
    Cbs,
}

/// The argument classes of a program in the order of `<program>_args.cpp` (the order of
/// the `m_Args` pushes).
fn class_order(program: Program) -> &'static [Class] {
    use Class::*;
    match program {
        Program::Blastn => &[
            Db, Generic, Nucl, Dmb, Filtering, Gapped, Hsp, Window, Formatting,
        ],
        Program::Blastp => &[
            Db, Generic, Filtering, Matrix, Threshold, Hsp, Window, Formatting, Gapped, Cbs,
        ],
        Program::Blastx => &[
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
        Program::Tblastn => &[
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
        Program::Tblastx => &[
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
#[derive(Clone, Debug, Default)]
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

    pub fn into_params(self) -> Vec<Param> {
        self.params
    }
}

/// The recorded `algorithm-options` of a run, without `EffectiveSearchSpace` (which each
/// archive appends).
pub fn algorithm_options(inputs: &OptionInputs) -> Vec<Param> {
    use ParamValue::{BigInteger, Boolean, Evalue, Integer, Real, Str};
    let program = inputs.program;
    let mut list = ParamList::default();
    list.set("EvalueThreshold", Evalue(10.0));
    for param in task_calls(&inputs.task) {
        list.set(param.name, param.value);
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
        program,
        Program::Blastp | Program::Tblastn | Program::Blastx
    );
    let mut threshold_set_by_word_size = false;
    for &class in class_order(program) {
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
                    Boolean(inputs.soft_masking.unwrap_or(program == Program::Blastn)),
                );
                if program != Program::Blastn {
                    let seg = inputs.seg.clone().unwrap_or(if program == Program::Blastp {
                        Filter::No
                    } else {
                        Filter::Seg {
                            window: 12,
                            locut: 2.2,
                            hicut: 2.5,
                        }
                    });
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
                        Some(ParamValue::Integer(value)) => value == default_threshold(program),
                        Some(_) => false,
                    };
                    if is_default {
                        let mut value = suggested_threshold(&matrix);
                        value += match program {
                            Program::Tblastn | Program::Tblastx => 2,
                            Program::Blastx => 1,
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
            Class::Window => {
                if let Some(size) = inputs.window_size {
                    list.set("WindowSize", Integer(size));
                } else if program != Program::Blastn {
                    list.set("WindowSize", Integer(suggested_window(&matrix)));
                }
            }
            Class::Gapped => list.set("UngappedMode", Boolean(inputs.ungapped)),
            Class::Intron => list.set("LongestIntronLength", Integer(inputs.max_intron_length)),
            Class::QueryGeneticCode => list.set("QueryGeneticCode", Integer(inputs.query_gencode)),
            Class::DbGeneticCode => list.set("DbGeneticCode", Integer(inputs.db_gencode)),
            Class::Formatting => list.set("HitlistSize", Integer(inputs.hitlist_size)),
            Class::Cbs => {
                list.set("CompositionBasedStats", Integer(inputs.comp_based_stats));
                if inputs.unified_p && inputs.comp_based_stats != 0 {
                    list.set("UnifiedP", Integer(1));
                }
                list.set("SmithWatermanMode", Boolean(inputs.use_sw_tback));
            }
        }
    }
    list.into_params()
}

/// The `Blast4-queue-search-request` service of a run: `megablast` for the megablast and
/// dc-megablast tasks of BLASTN, `plain` otherwise.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp (`x_Process_Query`
/// callers set the service from the options: `blastn` task `megablast` and `dc-megablast`
/// give "megablast", every other task "plain").
pub fn service(program: Program, task: &str) -> &'static str {
    if program == Program::Blastn && matches!(task, "megablast" | "dc-megablast") {
        "megablast"
    } else {
        "plain"
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
}
