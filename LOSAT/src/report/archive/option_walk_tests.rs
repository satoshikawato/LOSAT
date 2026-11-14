//! Every argument that the command line of psiblast, rpsblast and rpstblastn takes, against what
//! the search strategy's requests record of it (line RG F1b-1, review 3 R3-3): each application's
//! options (`option_inputs` of its own arguments, then `options::algorithm_options` of its
//! `RequestKind`) and the request's `HitlistSize` (`strategy::hitlist_size`) and `Web_StepNumber`.
//! The expectations follow NCBI's argument classes (`psiblast_args.cpp`, `rpsblast_args.cpp`,
//! `rpstblastn_args.cpp` and their `ExtractAlgorithmOptions` in `blast_args.cpp`) and the rows of
//! NCBI BLAST+ 2.17.0 of probe 9 (`F/f1b/p9_docs.txt`) and of review 3: an argument is recorded
//! where its class records it, or not at all. An argument of the command line without an
//! expectation fails the test, so a new argument needs its NCBI rule here.

use super::options::OptionInputs;
use super::strategy::{hitlist_size, HitlistArgs};
use super::ArchiveRequest;
use crate::cli::{Cli, Commands};

/// A change of the default request: the parameter in its place, a new one after another, or a
/// parameter gone.
#[derive(Clone, Copy, Debug)]
enum Edit {
    At(&'static str, &'static str),
    After(&'static str, &'static str, &'static str),
    Drop(&'static str),
}

use Edit::{After, At, Drop};

/// One run of the walk: the argument's id (clap's), the words added to the base command line,
/// and what changes in the request.
type Case = (&'static str, &'static [&'static str], &'static [Edit]);

struct App {
    program: &'static str,
    /// The base command line after the program's name.
    base: &'static [&'static str],
    /// The request's options of the base command line, `Name=Value` (`ParamValue`'s `Debug`).
    default: &'static [&'static str],
    cases: &'static [Case],
    /// The application's options and its `-db` from a parsed command line.
    inputs: fn(Commands) -> (OptionInputs, Option<String>),
}

fn parse(program: &str, words: &[&str]) -> Result<Commands, String> {
    let mut argv = vec!["losat", program];
    argv.extend_from_slice(words);
    crate::cli::try_parse_from::<Cli, _, _>(argv)
        .map(|cli| cli.command)
        .map_err(|error| crate::cli::render_message(&error))
}

/// The request's options, `HitlistSize` included, and `Web_StepNumber` last.
fn render(inputs: &OptionInputs, database: Option<&str>) -> Vec<String> {
    let request = ArchiveRequest::for_inputs(inputs, database);
    let mut listed: Vec<String> = request
        .algorithm
        .iter()
        .map(|param| format!("{}={:?}", param.name, param.value))
        .collect();
    if let Some(step) = request.step_number {
        listed.push(format!("Web_StepNumber=Integer({step})"));
    }
    listed
}

fn expected(default: &[&str], edits: &[Edit]) -> Vec<String> {
    let mut listed: Vec<String> = default.iter().map(|item| item.to_string()).collect();
    let at = |listed: &[String], name: &str| {
        listed
            .iter()
            .position(|item| item.split('=').next() == Some(name))
            .unwrap_or_else(|| panic!("{name} is not in {listed:?}"))
    };
    for edit in edits {
        match *edit {
            At(name, value) => {
                let index = at(&listed, name);
                listed[index] = format!("{name}={value}");
            }
            After(anchor, name, value) => {
                let index = at(&listed, anchor);
                listed.insert(index + 1, format!("{name}={value}"));
            }
            Drop(name) => {
                let index = at(&listed, name);
                listed.remove(index);
            }
        }
    }
    listed
}

fn walk(app: &App) {
    // Every argument of the program's command line has a case.
    let command = <Cli as clap::CommandFactory>::command();
    let subcommand = command
        .find_subcommand(app.program)
        .expect("the program's subcommand");
    let missing: Vec<String> = subcommand
        .get_arguments()
        .map(|arg| arg.get_id().as_str().to_string())
        .filter(|id| !matches!(id.as_str(), "help" | "ncbi_help"))
        .filter(|id| !app.cases.iter().any(|case| case.0 == id))
        .collect();
    assert!(
        missing.is_empty(),
        "{}: arguments without an expectation: {missing:?}",
        app.program
    );
    let run = |words: &[&str]| {
        let command = parse(app.program, words)
            .unwrap_or_else(|error| panic!("{} {words:?}: {error}", app.program));
        let (inputs, database) = (app.inputs)(command);
        render(&inputs, database.as_deref())
    };
    assert_eq!(run(app.base), app.default, "{} default", app.program);
    for &(id, words, edits) in app.cases {
        // `-subject` takes the place of `-db`, `-in_pssm` of `-query`.
        let mut line: Vec<&str> = Vec::new();
        let mut base = app.base.iter();
        while let Some(&word) = base.next() {
            let replaced = (word == "-db" && words.contains(&"-subject"))
                || (word == "-query" && words.contains(&"-in_pssm"));
            if replaced {
                base.next();
            } else {
                line.push(word);
            }
        }
        line.extend_from_slice(words);
        assert_eq!(
            run(&line),
            expected(app.default, edits),
            "{} -{id}: {words:?}",
            app.program
        );
    }
}

fn hitlist(
    outfmt: &str,
    max_target_seqs: Option<usize>,
    num_descriptions: Option<usize>,
    num_alignments: Option<usize>,
) -> i32 {
    hitlist_size(
        crate::blastinput::app::parse_formatting_string(outfmt)
            .expect("outfmt")
            .number,
        HitlistArgs {
            max_target_seqs,
            num_descriptions,
            num_alignments,
        },
    )
}

/// The database lists (`DbFilterArgs`): program options of the database (`OptionsHandle`), never
/// algorithm options.
const DB_FILTERS: &[Case] = &[
    ("seqidlist", &["-seqidlist", "s.txt"], &[]),
    ("negative_seqidlist", &["-negative_seqidlist", "s.txt"], &[]),
    ("taxids", &["-taxids", "9606"], &[]),
    ("negative_taxids", &["-negative_taxids", "9606"], &[]),
    ("taxidlist", &["-taxidlist", "t.txt"], &[]),
    ("negative_taxidlist", &["-negative_taxidlist", "t.txt"], &[]),
    (
        "no_taxid_expansion",
        &["-taxids", "9606", "-no_taxid_expansion"],
        &[],
    ),
];

/// The formatting arguments (`FormattingArgs`): the report's only.
const FORMATTING: &[Case] = &[
    ("line_length", &["-line_length", "80"], &[]),
    ("sort_hits", &["-sorthits", "1"], &[]),
    ("sort_hsps", &["-sorthsps", "1"], &[]),
];

fn with(cases: &[&'static [Case]]) -> &'static [Case] {
    Box::leak(cases.concat().into_boxed_slice())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:44-121 (the classes
// and their order, `options::class_order`); the defaults of NCBI BLAST+ 2.17.0 `p9.ps.q.db.itdef.f6`.
#[test]
fn every_psiblast_argument_is_recorded_as_ncbi_records_it() {
    const PSIBLAST: &[Case] = &[
        ("query", &[], &[]),
        ("subject", &["-subject", "s.faa"], &[Drop("Web_StepNumber")]),
        ("db", &[], &[]),
        (
            "dbsize",
            &["-dbsize", "1000000"],
            &[After(
                "PseudoCountWeight",
                "DbLength",
                "BigInteger(1000000)",
            )],
        ),
        ("show_gis", &["-show_gis"], &[]),
        ("parse_deflines", &["-parse_deflines"], &[]),
        (
            "searchsp",
            &["-searchsp", "1000000"],
            &[After(
                "PseudoCountWeight",
                "EffectiveSearchSpace",
                "BigInteger(1000000)",
            )],
        ),
        (
            "sum_stats",
            &["-sum_stats", "true"],
            &[After("PseudoCountWeight", "SumStatistics", "Boolean(true)")],
        ),
        (
            "evalue",
            &["-evalue", "0.001"],
            &[At("EvalueThreshold", "Evalue(0.001)")],
        ),
        (
            "threshold",
            &["-threshold", "12"],
            &[At("WordThreshold", "Integer(12)")],
        ),
        (
            "word_size",
            &["-word_size", "2"],
            &[After("PseudoCountWeight", "WordSize", "Integer(2)")],
        ),
        ("num_threads", &["-num_threads", "2"], &[]),
        ("out", &["-out", "o.txt"], &[]),
        (
            "export_search_strategy",
            &["-export_search_strategy", "s.asn"],
            &[],
        ),
        (
            "max_target_seqs",
            &["-max_target_seqs", "5"],
            &[At("HitlistSize", "Integer(5)")],
        ),
        (
            "num_descriptions",
            &["-num_descriptions", "3", "-num_alignments", "2"],
            &[At("HitlistSize", "Integer(3)")],
        ),
        (
            "num_alignments",
            &["-num_alignments", "600"],
            &[At("HitlistSize", "Integer(600)")],
        ),
        (
            "max_hsps_per_subject",
            &["-max_hsps", "3"],
            &[After(
                "PseudoCountWeight",
                "MaxHspsPerSubject",
                "Integer(3)",
            )],
        ),
        (
            "qcov_hsp_perc",
            &["-qcov_hsp_perc", "10"],
            &[After("PseudoCountWeight", "QueryCovHspPerc", "Real(10.0)")],
        ),
        (
            "culling_limit",
            &["-culling_limit", "2"],
            &[After("WordThreshold", "Culling", "Integer(2)")],
        ),
        // `SetXDropoff` has no remote option (blast_options_cxx.cpp:588-590).
        ("xdrop_ungap", &["-xdrop_ungap", "8"], &[]),
        (
            "xdrop_gap",
            &["-xdrop_gap", "16"],
            &[After("PseudoCountWeight", "GapXDropoff", "Real(16.0)")],
        ),
        (
            "xdrop_gap_final",
            &["-xdrop_gap_final", "26"],
            &[After("PseudoCountWeight", "GapXDropoffFinal", "Real(26.0)")],
        ),
        (
            "gap_trigger",
            &["-gap_trigger", "20"],
            &[At("GapTrigger", "Real(20.0)")],
        ),
        (
            "window_size",
            &["-window_size", "30"],
            &[At("WindowSize", "Integer(30)")],
        ),
        // The matrix's gap costs (generic class), its name (matrix class), its threshold and window
        // suggestions in place (BLOSUM45: 14/2, `matrix_gaps`; `p9.ps.q.db.opts4.f6`: 14 and 60).
        (
            "matrix",
            &["-matrix", "BLOSUM45"],
            &[
                After("PseudoCountWeight", "GapOpeningCost", "Integer(14)"),
                After("GapOpeningCost", "GapExtensionCost", "Integer(2)"),
                After("SegFiltering", "MatrixName", "Str(\"BLOSUM45\")"),
                At("WordThreshold", "Integer(14)"),
                At("WindowSize", "Integer(60)"),
            ],
        ),
        (
            "gap_open",
            &["-gapopen", "9"],
            &[After("PseudoCountWeight", "GapOpeningCost", "Integer(9)")],
        ),
        (
            "gap_extend",
            &["-gapextend", "2"],
            &[After("PseudoCountWeight", "GapExtensionCost", "Integer(2)")],
        ),
        // No `UnifiedP`: NCBI sets it for `eBlastp` only (blast_args.cpp:880-884; review 3 R3-1).
        (
            "comp_based_stats",
            &["-comp_based_stats", "3u"],
            &[At("CompositionBasedStats", "Integer(3)")],
        ),
        (
            "comp_based_stats",
            &["-comp_based_stats", "0"],
            &[At("CompositionBasedStats", "Integer(0)")],
        ),
        (
            "seg",
            &["-seg", "yes"],
            &[At("SegFiltering", "Boolean(true)")],
        ),
        (
            "seg",
            &["-seg", "1 1.5 1.8"],
            &[
                After("MaskAtHash", "SegFilteringWindow", "Integer(1)"),
                After("SegFilteringWindow", "SegFilteringLocut", "Real(1.5)"),
                After("SegFilteringLocut", "SegFilteringHicut", "Real(1.8)"),
                Drop("SegFiltering"),
            ],
        ),
        (
            "soft_masking",
            &["-soft_masking", "true"],
            &[At("MaskAtHash", "Boolean(true)")],
        ),
        ("lcase_masking", &["-lcase_masking"], &[]),
        (
            "use_sw_tback",
            &["-use_sw_tback"],
            &[At("SmithWatermanMode", "Boolean(true)")],
        ),
        ("query_loc", &["-query_loc", "1-10"], &[]),
        (
            "subject_loc",
            &["-subject", "s.faa", "-subject_loc", "1-10"],
            &[Drop("Web_StepNumber")],
        ),
        ("outfmt", &["-outfmt", "6"], &[]),
        (
            "num_iterations",
            &["-num_iterations", "2"],
            &[At("Web_StepNumber", "Integer(2)")],
        ),
        (
            "num_iterations",
            &["-num_iterations", "0"],
            &[Drop("Web_StepNumber")],
        ),
        (
            "pseudocount",
            &["-pseudocount", "2"],
            &[At("PseudoCountWeight", "Integer(2)")],
        ),
        (
            "inclusion_ethresh",
            &["-inclusion_ethresh", "0.01"],
            &[At("InclusionThreshold", "Real(0.01)")],
        ),
        ("out_pssm", &["-out_pssm", "p.asn"], &[]),
        ("out_ascii_pssm", &["-out_ascii_pssm", "p.txt"], &[]),
        (
            "save_pssm_after_last_round",
            &["-save_pssm_after_last_round"],
            &[],
        ),
        ("save_each_pssm", &["-save_each_pssm"], &[]),
        // A PSSM query's request is item F1b-2's (refused with the strategy until then).
        ("in_pssm", &["-in_pssm", "p.asn"], &[]),
    ];
    walk(&App {
        program: "psiblast",
        base: &["-query", "q.faa", "-db", "d"],
        default: &[
            "EvalueThreshold=Evalue(10.0)",
            "InclusionThreshold=Real(0.002)",
            "PseudoCountWeight=Integer(0)",
            "MaskAtHash=Boolean(false)",
            "SegFiltering=Boolean(false)",
            "WordThreshold=Integer(11)",
            "WindowSize=Integer(40)",
            "HitlistSize=Integer(500)",
            "CompositionBasedStats=Integer(2)",
            "SmithWatermanMode=Boolean(false)",
            "GapTrigger=Real(22.0)",
            "IgnoreMsaMaster=Boolean(false)",
            "Web_StepNumber=Integer(1)",
        ],
        cases: with(&[PSIBLAST, DB_FILTERS, FORMATTING]),
        inputs: |command| match command {
            Commands::Psiblast(args) => {
                let mut inputs = args.option_inputs();
                inputs.hitlist_size = hitlist(
                    &args.outfmt,
                    args.max_target_seqs,
                    args.num_descriptions,
                    args.num_alignments,
                );
                (inputs, args.db.clone())
            }
            _ => unreachable!(),
        },
    });
}

/// The RPS programs' generic and window arguments (`CGenericSearchArgs` of an RPS search: no word
/// size, threshold or gap costs; `CWindowSizeArg` without a suggestion), after `EvalueThreshold`.
const RPS_COMMON: &[Case] = &[
    ("query", &[], &[]),
    ("db", &[], &[]),
    ("out", &["-out", "o.txt"], &[]),
    (
        "export_search_strategy",
        &["-export_search_strategy", "s.asn"],
        &[],
    ),
    (
        "evalue",
        &["-evalue", "0.01"],
        &[At("EvalueThreshold", "Evalue(0.01)")],
    ),
    (
        "dbsize",
        &["-dbsize", "100000"],
        &[After("EvalueThreshold", "DbLength", "BigInteger(100000)")],
    ),
    (
        "searchsp",
        &["-searchsp", "1000000"],
        &[After(
            "EvalueThreshold",
            "EffectiveSearchSpace",
            "BigInteger(1000000)",
        )],
    ),
    (
        "comp_based_stats",
        &["-comp_based_stats", "0"],
        &[At("CompositionBasedStats", "Integer(0)")],
    ),
    ("outfmt", &["-outfmt", "6"], &[]),
    ("show_gis", &["-show_gis"], &[]),
    (
        "num_descriptions",
        &["-num_descriptions", "3", "-num_alignments", "2"],
        &[At("HitlistSize", "Integer(3)")],
    ),
    (
        "num_alignments",
        &["-num_alignments", "600"],
        &[At("HitlistSize", "Integer(600)")],
    ),
    (
        "max_target_seqs",
        &["-max_target_seqs", "3"],
        &[At("HitlistSize", "Integer(3)")],
    ),
    (
        "max_hsps_per_subject",
        &["-max_hsps", "2"],
        &[After("EvalueThreshold", "MaxHspsPerSubject", "Integer(2)")],
    ),
    (
        "qcov_hsp_perc",
        &["-qcov_hsp_perc", "10"],
        &[After("EvalueThreshold", "QueryCovHspPerc", "Real(10.0)")],
    ),
    ("xdrop_ungap", &["-xdrop_ungap", "10"], &[]),
    (
        "xdrop_gap",
        &["-xdrop_gap", "20"],
        &[After("EvalueThreshold", "GapXDropoff", "Real(20.0)")],
    ),
    (
        "xdrop_gap_final",
        &["-xdrop_gap_final", "30"],
        &[After("EvalueThreshold", "GapXDropoffFinal", "Real(30.0)")],
    ),
    (
        "sum_stats",
        &["-sum_stats", "false"],
        &[After("EvalueThreshold", "SumStatistics", "Boolean(false)")],
    ),
    (
        "seg",
        &["-seg", "yes"],
        &[At("SegFiltering", "Boolean(true)")],
    ),
    (
        "soft_masking",
        &["-soft_masking", "true"],
        &[At("MaskAtHash", "Boolean(true)")],
    ),
    (
        "window_size",
        &["-window_size", "30"],
        &[After("SegFiltering", "WindowSize", "Integer(30)")],
    ),
    ("lcase_masking", &["-lcase_masking"], &[]),
    ("query_loc", &["-query_loc", "1-10"], &[]),
    ("parse_deflines", &["-parse_deflines"], &[]),
    ("num_threads", &["-num_threads", "2"], &[]),
    ("mt_mode", &["-mt_mode", "1", "-num_threads", "2"], &[]),
];

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:70-127; the defaults
// of NCBI BLAST+ 2.17.0 `p9.rps.tk.f6`, the options of `p9.rps.q40.opts.f6`.
#[test]
fn every_rpsblast_argument_is_recorded_as_ncbi_records_it() {
    const RPSBLAST: &[Case] = &[
        // RPS-BLAST's reading of the mode drops the `u` (`rps_composition_mode`); no `UnifiedP`.
        (
            "comp_based_stats",
            &["-comp_based_stats", "2u"],
            &[At("CompositionBasedStats", "Integer(2)")],
        ),
    ];
    walk(&App {
        program: "rpsblast",
        base: &["-query", "q.faa", "-db", "d"],
        default: &[
            "EvalueThreshold=Evalue(10.0)",
            "MaskAtHash=Boolean(false)",
            "SegFiltering=Boolean(false)",
            "HitlistSize=Integer(500)",
            "CompositionBasedStats=Integer(1)",
            "SmithWatermanMode=Boolean(false)",
        ],
        cases: with(&[RPS_COMMON, RPSBLAST, FORMATTING]),
        inputs: |command| match command {
            Commands::Rpsblast(args) => {
                let mut inputs = args.option_inputs();
                inputs.hitlist_size = hitlist(
                    &args.outfmt,
                    args.max_target_seqs,
                    args.num_descriptions,
                    args.num_alignments,
                );
                (inputs, args.db.clone())
            }
            _ => unreachable!(),
        },
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:45-118; the defaults
// of NCBI BLAST+ 2.17.0 `p9.rtn.nq4.f6`, the options of `p9.rtn.c0.nq.opts.f6` (`-strand` is never
// recorded: `p9.rtn.nq4.minus.f6`).
#[test]
fn every_rpstblastn_argument_is_recorded_as_ncbi_records_it() {
    const RPSTBLASTN: &[Case] = &[
        (
            "query_gencode",
            &["-query_gencode", "2"],
            &[At("QueryGeneticCode", "Integer(2)")],
        ),
        ("strand", &["-strand", "minus"], &[]),
        (
            "ungapped",
            &["-ungapped", "-comp_based_stats", "0"],
            &[
                At("UngappedMode", "Boolean(true)"),
                At("CompositionBasedStats", "Integer(0)"),
            ],
        ),
        // `SetSmithWatermanMode` in every mode, 0 included (blast_args.cpp:885; review 3 R3-2).
        (
            "use_sw_tback",
            &["-use_sw_tback", "-comp_based_stats", "0"],
            &[
                At("CompositionBasedStats", "Integer(0)"),
                At("SmithWatermanMode", "Boolean(true)"),
            ],
        ),
    ];
    walk(&App {
        program: "rpstblastn",
        base: &["-query", "q.fna", "-db", "d"],
        default: &[
            "EvalueThreshold=Evalue(10.0)",
            "QueryGeneticCode=Integer(1)",
            "MaskAtHash=Boolean(false)",
            "SegFiltering=Boolean(false)",
            "UngappedMode=Boolean(false)",
            "HitlistSize=Integer(500)",
            "CompositionBasedStats=Integer(1)",
            "SmithWatermanMode=Boolean(false)",
        ],
        cases: with(&[RPS_COMMON, RPSTBLASTN, FORMATTING]),
        inputs: |command| match command {
            Commands::Rpstblastn(args) => {
                let mut inputs = args.option_inputs();
                inputs.hitlist_size = hitlist(
                    &args.outfmt,
                    args.max_target_seqs,
                    args.num_descriptions,
                    args.num_alignments,
                );
                (inputs, args.db.clone())
            }
            _ => unreachable!(),
        },
    });
}
