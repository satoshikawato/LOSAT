//! Every argument that the command line of psiblast, rpsblast and rpstblastn takes, and those of
//! psiblast and tblastn with a PSSM query (`-in_pssm`), against what the search strategy's
//! requests record of it (line RG F1b-1, review 3 R3-3; F1b-2): each application's options
//! (`option_inputs` of its own arguments, then `options::algorithm_options` of its `RequestKind`)
//! and the request's `HitlistSize` (`strategy::hitlist_size`) and `Web_StepNumber`.
//! The expectations follow NCBI's argument classes (`psiblast_args.cpp`, `rpsblast_args.cpp`,
//! `rpstblastn_args.cpp`, `tblastn_args.cpp` and their `ExtractAlgorithmOptions` in
//! `blast_args.cpp`) and the rows of NCBI BLAST+ 2.17.0 of probe 9 (`F/f1b/p9_docs.txt`), of
//! review 3 and of probe 7 (`F/f1d/p7ncbi/tbn.pssm.db.ss`): an argument is recorded where its class
//! records it, or not at all. An argument of the command line without an expectation fails the
//! test, so a new argument needs its NCBI rule here.

use super::options::OptionInputs;
use super::strategy::{hitlist_size, HitlistArgs};
use super::{ArchiveEncoding, ArchiveRequest};
use crate::cli::{Cli, Commands};

/// A change of the default request: the parameter in its place, a new one after another, or a
/// parameter gone; or the parser refuses the words with the base command line (`Conflicts`: an
/// argument that `-in_pssm` excludes).
#[derive(Clone, Copy, Debug)]
enum Edit {
    At(&'static str, &'static str),
    After(&'static str, &'static str, &'static str),
    Drop(&'static str),
    Conflicts,
}

use Edit::{After, At, Conflicts, Drop};

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
    let request =
        ArchiveRequest::for_inputs(ArchiveEncoding::Text, inputs, database, "TEST", "-x").unwrap();
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
            Conflicts => unreachable!("a case the parser refuses has no request"),
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
        // `-subject` takes the place of `-db`, `-in_pssm` of `-query`, and an argument of the base
        // given again its own place (the base is pairs of an argument and its value).
        let mut line: Vec<&str> = Vec::new();
        for pair in app.base.chunks(2) {
            let flag = pair[0];
            let replaced = (flag == "-db" && words.contains(&"-subject"))
                || (flag == "-query" && words.contains(&"-in_pssm"))
                || words.contains(&flag);
            if !replaced {
                line.extend_from_slice(pair);
            }
        }
        line.extend_from_slice(words);
        if matches!(edits, [Conflicts]) {
            assert!(
                parse(app.program, &line).is_err(),
                "{} -{id}: {words:?} is refused with {:?}",
                app.program,
                app.base
            );
            continue;
        }
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

/// psiblast's arguments (with a FASTA query or a PSSM query: one handle, one class list).
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
    // A PSSM query's request records the same options (`p9.pp.db.itdef.f6`; the PSSM walk).
    ("in_pssm", &["-in_pssm", "p.asn"], &[]),
];

/// psiblast's request of the base command line (`p9.ps.q.db.itdef.f6`, `p9.pp.db.itdef.f6`).
const PSIBLAST_DEFAULT: &[&str] = &[
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
];

fn psiblast_inputs(command: Commands) -> (OptionInputs, Option<String>) {
    match command {
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
    }
}

/// `cases` with the arguments `excluded` (their ids and words) refused by the parser
/// (`Conflicts`).
fn excluding(cases: &[Case], excluded: &[(&str, &'static [&'static str])]) -> &'static [Case] {
    Box::leak(
        cases
            .iter()
            .map(
                |&(id, words, edits)| match excluded.iter().find(|(name, _)| *name == id) {
                    Some(&(_, words)) => (id, words, &[Conflicts][..]),
                    None => (id, words, edits),
                },
            )
            .collect::<Vec<Case>>()
            .into_boxed_slice(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/psiblast_args.cpp:44-121 (the classes
// and their order, `options::class_order`); the defaults of NCBI BLAST+ 2.17.0 `p9.ps.q.db.itdef.f6`.
#[test]
fn every_psiblast_argument_is_recorded_as_ncbi_records_it() {
    walk(&App {
        program: "psiblast",
        base: &["-query", "q.faa", "-db", "d"],
        default: PSIBLAST_DEFAULT,
        cases: with(&[PSIBLAST, DB_FILTERS, FORMATTING]),
        inputs: psiblast_inputs,
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1241-1247 (`-in_pssm`
// excludes `-query` and `-query_loc`)
// ```c++
//     arg_desc.SetDependency(kArgPSIInputChkPntFile, CArgDescriptions::eExcludes,
//                            kArgQuery);
// ```
// A PSSM query's request records psiblast's options as a FASTA query's does (the handle and the
// classes are the same; NCBI BLAST+ 2.17.0 `p9.pp.db.itdef.f6` = `p9.ps.q.db.itdef.f6`, `pp.db.cbs*`,
// `pp.db.it*`, `pp.db.ssp.f6`): every psiblast argument again, with `-in_pssm` in the base.
#[test]
fn every_psiblast_argument_with_a_pssm_query_is_recorded_as_ncbi_records_it() {
    walk(&App {
        program: "psiblast",
        base: &["-in_pssm", "p.asn", "-db", "d"],
        default: PSIBLAST_DEFAULT,
        cases: excluding(
            with(&[PSIBLAST, DB_FILTERS, FORMATTING]),
            &[
                ("query", &["-query", "q.faa"]),
                ("query_loc", &["-query_loc", "1-10"]),
            ],
        ),
        inputs: psiblast_inputs,
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:44-131 (tblastn's
// classes and their order, `options::class_order`), 138-149 (PSI-TBLASTN's handle, whatever
// `-task`); psiblast_options.cpp:69-75 (`SetPSITblastnDefaults`); blast_args.cpp:866-868 (mode 0)
// and 1050-1055 (`-db_gencode` not for `ePSITblastn`); blast_options.c:1202-1205 (the threshold
// suggestion + 2). The defaults of NCBI BLAST+ 2.17.0 `tbn.pssm.db.ss` (line RG F1d probe 7).
#[test]
fn every_tblastn_argument_with_a_pssm_query_is_recorded_as_ncbi_records_it() {
    const TBLASTN: &[Case] = &[
        ("query", &["-query", "q.faa"], &[Conflicts]),
        ("query_loc", &["-query_loc", "1-10"], &[Conflicts]),
        ("in_pssm", &["-in_pssm", "p.asn"], &[]),
        // NCBI stops before the request with `-subject` ("Remote PSI-BL2SEQ is not supported"); the
        // options would be the same.
        ("subject", &["-subject", "s.fna"], &[]),
        (
            "subject_loc",
            &["-subject", "s.fna", "-subject_loc", "1-10"],
            &[],
        ),
        ("db", &[], &[]),
        (
            "dbsize",
            &["-dbsize", "1000000"],
            &[After("WordThreshold", "DbLength", "BigInteger(1000000)")],
        ),
        ("show_gis", &["-show_gis"], &[]),
        ("parse_deflines", &["-parse_deflines"], &[]),
        (
            "searchsp",
            &["-searchsp", "1000000"],
            &[After(
                "WordThreshold",
                "EffectiveSearchSpace",
                "BigInteger(1000000)",
            )],
        ),
        // PSI-TBLASTN's handle is not the task's: no word size 5 or threshold 20.
        ("task", &["-task", "tblastn-fast"], &[]),
        ("out", &["-out", "o.txt"], &[]),
        (
            "export_search_strategy",
            &["-export_search_strategy", "s.asn"],
            &[],
        ),
        (
            "evalue",
            &["-evalue", "0.001"],
            &[At("EvalueThreshold", "Evalue(0.001)")],
        ),
        (
            "word_size",
            &["-word_size", "4"],
            &[After("WordThreshold", "WordSize", "Integer(4)")],
        ),
        (
            "gap_open",
            &["-gapopen", "9"],
            &[After("WordThreshold", "GapOpeningCost", "Integer(9)")],
        ),
        (
            "gap_extend",
            &["-gapextend", "2"],
            &[After("WordThreshold", "GapExtensionCost", "Integer(2)")],
        ),
        // The constructor's code 1 stays (`ePSITblastn` is not among the class's programs).
        ("db_gencode", &["-db_gencode", "4"], &[]),
        (
            "max_intron_length",
            &["-max_intron_length", "100"],
            &[At("LongestIntronLength", "Integer(100)")],
        ),
        // BLOSUM45: 14/2 (generic class), the name (matrix class), the threshold 14 + 2 and the
        // window 60 in place.
        (
            "matrix",
            &["-matrix", "BLOSUM45"],
            &[
                After("WordThreshold", "GapOpeningCost", "Integer(14)"),
                After("GapOpeningCost", "GapExtensionCost", "Integer(2)"),
                After("SegFilteringHicut", "MatrixName", "Str(\"BLOSUM45\")"),
                At("WordThreshold", "Integer(16)"),
                At("WindowSize", "Integer(60)"),
            ],
        ),
        (
            "threshold",
            &["-threshold", "12"],
            &[At("WordThreshold", "Integer(12)")],
        ),
        ("xdrop_ungap", &["-xdrop_ungap", "8"], &[]),
        (
            "xdrop_gap",
            &["-xdrop_gap", "16"],
            &[After("WordThreshold", "GapXDropoff", "Real(16.0)")],
        ),
        (
            "xdrop_gap_final",
            &["-xdrop_gap_final", "26"],
            &[After("WordThreshold", "GapXDropoffFinal", "Real(26.0)")],
        ),
        // Mode 0 whatever is given.
        ("comp_based_stats", &["-comp_based_stats", "2"], &[]),
        ("outfmt", &["-outfmt", "6"], &[]),
        (
            "seg",
            &["-seg", "yes"],
            &[
                After("MaskAtHash", "SegFiltering", "Boolean(true)"),
                Drop("SegFilteringWindow"),
                Drop("SegFilteringLocut"),
                Drop("SegFilteringHicut"),
            ],
        ),
        ("lcase_masking", &["-lcase_masking"], &[]),
        (
            "soft_masking",
            &["-soft_masking", "true"],
            &[At("MaskAtHash", "Boolean(true)")],
        ),
        (
            "sum_stats",
            &["-sum_stats", "true"],
            &[After("WordThreshold", "SumStatistics", "Boolean(true)")],
        ),
        (
            "window_size",
            &["-window_size", "30"],
            &[At("WindowSize", "Integer(30)")],
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
            "max_hsps",
            &["-max_hsps", "3"],
            &[After("WordThreshold", "MaxHspsPerSubject", "Integer(3)")],
        ),
        (
            "qcov_hsp_perc",
            &["-qcov_hsp_perc", "10"],
            &[After("WordThreshold", "QueryCovHspPerc", "Real(10.0)")],
        ),
        (
            "culling_limit",
            &["-culling_limit", "2"],
            &[After("SegFilteringHicut", "Culling", "Integer(2)")],
        ),
        ("num_threads", &["-num_threads", "2"], &[]),
        ("mt_mode", &["-mt_mode", "1", "-num_threads", "2"], &[]),
        (
            "ungapped",
            &["-ungapped"],
            &[At("UngappedMode", "Boolean(true)")],
        ),
    ];
    walk(&App {
        program: "tblastn",
        base: &["-in_pssm", "p.asn", "-db", "d"],
        default: &[
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
        ],
        cases: with(&[TBLASTN, DB_FILTERS, FORMATTING]),
        inputs: |command| match command {
            Commands::Tblastn(args) => {
                let mut inputs = crate::algorithm::tblastn::archive::option_inputs(&args);
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
