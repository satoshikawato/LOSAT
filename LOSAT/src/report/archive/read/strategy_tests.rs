//! The option builder (`options`) and the queries' writer (`queries`) on strategy files read by
//! the text reader (`read_document`, `request::Request`): every NCBI strategy file of the writer's
//! fixtures read back into the setters of what it records, its queries written again (decision D2 of `F/f2/SURVEY_F2.md`: the writer and the builder share
//! only the field registry, and this test keeps them consistent), and the crafted files of the
//! survey's section 2.2 with what NCBI BLAST+ 2.17.0 did with each.

use std::path::Path;

use super::options::{import_request, EProgram, SeqRange, Setter};
use super::queries::{render_queries, Text};
use super::request::{ParamSet, Queries, Request};
use crate::report::archive::fields;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::{read_document, Value};

/// The tree that the probe build's exception reports name.
const PROBE_BUILD_ROOT: &str = "/opt/conda/conda-bld/blast_1787399212531/work/c++";

fn fixtures(dir: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(dir)
}

/// The request of a strategy file's text (`Blast4-request`), read as the import reads its first
/// object; `None` when the text reader fails.
fn request(text: &[u8]) -> Option<Request> {
    let (value, posts) = read_document(text, named::BLAST4_REQUEST, Names::STRATEGY);
    assert!(posts.is_empty(), "{posts:?}");
    Some(Request::from_strategy(value.ok()?).unwrap_or_else(|stop| panic!("{stop:?}")))
}

/// A value as the writer records it and a setter takes it.
#[derive(Debug, PartialEq)]
enum Recorded {
    Int(i64),
    Real(f64),
    Bool(bool),
    Text(String),
}

/// The value of a stored parameter (`Blast4-value`) of a scalar type, or of an e-value cutoff.
fn recorded_value(value: V<'_>) -> Recorded {
    let (variant, v) = value.choice();
    let v = match variant {
        "cutoff" => match v.choice() {
            ("e-value", e) => e,
            other => panic!("{other:?}"),
        },
        _ => v,
    };
    match v.value {
        Value::Int(n) => Recorded::Int(*n),
        Value::Enum(n) => Recorded::Int(i64::from(*n)),
        Value::Real(r) => Recorded::Real(*r),
        Value::Bool(b) => Recorded::Bool(*b),
        Value::Str(s) => Recorded::Text(String::from_utf8(s.clone()).unwrap()),
        other => panic!("{variant}: not a scalar: {other:?}"),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:140-545
/// The parameter whose value a setter takes, and that value: the field registry's name of each
/// setter of the builder's switch.
fn recorded(setter: &Setter) -> (&'static str, Recorded) {
    use Recorded::{Bool, Int, Real, Text};
    let int = |v: &i32| Int(i64::from(*v));
    match setter {
        Setter::BestHitScoreEdge(v) => (fields::BEST_HIT_SCORE_EDGE.name(), Real(*v)),
        Setter::BestHitOverhang(v) => (fields::BEST_HIT_OVERHANG.name(), Real(*v)),
        Setter::CompositionBasedStats(v) => (fields::COMPOSITION_BASED_STATS.name(), int(v)),
        Setter::CullingLimit(v) => (fields::CULLING_LIMIT.name(), int(v)),
        Setter::DbGeneticCode(v) => (fields::DB_GENETIC_CODE.name(), int(v)),
        Setter::DbLength(v) => (fields::DB_LENGTH.name(), Int(*v)),
        Setter::DustFiltering(v) => (fields::DUST_FILTERING.name(), Bool(*v)),
        Setter::DustFilteringLevel(v) => (fields::DUST_FILTERING_LEVEL.name(), int(v)),
        Setter::DustFilteringWindow(v) => (fields::DUST_FILTERING_WINDOW.name(), int(v)),
        Setter::DustFilteringLinker(v) => (fields::DUST_FILTERING_LINKER.name(), int(v)),
        Setter::DomainInclusionThreshold(v) => {
            (fields::DOMAIN_INCLUSION_THRESHOLD.name(), Real(*v))
        }
        Setter::EffectiveSearchSpace(v) => (fields::EFFECTIVE_SEARCH_SPACE.name(), Int(*v)),
        Setter::EvalueThreshold(v) => (fields::EVALUE_THRESHOLD.name(), Real(*v)),
        Setter::FilterString(v) => (fields::FILTER_STRING.name(), Text(v.clone())),
        Setter::GapExtensionCost(v) => (fields::GAP_EXTENSION_COST.name(), int(v)),
        Setter::GapOpeningCost(v) => (fields::GAP_OPENING_COST.name(), int(v)),
        Setter::GapTracebackAlgorithm(v) => (fields::GAP_TRACEBACK_ALGORITHM.name(), int(v)),
        Setter::GapTrigger(v) => (fields::GAP_TRIGGER.name(), Real(*v)),
        Setter::GapXDropoff(v) => (fields::GAP_X_DROPOFF.name(), Real(*v)),
        Setter::GapXDropoffFinal(v) => (fields::GAP_X_DROPOFF_FINAL.name(), Real(*v)),
        Setter::GapExtnAlgorithm(v) => (fields::GAP_EXTN_ALGORITHM.name(), int(v)),
        Setter::HitlistSize(v) => (fields::HITLIST_SIZE.name(), int(v)),
        Setter::InclusionThreshold(v) => (fields::INCLUSION_THRESHOLD.name(), Real(*v)),
        Setter::IgnoreMsaMaster(v) => (fields::IGNORE_MSA_MASTER.name(), Bool(*v)),
        Setter::LongestIntronLength(v) => (fields::LONGEST_INTRON_LENGTH.name(), int(v)),
        Setter::MbTemplateLength(v) => (fields::MB_TEMPLATE_LENGTH.name(), Int(i64::from(*v))),
        Setter::MbTemplateType(v) => (fields::MB_TEMPLATE_TYPE.name(), Int(i64::from(*v))),
        Setter::MatchReward(v) => (fields::MATCH_REWARD.name(), int(v)),
        Setter::MatrixName(v) => (fields::MATRIX_NAME.name(), Text(v.clone())),
        Setter::MismatchPenalty(v) => (fields::MISMATCH_PENALTY.name(), int(v)),
        Setter::MaskAtHash(v) => (fields::MASK_AT_HASH.name(), Bool(*v)),
        Setter::MaxHspsPerSubject(v) => (fields::MAX_HSPS_PER_SUBJECT.name(), int(v)),
        Setter::OutOfFrameMode(v) => (fields::OUT_OF_FRAME_MODE.name(), Bool(*v)),
        Setter::PhiPattern { pattern, .. } => (fields::PHI_PATTERN.name(), Text(pattern.clone())),
        Setter::PercentIdentity(v) => (fields::PERCENT_IDENTITY.name(), Real(*v)),
        Setter::PseudoCount(v) => (fields::PSEUDO_COUNT.name(), int(v)),
        Setter::QueryGeneticCode(v) => (fields::QUERY_GENETIC_CODE.name(), int(v)),
        Setter::QueryCovHspPerc(v) => (fields::QUERY_COV_HSP_PERC.name(), Real(*v)),
        Setter::RepeatFiltering(v) => (fields::REPEAT_FILTERING.name(), Bool(*v)),
        Setter::RepeatFilteringDb(v) => (fields::REPEAT_FILTERING_DB.name(), Text(v.clone())),
        Setter::StrandOption(v) => (fields::STRAND_OPTION.name(), int(v)),
        Setter::SegFiltering(v) => (fields::SEG_FILTERING.name(), Bool(*v)),
        Setter::SegFilteringWindow(v) => (fields::SEG_FILTERING_WINDOW.name(), int(v)),
        Setter::SegFilteringLocut(v) => (fields::SEG_FILTERING_LOCUT.name(), Real(*v)),
        Setter::SegFilteringHicut(v) => (fields::SEG_FILTERING_HICUT.name(), Real(*v)),
        Setter::SumStatisticsMode(v) => (fields::SUM_STATISTICS_MODE.name(), Bool(*v)),
        Setter::SmithWatermanMode(v) => (fields::SMITH_WATERMAN_MODE.name(), Bool(*v)),
        // `SetGappedMode` with the negation of `UngappedMode`'s value.
        Setter::GappedMode(v) => (fields::GAPPED_MODE.name(), Bool(!*v)),
        Setter::UnifiedP(v) => (fields::UNIFIED_P.name(), int(v)),
        Setter::SubjectBestHit => (fields::SUBJECT_BEST_HIT.name(), Bool(true)),
        Setter::WindowMaskerTaxId(v) => (fields::WINDOW_MASKER_TAX_ID.name(), int(v)),
        Setter::WindowSize(v) => (fields::WINDOW_SIZE.name(), int(v)),
        Setter::WordSize(v) => (fields::WORD_SIZE.name(), int(v)),
        Setter::WordThreshold(v) => (fields::WORD_THRESHOLD.name(), int(v)),
        Setter::UseIndex { index_name } => (fields::MB_INDEX_NAME.name(), Text(index_name.clone())),
    }
}

/// The strategy files of the writer's fixtures (`tests/fixtures/strategy`, NCBI's own exports,
/// equal to LOSAT's at every thread count) that hold a request: what NCBI's builder makes of
/// each is what the file records. Each parameter of the algorithm and program options is one
/// setter with its value, in the file's order; the ones the builder keeps aside are checked
/// where it keeps them; a file that holds `WindowMaskerDatabase`, which NCBI writes and its
/// builder has no arm for, fails as NCBI's import of it does (`I.blastn.wmdb.alone`: rc 3
/// `BLAST engine error: Internal: Error processing option [WindowMaskerDatabase] type [12].`).
#[test]
fn every_ncbi_strategy_file_builds_the_setters_it_records() {
    let dir = fixtures("strategy");
    let mut paths: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|e| e == "ss"))
        .collect();
    paths.sort();
    let (mut built, mut refused, mut rendered) = (0, 0, 0);
    for path in paths {
        let case = path.file_name().unwrap().to_string_lossy().into_owned();
        let text = std::fs::read(&path).unwrap();
        // What a failed run of NCBI's left (nothing, or a file the run did not write).
        if !text.starts_with(b"Blast4-request ::= {") {
            continue;
        }
        let request = request(&text).unwrap_or_else(|| panic!("{case}: not read"));
        // The queries NCBI's FASTA reader made (`CFastaReader` titles and user objects, local
        // ids, every coding it stores) are each written again.
        if let Queries::BioseqSet(set) = &request.queries {
            render_queries(V::new(named::BIOSEQ_SET, set), Text::held())
                .unwrap_or_else(|error| panic!("{case}: {error:?}"));
            rendered += 1;
        }
        let imported = import_request(&request);
        if request.params().iter().any(|param| {
            param.set == ParamSet::Algorithm
                && param.name == fields::WINDOW_MASKER_DATABASE.name().as_bytes()
        }) {
            let error = imported.expect_err(&case);
            assert_eq!(
                error.get_msg(),
                b"Internal: Error processing option [WindowMaskerDatabase] type [12].",
                "{case}"
            );
            refused += 1;
            continue;
        }
        let imported = imported.unwrap_or_else(|error| panic!("{case}: {error:?}"));
        let options = &imported.options;
        let ints = |value: V<'_>| -> Vec<i32> {
            value
                .choice()
                .1
                .ints()
                .map(|n| i32::try_from(n).unwrap())
                .collect()
        };
        let mut expected = Vec::new();
        for param in request.params() {
            let name = std::str::from_utf8(&param.name).unwrap();
            let value = param.value();
            let int = || match recorded_value(value) {
                Recorded::Int(n) => n,
                other => panic!("{case} {name}: {other:?}"),
            };
            match (param.set, name) {
                (ParamSet::Format, "Web_StepNumber") => {
                    assert_eq!(i64::from(imported.psi_num_of_iterations), int(), "{case}");
                    assert_eq!(options.program, EProgram::PsiBlast, "{case}");
                }
                (ParamSet::Format, _) => panic!("{case}: format option {name}"),
                (_, "LCaseMask") => assert!(options.query_masks.is_some(), "{case}"),
                (_, "TaxidList") => {
                    assert_eq!(options.taxid_list, Some(ints(value)), "{case}")
                }
                (_, "NegativeTaxidList") => {
                    assert_eq!(options.negative_taxid_list, Some(ints(value)), "{case}")
                }
                (_, "DbFilteringAlgorithmId") => {
                    assert_eq!(i64::from(imported.filtering_id), int(), "{case}")
                }
                (_, "DbFilteringAlgorithmKey") => assert_eq!(
                    Recorded::Text(imported.filtering_key.clone()),
                    recorded_value(value),
                    "{case}"
                ),
                (_, "SubjectMaskingType") => {
                    assert_eq!(i64::from(imported.subject_masking_type), int(), "{case}")
                }
                (_, "RequiredStart") => {
                    assert_eq!(i64::from(options.query_range.from), int(), "{case}")
                }
                (_, "RequiredEnd") => {
                    assert_eq!(i64::from(options.query_range.to_open), int(), "{case}")
                }
                _ => expected.push((name, recorded_value(value))),
            }
        }
        let got: Vec<_> = options.setters.iter().map(recorded).collect();
        assert_eq!(got, expected, "{case}");
        built += 1;
    }
    assert!(
        built >= 400 && refused >= 5 && rendered >= 370,
        "{built} built, {refused} refused, {rendered} query sets rendered"
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:62-75
/// An element stated as `NULL` in a list of the search (an algorithm option, a query's id, a
/// descriptor) is read as NCBI's reader reads it, the empty reference of its list, and the request
/// refuses it before the import reads any list (`request::Stop::NullElement`); a body other than a
/// search keeps `CImportStrategy`'s exception, which comes first.
#[test]
fn a_null_element_of_a_strategy_list_is_refused() {
    let text = std::fs::read(fixtures("strategy").join("blastp.db.batches.ss")).unwrap();
    let text = String::from_utf8(text).unwrap();
    let read = |text: &str| {
        let (value, posts) = read_document(text.as_bytes(), named::BLAST4_REQUEST, Names::STRATEGY);
        assert!(posts.is_empty(), "{posts:?}");
        value.unwrap()
    };
    assert!(Request::from_strategy(read(&text)).is_ok());
    for (at, null) in [
        ("    algorithm-options {\n", "      NULL,\n"),
        ("          id {\n", "            NULL,\n"),
        ("          descr {\n", "            NULL,\n"),
    ] {
        assert!(text.contains(at), "{at}");
        let crafted = text.replacen(at, &format!("{at}{null}"), 1);
        let value = read(&crafted);
        assert!(
            V::new(named::BLAST4_REQUEST, &value).holds_null_element(),
            "{at}"
        );
        let stop = Request::from_strategy(value).unwrap_err();
        assert!(
            matches!(stop, super::request::Stop::NullElement),
            "{at}: {stop:?}"
        );
    }
    let other = r#"Blast4-request ::= {
  body get-sequences {
    database { name "nr", type protein },
    seq-ids { NULL, local str "x" }
  }
}
"#;
    let value = read(other);
    assert!(V::new(named::BLAST4_REQUEST, &value).holds_null_element());
    let stop = Request::from_strategy(value).unwrap_err();
    let super::request::Stop::Exception(exception) = stop else {
        panic!("{stop:?}");
    };
    assert_eq!(exception.get_msg(), b"No body in CBlast4_request");
}

/// The crafted files of survey section 2.2 (`tests/fixtures/strategy_crafted`, the files NCBI
/// read): those NCBI's text reader refuses, those its builder refuses (NCBI's stderr less its
/// warnings), and those it builds, with the value each crafted change sets.
#[test]
fn the_crafted_files_of_the_survey_end_as_ncbis() {
    let dir = fixtures("strategy_crafted");
    let manifest = std::fs::read_to_string(dir.join("cases.tsv")).unwrap();
    let mut checked = 0;
    for line in manifest.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.split('\t');
        let (case, kind) = (fields.next().unwrap(), fields.next().unwrap());
        let text = std::fs::read(dir.join(format!("{case}.ss"))).unwrap();
        let (program, change) = case.split_once('.').unwrap();
        if kind == "read" {
            let (value, _) = read_document(&text, named::BLAST4_REQUEST, Names::STRATEGY);
            assert!(value.is_err(), "{case}");
            checked += 1;
            continue;
        }
        let request = request(&text).unwrap_or_else(|| panic!("{case}: not read"));
        let imported = import_request(&request);
        if kind == "builder" {
            let error = imported.expect_err(case);
            let ncbi = std::fs::read_to_string(dir.join(format!("{case}.err"))).unwrap();
            let report = if ncbi.starts_with("BLAST engine error: ") {
                format!(
                    "BLAST engine error: {}\n",
                    String::from_utf8_lossy(error.get_msg())
                )
            } else {
                format!(
                    "Error: {}\n",
                    String::from_utf8_lossy(&error.what_in(PROBE_BUILD_ROOT))
                )
            };
            assert_eq!(
                crate::wording::losat_text(&report),
                crate::wording::losat_text(&ncbi),
                "{case}"
            );
            checked += 1;
            continue;
        }
        assert_eq!(kind, "built", "{case}");
        let imported = imported.unwrap_or_else(|error| panic!("{case}: {error:?}"));
        let options = &imported.options;
        let evalues: Vec<f64> = options
            .setters
            .iter()
            .filter_map(|setter| match setter {
                Setter::EvalueThreshold(e) => Some(*e),
                _ => None,
            })
            .collect();
        match change {
            // The crafted e-value comes first in the algorithm options and the export's own
            // after it (the last call wins: the export's 10), or alone in the program options,
            // which come after (1e-30 wins).
            "evalue_cutoff" | "dupparam" => assert_eq!(evalues, [1e-30, 10.0], "{case}"),
            "algo_in_prog" => assert_eq!(evalues, [10.0, 1e-30], "{case}"),
            "web_step3" => {
                assert_eq!(options.program, EProgram::PsiBlast, "{case}");
                assert_eq!(imported.psi_num_of_iterations, 3, "{case}");
            }
            "phi" => {
                assert_eq!((program, options.program), ("blastp", EProgram::PhiBlastp))
            }
            "required_start_only" => assert_eq!(
                options.query_range,
                SeqRange {
                    from: 10,
                    to_open: u32::MAX
                }
            ),
            "qrange_out" => assert_eq!(
                options.query_range,
                SeqRange {
                    from: 5000,
                    to_open: 6000
                }
            ),
            "hsprange_culling" => assert!(options.setters.contains(&Setter::CullingLimit(1))),
            "taxid_empty" => assert_eq!(options.taxid_list.as_deref(), Some(&[][..])),
            "neg_taxid" => assert_eq!(options.negative_taxid_list.as_deref(), Some(&[9606][..])),
            "gilist" => assert_eq!(options.gi_list.as_deref(), Some(&[1, 2][..])),
            "entrez" => assert_eq!(options.entrez_query.as_deref(), Some("human[orgn]")),
            "firstdbseq" => assert_eq!(
                (options.first_db_seq, options.final_db_seq),
                (Some(2), Some(3))
            ),
            // The rest change what the builder does not read (the queries, the subject, the
            // ident, the objects after the first): built is all.
            _ => {}
        }
        checked += 1;
    }
    assert!(checked >= 60, "{checked} cases");
}

/// F2-Q1 (Owner 2026-11-20 「Q1=C」, PD-LOSAT-NCBI-DEFECTS exception 32): the import keeps the
/// command line's `-num_threads` for the strategy's subjects, which the run searches once with
/// them (`tests/import_search_strategy.rs` checks that the bytes are the `-num_threads 1` ones),
/// also where blastn runs a TBLASTX strategy.
#[test]
fn a_strategys_subjects_keep_the_threads_of_the_command_line() {
    use crate::cli::Commands;
    let dir = fixtures("strategy/import/ss");
    for (program, file) in [
        ("blastn", "blastn.sub.ss"),
        ("blastp", "blastp.sub.ss"),
        ("blastx", "blastx.sub.ss"),
        ("tblastn", "tblastn.sub.ss"),
        ("tblastx", "tblastx.sub.ss"),
        ("blastn", "tblastx.sub.ss"),
    ] {
        let path = dir.join(file);
        let argv = [
            "losat",
            program,
            "-import_search_strategy",
            path.to_str().unwrap(),
            "-num_threads",
            "4",
        ];
        let cli = crate::cli::try_parse_from(argv)
            .unwrap_or_else(|error| panic!("{program} {file}: {error}"));
        let mut posted = Vec::new();
        let command = super::import::apply(cli.command, None, &mut posted)
            .unwrap_or_else(|error| panic!("{program} {file}: {error:#}"));
        let threads = match &command {
            Commands::Blastn(args) => args.num_threads,
            Commands::Blastp(args) => args.num_threads,
            Commands::Blastx(args) => usize::try_from(args.num_threads).unwrap(),
            Commands::Tblastn(args) => args.num_threads,
            Commands::Tblastx(args) => args.num_threads,
            _ => panic!("{program} {file}: not a search"),
        };
        assert_eq!(threads, 4, "{program} {file}");
    }
}

/// One import of the strategy `file` by the running program `running` with `options`: the command.
fn imported(running: &str, file: &str, options: &[&str]) -> crate::cli::Commands {
    let path = fixtures("strategy/import/ss").join(file);
    let mut argv = vec![
        "losat",
        running,
        "-import_search_strategy",
        path.to_str().unwrap(),
    ];
    argv.extend_from_slice(options);
    let cli = crate::cli::try_parse_from(argv)
        .unwrap_or_else(|error| panic!("{running} {file}: {error}"));
    let mut posted = Vec::new();
    super::import::apply(cli.command, None, &mut posted)
        .unwrap_or_else(|error| panic!("{running} {file}: {error:#}"))
}

/// R2-B1 (review 2 of F2): a strategy of another program runs inside the running application,
/// whose identity the one typed value `ImportRun::app` carries into the run: the engine of the
/// strategy's program reads it (`app()`), not its own program's.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:543-549
#[test]
fn a_cross_program_import_carries_the_running_application() {
    use crate::blastinput::app::Program;
    use crate::cli::Commands;
    for (running, file, expected) in [
        ("blastn", "tblastx.db.ss", Program::Blastn),
        ("blastx", "tblastx.db.ss", Program::Blastx),
        ("tblastx", "blastn.db.ss", Program::Tblastx),
        ("tblastx", "blastx.db.ss", Program::Tblastx),
        ("blastn", "blastn.db.ss", Program::Blastn),
        ("tblastx", "tblastx.db.ss", Program::Tblastx),
    ] {
        let (app, own, import) = match imported(running, file, &[]) {
            Commands::Blastn(args) => (args.app(), Program::Blastn, args.import),
            Commands::Blastx(args) => (args.app(), Program::Blastx, args.import),
            Commands::Tblastx(args) => (args.app(), Program::Tblastx, args.import),
            _ => panic!("{running} {file}: not a search"),
        };
        assert_eq!(app, expected, "{running} {file}: the running application");
        assert_eq!(
            import.expect("the import's run state").app,
            expected,
            "{running} {file}"
        );
        assert_eq!(
            (own == expected),
            file.starts_with(running),
            "{running} {file}: whether the engine is the running program's"
        );
    }
}

/// R2-B1: the formatting handler's warnings of the run carry the running application's name, the
/// BLASTN engine's (tblastx running a blastn strategy) as the TBLASTX engine's
/// (`run_impl_tests`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2947-2989
#[test]
fn the_warnings_of_a_cross_program_run_name_the_running_application() {
    use crate::cli::Commands;
    let Commands::Blastn(args) = imported(
        "tblastx",
        "blastn.db.ss",
        &[
            "-outfmt",
            "6",
            "-line_length",
            "30",
            "-max_target_seqs",
            "3",
        ],
    ) else {
        panic!("tblastx running a blastn strategy");
    };
    let (_, warnings) = crate::algorithm::blastn::scoring::report_counts_of(&args, &[6]).unwrap();
    assert_eq!(
        String::from_utf8(warnings).unwrap(),
        "Warning: [tblastx] The parameter -line_length is not applicable for output formats > 4 .\n\
         Warning: [tblastx] Examining 5 or more matches is recommended\n"
    );
}

/// R2-B1: BLAST XML and JSON cite the megablast paper by the options' task, the pairwise report
/// by the formatter of the running application (BLASTN's makes it with `-task megablast`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:653-668
#[test]
fn the_xml_reference_follows_the_options_and_the_pairwise_one_the_application() {
    use crate::cli::Commands;
    for (running, pairwise) in [("blastn", true), ("tblastx", false)] {
        let Commands::Blastn(args) = imported(running, "blastn.db.ss", &[]) else {
            panic!("{running} running a blastn strategy");
        };
        assert_eq!(args.task, "megablast", "{running}");
        assert!(
            args.structured_megablast_reference(),
            "{running}: XML, JSON"
        );
        assert_eq!(
            args.pairwise_megablast_reference(),
            pairwise,
            "{running}: outfmt 0"
        );
    }
}
