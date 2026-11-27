use super::*;
use crate::blastinput::app::Program;
use crate::report::archive::bioseq::BioseqInput;
use crate::report::archive::masks::{Frame, MaskEntry};
use crate::report::archive::ArchiveHandle;

/// The strategy file's bytes, which the test reads while the export writes them.
#[derive(Clone, Default)]
struct Sink(Arc<Mutex<Vec<u8>>>);

impl StrategySink for Sink {
    fn write_request(
        &mut self,
        bytes: &[u8],
        _pieces: &[BufferPiece],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(())
    }

    fn write_partial(
        &mut self,
        bytes: &[u8],
        posts: &[(u64, Vec<u8>)],
        diagnostics: &mut dyn Write,
        archive: Option<&ArchiveHandle>,
    ) -> std::io::Result<()> {
        let _ = emit_strategy_posts(posts, None, diagnostics, archive)?;
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(())
    }
}

impl Sink {
    fn text(&self) -> String {
        String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
    }
}

/// A searched batch whose context 0 has the search space `space` (a query of 142 letters).
fn searched(space: i64) -> Option<crate::report::archive::SearchedContext0> {
    Some(crate::report::archive::SearchedContext0 {
        space,
        query_length: 142,
    })
}

fn request<'a>(queries: &'a [BioseqInput<'a>], lcase_masks: &'a [MaskEntry]) -> RequestInputs<'a> {
    RequestInputs {
        queries: crate::report::archive::writer::RequestQueries::Bioseqs(queries),
        subjects: &[],
        required: None,
        lcase_masks,
    }
}

fn protein<'a>(id: &'a str, letters: &'a [u8]) -> BioseqInput<'a> {
    BioseqInput {
        id,
        ids: None,
        seq_loc: None,
        defline: std::borrow::Cow::Borrowed(b">q a query"),
        title: b"q a query",
        letters,
        protein: true,
    }
}

/// A BLASTP export into a sink, with its own options handle (`searchsp`: `-searchsp`).
fn export(searchsp: Option<i64>) -> (StrategyExport, Sink) {
    let sink = Sink::default();
    let mut inputs = OptionInputs::new(Program::Blastp, "blastp");
    inputs.searchsp = searchsp;
    let export = StrategyExport::new(
        Box::new(sink.clone()),
        &inputs,
        Some(b"in/db/bp"),
        OptionsHandle::default(),
    );
    (export, sink)
}

// NCBI BLAST+ 2.17.0, `blastp -query q.faa -db in/db/bp -export_search_strategy f` with a
// query `>q a query` / `MWCHWCMWCHWC` (the layout of the RG-F1 probe's blastp.d.plain.ss).
const FIRST: &str = "Blast4-request ::= {\n  body queue-search {\n    program \"blastp\",\n    service \"plain\",\n    queries bioseq-set {\n      seq-set {\n        seq {\n          id {\n            local str \"Query_1\"\n          },\n          descr {\n            user {\n              type str \"CFastaReader\",\n              data {\n                {\n                  label str \"DefLine\",\n                  data str \">q a query\"\n                }\n              }\n            },\n            title \"q a query\"\n          },\n          inst {\n            repr raw,\n            mol aa,\n            length 12,\n            seq-data ncbieaa \"MWCHWCMWCHWC\"\n          }\n        }\n      }\n    },\n    subject database \"in/db/bp\",\n    algorithm-options {\n      {\n        name \"EvalueThreshold\",\n        value cutoff e-value { 1, 10, 1 }\n      },\n      {\n        name \"MaskAtHash\",\n        value boolean FALSE\n      },\n      {\n        name \"SegFiltering\",\n        value boolean FALSE\n      },\n      {\n        name \"WordThreshold\",\n        value integer 11\n      },\n      {\n        name \"WindowSize\",\n        value integer 40\n      },\n      {\n        name \"HitlistSize\",\n        value integer 500\n      },\n      {\n        name \"UngappedMode\",\n        value boolean FALSE\n      },\n      {\n        name \"CompositionBasedStats\",\n        value integer 2\n      },\n      {\n        name \"SmithWatermanMode\",\n        value boolean FALSE\n      }\n    }\n  }\n}\n";

fn with_space(space: i64) -> String {
    FIRST.replace(
        "        value boolean FALSE\n      }\n    }\n  }\n}\n",
        &format!(
            "        value boolean FALSE\n      }},\n      {{\n        name \"EffectiveSearchSpace\",\n        value big-integer {space}\n      }}\n    }}\n  }}\n}}\n"
        ),
    )
}

/// One request per batch, without `ident`; the `EffectiveSearchSpace` option is absent in
/// the first request and is the batch before's in the next ones, also after a batch that
/// was not searched.
#[test]
fn each_batch_is_a_request_with_the_search_space_of_the_batch_before() {
    let (export, sink) = export(None);
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    assert!(export.write_request(0..1, &request(&queries, &[])).unwrap());
    export.write_posts(&mut std::io::sink(), None).unwrap();
    assert_eq!(sink.text(), FIRST);
    export.note_searched(0..1, searched(817_380)).unwrap();
    for (at, searched) in [(1, None), (2, searched(5)), (3, None)] {
        assert!(export
            .write_request(at..at + 1, &request(&queries, &[]))
            .unwrap());
        export.note_searched(at..at + 1, searched).unwrap();
    }
    end_run(Some(&export), Ok(())).unwrap();
    assert_eq!(
        sink.text(),
        [
            FIRST.to_string(),
            with_space(817_380),
            with_space(817_380),
            with_space(5)
        ]
        .concat()
    );
}

/// `-searchsp` is in the first request as typed and keeps its place in the next ones, with
/// the batch before's value; a database list that keeps no sequence makes that value 0.
#[test]
fn searchsp_keeps_its_place_and_an_empty_database_gives_zero() {
    let (export, sink) = export(Some(1_000_000));
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    let mask = MaskEntry {
        id: "Query_1".to_string(),
        best: None,
        frame: Frame::NotSet,
        intervals: vec![(2, 5)],
    };
    export
        .write_request(0..1, &request(&queries, std::slice::from_ref(&mask)))
        .unwrap();
    export.note_searched(0..1, searched(1_000_000)).unwrap();
    export
        .options()
        .set_database_filter(
            Vec::new(),
            crate::report::archive::RequestLength::Lengths(None),
        )
        .unwrap();
    export.write_request(1..2, &request(&queries, &[])).unwrap();
    export.note_searched(1..2, searched(77)).unwrap();
    export.write_request(2..3, &request(&queries, &[])).unwrap();
    export.note_searched(2..3, None).unwrap();
    end_run(Some(&export), Ok(())).unwrap();
    let all = sink.text();
    let docs: Vec<&str> = all.split("Blast4-request ::= ").skip(1).collect();
    assert_eq!(docs.len(), 3);
    let space_at = |doc: &str| {
        let at = doc.find("name \"EffectiveSearchSpace\"").unwrap();
        let value = &doc[at..];
        let value = &value[value.find("big-integer ").unwrap() + 12..];
        (
            at,
            value[..value.find('\n').unwrap()].parse::<i64>().unwrap(),
        )
    };
    let places: Vec<(usize, i64)> = docs.iter().map(|doc| space_at(doc)).collect();
    assert_eq!(places[0].1, 1_000_000);
    assert_eq!(places[1], places[0]);
    assert_eq!(places[2], (places[0].0, 0));
    // The first request has the `LCaseMask` parameter, the others no `program-options`.
    assert!(docs[0].contains("program-options {\n      {\n        name \"LCaseMask\","));
    assert!(!docs[1].contains("program-options"));
    assert!(!docs[0].contains("ident"));
}

/// PD-LOSAT-NCBI-DEFECTS approved exception 17, Int8 -2^63 written as digits (review C of
/// line RG D1a-4, SF1): a big integer of -2^63 in a request, which NCBI writes as `big-integer -0` (probes q1-q12 of
/// review C), is written as its digits: BLASTP's `-dbsize` in every request, RPS-BLAST's
/// `DbLength` (`x_AdjustDbSize`, `rps_search::volumes::search_sets`) in the requests after the
/// first.
#[test]
fn a_big_integer_of_minus_two_to_the_63_is_written_as_its_digits() {
    let digits = "name \"DbLength\",\n        value big-integer -9223372036854775808\n";
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    let sink = Sink::default();
    let mut inputs = OptionInputs::new(Program::Blastp, "blastp");
    inputs.dbsize = Some(i64::MIN);
    let dbsize = StrategyExport::new(
        Box::new(sink.clone()),
        &inputs,
        Some(b"in/db/bp"),
        OptionsHandle::default(),
    );
    assert!(dbsize.write_request(0..1, &request(&queries, &[])).unwrap());
    end_run(Some(&dbsize), Ok(())).unwrap();
    let text = sink.text();
    assert!(text.contains(digits), "{text}");
    let (rps, sink) = export(None);
    rps.options()
        .set_search_sets(
            Vec::new(),
            crate::algorithm::rps_search::volumes::search_sets(Some(i64::MIN), None),
            Vec::new(),
        )
        .unwrap();
    for at in 0..2 {
        assert!(rps
            .write_request(at..at + 1, &request(&queries, &[]))
            .unwrap());
        rps.note_searched(at..at + 1, searched(5)).unwrap();
    }
    end_run(Some(&rps), Ok(())).unwrap();
    let all = sink.text();
    let docs: Vec<&str> = all.split("Blast4-request ::= ").skip(1).collect();
    assert_eq!(docs.len(), 2);
    assert!(!docs[0].contains("DbLength"), "{all}");
    assert!(docs[1].contains(digits), "{all}");
}

/// A batch searched again (a later part of a batch, a call that searches the batches before
/// again) writes no second request and does not move the search space again.
#[test]
fn a_request_written_once_is_not_written_again() {
    let (export, sink) = export(None);
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    assert!(export.takes(0..2));
    assert!(export.write_request(0..2, &request(&queries, &[])).unwrap());
    export.note_searched(0..2, searched(5)).unwrap();
    // The same batch, and a part of it: not taken, nothing written.
    assert!(!export.takes(0..2));
    assert!(!export.takes(1..2));
    assert!(!export.write_request(0..2, &request(&queries, &[])).unwrap());
    export.note_searched(0..2, searched(9)).unwrap();
    export.write_posts(&mut std::io::sink(), None).unwrap();
    assert_eq!(sink.text(), FIRST);
    assert!(export.write_request(2..3, &request(&queries, &[])).unwrap());
    end_run(Some(&export), Ok(())).unwrap();
    assert_eq!(sink.text(), [FIRST.to_string(), with_space(5)].concat());
}

/// NCBI's stop when it formats a batch: the requests of the batches after it are not written.
#[test]
fn no_request_follows_the_stop() {
    let (export, sink) = export(None);
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    export.stop_after(1);
    assert!(export.write_request(0..1, &request(&queries, &[])).unwrap());
    export.note_searched(0..1, searched(5)).unwrap();
    assert!(!export.takes(1..2));
    assert!(!export.write_request(1..2, &request(&queries, &[])).unwrap());
    export.note_searched(1..2, searched(6)).unwrap();
    end_run(Some(&export), Ok(())).unwrap();
    assert_eq!(sink.text(), FIRST);
}

/// A request that does not start where the one before ended is LOSAT's internal error (a batch
/// whose request was never written), in every build.
#[test]
fn a_gap_between_requests_is_an_error() {
    let (export, _) = export(None);
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    export.write_request(0..1, &request(&queries, &[])).unwrap();
    export.note_searched(0..1, None).unwrap();
    let gap = export
        .write_request(2..3, &request(&queries, &[]))
        .unwrap_err();
    assert!(gap.to_string().contains("LOSAT internal error"), "{gap}");
    // A batch noted out of order too.
    assert!(export.note_searched(3..4, None).is_err());
}

/// The archive and the search strategy of one run share the options handle: the database's
/// list reaches both, each with its own search space carry.
#[test]
fn one_options_handle_serves_the_archive_and_the_strategy() {
    let handle = ArchiveHandle::new(
        crate::report::archive::ArchiveEncoding::Text,
        ArchiveRequest {
            program: "blastp",
            application: "blastp",
            service: "plain",
            algorithm: Vec::new(),
            database: None,
            step_number: None,
        },
        OptionsHandle::default(),
    );
    let sink = Sink::default();
    let export = StrategyExport::new(
        Box::new(sink.clone()),
        &OptionInputs::new(Program::Blastp, "blastp"),
        Some(b"in/db/bp"),
        handle.options().clone(),
    );
    assert!(export.options().same(handle.options()));
    let list = crate::report::archive::options::Param::new(
        crate::report::archive::fields::TAXID_LIST,
        crate::report::archive::options::ParamValue::IntegerList(vec![9606]),
    );
    handle
        .options()
        .set_database_filter(
            vec![list],
            crate::report::archive::RequestLength::Lengths(None),
        )
        .unwrap();
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    export.write_request(0..1, &request(&queries, &[])).unwrap();
    export.note_searched(0..1, searched(5)).unwrap();
    export.write_request(1..2, &request(&queries, &[])).unwrap();
    end_run(Some(&export), Ok(())).unwrap();
    let all = sink.text();
    assert_eq!(all.matches("name \"TaxidList\"").count(), 2, "{all}");
    assert!(all.contains("value big-integer 0\n"), "{all}");
    // The archive's carry is its own: its first archive's option is its batch's.
    let space =
        |searched| crate::report::archive::tests::search_space(&handle.request_sets(searched));
    assert_eq!(space(None), None);
    assert_eq!(space(searched(5)), Some(0));
}

/// A post of a request that no batch start wrote is LOSAT's internal error at the end of a run
/// that ends as NCBI's (`Ok` or NCBI's error); a refusal of LOSAT's keeps its text. The end lets
/// the posts go (each run ends once).
#[test]
fn a_post_left_at_the_end_of_a_run_is_an_internal_error() {
    let queries = [BioseqInput {
        title: b"q caf\xe9",
        defline: std::borrow::Cow::Borrowed(b">q caf\xe9"),
        ..protein("Query_1", b"MWCHWCMWCHWC")
    }];
    let written = || {
        let (export, _) = export(None);
        export.write_request(0..1, &request(&queries, &[])).unwrap();
        export
    };
    let native = || -> Result<()> {
        Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: x\n".into(),
        }
        .into())
    };
    for result in [Ok(()), native()] {
        let export = written();
        let error = end_run(Some(&export), result).unwrap_err();
        assert!(
            error.to_string().starts_with("LOSAT internal error"),
            "{error}"
        );
        assert!(export.take_all_posts().is_empty());
    }
    let export = written();
    assert!(!export.take_all_posts().is_empty());
    end_run(Some(&export), Ok(())).unwrap();
    let export = written();
    let refusal = end_run(Some(&export), Err(anyhow::anyhow!("not supported"))).unwrap_err();
    assert_eq!(refusal.to_string(), "not supported");
    assert!(export.take_all_posts().is_empty());
    let export = written();
    export.take_all_posts();
    let error = end_run(Some(&export), native()).unwrap_err();
    assert!(error.is::<crate::cli::NativeError>());
    end_run(None, Ok(())).unwrap();
}

/// A request reaches the sink where NCBI writes it (`Posts::request` released at each batch's
/// start): at NCBI's error only the requests reached are written; at `Ok` or a refusal of LOSAT's
/// every request is written.
#[test]
fn requests_after_the_last_batch_start_reached_are_dropped_at_ncbis_error() {
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    let three = |export: &StrategyExport| {
        for at in 0..3 {
            export
                .write_request(at..at + 1, &request(&queries, &[]))
                .unwrap();
            export.note_searched(at..at + 1, searched(5)).unwrap();
        }
    };
    let native = || -> Result<()> {
        Err(crate::cli::NativeError {
            exit: 255,
            message: "Error: x\n".into(),
        }
        .into())
    };
    let all = [FIRST.to_string(), with_space(5), with_space(5)].concat();
    let (strategy, sink) = export(None);
    three(&strategy);
    assert_eq!(sink.text(), "");
    let marks = strategy.take_posts(usize::MAX);
    assert_eq!(marks.len(), 3);
    // The first batch's start writes the first request; the second reaches the first two.
    marks[0].1.write(&mut std::io::sink(), None).unwrap();
    assert_eq!(sink.text(), FIRST);
    marks[1].1.write(&mut std::io::sink(), None).unwrap();
    let error = end_run(Some(&strategy), native()).unwrap_err();
    assert!(error.is::<crate::cli::NativeError>());
    assert_eq!(sink.text(), [FIRST.to_string(), with_space(5)].concat());
    // A run without an error writes every request, and so does a refusal of LOSAT's.
    for result in [Ok(()), Err(anyhow::anyhow!("not supported"))] {
        let (strategy, sink) = export(None);
        three(&strategy);
        let _ = end_run(Some(&strategy), result);
        assert_eq!(sink.text(), all);
    }
}

/// The request of a PSSM query (`RequestQueries::Pssm`) as NCBI BLAST+ 2.17.0 writes it: the
/// `PssmWithParameters` as read, at the depth of `queries pssm` (its strings and octets wrapped
/// there), with psiblast's options and `Web_StepNumber` (`psiblast -in_pssm <rnd> -db
/// tests/db_search/paja`, fixture row `psiblast.pssm.db.itdef`) and PSI-TBLASTN's (`tblastn
/// -in_pssm <rnd> -db tests/db_search/tslice`, row `tblastn.db.in_pssm`).
#[test]
fn a_pssm_query_is_the_pssm_as_read() {
    let fixture =
        |path: &str| std::fs::read(format!("{}/tests/{path}", env!("CARGO_MANIFEST_DIR"))).unwrap();
    let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
        "fixtures/db_search/db.ps.asn.rnd.6.asn",
    ))
    .unwrap();
    for (kind, database, frozen) in [
        (
            crate::report::archive::options::RequestKind::Psiblast,
            "tests/db_search/paja",
            "psiblast.pssm.db.itdef.ss",
        ),
        (
            crate::report::archive::options::RequestKind::PsiTblastn,
            "tests/db_search/tslice",
            "tblastn.db.in_pssm.ss",
        ),
    ] {
        let sink = Sink::default();
        let mut inputs = OptionInputs::new(kind, "blastp");
        if kind == crate::report::archive::options::RequestKind::Psiblast {
            inputs.num_iterations = Some(1);
        }
        let export = StrategyExport::new(
            Box::new(sink.clone()),
            &inputs,
            Some(database.as_bytes()),
            OptionsHandle::default(),
        );
        let request = RequestInputs {
            queries: crate::report::archive::writer::RequestQueries::Pssm(&pssm),
            subjects: &[],
            required: None,
            lcase_masks: &[],
        };
        assert!(export.write_request(0..1, &request).unwrap());
        end_run(Some(&export), Ok(())).unwrap();
        assert_eq!(
            sink.text(),
            String::from_utf8(fixture(&format!("fixtures/strategy/{frozen}"))).unwrap(),
            "{kind:?}"
        );
    }
}

/// `-` is refused at `check_report_stream` / `write_with`; a file (including the same path as
/// `-out`) is created (emptied) with no refusal; a path that cannot be created is NCBI's argument
/// error.
#[test]
fn the_file_is_created_and_its_refusal_kept() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-strategy-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("s.asn");
    std::fs::write(&path, b"old").unwrap();
    let made = create_strategy_file(&path, "BLASTP").unwrap();
    assert!(matches!(made.target, StrategyTarget::Disk(_)) && made.refusal.is_none());
    assert_eq!(std::fs::read(&path).unwrap(), b"");
    let dash = create_strategy_file(Path::new("-"), "BLASTP").unwrap();
    assert!(matches!(dash.target, StrategyTarget::None));
    assert!(dash.refusal.as_deref().unwrap().contains("standard output"));
    // The refusal comes where the run asks for it, and before any request.
    let export = StrategyExport::for_run(
        create_strategy_file(Path::new("-"), "BLASTP").unwrap(),
        OptionInputs::new(Program::Blastp, "blastp"),
        Some(6),
        HitlistArgs::default(),
        None,
        &OptionsHandle::default(),
    );
    let refusal = check_report_stream(Some(&export)).unwrap_err();
    assert!(refusal.to_string().contains("standard output"), "{refusal}");
    assert!(export.write_with(0..1, |_| Ok(())).is_err());
    assert!(check_report_stream(None).is_ok());
    let missing = dir.join("no_such_dir").join("x");
    let error = create_strategy_file(&missing, "BLASTP").unwrap_err();
    let native = error.downcast_ref::<crate::cli::NativeError>().unwrap();
    assert_eq!(native.exit, 1);
    assert_eq!(
        native.message,
        format!(
            "Command line argument error: Argument \"export_search_strategy\". File is not accessible:  `{}'\n",
            missing.display()
        )
        .as_bytes()
    );
    // The null device as the strategy file: nothing written there is kept.
    #[cfg(unix)]
    drop(create_strategy_file(Path::new("/dev/null"), "BLASTP").unwrap());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// A run refused at NCBI's first request (`-`) that ends before it with an exception of NCBI's
/// that a report's writer carries (the `-html` links' registry, `NcbiException::into_io_error`):
/// the end takes it for NCBI's error, as the command line's exit does (`cli::NativeError::caught`,
/// exit status 255), and keeps its text, as it keeps a refusal of LOSAT's (the report stream
/// writes NCBI's prolog before either as NCBI's stream does: `check_report_stream`).
#[test]
fn an_exception_of_ncbis_from_a_writer_ends_the_run_as_ncbis_error() {
    let exception = || {
        anyhow::Error::from(
            crate::ncbi_exception::NcbiException::new(
                crate::ncbi_exception::ThrowSite {
                    file: "src/corelib/ncbireg.cpp",
                    line: 722,
                    kind: "CRegistryException::eEntry",
                    function: None,
                },
                "Invalid registry entry format",
            )
            .into_io_error(),
        )
    };
    let error = exception();
    let native = crate::cli::NativeError::caught(&error).expect("NCBI's error");
    assert_eq!(native.exit, 255);
    for error in [exception(), anyhow::anyhow!("a refusal of LOSAT's")] {
        let export = StrategyExport::for_run(
            create_strategy_file(Path::new("-"), "BLASTX").unwrap(),
            OptionInputs::new(Program::Blastx, "blastx"),
            Some(0),
            HitlistArgs::default(),
            None,
            &OptionsHandle::default(),
        );
        let message = error.to_string();
        let ended = end_run(Some(&export), Err(error)).unwrap_err();
        assert_eq!(ended.to_string(), message);
    }
}

/// The files that a receiver got (LOSAT Web's run files).
#[derive(Default)]
struct Received(Mutex<Vec<(Vec<u8>, Vec<u8>)>>);

impl FileReceiver for Received {
    fn receive(&self, name: &[u8], bytes: &[u8]) -> std::io::Result<()> {
        self.0.lock().unwrap().push((name.to_vec(), bytes.to_vec()));
        Ok(())
    }
}

impl Received {
    fn files(&self) -> Vec<(Vec<u8>, Vec<u8>)> {
        self.0.lock().unwrap().clone()
    }
}

/// LOSAT Web's strategy file (`RunFiles::Receiver`, line RG's F1w): the receiver gets it under
/// the name as typed where the run ends (`end_run`, whatever still holds the export), with the bytes
/// that the command line writes to the file;
/// a made file dropped before its export (an error after it is made) is the empty file that the
/// command line leaves; NCBI's error cuts it back as it cuts the command line's file
/// (`end_run`).
#[test]
fn a_receiver_gets_the_file_when_the_export_ends() {
    let received = Arc::new(Received::default());
    let files = RunFiles::Receiver(received.clone());
    let open = |name: &str| create_strategy_file_in(&files, Path::new(name), "BLASTP", &[]);
    let export = |file| {
        StrategyExport::for_run(
            file,
            OptionInputs::new(Program::Blastp, "blastp"),
            None,
            HitlistArgs::default(),
            Some(b"in/db/bp"),
            &OptionsHandle::default(),
        )
    };
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    let made = export(open("out/s.asn").unwrap());
    assert!(!made.refused());
    assert!(made.write_request(0..1, &request(&queries, &[])).unwrap());
    // The run's end hands the file over, while a clone of the export is still held; the drop
    // gives it no second time.
    let kept = made.clone();
    end_run(Some(&made), Ok(())).unwrap();
    assert_eq!(
        received.files(),
        [(b"out/s.asn".to_vec(), FIRST.as_bytes().to_vec())]
    );
    drop((made, kept));
    assert_eq!(received.files().len(), 1);
    // An error after the file is made: the file, empty.
    received.0.lock().unwrap().clear();
    drop(open("s.asn").unwrap());
    assert_eq!(received.files(), [(b"s.asn".to_vec(), Vec::new())]);
    // NCBI's error before the batch's start was reached: the request is cut away.
    received.0.lock().unwrap().clear();
    let cut = export(open("s.asn").unwrap());
    assert!(cut.write_request(0..1, &request(&queries, &[])).unwrap());
    let native = crate::cli::inaccessible("query", Path::new("q"));
    assert!(end_run(Some(&cut), Err(native)).is_err());
    drop(cut);
    assert_eq!(received.files(), [(b"s.asn".to_vec(), Vec::new())]);
    // A run that ends with the posts of a request not taken (an error between the request and its
    // batch's start) still gives the file at its end.
    received.0.lock().unwrap().clear();
    let early = export(open("s.asn").unwrap());
    assert!(early.write_request(0..1, &request(&queries, &[])).unwrap());
    let refusal = anyhow::anyhow!("a refusal of LOSAT's");
    assert!(end_run(Some(&early), Err(refusal)).is_err());
    drop(early);
    assert_eq!(
        received.files(),
        [(b"s.asn".to_vec(), FIRST.as_bytes().to_vec())]
    );
    // A file that does not reach the run's files is the run's internal error, not a run without
    // it; a refusal of LOSAT's keeps its own text.
    struct Closed;
    impl FileReceiver for Closed {
        fn receive(&self, _name: &[u8], _bytes: &[u8]) -> std::io::Result<()> {
            Err(std::io::Error::other("the archive was closed"))
        }
    }
    let closed = RunFiles::Receiver(Arc::new(Closed));
    let lost = |result: Result<()>| {
        let file = create_strategy_file_in(&closed, Path::new("s.asn"), "BLASTP", &[]);
        end_run(Some(&export(file.unwrap())), result)
            .unwrap_err()
            .to_string()
    };
    let error = lost(Ok(()));
    assert!(
        error.starts_with("LOSAT internal error: the -export_search_strategy file did not reach")
            && error.contains("the archive was closed"),
        "{error}"
    );
    assert_eq!(
        lost(Err(anyhow::anyhow!("a refusal of LOSAT's"))),
        "a refusal of LOSAT's"
    );
}

/// The names that LOSAT Web reads otherwise than as a file of the run (`RunFiles::strategy_target`):
/// NCBI's error at once where the command line cannot open the name; the command line's refusals of
/// `-` and of the standard output, raised at NCBI's first request with the report held; no file for
/// the null device; LOSAT Web's refusal, at the same point, of the standard error and of the file of
/// another option of the run, which names that option.
#[test]
fn a_receiver_refuses_what_the_command_line_refuses_or_has_no_file_for() {
    let received = Arc::new(Received::default());
    let files = RunFiles::Receiver(received.clone());
    let query = [RunFileName {
        option: "query",
        name: std::ffi::OsStr::new("q.fa"),
    }];
    let open = |name: &str| create_strategy_file_in(&files, Path::new(name), "BLASTX", &query);
    let error = open("/etc/s.asn").unwrap_err();
    let native = error.downcast_ref::<crate::cli::NativeError>().unwrap();
    assert_eq!(native.exit, 1);
    assert_eq!(
        native.message,
        b"Command line argument error: Argument \"export_search_strategy\". File is not accessible:  `/etc/s.asn'\n"
    );
    for (name, refusal) in [
        ("-", "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's BLASTX; give a file name"),
        ("/dev/stdout", "an -export_search_strategy file that is also the output (-out or the standard output) is not supported by LOSAT's BLASTX: NCBI writes both through separate streams and the file's bytes follow the order of their flushes"),
        ("/dev/stderr", "-export_search_strategy /dev/stderr (a device or a file of the system, where the run has no file) is not supported in the browser; give a file name of the run's own"),
        ("./q.fa", "-export_search_strategy ./q.fa (the file of -query too, whose bytes the command line's strategy file replaces or mixes with its own) is not supported in the browser; give a file name of the run's own"),
    ] {
        let file = open(name).unwrap();
        assert_eq!(file.refusal.as_deref(), Some(refusal), "{name}");
        let export = StrategyExport::for_run(
            file,
            OptionInputs::new(Program::Blastx, "blastx"),
            Some(0),
            HitlistArgs::default(),
            None,
            &OptionsHandle::default(),
        );
        assert!(export.refused(), "{name}");
        let error = check_report_stream(Some(&export)).unwrap_err();
        assert_eq!(error.to_string(), refusal);
    }
    let null = open("/dev/null").unwrap();
    assert!(matches!(null.target, StrategyTarget::None) && null.refusal.is_none());
    drop(null);
    assert!(received.files().is_empty());
}

/// A panic (LOSAT's own failure; the release profile aborts) ends the run's strategy file as a
/// refusal of LOSAT's ends it (`install_panic_hook`): every live export is kept for the hook
/// (`HOLDING`), whose writing gives its sink every request made, even while a request is being
/// made.
#[test]
fn a_panic_writes_the_requests_held() {
    let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
    let made = || {
        let (export, sink) = export(None);
        for at in 0..3 {
            export
                .write_request(at..at + 1, &request(&queries, &[]))
                .unwrap();
            export.note_searched(at..at + 1, searched(5)).unwrap();
        }
        (export, sink)
    };
    let kept = |export: &StrategyExport| {
        HOLDING
            .lock()
            .unwrap()
            .iter()
            .any(|shared| shared.upgrade().is_some_and(|s| Arc::ptr_eq(&s, &export.0)))
    };
    let (strategy, sink) = made();
    assert!(kept(&strategy));
    assert_eq!(sink.text(), "");
    // A panic while a request is made (the export's state locked): the requests made before it
    // are written (the hook does not take that lock), the one being made is not held yet. Only
    // this export: other tests' exports of the process stay as they are.
    let all = [FIRST.to_string(), with_space(5), with_space(5)].concat();
    strategy
        .write_with(3..4, |write| {
            write_held_of(&[Arc::downgrade(&strategy.0)]);
            assert_eq!(sink.text(), all);
            write(&request(&queries, &[]))
        })
        .unwrap();
    assert_eq!(sink.text(), all);
    // The run's end after it writes the request made since, and nothing twice.
    let _ = end_run(Some(&strategy), Ok(()));
    let all = [all, with_space(5)].concat();
    assert_eq!(sink.text(), all);
    // Every lock that the hook takes was given back (`hook_lock`): a panic after it waits for none
    // that its own thread holds.
    assert_eq!(HOOK_LOCKS.get(), 0);
}

/// The panic hook of a thread that holds one of the hook's locks (a panic inside such a step)
/// returns at once instead of waiting for it (`write_held_at_panic`, `hook_lock`), and the count
/// of such locks goes back to none when the lock is given back.
#[test]
fn the_panic_hook_does_not_wait_for_a_lock_its_own_thread_holds() {
    let lock = Mutex::new(());
    {
        let _held = hook_lock(&lock);
        assert_eq!(HOOK_LOCKS.get(), 1);
        write_held_at_panic();
    }
    assert_eq!(HOOK_LOCKS.get(), 0);
}

/// `RunFiles::strategy_target` checks `-logfile` against the run's real files (line WQA D19 round
/// R1): when `out` does not exist, `-logfile out/x.asn` stays on stderr (`Unopened`), so
/// `-export_search_strategy out/x.asn` is not refused as the file of `-logfile`; when `out` exists
/// in the run's files, `-logfile out/x.asn` opens that file and the strategy is refused.
#[test]
fn strategy_target_reads_logfile_against_runs_files() {
    struct WithRegistered(Vec<Vec<u8>>);
    impl crate::report::run_files::FileReceiver for WithRegistered {
        fn receive(&self, _name: &[u8], _bytes: &[u8]) -> std::io::Result<()> {
            Ok(())
        }
        fn registered_files(&self) -> &[Vec<u8>] {
            &self.0
        }
    }
    let others = [RunFileName {
        option: "logfile",
        name: std::ffi::OsStr::new("out/x.asn"),
    }];
    let without_out = RunFiles::Receiver(Arc::new(WithRegistered(Vec::new())));
    let file =
        create_strategy_file_in(&without_out, Path::new("out/x.asn"), "TBLASTX", &others).unwrap();
    assert!(
        file.refusal.is_none(),
        "without folder `out`, -logfile out/x.asn is unopened so strategy is a file: {:?}",
        file.refusal
    );
    assert!(matches!(file.target, StrategyTarget::Received(_)));

    let with_out = RunFiles::Receiver(Arc::new(WithRegistered(vec![b"out/keep".to_vec()])));
    let refused =
        create_strategy_file_in(&with_out, Path::new("out/x.asn"), "TBLASTX", &others).unwrap();
    assert_eq!(
        refused.refusal.as_deref(),
        Some(
            "-export_search_strategy out/x.asn (the file of -logfile too, whose bytes the \
             command line's strategy file replaces or mixes with its own) is not supported in \
             the browser; give a file name of the run's own"
        )
    );
}
