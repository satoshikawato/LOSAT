//! The search strategy file (`-export_search_strategy`): the `Blast4-request` of each query
//! batch, in ASN.1 text.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:636-651
//! ```c++
//! void
//! SaveSearchStrategy(const CArgs& args,
//!                    blast::CBlastAppArgs* cmdline_args,
//!                    CRef<blast::IQueryFactory> queries,
//!                    CRef<blast::CBlastOptionsHandle> opts_hndl,
//!                    CRef<objects::CPssmWithParameters> pssm
//!                      /* = CRef<objects::CPssmWithParameters>() */,
//!                     unsigned int num_iters /* =0 */)
//! {
//!     CNcbiOstream* out = cmdline_args->GetExportSearchStrategyStream(args);
//!     if ( !out ) {
//!         return;
//!     }
//!
//!     s_ExportSearchStrategy(out, queries, opts_hndl,
//!                            cmdline_args->GetBlastDatabaseArgs(),
//!                            pssm, num_iters);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:366-369
//! ```c++
//! void CExportStrategy::ExportSearchStrategy_ASN1(CNcbiOstream* out)
//! {
//! 	*out << MSerial_AsnText << *GetSearchStrategy();
//! }
//! ```
//! The applications call `SaveSearchStrategy` for every query batch, before the batch's
//! search (blastp_app.cpp:262, blastn_app.cpp:282, blastx_app.cpp:264, tblastn_app.cpp:280,
//! tblastx_app.cpp:181): the file holds one `Blast4-request` per batch. `CExportStrategy` also
//! builds the request of the archive (outfmt 11, build_archive.cpp:156-200), after the search.
//! The strategy's request therefore differs from the archive's of the same batch in what the
//! search changes (NCBI BLAST+ 2.17.0, 418 runs of the RG-F1 probe):
//!
//! - no `ident` (the applications pass no client id, blast_app_util.cpp:335-337);
//! - `EffectiveSearchSpace` is the option as it stands before the batch: `-searchsp` as typed in
//!   the first request, else absent there, and the value of the batch before in the later ones
//!   (the options handle is shared by the batches);
//! - blastn's `LCaseMask` parameters hold the lower-case runs of `-lcase_masking` only; the
//!   archive's also hold the DUST intervals that the search added to the query's masks.
//!
//! `ARCHIVE_FORMAT` does not change the file: it is always ASN.1 text.
//!
//! Which requests the file holds when a run stops follows from where they are written: the
//! programs that take their batches one at a time (BLASTN, TBLASTX, TBLASTN, BLASTX) write a
//! batch's request where NCBI does, once the batch is read and is not empty and before any error
//! of its set-up, search or reports is acted on (TBLASTN, TBLASTX and BLASTX may have searched the
//! batch ahead; only the order of the file and of what the run writes matters), so a batch that
//! fails after that point has its request and the batches after it have none; BLASTP, which
//! searches the batches of a call together, releases them where the call ends, up to the batch
//! that ends it. One cursor of the run (`State::next_query`), keyed by the run's index of each
//! batch's first query, keeps every request written once and in order (a batch searched again
//! writes none), and the batch's searched space moves the next request's `EffectiveSearchSpace`
//! option (`note_searched`). A broken order is LOSAT's internal error, not a file that differs.

use std::io::Write;
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, Mutex};

use anyhow::Result;

use super::options::OptionInputs;
use super::writer::{text_out, write_strategy_document, RequestInputs};
use super::{ArchiveEncoding, ArchiveRequest, OptionsHandle, SearchSpaceCarry};

struct State {
    sink: Box<dyn Write + Send>,
    /// The `EffectiveSearchSpace` option after the batches noted so far (`note_searched`).
    carry: SearchSpaceCarry,
    /// The run's searched queries whose batches' requests are written: the next request starts
    /// here.
    next_query: usize,
    /// The run's searched queries whose batches' searched spaces moved `carry`.
    noted_through: usize,
    /// No request starts at or after this query (`stop_after`; `usize::MAX` without a stop).
    stop: usize,
}

struct Shared {
    request: ArchiveRequest,
    /// The run's options handle, shared with the run's archive (`OptionsHandle`).
    options: OptionsHandle,
    state: Mutex<State>,
}

/// What a run with `-export_search_strategy` holds: the request of its options and the file.
#[derive(Clone)]
pub struct StrategyExport(Arc<Shared>);

impl std::fmt::Debug for StrategyExport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StrategyExport").finish_non_exhaustive()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3560-3569
// ```c++
// CNcbiOstream*
// CSearchStrategyArgs::GetExportStream(const CArgs& args) const
// {
//     CNcbiOstream* retval = NULL;
//     if (args.Exist(kArgOutputSearchStrategy) &&
//         args[kArgOutputSearchStrategy].HasValue()) {
//         retval = &args[kArgOutputSearchStrategy].AsOutputFile();
//     }
//     return retval;
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
// ```c++
//     if (GetExportSearchStrategyStream(args) ||
//            m_FormattingArgs->ArchiveFormatRequested(args)) {
//         locality = CBlastOptions::eBoth;
//     }
// ```
/// Creates (and empties) the strategy file where NCBI's `SetOptions` does: after the argument
/// parser accepted the command line and before any option handler runs, so before the
/// `-outfmt` text is read, the database list is resolved, and `-subject`, `-query` and `-out`
/// are opened. A file that cannot be created is NCBI's argument error (exit status 1).
///
/// `out` is the `-out` argument. LOSAT refuses a strategy file that is the standard output
/// (`-`) or the same file as the standard output or `-out`: NCBI writes each request to the
/// output's stream before the batch's report (`-`), or to a second stream on the same file,
/// whose bytes then replace or are replaced by the report's as the two buffers are flushed.
pub fn create_strategy_file(path: &Path, out: Option<&Path>, label: &str) -> Result<std::fs::File> {
    if path.as_os_str() == "-" {
        anyhow::bail!(
            "-export_search_strategy - (the search strategy on the standard output, between the parts of the report) is not supported by LOSAT's {label}; give a file name"
        );
    }
    let file = std::fs::File::create(path)
        .map_err(|_| crate::cli::inaccessible("export_search_strategy", path))?;
    if is_same_as_output(&file, path, out) {
        anyhow::bail!(
            "an -export_search_strategy file that is also the output (-out or the standard output) is not supported by LOSAT's {label}: NCBI writes both through separate streams and the file's bytes follow the order of their flushes"
        );
    }
    Ok(file)
}

/// Whether the created strategy `file` is the file of `-out` (which exists then: NCBI creates
/// it later) or of the standard output.
#[cfg(unix)]
fn is_same_as_output(file: &std::fs::File, _path: &Path, out: Option<&Path>) -> bool {
    use std::os::fd::AsFd;
    use std::os::unix::fs::{FileTypeExt, MetadataExt};
    let Ok(own) = file.metadata() else {
        return false;
    };
    // The null device: what NCBI writes there is lost whatever the order of the writes.
    if own.file_type().is_char_device()
        && std::fs::metadata("/dev/null").is_ok_and(|null| null.rdev() == own.rdev())
    {
        return false;
    }
    let same = |other: std::io::Result<std::fs::Metadata>| {
        other.is_ok_and(|other| other.dev() == own.dev() && other.ino() == own.ino())
    };
    let stdout_is_report = out.is_none_or(|out| out.as_os_str() == "-");
    match out.filter(|out| out.as_os_str() != "-") {
        Some(out) if same(std::fs::metadata(out)) => return true,
        _ => {}
    }
    stdout_is_report
        && std::io::stdout()
            .as_fd()
            .try_clone_to_owned()
            .is_ok_and(|fd| same(std::fs::File::from(fd).metadata()))
}

#[cfg(not(unix))]
fn is_same_as_output(_file: &std::fs::File, path: &Path, out: Option<&Path>) -> bool {
    let (Some(out), Ok(own)) = (
        out.filter(|out| out.as_os_str() != "-"),
        std::fs::canonicalize(path),
    ) else {
        return false;
    };
    std::fs::canonicalize(out).is_ok_and(|out| out == own)
}

/// The arguments that make the `HitlistSize` option of a run's requests besides its format
/// (`StrategyExport::for_run`).
#[derive(Clone, Copy, Debug, Default)]
pub struct HitlistArgs {
    pub max_target_seqs: Option<usize>,
    pub num_descriptions: Option<usize>,
    pub num_alignments: Option<usize>,
}

/// The `HitlistSize` option of a run of format `outfmt` (`report_counts`): the archive's is that
/// of outfmt 11, while outfmt 0 to 4 take the larger of `-num_descriptions` and `-num_alignments`
/// (blast_args.cpp:2910-2936, probe `nd` of RG-F1).
fn hitlist_size(outfmt: i32, counts: HitlistArgs) -> i32 {
    crate::blastinput::report_counts::report_counts(
        outfmt,
        counts.max_target_seqs,
        counts.num_descriptions,
        counts.num_alignments,
    )
    .hitlist_size as i32
}

/// The strategy file of a run's `-export_search_strategy` argument (`path`), created where NCBI
/// creates it (`create_strategy_file`), or `None` without the argument. A file name that is not
/// UTF-8 is refused first (`check_utf8_file_name`). `out` is the `-out` argument; `label` names
/// the program in a refusal.
pub fn open_strategy_file(
    path: Option<&Path>,
    out: Option<&Path>,
    label: &str,
) -> Result<Option<std::fs::File>> {
    path.map(|path| {
        crate::blastinput::input_files::check_utf8_file_name(
            path,
            "export_search_strategy",
            label,
        )?;
        create_strategy_file(path, out, label)
    })
    .transpose()
}

impl StrategyExport {
    /// The export of a run of `inputs.program` into `sink` (`create_strategy_file`). `database`
    /// is the `-db` argument as NCBI keeps it (`None` for `-subject`); `label` names the program
    /// in a refusal; `options` is the run's options handle, which its archive shares.
    /// The request is the archive's (`ArchiveRequest::for_inputs`), with the same refusal of a
    /// `-db` name or `-window_masker_db` path that ASN.1 text cannot hold.
    pub fn new(
        sink: Box<dyn Write + Send>,
        inputs: &OptionInputs,
        database: Option<&str>,
        label: &str,
        options: OptionsHandle,
    ) -> Result<Self> {
        let request = ArchiveRequest::for_inputs(
            ArchiveEncoding::Text,
            inputs,
            database,
            label,
            "-export_search_strategy",
        )?;
        Ok(Self(Arc::new(Shared {
            request,
            options,
            state: Mutex::new(State {
                sink,
                carry: SearchSpaceCarry::default(),
                next_query: 0,
                noted_through: 0,
                stop: usize::MAX,
            }),
        })))
    }

    /// The export of a run into its strategy file (`open_strategy_file`), buffered: the program's
    /// recorded options `inputs` with the `HitlistSize` of the run's format `outfmt` and `counts`
    /// (the archive's is outfmt 11's), and the run's options handle `options` (made once per run
    /// and given to its archive too, `ArchiveHandle::for_program`).
    pub fn for_run(
        file: std::fs::File,
        mut inputs: OptionInputs,
        outfmt: i32,
        counts: HitlistArgs,
        database: Option<&str>,
        label: &str,
        options: &OptionsHandle,
    ) -> Result<Self> {
        inputs.hitlist_size = hitlist_size(outfmt, counts);
        Self::new(
            Box::new(std::io::BufWriter::new(file)),
            &inputs,
            database,
            label,
            options.clone(),
        )
    }

    /// The run's options handle.
    pub fn options(&self) -> &OptionsHandle {
        &self.0.options
    }

    fn with_state<R>(&self, f: impl FnOnce(&mut State) -> R) -> R {
        let mut state = self
            .0
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        f(&mut state)
    }

    /// Whether the request of the batch of the run's queries `queries` is still to be written:
    /// not one already written (a batch that a run searches again, or a later part of a batch the
    /// search split) and none at or after the run's stop (`stop_after`).
    pub fn takes(&self, queries: Range<usize>) -> bool {
        self.with_state(|state| state.takes(&queries))
    }

    /// The run's searched queries whose batches' requests are written.
    pub fn written_through(&self) -> usize {
        self.with_state(|state| state.next_query)
    }

    /// Writes the request of the batch of the run's searched queries `queries`, where NCBI writes
    /// it: after the batch is read and made a query factory, before its search. `make` builds the
    /// request (from the batch's records, `requests::RunRequests`; it may refuse first) and hands
    /// it to its argument, which writes it with the `EffectiveSearchSpace` option the batches
    /// before left (`note_searched`), all under one lock. Returns false, calling nothing, when the
    /// request is not to be written (`takes`). The requests are written in the run's order, each
    /// starting where the one before ended, after the batch before was noted; anything else is
    /// LOSAT's internal error.
    pub fn write_with(
        &self,
        queries: Range<usize>,
        make: impl FnOnce(&mut dyn FnMut(&RequestInputs<'_>) -> Result<()>) -> Result<()>,
    ) -> Result<bool> {
        self.with_state(|state| {
            if !state.takes(&queries) {
                return Ok(false);
            }
            anyhow::ensure!(
                queries.start == state.next_query && state.noted_through == state.next_query,
                "LOSAT internal error: the search strategy's request of the queries {queries:?} does not follow the requests written ({}) and noted ({})",
                state.next_query,
                state.noted_through
            );
            let shared = &self.0;
            make(&mut |request| {
                let mut buf = Vec::new();
                let mut out = text_out(&mut buf, &mut *state.sink);
                write_strategy_document(
                    &mut out,
                    &shared.request,
                    request,
                    state.carry.option(),
                    &shared.options.database_params(),
                );
                let written = out.finish();
                // The file is complete after every request: the command line's exit does not run
                // destructors.
                let flushed = state.sink.flush();
                state.next_query = queries.end;
                written?;
                flushed?;
                Ok(())
            })?;
            Ok(true)
        })
    }

    /// `write_with` of a request already made.
    pub fn write_request(
        &self,
        queries: Range<usize>,
        request: &RequestInputs<'_>,
    ) -> Result<bool> {
        self.write_with(queries, |write| write(request))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:177-181
    // ```c++
    //     vector<Int8> eff_searchsp;
    //     for (size_t index = 0; index <= (size_t)qinfo->last_context; index++) {
    //         eff_searchsp.push_back(calc.GetEffSearchSpaceForContext(index));
    //     }
    //     options->SetEffectiveSearchSpace(eff_searchsp);
    // ```
    /// Moves the `EffectiveSearchSpace` option of the requests after the batch of the run's
    /// searched queries `queries` past that batch, once the batch's set-up has written it into the
    /// options handle (above): `searched` is the batch's option when NCBI searched the batch (the
    /// archive's `BatchArchive::searched_space`, from the same values), `None` for a batch NCBI does
    /// not search (`OptionsHandle::pass`). A batch noted already (a search of the batches again) is
    /// not passed again; a batch noted out of the run's order is LOSAT's internal error.
    pub fn note_searched(&self, queries: Range<usize>, searched: Option<i64>) -> Result<()> {
        self.with_state(|state| {
            if queries.end <= state.noted_through {
                return Ok(());
            }
            anyhow::ensure!(
                queries.start == state.noted_through,
                "LOSAT internal error: the search strategy's batch of the queries {queries:?} is noted after the batches through {}",
                state.noted_through
            );
            self.0.options.pass(&mut state.carry, searched);
            state.noted_through = queries.end;
            Ok(())
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
    // ```c++
    //     else {
    //         string prog_name = Blast_ProgramNameFromType(
    //                                            m_Options.GetProgramType());
    //         if (prog_name != "blastn" && prog_name != "megablast") {
    //             NCBI_THROW(blast::CBlastException, eInvalidArgument,
    //                         "unsupported score matrix");
    // ```
    /// NCBI stops the run when it formats the batch whose queries end at `queries_end`, after
    /// that batch's request (the BLAST XML report data of a matrix name it does not know,
    /// `blast_xml::matrix::xml_formatting_fails`): no request of a later batch is written.
    pub fn stop_after(&self, queries_end: usize) {
        self.with_state(|state| state.stop = state.stop.min(queries_end));
    }
}

impl State {
    /// `StrategyExport::takes`.
    fn takes(&self, queries: &Range<usize>) -> bool {
        queries.end > self.next_query && queries.start < self.stop
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::app::Program;
    use crate::report::archive::bioseq::BioseqInput;
    use crate::report::archive::masks::{Frame, MaskEntry};
    use crate::report::archive::ArchiveHandle;

    /// The strategy file's bytes, which the test reads while the export writes them.
    #[derive(Clone, Default)]
    struct Sink(Arc<Mutex<Vec<u8>>>);

    impl Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    impl Sink {
        fn text(&self) -> String {
            String::from_utf8(self.0.lock().unwrap().clone()).unwrap()
        }
    }

    fn request<'a>(
        queries: &'a [BioseqInput<'a>],
        lcase_masks: &'a [MaskEntry],
    ) -> RequestInputs<'a> {
        RequestInputs {
            queries,
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
            Some("in/db/bp"),
            "BLASTP",
            OptionsHandle::default(),
        )
        .unwrap();
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
        assert_eq!(sink.text(), FIRST);
        export.note_searched(0..1, Some(817_380)).unwrap();
        for (at, searched) in [(1, None), (2, Some(5)), (3, None)] {
            assert!(export
                .write_request(at..at + 1, &request(&queries, &[]))
                .unwrap());
            export.note_searched(at..at + 1, searched).unwrap();
        }
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
        export.note_searched(0..1, Some(1_000_000)).unwrap();
        export.options().set_database_filter(None, true).unwrap();
        export.write_request(1..2, &request(&queries, &[])).unwrap();
        export.note_searched(1..2, Some(77)).unwrap();
        export.write_request(2..3, &request(&queries, &[])).unwrap();
        export.note_searched(2..3, None).unwrap();
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

    /// A batch searched again (a later part of a batch, a call that searches the batches before
    /// again) writes no second request and does not move the search space again.
    #[test]
    fn a_request_written_once_is_not_written_again() {
        let (export, sink) = export(None);
        let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
        assert!(export.takes(0..2));
        assert!(export.write_request(0..2, &request(&queries, &[])).unwrap());
        export.note_searched(0..2, Some(5)).unwrap();
        // The same batch, and a part of it: not taken, nothing written.
        assert!(!export.takes(0..2));
        assert!(!export.takes(1..2));
        assert!(!export.write_request(0..2, &request(&queries, &[])).unwrap());
        export.note_searched(0..2, Some(9)).unwrap();
        assert_eq!(sink.text(), FIRST);
        assert!(export.write_request(2..3, &request(&queries, &[])).unwrap());
        assert_eq!(sink.text(), [FIRST.to_string(), with_space(5)].concat());
    }

    /// NCBI's stop when it formats a batch: the requests of the batches after it are not written.
    #[test]
    fn no_request_follows_the_stop() {
        let (export, sink) = export(None);
        let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
        export.stop_after(1);
        assert!(export.write_request(0..1, &request(&queries, &[])).unwrap());
        export.note_searched(0..1, Some(5)).unwrap();
        assert!(!export.takes(1..2));
        assert!(!export.write_request(1..2, &request(&queries, &[])).unwrap());
        export.note_searched(1..2, Some(6)).unwrap();
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
            ArchiveEncoding::Text,
            ArchiveRequest {
                program: "blastp",
                service: "plain",
                algorithm: Vec::new(),
                database: None,
            },
            OptionsHandle::default(),
        );
        let sink = Sink::default();
        let export = StrategyExport::new(
            Box::new(sink.clone()),
            &OptionInputs::new(Program::Blastp, "blastp"),
            Some("in/db/bp"),
            "BLASTP",
            handle.options().clone(),
        )
        .unwrap();
        assert!(export.options().same(handle.options()));
        let list = crate::report::archive::options::Param {
            name: "TaxidList",
            value: crate::report::archive::options::ParamValue::IntegerList(vec![9606]),
        };
        handle
            .options()
            .set_database_filter(Some(list), true)
            .unwrap();
        let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
        export.write_request(0..1, &request(&queries, &[])).unwrap();
        export.note_searched(0..1, Some(5)).unwrap();
        export.write_request(1..2, &request(&queries, &[])).unwrap();
        let all = sink.text();
        assert_eq!(all.matches("name \"TaxidList\"").count(), 2, "{all}");
        assert!(all.contains("value big-integer 0\n"), "{all}");
        // The archive's carry is its own: its first archive's option is its batch's.
        assert_eq!(handle.search_space_option(None), None);
        assert_eq!(handle.search_space_option(Some(5)), Some(0));
    }

    /// `-` and a file that is also `-out` are refused; a path that cannot be created is NCBI's
    /// argument error.
    #[test]
    fn the_file_is_created_or_refused() {
        let dir = std::env::temp_dir().join(format!("losat-strategy-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("s.asn");
        std::fs::write(&path, b"old").unwrap();
        drop(create_strategy_file(&path, None, "BLASTP").unwrap());
        assert_eq!(std::fs::read(&path).unwrap(), b"");
        let refusal = create_strategy_file(Path::new("-"), None, "BLASTP").unwrap_err();
        assert!(refusal.to_string().contains("standard output"), "{refusal}");
        let same = create_strategy_file(&path, Some(&path), "BLASTP").unwrap_err();
        assert!(same.to_string().contains("also the output"), "{same}");
        let other = dir.join("o.txt");
        drop(create_strategy_file(&path, Some(&other), "BLASTP").unwrap());
        let missing = dir.join("no_such_dir").join("x");
        let error = create_strategy_file(&missing, None, "BLASTP").unwrap_err();
        let native = error.downcast_ref::<crate::cli::NativeError>().unwrap();
        assert_eq!(native.exit, 1);
        assert_eq!(
            native.message,
            format!(
                "Command line argument error: Argument \"export_search_strategy\". File is not accessible:  `{}'\n",
                missing.display()
            )
        );
        // The null device as the strategy file and as `-out`: nothing written there is kept.
        #[cfg(unix)]
        drop(
            create_strategy_file(
                Path::new("/dev/null"),
                Some(Path::new("/dev/null")),
                "BLASTP",
            )
            .unwrap(),
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
