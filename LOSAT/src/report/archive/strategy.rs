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

use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};

use anyhow::Result;

use super::masks::MaskEntry;
use super::options::{OptionInputs, Param};
use super::writer::{text_out, write_strategy_document, BatchArchive};
use super::{ArchiveEncoding, ArchiveRequest};

struct State {
    sink: Box<dyn Write + Send>,
    /// The `EffectiveSearchSpace` option after the last written batch
    /// (`ArchiveHandle::search_space_option`).
    last_search_space: Option<i64>,
    /// The `TaxidList` or `NegativeTaxidList` program option (`set_database_filter`).
    database_list: Option<Param>,
    /// A list kept no sequence and `-searchsp` is not given (`set_database_filter`).
    empty_database: bool,
}

struct Shared {
    request: ArchiveRequest,
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
    use std::os::unix::fs::MetadataExt;
    let Ok(own) = file.metadata() else {
        return false;
    };
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

impl StrategyExport {
    /// The export of a run of `inputs.program` into `sink` (`create_strategy_file`). `database`
    /// is the `-db` argument as NCBI keeps it (`None` for `-subject`); `label` names the program
    /// in a refusal. The request is the archive's (`ArchiveRequest::for_inputs`), with the same
    /// refusal of a `-db` name or `-window_masker_db` path that ASN.1 text cannot hold.
    pub fn new(
        sink: Box<dyn Write + Send>,
        inputs: &OptionInputs,
        database: Option<&str>,
        label: &str,
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
            state: Mutex::new(State {
                sink,
                last_search_space: None,
                database_list: None,
                empty_database: false,
            }),
        })))
    }

    /// `new` into the created strategy file (`create_strategy_file`), buffered.
    pub fn for_file(
        file: std::fs::File,
        inputs: &OptionInputs,
        database: Option<&str>,
        label: &str,
    ) -> Result<Self> {
        Self::new(
            Box::new(std::io::BufWriter::new(file)),
            inputs,
            database,
            label,
        )
    }

    fn with_state<R>(&self, f: impl FnOnce(&mut State) -> R) -> R {
        let mut state = self
            .0
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        f(&mut state)
    }

    /// `ArchiveHandle::set_database_filter` for the strategy's requests.
    pub fn set_database_filter(&self, list: Option<Param>, empty_database: bool) {
        self.with_state(|state| {
            state.database_list = list;
            state.empty_database = empty_database;
        });
    }

    /// Writes the request of one query batch: what NCBI writes before the batch's search.
    /// `batch` is the batch as its archive takes it (the results' parts are not read);
    /// `lcase_masks` are the `LCaseMask` parameters before the search (the lower-case runs
    /// of `-lcase_masking`; BLASTN's archive adds DUST intervals to `batch.lcase_masks`).
    /// `batch.searched_space` then becomes the `EffectiveSearchSpace` option of the requests
    /// after this one, as the archive's handle keeps it (`ArchiveHandle::search_space_option`).
    pub fn write_batch(&self, batch: &BatchArchive<'_>, lcase_masks: &[MaskEntry]) -> Result<()> {
        self.with_state(|state| {
            let mut buf = Vec::new();
            let mut out = text_out(&mut buf, &mut *state.sink);
            write_strategy_document(
                &mut out,
                &self.0.request,
                batch,
                lcase_masks,
                state.last_search_space,
                state.database_list.as_ref(),
            );
            let written = out.finish();
            // The file is complete after every request: the command line's exit does not run
            // destructors.
            let flushed = state.sink.flush();
            if batch.searched_space.is_some() {
                state.last_search_space = if state.empty_database {
                    Some(0)
                } else {
                    batch.searched_space
                };
            }
            written?;
            flushed?;
            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::app::Program;
    use crate::report::archive::bioseq::BioseqInput;
    use crate::report::archive::masks::Frame;
    use crate::report::archive::KaBlock;

    #[derive(Clone, Default)]
    struct Shared(Arc<Mutex<Vec<u8>>>);

    impl Write for Shared {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(bytes);
            Ok(bytes.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    fn batch<'a>(
        queries: &'a [BioseqInput<'a>],
        subjects: &'a [BioseqInput<'a>],
        searched_space: Option<i64>,
    ) -> BatchArchive<'a> {
        BatchArchive {
            queries,
            subjects,
            required: None,
            lcase_masks: &[],
            masks: &[],
            ka_blocks: [KaBlock::missing(false), KaBlock::missing(true)],
            search_space: 0,
            length_adjustment: 0,
            searched_space,
        }
    }

    fn protein<'a>(id: &'a str, letters: &'a [u8]) -> BioseqInput<'a> {
        BioseqInput {
            id,
            ids: None,
            seq_loc: None,
            defline: b">q a query",
            title: b"q a query",
            letters,
            protein: true,
        }
    }

    // NCBI BLAST+ 2.17.0, `blastp -query q.faa -db in/db/bp -export_search_strategy f` with a
    // query `>q a query` / `MWCHWCMWCHWC` (the layout of the RG-F1 probe's blastp.d.plain.ss).
    const FIRST: &str = "Blast4-request ::= {\n  body queue-search {\n    program \"blastp\",\n    service \"plain\",\n    queries bioseq-set {\n      seq-set {\n        seq {\n          id {\n            local str \"Query_1\"\n          },\n          descr {\n            user {\n              type str \"CFastaReader\",\n              data {\n                {\n                  label str \"DefLine\",\n                  data str \">q a query\"\n                }\n              }\n            },\n            title \"q a query\"\n          },\n          inst {\n            repr raw,\n            mol aa,\n            length 12,\n            seq-data ncbieaa \"MWCHWCMWCHWC\"\n          }\n        }\n      }\n    },\n    subject database \"in/db/bp\",\n    algorithm-options {\n      {\n        name \"EvalueThreshold\",\n        value cutoff e-value { 1, 10, 1 }\n      },\n      {\n        name \"MaskAtHash\",\n        value boolean FALSE\n      },\n      {\n        name \"SegFiltering\",\n        value boolean FALSE\n      },\n      {\n        name \"WordThreshold\",\n        value integer 11\n      },\n      {\n        name \"WindowSize\",\n        value integer 40\n      },\n      {\n        name \"HitlistSize\",\n        value integer 500\n      },\n      {\n        name \"UngappedMode\",\n        value boolean FALSE\n      },\n      {\n        name \"CompositionBasedStats\",\n        value integer 2\n      },\n      {\n        name \"SmithWatermanMode\",\n        value boolean FALSE\n      }\n    }\n  }\n}\n";

    /// One request per batch, without `ident`; the `EffectiveSearchSpace` option is absent in
    /// the first request and is the batch before's in the next ones, also after a batch that
    /// was not searched.
    #[test]
    fn each_batch_is_a_request_with_the_search_space_of_the_batch_before() {
        let sink = Shared::default();
        let inputs = OptionInputs::new(Program::Blastp, "blastp");
        let export =
            StrategyExport::new(Box::new(sink.clone()), &inputs, Some("in/db/bp"), "BLASTP")
                .unwrap();
        let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
        export
            .write_batch(&batch(&queries, &[], Some(817_380)), &[])
            .unwrap();
        let first = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
        assert_eq!(first, FIRST);

        export
            .write_batch(&batch(&queries, &[], None), &[])
            .unwrap();
        export
            .write_batch(&batch(&queries, &[], Some(5)), &[])
            .unwrap();
        export
            .write_batch(&batch(&queries, &[], None), &[])
            .unwrap();
        let all = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
        let with_space = |space: i64| {
            FIRST.replace(
                "        value boolean FALSE\n      }\n    }\n  }\n}\n",
                &format!(
                    "        value boolean FALSE\n      }},\n      {{\n        name \"EffectiveSearchSpace\",\n        value big-integer {space}\n      }}\n    }}\n  }}\n}}\n"
                ),
            )
        };
        assert_eq!(
            all,
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
        let sink = Shared::default();
        let mut inputs = OptionInputs::new(Program::Blastp, "blastp");
        inputs.searchsp = Some(1_000_000);
        let export =
            StrategyExport::new(Box::new(sink.clone()), &inputs, Some("in/db/bp"), "BLASTP")
                .unwrap();
        let queries = [protein("Query_1", b"MWCHWCMWCHWC")];
        let mask = MaskEntry {
            id: "Query_1".to_string(),
            best: None,
            frame: Frame::NotSet,
            intervals: vec![(2, 5)],
        };
        export
            .write_batch(
                &batch(&queries, &[], Some(1_000_000)),
                std::slice::from_ref(&mask),
            )
            .unwrap();
        export.set_database_filter(None, true);
        export
            .write_batch(&batch(&queries, &[], Some(77)), &[])
            .unwrap();
        export
            .write_batch(&batch(&queries, &[], None), &[])
            .unwrap();
        let all = String::from_utf8(sink.0.lock().unwrap().clone()).unwrap();
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
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
