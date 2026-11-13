//! TBLASTN's application (NCBI c++/src/app/blast/tblastn_app.cpp:114-361, `Run` and
//! `x_RunMTBySplitDB`): the command line's `TblastnArgs::run`, web ABI v2's `run_local`, and the
//! search over the query batches that both reach. The options are in `args.rs`.
use anyhow::{Context, Result};
use std::io::Write;

use super::args::{search_settings, ResolvedTblastnArgs, SearchSettings, TblastnArgs};
use super::database::TblastnDatabase;
use super::stage_d_pipeline::SetupEngineError;
use super::stage_d_pipeline::{
    run_local_for_report_threads, LocalStageDProfile, SubjectSetNumbers,
};
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{render, ShownSubjects};
use super::x_subject_store::SubjectLetters;
use crate::algorithm::blastn::blast_engine::QueryReading;
use crate::algorithm::common::database::{HeaderReadError, PreliminaryReadError};
use crate::algorithm::common::x_db_stream::XReadError;
use crate::api::local_blast::{OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::tabular::SubjectFetch;

/// The failure of a batch whose search failed where NCBI's fails it with its engine error: a
/// read of the database (`common::x_db_stream::XReadError`), as the check's error fails the first
/// batch with a valid context, when the reader's error is NCBI's own, and the batch's setup
/// (`stage_d_pipeline::SetupEngineError`: a PSSM query whose `kbp_psi` is not valid, the average
/// subject length). LOSAT's refusal (and any other
/// error) is returned as `Err` to stop the run.
fn batch_failure(error: anyhow::Error) -> std::result::Result<anyhow::Error, anyhow::Error> {
    if let Some(SetupEngineError(text)) = error.downcast_ref::<SetupEngineError>() {
        return Ok(crate::blastinput::app::engine_error(text));
    }
    let error = match error.downcast::<XReadError>() {
        Ok(XReadError(read)) => read,
        Err(error) => return Err(error),
    };
    match error.downcast::<PreliminaryReadError>() {
        Ok(read_error) => Ok(crate::blastinput::app::engine_error(read_error.0)),
        Err(error) => Err(error),
    }
}

impl TblastnArgs {
    // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
    // results = lcl_blast.Run();
    // formatter.PrintOneResultSet(**result, query);
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1458:
    // PrintOneResultSet dispatches the complete result to outfmt 0/6/7.
    /// `settings` are the NCBI application settings that LOSAT reproduces
    /// (`ncbi_environment::check_ncbi_application_settings`): the readers of the subjects
    /// and of the queries (`fasta_reader`) use them.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
    /// ```c++
    ///     CNcbiApplication* app = CNcbiApplication::Instance();
    ///     if (app) {
    ///         const CNcbiRegistry& registry = app->GetConfig();
    ///         x_LoadDataLoadersConfig(registry);
    ///         x_LoadBlastDbDataLoaderConfig(registry);
    ///     }
    /// ```
    pub fn run(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    ) -> Result<()> {
        self.run_io(settings, RunIo::Files)
    }

    /// `run` over the inputs of web ABI v1 (`web_api::v1_tblastn`): the query and the subjects
    /// are the call's bytes, read as the regular files with the same bytes, and the report is
    /// returned instead of written to standard output. The `query` and `subject` of `self` are
    /// only the names the report shows for the inputs; `out` is `None`. The diagnostics go to
    /// standard error, as the command line's.
    #[cfg(any(target_arch = "wasm32", test))]
    pub(crate) fn run_in_memory(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
        query: &[u8],
        subject: &[u8],
    ) -> Result<Vec<u8>> {
        let report = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        self.run_io(
            settings,
            RunIo::Memory {
                query,
                subject,
                report: report.clone(),
            },
        )?;
        let bytes = std::mem::take(&mut *report.lock().expect("report mutex poisoned"));
        Ok(bytes)
    }

    fn run_io(
        self,
        settings: crate::blastinput::ncbi_environment::ApplicationSettings,
        io: RunIo<'_>,
    ) -> Result<()> {
        use crate::blastinput::app;
        use crate::blastinput::fasta_reader::{FastaInputSource, ReaderConfig};
        use crate::blastinput::input_files;
        use crate::blastinput::seq_range;
        let mut args = self;
        args.id_scope = args.parse_deflines.then(Default::default);
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
        // ```c
        //     if (GetExportSearchStrategyStream(args) ||
        //            m_FormattingArgs->ArchiveFormatRequested(args)) {
        //         locality = CBlastOptions::eBoth;
        //     }
        // ```
        // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
        // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
        // there too. `GetExportSearchStrategyStream` comes first: it creates the search strategy
        // file (`strategy::open_strategy_file`), which an error after it leaves empty.
        if args.export_search_strategy.is_some() {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:281-287,298-302
            // ```c++
            //         if (pssm.NotEmpty())
            //         {
            //             _ASSERT(queries.Empty());
            //             if(num_iters != 0)
            //             	retval.Reset(new blast::CExportStrategy(pssm, opts_hndl, search_db, client_id, num_iters));
            //             else
            //             	retval.Reset(new blast::CExportStrategy(pssm, opts_hndl, search_db, client_id));
            // ...
            //         if (pssm.NotEmpty())
            //         {
            //             NCBI_THROW(CInputException, eInvalidInput,
            //                        "Remote PSI-BL2SEQ is not supported");
            // ```
            // A PSSM search's request is the PSSM one (program `blastp`, service `psi`, the
            // whole `PssmWithParameters`; with `-subject` NCBI stops with that error after
            // the prolog): item F1b of line RG, with psiblast's.
            if args.in_pssm.is_some() {
                anyhow::bail!(
                        "-export_search_strategy with -in_pssm is not supported by LOSAT's TBLASTN yet (NCBI writes a PSI-BLAST request with -db and stops with \"Remote PSI-BL2SEQ is not supported\" with -subject)"
                    );
            }
        }
        let strategy_file = crate::report::archive::strategy::open_strategy_file(
            args.export_search_strategy.as_deref(),
            args.out.as_deref(),
            "TBLASTN",
        )?;
        let choice = app::parse_formatting_string(&args.outfmt)?;
        let format = app::report_format(&choice, app::Program::Tblastn, args.out.as_deref())?;
        // The BLAST archive (outfmt 11): its request, and the refusals of what LOSAT cannot
        // write, before anything is read or written.
        // The error archives follow the `-outfmt` text (`writes_error_archive`).
        // The run's options handle, which its archive and its search strategy share.
        let options = crate::report::archive::OptionsHandle::default();
        args.archive = super::archive::handle_for(
            &args,
            format == Some(app::ReportFormat::ArchiveFormat),
            &options,
        )?
        .map(|handle| handle.with_outfmt_text(&args.outfmt));
        // The search strategy's request, with the refusals of the archive's, and the archive's
        // options handle.
        args.strategy = strategy_file
            .map(|file| {
                crate::report::archive::strategy::StrategyExport::for_run(
                    file,
                    super::archive::option_inputs(&args),
                    choice.number,
                    crate::report::archive::strategy::HitlistArgs {
                        max_target_seqs: args.max_target_seqs,
                        num_descriptions: args.num_descriptions,
                        num_alignments: args.num_alignments,
                    },
                    args.database_name(),
                    "TBLASTN",
                    &options,
                )
            })
            .transpose()?;
        let archive = args.archive.clone();
        // LOSAT's thread capability is checked before any input or output, as for BLASTN.
        crate::utils::threading::validate_threads(args.num_threads)?;
        // The diagnostics stream: with the archive it also saves what it writes as the
        // messages of the next archive.
        let mut stderr = std::io::stderr();
        let mut archive_tee;
        let diagnostics: &mut (dyn Write + Send) = match &archive {
            Some(handle) => {
                archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
                &mut archive_tee
            }
            None => &mut stderr,
        };
        // The taxid or seqid list (line DI-8): NCBI's handler of the database arguments
        // (`s_GetTaxIDList`, `CSeqDBFileGiList`) comes before the other option handlers
        // (`m_BlastDbArgs` is among the first of `CBlastpAppArgs`' `m_Args`, blastp_args.cpp:61-65)
        // and before `CStdCmdLineArgs`, which opens the query and the output: its warning before
        // theirs (probe F04 of line DI-8), and a list error leaves no `-out`, only the error archive
        // that NCBI's `PrintErrorArchive` writes for outfmt 11 (probe E of line DI-8).
        let list = match args.database_name() {
            Some(_) => match args.db_filter.resolve(&settings, &mut *diagnostics) {
                Ok(list) => list,
                Err(error) => {
                    write_run_error_archive(&args, &io, archive.as_ref(), &error);
                    return Err(error);
                }
            },
            None => None,
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2553-2562
        // ```c
        //         CRef<blast::CBlastQueryVector> subjects;
        //         m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(),
        //                                        subj_range, parse_deflines,
        //                                        use_lcase_masks, subjects, m_IsMapper);
        //         m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects));
        //
        //     } else if (!m_IsIgBlast){
        //         // IgBlast permits use of germline database
        //         NCBI_THROW(CInputException, eInvalidInput,
        //            "Either a BLAST database or subject sequence(s) must be specified");
        //     }
        // ```
        // The handler of the database arguments reads the subjects (an empty subject set
        // fails there), before the query and the output are opened.
        // `-db`: the database is opened with the options (`search_cli`), as NCBI's
        // `InitializeSubject` opens it after `SetOptions`.
        let subjects_read = (|| -> Result<(
            Vec<FastaRecord>,
            seq_range::Placements,
            Option<Vec<crate::sequence::packed_subject::PackedNucleotide>>,
        )> {
            Ok(if args.db.is_some() {
            (Vec::new(), seq_range::Placements::default(), None)
        } else {
            let Some(subject_path) = args.subject.clone() else {
                return Err(app::missing_subject_error());
            };
            // Web ABI v1's subjects are the call's bytes: no file to name or open.
            let subject_input = match &io {
                RunIo::Files => {
                    input_files::check_utf8_file_name(&subject_path, "subject", "TBLASTN")?;
                    Input::File(input_files::open_input(&subject_path, "subject", "TBLASTN")?)
                }
                RunIo::Memory { subject, .. } => Input::Bytes(subject),
            };
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
            // ```c++
            //             subj_input_stream = &args[kArgSubject].AsInputFile();
            //         }
            //
            //         TSeqRange subj_range;
            //         if (args.Exist(kArgSubjectLocation) && args[kArgSubjectLocation]) {
            //             subj_range =
            //                 ParseSequenceRange(args[kArgSubjectLocation].AsString(),
            //                             "Invalid specification of subject location");
            //         }
            // ```
            // The subject range is read after the subject file is opened, before it is read.
            let subject_range = seq_range::parse_optional_range(
                args.subject_loc.as_deref(),
                seq_range::RangeRole::Subject,
                "TBLASTN",
            )?;
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:230-246
            // ```c++
            //     SDataLoaderConfig dlconfig(read_proteins);
            //     dlconfig.OptimizeForWholeLargeSequenceRetrieval();
            //
            //     CBlastInputSourceConfig iconfig(dlconfig);
            //     iconfig.SetRange(range);
            //     iconfig.SetBelieveDeflines(parse_deflines);
            //     iconfig.SetLowercaseMask(use_lcase_masking);
            //     iconfig.SetSubjectLocalIdMode();
            //     if (!read_proteins && gaps_to_Ns) {
            //         iconfig.SetConvertGapsToNs(true);
            //     }
            //
            //     CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
            //     CRef<CBlastInput> input(new CBlastInput(fasta));
            //     CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
            //     sequences = input->GetAllSeqs(*scope);
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2466-2470
            // ```c++
            //     EMoleculeType mol_type = Blast_SubjectIsNucleotide(opts.GetProgramType())
            //         ? CSearchDatabase::eBlastDbIsNucleotide
            //         : CSearchDatabase::eBlastDbIsProtein;
            //     m_IsProtein = (mol_type == CSearchDatabase::eBlastDbIsProtein);
            // ```
            // TBLASTN's subjects are nucleotide (`read_proteins` = `IsProtein()` is false, so the
            // reader assumes nucleotides, `fAssumeNuc`). NCBI reads the records one at a time, writes
            // each one's messages as it reads it, and checks each one's range after reading it: a
            // range that starts past the end of a record stops the reading there
            // (`read_subjects`). Records without residues (and intervals without letters) are
            // read without a message and reported when the formatter sets up the subjects
            // (`before_prolog`).
            let config = ReaderConfig::subject("TBLASTN", false, settings.data_loaders)
                .with_defline(archive.is_some() || args.strategy.is_some())
                .with_parse_deflines(args.parse_deflines);
            // EXPERIMENT (LOSAT_X_TBNSTORE): each record's letters packed as it is read
            // (`x_read_subject_stores`), with a subject range after the record is cut to its
            // interval (the record `cut_subjects` makes).
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
            // ```c
            //     virtual SBlastSequence GetCompressedPlusStrand() {
            //         SBlastSequence retval(size());
            //         string ncbi4na = kEmptyStr;
            //         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
            //         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
            // ```
            // Dispatch point: NCBI keeps the subjects in its scope and makes each one's search forms
            // from it. With the switch the records keep their deflines and the letters are held packed,
            // made again when a subject is searched or reported; without it the records keep them.
            // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs);
            // with LOSAT_X_TBNSTORE the packing of each record is in this interval too. No NCBI
            // counterpart.
            let t_read_fasta = crate::utils::xtiming::clock();
            // With the switch the records come back cut to the range (`Some` placements).
            let (read_subject_records, x_subject_stores, x_placements) = {
                let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
                // Both sides of the stores are recorded in x_subject_store.rs. The archive of
                // `-subject` writes every subject's whole letters (and its deflines) in every
                // batch, which the packed stores do not keep (they hold the interval of a
                // `-subject_loc` and no letters), so a run that writes it (or the search
                // strategy's requests, which hold the same Bioseqs) reads the records as the
                // switch's `=0` side does.
                let x_packed = super::x_subject_store::x_tbn_store()
                    && archive.is_none()
                    && args.strategy.is_none();
                match subject_input {
                    Input::File(file) => read_subject_input(
                        FastaInputSource::from_argument(&subject_path, file, config),
                        &args,
                        subject_range.as_ref(),
                        x_packed,
                        &mut warn,
                    )?,
                    Input::Bytes(bytes) => read_subject_input(
                        FastaInputSource::from_bytes(bytes, config),
                        &args,
                        subject_range.as_ref(),
                        x_packed,
                        &mut warn,
                    )?,
                }
            };
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
            let (subjects, subject_placements) = match x_placements {
                Some(placements) => (read_subject_records, placements),
                None => {
                    match seq_range::cut_subjects(&read_subject_records, subject_range.as_ref())? {
                        Some((cut, placements)) => (cut, placements),
                        None => (read_subject_records, seq_range::Placements::default()),
                    }
                }
            };
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
            // ```c
            // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
            //     : m_QueryVector(& queries)
            // {
            //     if (queries.Empty()) {
            //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
            //     }
            // ```
            // `m_Subjects.Reset(new blast::CObjMgr_QueryFactory(*subjects))` (blast_args.cpp:2557)
            // fails for a subject input without records.
            if subjects.is_empty() {
                return Err(app::empty_subjects_error());
            }
            (subjects, subject_placements, x_subject_stores)
            })
        })();
        let (subjects, subject_placements, x_subject_stores) = match subjects_read {
            Ok(read) => read,
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, &io, archive.as_ref(), &error);
                return Err(error);
            }
        };
        // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
        // ```c
        //     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue() &&
        //         m_InputStream == NULL) {
        // ...
        //         else {
        //             m_InputStream = &args[kArgQuery].AsInputFile();
        //         }
        //     }
        // ...
        //     else {
        //         m_OutputStream = &args[kArgOutput].AsOutputFile();
        //     }
        // ```
        // The query is opened and the output file created before the options are checked;
        // `-` is standard input and standard output. Web ABI v1's query is the call's bytes.
        let query_input = match &io {
            // `-in_pssm`: the query is the PSSM, which the options' handler reads
            // (`check_options`); NCBI opens the default `-query -` and never reads it.
            RunIo::Files if args.in_pssm.is_some() => Input::Bytes(&[]),
            RunIo::Files => {
                input_files::check_utf8_file_name(&args.query, "query", "TBLASTN")?;
                // NCBI opens the query file inside `Run` too: its error is saved and an error
                // archive follows.
                match input_files::open_input(&args.query, "query", "TBLASTN") {
                    Ok(file) => Input::File(file),
                    Err(error) => {
                        write_run_error_archive(&args, &io, archive.as_ref(), &error);
                        return Err(error);
                    }
                }
            }
            RunIo::Memory { query, .. } => Input::Bytes(query),
        };
        let mut stream = report_stream(&args, &io)?;
        let pairwise = format == Some(app::ReportFormat::Pairwise);
        // The format as the reports read it: the number and the custom fields, or the value as
        // given when it has a delimiter, which `normalized` drops and the custom tabular report
        // reads (`tabular::tabular_spec`, which parses it again as NCBI does).
        let outfmt = if choice.delimiter.is_empty() {
            choice.normalized()
        } else {
            args.outfmt.clone()
        };
        let result = {
            let mut outputs =
                ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut *diagnostics);
            search_cli(
                args,
                list,
                query_input,
                &subjects,
                &subject_placements,
                x_subject_stores.as_deref(),
                &mut outputs,
                &choice,
                settings,
            )
        };
        // NCBI writes an error archive after the archives of the batches before when the
        // application stops with one of its own errors, and when the run ends with messages that
        // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status).
        if let Some(handle) = &archive {
            match &result {
                Err(error) => {
                    if let Some(lines) = super::archive::error_lines(error, false) {
                        let _ = write_error_archive(&io, handle, error, &lines, &mut stream);
                    }
                }
                Ok(()) => {
                    if handle.has_messages() {
                        let _ =
                            crate::report::archive::write_error_archive(handle, &[], &mut stream);
                    }
                }
            }
        }
        let flushed = std::io::Write::flush(&mut stream);
        // NCBI reference: ncbi-blast/c++/src/app/blast/blast_app_util.hpp:252-255
        // ```c
        //     catch (const std::ios::failure&) {                                      \
        //         LOG_POST(Error << "BLAST failed to write output");                  \
        //         exit_code = BLAST_OUTPUT_ERROR;                                     \
        //     }                                                                       \
        // ```
        // The outfmt 0 formatter's stream throws when a write fails; for outfmt 6/7 NCBI
        // aborts and LOSAT reports the error (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES exception 3).
        if stream.failed && pairwise {
            return Err(crate::cli::NativeError {
                exit: 6,
                message: "BLAST failed to write output\n".to_string(),
            }
            .into());
        }
        result?;
        flushed.context("failed to write the output")
    }
}

/// Where `run` reads its inputs and writes its report.
enum RunIo<'a> {
    /// The command line: the files of `-query` and `-subject` (`-`: standard input), and the
    /// `-out` file or standard output.
    Files,
    /// Web ABI v1 (`run_in_memory`): the call's bytes, and the report kept in memory.
    #[cfg_attr(
        not(any(target_arch = "wasm32", test)),
        expect(
            dead_code,
            reason = "made only by `run_in_memory`, for the web ABI (wasm32) and tests"
        )
    )]
    Memory {
        query: &'a [u8],
        subject: &'a [u8],
        report: std::sync::Arc<std::sync::Mutex<Vec<u8>>>,
    },
}

/// An input of `run`: the opened file of `-query` or `-subject`, or web ABI v1's bytes.
enum Input<'a> {
    File(std::fs::File),
    Bytes(&'a [u8]),
}

/// The reader of `search_cli`'s queries over an `Input`.
enum QuerySource<'a> {
    File(crate::blastinput::fasta_reader::FastaInputSource<std::fs::File>),
    Bytes(crate::blastinput::fasta_reader::FastaInputSource<&'a [u8]>),
}

/// The report kept in memory (`RunIo::Memory`) as the report's stream.
struct MemoryReport(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);

impl Write for MemoryReport {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0
            .lock()
            .expect("report mutex poisoned")
            .extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// The subjects of `source` as `run` reads them (`read_subjects`), or with LOSAT_X_TBNSTORE's
/// packed stores (`packed`, see `run`), with the reader's messages to `warn`.
fn read_subject_input<R: std::io::Read>(
    mut source: crate::blastinput::fasta_reader::FastaInputSource<R>,
    args: &TblastnArgs,
    range: Option<&crate::blastinput::seq_range::SequenceRange>,
    packed: bool,
    warn: &mut dyn FnMut(&[u8]) -> std::io::Result<()>,
) -> Result<(
    Vec<FastaRecord>,
    Option<Vec<crate::sequence::packed_subject::PackedNucleotide>>,
    Option<crate::blastinput::seq_range::Placements>,
)> {
    source.set_id_scope(args.id_scope.clone());
    if packed {
        let (records, stores, placements) = super::x_subject_store::x_read_subject_stores(
            &mut source,
            range,
            warn,
            args.num_threads,
        )?;
        Ok((records, stores, Some(placements)))
    } else {
        Ok((
            crate::blastinput::fasta_reader::read_subjects(&mut source, range, warn)?,
            None,
            None,
        ))
    }
}

/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &TblastnArgs,
    io: &RunIo<'_>,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    // `-outfmt "11 qseqid"`: NCBI's `PrintErrorArchive` stops at the number before it opens `-out`.
    if !handle.error_archive_enabled() {
        return;
    }
    let Some(lines) = super::archive::error_lines(error, true) else {
        return;
    };
    if let Ok(mut report) = report_stream(args, io) {
        let _ = write_error_archive(io, handle, error, &lines, &mut report);
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The error archive of a run that stopped with `error` (`lines`), into `report`: on the command
/// line after the error's text (`write_error_archive_after_post`); web ABI v1 (`RunIo::Memory`)
/// reports the error itself, so there only the archive.
fn write_error_archive(
    io: &RunIo<'_>,
    handle: &crate::report::archive::ArchiveHandle,
    error: &anyhow::Error,
    lines: &[crate::report::diag_log::Message],
    report: &mut crate::cli::ReportStream,
) -> Result<()> {
    match io {
        RunIo::Files => {
            crate::report::archive::write_error_archive_after_post(handle, error, lines, report)
        }
        RunIo::Memory { .. } => crate::report::archive::write_error_archive(handle, lines, report),
    }
}

/// The report stream of a run: the `-out` file (created) or standard output, or the report
/// in memory (`RunIo::Memory`).
fn report_stream(args: &TblastnArgs, io: &RunIo<'_>) -> Result<crate::cli::ReportStream> {
    Ok(match io {
        RunIo::Files => cli_report_stream(open_out_file(args)?),
        RunIo::Memory { report, .. } => crate::cli::ReportStream {
            inner: Box::new(MemoryReport(report.clone())),
            failed: false,
        },
    })
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &TblastnArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    if let Some(path) = args.out.as_deref() {
        crate::blastinput::input_files::check_utf8_file_name(path, "out", "TBLASTN")?;
    }
    Ok(
        match args.out.as_deref().filter(|path| path.as_os_str() != "-") {
            Some(path) => Some(std::io::BufWriter::new(
                std::fs::File::create(path).map_err(|_| crate::cli::inaccessible("out", path))?,
            )),
            None => None,
        },
    )
}

/// The report stream of the command line: the `-out` file, or standard output.
fn cli_report_stream(
    out_file: Option<std::io::BufWriter<std::fs::File>>,
) -> crate::cli::ReportStream {
    crate::cli::ReportStream {
        inner: match out_file {
            Some(file) => Box::new(file) as Box<dyn std::io::Write + Send>,
            None => crate::cli::report_standard_output(),
        },
        failed: false,
    }
}

/// What `search_cli` and `run_local_database` do after the options are processed and before the
/// query: the opening of the database of `-db` (`list`: its taxid or seqid list, `None` for a
/// search without `-db`), and what the requested reports ask of it.
fn open_search_database(
    args: &mut TblastnArgs,
    list: Option<&crate::seqdb::DbList>,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<Option<TblastnDatabase>> {
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
    // ```c++
    //         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
    //         CRef<CLocalDbAdapter> db_adapter;
    //         CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    //         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
    //                          db_adapter, scope);
    // ```
    // The database of `-db` is opened after the options are processed, before
    // `Query is Empty!` (as BLASTN's).
    let database = match args.database_name() {
        Some(name) => {
            let database = super::database::TblastnDatabase::open_restricted(
                name,
                settings,
                list,
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis);
            crate::algorithm::common::database::note_archive_database_filter(
                args.archive.as_ref(),
                args.strategy.as_ref(),
                &args.db_filter,
                database.list_taxids(),
                database.statistics_length_is_empty(),
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "TBLASTN",
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1093-1104
    // ```c++
    // void CBlastTabularInfo::x_CheckTaxDB()
    // {
    //    	if( x_IsFieldRequested(eSubjectSciNames) ||
    //     ...
    //    		string resolved = SeqDB_ResolveDbPath("taxdb.bti");
    // ```
    // A tabular field list that shows a name of the taxonomy database asks whether `taxdb.bti`
    // is on the search path: its formatters warn without it (the plural names), and LOSAT
    // refuses a name of a subject with a taxid with it (`tabular`), as BLASTN.
    if outputs
        .formats
        .iter()
        .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
        .any(|spec| spec.needs().taxonomy_names)
    {
        args.tabular_taxdb = Some(match database.as_ref() {
            Some(database) => database.taxdb_found(),
            None => crate::algorithm::common::database::subject_taxdb_found(settings, "TBLASTN")?,
        });
    }
    // LOSAT's refusals of a database for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for format in &outputs.formats {
            let choice = crate::blastinput::app::parse_formatting_string(format.outfmt)?;
            if let Some(serial) = crate::blastinput::app::ReportFormat::from_number(choice.number)
                .and_then(crate::report::serial::SerialFormat::of)
            {
                crate::report::serial::seqalign::check_database(
                    serial,
                    database.name(),
                    &database.pairwise_database().definition,
                    crate::blastinput::app::Program::Tblastn,
                )?;
            }
        }
    }
    Ok(database)
}

/// The part of `run` after the files are opened: NCBI's processing of the options,
/// `Query is Empty!`, the query batch size and the subjects without letters, the reading of
/// the queries and the search.
#[allow(clippy::too_many_arguments)]
fn search_cli(
    mut args: TblastnArgs,
    list: Option<crate::seqdb::DbList>,
    query_input: Input<'_>,
    subjects: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    outputs: &mut ReportOutputs<'_>,
    choice: &crate::blastinput::app::FormatChoice,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    use crate::blastinput::fasta_reader::{
        read_queries, FastaInputSource, QueryEnd, QueryRecords, ReaderConfig,
    };
    let resolved = args.check_options(outputs.diagnostics)?;
    let database = open_search_database(&mut args, list.as_ref(), &settings, outputs)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:194-196,225-238
    // ```c++
    //         /*** Get the query sequence(s) or PSSM (these two options are mutually
    //          * exclusive) ***/
    //         CRef<CPssmWithParameters> pssm = m_CmdLineArgs->GetInputPssm();
    // ...
    //         } else {
    //             _ASSERT(pssm->HasQuery());
    //             _ASSERT(pssm->GetQuery().IsSeq());
    //             query.Reset(new CBlastQueryVector());
    //             scope->AddTopLevelSeqEntry(pssm->SetQuery());
    //
    //             CRef<CSeq_loc> sl(new CSeq_loc());
    //             sl->SetWhole().Assign(*pssm->GetQuery().GetSeq().GetFirstId());
    //
    //             CRef<CBlastSearchQuery> q(new CBlastSearchQuery(*sl, *scope));
    //             query->AddQuery(q);
    // ```
    // `-in_pssm` (PSI-TBLASTN): the PSSM's Bioseq is the one query, read with the options; no
    // FASTA input, no `Query is Empty!` and no query batches (one search of the PSSM,
    // tblastn_app.cpp:305-327).
    if let Some(psi) = resolved.pssm.clone() {
        super::in_pssm::check_limits(&args, &resolved, &psi, choice.number)?;
        drop(query_input);
        crate::blastinput::input_files::write_empty_subject_warnings_by_length(
            subjects,
            |index| {
                x_subject_stores.map_or(subjects[index].seq().len(), |stores| stores[index].len())
            },
            "tblastn",
            outputs.diagnostics,
        )?;
        crate::blastinput::app::xinclude_check_for(choice, args.out.as_deref())?;
        // The PSSM's Bioseq as the reports show it (`PsiSearch::query_record`); with
        // `-parse_deflines` its ids join the scope of the `-subject` records as a FASTA query's do
        // (`IdScope::add_record`; NCBI BLAST+ 2.17.0 fails on a subject with the query's id), and
        // the database's sequences are checked as for a FASTA query (`QueryScopes`).
        let records =
            vec![psi.query_record_in(args.parse_deflines, args.id_scope.as_deref(), "TBLASTN")?];
        return search(
            &args,
            resolved,
            1,
            &records,
            QueryReading::End,
            subjects,
            subject_placements,
            x_subject_stores,
            database.as_ref(),
            outputs,
        );
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:200-221
    // ```c++
    //         SDataLoaderConfig dlconfig =
    //             InitializeQueryDataLoaderConfiguration(query_opts->QueryIsProtein(),
    //                                                    db_adapter);
    // ...
    //         if (pssm.Empty()) {
    //             CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
    //                                          query_opts->UseLowercaseMasks(),
    //                                          query_opts->GetParseDeflines(),
    //                                          query_opts->GetRange());
    //             if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
    //                	ERR_POST(Warning << "Query is Empty!");
    //                	return BLAST_EXIT_SUCCESS;
    //             }
    //             fasta.Reset(new CBlastFastaInputSource(
    //                                          m_CmdLineArgs->GetInputStream(),
    //                                          iconfig));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:51,106
    // ```c++
    //     const bool kQueryIsProtein = true;
    // ...
    //     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:328-330
    // ```c++
    //     flags += (m_ReadProteins
    //               ? CFastaReader::fAssumeProt
    //               : CFastaReader::fAssumeNuc);
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // ```c++
    // 	char c;
    // 	CNcbiStreampos orig_p = in.tellg();
    // 	// Piped input
    // 	if(orig_p < 0)
    // 		return false;
    // ```
    // The query of tblastn is a protein (`QueryIsProtein()`, `fAssumeProt`) and its data
    // loaders are those of the subjects (`SDataLoaderConfig`). The position is taken on the
    // opened file, before it is read (`stream_is_empty`). When the subjects were read from
    // standard input too, `cin` has reached its end (a failed stream), so NCBI gets no
    // position for the query either, and its reader is at its end: no query batch is read.
    let config = ReaderConfig::query("TBLASTN", true, settings.data_loaders)
        .with_defline(args.archive.is_some() || args.strategy.is_some())
        .with_parse_deflines(args.parse_deflines);
    let (mut query_source, query_empty) = match query_input {
        Input::File(query_file) => {
            let shared_standard_input = args.query.as_os_str() == "-"
                && args
                    .subject
                    .as_deref()
                    .is_some_and(|path| path.as_os_str() == "-");
            let mut source = FastaInputSource::from_argument(&args.query, query_file, config);
            source.set_id_scope(args.id_scope.clone());
            let empty = !shared_standard_input && source.stream_is_empty();
            (QuerySource::File(source), empty)
        }
        // Web ABI v1's bytes, read as the regular file with the same bytes: `stream_is_empty`
        // finds such a file empty when it holds white space only.
        Input::Bytes(bytes) => {
            let mut source = FastaInputSource::from_bytes(bytes, config);
            source.set_id_scope(args.id_scope.clone());
            let empty = bytes
                .iter()
                .all(|&byte| crate::blastinput::fasta_reader::stream::input_space(byte));
            (QuerySource::Bytes(source), empty)
        }
    };
    if query_empty {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    let batch_size = before_prolog(subjects, x_subject_stores, outputs.diagnostics)?;
    // NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:267-270
    // ```c
    // 	if(UseXInclude(*fmt_args, args[kArgOutput].AsString())) {
    //         	formatter.SetBaseFile(args[kArgOutput].AsString());
    //         }
    //         formatter.PrintProlog();
    // ```
    // The XInclude formats written to standard output fail after the formatter is made.
    crate::blastinput::app::xinclude_check_for(choice, args.out.as_deref())?;
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-277
    // ```c++
    //             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //                 query =  input->GetNextSeqBatch(*scope);
    // ```
    // The records are read here, and their messages are written batch by batch
    // (`read_query_batches`). What ends the reading (the end of the input, an `eEOF`, a
    // reader error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads
    // it in, after the reports of the batches before.
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let QueryRecords { records, end } = match &mut query_source {
        QuerySource::File(source) => read_queries(source),
        QuerySource::Bytes(source) => read_queries(source),
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(
        &args,
        resolved,
        batch_size,
        &records,
        reading,
        subjects,
        subject_placements,
        x_subject_stores,
        database.as_ref(),
        outputs,
    )
}

/// NCBI's steps of tblastn between `Query is Empty!` and the outfmt 0 prolog that concern
/// the subjects: the query batch size, then the formatter, which sets up the subjects and
/// warns about each one without letters (a record without residues, or an interval that
/// starts just past its record's end). Such a subject stays in the database statistics
/// with no letters and is never searched (`x_scan_plan` skips subjects
/// shorter than a word). Returns the query batch size.
///
/// NCBI reference: ncbi-blast/c++/src/app/blast/tblastn_app.cpp:217-221
/// ```c
///             fasta.Reset(new CBlastFastaInputSource(
///                                          m_CmdLineArgs->GetInputStream(),
///                                          iconfig));
///             input.Reset(new CBlastInput(&*fasta,
///                                         m_CmdLineArgs->GetQueryBatchSize()));
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:129-136
/// ```c
///     if(m_IsDbScan) {
/// 	int num_seqs=0;
///         int total_length=0;
/// 	if (!is_remote_search)
///         {
///                 BlastSeqSrc* seqsrc = db_adapter.MakeSeqSrc();
///                 num_seqs=BlastSeqSrcGetNumSeqs(seqsrc);
///                 total_length=static_cast<int>(BlastSeqSrcGetTotLen(seqsrc));
/// ```
/// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
/// ```c
///         catch(CBlastException & e ) {
///         	// Skip bad subject sequence
///         	if(e.GetErrCode() == CBlastException::eInvalidArgument) {
///         		seqblk_vec->push_back(subj);
///         ...
///         		warning += "Subject sequence contains no data";
///         		ERR_POST(Warning << warning);
///         		continue;
/// ```
/// NCBI's gapped tblastn with the Gumbel block never stops for the average subject length,
/// also when every subject is without letters: its shortest subject is raised to
/// `BLAST_SEQSRC_MINLENGTH` before the division by 3 (`run_local_search_with_pool` ports
/// blast_setup.c:969-980 and seqsrc_multiseq.cpp:226-240). Without it (`-searchsp` or
/// `OLD_FSC`; LOSAT rejects `-ungapped`) the average is used, and below one letter the batch
/// fails (`stage_d_pipeline::SetupEngineError::SUBJECT_LENGTH_INVALID`).
fn before_prolog(
    subjects: &[FastaRecord],
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    diagnostics: &mut dyn Write,
) -> Result<u32> {
    let batch_size = crate::blastinput::app::query_batch_size("TBLASTN", 20000)?;
    // EXPERIMENT (LOSAT_X_TBNSTORE): a packed subject's length is its store's.
    crate::blastinput::input_files::write_empty_subject_warnings_by_length(
        subjects,
        |index| x_subject_stores.map_or(subjects[index].seq().len(), |stores| stores[index].len()),
        "tblastn",
        diagnostics,
    )?;
    Ok(batch_size)
}

// NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
// results = lcl_blast.Run();
// formatter.PrintOneResultSet(**result, query);
// NCBI c++/src/app/blast/blast_formatter.cpp:429-467:
// CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
// formatter.PrintProlog();
// ...
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one TBLASTN search over already parsed records and writes every requested
/// output format from the same result (the entry of web ABI v2; the command line reaches
/// `search` through `search_cli`).
///
/// Each requested `-outfmt` is read as on the command line, and NCBI's checks of the
/// options and of the subjects come before the search. The `-query` and `-subject` values
/// of `args` are used only as display names.
///
/// The records are those of NCBI's reader (`fasta_reader`): the subjects as NCBI's
/// nucleotide reader reads them, the queries as its protein reader reads them. Queries
/// without a record are NCBI's `Query is Empty!`.
pub fn run_local(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::TBLASTN_TASK.check(&mut args.task)?;
    super::in_pssm::check_records_only(&args)?;
    // The host has no database files; LOSAT Web searches the records it gives.
    if args.db.is_some() {
        anyhow::bail!("-db (a database search) is not supported in the browser yet");
    }
    // The BLAST archive (outfmt 11) among the formats: its request, the refusals of what LOSAT
    // cannot write, and a diagnostics stream that saves the messages. The records' deflines are
    // those the caller's reader kept; records without one get `>` and their title.
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| super::stage_e_report::format_number(format.outfmt) == 11);
    if archive_requested {
        args.archive = super::archive::handle_for(
            &args,
            true,
            &crate::report::archive::OptionsHandle::default(),
        )?
        .map(|handle| outputs.with_archive_outfmt(handle));
        if let Some(handle) = args.archive.clone() {
            let result = outputs.with_archive_diagnostics(handle.clone(), |view| {
                run_local_formats(args, query_records, subject_records, view)
            });
            return outputs.finish_archive(&handle, query_records.is_empty(), result);
        }
    }
    run_local_formats(args, query_records, subject_records, outputs)
}

/// The numbers of the requested formats, each read as on the command line. `out`: the file of
/// `api::run_to_file` (a format with the caller's receiver of side files takes that receiver's,
/// `FormatOutput::side_out`).
fn read_formats(outputs: &ReportOutputs<'_>, out: Option<&std::path::Path>) -> Result<Vec<i32>> {
    use crate::blastinput::app;
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        if app::report_format(&choice, app::Program::Tblastn, format.side_out(out))?.is_none() {
            app::formatting_handler_check(&choice, app::Program::Tblastn)?;
            app::xinclude_check(&choice)?;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2887
        // ```c++
        //     s_ValidateCustomDelim(m_CustomOutputFormatSpec,m_CustomDelim);
        // ```
        // The formatting handler checks the custom delimiter of each format (the command line's
        // in `check_options`).
        app::formatting_handler_check(&choice, app::Program::Tblastn)?;
        outfmts.push(choice.number);
    }
    Ok(outfmts)
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:189-193
// ```c++
//         CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
//         CRef<CLocalDbAdapter> db_adapter;
//         CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
//         InitializeSubject(db_args, m_OptsHndl, m_CmdLineArgs->ExecuteRemotely(),
//                          db_adapter, scope);
// ```
/// Runs one TBLASTN search of the records `query_records` against the database of `args` and
/// writes every requested output format from the same result (the library API's entry; the
/// command line's is `run`). `settings` are the NCBI application settings, which find the
/// database. The steps are `run`'s and `search_cli`'s, in their order; what differs is what they
/// read: the records are read already, and the outputs are the caller's. An archive (outfmt 11)
/// ends as the command line's does (`ReportOutputs::finish_archive`).
pub fn run_local_database(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::TBLASTN_TASK.check(&mut args.task)?;
    super::in_pssm::check_records_only(&args)?;
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| super::stage_e_report::format_number(format.outfmt) == 11);
    args.archive = super::archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::OptionsHandle::default(),
    )?
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads)?;
    match args.archive.clone() {
        Some(handle) => {
            let result = outputs.with_archive_diagnostics(handle.clone(), |view| {
                search_database(args, query_records, settings, view)
            });
            outputs.finish_archive(&handle, query_records.is_empty(), result)
        }
        None => search_database(args, query_records, settings, outputs),
    }
}

/// `run_local_database` from the taxid or seqid list on.
fn search_database(
    mut args: TblastnArgs,
    query_records: &[FastaRecord],
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The handler of the database arguments reads the list before the other handlers (`run`).
    let list = args
        .db_filter
        .resolve(&settings, &mut *outputs.diagnostics)?;
    let outfmts = read_formats(outputs, args.out.as_deref())?;
    let resolved = args.check_options_for_formats(&outfmts, outputs.diagnostics)?;
    let database = open_search_database(&mut args, list.as_ref(), &settings, outputs)?
        .ok_or_else(crate::blastinput::app::missing_subject_error)?;
    if query_records.is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    let batch_size = before_prolog(&[], None, outputs.diagnostics)?;
    search(
        &args,
        resolved,
        batch_size,
        query_records,
        QueryReading::Records,
        &[],
        &crate::blastinput::seq_range::Placements::default(),
        None,
        Some(&database),
        outputs,
    )
}

fn run_local_formats(
    args: TblastnArgs,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    let outfmts = read_formats(outputs, args.out.as_deref())?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        seq_range::RangeRole::Subject,
        "TBLASTN",
    )?;
    for record in
        &subject_records[..seq_range::subjects_read(subject_records, subject_range.as_ref())]
    {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = seq_range::cut_subjects(subject_records, subject_range.as_ref())?;
    if subject_records.is_empty() {
        return Err(app::empty_subjects_error());
    }
    let resolved = args.check_options_for_formats(&outfmts, outputs.diagnostics)?;
    if query_records.is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [tblastn] Query is Empty!\n")?;
        return Ok(());
    }
    let whole_subjects = seq_range::Placements::default();
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((cut, placements)) => (cut.as_slice(), placements),
        None => (subject_records, &whole_subjects),
    };
    // The caller owns the records: their letters are read in place.
    let batch_size = before_prolog(searched_subjects, None, outputs.diagnostics)?;
    search(
        &args,
        resolved,
        batch_size,
        query_records,
        QueryReading::Records,
        searched_subjects,
        subject_placements,
        None,
        None,
        outputs,
    )
}

/// The outfmt 0 prologs of a run that fails in its first query batch: NCBI writes the prolog
/// (`formatter.PrintProlog()`) before it reads the first batch. `subject_path` names a
/// FASTA subject set of `num_subjects` records (`subject_letters`); a database search shows
/// its `Database:` lines (`database`).
fn write_first_batch_prologs(
    outputs: &mut ReportOutputs<'_>,
    subject_path: Option<&std::path::Path>,
    num_subjects: usize,
    subject_letters: SubjectLetters<'_>,
    database: Option<&crate::report::pairwise::PairwiseDatabase>,
    taxdb_found: bool,
    // The program line of outfmt 0 (`TBLASTN`; `PSITBLASTN` with `-in_pssm`).
    program: &str,
) -> Result<()> {
    let name = subject_path.map_or_else(String::new, |path| {
        format!("User specified sequence set (Input: {})", path.display())
    });
    let total = match database {
        Some(_) => 0,
        None => subject_letters.lengths().sum(),
    };
    for format in outputs.formats.iter_mut() {
        if crate::blastinput::app::parse_formatting_string(format.outfmt)?.number == 0 {
            let mut writer = format.sink.open()?;
            crate::report::pairwise::write_tblastn_pairwise_prolog(
                &mut writer,
                program,
                crate::NCBI_BLAST_VERSION,
                &name,
                num_subjects,
                total,
                database,
            )?;
            writer.flush()?;
        } else if let Some(kind) = super::stage_e_report::xml2_kind(format.outfmt) {
            // The XML2 and JSON documents start with their prolog too.
            let mut writer = format.sink.open()?;
            writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
            writer.flush()?;
        } else if let Some(spec) = super::tabular::tabular_spec(format.outfmt) {
            // Outfmt 20's line of the field names and its formatter's warning (nothing for the
            // other tabular kinds).
            let mut writer = format.sink.open()?;
            super::tabular::write_tabular_prolog(
                &spec,
                taxdb_found,
                &mut writer,
                &mut *outputs.diagnostics,
            )?;
            writer.flush()?;
        }
    }
    Ok(())
}

/// The search of `run_local` and of the CLI, after NCBI's checks of the options, the
/// reading of the subjects and the steps before the prolog (`before_prolog`): LOSAT's
/// limits and the environment, then the query batches (`read_query_batches`), the search of
/// the batches that NCBI searches and their reports.
#[allow(clippy::too_many_arguments)]
fn search(
    args: &TblastnArgs,
    resolved: ResolvedTblastnArgs,
    batch_size: u32,
    query_records: &[FastaRecord],
    reading: QueryReading,
    subject_records: &[FastaRecord],
    subject_placements: &crate::blastinput::seq_range::Placements,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the subjects' letters, packed (the records keep none).
    x_subject_stores: Option<&[crate::sequence::packed_subject::PackedNucleotide]>,
    // `-db`: the database (the subject records are then none).
    database: Option<&TblastnDatabase>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    use crate::blastinput::app;
    use crate::blastinput::seq_range;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // ```c++
    //         formatter.PrintProlog();
    //         ...
    //         /*** Process the input ***/
    //         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //         ...
    //             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
    // ```
    // What the run writes before its query batches ends here.
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.queries_begin();
    }
    // The records are NCBI's reader's records: the subjects of its nucleotide reader (a
    // subject without letters stays in the database, its warning written before the prolog,
    // `before_prolog`) and the queries of its protein reader. The query range applies to
    // every query record as NCBI's batch reader reads it (`seq_range::cut_queries`): a record
    // is searched cut, and one whose interval starts just past its end is a query without
    // data (`Sequence contains no data`, as a record without residues).
    let input_records = query_records;
    let (searched_queries, query_input) = match resolved.query_range.as_ref() {
        Some(range) => {
            let ranged = seq_range::cut_queries(query_records, range);
            (std::borrow::Cow::Owned(ranged.records), ranged.input)
        }
        None => (
            std::borrow::Cow::Borrowed(query_records),
            seq_range::QueryInput::whole(query_records),
        ),
    };
    let SearchSettings {
        compo_adjust_mode,
        identity_cbs_reset,
        scoring,
    } = search_settings(&resolved)?;
    // `-in_pssm` (PSI-TBLASTN): the query's PSSM, which the batch's search and reports read.
    let pssm = resolved.pssm.clone();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_memento_priv.hpp:64-74
    // ```c++
    //     CBlastOptionsMemento(CBlastOptionsLocal* local_opts)
    //     {
    //         m_ProgramType = local_opts->GetProgramType();
    //         m_QueryOpts = local_opts->m_QueryOpts;
    // ...
    //         m_ExtnOpts = local_opts->m_ExtnOpts;
    // ```
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:288
    // ```c++
    //                 CLocalBlast lcl_blast(query_factory, opt_handle, db_adapter);
    // ```
    // IDENTITY's reset (`search_settings`) writes through the memento's pointers into the options
    // that every query batch reuses: the first batch resets the mode and warns its queries, and
    // the later batches start in mode 0 without the warning. These are the searched queries of
    // the first batch.
    let identity_warned = if identity_cbs_reset {
        query_input
            .batches(batch_size)
            .first()
            .map_or(0, |batch| batch.searched.end)
    } else {
        0
    };
    app::check_search_environment("TBLASTN", database.is_none(), outputs.only_archive())?;
    app::check_query_split_environment("TBLASTN", false, true)??;
    let subjects = subject_records;
    // The search reads the residues of the subject records in place (NCBI's subjects stay in the
    // reader's records, `CBlastQueryVector`); they are not copied for it. EXPERIMENT
    // (LOSAT_X_TBNSTORE): or from their packed form.
    let subject_seqs: Vec<&[u8]> = match x_subject_stores {
        Some(_) => Vec::new(),
        None => subjects.iter().map(FastaRecord::seq).collect(),
    };
    let subject_letters = match x_subject_stores {
        Some(stores) => SubjectLetters::Packed(stores),
        None => SubjectLetters::Records(&subject_seqs),
    };
    // What the requested formats need beyond outfmt 0, 6 and 7, decided once for the search
    // (`ReportNeeds`): with the Seq-aligns (outfmt 8, 12, 11), the engine positives of the
    // composition-adjusted matches (`LocalStageDProfile::engine_positives`).
    // The archive's Seq-aligns are those of outfmt 8.
    let engine_positives =
        crate::report::record::ReportNeeds::of(outputs.formats.iter().filter_map(|format| {
            app::ReportFormat::from_number(super::stage_e_report::format_number(format.outfmt))
        }))
        .seqalign
            || args.archive.is_some();
    // A database's numbers and `-dbsize` (`SubjectSetNumbers`), and the `Database:` lines of
    // a database search.
    let subject_set = SubjectSetNumbers {
        database: database.map(TblastnDatabase::numbers),
        dbsize: args.dbsize,
        searchsp: args.searchsp,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:922-923
        // ```c
        //     use_old_fsc = getenv("OLD_FSC");
        //     if (!use_old_fsc) sbp->gbp = s_BlastGumbelBlkNew();
        // ```
        // `-searchsp` sets OLD_FSC (blast_args.cpp:303-307).
        old_fsc: args.searchsp.is_some() || crate::blastinput::app::old_fsc_environment(),
    };
    // Without the Gumbel block NCBI's shortest subject length is the statistics' average
    // (blast_setup.c:975), which divides by their number of sequences.
    if let Some(database) = database.filter(|_| subject_set.old_fsc) {
        database.check_statistics_count()?;
    }
    let pairwise_database = database.map(TblastnDatabase::pairwise_database);
    // BLAST XML (5, and the XML2 documents of 16 and 15) cannot carry a database name outside
    // US-ASCII (NCBI's stream throws): refused before the search, and before anything is
    // written.
    if let Some(database) = database {
        if outputs
            .formats
            .iter()
            .any(|format| super::stage_e_report::is_structured_xml(format.outfmt))
        {
            crate::report::blast_xml::names::check_database_name("TBLASTN", database.name())?;
        }
    }
    let subject_path = match database {
        Some(_) => None,
        None => Some(
            args.subject
                .as_ref()
                .context("TBLASTN subject file missing")?,
        ),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-278
    // ```c++
    //             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
    //
    //                 query =  input->GetNextSeqBatch(*scope);
    //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
    // ```
    // NCBI reads the query batches one at a time after the outfmt 0 prolog. The run stops at
    // the first batch that NCBI fails (a reader error, a batch without a query, a batch of
    // queries without letters), after the reports of the batches before (without the
    // epilog) and the messages of that batch's records (`read_query_batches`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:59-68
    // ```c++
    // CPsiBlastImpl::CPsiBlastImpl(CRef<objects::CPssmWithParameters> pssm,
    //                              CRef<CLocalDbAdapter> subject,
    //                              CConstRef<CPSIBlastOptionsHandle> options)
    // : m_Pssm(pssm), m_Query(0), m_Subject(subject), m_OptsHandle(options),
    //     m_ResultType(eDatabaseSearch)
    // {
    //     x_Validate();
    //     x_ExtractQueryFromPssm();
    //     x_CreatePssmScoresFromFrequencyRatios();
    // }
    // ```
    // `-in_pssm`: NCBI's checks of the PSSM (`CPsiBlastValidate::Pssm`) and the PSSM engine's
    // scores of a PSSM of frequency ratios only (`PsiBlastComputePssmScores`) fail when the search
    // object is made, after the outfmt 0 prolog (tblastn_app.cpp:270, 322) and before any
    // subject is read.
    if let Some(error) = pssm
        .as_deref()
        .and_then(|psi| psi.validation_error.as_deref())
    {
        write_first_batch_prologs(
            outputs,
            subject_path.map(|path| path.as_path()),
            subjects.len(),
            subject_letters,
            pairwise_database.as_ref(),
            args.tabular_taxdb.unwrap_or(false),
            super::in_pssm::PROGRAM,
        )?;
        return Err(app::engine_error(error));
    }
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`); with
    // `-in_pssm`, of the PSSM's query when the reports show its ids (`PsiSearch::query_record`).
    let believed_queries = crate::blastinput::query_batch::believed_queries(
        args.parse_deflines,
        resolved.pssm.is_some(),
        input_records,
    );
    let mut batches = crate::blastinput::query_batch::read_query_batches(
        input_records,
        &searched_queries,
        &query_input,
        reading,
        batch_size,
        crate::blastinput::query_batch::QueryScopes::of_search(
            believed_queries,
            database.is_none(),
            believed_queries
                && crate::blastinput::query_batch::outfmts_have_xml(
                    outputs.formats.iter().map(|format| format.outfmt),
                ),
            "TBLASTN",
        ),
    );
    // LOSAT's refusal of two queries of one scope with one id comes before any report.
    if let Some(refusal) = batches.query_refusal.take() {
        return Err(refusal);
    }
    // NCBI crashes with IDENTITY's compressed lookup table once it searches a query batch: a run
    // that stops at its first batch gives NCBI's error.
    if batches.searched_done > 0 {
        super::args::check_identity_compressed_lookup_of(&resolved)?;
    }
    // The search strategy's requests of the batches (`-export_search_strategy`), with every
    // `-subject` record: written in the loop below, and for the batch that stops the run (every
    // query without letters) after the batches before it.
    let requests = args.strategy.as_ref().map(|_| {
        super::archive::requests(
            database.is_none().then_some(subject_records),
            resolved.lcase_masking,
        )
    });
    // The titles that the BLAST archive (outfmt 11) cannot write are refused after the batches
    // are searched (below): the batch that stops the run here, and one whose search fails, write
    // no archive.
    if batches.searched_done == 0 && batches.stop.is_some() {
        // A first batch whose set-up fails has its search strategy request, written before the
        // prolog (a refusal of its titles follows no byte).
        if let (Some(strategy), Some(requests), Some(span)) = (
            args.strategy.as_ref(),
            requests.as_ref(),
            batches.stop_batch.clone(),
        ) {
            requests.write(
                strategy,
                span.start,
                &searched_queries[span.clone()],
                query_input.placements.get(span.start),
            )?;
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // ```c
        // CObjMgr_QueryFactory::CObjMgr_QueryFactory(CBlastQueryVector & queries)
        //     : m_QueryVector(& queries)
        // {
        //     if (queries.Empty()) {
        //         NCBI_THROW(CBlastException, eInvalidArgument, "Empty CBlastQueryVector");
        //     }
        // ```
        // The first batch fails after the outfmt 0 prolog (`formatter.PrintProlog()`).
        write_first_batch_prologs(
            outputs,
            subject_path.map(|path| path.as_path()),
            subjects.len(),
            subject_letters,
            pairwise_database.as_ref(),
            args.tabular_taxdb.unwrap_or(false),
            "TBLASTN",
        )?;
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1407-1426
        // ```c
        //     db_length = BlastSeqSrcGetTotLen(seq_src);
        //
        //     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
        //
        //     /* iterate over all subject sequences */
        //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
        //            != BLAST_SEQSRC_EOF) {
        //        ...
        //        if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
        // ```
        // `-db`: the database's sequences, checked once for every batch (`XDbStream::new`)
        // where NCBI reads them first, in the preliminary search of the first batch with a
        // valid context (a batch without one is not searched, `CheckInternalData`). A database
        // without sequences (BLASTN's `Database memory map file error`,
        // blastn/blast_engine/run.rs) or a sequence file that cannot be read fails that batch,
        // after the reports of the batches before, which are searched without sequences; a run
        // without such a batch never reads the file. The
        // database is checked as NCBI's first read would fail, with the same error, and each
        // sequence is read where it is searched (`common/x_db_stream.rs`).
        let mut x_stream = None;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-980
        // ```c
        //    if (sbp->gbp) {
        //        min_subject_length = BlastSeqSrcGetMinSeqLen(seq_src);
        //        if (Blast_SubjectIsTranslated(program_number)) {
        //            min_subject_length/=3;
        //        }
        //    ...
        //    if(min_subject_length <=0) {
        // 	   return BLASTERR_SUBJECT_LENGTH_INVALID;
        //    }
        // ```
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:215-217
        // ```c
        //     case BLASTERR_SUBJECT_LENGTH_INVALID:
        //         new_msg->message = strdup("The average subject length is too short");
        // ```
        // With the Gumbel block a database's `GetMinLength` below 3 (an alias file's
        // `MIN_SEQ_LENGTH` 1 or 2, a list's shortest sequence) leaves no residue
        // (`min_subject_length` is NCBI's `Uint4`): the setup of the first batch with a valid
        // context fails after the database's first OID is found and before any sequence is read.
        // NCBI 2.17.0 (probe of line RG D1a, `db.ak.tbn.min1.6`): `BLAST engine error: The
        // average subject length is too short`, exit status 3.
        let (database_stores, mut stores_error) = match database {
            Some(database) => {
                let read = if database.num_oids() == 0 {
                    Err(PreliminaryReadError("Database memory map file error").into())
                } else if !subject_set.old_fsc && database.numbers().min_length as u32 / 3 == 0 {
                    Err(PreliminaryReadError("The average subject length is too short").into())
                } else {
                    database.x_stream(resolved.num_threads).map(|stream| {
                        x_stream = Some(stream);
                        Vec::new()
                    })
                };
                match read {
                    Ok(stores) => (stores, None),
                    Err(error) => (Vec::new(), Some(error)),
                }
            }
            None => (Vec::new(), None),
        };
        let search_letters = match (database, &x_stream) {
            (Some(_), Some(stream)) => SubjectLetters::Stream(stream),
            (Some(_), None) => SubjectLetters::Packed(&database_stores),
            (None, _) => subject_letters,
        };
        // The batches before the stop are searched and reported; a reader at the end of its
        // input reads no batch (`!input->End()`), and the reports have no query (an empty
        // pipe, or `-query -` after `-subject -`).
        let searched_done = batches.searched_done;
        let queries = &searched_queries[..searched_done];
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218 (the BLAST
        // XML report data of the first result NCBI formats: the matrix of its formatting table,
        // found by the name as typed)
        // ```c++
        //     else {
        //         string prog_name = Blast_ProgramNameFromType(
        //                                            m_Options.GetProgramType());
        //         if (prog_name != "blastn" && prog_name != "megablast") {
        //             NCBI_THROW(blast::CBlastException, eInvalidArgument,
        //                         "unsupported score matrix");
        // ```
        // NCBI fails when it formats the first batch in BLAST XML with a matrix name its XML data
        // does not know (`stage_e_report::render`, the same decision): the search strategy's file
        // holds the first batch's request alone, whatever the number of batches (probe 7 of RG-F1,
        // `-matrix blosum62` with outfmt 5 and 16).
        if let Some(strategy) = args.strategy.as_ref() {
            if crate::report::blast_xml::matrix::xml_formatting_fails(
                outputs.formats.iter().map(|format| format.outfmt),
                &resolved.matrix_name,
            ) {
                if let Some(end) = batches.first_request_end() {
                    strategy.stop_after(end);
                }
            }
        }
        // Each searched query's whole record and where its letters start in it (`-query_loc`).
        let whole_records: Vec<(&[u8], usize)> = (0..queries.len())
            .map(|index| {
                (
                    input_records[query_input.ordinal(index)].seq(),
                    query_input.placements.offset(index),
                )
            })
            .collect();
        let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
        let seg = resolved.seg.params();
        // NCBI c++/src/algo/blast/api/prelim_stage.cpp:172-188:
        // (*thread)->Run(); (*thread)->Join(&result);
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1411-1417:
        // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
        //                         CConstRef<blast::CBlastQueryVector> queries,
        // NCBI c++/src/app/blast/tblastn_app.cpp:251-252,275-301:
        // input.Reset(new CBlastInput(&*fasta, GetQueryBatchSize()));
        // for (; !input->End(); ...) {
        //     query = input->GetNextSeqBatch(*scope);
        //     CLocalBlast lcl_blast(query_factory, ..., db_adapter);
        //     results = lcl_blast.Run();
        // }
        // NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:70-127:
        // case eTblastn: retval = 20000;
        // Search each accumulated batch independently; preliminary linking
        // and sum-statistics depend on its concatenated query contexts.
        let mut results = Vec::with_capacity(queries.len());
        let mut lengths: Option<LocalSubjectParameters> = None;
        let mut ungapped_karlin = Vec::with_capacity(queries.len());
        let mut query_validity = Vec::with_capacity(queries.len());
        let mut query_batch_skipped = Vec::with_capacity(queries.len());
        // `-db`: the stored titles read so far, and the batch that a header file or the
        // sequence file failed with its error.
        let mut titles = super::database::DatabaseTitles::new();
        let mut failed: Option<(std::ops::Range<usize>, anyhow::Error)> = None;
        // Graduated 2026-10-22 (was LOSAT_X_TBNSSIDE, LOSAT_X_TBNBATCH): the batches of the
        // loop below searched ahead (`x_search_batches`), each with its range; `None` when there
        // is none. A batch whose range does not match is searched by the loop.
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:275-291
        // ```c
        //             for (; !input->End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
        // ...
        //                 query =  input->GetNextSeqBatch(*scope);
        //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
        // ...
        //                     CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
        //                     lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
        //                     results = lcl_blast.Run();
        // ```
        // Dispatch point: NCBI searches the query batches one after another (one CLocalBlast per
        // batch). x_search_batches searches them ahead (on one shared pool, or side by side) and
        // returns the reports in batch order; the loop below takes them from it.
        // A database whose sequences could not be read searches its batches one by one below.
        let x_batches = match stores_error {
            Some(_) => Vec::new(),
            None => query_input
                .batches(batch_size)
                .into_iter()
                .map(|batch| batch.searched)
                .filter(|searched| !searched.is_empty() && searched.end <= searched_done)
                .collect(),
        };
        // A read of the database that failed while the batches were searched ahead, or a batch
        // whose setup failed (`SetupEngineError`): they are searched again one at a time below, so
        // that the batch that fails does after the reports of the batches before.
        let x_searched_ahead = super::stage_d_pipeline::x_search_batches(
            &query_seqs,
            search_letters,
            x_batches,
            |range| LocalStageDProfile {
                seg: seg.as_ref(),
                soft_masking: resolved.soft_masking,
                mask_lowercase: resolved.lcase_masking,
                genetic_code: resolved.db_gencode,
                expect_value: resolved.evalue,
                max_target_seqs: resolved.hitlist_size,
                query_records: &whole_records[range],
                subject_set,
                engine_positives,
                culling_limit: resolved.culling_limit,
                pssm: pssm.as_deref(),
            },
            compo_adjust_mode,
            resolved.sum_stats,
            scoring,
            resolved.num_threads,
        );
        let mut x_searched = match x_searched_ahead {
            Ok(searched) => searched.map(Vec::into_iter),
            Err(error) if error.is::<XReadError>() || error.is::<SetupEngineError>() => None,
            Err(error) => return Err(error),
        };
        for range in query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty() && searched.end <= searched_done)
        {
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-280
            // ```c++
            //                 query =  input->GetNextSeqBatch(*scope);
            //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
            //
            //                 SaveSearchStrategy(args, m_CmdLineArgs, query_factory, m_OptsHndl);
            // ```
            // The batch's search strategy request (`-export_search_strategy`): once its query
            // factory is made, before any error of the batch is acted on (`x_search_batches` may
            // have searched it ahead; only the order of the file and of what the run writes
            // matters, and the file is apart from the reports). Its `EffectiveSearchSpace` option
            // is the one the batch before left (`note_searched` below).
            if let (Some(strategy), Some(requests)) = (args.strategy.as_ref(), requests.as_ref()) {
                requests.write(
                    strategy,
                    range.start,
                    &queries[range.clone()],
                    query_input.placements.get(range.start),
                )?;
            }
            let (mut batch_results, batch_lengths, batch_karlin, batch_validity) = {
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-291
                // ```c
                //                 query =  input->GetNextSeqBatch(*scope);
                //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
                // ...
                //                     CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
                //                     lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
                //                     results = lcl_blast.Run();
                // ```
                // Dispatch point: a report searched ahead for this very range is used; otherwise
                // run_local_for_report_threads searches the batch here, which is the port of
                // CLocalBlast::Run for one batch.
                let x_searched_here = match x_searched.as_mut().and_then(Iterator::next) {
                    Some((x_range, report)) if x_range == range => Ok(report),
                    _ => run_local_for_report_threads(
                        &query_seqs[range.clone()],
                        search_letters,
                        LocalStageDProfile {
                            seg: seg.as_ref(),
                            soft_masking: resolved.soft_masking,
                            mask_lowercase: resolved.lcase_masking,
                            genetic_code: resolved.db_gencode,
                            expect_value: resolved.evalue,
                            max_target_seqs: resolved.hitlist_size,
                            query_records: &whole_records[range.clone()],
                            subject_set,
                            engine_positives,
                            culling_limit: resolved.culling_limit,
                            pssm: pssm.as_deref(),
                        },
                        compo_adjust_mode,
                        resolved.sum_stats,
                        scoring,
                        resolved.num_threads,
                    ),
                };
                // A read of the database that failed in this
                // batch fails it as the check's error would (below); so does NCBI's engine error
                // of its setup.
                match x_searched_here {
                    Ok(report) => {
                        // The batch's set-up has written its `EffectiveSearchSpace` option, which
                        // the next batch's search strategy request holds.
                        if let Some(strategy) = args.strategy.as_ref() {
                            let (_, lengths, _, validity) = &report;
                            strategy.note_searched(
                                range.clone(),
                                super::archive::batch_searched_space(
                                    validity.iter().any(|&valid| valid),
                                    &lengths.lengths[0],
                                ),
                            )?;
                        }
                        report
                    }
                    Err(error) => match batch_failure(error) {
                        Ok(failure) => {
                            failed = Some((range, failure));
                            break;
                        }
                        Err(error) => return Err(error),
                    },
                }
            };
            // `-db`: the first batch with a valid context reads the sequences (above); a
            // reader error of NCBI's own fails it, LOSAT's refusal stops the run.
            if batch_validity.iter().any(|&valid| valid) {
                if let Some(error) = stores_error.take() {
                    match error.downcast::<PreliminaryReadError>() {
                        Ok(read_error) => {
                            failed = Some((range, app::engine_error(read_error.0)));
                            break;
                        }
                        Err(error) => return Err(error),
                    }
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:283-293
            // ```c++
            //     m_TbackSearch.Reset(new CBlastTracebackSearch(m_QueryFactory,
            //                                                   m_InternalData,
            //                                                   m_Opts,
            //                                                   seqinfo_src,
            //                                                   search_msgs));
            // ```
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1763-1777
            // ```c++
            //     if(hit_params->options->query_cov_hsp_perc > 0 || hit_params->options->max_hsps_per_subject > 0 ||
            //        (hit_params->options->hsp_filt_opt != NULL && hit_params->options->hsp_filt_opt->subject_besthit_opts != NULL)) {
            //     	s_FilterBlastResults(results, hit_params->options, query_info, program_number);
            //     }
            //
            //     /* Re-sort the hit lists according to their best e-values, because they
            //        could have changed. Only do this for a database search. */
            //     if (BlastSeqSrcGetTotLen(seq_src) > 0) {
            //         Blast_HSPResultsSortByEvalue(results);
            //     }
            // ```
            // `-max_hsps` and `-qcov_hsp_perc` filter the batch's final hit lists (each query's
            // context is its searched letters), which the search has sorted by e-value; the
            // sort is a total order (best e-value, first score, OID), so sorting again after
            // the filter gives NCBI's order. A search with hits has subject letters.
            // `-culling_limit` runs NCBI's culling pipe first on each query's hit list
            // (`BlastHSPStreamTBackClose`, before the filters of `s_FilterBlastResults`;
            // hspfilter_culling.c:699-733); its lists leave in the order of the tree, and the
            // sort by e-value below is the one NCBI runs unconditionally after the pipes.
            if resolved.max_hsps > 0
                || resolved.query_cov_hsp_perc > 0.0
                || resolved.culling_limit > 0
            {
                for (hitlist, query) in batch_results.iter_mut().zip(&query_seqs[range.clone()]) {
                    if resolved.culling_limit > 0 {
                        hitlist.cull(resolved.culling_limit, query.len())?;
                    }
                    if resolved.max_hsps > 0 || resolved.query_cov_hsp_perc > 0.0 {
                        hitlist.filter_results(
                            resolved.max_hsps,
                            resolved.query_cov_hsp_perc,
                            query.len(),
                        );
                    }
                    hitlist.sort_by_evalue_for_positive_totlen();
                }
            }
            // `-db`: the definition lines of the batch's final hit lists are read where NCBI
            // makes the batch's Seq-aligns (`database::read_titles`); a header file that NCBI
            // cannot read fails this batch, after the reports of the batches before.
            if let Some(database) = database {
                if let Err(error) =
                    super::database::read_titles(database, &batch_results, range.start, &mut titles)
                {
                    match error.downcast::<HeaderReadError>() {
                        Ok(failure) => {
                            failed = Some((range, failure.error));
                            break;
                        }
                        Err(error) => return Err(error),
                    }
                }
                if let Some(scopes) = batches.query_scopes.as_ref() {
                    crate::algorithm::common::database::refuse_reported_ids_of_queries(
                        &titles,
                        batch_results
                            .iter()
                            .enumerate()
                            .flat_map(|(index, hitlist)| {
                                hitlist.lists().iter().filter_map(move |list| {
                                    Some((range.start + index, u32::try_from(list.oid).ok()?))
                                })
                            }),
                        scopes,
                        "TBLASTN",
                    )?;
                }
            }
            results.append(&mut batch_results);
            if let Some(all_lengths) = &mut lengths {
                // NCBI tblastn_app.cpp:288-301 formats each batch's own
                // query-context lengths; retain those values in input order.
                all_lengths.lengths.extend(batch_lengths.lengths);
                all_lengths.cutoffs.extend(batch_lengths.cutoffs);
            } else {
                lengths = Some(batch_lengths);
            }
            ungapped_karlin.extend(batch_karlin);
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
            // if (CheckInternalData() != 0) each result in this Run() batch
            // receives an unsearched ancillary state and null align set.
            let skipped = batch_validity.iter().all(|&valid| !valid);
            query_batch_skipped.extend(std::iter::repeat_n(skipped, batch_validity.len()));
            query_validity.extend(batch_validity);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2777-2792
        // ```c
        //       sbp->kbp_std[context] = kbp = Blast_KarlinBlkNew();
        //       loop_status = Blast_KarlinBlkUngappedCalc(kbp, sbp->sfp[context]);
        //       if (loop_status) {
        //           contexts[context].is_valid = FALSE;
        // ```
        // A PSSM query whose letters give no Karlin-Altschul parameters (one that `-seg` masks
        // completely) makes NCBI BLAST+ 2.17.0's tblastn crash (SIGSEGV, probe `gs.masked`): LOSAT
        // refuses it (DECISIONS 35 Q5 of line PS, as psiblast's, DECISIONS 12 S4), before
        // anything is written.
        if pssm.is_some() && query_validity.iter().any(|&valid| !valid) {
            anyhow::bail!(
                "a PSSM query whose letters give no Karlin-Altschul parameters (one that -seg masks completely), with which NCBI BLAST+ 2.17.0's tblastn crashes, is not supported by LOSAT's TBLASTN"
            );
        }
        // A run without a searched query has no parameters of a query set.
        let mut lengths = lengths.unwrap_or(LocalSubjectParameters {
            lengths: Vec::new(),
            cutoffs: Vec::new(),
            link: None,
            hit_cutoff_min: 0,
            word_cutoff_min: 0,
            prelim_evalue: 0.0,
        });
        // A batch that a database file failed is not reported, nor are the batches after it.
        let queries = match &failed {
            Some((range, _)) => {
                results.truncate(range.start);
                lengths.lengths.truncate(range.start);
                lengths.cutoffs.truncate(range.start);
                ungapped_karlin.truncate(range.start);
                query_validity.truncate(range.start);
                query_batch_skipped.truncate(range.start);
                &queries[..range.start]
            }
            None => queries,
        };
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:277-294
        // ```c++
        //                 query =  input->GetNextSeqBatch(*scope);
        //                 query_factory.Reset(new CObjMgr_QueryFactory(*query));
        //     ...
        //                     results = lcl_blast.Run();
        //     ...
        //                 if (fmt_args->ArchiveFormatRequested(args)) {
        //                     formatter.WriteArchive(*query_factory, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
        // ```
        // The archive of each batch that NCBI completes holds its queries and, with `-subject`,
        // every subject. The batch that stops the run writes neither: a reader error, `Empty
        // CBlastQueryVector` or a batch of queries without letters (`read_query_batches`), and a
        // batch whose search fails (`failed`: the sequence file or a header file that NCBI cannot
        // read); the error archive holds no sequence. So a title the archive cannot write is
        // refused here, where the batches have been searched and their definition lines read, and
        // before the first byte of any report.
        // (The search strategy's requests check their titles where they are written, in the
        // loop above.)
        if args.archive.is_some() && !queries.is_empty() {
            super::archive::check_titles(queries, "query")?;
            super::archive::check_titles(subject_records, "subject")?;
        }
        // The request of the batch that stops the run (every query without letters) comes after
        // these reports (below): its titles are refused before them, when the reports end with
        // that batch (no batch failed, and NCBI's formats do not stop first).
        if let (Some(strategy), Some(requests), Some(span), None) = (
            args.strategy.as_ref(),
            requests.as_ref(),
            batches.stop_batch.clone(),
            &failed,
        ) {
            if strategy.takes(span.clone()) {
                requests.check_titles(&searched_queries[span])?;
            }
        }
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1443-1451:
        // if (results.HasWarnings()) ERR_POST(Warning << results.GetWarningStrings());
        // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
        // ```c
        //         } catch (const CException& e) {
        //             ...
        //             CRef<CSearchMessage> m
        //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
        //             messages[index].push_back(m);
        //             s_InvalidateQueryContexts(qinfo, index);
        //         }
        // ```
        // The warnings of each query are written before its report, once (`QueryWarnings`):
        // the reader's messages of the query's batch (before its first query,
        // `read_query_batches`), then the query's own, on one line in NCBI's order
        // (`query_warning`): the O characters read as X, and the message of its set-up (a
        // query without letters) or of its Karlin-Altschul parameters (an invalid query
        // with letters).
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1452
        // ```c++
        //     if (m_FormatType == CFormattingArgs::eAsnText
        //       ...
        //       || m_FormatType == CFormattingArgs::eJsonSeqalign
        //       || m_FormatType == CFormattingArgs::eSAM)
        //     {
        //         x_PrintStructuredReport(results, queries);
        //         return;
        //     }
        //     ...
        //     if (results.HasWarnings()) {
        //         ERR_POST(Warning << results.GetWarningStrings());
        //     }
        // ```
        // A structured format (outfmt 8, 12) writes none of a query's own warnings
        // (`ReportFormat::writes_query_warnings`); the warnings go with the first format.
        let result_warnings = outputs.formats.first().is_none_or(|format| {
            crate::blastinput::app::ReportFormat::from_number(super::stage_e_report::format_number(
                format.outfmt,
            ))
            .is_none_or(crate::blastinput::app::ReportFormat::writes_query_warnings)
        });
        let filtered_empty = database.is_some_and(|database| database.filtered_empty());
        // The texts of the subjects that the reports of the queries whose batch's prefetch is
        // off read one by one, by query, where the first format places them (`places_reads`).
        let mut subject_reads = crate::report::query_warnings::SubjectReads::default();
        let mut query_warning_lines: Vec<Vec<u8>> = queries
            .iter()
            .zip(&query_validity)
            .enumerate()
            .map(|(index, (query, valid))| {
                let mut lines = batches.warnings[index].clone();
                if !result_warnings {
                    return lines;
                }
                let messages = crate::report::query_warnings::query_messages(
                    query.seq(),
                    *valid,
                    index < identity_warned,
                    filtered_empty && !query_batch_skipped[index],
                );
                lines.extend(match pssm.as_deref() {
                    // A PSSM query's label is its Seq-id and the title of NCBI's query source of
                    // the PSSM's Bioseq (`PsiSearch::warning_title`), as psiblast's
                    // (`blastp::blast_engine::pssm_hooks::pssm_query_warning`).
                    Some(psi) => crate::report::query_warnings::query_warning_with_seq_id(
                        "tblastn",
                        psi.query.warning_id(),
                        psi.warning_title.as_deref().unwrap_or(&psi.query.title),
                        &messages,
                    ),
                    None => crate::report::query_warnings::query_warning(
                        "tblastn",
                        args.first_query_ordinal + query_input.ordinal(index),
                        query,
                        &messages,
                    ),
                });
                lines
            })
            .collect();
        // The query batches that were read and searched, by their searched queries (the batches
        // before a database failure or a reader stop): each is reported.
        let batch_spans: Vec<std::ops::Range<usize>> = query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty() && searched.end <= queries.len())
            .collect();
        // `-db`: the reported sequences by rank, with their records and lengths, and their packed
        // subjects where a report reads their letters; the database's other sequences are dropped
        // before the reports.
        let reported = match database {
            Some(database) => Some(
                super::database::reported_subjects(
                    database,
                    &titles,
                    &mut results,
                    search_letters,
                    super::stage_e_report::reads_subject_letters(
                        outputs.hits.is_some(),
                        outputs.formats.iter().map(|format| format.outfmt),
                    ),
                    parsed_shown(&outputs.formats),
                    &batch_spans,
                )
                // A read of a reported sequence that
                // failed: the run's error, as the check's error would be.
                .map_err(|error| batch_failure(error).unwrap_or_else(|error| error))?,
            ),
            None => None,
        };
        // The library's records name a database's sequences through this event.
        if let (Some(reported), Some(database), true, Some(observer)) = (
            &reported,
            database,
            outputs.hits.is_some(),
            outputs.observer.as_deref_mut(),
        ) {
            observer.database_subjects(reported.library_subjects(database));
        }
        drop(database_stores);
        // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:298-302
        // ```c++
        //                     BlastFormatter_PreFetchSequenceData(*results, scope,
        //                     		                            fmt_args->GetFormattedOutputChoice());
        //                     ITERATE(CSearchResultSet, result, *results) {
        //                         formatter.PrintOneResultSet(**result, query);
        //                     }
        // ```
        // Before the reports of a batch NCBI fetches every subject that the batch reports, and
        // with `taxdb.bti` on the search path the taxonomy database writes its messages there
        // (`TblastnDatabase::prefetch_messages`: the texts of opening it at the run's first
        // lookup, a volume's lookup of a taxid once per run): after the reader's messages of the
        // batch and before the warnings of its first query. They go with the warnings, where the
        // first format is outfmt 0, 5, 6 or 7, or with `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 or
        // 12, which write the reader's messages only (`database::prefetches_subjects`); XML2 and JSON
        // write theirs at each `Hit` (`TaxdbPrinter`). They are read here, after the refusals of
        // the reported subjects and before any format is written, so a refusal (a taxonomy
        // database or a taxid that LOSAT cannot read as NCBI does) comes before any byte, as
        // those refusals do; without `taxdb.bti` nothing is read. A failed batch is not reported
        // and prefetches nothing (`failed`, a header file or the sequence file that NCBI cannot
        // read during the batch's search).
        if let (Some(reported), Some(database)) = (&reported, database) {
            // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:716-725
            // ```c++
            // s_IsPrefetchFormat(blast::CFormattingArgs::EOutputFormat format_type)
            // {
            // 	if ((format_type == CFormattingArgs::eAsnText) ||
            // 	    (format_type == CFormattingArgs::eAsnBinary) ||
            // 	    (format_type == CFormattingArgs::eArchiveFormat)||
            // 	    (format_type == CFormattingArgs::eJsonSeqalign)) {
            // 		return false;
            // 	}
            // 	return true;
            // ```
            // Outfmt 10 and 20 prefetch as the tabular 6 and 7 do (a name of the taxonomy
            // database that the prefetch would look up is LOSAT's refusal, `tabular`).
            // Where the prefetch fails at a subject with parsed Seq-ids, NCBI writes its warning
            // after those messages (`TblastnDatabase::batch_prefetch`), with or without
            // `taxdb.bti`.
            let taxdb_found = database.taxdb_found();
            if outputs.formats.first().is_some_and(|format| {
                crate::algorithm::common::database::prefetches_subjects(Some(
                    match super::stage_e_report::format_number(format.outfmt) {
                        10 | 20 => 6,
                        number => number,
                    },
                ))
            }) && (taxdb_found || database.parsed_seqids())
            {
                // A custom tabular report with a plural name looks every taxid of its rows up,
                // taxid 0 too, which the prefetch does not: the rows of a batch with hits may open
                // the taxonomy database before a later batch's prefetch
                // (`TaxdbRun::open_for_rows`).
                let first_spec = outputs
                    .formats
                    .first()
                    .and_then(|format| super::tabular::tabular_spec(format.outfmt));
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
                // ```c++
                //        		const CBioseq_Handle& subject_bh =
                //                 scope.GetBioseqHandle(align.GetSeq_id(1));
                // ```
                // Whether the first format's writer writes the texts of the subjects that a
                // query's report reads one by one where it reads each
                // (`QueryWarnings::subject_read`): outfmt 0 (description lines and alignment
                // headings), 5 (`Hit`s) and the tabular rows.
                let places_reads = outputs.formats.first().is_some_and(|format| {
                    matches!(
                        super::stage_e_report::format_number(format.outfmt),
                        0 | 5 | 6 | 7 | 10 | 20
                    )
                });
                let plural_names =
                    taxdb_found && first_spec.as_ref().is_some_and(|spec| spec.checks_taxdb());
                let rows_open = |hitlists: &[super::stage_d_results::KappaResultHitList]| {
                    if plural_names
                        && hitlists
                            .iter()
                            .flat_map(super::stage_d_results::KappaResultHitList::lists)
                            .any(|list| !list.hsps.hsps.is_empty())
                    {
                        database.open_for_tabular_rows()?;
                    }
                    Ok::<(), anyhow::Error>(())
                };
                for span in &batch_spans {
                    let hitlists = &results[span.clone()];
                    // The batch's Seq-aligns, one for each HSP (`database::seq_aligns_of`).
                    let seq_aligns = super::database::seq_aligns_of(hitlists);
                    // Where `PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch off
                    // (`shared::prefetch_off`), each query's report fetches the subjects it
                    // shows, after the query's own warnings, where the first format's rows read
                    // their Bioseqs (`TabularSpec::subject_fetch`; for `nident` alone the rows
                    // with identities, `database::reads_subject`). NCBI looks the taxonomy up
                    // only where it fetches a subject: in such a report, or in the prefetch of a
                    // batch with Seq-aligns (one without them loads no id, below). Every batch's
                    // lookups are made here before any report is written, so the run's refusal
                    // of NCBI's taxonomy cache (`shared::TaxdbRun::cache_collision`) is asked
                    // only there: a later batch that fetches still refuses before any byte.
                    if crate::algorithm::common::database::prefetch_off(seq_aligns) {
                        // No prefetch; without `taxdb.bti` the reports' fetches look nothing up.
                        if !taxdb_found {
                            continue;
                        }
                        let rows = first_spec.as_ref().map_or(SubjectFetch::Every, |spec| {
                            spec.subject_fetch(crate::blastinput::app::Program::Tblastn)
                        });
                        for (query, hitlist) in span.clone().zip(hitlists) {
                            let shown = std::slice::from_ref(hitlist);
                            if hitlist
                                .lists()
                                .iter()
                                .any(|list| super::database::reads_subject(list, rows))
                            {
                                // The first format writes the texts where it reads each subject
                                // (`places_subject_reads`), else before the report.
                                if places_reads {
                                    subject_reads.set(
                                        query,
                                        database.subject_read_texts(
                                            reported,
                                            shown,
                                            rows,
                                            plural_names,
                                        )?,
                                    );
                                } else {
                                    let messages =
                                        database.prefetch_messages(reported, shown, rows)?;
                                    query_warning_lines[query].extend(messages);
                                }
                            }
                            rows_open(shown)?;
                        }
                        continue;
                    }
                    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-772
                    // ```c++
                    //        CScope::TIds ids;
                    //        vector<TSeqRange> ranges;
                    //        s_ExtractSeqidsAndRanges(results, ids, ranges);
                    //        _TRACE("Prefetching " << ids.size() << " sequence lengths");
                    //        LoadSequencesToScope(ids, ranges, scope);
                    // ```
                    // A batch without Seq-aligns prefetches no sequence.
                    if seq_aligns == 0 {
                        continue;
                    }
                    let messages = database.batch_prefetch(reported, hitlists)?;
                    rows_open(hitlists)?;
                    if messages.is_empty() {
                        continue;
                    }
                    let reader = batches.warnings[span.start].len();
                    query_warning_lines[span.start].splice(reader..reader, messages);
                }
            }
        }
        // `-db`, XML2 and JSON: the taxid of each reported sequence, and what the taxonomy
        // database says about it (the refusals of both come before any byte).
        let xml_taxa = match (&reported, database) {
            (Some(reported), Some(database))
                if outputs
                    .formats
                    .iter()
                    .any(|format| super::stage_e_report::xml2_kind(format.outfmt).is_some()) =>
            {
                database.read_taxids(reported, &results)?
            }
            _ => Default::default(),
        };
        // The custom tabular reports (`tabular`): whether `taxdb.bti` was found, the taxids of the
        // reported database sequences where a field shows them and the taxonomy database's names
        // that the rows look up (read here, before any byte, so that their refusals come before
        // it), and the `-query_loc` range's length (`qlen`, `qcovs`).
        let tabular_specs: Vec<_> = outputs
            .formats
            .iter()
            .filter_map(|format| super::tabular::tabular_spec(format.outfmt))
            .collect();
        let taxonomy = match (&reported, database) {
            (Some(reported), Some(database))
                if tabular_specs.iter().any(|spec| spec.needs().taxonomy) =>
            {
                database.tabular_taxonomy(reported, &results)?
            }
            _ => Vec::new(),
        };
        let names = match (tabular_specs.first(), database) {
            (Some(spec), Some(database)) => database.tabular_names(spec, &taxonomy)?,
            _ => None,
        };
        let tabular_inputs = super::tabular::TabularInputs {
            taxdb_found: args.tabular_taxdb.unwrap_or(false),
            taxonomy,
            names,
            range_length: resolved
                .query_range
                .as_ref()
                .map(crate::blastinput::seq_range::SequenceRange::length),
        };
        // The records of a subject set whose ids `-parse_deflines` read, by subject index.
        let records_parsed = if args.parse_deflines && reported.is_none() {
            crate::report::db_subject::parsed_records(
                subjects
                    .iter()
                    .map(|record| (record.believed(), record.title.as_slice())),
                args.show_gis,
                false,
                crate::blastinput::app::parse_formatting_string(&args.outfmt)?.number == 0,
                "TBLASTN",
            )?
        } else {
            Vec::new()
        };
        let (report_records, report_letters, shown) = match (&reported, database) {
            (Some(reported), Some(database)) => (
                reported.records.as_slice(),
                reported.letters(),
                ShownSubjects::Database {
                    name: database.name(),
                    pairwise: pairwise_database
                        .as_ref()
                        .expect("a database search has its Database: lines"),
                    reported,
                },
            ),
            _ => (
                subjects,
                subject_letters,
                ShownSubjects::Fasta {
                    path: subject_path.expect("a search of -subject has its path"),
                    parsed: &records_parsed,
                },
            ),
        };
        // The BLAST archive (outfmt 11): one document per query batch that was read and searched
        // (`batch_spans`).
        let archive_render = args
            .archive
            .as_ref()
            .map(|handle| super::archive::ArchiveRender {
                handle,
                spans: &batch_spans,
            });
        // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); it is
        // written after all batches are searched. No NCBI counterpart.
        let t_output = crate::utils::xtiming::clock();
        render(
            outputs,
            queries,
            report_records,
            report_letters,
            &shown,
            subject_set,
            &mut results,
            &lengths,
            &ungapped_karlin,
            &query_validity,
            &query_batch_skipped,
            scoring,
            &resolved.matrix_name,
            crate::blastinput::report_counts::ReportCounts {
                hitlist_size: resolved.hitlist_size,
                num_descriptions: resolved.num_descriptions,
                num_alignments: resolved.num_alignments,
            },
            resolved.formatting,
            resolved.db_gencode,
            seg.as_ref(),
            resolved.lcase_masking,
            &query_warning_lines,
            subject_reads.queries(),
            &super::stage_e_report::ReportRanges {
                queries: &query_input.placements,
                subjects: subject_placements,
                // NCBI reference: c++/src/app/blast/tblastn_app.cpp:342
                // ```c++
                //         formatter.PrintEpilog(opt);
                // ```
                // A run that stops at a query batch writes no epilog.
                epilog: batches.stop.is_none() && failed.is_none(),
            },
            super::blast_xml::XmlValues {
                expect: resolved.evalue,
                soft_masking: resolved.soft_masking,
                gapped: resolved.gapped,
                // IDENTITY's reset is NCBI's option of the search (`search_settings`).
                cbs: if identity_cbs_reset {
                    0
                } else {
                    resolved.composition_mode as i32
                },
                db_gencode: i32::from(resolved.db_gencode),
                identity_warned,
                filtered_empty,
            },
            xml_taxa,
            args.out.as_deref(),
            archive_render.as_ref(),
            &tabular_inputs,
            pssm.as_deref(),
        )?;
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
        // The failed batch's reading messages (its search strategy request was written before its
        // search), then its error.
        if let Some((range, error)) = failed {
            for warnings in &batches.warnings[range] {
                outputs.diagnostics.write_all(warnings)?;
            }
            return Err(error);
        }
    }
    // A batch whose set-up fails (every query without letters) has its search strategy request:
    // NCBI made its query factory before the search (`QueryBatches::stop_batch`), after the
    // batches above (a first batch's was written before the prolog).
    if let (Some(strategy), Some(requests), Some(span)) = (
        args.strategy.as_ref(),
        requests.as_ref(),
        batches.stop_batch.clone(),
    ) {
        requests.write(
            strategy,
            span.start,
            &searched_queries[span.clone()],
            query_input.placements.get(span.start),
        )?;
    }
    if let Some((warnings, error)) = batches.stop {
        outputs.diagnostics.write_all(&warnings)?;
        return Err(error);
    }
    Ok(())
}

/// What the reports `formats` show of the reported database sequences, for LOSAT's refusals of
/// those with parsed Seq-ids (`common::database::ParsedShown`): the tabular formats (6, 7, 10,
/// 20) whose fields read the subject (`report::tabular::reads_subject`), outfmt 0, the
/// structured documents (5, 13 to 16), the formats that join the alignments of following hits
/// whose Seq-align ids match (0, 5, 13 to 16), those that write the sequences' taxids (13 to 16,
/// the tabular reports of taxonomy fields), and whether the first format's batches prefetch the
/// subjects.
fn parsed_shown(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
) -> crate::algorithm::common::database::ParsedShown {
    let numbers: Vec<i32> = formats
        .iter()
        .map(|format| super::stage_e_report::format_number(format.outfmt))
        .collect();
    let pairwise = numbers.contains(&0);
    crate::algorithm::common::database::ParsedShown {
        tabular: formats.iter().any(|format| {
            matches!(
                super::stage_e_report::format_number(format.outfmt),
                6 | 7 | 10 | 20
            ) && crate::report::tabular::reads_subject(
                super::tabular::tabular_spec(format.outfmt).as_ref(),
                crate::blastinput::app::Program::Tblastn,
            )
        }),
        pairwise,
        structured: formats
            .iter()
            .any(|format| super::stage_e_report::is_structured_xml(format.outfmt)),
        joined: pairwise
            || formats.iter().any(|format| {
                super::stage_e_report::is_structured_xml(format.outfmt)
                    || super::tabular::tabular_spec(format.outfmt)
                        .is_some_and(|spec| spec.needs().joins_subjects())
            }),
        sam: false,
        taxonomy: formats.iter().any(|format| {
            super::stage_e_report::xml2_kind(format.outfmt).is_some()
                || super::tabular::tabular_spec(format.outfmt)
                    .is_some_and(|spec| spec.needs().taxonomy)
        }),
        words: formats
            .iter()
            .filter(|format| {
                matches!(
                    super::stage_e_report::format_number(format.outfmt),
                    6 | 7 | 10 | 20
                )
            })
            .fold(
                crate::report::tabular::PrintedIds::default(),
                |words, format| {
                    words.or(crate::report::tabular::PrintedIds::of(
                        super::tabular::tabular_spec(format.outfmt).as_ref(),
                    ))
                },
            ),
        // NCBI's `s_IsPrefetchFormat` prefetches for every tabular format (quoted in `run`).
        prefetched: crate::algorithm::common::database::prefetches_subjects(numbers.first().map(
            |&number| match number {
                10 | 20 => 6,
                number => number,
            },
        )) && !crate::algorithm::common::database::prefetch_always_off(),
        one_by_one: crate::algorithm::common::database::prefetch_may_be_off(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A read that failed during the search fails its batch with
    // NCBI's engine error (exit status 3) when the reader's error is NCBI's own, and stops the run
    // with LOSAT's refusal or any other error, as the check's error does.
    #[test]
    fn a_read_failed_in_the_search_fails_as_the_check_s_error() {
        let engine = batch_failure(
            XReadError(PreliminaryReadError("Database memory map file error").into()).into(),
        )
        .expect("the engine error");
        let native = engine
            .downcast::<crate::cli::NativeError>()
            .expect("a NativeError");
        assert_eq!(native.exit, 3);
        assert_eq!(
            native.message,
            "BLAST engine error: Database memory map file error\n"
        );
        let refusal = batch_failure(XReadError(anyhow::anyhow!("refused")).into()).unwrap_err();
        assert_eq!(refusal.to_string(), "refused");
        let other = batch_failure(anyhow::anyhow!("other")).unwrap_err();
        assert_eq!(other.to_string(), "other");
    }

    // NCBI 2.17.0 tblastn with `-searchsp` or `OLD_FSC` and `-subject` of three empty records
    // (env parity probe f_probe of 2026-11-12): the batch fails with NCBI's engine error.
    #[test]
    fn a_short_average_subject_length_fails_the_batch_with_the_engine_error() {
        let native = batch_failure(SetupEngineError::SUBJECT_LENGTH_INVALID.into())
            .expect("the engine error")
            .downcast::<crate::cli::NativeError>()
            .expect("a NativeError");
        assert_eq!(native.exit, 3);
        assert_eq!(
            native.message,
            "BLAST engine error: The average subject length is too short\n"
        );
    }
}
